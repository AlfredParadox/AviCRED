# ============================================================================================
# PROJECT: AVICRED
# OBJECTIVE: OBSERVER RELIABILITY AND CHECKLIST-DATA QUALITY ASSESSMENT OF eBIRD DATA
#            VIA MULTI-METHOD FRAMEWORK
#
# RUN THE R SCRIPT TO FILTER AND PROCESS EBIRD DATA AND ASSESS CHECKLIST DATA-QUALITY WITH
# OBSERVER CREDIBILITY BY COMPUTING AVICRED SCORE. SUBSEQUENTLY, VALIDATE RESULTS WITH
# MACHINE-LEARNING ALGORITHMS TO UNDERSAND OBSERVER AND CHECKLIST BEHAVIOURAL PROFILES
# ============================================================================================

############# CALL LIBRARIES ##############   

##### INSTALL LIBRARIES FIRST FROM install_dependencies.R file in the folder

library(dplyr)
library(Hmisc)
library(data.table)
library(DBI)
library(duckdb)
library(ggplot2)
library(tidyr)
library(scales)
library(extrafont)
library(isotree)
library(ComplexUpset)
library(ggalluvial)
library(sf)
library(hexbin)
library(treemapify)
library(terra)
library(patchwork)
library(ggridges)
library(caret)
library(yardstick)
library(viridis)
library(purrr)

loadfonts(device = "win")
windowsFonts(
  Times = windowsFont("Times New Roman")
)
setwd("file path")

ebd_file <- "file path/ebd raw.txt"

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = "file path/avicred.duckdb",
  read_only = FALSE
)

################ YEAR-WISE CHECKLISTS GRAPH ############

dbExecute(con, "
CREATE OR REPLACE VIEW ebd_file AS
SELECT *
FROM read_csv_auto(
'file path/ebd raw.txt',
delim = '\t',
header = TRUE,
ignore_errors = TRUE
)
")

yearwise_checklists <- dbGetQuery(con, "
SELECT
YEAR(\"OBSERVATION DATE\") AS Year,
COUNT(DISTINCT \"SAMPLING EVENT IDENTIFIER\") AS Checklists
FROM ebd_file
WHERE
\"ALL SPECIES REPORTED\" = 1
AND COUNTRY = 'India'
GROUP BY Year
ORDER BY Year
")

ggplot(
  yearwise_checklists,
  aes(
    x = Year,
    y = Checklists,
    group = 1
  )
) +
  geom_line(
    colour = "#006D6F",
    linewidth = 1
  ) +
  geom_point(
    colour = "#006D6F",
    fill = "#B2DFDB",
    shape = 21,
    size = 2.8,
    stroke = 0.8
  ) +
  geom_text(
    aes(label = comma(Checklists)),
    vjust = -0.7,
    family = "Times New Roman",
    size = 3.2
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Year-wise Complete eBird Checklists in India (2000–2025)",
    x = "Year",
    y = "Number of Complete Checklists"
  ) +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.grid.minor = element_blank()
  )

############### SPLIT INTO 2010-25 ################

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_2010_25 AS
SELECT *
FROM ebd_file
WHERE YEAR(\"OBSERVATION DATE\") BETWEEN 2010 AND 2025
")

################# EFFORT FILTERING AS PER EBIRD RECOMMENDATION ###########

dbExecute(con, "
CREATE OR REPLACE VIEW ebd_2010_25_view AS
SELECT *
FROM ebd_2010_25
")

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_filtered AS
SELECT *
FROM ebd_2010_25
WHERE
  \"ALL SPECIES REPORTED\" = 1 
  AND \"PROTOCOL NAME\" IN ('Traveling', 'Stationary')
  AND \"DURATION MINUTES\" BETWEEN 10 AND 300
  AND \"NUMBER OBSERVERS\" <= 10
  AND (
        \"PROTOCOL NAME\" = 'Stationary'
        OR COALESCE(\"EFFORT DISTANCE KM\",0) <= 8
      )
  AND APPROVED = 1
")

dbGetQuery(con, "
SELECT COUNT(*) AS records_after_filter
FROM ebd_filtered
")

dbGetQuery(con, "
SELECT COUNT(DISTINCT \"SAMPLING EVENT IDENTIFIER\") AS checklists_after_filter
FROM ebd_filtered
")

dbExecute(con, "
COPY ebd_filtered
TO 'D:/ebd_filtered.txt'
(
  FORMAT CSV,
  DELIMITER '\t',
  HEADER TRUE
)
")

############ LOAD FILE IF SAVED PREVIOUSLY ############

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_filtered AS
SELECT *
FROM read_csv_auto(
  'FILE PATH/ebd_filtered.csv',
  header = TRUE
)
")

dbGetQuery(con, '
SELECT
    COUNT(DISTINCT "SAMPLING EVENT IDENTIFIER") AS n_checklists,
    COUNT(DISTINCT "OBSERVER ID") AS n_observers
FROM ebd_filtered
')

############## MASTER CHECKLIST ##############

dbExecute(con, "
CREATE OR REPLACE TABLE checklist_master AS
SELECT
  -- identifiers
  \"SAMPLING EVENT IDENTIFIER\" AS checklist_id,
  MAX(\"OBSERVER ID\") AS observer_id,
  -- date
  MAX(\"OBSERVATION DATE\") AS obs_date,
  -- location
  MAX(\"LATITUDE\") AS latitude,
  MAX(\"LONGITUDE\") AS longitude,
  MAX(STATE) AS state,
  MAX(COUNTY) AS county,
  MAX(\"LOCALITY ID\") AS locality_id,
  MAX(LOCALITY) AS locality,
  -- effort
  MAX(\"PROTOCOL NAME\") AS protocol,
  COALESCE(MAX(\"EFFORT DISTANCE KM\"), 0) AS distance_km,
  MAX(\"DURATION MINUTES\") AS duration_min,
  MAX(\"NUMBER OBSERVERS\") AS n_observers,
  -- checklist completeness
  MAX(\"ALL SPECIES REPORTED\") AS complete_flag,
  -- richness
  COUNT(DISTINCT \"COMMON NAME\") AS richness,
  -- abundance
  SUM(
    CASE
      WHEN regexp_matches(\"OBSERVATION COUNT\", '^[0-9]+$')
      THEN CAST(\"OBSERVATION COUNT\" AS INTEGER)
      ELSE 0
    END
  ) AS abundance_per_checklist
FROM ebd_filtered
GROUP BY \"SAMPLING EVENT IDENTIFIER\"
")

dbExecute(con, "
COPY checklist_master
TO 'file path/checklist_master.txt'
(
  FORMAT CSV,
  DELIMITER '\t',
  HEADER TRUE
)
")

dbExecute(con, "
COPY checklist_master
TO 'file path/checklist_master.csv'
(
  FORMAT CSV,
  DELIMITER ',',
  HEADER TRUE
)
")

############ LOAD FILE IF SAVED PREVIOUSLY ############

dbExecute(con, "
CREATE OR REPLACE TABLE checklist_master AS
SELECT *
FROM read_csv_auto(
  'file path/checklist_master.txt',
  delim = '\t',
  header = TRUE
)
")

############## MASTER CHECKLIST STATS ##############

dbGetQuery(con, "
SELECT COUNT(DISTINCT observer_id) AS n_observers
FROM checklist_master
")

dbGetQuery(con, "
SELECT
  MIN(obs_date) AS first_date,
  MAX(obs_date) AS last_date
FROM checklist_master
")

dbGetQuery(con, "
SELECT
  MIN(richness) AS min_rich,
  AVG(richness) AS mean_rich,
  MAX(richness) AS max_rich
FROM checklist_master
")

dbGetQuery(con, "
SELECT
  MIN(abundance_per_checklist) AS min_abund,
  AVG(abundance_per_checklist) AS mean_abund,
  MAX(abundance_per_checklist) AS max_abund
FROM checklist_master
")

#================ CALCULATE ORI METRICS ================#

dbExecute(con, "
CREATE OR REPLACE TABLE observer_metrics AS
SELECT
    observer_id,
    COUNT(*) AS n_checklists,
    COUNT(DISTINCT obs_date) AS n_birding_days,
    DATE_DIFF(
        'year',
        MIN(obs_date),
        MAX(obs_date)
    ) + 1 AS years_active
FROM checklist_master
GROUP BY observer_id
")

############# CALCULATE SPECIES OBSERVED PER OBSERVER #########

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_filtered_species AS
SELECT *
FROM ebd_filtered
WHERE
  \"ALL SPECIES REPORTED\" = 1
  AND \"PROTOCOL NAME\" IN ('Traveling', 'Stationary')
  AND \"DURATION MINUTES\" BETWEEN 10 AND 300
  AND \"NUMBER OBSERVERS\" <= 10
  AND (
        \"PROTOCOL NAME\" = 'Stationary'
        OR COALESCE(\"EFFORT DISTANCE KM\",0) <= 8
      )
  AND APPROVED = 1
  AND \"OBSERVATION DATE\" BETWEEN DATE '2010-01-01'
                              AND DATE '2025-12-31'
")


dbExecute(con, "
CREATE OR REPLACE TABLE observer_species AS
SELECT
  \"OBSERVER ID\" AS observer_id,
  COUNT(DISTINCT \"COMMON NAME\") AS total_species
FROM ebd_filtered_species
GROUP BY \"OBSERVER ID\"
")

########### CREATE FINAL OBSERVER TABLE #############

dbExecute(con, "
CREATE OR REPLACE TABLE observer_master AS
SELECT
    m.observer_id,
    m.years_active,
    m.n_checklists,
    s.total_species
FROM observer_metrics m
LEFT JOIN observer_species s
ON m.observer_id = s.observer_id
")

observer_master <- dbGetQuery(
  con,
  "SELECT * FROM observer_master"
)

############ CHECK RAW ORI VARIABLES ############

ori <- observer_master[, c(
  "years_active",
  "n_checklists",
  "total_species"
)]

cor(
  ori,
  use = "complete.obs",
  method = "spearman"
)

############ STANDARDIZE ORI PARAMETERS ############
# NEW ORI = years active + checklist count + total species
# Birding days excluded because of high redundancy with checklist count

observer_master$z_years_active <- as.numeric(
  scale(observer_master$years_active)
)

observer_master$z_checklists <- as.numeric(
  scale(log1p(observer_master$n_checklists))
)

observer_master$z_total_species <- as.numeric(
  scale(log1p(observer_master$total_species))
)

vif(
  lm(
    z_years_active ~ z_checklists + z_total_species,
    data = observer_master
  )
)

############ CALCULATE NEW ORI ############
# ORI3a

observer_master$ORI <- rowMeans(
  observer_master[, c(
    "z_years_active",
    "z_checklists",
    "z_total_species"
  )],
  na.rm = FALSE
)

############ BOUND ORI TO 0–1 ############

observer_master$ORI_bounded <- pnorm(
  observer_master$ORI
)

############ CHECK ORI COMPONENTS ############

sapply(
  observer_master[, c(
    "z_years_active",
    "z_checklists",
    "z_total_species"
  )],
  function(x) c(
    mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE)
  )
)

############ CHECK NEW ORI ############

summary(observer_master$ORI)

range(
  observer_master$ORI,
  na.rm = TRUE
)

############ CHECK BOUNDED ORI ############

summary(observer_master$ORI_bounded)

range(
  observer_master$ORI_bounded,
  na.rm = TRUE
)

############ VIEW FINAL TABLE ############

View(observer_master)

############ WRITE UPDATED TABLE BACK TO DUCKDB ############

DBI::dbWriteTable(
  con,
  "observer_master",
  observer_master,
  overwrite = TRUE
)

################# ORI STATS ###########

dbGetQuery(con, "
SELECT
  quantile_cont(n_checklists, 0.5) AS median,
  quantile_cont(n_checklists, 0.9) AS p90,
  quantile_cont(n_checklists, 0.99) AS p99,
  MAX(n_checklists) AS max_checklists
FROM observer_master
")

dbGetQuery(con, "
SELECT
  quantile_cont(years_active, 0.5) AS median,
  quantile_cont(years_active, 0.9) AS p90,
  quantile_cont(years_active, 0.99) AS p99,
  MAX(years_active) AS max_hours
FROM observer_master
")

dbGetQuery(con, "
SELECT
  quantile_cont(total_species, 0.5) AS median,
  quantile_cont(total_species, 0.9) AS p90,
  quantile_cont(total_species, 0.99) AS p99,
  MAX(total_species) AS max_species
FROM observer_master
")

dbGetQuery(con, "
SELECT
  quantile_cont(n_birding_days, 0.5) AS median,
  quantile_cont(n_birding_days, 0.9) AS p90,
  quantile_cont(n_birding_days, 0.99) AS p99,
  MAX(n_birding_days) AS max_years
FROM observer_master
")

#=============== RICHNESS AND ABUNDANCE CAS SCORES =============#

# ============================================================
# STEP 1. CHECK CHECKLIST MASTER UNIQUENESS
# ============================================================

dbGetQuery(con, "
    SELECT
        COUNT(*) AS total_rows,
        COUNT(DISTINCT checklist_id) AS unique_checklists,
        SUM(CASE WHEN checklist_id IS NULL THEN 1 ELSE 0 END)
            AS missing_checklist_id
    FROM checklist_master
")

# ============================================================
# STEP 2. BUILD STATE-CODE LOOKUP FROM EBD_FILTERED
# ============================================================

dbExecute(con, "
    CREATE OR REPLACE TABLE cas_state_lookup AS
    SELECT
        \"SAMPLING EVENT IDENTIFIER\" AS checklist_id,
        MIN(\"STATE CODE\") AS state_code,
        COUNT(DISTINCT \"STATE CODE\") AS n_state_codes
    FROM ebd_filtered
    WHERE
        \"SAMPLING EVENT IDENTIFIER\" IS NOT NULL
        AND \"STATE CODE\" IS NOT NULL
    GROUP BY \"SAMPLING EVENT IDENTIFIER\"
")

# Check whether any checklist has conflicting state codes
dbGetQuery(con, "
    SELECT COUNT(*) AS checklists_with_multiple_state_codes
    FROM cas_state_lookup
    WHERE n_state_codes > 1
")

# ============================================================
# STEP 3. DERIVE CHECKLIST-LEVEL CAS INPUT TABLE
# ============================================================

dbExecute(con, "
    CREATE OR REPLACE TABLE cas_final_checklists AS

    SELECT
        c.checklist_id,
        c.observer_id,
        c.protocol,
        s.state_code,
        c.obs_date,

        CASE
            WHEN MONTH(c.obs_date) IN (11, 12, 1, 2)
                THEN 'Winter'
            WHEN MONTH(c.obs_date) IN (3, 4, 5)
                THEN 'Pre-monsoon'
            WHEN MONTH(c.obs_date) IN (6, 7, 8, 9)
                THEN 'Monsoon'
            WHEN MONTH(c.obs_date) = 10
                THEN 'Post-monsoon'
        END AS season,

        c.richness,
        c.abundance_per_checklist AS abundance,
        c.duration_min,
        c.distance_km,

        LN(1 + c.duration_min) AS log_duration,

        CASE
            WHEN c.protocol = 'Traveling'
                 AND c.distance_km IS NOT NULL
                 AND c.distance_km >= 0
                THEN LN(1 + c.distance_km)
            ELSE 0
        END AS log_distance

    FROM checklist_master c

    INNER JOIN cas_state_lookup s
        ON c.checklist_id = s.checklist_id

    WHERE
        c.checklist_id IS NOT NULL
        AND c.richness IS NOT NULL
        AND c.abundance_per_checklist IS NOT NULL
        AND c.duration_min IS NOT NULL
        AND c.duration_min >= 0
        AND c.protocol IS NOT NULL
        AND c.obs_date IS NOT NULL
        AND s.state_code IS NOT NULL
")

# ============================================================
# STEP 4. VALIDATE THE RECONSTRUCTED TABLE
# ============================================================

dbGetQuery(con, "
    SELECT
        COUNT(*) AS n_rows,
        COUNT(DISTINCT checklist_id) AS unique_checklists,
        SUM(CASE WHEN richness IS NULL THEN 1 ELSE 0 END)
            AS missing_richness,
        SUM(CASE WHEN abundance IS NULL THEN 1 ELSE 0 END)
            AS missing_abundance,
        SUM(CASE WHEN state_code IS NULL THEN 1 ELSE 0 END)
            AS missing_state_code,
        SUM(CASE WHEN season IS NULL THEN 1 ELSE 0 END)
            AS missing_season,
        SUM(CASE WHEN log_duration IS NULL THEN 1 ELSE 0 END)
            AS missing_log_duration,
        SUM(CASE WHEN log_distance IS NULL THEN 1 ELSE 0 END)
            AS missing_log_distance
    FROM cas_final_checklists
")

# Protocol diagnostics
dbGetQuery(con, "
    SELECT
        protocol,
        COUNT(*) AS n,
        SUM(CASE WHEN distance_km IS NULL THEN 1 ELSE 0 END)
            AS missing_distance
    FROM cas_final_checklists
    GROUP BY protocol
    ORDER BY protocol
")

# Inspect seasons and states
dbGetQuery(con, "
    SELECT season, COUNT(*) AS n
    FROM cas_final_checklists
    GROUP BY season
    ORDER BY season
")

dbGetQuery(con, "
    SELECT state_code, COUNT(*) AS n
    FROM cas_final_checklists
    GROUP BY state_code
    ORDER BY n DESC
    LIMIT 20
")

# ============================================================
# STEP 5. LOAD CHECKLIST DATA FROM DUCKDB
# ============================================================

cas_final <- dbGetQuery(con, "
    SELECT *
    FROM cas_final_checklists
")

# Set factor variables
cas_final$protocol <- factor(cas_final$protocol)
cas_final$state_code <- factor(cas_final$state_code)
cas_final$season <- factor(
  cas_final$season,
  levels = c(
    'Winter',
    'Pre-monsoon',
    'Monsoon',
    'Post-monsoon'
  )
)

# ============================================================
# STEP 6. FIT POOLED RICHNESS/ ABUNDANCE MODEL
# ============================================================

pooled_richness_final <- glm(
  richness ~ protocol + state_code + season +
    log_duration + log_distance,
  family = quasipoisson(link = "log"),
  data = cas_final
)

pooled_abundance_final <- glm(
  abundance ~ protocol + state_code + season +
    log_duration + log_distance,
  family = quasipoisson(link = "log"),
  data = cas_final
)

# ============================================================
# STEP 7. PREDICT EXPECTED RICHNESS AND ABUNDANCE
# ============================================================

cas_final$expected_richness <- predict(
  pooled_richness_final,
  newdata = cas_final,
  type = "response"
)

cas_final$expected_abundance <- predict(
  pooled_abundance_final,
  newdata = cas_final,
  type = "response"
)

# ============================================================
# STEP 8. CALCULATE RESIDUALS
# ============================================================

cas_final$richness_residual <-
  cas_final$richness - cas_final$expected_richness

cas_final$abundance_residual <-
  cas_final$abundance - cas_final$expected_abundance

# ============================================================
# STEP 9. ROBUST Z-SCORES USING MEDIAN AND MAD
# ============================================================

richness_median <- median(
  cas_final$richness_residual,
  na.rm = TRUE
)

richness_mad <- mad(
  cas_final$richness_residual,
  na.rm = TRUE
)

abundance_median <- median(
  cas_final$abundance_residual,
  na.rm = TRUE
)

abundance_mad <- mad(
  cas_final$abundance_residual,
  na.rm = TRUE
)

if (richness_mad == 0 || abundance_mad == 0) {
  stop("A residual MAD is zero; inspect model residuals before scoring.")
}

cas_final$richness_robust_z <-
  (cas_final$richness_residual - richness_median) /
  (1.4826 * richness_mad)

cas_final$abundance_robust_z <-
  (cas_final$abundance_residual - abundance_median) /
  (1.4826 * abundance_mad)

# ============================================================
# STEP 10. CALCULATE COMPONENT SCORES
# ============================================================

cas_final$richness_score <- ifelse(
  cas_final$richness_robust_z <= 0,
  1,
  exp(-cas_final$richness_robust_z / 3.5)
)

cas_final$abundance_score <- ifelse(
  cas_final$abundance_robust_z <= 0,
  1,
  exp(-cas_final$abundance_robust_z / 3.5)
)

# ============================================================
# STEP 11. FINAL CHECKLIST AUTHENTICITY SCORE
# ============================================================

cas_final$CAS <- (
  cas_final$richness_score +
    cas_final$abundance_score
) / 2

# ============================================================
# STEP 12. SUMMARY DIAGNOSTICS
# ============================================================

summary_cas <- data.frame(
  n = nrow(cas_final),
  min_CAS = min(cas_final$CAS, na.rm = TRUE),
  max_CAS = max(cas_final$CAS, na.rm = TRUE),
  mean_CAS = mean(cas_final$CAS, na.rm = TRUE),
  median_CAS = median(cas_final$CAS, na.rm = TRUE),
  missing_CAS = sum(is.na(cas_final$CAS))
)

print(summary_cas)

# Protocol-level diagnostics
print(
  aggregate(
    cbind(
      CAS,
      richness_score,
      abundance_score
    ) ~ protocol,
    data = cas_final,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
)

# ============================================================
# STEP 13. SAVE RESULTS BACK TO DUCKDB
# ============================================================

dbWriteTable(
  con,
  "cas_final",
  cas_final,
  overwrite = TRUE
)

# ============================================================
# STEP 14. EXPORT RESULTS
# ============================================================

dbExecute(con, "
    COPY cas_final
    TO 'file path/CAS_FINAL.csv'
    WITH (HEADER, DELIMITER ',')
")

dbExecute(con, "
    COPY cas_final
    TO 'file path/CAS_FINAL.txt'
    WITH (HEADER, DELIMITER '\t')
")

# ============================================================
# FLAG ANOMALOUS CHECKLISTS USING EXISTING CAS RESULTS
# ============================================================

# 1. Individual anomaly flags
cas_final$richness_flag <- as.integer(
  cas_final$richness_robust_z >= 3.5
)

cas_final$abundance_flag <- as.integer(
  cas_final$abundance_robust_z >= 3.5
)

# 2. Both richness and abundance anomalous
cas_final$both_flags <- as.integer(
  cas_final$richness_flag == 1 &
    cas_final$abundance_flag == 1
)

# 3. Mutually exclusive anomaly categories
cas_final$anomaly_category <- ifelse(
  cas_final$both_flags == 1,
  "Both richness and abundance",
  ifelse(
    cas_final$richness_flag == 1,
    "Richness only",
    ifelse(
      cas_final$abundance_flag == 1,
      "Abundance only",
      "Not flagged"
    )
  )
)

# ============================================================
# 1. SUMMARISE FLAGS
# ============================================================

flag_summary <- data.frame(
  flag = c(
    "Richness anomalous",
    "Abundance anomalous",
    "Both anomalous"
  ),
  checklists = c(
    sum(cas_final$richness_flag, na.rm = TRUE),
    sum(cas_final$abundance_flag, na.rm = TRUE),
    sum(cas_final$both_flags, na.rm = TRUE)
  )
)

flag_summary$percentage <- round(
  100 * flag_summary$checklists / nrow(cas_final),
  3
)

print(flag_summary)

cat("\nAnomaly category counts:\n")
print(table(cas_final$anomaly_category, useNA = "ifany"))

# ============================================================
# 2. SAVE RESULTS TO DUCKDB
# ============================================================

dbWriteTable(
  con,
  "cas_final",
  cas_final,
  overwrite = TRUE
)

# ============================================================
# 3. EXPORT COMPLETE RESULTS
# ============================================================

dbExecute(con, "
    COPY cas_final
    TO 'file path/CAS_FINAL_WITH_FLAGS.csv'
    WITH (HEADER, DELIMITER ',')
")

# ============================================================
# 4. MEAN CAS PER OBSERVER
# ============================================================

dbExecute(con, "
CREATE OR REPLACE TABLE observer_CAS AS

SELECT

    observer_id,

    COUNT(*) AS n_checklists,

    AVG(CAS) AS mean_CAS,

    MEDIAN(CAS) AS median_CAS,

    STDDEV_SAMP(CAS) AS sd_CAS,

    MIN(CAS) AS min_CAS,

    MAX(CAS) AS max_CAS

FROM cas_final

WHERE
    CAS IS NOT NULL
    AND observer_id IS NOT NULL

GROUP BY
    observer_id
")

dbGetQuery(con, "
SELECT *
FROM observer_CAS
LIMIT 10
")

dbGetQuery(con, "
SELECT
    MIN(mean_CAS) AS min,
    QUANTILE(mean_CAS, 0.25) AS Q1,
    MEDIAN(mean_CAS) AS median,
    AVG(mean_CAS) AS mean,
    QUANTILE(mean_CAS, 0.75) AS Q3,
    MAX(mean_CAS) AS max
FROM observer_CAS
")

dbExecute(con, "
COPY observer_CAS
TO 'file path/OBSERVER_mean CAS.csv'
WITH (
    HEADER,
    DELIMITER ','
)
")

dbExecute(con, "
COPY observer_CAS
TO 'file path/OBSERVER_mean CAS.txt'
WITH (
    HEADER,
    DELIMITER '\t'
)
")

#============== OBSERVER-LEVEL CHECKLIST ANOMALY FLAGS ============#
# ============================================================
# 1. IMPORT COMPLETE CAS RESULTS WITH ANOMALY FLAGS
# ============================================================

dbExecute(con, "
CREATE OR REPLACE TABLE cas_final_import AS
SELECT *
FROM read_csv(
    'file path/CAS_FINAL_WITH_FLAGS.csv',
    header = TRUE,
    strict_mode = FALSE,
    null_padding = TRUE
)
")

# Verify the imported columns
print(dbGetQuery(con, "DESCRIBE cas_final_import"))

# Check record counts
print(dbGetQuery(con, "
SELECT
    COUNT(*) AS n_rows,
    COUNT(DISTINCT checklist_id) AS n_checklists,
    COUNT(DISTINCT observer_id) AS n_observers
FROM cas_final_import
"))

# ============================================================
# 2. CREATE OBSERVER-LEVEL CAS AND ANOMALY SUMMARY
# ============================================================

dbExecute(con, "
CREATE OR REPLACE TABLE observer_CAS_anomalies AS
SELECT
    observer_id,

    COUNT(*) AS n_checklists,

    COUNT(richness_robust_z) AS n_richness_evaluable,
    COUNT(abundance_robust_z) AS n_abundance_evaluable,

    SUM(richness_flag) AS richness_flag_count,
    SUM(abundance_flag) AS abundance_flag_count,
    SUM(both_flags) AS both_flag_count,

    AVG(richness_flag) AS richness_anomaly_rate,
    AVG(abundance_flag) AS abundance_anomaly_rate,
    AVG(both_flags) AS both_anomaly_rate,

    AVG(richness_score) AS mean_richness_score,
    AVG(abundance_score) AS mean_abundance_score,

    AVG(CAS) AS mean_CAS,
    MEDIAN(CAS) AS median_CAS,
    STDDEV_SAMP(CAS) AS sd_CAS,
    MIN(CAS) AS min_CAS,
    MAX(CAS) AS max_CAS

FROM cas_final_import
WHERE observer_id IS NOT NULL
GROUP BY observer_id
")

# ============================================================
# 3. CREATE OVERALL CHECKLIST-LEVEL ANOMALY SUMMARY
# ============================================================

dbExecute(con, "
CREATE OR REPLACE TABLE overall_CAS_anomalies AS
SELECT
    COUNT(*) AS n_checklists,

    COUNT(richness_robust_z) AS n_richness_evaluable,
    COUNT(abundance_robust_z) AS n_abundance_evaluable,

    SUM(richness_flag) AS richness_flagged,
    SUM(abundance_flag) AS abundance_flagged,
    SUM(both_flags) AS both_flagged,

    AVG(richness_flag) AS richness_anomaly_rate,
    AVG(abundance_flag) AS abundance_anomaly_rate,
    AVG(both_flags) AS both_anomaly_rate,

    AVG(richness_score) AS mean_richness_score,
    AVG(abundance_score) AS mean_abundance_score,
    AVG(CAS) AS mean_CAS,
    MEDIAN(CAS) AS median_CAS

FROM cas_final_import
")

# ============================================================
# 4. CONSOLIDATE CHECKLIST, OBSERVER AND OVERALL RESULTS
# ============================================================

# Use the checklist table's columns as the shared schema.
# Observer-level and overall records populate their relevant fields;
# non-applicable fields remain NULL.

dbExecute(con, "
CREATE OR REPLACE TABLE CAS_CONSOLIDATED AS

-- CHECKLIST-LEVEL RECORDS
SELECT
    'CHECKLIST' AS record_type,
    checklist_id,
    observer_id,
    protocol,
    state_code,
    obs_date,
    season,

    richness,
    abundance,
    duration_min,
    distance_km,

    expected_richness,
    expected_abundance,
    richness_residual,
    abundance_residual,
    richness_robust_z,
    abundance_robust_z,

    richness_score,
    abundance_score,
    CAS,

    richness_flag,
    abundance_flag,
    both_flags,
    anomaly_category,

    NULL::BIGINT AS n_checklists,
    NULL::BIGINT AS n_richness_evaluable,
    NULL::BIGINT AS n_abundance_evaluable,

    NULL::BIGINT AS richness_flag_count,
    NULL::BIGINT AS abundance_flag_count,
    NULL::BIGINT AS both_flag_count,

    NULL::DOUBLE AS richness_anomaly_rate,
    NULL::DOUBLE AS abundance_anomaly_rate,
    NULL::DOUBLE AS both_anomaly_rate,

    NULL::DOUBLE AS mean_richness_score,
    NULL::DOUBLE AS mean_abundance_score,
    NULL::DOUBLE AS mean_CAS,
    NULL::DOUBLE AS median_CAS,
    NULL::DOUBLE AS sd_CAS,
    NULL::DOUBLE AS min_CAS,
    NULL::DOUBLE AS max_CAS

FROM cas_final_import

UNION ALL

-- OBSERVER-LEVEL RECORDS
SELECT
    'OBSERVER' AS record_type,
    NULL::VARCHAR AS checklist_id,
    observer_id,
    NULL::VARCHAR AS protocol,
    NULL::VARCHAR AS state_code,
    NULL AS obs_date,
    NULL::VARCHAR AS season,

    NULL::BIGINT AS richness,
    NULL::DOUBLE AS abundance,
    NULL::DOUBLE AS duration_min,
    NULL::DOUBLE AS distance_km,

    NULL::DOUBLE AS expected_richness,
    NULL::DOUBLE AS expected_abundance,
    NULL::DOUBLE AS richness_residual,
    NULL::DOUBLE AS abundance_residual,
    NULL::DOUBLE AS richness_robust_z,
    NULL::DOUBLE AS abundance_robust_z,

    NULL::DOUBLE AS richness_score,
    NULL::DOUBLE AS abundance_score,
    NULL::DOUBLE AS CAS,

    NULL::INTEGER AS richness_flag,
    NULL::INTEGER AS abundance_flag,
    NULL::INTEGER AS both_flags,
    NULL::VARCHAR AS anomaly_category,

    n_checklists,
    n_richness_evaluable,
    n_abundance_evaluable,

    richness_flag_count,
    abundance_flag_count,
    both_flag_count,

    richness_anomaly_rate,
    abundance_anomaly_rate,
    both_anomaly_rate,

    mean_richness_score,
    mean_abundance_score,
    mean_CAS,
    median_CAS,
    sd_CAS,
    min_CAS,
    max_CAS

FROM observer_CAS_anomalies

UNION ALL

-- OVERALL RECORD
SELECT
    'OVERALL' AS record_type,
    NULL::VARCHAR AS checklist_id,
    NULL::VARCHAR AS observer_id,
    NULL::VARCHAR AS protocol,
    NULL::VARCHAR AS state_code,
    NULL AS obs_date,
    NULL::VARCHAR AS season,

    NULL::BIGINT AS richness,
    NULL::DOUBLE AS abundance,
    NULL::DOUBLE AS duration_min,
    NULL::DOUBLE AS distance_km,

    NULL::DOUBLE AS expected_richness,
    NULL::DOUBLE AS expected_abundance,
    NULL::DOUBLE AS richness_residual,
    NULL::DOUBLE AS abundance_residual,
    NULL::DOUBLE AS richness_robust_z,
    NULL::DOUBLE AS abundance_robust_z,

    NULL::DOUBLE AS richness_score,
    NULL::DOUBLE AS abundance_score,
    mean_CAS AS CAS,

    NULL::INTEGER AS richness_flag,
    NULL::INTEGER AS abundance_flag,
    NULL::INTEGER AS both_flags,
    NULL::VARCHAR AS anomaly_category,

    n_checklists,
    n_richness_evaluable,
    n_abundance_evaluable,

    richness_flagged AS richness_flag_count,
    abundance_flagged AS abundance_flag_count,
    both_flagged AS both_flag_count,

    richness_anomaly_rate,
    abundance_anomaly_rate,
    both_anomaly_rate,

    mean_richness_score,
    mean_abundance_score,
    mean_CAS,
    median_CAS,
    NULL::DOUBLE AS sd_CAS,
    NULL::DOUBLE AS min_CAS,
    NULL::DOUBLE AS max_CAS

FROM overall_CAS_anomalies
")

dbExecute(con, "
COPY CAS_CONSOLIDATED
TO 'file path/CAS_CONSOLIDATED_RESULTS.csv'
WITH (HEADER, DELIMITER ',')
")

#================== COMPLETENESS INDEX ==================#
          
          dbExecute(con, "
CREATE OR REPLACE TABLE ebd_raw AS
SELECT *
FROM read_csv(
  'file path/ebd raw.txt',
  delim='\t',
  header=TRUE,
  quote='',
  strict_mode=false
)
")

################# COMPLETENESS SCORE (CS) #################

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_raw AS
SELECT *
FROM read_csv(
  'file path/ebd raw.txt',
  delim='\t',
  header=TRUE,
  quote='',
  strict_mode=false
)
")

dbExecute(con, "
CREATE OR REPLACE TABLE completeness_checklist AS
SELECT
  \"OBSERVER ID\" AS observer_id,

  COUNT(DISTINCT \"SAMPLING EVENT IDENTIFIER\") AS total_checklists,

  COUNT(DISTINCT CASE
      WHEN \"ALL SPECIES REPORTED\" = 1
      THEN \"SAMPLING EVENT IDENTIFIER\"
  END) AS complete_checklists,

  COUNT(DISTINCT CASE
      WHEN \"ALL SPECIES REPORTED\" = 0
      THEN \"SAMPLING EVENT IDENTIFIER\"
  END) AS incomplete_checklists,

  COUNT(DISTINCT CASE
      WHEN \"ALL SPECIES REPORTED\" = 1
      THEN \"SAMPLING EVENT IDENTIFIER\"
  END) * 1.0
  /
  COUNT(DISTINCT \"SAMPLING EVENT IDENTIFIER\")
  AS completeness_score

FROM ebd_raw

WHERE
  \"OBSERVATION DATE\" BETWEEN DATE '2010-01-01' AND DATE '2025-12-31'
  AND APPROVED = 1
  AND \"PROTOCOL NAME\" IN ('Traveling','Stationary')
  AND \"DURATION MINUTES\" BETWEEN 10 AND 300
  AND \"NUMBER OBSERVERS\" <= 10
  AND (
        \"PROTOCOL NAME\" = 'Stationary'
        OR COALESCE(\"EFFORT DISTANCE KM\",0) <= 8
      )

GROUP BY \"OBSERVER ID\"
")

completeness_checklist <- dbGetQuery(con,"
SELECT *
FROM completeness_checklist
")

############# CHECK CHECKLIST DISTRIBUTION ACROSS PERCENTILES AND MEDIAN ###########

dbGetQuery(con, "
SELECT
    MIN(total_checklists) AS min_chk,
    PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY total_checklists) AS q1,
    MEDIAN(total_checklists) AS median_chk,
    PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY total_checklists) AS q3,
    MAX(total_checklists) AS max_chk
FROM completeness_checklist
")

############### SENSITIVITY ANALYSIS #################

comp <- completeness_checklist %>%
  select(
    observer_id,
    completeness_score,
    total_checklists
  )

taus <- c(10,20,30,40,50)

for(t in taus){       # Calculate weighted scores for all tau values
  
  comp[[paste0("weight_",t)]] <-
    comp$completeness_score *
    (1-exp(-comp$total_checklists/t))
  
}

head(comp)

weights <- comp %>%
  select(starts_with("weight_"))

cor(
  weights,
  method = "spearman"
)

############ PAIRWISE SIGNIFICANCE TESTS ############

results <- data.frame(
  Tau_1 = character(),
  Tau_2 = character(),
  Spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(i in 1:(length(taus)-1)){
  for(j in (i+1):length(taus)){
    
    test <- cor.test(
      
      comp[[paste0("weight_",taus[i])]],
      comp[[paste0("weight_",taus[j])]],
      
      method = "spearman"
    )
    results <- rbind(
      results,
      data.frame(
        Tau_1 = paste0("Tau ",taus[i]),
        Tau_2 = paste0("Tau ",taus[j]),
        Spearman_rho = round(unname(test$estimate),4),
        p_value = signif(test$p.value,3)
      )
    )
  }
}

print(results)

############ FINAL WEIGHTED CHECKLIST COMPLETENESS INDEX (CCI) ############

completeness_checklist$weighted_completeness_score <-   # Tau value of 20 checklists was selected based on sensitivity analysis and median checklist number (5 checklists)
  
  completeness_checklist$completeness_score *
  (
    1-exp(
      -completeness_checklist$total_checklists/20
    )
  )

head(completeness_checklist)

write.csv(
  completeness_checklist,
  "file path/completeness new.csv",
  row.names = FALSE
)

#==================  CHECKLIST CONSISTENCY =================#

dbExecute(con, "
CREATE OR REPLACE TABLE checklist_gaps AS
SELECT
  observer_id,
  checklist_id,
  obs_date,
  LAG(obs_date) OVER (
    PARTITION BY observer_id
    ORDER BY obs_date
  ) AS prev_date,
  DATE_DIFF(
    'day',
    LAG(obs_date) OVER (
      PARTITION BY observer_id
      ORDER BY obs_date
    ),
    obs_date
  ) AS gap_days
FROM checklist_master
")

############# CALCULATE ROBUST STATISICS AND CONSISTENCY SCORE ############

dbExecute(con, "
CREATE OR REPLACE TABLE checklist_regularity AS
WITH observer_gaps AS (

SELECT
    observer_id,
    median(gap_days) AS median_gap_days
FROM checklist_gaps
WHERE gap_days IS NOT NULL
GROUP BY observer_id
),

observer_mad AS (
SELECT
    g.observer_id,
    median(abs(g.gap_days - o.median_gap_days)) AS mad_gap_days
FROM checklist_gaps g
JOIN observer_gaps o
ON g.observer_id = o.observer_id
WHERE g.gap_days IS NOT NULL
GROUP BY g.observer_id
)

SELECT
o.observer_id,
o.median_gap_days,
m.mad_gap_days,
(1.4826 * m.mad_gap_days) /
(o.median_gap_days + 1) AS robust_cv,

EXP(
-
(
(1.4826 * m.mad_gap_days) /
(o.median_gap_days + 1)
)
) AS consistency_score

FROM observer_gaps o
JOIN observer_mad m
ON o.observer_id = m.observer_id
")

regularity <- dbGetQuery(con, "
SELECT consistency_score
FROM checklist_regularity
")

checklist_regularity <- dbGetQuery(con, "
SELECT *
FROM checklist_regularity
")
View(checklist_regularity)

summary(regularity$consistency_score)

write.csv(
  checklist_regularity,
  "file path/consistency new.csv",
  row.names = FALSE
)

############### SAMPLING BEHAVIOUR INDEX ################

consistency <- fread("file path/consistency new.csv")
completeness <- fread("file path/completeness new.csv")

consistency[, observer_id := as.character(observer_id)]
completeness[, observer_id := as.character(observer_id)]

# Merge observer-level components
sbi_master <- merge(
  consistency,
  completeness,
  by = "observer_id",
  all = TRUE
)

# Calculate SBI only when both components are available
sbi_master[, SBI := fifelse(
  !is.na(consistency_score) &
    !is.na(weighted_completeness_score),
  (consistency_score +
     weighted_completeness_score) / 2,
  NA_real_
)]

# Check results
summary(sbi_master$SBI)
range(sbi_master$SBI, na.rm = TRUE)

head(sbi_master[, .(
  observer_id,
  consistency_score,
  weighted_completeness_score,
  SBI
)])

# Export
fwrite(
  sbi_master,
  "file path/sbi_master.csv"
)

############### OBSERVER MASTER TABLE AND CHECKLIST MASTER TABLE JOIN ##########

ori_master <- fread("file path/ori_master.csv")
cas_master <- fread("file path/OBSERVER_mean CAS.csv")
sbi_master <- fread("file path/sbi_master.csv")

# ============================================================
# STANDARDIZE OBSERVER IDs
# ============================================================

ori_master[, observer_id := trimws(as.character(observer_id))]
cas_master[, observer_id := trimws(as.character(observer_id))]
sbi_master[, observer_id := trimws(as.character(observer_id))]

# ============================================================
# CHECK NUMBER OF OBSERVERS IN EACH TABLE
# ============================================================

cat("\nNumber of unique observers:\n")

cat("ORI :", uniqueN(ori_master$observer_id), "\n")
cat("CAS :", uniqueN(cas_master$observer_id), "\n")
cat("SBI :", uniqueN(sbi_master$observer_id), "\n")

# ============================================================
# CHECK WHETHER THERE ARE DUPLICATE OBSERVER IDs
# ============================================================

cat("\nDuplicate observer IDs:\n")

cat("ORI :", ori_master[, anyDuplicated(observer_id) > 0], "\n")
cat("CAS :", cas_master[, anyDuplicated(observer_id) > 0], "\n")
cat("SBI :", sbi_master[, anyDuplicated(observer_id) > 0], "\n")

# ============================================================
# FIND OBSERVERS PRESENT IN ALL THREE TABLES
# ============================================================

ori_ids <- unique(ori_master$observer_id)
cas_ids <- unique(cas_master$observer_id)
sbi_ids <- unique(sbi_master$observer_id)

common_ids <- Reduce(
  intersect,
  list(ori_ids, cas_ids, sbi_ids)
)

cat("\nObservers present in ALL THREE:\n")
cat(length(common_ids), "\n")

# ============================================================
# OBSERVERS MISSING FROM EACH TABLE
# ============================================================

missing_cas <- setdiff(ori_ids, cas_ids)
missing_sbi <- setdiff(ori_ids, sbi_ids)

missing_ori_from_cas <- setdiff(cas_ids, ori_ids)
missing_ori_from_sbi <- setdiff(sbi_ids, ori_ids)

cat("\nObservers in ORI but NOT CAS :", length(missing_cas), "\n")
cat("Observers in ORI but NOT SBI :", length(missing_sbi), "\n")

cat("\nObservers in CAS but NOT ORI :", length(missing_ori_from_cas), "\n")
cat("Observers in SBI but NOT ORI :", length(missing_ori_from_sbi), "\n")

# ============================================================
# CREATE MATCHING TABLE
# ============================================================

observer_matching <- data.table(
  observer_id = sort(
    unique(
      c(ori_ids, cas_ids, sbi_ids)
    )
  )
)

observer_matching[, in_ORI := observer_id %in% ori_ids]
observer_matching[, in_CAS := observer_id %in% cas_ids]
observer_matching[, in_SBI := observer_id %in% sbi_ids]

observer_matching[, in_all_three :=
                    in_ORI & in_CAS & in_SBI
]

# ============================================================
# SUMMARY OF MATCHING
# ============================================================

cat("\nObserver matching summary:\n")

print(
  observer_matching[
    ,
    .(
      total_observers = .N,
      ORI = sum(in_ORI),
      CAS = sum(in_CAS),
      SBI = sum(in_SBI),
      all_three = sum(in_all_three)
    )
  ]
)

# ============================================================
# SHOW OBSERVERS NOT PRESENT IN ALL THREE
# ============================================================

observers_left_out <- observer_matching[
  in_all_three == FALSE
]

cat("\nObservers NOT present in all three datasets:\n")
cat(nrow(observers_left_out), "\n")

if (nrow(observers_left_out) > 0) {
  print(observers_left_out)
}

# ============================================================
# KEEP ONLY REQUIRED COLUMNS FROM ORI
# ============================================================

ori_final <- ori_master[
  ,
  .(
    observer_id,
    years_active,
    n_checklists,
    total_species,
    ORI,
    ORI_bounded
  )
]

# ============================================================
# KEEP ONLY CAS
# ============================================================

cas_final <- cas_master[
  ,
  .(
    observer_id,
    mean_CAS
  )
]

# ============================================================
# KEEP SBI AND ITS COMPONENTS
# ============================================================

sbi_final <- sbi_master[
  ,
  .(
    observer_id,
    consistency_score,
    weighted_completeness_score,
    SBI
  )
]

# ============================================================
# MERGE ORI + CAS + SBI
# ============================================================

observer_final <- merge(
  ori_final,
  cas_final,
  by = "observer_id",
  all = FALSE
)

observer_final <- merge(
  observer_final,
  sbi_final[!is.na(SBI), ],
  by = "observer_id",
  all = FALSE
)

cat("Final observers:", nrow(observer_final), "\n")

# Verify complete three-component data
stopifnot(
  all(!is.na(observer_final$ORI_bounded)),
  all(!is.na(observer_final$mean_CAS)),
  all(!is.na(observer_final$SBI))
)

# ============================================================
# CHECK FOR NA AFTER MERGING
# ============================================================

cat("\nNA counts after merging:\n")

print(
  observer_final[
    ,
    .(
      ORI_NA = sum(is.na(ORI_bounded)),
      CAS_NA = sum(is.na(mean_CAS)),
      SBI_NA = sum(is.na(SBI)),
      Consistency_NA = sum(is.na(consistency_score)),
      Completeness_NA = sum(is.na(weighted_completeness_score))
    )
  ]
)

# ============================================================
# KEEP ONLY OBSERVERS WITH ALL THREE FINAL COMPONENTS
# ============================================================

observer_final <- observer_final[
  !is.na(ORI_bounded) &
    !is.na(mean_CAS) &
    !is.na(SBI)
]

# ============================================================
# CALCULATE NEW AVICRED
# ============================================================

observer_final[
  ,
  AviCRED := (ORI_bounded + mean_CAS + SBI) / 3
]

# ============================================================
# CHECK FINAL AVICRED
# ============================================================

summary(observer_final$AviCRED)

range(
  observer_final$AviCRED,
  na.rm = TRUE
)

cat("\nFinal number of observers:\n")
cat(nrow(observer_final), "\n")

# ============================================================
# WRITE TO DUCKDB
# ============================================================

dbWriteTable(
  con,
  "observer_final",
  observer_final,
  overwrite = TRUE
)

# Classify observers by checklist history
observer_final <- observer_final %>%
  mutate(
    checklist_history_status = if_else(
      n_checklists >= 50,
      "Sufficient checklist history",
      "Insufficient checklist history"
    )
  )

# Check classification counts
cat("\nChecklist history classification:\n")

print(
  observer_final %>%
    count(checklist_history_status)
)

# Summarise checklist history by classification
print(
  observer_final %>%
    group_by(checklist_history_status) %>%
    summarise(
      observers = n(),
      min_checklists = min(n_checklists, na.rm = TRUE),
      median_checklists = median(n_checklists, na.rm = TRUE),
      max_checklists = max(n_checklists, na.rm = TRUE),
      .groups = "drop"
    )
)

write_csv(
  observer_final,
  "file path/Avicred new.csv"
)

############ SCORE STATS ###########

dbGetQuery(con, "
SELECT

    --------------------------------------------------
    -- ORI
    --------------------------------------------------

    MIN(ORI_bounded) AS min_ori,
    quantile_cont(ORI_bounded, 0.5) AS median_ori,
    MAX(ORI_bounded) AS max_ori,

    --------------------------------------------------
    -- CAS
    --------------------------------------------------

    MIN(mean_CAS) AS min_cas,
    quantile_cont(mean_CAS, 0.5) AS median_cas,
    MAX(mean_CAS) AS max_cas,

    --------------------------------------------------
    -- SBI
    --------------------------------------------------

    MIN(SBI) AS min_sbi,
    quantile_cont(SBI, 0.5) AS median_sbi,
    MAX(SBI) AS max_sbi,

    --------------------------------------------------
    -- FINAL AVICRED
    --------------------------------------------------

    MIN(AviCRED) AS min_avicred,
    quantile_cont(AviCRED, 0.5) AS median_avicred,
    MAX(AviCRED) AS max_avicred

FROM observer_final
")

############# NEW AVICRED SENSITIVITY ANALYSIS #############

observer_final <- read.csv(
  "file path/Observer final new.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

names(observer_final)

analysis_components <- observer_final[, c(
  "observer_id",
  "ORI_bounded",
  "mean_CAS",
  "consistency_score",
  "weighted_completeness_score",
  "SBI",
  "AviCRED"
)]

setnames(
  analysis_components,
  old = c(
    "ORI_bounded",
    "mean_CAS",
    "consistency_score",
    "weighted_completeness_score"
  ),
  new = c(
    "ORI",
    "CAS",
    "Consistency",
    "Completeness"
  )
)

# Select the actual columns and assign clear labels
analysis_components_cor <- data.frame(
  ORI          = observer_final$ORI_bounded,
  CAS          = observer_final$mean_CAS,
  SBI          = observer_final$SBI,
  Consistency  = observer_final$consistency_score,
  Completeness = observer_final$weighted_completeness_score
)

# Spearman correlation matrix
cor_matrix <- cor(
  analysis_components_cor,
  use = "pairwise.complete.obs",
  method = "spearman"
)

cor_results <- rcorr(
  as.matrix(analysis_components_cor),
  type = "spearman"
)

round(cor_results$r, 3)  # Correlation coefficients
round(cor_results$P, 4)  # P-values

#========= AVICRED COMPONENT CONTRIBUTION, RANK SENSITIVITY AND BOOTSTRAP =======#

output_dir <- "file path/ABLATION"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(20261017)
B <- 2000

# ------------------------------------------------------------
# 1. SELECT COMPLETE OBSERVER RECORDS
# ------------------------------------------------------------

required_vars <- c(
  "observer_id",
  "years_active",
  "n_checklists",
  "total_species",
  "consistency_score",
  "weighted_completeness_score",
  "mean_CAS",
  "ORI_bounded",
  "SBI",
  "AviCRED"
)

if (!all(required_vars %in% names(observer_final))) {
  stop(
    "Missing columns: ",
    paste(
      setdiff(required_vars, names(observer_final)),
      collapse = ", "
    )
  )
}

av <- observer_final[
  complete.cases(observer_final[, required_vars]),
  required_vars
]

if (anyDuplicated(av$observer_id)) {
  stop("Duplicate observer IDs detected.")
}

if (nrow(av) < 20) {
  stop("Too few complete observer records for this analysis.")
}

cat("Complete observers:", nrow(av), "\n")

# Six underlying variables
basic_vars <- c(
  "years_active",
  "n_checklists",
  "total_species",
  "consistency_score",
  "weighted_completeness_score",
  "mean_CAS"
)

# ------------------------------------------------------------
# 2. VERIFY ORIGINAL AVICRED FORMULA
# ------------------------------------------------------------

av$AviCRED_reconstructed <- (
  av$ORI_bounded +
    av$mean_CAS +
    av$SBI
) / 3

max_discrepancy <- max(
  abs(av$AviCRED - av$AviCRED_reconstructed)
)

cat(
  "Maximum AviCRED reconstruction discrepancy:",
  format(max_discrepancy, scientific = TRUE),
  "\n"
)

if (max_discrepancy > 1e-8) {
  warning(
    "Reconstructed AviCRED differs from the stored AviCRED score. ",
    "Check the original formula and any rounding."
  )
}

# ------------------------------------------------------------
# 3. DESCRIPTIVE STATISTICS FOR SIX BASIC VARIABLES
# ------------------------------------------------------------

basic_summary <- do.call(
  rbind,
  lapply(basic_vars, function(v) {
    
    x <- av[[v]]
    
    data.frame(
      variable = v,
      n = length(x),
      mean = mean(x),
      SD = sd(x),
      median = median(x),
      minimum = min(x),
      maximum = max(x)
    )
  })
)

rownames(basic_summary) <- NULL

# ------------------------------------------------------------
# 4. CORRELATIONS AMONG BASIC VARIABLES
# ------------------------------------------------------------

basic_cor_pearson <- cor(
  av[, basic_vars],
  method = "pearson"
)

basic_cor_spearman <- cor(
  av[, basic_vars],
  method = "spearman"
)

# ------------------------------------------------------------
# 5. COMPONENT CONTRIBUTIONS TO AVICRED
# ------------------------------------------------------------

# Each component receives an equal weight of 1/3.

av$ORI_contribution <- av$ORI_bounded / 3
av$CAS_contribution <- av$mean_CAS / 3
av$SBI_contribution <- av$SBI / 3

# Verify that the contributions sum to AviCRED.
contribution_discrepancy <- max(
  abs(
    av$ORI_contribution +
      av$CAS_contribution +
      av$SBI_contribution -
      av$AviCRED
  )
)

cat(
  "Maximum contribution-sum discrepancy:",
  format(contribution_discrepancy, scientific = TRUE),
  "\n"
)

contribution_vars <- c(
  "ORI_contribution",
  "CAS_contribution",
  "SBI_contribution"
)

contribution_summary <- do.call(
  rbind,
  lapply(contribution_vars, function(v) {
    
    x <- av[[v]]
    
    data.frame(
      component = v,
      n = length(x),
      mean = mean(x),
      SD = sd(x),
      median = median(x),
      minimum = min(x),
      maximum = max(x)
    )
  })
)

rownames(contribution_summary) <- NULL

# ------------------------------------------------------------
# 6. ORIGINAL AVICRED COMPONENT RANKING SENSITIVITY
# ------------------------------------------------------------

# The denominator remains 3 after component removal.
# This is a fixed-weight ablation: the removed component
# contributes zero, rather than reweighting the other components.

av$no_ORI <- (
  av$mean_CAS + av$SBI
) / 3

av$no_CAS <- (
  av$ORI_bounded + av$SBI
) / 3

av$no_SBI <- (
  av$ORI_bounded + av$mean_CAS
) / 3


rank_metrics <- function(full, reduced, component) {
  
  r_full <- rank(-full, ties.method = "average")
  r_reduced <- rank(-reduced, ties.method = "average")
  
  top10_full <- order(-full)[seq_len(10)]
  top10_reduced <- order(-reduced)[seq_len(10)]
  
  top20_full <- order(-full)[seq_len(20)]
  top20_reduced <- order(-reduced)[seq_len(20)]
  
  data.frame(
    component_removed = component,
    Pearson = cor(full, reduced, method = "pearson"),
    Spearman = cor(full, reduced, method = "spearman"),
    mean_rank_change = mean(abs(r_full - r_reduced)),
    median_rank_change = median(abs(r_full - r_reduced)),
    top10_overlap = length(
      intersect(top10_full, top10_reduced)
    ) / 10,
    top20_overlap = length(
      intersect(top20_full, top20_reduced)
    ) / 20
  )
}


rank_results <- rbind(
  rank_metrics(av$AviCRED, av$no_ORI, "ORI"),
  rank_metrics(av$AviCRED, av$no_CAS, "CAS"),
  rank_metrics(av$AviCRED, av$no_SBI, "SBI")
)

rownames(rank_results) <- NULL

# ------------------------------------------------------------
# 7. BOOTSTRAP MEANS AND SDs OF THE SIX BASIC VARIABLES
# ------------------------------------------------------------

# Observer-level non-parametric bootstrap.
# Each replicate samples observers with replacement.
# The original variable scales are retained for this analysis.

boot_basic <- function(d) {
  
  unlist(
    lapply(basic_vars, function(v) {
      
      x <- d[[v]]
      
      c(
        setNames(mean(x), paste0(v, "_mean")),
        setNames(sd(x), paste0(v, "_SD"))
      )
    })
  )
}

n_obs <- nrow(av)

observed <- boot_basic(av)

boot_results <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(observed),
  dimnames = list(NULL, names(observed))
)

for (b in seq_len(B)) {
  
  idx <- sample.int(
    n_obs,
    size = n_obs,
    replace = TRUE
  )
  
  boot_results[b, ] <- boot_basic(av[idx, , drop = FALSE])
}


bootstrap_summary <- do.call(
  rbind,
  lapply(names(observed), function(stat_name) {
    
    x <- boot_results[, stat_name]
    x <- x[is.finite(x)]
    
    ci <- quantile(
      x,
      probs = c(0.025, 0.975),
      names = FALSE
    )
    
    data.frame(
      statistic = stat_name,
      estimate = unname(observed[stat_name]),
      bootstrap_mean = mean(x),
      lower_95_CI = ci[1],
      upper_95_CI = ci[2],
      bootstrap_SD = sd(x),
      successful_replicates = length(x)
    )
  })
)

rownames(bootstrap_summary) <- NULL

# ------------------------------------------------------------
# 8. BUILD ONE CONSOLIDATED RESULTS TABLE
# ------------------------------------------------------------

make_rows <- function(
    analysis,
    metric,
    variable = NA_character_,
    variable_2 = NA_character_,
    n = NA_integer_,
    estimate = NA_real_,
    bootstrap_mean = NA_real_,
    lower_95_CI = NA_real_,
    upper_95_CI = NA_real_,
    bootstrap_SD = NA_real_,
    pearson = NA_real_,
    spearman = NA_real_,
    mean_rank_change = NA_real_,
    median_rank_change = NA_real_,
    top10_overlap = NA_real_,
    top20_overlap = NA_real_) {
  
  data.frame(
    analysis = analysis,
    metric = metric,
    variable = variable,
    variable_2 = variable_2,
    n = n,
    estimate = estimate,
    bootstrap_mean = bootstrap_mean,
    lower_95_CI = lower_95_CI,
    upper_95_CI = upper_95_CI,
    bootstrap_SD = bootstrap_SD,
    pearson = pearson,
    spearman = spearman,
    mean_rank_change = mean_rank_change,
    median_rank_change = median_rank_change,
    top10_overlap = top10_overlap,
    top20_overlap = top20_overlap,
    stringsAsFactors = FALSE
  )
}

all_results <- list()

# A. Basic-variable descriptive statistics

for (i in seq_len(nrow(basic_summary))) {
  
  v <- basic_summary$variable[i]
  
  for (stat in c("mean", "SD", "median", "minimum", "maximum")) {
    
    all_results[[length(all_results) + 1]] <- make_rows(
      analysis = "Basic-variable descriptives",
      metric = stat,
      variable = v,
      n = basic_summary$n[i],
      estimate = basic_summary[[stat]][i]
    )
  }
}

# B. Component-contribution descriptive statistics

for (i in seq_len(nrow(contribution_summary))) {
  
  v <- contribution_summary$component[i]
  
  for (stat in c("mean", "SD", "median", "minimum", "maximum")) {
    
    all_results[[length(all_results) + 1]] <- make_rows(
      analysis = "AviCRED component contributions",
      metric = stat,
      variable = v,
      n = contribution_summary$n[i],
      estimate = contribution_summary[[stat]][i]
    )
  }
}

# C. Pearson and Spearman correlations

for (method in c("Pearson", "Spearman")) {
  
  correlation_matrix <- if (method == "Pearson") {
    basic_cor_pearson
  } else {
    basic_cor_spearman
  }
  
  for (i in seq_along(basic_vars)) {
    for (j in seq_along(basic_vars)) {
      
      if (i < j) {
        
        all_results[[length(all_results) + 1]] <- make_rows(
          analysis = "Basic-variable correlations",
          metric = method,
          variable = basic_vars[i],
          variable_2 = basic_vars[j],
          n = n_obs,
          estimate = correlation_matrix[i, j]
        )
      }
    }
  }
}

# D. Original AviCRED component sensitivity

for (i in seq_len(nrow(rank_results))) {
  
  x <- rank_results[i, ]
  
  all_results[[length(all_results) + 1]] <- make_rows(
    analysis = "Original AviCRED sensitivity",
    metric = paste("Removed", x$component_removed),
    variable = x$component_removed,
    n = n_obs,
    pearson = x$Pearson,
    spearman = x$Spearman,
    mean_rank_change = x$mean_rank_change,
    median_rank_change = x$median_rank_change,
    top10_overlap = x$top10_overlap,
    top20_overlap = x$top20_overlap
  )
}

# E. Bootstrap descriptive statistics and confidence intervals

for (i in seq_len(nrow(bootstrap_summary))) {
  
  x <- bootstrap_summary[i, ]
  
  all_results[[length(all_results) + 1]] <- make_rows(
    analysis = "Bootstrap descriptive statistics",
    metric = x$statistic,
    variable = sub("_(mean|SD)$", "", x$statistic),
    n = n_obs,
    estimate = x$estimate,
    bootstrap_mean = x$bootstrap_mean,
    lower_95_CI = x$lower_95_CI,
    upper_95_CI = x$upper_95_CI,
    bootstrap_SD = x$bootstrap_SD
  )
}

# F. AviCRED formula verification

all_results[[length(all_results) + 1]] <- make_rows(
  analysis = "Formula verification",
  metric = "Maximum absolute reconstruction discrepancy",
  n = n_obs,
  estimate = max_discrepancy
)

# Combine results into one table

combined_results <- do.call(rbind, all_results)

rownames(combined_results) <- NULL

# ------------------------------------------------------------
# 9. SAVE FINAL OUTPUTS
# ------------------------------------------------------------

# Main consolidated results file
write.csv(
  combined_results,
  file.path(
    output_dir,
    "AviCRED_ALL_ANALYSES_SUMMARY.csv"
  ),
  row.names = FALSE,
  na = ""
)

# Observer-level scores and component contributions
write.csv(
  av,
  file.path(
    output_dir,
    "AviCRED_observer_analysis.csv"
  ),
  row.names = FALSE,
  na = ""
)

# Bootstrap replicates retained for reproducibility
write.csv(
  boot_results,
  file.path(
    output_dir,
    "component_variable_bootstrap_replicates.csv"
  ),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 10. DISPLAY RESULTS
# ------------------------------------------------------------

cat("\n============================================\n")
cat("BASIC VARIABLE DESCRIPTIVE STATISTICS\n")
cat("============================================\n")
print(basic_summary, row.names = FALSE, digits = 5)

cat("\n============================================\n")
cat("PEARSON CORRELATIONS\n")
cat("============================================\n")
print(round(basic_cor_pearson, 3))

cat("\n============================================\n")
cat("SPEARMAN CORRELATIONS\n")
cat("============================================\n")
print(round(basic_cor_spearman, 3))

cat("\n============================================\n")
cat("AVICRED COMPONENT CONTRIBUTIONS\n")
cat("============================================\n")
print(contribution_summary, row.names = FALSE, digits = 5)

cat("\n============================================\n")
cat("ORI-CAS-SBI RANK SENSITIVITY\n")
cat("============================================\n")
print(rank_results, row.names = FALSE, digits = 5)

cat("\n============================================\n")
cat("BOOTSTRAP MEANS, SDs AND 95% CIs\n")
cat("============================================\n")
print(bootstrap_summary, row.names = FALSE, digits = 5)

cat("\n============================================\n")
cat("CONSOLIDATED RESULTS\n")
cat("============================================\n")
cat("Total result rows:", nrow(combined_results), "\n")

# ============================================================
# AVICRED: BASIC-COMPONENT ABLATION 
# ============================================================

# 1. CHECK REQUIRED COLUMNS
required_cols <- c(
  "years_active",
  "n_checklists",
  "total_species",
  "ORI",
  "ORI_bounded",
  "mean_CAS",
  "consistency_score",
  "weighted_completeness_score",
  "SBI",
  "AviCRED"
)

missing_cols <- setdiff(required_cols, names(observer_final))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns in observer_final: ",
    paste(missing_cols, collapse = ", ")
  )
}

# 2. PREPARE ANALYSIS DATA
# Use identical observers across all ablations.
dat <- observer_final %>%
  select(all_of(required_cols)) %>%
  filter(if_all(everything(), ~ is.finite(.x)))

cat("Observers included:", nrow(dat), "\n")

# 3. VERIFY THE ORIGINAL COMPOSITE
avi_reconstructed <- (
  dat$ORI_bounded +
    dat$mean_CAS +
    dat$SBI
) / 3

cat(
  "Maximum AviCRED reconstruction difference:",
  max(abs(dat$AviCRED - avi_reconstructed)),
  "\n"
)

# Verify SBI
SBI_reconstructed <- (
  dat$consistency_score +
    dat$weighted_completeness_score
) / 2

cat(
  "Maximum SBI reconstruction difference:",
  max(abs(dat$SBI - SBI_reconstructed)),
  "\n"
)

if (max(abs(dat$SBI - SBI_reconstructed)) > 1e-6) {
  stop("Stored SBI does not match its assumed formula.")
}

# 4. STANDARDISE ORI INPUTS
z_years <- as.numeric(scale(dat$years_active))

z_checklists <- as.numeric(
  scale(log1p(dat$n_checklists))
)

z_species <- as.numeric(
  scale(log1p(dat$total_species))
)

# Verify ORI_bounded
ORI_reconstructed <- pnorm(dat$ORI)

cat(
  "Maximum ORI_bounded reconstruction difference:",
  max(abs(dat$ORI_bounded - ORI_reconstructed)),
  "\n"
)

if (max(abs(dat$ORI_bounded - ORI_reconstructed)) > 1e-6) {
  stop("ORI_bounded does not match pnorm(ORI).")
}

# 5. ORIGINAL SCORES AND RANKS
original_score <- dat$AviCRED

original_rank <- rank(
  -original_score,
  ties.method = "average"
)

original_order <- order(
  -original_score,
  seq_along(original_score)
)

original_top10 <- original_order[seq_len(10)]
original_top20 <- original_order[seq_len(20)]

# 6. ABLATION EVALUATION FUNCTION WITH P-VALUES

evaluate_ablation <- function(
    modified_score,
    component,
    variable_removed) {
  
  modified_rank <- rank(
    -modified_score,
    ties.method = "average"
  )
  
  modified_order <- order(
    -modified_score,
    seq_along(modified_score)
  )
  
  modified_top10 <- modified_order[seq_len(10)]
  modified_top20 <- modified_order[seq_len(20)]
  
  rank_change <- abs(modified_rank - original_rank)
  
  # Pearson correlation test
  pearson_test <- cor.test(
    original_score,
    modified_score,
    method = "pearson"
  )
  
  # Spearman correlation test
  spearman_test <- cor.test(
    original_score,
    modified_score,
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    component = component,
    variable_removed = variable_removed,
    N = length(original_score),
    
    Pearson_r = unname(pearson_test$estimate),
    Pearson_p = pearson_test$p.value,
    
    Spearman_rho = unname(spearman_test$estimate),
    Spearman_p = spearman_test$p.value,
    
    MAE = mean(
      abs(original_score - modified_score)
    ),
    
    RMSE = sqrt(
      mean((original_score - modified_score)^2)
    ),
    
    mean_rank_change = mean(rank_change),
    
    median_rank_change = median(rank_change),
    
    mean_rank_change_percent =
      100 * mean(rank_change) / length(original_score),
    
    top10_overlap = length(
      intersect(original_top10, modified_top10)
    ) / 10,
    
    top20_overlap = length(
      intersect(original_top20, modified_top20)
    ) / 20
  )
}

# 7. RUN ALL SIX VARIABLE-LEVEL ABLATIONS
results <- list()

# ------------------------------------------------------------
# ORI ABLATIONS
# Reconstruct ORI from the mean of the two remaining
# standardised inputs, then apply the standard normal CDF.
# ------------------------------------------------------------

# 7.1 Remove years_active
ORI_no_years <- pnorm(
  (z_checklists + z_species) / 2
)

results[["years_active"]] <- evaluate_ablation(
  (ORI_no_years + dat$mean_CAS + dat$SBI) / 3,
  "ORI",
  "years_active"
)

# 7.2 Remove n_checklists
ORI_no_checklists <- pnorm(
  (z_years + z_species) / 2
)

results[["n_checklists"]] <- evaluate_ablation(
  (ORI_no_checklists + dat$mean_CAS + dat$SBI) / 3,
  "ORI",
  "n_checklists"
)

# 7.3 Remove total_species
ORI_no_species <- pnorm(
  (z_years + z_checklists) / 2
)

results[["total_species"]] <- evaluate_ablation(
  (ORI_no_species + dat$mean_CAS + dat$SBI) / 3,
  "ORI",
  "total_species"
)

# ------------------------------------------------------------
# CAS ABLATION
# Remove mean_CAS; retain the original denominator of 3.
# ------------------------------------------------------------

# 7.4 Remove mean_CAS
results[["mean_CAS"]] <- evaluate_ablation(
  (dat$ORI_bounded + 0 + dat$SBI) / 3,
  "CAS",
  "mean_CAS"
)

# ------------------------------------------------------------
# SBI ABLATIONS
# This reproduces your existing method: substitute the
# remaining SBI input directly, without averaging it with zero.
# ------------------------------------------------------------

# 7.5 Remove consistency_score
SBI_no_consistency <- dat$weighted_completeness_score

results[["consistency_score"]] <- evaluate_ablation(
  (dat$ORI_bounded + dat$mean_CAS +
     SBI_no_consistency) / 3,
  "SBI",
  "consistency_score"
)

# 7.6 Remove weighted_completeness_score
SBI_no_completeness <- dat$consistency_score

results[["weighted_completeness_score"]] <- evaluate_ablation(
  (dat$ORI_bounded + dat$mean_CAS +
     SBI_no_completeness) / 3,
  "SBI",
  "weighted_completeness_score"
)

# 8. COMBINE RESULTS
basic_ablation_results <- bind_rows(results) %>%
  mutate(
    top10_agreement_percent = 100 * top10_overlap,
    top20_agreement_percent = 100 * top20_overlap
  )

# 9. DISPLAY RESULTS
cat("\nBASIC-VARIABLE ABLATION SENSITIVITY RESULTS\n")
cat("====================================================\n")

print(
  basic_ablation_results,
  row.names = FALSE,
  digits = 5
)

# 10. SAVE RESULTS
write.csv(
  basic_ablation_results,
  "FILE PATH/AviCRED_basic_variable_ablation.csv",
  row.names = FALSE
)

# ============================================================
# AVICRED CHECKLIST-HISTORY THRESHOLDING ANALYSIS
# INTEGRATED SCORE: ORI + MEAN CAS + SBI
#
# Method:
#   - Calculate checklist CAS from State x Season contextual
#     medians and MADs using the established CAS formula
#   - Randomly sample n checklists per eligible observer
#   - Recalculate sampled mean CAS and integrated AviCRED
#   - Compare sampled AviCRED against full-history AviCRED
#   - Repeat each threshold 30 times
#   - Calculate MAE, RMSE, Spearman rho, four-class agreement,
#     quadratic weighted kappa, retention and diminishing returns
# ============================================================

AVICRED_FILE <- "FILE PATH/Avicred new.csv"

CHECKLIST_FILE <- file.choose("FILE PATH/checlist_master.csv")

OUT_FILE <- "FILE PATH/AviCRED_threshold_all_results.csv"

N_REPLICATES <- 30L

THRESHOLDS <- c(
  2L, 3L, 5L, 10L, 20L, 30L, 50L, 75L, 100L
)

RANDOM_SEED <- 20261011L

# CAS robust-Z multiplier and decay constant.
# Robust Z = 0.6745 * (x - median) / MAD
# Score = 1 for robust Z <= 0; exp(-Z/3.5) otherwise.
ROBUST_Z_MULTIPLIER <- 0.6745
CAS_DECAY_CONSTANT <- 3.5

set.seed(RANDOM_SEED)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table first: install.packages('data.table')")
}

library(data.table)

if (!file.exists(AVICRED_FILE)) {
  stop(
    "AviCRED CSV not found: ",
    AVICRED_FILE,
    "\nCorrect AVICRED_FILE before continuing."
  )
}

if (!file.exists(CHECKLIST_FILE)) {
  stop("Checklist-master CSV not found: ", CHECKLIST_FILE)
}

avicred_reference <- data.table::fread(
  AVICRED_FILE,
  na.strings = c("", "NA", "NaN", "NULL")
)

checklist_data <- data.table::fread(
  CHECKLIST_FILE,
  na.strings = c("", "NA", "NaN", "NULL")
)

cat("\nAviCRED reference rows:", nrow(avicred_reference), "\n")
cat("Checklist-master rows:", nrow(checklist_data), "\n")

cat("\nAviCRED CSV columns:\n")
print(names(avicred_reference))

cat("\nChecklist-master CSV columns:\n")
print(names(checklist_data))

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

normalise_names <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

find_column <- function(data, candidates, description,
                        required = TRUE) {
  
  actual_names <- names(data)
  actual_normalised <- normalise_names(actual_names)
  candidate_normalised <- normalise_names(candidates)
  
  idx <- match(candidate_normalised, actual_normalised)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) > 0L) {
    return(actual_names[idx[1L]])
  }
  
  if (required) {
    stop(
      "\nCould not identify the required column: ",
      description,
      "\nTried these names: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(actual_names, collapse = ", "),
      "\nNo substitute column has been guessed."
    )
  }
  
  return(NULL)
}

to_numeric_safe <- function(x, column_name) {
  
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  original <- as.character(x)
  
  result <- suppressWarnings(as.numeric(original))
  
  failed <- !is.na(original) &
    nzchar(trimws(original)) &
    is.na(result)
  
  if (any(failed)) {
    warning(
      sum(failed),
      " values in column '",
      column_name,
      "' could not be converted to numeric."
    )
  }
  
  result
}

robust_z <- function(x, centre, mad_value) {
  
  if (!is.finite(centre) || !is.finite(mad_value) ||
      mad_value <= 0) {
    return(rep(0, length(x)))
  }
  
  ROBUST_Z_MULTIPLIER * (x - centre) / mad_value
}


cas_component_score <- function(z) {
  
  out <- rep(NA_real_, length(z))
  
  valid <- is.finite(z)
  
  out[valid] <- 1
  
  positive <- valid & z > 0
  
  out[positive] <- exp(
    -z[positive] / CAS_DECAY_CONSTANT
  )
  
  out
}

spearman_safe <- function(x, y) {
  
  keep <- is.finite(x) & is.finite(y)
  
  if (sum(keep) < 3L ||
      length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  
  suppressWarnings(
    cor(x[keep], y[keep], method = "spearman")
  )
}


mae_safe <- function(observed, estimated) {
  
  keep <- is.finite(observed) & is.finite(estimated)
  
  if (!any(keep)) return(NA_real_)
  
  mean(abs(estimated[keep] - observed[keep]))
}


rmse_safe <- function(observed, estimated) {
  
  keep <- is.finite(observed) & is.finite(estimated)
  
  if (!any(keep)) return(NA_real_)
  
  sqrt(mean((estimated[keep] - observed[keep])^2))
}


# Quadratic weighted Cohen's kappa.
# Categories must be integers from 1 to K.

quadratic_weighted_kappa <- function(reference_class,
                                     estimated_class,
                                     K = 4L) {
  
  keep <- !is.na(reference_class) &
    !is.na(estimated_class) &
    reference_class >= 1L &
    reference_class <= K &
    estimated_class >= 1L &
    estimated_class <= K
  
  reference_class <- as.integer(reference_class[keep])
  estimated_class <- as.integer(estimated_class[keep])
  
  n <- length(reference_class)
  
  if (n < 2L) return(NA_real_)
  
  observed <- table(
    factor(reference_class, levels = seq_len(K)),
    factor(estimated_class, levels = seq_len(K))
  )
  
  observed <- as.matrix(observed) / n
  
  ref_marginal <- rowSums(observed)
  est_marginal <- colSums(observed)
  
  expected <- outer(ref_marginal, est_marginal)
  
  weights <- outer(
    seq_len(K),
    seq_len(K),
    function(i, j) ((i - j)^2) / ((K - 1)^2)
  )
  
  observed_disagreement <- sum(weights * observed)
  expected_disagreement <- sum(weights * expected)
  
  if (!is.finite(expected_disagreement) ||
      expected_disagreement == 0) {
    return(NA_real_)
  }
  
  1 - observed_disagreement / expected_disagreement
}


classify_by_fixed_breaks <- function(x, breaks) {
  
  out <- rep(NA_integer_, length(x))
  
  valid <- is.finite(x)
  
  out[valid] <- as.integer(
    cut(
      x[valid],
      breaks = breaks,
      include.lowest = TRUE,
      labels = FALSE
    )
  )
  
  out
}

# ============================================================
# 2. IDENTIFY ACTUAL CSV COLUMNS
# ============================================================

# ---- AviCRED reference columns ----

REF_OBSERVER_COL <- find_column(
  avicred_reference,
  c("observer_id", "OBSERVER ID", "observerid"),
  "observer ID in AviCRED CSV"
)

REF_ORI_COL <- find_column(
  avicred_reference,
  c("ORI_bounded", "ori_bounded", "ORI bounded"),
  "bounded ORI score in AviCRED CSV"
)

REF_CAS_COL <- find_column(
  avicred_reference,
  c("mean_CAS", "mean CAS", "mean_cas"),
  "full-history mean CAS in AviCRED CSV"
)

REF_SBI_COL <- find_column(
  avicred_reference,
  c("SBI", "sbi"),
  "SBI score in AviCRED CSV"
)

REF_AVICRED_COL <- find_column(
  avicred_reference,
  c("AviCRED", "avicred", "AviCRED_score"),
  "full-history integrated AviCRED score in AviCRED CSV"
)

# ---- Checklist-master columns ----

CHECK_OBSERVER_COL <- find_column(
  checklist_data,
  c("observer_id", "OBSERVER ID", "observerid"),
  "observer ID in checklist-master CSV"
)

CHECK_ID_COL <- find_column(
  checklist_data,
  c("checklist_id", "SAMPLING EVENT IDENTIFIER",
    "sampling_event_identifier", "event_id"),
  "unique checklist ID in checklist-master CSV"
)

CHECK_RICHNESS_COL <- find_column(
  checklist_data,
  c("richness", "species_richness", "n_species"),
  "checklist richness in checklist-master CSV"
)

CHECK_ABUNDANCE_COL <- find_column(
  checklist_data,
  c("abundance", "abundance_per_checklist",
    "total_abundance", "checklist_abundance"),
  "total checklist abundance in checklist-master CSV"
)

CHECK_DATE_COL <- find_column(
  checklist_data,
  c("obs_date", "observation_date", "date",
    "observationdate"),
  "checklist observation date in checklist-master CSV"
)

CHECK_STATE_COL <- find_column(
  checklist_data,
  c("state", "state_province", "stateprovince",
    "admin1", "state_name"),
  "state in checklist-master CSV"
)


cat("\nColumn mapping used:\n")

print(
  data.frame(
    input = c(
      "AviCRED observer ID",
      "ORI_bounded",
      "mean_CAS",
      "SBI",
      "AviCRED",
      "Checklist observer ID",
      "Checklist ID",
      "Richness",
      "Abundance",
      "Date",
      "State"
    ),
    column = c(
      REF_OBSERVER_COL,
      REF_ORI_COL,
      REF_CAS_COL,
      REF_SBI_COL,
      REF_AVICRED_COL,
      CHECK_OBSERVER_COL,
      CHECK_ID_COL,
      CHECK_RICHNESS_COL,
      CHECK_ABUNDANCE_COL,
      CHECK_DATE_COL,
      CHECK_STATE_COL
    )
  )
)

# ============================================================
# 3. PREPARE FULL-HISTORY AVICRED REFERENCE
# ============================================================

reference <- data.table(
  observer_id = as.character(
    avicred_reference[[REF_OBSERVER_COL]]
  ),
  ORI_bounded = to_numeric_safe(
    avicred_reference[[REF_ORI_COL]],
    REF_ORI_COL
  ),
  mean_CAS_reference = to_numeric_safe(
    avicred_reference[[REF_CAS_COL]],
    REF_CAS_COL
  ),
  SBI = to_numeric_safe(
    avicred_reference[[REF_SBI_COL]],
    REF_SBI_COL
  ),
  AviCRED_reference = to_numeric_safe(
    avicred_reference[[REF_AVICRED_COL]],
    REF_AVICRED_COL
  )
)

reference <- reference[
  !is.na(observer_id) &
    nzchar(observer_id)
]

# Resolve duplicate observer rows only if they are exact duplicates
# across the required reference fields. Conflicting duplicates stop
# the analysis rather than being silently averaged.

reference_counts <- reference[
  ,
  .(
    n_rows = .N,
    n_ORI = uniqueN(ORI_bounded),
    n_CAS = uniqueN(mean_CAS_reference),
    n_SBI = uniqueN(SBI),
    n_AviCRED = uniqueN(AviCRED_reference)
  ),
  by = observer_id
]

conflicting_reference <- reference_counts[
  n_ORI > 1L |
    n_CAS > 1L |
    n_SBI > 1L |
    n_AviCRED > 1L
]

if (nrow(conflicting_reference) > 0L) {
  stop(
    "\nThe AviCRED CSV contains conflicting rows for ",
    nrow(conflicting_reference),
    " observer IDs. Resolve duplicate observer records before ",
    "running the threshold analysis."
  )
}

reference <- unique(reference, by = "observer_id")

reference <- reference[
  is.finite(ORI_bounded) &
    is.finite(mean_CAS_reference) &
    is.finite(SBI) &
    is.finite(AviCRED_reference)
]

if (nrow(reference) == 0L) {
  stop("No complete observer rows remain in the AviCRED reference.")
}

# Verify that the integrated score in the CSV is consistent with
# the current three-component definition:
# AviCRED = mean(ORI_bounded, full-history mean CAS, SBI).

reference[, AviCRED_recalculated := rowMeans(
  cbind(
    ORI_bounded,
    mean_CAS_reference,
    SBI
  )
)]

reference[, absolute_reference_difference := abs(
  AviCRED_reference - AviCRED_recalculated
)]

max_reference_difference <- max(
  reference$absolute_reference_difference,
  na.rm = TRUE
)

cat(
  "\nMaximum difference between reference AviCRED and ",
  "mean(ORI_bounded, mean_CAS, SBI): ",
  format(max_reference_difference, scientific = TRUE),
  "\n",
  sep = ""
)

if (max_reference_difference > 1e-5) {
  stop(
    "\nThe AviCRED column does not match the current ",
    "ORI_bounded + mean_CAS + SBI equal-weight formula. ",
    "Do not continue until the reference score definition ",
    "is reconciled."
  )
}

reference[, AviCRED_reference := AviCRED_recalculated]

reference[, absolute_reference_difference := NULL]
reference[, AviCRED_recalculated := NULL]

# ============================================================
# 4. PREPARE CHECKLIST-LEVEL DATA
# ============================================================

checklists <- data.table(
  observer_id = as.character(
    checklist_data[[CHECK_OBSERVER_COL]]
  ),
  checklist_id = as.character(
    checklist_data[[CHECK_ID_COL]]
  ),
  richness = to_numeric_safe(
    checklist_data[[CHECK_RICHNESS_COL]],
    CHECK_RICHNESS_COL
  ),
  abundance = to_numeric_safe(
    checklist_data[[CHECK_ABUNDANCE_COL]],
    CHECK_ABUNDANCE_COL
  ),
  obs_date = as.Date(
    as.character(checklist_data[[CHECK_DATE_COL]])
  ),
  state = trimws(
    as.character(checklist_data[[CHECK_STATE_COL]])
  )
)

# Remove records missing identifiers or the fields required for CAS.

checklists <- checklists[
  !is.na(observer_id) &
    nzchar(observer_id) &
    !is.na(checklist_id) &
    nzchar(checklist_id) &
    is.finite(richness) &
    is.finite(abundance) &
    !is.na(obs_date) &
    !is.na(state) &
    nzchar(state)
]

# The sampling unit must be one checklist per observer.
# Exact duplicate checklist-observer rows are collapsed.
# Conflicting duplicates are not allowed.

duplicate_check <- checklists[
  ,
  .(
    n_rows = .N,
    n_richness = uniqueN(richness),
    n_abundance = uniqueN(abundance),
    n_date = uniqueN(obs_date),
    n_state = uniqueN(state)
  ),
  by = .(observer_id, checklist_id)
]

conflicting_checklists <- duplicate_check[
  n_richness > 1L |
    n_abundance > 1L |
    n_date > 1L |
    n_state > 1L
]

if (nrow(conflicting_checklists) > 0L) {
  stop(
    "\nThere are ",
    nrow(conflicting_checklists),
    " observer-checklist pairs with conflicting values. ",
    "Resolve them before sampling."
  )
}

checklists <- unique(
  checklists,
  by = c("observer_id", "checklist_id")
)

# Derive four seasons:
# Winter = Dec-Feb
# Summer = Mar-May
# Monsoon = Jun-Sep
# Post-Monsoon = Oct-Nov

month_number <- as.integer(format(checklists$obs_date, "%m"))

checklists[, season := fifelse(
  month_number %in% c(12L, 1L, 2L),
  "Winter",
  fifelse(
    month_number %in% 3L:5L,
    "Summer",
    fifelse(
      month_number %in% 6L:9L,
      "Monsoon",
      "Post-Monsoon"
    )
  )
)]

# CAS uses log-transformed abundance to reduce right skewness.

checklists[, log_abundance := log1p(
  pmax(abundance, 0)
)]

cat("\nUsable checklist records:", nrow(checklists), "\n")
cat(
  "Unique observers in checklist data:",
  uniqueN(checklists$observer_id),
  "\n"
)
cat(
  "Unique observer-checklist pairs:",
  nrow(checklists),
  "\n"
)

# ============================================================
# 5. CALCULATE CONTEXTUAL CAS BASELINES
# ============================================================

# Baselines are calculated once from the complete checklist dataset.
# The same baseline is then held fixed for all threshold replicates.
# This isolates the effect of reducing observer checklist history
# rather than changing the ecological reference distribution.

context_baselines <- checklists[
  ,
  .(
    richness_median = median(richness, na.rm = TRUE),
    richness_mad = mad(
      richness,
      center = median(richness, na.rm = TRUE),
      constant = 1,
      na.rm = TRUE
    ),
    abundance_median = median(log_abundance, na.rm = TRUE),
    abundance_mad = mad(
      log_abundance,
      center = median(log_abundance, na.rm = TRUE),
      constant = 1,
      na.rm = TRUE
    ),
    n_context_checklists = .N
  ),
  by = .(state, season)
]

# Fallback for contexts with zero MAD:
# use IQR / 1.349; if that is also zero, use SD.
# If all dispersion measures are zero, use a small-context fallback
# that assigns zero robust Z-scores rather than dividing by zero.

context_baselines[
  ,
  richness_iqr := NA_real_
]

context_baselines[
  ,
  abundance_iqr := NA_real_
]

for (i in seq_len(nrow(context_baselines))) {
  
  group_state <- context_baselines$state[i]
  group_season <- context_baselines$season[i]
  
  idx <- checklists$state == group_state &
    checklists$season == group_season
  
  x_richness <- checklists$richness[idx]
  x_abundance <- checklists$log_abundance[idx]
  
  if (!is.finite(context_baselines$richness_mad[i]) ||
      context_baselines$richness_mad[i] <= 0) {
    
    q <- IQR(x_richness, na.rm = TRUE) / 1.349
    
    if (is.finite(q) && q > 0) {
      context_baselines$richness_mad[i] <- q
    } else {
      s <- sd(x_richness, na.rm = TRUE)
      
      if (is.finite(s) && s > 0) {
        context_baselines$richness_mad[i] <- s
      } else {
        context_baselines$richness_mad[i] <- NA_real_
      }
    }
  }
  
  if (!is.finite(context_baselines$abundance_mad[i]) ||
      context_baselines$abundance_mad[i] <= 0) {
    
    q <- IQR(x_abundance, na.rm = TRUE) / 1.349
    
    if (is.finite(q) && q > 0) {
      context_baselines$abundance_mad[i] <- q
    } else {
      s <- sd(x_abundance, na.rm = TRUE)
      
      if (is.finite(s) && s > 0) {
        context_baselines$abundance_mad[i] <- s
      } else {
        context_baselines$abundance_mad[i] <- NA_real_
      }
    }
  }
}

checklists <- merge(
  checklists,
  context_baselines[
    ,
    .(
      state,
      season,
      richness_median,
      richness_mad,
      abundance_median,
      abundance_mad
    )
  ],
  by = c("state", "season"),
  all.x = TRUE,
  sort = FALSE
)

checklists[
  ,
  richness_z := ifelse(
    is.finite(richness_mad) & richness_mad > 0,
    ROBUST_Z_MULTIPLIER *
      (richness - richness_median) / richness_mad,
    0
  )
]

checklists[
  ,
  abundance_z := ifelse(
    is.finite(abundance_mad) & abundance_mad > 0,
    ROBUST_Z_MULTIPLIER *
      (log_abundance - abundance_median) / abundance_mad,
    0
  )
]

checklists[
  ,
  richness_score := cas_component_score(richness_z)
]

checklists[
  ,
  abundance_score := cas_component_score(abundance_z)
]

checklists[
  ,
  CAS := rowMeans(
    cbind(richness_score, abundance_score),
    na.rm = FALSE
  )
]

checklists <- checklists[
  is.finite(CAS)
]

if (nrow(checklists) == 0L) {
  stop("No finite checklist-level CAS scores could be calculated.")
}

cat("\nChecklist-level CAS calculated.\n")
cat("Checklist rows with finite CAS:", nrow(checklists), "\n")


# ============================================================
# 6. CHECK FULL-HISTORY CAS AGAINST THE AVICRED CSV
# ============================================================

# This check determines whether the newly calculated checklist CAS
# agrees with the observer-level mean_CAS in the reference CSV.
# A discrepancy is reported because it can indicate differences in
# filtering, reference groups, abundance definition, or CAS history.

full_history_cas <- checklists[
  ,
  .(
    mean_CAS_calculated = mean(CAS, na.rm = TRUE),
    n_checklists_available = .N
  ),
  by = observer_id
]

reference_check <- merge(
  reference[
    ,
    .(
      observer_id,
      mean_CAS_reference,
      AviCRED_reference
    )
  ],
  full_history_cas,
  by = "observer_id",
  all = FALSE
)

reference_check <- reference_check[
  is.finite(mean_CAS_reference) &
    is.finite(mean_CAS_calculated)
]

if (nrow(reference_check) == 0L) {
  stop(
    "No overlapping observers between the AviCRED CSV ",
    "and checklist-master CSV."
  )
}

CAS_reference_MAE <- mae_safe(
  reference_check$mean_CAS_reference,
  reference_check$mean_CAS_calculated
)

CAS_reference_RMSE <- rmse_safe(
  reference_check$mean_CAS_reference,
  reference_check$mean_CAS_calculated
)

CAS_reference_rho <- spearman_safe(
  reference_check$mean_CAS_reference,
  reference_check$mean_CAS_calculated
)

cat("\nFull-history CAS agreement diagnostics:\n")

print(
  data.frame(
    overlapping_observers = nrow(reference_check),
    CAS_MAE = CAS_reference_MAE,
    CAS_RMSE = CAS_reference_RMSE,
    CAS_Spearman_rho = CAS_reference_rho
  )
)

# A very low correlation or large discrepancy means the CSV and
# checklist data do not reproduce the same CAS reference. Stop rather
# than report a misleading threshold recommendation.

if (!is.finite(CAS_reference_rho) ||
    CAS_reference_rho < 0.90) {
  
  stop(
    "\nCalculated full-history CAS has Spearman rho < 0.90 ",
    "against mean_CAS in the AviCRED reference CSV. ",
    "The checklist data and reference do not reproduce the same ",
    "CAS definition or checklist universe. Resolve this mismatch ",
    "before interpreting threshold results."
  )
}

# ============================================================
# 7. BUILD FULL-HISTORY OBSERVER REFERENCE
# ============================================================

# Only observers with both:
#   a) a valid integrated AviCRED reference score, and
#   b) usable checklist-level CAS history
# can enter the threshold experiment.

observer_history <- checklists[
  ,
  .(
    n_checklists = .N
  ),
  by = observer_id
]

analysis_data <- merge(
  reference,
  observer_history,
  by = "observer_id",
  all = FALSE
)

analysis_data <- merge(
  analysis_data,
  full_history_cas[
    ,
    .(
      observer_id,
      mean_CAS_calculated,
      n_checklists_available
    )
  ],
  by = "observer_id",
  all = FALSE
)

# Retain the CSV's original mean CAS as the full-history reference.
# The recalculated mean CAS is used for sampled-history estimation.
analysis_data <- analysis_data[
  is.finite(ORI_bounded) &
    is.finite(SBI) &
    is.finite(AviCRED_reference) &
    is.finite(mean_CAS_reference) &
    n_checklists >= min(THRESHOLDS)
]

if (nrow(analysis_data) < 10L) {
  stop(
    "Fewer than 10 observers meet the minimum threshold ",
    "and have complete integrated AviCRED components."
  )
}

cat(
  "\nObservers eligible for at least one threshold:",
  nrow(analysis_data),
  "\n"
)

# ============================================================
# 8. DEFINE FIXED FOUR-CLASS REFERENCE BOUNDARIES
# ============================================================

# Use quartiles of the full-history AviCRED reference scores.
# Apply these same boundaries to every threshold and replicate.
# This avoids changing class definitions between replicates.

class_breaks <- as.numeric(
  quantile(
    analysis_data$AviCRED_reference,
    probs = c(0, 0.25, 0.50, 0.75, 1),
    na.rm = TRUE,
    type = 7,
    names = FALSE
  )
)

# If ties cause repeated quantile breaks, stop because four distinct
# ordinal classes cannot be formed reliably using these cutpoints.

if (any(diff(class_breaks) <= 0)) {
  stop(
    "The reference AviCRED quartiles contain tied boundaries. ",
    "A four-class comparison using fixed quartile cutpoints ",
    "is not identifiable. Define and justify fixed class ",
    "boundaries before continuing."
  )
}

analysis_data[
  ,
  reference_class := classify_by_fixed_breaks(
    AviCRED_reference,
    class_breaks
  )
]

# ============================================================
# 9. RANDOM-SUBSAMPLING THRESHOLD ANALYSIS
# ============================================================

replicate_results <- list()
result_counter <- 0L

for (threshold_n in THRESHOLDS) {
  
  cat(
    "\nRunning threshold:",
    threshold_n,
    "checklists per observer\n"
  )
  
  eligible_ids <- analysis_data[
    n_checklists >= threshold_n,
    observer_id
  ]
  
  n_eligible <- length(eligible_ids)
  
  if (n_eligible < 10L) {
    warning(
      "Only ",
      n_eligible,
      " observers meet threshold ",
      threshold_n,
      "."
    )
  }
  
  eligible_reference <- analysis_data[
    observer_id %in% eligible_ids
  ]
  
  # Reference metrics are calculated only for observers eligible at
  # the current threshold. This ensures a fair comparison of scores
  # for the same observer set within each replicate.
  
  for (replicate_id in seq_len(N_REPLICATES)) {
    
    set.seed(
      RANDOM_SEED +
        threshold_n * 1000L +
        replicate_id
    )
    
    sampled_rows <- vector(
      mode = "list",
      length = n_eligible
    )
    
    names(sampled_rows) <- eligible_ids
    
    # Sample exactly n checklists without replacement per observer.
    
    for (j in seq_along(eligible_ids)) {
      
      observer <- eligible_ids[j]
      
      observer_checklists <- checklists[
        observer_id == observer
      ]
      
      sampled_index <- sample.int(
        n = nrow(observer_checklists),
        size = threshold_n,
        replace = FALSE
      )
      
      sampled_rows[[j]] <- data.table(
        observer_id = observer,
        sampled_mean_CAS = mean(
          observer_checklists$CAS[sampled_index],
          na.rm = TRUE
        )
      )
    }
    
    sampled_CAS <- rbindlist(sampled_rows)
    
    replicate_data <- merge(
      eligible_reference,
      sampled_CAS,
      by = "observer_id",
      all = FALSE
    )
    
    # Integrated score uses the three current components with equal
    # weights. ORI and SBI remain fixed; mean CAS is estimated from
    # the sampled checklist history.
    
    replicate_data[
      ,
      AviCRED_sampled := rowMeans(
        cbind(
          ORI_bounded,
          sampled_mean_CAS,
          SBI
        )
      )
    ]
    
    # Classify sampled scores using the fixed full-history quartiles.
    
    replicate_data[
      ,
      sampled_class := classify_by_fixed_breaks(
        AviCRED_sampled,
        class_breaks
      )
    ]
    
    replicate_data <- replicate_data[
      is.finite(AviCRED_reference) &
        is.finite(AviCRED_sampled) &
        !is.na(reference_class) &
        !is.na(sampled_class)
    ]
    
    n_compared <- nrow(replicate_data)
    
    result_counter <- result_counter + 1L
    
    replicate_results[[result_counter]] <- data.table(
      threshold_n = threshold_n,
      replicate = replicate_id,
      eligible_observers = n_eligible,
      compared_observers = n_compared,
      retained_percent_of_reference =
        100 * n_eligible / nrow(reference),
      CAS_MAE = mae_safe(
        replicate_data$mean_CAS_reference,
        replicate_data$sampled_mean_CAS
      ),
      CAS_RMSE = rmse_safe(
        replicate_data$mean_CAS_reference,
        replicate_data$sampled_mean_CAS
      ),
      CAS_Spearman_rho = spearman_safe(
        replicate_data$mean_CAS_reference,
        replicate_data$sampled_mean_CAS
      ),
      AviCRED_MAE = mae_safe(
        replicate_data$AviCRED_reference,
        replicate_data$AviCRED_sampled
      ),
      AviCRED_RMSE = rmse_safe(
        replicate_data$AviCRED_reference,
        replicate_data$AviCRED_sampled
      ),
      AviCRED_Spearman_rho = spearman_safe(
        replicate_data$AviCRED_reference,
        replicate_data$AviCRED_sampled
      ),
      four_class_exact_agreement = mean(
        replicate_data$reference_class ==
          replicate_data$sampled_class
      ),
      four_class_weighted_kappa =
        quadratic_weighted_kappa(
          replicate_data$reference_class,
          replicate_data$sampled_class,
          K = 4L
        )
    )
  }
}

replicate_results <- rbindlist(
  replicate_results,
  use.names = TRUE,
  fill = TRUE
)

# ============================================================
# 10. SUMMARISE REPLICATES FOR EACH THRESHOLD
# ============================================================

metric_columns <- c(
  "CAS_MAE",
  "CAS_RMSE",
  "CAS_Spearman_rho",
  "AviCRED_MAE",
  "AviCRED_RMSE",
  "AviCRED_Spearman_rho",
  "four_class_exact_agreement",
  "four_class_weighted_kappa"
)

summary_results <- replicate_results[
  ,
  {
    out <- list(
      threshold_n = first(threshold_n),
      n_replicates = .N,
      eligible_observers_mean = mean(
        eligible_observers,
        na.rm = TRUE
      ),
      compared_observers_mean = mean(
        compared_observers,
        na.rm = TRUE
      ),
      retained_percent_mean = mean(
        retained_percent_of_reference,
        na.rm = TRUE
      ),
      retained_percent_sd = sd(
        retained_percent_of_reference,
        na.rm = TRUE
      )
    )
    
    for (metric in metric_columns) {
      
      out[[paste0(metric, "_mean")]] <- mean(
        get(metric),
        na.rm = TRUE
      )
      
      out[[paste0(metric, "_sd")]] <- sd(
        get(metric),
        na.rm = TRUE
      )
    }
    
    out
  },
  by = threshold_n
]

setorder(summary_results, threshold_n)

# ============================================================
# 11. RETENTION AND DIMINISHING-RETURNS METRICS
# ============================================================

summary_results[
  ,
  `:=`(
    observer_retention_percent = retained_percent_mean,
    observer_loss_percent = 100 - retained_percent_mean
  )
]

summary_results[
  ,
  previous_threshold_n := shift(threshold_n)
]

summary_results[
  ,
  previous_AviCRED_rho := shift(
    AviCRED_Spearman_rho_mean
  )
]

summary_results[
  ,
  previous_AviCRED_MAE := shift(
    AviCRED_MAE_mean
  )
]

summary_results[
  ,
  rho_gain_from_previous_threshold :=
    AviCRED_Spearman_rho_mean -
    previous_AviCRED_rho
]

summary_results[
  ,
  MAE_reduction_from_previous_threshold :=
    previous_AviCRED_MAE -
    AviCRED_MAE_mean
]

summary_results[
  ,
  incremental_retention_loss_percentage_points :=
    shift(observer_retention_percent) -
    observer_retention_percent
]

summary_results[
  ,
  incremental_retention_loss_percentage_points :=
    ifelse(
      is.na(incremental_retention_loss_percentage_points),
      NA_real_,
      incremental_retention_loss_percentage_points
    )
]

# Efficiency of each increase in threshold:
# improvement in AviCRED rank agreement per percentage point of
# additional observer-retention loss.

summary_results[
  ,
  rho_gain_per_retention_point :=
    ifelse(
      is.finite(incremental_retention_loss_percentage_points) &
        incremental_retention_loss_percentage_points > 0,
      rho_gain_from_previous_threshold /
        incremental_retention_loss_percentage_points,
      NA_real_
    )
]

summary_results[
  ,
  previous_threshold_n := NULL
]

summary_results[
  ,
  previous_AviCRED_rho := NULL
]

summary_results[
  ,
  previous_AviCRED_MAE := NULL
]

# ============================================================
# 12. EXPORT ONE COMBINED CSV
# ============================================================

# One row per candidate threshold, with mean and SD across replicates.
# Replicate-level data are retained in memory as replicate_results.

summary_results[
  ,
  CAS_reference_MAE := CAS_reference_MAE
]

summary_results[
  ,
  CAS_reference_RMSE := CAS_reference_RMSE
]

summary_results[
  ,
  CAS_reference_Spearman_rho := CAS_reference_rho
]

summary_results[
  ,
  full_reference_observers := nrow(reference)
]

summary_results[
  ,
  checklist_data_observers := uniqueN(checklists$observer_id)
]

summary_results[
  ,
  class_1_lower_break := class_breaks[1]
]

summary_results[
  ,
  class_2_break := class_breaks[2]
]

summary_results[
  ,
  class_3_break := class_breaks[3]
]

summary_results[
  ,
  class_4_break := class_breaks[4]
]

summary_results[
  ,
  class_5_upper_break := class_breaks[5]
]

data.table::fwrite(
  summary_results,
  OUT_FILE,
  na = ""
)

cat("\n============================================================\n")
cat("THRESHOLD ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("\nOutput CSV:\n", OUT_FILE, "\n", sep = "")

cat("\nThreshold summary:\n")

print(
  summary_results[
    ,
    .(
      threshold_n,
      eligible_observers_mean,
      observer_retention_percent,
      CAS_MAE_mean,
      CAS_Spearman_rho_mean,
      AviCRED_MAE_mean,
      AviCRED_RMSE_mean,
      AviCRED_Spearman_rho_mean,
      four_class_exact_agreement_mean,
      four_class_weighted_kappa_mean,
      rho_gain_from_previous_threshold
    )
  ]
)

####====== AFTETR THRESHOLDING 50 IS DECIDED TO BE CUTOFF CHECKLIST THRESHOLD =======##

# ==================== CHECKLIST ISOLATION FOREST ========================#

checklist_master <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
)

print(dim(checklist_master))
print(names(checklist_master))
print(sapply(checklist_master, class))

# ------------------------------------------------------------
# A. CHECKLIST MASTER: MISSINGNESS AND UNIQUE VALUES
# ------------------------------------------------------------

checklist_variables <- c(
  "richness",
  "abundance_per_checklist",
  "duration_min",
  "distance_km",
  "n_observers",
  "complete_flag",
  "protocol",
  "state",
  "county",
  "locality_id",
  "latitude",
  "longitude",
  "obs_date"
)

checklist_diagnostics <- data.frame(
  variable = checklist_variables,
  missing_n = sapply(
    checklist_master[checklist_variables],
    function(x) sum(is.na(x))
  ),
  unique_n = sapply(
    checklist_master[checklist_variables],
    function(x) length(unique(x[!is.na(x)]))
  ),
  row.names = NULL
)

print(checklist_diagnostics)

# ------------------------------------------------------------
# B. CHECKLIST ID UNIQUENESS
# ------------------------------------------------------------

cat(
  "\nChecklist rows:",
  nrow(checklist_master),
  "\nUnique checklist IDs:",
  dplyr::n_distinct(checklist_master$checklist_id),
  "\nDuplicated checklist IDs:",
  sum(duplicated(checklist_master$checklist_id)),
  "\n"
)

# ------------------------------------------------------------
# C. OBSERVATION DISTRIBUTIONS
# ------------------------------------------------------------

summary(
  checklist_master[, c(
    "richness",
    "abundance_per_checklist",
    "duration_min",
    "distance_km",
    "n_observers"
  )]
)

cat("\nProtocol frequencies:\n")
print(table(checklist_master$protocol, useNA = "ifany"))

cat("\nCompleteness frequencies:\n")
print(table(checklist_master$complete_flag, useNA = "ifany"))

# ------------------------------------------------------------
# D. DERIVE MONTH AND METEOROLOGICAL SEASON
# ------------------------------------------------------------

checklist_master$obs_date <- as.Date(
  checklist_master$obs_date
)

checklist_master$month <- as.integer(
  format(checklist_master$obs_date, "%m")
)

checklist_master$season <- case_when(
  checklist_master$month %in% c(12, 1, 2) ~ "Winter",
  checklist_master$month %in% 3:5 ~ "Pre-monsoon",
  checklist_master$month %in% 6:9 ~ "Monsoon",
  checklist_master$month %in% 10:11 ~ "Post-monsoon",
  TRUE ~ NA_character_
)

print(table(checklist_master$season, useNA = "ifany"))

# ============================================================
# 1. PREPARE CHECKLIST DATA
# ============================================================

set.seed(20261011)

checklists <- checklist_master %>%
  mutate(
    obs_date = as.Date(obs_date),
    month = as.integer(format(obs_date, "%m")),
    
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% 3:5 ~ "Pre-monsoon",
      month %in% 6:9 ~ "Monsoon",
      month %in% 10:11 ~ "Post-monsoon",
      TRUE ~ NA_character_
    ),
    
    richness = as.numeric(richness),
    abundance_per_checklist =
      as.numeric(abundance_per_checklist),
    duration_min = as.numeric(duration_min),
    distance_km = as.numeric(distance_km),

    log_abundance = log1p(abundance_per_checklist),
    log_duration = log1p(duration_min),
    log_distance = log1p(distance_km),
    
    group = paste(state, season, sep = "_")
  )

stopifnot(!anyDuplicated(checklists$checklist_id))

# ============================================================
# 2. DEFINE EXPANDED PREDICTORS
# ============================================================

numeric_features <- c(
  "richness",
  "log_abundance",
  "log_duration",
  "log_distance"
)

checklists$if_anomaly_score <- NA_real_
checklists$if_typicality <- NA_real_
checklists$if_anomaly_5pct <- NA_integer_
checklists$if_anomaly_10pct <- NA_integer_

valid <- complete.cases(
  checklists[, c(
    numeric_features,
    "protocol",
    "state",
    "season"
  )]
)

valid <- valid & apply(
  as.matrix(checklists[, numeric_features, drop = FALSE]),
  1,
  function(x) all(is.finite(x))
)

# Continue with the existing State x Season forest loop,
# score transformation, and coverage reporting.

# ============================================================
# 4. FIT STATE x SEASON FORESTS
# ============================================================

groups <- unique(checklists$group[valid])

group_summary <- data.frame(
  group = character(),
  n_checklists = integer(),
  status = character(),
  stringsAsFactors = FALSE
)

for (g in groups) {
  
  idx <- which(valid & checklists$group == g)
  n <- length(idx)
  
  # Retain the original minimum group size
  if (n < 100) {
    group_summary <- rbind(
      group_summary,
      data.frame(
        group = g,
        n_checklists = n,
        status = "Skipped: fewer than 100"
      )
    )
    next
  }
  
  # Encode protocol as categorical dummy variables
  protocol_matrix <- model.matrix(
    ~ protocol - 1,
    data = checklists[idx, , drop = FALSE]
  )
  
  X <- cbind(
    as.matrix(checklists[idx, numeric_features]),
    protocol_matrix
  )
  
  # Remove constant features within each group
  keep <- apply(
    X,
    2,
    function(x) length(unique(x)) > 1
  )
  
  X <- X[, keep, drop = FALSE]
  
  if (ncol(X) == 0) {
    group_summary <- rbind(
      group_summary,
      data.frame(
        group = g,
        n_checklists = n,
        status = "Skipped: no varying predictors"
      )
    )
    next
  }
  
  model <- isolation.forest(
    X,
    ntrees = 300,
    sample_size = min(256, n),
    ndim = 1,
    nthreads = 2
  )
  
  # Score in batches to limit prediction memory
  batch_size <- 50000
  group_scores <- numeric(n)
  
  for (start in seq(1, n, by = batch_size)) {
    
    end <- min(start + batch_size - 1, n)
    rows <- start:end
    
    group_scores[rows] <- as.numeric(
      predict(
        model,
        X[rows, , drop = FALSE],
        type = "score"
      )
    )
  }
  
  checklists$if_anomaly_score[idx] <- group_scores
  
  group_summary <- rbind(
    group_summary,
    data.frame(
      group = g,
      n_checklists = n,
      status = "Scored"
    )
  )
  
  rm(model, X, protocol_matrix, group_scores)
  gc(verbose = FALSE)
}

# ============================================================
# 5. SCORE DIRECTION AND IF TYPICALITY TRANSFORMATION
# ============================================================

summary(checklists$if_anomaly_score)
range(checklists$if_anomaly_score, na.rm = TRUE)

# Assuming the installed isotree version returns anomaly
# scores where higher values indicate greater atypicality:
#
# Higher anomaly score  = more atypical
# Higher IF typicality  = more typical
#
# This transformation reverses the score direction.
# It does NOT establish genuine checklist authenticity.

scored <- is.finite(checklists$if_anomaly_score)

checklists$if_typicality <- NA_real_

checklists$if_typicality[scored] <-
  1 - checklists$if_anomaly_score[scored]

cat("\nIF anomaly score summary:\n")
print(summary(checklists$if_anomaly_score))

cat("\nIF typicality score summary:\n")
print(summary(checklists$if_typicality))

cat("\nIF typicality score range:\n")
print(range(checklists$if_typicality, na.rm = TRUE))

# ============================================================
# 6. SUMMARISE IF ANOMALY AND TYPICALITY SCORES
# ============================================================

cat("\nIF anomaly score summary:\n")
print(summary(checklists$if_anomaly_score))

cat("\nIF typicality summary:\n")
print(summary(checklists$if_typicality))


# ============================================================
# 7. FLAG TOP 5% AND TOP 10% ANOMALOUS / TYPICAL CHECKLISTS
# ============================================================

# Anomaly flags: highest anomaly scores
q95 <- quantile(
  checklists$if_anomaly_score[scored],
  0.95,
  na.rm = TRUE
)

q90 <- quantile(
  checklists$if_anomaly_score[scored],
  0.90,
  na.rm = TRUE
)

checklists$if_anomaly_5pct <- NA_integer_
checklists$if_anomaly_10pct <- NA_integer_

checklists$if_anomaly_5pct[scored] <-
  as.integer(checklists$if_anomaly_score[scored] >= q95)

checklists$if_anomaly_10pct[scored] <-
  as.integer(checklists$if_anomaly_score[scored] >= q90)

# Typicality flags: highest IF typicality scores
q_typical_95 <- quantile(
  checklists$if_typicality[scored],
  0.95,
  na.rm = TRUE
)

q_typical_90 <- quantile(
  checklists$if_typicality[scored],
  0.90,
  na.rm = TRUE
)

checklists$if_top5_typicality <- NA_integer_
checklists$if_top10_typicality <- NA_integer_

checklists$if_top5_typicality[scored] <-
  as.integer(
    checklists$if_typicality[scored] >= q_typical_95
  )

checklists$if_top10_typicality[scored] <-
  as.integer(
    checklists$if_typicality[scored] >= q_typical_90
  )

# Report thresholds and counts
cat("\nIF anomaly thresholds:\n")
print(c(
  top_5pct_anomaly = q95,
  top_10pct_anomaly = q90
))

cat("\nIF typicality thresholds:\n")
print(c(
  top_5pct_typicality = q_typical_95,
  top_10pct_typicality = q_typical_90
))

cat("\nAnomaly flag counts:\n")
print(table(checklists$if_anomaly_5pct, useNA = "ifany"))
print(table(checklists$if_anomaly_10pct, useNA = "ifany"))

cat("\nTypicality flag counts:\n")
print(table(checklists$if_top5_typicality, useNA = "ifany"))
print(table(checklists$if_top10_typicality, useNA = "ifany"))

# ============================================================
# 8. GROUP COVERAGE AND EXPORT
# ============================================================

print(table(group_summary$status))

print(
  head(
    group_summary[
      order(-group_summary$n_checklists),
    ],
    20
  )
)

cat("\nTotal checklists:", nrow(checklists), "\n")
cat("Scored checklists:", sum(scored), "\n")
cat("Unscored checklists:", sum(!scored), "\n")

cat(
  "Percentage scored:",
  round(100 * mean(scored), 2),
  "%\n"
)

write.csv(
  checklists,
  "Isolation checklists new.csv",
  row.names = FALSE
)

# ============================================================
# 9. LOAD CAS DATA
# ============================================================

cas_final <- read.csv(
  "CAS_FINAL.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("CAS file rows:", nrow(cas_final), "\n")
print(names(cas_final))

# Identify the expected columns
required_cas_cols <- c("checklist_id", "CAS")

if (!all(required_cas_cols %in% names(cas_final))) {
  stop(
    "cas_final.csv must contain checklist_id and CAS. ",
    "Check the column names printed above."
  )
}

cas_final$checklist_id <- as.character(
  cas_final$checklist_id
)

cas_final$CAS <- as.numeric(cas_final$CAS)

# Remove missing or non-finite CAS values
cas_final <- cas_final[
  !is.na(cas_final$checklist_id) &
    nzchar(cas_final$checklist_id) &
    is.finite(cas_final$CAS),
  c("checklist_id", "CAS"),
  drop = FALSE
]

# Ensure one CAS value per checklist
if (anyDuplicated(cas_final$checklist_id)) {
  stop(
    "Duplicate checklist_id values detected in cas_final.csv. ",
    "Resolve duplicates before validation."
  )
}

# ============================================================
# 10. MATCH CAS TO CHECKLIST IF RESULTS
# ============================================================

checklists$checklist_id <- as.character(
  checklists$checklist_id
)

if ("CAS" %in% names(checklists)) {
  checklists$CAS <- NULL
}

cas_index <- match(
  checklists$checklist_id,
  cas_final$checklist_id
)

checklists$CAS <- cas_final$CAS[cas_index]

cat("\nChecklist-level CAS matching\n")
cat("Total IF checklists:", nrow(checklists), "\n")
cat(
  "Matched CAS checklists:",
  sum(is.finite(checklists$CAS)),
  "\n"
)
cat(
  "Unmatched checklists:",
  sum(!is.finite(checklists$CAS)),
  "\n"
)

validation_data <- checklists[
  is.finite(checklists$CAS) &
    is.finite(checklists$if_anomaly_score) &
    is.finite(checklists$if_typicality),
  ,
  drop = FALSE
]

cat(
  "Checklists used for validation:",
  nrow(validation_data),
  "\n"
)

if (nrow(validation_data) < 3) {
  stop("Insufficient matched checklists for validation.")
}

# ============================================================
# 11. CAS CLASS COMPARISON
# ============================================================

validation_data$CAS_class <- cut(
  validation_data$CAS,
  breaks = c(-Inf, 0.25, 0.50, 0.75, 0.90, Inf),
  labels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  include.lowest = TRUE
)

print(table(
  validation_data$CAS_class,
  useNA = "ifany"
))

print(
  aggregate(
    if_typicality ~ CAS_class,
    data = validation_data,
    FUN = mean
  )
)

print(
  aggregate(
    if_typicality ~ CAS_class,
    data = validation_data,
    FUN = median
  )
)

# ============================================================
# 12. TOP 10% ANOMALIES VS CAS
# ============================================================

cutoff_10 <- quantile(
  validation_data$if_anomaly_score,
  probs = 0.90,
  na.rm = TRUE
)

validation_data$IF_top10 <- as.integer(
  validation_data$if_anomaly_score >= cutoff_10
)

cat("\nMean CAS by IF anomaly group:\n")

print(
  aggregate(
    CAS ~ IF_top10,
    data = validation_data,
    FUN = mean
  )
)

cat("\nCAS classes by IF top-10% anomaly flag:\n")

print(
  table(
    validation_data$CAS_class,
    validation_data$IF_top10
  )
)

# ============================================================
# 14. HIGH CAS VS IF TOP 10% ANOMALY OVERLAP
# ============================================================

validation_data$CAS_high <- as.integer(
  validation_data$CAS >= 0.75
)

comparison <- table(
  CAS_high = factor(
    validation_data$CAS_high,
    levels = c(0, 1)
  ),
  IF_top10 = factor(
    validation_data$IF_top10,
    levels = c(0, 1)
  )
)

cat("\nHigh-CAS vs IF top-10% anomaly table:\n")
print(comparison)

n_high_CAS <- sum(validation_data$CAS_high == 1)
n_IF_top10 <- sum(validation_data$IF_top10 == 1)

n_overlap <- sum(
  validation_data$CAS_high == 1 &
    validation_data$IF_top10 == 1
)

cat("\nHigh-CAS checklists:", n_high_CAS, "\n")
cat("IF top-10% anomalous checklists:", n_IF_top10, "\n")
cat("Overlap:", n_overlap, "\n")

cat(
  "Percentage of high-CAS checklists flagged as anomalous:",
  if (n_high_CAS > 0) {
    round(100 * n_overlap / n_high_CAS, 2)
  } else {
    NA_real_
  },
  "%\n"
)

cat(
  "Percentage of IF top-10% anomalies with high CAS:",
  if (n_IF_top10 > 0) {
    round(100 * n_overlap / n_IF_top10, 2)
  } else {
    NA_real_
  },
  "%\n"
)

# ============================================================
# 15. EXPORT VALIDATION RESULTS
# ============================================================

write.csv(
  validation_data,
  "Checklists final new.csv",
  row.names = FALSE,
  na = ""
)

#================ OBSERVER ISOLATION FOREST =================#

set.seed(20261011)

observer_final <- read.csv(
  "AVICRED NEW.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================================
# 1. DEFINE BASE-COMPONENT PREDICTORS
# ============================================================

feature_cols <- c(
  "years_active",
  "n_checklists",
  "total_species",
  "consistency_score",
  "weighted_completeness_score"
)

required_cols <- c(
  "observer_id",
  feature_cols,
  "AviCRED"
)

missing_cols <- setdiff(
  required_cols,
  names(observer_final)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Convert predictors and outcome to numeric
for (v in c(feature_cols, "AviCRED")) {
  observer_final[[v]] <- as.numeric(
    observer_final[[v]]
  )
}

# ============================================================
# 2. RETAIN COMPLETE, FINITE RECORDS
# ============================================================

complete_rows <- complete.cases(
  observer_final[, c(feature_cols, "AviCRED"),
                 drop = FALSE]
)

finite_rows <- apply(
  as.matrix(
    observer_final[, feature_cols, drop = FALSE]
  ),
  1,
  function(x) all(is.finite(x))
) & is.finite(observer_final$AviCRED)

ml_complete <- observer_final[
  complete_rows & finite_rows,
  ,
  drop = FALSE
]

if (anyDuplicated(ml_complete$observer_id)) {
  stop("Duplicate observer IDs detected.")
}

cat("Total observers:", nrow(observer_final), "\n")
cat("Complete observers:", nrow(ml_complete), "\n")
cat(
  "Excluded observers:",
  nrow(observer_final) - nrow(ml_complete),
  "\n"
)

if (nrow(ml_complete) < 2) {
  stop("Insufficient complete observer records.")
}

# ============================================================
# 3. STANDARDISE BASE-COMPONENT PREDICTORS
# ============================================================

features <- scale(
  ml_complete[, feature_cols, drop = FALSE]
)

# Remove constant predictors, if any
varying <- apply(
  features,
  2,
  function(x) length(unique(x)) > 1
)

features <- features[, varying, drop = FALSE]

if (ncol(features) == 0) {
  stop("No varying predictors available.")
}

if (any(!is.finite(features))) {
  stop("Non-finite values detected after scaling.")
}

# ============================================================
# 4. FIT ISOLATION FOREST
# ============================================================

iso_model <- isolation.forest(
  features,
  ntrees = 500,
  sample_size = min(256, nrow(features)),
  ndim = 1,
  nthreads = 2
)

ml_complete$ml_score <- as.numeric(
  predict(
    iso_model,
    features,
    type = "score"
  )
)

cat("\nIF anomaly score summary:\n")
print(summary(ml_complete$ml_score))

# ============================================================
# 5. TOP 10% ANOMALY FLAGS
# ============================================================

cutoff_10 <- quantile(
  ml_complete$ml_score,
  probs = 0.90,
  na.rm = TRUE
)

ml_complete$ml_10pct <- as.integer(
  ml_complete$ml_score >= cutoff_10
)

print(table(ml_complete$ml_10pct))

# ============================================================
# 6. AVICRED QUINTILE CLASSES
# ============================================================

ml_complete$avicred_quintile <- dplyr::ntile(
  ml_complete$AviCRED,
  5
)

ml_complete$avicred_class <- factor(
  ml_complete$avicred_quintile,
  levels = 1:5,
  labels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  ordered = TRUE
)

print(table(ml_complete$avicred_class))

# ============================================================
# 7. INDEPENDENT VALIDATION
# ============================================================

cat("\nSpearman correlation: IF score vs AviCRED\n")

print(
  cor.test(
    ml_complete$ml_score,
    ml_complete$AviCRED,
    method = "spearman",
    exact = FALSE
  )
)

cat("\nMean AviCRED by IF anomaly flag:\n")

print(
  aggregate(
    AviCRED ~ ml_10pct,
    data = ml_complete,
    FUN = mean
  )
)

cat("\nAviCRED classes by IF top-10% flag:\n")

print(
  table(
    ml_complete$avicred_class,
    ml_complete$ml_10pct
  )
)

# ============================================================
# 8. EXPORT RESULTS
# ============================================================

write.csv(
  ml_complete,
  "Observer final new.csv",
  row.names = FALSE,
  na = ""
)

#================= ISOLATION FOREST VS AVICRED/CAS CORRELATION =================#

############ OBSERVER LEVEL #################

observer_final_new <- read.csv(
  file.choose(),          # Select observer_final_new.csv
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(observer_final_new)

# 1. Prepare the analysis dataset
analysis_data <- observer_final_new %>%
  transmute(
    observer_id,
    ml_score = as.numeric(ml_score),
    AviCRED = as.numeric(AviCRED),
    ORI_bounded = as.numeric(ORI_bounded),
    mean_CAS = as.numeric(mean_CAS),
    SBI = as.numeric(SBI)
  ) %>%
  filter(
    is.finite(ml_score)
  )

metrics <- c("AviCRED", "ORI_bounded", "mean_CAS", "SBI")

# 2. Calculate correlations using pairwise-complete observations
stats_list <- lapply(metrics, function(metric) {
  
  d <- analysis_data %>%
    select(ml_score, all_of(metric)) %>%
    filter(
      is.finite(ml_score),
      is.finite(.data[[metric]])
    )
  
  pearson_test <- cor.test(
    d$ml_score,
    d[[metric]],
    method = "pearson"
  )
  
  spearman_test <- cor.test(
    d$ml_score,
    d[[metric]],
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    Metric = metric,
    N = nrow(d),
    Pearson_r = unname(pearson_test$estimate),
    Pearson_CI_lower = pearson_test$conf.int[1],
    Pearson_CI_upper = pearson_test$conf.int[2],
    Pearson_p = pearson_test$p.value,
    Spearman_rho = unname(spearman_test$estimate),
    Spearman_p = spearman_test$p.value,
    R2 = unname(pearson_test$estimate)^2
  )
})

stats_df <- bind_rows(stats_list)

cat("\nNEW MODEL: ISOLATION FOREST ASSOCIATIONS\n")
print(stats_df, row.names = FALSE)

write.csv(
  stats_df,
  "New_AviCRED_IF_Correlation_Results.csv",
  row.names = FALSE
)

# 3. Create annotation labels
stats_df$label <- paste0(
  "N = ", format(stats_df$N, big.mark = ","),
  "\nPearson r = ", sprintf("%.3f", stats_df$Pearson_r),
  "\nSpearman \u03c1 = ", sprintf("%.3f", stats_df$Spearman_rho),
  "\nR\u00b2 = ", sprintf("%.3f", stats_df$R2)
)

# 4. Long-format data for plotting
corr_df <- analysis_data %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(
    is.finite(ml_score),
    is.finite(Value)
  )

corr_df$Metric <- factor(
  corr_df$Metric,
  levels = metrics
)

stats_df$Metric <- factor(
  stats_df$Metric,
  levels = metrics
)

# 5. Plot
p <- ggplot(
  corr_df,
  aes(x = Value, y = ml_score)
) +
  geom_point(
    alpha = 0.20,
    size = 0.8
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "blue",
    fill = "lightblue",
    alpha = 0.25,
    linewidth = 1.2
  ) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    scales = "free_x"
  ) +
  geom_text(
    data = stats_df,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.1,
    family = "serif",
    size = 3.5
  ) +
  labs(
    x = "Observer-level metric",
    y = "Isolation Forest score"
  ) +
  theme_classic(
    base_family = "serif"
  ) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(size = 10),
    strip.text = element_text(
      size = 11,
      face = "bold"
    )
  )

print(p)

ggsave(
  "New_AviCRED_IF_Comparison.jpg",
  plot = p,
  width = 15,
  height = 5.5,
  dpi = 400,
  bg = "white"
)

# ============================================================
# SIMPLE LINEAR REGRESSION: IF SCORE AGAINST AviCRED COMPONENTS
# ============================================================

library(broom)

# 1. Prepare observer-level data
regression_data <- observer_final_new %>%
  transmute(
    ml_score   = as.numeric(ml_score),
    ORI        = as.numeric(ORI_bounded),
    SBI        = as.numeric(SBI),
    CAS        = as.numeric(mean_CAS),
    AviCRED    = as.numeric(AviCRED)
  ) %>%
  filter(
    is.finite(ml_score),
    is.finite(ORI),
    is.finite(SBI),
    is.finite(CAS),
    is.finite(AviCRED)
  )

cat("Observers included:", nrow(regression_data), "\n")

# 2. Fit four separate simple linear regression models
model_ORI <- lm(ml_score ~ ORI, data = regression_data)
model_SBI <- lm(ml_score ~ SBI, data = regression_data)
model_CAS <- lm(ml_score ~ CAS, data = regression_data)
model_AviCRED <- lm(ml_score ~ AviCRED, data = regression_data)

models <- list(
  ORI_bounded = model_ORI,
  SBI         = model_SBI,
  mean_CAS    = model_CAS,
  AviCRED     = model_AviCRED
)

# 3. Extract model-level statistics
model_summary <- bind_rows(
  lapply(names(models), function(nm) {
    m <- models[[nm]]
    s <- summary(m)
    
    data.frame(
      Metric = nm,
      N = nobs(m),
      Intercept = unname(coef(m)[1]),
      Slope = unname(coef(m)[2]),
      R2 = s$r.squared,
      Adjusted_R2 = s$adj.r.squared,
      Model_p = pf(
        s$fstatistic[1],
        s$fstatistic[2],
        s$fstatistic[3],
        lower.tail = FALSE
      )
    )
  })
)

# 4. Extract slope confidence intervals and coefficient p-values
slope_results <- bind_rows(
  lapply(names(models), function(nm) {
    m <- models[[nm]]
    
    tidy(m, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      transmute(
        Metric = nm,
        Slope_CI_lower = conf.low,
        Slope_CI_upper = conf.high,
        Slope_p = p.value
      )
  })
)

# 5. Combine results
regression_results <- left_join(
  model_summary,
  slope_results,
  by = "Metric"
)

print(regression_results, row.names = FALSE)

# 6. Export results
write.csv(
  regression_results,
  "AviCRED_IF_Simple_Linear_Regression.csv",
  row.names = FALSE
)

# 7. Full diagnostic output for each model
summary(model_ORI)
summary(model_SBI)
summary(model_CAS)
summary(model_AviCRED)

# ============================================================
# CAS vs ISOLATION FOREST TYPICALITY
# REGRESSION, CLUSTERED INFERENCE, DIAGNOSTICS, AND CSV EXPORTS
# ============================================================

library(sandwich)
library(lmtest)
library(splines)

# ------------------------------------------------------------
# 0. OUTPUT DIRECTORY
# ------------------------------------------------------------

output_dir <- "file path/CAS regression"

save_csv <- function(x, filename) {
  write.csv(
    x,
    file.path(output_dir, filename),
    row.names = FALSE,
    na = "")
}

# ============================================================
# CAS vs ISOLATION FOREST TYPICALITY
# CORRELATION, REGRESSION, CLUSTERED INFERENCE AND VISUALISATION
# ============================================================

library(sandwich)

set.seed(20261017)

# ============================================================
# 2. IMPORT AND PREPARE DATA
# ============================================================

checklists <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("CAS", "if_typicality", "observer_id")

if (!all(required_cols %in% names(checklists))) {
  stop(
    "Missing required columns: ",
    paste(
      setdiff(required_cols, names(checklists)),
      collapse = ", "
    )
  )
}

analysis_data <- data.frame(
  CAS = suppressWarnings(as.numeric(checklists$CAS)),
  IF_typicality = suppressWarnings(
    as.numeric(checklists$if_typicality)
  ),
  observer_id = trimws(as.character(checklists$observer_id))
) %>%
  filter(
    is.finite(CAS),
    is.finite(IF_typicality),
    !is.na(observer_id),
    observer_id != "",
    observer_id != "NA"
  )

# Remove unused factor levels if applicable
analysis_data$observer_id <- factor(analysis_data$observer_id)

if (nrow(analysis_data) < 3) {
  stop("Insufficient valid observations.")
}

if (dplyr::n_distinct(analysis_data$CAS) < 2 ||
    dplyr::n_distinct(analysis_data$IF_typicality) < 2) {
  stop("CAS and IF Typicality must both vary.")
}

n_checklists <- nrow(analysis_data)
n_observers <- dplyr::n_distinct(analysis_data$observer_id)

cat("\n========== DATA SUMMARY ==========\n")
cat("Valid checklists:", format(n_checklists, big.mark = ","), "\n")
cat("Unique observers:", format(n_observers, big.mark = ","), "\n")
cat("Excluded rows:", nrow(checklists) - n_checklists, "\n")

cat("\nCAS summary:\n")
print(summary(analysis_data$CAS))

cat("\nIF Typicality summary:\n")
print(summary(analysis_data$IF_typicality))

# ============================================================
# 3. CORRELATION ANALYSES
# ============================================================

pearson <- cor.test(
  analysis_data$CAS,
  analysis_data$IF_typicality,
  method = "pearson"
)

spearman <- cor.test(
  analysis_data$CAS,
  analysis_data$IF_typicality,
  method = "spearman",
  exact = FALSE
)

cat("\n========== PEARSON CORRELATION ==========\n")
print(pearson)

cat("\n========== SPEARMAN CORRELATION ==========\n")
print(spearman)

# IMPORTANT:
# These are ordinary correlation tests.
# Their p-values and Pearson CI do NOT account for repeated
# checklists contributed by the same observer.
# Cluster-adjusted regression inference is performed below.

# ============================================================
# 4. CAS DECILE SUMMARY
# ============================================================

cas_bins <- analysis_data %>%
  mutate(
    CAS_decile = dplyr::ntile(CAS, 10)
  ) %>%
  group_by(CAS_decile) %>%
  summarise(
    n = n(),
    CAS_median = median(CAS),
    IF_median = median(IF_typicality),
    IF_mean = mean(IF_typicality),
    .groups = "drop"
  )

cat("\n========== CAS DECILE SUMMARY ==========\n")
print(cas_bins)

# ============================================================
# 5. FIT LINEAR AND QUADRATIC MODELS
# ============================================================

lm_model <- lm(
  IF_typicality ~ CAS,
  data = analysis_data
)

quad_model <- lm(
  IF_typicality ~ CAS + I(CAS^2),
  data = analysis_data
)

cat("\n========== LINEAR MODEL ==========\n")
print(summary(lm_model))

cat("\n========== QUADRATIC MODEL ==========\n")
print(summary(quad_model))

# AIC is a supplementary fit comparison.
# It does not account for observer clustering.
cat("\n========== AIC COMPARISON ==========\n")
print(AIC(lm_model, quad_model))

# ============================================================
# 6. OBSERVER-CLUSTERED INFERENCE
# ============================================================

# Account for multiple checklists from the same observer.
# HC1 cluster-robust covariance estimator.

cluster_vcov <- sandwich::vcovCL(
  quad_model,
  cluster = analysis_data$observer_id,
  type = "HC1"
)

cluster_est <- coef(quad_model)
cluster_se <- sqrt(diag(cluster_vcov))

# Use a t approximation with number of observer clusters - 1
# degrees of freedom for confidence intervals and p-values.
cluster_df <- n_observers - 1

if (cluster_df < 1) {
  stop("Insufficient observer clusters for clustered inference.")
}

cluster_t <- cluster_est / cluster_se

cluster_p <- 2 * pt(
  abs(cluster_t),
  df = cluster_df,
  lower.tail = FALSE
)

critical_t <- qt(0.975, df = cluster_df)

cluster_results <- data.frame(
  Term = names(cluster_est),
  Estimate = unname(cluster_est),
  Cluster_SE = unname(cluster_se),
  CI_Lower = unname(cluster_est - critical_t * cluster_se),
  CI_Upper = unname(cluster_est + critical_t * cluster_se),
  t_value = unname(cluster_t),
  df = cluster_df,
  p_value = unname(cluster_p),
  row.names = NULL
)

cat("\n========== QUADRATIC MODEL: CLUSTER-ROBUST RESULTS ==========\n")
print(cluster_results)

# Test whether the quadratic term contributes evidence of curvature.
quadratic_term <- "I(CAS^2)"

quadratic_result <- cluster_results[
  cluster_results$Term == quadratic_term,
  ,
  drop = FALSE
]

cat("\n========== CLUSTER-ROBUST TEST OF CURVATURE ==========\n")
print(quadratic_result)

if (nrow(quadratic_result) == 0) {
  stop("Could not locate the quadratic coefficient.")
}

# ============================================================
# 7. QUADRATIC TURNING POINT
# ============================================================

b1 <- unname(coef(quad_model)["CAS"])
b2 <- unname(coef(quad_model)[quadratic_term])

cas_min <- min(analysis_data$CAS)
cas_max <- max(analysis_data$CAS)

if (is.finite(b2) &&
    abs(b2) > .Machine$double.eps) {
  
  turning_point <- -b1 / (2 * b2)
  
  cat("\n========== QUADRATIC TURNING POINT ==========\n")
  cat("Turning point:", turning_point, "\n")
  cat("Observed CAS range:", cas_min, "to", cas_max, "\n")
  
  if (is.finite(turning_point) &&
      turning_point >= cas_min &&
      turning_point <= cas_max) {
    cat("The turning point lies within the observed CAS range.\n")
  } else {
    cat("Caution: the turning point lies outside the observed CAS range.\n")
  }
  
  if (quadratic_result$p_value < 0.05) {
    cat("The quadratic coefficient is significant at the 0.05 level ",
        "under the cluster-based t approximation.\n", sep = "")
  } else {
    cat("The quadratic coefficient is not significant at the 0.05 level ",
        "under the cluster-based t approximation.\n", sep = "")
  }
  
} else {
  turning_point <- NA_real_
  cat("\nTurning point not calculated: quadratic coefficient is effectively zero.\n")
}

# ============================================================
# 8. PREDICTION CURVES
# ============================================================

prediction_grid <- data.frame(
  CAS = seq(cas_min, cas_max, length.out = 300)
)

prediction_grid$Linear <- predict(
  lm_model,
  newdata = prediction_grid
)

prediction_grid$Quadratic <- predict(
  quad_model,
  newdata = prediction_grid
)

# Sample points for visualisation only.
# All retained observations are still used in model fitting.
set.seed(20261017)

plot_n <- min(30000, n_checklists)

plot_data <- analysis_data[
  sample.int(n_checklists, plot_n),
  ,
  drop = FALSE
]

# ============================================================
# 9. FINAL PUBLICATION FIGURE
# ============================================================

p_cas_if <- ggplot() +
  
  geom_point(
    data = plot_data,
    aes(x = CAS, y = IF_typicality),
    alpha = 0.10,
    size = 0.4,
    colour = "steelblue"
  ) +
  
  geom_line(
    data = prediction_grid,
    aes(x = CAS, y = Linear),
    colour = "red",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  geom_line(
    data = prediction_grid,
    aes(x = CAS, y = Quadratic),
    colour = "darkgreen",
    linewidth = 1.1
  ) +
  
  geom_line(
    data = cas_bins,
    aes(x = CAS_median, y = IF_median),
    colour = "black",
    linewidth = 0.8
  ) +
  
  geom_point(
    data = cas_bins,
    aes(x = CAS_median, y = IF_median),
    colour = "black",
    size = 2.5
  ) +
  
  labs(
    title = "Checklist Atypicality vs Isolation Forest Typicality",
    subtitle = paste0(
      "Linear fit (red dashed); quadratic fit (green); ",
      "decile medians (black); plotted points = ",
      format(plot_n, big.mark = ","),
      "; analysed checklists = ",
      format(n_checklists, big.mark = ",")
    ),
    x = "Checklist Atypicality Score (CAS)",
    y = "Isolation Forest Typicality Score"
  ) +
  
  theme_classic(
    base_family = "Times New Roman",
    base_size = 14
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(size = 12)
  )

print(p_cas_if)

ggplot2::ggsave(
  filename = file.path(
    getwd(),
    "CAS_IF_Regression.jpg"
  ),
  plot = p_cas_if,
  device = "jpeg",
  dpi = 400,
  width = 12,
  height = 8,
  units = "in",
  quality = 100,
  bg = "white"
)

cat("\nFigure saved to:\n")
cat(file.path(getwd(), "CAS_IF_Regression.jpg"), "\n")

# ============================================================
# 10. RESIDUAL DIAGNOSTICS
# ============================================================

# Diagnostics based on the full-data quadratic model.
# Sample rows for plotting to avoid excessive visual density.

set.seed(20261017)

diag_n <- min(100000, n_checklists)

diag_rows <- sample.int(n_checklists, diag_n)

diag_data <- data.frame(
  Fitted = fitted(quad_model)[diag_rows],
  Residual = residuals(quad_model)[diag_rows]
)

# Residuals vs fitted
p_residuals <- ggplot(
  diag_data,
  aes(x = Fitted, y = Residual)
) +
  geom_point(
    alpha = 0.15,
    size = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    colour = "blue"
  ) +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted IF Typicality",
    y = "Residuals"
  ) +
  theme_classic(
    base_family = "Times New Roman",
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

print(p_residuals)

ggplot2::ggsave(
  filename = file.path(
    getwd(),
    "CAS_IF_Residual_Diagnostics.jpg"
  ),
  plot = p_residuals,
  device = "jpeg",
  dpi = 400,
  width = 10,
  height = 7,
  units = "in",
  quality = 100,
  bg = "white"
)

# ============================================================
# AVICRED COMPONENT CORRELATIONS
# SBI vs CAS | ORI vs CAS | ORI vs SBI
# Points filled by AviCRED credibility class
# Pearson and Spearman correlations shown in each panel
# ============================================================

ml_complete <- read.csv(
  "file path/Observer final new.csv",
  stringsAsFactors = FALSE
)

required_columns <- c(
  "SBI",
  "mean_CAS",
  "ORI_bounded",
  "AviCRED"
)

missing_columns <- setdiff(
  required_columns,
  names(ml_complete)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Ensure analysis variables are numeric
numeric_columns <- c(
  "SBI",
  "mean_CAS",
  "ORI_bounded",
  "AviCRED"
)

for (v in numeric_columns) {
  ml_complete[[v]] <- as.numeric(
    as.character(ml_complete[[v]])
  )
}

# ------------------------------------------------------------
# 2. DEFINE AVICRED CREDIBILITY CLASSES
# ------------------------------------------------------------

# Retain these only if they match the intended classification.
class_breaks <- c(
  -Inf, 0.20, 0.40, 0.60, 0.80, Inf
)

class_labels <- c(
  "Very Low",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

ml_complete$AviCRED_class <- cut(
  ml_complete$AviCRED,
  breaks = class_breaks,
  labels = class_labels,
  include.lowest = TRUE,
  ordered_result = TRUE
)

# Check class distribution
cat("\nAVICRED CLASS COUNTS\n\n")

print(
  table(
    ml_complete$AviCRED_class,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# 3. DEFINE DARKER, HIGH-CONTRAST COLOUR PALETTE
# ------------------------------------------------------------

avicred_palette <- c(
  "Very Low"  = "#A50F15",  # Deep crimson
  "Low"       = "#E87500",  # Saturated orange
  "Moderate"  = "#7561A8",  # Purple
  "High"      = "#16817A",  # Deep teal
  "Very High" = "#225EA8"   # Strong blue
)

# ------------------------------------------------------------
# 4. DEFINE CORRELATION PAIRS
# ------------------------------------------------------------

correlation_pairs <- list(
  "A. SBI vs CAS" = c("SBI", "mean_CAS"),
  "B. CAS vs ORI" = c("mean_CAS", "ORI_bounded"),
  "C. ORI vs SBI" = c("ORI_bounded", "SBI")
)

# ------------------------------------------------------------
# 5. CALCULATE PEARSON AND SPEARMAN CORRELATIONS
# ------------------------------------------------------------

calculate_correlations <- function(
    data, x_var, y_var, pair_name) {
  
  valid <- is.finite(data[[x_var]]) &
    is.finite(data[[y_var]])
  
  x <- data[[x_var]][valid]
  y <- data[[y_var]][valid]
  
  n <- length(x)
  
  if (n < 3) {
    stop(
      "Insufficient observations for ",
      pair_name
    )
  }
  
  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    stop(
      "No variation in one or both variables for ",
      pair_name
    )
  }
  
  pearson_test <- cor.test(
    x,
    y,
    method = "pearson"
  )
  
  spearman_test <- cor.test(
    x,
    y,
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    Pair = pair_name,
    Variable_X = x_var,
    Variable_Y = y_var,
    N = n,
    Pearson_r = unname(
      pearson_test$estimate
    ),
    Pearson_CI_lower = pearson_test$conf.int[1],
    Pearson_CI_upper = pearson_test$conf.int[2],
    Pearson_p = pearson_test$p.value,
    Spearman_rho = unname(
      spearman_test$estimate
    ),
    Spearman_p = spearman_test$p.value,
    stringsAsFactors = FALSE
  )
}

correlation_results <- do.call(
  rbind,
  lapply(
    names(correlation_pairs),
    function(pair_name) {
      
      vars <- correlation_pairs[[pair_name]]
      
      calculate_correlations(
        data = ml_complete,
        x_var = vars[1],
        y_var = vars[2],
        pair_name = pair_name
      )
    }
  )
)

print(
  correlation_results,
  row.names = FALSE,
  digits = 5
)

write.csv(
  correlation_results,
  "AviCRED_Component_Correlations.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. CREATE EACH CORRELATION PANEL
# ------------------------------------------------------------

make_correlation_plot <- function(
    data,
    x_var,
    y_var,
    panel_title,
    correlation_row) {
  
  # Prepare panel data
  plot_df <- data.frame(
    x = data[[x_var]],
    y = data[[y_var]],
    AviCRED_class = data$AviCRED_class
  )
  
  # Retain valid observations
  valid <- is.finite(plot_df$x) &
    is.finite(plot_df$y) &
    !is.na(plot_df$AviCRED_class)
  
  plot_df <- plot_df[
    valid,
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_df) < 3) {
    stop(
      "Insufficient valid plotting observations for ",
      panel_title
    )
  }
  
  # Preserve ordered credibility classes
  plot_df$AviCRED_class <- factor(
    plot_df$AviCRED_class,
    levels = class_labels,
    ordered = TRUE
  )
  
  # Statistics annotation
  stats_label <- sprintf(
    "Pearson r = %.3f\nSpearman \u03C1 = %.3f\nN = %s",
    correlation_row$Pearson_r,
    correlation_row$Spearman_rho,
    format(
      correlation_row$N,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  )
  
  # Construct plot
  ggplot(
    plot_df,
    aes(
      x = x,
      y = y,
      fill = AviCRED_class
    )
  ) +
    
    # Observations
    geom_point(
      shape = 21,
      colour = "grey15",
      stroke = 0.12,
      size = 1.8,
      alpha = 0.70
    ) +
    
    # Overall linear regression line and 95% CI
    geom_smooth(
      aes(
        x = x,
        y = y
      ),
      inherit.aes = FALSE,
      data = plot_df,
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = "#D7191C",
      fill = "#D7191C",
      linewidth = 1.05,
      alpha = 0.13
    ) +
    
    # Credibility class colours
    scale_fill_manual(
      values = avicred_palette,
      breaks = class_labels,
      drop = FALSE,
      name = "AviCRED class"
    ) +
    
    # Correlation statistics
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = stats_label,
      hjust = 1.05,
      vjust = 1.10,
      size = 3.2,
      family = "serif",
      label.size = 0.2,
      fill = "white",
      alpha = 0.92
    ) +
    
    # Axis labels and title
    labs(
      title = panel_title,
      x = x_var,
      y = y_var
    ) +
    
    # Publication-style theme
    theme_classic(
      base_size = 12,
      base_family = "serif"
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12
      ),
      
      axis.title = element_text(
        face = "bold",
        size = 11
      ),
      
      axis.text = element_text(
        colour = "black",
        size = 9
      ),
      
      legend.title = element_text(
        face = "bold",
        size = 10
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      legend.position = "right",
      
      plot.margin = margin(
        8, 10, 8, 8
      )
    )
}

# ------------------------------------------------------------
# 7. GENERATE THE THREE PANELS
# ------------------------------------------------------------

p1 <- make_correlation_plot(
  ml_complete,
  "SBI",
  "mean_CAS",
  "A. SBI vs CAS",
  correlation_results[1, ]
)

p2 <- make_correlation_plot(
  ml_complete,
  "mean_CAS",
  "ORI_bounded",
  "B. CAS vs ORI",
  correlation_results[2, ]
)

p3 <- make_correlation_plot(
  ml_complete,
  "ORI_bounded",
  "SBI",
  "C. ORI vs SBI",
  correlation_results[3, ]
)

# ------------------------------------------------------------
# 8. COMBINE PANELS WITH SHARED LEGEND
# ------------------------------------------------------------

avicred_correlation_figure <-
  (p1 | p2 | p3) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

print(
  avicred_correlation_figure
)

# ------------------------------------------------------------
# 9. EXPORT PUBLICATION-QUALITY FIGURES
# ------------------------------------------------------------

ggsave(
  filename = "AviCRED_Component_Correlations.png",
  plot = avicred_correlation_figure,
  width = 15,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# # Close any open graphics devices
while (grDevices::dev.cur() > 1) {
  grDevices::dev.off()
}

# ------------------------------------------------------------
# 8. EXPORT FIGURES AND CORRELATION RESULTS
# ------------------------------------------------------------

output_dir <- "file path"

ggsave(
  filename = file.path(
    output_dir, "AviCRED_Component_Correlations.png"
  ),
  plot = avicred_correlation_figure,
  width = 15,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir, "AviCRED_Component_Correlations.pdf"
  ),
  plot = avicred_correlation_figure,
  width = 15,
  height = 5.8,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

write.csv(
  correlation_results,
  file.path(output_dir, "AviCRED_Component_Correlations.csv"),
  row.names = FALSE
)

# ============================================================
# REGRESSION ANALYSIS OF AVICRED COMPONENT RELATIONSHIPS
# Linear vs quadratic vs spline regression
# ============================================================

# ------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------

df <- read.csv(
  "file path/Observer final new.csv",
  stringsAsFactors = FALSE
)

output_dir <- "file path/AviCRED_regression_results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Use the existing component columns only.
required_columns <- c(
  "ORI_bounded",
  "mean_CAS",
  "SBI"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. DEFINE THE THREE PAIRWISE RELATIONSHIPS
# ------------------------------------------------------------

# Each model estimates the second variable as a function
# of the first variable.

model_pairs <- list(
  list(
    label = "ORI_vs_CAS",
    x = "ORI_bounded",
    y = "mean_CAS"
  ),
  list(
    label = "CAS_vs_ORI",
    x = "mean_CAS",
    y = "ORI_bounded"
  ),
  list(
    label = "ORI_vs_SBI",
    x = "ORI_bounded",
    y = "SBI"
  )
)

# ------------------------------------------------------------
# 3. STORAGE FOR RESULTS
# ------------------------------------------------------------

model_fit_results <- list()
coefficient_results <- list()
comparison_results <- list()
diagnostic_results <- list()
plots <- list()

# ------------------------------------------------------------
# 4. FIT MODELS FOR EACH PAIR
# ------------------------------------------------------------

for (pair in model_pairs) {
  
  x <- pair$x
  y <- pair$y
  pair_name <- pair$label
  
  cat("\n\n================================================\n")
  cat("RELATIONSHIP:", pair_name, "\n")
  cat("Response:", y, "| Predictor:", x, "\n")
  cat("================================================\n")
  
  # Use complete observations for this pair.
  dat <- df %>%
    select(all_of(c(x, y))) %>%
    filter(
      is.finite(.data[[x]]),
      is.finite(.data[[y]])
    )
  
  if (nrow(dat) < 30) {
    warning("Too few complete observations for ", pair_name)
    next
  }
  
  # ----------------------------------------------------------
  # A. LINEAR MODEL
  # ----------------------------------------------------------
  
  formula_linear <- reformulate(x, response = y)
  
  model_linear <- lm(
    formula_linear,
    data = dat
  )
  
  # ----------------------------------------------------------
  # B. QUADRATIC MODEL
  # ----------------------------------------------------------
  
  formula_quadratic <- as.formula(
    paste0(y, " ~ ", x, " + I(", x, "^2)")
  )
  
  model_quadratic <- lm(
    formula_quadratic,
    data = dat
  )
  
  # ----------------------------------------------------------
  # C. NATURAL SPLINE MODEL
  # ----------------------------------------------------------
  
  formula_spline <- as.formula(
    paste0(y, " ~ splines::ns(", x, ", df = 3)")
  )
  
  model_spline <- lm(
    formula_spline,
    data = dat
  )
  
  # ----------------------------------------------------------
  # D. MODEL SUMMARIES
  # ----------------------------------------------------------
  
  models <- list(
    Linear = model_linear,
    Quadratic = model_quadratic,
    Spline = model_spline
  )
  
  cat("\nLINEAR MODEL:\n")
  print(summary(model_linear))
  
  cat("\nQUADRATIC MODEL:\n")
  print(summary(model_quadratic))
  
  cat("\nSPLINE MODEL:\n")
  print(summary(model_spline))
  
  # ----------------------------------------------------------
  # E. COMPARE MODEL FIT
  # ----------------------------------------------------------
  
  fit_table <- bind_rows(
    lapply(names(models), function(model_name) {
      
      m <- models[[model_name]]
      sm <- summary(m)
      
      data.frame(
        Pair = pair_name,
        Response = y,
        Predictor = x,
        Model = model_name,
        N = nobs(m),
        R_squared = sm$r.squared,
        Adjusted_R_squared = sm$adj.r.squared,
        AIC = AIC(m),
        BIC = BIC(m),
        Residual_SE = sm$sigma,
        stringsAsFactors = FALSE
      )
    })
  )
  
  model_fit_results[[pair_name]] <- fit_table
  
  cat("\nMODEL FIT COMPARISON:\n")
  print(fit_table)
  
  # ----------------------------------------------------------
  # F. TEST NONLINEARITY
  # ----------------------------------------------------------
  
  # Does adding a quadratic term improve the linear model?
  quadratic_comparison <- anova(
    model_linear,
    model_quadratic
  )
  
  # Does the spline improve fit over a linear relationship?
  spline_comparison <- anova(
    model_linear,
    model_spline
  )
  
  cat("\nLINEAR VS QUADRATIC:\n")
  print(quadratic_comparison)
  
  cat("\nLINEAR VS SPLINE:\n")
  print(spline_comparison)
  
  comparison_results[[pair_name]] <- data.frame(
    Pair = pair_name,
    Comparison = c(
      "Linear vs Quadratic",
      "Linear vs Spline"
    ),
    F_statistic = c(
      quadratic_comparison$F[2],
      spline_comparison$F[2]
    ),
    P_value = c(
      quadratic_comparison$`Pr(>F)`[2],
      spline_comparison$`Pr(>F)`[2]
    ),
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
  # G. COEFFICIENTS AND 95% CONFIDENCE INTERVALS
  # ----------------------------------------------------------
  
  for (model_name in names(models)) {
    
    m <- models[[model_name]]
    
    coefficient_table <- as.data.frame(
      summary(m)$coefficients
    )
    
    coefficient_table$Term <- rownames(coefficient_table)
    rownames(coefficient_table) <- NULL
    
    ci <- as.data.frame(confint(m, level = 0.95))
    ci$Term <- rownames(ci)
    rownames(ci) <- NULL
    
    names(coefficient_table) <- c(
      "Estimate",
      "Std_Error",
      "Statistic",
      "P_value",
      "Term"
    )
    
    names(ci) <- c(
      "CI_Lower",
      "CI_Upper",
      "Term"
    )
    
    coefficient_table <- left_join(
      coefficient_table,
      ci,
      by = "Term"
    ) %>%
      mutate(
        Pair = pair_name,
        Model = model_name,
        Response = y,
        Predictor = x,
        .before = 1
      )
    
    coefficient_results[[
      paste(pair_name, model_name, sep = "_")
    ]] <- coefficient_table
  }
  
  # ----------------------------------------------------------
  # H. RESIDUAL DIAGNOSTICS
  # ----------------------------------------------------------
  
  # Diagnostics for the spline model, used as a flexible
  # candidate model. Also inspect linear/quadratic diagnostics
  # before selecting the final model.
  
  residual_values <- residuals(model_spline)
  fitted_values <- fitted(model_spline)
  
  diagnostic_results[[pair_name]] <- data.frame(
    Pair = pair_name,
    N = nobs(model_spline),
    Shapiro_W = if (length(residual_values) >= 3 &&
                    length(residual_values) <= 5000) {
      unname(shapiro.test(residual_values)$statistic)
    } else {
      NA_real_
    },
    Shapiro_P = if (length(residual_values) >= 3 &&
                    length(residual_values) <= 5000) {
      shapiro.test(residual_values)$p.value
    } else {
      NA_real_
    },
    Residual_Mean = mean(residual_values),
    Residual_SD = sd(residual_values),
    stringsAsFactors = FALSE
  )
  
  # Save standard diagnostic plots for the spline model.
  jpeg(
    filename = file.path(
      output_dir,
      paste0(pair_name, "_spline_diagnostics.jpg")
    ),
    width = 2400,
    height = 2400,
    res = 300
  )
  
  par(
    mfrow = c(2, 2),
    family = "serif"
  )
  
  plot(model_spline)
  
  # ----------------------------------------------------------
  # I. PREDICTION CURVES WITH 95% CONFIDENCE BANDS
  # ----------------------------------------------------------
  
  prediction_data <- data.frame(
    predictor_value = seq(
      min(dat[[x]]),
      max(dat[[x]]),
      length.out = 300
    )
  )
  
  names(prediction_data) <- x
  
  prediction_data$Linear <- predict(
    model_linear,
    newdata = prediction_data
  )
  
  prediction_data$Quadratic <- predict(
    model_quadratic,
    newdata = prediction_data
  )
  
  spline_prediction <- predict(
    model_spline,
    newdata = prediction_data,
    interval = "confidence",
    level = 0.95
  )
  
  prediction_data$Spline <- spline_prediction[, "fit"]
  prediction_data$Spline_lower <- spline_prediction[, "lwr"]
  prediction_data$Spline_upper <- spline_prediction[, "upr"]
  
  p <- ggplot(dat, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(
      alpha = 0.12,
      size = 0.8
    ) +
    geom_ribbon(
      data = prediction_data,
      aes(
        x = .data[[x]],
        ymin = Spline_lower,
        ymax = Spline_upper
      ),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    geom_line(
      data = prediction_data,
      aes(x = .data[[x]], y = Spline),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    labs(
      title = paste(y, "against", x),
      subtitle = "Natural spline regression with 95% confidence band",
      x = x,
      y = y
    ) +
    theme_classic(base_family = "serif") +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 16
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11
      ),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
  
  plots[[pair_name]] <- p
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(pair_name, "_spline_regression.jpg")
    ),
    plot = p,
    width = 8,
    height = 6,
    dpi = 400
  )
}

# ------------------------------------------------------------
# 5. COMBINE AND SAVE ALL RESULTS
# ------------------------------------------------------------

all_model_fits <- bind_rows(model_fit_results)
all_coefficients <- bind_rows(coefficient_results)
all_comparisons <- bind_rows(comparison_results)
all_diagnostics <- bind_rows(diagnostic_results)

write.csv(
  all_model_fits,
  file.path(output_dir, "Model_fit_comparison.csv"),
  row.names = FALSE
)

write.csv(
  all_coefficients,
  file.path(output_dir, "Regression_coefficients_and_CI.csv"),
  row.names = FALSE
)

write.csv(
  all_comparisons,
  file.path(output_dir, "Nonlinearity_model_comparisons.csv"),
  row.names = FALSE
)

write.csv(
  all_diagnostics,
  file.path(output_dir, "Spline_model_diagnostics.csv"),
  row.names = FALSE
)

cat("\n\n============================================\n")
cat("REGRESSION ANALYSIS COMPLETE\n")
cat("Results saved to:\n", output_dir, "\n")
cat("============================================\n")

cat("\nModel fit comparison:\n")
print(all_model_fits)

cat("\nTests of nonlinearity:\n")
print(all_comparisons)

cat("\nDiagnostic summary:\n")
print(all_diagnostics)

# ============================================================
# AVICRED: PAIRWISE COMPONENT CORRELATION AND DISCORDANCE
# Components: ORI, CAS, SBI
# ============================================================

# Check required variables
cas_var <- intersect(c("CAS", "mean_CAS"), names(ml_complete))[1]
required_vars <- c("ORI", cas_var, "SBI", "AviCRED", "avicred_class")
missing_vars <- setdiff(required_vars, names(ml_complete))

if (length(missing_vars) > 0 || is.na(cas_var)) {
  cat("Missing required variables:", paste(missing_vars, collapse = ", "), "\n")
} else {
  cat("All required variables are present. CAS column:", cas_var, "\n")
  print(str(ml_complete[, required_vars]))
  print(colSums(is.na(ml_complete[, required_vars])))
}

analysis_data <- ml_complete %>%
  select(
    observer_id, n_checklists, years_active,
    ORI_bounded, mean_CAS, SBI,
    AviCRED, avicred_class
  ) %>%
  filter(
    is.finite(ORI_bounded),
    is.finite(mean_CAS),
    is.finite(SBI),
    is.finite(AviCRED),
    !is.na(avicred_class)
  )

cat("Observers included:", nrow(analysis_data), "\n")

# Standardise component names for plotting
component_data <- analysis_data %>%
  transmute(
    observer_id,
    ORI = ORI_bounded,
    CAS = mean_CAS,
    SBI,
    AviCRED,
    avicred_class
  )

# Convert all three pairwise relationships to long format
pairwise_scatter_data <- bind_rows(
  component_data %>%
    transmute(
      observer_id, AviCRED,
      comparison = "CAS vs ORI",
      x = CAS, y = ORI
    ),
  component_data %>%
    transmute(
      observer_id, AviCRED,
      comparison = "CAS vs SBI",
      x = CAS, y = SBI
    ),
  component_data %>%
    transmute(
      observer_id, AviCRED,
      comparison = "ORI vs SBI",
      x = ORI, y = SBI
    )
)

# CAS vs ORI
cas_ori <- component_data %>%
  mutate(
    discordance_group = case_when(
      CAS >= quantile(CAS, 0.75) &
        ORI <= quantile(ORI, 0.25) ~ "High CAS - Low ORI",
      CAS <= quantile(CAS, 0.25) &
        ORI >= quantile(ORI, 0.75) ~ "Low CAS - High ORI",
      TRUE ~ NA_character_
    ),
    comparison = "CAS vs ORI"
  )

# CAS vs SBI
cas_sbi <- component_data %>%
  mutate(
    discordance_group = case_when(
      CAS >= quantile(CAS, 0.75) &
        SBI <= quantile(SBI, 0.25) ~ "High CAS - Low SBI",
      CAS <= quantile(CAS, 0.25) &
        SBI >= quantile(SBI, 0.75) ~ "Low CAS - High SBI",
      TRUE ~ NA_character_
    ),
    comparison = "CAS vs SBI"
  )

# ORI vs SBI
ori_sbi <- component_data %>%
  mutate(
    discordance_group = case_when(
      ORI >= quantile(ORI, 0.75) &
        SBI <= quantile(SBI, 0.25) ~ "High ORI - Low SBI",
      ORI <= quantile(ORI, 0.25) &
        SBI >= quantile(SBI, 0.75) ~ "Low ORI - High SBI",
      TRUE ~ NA_character_
    ),
    comparison = "ORI vs SBI"
  )

discordant_all <- bind_rows(cas_ori, cas_sbi, ori_sbi) %>%
  filter(!is.na(discordance_group)) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("CAS vs ORI", "CAS vs SBI", "ORI vs SBI")
    )
  )

print(
  discordant_all %>%
    count(comparison, discordance_group)
)

p_all_boxplot <- ggplot(
  discordant_all,
  aes(
    x = interaction(
      comparison, discordance_group,
      sep = "\n", lex.order = TRUE
    ),
    y = AviCRED,
    fill = comparison
  )
) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 1.8,
    outlier.alpha = 0.7,
    alpha = 0.65,
    linewidth = 0.6
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.3f", after_stat(y))),
    colour = "black",
    size = 4,
    family = "Times New Roman",
    fontface = "bold",
    vjust = -0.7
  ) +
  labs(
    x = "Component Comparison and Discordance Group",
    y = "AviCRED",
    fill = "Component Comparison"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 14,
      colour = "black"
    ),
    plot.title = element_text(
      family = "Times New Roman",
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  )

print(p_all_boxplot)

ggsave(
  filename = "AviCRED_Component_Discordance_Boxplot.jpg",
  plot = p_all_boxplot,
  device = "jpeg",
  dpi = 400,
  width = 15,
  height = 8,
  units = "in",
  quality = 100,
  bg = "white"
)

# ================== ISOLATION FOREST VS AVICRED STACKED ================= #

# ============================================================
# CHECKLIST-LEVEL CAS STACKED
# ============================================================

library(tidyverse)
library(scales)

# 1. LOAD CHECKLIST CSV
checklists <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. CHECK REQUIRED COLUMNS
required_columns <- c(
  "checklist_id",
  "CAS_class",
  "if_top10_typicality"
)

missing_columns <- setdiff(
  required_columns,
  names(checklists)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# 3. CLEAN DATA AND DEFINE TYPICALITY GROUPS
cas_levels <- c(
  "Very Low",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

checklists <- checklists %>%
  filter(
    !is.na(CAS_class),
    !is.na(if_top10_typicality)
  ) %>%
  mutate(
    CAS_class = factor(
      CAS_class,
      levels = cas_levels,
      ordered = TRUE
    ),
    IF_Typicality = case_when(
      if_top10_typicality == 1 ~ "Top 10% typicality",
      if_top10_typicality == 0 ~ "Remaining 90%",
      TRUE ~ NA_character_
    ),
    IF_Typicality = factor(
      IF_Typicality,
      levels = c(
        "Top 10% typicality",
        "Remaining 90%"
      )
    )
  ) %>%
  filter(!is.na(IF_Typicality))

# 4. COUNT CAS CLASSES WITHIN EACH IF GROUP
plot_df <- checklists %>%
  count(
    IF_Typicality,
    CAS_class,
    .drop = FALSE
  ) %>%
  group_by(IF_Typicality) %>%
  mutate(
    total_group = sum(n),
    Proportion = n / total_group,
    Percentage = 100 * Proportion,
    label = ifelse(
      Proportion >= 0.035,
      paste0(round(Percentage, 1), "%"),
      ""
    )
  ) %>%
  ungroup() %>%
  filter(n > 0)

# 5. GROUP TOTALS FOR SAMPLE-SIZE LABELS
total_df <- checklists %>%
  count(IF_Typicality, name = "total")

print(total_df)

windowsFonts(
  Times = windowsFont("Times New Roman")
)

plot_df <- plot_df %>%
  mutate(
    IF_Typicality = factor(
      IF_Typicality,
      levels = c("Top 10% typicality", "Remaining 90%")
    ),
    CAS_class = factor(
      CAS_class,
      levels = c(
        "Very Low",
        "Low",
        "Moderate",
        "High",
        "Very High"
      ),
      ordered = TRUE
    ),
    label = ifelse(
      Proportion >= 0.035,
      scales::percent(Proportion, accuracy = 0.1),
      ""
    )
  )

total_df <- total_df %>%
  mutate(
    IF_Typicality = factor(
      IF_Typicality,
      levels = c("Top 10% typicality", "Remaining 90%")
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x = IF_Typicality,
    y = Proportion,
    fill = CAS_class
  )
) +
  
  # Stacked bars
  geom_col(
    colour = "black",
    linewidth = 0.25,
    width = 0.65
  ) +
  
  # Percentage labels inside stacked bars
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    family = "Times",
    fontface = "plain",
    size = 4,
    colour = "black"
  ) +
  
  # Total sample size above each bar
  geom_text(
    data = total_df,
    aes(
      x = IF_Typicality,
      y = 1.045,
      label = paste0("n = ", scales::comma(total))
    ),
    inherit.aes = FALSE,
    family = "Times",
    fontface = "bold",
    size = 4,
    colour = "black"
  ) +
  
  # CAS class colours
  scale_fill_manual(
    values = c(
      "Very Low"  = "#B2182B",
      "Low"       = "#EF8A62",
      "Moderate"  = "#FDDBC7",
      "High"      = "#A6DBA0",
      "Very High" = "#1B7837"
    ),
    drop = FALSE
  ) +
  
  # Y-axis
  scale_y_continuous(
    labels = scales::label_percent(),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_cartesian(
    ylim = c(0, 1.12),
    clip = "off"
  ) +
  
  # Titles and labels
  labs(
    x = NULL,
    y = "Percentage of Checklists",
    fill = "CAS Class"
  ) +
  
  # Base theme
  theme_minimal(base_family = "Times", base_size = 12) +
  
  # Explicitly set Times for every text element
  theme(
    text = element_text(
      family = "Times",
      colour = "black"
    ),
    
    plot.title = element_text(
      family = "Times",
      face = "bold",
      size = 14,
      hjust = 0.5,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      family = "Times",
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      family = "Times",
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    axis.text.x = element_text(
      family = "Times",
      face = "bold",
      size = 11,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      family = "Times",
      face = "plain",
      size = 11,
      colour = "black"
    ),
    
    legend.title = element_text(
      family = "Times",
      face = "bold",
      size = 11,
      colour = "black"
    ),
    
    legend.text = element_text(
      family = "Times",
      face = "plain",
      size = 10,
      colour = "black"
    ),
    
    legend.position = "bottom",
    
    plot.margin = margin(15, 10, 10, 10)
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

print(p)

jpeg(
  filename = "CAS_Class_vs_IF_Typicality.jpg",
  width = 4000,
  height = 2800,
  res = 400,
  quality = 100
)

#================= OBSERVER ISOLATION FOREST STACKED ===============#

ml_complete <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Loaded observer-level dataset\n")
cat("Rows:", nrow(ml_complete), "\n")
cat("Columns:", ncol(ml_complete), "\n")

required_columns <- c(
  "avicred_class",
  "ml_10pct"
)

missing_columns <- setdiff(
  required_columns,
  names(ml_complete)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. CLEAN DATA
# ------------------------------------------------------------

ml_complete <- ml_complete %>%
  filter(
    !is.na(avicred_class),
    !is.na(ml_10pct),
    ml_10pct %in% c(0, 1)
  )

ml_complete$avicred_class <- factor(
  ml_complete$avicred_class,
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  ordered = TRUE
)

ml_complete <- ml_complete %>%
  filter(!is.na(avicred_class))

# ------------------------------------------------------------
# 3. DEFINE TWO ISOLATION FOREST GROUPS
# ------------------------------------------------------------

ml_complete <- ml_complete %>%
  mutate(
    IF_Class = factor(
      if_else(
        ml_10pct == 1,
        "Top 10%",
        "Remaining 90%"
      ),
      levels = c(
        "Top 10%",
        "Remaining 90%"
      )
    )
  )

# ------------------------------------------------------------
# 4. CALCULATE CLASS PROPORTIONS
# ------------------------------------------------------------

if_df_obs <- ml_complete %>%
  count(IF_Class, avicred_class) %>%
  group_by(IF_Class) %>%
  mutate(
    Proportion = n / sum(n)
  ) %>%
  ungroup() %>%
  rename(Group = IF_Class)

# Sample sizes
if_totals_obs <- ml_complete %>%
  count(IF_Class, name = "total") %>%
  rename(Group = IF_Class)

print(if_totals_obs)

windowsFonts(
  Times = windowsFont("Times New Roman")
)

# ------------------------------------------------------------
# 5. CREATE PLOT
# ------------------------------------------------------------

p <- ggplot(
  if_df_obs,
  aes(
    x = Group,
    y = Proportion,
    fill = avicred_class
  )
) +
  
  geom_col(
    colour = "black",
    linewidth = 0.2,
    width = 0.70
  ) +
  
  # Percentage labels inside the stacked bars
  geom_text(
    aes(
      label = ifelse(
        Proportion >= 0.035,
        paste0(round(Proportion * 100, 1), "%"),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    family = "Times",
    size = 3.5,
    colour = "black"
  ) +
  
  # Sample sizes above the bars
  geom_text(
    data = if_totals_obs,
    aes(
      x = Group,
      y = 1.045,
      label = scales::comma(total)
    ),
    inherit.aes = FALSE,
    family = "Times",
    fontface = "bold",
    size = 4,
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Very Low"  = "#FDD0A2",
      "Low"       = "#FC9272",
      "Moderate"  = "#FB6A4A",
      "High"      = "#807DBA",
      "Very High" = "#41B6C4"
    ),
    drop = FALSE
  ) +
  
  scale_y_continuous(
    labels = scales::label_percent(),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_cartesian(
    ylim = c(0, 1.10),
    clip = "off"
  ) +
  
  labs(
    title = "Distribution of AviCRED Classes: Isolation Forest Top 10% vs Remaining 90%",
    x = NULL,
    y = "Percentage of Observers",
    fill = "AviCRED Class"
  ) +
  
  theme_minimal(
    base_family = "Times",
    base_size = 12
  ) +
  
  theme(
    text = element_text(
      family = "Times",
      colour = "black"
    ),
    
    plot.title = element_text(
      family = "Times",
      face = "bold",
      size = 14,
      hjust = 0.5,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      family = "Times",
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    axis.text.x = element_text(
      family = "Times",
      face = "bold",
      size = 11,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      family = "Times",
      size = 10,
      colour = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      family = "Times",
      face = "bold",
      size = 11,
      colour = "black"
    ),
    
    legend.text = element_text(
      family = "Times",
      size = 10,
      colour = "black"
    ),
    
    plot.margin = margin(15, 12, 10, 12)
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

print(p)

jpeg(
  filename = "AviCRED_Classes_IF_Top10_vs_Remaining90.jpg",
  width = 4200,
  height = 2800,
  res = 400,
  quality = 100
)

# ================= UPSET PLOT ===============#

# ============================================================
# CHECKLIST IF TOP 10% TYPICALITY VS CAS VERY HIGH
# ============================================================

library(scales)
library(grid)

# ------------------------------------------------------------
# 1. CALCULATE INTERSECTION COUNTS
# ------------------------------------------------------------

upset_data1 <- data.frame(
  IF_Top10 = checklists$if_top10_typicality == 1,
  CAS_VeryHigh = checklists$CAS_high == 1
) %>%
  filter(
    !is.na(IF_Top10),
    !is.na(CAS_VeryHigh)
  ) %>%
  mutate(
    Intersection = case_when(
      IF_Top10 & CAS_VeryHigh ~ "Both",
      IF_Top10 & !CAS_VeryHigh ~ "IF only",
      !IF_Top10 & CAS_VeryHigh ~ "CAS only",
      TRUE ~ "Neither"
    )
  )

# Traditional UpSet plots omit the "Neither" combination
intersection_counts <- upset_data1 %>%
  filter(Intersection != "Neither") %>%
  count(Intersection, name = "Count") %>%
  filter(Count > 0) %>%
  arrange(desc(Count))

intersection_order <- intersection_counts$Intersection

intersection_counts$Intersection <- factor(
  intersection_counts$Intersection,
  levels = intersection_order
)

# ------------------------------------------------------------
# 2. MEMBERSHIP MATRIX
# ------------------------------------------------------------

matrix_data <- expand.grid(
  Intersection = intersection_order,
  Set = c(
    "IF Top 10% typicality",
    "CAS Very High"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Present = case_when(
      Intersection == "Both" ~ TRUE,
      Intersection == "IF only" &
        Set == "IF Top 10% typicality" ~ TRUE,
      Intersection == "CAS only" &
        Set == "CAS Very High" ~ TRUE,
      TRUE ~ FALSE
    ),
    Intersection = factor(
      Intersection,
      levels = intersection_order
    ),
    Set = factor(
      Set,
      levels = c(
        "IF Top 10% typicality",
        "CAS Very High"
      )
    )
  )

windowsFonts(Times = windowsFont("Times New Roman"))

# Intersection-size bar plot
p_bars <- p_bars +
  theme_classic(base_family = "Times", base_size = 14) +
  theme(
    plot.title = element_text(
      family = "Times", size = 14, face = "bold",
      hjust = 0.5, colour = "black"
    ),
    axis.title = element_text(
      family = "Times", size = 14, face = "bold",
      colour = "black"
    ),
    axis.text = element_text(
      family = "Times", size = 14, colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

# Intersection matrix
p_matrix <- p_matrix +
  theme_classic(base_family = "Times", base_size = 14) +
  theme(
    axis.title = element_text(
      family = "Times", size = 14, face = "bold",
      colour = "black"
    ),
    axis.text = element_text(
      family = "Times", size = 14, face = "bold",
      colour = "black"
    ),
    axis.line = element_blank(),
    axis.ticks.x = element_blank()
  )

# Reassemble the UpSet plot
p_upset <- p_bars / p_matrix +
  patchwork::plot_layout(heights = c(3, 1.25))

print(p_upset)

# Save at publication resolution
ggsave(
  "Checklist upset.jpg",
  plot = p_upset,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

#=========== OBSERVER UPSET =================#

# ------------------------------------------------------------
# 1. PREPARE DATA (ml_complete IS ALREADY LOADED)
# ------------------------------------------------------------

upset_data2 <- data.frame(
  IF_Top5 = !is.na(ml_complete$ml_5pct) &
    ml_complete$ml_5pct == 1,
  
  AviCRED_VeryHigh = !is.na(ml_complete$avicred_class) &
    ml_complete$avicred_class == "Very High"
)

colnames(upset_data2) <- c(
  "Top 5% Isolation Forest Score",
  "AviCRED Class Very High"
)

upset_data2 <- upset_data2[
  rowSums(upset_data2) > 0,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 2. CALCULATE INTERSECTION COUNTS
# ------------------------------------------------------------

intersection_counts <- data.frame(
  Intersection = c(
    "Both",
    "IF only",
    "AviCRED only"
  ),
  Count = c(
    sum(upset_data2[[1]] & upset_data2[[2]]),
    sum(upset_data2[[1]] & !upset_data2[[2]]),
    sum(!upset_data2[[1]] & upset_data2[[2]])
  )
) %>%
  filter(Count > 0) %>%
  arrange(desc(Count))

intersection_order <- intersection_counts$Intersection

intersection_counts$Intersection <- factor(
  intersection_counts$Intersection,
  levels = intersection_order
)

# ------------------------------------------------------------
# 3. PREPARE INTERSECTION MATRIX
# ------------------------------------------------------------

matrix_data <- expand.grid(
  Intersection = intersection_order,
  Set = c(
    "Top 5% Isolation Forest Score",
    "AviCRED Class Very High"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Present = case_when(
      Intersection == "Both" ~ TRUE,
      Intersection == "IF only" &
        Set == "Top 5% Isolation Forest Score" ~ TRUE,
      Intersection == "AviCRED only" &
        Set == "AviCRED Class Very High" ~ TRUE,
      TRUE ~ FALSE
    ),
    Intersection = factor(
      Intersection,
      levels = intersection_order
    ),
    Set = factor(
      Set,
      levels = c(
        "AviCRED Class Very High",
        "Top 5% Isolation Forest Score"
      )
    )
  )

windowsFonts(Times = windowsFont("Times New Roman"))

# ------------------------------------------------------------
# 5. INTERSECTION-SIZE BARS
# ------------------------------------------------------------

p_bars <- ggplot(
  intersection_counts,
  aes(x = Intersection, y = Count)
) +
  geom_col(
    fill = "grey50",
    colour = "black",
    linewidth = 0.4,
    width = 0.68
  ) +
  geom_text(
    aes(label = scales::comma(Count)),
    family = "Times",
    size = 14 / ggplot2::.pt,
    colour = "black",
    vjust = -0.45
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    title = "Overlap between Isolation Forest and AviCRED",
    y = "Intersection size",
    x = NULL
  ) +
  theme_classic(
    base_family = "Times",
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      family = "Times",
      size = 14,
      face = "bold",
      hjust = 0.5,
      colour = "black"
    ),
    axis.title.y = element_text(
      family = "Times",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      family = "Times",
      size = 14,
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

# ------------------------------------------------------------
# 6. SET-INTERSECTION MATRIX
# ------------------------------------------------------------

p_matrix <- ggplot(
  matrix_data,
  aes(x = Intersection, y = Set)
) +
  geom_point(
    shape = 21,
    size = 4.5,
    stroke = 0.5,
    fill = "grey85",
    colour = "grey55"
  ) +
  geom_point(
    data = filter(matrix_data, Present),
    shape = 16,
    size = 3.5,
    colour = "black"
  ) +
  geom_segment(
    data = data.frame(
      Intersection = factor(
        "Both",
        levels = intersection_order
      ),
      x = match("Both", intersection_order)
    ) %>%
      filter(!is.na(x)),
    aes(
      x = x,
      xend = x,
      y = 1,
      yend = 2
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.6
  ) +
  scale_y_discrete(
    limits = c(
      "AviCRED Class Very High",
      "Top 5% Isolation Forest Score"
    )
  ) +
  labs(
    x = "Intersection",
    y = NULL
  ) +
  theme_classic(
    base_family = "Times",
    base_size = 14
  ) +
  theme(
    axis.title.x = element_text(
      family = "Times",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_text(
      family = "Times",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      family = "Times",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.line = element_blank(),
    axis.ticks.x = element_blank()
  )

# ------------------------------------------------------------
# 7. COMBINE INTO UPSET PLOT
# ------------------------------------------------------------

p_upset_observer <- p_bars / p_matrix +
  plot_layout(heights = c(3, 1.5))

print(p_upset_observer)

# ------------------------------------------------------------
# 8. SAVE PUBLICATION-QUALITY OUTPUTS
# ------------------------------------------------------------

ggsave(
  "Observer_UpSet_Monochrome.jpg",
  plot = p_upset_observer,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


#=================== OBSERVER METRICS MAP =================#

con <- dbConnect(duckdb())

ebd_gis <- dbGetQuery(
  con,
  "
SELECT DISTINCT
    \"SAMPLING EVENT IDENTIFIER\" AS checklist_id,
    \"OBSERVER ID\" AS observer_id,
    CAST(\"OBSERVATION DATE\" AS DATE) AS obs_date,
    \"LATITUDE\" AS latitude,
    \"LONGITUDE\" AS longitude,
    \"STATE\" AS state
FROM read_csv_auto(
    'FILE PATH/ebd_filtered.txt',
    delim='\t'
)
"
)

india_shp <- st_read(   #INSERT RESPECTIVE SHAPEFILE AFTER DOWNLOADING
  "file path/India_State_Boundary.shp"
)

dbDisconnect(con, shutdown = TRUE)

gis_data <- ebd_gis %>%
  left_join(
    ml_complete %>%
      select(
        observer_id,
        AviCRED,
        avicred_class,
        cluster,
        ml_score,
        ml_5pct,
        ml_10pct
      ),
    by = c("observer_id")
  )

View(gis_data)

gis_data <- gis_data %>%   #PREPARE MAPPING VARIABLES TO PLOTTED
  mutate(
    AviCRED_group = factor(
      avicred_class,
      levels = c("Very Low", "Low", "Moderate", "High", "Very High")
    ),
    if_group = case_when(
      ml_score >= quantile(ml_score, 0.95, na.rm = TRUE) ~ "Top 5%",
      ml_score >= quantile(ml_score, 0.90, na.rm = TRUE) ~ "Top 10%",
      TRUE ~ "Rest"
    ),
    cluster_group = factor(cluster)
  )

gis_sf <- st_as_sf(  #CREATE SPATIAL COORDINATES
  gis_data,
  coords = c("longitude", "latitude"),
  crs = 4326   #CRS EPSG:4326
)

india <- st_transform(india_shp, 3857)  #CHANGE CRS TO EPSG:3857 (distance measured in meters)

gis_sf <- st_transform(gis_sf, 3857) #CHANGE CRS TO EPSG:3857 (distance measured in meters)

grid25 <- st_make_grid(
  india,
  cellsize = 25000,   # MAKE 25 X 25KM GRIDS
  square = TRUE #SQUARE GRID CELLS
)

grid25 <- st_sf(  #TRANSFORM TOSF OBJECT AND ASSIGN UNIQUE ID TO GRIDS
  grid_id = seq_along(grid25),
  geometry = grid25
)

gis_join <- st_join(
  gis_sf,
  grid25
)

observer_grid <- gis_join %>%
  st_drop_geometry() %>%
  distinct(
    grid_id,
    observer_id,
    avicred_class,
    ml_5pct,
    ml_10pct,
    cluster
  )

obs_birdcred <- observer_grid %>%
  distinct(grid_id, #AVOIDS COUNTING DUPLICATE OBSERVERS PER GRID
           observer_id, 
           avicred_class) %>% 
  count(grid_id, avicred_class, name = "n") %>%  #DOMINANT AVICRED CATEGORY DISPLAY PER GRID
  group_by(grid_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

avicred_map <- left_join(grid25, obs_birdcred, by = "grid_id")

obs_ml <- observer_grid %>%
  distinct(   #AVOIDS COUNTING DUPLICATE OBSERVERS PER GRID
    grid_id,
    observer_id,
    ml_5pct,
    ml_10pct
  ) %>%
  mutate(
    ml_class = case_when(
      ml_5pct == 1 ~ "Top 5%",
      ml_10pct == 1 ~ "Top 10%",
      TRUE ~ "Remaining 90%"
    )
  ) %>%
  count(grid_id, ml_class, name = "n") %>%  #DOMINANT IF CLASSIFICATION DISPLAY PER GRID
  group_by(grid_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

ml_map <- left_join(grid25, obs_ml, by = "grid_id")

obs_cluster <- observer_grid %>%
  distinct(grid_id,  #AVOIDS COUNTING DUPLICATE OBSERVERS PER GRID 
           observer_id, 
           cluster
  ) %>% 
  count(grid_id, cluster, name = "n") %>% #DOMINANT CLUSTER CATEGORY DISPLAY PER GRID
  group_by(grid_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

cluster_map <- left_join(grid25, obs_cluster, by = "grid_id")

##################### AVICRED MAP ###############

ggplot() +
  geom_sf(
    data = avicred_map,
    aes(fill = avicred_class),
    colour = NA
  ) +
  geom_sf(
    data = india_shp,
    fill = NA,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Very Low" = "#EF2917",
      "Low" = "#A40E4C",
      "Moderate" = "#E6AF2E",
      "High" = "#307351",
      "Very High" = "#A5FFD6"
    ),
    na.translate = FALSE,
    name = "AviCRED Class"
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title = "AviCRED Class Distribution"
  ) +
  theme_bw() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    legend.title = element_text(
      face = "bold",
      family = "Times New Roman",
      size = 12
    ),
    legend.text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white"
    )
  )

################# ISOLATION MAP ###############

ggplot() +
  geom_sf(
    data = ml_map,
    aes(fill = ml_class),
    colour = NA
  ) +
  geom_sf(
    data = india_shp,
    fill = NA,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Top 5%" = "#132A13",
      "Top 10%" = "#E6AF2E",
      "Remaining 90%" = "#7FB685"
    ),
    na.translate = FALSE,
    name = "Isolation Forest"
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title = "Isolation Forest Distribution"
  ) +
  theme_bw() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    legend.title = element_text(
      face = "bold",
      family = "Times New Roman",
      size = 12
    ),
    legend.text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white"
    )
  )

############### CHECKLIST MAP ###############

gis_data_chk <- ebd_gis %>%
  left_join(
    checklists %>%
      select(
        checklist_id,
        CAS,
        CAS_class,
        cluster,
        if_authenticity,
        if_5pct,
        if_10pct
      ),
    by = "checklist_id"
  )

gis_data_chk <- gis_data_chk %>%  #PREPARE MAPPING VARIABLES TO BE PLOTTED
  mutate(
    CAS_group = factor(
      CAS_class,
      levels = c(
        "Very Low",
        "Low",
        "Moderate",
        "High",
        "Very High"
      ),
      ordered = TRUE
    ),
    if_group = case_when(
      if_5pct == 1 ~ "Top 5%",
      if_10pct == 1 ~ "Top 10%",
      TRUE ~ "Remaining 90%"
    ),
    if_group = factor(
      if_group,
      levels = c(
        "Top 5%",
        "Top 10%",
        "Remaining 90%"
      )
    ),
    cluster_group = factor(cluster)
  )

gis_sf <- st_as_sf(
  gis_data_chk,
  coords = c("longitude", "latitude"),
  crs = 4326  #CRS EPSG:4326
)

gis_sf <- st_transform(gis_sf, 3857) #CHANGE CRS TO EPSG:3857 (distance measured in meters)

gis_join_chk <- st_join(
  gis_sf,
  grid25
)

############## CAS MAP  #############

chk_cas <- gis_join_chk %>%
  st_drop_geometry() %>%
  distinct(  #AVOID COUNTING DUPLICATE CHECKLISTS PER GRID
    grid_id,
    checklist_id,
    CAS_class
  ) %>%
  count(   #DOMINANT CAS CATEGORY DISPLAY PER GRID
    grid_id,
    CAS_class,
    name = "n"
  ) %>%
  group_by(grid_id) %>%
  slice_max(
    n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

cas_map <- left_join(
  grid25,
  chk_cas,
  by = "grid_id"
)

cas_map <- cas_map %>%
  mutate(
    CAS_class = factor(
      CAS_class,
      levels = c(
        "Very Low",
        "Low",
        "Moderate",
        "High",
        "Very High"
      )
    )
  )

ggplot() +
  geom_sf(
    data = cas_map,
    aes(fill = CAS_class),
    colour = NA
  ) +
  geom_sf(
    data = india_shp,
    fill = NA,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Very Low" = "#020887",
      "Low" = "#B9D6F2",
      "Moderate" = "#F4B860",
      "High" = "#93E5AB",
      "Very High" = "#60992D"
    ),
    na.translate = FALSE,
    name = "CAS Class"
  ) +
  coord_sf(expand = FALSE) +
  labs(title = "Checklist CAS Distribution") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

################## IF SCORE MAP ###########

chk_ml <- gis_join_chk %>%
  st_drop_geometry() %>%
  distinct(  #AVOID COUNTING DUPLICATE CHECKLISTS PER GRID
    grid_id,
    checklist_id,
    if_5pct,
    if_10pct
  ) %>%
  mutate(
    ml_class = case_when(
      if_5pct == 1 ~ "Top 5%",
      if_10pct == 1 ~ "Top 10%",
      TRUE ~ "Remaining 90%"
    ),
    ml_class = factor(
      ml_class,
      levels = c(
        "Top 5%",
        "Top 10%",
        "Remaining 90%"
      )
    )
  ) %>%
  count(  #DOMINANT IF TYPICALITY CATEGORY DISPLAY PER GRID
    grid_id,
    ml_class,
    name = "n"
  ) %>%
  group_by(grid_id) %>%
  slice_max(
    n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

ml_map <- left_join(
  grid25,
  chk_ml,
  by = "grid_id"
)

ggplot() +
  geom_sf(
    data = ml_map,
    aes(fill = ml_class),
    colour = NA
  ) +
  geom_sf(
    data = india_shp,
    fill = NA,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Top 5%" = "#020887",
      "Top 10%" = "#345830",
      "Remaining 90%" = "#8CB369"
    ),
    na.translate = FALSE,
    name = "IF Authenticity Class"
  ) +
  coord_sf(expand = FALSE) +
  labs(title = "Spatial Distribution of IF Authenticity Classes") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


#============ FLAGGED RECORDS DISTRIBUTION ACROSS CHECKLISTS AND OBSERVER METRICS ===========#

df <- file.choose()

df <- read.csv(
  df,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

chk_df <- file.choose()

chk_df <- read.csv(
  chk_df,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(chk_df)

# ============================================================
# SET THE CORRECT ORDER OF AVICRED AND CAS CLASSES
# ============================================================

# AviCRED classes: lowest to highest credibility
df$avicred_class <- factor(
  df$avicred_class,
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  ordered = TRUE
)

# CAS classes: lowest to highest typicality
chk_df$CAS_class <- factor(
  chk_df$CAS_class,
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  ordered = TRUE
)

chk_df <- chk_df %>%  # Keep variables required for the anomaly analysis
  select(
    checklist_id,
    observer_id,
    CAS_class,
    IF_top10,
    richness_flag,
    abundance_flag
  )

chk_df <- chk_df %>%  # Classify each checklist according to its detected anomaly type
  mutate(
    anomaly_type = case_when(
      richness_flag == 0 & abundance_flag == 0 ~ "None",
      richness_flag == 1 & abundance_flag == 0 ~ "Richness",
      richness_flag == 0 & abundance_flag == 1 ~ "Abundance",
      TRUE ~ "Both"
    )
  )

chk_df$anomaly_type <- factor(  # Set the order of checklist anomaly categories
  chk_df$anomaly_type,
  levels = c(
    "None",
    "Richness",
    "Abundance",
    "Both"
  )
)

observer_flags <- chk_df %>%  # Count the different anomaly types for each observer
  group_by(observer_id) %>%
  summarise(
    n_none = sum(anomaly_type == "None"),
    n_richness = sum(anomaly_type == "Richness"),
    n_abundance = sum(anomaly_type == "Abundance"),
    n_both = sum(anomaly_type == "Both"),
    n_checklists = n(),
    .groups = "drop"
  )

df <- df %>%  # Add observer-level anomaly counts to the AviCRED dataset
  left_join(
    observer_flags,
    by = "observer_id"
  )

df <- df %>%  # Assign each observer an overall anomaly type based on their checklist history
  mutate(
    observer_anomaly_type = case_when(
      n_richness == 0 &
        n_abundance == 0 &
        n_both == 0 ~ "None",
      n_both >= pmax(
        n_richness,
        n_abundance
      ) &
        n_both > 0 ~ "Both",
      n_richness >= n_abundance &
        n_richness > 0 ~ "Richness",
      TRUE ~ "Abundance"
    )
  )

df$observer_anomaly_type <- factor(  # Set the order of observer anomaly categories
  df$observer_anomaly_type,
  levels = c(
    "None",
    "Richness",
    "Abundance",
    "Both"
  )
)

obs_avicred <- df %>%  # Calculate observer anomaly proportions within each AviCRED class
  filter(!is.na(avicred_class)) %>%
  count(avicred_class, observer_anomaly_type) %>%
  group_by(avicred_class) %>%
  mutate(prop = n / sum(n))

obs_if <- df %>%  # Classify observers according to Isolation Forest thresholds
  filter(
    !is.na(ml_10pct)
  ) %>%
  mutate(
    IF = case_when(
      ml_10pct == 1 ~ "Top 10%",
      TRUE ~ "Remaining 90%"
    ),
    IF = factor(
      IF,
      levels = c(
        "Top 10%",
        "Remaining 90%"
      )
    )
  ) %>%
  count(IF, observer_anomaly_type) %>%
  group_by(IF) %>%
  mutate(prop = n / sum(n))

cas_plot <- chk_df %>%  # Calculate checklist anomaly proportions within each CAS class
  count(CAS_class, anomaly_type) %>%
  group_by(CAS_class) %>%
  mutate(prop = n / sum(n))

checklists_if <- chk_df %>%  # Prepare checklist data for Isolation Forest classification
  filter(
    !is.na(IF_top10)
  )

chk_if <- checklists_if %>%   # Calculate checklist anomaly proportions within Isolation Forest groups
  mutate(
    IF = case_when(
      IF_top10 == 1 ~ "Top 10%",
      TRUE ~ "Remaining 90%"
    ),
    IF = factor(
      IF,
      levels = c(
        "Top 10%",
        "Remaining 90%"
      )
    )
  ) %>%
  count(IF, anomaly_type) %>%
  group_by(IF) %>%
  mutate(prop = n / sum(n))

cols <- c(  # Define colours for the four anomaly categories
  "None" = "grey80",
  "Richness" = "#1f78b4",
  "Abundance" = "#ff7f00",
  "Both" = "#d62728"
)

p1 <- ggplot(
  obs_avicred,
  aes(
    x = avicred_class,
    y = prop,
    fill = observer_anomaly_type
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = ifelse(
        prop > 0.03,
        scales::percent(prop, accuracy = 0.1),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    family = "Times",
    fontface = "bold",
    size = 3.2,
    colour = "black"
  ) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(
    values = cols,
    name = "Observer Anomaly Type"
  ) +
  labs(
    title = "AviCRED Classes",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(base_family = "Times") +
  theme(
    text = element_text(
      family = "Times"
    ),
    plot.title = element_text(
      family = "Times",
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      family = "Times",
      size = 10
    ),
    axis.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      family = "Times",
      size = 10
    ),
    legend.position = "bottom"
  )

p2 <- ggplot(
  obs_if,
  aes(
    x = IF,
    y = prop,
    fill = observer_anomaly_type
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = ifelse(
        prop > 0.03,
        scales::percent(prop, accuracy = 0.1),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    family = "Times",
    fontface = "bold",
    size = 3.2,
    colour = "black"
  ) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(
    values = cols,
    name = "Observer Anomaly Type"
  ) +
  labs(
    title = "Observer IF Anomaly Classes",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(base_family = "Times") +
  theme(
    text = element_text(
      family = "Times"
    ),
    plot.title = element_text(
      family = "Times",
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      family = "Times",
      size = 10
    ),
    axis.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      family = "Times",
      size = 10
    ),
    legend.position = "bottom"
  )

p3 <- ggplot(
  cas_plot,
  aes(
    x = CAS_class,
    y = prop,
    fill = anomaly_type
  )
) +
  geom_col() +
  geom_text(
    aes(label = ifelse(prop > 0.03,
                       scales::percent(prop, accuracy = 0.1),
                       ""
    )),
    position = position_stack(vjust = 0.5),
    family = "Times",
    fontface = "bold",
    size = 3.2,
    colour = "black"
  ) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(
    values = cols,
    name = "Checklist Anomaly Type"
  ) +
  labs(
    title = "CAS Classes",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(base_family = "Times") +
  theme(
    text = element_text(
      family = "Times"
    ),
    plot.title = element_text(
      family = "Times",
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      family = "Times",
      size = 10
    ),
    axis.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      family = "Times",
      size = 10
    ),
    legend.position = "bottom"
  )

p4 <- ggplot(
  chk_if,
  aes(
    x = IF,
    y = prop,
    fill = anomaly_type
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = ifelse(
        prop > 0.03,
        scales::percent(prop, accuracy = 0.1),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    family = "Times",
    fontface = "bold",
    size = 3.2,
    colour = "black"
  ) +
  scale_y_continuous(
    labels = percent
  ) +
  scale_fill_manual(
    values = cols,
    name = "Checklist Anomaly Type"
  ) +
  labs(
    title = "Checklist IF Typicality Classes",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(base_family = "Times") +
  theme(
    text = element_text(
      family = "Times"
    ),
    plot.title = element_text(
      family = "Times",
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      family = "Times",
      size = 10
    ),
    axis.title = element_text(
      family = "Times",
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      family = "Times",
      size = 10
    ),
    legend.position = "bottom"
  )

(p1 | p2) /
  (p3 | p4) 

windowsFonts(
  Times = windowsFont("Times New Roman")
)
ggsave(
  filename = "file path/avicred flags.jpg",
  plot = (p1 | p2) /
    (p3 | p4) ,
  device = grDevices::jpeg,
  type = "cairo",
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)

##################### MEAN CAS VS AVICRED CORRELATION ##############

plot_data <- ml_complete %>%
  dplyr::filter(
    is.finite(mean_CAS),
    is.finite(AviCRED)
  )

pearson <- cor.test(
  plot_data$mean_CAS,
  plot_data$AviCRED,
  method = "pearson"
)

spearman <- cor.test(
  plot_data$mean_CAS,
  plot_data$AviCRED,
  method = "spearman",
  exact = FALSE
)

# Print results
print(pearson)
print(spearman)

pearson_r <- round(unname(pearson$estimate), 3)
pearson_p <- signif(pearson$p.value, 3)

spearman_rho <- round(unname(spearman$estimate), 3)
spearman_p <- signif(spearman$p.value, 3)

p <- ggplot(
  plot_data,
  aes(x = mean_CAS, y = AviCRED)
) +
  geom_point(alpha = 0.4) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    subtitle = paste0(
      "Pearson r = ", pearson_r,
      " (p = ", pearson_p, ")   |   ",
      "Spearman \u03C1 = ", spearman_rho,
      " (p = ", spearman_p, ")"
    ),
    x = "Mean CAS",
    y = "AviCRED Score"
  ) +
  theme_bw(
    base_size = 14,
    base_family = "Times New Roman"
  ) +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 14
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.subtitle = element_text(size = 14)
  )

print(p)

ggsave(
  filename = "CAS_AviCRED_Correlation.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)

#============= CHECKLISTS SUBMITTED PER AVICRED CLASS ===========#

# ============================================================
# 1. CHECK REQUIRED COLUMNS
# ============================================================

stopifnot(
  all(c("observer_id", "avicred_class") %in% names(avicred_new)),
  all(c("observer_id", "checklist_id") %in% names(checklists_new))
)

# Standardise IDs before joining
avicred_new <- avicred_new %>%
  dplyr::mutate(observer_id = as.character(observer_id))

checklists_new <- checklists_new %>%
  dplyr::mutate(
    observer_id = as.character(observer_id),
    checklist_id = as.character(checklist_id)
  )

# Ensure each observer has only one AviCRED classification
avicred_classes <- avicred_new %>%
  dplyr::filter(!is.na(observer_id), !is.na(avicred_class)) %>%
  dplyr::distinct(observer_id, avicred_class)

# Stop if any observer has conflicting classifications
conflicting_observers <- avicred_classes %>%
  dplyr::count(observer_id, name = "n_classes") %>%
  dplyr::filter(n_classes > 1)

if (nrow(conflicting_observers) > 0) {
  stop("Some observers have multiple AviCRED classes. Resolve these before counting.")
}

avicred_classes <- avicred_classes %>%
  dplyr::distinct(observer_id, .keep_all = TRUE)

# ============================================================
# 2. ASSIGN EACH CHECKLIST TO ITS OBSERVER'S AVICRED CLASS
# ============================================================

checklist_class_data <- checklists_new %>%
  dplyr::filter(
    !is.na(checklist_id),
    !is.na(observer_id)
  ) %>%
  dplyr::distinct(checklist_id, .keep_all = TRUE) %>%
  dplyr::left_join(
    avicred_classes,
    by = "observer_id"
  )

# ============================================================
# 3. COUNT UNIQUE CHECKLISTS PER AVICRED CLASS
# ============================================================

treemap_class <- checklist_class_data %>%
  dplyr::filter(!is.na(avicred_class)) %>%
  dplyr::group_by(avicred_class) %>%
  dplyr::summarise(
    n_checklists = dplyr::n_distinct(checklist_id),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    avicred_class = factor(
      avicred_class,
      levels = c(
        "Very Low",
        "Low",
        "Moderate",
        "High",
        "Very High"
      ),
      ordered = TRUE
    )
  ) %>%
  dplyr::arrange(avicred_class)

# Inspect the counts
print(treemap_class)

ggplot(
  treemap_class,
  aes(
    area = n_checklists,
    fill = avicred_class,
    label = paste0(
      as.character(avicred_class),
      "\n",
      scales::comma(n_checklists),
      " checklists"
    )
  )
) +
  treemapify::geom_treemap(
    colour = "white",
    linewidth = 1
  ) +
  treemapify::geom_treemap_text(
    place = "centre",
    reflow = TRUE,
    colour = "black",
    fontface = "bold",
    family = "Times",
    size = 16,
    min.size = 0
  ) +
  scale_fill_manual(
    values = c(
      "Very Low" = "#E8F6EF",
      "Low" = "#ABEBC6",
      "Moderate" = "#52BE80",
      "High" = "#239B56",
      "Very High" = "#145A32"
    ),
    breaks = c(
      "Very Low",
      "Low",
      "Moderate",
      "High",
      "Very High"
    ),
    drop = FALSE,
    name = "AviCRED Class"
  ) +
  theme_void() +
  theme(
    text = element_text(
      family = "Times",
      colour = "black"
    ),
    plot.title = element_text(
      family = "Times",
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(
      family = "Times",
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      family = "Times",
      size = 14
    ),
    plot.margin = margin(
      t = 15,
      r = 35,
      b = 15,
      l = 35
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

#============= NUMBER OF ANOMALOUS CHECKLISTS PER AVICRED CLASS ===============#

dbWriteTable(  # Store checklist-level anomaly flags in the DuckDB database
  con,
  "chk_flags",
  checklists %>%
    select(
      checklist_id,
      richness_flag,
      abundance_flag
    ),
  overwrite = TRUE
)

dbWriteTable(  # Store observer-level AviCRED classifications in the database
  con,
  "avicred_tbl",
  ml_complete %>%
    select(
      observer_id,
      avicred_class
    ),
  overwrite = TRUE
)

# Count different types of checklist anomalies within each AviCRED class

anomaly_counts <- dbGetQuery(
  con,
  "
  SELECT
      a.avicred_class,
      CASE
          WHEN c.richness_flag = 1 AND c.abundance_flag = 0
              THEN 'Richness'
          WHEN c.richness_flag = 0 AND c.abundance_flag = 1
              THEN 'Abundance'
          WHEN c.richness_flag = 1 AND c.abundance_flag = 1
              THEN 'Both'
          ELSE 'None'
      END AS anomaly_type,
      COUNT(*) AS n_checklists

  FROM chk_flags AS c
  LEFT JOIN checklist_master AS cm
      ON c.checklist_id = cm.checklist_id
  LEFT JOIN avicred_tbl AS a
      ON cm.observer_id = a.observer_id
  WHERE a.avicred_class IS NOT NULL
  GROUP BY
      a.avicred_class,
      anomaly_type
  "
)

anomaly_counts$avicred_class <- factor(  # Arrange AviCRED classes in the desired order
  anomaly_counts$avicred_class,
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  )
)

anomaly_counts$anomaly_type <- factor(  # Set the desired order for anomaly types
  anomaly_counts$anomaly_type,
  levels = c(
    "Richness",
    "Abundance",
    "Both"
  )
)

anomaly_counts <- anomaly_counts %>%  # Remove checklists with no detected anomaly
  filter(anomaly_type != "None")

segment_labels <- anomaly_counts %>%  # Prepare labels for High and Very High AviCRED classes
  filter(avicred_class %in% c("High", "Very High"))

total_labels <- anomaly_counts %>%  # Calculate total anomalies for the remaining AviCRED classes
  filter(!(avicred_class %in% c("High", "Very High"))) %>%
  group_by(avicred_class) %>%
  summarise(
    total = sum(n_checklists),
    .groups = "drop"
  )

ggplot(
  anomaly_counts,
  aes(
    x = avicred_class,
    y = n_checklists,
    fill = anomaly_type
  )
) +
  geom_bar(
    stat = "identity",
    colour = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "Richness" = "#7FC97F",
      "Abundance" = "#FDB462",
      "Both" = "#BEAED4"
    )
  ) +
  scale_y_continuous(
    labels = comma
  ) +
  labs(
    x = "AviCRED Class",
    y = "Number of Flagged Checklists",
    fill = "Anomaly Type"
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Times",
      size = 14
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      colour = "black",
      size = 14
    ),
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.text = element_text(
      size = 14
    ),
    legend.position = "bottom",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  geom_text(
    data = segment_labels,
    aes(label = scales::comma(n_checklists)),
    position = position_stack(vjust = 0.5),
    family = "Times",
    fontface = "bold",
    size = 3.5,
    colour = "black"
  ) +
  geom_text(
    data = total_labels,
    aes(
      x = avicred_class,
      y = total,
      label = scales::comma(total)
    ),
    inherit.aes = FALSE,
    family = "Times",
    fontface = "bold",
    size = 4,
    vjust = -0.5
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

################# OBSERVER YEARS VS CHECKLISTS SUBMITTED #################

year_group <- ml_complete %>%           
  filter(!is.na(avicred_class)) %>%
  mutate(
    avicred_class = factor(
      avicred_class,
      levels = c(
        "Very Low",
        "Low",
        "Moderate",
        "High",
        "Very High"
      )
    )
  ) %>%
  group_by(years_active, avicred_class) %>%
  summarise(
    Total_Checklists = sum(n_checklists, na.rm = TRUE),
    .groups = "drop"
  )

year_totals <- year_group %>%
  group_by(years_active) %>%   # CALCULATE TOTAL CHECKLISTS SUBMITTED BY OBSERVRS IN EACH YEAR GROUP
  summarise(
    Total = sum(Total_Checklists),
    .groups = "drop"
  )

ggplot(
  year_group,
  aes(
    x = factor(years_active),
    y = Total_Checklists,
    fill = avicred_class
  )
) +
  geom_col(
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(  
    data = year_totals,
    aes(
      x = factor(years_active),
      y = Total,
      label = comma(Total)
    ),
    inherit.aes = FALSE,
    vjust = -0.3,
    size = 3,
    family = "Times",
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "Very Low"  = "#FDD0A2",
      "Low"       = "#FC9272",
      "Moderate"  = "#FB6A4A",
      "High"      = "#807DBA",
      "Very High" = "#41B6C4"
    ),
    breaks = c(
      "Very Low",
      "Low",
      "Moderate",
      "High",
      "Very High"
    )
  ) +
  
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.08))
  ) +
  
  labs(
    x = "Years Active",
    y = "Total Checklists Submitted",
    fill = "AviCRED Class"
  ) +
  
  theme_bw(base_family = "Times") +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

###################### CAS VS AVICRED CLASSES HEATMAP #####################

# Extract the checklist–observer relationship from the EBD

checklist_observer <- dbGetQuery(
  con,
  '
  SELECT DISTINCT
    "SAMPLING EVENT IDENTIFIER" AS checklist_id,
    "OBSERVER ID" AS observer_id
  FROM ebd_filtered
  '
)

cas_avicred <- checklists %>%
  dplyr::select(
    checklist_id,
    CAS_class
  ) %>%
  dplyr::left_join(
    checklist_observer,
    by = "checklist_id"
  ) %>%
  dplyr::left_join(
    ml_complete %>%
      dplyr::select(
        observer_id,
        avicred_class
      ),
    by = "observer_id"
  ) %>%
  dplyr::filter(
    !is.na(CAS_class),
    !is.na(avicred_class)
  )

cas_avicred_matrix <- cas_avicred %>%  # Create a frequency matrix of CAS and AviCRED classes
  count(
    CAS_class,
    avicred_class
  ) %>%
  pivot_wider(
    names_from = avicred_class,
    values_from = n,
    values_fill = 0
  )

cas_avicred_matrix

cas_avicred_percent <- cas_avicred %>%  # Calculate percentages within each CAS class
  count(
    CAS_class,
    avicred_class
  ) %>%
  group_by(CAS_class) %>%
  mutate(
    Percentage = round(100 * n / sum(n), 2)  # Percentage within each CAS class
  ) %>%
  pivot_wider(
    names_from = avicred_class,
    values_from = Percentage,
    values_fill = 0
  )

cas_avicred_percent

heatmap_df <- cas_avicred %>%  # Prepare data for the CAS–AviCRED heatmap
  count(CAS_class, avicred_class)

levels_order <- c(  # Set the desired order of credibility classes
  "Very Low",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

heatmap_df$CAS_class <- factor( # Apply the same order to both axes
  heatmap_df$CAS_class,
  levels = levels_order
)

heatmap_df$avicred_class <- factor(
  heatmap_df$avicred_class,
  levels = levels_order
)

ggplot(
  heatmap_df,
  aes(
    x = avicred_class,
    y = CAS_class,
    fill = n
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(label = scales::comma(n)),
    family = "Times",
    size = 3.5
  ) +
  
  scale_fill_viridis_c(
    option = "C",
    trans = "log10",
    labels = scales::comma,
    name = "No. of checklists"
  ) +
  
  labs(
    x = "AviCRED Class",
    y = "CAS Class"
  ) +
  
  theme_minimal() +
  
  theme(
    text = element_text(
      family = "Times",
      size = 14
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  )

###################### DISTRIBUTION OF ORI, SBI and CAS ACROSS AVICRED CLASSES #################

library(ggridges)

# 1. PREPARE DATA
plot_df <- ml_complete %>%
  dplyr::filter(
    !is.na(ORI_bounded),
    !is.na(mean_CAS),
    !is.na(SBI),
    !is.na(avicred_class)
  ) %>%
  dplyr::mutate(
    avicred_class = factor(
      avicred_class,
      levels = c(
        "Very Low", "Low", "Moderate", "High", "Very High"
      ),
      ordered = TRUE
    )
  )

# 2. CALCULATE COMPONENT MEDIANS
medians_ORI <- plot_df %>%
  dplyr::group_by(avicred_class) %>%
  dplyr::summarise(
    median_value = median(ORI_bounded, na.rm = TRUE),
    .groups = "drop"
  )

medians_CAS <- plot_df %>%
  dplyr::group_by(avicred_class) %>%
  dplyr::summarise(
    median_value = median(mean_CAS, na.rm = TRUE),
    .groups = "drop"
  )

medians_SBI <- plot_df %>%
  dplyr::group_by(avicred_class) %>%
  dplyr::summarise(
    median_value = median(SBI, na.rm = TRUE),
    .groups = "drop"
  )

# 3. AVICRED CLASS COLOURS
avicred_colours <- c(
  "Very Low"  = "#FDD0A2",
  "Low"       = "#FC9272",
  "Moderate"  = "#FB6A4A",
  "High"      = "#807DBA",
  "Very High" = "#41B6C4"
)

# 4. COMMON THEME
ridge_theme <- theme_ridges() +
  theme(
    text = element_text(
      family = "Times",
      size = 12,
      colour = "black"
    ),
    plot.title = element_text(
      family = "Times",
      size = 14,
      face = "bold",
      hjust = 0.5,
      colour = "black"
    ),
    axis.title = element_text(
      family = "Times",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text = element_text(
      family = "Times",
      size = 14,
      colour = "black"
    ),
    legend.position = "none",
    plot.margin = margin(8, 14, 8, 14)
  )

# 5. ORI_BOUNDED RIDGE PLOT
ORI_ridge <- ggplot(
  plot_df,
  aes(
    x = ORI_bounded,
    y = avicred_class,
    fill = avicred_class
  )
) +
  geom_density_ridges(
    alpha = 0.85,
    scale = 1.2,
    rel_min_height = 0.001,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    data = medians_ORI,
    aes(xintercept = median_value),
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = medians_ORI,
    aes(
      x = median_value,
      y = avicred_class,
      label = sprintf("%.2f", median_value)
    ),
    inherit.aes = FALSE,
    vjust = -0.8,
    size = 3.5,
    family = "Times",
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_manual(
    values = avicred_colours,
    drop = FALSE
  ) +
  labs(
    title = "Observer Reliability Index",
    x = "Bounded ORI",
    y = "AviCRED Class"
  ) +
  ridge_theme

# 6. MEAN CAS RIDGE PLOT
CAS_ridge <- ggplot(
  plot_df,
  aes(
    x = mean_CAS,
    y = avicred_class,
    fill = avicred_class
  )
) +
  geom_density_ridges(
    alpha = 0.85,
    scale = 1.2,
    rel_min_height = 0.001,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    data = medians_CAS,
    aes(xintercept = median_value),
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = medians_CAS,
    aes(
      x = median_value,
      y = avicred_class,
      label = sprintf("%.2f", median_value)
    ),
    inherit.aes = FALSE,
    vjust = -0.8,
    size = 3.5,
    family = "Times",
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_manual(
    values = avicred_colours,
    drop = FALSE
  ) +
  labs(
    title = "Checklist Atypicality Score",
    x = "Mean CAS",
    y = "AviCRED Class"
  ) +
  ridge_theme

# 7. SBI RIDGE PLOT
SBI_ridge <- ggplot(
  plot_df,
  aes(
    x = SBI,
    y = avicred_class,
    fill = avicred_class
  )
) +
  geom_density_ridges(
    alpha = 0.85,
    scale = 1.2,
    rel_min_height = 0.001,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    data = medians_SBI,
    aes(xintercept = median_value),
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = medians_SBI,
    aes(
      x = median_value,
      y = avicred_class,
      label = sprintf("%.2f", median_value)
    ),
    inherit.aes = FALSE,
    vjust = -0.8,
    size = 3.5,
    family = "Times",
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_manual(
    values = avicred_colours,
    drop = FALSE
  ) +
  labs(
    title = "Sampling Behaviour Index",
    x = "SBI",
    y = "AviCRED Class"
  ) +
  ridge_theme

# 8. COMBINE THE THREE PANELS
component_ridges <- (
  ORI_ridge | CAS_ridge | SBI_ridge
)

print(component_ridges)

ggsave(
  filename = "AviCRED_Component_Ridge_Plots.png",
  plot = component_ridges,
  width = 16,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

#==================CHECKLIST HISTORY COMPARISON: SUFFICIENT VS INSUFFICIENT =================#

# Prepare data
plot_data <- df %>%
  filter(!is.na(checklist_history)) %>%
  select(checklist_history, SBI, mean_CAS, ORI_bounded, AviCRED) %>%
  pivot_longer(
    cols = c(SBI, mean_CAS, ORI_bounded, AviCRED),
    names_to = "Component",
    values_to = "Score"
  ) %>%
  filter(!is.na(Score)) %>%
  mutate(
    checklist_history = factor(
      checklist_history,
      levels = c("Insufficient", "Sufficient")
    ),
    Component = factor(
      Component,
      levels = c("SBI", "mean_CAS", "ORI_bounded", "AviCRED"),
      labels = c("SBI", "Mean CAS", "ORI (bounded)", "AviCRED")
    )
  )

# Consistent colours
component_colours <- c(
  "SBI" = "#2878B5",
  "Mean CAS" = "#E68632",
  "ORI (bounded)" = "#39875B",
  "AviCRED" = "#9B59B6"
)

# Insufficient history panel
p_insufficient <- plot_data %>%
  filter(checklist_history == "Insufficient") %>%
  ggplot(aes(x = Component, y = Score, fill = Component)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = 16,
    outlier.size = 0.8,
    outlier.alpha = 0.35,
    outlier.colour = "black",
    colour = "black",
    linewidth = 0.4
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.3f", after_stat(y))),
    colour = "black",
    size = 3.5,
    fontface = "bold",
    vjust = -0.8
  )+ 
  scale_fill_manual(values = component_colours) +
  labs(
    title = "Insufficient Checklist history (<50 checklists)",
    x = NULL,
    y = "Score"
  ) +
  theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "none"
  )

# Sufficient history panel
p_sufficient <- plot_data %>%
  filter(checklist_history == "Sufficient") %>%
  ggplot(aes(x = Component, y = Score, fill = Component)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = 16,
    outlier.size = 0.8,
    outlier.alpha = 0.35,
    outlier.colour = "black",
    colour = "black",
    linewidth = 0.4
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.3f", after_stat(y))),
    colour = "black",
    size = 3.5,
    fontface = "bold",
    vjust = -0.8
  ) +
  scale_fill_manual(values = component_colours) +
  labs(
    title = "Sufficient Checklist history (≥50 checklists)",
    x = NULL,
    y = "Score"
  ) +
  theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "none"
  )

# Combine panels
combined_plot <- p_insufficient + p_sufficient 

print(combined_plot)

# Export to D: drive
ggsave(
  "file path/AVICRED CHECKLIST HISTORY COMPARISON.jpg",
  plot = combined_plot,
  width = 12,
  height = 6,
  units = "in",
  dpi = 400,
  bg = "white"
)

output_file <- "file path/Checklist_history_group_comparison_results.csv"

results_list <- list()

# ============================================================
# WILCOXON TESTS: INSUFFICIENT VS SUFFICIENT HISTORY
# ============================================================

components <- c("SBI", "mean_CAS", "ORI_bounded", "AviCRED")

for (v in components) {
  
  test_data <- df %>%
    filter(
      !is.na(checklist_history),
      !is.na(.data[[v]])
    )
  
  test_data$checklist_history <- factor(
    test_data$checklist_history,
    levels = c("Insufficient", "Sufficient")
  )
  
  test_data <- test_data %>%
    filter(!is.na(checklist_history))
  
  test_result <- wilcox.test(
    test_data[[v]] ~ test_data$checklist_history,
    exact = FALSE,
    conf.int = TRUE
  )
  
  group_summary <- test_data %>%
    group_by(checklist_history) %>%
    summarise(
      n = n(),
      Median = median(.data[[v]], na.rm = TRUE),
      Q1 = quantile(.data[[v]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[v]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  insufficient <- group_summary %>%
    filter(checklist_history == "Insufficient")
  
  sufficient <- group_summary %>%
    filter(checklist_history == "Sufficient")
  
  results_list[[paste0("Wilcoxon_", v)]] <- data.frame(
    Analysis = "Wilcoxon rank-sum test",
    Variable = v,
    Group_1 = "Insufficient",
    N_1 = insufficient$n,
    Median_1 = insufficient$Median,
    Q1_1 = insufficient$Q1,
    Q3_1 = insufficient$Q3,
    Group_2 = "Sufficient",
    N_2 = sufficient$n,
    Median_2 = sufficient$Median,
    Q1_2 = sufficient$Q1,
    Q3_2 = sufficient$Q3,
    Statistic = unname(test_result$statistic),
    Estimate = if (!is.null(test_result$estimate)) {
      unname(test_result$estimate[1])
    } else {
      NA_real_
    },
    CI_Lower = if (!is.null(test_result$conf.int)) {
      test_result$conf.int[1]
    } else {
      NA_real_
    },
    CI_Upper = if (!is.null(test_result$conf.int)) {
      test_result$conf.int[2]
    } else {
      NA_real_
    },
    P_value = test_result$p.value,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# CHI-SQUARE TEST: AVICRED CLASS BY HISTORY GROUP
# ============================================================

class_data <- df %>%
  filter(
    !is.na(checklist_history),
    !is.na(avicred_class)
  )

class_table <- table(
  class_data$checklist_history,
  class_data$avicred_class
)

class_test <- chisq.test(class_table)

class_percentages <- prop.table(class_table, margin = 1) * 100

for (i in seq_len(nrow(class_table))) {
  for (j in seq_len(ncol(class_table))) {
    
    results_list[[paste0("AviCRED_class_", i, "_", j)]] <-
      data.frame(
        Analysis = "Chi-square: AviCRED class by history",
        Variable = "avicred_class",
        Group_1 = rownames(class_table)[i],
        N_1 = sum(class_table[i, ]),
        Median_1 = NA_real_,
        Q1_1 = NA_real_,
        Q3_1 = NA_real_,
        Group_2 = colnames(class_table)[j],
        N_2 = as.numeric(class_table[i, j]),
        Median_2 = NA_real_,
        Q1_2 = NA_real_,
        Q3_2 = NA_real_,
        Statistic = unname(class_test$statistic),
        Estimate = as.numeric(class_percentages[i, j]),
        CI_Lower = NA_real_,
        CI_Upper = NA_real_,
        P_value = class_test$p.value,
        stringsAsFactors = FALSE
      )
  }
}

# ============================================================
# CHI-SQUARE TEST: ISOLATION FOREST CLASSIFICATION
# ============================================================

if_data <- df %>%
  filter(
    !is.na(checklist_history),
    !is.na(ml_10pct)
  )

if_table <- table(
  if_data$checklist_history,
  factor(
    if_data$ml_10pct,
    levels = c(0, 1),
    labels = c("Remaining 90%", "Top 10%")
  )
)

if_test <- chisq.test(if_table)

if_percentages <- prop.table(if_table, margin = 1) * 100

for (i in seq_len(nrow(if_table))) {
  for (j in seq_len(ncol(if_table))) {
    
    results_list[[paste0("IF_class_", i, "_", j)]] <-
      data.frame(
        Analysis = "Chi-square: IF classification by history",
        Variable = "ml_10pct",
        Group_1 = rownames(if_table)[i],
        N_1 = sum(if_table[i, ]),
        Median_1 = NA_real_,
        Q1_1 = NA_real_,
        Q3_1 = NA_real_,
        Group_2 = colnames(if_table)[j],
        N_2 = as.numeric(if_table[i, j]),
        Median_2 = NA_real_,
        Q1_2 = NA_real_,
        Q3_2 = NA_real_,
        Statistic = unname(if_test$statistic),
        Estimate = as.numeric(if_percentages[i, j]),
        CI_Lower = NA_real_,
        CI_Upper = NA_real_,
        P_value = if_test$p.value,
        stringsAsFactors = FALSE
      )
  }
}

combined_results <- bind_rows(results_list)

write.csv(
  combined_results,
  output_file,
  row.names = FALSE,
  na = ""
)

###################### VARIANCE OF ORI VS CHECKLISTS SUBMITTED AND NUMBER OF OBSERVERS ###############

plot_df <- ml_complete %>%
  dplyr::filter(
    !is.na(ORI_bounded),
    !is.na(n_checklists),
    is.finite(ORI_bounded),
    is.finite(n_checklists),
    n_checklists >= 0
  ) %>%
  dplyr::mutate(
    log_checklists = log10(n_checklists + 1)
  )

# Verify the data used in the plot
cat("Observers included:", nrow(plot_df), "\n")
cat("Unique ORI_bounded values:", dplyr::n_distinct(plot_df$ORI_bounded), "\n")
cat("Checklist count range:",
    min(plot_df$n_checklists), "to",
    max(plot_df$n_checklists), "\n")

ORI_effort_plot <- ggplot(
  plot_df,
  aes(
    x = ORI_bounded,
    y = log_checklists
  )
) +
  geom_hex(
    bins = 60
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    name = "Observer count",
    labels = scales::label_comma()
  ) +
  geom_smooth(
    method = "loess",
    span = 1,
    colour = "red",
    linewidth = 1.0,
    se = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.02)
  ) +
  labs(
    x = "Bounded Observer Reliability Index (ORI)",
    y = expression(Log[10] * "(Total checklists + 1)")
  ) +
  theme_classic() +
  theme(
    text = element_text(
      family = "Times",
      size = 14,
      colour = "black"
    ),
    axis.title = element_text(
      family = "Times",
      size = 14
    ),
    axis.text = element_text(
      family = "Times",
      size = 14,
      colour = "black"
    ),
    legend.title = element_text(
      family = "Times",
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      family = "Times",
      size = 14
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    plot.margin = margin(12, 15, 12, 15)
  )

print(ORI_effort_plot)

ggsave(
  filename = "ORI_Checklist_Effort_Hexbin.png",
  plot = ORI_effort_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

#============= IF ANOMALY SCORE VS AVICRED AGREEMENT ==============#

############# CONFUSION MATRIX ##############

observer_cm <- ml_complete %>% # Prepare observer-level data for comparing AviCRED with Isolation Forest
  filter(
    !is.na(avicred_class),   
    !is.na(ml_10pct)         
  ) %>%
  mutate(
    avicred_high = if_else(
      avicred_class %in% c("High", "Very High"),  # Classify High/Very High as positive
      1,
      0
    )
  )

cm_obs <- confusionMatrix(  # Generate confusion matrix for AviCRED vs Isolation Forest
  data = factor(observer_cm$avicred_high),
  reference = factor(observer_cm$ml_10pct),
  positive = "1"
)

cm_obs$table       # Confusion matrix counts
cm_obs$overall     # Overall performance metrics
cm_obs$byClass     # Class-specific performance metrics

checklist_cm <- checklists %>%  # Prepare checklist-level data for comparing CAS with Isolation Forest
  filter(
    !is.na(CAS_class),       
    !is.na(if_top10_typicality)         
  ) %>%
  mutate(
    cas_high = if_else(
      CAS_class %in% c("High", "Very High"),  # Classify High/Very High as positive
      1,
      0
    )
  )

cm_chk <- confusionMatrix(  # Generate confusion matrix for CAS vs Isolation Forest
  data = factor(checklist_cm$cas_high),
  reference = factor(checklist_cm$if_top10_typicality),
  positive = "1"
)

cm_chk$table       # Confusion matrix counts
cm_chk$overall     # Overall performance metrics
cm_chk$byClass     # Class-specific performance metrics

############# CHECKLIST MATRIX ##############

chk_df_plot <- as.data.frame(cm_chk$table)  # Convert the checklist confusion matrix to a data frame

chk_df_plot$Prediction <- factor(  # Rename prediction categories for plotting
  chk_df_plot$Prediction,
  levels = c("0", "1"),
  labels = c("Remaining", "High")
)

chk_df_plot$Reference <- factor(  # Rename reference categories for plotting
  chk_df_plot$Reference,
  levels = c("0", "1"),
  labels = c("Remaining", "Top 10%")
)

p1 <- ggplot(
  chk_df_plot,
  aes(
    x = Reference,
    y = Prediction,
    fill = Freq
  )
) +
  geom_tile(
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = scales::comma(Freq)),
    family = "Times",
    fontface = "bold",
    size = 4.5
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#F9C80E"  
  ) +
  labs(
    title = "Checklist-level",
    x = "IF Typicality Top 10%",
    y = "CAS High/Very High",
    fill = "Count"
  ) +
  theme_bw(base_family = "Times") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 14
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      size = 14
    ),
    axis.text.y = element_text(
      size = 14,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    ),
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.text = element_text(
      size = 14
    ),
    panel.grid = element_blank()
  )

################ OBSERVER MATRIX ####################

obs_df <- as.data.frame(cm_obs$table) # Convert the observer confusion matrix to a data frame

obs_df$Prediction <- factor(  # Rename prediction categories for plotting
  obs_df$Prediction,
  levels = c("0", "1"),
  labels = c("Remaining", "High")
)

obs_df$Reference <- factor(  # Rename reference categories for plotting
  obs_df$Reference,
  levels = c("0", "1"),
  labels = c("Remaining", "Top 10%")
)

max_count <- max(  # Find the largest cell count across both confusion matrices
  c(obs_df$Freq, chk_df_plot$Freq)
)

p2 <- ggplot(
  obs_df,
  aes(
    x = Reference,
    y = Prediction,
    fill = Freq
  )
) +
  geom_tile(
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = scales::comma(Freq)),
    family = "Times",
    fontface = "bold",
    size = 4.5
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#7678ED"
  ) +
  labs(
    title = "Observer-level",
    x = "IF Top 10%",
    y = "AviCRED High/Very High",
    fill = "Count"
  ) +
  theme_bw(base_family = "Times") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 14
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      size = 14
    ),
    axis.text.y = element_text(
      size = 14,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    ),
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.text = element_text(
      size = 14
    ),
    panel.grid = element_blank()
  )

combined_plot <- (p1 + p2) +
  plot_layout(guides = "collect") 

combined_plot

# ============================================================
# MCC: INDEPENDENT IF VALIDATION OF AVICRED AND CAS
# ============================================================

library(yardstick)

# ============================================================
# 1. OBSERVER-LEVEL: AVICRED VS IF
# ============================================================

observer_cm <- ml_complete %>%
  dplyr::filter(
    !is.na(avicred_class),
    !is.na(ml_10pct)
  ) %>%
  dplyr::mutate(
    avicred_high = dplyr::if_else(
      avicred_class %in% c("High", "Very High"),
      1L, 0L
    ),
    if_top10 = dplyr::if_else(
      ml_10pct == 1,
      1L, 0L
    )
  )

# Ensure one record per observer
if (anyDuplicated(observer_cm$observer_id) > 0) {
  stop("Duplicate observer IDs found in observer_cm.")
}

# MCC
observer_mcc <- observer_cm %>%
  dplyr::mutate(
    avicred_high = factor(avicred_high, levels = c(0, 1)),
    if_top10 = factor(if_top10, levels = c(0, 1))
  ) %>%
  yardstick::mcc(
    truth = avicred_high,
    estimate = if_top10
  )

# Counts and overlap
observer_results <- observer_cm %>%
  dplyr::summarise(
    High_Credibility = sum(avicred_high == 1),
    IF_Top10 = sum(if_top10 == 1),
    Overlap = sum(avicred_high == 1 & if_top10 == 1)
  ) %>%
  dplyr::mutate(
    `Overlap vs AviCRED (%)` = dplyr::if_else(
      High_Credibility > 0,
      100 * Overlap / High_Credibility,
      NA_real_
    ),
    `Overlap vs IF (%)` = dplyr::if_else(
      IF_Top10 > 0,
      100 * Overlap / IF_Top10,
      NA_real_
    ),
    MCC = observer_mcc$.estimate
  )

print(observer_results)

# ============================================================
# 2. CHECKLIST-LEVEL: CAS VS IF
# ============================================================

checklist_cm <- checklists %>%
  dplyr::filter(
    !is.na(CAS_class),
    !is.na(if_top10_typicality)
  ) %>%
  dplyr::mutate(
    cas_high = dplyr::if_else(
      CAS_class %in% c("High", "Very High"),
      1L, 0L
    ),
    if_top10 = dplyr::if_else(
      if_top10_typicality == 1,
      1L, 0L
    )
  )

# Ensure one record per checklist
if (anyDuplicated(checklist_cm$checklist_id) > 0) {
  stop(
    "Duplicate checklist IDs found. Resolve duplicates before calculating MCC."
  )
}

# MCC
checklist_mcc <- checklist_cm %>%
  dplyr::mutate(
    cas_high = factor(cas_high, levels = c(0, 1)),
    if_top10 = factor(if_top10, levels = c(0, 1))
  ) %>%
  yardstick::mcc(
    truth = cas_high,
    estimate = if_top10
  )

# Counts and overlap
checklist_results <- checklist_cm %>%
  dplyr::summarise(
    High_CAS = sum(cas_high == 1),
    IF_Top10 = sum(if_top10 == 1),
    Overlap = sum(cas_high == 1 & if_top10 == 1)
  ) %>%
  dplyr::mutate(
    `Overlap vs CAS (%)` = dplyr::if_else(
      High_CAS > 0,
      100 * Overlap / High_CAS,
      NA_real_
    ),
    `Overlap vs IF (%)` = dplyr::if_else(
      IF_Top10 > 0,
      100 * Overlap / IF_Top10,
      NA_real_
    ),
    MCC = checklist_mcc$.estimate
  )

print(checklist_results)

# ============================================================
# 3. COMBINE RESULTS INTO A SUMMARY TABLE
# ============================================================

summary_table <- tibble::tibble(
  Comparison = c(
    "Checklist (CAS vs IF)",
    "Observer (AviCRED vs IF)"
  ),
  
  `High Credibility Total` = c(
    checklist_results$High_CAS,
    observer_results$High_Credibility
  ),
  
  `IF Top 10% Total` = c(
    checklist_results$IF_Top10,
    observer_results$IF_Top10
  ),
  
  Overlap = c(
    checklist_results$Overlap,
    observer_results$Overlap
  ),
  
  `Overlap vs High (%)` = round(
    c(
      checklist_results$`Overlap vs CAS (%)`,
      observer_results$`Overlap vs AviCRED (%)`
    ),
    2
  ),
  
  `Overlap vs IF (%)` = round(
    c(
      checklist_results$`Overlap vs IF (%)`,
      observer_results$`Overlap vs IF (%)`
    ),
    2
  ),
  
  MCC = round(
    c(
      checklist_results$MCC,
      observer_results$MCC
    ),
    3
  )
)

print(summary_table)

# Export
write.csv(
  summary_table,
  "MCC Summary.csv",
  row.names = FALSE
)

#============ AGREEMENT BETWEEN ANOMALY FLAGGED OBSERVERS ACROSS IF AND AVICRED =========#

# ============================================================
# 1. RETRIEVE CHECKLIST-OBSERVER MAPPING
# ============================================================

checklist_obs <- DBI::dbGetQuery(
  con,
  "
  SELECT
      checklist_id,
      observer_id
  FROM checklist_master
  "
)

# Confirm that checklist IDs map uniquely to observers
duplicate_checklists <- checklist_obs %>%
  dplyr::count(checklist_id) %>%
  dplyr::filter(n > 1)

if (nrow(duplicate_checklists) > 0) {
  stop("Duplicate checklist IDs found in checklist_master.")
}

# ============================================================
# 2. JOIN OBSERVER IDs TO CHECKLIST-LEVEL ANOMALY DATA
# ============================================================

checklists_obs <- checklists %>%
  dplyr::select(-dplyr::any_of("observer_id")) %>%
  dplyr::left_join(
    checklist_obs,
    by = "checklist_id"
  )

stopifnot("observer_id" %in% names(checklists_obs))

cat("Checklist rows:", nrow(checklists_obs), "\n")
cat(
  "Rows without observer ID:",
  sum(is.na(checklists_obs$observer_id)),
  "\n"
)
cat(
  "Unique observers:",
  dplyr::n_distinct(checklists_obs$observer_id, na.rm = TRUE),
  "\n"
)

# ============================================================
# 3. CREATE OBSERVER-LEVEL ANOMALY FLAG
# ============================================================

observer_anomaly <- checklists_obs %>%
  dplyr::filter(!is.na(observer_id)) %>%
  dplyr::group_by(observer_id) %>%
  dplyr::summarise(
    richness_anomaly = as.integer(
      any(richness_flag == 1, na.rm = TRUE)
    ),
    abundance_anomaly = as.integer(
      any(abundance_flag == 1, na.rm = TRUE)
    ),
    anomaly_flag = as.integer(
      any(
        richness_flag == 1 | abundance_flag == 1,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

# Inspect observer-level anomaly counts
print(
  observer_anomaly %>%
    dplyr::count(anomaly_flag) %>%
    dplyr::mutate(
      percentage = round(100 * n / sum(n), 2)
    )
)

# ============================================================
# 4. COMBINE ANOMALY, AVICRED AND IF INFORMATION
# ============================================================

observer_cm <- ml_complete %>%
  dplyr::left_join(
    observer_anomaly,
    by = "observer_id"
  ) %>%
  dplyr::mutate(
    # Missing anomaly flags are set to zero only for observers
    # with no matched anomaly record; inspect unmatched IDs below.
    anomaly_flag = tidyr::replace_na(anomaly_flag, 0L),
    
    avicred_high = as.integer(
      avicred_class %in% c("High", "Very High")
    ),
    
    # Use the 10% threshold as the single IF comparison group.
    if_top10 = as.integer(ml_10pct == 1)
  )

# ============================================================
# 5. VALIDATE THE MERGED OBSERVER DATA
# ============================================================

cat("\nRows in ml_complete:", nrow(ml_complete), "\n")
cat("Rows in observer_cm:", nrow(observer_cm), "\n")
cat(
  "Unique observer IDs in observer_cm:",
  dplyr::n_distinct(observer_cm$observer_id),
  "\n"
)

cat("\nAnomaly flag counts:\n")
print(table(observer_cm$anomaly_flag, useNA = "ifany"))

cat("\nAviCRED high-class counts:\n")
print(table(observer_cm$avicred_high, useNA = "ifany"))

cat("\nIF top-10% counts:\n")
print(table(observer_cm$if_top10, useNA = "ifany"))

# Check for unmatched observers
unmatched_observers <- ml_complete %>%
  dplyr::anti_join(
    observer_anomaly,
    by = "observer_id"
  )

########### MATTHEW'S COEFFICIENT AND PERCENTAGE OVERLAP #############

calc_validation_metrics <- function(prediction, reference, model_name) {
  
  prediction <- factor(prediction, levels = c(0, 1))
  reference  <- factor(reference,  levels = c(0, 1))
  
  cm <- table(
    Prediction = prediction,
    Reference = reference
  )
  
  TN <- as.numeric(cm["0","0"])
  FP <- as.numeric(cm["1","0"])
  FN <- as.numeric(cm["0","1"])
  TP <- as.numeric(cm["1","1"])
  
  # Matthews Correlation Coefficient (MCC)
  denominator <- sqrt(
    as.numeric(TP + FP) *
      as.numeric(TP + FN) *
      as.numeric(TN + FP) *
      as.numeric(TN + FN)
  )
  
  MCC <- if (denominator == 0) {
    NA_real_
  } else {
    (TP * TN - FP * FN) / denominator
  }
  
  ## Overlap statistics
  prediction_overlap <- TP / (TP + FP) * 100
  anomaly_detection  <- TP / (TP + FN) * 100
  
  data.frame(
    Model = model_name,
    True_Negative = TN,
    False_Positive = FP,
    False_Negative = FN,
    True_Positive = TP,
    MCC = MCC,
    MCC_percent = MCC * 100,
    Prediction_Overlap_percent = prediction_overlap,
    Observer_Anomaly_Detection_percent = anomaly_detection
  )
}

validation_results <- rbind(
  
  calc_validation_metrics(
    observer_cm$avicred_high,
    observer_cm$anomaly_flag,
    "AviCRED vs Observer Anomaly"
  ),
  
  calc_validation_metrics(
    observer_cm$if_top10,
    observer_cm$anomaly_flag,
    "Isolation Forest vs Observer Anomaly"
  )
)

validation_results

write.csv(
  validation_results,
  "MCC_summary_flagged.csv",
  row.names = FALSE
)

########## CONFUSION MATRIX ###########

cm_avicred <- confusionMatrix(   # Compare AviCRED classification with observer-level anomalies
  data = factor(observer_cm$avicred_high, levels = c(0, 1)),
  reference = factor(observer_cm$anomaly_flag, levels = c(0, 1)),
  positive = "1"
)

cm_if <- confusionMatrix(  # Compare Isolation Forest classification with observer-level anomalies
  data = factor(observer_cm$if_top10, levels = c(0, 1)),
  reference = factor(observer_cm$anomaly_flag, levels = c(0, 1)),
  positive = "1"
)

cm1_df <- as.data.frame(cm_avicred$table) %>%  # Prepare AviCRED confusion matrix for plotting
  mutate(
    Prediction = factor(
      Prediction,
      levels = c("0", "1"),
      labels = c("Remaining", "AviCRED/IF Top 10%\nMetric")
    ),
    Reference = factor(
      Reference,
      levels = c("0", "1"),
      labels = c("No Observer\nAnomaly", "Observer\nAnomaly")
    ),
    Panel = "AviCRED vs Observer Anomaly"
  )

cm2_df <- as.data.frame(cm_if$table) %>% # Prepare Isolation Forest confusion matrix for plotting
  mutate(
    Prediction = factor(
      Prediction,
      levels = c("0", "1"),
      labels = c("Remaining", "AviCRED/IF Top 10%\nMetric")
    ),
    Reference = factor(
      Reference,
      levels = c("0", "1"),
      labels = c("No Observer\nAnomaly", "Observer\nAnomaly")
    ),
    Panel = "Isolation Forest vs Observer Anomaly"
  )

confusion_df <- bind_rows(cm1_df, cm2_df)  # Combine both confusion matrices into a single data frame

ggplot(
  confusion_df,
  aes(
    x = Reference,
    y = Prediction,
    fill = Freq
  )
) +
  geom_tile(
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = scales::comma(Freq)),
    family = "Times",
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(
    ~Panel,
    nrow = 1
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#93032E"
  ) +
  labs(
    x = "Observer Anomaly Status",
    y = "Predicted Classification",
    fill = "Count"
  ) +
  theme_bw(base_family = "Times") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(
      size = 14
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    legend.title = element_text(
      face = "bold"
    ),
    panel.grid = element_blank()
  )

#============= ANOMALOUS CHECKLISTS PER OBSERVER CLASSIFICATION ===========#

# ============================================================
# 1. CREATE CHECKLIST-LEVEL ANOMALY FLAGS
# ============================================================

checklist_anomalies <- checklists %>%
  dplyr::transmute(
    checklist_id,
    anomaly_flag = as.integer(
      richness_flag == 1 | abundance_flag == 1
    )
  ) %>%
  dplyr::left_join(
    checklist_obs,
    by = "checklist_id"
  )

# ============================================================
# 2. LINK CHECKLISTS TO OBSERVER-LEVEL AVICRED METRICS
# ============================================================

observer_checklists <- ml_complete %>%
  dplyr::select(
    observer_id,
    avicred_class,
    ml_10pct
  ) %>%
  dplyr::left_join(
    checklist_anomalies,
    by = "observer_id"
  )

# ============================================================
# 3. SUMMARISE CHECKLIST ANOMALIES BY AVICRED CLASS
# ============================================================

anomalies_by_avicred <- observer_checklists %>%
  dplyr::group_by(avicred_class) %>%
  dplyr::summarise(
    Observers = dplyr::n_distinct(observer_id),
    Total_Checklists = sum(!is.na(checklist_id)),
    Anomalous_Checklists = sum(
      anomaly_flag == 1,
      na.rm = TRUE
    ),
    Mean_Anomalous_Checklists_per_Observer =
      Anomalous_Checklists / Observers,
    Percent_Anomalous_Checklists = dplyr::if_else(
      Total_Checklists > 0,
      100 * Anomalous_Checklists / Total_Checklists,
      NA_real_
    ),
    .groups = "drop"
  )

print(anomalies_by_avicred)

# ============================================================
# 4. SUMMARISE CHECKLIST ANOMALIES BY IF TOP-10% GROUP
# ============================================================

anomalies_by_if <- observer_checklists %>%
  dplyr::mutate(
    IF_Group = dplyr::if_else(
      ml_10pct == 1,
      "Top 10%",
      "Remaining"
    )
  ) %>%
  dplyr::group_by(IF_Group) %>%
  dplyr::summarise(
    Observers = dplyr::n_distinct(observer_id),
    Total_Checklists = sum(!is.na(checklist_id)),
    Anomalous_Checklists = sum(
      anomaly_flag == 1,
      na.rm = TRUE
    ),
    Mean_Anomalous_Checklists_per_Observer =
      Anomalous_Checklists / Observers,
    Percent_Anomalous_Checklists = dplyr::if_else(
      Total_Checklists > 0,
      100 * Anomalous_Checklists / Total_Checklists,
      NA_real_
    ),
    .groups = "drop"
  )

print(anomalies_by_if)

################ NUMBER OF OBSERVERS WITH ANOMALOUS RECORDS PER CLASSIFCATION ############

# ============================================================
# 1. COMBINE OBSERVER ANOMALY INFORMATION WITH AVICRED
# ============================================================

observer_summary <- ml_complete %>%
  dplyr::left_join(
    observer_anomaly %>%
      dplyr::select(observer_id, anomaly_flag),
    by = "observer_id"
  ) %>%
  dplyr::mutate(
    anomaly_flag = tidyr::replace_na(anomaly_flag, 0L)
  )

# Check that each observer occurs only once
if (anyDuplicated(observer_summary$observer_id) > 0) {
  stop("Duplicate observer IDs found in observer_summary.")
}

# ============================================================
# 2. SUMMARISE ANOMALOUS OBSERVERS BY AVICRED CLASS
# ============================================================

anomalies_by_avicred <- observer_summary %>%
  dplyr::filter(!is.na(avicred_class)) %>%
  dplyr::group_by(avicred_class) %>%
  dplyr::summarise(
    Total_Observers = dplyr::n_distinct(observer_id),
    Observers_with_Anomalies = sum(anomaly_flag == 1, na.rm = TRUE),
    Percent = 100 * Observers_with_Anomalies / Total_Observers,
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(
      avicred_class,
      levels = c("Very Low", "Low", "Moderate", "High", "Very High")
    )
  )

print(anomalies_by_avicred)

# ============================================================
# 3. SUMMARISE ANOMALOUS OBSERVERS BY IF TOP-10% GROUP
# ============================================================

anomalies_by_if <- observer_summary %>%
  dplyr::mutate(
    IF_Class = dplyr::if_else(
      ml_10pct == 1,
      "Top 10%",
      "Remaining",
      missing = NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(IF_Class)) %>%
  dplyr::group_by(IF_Class) %>%
  dplyr::summarise(
    Total_Observers = dplyr::n_distinct(observer_id),
    Observers_with_Anomalies = sum(anomaly_flag == 1, na.rm = TRUE),
    Percent = 100 * Observers_with_Anomalies / Total_Observers,
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(IF_Class, levels = c("Top 10%", "Remaining"))
  )

print(anomalies_by_if)

############### ANOMALOUS CHECKLISTS PER CHECKLIST CLASSIFICATION ############

# ============================================================
# 1. SUMMARISE ANOMALOUS CHECKLISTS BY CAS CLASS
# ============================================================

checklist_anomalies <- checklists %>%
  dplyr::mutate(
    anomaly_flag = as.integer(
      richness_flag == 1 | abundance_flag == 1
    )
  )

anomalies_by_cas <- checklist_anomalies %>%
  dplyr::filter(!is.na(CAS_class)) %>%
  dplyr::group_by(CAS_class) %>%
  dplyr::summarise(
    Total_Checklists = dplyr::n(),
    Anomalous_Checklists = sum(anomaly_flag, na.rm = TRUE),
    Percentage_Anomalous = dplyr::if_else(
      Total_Checklists > 0,
      100 * Anomalous_Checklists / Total_Checklists,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(
      CAS_class,
      levels = c("Very Low", "Low", "Moderate", "High", "Very High")
    )
  )

print(anomalies_by_cas)

# ============================================================
# 2. SUMMARISE ANOMALOUS CHECKLISTS BY IF TOP-10% GROUP
# ============================================================

anomalies_by_if <- checklist_anomalies %>%
  dplyr::mutate(
    IF_Group = dplyr::case_when(
      if_top10_typicality == 1 ~ "IF Top 10%",
      if_top10_typicality == 0 ~ "Remaining 90%",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(IF_Group)) %>%
  dplyr::group_by(IF_Group) %>%
  dplyr::summarise(
    Total_Checklists = dplyr::n(),
    Anomalous_Checklists = sum(anomaly_flag, na.rm = TRUE),
    Percentage_Anomalous = dplyr::if_else(
      Total_Checklists > 0,
      100 * Anomalous_Checklists / Total_Checklists,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(IF_Group, levels = c("IF Top 10%", "Remaining 90%"))
  )

print(anomalies_by_if)








