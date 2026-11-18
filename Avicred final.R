# ============================================================================================
# PROJECT: AVICRED
# OBJECTIVE: CHECKLIST-DATA QUALITY ASSESSMENT OF eBIRD DATA VIA MULTI-METHOD FRAMEWORK
#
# RUN THE R SCRIPT TO FILTER AND PROCESS EBIRD DATA AND ASSESS CHECKLIST DATA-QUALITY WITH
# BY COMPUTING CAS SCORE. SUBSEQUENTLY, VALIDATE RESULTS WITH
# MACHINE-LEARNING ALGORITHMS TO UNDERSAND CHECKLIST PATTERNS AND TRENDS
# ============================================================================================

############# CALL LIBRARIES ##############   

#=== INSTALL LIBRARIES FIRST FROM install_dependencies.R file in the folder ====#

library(dplyr)
library(readr)
library(stringr)
library(Hmisc)
library(data.table)
library(DBI)
library(duckdb)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(scales)
library(extrafont)
library(isotree)
library(ComplexUpset)
library(sf)
library(treemapify)
library(terra)
library(patchwork)
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
  dbdir = "FILE PATH/avicred.duckdb",
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
TO 'FILE PATH/ebd_filtered.txt'
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

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_filtered AS
SELECT *
FROM read_csv(
  'FILE PATH/ebd_filtered.txt',
  header = TRUE,
  delim = '\t'
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

#=============== CAS SCORE FORMULATION =============#

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

###################### CAS POST-HOC VALIDATION ############

# ============================================================
# INDEPENDENT ISOLATION FOREST BENCHMARK
#
# 1. Fit a separate Isolation Forest for each state_season.
# 2. Calculate a continuous IF anomaly score for every
#    eligible checklist.
# 3. Calculate group-specific 10th and 90th percentile cutoffs.
# 4. Assign bottom 10%, top 10%, and combined binary flags.
# 5. Preserve all existing CAS columns and values.
# 6. Export the complete dataset and group-level diagnostics.
# ============================================================

# ------------------------------------------------------------
# 1. CHECK INPUT DATA
# ------------------------------------------------------------

# Read CSV file
csv_file <- "FILE PATH/CAS_FINAL_WITH_FLAGS.csv"

cas <- read.csv(
  csv_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check data
head(cas)
str(cas)
names(cas)

stopifnot(exists("cas"))

cas <- as.data.table(cas)

required_columns <- c(
  "checklist_id",
  "state_season",
  "protocol",
  "richness",
  "abundance",
  "n_observers",
  "log_duration",
  "log_distance",
  "start_time_sin",
  "start_time_cos",
  "CAS"
)

missing_columns <- setdiff(required_columns, names(cas))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(cas$checklist_id) > 0) {
  stop("Duplicate checklist IDs detected. Resolve these before IF modelling.")
}

cat("\nINPUT DATA\n")
cat("Rows:", nrow(cas), "\n")
cat("Columns:", ncol(cas), "\n")
cat("Unique checklists:", uniqueN(cas$checklist_id), "\n")
cat("State-season groups:", uniqueN(cas$state_season), "\n")

# ------------------------------------------------------------
# 3. DEFINE IF FEATURES
# ------------------------------------------------------------
#
# Ecological features:
#   richness
#   log1p(abundance)
#
# CAS effort/context predictors:
#   protocol
#   log_duration
#   log_distance
#   start_time_sin
#   start_time_cos
#   n_observers
#
# state_season is used for grouping, NOT as a predictor.
#
# CAS itself is NOT used as an IF predictor.
# CAS scores and flags are NOT used to train the IF model.
# ------------------------------------------------------------


# Convert protocol to a factor
cas[, protocol_factor := factor(protocol)]

# Create protocol dummy variables
protocol_matrix <- model.matrix(
  ~ protocol_factor,
  data = cas
)[, -1, drop = FALSE]

# Prepare numerical features
numeric_features <- data.table(
  richness = as.numeric(cas$richness),
  log_abundance = log1p(as.numeric(cas$abundance)),
  n_observers = as.numeric(cas$n_observers),
  log_duration = as.numeric(cas$log_duration),
  log_distance = as.numeric(cas$log_distance),
  start_time_sin = as.numeric(cas$start_time_sin),
  start_time_cos = as.numeric(cas$start_time_cos)
)

# Combine numerical and protocol features
X <- as.matrix(
  cbind(
    numeric_features,
    as.data.frame(protocol_matrix)
  )
)

storage.mode(X) <- "double"

cat("\nIF FEATURE MATRIX\n")
cat("Number of features:", ncol(X), "\n")
print(colnames(X))

# ------------------------------------------------------------
# 4. INITIALISE OUTPUT COLUMNS
# ------------------------------------------------------------

cas[, IF_anomaly_score := NA_real_]

cas[, IF_bottom10_flag := NA_integer_]
cas[, IF_top10_flag := NA_integer_]
cas[, IF_either10_flag := NA_integer_]

cas[, IF_bottom10_cutoff := NA_real_]
cas[, IF_top10_cutoff := NA_real_]

cas[, IF_scoring_status := "Not scored"]

# ------------------------------------------------------------
# 5. IDENTIFY STATE-SEASON GROUPS
# ------------------------------------------------------------

group_indices <- split(
  seq_len(nrow(cas)),
  as.character(cas$state_season)
)

cat("\nNumber of groups:", length(group_indices), "\n")

# ------------------------------------------------------------
# 6. ISOLATION FOREST SETTINGS
# ------------------------------------------------------------

N_TREES <- 300

SAMPLE_SIZE <- 256

N_THREADS <- 2

MIN_GROUP_SIZE <- 100

BASE_SEED <- 20261105

# ------------------------------------------------------------
# 7. CREATE GROUP-LEVEL DIAGNOSTICS TABLE
# ------------------------------------------------------------

IF_summary <- data.table(
  state_season = character(),
  total_checklists = integer(),
  eligible_checklists = integer(),
  scored_checklists = integer(),
  bottom10_cutoff = numeric(),
  top10_cutoff = numeric(),
  bottom10_n = integer(),
  top10_n = integer(),
  bottom10_percent = numeric(),
  top10_percent = numeric(),
  overlapping_flags = integer(),
  scoring_status = character()
)

# ------------------------------------------------------------
# 8. FIT ISOLATION FOREST SEPARATELY WITHIN EACH GROUP
# ------------------------------------------------------------

group_names <- names(group_indices)

for (i in seq_along(group_indices)) {
  
  current_group <- group_names[i]
  
  idx <- group_indices[[i]]
  
  cat(
    "\nProcessing group",
    i, "of", length(group_indices),
    ":", current_group,
    "\n"
  )
  
  # Extract feature matrix for this group
  X_group <- X[idx, , drop = FALSE]
  
  # Identify rows with complete, finite predictor values
  valid_rows <- complete.cases(X_group) &
    apply(X_group, 1, function(x) all(is.finite(x)))
  
  valid_idx <- idx[valid_rows]
  
  X_valid <- X_group[valid_rows, , drop = FALSE]
  
  n_total <- length(idx)
  n_valid <- nrow(X_valid)
  
  # ----------------------------------------------------------
  # Check minimum group size
  # ----------------------------------------------------------
  
  if (n_valid < MIN_GROUP_SIZE) {
    
    cas[idx, IF_scoring_status := "Insufficient group size"]
    
    IF_summary <- rbind(
      IF_summary,
      data.table(
        state_season = current_group,
        total_checklists = n_total,
        eligible_checklists = n_valid,
        scored_checklists = 0L,
        bottom10_cutoff = NA_real_,
        top10_cutoff = NA_real_,
        bottom10_n = NA_integer_,
        top10_n = NA_integer_,
        bottom10_percent = NA_real_,
        top10_percent = NA_real_,
        overlapping_flags = NA_integer_,
        scoring_status = "Insufficient group size"
      )
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # Remove features with no variation within this group
  # ----------------------------------------------------------
  
  variable_features <- vapply(
    seq_len(ncol(X_valid)),
    function(j) {
      length(unique(X_valid[, j])) > 1
    },
    logical(1)
  )
  
  X_valid <- X_valid[, variable_features, drop = FALSE]
  
  if (ncol(X_valid) < 2) {
    
    cas[idx, IF_scoring_status := "Insufficient feature variation"]
    
    IF_summary <- rbind(
      IF_summary,
      data.table(
        state_season = current_group,
        total_checklists = n_total,
        eligible_checklists = n_valid,
        scored_checklists = 0L,
        bottom10_cutoff = NA_real_,
        top10_cutoff = NA_real_,
        bottom10_n = NA_integer_,
        top10_n = NA_integer_,
        bottom10_percent = NA_real_,
        top10_percent = NA_real_,
        overlapping_flags = NA_integer_,
        scoring_status = "Insufficient feature variation"
      )
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # Fit genuine Isolation Forest
  # ----------------------------------------------------------
  
  set.seed(BASE_SEED + i)
  
  model_if <- isotree::isolation.forest(
    X_valid,
    ntrees = N_TREES,
    sample_size = min(SAMPLE_SIZE, nrow(X_valid)),
    ndim = 1,
    nthreads = N_THREADS
  )
  
  
  # ----------------------------------------------------------
  # Predict continuous anomaly scores
  # ----------------------------------------------------------
  
  scores <- as.numeric(
    predict(
      model_if,
      X_valid,
      type = "score"
    )
  )
  
  if (length(scores) != n_valid ||
      any(!is.finite(scores))) {
    
    stop(
      "Invalid IF scores in group: ",
      current_group
    )
  }
  
  # Save continuous scores
  cas[valid_idx, IF_anomaly_score := scores]
  
  cas[valid_idx, IF_scoring_status := "Scored"]
  
  
  # ----------------------------------------------------------
  # Calculate within-group percentile thresholds
  # ----------------------------------------------------------
  
  cutoffs <- quantile(
    scores,
    probs = c(0.10, 0.90),
    na.rm = TRUE,
    type = 7,
    names = FALSE
  )
  
  bottom_cutoff <- cutoffs[1]
  top_cutoff <- cutoffs[2]
  
  # Save group-specific cutoffs against each checklist
  cas[valid_idx, IF_bottom10_cutoff := bottom_cutoff]
  cas[valid_idx, IF_top10_cutoff := top_cutoff]
  
  
  # ----------------------------------------------------------
  # Assign binary flags
  # ----------------------------------------------------------
  
  # If the score distribution has no separable lower and
  # upper tails, retain NA flags rather than assigning
  # every checklist to both extreme categories.
  
  if (top_cutoff <= bottom_cutoff) {
    
    cas[valid_idx, IF_scoring_status :=
          "Scored - no separable percentile tails"]
    
    bottom_flag <- rep(NA_integer_, n_valid)
    top_flag <- rep(NA_integer_, n_valid)
    either_flag <- rep(NA_integer_, n_valid)
    
  } else {
    
    # Bottom 10% = lowest anomaly scores
    bottom_flag <- as.integer(
      scores <= bottom_cutoff
    )
    
    # Top 10% = highest anomaly scores
    top_flag <- as.integer(
      scores >= top_cutoff
    )
    
    # Either extreme
    either_flag <- as.integer(
      bottom_flag == 1 | top_flag == 1
    )
    
    cas[valid_idx, IF_bottom10_flag := bottom_flag]
    cas[valid_idx, IF_top10_flag := top_flag]
    cas[valid_idx, IF_either10_flag := either_flag]
  }
  
  
  # ----------------------------------------------------------
  # Group-level summary
  # ----------------------------------------------------------
  
  bottom_n <- if (all(is.na(bottom_flag))) {
    NA_integer_
  } else {
    sum(bottom_flag == 1, na.rm = TRUE)
  }
  
  top_n <- if (all(is.na(top_flag))) {
    NA_integer_
  } else {
    sum(top_flag == 1, na.rm = TRUE)
  }
  
  overlap_n <- if (all(is.na(bottom_flag))) {
    NA_integer_
  } else {
    sum(
      bottom_flag == 1 & top_flag == 1,
      na.rm = TRUE
    )
  }
  
  bottom_percent <- if (is.na(bottom_n)) {
    NA_real_
  } else {
    100 * bottom_n / n_valid
  }
  
  top_percent <- if (is.na(top_n)) {
    NA_real_
  } else {
    100 * top_n / n_valid
  }
  
  group_status <- if (is.na(bottom_n)) {
    "Scored - no separable percentile tails"
  } else {
    "Scored"
  }
  
  IF_summary <- rbind(
    IF_summary,
    data.table(
      state_season = current_group,
      total_checklists = n_total,
      eligible_checklists = n_valid,
      scored_checklists = n_valid,
      bottom10_cutoff = bottom_cutoff,
      top10_cutoff = top_cutoff,
      bottom10_n = bottom_n,
      top10_n = top_n,
      bottom10_percent = bottom_percent,
      top10_percent = top_percent,
      overlapping_flags = overlap_n,
      scoring_status = group_status
    )
  )
  
  
  # Release model memory before next group
  rm(model_if, X_group, X_valid, scores)
  
  gc(verbose = FALSE)
  
}

# ------------------------------------------------------------
# 9. REMOVE TEMPORARY COLUMNS
# ------------------------------------------------------------

cas[, protocol_factor := NULL]

# ------------------------------------------------------------
# 10. FINAL CHECKS
# ------------------------------------------------------------

cat("\n\n====================================================\n")
cat("ISOLATION FOREST PROCESSING COMPLETED\n")
cat("====================================================\n")

cat("\nTotal checklists:", nrow(cas), "\n")

cat(
  "Scored checklists:",
  sum(!is.na(cas$IF_anomaly_score)),
  "\n"
)

cat(
  "Unscored checklists:",
  sum(is.na(cas$IF_anomaly_score)),
  "\n"
)

cat(
  "Bottom 10% flagged:",
  sum(cas$IF_bottom10_flag == 1, na.rm = TRUE),
  "\n"
)

cat(
  "Top 10% flagged:",
  sum(cas$IF_top10_flag == 1, na.rm = TRUE),
  "\n"
)

cat(
  "Either extreme flagged:",
  sum(cas$IF_either10_flag == 1, na.rm = TRUE),
  "\n"
)

cat(
  "Overlapping flags:",
  sum(
    cas$IF_bottom10_flag == 1 &
      cas$IF_top10_flag == 1,
    na.rm = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 11. GROUP-LEVEL FLAG DISTRIBUTION
# ------------------------------------------------------------

cat("\nGROUP-LEVEL IF SUMMARY:\n")

print(IF_summary)

# ------------------------------------------------------------
# 12. CHECK FLAG DISTRIBUTIONS
# ------------------------------------------------------------

cat("\nBOTTOM 10% FLAG DISTRIBUTION:\n")
print(table(cas$IF_bottom10_flag, useNA = "ifany"))

cat("\nTOP 10% FLAG DISTRIBUTION:\n")
print(table(cas$IF_top10_flag, useNA = "ifany"))

cat("\nCOMBINED FLAG DISTRIBUTION:\n")
print(table(cas$IF_either10_flag, useNA = "ifany"))

# ============================================================
# SAVE IF OUTPUTS TO THE REQUIRED FOLDER
# ============================================================

output_directory <- "FILE PATH"

# Create folder if it does not exist
if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

# Save complete CAS dataset with IF results
data.table::fwrite(
  cas,
  file.path(output_directory, "CAS_LATEST_with_IF.csv")
)

# Save state-season IF summary
data.table::fwrite(
  IF_summary,
  file.path(output_directory, "IF_state_season_summary.csv")
)

# Save unscored checklists, if object exists
if (exists("IF_unscored")) {
  data.table::fwrite(
    IF_unscored,
    file.path(output_directory, "IF_unscored_checklists.csv")
  )
}

cat("Files saved successfully in:\n", output_directory, "\n")

# ============================================================
# AVICRED - ANALYSIS 1
# AGREEMENT BETWEEN CAS AND ISOLATION FOREST
#
# PRIMARY ANALYSES:
#   1. Pearson correlation: CAS vs IF anomaly score
#   2. Spearman correlation: CAS vs IF anomaly score
#   3. ROC-AUC: IF top 10% vs IF bottom 10%
#   4. ROC-AUC: IF top 10% vs remaining 90%
#   5. Sensitivity, specificity, PPV, NPV, balanced accuracy
#   6. CAS class distribution across IF groups
#   7. CAS class-wise IF score distributions
# ============================================================

library(pROC)

# ------------------------------------------------------------
# 1. OUTPUT DIRECTORY
# ------------------------------------------------------------

output_directory <- "FILE PATH/ROC"

if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 2. LOAD DATA IF REQUIRED
# ------------------------------------------------------------

if (!exists("cas")) {
  
  cas <- fread(
    file.path(
      output_directory,
      "CAS_LATEST_with_IF.csv"
    )
  )
  
} else {
  
  cas <- as.data.table(cas)
  
}

# ------------------------------------------------------------
# 3. REMOVE OBSOLETE CAS_SCORE COLUMN
# ------------------------------------------------------------

if ("CAS_score" %in% names(cas)) {
  cas[, CAS_score := NULL]
}

stopifnot("CAS" %in% names(cas))

stopifnot("IF_anomaly_score" %in% names(cas))

stopifnot("IF_top10_flag" %in% names(cas))

stopifnot("IF_bottom10_flag" %in% names(cas))

stopifnot("IF_either10_flag" %in% names(cas))

# ------------------------------------------------------------
# 4. VALIDATE INPUT DATA
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("ANALYSIS 3: CAS AND ISOLATION FOREST AGREEMENT\n")
cat("====================================================\n")

cat("\nTotal checklists:", nrow(cas), "\n")

cat(
  "Scored IF checklists:",
  sum(!is.na(cas$IF_anomaly_score)),
  "\n"
)

cat(
  "Valid CAS values:",
  sum(is.finite(cas$CAS)),
  "\n"
)

cat(
  "CAS missing values:",
  sum(is.na(cas$CAS)),
  "\n"
)

cat(
  "IF missing values:",
  sum(is.na(cas$IF_anomaly_score)),
  "\n"
)

# ------------------------------------------------------------
# 5. CREATE ANALYSIS DATASET
# ------------------------------------------------------------

analysis_data <- cas[
  is.finite(CAS) &
    is.finite(IF_anomaly_score) &
    !is.na(IF_top10_flag) &
    !is.na(IF_bottom10_flag)
]

cat(
  "\nChecklists available for continuous agreement:",
  nrow(analysis_data),
  "\n"
)

# ------------------------------------------------------------
# 6. DEFINE IF GROUPS
# ------------------------------------------------------------

analysis_data[
  ,
  IF_group := fcase(
    IF_top10_flag == 1,
    "Top 10% - Most atypical",
    
    IF_bottom10_flag == 1,
    "Bottom 10% - Most typical",
    
    default = "Middle 80%"
  )
]

analysis_data[
  ,
  IF_group := factor(
    IF_group,
    levels = c(
      "Bottom 10% - Most typical",
      "Middle 80%",
      "Top 10% - Most atypical"
    )
  )
]

# ============================================================
# CONTINUOUS CORRELATION BETWEEN CAS AND IF
# ============================================================

# ------------------------------------------------------------
# 7. PEARSON CORRELATION
# ------------------------------------------------------------

pearson_test <- cor.test(
  analysis_data$CAS,
  analysis_data$IF_anomaly_score,
  method = "pearson",
  conf.level = 0.95
)

pearson_r <- unname(
  pearson_test$estimate
)

pearson_p <- pearson_test$p.value

pearson_ci <- pearson_test$conf.int

pearson_r2 <- pearson_r^2

# ------------------------------------------------------------
# 8. SPEARMAN CORRELATION
# ------------------------------------------------------------

spearman_test <- cor.test(
  analysis_data$CAS,
  analysis_data$IF_anomaly_score,
  method = "spearman",
  exact = FALSE
)

spearman_rho <- unname(
  spearman_test$estimate
)

spearman_p <- spearman_test$p.value

# Approximate 95% CI using Fisher transformation
# This is an asymptotic approximation for Spearman rho.

n_cor <- nrow(analysis_data)

spearman_z <- atanh(
  max(min(spearman_rho, 0.999999), -0.999999)
)

spearman_se <- 1 / sqrt(n_cor - 3)

spearman_ci <- tanh(
  spearman_z + c(-1, 1) * 1.96 * spearman_se
)

# ------------------------------------------------------------
# 9. CORRELATION RESULTS TABLE
# ------------------------------------------------------------

correlation_results <- data.table(
  
  Analysis = c(
    "Pearson correlation",
    "Spearman correlation"
  ),
  
  Correlation = c(
    pearson_r,
    spearman_rho
  ),
  
  CI_95_lower = c(
    pearson_ci[1],
    spearman_ci[1]
  ),
  
  CI_95_upper = c(
    pearson_ci[2],
    spearman_ci[2]
  ),
  
  P_value = c(
    pearson_p,
    spearman_p
  ),
  
  R_squared = c(
    pearson_r2,
    NA_real_
  ),
  
  Sample_size = c(
    n_cor,
    n_cor
  )
  
)

cat("\n\nCONTINUOUS CORRELATION RESULTS:\n")

print(correlation_results)

# ------------------------------------------------------------
# 10. LINEAR REGRESSION
# ------------------------------------------------------------

# IF anomaly score as response
# CAS as explanatory variable

linear_model <- lm(
  IF_anomaly_score ~ CAS,
  data = analysis_data
)

linear_summary <- summary(linear_model)

linear_results <- data.table(
  
  Term = names(coef(linear_model)),
  
  Estimate = as.numeric(
    coef(linear_model)
  ),
  
  Std_Error = as.numeric(
    coef(linear_summary)[, "Std. Error"]
  ),
  
  T_value = as.numeric(
    coef(linear_summary)[, "t value"]
  ),
  
  P_value = as.numeric(
    coef(linear_summary)[, "Pr(>|t|)"]
  )
  
)

linear_model_metrics <- data.table(
  Sample_size = nrow(analysis_data),
  R_squared = linear_summary$r.squared,
  Adjusted_R_squared = linear_summary$adj.r.squared,
  Residual_SE = linear_summary$sigma
)

cat("\n\nLINEAR REGRESSION RESULTS:\n")
print(linear_results)

cat("\nMODEL FIT:\n")
print(linear_model_metrics)

# ============================================================
# ANALYSIS 2
# ROC-AUC: TOP 10% VS BOTTOM 10%
# PRIMARY DISCRIMINATION TEST
# ============================================================

# ------------------------------------------------------------
# 11. CREATE EXTREME-GROUP DATASET
# ------------------------------------------------------------

extreme_data <- analysis_data[
  IF_top10_flag == 1 |
    IF_bottom10_flag == 1
]

extreme_data[
  ,
  IF_reference_group := factor(
    ifelse(
      IF_top10_flag == 1,
      "Top 10% - Most atypical",
      "Bottom 10% - Most typical"
    ),
    levels = c(
      "Bottom 10% - Most typical",
      "Top 10% - Most atypical"
    )
  )
]

cat("\n\nPRIMARY ROC DATASET:\n")

print(
  table(extreme_data$IF_reference_group)
)

# ------------------------------------------------------------
# 12. PRIMARY ROC CURVE
# ------------------------------------------------------------

# The case group is the most atypical IF group.
# Higher CAS is expected to predict this group.

roc_primary <- pROC::roc(
  response = extreme_data$IF_reference_group,
  predictor = extreme_data$CAS,
  levels = c(
    "Bottom 10% - Most typical",
    "Top 10% - Most atypical"
  ),
  direction = "<",
  quiet = TRUE
)

# ------------------------------------------------------------
# 13. PRIMARY AUC AND CONFIDENCE INTERVAL
# ------------------------------------------------------------

auc_primary <- as.numeric(
  pROC::auc(roc_primary)
)

auc_primary_ci <- pROC::ci.auc(
  roc_primary,
  conf.level = 0.95,
  method = "delong"
)

cat("\n\nPRIMARY ROC RESULTS:\n")

cat("AUC:", auc_primary, "\n")

cat(
  "95% CI:",
  as.numeric(auc_primary_ci[1]),
  "-",
  as.numeric(auc_primary_ci[3]),
  "\n"
)

# ------------------------------------------------------------
# 14. OPTIMAL CAS THRESHOLD
# ------------------------------------------------------------

primary_best <- pROC::coords(
  roc_primary,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  transpose = FALSE
)

primary_threshold <- as.numeric(
  primary_best["threshold"]
)

primary_sensitivity <- as.numeric(
  primary_best["sensitivity"]
)

primary_specificity <- as.numeric(
  primary_best["specificity"]
)

# ------------------------------------------------------------
# 15. PRIMARY CLASSIFICATION METRICS
# ------------------------------------------------------------

extreme_data[
  ,
  CAS_predicted_top := as.integer(
    CAS >= primary_threshold
  )
]

extreme_data[
  ,
  IF_actual_top := as.integer(
    IF_top10_flag == 1
  )
]

TP <- sum(
  extreme_data$CAS_predicted_top == 1 &
    extreme_data$IF_actual_top == 1
)

TN <- sum(
  extreme_data$CAS_predicted_top == 0 &
    extreme_data$IF_actual_top == 0
)

FP <- sum(
  extreme_data$CAS_predicted_top == 1 &
    extreme_data$IF_actual_top == 0
)

FN <- sum(
  extreme_data$CAS_predicted_top == 0 &
    extreme_data$IF_actual_top == 1
)

primary_precision <- ifelse(
  TP + FP > 0,
  TP / (TP + FP),
  NA_real_
)

primary_npv <- ifelse(
  TN + FN > 0,
  TN / (TN + FN),
  NA_real_
)

primary_balanced_accuracy <- mean(
  c(
    primary_sensitivity,
    primary_specificity
  )
)

primary_metrics <- data.table(
  
  Metric = c(
    "AUC",
    "AUC_95CI_lower",
    "AUC_95CI_upper",
    "Optimal_CAS_threshold",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NPV",
    "Balanced_accuracy",
    "True_positive",
    "True_negative",
    "False_positive",
    "False_negative"
  ),
  
  Value = c(
    auc_primary,
    as.numeric(auc_primary_ci[1]),
    as.numeric(auc_primary_ci[3]),
    primary_threshold,
    primary_sensitivity,
    primary_specificity,
    primary_precision,
    primary_npv,
    primary_balanced_accuracy,
    TP,
    TN,
    FP,
    FN
  )
  
)

cat("\nPRIMARY CLASSIFICATION METRICS:\n")

print(primary_metrics)

# ============================================================
# ANALYSIS 3
# ROC-AUC: TOP 10% VS REMAINING 90%
# SECONDARY DISCRIMINATION TEST
# ============================================================

# ------------------------------------------------------------
# 16. CREATE SECONDARY ROC DATASET
# ------------------------------------------------------------

secondary_data <- copy(
  analysis_data
)

secondary_data[
  ,
  IF_reference_group := factor(
    ifelse(
      IF_top10_flag == 1,
      "Top 10% - Most atypical",
      "Remaining 90%"
    ),
    levels = c(
      "Remaining 90%",
      "Top 10% - Most atypical"
    )
  )
]

cat("\n\nSECONDARY ROC DATASET:\n")

print(
  table(secondary_data$IF_reference_group)
)

# ------------------------------------------------------------
# 17. SECONDARY ROC CURVE
# ------------------------------------------------------------

roc_secondary <- pROC::roc(
  response = secondary_data$IF_reference_group,
  predictor = secondary_data$CAS,
  levels = c(
    "Remaining 90%",
    "Top 10% - Most atypical"
  ),
  direction = "<",
  quiet = TRUE
)

# ------------------------------------------------------------
# 18. SECONDARY AUC AND CI
# ------------------------------------------------------------

auc_secondary <- as.numeric(
  pROC::auc(roc_secondary)
)

auc_secondary_ci <- pROC::ci.auc(
  roc_secondary,
  conf.level = 0.95,
  method = "delong"
)

cat("\n\nSECONDARY ROC RESULTS:\n")

cat("AUC:", auc_secondary, "\n")

cat(
  "95% CI:",
  as.numeric(auc_secondary_ci[1]),
  "-",
  as.numeric(auc_secondary_ci[3]),
  "\n"
)

# ------------------------------------------------------------
# 19. SECONDARY OPTIMAL THRESHOLD
# ------------------------------------------------------------

secondary_best <- pROC::coords(
  roc_secondary,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  transpose = FALSE
)

secondary_threshold <- as.numeric(
  secondary_best["threshold"]
)

secondary_sensitivity <- as.numeric(
  secondary_best["sensitivity"]
)

secondary_specificity <- as.numeric(
  secondary_best["specificity"]
)

# ------------------------------------------------------------
# 20. SECONDARY CLASSIFICATION METRICS
# ------------------------------------------------------------

secondary_data[
  ,
  CAS_predicted_top := as.integer(
    CAS >= secondary_threshold
  )
]

secondary_data[
  ,
  IF_actual_top := as.integer(
    IF_top10_flag == 1
  )
]

TP2 <- sum(
  secondary_data$CAS_predicted_top == 1 &
    secondary_data$IF_actual_top == 1
)

TN2 <- sum(
  secondary_data$CAS_predicted_top == 0 &
    secondary_data$IF_actual_top == 0
)

FP2 <- sum(
  secondary_data$CAS_predicted_top == 1 &
    secondary_data$IF_actual_top == 0
)

FN2 <- sum(
  secondary_data$CAS_predicted_top == 0 &
    secondary_data$IF_actual_top == 1
)

secondary_precision <- ifelse(
  TP2 + FP2 > 0,
  TP2 / (TP2 + FP2),
  NA_real_
)

secondary_npv <- ifelse(
  TN2 + FN2 > 0,
  TN2 / (TN2 + FN2),
  NA_real_
)

secondary_balanced_accuracy <- mean(
  c(
    secondary_sensitivity,
    secondary_specificity
  )
)

secondary_metrics <- data.table(
  
  Metric = c(
    "AUC",
    "AUC_95CI_lower",
    "AUC_95CI_upper",
    "Optimal_CAS_threshold",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NPV",
    "Balanced_accuracy",
    "True_positive",
    "True_negative",
    "False_positive",
    "False_negative"
  ),
  
  Value = c(
    auc_secondary,
    as.numeric(auc_secondary_ci[1]),
    as.numeric(auc_secondary_ci[3]),
    secondary_threshold,
    secondary_sensitivity,
    secondary_specificity,
    secondary_precision,
    secondary_npv,
    secondary_balanced_accuracy,
    TP2,
    TN2,
    FP2,
    FN2
  )
  
)

cat("\nSECONDARY CLASSIFICATION METRICS:\n")

print(secondary_metrics)

# ============================================================
# ANALYSIS 4
# CAS CLASS DISTRIBUTION ACROSS IF GROUPS
# ============================================================

# ------------------------------------------------------------
# 21. CHECK CAS CLASS LEVELS
# ------------------------------------------------------------

cat("\n\nCAS CLASS DISTRIBUTION:\n")

print(
  table(cas$CAS_class, useNA = "ifany")
)

if (!"CAS_class" %in% names(analysis_data)) {
  stop("CAS_class column is missing.")
}

# Preserve the established CAS class order
cas_class_levels <- c(
  "Very Low",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

analysis_data[
  ,
  CAS_class := factor(
    as.character(CAS_class),
    levels = cas_class_levels,
    ordered = TRUE
  )
]

names(analysis_data)

# ------------------------------------------------------------
# 22. CAS CLASS COUNTS BY IF GROUP
# ------------------------------------------------------------

CAS_IF_counts <- analysis_data[
  !is.na(CAS_class),
  .(
    Number_of_checklists = .N
  ),
  by = .(
    IF_group,
    CAS_class
  )
]

CAS_IF_counts[
  ,
  Group_total := sum(Number_of_checklists),
  by = IF_group
]

CAS_IF_counts[
  ,
  Percentage := 100 *
    Number_of_checklists / Group_total
]

setorder(
  CAS_IF_counts,
  IF_group,
  CAS_class
)

cat("\nCAS CLASS COUNTS BY IF GROUP:\n")

print(CAS_IF_counts)

# ------------------------------------------------------------
# 23. PRIMARY EXTREME-GROUP CAS CLASS DISTRIBUTION
# ------------------------------------------------------------

CAS_extreme_distribution <- CAS_IF_counts[
  IF_group %in% c(
    "Bottom 10% - Most typical",
    "Top 10% - Most atypical"
  )
]

# ------------------------------------------------------------
# 24. STACKED BAR PLOT
# ------------------------------------------------------------

plot_CAS_IF <- ggplot(
  CAS_extreme_distribution,
  aes(
    x = IF_group,
    y = Percentage,
    fill = CAS_class
  )
) +
  
  geom_col(
    position = "fill",
    width = 0.65,
    colour = "black",
    linewidth = 0.2
  ) +
  
  scale_y_continuous(
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  scale_fill_brewer(
    palette = "YlGnBu",
    drop = FALSE
  ) +
  
  labs(
    x = "Isolation Forest reference group",
    y = "Checklists (%)",
    fill = "CAS class"
  ) +
  
  theme_classic(
    base_family = "serif",
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

print(plot_CAS_IF)

ggsave(
  file.path(
    output_directory,
    "Analysis3_CAS_IF_extreme_class_distribution.jpg"
  ),
  plot_CAS_IF,
  width = 8,
  height = 6,
  dpi = 600
)

# ============================================================
# ANALYSIS 5
# CAS CLASS DISTRIBUTION:
# TOP 10% IF ANOMALIES VS ALL REMAINING CHECKLISTS
# ============================================================

# ------------------------------------------------------------
# 25. CAS CLASS DISTRIBUTION:
# TOP 10% IF ANOMALY VS REMAINING 90%
# ------------------------------------------------------------

analysis_data[
  ,
  IF_top10_vs_remaining := factor(
    fifelse(
      IF_top10_flag == TRUE,
      "Top 10% IF anomaly",
      "Remaining 90% checklists"
    ),
    levels = c(
      "Top 10% IF anomaly",
      "Remaining 90% checklists"
    ),
    ordered = TRUE
  )
]

# ------------------------------------------------------------
# 26. DEFINE CAS CLASS ORDER
# ------------------------------------------------------------

cas_class_levels <- c(
  "Very Low",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

analysis_data[
  ,
  CAS_class := factor(
    as.character(CAS_class),
    levels = cas_class_levels,
    ordered = TRUE
  )
]

# ------------------------------------------------------------
# 27. CALCULATE CHECKLIST COUNTS
# ------------------------------------------------------------

CAS_IF_top10_counts <- analysis_data[
  !is.na(CAS_class) &
    !is.na(IF_top10_vs_remaining),
  .(
    Number_of_checklists = .N
  ),
  by = .(
    IF_top10_vs_remaining,
    CAS_class
  )
]

# ------------------------------------------------------------
# 28. CALCULATE GROUP TOTALS AND PERCENTAGES
# ------------------------------------------------------------

CAS_IF_top10_counts[
  ,
  Group_total := sum(Number_of_checklists),
  by = IF_top10_vs_remaining
]

CAS_IF_top10_counts[
  ,
  Percentage := (
    Number_of_checklists / Group_total
  ) * 100
]

# ------------------------------------------------------------
# 29. GROUP TOTAL LABELS
# ------------------------------------------------------------

CAS_IF_top10_totals <- unique(
  CAS_IF_top10_counts[
    ,
    .(
      IF_top10_vs_remaining,
      Group_total
    )
  ]
)

CAS_IF_top10_totals[
  ,
  Group_label := paste0(
    "n = ",
    format(
      Group_total,
      big.mark = ",",
      scientific = FALSE
    )
  )
]

# ------------------------------------------------------------
# 30. GENERATE 100% STACKED BAR PLOT
# ------------------------------------------------------------

plot_CAS_IF_top10 <- ggplot(
  CAS_IF_top10_counts,
  aes(
    x = IF_top10_vs_remaining,
    y = Number_of_checklists,
    fill = CAS_class
  )
) +
  
  geom_col(
    position = "fill",
    width = 0.65,
    colour = "black",
    linewidth = 0.25
  ) +
  
  # Group sample sizes above bars
  geom_text(
    data = CAS_IF_top10_totals,
    aes(
      x = IF_top10_vs_remaining,
      y = 1.035,
      label = Group_label
    ),
    inherit.aes = FALSE,
    family = "serif",
    fontface = "bold",
    size = 4.2
  ) +
  
  scale_fill_brewer(
    palette = "YlGnBu",
    name = "Checklist Atypicality Score"
  ) +
  
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_cartesian(
    ylim = c(0, 1.08),
    clip = "off"
  ) +
  
  labs(
    x = NULL,
    y = "Percentage of checklists"
  ) +
  
  theme_classic(
    base_family = "serif",
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      t = 20,
      r = 15,
      b = 10,
      l = 10
    ),
    
    legend.position = "right"
  )

# ------------------------------------------------------------
# 31. DISPLAY PLOT
# ------------------------------------------------------------

print(plot_CAS_IF_top10)

# ------------------------------------------------------------
# 32. SAVE HIGH-RESOLUTION JPG
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_directory,
    "CAS_Class_Distribution_IF_Top10_vs_Remaining90.jpg"
  ),
  plot = plot_CAS_IF_top10,
  width = 9,
  height = 6,
  units = "in",
  dpi = 400,
  bg = "white"
)

# ============================================================
# ANALYSIS 6
# IF SCORE DISTRIBUTION ACROSS CAS CLASSES
# ============================================================

# ------------------------------------------------------------
# 33. SUMMARISE IF SCORES BY CAS CLASS
# ------------------------------------------------------------

CAS_class_IF_summary <- analysis_data[
  !is.na(CAS_class),
  .(
    Number_of_checklists = .N,
    
    Mean_IF_score = mean(
      IF_anomaly_score,
      na.rm = TRUE
    ),
    
    SD_IF_score = sd(
      IF_anomaly_score,
      na.rm = TRUE
    ),
    
    Median_IF_score = median(
      IF_anomaly_score,
      na.rm = TRUE
    ),
    
    Q1_IF_score = quantile(
      IF_anomaly_score,
      0.25,
      na.rm = TRUE
    ),
    
    Q3_IF_score = quantile(
      IF_anomaly_score,
      0.75,
      na.rm = TRUE
    ),
    
    Minimum_IF_score = min(
      IF_anomaly_score,
      na.rm = TRUE
    ),
    
    Maximum_IF_score = max(
      IF_anomaly_score,
      na.rm = TRUE
    )
    
  ),
  by = CAS_class
]

cat("\nIF SCORE DISTRIBUTION BY CAS CLASS:\n")

print(CAS_class_IF_summary)

# ------------------------------------------------------------
# 34. BOX PLOT
# ------------------------------------------------------------

plot_IF_CAS <- ggplot(
  analysis_data[
    !is.na(CAS_class)
  ],
  aes(
    x = CAS_class,
    y = IF_anomaly_score,
    fill = CAS_class
  )
) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.8
  ) +
  
  # Display outlier observations
  geom_boxplot(
    fill = NA,
    colour = "black",
    outlier.shape = 16,
    outlier.size = 0.5,
    outlier.alpha = 0.35
  ) +
  
  stat_summary(
    fun = median,
    geom = "text",
    aes(
      label = sprintf("%.3f", after_stat(y))
    ),
    family = "Times New Roman",
    colour = "black",
    fontface = "bold",
    size = 3.5,
    vjust = -0.9
  ) +
  
  scale_fill_brewer(
    palette = "YlGnBu",
    drop = FALSE
  ) +
  
  labs(
    x = "CAS class",
    y = "Isolation Forest anomaly score"
  ) +
  
  theme_classic(
    base_family = "serif",
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    legend.position = "none"
  )

print(plot_IF_CAS)

ggsave(
  file.path(output_directory, "IF_anomaly_score_by_CAS_class.jpg"),
  plot = plot_IF_CAS,
  width = 9,
  height = 6,
  dpi = 600
)

# ============================================================
# ANALYSIS 7
# MONOTONIC TREND ACROSS CAS CLASSES
# ============================================================

# ------------------------------------------------------------
# 35. SPEARMAN CORRELATION BETWEEN CAS CLASS AND IF SCORE
# ------------------------------------------------------------

# Retain complete observations only.

class_trend_data <- analysis_data[
  !is.na(CAS_class) &
    !is.na(IF_anomaly_score) &
    is.finite(IF_anomaly_score)
]

# Ensure CAS classes have the correct ordinal order.

cas_class_levels <- c(
  "Very Low",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

class_trend_data[
  ,
  CAS_class := factor(
    as.character(CAS_class),
    levels = cas_class_levels,
    ordered = TRUE
  )
]

# Remove any unmatched class labels.

class_trend_data <- class_trend_data[
  !is.na(CAS_class)
]

# Spearman rank correlation.

class_trend_test <- cor.test(
  x = as.numeric(class_trend_data$CAS_class),
  y = class_trend_data$IF_anomaly_score,
  method = "spearman",
  exact = FALSE
)

class_trend_results <- data.table(
  
  Analysis = "CAS class ordinal vs IF anomaly score",
  
  Spearman_rho = unname(
    class_trend_test$estimate
  ),
  
  P_value = class_trend_test$p.value,
  
  Sample_size = nrow(class_trend_data)
  
)

cat("\nCAS CLASS MONOTONIC TREND:\n")

print(class_trend_results)

# ------------------------------------------------------------
# 36. SAVE SPEARMAN RESULTS
# ------------------------------------------------------------

data.table::fwrite(
  class_trend_results,
  file.path(
    output_directory,
    "CAS_IF_Spearman_Class_Trend.csv"
  )
)

# ============================================================
# ANALYSIS 8
# COMBINE ROC RESULTS
# ============================================================

# ------------------------------------------------------------
# 37. ROC SUMMARY TABLE
# ------------------------------------------------------------

ROC_summary <- data.table(
  
  Comparison = c(
    "IF top 10% vs bottom 10%",
    "IF top 10% vs remaining 90%"
  ),
  
  AUC = c(
    auc_primary,
    auc_secondary
  ),
  
  CI_95_lower = c(
    as.numeric(auc_primary_ci[1]),
    as.numeric(auc_secondary_ci[1])
  ),
  
  CI_95_upper = c(
    as.numeric(auc_primary_ci[3]),
    as.numeric(auc_secondary_ci[3])
  ),
  
  Optimal_CAS_threshold = c(
    primary_threshold,
    secondary_threshold
  ),
  
  Sensitivity = c(
    primary_sensitivity,
    secondary_sensitivity
  ),
  
  Specificity = c(
    primary_specificity,
    secondary_specificity
  ),
  
  PPV = c(
    primary_precision,
    secondary_precision
  ),
  
  NPV = c(
    primary_npv,
    secondary_npv
  ),
  
  Balanced_accuracy = c(
    primary_balanced_accuracy,
    secondary_balanced_accuracy
  ),
  
  Sample_size = c(
    length(roc_primary$response),
    length(roc_secondary$response)
  )
  
)

cat("\n\nFINAL ROC SUMMARY:\n")

print(ROC_summary)

# ============================================================
# ANALYSIS 8: ROC FIGURE
# ============================================================

# ------------------------------------------------------------
# 38. AUC LEGEND LABELS
# ------------------------------------------------------------

auc_label_primary <- sprintf(
  "Top 10%% vs Bottom 10%%: AUC = %.3f (95%% CI: %.3f–%.3f)",
  auc_primary,
  as.numeric(auc_primary_ci[1]),
  as.numeric(auc_primary_ci[3])
)

auc_label_secondary <- sprintf(
  "Top 10%% vs Remaining 90%%: AUC = %.3f (95%% CI: %.3f–%.3f)",
  auc_secondary,
  as.numeric(auc_secondary_ci[1]),
  as.numeric(auc_secondary_ci[3])
)

auc_legend_labels <- c(
  
  "Top 10% vs Bottom 10%" =
    auc_label_primary,
  
  "Top 10% vs Remaining 90%" =
    auc_label_secondary
  
)

# ------------------------------------------------------------
# 39. ROC COORDINATES
# ------------------------------------------------------------

roc_plot_data <- data.table::rbindlist(
  
  list(
    
    data.table(
      FPR = 1 - roc_primary$specificities,
      TPR = roc_primary$sensitivities,
      Comparison = "Top 10% vs Bottom 10%"
    ),
    
    data.table(
      FPR = 1 - roc_secondary$specificities,
      TPR = roc_secondary$sensitivities,
      Comparison = "Top 10% vs Remaining 90%"
    )
    
  )
  
)

roc_plot_data[
  ,
  Comparison := factor(
    Comparison,
    levels = c(
      "Top 10% vs Bottom 10%",
      "Top 10% vs Remaining 90%"
    )
  )
]

# ------------------------------------------------------------
# 40. ROC FIGURE
# ------------------------------------------------------------

plot_ROC <- ggplot(
  
  roc_plot_data,
  
  aes(
    x = FPR,
    y = TPR,
    colour = Comparison
  )
  
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.7
  ) +
  
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  
  scale_colour_manual(
    
    values = c(
      "Top 10% vs Bottom 10%" = "#F8766D",
      "Top 10% vs Remaining 90%" = "#00BFC4"
    ),
    
    breaks = c(
      "Top 10% vs Bottom 10%",
      "Top 10% vs Remaining 90%"
    ),
    
    labels = auc_legend_labels,
    
    name = "ROC comparison"
    
  ) +
  
  guides(
    
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.4
      )
    )
    
  ) +
  
  labs(
    
    x = "1 - Specificity",
    
    y = "Sensitivity"
    
  ) +
  
  theme_classic(
    
    base_family = "Times New Roman",
    base_size = 13
    
  ) +
  
  theme(
    
    text = element_text(
      family = "Times New Roman",
      face = "plain",
      colour = "black"
    ),
    
    axis.title = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      family = "Times New Roman",
      colour = "black",
      size = 12
    ),
    
    # Compact horizontal legend below the graph
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.title = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      family = "Times New Roman",
      size = 9.5,
      colour = "black"
    ),
    
    legend.key.width = unit(
      0.8,
      "cm"
    ),
    
    legend.key.height = unit(
      0.35,
      "cm"
    ),
    
    legend.spacing.x = unit(
      0.35,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.1,
      "cm"
    ),
    
    legend.margin = margin(
      t = 5,
      b = 0
    ),
    
    plot.margin = margin(
      8,
      10,
      5,
      5
    )
    
  )

# ------------------------------------------------------------
# 41. DISPLAY
# ------------------------------------------------------------

print(plot_ROC)

ggsave(
  
  filename = file.path(
    output_directory,
    "CAS_IF_ROC_comparison.jpg"
  ),
  
  plot = plot_ROC,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  device = "jpeg",
  bg = "white"
)

# ------------------------------------------------------------
# 42. SAVE COMBINED ROC SUMMARY
# ------------------------------------------------------------

data.table::fwrite(
  
  ROC_summary,
  
  file.path(
    output_directory,
    "CAS_IF_ROC_Comparison_Summary.csv"
  )
  
)

# ============================================================
# SAVE ALL RESULTS
# ============================================================

# ------------------------------------------------------------
# EXPORT TABLES
# ------------------------------------------------------------

fwrite(
  correlation_results,
  file.path(
    output_directory,
    "Analysis3_CAS_IF_correlations.csv"
  )
)

fwrite(
  linear_results,
  file.path(
    output_directory,
    "Analysis3_CAS_IF_linear_regression.csv"
  )
)

fwrite(
  linear_model_metrics,
  file.path(
    output_directory,
    "Analysis3_CAS_IF_regression_model_fit.csv"
  )
)

fwrite(
  primary_metrics,
  file.path(
    output_directory,
    "Analysis3_IF_top10_vs_bottom10_metrics.csv"
  )
)

fwrite(
  secondary_metrics,
  file.path(
    output_directory,
    "Analysis3_IF_top10_vs_remaining90_metrics.csv"
  )
)

fwrite(
  ROC_summary,
  file.path(
    output_directory,
    "Analysis3_CAS_IF_ROC_summary.csv"
  )
)

fwrite(
  CAS_IF_counts,
  file.path(
    output_directory,
    "Analysis3_CAS_class_by_IF_group.csv"
  )
)

fwrite(
  CAS_class_IF_summary,
  file.path(
    output_directory,
    "Analysis3_IF_score_by_CAS_class.csv"
  )
)

fwrite(
  class_trend_results,
  file.path(
    output_directory,
    "Analysis3_CAS_class_IF_trend.csv"
  )
)

# ------------------------------------------------------------
# SAVE ROC OBJECTS
# ------------------------------------------------------------

saveRDS(
  roc_primary,
  file.path(
    output_directory,
    "Analysis3_ROC_primary.rds"
  )
)

saveRDS(
  roc_secondary,
  file.path(
    output_directory,
    "Analysis3_ROC_secondary.rds"
  )
)

#================ PRE-SENSITIVITY/ SUBSAMPLING DATA CHECKS ==============#

library(data.table)

# ------------------------------------------------------------
# 1. LOAD CAS FILE
# ------------------------------------------------------------

CAS_data <- fread(    # USE THE LATEST CAS FILE WITH IF SCORES PRODUCED FROM ROC ANALYSIS
  file.choose(),
  showProgress = TRUE
)

# Work on a copy
CAS_data <- as.data.table(
  copy(CAS_data)
)

# ------------------------------------------------------------
# 2. REQUIRED VARIABLES
# ------------------------------------------------------------

required_columns <- c(
  "checklist_id",
  "observer_id",
  "state_season",
  "protocol",
  "richness",
  "abundance",
  "CAS",
  "CAS_class"
)

missing_columns <- setdiff(
  required_columns,
  names(CAS_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. RETAIN SUCCESSFULLY SCORED CHECKLISTS
# ------------------------------------------------------------

if ("CAS_scoring_status" %in% names(CAS_data)) {
  
  CAS_scored <- CAS_data[
    CAS_scoring_status == "Scored" &
      !is.na(CAS) &
      is.finite(CAS)
  ]
  
} else {
  
  CAS_scored <- CAS_data[
    !is.na(CAS) &
      is.finite(CAS)
  ]
  
}

# ------------------------------------------------------------
# 4. BASIC VALIDATION
# ------------------------------------------------------------

cat(
  "\nTotal extracted rows:",
  nrow(CAS_data)
)

cat(
  "\nSuccessfully scored checklists:",
  nrow(CAS_scored)
)

cat(
  "\nUnique observers:",
  uniqueN(CAS_scored$observer_id)
)

cat(
  "\nUnique state-season groups:",
  uniqueN(CAS_scored$state_season)
)

cat(
  "\nUnique protocols:\n"
)

print(
  CAS_scored[
    ,
    .N,
    by = protocol
  ]
)

# ------------------------------------------------------------
# 5. CHECK DUPLICATE CHECKLISTS
# ------------------------------------------------------------

duplicate_checklists <- CAS_scored[
  ,
  .N,
  by = checklist_id
][N > 1]

cat(
  "\nDuplicate checklist IDs:",
  nrow(duplicate_checklists)
)

if (nrow(duplicate_checklists) > 0) {
  
  stop(
    "Duplicate checklist IDs detected. ",
    "Resolve before analysis."
  )
  
}

# ------------------------------------------------------------
# 6. CHECK CAS VALIDITY
# ------------------------------------------------------------

cat("\nCAS range:\n")

print(
  CAS_scored[
    ,
    .(
      Number_of_checklists = uniqueN(checklist_id),
      Minimum_CAS = min(CAS),
      Maximum_CAS = max(CAS),
      Mean_CAS = mean(CAS),
      SD_CAS = sd(CAS),
      Median_CAS = median(CAS)
    )
  ]
)

if (any(
  CAS_scored$CAS < 0 |
  CAS_scored$CAS > 1
)) {
  
  stop(
    "CAS values outside the expected 0–1 range."
  )
  
}

# ------------------------------------------------------------
# 7. CHECK CAS CLASS DISTRIBUTION
# ------------------------------------------------------------

cat("\nCAS CLASS DISTRIBUTION:\n")

CAS_class_distribution <- CAS_scored[
  ,
  .(
    Number_of_checklists = uniqueN(checklist_id),
    Percentage = 100 * uniqueN(checklist_id) /
      uniqueN(CAS_scored$checklist_id)
  ),
  by = CAS_class
][order(CAS_class)]

print(CAS_class_distribution)

# ------------------------------------------------------------
# 8. OBSERVER CONTRIBUTION
# ------------------------------------------------------------

observer_summary <- CAS_scored[
  ,
  .(
    Number_of_checklists = .N
  ),
  by = observer_id
]

setorder(
  observer_summary,
  -Number_of_checklists
)

cat("\nObserver contribution summary:\n")

print(
  observer_summary[
    ,
    .(
      Number_of_observers = .N,
      Mean_checklists = mean(Number_of_checklists),
      Median_checklists = median(Number_of_checklists),
      Minimum_checklists = min(Number_of_checklists),
      Maximum_checklists = max(Number_of_checklists)
    )
  ]
)

# ------------------------------------------------------------
# 9. STATE-SEASON × PROTOCOL STRUCTURE
# ------------------------------------------------------------

subsampling_strata <- CAS_scored[
  ,
  .(
    Number_of_checklists = .N,
    Number_of_observers = uniqueN(observer_id),
    Mean_CAS = mean(CAS, na.rm = TRUE),
    SD_CAS = sd(CAS, na.rm = TRUE)
  ),
  by = .(
    state_season,
    protocol
  )
]

setorder(
  subsampling_strata,
  state_season,
  protocol
)

cat("\nState-season × protocol structure:\n")

print(subsampling_strata)

# ------------------------------------------------------------
# 10. EXPORT VALIDATION TABLES
# ------------------------------------------------------------

fwrite(
  observer_summary,
  file.path(
    output_directory,
    "CAS_Observer_Contribution.csv"
  )
)

fwrite(
  subsampling_strata,
  file.path(
    output_directory,
    "CAS_StateSeason_Protocol_Structure.csv"
  )
)

# ============================================================
# CONTEMPORARY CAS
# SUBSAMPLING STRATA SIZE DIAGNOSTIC
# ============================================================

# ------------------------------------------------------------
# 1. STRATA SIZE DISTRIBUTION
# ------------------------------------------------------------

strata_size_summary <- subsampling_strata[
  ,
  .(
    Number_of_strata = .N,
    Mean_checklists = mean(Number_of_checklists),
    Median_checklists = median(Number_of_checklists),
    Minimum_checklists = min(Number_of_checklists),
    Maximum_checklists = max(Number_of_checklists),
    Total_checklists = sum(Number_of_checklists)
  )
]

print(strata_size_summary)

# ------------------------------------------------------------
# 2. GROUP STRATA BY CHECKLIST COUNT
# ------------------------------------------------------------

strata_size_classes <- copy(subsampling_strata)

strata_size_classes[
  ,
  Strata_size_class := cut(
    Number_of_checklists,
    breaks = c(
      0, 1, 5, 10, 20, 50, 100, Inf
    ),
    labels = c(
      "1",
      "2-5",
      "6-10",
      "11-20",
      "21-50",
      "51-100",
      ">100"
    ),
    right = TRUE
  )
]

strata_size_classes <- strata_size_classes[
  ,
  .(
    Number_of_strata = .N,
    Total_checklists = sum(Number_of_checklists),
    Mean_checklists = mean(Number_of_checklists)
  ),
  by = Strata_size_class
][
  order(
    factor(
      Strata_size_class,
      levels = c(
        "1", "2-5", "6-10", "11-20",
        "21-50", "51-100", ">100"
      )
    )
  )
]

print(strata_size_classes)

# ------------------------------------------------------------
# 3. SINGLE-CHECKLIST AND SMALL STRATA
# ------------------------------------------------------------

single_checklist_strata <- subsampling_strata[
  Number_of_checklists == 1
]

print(single_checklist_strata)

small_strata_summary <- subsampling_strata[
  Number_of_checklists <= 5,
  .(
    Number_of_strata = .N,
    Total_checklists = sum(Number_of_checklists)
  )
]

print(small_strata_summary)

# ------------------------------------------------------------
# 4. OBSERVER CONTRIBUTION DISTRIBUTION
# ------------------------------------------------------------

observer_contribution_classes <- copy(observer_summary)

observer_contribution_classes[
  ,
  Contribution_class := cut(
    Number_of_checklists,
    breaks = c(
      0, 1, 2, 5, 10, 20, 50, 100, Inf
    ),
    labels = c(
      "1",
      "2",
      "3-5",
      "6-10",
      "11-20",
      "21-50",
      "51-100",
      ">100"
    ),
    right = TRUE
  )
]

observer_contribution_classes <- observer_contribution_classes[
  ,
  .(
    Number_of_observers = .N,
    Total_checklists = sum(Number_of_checklists),
    Mean_checklists = mean(Number_of_checklists),
    Median_checklists = median(Number_of_checklists)
  ),
  by = Contribution_class
][
  order(
    factor(
      Contribution_class,
      levels = c(
        "1", "2", "3-5", "6-10",
        "11-20", "21-50", "51-100", ">100"
      )
    )
  )
]

print(observer_contribution_classes)

# ------------------------------------------------------------
# 5. EXPORT RESULTS
# ------------------------------------------------------------

fwrite(
  strata_size_summary,
  file.path(
    output_directory,
    "CAS_Strata_Size_Summary.csv"
  )
)

fwrite(
  strata_size_classes,
  file.path(
    output_directory,
    "CAS_Strata_Size_Classes.csv"
  )
)

fwrite(
  single_checklist_strata,
  file.path(
    output_directory,
    "CAS_Single_Checklist_Strata.csv"
  )
)

fwrite(
  small_strata_summary,
  file.path(
    output_directory,
    "CAS_Small_Strata_Summary.csv"
  )
)

fwrite(
  observer_contribution_classes,
  file.path(
    output_directory,
    "CAS_Observer_Contribution_Classes.csv"
  )
)

# ------------------------------------------------------------
# 6. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CONTEMPORARY 10% CAS SAMPLE: STRATA DIAGNOSTICS COMPLETED\n")
cat("============================================================\n")

cat(
  "\nTotal state-season-protocol strata:",
  nrow(subsampling_strata)
)

cat(
  "\nTotal sampled checklists:",
  sum(subsampling_strata$Number_of_checklists)
)

cat(
  "\nTotal observers:",
  nrow(observer_summary)
)

cat(
  "\nStrata containing five or fewer checklists:",
  small_strata_summary$Number_of_strata
)

cat(
  "\nObservers contributing more than 100 checklists:",
  observer_contribution_classes[
    Contribution_class == ">100",
    Number_of_observers
  ]
)

cat("\n\nAll diagnostic tables exported to:\n")
cat(output_directory, "\n")

# ============================================================
# AVICRED
# CAS SENSITIVITY ANALYSIS
#
# Primary contemporary CAS reference:
#   Scaling parameter = 3.5
#   Richness weight   = 0.50
#   Abundance weight  = 0.50
#
# This analysis evaluates:
#   1. Scaling-threshold sensitivity
#   2. Component-weight sensitivity
#   3. Single-component sensitivity
# ============================================================

library(showtext)
library(sysfonts)

font_add(
  family = "Times New Roman",
  regular = "C:/Windows/Fonts/times.ttf",
  bold = "C:/Windows/Fonts/timesbd.ttf",
  italic = "C:/Windows/Fonts/timesi.ttf",
  bolditalic = "C:/Windows/Fonts/timesbi.ttf"
)

showtext_auto()

# ------------------------------------------------------------
# 1. OUTPUT DIRECTORY
# ------------------------------------------------------------

output_directory <- file.path(
  "FILE PATH",
  "sensitivity full"
)

if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 2. PREPARE CONTEMPORARY DATA
# ------------------------------------------------------------

if (!exists("CAS_scored")) {
  stop(
    "Object 'CAS_scored' was not found. ",
    "Load and prepare the contemporary 10% scored dataset first."
  )
}

sensitivity_data <- as.data.table(
  copy(CAS_scored)
)

# ------------------------------------------------------------
# 3. CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required_columns <- c(
  "checklist_id",
  "observer_id",
  "CAS",
  "CAS_class",
  "richness_robust_z",
  "abundance_robust_z"
)

missing_columns <- setdiff(
  required_columns,
  names(sensitivity_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 4. DATA VALIDATION
# ------------------------------------------------------------

if (anyDuplicated(sensitivity_data$checklist_id) > 0) {
  stop(
    "Duplicate checklist IDs detected. ",
    "The sensitivity analysis requires one row per checklist."
  )
}

if (anyNA(sensitivity_data$richness_robust_z) ||
    anyNA(sensitivity_data$abundance_robust_z) ||
    anyNA(sensitivity_data$CAS)) {
  
  stop(
    "Missing CAS or robust z-score values detected."
  )
}

if (any(!is.finite(sensitivity_data$richness_robust_z)) ||
    any(!is.finite(sensitivity_data$abundance_robust_z)) ||
    any(!is.finite(sensitivity_data$CAS))) {
  
  stop(
    "Non-finite CAS or robust z-score values detected."
  )
}

if (any(sensitivity_data$CAS < 0 |
        sensitivity_data$CAS > 1)) {
  
  stop(
    "CAS values outside the expected 0–1 range detected."
  )
}

# ------------------------------------------------------------
# 5. VERIFY PRIMARY CAS FORMULA
# ------------------------------------------------------------

calculated_primary_CAS <- 1 - (
  exp(
    -abs(sensitivity_data$richness_robust_z) / 3.5
  ) +
    exp(
      -abs(sensitivity_data$abundance_robust_z) / 3.5
    )
) / 2

if (!isTRUE(all.equal(
  sensitivity_data$CAS,
  calculated_primary_CAS,
  tolerance = 1e-8
))) {
  
  stop(
    "Stored CAS does not match the primary CAS formula. ",
    "Check the robust z-score columns and calculations."
  )
}

# ------------------------------------------------------------
# 6. VERIFY PRIMARY CLASSIFICATION
# ------------------------------------------------------------

classify_CAS <- function(x) {
  
  cut(
    x,
    breaks = c(
      -Inf,
      0.20,
      0.40,
      0.60,
      0.80,
      Inf
    ),
    labels = c(
      "Very Low",
      "Low",
      "Moderate",
      "High",
      "Very High"
    ),
    right = FALSE,
    ordered_result = TRUE
  )
}

primary_class <- classify_CAS(
  sensitivity_data$CAS
)

if (!all(
  as.character(sensitivity_data$CAS_class) ==
  as.character(primary_class)
)) {
  
  stop(
    "Stored CAS classes do not match the primary ",
    "CAS classification boundaries."
  )
  
} else {
  
  cat(
    "\nPrimary CAS classification verified successfully.",
    "\nClassification agreement: 100%\n"
  )
}

# ------------------------------------------------------------
# 7. SCALING-THRESHOLD SENSITIVITY
# ------------------------------------------------------------

scaling_values <- c(
  2.5,
  3.0,
  3.5,
  4.0,
  4.5
)

scaling_results <- rbindlist(
  lapply(
    scaling_values,
    function(k) {
      
      alternative_CAS <- 1 - (
        exp(
          -abs(sensitivity_data$richness_robust_z) / k
        ) +
          exp(
            -abs(sensitivity_data$abundance_robust_z) / k
          )
      ) / 2
      
      alternative_class <- classify_CAS(
        alternative_CAS
      )
      
      data.table(
        checklist_id = sensitivity_data$checklist_id,
        scaling_parameter = k,
        CAS = alternative_CAS,
        CAS_class = alternative_class
      )
    }
  )
)

# ------------------------------------------------------------
# 8. COMPONENT-WEIGHT SENSITIVITY
# ------------------------------------------------------------

weight_combinations <- data.table(
  richness_weight = c(
    0.40,
    0.50,
    0.60
  ),
  abundance_weight = c(
    0.60,
    0.50,
    0.40
  )
)

weight_results <- rbindlist(
  lapply(
    seq_len(nrow(weight_combinations)),
    function(i) {
      
      wr <- weight_combinations$richness_weight[i]
      wa <- weight_combinations$abundance_weight[i]
      
      alternative_CAS <- 1 - (
        wr * exp(
          -abs(sensitivity_data$richness_robust_z) / 3.5
        ) +
          wa * exp(
            -abs(sensitivity_data$abundance_robust_z) / 3.5
          )
      )
      
      data.table(
        checklist_id = sensitivity_data$checklist_id,
        richness_weight = wr,
        abundance_weight = wa,
        CAS = alternative_CAS,
        CAS_class = classify_CAS(alternative_CAS)
      )
    }
  )
)

# ------------------------------------------------------------
# 9. SINGLE-COMPONENT SENSITIVITY
# ------------------------------------------------------------

single_component_results <- rbindlist(
  list(
    
    data.table(
      checklist_id = sensitivity_data$checklist_id,
      component = "Richness only",
      CAS = 1 - exp(
        -abs(sensitivity_data$richness_robust_z) / 3.5
      )
    ),
    
    data.table(
      checklist_id = sensitivity_data$checklist_id,
      component = "Abundance only",
      CAS = 1 - exp(
        -abs(sensitivity_data$abundance_robust_z) / 3.5
      )
    )
    
  )
)

single_component_results[
  ,
  CAS_class := classify_CAS(CAS)
]

# ------------------------------------------------------------
# 10. SCALING SUMMARY
# ------------------------------------------------------------

scaling_summary <- scaling_results[
  ,
  .(
    Number_of_checklists = .N,
    Mean_CAS = mean(CAS),
    SD_CAS = sd(CAS),
    Median_CAS = median(CAS),
    Minimum_CAS = min(CAS),
    Maximum_CAS = max(CAS),
    Very_Low_pct = mean(CAS_class == "Very Low") * 100,
    Low_pct = mean(CAS_class == "Low") * 100,
    Moderate_pct = mean(CAS_class == "Moderate") * 100,
    High_pct = mean(CAS_class == "High") * 100,
    Very_High_pct = mean(CAS_class == "Very High") * 100
  ),
  by = scaling_parameter
]

# ------------------------------------------------------------
# 11. WEIGHT SUMMARY
# ------------------------------------------------------------

weight_summary <- weight_results[
  ,
  .(
    Number_of_checklists = .N,
    Mean_CAS = mean(CAS),
    SD_CAS = sd(CAS),
    Median_CAS = median(CAS),
    Minimum_CAS = min(CAS),
    Maximum_CAS = max(CAS),
    Very_Low_pct = mean(CAS_class == "Very Low") * 100,
    Low_pct = mean(CAS_class == "Low") * 100,
    Moderate_pct = mean(CAS_class == "Moderate") * 100,
    High_pct = mean(CAS_class == "High") * 100,
    Very_High_pct = mean(CAS_class == "Very High") * 100
  ),
  by = .(
    richness_weight,
    abundance_weight
  )
]

# ------------------------------------------------------------
# 12. SINGLE-COMPONENT SUMMARY
# ------------------------------------------------------------

single_component_summary <- single_component_results[
  ,
  .(
    Number_of_checklists = .N,
    Mean_CAS = mean(CAS),
    SD_CAS = sd(CAS),
    Median_CAS = median(CAS),
    Minimum_CAS = min(CAS),
    Maximum_CAS = max(CAS),
    Very_Low_pct = mean(CAS_class == "Very Low") * 100,
    Low_pct = mean(CAS_class == "Low") * 100,
    Moderate_pct = mean(CAS_class == "Moderate") * 100,
    High_pct = mean(CAS_class == "High") * 100,
    Very_High_pct = mean(CAS_class == "Very High") * 100
  ),
  by = component
]

# ------------------------------------------------------------
# 13. CLASSIFICATION AGREEMENT WITH PRIMARY CAS
# ------------------------------------------------------------

scaling_agreement <- scaling_results[
  ,
  .(
    Agreement_with_primary_pct = mean(
      as.character(CAS_class) ==
        as.character(primary_class)
    ) * 100
  ),
  by = scaling_parameter
]

weight_agreement <- weight_results[
  ,
  .(
    Agreement_with_primary_pct = mean(
      as.character(CAS_class) ==
        as.character(primary_class)
    ) * 100
  ),
  by = .(
    richness_weight,
    abundance_weight
  )
]

# ------------------------------------------------------------
# 14. CORRELATION WITH PRIMARY CAS
# ------------------------------------------------------------

scaling_correlation <- scaling_results[
  ,
  .(
    Pearson_r = cor(
      CAS,
      sensitivity_data$CAS,
      method = "pearson"
    ),
    Spearman_rho = cor(
      CAS,
      sensitivity_data$CAS,
      method = "spearman"
    )
  ),
  by = scaling_parameter
]

weight_correlation <- weight_results[
  ,
  .(
    Pearson_r = cor(
      CAS,
      sensitivity_data$CAS,
      method = "pearson"
    ),
    Spearman_rho = cor(
      CAS,
      sensitivity_data$CAS,
      method = "spearman"
    )
  ),
  by = .(
    richness_weight,
    abundance_weight
  )
]

# ------------------------------------------------------------
# 15. EXPORT CHECKLIST-LEVEL RESULTS
# ------------------------------------------------------------

fwrite(
  scaling_results,
  file.path(
    output_directory,
    "CAS_Sensitivity_Scaling_Checklist_Level.csv"
  )
)

fwrite(
  weight_results,
  file.path(
    output_directory,
    "CAS_Sensitivity_Weights_Checklist_Level.csv"
  )
)

fwrite(
  single_component_results,
  file.path(
    output_directory,
    "CAS_Single_Component_Checklist_Level.csv"
  )
)

# ------------------------------------------------------------
# 16. EXPORT SUMMARY RESULTS
# ------------------------------------------------------------

fwrite(
  scaling_summary,
  file.path(
    output_directory,
    "CAS_Sensitivity_Scaling_Summary.csv"
  )
)

fwrite(
  weight_summary,
  file.path(
    output_directory,
    "CAS_Sensitivity_Weights_Summary.csv"
  )
)

fwrite(
  single_component_summary,
  file.path(
    output_directory,
    "CAS_Single_Component_Summary.csv"
  )
)

# ------------------------------------------------------------
# 17. EXPORT CLASSIFICATION AGREEMENT
# ------------------------------------------------------------

fwrite(
  scaling_agreement,
  file.path(
    output_directory,
    "CAS_Scaling_Class_Agreement.csv"
  )
)

fwrite(
  weight_agreement,
  file.path(
    output_directory,
    "CAS_Weight_Class_Agreement.csv"
  )
)

# ------------------------------------------------------------
# 18. EXPORT CORRELATION RESULTS
# ------------------------------------------------------------

fwrite(
  scaling_correlation,
  file.path(
    output_directory,
    "CAS_Scaling_Correlation.csv"
  )
)

fwrite(
  weight_correlation,
  file.path(
    output_directory,
    "CAS_Weight_Correlation.csv"
  )
)

# ------------------------------------------------------------
# 19. PRINT RESULTS
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("CONTEMPORARY CAS SENSITIVITY ANALYSIS\n")
cat("====================================================\n")

cat(
  "\nNumber of scored checklists:",
  nrow(sensitivity_data),
  "\n"
)

cat(
  "Number of unique observers:",
  uniqueN(sensitivity_data$observer_id),
  "\n"
)

cat(
  "Primary mean CAS:",
  mean(sensitivity_data$CAS),
  "\n"
)

cat(
  "Primary SD CAS:",
  sd(sensitivity_data$CAS),
  "\n"
)

cat(
  "Primary median CAS:",
  median(sensitivity_data$CAS),
  "\n"
)

cat("\n============================================\n")
cat("SCALING-THRESHOLD SENSITIVITY\n")
cat("============================================\n")

print(scaling_summary)

cat("\nCLASSIFICATION AGREEMENT:\n")
print(scaling_agreement)

cat("\nCORRELATION WITH PRIMARY CAS:\n")
print(scaling_correlation)


cat("\n============================================\n")
cat("COMPONENT-WEIGHT SENSITIVITY\n")
cat("============================================\n")

print(weight_summary)

cat("\nCLASSIFICATION AGREEMENT:\n")
print(weight_agreement)

cat("\nCORRELATION WITH PRIMARY CAS:\n")
print(weight_correlation)

cat("\n============================================\n")
cat("SINGLE-COMPONENT SENSITIVITY\n")
cat("============================================\n")

print(single_component_summary)

cat("\n====================================================\n")
cat("CONTEMPORARY 10% CAS SENSITIVITY ANALYSIS COMPLETED\n")
cat("====================================================\n")

cat("\nOutput directory:\n")
cat(output_directory, "\n")

# ============================================================
# AVICRED
# CONTEMPORARY CAS ROBUSTNESS ANALYSIS
# CHECKLIST-LEVEL AND OBSERVER-LEVEL SUBSAMPLING
# HIGH-CONTRIBUTOR EXCLUSION
#
# Analyses:
# 1. Checklist-level random subsampling
# 2. Observer-level cluster subsampling
# 3. High-contributor exclusion
#
# 1,000 repetitions per subsampling level
# Original CAS values are retained without recalculation.
#
# IMPORTANT:
# This analysis examines the stability of CAS distributions
#
# Outputs:
# - Baseline CAS summary
# - Full replicate-level results
# - Subsampling summary tables
# - Observer contribution table
# - High-contributor exclusion results
# - Baseline comparison tables
# ============================================================

# ============================================================
# 1. SETTINGS
# ============================================================

set.seed(20261105)

number_of_repetitions <- 1000

sampling_fractions <- c(
  0.90,
  0.70,
  0.50,
  0.30,
  0.10
)

high_contributor_percentages <- c(
  1,
  2.5,
  5
)

output_directory <- "FILE PATH"

if (!dir.exists(output_directory)) {
  
  dir.create(
    output_directory,
    recursive = TRUE
  )
  
}

# ============================================================
# 2. VALIDATE INPUT DATA
# ============================================================

if (!exists("CAS_scored")) {
  
  stop(
    "CAS_scored does not exist. ",
    "Load the scored contemporary 10% sample first."
  )
  
}

if (!is.data.table(CAS_scored)) {
  
  CAS_scored <- as.data.table(
    CAS_scored
  )
  
}

required_columns <- c(
  "checklist_id",
  "observer_id",
  "CAS"
)

missing_columns <- setdiff(
  required_columns,
  names(CAS_scored)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
  
}

# ============================================================
# 3. PREPARE ANALYSIS DATA
# ============================================================

subsampling_data <- copy(
  CAS_scored
)

# Retain only valid, scored CAS records.
# Original CAS values are not recalculated.

subsampling_data <- subsampling_data[
  !is.na(CAS) &
    is.finite(CAS) &
    !is.na(observer_id) &
    observer_id != "" &
    !is.na(checklist_id) &
    checklist_id != ""
]

# ------------------------------------------------------------
# 3A. Validate CAS range
# ------------------------------------------------------------

if (any(
  subsampling_data$CAS < 0 |
  subsampling_data$CAS > 1
)) {
  
  stop(
    "CAS values outside the expected 0–1 range were detected."
  )
  
}

# ------------------------------------------------------------
# 3B. Check checklist uniqueness
# ------------------------------------------------------------

if (anyDuplicated(subsampling_data$checklist_id) > 0) {
  
  stop(
    "Duplicate checklist IDs detected. ",
    "The analysis requires one row per checklist."
  )
  
}

# ------------------------------------------------------------
# 3C. Define contemporary CAS classes
# ------------------------------------------------------------

CAS_levels <- c(
  "Very Low",
  "Low",
  "Moderate",
  "High",
  "Very High"
)

classify_CAS <- function(x) {
  
  cut(
    x,
    breaks = c(
      -Inf,
      0.20,
      0.40,
      0.60,
      0.80,
      Inf
    ),
    labels = CAS_levels,
    right = FALSE,
    ordered_result = TRUE
  )
  
}

# Recalculate classification labels only.
# The CAS scores themselves remain unchanged.

recomputed_class <- classify_CAS(
  subsampling_data$CAS
)

# ------------------------------------------------------------
# 3D. Verify stored classification, if available
# ------------------------------------------------------------

if ("CAS_class" %in% names(subsampling_data)) {
  
  classification_agreement <- mean(
    as.character(subsampling_data$CAS_class) ==
      as.character(recomputed_class)
  ) * 100
  
  cat(
    "\nStored versus recomputed classification agreement:",
    round(classification_agreement, 5),
    "%\n"
  )
  
  if (classification_agreement < 100) {
    
    stop(
      "Stored CAS_class does not match the contemporary ",
      "classification boundaries. Investigate before proceeding."
    )
    
  }
  
}

subsampling_data[
  ,
  CAS_class := recomputed_class
]

# Confirm that every checklist has a class.

if (anyNA(subsampling_data$CAS_class)) {
  
  stop(
    "Some CAS values could not be assigned to a classification."
  )
  
}


setkey(
  subsampling_data,
  checklist_id
)

# ============================================================
# 4. BASELINE DATASET INFORMATION
# ============================================================

baseline_checklists <- nrow(
  subsampling_data
)

baseline_observers <- uniqueN(
  subsampling_data$observer_id
)

cat("\n")
cat("============================================================\n")
cat("BASELINE CONTEMPORARY 10% CAS DATASET\n")
cat("============================================================\n")

cat(
  "Number of scored checklists:",
  baseline_checklists,
  "\n"
)

cat(
  "Number of unique observers:",
  baseline_observers,
  "\n"
)

cat(
  "Mean CAS:",
  mean(subsampling_data$CAS),
  "\n"
)

cat(
  "SD CAS:",
  sd(subsampling_data$CAS),
  "\n"
)

cat(
  "Median CAS:",
  median(subsampling_data$CAS),
  "\n"
)

# ============================================================
# 5. FUNCTION TO CALCULATE CAS DISTRIBUTION METRICS
# ============================================================

calculate_CAS_metrics <- function(data) {
  
  data <- as.data.table(data)
  
  n_checklists <- nrow(data)
  
  if (n_checklists == 0) {
    stop("Cannot calculate CAS metrics for an empty dataset.")
  }
  
  n_observers <- uniqueN(
    data$observer_id
  )
  
  class_counts <- table(
    factor(
      data$CAS_class,
      levels = CAS_levels
    )
  )
  
  class_percentages <- (
    as.numeric(class_counts) /
      n_checklists
  ) * 100
  
  result <- data.table(
    
    Number_of_checklists = n_checklists,
    
    Number_of_observers = n_observers,
    
    Mean_CAS = mean(data$CAS),
    
    SD_CAS = sd(data$CAS),
    
    Median_CAS = median(data$CAS),
    
    Minimum_CAS = min(data$CAS),
    
    Maximum_CAS = max(data$CAS),
    
    Very_Low_pct = class_percentages[1],
    
    Low_pct = class_percentages[2],
    
    Moderate_pct = class_percentages[3],
    
    High_pct = class_percentages[4],
    
    Very_High_pct = class_percentages[5]
    
  )
  
  return(result)
  
}

# ============================================================
# 6. BASELINE CAS SUMMARY
# ============================================================

baseline_summary <- calculate_CAS_metrics(
  subsampling_data
)

baseline_summary[
  ,
  Analysis := "Original sample"
]

setcolorder(
  baseline_summary,
  c(
    "Analysis",
    setdiff(
      names(baseline_summary),
      "Analysis"
    )
  )
)

# Store baseline statistics explicitly for later comparisons.

baseline_mean <- baseline_summary$Mean_CAS

baseline_sd <- baseline_summary$SD_CAS

baseline_median <- baseline_summary$Median_CAS

cat("\n")
cat("============================================================\n")
cat("BASELINE CAS SUMMARY\n")
cat("============================================================\n")

print(
  baseline_summary
)

fwrite(
  baseline_summary,
  file.path(
    output_directory,
    "CAS_Baseline_Summary.csv"
  )
)

# ============================================================
# 7. CHECKLIST-LEVEL RANDOM SUBSAMPLING
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CHECKLIST-LEVEL RANDOM SUBSAMPLING\n")
cat("============================================================\n")

checklist_replicate_results <- list()

replicate_counter <- 1

for (fraction in sampling_fractions) {
  
  cat(
    "\nSampling fraction:",
    fraction * 100,
    "%\n"
  )
  
  sample_size <- floor(
    baseline_checklists * fraction
  )
  
  for (rep in seq_len(number_of_repetitions)) {
    
    # Randomly select checklist rows.
    
    selected_rows <- sample.int(
      n = baseline_checklists,
      size = sample_size,
      replace = FALSE
    )
    
    sampled_data <- subsampling_data[
      selected_rows
    ]
    
    # Calculate metrics using original CAS values.
    
    metrics <- calculate_CAS_metrics(
      sampled_data
    )
    
    metrics[
      ,
      `:=`(
        Sampling_method = "Checklist-level",
        Sampling_fraction = fraction,
        Sampling_percentage = fraction * 100,
        Repetition = rep
      )
    ]
    
    checklist_replicate_results[[replicate_counter]] <-
      metrics
    
    replicate_counter <- replicate_counter + 1
    
  }
  
  cat(
    "Completed",
    number_of_repetitions,
    "repetitions.\n"
  )
  
}

checklist_replicate_results <- rbindlist(
  checklist_replicate_results
)

setcolorder(
  checklist_replicate_results,
  c(
    "Sampling_method",
    "Sampling_fraction",
    "Sampling_percentage",
    "Repetition",
    setdiff(
      names(checklist_replicate_results),
      c(
        "Sampling_method",
        "Sampling_fraction",
        "Sampling_percentage",
        "Repetition"
      )
    )
  )
)

# ============================================================
# 8. OBSERVER-LEVEL CLUSTER SUBSAMPLING
# ============================================================

cat("\n")
cat("============================================================\n")
cat("OBSERVER-LEVEL CLUSTER SUBSAMPLING\n")
cat("============================================================\n")

all_observers <- unique(
  subsampling_data$observer_id
)

observer_replicate_results <- list()

replicate_counter <- 1

for (fraction in sampling_fractions) {
  
  cat(
    "\nObserver sampling fraction:",
    fraction * 100,
    "%\n"
  )
  
  observer_sample_size <- floor(
    length(all_observers) * fraction
  )
  
  for (rep in seq_len(number_of_repetitions)) {
    
    # Randomly select observers.
    
    selected_observer_indices <- sample.int(
      n = length(all_observers),
      size = observer_sample_size,
      replace = FALSE
    )
    
    selected_observers <- all_observers[
      selected_observer_indices
    ]
    
    # Retain every checklist belonging to selected observers.
    
    sampled_data <- subsampling_data[
      observer_id %chin% selected_observers
    ]
    
    # Calculate metrics using original CAS values.
    
    metrics <- calculate_CAS_metrics(
      sampled_data
    )
    
    metrics[
      ,
      `:=`(
        Sampling_method = "Observer-level",
        Sampling_fraction = fraction,
        Sampling_percentage = fraction * 100,
        Repetition = rep
      )
    ]
    
    observer_replicate_results[[replicate_counter]] <-
      metrics
    
    replicate_counter <- replicate_counter + 1
    
  }
  
  cat(
    "Completed",
    number_of_repetitions,
    "repetitions.\n"
  )
  
}

observer_replicate_results <- rbindlist(
  observer_replicate_results
)

setcolorder(
  observer_replicate_results,
  c(
    "Sampling_method",
    "Sampling_fraction",
    "Sampling_percentage",
    "Repetition",
    setdiff(
      names(observer_replicate_results),
      c(
        "Sampling_method",
        "Sampling_fraction",
        "Sampling_percentage",
        "Repetition"
      )
    )
  )
)

# ============================================================
# 9. COMBINE ALL SUBSAMPLING REPLICATES
# ============================================================

all_subsampling_replicates <- rbindlist(
  list(
    checklist_replicate_results,
    observer_replicate_results
  ),
  use.names = TRUE
)

fwrite(
  all_subsampling_replicates,
  file.path(
    output_directory,
    "CAS_All_Subsampling_Replicates.csv"
  )
)

# ============================================================
# 10. SUMMARISE SUBSAMPLING RESULTS
# ============================================================

metric_columns <- c(
  
  "Number_of_checklists",
  
  "Number_of_observers",
  
  "Mean_CAS",
  
  "SD_CAS",
  
  "Median_CAS",
  
  "Minimum_CAS",
  
  "Maximum_CAS",
  
  "Very_Low_pct",
  
  "Low_pct",
  
  "Moderate_pct",
  
  "High_pct",
  
  "Very_High_pct"
  
)


# Explicit summary function to ensure stable column names.

summarise_replicates <- function(data) {
  
  summary_values <- list()
  
  for (metric in metric_columns) {
    
    x <- data[[metric]]
    
    summary_values[[paste0(metric, "_Mean")]] <-
      mean(x, na.rm = TRUE)
    
    summary_values[[paste0(metric, "_SD_between_repetitions")]] <-
      sd(x, na.rm = TRUE)
    
    summary_values[[paste0(metric, "_Lower_95")]] <-
      as.numeric(
        quantile(
          x,
          0.025,
          na.rm = TRUE,
          names = FALSE
        )
      )
    
    summary_values[[paste0(metric, "_Upper_95")]] <-
      as.numeric(
        quantile(
          x,
          0.975,
          na.rm = TRUE,
          names = FALSE
        )
      )
    
  }
  
  return(
    as.data.table(summary_values)
  )
  
}


subsampling_groups <- unique(
  all_subsampling_replicates[
    ,
    .(
      Sampling_method,
      Sampling_fraction,
      Sampling_percentage
    )
  ]
)

subsampling_summary_list <- vector(
  "list",
  nrow(subsampling_groups)
)

for (i in seq_len(nrow(subsampling_groups))) {
  
  group_method <- subsampling_groups$Sampling_method[i]
  
  group_fraction <- subsampling_groups$Sampling_fraction[i]
  
  group_data <- all_subsampling_replicates[
    Sampling_method == group_method &
      Sampling_fraction == group_fraction
  ]
  
  group_summary <- summarise_replicates(
    group_data
  )
  
  group_summary[
    ,
    `:=`(
      Sampling_method = group_method,
      Sampling_fraction = group_fraction,
      Sampling_percentage =
        subsampling_groups$Sampling_percentage[i]
    )
  ]
  
  subsampling_summary_list[[i]] <- group_summary
  
}


subsampling_summary <- rbindlist(
  subsampling_summary_list,
  fill = TRUE
)

setcolorder(
  subsampling_summary,
  c(
    "Sampling_method",
    "Sampling_fraction",
    "Sampling_percentage",
    setdiff(
      names(subsampling_summary),
      c(
        "Sampling_method",
        "Sampling_fraction",
        "Sampling_percentage"
      )
    )
  )
)

cat("\n")
cat("============================================================\n")
cat("SUBSAMPLING SUMMARY\n")
cat("============================================================\n")

print(
  subsampling_summary
)

fwrite(
  subsampling_summary,
  file.path(
    output_directory,
    "CAS_Subsampling_Summary.csv"
  )
)

# ============================================================
# 11. OBSERVER CONTRIBUTION TABLE
# ============================================================

observer_contribution <- subsampling_data[
  ,
  .(
    Number_of_checklists = .N
  ),
  by = observer_id
]

# Deterministic ordering, including ties.

setorder(
  observer_contribution,
  -Number_of_checklists,
  observer_id
)

observer_contribution[
  ,
  Observer_rank := seq_len(.N)
]

observer_contribution[
  ,
  Checklist_contribution_pct :=
    (Number_of_checklists / baseline_checklists) * 100
]

observer_contribution[
  ,
  Cumulative_checklist_contribution_pct :=
    cumsum(Number_of_checklists) /
    baseline_checklists * 100
]

cat("\n")
cat("============================================================\n")
cat("OBSERVER CONTRIBUTION SUMMARY\n")
cat("============================================================\n")

print(
  observer_contribution[
    ,
    .(
      Number_of_observers = .N,
      Total_checklists = sum(Number_of_checklists),
      Mean_checklists = mean(Number_of_checklists),
      Median_checklists = median(Number_of_checklists),
      Minimum_checklists = min(Number_of_checklists),
      Maximum_checklists = max(Number_of_checklists)
    )
  ]
)

fwrite(
  observer_contribution,
  file.path(
    output_directory,
    "CAS_Observer_Contribution.csv"
  )
)

# ============================================================
# 12. HIGH-CONTRIBUTOR EXCLUSION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("HIGH-CONTRIBUTOR EXCLUSION ANALYSIS\n")
cat("============================================================\n")

high_contributor_results <- list()

high_contributor_observers <- list()

for (percentage in high_contributor_percentages) {
  
  # Number of observers to exclude.
  
  number_to_exclude <- ceiling(
    baseline_observers * percentage / 100
  )
  
  # Identify the highest-contributing observers.
  
  excluded_observers <- observer_contribution[
    seq_len(number_to_exclude),
    observer_id
  ]
  
  # Retain all checklists from the remaining observers.
  
  retained_data <- subsampling_data[
    !observer_id %chin% excluded_observers
  ]
  
  # Calculate the CAS distribution.
  
  metrics <- calculate_CAS_metrics(
    retained_data
  )
  
  metrics[
    ,
    `:=`(
      Exclusion_percentage = percentage,
      Number_of_observers_excluded = number_to_exclude,
      Number_of_checklists_excluded =
        baseline_checklists - nrow(retained_data),
      Percentage_of_checklists_excluded =
        (baseline_checklists - nrow(retained_data)) /
        baseline_checklists * 100
    )
  ]
  
  high_contributor_results[[as.character(percentage)]] <-
    metrics
  
  # Save identities and contribution statistics
  # of excluded observers.
  
  high_contributor_observers[[as.character(percentage)]] <-
    observer_contribution[
      observer_id %chin% excluded_observers
    ][
      ,
      Exclusion_percentage := percentage
    ]
  
  cat("\n")
  cat(
    "Top",
    percentage,
    "% of observers excluded.\n"
  )
  
  cat(
    "Observers excluded:",
    number_to_exclude,
    "\n"
  )
  
  cat(
    "Checklists excluded:",
    baseline_checklists - nrow(retained_data),
    "\n"
  )
  
  cat(
    "Percentage of checklists excluded:",
    round(
      (baseline_checklists - nrow(retained_data)) /
        baseline_checklists * 100,
      2
    ),
    "%\n"
  )
  
}

high_contributor_results <- rbindlist(
  high_contributor_results,
  fill = TRUE
)

high_contributor_observers <- rbindlist(
  high_contributor_observers,
  fill = TRUE
)


# Add original baseline for direct comparison.

baseline_for_exclusion <- copy(
  baseline_summary
)

baseline_for_exclusion[
  ,
  `:=`(
    Exclusion_percentage = 0,
    Number_of_observers_excluded = 0,
    Number_of_checklists_excluded = 0,
    Percentage_of_checklists_excluded = 0
  )
]

high_contributor_results <- rbind(
  baseline_for_exclusion,
  high_contributor_results,
  fill = TRUE
)

setorder(
  high_contributor_results,
  Exclusion_percentage
)

cat("\n")
cat("============================================================\n")
cat("HIGH-CONTRIBUTOR EXCLUSION RESULTS\n")
cat("============================================================\n")

print(
  high_contributor_results
)

fwrite(
  high_contributor_results,
  file.path(
    output_directory,
    "CAS_High_Contributor_Exclusion.csv"
  )
)

fwrite(
  high_contributor_observers,
  file.path(
    output_directory,
    "CAS_Excluded_Observers.csv"
  )
)

# ============================================================
# 13. COMPARISON WITH BASELINE CAS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUBSAMPLING COMPARISON WITH BASELINE\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 13A. Calculate deviations from baseline
# ------------------------------------------------------------

subsampling_deviation <- copy(
  subsampling_summary
)

subsampling_deviation[
  ,
  Mean_CAS_difference_from_baseline :=
    Mean_CAS_Mean - baseline_mean
]

subsampling_deviation[
  ,
  SD_CAS_difference_from_baseline :=
    SD_CAS_Mean - baseline_sd
]

subsampling_deviation[
  ,
  Median_CAS_difference_from_baseline :=
    Median_CAS_Mean - baseline_median
]

# ------------------------------------------------------------
# 13B. Calculate class percentage deviations
# ------------------------------------------------------------

class_names <- c(
  "Very_Low_pct",
  "Low_pct",
  "Moderate_pct",
  "High_pct",
  "Very_High_pct"
)

for (class_name in class_names) {
  
  summary_column <- paste0(
    class_name,
    "_Mean"
  )
  
  deviation_column <- paste0(
    class_name,
    "_difference_from_baseline"
  )
  
  baseline_value <- baseline_summary[[class_name]]
  
  subsampling_deviation[
    ,
    (deviation_column) :=
      get(summary_column) - baseline_value
  ]
  
}

# ------------------------------------------------------------
# 13C. Export subsampling comparison
# ------------------------------------------------------------

fwrite(
  subsampling_deviation,
  file.path(
    output_directory,
    "CAS_Subsampling_Baseline_Comparison.csv"
  )
)

cat("\nSUBSAMPLING BASELINE COMPARISON:\n")

print(
  subsampling_deviation
)


# ----------------------------------------------------
# 13D. HIGH-CONTRIBUTOR EXCLUSION COMPARISON
# ----------------------------------------------------

high_contributor_deviation <- copy(
  high_contributor_results
)

high_contributor_deviation[
  ,
  Mean_CAS_difference_from_baseline :=
    Mean_CAS - baseline_mean
]

high_contributor_deviation[
  ,
  SD_CAS_difference_from_baseline :=
    SD_CAS - baseline_sd
]

high_contributor_deviation[
  ,
  Median_CAS_difference_from_baseline :=
    Median_CAS - baseline_median
]


for (class_name in class_names) {
  
  deviation_column <- paste0(
    class_name,
    "_difference_from_baseline"
  )
  
  baseline_value <- baseline_summary[[class_name]]
  
  high_contributor_deviation[
    ,
    (deviation_column) :=
      get(class_name) - baseline_value
  ]
  
}


# ------------------------------------------------------------
# 13E. Export high-contributor comparison
# ------------------------------------------------------------

fwrite(
  high_contributor_deviation,
  file.path(
    output_directory,
    "CAS_High_Contributor_Baseline_Comparison.csv"
  )
)

cat("\nHIGH-CONTRIBUTOR BASELINE COMPARISON:\n")

print(
  high_contributor_deviation
)


# ============================================================
# 14. FINAL CONSOLE OUTPUT
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("CONTEMPORARY 10% CAS SUBSAMPLING ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "\nOriginal scored checklists:",
  baseline_checklists,
  "\n"
)

cat(
  "Original observers:",
  baseline_observers,
  "\n"
)

cat(
  "Checklist-level repetitions:",
  nrow(checklist_replicate_results),
  "\n"
)

cat(
  "Observer-level repetitions:",
  nrow(observer_replicate_results),
  "\n"
)

cat(
  "Total subsampling replicates:",
  nrow(all_subsampling_replicates),
  "\n"
)

cat(
  "High-contributor exclusion levels:",
  length(high_contributor_percentages),
  "\n"
)

cat(
  "\nOutput directory:\n",
  output_directory,
  "\n"
)

# ============================================================
# HIGH-CONTRIBUTOR CAS DIAGNOSTICS
# Contemporary 10% sample
# Observer-level CAS distribution and contribution
# ============================================================

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

output_directory <- 
  "FILE PATH"

# Create output directory if it does not already exist
if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(output_directory)) {
  stop(
    paste(
      "Output directory does not exist and could not be created:",
      output_directory
    )
  )
}

cat(
  "\nOutput directory confirmed:",
  output_directory,
  "\n"
)

# High-contributor exclusion thresholds
high_contributor_percentages <- c(
  1,
  2.5,
  5
)

# ------------------------------------------------------------
# 2. Validate input
# ------------------------------------------------------------

required_columns <- c(
  "observer_id",
  "CAS"
)

missing_columns <- setdiff(
  required_columns,
  names(subsampling_data)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "The following required columns are missing from subsampling_data:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 3. Apply contemporary CAS class boundaries
# ------------------------------------------------------------
# Contemporary CAS boundaries:
#
# Very Low  = 0.00–<0.20
# Low       = 0.20–<0.40
# Moderate  = 0.40–<0.60
# High      = 0.60–<0.80
# Very High = >=0.80

subsampling_data[
  ,
  CAS_class := cut(
    CAS,
    breaks = c(
      -Inf,
      0.20,
      0.40,
      0.60,
      0.80,
      Inf
    ),
    labels = c(
      "Very Low",
      "Low",
      "Moderate",
      "High",
      "Very High"
    ),
    right = FALSE
  )
]

# ------------------------------------------------------------
# 4. Calculate observer-level CAS statistics
# ------------------------------------------------------------

observer_CAS_diagnostic <- subsampling_data[
  ,
  .(
    Number_of_checklists = .N,
    Mean_CAS = mean(CAS, na.rm = TRUE),
    SD_CAS = sd(CAS, na.rm = TRUE),
    Median_CAS = median(CAS, na.rm = TRUE),
    
    Very_Low_pct =
      mean(CAS_class == "Very Low", na.rm = TRUE) * 100,
    
    Low_pct =
      mean(CAS_class == "Low", na.rm = TRUE) * 100,
    
    Moderate_pct =
      mean(CAS_class == "Moderate", na.rm = TRUE) * 100,
    
    High_pct =
      mean(CAS_class == "High", na.rm = TRUE) * 100,
    
    Very_High_pct =
      mean(CAS_class == "Very High", na.rm = TRUE) * 100
  ),
  by = observer_id
]

# ------------------------------------------------------------
# 5. Create observer contribution ranking
# ------------------------------------------------------------
# Observers are ranked by the number of contemporary
# checklists contributed to the 10% sample.
#
# Rank 1 = observer with the largest number of checklists.

setorder(
  observer_CAS_diagnostic,
  -Number_of_checklists,
  observer_id
)

observer_CAS_diagnostic[
  ,
  Observer_rank := seq_len(.N)
]

total_observers <- nrow(observer_CAS_diagnostic)

# ------------------------------------------------------------
# 6. Calculate observer contribution percentages
# ------------------------------------------------------------

total_checklists <- sum(
  observer_CAS_diagnostic$Number_of_checklists
)

observer_CAS_diagnostic[
  ,
  Checklist_contribution_pct :=
    (Number_of_checklists / total_checklists) * 100
]

observer_CAS_diagnostic[
  ,
  Cumulative_checklist_contribution_pct :=
    cumsum(Checklist_contribution_pct)
]

# ------------------------------------------------------------
# 7. Define contribution groups
# ------------------------------------------------------------
# The groups are based on percentages of observers:
#
# Top 1%
# Top 2.5%, excluding the top 1%
# Top 5%, excluding the top 2.5%
# Remaining observers
#
# For 23,253 observers:
#
# 1%   ≈ 233 observers
# 2.5% ≈ 582 observers
# 5%   ≈ 1,163 observers

top_1_cutoff <- ceiling(
  total_observers * 0.01
)

top_2_5_cutoff <- ceiling(
  total_observers * 0.025
)

top_5_cutoff <- ceiling(
  total_observers * 0.05
)

observer_CAS_diagnostic[
  ,
  Contribution_group := fcase(
    Observer_rank <= top_1_cutoff,
    "Top 1%",
    
    Observer_rank <= top_2_5_cutoff,
    "Top 2.5% (excluding top 1%)",
    
    Observer_rank <= top_5_cutoff,
    "Top 5% (excluding top 2.5%)",
    
    default = "Remaining observers"
  )
]

# ------------------------------------------------------------
# 8. Summarise observer-level CAS characteristics
# ------------------------------------------------------------

observer_group_summary <- observer_CAS_diagnostic[
  ,
  .(
    Number_of_observers = .N,
    
    Total_checklists =
      sum(Number_of_checklists),
    
    Mean_checklists_per_observer =
      mean(Number_of_checklists),
    
    Median_checklists_per_observer =
      median(Number_of_checklists),
    
    Mean_observer_CAS =
      mean(Mean_CAS, na.rm = TRUE),
    
    SD_observer_CAS =
      sd(Mean_CAS, na.rm = TRUE),
    
    Median_observer_CAS =
      median(Mean_CAS, na.rm = TRUE),
    
    Minimum_observer_CAS =
      min(Mean_CAS, na.rm = TRUE),
    
    Maximum_observer_CAS =
      max(Mean_CAS, na.rm = TRUE)
  ),
  by = Contribution_group
]

# Order groups logically
group_order <- c(
  "Top 1%",
  "Top 2.5% (excluding top 1%)",
  "Top 5% (excluding top 2.5%)",
  "Remaining observers"
)

observer_group_summary[
  ,
  Contribution_group :=
    factor(
      Contribution_group,
      levels = group_order
    )
]

setorder(
  observer_group_summary,
  Contribution_group
)

observer_group_summary[
  ,
  Contribution_group :=
    as.character(Contribution_group)
]

cat("\n============================================================\n")
cat("CONTEMPORARY OBSERVER CONTRIBUTION GROUP SUMMARY\n")
cat("============================================================\n")

print(observer_group_summary)

# ------------------------------------------------------------
# 9. Checklist-weighted CAS by contribution group
# ------------------------------------------------------------
# This weights each observer according to the number of
# checklists they contributed.
#
# This is particularly important because an observer with
# 1 checklist should not have the same influence on the
# overall CAS distribution as an observer with 1,000 checklists.

weighted_group_summary <- observer_CAS_diagnostic[
  ,
  .(
    Number_of_observers = .N,
    
    Number_of_checklists =
      sum(Number_of_checklists),
    
    Checklist_contribution_pct =
      sum(Number_of_checklists) /
      total_checklists * 100,
    
    Checklist_weighted_mean_CAS =
      weighted.mean(
        Mean_CAS,
        w = Number_of_checklists,
        na.rm = TRUE
      ),
    
    Checklist_weighted_median_CAS = {
      expanded_CAS <- rep(
        Mean_CAS,
        Number_of_checklists
      )
      
      median(expanded_CAS, na.rm = TRUE)
    }
  ),
  by = Contribution_group
]

weighted_group_summary[
  ,
  Contribution_group :=
    factor(
      Contribution_group,
      levels = group_order
    )
]

setorder(
  weighted_group_summary,
  Contribution_group
)

weighted_group_summary[
  ,
  Contribution_group :=
    as.character(Contribution_group)
]

cat("\n============================================================\n")
cat("CHECKLIST-WEIGHTED CAS BY CONTRIBUTION GROUP\n")
cat("============================================================\n")

print(weighted_group_summary)

# ------------------------------------------------------------
# 10. Overall observer activity and CAS association
# ------------------------------------------------------------

observer_activity_correlation <- data.table(
  Number_of_observers = total_observers,
  
  Spearman_rho = cor(
    observer_CAS_diagnostic$Number_of_checklists,
    observer_CAS_diagnostic$Mean_CAS,
    method = "spearman",
    use = "complete.obs"
  ),
  
  Pearson_r = cor(
    observer_CAS_diagnostic$Number_of_checklists,
    observer_CAS_diagnostic$Mean_CAS,
    method = "pearson",
    use = "complete.obs"
  )
)

cat("\n============================================================\n")
cat("OBSERVER ACTIVITY AND MEAN CAS ASSOCIATION\n")
cat("============================================================\n")

print(observer_activity_correlation)

# ------------------------------------------------------------
# 11. Calculate exact high-contributor exclusion sizes
# ------------------------------------------------------------

high_contributor_exclusion_summary <- rbindlist(
  lapply(
    high_contributor_percentages,
    function(pct) {
      
      n_excluded <- ceiling(
        total_observers * pct / 100
      )
      
      excluded_observers <- observer_CAS_diagnostic[
        Observer_rank <= n_excluded,
        observer_id
      ]
      
      remaining_data <- subsampling_data[
        !observer_id %in% excluded_observers
      ]
      
      data.table(
        Exclusion_percentage = pct,
        Observers_excluded = n_excluded,
        Observers_remaining =
          total_observers - n_excluded,
        
        Checklists_excluded =
          nrow(subsampling_data) -
          nrow(remaining_data),
        
        Checklists_remaining =
          nrow(remaining_data),
        
        Percentage_checklists_excluded =
          (
            nrow(subsampling_data) -
              nrow(remaining_data)
          ) /
          nrow(subsampling_data) * 100,
        
        Remaining_mean_CAS =
          mean(
            remaining_data$CAS,
            na.rm = TRUE
          ),
        
        Remaining_SD_CAS =
          sd(
            remaining_data$CAS,
            na.rm = TRUE
          ),
        
        Remaining_median_CAS =
          median(
            remaining_data$CAS,
            na.rm = TRUE
          )
      )
    }
  )
)

cat("\n============================================================\n")
cat("HIGH-CONTRIBUTOR EXCLUSION SUMMARY\n")
cat("============================================================\n")

print(
  high_contributor_exclusion_summary
)

# ------------------------------------------------------------
# 12. Export observer-level diagnostic table
# ------------------------------------------------------------

fwrite(
  observer_CAS_diagnostic,
  file.path(
    output_directory,
    "CAS_Observer_Level_Diagnostic.csv"
  )
)

# ------------------------------------------------------------
# 13. Export observer contribution group summary
# ------------------------------------------------------------

fwrite(
  observer_group_summary,
  file.path(
    output_directory,
    "CAS_Observer_Group_Summary.csv"
  )
)

# ------------------------------------------------------------
# 14. Export checklist-weighted contribution summary
# ------------------------------------------------------------

fwrite(
  weighted_group_summary,
  file.path(
    output_directory,
    "CAS_Contribution_Weighted_Summary.csv"
  )
)

# ------------------------------------------------------------
# 15. Export observer activity correlation
# ------------------------------------------------------------

fwrite(
  observer_activity_correlation,
  file.path(
    output_directory,
    "CAS_Observer_Activity_Correlation.csv"
  )
)

# ------------------------------------------------------------
# 16. Export high-contributor exclusion results
# ------------------------------------------------------------

fwrite(
  high_contributor_exclusion_summary,
  file.path(
    output_directory,
    "CAS_High_Contributor_Exclusion_Summary.csv"
  )
)

# ------------------------------------------------------------
# 17. Final diagnostic information
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CONTEMPORARY OBSERVER-LEVEL DIAGNOSTICS COMPLETED\n")
cat("============================================================\n")

cat(
  "\nTotal contemporary 10% sample checklists:",
  nrow(subsampling_data),
  "\n"
)

cat(
  "Total observers:",
  total_observers,
  "\n"
)

cat(
  "Top 1% observer cutoff:",
  top_1_cutoff,
  "observers\n"
)

cat(
  "Top 2.5% observer cutoff:",
  top_2_5_cutoff,
  "observers\n"
)

cat(
  "Top 5% observer cutoff:",
  top_5_cutoff,
  "observers\n"
)

cat(
  "\nObserver-level diagnostic files saved to:\n",
  output_directory,
  "\n"
)

# ============================================================
# CAS AND ISOLATION FOREST ANOMALY CHECKLISTS SUMMARY
# ============================================================

# ------------------------------------------------------------
# 1. PREPARE ANALYSIS DATA
# ------------------------------------------------------------

dt <- copy(as.data.table(analysis_data))

dt[, CAS_class := factor(
  as.character(CAS_class),
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  ordered = TRUE
)]

# Retain valid CAS and IF classifications
dt_valid <- dt[
  !is.na(CAS_class) &
    !is.na(either_anomaly_flag) &
    !is.na(IF_top10_flag) &
    !is.na(IF_bottom10_flag)
]

# ------------------------------------------------------------
# 2. TABLE 1: CHECKLIST COUNTS AND ANOMALIES BY CAS CLASS
# ------------------------------------------------------------

CAS_class_summary <- dt_valid[
  ,
  .(
    Total_checklists = .N,
    
    # CAS component-level anomaly flags
    CAS_anomalous_checklists = sum(
      either_anomaly_flag == 1
    ),
    
    CAS_anomalous_percentage = mean(
      either_anomaly_flag == 1
    ) * 100,
    
    # IF top 10%
    IF_top10_checklists = sum(
      IF_top10_flag == 1
    ),
    
    IF_top10_percentage = mean(
      IF_top10_flag == 1
    ) * 100,
    
    # IF bottom 10%
    IF_bottom10_checklists = sum(
      IF_bottom10_flag == 1
    ),
    
    IF_bottom10_percentage = mean(
      IF_bottom10_flag == 1
    ) * 100
  ),
  by = CAS_class
]

CAS_class_summary[
  ,
  Percentage_of_total_checklists :=
    Total_checklists / sum(Total_checklists) * 100
]

setcolorder(
  CAS_class_summary,
  c(
    "CAS_class",
    "Total_checklists",
    "Percentage_of_total_checklists",
    "CAS_anomalous_checklists",
    "CAS_anomalous_percentage",
    "IF_top10_checklists",
    "IF_top10_percentage",
    "IF_bottom10_checklists",
    "IF_bottom10_percentage"
  )
)

print(CAS_class_summary)

# ------------------------------------------------------------
# 3. TABLE 2: CHECKLIST COUNTS AND CAS ANOMALIES BY IF GROUP
# ------------------------------------------------------------

dt_valid[
  ,
  IF_group_summary := fcase(
    IF_top10_flag == 1, "IF Top 10% Anomalies",
    IF_bottom10_flag == 1, "IF Bottom 10%",
    default = "IF Middle 80%"
  )
]

dt_valid[
  ,
  IF_group_summary := factor(
    IF_group_summary,
    levels = c(
      "IF Bottom 10%",
      "IF Middle 80%",
      "IF Top 10% Anomalies"
    )
  )
]

IF_group_summary <- dt_valid[
  ,
  .(
    Total_checklists = .N,
    
    CAS_anomalous_checklists = sum(
      either_anomaly_flag == 1
    ),
    
    CAS_anomalous_percentage = mean(
      either_anomaly_flag == 1
    ) * 100,
    
    CAS_non_anomalous_checklists = sum(
      either_anomaly_flag == 0
    ),
    
    CAS_non_anomalous_percentage = mean(
      either_anomaly_flag == 0
    ) * 100
  ),
  by = IF_group_summary
]

IF_group_summary[
  ,
  Percentage_of_total_checklists :=
    Total_checklists / sum(Total_checklists) * 100
]

setcolorder(
  IF_group_summary,
  c(
    "IF_group_summary",
    "Total_checklists",
    "Percentage_of_total_checklists",
    "CAS_anomalous_checklists",
    "CAS_anomalous_percentage",
    "CAS_non_anomalous_checklists",
    "CAS_non_anomalous_percentage"
  )
)

print(IF_group_summary)

# ------------------------------------------------------------
# 4. SAVE BOTH TABLES
# ------------------------------------------------------------

output_directory <- "FILE PATH"

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  CAS_class_summary,
  file.path(
    output_directory,
    "CAS_Class_Anomaly_Summary.csv"
  )
)

fwrite(
  IF_group_summary,
  file.path(
    output_directory,
    "IF_Group_CAS_Anomaly_Summary.csv"
  )
)

# ====================================================================
# UPSET PLOT: CHECKLIST IF TOP 10% VS CAS VERY HIGH INTERSECTION
# ====================================================================

library(grid)

# ------------------------------------------------------------
# 1. CREATE CHECKLIST-LEVEL MEMBERSHIP DATA
# ------------------------------------------------------------

upset_data1 <- checklists %>%
  select(
    checklist_id,
    CAS_class,
    IF_top10_flag
  ) %>%
  distinct(checklist_id, .keep_all = TRUE) %>%
  filter(
    !is.na(CAS_class),
    !is.na(IF_top10_flag)
  ) %>%
  mutate(
    IF_Top10 = IF_top10_flag == 1,
    CAS_VeryHigh = CAS_class == "Very High"
  ) %>%
  mutate(
    Intersection = case_when(
      IF_Top10 & CAS_VeryHigh ~ "Both",
      IF_Top10 & !CAS_VeryHigh ~ "IF only",
      !IF_Top10 & CAS_VeryHigh ~ "CAS only",
      TRUE ~ "Neither"
    )
  )

# ------------------------------------------------------------
# 2. CHECK THE ACTUAL TOTALS
# ------------------------------------------------------------

cat("\nCAS Very High total:\n")
print(sum(upset_data1$CAS_VeryHigh))

cat("\nIF Top 10% total:\n")
print(sum(upset_data1$IF_Top10))

cat("\nChecklist-level intersection counts:\n")
print(
  upset_data1 %>%
    count(Intersection, name = "Count")
)

# ------------------------------------------------------------
# 3. INTERSECTION COUNTS
# ------------------------------------------------------------

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
# 4. MEMBERSHIP MATRIX
# ------------------------------------------------------------

matrix_data <- expand.grid(
  Intersection = intersection_order,
  Set = c(
    "IF Top 10%",
    "CAS Very High"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Present = case_when(
      Intersection == "Both" ~ TRUE,
      
      Intersection == "IF only" &
        Set == "IF Top 10%" ~ TRUE,
      
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
        "IF Top 10%",
        "CAS Very High"
      )
    )
  )

# ------------------------------------------------------------
# 5. TIMES NEW ROMAN
# ------------------------------------------------------------

windowsFonts(
  TimesNewRoman = windowsFont("Times New Roman")
)

# ------------------------------------------------------------
# 6. INTERSECTION-SIZE BAR PLOT
# ------------------------------------------------------------

p_bars <- ggplot(
  intersection_counts,
  aes(
    x = Intersection,
    y = Count
  )
) +
  geom_col(
    fill = "grey50",
    colour = "black",
    linewidth = 0.4,
    width = 0.68
  ) +
  geom_text(
    aes(
      label = scales::comma(Count)
    ),
    family = "TimesNewRoman",
    size = 14 / ggplot2::.pt,
    colour = "black",
    vjust = -0.45
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    y = "Intersection size",
    x = NULL
  ) +
  theme_classic(
    base_family = "TimesNewRoman",
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      family = "TimesNewRoman",
      size = 14,
      face = "bold",
      hjust = 0.5,
      colour = "black"
    ),
    axis.title.y = element_text(
      family = "TimesNewRoman",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      family = "TimesNewRoman",
      size = 14,
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

# ------------------------------------------------------------
# 7. INTERSECTION MATRIX
# ------------------------------------------------------------

p_matrix <- ggplot(
  matrix_data,
  aes(
    x = Intersection,
    y = Set
  )
) +
  geom_point(
    aes(fill = Present),
    shape = 21,
    size = 5,
    colour = "black",
    stroke = 0.7
  ) +
  geom_line(
    data = matrix_data %>%
      group_by(Intersection) %>%
      filter(sum(Present) > 1),
    aes(
      group = Intersection
    ),
    linewidth = 0.8,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "black",
      "FALSE" = "white"
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_family = "TimesNewRoman",
    base_size = 14
  ) +
  theme(
    axis.title = element_text(
      family = "TimesNewRoman",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      family = "TimesNewRoman",
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank()
  )

# ------------------------------------------------------------
# 8. REASSEMBLE UPSET PLOT
# ------------------------------------------------------------

p_upset <- p_bars / p_matrix +
  patchwork::plot_layout(
    heights = c(3, 1.25)
  )

print(p_upset)

ggsave(
  "Checklist upset.jpg",
  plot = p_upset,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ############### CHECKLIST MAP ###############

library(ggspatial)

############################################################
############## PREPARE MAPPING DATA ########################
############################################################

checklists <- read.csv(
  "FILE PATH/CAS_LATEST_with_IF.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

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
  FROM read_csv(
    'FILE PATH/ebd_filtered.txt',
    delim = '\t',
    header = true,
    ignore_errors = true
  )
  "
)

gis_data_chk <- ebd_gis %>%
  left_join(
    checklists %>%
      select(
        checklist_id,
        CAS,
        CAS_class,
        IF_anomaly_score,
        IF_bottom10_flag,
        IF_top10_flag
      ),
    by = "checklist_id"
  )

############## PREPARE MAPPING VARIABLES ###################

gis_data_chk <- gis_data_chk %>%
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
    
    IF_group = case_when(
      IF_top10_flag == 1 ~ "Top 10%",
      IF_bottom10_flag == 1 ~ "Bottom 10%",
      TRUE ~ "Remaining 80%"
    ),
    
    IF_group = factor(
      IF_group,
      levels = c(
        "Top 10%",
        "Bottom 10%",
        "Remaining 80%"
      )
    )
    
  )

############################################################
############## CONVERT TO SPATIAL DATA #####################
############################################################

gis_sf <- st_as_sf(
  gis_data_chk,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

india_shp <- st_read(   #INSERT RESPECTIVE SHAPEFILE AFTER DOWNLOADING
  "FILE PATH/India_State_Boundary.shp"
)

# Transform to EPSG:3857

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

# Spatially join checklist locations to grid cells

gis_join_chk <- st_join(
  gis_sf,
  grid25
)

############################################################
###################### CAS MAP #############################
############################################################

chk_cas <- gis_join_chk %>%
  st_drop_geometry() %>%
  distinct(
    grid_id,
    checklist_id,
    CAS_class
  ) %>%
  filter(!is.na(CAS_class)) %>%
  count(
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


# Join dominant CAS category to grid

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
      ),
      ordered = TRUE
    )
  )


############################################################
###################### IF MAP ##############################
############################################################

chk_if <- gis_join_chk %>%
  st_drop_geometry() %>%
  distinct(
    grid_id,
    checklist_id,
    IF_top10_flag,
    IF_bottom10_flag
  ) %>%
  filter(
    !is.na(IF_top10_flag),
    !is.na(IF_bottom10_flag)
  ) %>%
  mutate(
    
    IF_class = case_when(
      IF_top10_flag == 1 ~ "Top 10%",
      IF_bottom10_flag == 1 ~ "Bottom 10%",
      TRUE ~ "Remaining 80%"
    ),
    
    IF_class = factor(
      IF_class,
      levels = c(
        "Top 10%",
        "Bottom 10%",
        "Remaining 80%"
      )
    )
    
  ) %>%
  count(
    grid_id,
    IF_class,
    name = "n"
  ) %>%
  group_by(grid_id) %>%
  slice_max(
    n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()


# Join dominant IF category to grid

if_map <- left_join(
  grid25,
  chk_if,
  by = "grid_id"
)

############################################################
############## COMMON MAP THEME AND ELEMENTS ################
############################################################

map_theme <- theme_bw(
  base_family = "Times New Roman"
) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 1
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 15,
      margin = margin(b = 10)
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    legend.key.size = unit(0.9, "cm"),
    
    legend.key.width = unit(1.0, "cm"),
    
    legend.key.height = unit(0.9, "cm"),
    
    legend.spacing.y = unit(0.3, "cm"),
    
    legend.margin = margin(
      15, 18, 15, 18
    ),
    
    legend.box.margin = margin(8, 8, 8, 8),
    
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      colour = "black"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    plot.margin = margin(6, 6, 6, 6)
  )

############################################################
###################### CAS MAP PLOT ########################
############################################################

p_cas <- ggplot() +
  
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
    name = "CAS Class",
    drop = FALSE
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.30,
    text_family = "Times New Roman",
    text_cex = 1.1,
    line_width = 0.9,
    height = unit(0.3, "cm")
  ) +
  
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.4, "cm"),
    width = unit(1.4, "cm")
  ) +
  
  coord_sf(
    crs = st_crs(3857),
    datum = st_crs(4326),
    expand = FALSE,
    label_graticule = "SW"
  ) +
  
  labs(
    title = "(A) Checklist CAS Distribution"
  ) +
  
  map_theme


############################################################
###################### IF MAP PLOT #########################
############################################################

p_if <- ggplot() +
  
  geom_sf(
    data = if_map,
    aes(fill = IF_class),
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
      "Top 10%" = "#020887",
      "Bottom 10%" = "#D1495B",
      "Remaining 80%" = "#8CB369"
    ),
    na.translate = FALSE,
    name = "IF Anomaly Class",
    drop = FALSE
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.30,
    text_family = "Times New Roman",
    text_cex = 1.1,
    line_width = 0.9,
    height = unit(0.3, "cm")
  ) +
  
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.4, "cm"),
    width = unit(1.4, "cm")
  ) +
  
  coord_sf(
    crs = st_crs(3857),
    datum = st_crs(4326),
    expand = FALSE,
    label_graticule = "SW"
  ) +
  
  labs(
    title = "(B) Spatial Distribution of IF Anomaly Classes"
  ) +
  
  map_theme

############################################################
################## COMBINE BOTH MAPS #######################
############################################################

combined_map <- p_cas | p_if

print(combined_map)

############################################################
#################### SAVE FINAL MAP ########################
############################################################

ggsave(
  filename = "CAS_IF_Spatial_Distribution.jpg",
  plot = combined_map,
  width = 20,
  height = 10,
  units = "in",
  dpi = 400,
  device = "jpeg",
  bg = "white"
)

# ============================================================
# RICHNESS AND ABUNDANCE ANOMALY MAPS
# ============================================================

# ============================================================
# 1. ADD RICHNESS AND ABUNDANCE FLAGS
#    TO THE EXISTING GIS DATA
# ============================================================

gis_data_flag <- ebd_gis %>%
  left_join(
    checklists %>%
      select(
        checklist_id,
        richness_anomaly_flag,
        abundance_anomaly_flag
      ),
    by = "checklist_id"
  )

# Convert to sf
gis_sf_flag <- st_as_sf(
  gis_data_flag,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

# Transform to the same CRS already used for grid25
gis_sf_flag <- st_transform(gis_sf_flag, 3857)

# ============================================================
# 2. JOIN CHECKLISTS TO EXISTING 25-KM GRID
# ============================================================

gis_join_flag <- st_join(
  gis_sf_flag,
  grid25
)

# ============================================================
# 3. COUNT UNIQUE RICHNESS-ANOMALOUS CHECKLISTS PER GRID
# ============================================================

richness_count <- gis_join_flag %>%
  filter(
    !is.na(checklist_id),
    !is.na(richness_anomaly_flag),
    richness_anomaly_flag == 1
  ) %>%
  distinct(grid_id, checklist_id) %>%
  count(grid_id, name = "richness_n")


richness_map <- grid25 %>%
  left_join(richness_count, by = "grid_id") %>%
  mutate(
    richness_n = coalesce(richness_n, 0L)
  )

# ============================================================
# 4. COUNT UNIQUE ABUNDANCE-ANOMALOUS CHECKLISTS PER GRID
# ============================================================

abundance_count <- gis_join_flag %>%
  filter(
    !is.na(checklist_id),
    !is.na(abundance_anomaly_flag),
    abundance_anomaly_flag == 1
  ) %>%
  distinct(grid_id, checklist_id) %>%
  count(grid_id, name = "abundance_n")

abundance_map <- grid25 %>%
  left_join(abundance_count, by = "grid_id") %>%
  mutate(
    abundance_n = coalesce(abundance_n, 0L)
  )

# ============================================================
# 5. CHECK ACTUAL RANGES
# ============================================================

cat(
  "\nMaximum richness-anomalous checklists in a grid:",
  max(richness_map$richness_n, na.rm = TRUE),
  "\n"
)

cat(
  "Maximum abundance-anomalous checklists in a grid:",
  max(abundance_map$abundance_n, na.rm = TRUE),
  "\n"
)

# ============================================================
# 6. CREATE SEPARATE BREAKS FOR EACH MAP
# ============================================================

richness_max <- max(
  richness_map$richness_n,
  na.rm = TRUE
)

abundance_max <- max(
  abundance_map$abundance_n,
  na.rm = TRUE
)

# ---- Richness breaks ----

richness_breaks <- c(
  0,
  25,
  50,
  100,
  250,
  richness_max
)

# Remove duplicate values if the maximum is exactly 25, 50, 100, or 250
richness_breaks <- unique(richness_breaks)


# ---- Abundance breaks ----

abundance_breaks <- c(
  0,
  25,
  50,
  100,
  250,
  500,
  abundance_max
)

# Remove duplicate values if the maximum equals one of the break values
abundance_breaks <- unique(abundance_breaks)

# ============================================================
# 7. CONVERT COUNTS TO FACTOR CLASSES
# ============================================================

# ---- Richness classes ----

richness_map <- richness_map %>%
  mutate(
    richness_class = cut(
      richness_n,
      breaks = c(0, 25, 50, 100, 250, Inf),
      labels = c(
        "1–25",
        "25–50",
        "50–100",
        "100–250",
        "250–highest"
      ),
      include.lowest = TRUE,
      right = TRUE
    )
  ) %>%
  mutate(
    richness_class = factor(
      richness_class,
      levels = c(
        "1–25",
        "25–50",
        "50–100",
        "100–250",
        "250–highest"
      )
    )
  )


# ---- Abundance classes ----

abundance_map <- abundance_map %>%
  mutate(
    abundance_class = cut(
      abundance_n,
      breaks = c(0, 25, 50, 100, 250, 500, Inf),
      labels = c(
        "1–25",
        "25–50",
        "50–100",
        "100–250",
        "250–500",
        "500–highest"
      ),
      include.lowest = TRUE,
      right = TRUE
    )
  ) %>%
  mutate(
    abundance_class = factor(
      abundance_class,
      levels = c(
        "1–25",
        "25–50",
        "50–100",
        "100–250",
        "250–500",
        "500–highest"
      )
    )
  )

# ============================================================
# 8. CLIP GRID TO INDIA
# ============================================================

# Transform India boundary to the same CRS
india_3857 <- st_transform(
  india_shp,
  3857
)

# Union all state polygons into one India polygon
india_union <- st_union(india_3857)

# Keep only the portion of each grid cell that falls within India
richness_map <- st_intersection(
  richness_map,
  india_union
)

abundance_map <- st_intersection(
  abundance_map,
  india_union
)

# ============================================================
# 9. COLOUR PALETTES
# ============================================================

richness_colours <- c(
  "#DCEAF7",
  "#A9CDE5",
  "#6BAED6",
  "#2171B5",
  "#3F007D"
)


abundance_colours <- c(
  "#FFF7BC",
  "#FEE391",
  "#FEC44F",
  "#FE9929",
  "#CC4C02",
  "#8C2D04"
)

# ============================================================
# 10. RICHNESS MAP
# ============================================================

p_richness <- ggplot() +

  geom_sf(
    data = richness_map,
    aes(fill = richness_class),
    colour = NA
  ) +

  geom_sf(
    data = india_3857,
    fill = NA,
    colour = "black",
    linewidth = 0.35
  ) +

  scale_fill_manual(
    values = setNames(richness_colours, levels(richness_map$richness_class)),
    drop = TRUE,
    name = "Richness-anomalous\nchecklists"
  ) +

  annotation_scale(
    location = "bl",
    width_hint = 0.3,
    text_family = "Times New Roman"
  ) +

  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering(
      text_family = "Times New Roman"
    )
  ) +

  coord_sf(
    datum = st_crs(4326),
    expand = FALSE
  ) +

  labs(
    title = "Spatial distribution of richness-anomalous checklists"
  ) +

  theme_minimal(base_family = "Times New Roman") +

  theme(
    # Title formatting
    plot.title = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 14,
      hjust = 0
    ),

    legend.position = c(0.88, 0.90),

    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.85, "cm"),

    legend.spacing.y = unit(0.35, "cm"),

    legend.margin = margin(10, 12, 10, 12),

    legend.title = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 15,
      margin = margin(b = 9)
    ),

    legend.text = element_text(
      family = "Times New Roman",
      size = 13
    ),

    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.6
    )
  )

# ============================================================
# 11. ABUNDANCE MAP
# ============================================================

p_abundance <- ggplot() +

  geom_sf(
    data = abundance_map,
    aes(fill = abundance_class),
    colour = NA
  ) +

  geom_sf(
    data = india_3857,
    fill = NA,
    colour = "black",
    linewidth = 0.35
  ) +

  scale_fill_manual(
    values = setNames(abundance_colours, levels(abundance_map$abundance_class)),
    drop = TRUE,
    name = "Abundance-anomalous\nchecklists"
  ) +

  annotation_scale(
    location = "bl",
    width_hint = 0.3,
    text_family = "Times New Roman"
  ) +

  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering(
      text_family = "Times New Roman"
    )
  ) +

  coord_sf(
    datum = st_crs(4326),
    expand = FALSE
  ) +

  labs(
    title = "Spatial distribution of abundance-anomalous checklists"
  ) +

  theme_minimal(base_family = "Times New Roman") +

  theme(
    # Title formatting
    plot.title = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 14,
      hjust = 0
    ),

    legend.position = c(0.88, 0.90),

    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.85, "cm"),

    legend.spacing.y = unit(0.35, "cm"),

    legend.margin = margin(10, 12, 10, 12),

    legend.title = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 15,
      margin = margin(b = 9)
    ),

    legend.text = element_text(
      family = "Times New Roman",
      size = 13
    ),

    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.6
    )
  )


# ============================================================
# 12. DISPLAY SIDE-BY-SIDE
# ============================================================

combined_map <- p_richness | p_abundance

combined_map

# ============================================================
# 13. SAVE
# ============================================================

ggsave(
  "FILE PATH/Richness_Abundance_Anomaly_Maps.jpg",
  combined_map,
  width = 20,
  height = 10,
  units = "in",
  dpi = 400
)

# ============================================================
# SUMMARY OF CHECKLISTS AND ANOMALY FLAGS
# SEPARATELY FOR CAS AND IF CLASSIFICATIONS
# ============================================================

# ------------------------------------------------------------
# 1. CAS CLASSIFICATION SUMMARY
# ------------------------------------------------------------

CAS_summary <- checklists %>%
  filter(!is.na(CAS_class)) %>%
  group_by(CAS_class) %>%
  summarise(
    Number_of_checklists = n(),
    
    Richness_flagged = sum(
      richness_anomaly_flag == 1,
      na.rm = TRUE
    ),
    
    Abundance_flagged = sum(
      abundance_anomaly_flag == 1,
      na.rm = TRUE
    ),
    
    Either_flagged = sum(
      either_anomaly_flag == 1,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
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
  ) %>%
  arrange(CAS_class)

print(CAS_summary)

# ------------------------------------------------------------
# 2. IF CLASSIFICATION SUMMARY
# ------------------------------------------------------------

IF_summary <- checklists %>%
  filter(
    !is.na(IF_top10_flag),
    !is.na(IF_bottom10_flag)
  ) %>%
  mutate(
    IF_classification = case_when(
      IF_top10_flag == 1 ~ "Top 10%",
      IF_bottom10_flag == 1 ~ "Bottom 10%",
      IF_top10_flag == 0 &
        IF_bottom10_flag == 0 ~ "Remaining 80%"
    )
  ) %>%
  group_by(IF_classification) %>%
  summarise(
    Number_of_checklists = n(),
    Richness_flagged = sum(richness_anomaly_flag == 1, na.rm = TRUE),
    Abundance_flagged = sum(abundance_anomaly_flag == 1, na.rm = TRUE),
    Either_flagged = sum(either_anomaly_flag == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    IF_classification = factor(
      IF_classification,
      levels = c("Top 10%", "Remaining 80%", "Bottom 10%")
    )
  ) %>%
  arrange(IF_classification)

print(IF_summary)

# ------------------------------------------------------------
# 3. SAVE BOTH TABLES
# ------------------------------------------------------------

write.csv(
  CAS_summary,
  "CAS_class_anomaly_summary.csv",
  row.names = FALSE
)

write.csv(
  IF_summary,
  "IF_classification_anomaly_summary.csv",
  row.names = FALSE
)

################ CONFUSION MATRIX CAS VS IF ###############

# Prepare checklist-level comparison
confusion_data <- checklists %>%
  distinct(checklist_id, .keep_all = TRUE) %>%
  filter(
    !is.na(IF_top10_flag),
    !is.na(CAS_class)
  ) %>%
  mutate(
    IF_reference = factor(
      ifelse(IF_top10_flag == 1, "Top 10%", "Not Top 10%"),
      levels = c("Not Top 10%", "Top 10%")
    ),
    
    CAS_prediction = factor(
      ifelse(
        CAS_class %in% c("High", "Very High"),
        "High/Very High",
        "Not High"
      ),
      levels = c("Not High", "High/Very High")
    )
  )

# Confusion matrix
confusion_matrix <- table(
  IF = confusion_data$IF_reference,
  CAS = confusion_data$CAS_prediction
)

print(confusion_matrix)

# Extract TP, TN, FP and FN
TN <- confusion_matrix["Not Top 10%", "Not High"]
FP <- confusion_matrix["Not Top 10%", "High/Very High"]
FN <- confusion_matrix["Top 10%", "Not High"]
TP <- confusion_matrix["Top 10%", "High/Very High"]

# Convert confusion matrix counts to numeric
TP <- as.numeric(TP)
TN <- as.numeric(TN)
FP <- as.numeric(FP)
FN <- as.numeric(FN)

# Calculate Matthews Correlation Coefficient
MCC <- (TP * TN - FP * FN) /
  sqrt(
    (TP + FP) *
      (TP + FN) *
      (TN + FP) *
      (TN + FN)
  )

MCC

# Additional classification metrics
accuracy <- (TP + TN) / sum(confusion_matrix)

sensitivity <- TP / (TP + FN)

specificity <- TN / (TN + FP)

precision <- TP / (TP + FP)

# Summary
confusion_summary <- data.frame(
  Metric = c(
    "True Positive",
    "True Negative",
    "False Positive",
    "False Negative",
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Precision",
    "Matthews Correlation Coefficient"
  ),
  Value = c(
    TP, TN, FP, FN,
    accuracy, sensitivity, specificity,
    precision, MCC
  )
)

print(confusion_summary)

# Export results
write.csv(
  as.data.frame.matrix(confusion_matrix),
  "CAS_IF_confusion_matrix.csv"
)

write.csv(
  confusion_summary,
  "CAS_IF_confusion_metrics.csv",
  row.names = FALSE
)

# Prepare binary classifications
confusion_data <- checklists %>%
  filter(
    !is.na(IF_top10_flag),
    !is.na(CAS_class)
  ) %>%
  mutate(
    IF_class = ifelse(IF_top10_flag == 1, "Positive", "Negative"),
    CAS_class_binary = ifelse(
      CAS_class %in% c("High", "Very High"),
      "Positive", "Negative"
    )
  )

# Generate confusion matrix
confusion_matrix <- table(
  CAS = factor(
    confusion_data$CAS_class_binary,
    levels = c("Positive", "Negative")
  ),
  IF = factor(
    confusion_data$IF_class,
    levels = c("Positive", "Negative")
  )
)

print(confusion_matrix)

# Convert to plotting data
cm_plot_data <- as.data.frame(confusion_matrix)

# Calculate percentages
cm_plot_data <- cm_plot_data %>%
  mutate(
    Percentage = Freq / sum(Freq) * 100,
    Label = paste0(
      comma(Freq),
      "\n(",
      sprintf("%.2f", Percentage),
      "%)"
    )
  )

# Plot confusion matrix
p_cm <- ggplot(
  cm_plot_data,
  aes(x = IF, y = CAS, fill = Freq)
) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(
    aes(label = Label),
    family = "Times New Roman",
    fontface = "bold",
    size = 5
  ) +
  scale_fill_gradient(
    low = "#F0F4E8",
    high = "#287A55",
    labels = comma,
    name = "Checklists"
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(cm_plot_data$CAS))) +
  labs(
    x = "Isolation Forest (Top 10% Anomaly)",
    y = "Checklist Atypicality Score (CAS)"
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(p_cm)

# Save at 600 dpi
ggsave(
  "CAS_IF_Confusion_Matrix.jpg",
  plot = p_cm,
  width = 8,
  height = 6,
  dpi = 600
)

################# TREEMAP OF ANOMALOUS RECORDS #############

# ============================================================
# 1. CHECK REQUIRED COLUMNS
# ============================================================

stopifnot(
  all(c(
    "checklist_id",
    "CAS_class",
    "richness_anomaly_flag",
    "abundance_anomaly_flag",
    "IF_top10_flag",
    "IF_bottom10_flag"
  ) %in% names(checklists))
)

# ============================================================
# 2. FILTER NA VALUES AND KEEP ONE RECORD PER CHECKLIST
# ============================================================

checklist_data <- checklists %>%
  dplyr::filter(
    !is.na(checklist_id),
    !is.na(CAS_class),
    !is.na(richness_anomaly_flag),
    !is.na(abundance_anomaly_flag),
    !is.na(IF_top10_flag),
    !is.na(IF_bottom10_flag)
  ) %>%
  dplyr::distinct(
    checklist_id,
    .keep_all = TRUE
  )

# ============================================================
# 3. DEFINE A SINGLE ANOMALOUS-CHECKLIST FLAG
# ============================================================

# A checklist is anomalous if:
#   richness flag = 1 OR abundance flag = 1
#
# Each checklist is counted ONLY ONCE.

checklist_data <- checklist_data %>%
  dplyr::mutate(
    anomalous = dplyr::if_else(
      richness_anomaly_flag == 1 |
        abundance_anomaly_flag == 1,
      1L,
      0L
    )
  )

# ============================================================
# 4. CAS ANOMALOUS CHECKLIST COUNTS
# ============================================================

CAS_treemap <- checklist_data %>%
  dplyr::filter(anomalous == 1) %>%
  dplyr::group_by(CAS_class) %>%
  dplyr::summarise(
    n_anomalous = dplyr::n_distinct(checklist_id),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
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
    )
  ) %>%
  dplyr::arrange(CAS_class)

# ============================================================
# 5. CREATE IF CLASSIFICATION
# ============================================================

checklist_data <- checklist_data %>%
  dplyr::mutate(
    IF_classification = dplyr::case_when(
      IF_top10_flag == 1 ~ "Top 10%",
      IF_bottom10_flag == 1 ~ "Bottom 10%",
      IF_top10_flag == 0 &
        IF_bottom10_flag == 0 ~ "Remaining 80%",
      TRUE ~ NA_character_
    )
  )

# ============================================================
# 6. IF ANOMALOUS CHECKLIST COUNTS
# ============================================================

IF_treemap <- checklist_data %>%
  dplyr::filter(
    anomalous == 1,
    !is.na(IF_classification)
  ) %>%
  dplyr::group_by(IF_classification) %>%
  dplyr::summarise(
    n_anomalous = dplyr::n_distinct(checklist_id),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    IF_classification = factor(
      IF_classification,
      levels = c(
        "Top 10%",
        "Bottom 10%",
        "Remaining 80%"
      ),
      ordered = TRUE
    )
  ) %>%
  dplyr::arrange(IF_classification)

# ============================================================
# 7. PRINT THE COUNTS
# ============================================================

print(CAS_treemap)
print(IF_treemap)

# ============================================================
# 8. CAS TREEMAP
# ============================================================

p_cas <- ggplot(
  CAS_treemap,
  aes(
    area = n_anomalous,
    fill = CAS_class,
    label = paste0(
      as.character(CAS_class),
      "\n",
      scales::comma(n_anomalous),
      " anomalous"
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
    family = "Times New Roman",
    size = 15,
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
    name = "CAS Class"
  ) +
  labs(
    title = "CAS classes",
    subtitle = "Anomalous checklists"
  ) +
  theme_void() +
  theme(
    text = element_text(
      family = "Times New Roman",
      colour = "black"
    ),
    plot.title = element_text(
      family = "Times New Roman",
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      family = "Times New Roman",
      size = 13,
      hjust = 0.5
    ),
    legend.position = "bottom",
    legend.title = element_text(
      family = "Times New Roman",
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(
      family = "Times New Roman",
      size = 12
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

# ============================================================
# 9. IF TREEMAP
# ============================================================

p_if <- ggplot(
  IF_treemap,
  aes(
    area = n_anomalous,
    fill = IF_classification,
    label = paste0(
      as.character(IF_classification),
      "\n",
      scales::comma(n_anomalous),
      " anomalous"
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
    family = "Times New Roman",
    size = 15,
    min.size = 0
  ) +
  scale_fill_manual(
    values = c(
      "Top 10%" = "#922B21",
      "Bottom 10%" = "#F5B041",
      "Remaining 80%" = "#D5D8DC"
    ),
    breaks = c(
      "Top 10%",
      "Bottom 10%",
      "Remaining 80%"
    ),
    drop = FALSE,
    name = "IF Class"
  ) +
  labs(
    title = "Isolation Forest classes",
    subtitle = "Anomalous checklists"
  ) +
  theme_void() +
  theme(
    text = element_text(
      family = "Times New Roman",
      colour = "black"
    ),
    plot.title = element_text(
      family = "Times New Roman",
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      family = "Times New Roman",
      size = 13,
      hjust = 0.5
    ),
    legend.position = "bottom",
    legend.title = element_text(
      family = "Times New Roman",
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(
      family = "Times New Roman",
      size = 12
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

# ============================================================
# 10. COMBINE THE TWO TREEMAPS WITH MORE SPACE BETWEEN THEM
# ============================================================

combined_treemaps <- p_cas +
  patchwork::plot_spacer() +
  p_if +
  patchwork::plot_layout(
    ncol = 3,
    widths = c(1, 0.02, 1)
  )

print(combined_treemaps)

# ============================================================
# 11. SAVE 400 DPI JPG
# ============================================================

ggsave(
  filename = "CAS_IF_anomalous_checklist_treemaps.jpg",
  plot = combined_treemaps,
  width = 18,
  height = 8,
  units = "in",
  dpi = 400,
  device = "jpeg"
)

# ============================================================
# TOP 1% OBSERVER CONTRIBUTOR PROFILE
# ============================================================

# ============================================================
# 1. FILE PATHS
# ============================================================

input_csv <- "FILE PATH"

output_csv <- "FILE PATH"


# ============================================================
# 2. READ CHECKLIST CSV
# ============================================================

checklists <- read_csv(
  input_csv,
  show_col_types = FALSE
)

cat("\nCSV loaded successfully.\n")
cat("Rows:", nrow(checklists), "\n")
cat("Columns:", ncol(checklists), "\n")


# ============================================================
# 3. STANDARDISE ONLY THE VARIABLES WE NEED
# ============================================================
# Using the exact column names from your CSV:
#
# checklist_id
# observer_id
# obs_date
# CAS_class
# richness_anomaly_flag
# abundance_anomaly_flag
# either_anomaly_flag


checklists <- checklists %>%
  mutate(
    
    observer_id = str_trim(
      as.character(observer_id)
    ),
    
    checklist_id = str_trim(
      as.character(checklist_id)
    ),
    
    obs_date = as.Date(obs_date),
    
    active_year = as.integer(
      format(obs_date, "%Y")
    )
  ) %>%
  
  filter(
    !is.na(observer_id),
    observer_id != "",
    !is.na(checklist_id),
    checklist_id != ""
  )


# ============================================================
# 4. STANDARDISE ANOMALY FLAGS
# ============================================================
# Your CSV already contains:
#
# richness_anomaly_flag
# abundance_anomaly_flag
# either_anomaly_flag
#
# This converts them safely to TRUE/FALSE.

flag_to_logical <- function(x) {
  
  if (is.logical(x)) {
    return(
      replace(x, is.na(x), FALSE)
    )
  }
  
  if (is.numeric(x)) {
    return(
      replace(x != 0, is.na(x), FALSE)
    )
  }
  
  x <- str_to_lower(
    str_trim(
      as.character(x)
    )
  )
  
  x %in% c(
    "1",
    "true",
    "t",
    "yes",
    "y"
  )
}


checklists <- checklists %>%
  mutate(
    
    richness_anomaly_flag =
      flag_to_logical(
        richness_anomaly_flag
      ),
    
    abundance_anomaly_flag =
      flag_to_logical(
        abundance_anomaly_flag
      ),
    
    either_anomaly_flag =
      flag_to_logical(
        either_anomaly_flag
      )
  )

# ============================================================
# 5. CHECK FOR DUPLICATE CHECKLIST IDs
# ============================================================

duplicate_checklists <- checklists %>%
  count(
    checklist_id,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1
  )

cat(
  "\nDuplicate checklist IDs:",
  nrow(duplicate_checklists),
  "\n"
)

# ============================================================
# 6. CREATE ONE RECORD PER CHECKLIST
# ============================================================
#
# We use checklist_id as the unit of analysis.
#
# If a checklist appears more than once:
#
# - CAS_class = first non-missing CAS class
# - richness_anomaly_flag = TRUE if any row is TRUE
# - abundance_anomaly_flag = TRUE if any row is TRUE
# - either_anomaly_flag = TRUE if any row is TRUE
# - obs_date = first non-missing date
#
# This prevents duplicate rows from inflating checklist counts.


checklist_level <- checklists %>%
  
  group_by(
    observer_id,
    checklist_id
  ) %>%
  
  summarise(
    
    CAS_class =
      first(
        CAS_class[
          !is.na(CAS_class)
        ],
        default = NA_character_
      ),
    
    richness_anomaly_flag =
      any(
        richness_anomaly_flag,
        na.rm = TRUE
      ),
    
    abundance_anomaly_flag =
      any(
        abundance_anomaly_flag,
        na.rm = TRUE
      ),
    
    either_anomaly_flag =
      any(
        either_anomaly_flag,
        na.rm = TRUE
      ),
    
    obs_date =
      first(
        obs_date[
          !is.na(obs_date)
        ],
        default = as.Date(NA)
      ),
    
    active_year =
      first(
        active_year[
          !is.na(active_year)
        ],
        default = NA_integer_
      ),
    
    .groups = "drop"
  )


# ============================================================
# 7. COUNT CHECKLISTS PER OBSERVER
# ============================================================
#
# This is the definition of contributor status.
#
# TOP 1% =
# the 1% of observers with the highest number of
# submitted checklists.
#
# We are NOT selecting the top 1% of checklists.


observer_contribution <- checklist_level %>%
  
  group_by(
    observer_id
  ) %>%
  
  summarise(
    
    total_checklists =
      n_distinct(checklist_id),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    desc(total_checklists),
    observer_id
  )


# ============================================================
# 8. DETERMINE NUMBER OF TOP 1% OBSERVERS
# ============================================================

n_observers <- nrow(
  observer_contribution
)

n_top1 <- ceiling(
  n_observers * 0.01
)

cat("\n============================================\n")
cat("OBSERVER CONTRIBUTION\n")
cat("============================================\n")

cat(
  "Total observers:",
  n_observers,
  "\n"
)

cat(
  "Top 1% observers:",
  n_top1,
  "\n"
)

# ============================================================
# 9. SELECT TOP 1% CONTRIBUTORS
# ============================================================

top1_observers <- observer_contribution %>%
  
  arrange(
    desc(total_checklists),
    observer_id
  ) %>%
  
  slice_head(
    n = n_top1
  ) %>%
  
  mutate(
    contributor_rank = row_number()
  )

# ============================================================
# 10. EXTRACT ALL CHECKLISTS FROM TOP 1%
# ============================================================

top1_checklists <- checklist_level %>%
  
  inner_join(
    top1_observers %>%
      select(
        observer_id,
        contributor_rank,
        total_checklists
      ),
    by = "observer_id"
  )


cat(
  "\nChecklists from top 1% observers:",
  nrow(top1_checklists),
  "\n"
)

# ============================================================
# 11. COUNT CAS CLASSES PER OBSERVER
# ============================================================
#
# Your actual CAS classes are:
#
# Very Low
# Low
# Moderate
# High
# Very High
#
# Standardise them explicitly so that blank or unexpected
# values cannot become column names in pivot_wider().
# ============================================================

cas_counts <- top1_checklists %>%
  
  mutate(
    
    CAS_class = str_to_lower(
      str_trim(
        as.character(CAS_class)
      )
    ),
    
    CAS_class = case_when(
      
      CAS_class == "very low" ~ "CAS_Very_Low",
      
      CAS_class == "low" ~ "CAS_Low",
      
      CAS_class == "moderate" ~ "CAS_Moderate",
      
      CAS_class == "high" ~ "CAS_High",
      
      CAS_class == "very high" ~ "CAS_Very_High",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  # Remove missing / blank / unrecognised CAS classes
  filter(
    !is.na(CAS_class),
    CAS_class != ""
  ) %>%
  
  count(
    observer_id,
    CAS_class,
    name = "n_checklists"
  ) %>%
  
  pivot_wider(
    names_from = CAS_class,
    values_from = n_checklists,
    values_fill = 0
  )

# ============================================================
# 12. ENSURE ALL FIVE CAS CLASS COLUMNS EXIST
# ============================================================

required_CAS_columns <- c(
  "CAS_Very_Low",
  "CAS_Low",
  "CAS_Moderate",
  "CAS_High",
  "CAS_Very_High"
)

for (col in required_CAS_columns) {
  
  if (!col %in% names(cas_counts)) {
    
    cas_counts[[col]] <- 0L
    
  }
}

# ============================================================
# 13. CHECK CAS COUNTS
# ============================================================

cat("\n============================================\n")
cat("CAS CLASS COUNTS IN TOP 1% OBSERVERS\n")
cat("============================================\n")

print(
  cas_counts %>%
    select(
      observer_id,
      all_of(required_CAS_columns)
    ) %>%
    head(10)
)

# ============================================================
# 14. COUNT ANOMALOUS CHECKLISTS PER OBSERVER
# ============================================================

anomaly_counts <- top1_checklists %>%
  
  group_by(
    observer_id
  ) %>%
  
  summarise(
    
    richness_anomalous_checklists =
      sum(
        richness_anomaly_flag,
        na.rm = TRUE
      ),
    
    abundance_anomalous_checklists =
      sum(
        abundance_anomaly_flag,
        na.rm = TRUE
      ),
    
    either_anomalous_checklists =
      sum(
        either_anomaly_flag,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

# ============================================================
# 15. CALCULATE ACTIVE YEARS
# ============================================================
#
# Active years =
# number of DISTINCT calendar years in which the observer
# submitted at least one checklist.
#
# Example:
#
# 2012, 2013, 2015, 2019, 2020
#
# = 5 active years.
#
# ============================================================

active_years <- top1_checklists %>%
  
  filter(
    !is.na(active_year)
  ) %>%
  
  group_by(
    observer_id
  ) %>%
  
  summarise(
    
    active_years =
      n_distinct(active_year),
    
    first_active_year =
      min(active_year),
    
    last_active_year =
      max(active_year),
    
    .groups = "drop"
  )

# ============================================================
# 16. COMBINE OBSERVER-LEVEL RESULTS
# ============================================================

top1_observer_profile <- top1_observers %>%
  
  left_join(
    cas_counts,
    by = "observer_id"
  ) %>%
  
  left_join(
    anomaly_counts,
    by = "observer_id"
  ) %>%
  
  left_join(
    active_years,
    by = "observer_id"
  ) %>%
  
  mutate(
    
    CAS_Very_Low =
      coalesce(
        CAS_Very_Low,
        0L
      ),
    
    CAS_Low =
      coalesce(
        CAS_Low,
        0L
      ),
    
    CAS_Moderate =
      coalesce(
        CAS_Moderate,
        0L
      ),
    
    CAS_High =
      coalesce(
        CAS_High,
        0L
      ),
    
    CAS_Very_High =
      coalesce(
        CAS_Very_High,
        0L
      ),
    
    richness_anomalous_checklists =
      coalesce(
        richness_anomalous_checklists,
        0L
      ),
    
    abundance_anomalous_checklists =
      coalesce(
        abundance_anomalous_checklists,
        0L
      ),
    
    either_anomalous_checklists =
      coalesce(
        either_anomalous_checklists,
        0L
      ),
    
    active_years =
      coalesce(
        active_years,
        0L
      )
  ) %>%
  
  select(
    
    contributor_rank,
    observer_id,
    total_checklists,
    
    CAS_Very_Low,
    CAS_Low,
    CAS_Moderate,
    CAS_High,
    CAS_Very_High,
    
    richness_anomalous_checklists,
    abundance_anomalous_checklists,
    either_anomalous_checklists,
    
    active_years,
    first_active_year,
    last_active_year
  ) %>%
  
  arrange(
    contributor_rank
  )

# ============================================================
# 17. FINAL CHECKS
# ============================================================

cat("\n============================================\n")
cat("FINAL TOP 1% OBSERVER PROFILE\n")
cat("============================================\n")

cat(
  "Top 1% observers:",
  nrow(top1_observer_profile),
  "\n"
)

cat(
  "Total checklists represented:",
  sum(
    top1_observer_profile$total_checklists
  ),
  "\n"
)

cat("\nFirst 10 observers:\n\n")

print(
  head(
    top1_observer_profile,
    10
  )
)

# ============================================================
# 18. CHECK CAS COUNTS
# ============================================================

top1_observer_profile <- top1_observer_profile %>%
  
  mutate(
    
    CAS_class_total =
      CAS_Very_Low +
      CAS_Low +
      CAS_Moderate +
      CAS_High +
      CAS_Very_High
  )

cat("\n============================================\n")
cat("CAS COUNT CHECK\n")
cat("============================================\n")

print(
  summary(
    top1_observer_profile$CAS_class_total
  )
)

# ============================================================
# 19. SAVE FINAL CSV
# ============================================================

write_csv(
  top1_observer_profile,
  output_csv
)

names(top1_observer_profile)

cat("\n============================================\n")
cat("OUTPUT SAVED\n")
cat("============================================\n")

cat(
  output_csv,
  "\n"
)

#================== CAS CLASSES ACORSS OBSERVER ACTIVE YEARS ================#

# ============================================================
# 1. READ DATA
# ============================================================

top1_observer_profile <- read_csv(
  "FILE PATH/top_1_percent_observer_profile.csv",
  show_col_types = FALSE
)

# ============================================================
# 2. CALCULATE ANOMALY RATES AND CAS PROPORTIONS
# ============================================================

top1_observer_profile <- top1_observer_profile %>%
  mutate(
    either_anomaly_rate =
      either_anomalous_checklists / total_checklists,
    
    richness_anomaly_rate =
      richness_anomalous_checklists / total_checklists,
    
    abundance_anomaly_rate =
      abundance_anomalous_checklists / total_checklists,
    
    CAS_High_Very_High =
      CAS_High + CAS_Very_High,
    
    CAS_High_Very_High_rate =
      CAS_High_Very_High / total_checklists
  )

# ============================================================
# 3. DEFINE EXPERIENCE GROUPS
# ============================================================

top1_observer_profile <- top1_observer_profile %>%
  mutate(
    experience_group = case_when(
      active_years >= 1 & active_years <= 5 ~ "1–5 years",
      active_years >= 6 & active_years <= 10 ~ "6–10 years",
      active_years >= 11 & active_years <= 15 ~ "11–15 years",
      active_years >= 16 & active_years <= 20 ~ "16–20 years",
      active_years > 20 ~ ">20 years",
      TRUE ~ NA_character_
    ),
    
    experience_group = factor(
      experience_group,
      levels = c(
        "1–5 years",
        "6–10 years",
        "11–15 years",
        "16–20 years",
        ">20 years"
      )
    )
  )

# ============================================================
# 4. CONVERT CAS COUNTS TO LONG FORMAT
# ============================================================

cas_long <- top1_observer_profile %>%
  filter(!is.na(experience_group)) %>%
  group_by(experience_group) %>%
  summarise(
    Very_Low = sum(CAS_Very_Low, na.rm = TRUE),
    Low = sum(CAS_Low, na.rm = TRUE),
    Moderate = sum(CAS_Moderate, na.rm = TRUE),
    High = sum(CAS_High, na.rm = TRUE),
    Very_High = sum(CAS_Very_High, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Very_Low, Low, Moderate, High, Very_High),
    names_to = "CAS_class",
    values_to = "checklist_count"
  ) %>%
  group_by(experience_group) %>%
  mutate(
    proportion = checklist_count / sum(checklist_count) * 100
  ) %>%
  ungroup()

cas_long$CAS_class <- factor(
  cas_long$CAS_class,
  levels = c(
    "Very_Low", "Low", "Moderate", "High", "Very_High"
  ),
  labels = c(
    "Very Low", "Low", "Moderate", "High", "Very High"
  )
)

# ============================================================
# 5. STACKED BAR PLOT
# ============================================================

p_stacked <- ggplot(
  cas_long,
  aes(
    x = experience_group,
    y = proportion,
    fill = CAS_class
  )
) +
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(
      label = paste0(round(proportion, 1), "%")
    ),
    position = position_stack(vjust = 0.5),
    size = 3.8,
    fontface = "bold",
    family = "Times New Roman",
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Very Low" = "#2166AC",
      "Low" = "#67A9CF",
      "Moderate" = "#D1E5F0",
      "High" = "#EF8A62",
      "Very High" = "#B2182B"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Observer experience (active years)",
    y = "Percentage of checklists",
    fill = "CAS class"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    axis.title.x = element_text(
      size = 15,
      face = "bold",
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 15,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      colour = "black"
    ),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      size = 12
    ),
    legend.position = "right",
    legend.key.size = unit(0.7, "cm")
  )

p_stacked

# ============================================================
# 6. SAVE STACKED BAR PLOT
# ============================================================

ggsave(
  "FILE PATH/CAS_composition_experience.jpg",
  plot = p_stacked,
  width = 9,
  height = 6,
  dpi = 400
)

#============= SCATTER PLOT ANOMALY RATE VS CHECKLISTS/YEARS =============#

library(mgcv)

# ============================================================
# 1. CORRELATION TESTS
# ============================================================

# Panel A: Active years vs anomaly rate

pearson_A <- cor.test(
  top1_observer_profile$active_years,
  top1_observer_profile$either_anomaly_rate,
  method = "pearson"
)

spearman_A <- cor.test(
  top1_observer_profile$active_years,
  top1_observer_profile$either_anomaly_rate,
  method = "spearman",
  exact = FALSE
)


# Panel B: Log total checklists vs anomaly rate

pearson_B <- cor.test(
  log10(top1_observer_profile$total_checklists),
  top1_observer_profile$either_anomaly_rate,
  method = "pearson"
)

spearman_B <- cor.test(
  top1_observer_profile$total_checklists,
  top1_observer_profile$either_anomaly_rate,
  method = "spearman",
  exact = FALSE
)

# ============================================================
# 2. FUNCTION FOR SIGNIFICANCE SYMBOLS
# ============================================================

significance_symbol <- function(p) {
  if (p < 0.05) {
    "*"
  } else {
    "ns"
  }
}

# ============================================================
# 3. FORMAT CORRELATION RESULTS
# ============================================================

correlation_label <- function(pearson, spearman) {
  
  paste0(
    "Pearson r = ", round(unname(pearson$estimate), 3),
    ", p = ", format.pval(pearson$p.value, digits = 3, eps = 0.001),
    " ", significance_symbol(pearson$p.value),
    
    "\n",
    
    "Spearman \u03C1 = ", round(unname(spearman$estimate), 3),
    ", p = ", format.pval(spearman$p.value, digits = 3, eps = 0.001),
    " ", significance_symbol(spearman$p.value)
  )
}

label_A <- correlation_label(pearson_A, spearman_A)
label_B <- correlation_label(pearson_B, spearman_B)

# ============================================================
# 4. PANEL A: ACTIVE YEARS VS ANOMALY RATE
# ============================================================

p1 <- ggplot(
  top1_observer_profile,
  aes(
    x = active_years,
    y = either_anomaly_rate
  )
) +
  geom_point(
    alpha = 0.45,
    size = 2,
    colour = "#2166AC"
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    colour = "#B2182B",
    fill = "#EF8A62",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_A,
    hjust = 1.05,
    vjust = 1.3,
    size = 4.5,
    family = "Times New Roman",
    fontface = "bold",
    lineheight = 1.3
  ) +
  scale_y_continuous(
    limits = c(0, 0.30),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "(A) Observer experience",
    x = "Active years",
    y = "Either-anomaly rate"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 11,
      colour = "black"
    )
  )

# ============================================================
# 5. PANEL B: TOTAL CHECKLISTS VS ANOMALY RATE
# ============================================================

p2 <- ggplot(
  top1_observer_profile,
  aes(
    x = total_checklists,
    y = either_anomaly_rate
  )
) +
  geom_point(
    alpha = 0.45,
    size = 2,
    colour = "#2166AC"
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    colour = "#B2182B",
    fill = "#EF8A62",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_B,
    hjust = 1.05,
    vjust = 1.3,
    size = 4.5,
    family = "Times New Roman",
    fontface = "bold",
    lineheight = 1.3
  ) +
  scale_x_log10(
    breaks = c(1000, 3000, 10000, 30000),
    labels = comma
  ) +
  scale_y_continuous(
    limits = c(0, 0.30),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "(B) Observer contribution",
    x = "Total checklists submitted (log scale)",
    y = "Either-anomaly rate"
  ) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 11,
      colour = "black"
    )
  )

# ============================================================
# 6. COMBINE PANELS
# ============================================================

p_scatter <- p1 + p2 +
  plot_layout(ncol = 2)

p_scatter

# ============================================================
# 7. SAVE FIGURE
# ============================================================

ggsave(
  "FILE PATH/Observer_experience_anomaly_scatter.jpg",
  plot = p_scatter,
  width = 12,
  height = 6,
  dpi = 400
)

# ============================================================
# 8. PRINT CORRELATION RESULTS
# ============================================================

cat("\nPANEL A: ACTIVE YEARS VS ANOMALY RATE\n")
print(pearson_A)
print(spearman_A)

cat("\nPANEL B: TOTAL CHECKLISTS VS ANOMALY RATE\n")
print(pearson_B)
print(spearman_B)

# ============================================================
# BETA-BINOMIAL REGRESSION RESULTS
# ============================================================

library(glmmTMB)

# ------------------------------------------------------------
# 1. FIT BETA-BINOMIAL MODEL
# ------------------------------------------------------------

model_bb <- glmmTMB(
  cbind(
    either_anomalous_checklists,
    total_checklists - either_anomalous_checklists
  ) ~ active_years + log_checklists,
  family = betabinomial(link = "logit"),
  data = glm_data
)

# ------------------------------------------------------------
# 2. EXTRACT CONDITIONAL MODEL COEFFICIENTS
# ------------------------------------------------------------

coef_bb <- summary(model_bb)$coefficients$cond

# ------------------------------------------------------------
# 3. EXTRACT 95% CONFIDENCE INTERVALS
# ------------------------------------------------------------

ci_bb <- confint(
  model_bb,
  parm = "beta_",
  method = "Wald"
)

# ------------------------------------------------------------
# 4. COMBINE RESULTS
# ------------------------------------------------------------

bb_results <- data.frame(
  term = rownames(coef_bb),
  estimate = coef_bb[, "Estimate"],
  std_error = coef_bb[, "Std. Error"],
  z_value = coef_bb[, "z value"],
  p_value = coef_bb[, "Pr(>|z|)"],
  OR = exp(coef_bb[, "Estimate"]),
  OR_2.5 = exp(ci_bb[, "2.5 %"]),
  OR_97.5 = exp(ci_bb[, "97.5 %"]),
  row.names = NULL
)

# ------------------------------------------------------------
# 5. VIEW RESULTS
# ------------------------------------------------------------

print(bb_results)

# ------------------------------------------------------------
# 6. SAVE RESULTS AS CSV
# ------------------------------------------------------------

write_csv(
  bb_results,
  "FILE PATH/beta_binomial_model_results.csv"
)

################### TOP 1% CONTRIBUTOR CHECKLISTS ANOMALY ANALYSES ###########

library(DescTools)

# ============================================================
# 1. PREPARE TOP 1% CHECKLIST DATA
# ============================================================

cas_test_data <- top1_checklists %>%
  filter(
    !is.na(checklist_id),
    !is.na(observer_id),
    !is.na(CAS_class)
  ) %>%
  mutate(
    CAS_class = str_to_lower(
      str_trim(as.character(CAS_class))
    ),
    CAS_class = recode(
      CAS_class,
      "very low" = "Very Low",
      "low" = "Low",
      "moderate" = "Moderate",
      "high" = "High",
      "very high" = "Very High"
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
    )
  ) %>%
  filter(!is.na(CAS_class)) %>%
  distinct(checklist_id, .keep_all = TRUE) %>%
  mutate(
    richness_flag = as.integer(richness_anomaly_flag),
    abundance_flag = as.integer(abundance_anomaly_flag),
    either_flag = as.integer(either_anomaly_flag)
  )

# Check dataset

cat("\nTotal top 1% checklists:", nrow(cas_test_data), "\n")

cat(
  "Unique observers:",
  n_distinct(cas_test_data$observer_id),
  "\n"
)

print(
  table(cas_test_data$CAS_class)
)

# ============================================================
# 2. DEFINE ANOMALY TYPES
# ============================================================

flag_variables <- c(
  "Richness" = "richness_flag",
  "Abundance" = "abundance_flag",
  "Either anomaly" = "either_flag"
)

# ============================================================
# 3. CALCULATE ANOMALY PREVALENCE
#    WITH EXACT BINOMIAL 95% CONFIDENCE INTERVALS
# ============================================================

prevalence_summary <- bind_rows(
  lapply(names(flag_variables), function(type) {
    
    flag_col <- flag_variables[[type]]
    
    cas_test_data %>%
      group_by(CAS_class) %>%
      summarise(
        total_checklists = sum(!is.na(.data[[flag_col]])),
        flagged_checklists = sum(.data[[flag_col]] == 1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        anomaly_type = type,
        prevalence = 100 * flagged_checklists / total_checklists
      ) %>%
      rowwise() %>%
      mutate(
        lower_CI = 100 * binom.test(
          flagged_checklists,
          total_checklists
        )$conf.int[1],
        upper_CI = 100 * binom.test(
          flagged_checklists,
          total_checklists
        )$conf.int[2]
      ) %>%
      ungroup()
    
  })
)

prevalence_summary <- prevalence_summary %>%
  mutate(
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
    anomaly_type = factor(
      anomaly_type,
      levels = c(
        "Richness",
        "Abundance",
        "Either anomaly"
      )
    )
  )

print(prevalence_summary)

# ============================================================
# 4. PREVALENCE PLOT
# ============================================================

p_prevalence <- ggplot(
  prevalence_summary,
  aes(
    x = CAS_class,
    y = prevalence,
    colour = anomaly_type,
    group = anomaly_type
  )
) +
  
  geom_line(linewidth = 1) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(
      ymin = lower_CI,
      ymax = upper_CI
    ),
    width = 0.12,
    linewidth = 0.7
  ) +
  
  scale_colour_manual(
    values = c(
      "Richness" = "#2171B5",
      "Abundance" = "#D95F0E",
      "Either anomaly" = "#238B45"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = "CAS class",
    y = "Anomaly-flagged checklists (%)",
    colour = "Anomaly type"
  ) +
  
  theme_classic(
    base_family = "Times New Roman"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(size = 11),
    
    legend.position = "bottom"
  )

print(p_prevalence)

# ============================================================
# 5. COCHRAN-ARMITAGE LINEAR TREND TEST
# ============================================================
# Tests whether anomaly prevalence changes linearly
# across the five ordered CAS classes.
#
# The trend statistic is calculated from the 2 x 5
# contingency table using scores 1, 2, 3, 4, 5.
# ============================================================

cochran_armitage <- function(flag, group) {
  
  tab <- table(
    factor(flag, levels = c(0, 1)),
    factor(
      group,
      levels = c(
        "Very Low",
        "Low",
        "Moderate",
        "High",
        "Very High"
      )
    )
  )
  
  n <- colSums(tab)
  events <- tab["1", ]
  
  scores <- 1:5
  
  N <- sum(n)
  A <- sum(events)
  non_events <- N - A
  
  if (N <= 1 || A == 0 || non_events == 0) {
    return(
      data.frame(
        Chi_square = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  mean_score <- sum(n * scores) / N
  
  numerator <- sum(
    scores * (events - n * A / N)
  )
  
  denominator <- sqrt(
    (A * non_events / (N * (N - 1))) *
      sum(n * (scores - mean_score)^2)
  )
  
  Z <- numerator / denominator
  
  data.frame(
    Chi_square = Z^2,
    Z_statistic = Z,
    p_value = pchisq(
      Z^2,
      df = 1,
      lower.tail = FALSE
    )
  )
}


cochran_results <- bind_rows(
  lapply(names(flag_variables), function(type) {
    
    flag_col <- flag_variables[[type]]
    
    result <- cochran_armitage(
      cas_test_data[[flag_col]],
      cas_test_data$CAS_class
    )
    
    data.frame(
      Anomaly_type = type,
      result
    )
  })
)

print(cochran_results)

# ============================================================
# CHI-SQUARE TEST AND CRAMER'S V
# ============================================================

# Define anomaly variables
flag_variables <- c(
  "Richness anomaly" = "richness_flag",
  "Abundance anomaly" = "abundance_flag",
  "Either anomaly" = "either_flag"
)

# Check required objects and columns
stopifnot(
  exists("cas_test_data"),
  all(c("CAS_class", unname(flag_variables)) %in%
        names(cas_test_data))
)

# Store results
chi_test_objects <- list()
contingency_tables <- list()

# Run analyses
chi_results <- bind_rows(
  lapply(names(flag_variables), function(type) {
    
    flag_col <- flag_variables[[type]]
    
    # Remove missing observations for this particular test
    test_data <- cas_test_data %>%
      filter(
        !is.na(CAS_class),
        !is.na(.data[[flag_col]])
      )
    
    # Construct contingency table
    tab <- table(
      CAS_class = test_data$CAS_class,
      Anomaly = factor(
        test_data[[flag_col]],
        levels = c(0, 1)
      )
    )
    
    # Chi-square test
    chi <- chisq.test(
      tab,
      correct = FALSE
    )
    
    # Cramer's V
    cramers_v <- sqrt(
      as.numeric(chi$statistic) /
        (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1))
    )
    
    # Store individual results
    chi_test_objects[[type]] <<- chi
    contingency_tables[[type]] <<- tab
    
    # Return summary
    data.frame(
      Anomaly_type = type,
      N = sum(tab),
      Chi_square = as.numeric(chi$statistic),
      df = as.numeric(chi$parameter),
      p_value = chi$p.value,
      Cramers_V = cramers_v
    )
  })
)

# Display contingency tables
contingency_tables

# Display individual chi-square tests
chi_test_objects

# Display combined results
print(chi_results)

# ============================================================
# 7. SAVE RESULTS
# ============================================================

output_dir <- "FILE PATH"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

write.csv(
  prevalence_summary,
  file.path(
    output_dir,
    "Top1_CAS_anomaly_prevalence.csv"
  ),
  row.names = FALSE
)

write.csv(
  cochran_results,
  file.path(
    output_dir,
    "Top1_Cochran_Armitage_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  chi_results,
  file.path(
    output_dir,
    "Top1_Chi_square_CramersV_results.csv"
  ),
  row.names = FALSE
)

ggsave(
  file.path(
    output_dir,
    "Top1_CAS_anomaly_prevalence.jpg"
  ),
  plot = p_prevalence,
  width = 10,
  height = 7,
  dpi = 400
)

cat("\nAll three analyses completed and saved.\n")

# Total records
nrow(cas_test_data)

# Unique checklist IDs
dplyr::n_distinct(cas_test_data$checklist_id)

# Number of duplicated checklist IDs
sum(duplicated(cas_test_data$checklist_id))

# Check whether a checklist ID is associated with multiple observers
cas_test_data %>%
  group_by(checklist_id) %>%
  summarise(
    n_observers = n_distinct(observer_id),
    n_CAS_classes = n_distinct(CAS_class),
    .groups = "drop"
  ) %>%
  filter(n_observers > 1 | n_CAS_classes > 1)

# Export top 1% contributor checklists
write.csv(
  top1_checklists,
  "FILE PATH/top1_checklists.csv",
  row.names = FALSE,
  na = ""
)
