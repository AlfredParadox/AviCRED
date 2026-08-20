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

loadfonts(device = "win")
windowsFonts(
  Times = windowsFont("Times New Roman")
)
setwd("file path")

ebd_file <- "ebd raw.txt"

con <- dbConnect(
  duckdb(),
  dbdir = "avicred.duckdb"
)

################ YEAR-WISE CHECKLISTS GRAPH ############

dbExecute(con, "
CREATE OR REPLACE VIEW ebd_file AS
SELECT *
FROM read_csv_auto(
'ebd raw.txt',
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

############### SPLIT 2000-09 AND 2010-25 ################

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_2000_09 AS
SELECT *
FROM ebd
WHERE YEAR(\"OBSERVATION DATE\") BETWEEN 2000 AND 2009
")

dbExecute(con, "
CREATE OR REPLACE TABLE ebd_2010_25 AS
SELECT *
FROM ebd
WHERE YEAR(\"OBSERVATION DATE\") BETWEEN 2010 AND 2025
")

dbGetQuery(con, "
  SELECT
    '2000_2009' AS period,
    COUNT(DISTINCT \"SAMPLING EVENT IDENTIFIER\") AS n_checklists
  FROM ebd_2000_09

  UNION ALL

  SELECT
    '2010_2025' AS period,
    COUNT(DISTINCT \"SAMPLING EVENT IDENTIFIER\") AS n_checklists
  FROM ebd_2010_25
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
TO 'ebd_filtered.txt'
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
  'ebd_filtered.txt',
  delim = '\t',
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
TO 'checklist_master.txt'
(
  FORMAT CSV,
  DELIMITER '\t',
  HEADER TRUE
)
")

dbGetQuery(con, "
SELECT column_name
FROM information_schema.columns
WHERE table_name = 'checklist_master'
")

############ LOAD FILE IF SAVED PREVIOUSLY ############

dbExecute(con, "
CREATE OR REPLACE TABLE checklist_master AS
SELECT *
FROM read_csv_auto(
  'checklist_master.txt',
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

#================ ORI METRICS ================

############ CALCULATE OBSERVER METRICS #########

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
FROM ebd_raw
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

########### CREATE FINAL TABLE #############

dbExecute(con, "
CREATE OR REPLACE TABLE observer_master AS
SELECT
    m.observer_id,
    m.years_active,
    m.n_birding_days,
    m.n_checklists,
    s.total_species
FROM observer_metrics m
LEFT JOIN observer_species s
ON m.observer_id = s.observer_id
")

observer_master <- dbGetQuery(con, "
SELECT *
FROM observer_master
")

ori <- observer_master[, c(
  "years_active",
  "n_birding_days",
  "n_checklists",
  "total_species"
)]

cor(ori, use = "complete.obs")

########### STANDARDIZE PARAMETERS ############

ori_scores <- data.frame(
  years_active   = scale(observer_master$years_active),
  n_birding_days = scale(log1p(observer_master$n_birding_days)),
  n_checklists   = scale(log1p(observer_master$n_checklists)),
  total_species  = scale(log1p(observer_master$total_species))
)

observer_master$ORI <- rowMeans(ori_scores)

sapply(ori_scores, max, na.rm = TRUE)

DBI::dbWriteTable(   # WRITE UPDATED TABLE BACK TO DUCKDB
  con,
  "observer_master",
  observer_master,
  overwrite = TRUE
)

write.csv(
  observer_master,
  "file path/ori_master.csv",
  row.names = FALSE
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

#=============== RICHNESS AND ABUNDANCE CAS SCORES =============

################## MAKE STATE X SEASON GROUPINGS #############

dbExecute(con, "
CREATE OR REPLACE TABLE checklist_master_season AS
SELECT
    *,
    CASE
        WHEN EXTRACT(MONTH FROM obs_date) IN (12,1,2) THEN 'Winter'
        WHEN EXTRACT(MONTH FROM obs_date) IN (3,4,5) THEN 'Summer'
        WHEN EXTRACT(MONTH FROM obs_date) IN (6,7,8,9) THEN 'Monsoon'
        WHEN EXTRACT(MONTH FROM obs_date) IN (10,11) THEN 'Post-monsoon'
    END AS season
FROM checklist_master
")

################### CREATE CONTEXTUAL AUTHENTICITY SCORES ############

dbExecute(con, "
CREATE OR REPLACE TABLE contextual_stats AS
WITH med AS (
SELECT
state,
season,
median(richness) AS median_richness,
median(abundance_per_checklist) AS median_abundance
FROM checklist_master_season

GROUP BY
state,
season
)

SELECT
m.state,
m.season,
m.median_richness,
m.median_abundance,
median(abs(c.richness-m.median_richness)) AS mad_richness,
median(abs(c.abundance_per_checklist-m.median_abundance))
AS mad_abundance

FROM med m
JOIN checklist_master_season c
ON m.state=c.state
AND m.season=c.season

GROUP BY

m.state,
m.season,
m.median_richness,
m.median_abundance
")

############## CALCULATE FLAGS, RCHNESS/ABUNDANCE SCORES AND ROBUST Z-SCORE ############

dbExecute(con, "
CREATE OR REPLACE TABLE checklist_authenticity AS
SELECT
c.*,
s.median_richness,
s.median_abundance,
s.mad_richness,
s.mad_abundance,

CASE
WHEN s.mad_richness=0 THEN NULL
ELSE
0.6745*(c.richness-s.median_richness)/s.mad_richness
END AS richness_robust_z,

CASE
WHEN s.mad_abundance=0 THEN NULL
ELSE
0.6745*(c.abundance_per_checklist-s.median_abundance)/s.mad_abundance
END AS abundance_robust_z,

CASE
WHEN s.mad_richness=0 THEN 0
WHEN
0.6745*(c.richness-s.median_richness)/s.mad_richness>3.5
THEN 1
ELSE 0
END AS richness_flag,

CASE
WHEN s.mad_abundance=0 THEN 0
WHEN
0.6745*(c.abundance_per_checklist-s.median_abundance)/s.mad_abundance>3.5
THEN 1
ELSE 0
END AS abundance_flag,

CASE
    WHEN richness_robust_z <= 0 THEN 1.0
    ELSE EXP(-richness_robust_z / 3.5)
END AS richness_score,

CASE
    WHEN abundance_robust_z <= 0 THEN 1.0
    ELSE EXP(-abundance_robust_z / 3.5)
END AS abundance_score,

(
    (
        CASE
            WHEN richness_robust_z <= 0 THEN 1.0
            ELSE EXP(-richness_robust_z / 3.5)
        END
    ) +
    (
        CASE
            WHEN abundance_robust_z <= 0 THEN 1.0
            ELSE EXP(-abundance_robust_z / 3.5)
        END
    )
) / 2 AS CAS

FROM checklist_master_season c
JOIN contextual_stats s
ON
c.state=s.state
AND
c.season=s.season
")

############## CALCULATE MEAN CAS PER OBSERVER, ANMALY RATES, FLAG COUNTS ###############

dbExecute(con, "
CREATE OR REPLACE TABLE observer_authenticity AS
SELECT

    observer_id,

    -- Number of submitted checklists
    COUNT(*) AS n_checklists,

    -- Number of flagged checklists
    SUM(richness_flag) AS richness_flag_count,
    SUM(abundance_flag) AS abundance_flag_count,

    -- Number flagged for both richness and abundance
    SUM(
        CASE
            WHEN richness_flag = 1
             AND abundance_flag = 1
            THEN 1
            ELSE 0
        END
    ) AS both_flag_count,

    -- Proportion (anomaly rate)
    AVG(richness_flag) AS richness_anomaly_rate,
    AVG(abundance_flag) AS abundance_anomaly_rate,

    AVG(
        CASE
            WHEN richness_flag = 1
             AND abundance_flag = 1
            THEN 1
            ELSE 0
        END
    ) AS both_anomaly_rate,

    -- Mean checklist scores
    AVG(richness_score) AS richness_score,
    AVG(abundance_score) AS abundance_score,

    -- Observer Checklist Authenticity Score
    AVG(CAS) AS mean_CAS

FROM checklist_authenticity

GROUP BY observer_id
")

checklist_authenticity <- dbGetQuery(con, "
SELECT *
FROM checklist_authenticity
")

observer_authenticity <- dbGetQuery(con, "
SELECT *
FROM observer_authenticity
")

write.csv(
  observer_authenticity,
  "file path/Observer mean CAS.csv",
  row.names = FALSE
)

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

#==================  CHECKLIST CONSISTENCY =================#

############# SUCCESSIVE CALCULATE CHECKLIST GAPS ##########

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

summary(regularity$consistency_score)

#================ CHECKLIST SPATIAL COVERAGE SCORE (SCS) ================#

############ CALCULATE SHANNON ENTROPY, EFFECTIVE LOCALITIES, SCS SCORE #########

dbExecute(con, "
CREATE OR REPLACE TABLE spatial_coverage AS
WITH loc_counts AS (

SELECT
    observer_id,
    locality_id,
    COUNT(*) AS n_checklists
FROM checklist_master
GROUP BY observer_id, locality_id

),

proportions AS (
SELECT
    observer_id,
    locality_id,
    n_checklists,
    n_checklists * 1.0 /
    SUM(n_checklists) OVER (PARTITION BY observer_id) AS p_i
FROM loc_counts

),

entropy AS (
SELECT

observer_id,
-SUM(p_i * LN(p_i)) AS shannon_spatial_coverage,
COUNT(DISTINCT locality_id) AS n_localities,
EXP(-SUM(p_i * LN(p_i))) AS effective_localities

FROM proportions
GROUP BY observer_id
)

SELECT

observer_id,
shannon_spatial_coverage,
n_localities,
effective_localities,

CASE
WHEN n_localities <= 1 THEN 0
ELSE shannon_spatial_coverage / LN(n_localities)
END AS spatial_coverage_score

FROM entropy
")

#=============== CHECKLIST ATTRIBUTE INDEX (CAI) MASTER TABLE =================#

########## JOIN ALL PARAMETERS TO MAKE MASTER CAI TABLE #########

dbExecute(con, "
CREATE OR REPLACE TABLE cai_master AS

SELECT

    os.observer_id,

    --------------------------------------------------
    -- Checklist Authenticity Score (CAS)
    --------------------------------------------------

    os.richness_anomaly_rate,
    os.abundance_anomaly_rate,
    os.both_anomaly_rate,
    os.richness_score,
    os.abundance_score,
    os.mean_CAS,
    os.richness_flag_count,
    os.abundance_flag_count,
    os.both_flag_count,

    --------------------------------------------------
    -- Checklist completeness
    --------------------------------------------------

    cc.total_checklists,
    cc.complete_checklists,
    cc.incomplete_checklists,
    cc.completeness_score,
    cc.weighted_completeness_score,

    --------------------------------------------------
    -- Checklist consistency (Robust)
    --------------------------------------------------

    kr.median_gap_days,
    kr.mad_gap_days,
    kr.robust_cv,
    kr.consistency_score,

    --------------------------------------------------
    -- Spatial coverage
    --------------------------------------------------

    sc.n_localities,
    sc.shannon_spatial_coverage,
    sc.spatial_coverage_score,
    sc.effective_localities

FROM observer_authenticity AS os

LEFT JOIN completeness_checklist AS cc
ON os.observer_id = cc.observer_id

LEFT JOIN checklist_regularity AS kr
ON os.observer_id = kr.observer_id

LEFT JOIN spatial_coverage AS sc
ON os.observer_id = sc.observer_id
")

cai_master <- dbGetQuery(con, "
SELECT *
FROM cai_master
")

############# STANDARDIZE PARAMETERS ############

scores <- data.frame(
  CAS          = cai_master$mean_CAS,
  completeness = cai_master$weighted_completeness_score,
  consistency  = cai_master$consistency_score,
  spatial      = cai_master$spatial_coverage_score
)

cai_master$CAI <- rowMeans(
  scores,
  na.rm = TRUE
)

cai_master$CAI[is.nan(cai_master$CAI)] <- NA

DBI::dbWriteTable(  # WIRTE UPDATED TABLE BACK TO DUCKDB
  con,
  "cai_master",
  cai_master,
  overwrite = TRUE
)

View(cai_master)

write.csv(
  cai_master,
  "file path/cai_master.csv",
  row.names = FALSE
)

############### OBSERVER MASTER TABLE AND CHECKLIST MASTER TABLE JOIN ##########

ori_master <- fread("ori_master.csv")
cai_master <- fread("cai_master.csv")

ori_master$observer_id <- as.character(observer_master$observer_id)
cai_master$observer_id <- as.character(cai_master$observer_id)

dbWriteTable(
  con,
  "observer_master",
  ori_master,
  overwrite = TRUE
)

dbWriteTable(
  con,
  "cai_master",
  cai_master,
  overwrite = TRUE
)

#=============== AVICRED ==================

dbExecute(con, "
CREATE OR REPLACE TABLE observer_final AS
SELECT

    o.observer_id,

    --------------------------------------------------
    -- ORI
    --------------------------------------------------

    o.years_active,
    o.n_birding_days,
    o.n_checklists,
    o.total_species,
    o.ORI,

    --------------------------------------------------
    -- CAI Components
    --------------------------------------------------

    c.mean_CAS,
    c.richness_score,
    c.abundance_score,
    c.weighted_completeness_score,
    c.consistency_score,
    c.spatial_coverage_score,
    c.CAI,

    --------------------------------------------------
    -- Descriptive stats
    --------------------------------------------------

    c.n_localities,
    c.shannon_spatial_coverage,
    c.effective_localities,
    c.richness_anomaly_rate,
    c.abundance_anomaly_rate,
    c.both_anomaly_rate,
    c.abundance_flag_count,
    c.richness_flag_count,
    c.both_flag_count,

FROM observer_master o

LEFT JOIN cai_master c
ON TRIM(o.observer_id) = TRIM(c.observer_id)
")

observer_final <- dbGetQuery(con, "
SELECT *
FROM observer_final
")

########### CALCULATE AVICRED RELIABILITY SCORE ########

observer_final$AviCRED <-
  (observer_final$ORI + observer_final$CAI) / 2

DBI::dbWriteTable(   #WRITE UPDATED TABLE BACK TO DUCKDB
  con,
  "observer_final",
  observer_final,
  overwrite = TRUE
)

############ SCORE STATS ###########

dbGetQuery(con, "
SELECT

    MIN(ORI) AS min_ori,
    quantile_cont(ORI, 0.5) AS median_ori,
    MAX(ORI) AS max_ori,

    MIN(CAI) AS min_cai,
    quantile_cont(CAI, 0.5) AS median_cai,
    MAX(CAI) AS max_cai,

    MIN(AviCRED) AS min_avicred,
    quantile_cont(AviCRED, 0.5) AS median_avicred,
    MAX(AviCRED) AS max_avicred

FROM observer_final
")

# ==================== ISOLATION FOREST ========================#

################## CHECKLIST ###############

checklists <- dbGetQuery(con, "
SELECT
    checklist_id,
    state,
    season,
    richness,
    abundance_per_checklist,
    CAS
FROM checklist_authenticity
")

checklists <- checklists %>%   #CHECK FOR MISSING VALUES AND OMIT
  filter(
    !is.na(richness),
    !is.na(abundance_per_checklist),
    !is.na(CAS)
  )

checklists$log_abundance <-
  log1p(checklists$abundance_per_checklist)  #LOG TRANSFORM ABUNDANCE TO REDUCE POSITVE SKEWNESS

checklists$group <-   #CREATE STATE X SEASON GROUPS FOR ISOLATION FOREST
  paste(
    checklists$state,
    checklists$season,
    sep = "_"
  )

checklists$iforest_score <- NA

groups <- unique(checklists$group)

for (g in groups) {
  idx <- which(checklists$group == g)

  if (length(idx) < 100) {
    next
  }

  features <- data.frame(
    richness = checklists$richness[idx],
    log_abundance = checklists$log_abundance[idx]
  )

  model <- isolation.forest(
    features,
    ntrees = 300,  #SET NUMBER OF TREES BASED ON DATA COMPLEXITY, HETEROGENITY AND NUMBER OF VARIABLES
    sample_size = min(256, nrow(features)),
    ndim = 1,
    nthreads = 0
  )

  checklists$iforest_score[idx] <-
    predict(
      model,
      features,
      type = "score"
    )
}

checklists$if_authenticity <-  #TRANFORM IF ANOMALY SCORE TO IF AUTHENTICITY SCORE FOR CHECKLISTS
  1 - checklists$iforest_score

summary(
  checklists[, c(
    "richness",
    "abundance_per_checklist"
  )]
)

summary(checklists$if_authenticity)

cor(                  #CHECK CHECKLIST CAS AND IF AUTHENTICITY CORRELATION
  checklists$CAS,
  checklists$if_authenticity,
  method = "pearson",
  use = "complete.obs"
)

######## RETRIEVE CHECKLIS FLAGS AND CREATE CAS CLASSES ########

flags <- dbGetQuery(con, "  
SELECT
    checklist_id,
    richness_flag,
    abundance_flag
FROM checklist_anomalies
")

checklists <- merge(
  checklists,
  flags,
  by = "checklist_id",
  all.x = TRUE
)

checklists$CAS_class <- cut(    #CREATE CAS CLASS FIXED THRESHOLDS
  checklists$CAS,
  breaks = c(
    -Inf,
    0.25,
    0.50,
    0.75,
    0.90,
    Inf
  ),
  labels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  include.lowest = TRUE
)

aggregate(
  if_authenticity ~ CAS_class,
  checklists,
  mean
)

checklists$if_5pct <- ifelse(  # TOP 5% IF AUTHENTICITY SCORES
  checklists$if_authenticity >=
    quantile(
      checklists$if_authenticity,
      0.95,
      na.rm = TRUE
    ),
  1, 0
)

checklists$if_10pct <- ifelse(  # TOP 10% IF AUTHENTICITY SCORES
  checklists$if_authenticity >=
    quantile(
      checklists$if_authenticity,
      0.90,
      na.rm = TRUE
    ),
  1, 0
)

################# OBSERVER ####################

feature_cols <- c(
  
  # ---------- ORI components ----------
  "years_active",
  "n_birding_days",
  "n_checklists",
  "total_species",

  # ---------- CAI components ----------
  "mean_CAS",
  "weighted_completeness_score",
  "consistency_score",
  "spatial_coverage_score"
)

complete_rows <- complete.cases(  #REMOVE NA AND MISSING VALUES
  observer_final[, feature_cols]
)

ml_complete <- observer_final[complete_rows, ]

features <- scale(    #STANDARDIZE FEATURES PRIOR TO ISOLATION FOREST AS OBSERVER FEATURES ARE HIGHLY HETEROGENOUS
  ml_complete[, feature_cols]
)

iso_model <- isolation.forest(
  features,
  ntrees = 500,  #SET NUMBER OF TREES BASED ON DATA COMPLEXITY, HETEROGENITY AND NUMBER OF VARIABLES
  sample_size = 256
)

ml_complete$ml_score <- predict(  #CALCULATE OBSERVER ISOLATION FOREST SCORE OR ML SCORE
  iso_model,
  features,
  type = "score"
)

avicred_quintiles <- quantile(  #CREATE QUINTILE BASED AVICRED CLASSES
  ml_complete$AviCRED,
  probs = c(0, 0.20, 0.40, 0.60, 0.80, 1),
  na.rm = TRUE
)

ml_complete$avicred_class <- cut(
  ml_complete$AviCRED,
  breaks = avicred_quintiles,
  include.lowest = TRUE,
  labels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  )
)

table(ml_complete$AviCRED)

cutoff_5 <- quantile(  #TOP 5% IF ANOMALY SCORERS
  ml_complete$ml_score,
  probs = 0.95,
  na.rm = TRUE
)

ml_complete$ml_5pct <- as.integer(
  ml_complete$ml_score >= cutoff_5
)

cutoff_10 <- quantile( #TOP 10% IF ANOMALY SCORERS
  ml_complete$ml_score,
  probs = 0.90,
  na.rm = TRUE
)

ml_complete$ml_10pct <- as.integer(
  ml_complete$ml_score >= cutoff_10
)

# ======================== K-MEANS CLUSTERING ======================#

############## ELBOW METHOD FOR OPTIMAL CLUSTER DETERMINATION IN CHECKLISTS ANALYSIS ##############

chk_features <- scale(         #STANDARDIZE PRIOR TO ELBOW METHOD
  checklists[, c("richness", "log_abundance")]
)

set.seed(123)

wss <- numeric(10)

for (k in 1:10) {  
  km <- kmeans(
    chk_features,
    centers = k,
    nstart = 25,
    iter.max = 100
  )
  wss[k] <- km$tot.withinss
}

plot(    #DETERMINE OPTIMAL CLUSTER FROM ELBOW PLOT (2 IN THIS CASE)
  1:10,
  wss,
  type = "b",
  pch = 19,
  xlab = "Number of clusters (k)",
  ylab = "Total Within-Cluster Sum of Squares"
)

groups <- unique(checklists$group)

###### RUN ACTUAL K-MEANS WITH OPTIMAL CLUSTER COUNT ######

for (g in groups) {
  idx <- which(checklists$group == g)

  if (length(idx) < 100) {
    next
  }

  group_features <- scale(        #STANDARDIZE PRIOR TO ACTUAL CLUSTERING
    checklists[idx, c("richness", "log_abundance")]
  )
  km <- kmeans(
    group_features,
    centers = 2,  # OPTIMAL CLUSTERS DETERMINED TO BE 2 FROM ELBOW PLOT
    nstart = 50
  )
  checklists$cluster[idx] <- km$cluster
}

write.csv(
  checklists,
  "file path/Checklists FINAL.csv",
  row.names = FALSE
)

########### ELBOW METHOD FOR OPTIMAL CLUSTER DETERMINATION FOR OBSERVER ANALYSIS ###############

obs_features <- ml_complete %>%
  dplyr::select(n_checklists, years_active, n_birding_days, total_species, mean_CAS, weighted_completeness_score, consistency_score, spatial_coverage_score) %>%
  na.omit()

obs_features <- scale(obs_features) #SCALE FEATURES PRIOR TO ELBOW METHOD 

set.seed(123)

wss <- numeric(10)
for (k in 1:10) {
  km <- kmeans(
    obs_features,
    centers = k,
    nstart = 25,
    iter.max = 100
  )
  wss[k] <- km$tot.withinss
}

plot(   #DETERMINE OPTIMAL CLUSTER BY ELBOW PLOT (3 IN THIS CASE)
  1:10,
  wss,
  type = "b",
  pch = 19,
  xlab = "Number of clusters (k)",
  ylab = "Total Within-Cluster Sum of Squares"
)

####### RUN ACTUAL K-MEANS WITH OPTIMAL CLUSTER COUNT #####

features_scaled <- scale(   #SCALE FEATURES PRIOR TO ACTUAL CLUSTERING
  ml_complete[, c(
    "years_active",
    "n_birding_days",
    "n_checklists",
    "total_species",
    "mean_CAS",
    "weighted_completeness_score",
    "consistency_score",
    "spatial_coverage_score"
  )]
)

km <- kmeans(
  features_scaled,
  centers = 3,  #OPTIMAL CLUSTER COUNT DETERMINED TO BE 3 FROM ELBOW PLOT
  nstart = 50,
  iter.max = 100
)

ml_complete$cluster <- km$cluster

table(
  ml_complete$cluster,
  ml_complete$avicred_class
)

table(
  ml_complete$cluster,
  ml_complete$ml_5pct
)

write.csv(
  ml_complete,
  "file path/AviCRED FINAL.csv",
  row.names = FALSE
)

#=================== ORI vs CAI CORRELATION AND DISCORDANCE  ==============#

######### CREATE ORI-CAI GROUPS FOR DISCORDANCE ANALYSIS #############

discordance_df <- ml_complete %>%
  filter(
    !is.na(ORI),
    !is.na(CAI),
    !is.na(AviCRED),
    !is.na(avicred_class)
  ) %>%
  mutate(
    ORI_group = case_when(
      ORI >= quantile(ORI, 0.75) ~ "High ORI",
      ORI <= quantile(ORI, 0.25) ~ "Low ORI",
      TRUE ~ "Medium ORI"
    ),
    CAI_group = case_when(
      CAI >= quantile(CAI, 0.75) ~ "High CAI",
      CAI <= quantile(CAI, 0.25) ~ "Low CAI",
      TRUE ~ "Medium CAI"
    ),
    discordance_group = case_when(
      ORI_group == "High ORI" &
        CAI_group == "Low CAI" ~ "High ORI - Low CAI",
      
      ORI_group == "Low ORI" &
        CAI_group == "High CAI" ~ "Low ORI - High CAI",
      
      TRUE ~ "Other"
    )
  )

pearson_test <- cor.test(   # Pearson correlation
  discordance_df$ORI,
  discordance_df$CAI,
  method = "pearson"
)

spearman_test <- cor.test(  # Spearman correlation
  discordance_df$ORI,
  discordance_df$CAI,
  method = "spearman"
)

pearson_r <- round(pearson_test$estimate, 3)
spearman_r <- round(spearman_test$estimate, 3)

############## ORI-CAI CORRELATION ################

ggplot(
  discordance_df,
  aes(
    x = ORI,
    y = CAI,
    colour = AviCRED
  )
)+
  geom_point(
    alpha = 0.5,
    size = 1.5
  )+
  scale_colour_viridis_c(
    name = "AviCRED"
  )+
  geom_smooth(
    method = "lm",
    colour = "black",
    linewidth = 1,
    se = FALSE
  )+
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.5,
    label = paste0(
      "Pearson r = ", round(pearson_r, 3),
      "\nSpearman \u03c1 = ", round(spearman_r, 3)
    ),
    family = "Times New Roman",
    size = 4
  ) +
  
  labs(
    title = "ORI–CAI Relationship with AviCRED Distribution",
    x = "Observer Reliability Index (ORI)",
    y = "Checklist Attribute Index (CAI)",
    colour = "AviCRED"
  )+
  
  theme_classic()+
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
    axis.title = element_text(
      face = "bold"
    )
  )

############## ORI-CAI DISCORDANCE BOX PLOT ##############

discordance_box <- discordance_df %>%
  mutate(
    group = case_when(
      ORI_group == "High ORI" &
        CAI_group == "Low CAI" ~ "High ORI - Low CAI",
      ORI_group == "Low ORI" &
        CAI_group == "High CAI" ~ "Low ORI - High CAI",
      TRUE ~ "Other"
    )
  ) %>%
  filter(group != "Other")

box_labels <- discordance_box %>%  #CALCULATE MIN, MEDIAN AND MAX
  group_by(group) %>%
  summarise(
    Min = min(AviCRED),
    Q2 = median(AviCRED),
    Max = max(AviCRED)
  ) %>%
  pivot_longer(
    cols = c(Min, Q2, Max),
    names_to = "label",
    values_to = "value"
  )

ggplot(
  discordance_box,
  aes(
    x = group,
    y = AviCRED,
    fill = group
  )
)+
  geom_boxplot(
    width = 0.6,
    outlier.alpha = 0.3
  )+
  geom_text(
    data = box_labels,
    aes(
      x = group,
      y = value,
      label = round(value, 3)
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.5,
    family = "Times New Roman",
    vjust = -0.5
  )+
  labs(
    title = "AviCRED Scores Across ORI–CAI Discordant Observer Groups",
    x = "Observer Group",
    y = "AviCRED Score"
  )+
  theme_classic()+
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
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "none"
  )

#================= ISOLATION FOREST VS AVICRED/CAS CORRELATION =================#

############# CHECKLIST LEVEL ##############

pearson <- cor.test(
  checklists$CAS,
  checklists$if_authenticity,
  method = "pearson"
)

spearman <- cor.test(
  checklists$CAS,
  checklists$if_authenticity,
  method = "spearman"
)

# Plot
ggplot(
  checklists,
  aes(
    x = CAS,
    y = if_authenticity
  )
) +
  
  geom_point(
    alpha = 0.08,
    size = 0.4,
    colour = "steelblue"
  ) +
  
  geom_smooth(
    method = "lm",
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic(base_family = "Times New Roman") +
  
  labs(
    x = "Checklist Authenticity Score (CAS)",
    y = "Isolation Forest Score"
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.5,
    size = 5,
    family = "Times New Roman",
    label = paste0(
      "Pearson r = ",
      round(pearson$estimate, 3),
      "\n",
      "Spearman \u03C1 = ",
      round(spearman$estimate, 3),
      "\n"
    )
  )

############ OBSERVER LEVEL #################

stats_df <- data.frame(
  Metric = c("AviCRED", "ORI", "CAI"),
  Pearson = c(
    cor(ml_complete$ml_score,
        ml_complete$AviCRED,
        method = "pearson",
        use = "complete.obs"),
    
    cor(ml_complete$ml_score,
        ml_complete$ORI,
        method = "pearson",
        use = "complete.obs"),
    
    cor(ml_complete$ml_score,
        ml_complete$CAI,
        method = "pearson",
        use = "complete.obs")
  )
)

stats_df$Spearman <- c(
  cor(ml_complete$ml_score,
      ml_complete$AviCRED,
      method = "spearman",
      use = "complete.obs"),
  
  cor(ml_complete$ml_score,
      ml_complete$ORI,
      method = "spearman",
      use = "complete.obs"),
  
  cor(ml_complete$ml_score,
      ml_complete$CAI,
      method = "spearman",
      use = "complete.obs")
)

stats_df$R2 <- stats_df$Pearson^2

stats_df

stats_df$label <- paste0(
  "Pearson r = ", round(stats_df$Pearson, 3),
  "\nSpearman ρ = ", round(stats_df$Spearman, 3)
)

corr_df <- ml_complete %>%
  filter(!is.na(ml_score)) %>%
  select(
    ml_score,
    AviCRED,
    ORI,
    CAI
  ) %>%
  pivot_longer(
    cols = c(AviCRED, ORI, CAI),
    names_to = "Metric",
    values_to = "Value"
  )

ggplot(
  corr_df,
  aes(
    x = Value,
    y = ml_score
  )
) +
  geom_point(
    alpha = 0.2,
    size = 1
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
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
    hjust = -0.1,
    vjust = 1.1,
    family = "Times New Roman",
    size = 4
  ) +
  
  labs(
    title = "Relationship between Isolation Forest Score and AviCRED Metrics",
    x = "Metric Value",
    y = "Isolation Forest Score"
  ) +
  
  theme_classic(base_family = "Times New Roman") +
  
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
    axis.text = element_text(
      size = 12
    ),
    strip.text = element_text(
      size = 12,
      face = "bold"
    )
  )

# ================== ISOLATION FOREST/CLUSTERS VS AVICRED STACKED ================= #

################ CHECKLIST ML ###############

checklists <- checklists %>%  #REMOVE NA's AND MISSING VALUES 
  filter(
    !is.na(CAS_class),
    !is.na(cluster),
    !is.na(if_10pct),
    !is.na(if_5pct)
  )

checklists$CAS_class <- factor(
  checklists$CAS_class,
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  ordered = TRUE
)

checklists$IF_Class <- case_when( #CREATE IF CLASSES FROM TOP 5% AND 10% CLASSIFICATIONS
  checklists$if_5pct == 1 ~ "Top 5%",
  checklists$if_10pct == 1 ~ "5–10%",
  TRUE ~ "Remaining 90%"
)

checklists$IF_Class <- factor(
  checklists$IF_Class,
  levels = c(
    "Top 5%",
    "5–10%",
    "Remaining 90%"
  )
)

if_df <- checklists %>% #SEGREGATE CAS CLASSES AS PER IF AUTHENTICITY CLASSIFICATIONS
  count(
    IF_Class,
    CAS_class
  ) %>%
  group_by(IF_Class) %>%
  mutate(
    Proportion = n / sum(n),
    Panel = "Isolation Forest"
  ) %>%
  ungroup() %>%
  rename(
    Group = IF_Class
  )

if_totals <- checklists %>%
  count(IF_Class) %>%
  mutate(
    Panel = "Isolation Forest"
  ) %>%
  rename(
    Group = IF_Class,
    total = n
  )

cluster_df <- checklists %>%  #SEGREGATE CAS CLASSES AS PER CLUSTER CLASSIFICATIONS
  mutate(
    cluster_label = factor(
      cluster,
      levels = c(1, 2),
      labels = c(
        "Cluster 1",
        "Cluster 2"
      )
    )
  ) %>%
  count(
    cluster_label,
    CAS_class
  ) %>%
  group_by(cluster_label) %>%
  mutate(
    Proportion = n / sum(n),
    Panel = "Checklist Clusters"
  ) %>%
  ungroup() %>%
  rename(
    Group = cluster_label
  )

cluster_totals <- checklists %>%
  mutate(
    cluster_label = factor(
      cluster,
      levels = c(1, 2),
      labels = c(
        "Cluster 1",
        "Cluster 2"
      )
    )
  ) %>%
  count(cluster_label) %>%
  mutate(
    Panel = "Checklist Clusters"
  ) %>%
  rename(
    Group = cluster_label,
    total = n
  )

plot_df <- bind_rows(
  if_df,
  cluster_df
)

total_df <- bind_rows(
  if_totals,
  cluster_totals
)

ggplot(
  plot_df,
  aes(
    x = Group,
    y = Proportion,
    fill = CAS_class
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill",
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(
      label = paste0(
        round(
          Proportion * 100,
          1
        ),
        "%"
      )
    ),
    position = position_fill(
      vjust = 0.5
    ),
    family = "Times New Roman",
    size = 3,
    colour = "black"
  ) +
  geom_text(
    data = total_df,
    aes(
      x = Group,
      y = 1.04,
      label = comma(total)
    ),
    inherit.aes = FALSE,
    family = "Times New Roman",
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(
    ~Panel,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "Very Low"  = "#B2182B",
      "Low"       = "#EF8A62",
      "Moderate"  = "#FDDBC7",
      "High"      = "#A6DBA0",
      "Very High" = "#1B7837"
    )
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(
      0,
      1.1
    )
  ) +
  labs(
    title = "Distribution of CAS Classes across Checklist Metrics",
    x = NULL,
    y = "Percentage of Checklists",
    fill = "CAS Class"
  ) +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

############# OBSERVER ML ###############

ml_complete <- ml_complete %>% #REMOVE NA's AND MISSING VALUES
  filter(
    !is.na(avicred_class),
    !is.na(cluster),
    !is.na(ml_10pct),
    !is.na(ml_5pct)
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

ml_complete$IF_Class <- case_when(  #CREATE IF CLASSES FROM TOP 5% AND 10% CLASSIFICATIONS
  ml_complete$ml_5pct == 1 ~ "Top 5%",
  ml_complete$ml_10pct == 1 ~ "5–10%",
  TRUE ~ "Remaining 90%"
)

ml_complete$IF_Class <- factor(
  ml_complete$IF_Class,
  levels = c(
    "Top 5%",
    "5–10%",
    "Remaining 90%"
  )
)

if_df_obs <- ml_complete %>%    #SEGREGATE AVICRED CLASSES AS PER IF ANOMALY CLASSIFICATIONS
  count(IF_Class, avicred_class) %>%
  group_by(IF_Class) %>%
  mutate(
    Proportion = n / sum(n),
    Panel = "Isolation Forest"
  ) %>%
  ungroup()

if_totals_obs <- ml_complete %>%
  count(IF_Class) %>%
  mutate(
    Panel = "Isolation Forest"
  )

cluster_df_obs <- ml_complete %>%  #SEGREGATE AVICRED CLASSES AS PER CLUSTER CLASSIFICATIONS
  mutate(
    cluster_label = factor(
      cluster,
      levels = c(1, 2, 3),
      labels = c(
        "Cluster 1",
        "Cluster 2",
        "Cluster 3"
      )
    )
  ) %>%
  count(cluster_label, avicred_class) %>%
  group_by(cluster_label) %>%
  mutate(
    Proportion = n / sum(n),
    Panel = "Clusters"
  ) %>%
  ungroup() %>%
  rename(Group = cluster_label)

cluster_totals_obs <- ml_complete %>%
  mutate(
    cluster_label = factor(
      cluster,
      levels = c(1, 2, 3),
      labels = c(
        "Cluster 1",
        "Cluster 2",
        "Cluster 3"
      )
    )
  ) %>%
  count(cluster_label) %>%
  mutate(
    Panel = "Clusters"
  ) %>%
  rename(
    Group = cluster_label,
    total = n
  )

if_df_obs <- if_df_obs %>%
  rename(
    Group = IF_Class
  )

if_totals_obs <- if_totals_obs %>%
  rename(
    Group = IF_Class,
    total = n
  )

if_df_obs$Group <- factor(
  if_df_obs$Group,
  levels = c(
    "Top 5%",
    "5–10%",
    "Remaining 90%"
  )
)

if_totals_obs$Group <- factor(
  if_totals_obs$Group,
  levels = c(
    "Top 5%",
    "5–10%",
    "Remaining 90%"
  )
)

plot_df_obs <- bind_rows(
  if_df_obs,
  cluster_df_obs
)

total_df_obs <- bind_rows(
  if_totals_obs,
  cluster_totals_obs
)

ggplot(
  plot_df_obs,
  aes(
    x = Group,
    y = Proportion,
    fill = avicred_class
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill",
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(
      label = paste0(
        round(Proportion * 100, 1),
        "%"
      )
    ),
    position = position_fill(vjust = 0.5),
    family = "Times New Roman",
    size = 3,
    colour = "black"
  ) +
  geom_text(
    data = total_df_obs,
    aes(
      x = Group,
      y = 1.04,
      label = scales::comma(total)
    ),
    inherit.aes = FALSE,
    family = "Times New Roman",
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(
    ~Panel,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "Very Low"  = "#FDD0A2",
      "Low"       = "#FC9272",
      "Moderate"  = "#FB6A4A",
      "High"      = "#807DBA",
      "Very High" = "#41B6C4"
    )
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1.1)
  ) +
  labs(
    title = "Distribution of AviCRED Classes across Observer Metrics",
    x = NULL,
    y = "Percentage of Observers",
    fill = "AviCRED Class"
  ) +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

#=============== RELATIVE IMPORTANCE OF CLUSTER ===============#

########## CHECKLIST LEVEL ##############

checklists_plot <- checklists %>% #REMOVE MISSING VALUES
  filter(
    !is.na(cluster),
    !is.na(if_authenticity)
  )

checklists_plot <- checklists_plot %>%
  group_by(cluster) %>%
  summarise(
    Richness = mean(richness, na.rm = TRUE),
    Abundance = mean(abundance_per_checklist, na.rm = TRUE),
    Log_Abundance = mean(log_abundance, na.rm = TRUE),
    .groups = "drop"
  )

checklists_plot

profile_long <- checklists_plot %>%
  pivot_longer(
    cols = -cluster,
    names_to = "Variable",
    values_to = "Value"
  )

ggplot(
  profile_long,
  aes(
    x = reorder(Variable, Value),
    y = Value,
    fill = factor(cluster)
  )
) +
  geom_col(
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    aes(
      label = round(Value, 2)
    ),
    position = position_dodge(width = 0.9),
    hjust = ifelse(profile_long$Value > 0, -0.2, 1.2),
    family = "Times New Roman",
    size = 3.5
  ) +
  coord_flip() +
  labs(
    title = "Checklist Cluster Characteristics",
    x = "",
    y = "Values",
    fill = "Cluster"
  ) +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

############## OBSERVER LEVEL ###################

cluster_profile <- ml_complete %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>%
  summarise(
    # ORI
    Years_Active = mean(years_active, na.rm = TRUE),
    Birding_Days = mean(n_birding_days, na.rm = TRUE),
    Checklists = mean(n_checklists, na.rm = TRUE),
    Total_Species = mean(total_species, na.rm = TRUE),

    # CAI
    Mean_CAS = mean(mean_CAS, na.rm = TRUE),
    Completeness = mean(weighted_completeness_score, na.rm = TRUE),
    Consistency = mean(consistency_score, na.rm = TRUE),
    Spatial_Coverage = mean(spatial_coverage_score, na.rm = TRUE),
    .groups = "drop"
  )

cluster_profile

profile_z <- cluster_profile

profile_z[, -1] <- scale( #STANDARDIZE BEFORE RELATIVE IMPORTANCE ANALYSIS
  profile_z[, -1]
)

profile_z

profile_long <- profile_z %>%
  pivot_longer(
    cols = -cluster,
    names_to = "Variable",
    values_to = "ZScore"
  )

ggplot(
  profile_long,
  aes(
    x = reorder(Variable, ZScore),
    y = ZScore,
    fill = factor(cluster)
  )
) +
  geom_col(
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    aes(
      label = round(ZScore, 2)
    ),
    position = position_dodge(width = 0.9),
    hjust = ifelse(profile_long$ZScore > 0, -0.2, 1.2),
    family = "Times New Roman",
    size = 3.5
  ) +
  coord_flip() +
  labs(
    title = "Relative Importance of Variables Across Clusters",
    x = "",
    y = "Standardized Mean (Z-score)",
    fill = "Cluster"
  ) +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

# ================= UPSET PLOT ===============#

############# checklist #####################

upset_data1 <- data.frame(
  cluster1_flag = checklists$cluster == 1,
  cluster2_flag = checklists$cluster == 2,
  if_top5_flag = checklists$if_5pct == 1,
  CAS_high_flag = checklists$CAS_class == "Very High"
)

colnames(upset_data1) <- c(
  "K-Means Cluster 1",
  "K-Means Cluster 2",
  "Top 5% Isolation Score",
  "CAS Class Very High"
)

upset_data1 <- upset_data1[
  rowSums(upset_data1) > 0,
]

ComplexUpset::upset(
  upset_data1,
  c(
    "K-Means Cluster 1",
    "K-Means Cluster 2",
    "Top 5% Isolation Score",
    "CAS Class Very High"
  )
) +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

############# OBSERVER ################

upset_data2 <- data.frame(
  cluster1_flag      = ml_complete$cluster == 1,
  cluster2_flag      = ml_complete$cluster == 2,
  cluster3_flag      = ml_complete$cluster == 3,
  ml_top5_flag       = ml_complete$ml_5pct == 1,
  avicred_high_flag  = ml_complete$avicred_class == "Very High"
)

colnames(upset_data2) <- c(
  "K-Means Cluster 1",
  "K-Means Cluster 2",
  "K-Means Cluster 3",
  "Top 5% Isolation Score",
  "AviCRED Class Very High"
)

upset_data2 <- upset_data2[rowSums(upset_data2) > 0, ]

upset_data2 <- subset(
  upset_data2,
  `K-Means Cluster 1` |
    `K-Means Cluster 2` |
    `K-Means Cluster 3` |
    `Top 5% Isolation Score` |
    `AviCRED Class Very High`
)

ComplexUpset::upset(
  upset_data2,
  c(
    "K-Means Cluster 1",
    "K-Means Cluster 2",
    "K-Means Cluster 3",
    "Top 5% Isolation Score",
    "AviCRED Class Very High"
  )
) +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Times New Roman",
      size = 12
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

################## EFFECTIVE LOCALITIES vs MEAN CAS/AVICRED CORRELATION ############

ml_complete %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>%
  summarise(
    Median_n_localities = median(n_localities, na.rm = TRUE),
    Mean_n_localities = mean(n_localities, na.rm = TRUE),
    Median_effective = median(effective_localities, na.rm = TRUE),
    Mean_effective = mean(effective_localities, na.rm = TRUE)
  )

p1_pearson <- cor.test(
  ml_complete$effective_localities,
  ml_complete$mean_CAS,
  method = "pearson"
)

p1_spearman <- cor.test(
  ml_complete$effective_localities,
  ml_complete$mean_CAS,
  method = "spearman"
)

p2_pearson <- cor.test(
  ml_complete$effective_localities,
  ml_complete$AviCRED,
  method = "pearson"
)

p2_spearman <- cor.test(             
  ml_complete$effective_localities,
  ml_complete$AviCRED,
  method = "spearman"
)

plot1 <- ggplot(
  ml_complete,
  aes(
    x = effective_localities,
    y = mean_CAS
  )
) +
  geom_point(alpha = 0.4) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    title = "(A) Mean CAS",
    x = "Effective Localities (e^H)",
    y = "Mean CAS",
    subtitle = paste0(
      "r = ", round(p1_pearson$estimate, 3),
      " | \u03C1 = ", round(p1_spearman$estimate, 3)
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(
      family = "Times New Roman"
    )
  )

plot2 <- ggplot(
  ml_complete,
  aes(
    effective_localities,
    AviCRED
  )
) +
  geom_point(alpha = .3) +
  scale_x_log10() +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    title = "(B) AviCRED Score",
    x = "Effective Localities (e^H)",
    y = "AviCRED Score",
    subtitle = paste0(
      "r = ", round(p2_pearson$estimate, 3),
      " | \u03C1 = ", round(p2_spearman$estimate, 3)
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(
      family = "Times New Roman"
    )
  )

(plot1 + plot2)

#=============== ALLUVIAL PLOT ===============#

############# CHECKLIST ###############

checklists$ml_class <- "Remaining 90%"

checklists$ml_class[
  checklists$if_authenticity >=
    quantile(
      checklists$if_authenticity,
      0.90,
      na.rm = TRUE
    )
] <- "Top 10%"

checklists$ml_class[
  checklists$if_authenticity >=
    quantile(
      checklists$if_authenticity,
      0.95,
      na.rm = TRUE
    )
] <- "Top 5%"

checklists$ml_class <- factor(
  checklists$ml_class,
  levels = c(
    "Top 5%",
    "Top 10%",
    "Remaining 90%"
  )
)

checklists$CAS_class <- factor(
  checklists$CAS_class,
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  ),
  ordered = TRUE
)

plot_data_chk <- checklists %>%
  filter(
    !is.na(cluster),
    !is.na(CAS_class),
    !is.na(ml_class)
  ) %>%
  mutate(
    cluster = factor(cluster)
  ) %>%
  group_by(
    cluster,
    CAS_class,
    ml_class
  ) %>%
  summarise(
    Freq = n(),
    .groups = "drop"
  )

ggplot(
  plot_data_chk,
  aes(
    axis1 = cluster,
    axis2 = CAS_class,
    axis3 = ml_class,
    y = Freq
  )
) +
  geom_alluvium(
    aes(fill = cluster),
    width = 0.6,
    alpha = 1,
    curve_type = "xspline",
    segments = 100
  ) +
  geom_stratum(
    width = 0.6,
    fill = "grey90",
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "1" = "#E3B505",
      "2" = "#610345"
    )
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    family = "Times New Roman",
    size = 4
  ) +
  scale_x_discrete(
    limits = c(
      "Cluster",
      "CAS",
      "Isolation Forest"
    ),
    expand = c(0.05, 0.05)
  ) +
  labs(
    title = "Checklist Classification Flow Across Methods",
    x = NULL,
    y = "Number of Checklists",
    fill = "Cluster"
  ) +
  theme_minimal(
    base_family = "Times New Roman"
  ) +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

############# OBSERVER ##############

plot_data <- ml_complete %>%
  mutate(
    cluster = factor(
      cluster,
      levels = c(2, 1, 3)
    ),
    AviCRED_class = factor(
      avicred_class,
      levels = c(
        "Very Low",
        "Low",
        "Moderate",
        "High",
        "Very High"
      )
    ),
    ml_class = case_when(
      ml_5pct == 1 ~ "Top 5%",
      ml_10pct == 1 ~ "Top 10%",
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
  group_by(
    cluster,
    AviCRED_class,
    ml_class
  ) %>%
  summarise(
    Freq = n(),
    .groups = "drop"
  ) %>%
  arrange(cluster)

ggplot(
  plot_data,
  aes(
    axis1 = cluster,
    axis2 = AviCRED_class,
    axis3 = ml_class,
    y = Freq
  )
) +
  scale_fill_manual(values = c(
    "2" = "#68A691",
    "1" = "#004643",
    "3" = "#E6AF2E"
  )) +
  geom_alluvium(
    aes(fill = cluster),
    width = 0.6,
    alpha = 1,
    curve_type = "xspline",
    segments = 100
  ) +
  geom_stratum(
    width = 0.6,
    fill = "grey90",
    colour = "black"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    family = "Times New Roman",
    size = 4
  ) +
  scale_x_discrete(
    limits = c("Cluster", "AviCRED", "Isolation Forest"),
    expand = c(0.05, 0.05)
  ) +
  labs(
    title = "Observer Classification Flow Across Methods",
    x = NULL,
    y = "Number of Observers",
    fill = "Cluster"
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

###################### OBSERVER MAP #####################

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
    'ebd_filtered.txt',
    delim='\t'
)
"
)

india_shp <- st_read(   #INSERT RESPECTIVE SHAPEFILE AFTER DOWNLOADING
  "file path/India_State_Boundary"
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

################## CLUSTER MAP ##############

ggplot() +
  geom_sf(
    data = cluster_map,
    aes(fill = factor(cluster)),
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
      "1" = "#DEEFB7",
      "2" = "#4F772D",
      "3" = "#2B193D"
    ),
    na.translate = FALSE,
    name = "K-Means Cluster"
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title = "K-Means Cluster Distribution"
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
  count(  #DOMINANT IF AUTHENTICITY CATEGORY DISPLAY PER GRID
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

################ CLUSTER MAP ##############

chk_cluster <- gis_join_chk %>%
  st_drop_geometry() %>%
  distinct(      #AVOID COUNTING DUPLICATE CHECKLISTS PER GRID
    grid_id,
    checklist_id,
    cluster
  ) %>%
  count(  #DOMINANT CLUSTER CATEGORY DISPLAY PER GRID
    grid_id,
    cluster,
    name = "n"
  ) %>%
  group_by(grid_id) %>%
  slice_max(
    n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    cluster = factor(
      cluster,
      levels = c(1, 2)
    )
  )

cluster_map <- left_join(
  grid25,
  chk_cluster,
  by = "grid_id"
)

ggplot() +
  geom_sf(
    data = cluster_map,
    aes(fill = factor(cluster)),
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
      "2" = "#307351",
      "1" = "#E3B505"
    ),
    breaks = c("2", "1"),
    na.translate = FALSE,
    name = "K-Means Cluster"
  ) +
  coord_sf(expand = FALSE) +
  labs(title = "K-Means Cluster Checklist Distribution") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

############## FLAGGED RECORDS DISTRIBUTION ACROSS CHECKLISTS AND OBSERVER METRICS #############

df <- ml_complete  # Create working copies of observer and checklist-level datasets
chk_df <- checklists

checklist_obs <- dbGetQuery(  # Get the observer associated with each checklist
  con,
  "
  SELECT
      checklist_id,
      observer_id
  FROM checklist_master
  "
)

chk_df <- chk_df %>%  # Add observer IDs to the checklist dataset
  left_join(
    checklist_obs,
    by = "checklist_id"
  )

chk_df <- chk_df %>%  # Keep variables required for the anomaly analysis
  select(
    checklist_id,
    observer_id,
    CAS_class,
    cluster,
    if_5pct,
    if_10pct,
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

obs_cluster <- df %>%  # Calculate observer anomaly proportions within each cluster
  filter(
    !is.na(cluster),
    !is.na(ml_5pct),
    !is.na(ml_10pct)
  ) %>%
  count(cluster, observer_anomaly_type) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

obs_if <- df %>%  # Classify observers according to Isolation Forest thresholds
  filter(
    !is.na(cluster),
    !is.na(ml_5pct),
    !is.na(ml_10pct)
  ) %>%
  mutate(
    IF = case_when(
      ml_5pct == 1 ~ "Top 5%",
      ml_10pct == 1 ~ "Top 10%",
      TRUE ~ "Remaining 90%"
    ),
    IF = factor(
      IF,
      levels = c(
        "Top 5%",
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

checklists_cluster <- chk_df %>%  # Calculate checklist anomaly proportions within each cluster
  filter(!is.na(cluster))

chk_cluster <- checklists_cluster %>%
  count(cluster, anomaly_type) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

checklists_if <- chk_df %>%  # Prepare checklist data for Isolation Forest classification
  filter(
    !is.na(if_5pct),
    !is.na(if_10pct)
  )

chk_if <- checklists_if %>%   # Calculate checklist anomaly proportions within Isolation Forest groups
  mutate(
    IF = case_when(
      if_5pct == 1 ~ "Top 5%",
      if_10pct == 1 ~ "Top 10%",
      TRUE ~ "Remaining 90%"
    ),
    IF = factor(
      IF,
      levels = c(
        "Top 5%",
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
    family = "Times New Roman",
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
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(size = 10)
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
    family = "Times New Roman",
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
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(size = 10)
  )

p3 <- ggplot(
  obs_cluster,
  aes(
    x = factor(cluster),
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
    family = "Times New Roman",
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
    title = "Observer Cluster",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(size = 10)
  )

p4 <- ggplot(
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
    family = "Times New Roman",
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
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10)
  )

p5 <- ggplot(
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
    family = "Times New Roman",
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
    title = "Checklist IF Authenticity Classes",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(
    base_family = "Times New Roman"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.position = "bottom",
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      size = 10
    )
  )

p6 <- ggplot(
  chk_cluster,
  aes(
    x = factor(cluster),
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
    family = "Times New Roman",
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
    title = "Checklist Cluster",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(
    base_family = "Times New Roman"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.position = "bottom",
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      size = 10
    )
  )

(p1 | p2 | p3) /
  (p4 | p5 | p6) 

############## SAMPLING EFFORT HEXBIN MAP ###############

checklist_coords <- dbGetQuery(
  con,
  "
SELECT DISTINCT
    \"SAMPLING EVENT IDENTIFIER\" AS checklist_id,
    \"LATITUDE\" AS latitude,
    \"LONGITUDE\" AS longitude
FROM ebd_filtered
"
)

checklist_coords <- checklist_coords %>%
  filter(
    !is.na(latitude),
    !is.na(longitude)
  )

checklist_sf <- st_as_sf(   
  checklist_coords,
  coords = c("longitude","latitude"),
  crs = 4326
)

checklist_sf <- st_transform( #CONVERT CHECKLISTS LAT LONG TO EPSG:3857 from EPSG:4326
  checklist_sf,
  3857
)

coords <- st_coordinates(checklist_sf)

checklist_sf$x <- coords[,1]
checklist_sf$y <- coords[,2]

ggplot() +
  geom_sf(
    data = india,
    fill = "grey95",
    colour = "black",
    linewidth = 0.35
  ) +
  coord_sf(expand = FALSE) +
  geom_hex(
    data = checklist_sf,
    aes(
      x = x,
      y = y
    ),
    bins = 200  #SET HIGHER NUMBER OF BINS TO ENHANCE SPATIAL RESOLUTION
  ) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "Number of\nchecklists"
  ) +
  coord_sf() +
  labs(
    title = "Spatial distribution of eBird sampling effort across India",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      size = 10
    ),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank()
  )

##################### MEAN CAS VS AVICRED CORRELATION ##############

pearson <- cor.test(
  df$mean_CAS,
  df$AviCRED,
  method = "pearson"
)


spearman <- cor.test(
  df$mean_CAS,
  df$AviCRED,
  method = "spearman"
)

pearson_r <- round(pearson$estimate, 3)
pearson_p <- signif(pearson$p.value, 3)

spearman_rho <- round(spearman$estimate, 3)
spearman_p <- signif(spearman$p.value, 3)

ggplot(df, aes(x = mean_CAS, y = AviCRED)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Relationship between Mean CAS and AVICRED Score",
    subtitle = paste0(
      "Pearson r = ", pearson_r,
      " (p = ", pearson_p, ")   |   ",
      "Spearman \u03C1 = ", spearman_rho,
      " (p = ", spearman_p, ")"
    ),
    x = "Mean CAS",
    y = "AviCRED Score"
  ) +
  theme_bw()

################### CHECKLISTS SUBMITTED PER AVICRED CLASS ############

# Count the number of checklists associated with each AviCRED class

treemap_class <- dbGetQuery(
  con,
  "
  SELECT
      a.avicred_class,
      COUNT(cm.checklist_id) AS n_checklists

  FROM chk_ids AS ch

  LEFT JOIN checklist_master AS cm
      ON ch.checklist_id = cm.checklist_id

  LEFT JOIN avicred_tbl AS a
      ON cm.observer_id = a.observer_id

  WHERE a.avicred_class IS NOT NULL

  GROUP BY a.avicred_class
  "
)

treemap_class$avicred_class <- factor(  # Arrange AviCRED classes in the desired order
  treemap_class$avicred_class,
  levels = c(
    "Very Low",
    "Low",
    "Moderate",
    "High",
    "Very High"
  )
)

ggplot(
  treemap_class,
  aes(
    area = n_checklists,
    fill = avicred_class,
    label = paste0(
      avicred_class,
      "\n",
      comma(n_checklists)
    )
  )
) +
  geom_treemap(
    colour = "white",
    linewidth = 1
  ) +
  geom_treemap_text(
    place = "centre",
    reflow = TRUE,
    colour = "black",
    fontface = "bold",
    family = "Times New Roman",
    size = 16,
    min.size = 0
  ) +
  scale_fill_manual(
    values = c(
      "#E8F6EF",
      "#ABEBC6",
      "#52BE80",
      "#239B56",
      "#145A32"
    )
  ) +
  labs(
    title = "Checklist Contributions Across AviCRED Classes",
    fill = "AviCRED Class"
  ) +
  theme_void() +
  theme(
    text = element_text(
      family = "Times New Roman"
    ),
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(
      size = 12
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

############### NUMBER OF ANOMALOUS CHECKLISTS PER AVICRED CLASS ##################

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
      family = "Times New Roman",
      size = 12
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      colour = "black",
      size = 12
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 12
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
    family = "Times New Roman",
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
    family = "Times New Roman",
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
    family = "Times New Roman",
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
    title = "Checklist Contributions Across Observer Experience and AviCRED Classes",
    x = "Years Active",
    y = "Total Checklists Submitted",
    fill = "AviCRED Class"
  ) +
  
  theme_bw(base_family = "Times New Roman") +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
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

cas_avicred <- checklists %>%  # Combine CAS classifications with the corresponding observer's AviCRED class
  select(
    checklist_id,
    CAS_class
  ) %>%
  left_join(
    checklist_observer,
    by = "checklist_id"
  ) %>%
  left_join(
    df %>%
      select(
        observer_id,
        avicred_class
      ),
    by = "observer_id"
  ) %>%
  filter(
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
    family = "Times New Roman",
    size = 3.5
  ) +
  
  scale_fill_viridis_c(
    option = "C",
    trans = "log10",
    labels = scales::comma,
    name = "No. of checklists"
  ) +
  
  labs(
    title = "Checklist Counts Across CAS and AviCRED Classes",
    x = "AviCRED Class",
    y = "CAS Class"
  ) +
  
  theme_minimal() +
  
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
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank()
  )

###################### DISTRIBUTION OF ORI/CAI ACROSS AVICRED CLASSES #################

plot_df <- ml_complete %>%   #REMOVE NAs
  filter(
    !is.na(ORI),
    !is.na(CAI),
    !is.na(avicred_class)
  )

plot_df$avicred_class <- factor(  #ORDER AVICRED CLASSES
  plot_df$avicred_class,
  levels = c("Very Low", "Low", "Moderate", "High", "Very High")
)

medians1 <- plot_df %>%   #CALCULATE CAI MEDIANS
  group_by(avicred_class) %>%
  summarise(
    median_CAI = median(CAI),
    .groups = "drop"
  )

medians2 <- plot_df %>%   #CALCULATE ORI MEDIANS
  group_by(avicred_class) %>%
  summarise(
    median_ORI = median(ORI),
    .groups = "drop"
  )

CAI_ridge <- ggplot(
  plot_df,
  aes(
    x = CAI,
    y = avicred_class,
    fill = avicred_class
  )
)+
  scale_fill_manual(
    values = c(
      "Very Low"  = "#FDD0A2",
      "Low"       = "#FC9272",
      "Moderate"  = "#FB6A4A",
      "High"      = "#807DBA",
      "Very High" = "#41B6C4"
    )
  ) +
  
  geom_density_ridges(
    alpha = 0.85,
    scale = 1.2,
    rel_min_height = 0.001,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_vline(        #MEDIAN LINES
    data = medians1,
    aes(
      xintercept = median_CAI
    ),
    linetype = "dashed",
    linewidth = 0.3,
    inherit.aes = FALSE
  )+
  
  geom_text(
    data = medians1,
    aes(
      x = median_CAI,
      y = avicred_class,
      label = round(median_CAI, 2)
    ),
    inherit.aes = FALSE,
    vjust = -0.8,
    size = 3.5,
    family = "Times New Roman",
    fontface = "bold"
  )+
  
  labs(
    title = "Distribution of CAI Across AviCRED Classes",
    x = "Checklist Attribute Index (CAI)",
    y = "AviCRED Class"
  )+
  
  theme_ridges()+
  
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
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.position = "none"
  )

ORI_ridge <- ggplot(
  plot_df,
  aes(
    x = ORI,
    y = avicred_class,
    fill = avicred_class
  )
)+   
  scale_fill_manual(
  values = c(
    "Very Low"  = "#FDD0A2",
    "Low"       = "#FC9272",
    "Moderate"  = "#FB6A4A",
    "High"      = "#807DBA",
    "Very High" = "#41B6C4"
  )
)  +
  
  geom_density_ridges(
    alpha = 0.85,
    scale = 1.2,
    rel_min_height = 0.001,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_vline(     #MEDIAN LINES
    data = medians2,
    aes(
      xintercept = median_ORI
    ),
    linetype = "dashed",
    linewidth = 0.3,
    inherit.aes = FALSE
  )+
  
  geom_text(
    data = medians2,
    aes(
      x = median_ORI,
      y = avicred_class,
      label = round(median_ORI, 2)
    ),
    inherit.aes = FALSE,
    vjust = -0.8,
    size = 3.5,
    family = "Times New Roman",
    fontface = "bold"
  )+
  
  labs(
    title = "Distribution of ORI Across AviCRED Classes",
    x = "Observer Reliability Index (ORI)",
    y = "AviCRED Class"
  )+
  
  theme_ridges()+
  
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
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.position = "none"
  )

(CAI_ridge + ORI_ridge)

###################### VARIANCE OF ORI VS CHECKLISTS SUBMITTED AND NUMBER OF OBSERVERS ###############

plot_df <- ml_complete %>% #REMOVE MISSING VALUES
  filter(
    !is.na(ORI),
    !is.na(n_checklists)
  ) %>%
  mutate(
    log_checklists = log10(n_checklists + 1)
  )

ggplot(
  plot_df,
  aes(
    x = ORI,
    y = log_checklists
  )
)+
  geom_hex(
    bins = 60
  )+
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    name = "Density of Observers"
  )+
  geom_smooth(
    method = "loess",
    span = 1,
    colour = "red",
    linewidth = 1.2,
    se = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      min(plot_df$ORI),
      max(plot_df$ORI)
    ),
    expand = c(0.02,0.02)
  )+
  scale_y_continuous(
    expand = c(0.02,0.02)
  )+
  labs(
    title = "Observer Reliability and Checklist Contributions",
    x = "Observer Reliability Index (ORI)",
    y = expression(Log[10]*"(Total Checklists + 1)")
  )+
  theme_classic()+
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
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    legend.title = element_text(
      face = "bold"
    )
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
    !is.na(if_10pct)         
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
  reference = factor(checklist_cm$if_10pct),
  positive = "1"
)

cm_chk$table       # Confusion matrix counts
cm_chk$overall     # Overall performance metrics
cm_chk$byClass     # Class-specific performance metrics

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

p1 <- ggplot(
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
    family = "Times New Roman",
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
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(
      size = 12
    ),
    axis.text.y = element_text(
      size = 12,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 11
    ),
    panel.grid = element_blank()
  )

p1

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

p2 <- ggplot(
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
    family = "Times New Roman",
    fontface = "bold",
    size = 4.5
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#F9C80E"  
  ) +
  labs(
    title = "Checklist-level",
    x = "IF Top 10%",
    y = "CAS High/Very High",
    fill = "Count"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(
      size = 12
    ),
    axis.text.y = element_text(
      size = 12,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 11
    ),
    panel.grid = element_blank()
  )

combined_plot <- (p1 + p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Agreement Between AviCRED/CAS and Isolation Forest Classifications",
    theme = theme(
      plot.title = element_text(
        family = "Times New Roman",
        face = "bold",
        size = 15,
        hjust = 0.5
      )
    )
  ) &
  theme(
    legend.position = "right"
  )

combined_plot

#============== MATTHEW'S CORRELATION COEFFICIENT BETWEEN IF AND AVICRED =================#

############# OBSERVER-LEVEL ##############

observer_cm <- ml_complete %>%
  filter(
    !is.na(avicred_class),
    !is.na(ml_5pct),
    !is.na(ml_10pct)
  ) %>%
  mutate(
    avicred_high = if_else(
      avicred_class %in% c("High", "Very High"),
      1,
      0
    ),
    
    if_top10 = if_else(
      ml_5pct == 1 | ml_10pct == 1,
      1,
      0
    )
  )

observer_mcc <- observer_cm %>%  #CALCULATE MCC
  mutate(
    avicred_high = factor(avicred_high),
    if_top10 = factor(if_top10)
  ) %>%
  mcc(
    truth = avicred_high,
    estimate = if_top10
  )

observer_results <- observer_cm %>%  #CHECK OVERLAP BETWEEN IF AND AVICRED
  summarise(
    High_Credibility = sum(avicred_high),
    IF_Top10 = sum(if_top10),
    Overlap = sum(
      avicred_high == 1 &
        if_top10 == 1
    )
  ) %>%
  mutate(
    `Overlap vs AviCRED (%)` =
      100 * Overlap / High_Credibility,
    
    `Overlap vs IF (%)` =
      100 * Overlap / IF_Top10,
    
    MCC = observer_mcc$.estimate
  )

observer_results

########### CHECKLIST-LEVEL ###########

checklist_cm <- checklists %>%
  filter(
    !is.na(CAS_class),
    !is.na(if_5pct),
    !is.na(if_10pct)
  ) %>%
  mutate(
    cas_high = if_else(
      CAS_class %in% c("High", "Very High"),
      1,
      0
    ),
    
    if_top10 = if_else(
      if_5pct == 1 | if_10pct == 1,
      1,
      0
    )
  )

checklist_mcc <- checklist_cm %>%  #CALCULATE MCC
  mutate(
    cas_high = factor(cas_high),
    if_top10 = factor(if_top10)
  ) %>%
  mcc(
    truth = cas_high,
    estimate = if_top10
  )

checklist_results <- checklist_cm %>%   #CHECK OVERLAP BEWTEEN CAS AND IF
  summarise(
    High_Credibility = sum(cas_high),
    IF_Top10 = sum(if_top10),
    Overlap = sum(
      cas_high == 1 &
        if_top10 == 1
    )
  ) %>%
  mutate(
    `Overlap vs CAS (%)` =
      100 * Overlap / High_Credibility,
    
    `Overlap vs IF (%)` =
      100 * Overlap / IF_Top10,
    
    MCC = checklist_mcc$.estimate
  )

checklist_results

########### SUMMARY TABLE ###########

summary_table <- tibble(
  
  Comparison = c(
    "Checklist (CAS vs IF)",
    "Observer (AviCRED vs IF)"
  ),
  
  `High Credibility Total` = c(
    checklist_results$High_Credibility,
    observer_results$High_Credibility
  ),
  
  `Combined IF Top 10%` = c(
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

summary_table

write.csv(
  summary_table,
  "MCC Summary.csv",
  row.names = FALSE
)

#============ AGREEMENT BETWEEN ANOMALY FLAGGED OBSERVERS ACROSS IF AND AVICRED =========#

# Get checklist-observer relationships from the database
checklist_obs <- dbGetQuery(
  con,
  "
  SELECT
      checklist_id,
      observer_id
  FROM checklist_master
  "
)

checklists_obs <- checklists %>%   # Add observer IDs to the checklist-level dataset
  left_join(
    checklist_obs,
    by = "checklist_id"
  )

observer_anomaly <- checklists_obs %>%   # Flag observers who have at least one anomalous checklist
  group_by(observer_id) %>%
  summarise(
    anomaly_flag = if_else(
      any(richness_flag == 1 | abundance_flag == 1),  
      1,
      0
    ),
    .groups = "drop"
  )

observer_cm <- ml_complete %>%  # Combine anomaly, AviCRED and Isolation Forest information
  left_join(
    observer_anomaly,
    by = "observer_id"
  ) %>%
  mutate(
    anomaly_flag = replace_na(anomaly_flag, 0),  
    
    avicred_high = if_else(
      avicred_class %in% c("High", "Very High"),  
      1,
      0
    ),
    
    if_top10 = if_else(
      ml_5pct == 1 | ml_10pct == 1,  
      1,
      0
    )
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
    family = "Times New Roman",
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
    title = "Observer-level Confusion Matrices",
    x = "Observer Anomaly Status",
    y = "Predicted Classification",
    fill = "Count"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      size = 12
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.title = element_text(
      face = "bold"
    ),
    panel.grid = element_blank()
  )

################# HIGHEST/LOWEST OBSERVERS ACROSS AVICRED CATEGORY ############

observer_examples <- ml_complete %>%
  filter(
    !is.na(avicred_class),
    avicred_class %in% c("Very High", "Moderate", "Very Low")  #SELECT ANY AVICRED CATEGORY
  ) %>%
  group_by(avicred_class) %>%
  arrange(desc(AviCRED)) %>%
  slice(c(1:5, (n()-4):n())) %>%  #SELECT TOP 5 AND BOTTOM 5 OBSERVERS IN EACH CATEGORY
  ungroup() %>%
  select(
    observer_id,
    avicred_class,
    AviCRED,
    ORI,
    CAI,
    n_checklists,
    years_active,
    n_birding_days,
    total_species,
    mean_CAS,
    weighted_completeness_score,
    consistency_score,
    spatial_coverage_score
  )

write.csv(
  observer_examples,
  "AviCRED_representative_observers.csv",
  row.names = FALSE
)

############# ANOMALOUS RICHNESS/ABUNDANCE CHECKLISTS #############

top10_richness <- checklists %>%   #TOP 10 CHECKLISTS WITH HIGHEST RICHNESS
  arrange(desc(richness)) %>%
  slice_head(n = 10) %>%
  mutate(
    anomaly_type = "High richness"
  )

top10_abundance <- checklists %>%  #TOP 10 CHECKLISTS WITH HIGHEST ABUNDANCE
  arrange(desc(abundance_per_checklist)) %>%
  slice_head(n = 10) %>%
  mutate(
    anomaly_type = "High abundance"
  )

anomaly_checklists <- bind_rows(
  top10_richness,
  top10_abundance,
) %>%
  select(
    checklist_id,
    anomaly_type,
    state,
    season,
    richness,
    abundance_per_checklist,
    CAS
  )

anomaly_checklists

duckdb_register(
  con,
  "anomaly_checklists",
  anomaly_checklists
)

# RETRIEVE OBSERVATION RECORDS FOR EACH CHECKLISTS

species_anomaly_details <- dbGetQuery(
  con,
  '
  SELECT

      a.anomaly_type,
      a.checklist_id,
      a.state,
      a.season,
      a.richness,
      a.abundance_per_checklist,
      a.CAS,

      e."COMMON NAME" AS common_name,
      e."SCIENTIFIC NAME" AS scientific_name,

      TRY_CAST(
          e."OBSERVATION COUNT"
          AS DOUBLE
      ) AS species_abundance


  FROM anomaly_checklists a
  INNER JOIN ebd_filtered e
  ON a.checklist_id =
     e."SAMPLING EVENT IDENTIFIER"

  ORDER BY
      a.anomaly_type,
      a.checklist_id,
      species_abundance DESC
  '
)

species_anomaly_details %>%
  distinct(
    anomaly_type,
    checklist_id
  ) %>%
  count(anomaly_type)

head(species_anomaly_details)

write.csv(
  species_anomaly_details,
  "AviCRED_species_level_anomaly.csv",
  row.names = FALSE
)

################# ANOMALOUS CHECKLISTS PER OBSERVER CLASSIFICATION ############

checklist_anomalies <- checklists %>%  # Create checklist-level anomaly flag
  transmute(
    checklist_id,
    anomaly_flag = if_else(
      richness_flag == 1 | abundance_flag == 1,  # Flag richness or abundance anomalies
      1,
      0
    )
  )

checklist_anomalies <- checklist_anomalies %>%   # Link checklist anomalies to their respective observers
  left_join(
    checklist_obs,
    by = "checklist_id"
  )

observer_checklists <- ml_complete %>%  # Combine observer-level metrics with checklist anomaly information
  select(
    observer_id,
    avicred_class,
    cluster,
    ml_5pct,
    ml_10pct
  ) %>%
  left_join(
    checklist_anomalies,
    by = "observer_id"
  )

observer_checklists %>% # Summarise checklist anomalies by AviCRED class
  group_by(avicred_class) %>%
  summarise(
    Observers = n_distinct(observer_id),
    Total_Checklists = n(),
    Anomalous_Checklists = sum(anomaly_flag),
    Mean_Anomalous_Checklists_per_Observer =
      Anomalous_Checklists / Observers,
    .groups = "drop"
  )

observer_checklists %>% # Summarise checklist anomalies by cluster
  group_by(cluster) %>%
  summarise(
    Observers = n_distinct(observer_id),
    Total_Checklists = n(),
    Anomalous_Checklists = sum(anomaly_flag),
    Mean_Anomalous_Checklists_per_Observer =
      Anomalous_Checklists / Observers,
    .groups = "drop"
  )

observer_checklists %>%  # Classify observers according to their Isolation Forest thresholds
  mutate(
    IF_Group = case_when(
      ml_5pct == 1 ~ "Top 5%",       
      ml_10pct == 1 ~ "Top 10%",     
      TRUE ~ "Remaining"             
    )
  ) %>%
  group_by(IF_Group) %>%
  summarise(
    Observers = n_distinct(observer_id),
    Total_Checklists = n(),
    Anomalous_Checklists = sum(anomaly_flag),
    Mean_Anomalous_Checklists_per_Observer =
      Anomalous_Checklists / Observers,
    Percent_Anomalous_Checklists =
      100 * Anomalous_Checklists / Total_Checklists,
    .groups = "drop"
  )

################ NUMBER OF OBSERVERS WITH ANOMALOUS RECORDS PER CLASSIFCATION ############

observer_summary <- ml_complete %>%  # Combine observer-level anomaly information with the AviCRED dataset
  left_join(observer_anomaly, by = "observer_id") %>%
  mutate(
    anomaly_flag = tidyr::replace_na(anomaly_flag, 0)  # No anomaly = 0
  )

observer_summary %>%   # Summarise observers with anomalies by AviCRED class
  group_by(avicred_class) %>%
  summarise(
    Total_Observers = n(),
    Observers_with_Anomalies = sum(anomaly_flag),  # Number of anomalous observers
    Percent = 100 * Observers_with_Anomalies / Total_Observers,
    .groups = "drop"
  ) %>%
  arrange(
    factor(
      avicred_class,
      levels = c("Very Low", "Low", "Moderate", "High", "Very High")  # Order classes
    )
  )

observer_summary %>%   # Summarise observers with anomalies by Isolation Forest group
  mutate(
    IF_Class = case_when(
      ml_5pct == 1 ~ "Top 5%",       
      ml_10pct == 1 ~ "Top 10%",    
      TRUE ~ "Remaining"
    )
  ) %>%
  group_by(IF_Class) %>%
  summarise(
    Total_Observers = n(),
    Observers_with_Anomalies = sum(anomaly_flag),
    Percent = 100 * Observers_with_Anomalies / Total_Observers,
    .groups = "drop"
  )

observer_summary %>%  # Summarize observers with anomalies by cluster
  group_by(cluster) %>%
  summarise(
    Total_Observers = n(),
    Observers_with_Anomalies = sum(anomaly_flag),
    Percent = 100 * Observers_with_Anomalies / Total_Observers,
    .groups = "drop"
  ) %>%
  arrange(cluster)  

############### ANOMALOUS CHECKLISTS PER CHECKLIST CLASSIFICATION ############

checklists %>% # Summarise anomalous checklists by CAS class
  mutate(
    anomaly_flag = if_else(
      richness_flag == 1 | abundance_flag == 1,  # Flag richness or abundance anomalies
      1, 0
    )
  ) %>%
  group_by(CAS_class) %>%
  summarise(
    Total_Checklists = n(),
    Anomalous_Checklists = sum(anomaly_flag),
    Percentage_Anomalous = 100 * Anomalous_Checklists / Total_Checklists
  )

checklists %>%  # Summarise anomalous checklists by cluster
  mutate(
    anomaly_flag = if_else(
      richness_flag == 1 | abundance_flag == 1,  # Flag richness or abundance anomalies
      1, 0
    )
  ) %>%
  group_by(cluster) %>%
  summarise(
    Total_Checklists = n(),
    Anomalous_Checklists = sum(anomaly_flag),
    Percentage_Anomalous = 100 * Anomalous_Checklists / Total_Checklists
  )

checklists %>%  # Compare checklist anomalies across Isolation Forest groups
  mutate(
    anomaly_flag = if_else(
      richness_flag == 1 | abundance_flag == 1,  # Flag richness or abundance anomalies
      1, 0
    )
  ) %>%
  summarise(
    
    ## Top 5% of Isolation Forest anomalies
    `IF Top 5% Total` = sum(if_5pct == 1, na.rm = TRUE),
    `IF Top 5% Anomalous` = sum(if_5pct == 1 & anomaly_flag == 1, na.rm = TRUE),
    `IF Top 5% % Anomalous` = 100 * `IF Top 5% Anomalous` / `IF Top 5% Total`,
    
    ## Top 10% of Isolation Forest anomalies
    `IF Top 10% Total` = sum(if_10pct == 1, na.rm = TRUE),
    `IF Top 10% Anomalous` = sum(if_10pct == 1 & anomaly_flag == 1, na.rm = TRUE),
    `IF Top 10% % Anomalous` = 100 * `IF Top 10% Anomalous` / `IF Top 10% Total`,
    
    ## Checklists outside the top 10% anomaly group
    `Remaining 90% Total` = sum(if_10pct == 0, na.rm = TRUE),
    `Remaining 90% Anomalous` = sum(if_10pct == 0 & anomaly_flag == 1, na.rm = TRUE),
    `Remaining 90% % Anomalous` = 100 * `Remaining 90% Anomalous` / `Remaining 90% Total`
  )
