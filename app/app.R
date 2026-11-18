
# ============================================================
# AVICRED
# AVIAN CHECKLIST RELIABILITY EVALUATION DATABASE
# Complete Shiny application with CSV-based CAS reference engine
# ============================================================

library(shiny)
library(bslib)
library(shinyjs)
library(DT)
library(dplyr)
library(readr)
library(stringr)
library(scales)

# ============================================================
# 1. LOAD CAS REFERENCE CSV FILES
# ============================================================

reference_folder <- "FILE PATH/AVICRED app/CAS_reference"

cas_expected_file <- file.path(
  reference_folder,
  "CAS_expected_value_variation.csv"
)

cas_protocol_file <- file.path(
  reference_folder,
  "CAS_protocol_expected_variation.csv"
)

cas_residual_file <- file.path(
  reference_folder,
  "CAS_reference_residual_statistics.csv"
)

cas_counts_file <- file.path(
  reference_folder,
  "CAS_state_season_counts.csv"
)

# Check that all reference files exist

reference_files <- c(
  cas_expected_file,
  cas_protocol_file,
  cas_residual_file,
  cas_counts_file
)

missing_reference_files <- reference_files[
  !file.exists(reference_files)
]

if (length(missing_reference_files) > 0) {
  stop(
    paste(
      "The following CAS reference files are missing:",
      paste(missing_reference_files, collapse = "\n")
    )
  )
}

# Read reference files

cas_expected <- read_csv(
  cas_expected_file,
  show_col_types = FALSE
)

cas_protocol <- read_csv(
  cas_protocol_file,
  show_col_types = FALSE
)

cas_residual_stats <- read_csv(
  cas_residual_file,
  show_col_types = FALSE
)

cas_counts <- read_csv(
  cas_counts_file,
  show_col_types = FALSE
)

# Validate required columns

required_expected <- c(
  "state_code",
  "season",
  "state_season",
  "richness_expected_mean",
  "abundance_expected_mean"
)

required_residual <- c(
  "state_season",
  "richness_residual_median",
  "richness_residual_MAD",
  "abundance_residual_median",
  "abundance_residual_MAD"
)

if (!all(required_expected %in% names(cas_expected))) {
  stop("CAS expected-value CSV has missing required columns.")
}

if (!all(required_residual %in% names(cas_residual_stats))) {
  stop("CAS residual-statistics CSV has missing required columns.")
}

# Clean reference identifiers

cas_expected <- cas_expected %>%
  mutate(
    state_code = as.character(state_code),
    season = as.character(season),
    state_season = as.character(state_season)
  )

cas_residual_stats <- cas_residual_stats %>%
  mutate(
    state_season = as.character(state_season)
  )

cas_protocol <- cas_protocol %>%
  mutate(
    state_code = as.character(state_code),
    season = as.character(season),
    protocol = as.character(protocol)
  )

# Check duplicate state-season reference groups

if (anyDuplicated(cas_expected$state_season)) {
  stop("Duplicate state-season groups in CAS expected-value reference.")
}

if (anyDuplicated(cas_residual_stats$state_season)) {
  stop("Duplicate state-season groups in CAS residual reference.")
}

# Combine expected values and residual statistics

cas_reference <- cas_expected %>%
  select(
    state_season,
    reference_n = n,
    expected_richness = richness_expected_mean,
    expected_abundance = abundance_expected_mean
  ) %>%
  left_join(
    cas_residual_stats %>%
      select(
        state_season,
        richness_residual_median,
        richness_residual_MAD,
        abundance_residual_median,
        abundance_residual_MAD
      ),
    by = "state_season"
  )

# Check reference uniqueness

if (anyDuplicated(cas_reference$state_season)) {
  stop("Duplicate state-season groups in combined CAS reference.")
}

# Check residual calibration

if (
  any(
    is.na(cas_reference$richness_residual_MAD) |
    cas_reference$richness_residual_MAD <= 0
  ) ||
  any(
    is.na(cas_reference$abundance_residual_MAD) |
    cas_reference$abundance_residual_MAD <= 0
  )
) {
  warning(
    "Some state-season groups have missing or invalid residual MAD values."
  )
}

cat("\nCAS REFERENCE DATABASE LOADED\n")
cat("Expected-value groups:", nrow(cas_expected), "\n")
cat("Protocol reference rows:", nrow(cas_protocol), "\n")
cat("Residual reference groups:", nrow(cas_residual_stats), "\n")
cat("Checklist count groups:", nrow(cas_counts), "\n")

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# SAFE NUMERIC CONVERSION
# ------------------------------------------------------------

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# ------------------------------------------------------------
# STANDARDISE COLUMN NAMES
# ------------------------------------------------------------

standardise_names <- function(data) {
  
  names(data) <- tolower(names(data))
  
  names(data) <- gsub(
    "[^a-z0-9]+",
    "_",
    names(data)
  )
  
  names(data) <- gsub("^_|_$", "", names(data))
  
  aliases <- c(
    
    "sampling_event_identifier" = "checklist_id",
    "observation_date" = "obs_date",
    "date" = "obs_date",
    
    "time_observations_started" = "start_time",
    "start_hour" = "start_time",
    
    "state" = "state_code",
    "county" = "county_code",
    
    "observation_type" = "protocol",
    "protocol_type" = "protocol",
    
    "duration_minutes" = "duration_min",
    "effort_distance_km" = "distance_km",
    "number_observers" = "n_observers",
    
    "all_species_reported" = "complete_flag",
    
    "scientific_name" = "species_name",
    "observation_count" = "count",
    
    "richness" = "observed_richness",
    "abundance" = "observed_abundance",
    "abundance_per_checklist" = "observed_abundance",
    
    "expected_species_richness" = "expected_richness"
    
  )
  
  for (source_name in names(aliases)) {
    
    target_name <- unname(aliases[source_name])
    
    if (
      source_name %in% names(data) &&
      !target_name %in% names(data)
    ) {
      
      names(data)[
        names(data) == source_name
      ] <- target_name
      
    }
  }
  
  data
}

# ------------------------------------------------------------
# FIRST NON-MISSING VALUE
# ------------------------------------------------------------

first_nonmissing <- function(x) {
  
  valid <- x[
    !is.na(x) &
      trimws(as.character(x)) != ""
  ]
  
  if (length(valid) == 0) {
    return(NA)
  }
  
  valid[[1]]
}

# ------------------------------------------------------------
# CHECK WHETHER CHECKLIST METADATA ARE CONSISTENT
# ------------------------------------------------------------

consistent_value <- function(x) {
  
  valid <- unique(
    trimws(
      as.character(
        x[
          !is.na(x) &
            trimws(as.character(x)) != ""
        ]
      )
    )
  )
  
  if (length(valid) == 0) {
    return(NA_character_)
  }
  
  if (length(valid) > 1) {
    return(NA_character_)
  }
  
  valid[[1]]
}

# ------------------------------------------------------------
# STANDARDISE STATE CODES
# ------------------------------------------------------------

normalise_state_code <- function(x) {
  
  x <- toupper(trimws(as.character(x)))
  
  x[x == ""] <- NA_character_
  
  x
}

# ------------------------------------------------------------
# STANDARDISE LOCALITY IDS
# ------------------------------------------------------------

normalise_locality_id <- function(x) {
  
  x <- trimws(as.character(x))
  
  x[x == ""] <- NA_character_
  
  x
}

# ------------------------------------------------------------
# PARSE EBIRD OBSERVATION DATES
# ------------------------------------------------------------

parse_ebird_date <- function(x) {
  
  if (inherits(x, "Date")) {
    return(x)
  }
  
  x <- trimws(as.character(x))
  
  x[x == ""] <- NA_character_
  
  parsed <- as.Date(
    rep(NA_character_, length(x))
  )
  
  date_formats <- c(
    "%Y-%m-%d",
    "%m/%d/%Y",
    "%d/%m/%Y",
    "%Y/%m/%d",
    "%d-%m-%Y",
    "%m-%d-%Y"
  )
  
  for (fmt in date_formats) {
    
    missing <- is.na(parsed) & !is.na(x)
    
    if (!any(missing)) {
      break
    }
    
    parsed[missing] <- suppressWarnings(
      as.Date(x[missing], format = fmt)
    )
  }
  
  parsed
}

# ------------------------------------------------------------
# DERIVE SEASON FROM OBSERVATION DATE
# ------------------------------------------------------------

derive_season <- function(date_value) {
  
  parsed_date <- parse_ebird_date(date_value)
  
  month_value <- as.integer(
    format(parsed_date, "%m")
  )
  
  dplyr::case_when(
    
    month_value %in% c(11, 12, 1, 2) ~ "Winter",
    
    month_value %in% 3:5 ~ "Pre-monsoon",
    
    month_value %in% 6:9 ~ "Monsoon",
    
    month_value == 10 ~ "Post-monsoon",
    
    TRUE ~ NA_character_
    
  )
}

# ------------------------------------------------------------
# BUILD THE STATE-SEASON KEY USING THE ACTUAL REFERENCE CSV
# This avoids assuming whether the CSV uses "_", " | ", or
# another separator. The reference key is the authority.
# ------------------------------------------------------------

lookup_state_season_key <- function(state_code, season) {
  
  state_code <- normalise_state_code(state_code)
  season <- trimws(as.character(season))
  season[season == ""] <- NA_character_
  
  ref <- cas_expected %>%
    transmute(
      .state = normalise_state_code(state_code),
      .season = tolower(trimws(as.character(season))),
      .key = trimws(as.character(state_season))
    ) %>%
    distinct()
  
  lookup <- paste(ref$.state, ref$.season, sep = "\r")
  query <- paste(state_code, tolower(season), sep = "\r")
  
  result <- ref$.key[match(query, lookup)]
  result[is.na(state_code) | is.na(season)] <- NA_character_
  result
}

# ------------------------------------------------------------
# ANOMALY FLAGS
# ------------------------------------------------------------

make_flag <- function(z) {
  
  ifelse(
    is.na(z),
    NA_character_,
    ifelse(abs(z) >= 3.5, "Yes", "No")
  )
}

# ============================================================
# 2A. VALIDATE LOCATION AND SEASON
# ============================================================

validate_location_season <- function(data) {
  
  data <- standardise_names(data)
  
  # Ensure metadata columns exist
  
  if (!"state_code" %in% names(data)) {
    data$state_code <- NA_character_
  }
  
  if (!"state_season" %in% names(data)) {
    data$state_season <- NA_character_
  }
  
  if (!"locality_id" %in% names(data)) {
    data$locality_id <- NA_character_
  }
  
  if (!"obs_date" %in% names(data)) {
    data$obs_date <- NA_character_
  }
  
  if (!"season" %in% names(data)) {
    data$season <- NA_character_
  }
  
  # Preserve the originally supplied state-season key
  
  supplied_state_season <- trimws(
    as.character(data$state_season)
  )
  
  supplied_state_season[
    supplied_state_season == ""
  ] <- NA_character_
  
  # Extract state code from existing state-season key
  # Expected format: STATE CODE | Season
  
  key_parts <- strsplit(
    supplied_state_season,
    " | ",
    fixed = TRUE
  )
  
  key_state <- vapply(
    key_parts,
    function(x) {
      if (length(x) >= 2) {
        x[1]
      } else {
        NA_character_
      }
    },
    character(1)
  )
  
  key_season <- vapply(
    key_parts,
    function(x) {
      if (length(x) >= 2) {
        x[length(x)]
      } else {
        NA_character_
      }
    },
    character(1)
  )
  
  # Use state_code if supplied; otherwise extract it from the existing state_season key
  
  supplied_state_code <- normalise_state_code(
    data$state_code
  )
  
  use_key_state <- is.na(supplied_state_code)
  
  supplied_state_code[use_key_state] <- normalise_state_code(
    key_state[use_key_state]
  )
  
  # Parse date and derive season
  
  parsed_date <- parse_ebird_date(data$obs_date)
  
  date_season <- derive_season(parsed_date)
  
  # Date-derived season is authoritative
  
  supplied_season <- trimws(
    as.character(data$season)
  )
  
  supplied_season[supplied_season == ""] <- NA_character_
  
  # Compare any supplied season with date-derived season
  
  season_conflict <- dplyr::case_when(
    
    is.na(date_season) ~ "Cannot verify supplied season",
    
    !is.na(supplied_season) &
      tolower(supplied_season) != tolower(date_season) ~ "Yes",
    
    is.na(supplied_season) &
      !is.na(key_season) &
      tolower(trimws(key_season)) != tolower(date_season) ~ "Yes",
    
    TRUE ~ "No"
    
  )
  
  canonical_state_season <- lookup_state_season_key(
    supplied_state_code,
    date_season
  )
  
  existing_key_conflict <- !is.na(supplied_state_season) &
    !is.na(canonical_state_season) &
    tolower(supplied_state_season) !=
    tolower(canonical_state_season)
  
  # Validate state code against expected-value reference
  
  valid_states <- unique(
    normalise_state_code(cas_expected$state_code)
  )
  
  valid_states <- valid_states[
    !is.na(valid_states)
  ]
  
  data <- data %>%
    mutate(
      
      state_code = supplied_state_code,
      
      obs_date_parsed = parsed_date,
      
      season_from_date = date_season,
      
      season = date_season,
      
      season_conflict = season_conflict,
      
      state_validation_status = case_when(
        
        is.na(state_code) ~ "Missing state code",
        
        state_code %in% valid_states ~ "Valid state code",
        
        TRUE ~ "State code not found in CAS reference"
        
      ),
      
      locality_id = normalise_locality_id(locality_id),
      
      locality_validation_status = case_when(
        
        is.na(locality_id) ~ "Missing locality ID",
        
        TRUE ~ "Locality ID supplied"
        
      ),
      
      date_validation_status = case_when(
        
        is.na(obs_date_parsed) ~
          "Invalid or missing observation date",
        
        TRUE ~
          "Valid observation date"
        
      ),
      
      season_validation_status = case_when(
        
        is.na(season) ~
          "Season unavailable: date missing or invalid",
        
        season_conflict == "Yes" ~
          "Supplied season conflicts with observation date",
        
        TRUE ~
          "Season derived from observation date"
        
      ),
      
      state_season = ifelse(
        
        state_validation_status == "Valid state code" &
          !is.na(season) &
          season_conflict != "Yes" &
          !existing_key_conflict,
        
        canonical_state_season,
        
        NA_character_
        
      )
      
    )
  
  # Check exact state-season coverage
  
  valid_state_seasons <- unique(
    cas_reference$state_season
  )
  
  data <- data %>%
    mutate(
      
      reference_match_status = case_when(
        
        is.na(state_season) ~
          "Not matched: state or season unavailable",
        
        state_season %in% valid_state_seasons ~
          "Matched to CAS reference",
        
        TRUE ~
          "Not matched: state-season reference unavailable"
        
      ),
      
      location_season_validation = case_when(
        
        state_validation_status != "Valid state code" ~
          state_validation_status,
        
        date_validation_status != "Valid observation date" ~
          date_validation_status,
        
        season_conflict == "Yes" ~
          "Supplied season conflicts with observation date",
        
        existing_key_conflict ~
          "Existing state-season key conflicts with state/date",
        
        reference_match_status != "Matched to CAS reference" ~
          reference_match_status,
        
        TRUE ~
          "Valid state-season reference assignment"
        
      )
      
    )
  
  data
}

# ============================================================
# 3. CAS SCORING FUNCTION — REVISED
#    Derives state_season from state_code + obs_date
#    Matches checklist data to the fixed CAS reference table
# ============================================================

calculate_CAS <- function(data) {
  
  data <- standardise_names(data)
  
  if (!"state_code" %in% names(data)) data$state_code <- NA_character_
  if (!"obs_date" %in% names(data)) data$obs_date <- NA_character_
  if (!"state_season" %in% names(data)) data$state_season <- NA_character_
  if (!"season" %in% names(data)) data$season <- NA_character_
  if (!"observed_richness" %in% names(data)) data$observed_richness <- NA_real_
  if (!"observed_abundance" %in% names(data)) data$observed_abundance <- NA_real_
  
  data <- data %>%
    mutate(
      state_code = normalise_state_code(state_code),
      observed_richness = safe_numeric(observed_richness),
      observed_abundance = safe_numeric(observed_abundance)
    )
  
  # Derive season and state_season through the same validator used
  # by species-level CSV processing. 
  data <- validate_location_season(data)
  
  # Remove any old score/reference columns before joining afresh.
  old_columns <- intersect(
    c(
      "reference_n", "expected_richness", "expected_abundance",
      "richness_residual_median", "abundance_residual_median",
      "richness_residual_MAD", "abundance_residual_MAD",
      "richness_residual", "abundance_residual",
      "richness_robust_z", "abundance_robust_z",
      "richness_CAS_component", "abundance_CAS_component",
      "proto_CAS", "CAS", "CAS_class", "CAS_scoring_status",
      "richness_anomaly_flag", "abundance_anomaly_flag",
      "either_anomaly_flag", "richness_anomaly_direction",
      "abundance_anomaly_direction"
    ),
    names(data)
  )
  if (length(old_columns)) data <- data %>% select(-all_of(old_columns))
  
  required_reference_columns <- c(
    "state_season", "reference_n", "expected_richness",
    "expected_abundance", "richness_residual_median",
    "abundance_residual_median", "richness_residual_MAD",
    "abundance_residual_MAD"
  )
  missing_reference_columns <- setdiff(required_reference_columns, names(cas_reference))
  if (length(missing_reference_columns)) {
    stop(
      "Missing columns in cas_reference: ",
      paste(missing_reference_columns, collapse = ", ")
    )
  }
  if (anyDuplicated(cas_reference$state_season)) {
    stop("Duplicate state_season entries found in cas_reference.")
  }
  
  # The validator has a diagnostic match status of its own. Remove it
  # before the join to avoid .x/.y columns and calculate one definitive status.
  data <- data %>% select(-any_of("reference_match_status"))
  
  data <- data %>%
    left_join(
      cas_reference %>%
        select(all_of(required_reference_columns)),
      by = "state_season"
    ) %>%
    mutate(
      reference_match_status = case_when(
        is.na(state_season) | state_season == "" ~ "Missing state-season",
        is.na(expected_richness) | is.na(expected_abundance) ~
          "State-season reference unavailable",
        TRUE ~ "Matched to CAS reference"
      ),
      richness_residual = observed_richness - expected_richness,
      abundance_residual = observed_abundance - expected_abundance,
      richness_robust_z = ifelse(
        !is.na(richness_residual) & !is.na(richness_residual_median) &
          !is.na(richness_residual_MAD) & richness_residual_MAD > 0,
        (richness_residual - richness_residual_median) / richness_residual_MAD,
        NA_real_
      ),
      abundance_robust_z = ifelse(
        !is.na(abundance_residual) & !is.na(abundance_residual_median) &
          !is.na(abundance_residual_MAD) & abundance_residual_MAD > 0,
        (abundance_residual - abundance_residual_median) / abundance_residual_MAD,
        NA_real_
      ),
      richness_CAS_component = exp(-abs(richness_robust_z) / 3.5),
      abundance_CAS_component = exp(-abs(abundance_robust_z) / 3.5),
      proto_CAS = (richness_CAS_component + abundance_CAS_component) / 2,
      CAS = 1 - proto_CAS,
      richness_anomaly_flag = make_flag(richness_robust_z),
      abundance_anomaly_flag = make_flag(abundance_robust_z),
      either_anomaly_flag = case_when(
        richness_anomaly_flag == "Yes" | abundance_anomaly_flag == "Yes" ~ "Yes",
        richness_anomaly_flag == "No" & abundance_anomaly_flag == "No" ~ "No",
        TRUE ~ NA_character_
      ),
      richness_anomaly_direction = case_when(
        is.na(richness_robust_z) ~ NA_character_,
        richness_robust_z >= 3.5 ~ "Higher than expected",
        richness_robust_z <= -3.5 ~ "Lower than expected",
        TRUE ~ "Within threshold"
      ),
      abundance_anomaly_direction = case_when(
        is.na(abundance_robust_z) ~ NA_character_,
        abundance_robust_z >= 3.5 ~ "Higher than expected",
        abundance_robust_z <= -3.5 ~ "Lower than expected",
        TRUE ~ "Within threshold"
      ),
      CAS_scoring_status = case_when(
        is.na(state_code) | state_code == "" ~ "Not scored: Missing state code",
        date_validation_status != "Valid observation date" ~
          "Not scored: Invalid or missing observation date",
        reference_match_status != "Matched to CAS reference" ~
          "Not scored: State-season reference unavailable",
        is.na(richness_residual_MAD) | is.na(abundance_residual_MAD) |
          richness_residual_MAD <= 0 | abundance_residual_MAD <= 0 ~
          "Not scored: Invalid residual calibration",
        is.na(observed_richness) | is.na(observed_abundance) ~
          "Not scored: Observed richness or abundance unavailable",
        is.na(richness_robust_z) | is.na(abundance_robust_z) ~
          "Not scored: Residual calculation unavailable",
        TRUE ~ "Scored"
      ),
      CAS_class = ifelse(
        CAS_scoring_status == "Scored",
        as.character(cut(
          CAS,
          breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
          labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
          right = FALSE
        )),
        NA_character_
      )
    )
  
  score_columns <- c(
    "richness_residual", "abundance_residual", "richness_robust_z",
    "abundance_robust_z", "richness_CAS_component",
    "abundance_CAS_component", "proto_CAS", "CAS"
  )
  data <- data %>%
    mutate(across(
      all_of(score_columns),
      ~ ifelse(CAS_scoring_status == "Scored", .x, NA_real_)
    ))
  
  data
}

# ------------------------------------------------------------
# CAS FAILURE FALLBACK
# Preserves all checklist rows and diagnostics.
# ------------------------------------------------------------

cas_fallback <- function(data, reason = "Unknown CAS error") {
  
  # Make a best-effort attempt to retain date-derived season and the
  # exact reference key even when the scoring/join stage has failed.
  data <- tryCatch(
    validate_location_season(data),
    error = function(e) data
  )
  
  score_columns <- c(
    "expected_richness", "expected_abundance",
    "richness_residual_median", "abundance_residual_median",
    "richness_residual_MAD", "abundance_residual_MAD",
    "richness_residual", "abundance_residual",
    "richness_robust_z", "abundance_robust_z",
    "richness_CAS_component", "abundance_CAS_component",
    "proto_CAS", "CAS"
  )
  for (nm in score_columns) data[[nm]] <- NA_real_
  data$CAS_class <- NA_character_
  data$richness_anomaly_flag <- NA_character_
  data$abundance_anomaly_flag <- NA_character_
  data$either_anomaly_flag <- NA_character_
  data$richness_anomaly_direction <- NA_character_
  data$abundance_anomaly_direction <- NA_character_
  data$reference_match_status <- if ("state_season" %in% names(data)) {
    ifelse(is.na(data$state_season), "Not matched: state-season unavailable",
           "CAS engine unavailable")
  } else {
    rep("CAS engine unavailable", nrow(data))
  }
  data$CAS_scoring_status <- paste0("CAS unavailable: ", reason)
  data$CAS_engine_error <- as.character(reason)
  data
}

# ============================================================
# 4. PROCESS CHECKLIST-LEVEL DATA
# ============================================================

process_checklist_data <- function(data) {
  
  data <- standardise_names(data)
  
  required <- c("checklist_id", "observer_id", "duration_min", "distance_km",
                "observed_richness", "observed_abundance")
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # API and checklist CSVs may omit these fields. Keep the records;
  # CAS will mark them unscored rather than aborting the whole upload.
  if (!"state_code" %in% names(data)) data$state_code <- NA_character_
  if (!"obs_date" %in% names(data)) data$obs_date <- NA_character_
  if (!"state_season" %in% names(data)) data$state_season <- NA_character_
  if (!"season" %in% names(data)) data$season <- NA_character_
  if (!"locality_id" %in% names(data)) data$locality_id <- NA_character_
  if (!"complete_flag" %in% names(data)) data$complete_flag <- NA
  
  numeric_columns <- c("duration_min", "distance_km",
                       "observed_richness", "observed_abundance")
  data <- data %>% mutate(across(all_of(numeric_columns), safe_numeric))
  
  # Exactly one CAS attempt. On engine/reference failure, retain rows
  # and return explicit unavailable status without assigning a score.
  data <- tryCatch(
    calculate_CAS(data),
    error = function(e) cas_fallback(data, conditionMessage(e))
  )
  
  data <- data %>%
    mutate(
      duration_effort_flag = case_when(
        is.na(duration_min) ~ NA_character_,
        duration_min >= 150 ~ "Yes",
        TRUE ~ "No"
      ),
      distance_effort_flag = case_when(
        is.na(distance_km) ~ NA_character_,
        distance_km >= 4 ~ "Yes",
        TRUE ~ "No"
      ),
      combined_effort_flag = case_when(
        duration_effort_flag == "Yes" | distance_effort_flag == "Yes" ~ "Yes",
        duration_effort_flag == "No" & distance_effort_flag == "No" ~ "No",
        TRUE ~ NA_character_
      )
    )
  
  data
}

# ============================================================
# 5. PROCESS SPECIES-LEVEL EBD DATA
# ============================================================

process_species_data <- function(data) {
  
  data <- standardise_names(data)
  
  required <- c(
    "checklist_id",
    "observer_id",
    "species_name",
    "count"
  )
  
  missing_columns <- setdiff(
    required,
    names(data)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "Missing required EBD columns:",
        paste(missing_columns, collapse = ", ")
      )
    )
    
  }
  
  # Validate state-season information before aggregation
  
  data <- validate_location_season(data)
  
  # Ensure metadata columns exist for later summarisation
  
  if (!"locality_id" %in% names(data)) {
    data$locality_id <- NA_character_
  }
  
  if (!"obs_date" %in% names(data)) {
    data$obs_date <- NA_character_
  }
  
  if (!"season" %in% names(data)) {
    data$season <- NA_character_
  }
  
  # ----------------------------------------------------------
  # 2. CLEAN IDENTIFIERS AND SPECIES NAMES
  # ----------------------------------------------------------
  
  data <- data %>%
    mutate(
      
      checklist_id = as.character(checklist_id),
      
      observer_id = as.character(observer_id),
      
      species_name = trimws(
        as.character(species_name)
      ),
      
      state_season = trimws(
        as.character(state_season)
      )
      
    )
  
  if (any(
    is.na(data$checklist_id) |
    data$checklist_id == ""
  )) {
    
    stop("Some records have missing checklist identifiers.")
    
  }
  
  if (any(
    is.na(data$species_name) |
    data$species_name == ""
  )) {
    
    stop("Some records have missing scientific names.")
    
  }
  
  # ----------------------------------------------------------
  # 3. INTERPRET SPECIES COUNTS
  # Treat "X" as a minimum of 1 individual
  # ----------------------------------------------------------
  
  original_count <- trimws(
    as.character(data$count)
  )
  
  numeric_count <- ifelse(
    !is.na(original_count) &
      toupper(original_count) == "X",
    1,
    safe_numeric(original_count)
  )
  
  invalid_count <- !is.na(original_count) &
    original_count != "" &
    is.na(numeric_count)
  
  if (any(invalid_count)) {
    
    stop("Unexpected non-numeric values in OBSERVATION COUNT.")
    
  }
  
  if (any(numeric_count < 0, na.rm = TRUE)) {
    
    stop("Species counts cannot be negative.")
    
  }
  
  data$count_numeric <- numeric_count
  
  # X is now treated as a numeric count of 1
  data$presence_only <- FALSE
  
  # ----------------------------------------------------------
  # A. SPECIES-LEVEL ABUNDANCE
  # ----------------------------------------------------------
  
  species_data <- data %>%
    group_by(checklist_id, species_name) %>%
    summarise(
      
      species_abundance = if (
        all(is.na(count_numeric))
      ) {
        
        NA_real_
        
      } else {
        
        sum(count_numeric, na.rm = TRUE)
        
      },
      
      abundance_incomplete = if (
        all(is.na(count_numeric))
      ) {
        
        "Yes"
        
      } else {
        
        "No"
        
      },
      
      .groups = "drop"
      
    )
  
  # ----------------------------------------------------------
  # A2. IDENTIFY HIGHEST-ABUNDANCE SPECIES
  # ----------------------------------------------------------
  
  highest_abundance_species <- species_data %>%
    group_by(checklist_id) %>%
    summarise(
      
      highest_species_abundance = if (
        all(is.na(species_abundance))
      ) {
        
        NA_real_
        
      } else {
        
        max(species_abundance, na.rm = TRUE)
        
      },
      
      highest_abundance_species = if (
        all(is.na(species_abundance))
      ) {
        
        NA_character_
        
      } else {
        
        paste(
          species_name[
            !is.na(species_abundance) &
              species_abundance == max(
                species_abundance,
                na.rm = TRUE
              )
          ],
          collapse = "; "
        )
        
      },
      
      .groups = "drop"
      
    )
  
  # ----------------------------------------------------------
  # B. CHECKLIST-LEVEL METRICS
  # ----------------------------------------------------------
  
  checklist_metrics <- species_data %>%
    group_by(checklist_id) %>%
    summarise(
      
      observed_richness = n_distinct(species_name),
      
      observed_abundance = if (
        any(is.na(species_abundance))
      ) {
        
        NA_real_
        
      } else {
        
        sum(species_abundance)
        
      },
      
      abundance_incomplete = if (
        any(is.na(species_abundance))
      ) {
        
        "Yes"
        
      } else {
        
        "No"
        
      },
      
      # Shannon-Wiener diversity
      
      shannon_wiener = {
        
        abund <- species_abundance
        
        if (
          any(is.na(abund)) ||
          sum(abund) <= 0
        ) {
          
          NA_real_
          
        } else {
          
          p <- abund[abund > 0] / sum(abund)
          
          -sum(p * log(p))
          
        }
        
      },
      
      # Simpson diversity
      
      simpson_diversity = {
        
        abund <- species_abundance
        
        if (
          any(is.na(abund)) ||
          sum(abund) <= 0
        ) {
          
          NA_real_
          
        } else {
          
          p <- abund[abund > 0] / sum(abund)
          
          1 - sum(p^2)
          
        }
        
      },
      
      # Pielou evenness
      
      pielou_evenness = {
        
        abund <- species_abundance
        
        if (
          any(is.na(abund)) ||
          sum(abund) <= 0
        ) {
          
          NA_real_
          
        } else {
          
          p <- abund[abund > 0] / sum(abund)
          
          H <- -sum(p * log(p))
          
          S <- length(p)
          
          if (S <= 1) {
            
            NA_real_
            
          } else {
            
            H / log(S)
            
          }
          
        }
        
      },
      
      # Berger-Parker dominance
      
      berger_parker = {
        
        abund <- species_abundance
        
        if (
          any(is.na(abund)) ||
          sum(abund) <= 0
        ) {
          
          NA_real_
          
        } else {
          
          max(abund) / sum(abund)
          
        }
        
      },
      
      # Bias-corrected Chao1
      
      chao1 = {
        
        abund <- species_abundance
        
        if (
          any(is.na(abund)) ||
          sum(abund) <= 0
        ) {
          
          NA_real_
          
        } else {
          
          abund <- abund[abund > 0]
          
          Sobs <- length(abund)
          
          F1 <- sum(abund == 1)
          
          F2 <- sum(abund == 2)
          
          Sobs +
            (F1 * (F1 - 1)) /
            (2 * (F2 + 1))
          
        }
        
      },
      
      .groups = "drop"
      
    )
  
  checklist_metrics <- checklist_metrics %>%
    left_join(
      highest_abundance_species,
      by = "checklist_id"
    )
  
  # ----------------------------------------------------------
  # C. RETAIN CHECKLIST METADATA
  # ----------------------------------------------------------
  
  metadata_columns <- intersect(
    
    c(
      "obs_date",
      "start_time",
      "state_season",
      "season",
      "county_code",
      "locality_id",
      "locality",
      "latitude",
      "longitude",
      "protocol",
      "protocol_name",
      "duration_min",
      "distance_km",
      "n_observers",
      "complete_flag",
      
      # Retain validation information
      
      "state_code",
      "state_validation_status",
      "locality_validation_status",
      "date_validation_status",
      "season_validation_status",
      "reference_match_status",
      "location_season_validation"
      
    ),
    
    names(data)
    
  )
  
  # ----------------------------------------------------------
  # D. SUMMARISE CHECKLIST METADATA
  # ----------------------------------------------------------
  
  checklist_metadata <- data %>%
    group_by(checklist_id) %>%
    summarise(
      
      observer_id = consistent_value(observer_id),
      
      across(
        all_of(metadata_columns),
        consistent_value
      ),
      
      # Detect state-season conflicts
      
      state_metadata_conflict =
        
        n_distinct(
          na.omit(
            trimws(as.character(state_season))
          )
        ) > 1,
      
      # Detect locality conflicts
      
      locality_metadata_conflict =
        
        n_distinct(
          na.omit(
            normalise_locality_id(locality_id)
          )
        ) > 1,
      
      # Detect observation date conflicts
      
      date_metadata_conflict =
        
        n_distinct(
          na.omit(
            as.character(
              parse_ebird_date(obs_date)
            )
          )
        ) > 1,
      
      .groups = "drop"
      
    )
  
  # ----------------------------------------------------------
  # E. COMBINE CHECKLIST METRICS AND METADATA
  # ----------------------------------------------------------
  
  checklist_data <- checklist_metrics %>%
    left_join(
      checklist_metadata,
      by = "checklist_id"
    )
  
  # ----------------------------------------------------------
  # F. ENSURE REQUIRED EFFORT AND METADATA COLUMNS EXIST
  # ----------------------------------------------------------
  
  if (!"duration_min" %in% names(checklist_data)) {
    
    checklist_data$duration_min <- NA_real_
    
  }
  
  if (!"distance_km" %in% names(checklist_data)) {
    
    checklist_data$distance_km <- NA_real_
    
  }
  
  if (!"state_season" %in% names(checklist_data)) {
    
    checklist_data$state_season <- NA_character_
    
  }
  
  if (!"locality_id" %in% names(checklist_data)) {
    
    checklist_data$locality_id <- NA_character_
    
  }
  
  if (!"obs_date" %in% names(checklist_data)) {
    
    checklist_data$obs_date <- NA_character_
    
  }
  
  # ----------------------------------------------------------
  # G. IDENTIFY CHECKLIST METADATA CONFLICTS
  # ----------------------------------------------------------
  
  checklist_data <- checklist_data %>%
    mutate(
      
      metadata_conflict = case_when(
        
        state_metadata_conflict ~
          "State-season conflict",
        
        locality_metadata_conflict ~
          "Locality ID conflict",
        
        date_metadata_conflict ~
          "Observation date conflict",
        
        TRUE ~
          "No metadata conflict"
        
      )
      
    )
  
  # ----------------------------------------------------------
  # H. PROCESS CHECKLIST METRICS AND CALCULATE CAS
  # ----------------------------------------------------------
  
  checklist_data <- process_checklist_data(checklist_data)
  
  # ----------------------------------------------------------
  # I. OVERRIDE SCORING STATUS FOR CONFLICTING METADATA
  # ----------------------------------------------------------
  
  checklist_data <- checklist_data %>%
    mutate(
      
      CAS_scoring_status = ifelse(
        
        metadata_conflict != "No metadata conflict",
        
        paste(
          "Not scored:",
          metadata_conflict
        ),
        
        CAS_scoring_status
        
      ),
      
      CAS = ifelse(
        
        metadata_conflict == "No metadata conflict",
        
        CAS,
        
        NA_real_
        
      ),
      
      CAS_class = ifelse(
        
        metadata_conflict == "No metadata conflict",
        
        CAS_class,
        
        NA_character_
        
      )
      
    )
  
  # ----------------------------------------------------------
  # J. RETURN CHECKLIST AND SPECIES DATA
  # ----------------------------------------------------------
  
  return(
    
    list(
      
      checklists = checklist_data,
      
      species = species_data
      
    )
    
  )
  
}

# ============================================================
# 6. READ CSV / TSV / TXT FILES
# ============================================================

read_ebird_file <- function(file_path, file_name) {
  
  extension <- tolower(
    tools::file_ext(file_name)
  )
  
  if (extension %in% c("txt", "tsv")) {
    
    data <- read_delim(
      file_path,
      delim = "\t",
      show_col_types = FALSE,
      name_repair = "minimal",
      progress = FALSE
    )
    
  } else if (extension == "csv") {
    
    data <- read_csv(
      file_path,
      show_col_types = FALSE,
      name_repair = "minimal",
      progress = FALSE
    )
    
  } else {
    
    stop(
      "Unsupported file format. Upload CSV, TSV or TXT."
    )
  }
  
  data
}

# ============================================================
# 7. USER INTERFACE
# ============================================================

ui <- page_navbar(
  
  title = div(
    class = "brand",
    icon("dove"),
    span("AviCRED")
  ),
  
  id = "main_nav",
  
  theme = bs_theme(
    
    version = 5,
    
    primary = "#245B45",
    secondary = "#8DAF91",
    success = "#477A59",
    info = "#A9C5A5",
    warning = "#D8C49A",
    danger = "#A85C4B",
    
    bg = "#F8F7F0",
    fg = "#263B30",
    
    base_font = "Times New Roman, Times, serif",
    heading_font = "Times New Roman, Times, serif"
  ),
  
  header = tags$head(
    
    tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = "style.css"
    )
  ),
  
  # ==========================================================
  # HOME PAGE
  # ==========================================================
  
  nav_panel(
    "Home",
    
    div(
      class = "home-page",
      
      div(
        class = "hero",
        
        div(
          class = "hero-content",
          
          p(
            class = "hero-kicker",
            "CITIZEN SCIENCE  |  BIRD DATA  |  RESEARCH"
          ),
          
          h1(
            class = "hero-title",
            "AVICRED"
          ),
          
          h2(
            class = "hero-full-name",
            "Avian Checklist Reliability Evaluation Database"
          ),
          
          div(
            class = "hero-divider"
          ),
          
          h3(
            "Making bird observation data more reliable."
          ),
          
          p(
            class = "hero-description",
            paste(
              "A transparent and interpretable framework for",
              "checklist-level assessment of bird observation data."
            )
          ),
          
          p(
            class = "hero-description",
            paste(
              "Explore checklist characteristics, species richness,",
              "abundance patterns, ecological indices and sampling effort."
            )
          ),
          
          actionButton(
            "start_assessment",
            
            tagList(
              icon("clipboard-check"),
              " Start assessment"
            ),
            
            class = "btn btn-light btn-lg hero-button"
          )
        )
      ),
      
      div(
        class = "home-intro",
        
        h2(
          "Understanding the quality of bird observation data"
        ),
        
        p(
          paste(
            "AviCRED is designed to support birders, researchers and",
            "data curators in examining the statistical characteristics",
            "of bird observation checklists. It brings together checklist",
            "metrics, anomaly diagnostics and ecological information in",
            "a single assessment interface."
          )
        )
      ),
      
      layout_columns(
        
        card(
          class = "feature-card feature-card-green",
          
          card_body(
            
            icon(
              "clipboard-check",
              class = "feature-icon"
            ),
            
            h3("Checklist assessment"),
            
            p(
              "Examine individual bird observation checklists using",
              "richness, abundance and sampling-effort information."
            )
          )
        ),
        
        card(
          class = "feature-card feature-card-cream",
          
          card_body(
            
            icon(
              "chart-line",
              class = "feature-icon"
            ),
            
            h3("Anomaly diagnostics"),
            
            p(
              "Identify statistically unusual richness and abundance",
              "patterns while distinguishing atypical observations",
              "from confirmed data errors."
            )
          )
        ),
        
        card(
          class = "feature-card feature-card-sage",
          
          card_body(
            
            icon(
              "leaf",
              class = "feature-icon"
            ),
            
            h3("Ecological information"),
            
            p(
              "Explore species diversity, evenness, dominance and",
              "rank-abundance patterns alongside checklist metrics."
            )
          )
        ),
        
        col_widths = c(4, 4, 4)
      ),
      
      div(
        class = "home-footer",
        
        h3(
          "Transparent assessment. Meaningful interpretation."
        ),
        
        p(
          paste(
            "AviCRED is intended to support the evaluation of",
            "citizen-science bird data, not to automatically declare",
            "observations correct or incorrect."
          )
        )
      )
    )
  ),
  
  # ==========================================================
  # ASSESSMENT PAGE
  # ==========================================================
  
  nav_panel(
    "Assess data",
    
    h2("Checklist assessment"),
    
    p(
      "Upload checklist-level or species-level eBird data ",
      "or retrieve checklists directly from eBird."
    ),
    
    layout_sidebar(
      
      sidebar = sidebar(
        
        radioButtons(
          "input_method",
          "Select input method",
          
          choices = c(
            "Checklist-level CSV" = "checklist",
            "Species-level EBD" = "species",
            "eBird checklist URLs" = "url"
          ),
          
          selected = "checklist"
        ),
        
        # ================================================
        # CSV / EBD FILE INPUT
        # ================================================
        
        conditionalPanel(
          condition = "input.input_method != 'url'",
          
          fileInput(
            "csv_file",
            "Upload eBird EBD / CSV / TSV file",
            accept = c(".csv", ".txt", ".tsv")
          )
        ),
        
        # ================================================
        # EBIRD API INPUT
        # ================================================
        
        conditionalPanel(
          condition = "input.input_method == 'url'",
          
          passwordInput(
            "ebird_api_key",
            "Enter your eBird API key",
            placeholder = "Paste your personal eBird API key",
            width = "100%"
          ),
          
          textAreaInput(
            "checklist_urls",
            "Paste eBird checklist URLs",
            
            placeholder = paste(
              "https://ebird.org/checklist/S123456789",
              "https://ebird.org/checklist/S123456790",
              "One checklist URL per line"
            ),
            
            rows = 6,
            width = "100%"
          ),
          
          helpText(
            "Enter your API key first, followed by one or more checklist URLs."
          ),
          
          helpText(
            "A maximum of 50 checklists can be retrieved per assessment."
          )
        ),
        
        # ================================================
        # RUN ASSESSMENT
        # ================================================
        
        actionButton(
          "run_assessment",
          "Process data",
          class = "btn-primary w-100"
        )
      ),
      
      card(
        
        card_header("Input requirements"),
        
        uiOutput("input_instructions"),
        
        hr(),
        
        textOutput("processing_status")
      )
    )
  ),
  
  # ==========================================================
  # RESULTS PAGE
  # ==========================================================
  
  nav_panel(
    "Results",
    
    div(
      style = "
        width: 100%;
        min-height: 100vh;
        padding: 10px 15px;
        overflow-x: hidden;
        overflow-y: visible;
      ",
      
      h2("Assessment results"),
      
      p(
        "Checklist-level assessment, CAS scores and diagnostic results."
      ),
      
      # ======================================================
      # DESCRIPTIVE SUMMARY
      # ======================================================
      
      card(
        
        card_header("Descriptive summary"),
        
        h4("Overall dataset summary"),
        
        tableOutput("overall_summary_table"),
        
        hr(),
        
        # ----------------------------------------------------
        # INDIVIDUAL CHECKLIST EXPLORER
        # ----------------------------------------------------
        
        h4("Individual checklist explorer"),
        
        p(
          "Select a checklist ID to view its individual sampling and observation details."
        ),
        
        uiOutput("checklist_selector"),
        
        tableOutput("individual_checklist_table"),
        
        hr(),
        
        # ----------------------------------------------------
        # SEASON-WISE SUMMARY
        # ----------------------------------------------------
        
        h4("Season-wise summary"),
        
        tableOutput("season_summary_table"),
        
        uiOutput("validation_summary"),
        
        style = "
          width: 100%;
          margin-bottom: 12px;
        "
      ),
      
      # ======================================================
      # MAIN CHECKLIST RESULTS TABLE
      # ======================================================
      
      card(
        
        card_header("Checklist-level results"),
        
        downloadButton(
          "download_results",
          "Download results CSV"
        ),
        
        br(),
        br(),
        
        div(
          style = "
            width: 100%;
            overflow-x: auto;
            overflow-y: visible;
          ",
          
          DTOutput("results_table")
        ),
        
        style = "
          width: 100%;
          margin-bottom: 15px;
        "
      )
    )
  ),
  
  # ==========================================================
  # METHODOLOGY & USER GUIDE PAGE
  # ==========================================================
  
  nav_panel(
    
    "Methodology & User Guide",
    
    h2("AviCRED: Methodology & User Guide"),
    
    p(
      "AviCRED is a checklist-level framework for evaluating ",
      "the statistical consistency and atypicality of bird ",
      "observation checklists submitted through citizen-science ",
      "biodiversity databases."
    ),
    
    p(
      "The application integrates reference-based checklist ",
      "calibration, anomaly detection, sampling-effort diagnostics ",
      "and ecological diversity indices."
    ),
    
    tags$hr(),
    
    h3("1. Input Data and File Requirements"),
    
    h4("Supported input formats"),
    
    p(
      "AviCRED accepts checklist-level and species-level eBird ",
      "data files in CSV, TSV and TXT formats. Direct eBird URL ",
      "retrieval is available through the eBird API."
    ),
    
    h4("Checklist-level input"),
    
    p("One row per sampling event. Required information includes:"),
    
    tags$ul(
      tags$li("SAMPLING EVENT IDENTIFIER"),
      tags$li("OBSERVER ID"),
      tags$li("DURATION MINUTES"),
      tags$li("EFFORT DISTANCE KM"),
      tags$li("Observed richness"),
      tags$li("Observed abundance"),
      tags$li("STATE CODE"),
      tags$li("OBSERVATION DATE")
    ),
    
    p(
      "Season may alternatively be supplied directly if available."
    ),
    
    h4("Species-level input"),
    
    p(
      "Species-level eBird EBD files contain individual species ",
      "observations. AviCRED aggregates records by sampling event ",
      "to calculate checklist-level richness and abundance."
    ),
    
    tags$ul(
      tags$li("SAMPLING EVENT IDENTIFIER"),
      tags$li("OBSERVER ID"),
      tags$li("SCIENTIFIC NAME"),
      tags$li("OBSERVATION COUNT"),
      tags$li("DURATION MINUTES"),
      tags$li("EFFORT DISTANCE KM"),
      tags$li("STATE CODE"),
      tags$li("OBSERVATION DATE")
    ),
    
    p(
      "Species counts recorded as X are treated as a minimum of ",
      "one individual when calculating abundance and ecological indices."
    ),
    
    tags$hr(),
    
    h3("2. Checklist Atypicality Score (CAS)"),
    
    h4("Conceptual framework"),
    
    p(
      "CAS evaluates checklist-level statistical atypicality by ",
      "comparing observed richness and abundance against reference ",
      "expectations for the corresponding state-season combination."
    ),
    
    p(
      "The scoring framework uses two related measures: Proto-CAS ",
      "and CAS. Proto-CAS represents the averaged conformity score ",
      "derived from richness and abundance, whereas CAS is its ",
      "complement and represents checklist-level statistical atypicality."
    ),
    
    h4("Seasonal classification"),
    
    tags$table(
      class = "table table-bordered table-striped",
      
      tags$thead(
        tags$tr(
          tags$th("Season"),
          tags$th("Months")
        )
      ),
      
      tags$tbody(
        tags$tr(
          tags$td("Winter"),
          tags$td("November–February")
        ),
        
        tags$tr(
          tags$td("Pre-monsoon"),
          tags$td("March–May")
        ),
        
        tags$tr(
          tags$td("Monsoon"),
          tags$td("June–September")
        ),
        
        tags$tr(
          tags$td("Post-monsoon"),
          tags$td("October")
        )
      )
    ),
    
    h4("Expected values and residual calculation"),
    
    p(
      "Expected richness and abundance are retrieved from the ",
      "CAS expected-value reference file using the state-season ",
      "combination. The current scoring implementation uses ",
      "state-season expected means."
    ),
    
    p(
      "Richness residual = Observed richness − Expected richness"
    ),
    
    p(
      "Abundance residual = Observed abundance − Expected abundance"
    ),
    
    h4("Robust standardisation"),
    
    p(
      "Residuals are standardised using the corresponding ",
      "state-season reference residual median and median ",
      "absolute deviation (MAD)."
    ),
    
    p(
      "Robust Z = (Residual − Reference residual median) / MAD"
    ),
    
    h4("Proto-CAS and CAS calculation"),
    
    p(
      "The standardised deviations are converted into richness ",
      "and abundance conformity components using an exponential ",
      "decay function:"
    ),
    
    p(
      "Richness component = exp(−|Z richness| / 3.5)"
    ),
    
    p(
      "Abundance component = exp(−|Z abundance| / 3.5)"
    ),
    
    p(
      "Each component ranges from 0 to 1. Values closer to 1 ",
      "indicate closer agreement with the corresponding reference expectation."
    ),
    
    p(
      style = "font-weight:bold;text-align:center;",
      "Proto-CAS = (Richness component + Abundance component) / 2"
    ),
    
    p(
      "Proto-CAS represents the combined conformity score, ",
      "integrating the richness and abundance components with equal weight."
    ),
    
    p(
      style = "font-weight:bold;text-align:center;",
      "CAS = 1 − Proto-CAS"
    ),
    
    p(
      "CAS is the complement of Proto-CAS and represents statistical ",
      "atypicality. Higher Proto-CAS values indicate greater conformity ",
      "with reference expectations, whereas higher CAS values indicate ",
      "greater statistical atypicality."
    ),
    
    h4("CAS classification"),
    
    tags$table(
      class = "table table-bordered table-striped",
      
      tags$thead(
        tags$tr(
          tags$th("CAS range"),
          tags$th("Class"),
          tags$th("Interpretation")
        )
      ),
      
      tags$tbody(
        tags$tr(
          tags$td("0.00–<0.20"),
          tags$td("Very Low"),
          tags$td("Close agreement with reference expectations")
        ),
        
        tags$tr(
          tags$td("0.20–<0.40"),
          tags$td("Low"),
          tags$td("Relatively small deviation")
        ),
        
        tags$tr(
          tags$td("0.40–<0.60"),
          tags$td("Moderate"),
          tags$td("Intermediate deviation")
        ),
        
        tags$tr(
          tags$td("0.60–<0.80"),
          tags$td("High"),
          tags$td("Substantial deviation")
        ),
        
        tags$tr(
          tags$td("0.80–1.00"),
          tags$td("Very High"),
          tags$td("Strong deviation from reference expectations")
        )
      )
    ),
    
    tags$hr(),
    
    h3("3. Anomaly Detection"),
    
    p(
      "Richness and abundance anomalies are assessed independently ",
      "using their respective robust Z-scores."
    ),
    
    tags$ul(
      tags$li(
        tags$b("Richness anomaly: "),
        "absolute richness robust Z-score ≥ 3.5."
      ),
      
      tags$li(
        tags$b("Abundance anomaly: "),
        "absolute abundance robust Z-score ≥ 3.5."
      )
    ),
    
    p(
      "Both unusually high and unusually low values can be flagged. ",
      "A checklist may have a richness anomaly, an abundance ",
      "anomaly, both, or neither."
    ),
    
    p(
      "An anomaly indicates statistical atypicality and does not ",
      "automatically indicate an incorrect observation."
    ),
    
    tags$hr(),
    
    h3("4. Sampling Effort Diagnostics"),
    
    tags$ul(
      tags$li(
        tags$b("Duration effort flag: "),
        "duration ≥ 150 minutes."
      ),
      
      tags$li(
        tags$b("Distance effort flag: "),
        "Traveling checklist distance ≥ 4 km."
      )
    ),
    
    p(
      "Effort diagnostics are reported separately from CAS and ",
      "do not directly modify or penalise the score."
    ),
    
    tags$hr(),
    
    h3("5. Ecological Diversity Indices"),
    
    p(
      "For species-level input, AviCRED calculates additional ",
      "ecological indices from checklist-level species abundance."
    ),
    
    tags$ul(
      tags$li(
        tags$b("Shannon–Wiener diversity: "),
        "measures diversity considering richness and relative abundance."
      ),
      
      tags$li(
        tags$b("Simpson diversity (1 − D): "),
        "measures the probability that two randomly selected ",
        "individuals belong to different species."
      ),
      
      tags$li(
        tags$b("Pielou's evenness: "),
        "measures the distribution of individuals among species."
      ),
      
      tags$li(
        tags$b("Berger–Parker dominance: "),
        "measures the proportional abundance of the dominant species."
      ),
      
      tags$li(
        tags$b("Chao1: "),
        "estimates potential species richness using singleton ",
        "and doubleton frequencies."
      )
    ),
    
    p(
      "Abundance-dependent indices are unavailable when abundance ",
      "information is incomplete. These indices do not contribute ",
      "to CAS."
    ),
    
    tags$hr(),
    
    h3("6. CAS Reference Datasets"),
    
    tags$table(
      class = "table table-bordered table-striped",
      
      tags$thead(
        tags$tr(
          tags$th("Reference file"),
          tags$th("Purpose")
        )
      ),
      
      tags$tbody(
        tags$tr(
          tags$td("CAS_expected_value_variation.csv"),
          tags$td(
            "State-season expected richness and abundance summaries."
          )
        ),
        
        tags$tr(
          tags$td("CAS_protocol_expected_variation.csv"),
          tags$td(
            "Protocol-specific expected richness and abundance variation."
          )
        ),
        
        tags$tr(
          tags$td("CAS_reference_residual_statistics.csv"),
          tags$td(
            "State-season residual medians and MAD values."
          )
        ),
        
        tags$tr(
          tags$td("CAS_state_season_counts.csv"),
          tags$td(
            "Checklist sample sizes and reference coverage."
          )
        )
      )
    ),
    
    tags$hr(),
    
    h3("7. Interpretation and Limitations"),
    
    tags$ul(
      tags$li(
        "CAS measures statistical atypicality, not the probability ",
        "that a checklist contains erroneous data."
      ),
      
      tags$li(
        "Proto-CAS represents the averaged richness and abundance ",
        "conformity score, while CAS represents its complement."
      ),
      
      tags$li(
        "Unusual richness or abundance may reflect genuine ecological ",
        "events or unusual species assemblages."
      ),
      
      tags$li(
        "CAS depends on the availability and representativeness of ",
        "state-season reference values and residual statistics."
      ),
      
      tags$li(
        "Checklists without valid expected values or usable residual ",
        "dispersion estimates cannot be reliably scored."
      ),
      
      tags$li(
        "Incomplete abundance prevents complete two-component CAS calculation."
      ),
      
      tags$li(
        "Effort indicators are contextual information, not direct ",
        "evidence of checklist inauthenticity."
      )
    ),
    
    tags$hr(),
    
    h3("8. Responsible Use"),
    
    p(
      "AviCRED is an analytical support tool for evaluating ",
      "citizen-science bird-observation data. CAS values, anomaly ",
      "flags and effort diagnostics should be interpreted alongside ",
      "ecological knowledge, sampling context and independent ",
      "verification where necessary."
    ),
    
    p(
      style = "text-align:center;font-style:italic;",
      "AviCRED — Checklist-level assessment of citizen-science bird data."
    )
  )
)

# ============================================================
# eBIRD API 2.0 CHECKLIST RETRIEVAL
# Returns checklist-level data and attaches species-level data
# ============================================================

# ============================================================
# EXTRACT eBIRD CHECKLIST IDs
# ============================================================

extract_ebird_checklist_ids <- function(x) {
  
  if (is.null(x) || length(x) == 0) {
    return(character(0))
  }
  
  x <- as.character(x)
  
  # Combine all pasted lines
  x <- paste(x, collapse = "\n")
  
  # Split URLs / IDs separated by whitespace, commas,
  # semicolons, tabs, or line breaks
  tokens <- unlist(
    strsplit(
      x,
      "[\r\n,;[:space:]]+"
    )
  )
  
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  
  if (!length(tokens)) {
    return(character(0))
  }
  
  # eBird checklist IDs begin with S followed by digits
  ids <- vapply(
    tokens,
    function(z) {
      
      m <- regmatches(
        z,
        regexpr(
          "S[0-9]{5,}",
          z,
          perl = TRUE
        )
      )
      
      if (
        length(m) == 0 ||
        is.na(m) ||
        !nzchar(m)
      ) {
        NA_character_
      } else {
        m
      }
    },
    character(1)
  )
  
  ids <- ids[
    !is.na(ids) &
      nzchar(ids)
  ]
  
  unique(ids)
}

# ============================================================
# RETRIEVE ONE eBIRD CHECKLIST
# ============================================================

fetch_ebird_checklist <- function(
    sub_id,
    api_key
) {
  
  if (
    is.null(api_key) ||
    !nzchar(trimws(api_key))
  ) {
    stop(
      "Enter your personal eBird API key before retrieving checklist URLs."
    )
  }
  
  sub_id <- trimws(as.character(sub_id))
  api_key <- trimws(as.character(api_key))
  
  endpoint <- paste0(
    "https://api.ebird.org/v2/product/checklist/view/",
    sub_id
  )
  
  response <- httr::GET(
    endpoint,
    httr::add_headers(
      `X-eBirdApiToken` = api_key
    ),
    httr::user_agent(
      "AviCRED Shiny application"
    ),
    httr::timeout(45)
  )
  
  status <- httr::status_code(response)
  
  if (status != 200) {
    
    body <- tryCatch(
      httr::content(
        response,
        as = "text",
        encoding = "UTF-8"
      ),
      error = function(e) {
        ""
      }
    )
    
    stop(
      "eBird API returned HTTP ",
      status,
      " for ",
      sub_id,
      if (nzchar(body)) {
        paste0(
          ". ",
          substr(body, 1, 500)
        )
      } else {
        ""
      }
    )
  }
  
  txt <- httr::content(
    response,
    as = "text",
    encoding = "UTF-8"
  )
  
  jsonlite::fromJSON(
    txt,
    simplifyVector = FALSE
  )
}

fetch_ebird_checklists <- function(
    url_text,
    api_key,
    max_checklists = 50L
) {
  
  if (
    is.null(api_key) ||
    !nzchar(trimws(api_key))
  ) {
    stop(
      "Enter your personal eBird API key before retrieving checklist URLs."
    )
  }
  
  ids <- extract_ebird_checklist_ids(url_text)
  
  if (!length(ids)) {
    stop(
      "No valid eBird checklist IDs were found. ",
      "Paste URLs containing IDs such as S123456789."
    )
  }
  
  if (length(ids) > max_checklists) {
    stop(
      "A maximum of ", max_checklists,
      " checklists can be retrieved per run. You supplied ",
      length(ids), "."
    )
  }
  
  checklist_records <- vector("list", length(ids))
  species_records_list <- vector("list", length(ids))
  errors <- character(0)
  
  for (i in seq_along(ids)) {
    
    payload <- tryCatch(
      fetch_ebird_checklist(ids[[i]], api_key),
      error = function(e) {
        
        errors <<- c(
          errors,
          paste0(
            ids[[i]], ": ",
            conditionMessage(e)
          )
        )
        
        NULL
      }
    )
    
    if (is.null(payload)) next
    
    obs <- payload$obs
    
    if (is.null(obs) || length(obs) == 0) {
      
      errors <- c(
        errors,
        paste0(
          ids[[i]],
          ": checklist returned no observations"
        )
      )
      
      next
    }
    
    # ----------------------------------------------------------
    # EXTRACT CHECKLIST DATE AND START TIME FROM API obsDt
    # ----------------------------------------------------------
    
    obs_datetime <- if (
      is.null(payload$obsDt) ||
      length(payload$obsDt) == 0
    ) {
      NA_character_
    } else {
      trimws(as.character(payload$obsDt))
    }
    
    # Extract date separately
    obs_date_value <- NA_character_
    
    if (
      !is.na(obs_datetime) &&
      nzchar(obs_datetime)
    ) {
      date_match <- regmatches(
        obs_datetime,
        regexpr(
          "^\\d{4}-\\d{2}-\\d{2}",
          obs_datetime
        )
      )
      
      if (length(date_match) > 0 && nzchar(date_match)) {
        obs_date_value <- date_match
      }
    }
    
    # Extract start time only when eBird marks it as valid
    start_time_value <- NA_character_
    
    if (
      isTRUE(payload$obsTimeValid) &&
      !is.na(obs_datetime) &&
      nzchar(obs_datetime)
    ) {
      
      time_part <- sub(
        "^\\d{4}-\\d{2}-\\d{2}[T ]",
        "",
        obs_datetime
      )
      
      time_match <- regmatches(
        time_part,
        regexpr(
          "^\\d{1,2}:\\d{2}(:\\d{2})?",
          time_part
        )
      )
      
      if (length(time_match) > 0 && nzchar(time_match)) {
        
        parsed_time <- suppressWarnings(
          as.POSIXct(
            time_match,
            format = if (nchar(time_match) > 5) {
              "%H:%M:%S"
            } else {
              "%H:%M"
            },
            tz = "UTC"
          )
        )
        
        if (!is.na(parsed_time)) {
          start_time_value <- format(
            parsed_time,
            "%H:%M",
            tz = "UTC"
          )
        }
      }
    }
    
    # --------------------------------------------------------
    # SPECIES NAMES
    # --------------------------------------------------------
    
    species_names <- vapply(obs, function(z) {
      
      candidates <- c(
        z$sciName %||% NA_character_,
        z$comName %||% NA_character_,
        z$speciesCode %||% NA_character_
      )
      
      candidates <- trimws(as.character(candidates))
      
      valid <- candidates[
        !is.na(candidates) & candidates != ""
      ]
      
      if (length(valid) == 0) {
        NA_character_
      } else {
        valid[1]
      }
      
    }, character(1))
    
    species_codes <- vapply(obs, function(z) {
      
      if (is.null(z$speciesCode)) {
        NA_character_
      } else {
        as.character(z$speciesCode)
      }
      
    }, character(1))
    
    # --------------------------------------------------------
    # SPECIES COUNTS
    # X = minimum of one individual
    # --------------------------------------------------------
    
    count_strings <- vapply(obs, function(z) {
      
      value <- z$howManyStr
      
      if (
        is.null(value) ||
        length(value) == 0 ||
        is.na(value) ||
        trimws(as.character(value)) == ""
      ) {
        value <- z$howMany
      }
      
      if (
        is.null(value) ||
        length(value) == 0 ||
        is.na(value) ||
        trimws(as.character(value)) == ""
      ) {
        value <- "X"
      }
      
      as.character(value)
      
    }, character(1))
    
    count_numeric <- vapply(count_strings, function(x) {
      
      x <- trimws(as.character(x))
      
      if (toupper(x) == "X") {
        return(1)
      }
      
      value <- suppressWarnings(as.numeric(x))
      
      if (is.na(value) || value < 0) {
        return(NA_real_)
      }
      
      value
      
    }, numeric(1))
    
    # --------------------------------------------------------
    # SPECIES-LEVEL DATA FOR THIS CHECKLIST
    # --------------------------------------------------------
    
    species_raw <- data.frame(
      checklist_id = as.character(
        payload$subId %||% ids[[i]]
      ),
      species_name = species_names,
      species_code = species_codes,
      count = count_strings,
      count_numeric = count_numeric,
      stringsAsFactors = FALSE
    )
    
    species_raw <- species_raw %>%
      dplyr::filter(
        !is.na(species_name),
        trimws(species_name) != ""
      )
    
    species_by_checklist <- species_raw %>%
      dplyr::group_by(
        checklist_id,
        species_name
      ) %>%
      dplyr::summarise(
        
        species_abundance = if (
          all(is.na(count_numeric))
        ) {
          NA_real_
        } else {
          sum(count_numeric, na.rm = TRUE)
        },
        
        abundance_incomplete = if (
          all(is.na(count_numeric))
        ) {
          "Yes"
        } else {
          "No"
        },
        
        .groups = "drop"
      )
    
    species_records_list[[i]] <- species_by_checklist
    
    # --------------------------------------------------------
    # CHECKLIST RICHNESS AND TOTAL ABUNDANCE
    # --------------------------------------------------------
    
    observed_richness <- length(
      unique(stats::na.omit(species_codes))
    )
    
    observed_abundance <- if (
      any(is.na(count_numeric))
    ) {
      NA_real_
    } else {
      sum(count_numeric)
    }
    
    # --------------------------------------------------------
    # HIGHEST-ABUNDANCE SPECIES
    # --------------------------------------------------------
    
    valid_species <- species_by_checklist %>%
      dplyr::filter(
        !is.na(species_abundance)
      )
    
    if (nrow(valid_species) == 0) {
      
      highest_species_abundance <- NA_real_
      highest_abundance_species <- NA_character_
      
    } else {
      
      highest_species_abundance <- max(
        valid_species$species_abundance
      )
      
      highest_abundance_species <- paste(
        valid_species$species_name[
          valid_species$species_abundance ==
            highest_species_abundance
        ],
        collapse = "; "
      )
    }
    
    # --------------------------------------------------------
    # CHECKLIST METADATA
    # --------------------------------------------------------
    
    protocol_code <- as.character(
      payload$protocolId %||% NA_character_
    )
    
    protocol_name <- if (
      identical(protocol_code, "P21")
    ) {
      "Traveling"
    } else if (
      identical(protocol_code, "P22")
    ) {
      "Stationary"
    } else {
      protocol_code
    }
    
    duration_min <- if (
      is.null(payload$durationHrs)
    ) {
      NA_real_
    } else {
      as.numeric(payload$durationHrs) * 60
    }
    
    distance_km <- if (
      is.null(payload$effortDistanceKm)
    ) {
      NA_real_
    } else {
      as.numeric(payload$effortDistanceKm)
    }
    
    state_code <- if (
      is.null(payload$subnational1Code)
    ) {
      NA_character_
    } else {
      toupper(as.character(payload$subnational1Code))
    }
    
    display_name <- if (
      is.null(payload$userDisplayName)
    ) {
      NA_character_
    } else {
      as.character(payload$userDisplayName)
    }
    
    observer_proxy <- if (
      is.na(display_name) ||
      !nzchar(display_name)
    ) {
      paste0("API_USER_", ids[[i]])
    } else {
      paste0("API_USER_", display_name)
    }
    
    # --------------------------------------------------------
    # CHECKLIST-LEVEL RECORD
    # --------------------------------------------------------
    
    checklist_records[[i]] <- data.frame(
      
      checklist_id = as.character(
        payload$subId %||% ids[[i]]
      ),
      
      observer_id = observer_proxy,
      
      observer_name = display_name,
      
      obs_date = obs_date_value,
      
      start_time = start_time_value,
      
      state_code = state_code,
      
      protocol = protocol_name,
      
      duration_min = duration_min,
      
      distance_km = distance_km,
      
      complete_flag = isTRUE(
        payload$allObsReported
      ),
      
      observed_richness = observed_richness,
      
      observed_abundance = observed_abundance,
      
      highest_species_abundance =
        highest_species_abundance,
      
      highest_abundance_species =
        highest_abundance_species,
      
      api_species_codes = paste(
        stats::na.omit(species_codes),
        collapse = ";"
      ),
      
      api_species_counts = paste(
        count_strings,
        collapse = ";"
      ),
      
      stringsAsFactors = FALSE
    )
    
    Sys.sleep(0.15)
  }
  
  # ----------------------------------------------------------
  # COMBINE CHECKLIST RECORDS
  # ----------------------------------------------------------
  
  out <- dplyr::bind_rows(
    checklist_records
  )
  
  if (!nrow(out)) {
    stop(
      "No checklist could be retrieved. ",
      paste(errors, collapse = " | ")
    )
  }
  
  # ----------------------------------------------------------
  # COMBINE SPECIES RECORDS
  # ----------------------------------------------------------
  
  species_out <- dplyr::bind_rows(
    species_records_list
  )
  
  # Preserve both outputs without changing the existing
  # data-frame return structure expected by the server.
  attr(out, "api_errors") <- errors
  attr(out, "api_species_records") <- species_out
  
  return(out)
}

# ============================================================
# 8. SERVER
# ============================================================

server <- function(input, output, session) {
  
  results <- reactiveVal(NULL)
  
  species_records <- reactiveVal(NULL)
  
  # Keep the rest of your existing server code below.
  
  # ----------------------------------------------------------
  # LOCATION AND SEASON VALIDATION SUMMARY
  # ----------------------------------------------------------
  
  output$validation_summary <- renderUI({
    
    req(results())
    
    data <- results()
    
    total <- nrow(data)
    
    valid_state <- sum(
      data$state_validation_status == "Valid state code",
      na.rm = TRUE
    )
    
    valid_dates <- sum(
      data$date_validation_status == "Valid observation date",
      na.rm = TRUE
    )
    
    valid_localities <- sum(
      data$locality_validation_status == "Locality ID supplied",
      na.rm = TRUE
    )
    
    matched_reference <- sum(
      data$reference_match_status == "Matched to CAS reference",
      na.rm = TRUE
    )
    
    scored <- sum(
      data$CAS_scoring_status == "Scored",
      na.rm = TRUE
    )
    
    tagList(
      
      p(paste(
        "Total checklists:",
        format(total, big.mark = ",")
      )),
      
      p(paste(
        "Valid state codes:",
        format(valid_state, big.mark = ",")
      )),
      
      p(paste(
        "Locality IDs supplied:",
        format(valid_localities, big.mark = ",")
      )),
      
      p(paste(
        "Valid observation dates:",
        format(valid_dates, big.mark = ",")
      )),
      
      p(paste(
        "Matched state-season references:",
        format(matched_reference, big.mark = ",")
      )),
      
      p(paste(
        "Successfully scored:",
        format(scored, big.mark = ",")
      ))
      
    )
    
  })
  
  # ==========================================================
  # OVERALL DATASET SUMMARY
  # ==========================================================
  
  output$overall_summary_table <- renderTable({
    
    req(results())
    
    dat <- results()
    
    get_column <- function(data, candidates) {
      
      matched <- candidates[candidates %in% names(data)]
      
      if (length(matched) == 0) {
        return(NULL)
      }
      
      matched[1]
    }
    
    richness_col <- get_column(
      dat,
      c("observed_richness", "richness")
    )
    
    abundance_col <- get_column(
      dat,
      c("observed_abundance", "abundance")
    )
    
    duration_col <- get_column(
      dat,
      c("duration_min", "duration")
    )
    
    distance_col <- get_column(
      dat,
      c("distance_km", "distance")
    )
    
    checklist_col <- get_column(
      dat,
      c("checklist_id")
    )
    
    observer_col <- get_column(
      dat,
      c("observer_id")
    )
    
    total_checklists <- if (!is.null(checklist_col)) {
      dplyr::n_distinct(dat[[checklist_col]], na.rm = TRUE)
    } else {
      nrow(dat)
    }
    
    total_observers <- if (!is.null(observer_col)) {
      dplyr::n_distinct(dat[[observer_col]], na.rm = TRUE)
    } else {
      NA_integer_
    }
    
    mean_richness <- if (!is.null(richness_col)) {
      mean(dat[[richness_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    
    median_richness <- if (!is.null(richness_col)) {
      median(dat[[richness_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    
    mean_abundance <- if (!is.null(abundance_col)) {
      mean(dat[[abundance_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    
    median_abundance <- if (!is.null(abundance_col)) {
      median(dat[[abundance_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    
    mean_duration <- if (!is.null(duration_col)) {
      mean(dat[[duration_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    
    mean_distance <- if (!is.null(distance_col)) {
      mean(dat[[distance_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    
    data.frame(
      Metric = c(
        "Total checklists",
        "Unique observers",
        "Mean species richness",
        "Median species richness",
        "Mean total abundance",
        "Median total abundance",
        "Mean checklist duration (minutes)",
        "Mean distance travelled (km)"
      ),
      
      Value = c(
        total_checklists,
        total_observers,
        mean_richness,
        median_richness,
        mean_abundance,
        median_abundance,
        mean_duration,
        mean_distance
      )
    )
    
  }, digits = 3, striped = TRUE, bordered = TRUE)
  
  
  # ==========================================================
  # CHECKLIST DROPDOWN
  # ==========================================================
  
  output$checklist_selector <- renderUI({
    
    req(results())
    
    dat <- results()
    
    req("checklist_id" %in% names(dat))
    
    checklist_ids <- sort(
      unique(
        as.character(
          dat$checklist_id[
            !is.na(dat$checklist_id) &
              dat$checklist_id != ""
          ]
        )
      )
    )
    
    selectInput(
      inputId = "selected_checklist",
      label = "Select checklist ID",
      choices = checklist_ids,
      selected = if (length(checklist_ids) > 0) {
        checklist_ids[1]
      } else {
        NULL
      },
      width = "100%"
    )
    
  })
  
  # ==========================================================
  # INDIVIDUAL CHECKLIST DETAILS
  # ==========================================================
  
  output$individual_checklist_table <- renderTable({
    
    req(results())
    req(input$selected_checklist)
    
    dat <- results()
    
    selected <- dat %>%
      dplyr::filter(
        as.character(checklist_id) == input$selected_checklist
      )
    
    req(nrow(selected) > 0)
    
    selected <- selected[1, , drop = FALSE]
    
    # --------------------------------------------------------
    # HELPER: EXTRACT CHECKLIST VALUES
    # --------------------------------------------------------
    
    get_value <- function(data, candidates) {
      
      matched <- candidates[candidates %in% names(data)]
      
      if (length(matched) == 0) {
        return("Not available")
      }
      
      value <- data[[matched[1]]][1]
      
      if (length(value) == 0 || is.na(value)) {
        return("Not available")
      }
      
      value <- trimws(as.character(value))
      
      if (value == "") {
        return("Not available")
      }
      
      value
    }
    
    # --------------------------------------------------------
    # HELPER: FORMAT START TIME AS 24-HOUR TIME
    # --------------------------------------------------------
    
    format_24hr_time <- function(x) {
      
      if (
        is.null(x) ||
        length(x) == 0 ||
        is.na(x) ||
        trimws(as.character(x)) == ""
      ) {
        return("Not available")
      }
      
      # Convert input to numeric seconds
      seconds <- suppressWarnings(
        as.numeric(as.character(x))
      )
      
      # Handle eBird start times stored as seconds since midnight
      if (!is.na(seconds)) {
        
        if (seconds < 0 || seconds >= 86400) {
          return("Not available")
        }
        
        hours <- floor(seconds / 3600)
        minutes <- floor((seconds %% 3600) / 60)
        
        return(
          sprintf("%02d:%02d", hours, minutes)
        )
      }
      
      # Handle time values already stored as text
      x <- trimws(as.character(x))
      
      parsed_time <- suppressWarnings(
        as.POSIXct(
          x,
          format = "%I:%M:%S %p",
          tz = "UTC"
        )
      )
      
      if (is.na(parsed_time)) {
        parsed_time <- suppressWarnings(
          as.POSIXct(
            x,
            format = "%I:%M %p",
            tz = "UTC"
          )
        )
      }
      
      if (is.na(parsed_time)) {
        parsed_time <- suppressWarnings(
          as.POSIXct(
            x,
            format = "%H:%M:%S",
            tz = "UTC"
          )
        )
      }
      
      if (is.na(parsed_time)) {
        parsed_time <- suppressWarnings(
          as.POSIXct(
            x,
            format = "%H:%M",
            tz = "UTC"
          )
        )
      }
      
      if (is.na(parsed_time)) {
        return(x)
      }
      
      format(parsed_time, "%H:%M", tz = "UTC")
    }
    
    # --------------------------------------------------------
    # BASIC CHECKLIST DETAILS
    # --------------------------------------------------------
    
    checklist_id <- get_value(
      selected,
      c("checklist_id")
    )
    
    observation_date <- get_value(
      selected,
      c("obs_date", "observation_date")
    )
    
    start_time_raw <- get_value(
      selected,
      c("start_time", "start_time_raw", "start_hour")
    )
    
    start_time <- format_24hr_time(
      start_time_raw
    )
    
    # --------------------------------------------------------
    # STATE AND SEASON
    # --------------------------------------------------------
    
    state <- get_value(
      selected,
      c("state_code", "state", "state_name")
    )
    
    season <- get_value(
      selected,
      c("season")
    )
    
    # --------------------------------------------------------
    # SAMPLING DETAILS
    # --------------------------------------------------------
    
    protocol <- get_value(
      selected,
      c("protocol", "protocol_name")
    )
    
    richness <- get_value(
      selected,
      c("observed_richness", "richness")
    )
    
    abundance <- get_value(
      selected,
      c("observed_abundance", "abundance")
    )
    
    distance <- get_value(
      selected,
      c("distance_km", "distance")
    )
    
    duration <- get_value(
      selected,
      c("duration_min", "duration")
    )
    
    observers <- get_value(
      selected,
      c("n_observers", "number_observers", "observer_count")
    )
    
    # --------------------------------------------------------
    # OBSERVER COUNT FALLBACK
    # --------------------------------------------------------
    
    if (
      observers == "Not available" &&
      "observer_id" %in% names(selected)
    ) {
      observers <- "1"
    }
    
    # --------------------------------------------------------
    # HIGHEST SPECIES ABUNDANCE
    # --------------------------------------------------------
    
    highest_species <- get_value(
      selected,
      c("highest_abundance_species")
    )
    
    highest_species_count <- get_value(
      selected,
      c("highest_species_abundance")
    )
    
    # Use the checklist-level result first.
    # This is already calculated by process_species_data()
    # for species-level CSV/EBD input.
    
    if (highest_species != "Not available") {
      
      if (highest_species_count != "Not available") {
        
        highest_species <- paste0(
          highest_species,
          " (",
          highest_species_count,
          " individuals)"
        )
        
      }
      
    } else {
      
      # Fallback for API or other species-level records
      
      highest_species <- "Not available"
      
      if (
        !is.null(species_records()) &&
        nrow(species_records()) > 0
      ) {
        
        species_dat <- species_records()
        
        species_name_col <- c(
          "species_name",
          "common_name",
          "COMMON NAME"
        )
        
        abundance_col <- c(
          "species_abundance",
          "count_numeric",
          "count",
          "species_count",
          "abundance",
          "observation_count",
          "OBSERVATION COUNT"
        )
        
        species_name_col <- species_name_col[
          species_name_col %in% names(species_dat)
        ]
        
        abundance_col <- abundance_col[
          abundance_col %in% names(species_dat)
        ]
        
        if (
          length(species_name_col) > 0 &&
          length(abundance_col) > 0 &&
          "checklist_id" %in% names(species_dat)
        ) {
          
          selected_species <- species_dat %>%
            dplyr::filter(
              as.character(checklist_id) ==
                input$selected_checklist
            ) %>%
            dplyr::mutate(
              .abundance_numeric = suppressWarnings(
                as.numeric(
                  as.character(
                    .data[[abundance_col[1]]]
                  )
                )
              )
            ) %>%
            dplyr::filter(
              !is.na(.abundance_numeric)
            )
          
          if (nrow(selected_species) > 0) {
            
            highest_value <- max(
              selected_species$.abundance_numeric
            )
            
            highest_names <- selected_species %>%
              dplyr::filter(
                .abundance_numeric == highest_value
              ) %>%
              dplyr::pull(
                all_of(species_name_col[1])
              ) %>%
              unique()
            
            highest_species <- paste0(
              paste(highest_names, collapse = "; "),
              " (",
              highest_value,
              " individuals)"
            )
          }
        }
      }
    }
    
    # --------------------------------------------------------
    # OUTPUT TABLE
    # --------------------------------------------------------
    
    data.frame(
      
      Metric = c(
        "Checklist ID",
        "Observation date",
        "Start time (24-hour)",
        "State",
        "Season",
        "Protocol",
        "Species richness",
        "Total abundance",
        "Highest species abundance",
        "Distance travelled (km)",
        "Duration (minutes)",
        "Number of observers"
      ),
      
      Value = c(
        checklist_id,
        observation_date,
        start_time,
        state,
        season,
        protocol,
        richness,
        abundance,
        highest_species,
        distance,
        duration,
        observers
      ),
      
      stringsAsFactors = FALSE
    )
    
  }, digits = 3, striped = TRUE, bordered = TRUE)
  
  # ============================================================
  # SEASON SUMMARY TABLE
  # ============================================================
  
  output$season_summary_table <- renderTable({
    
    req(results())
    
    results() %>%
      count(season, name = "Number of checklists")
    
  })
  
  # ----------------------------------------------------------
  # NAVIGATION
  # ----------------------------------------------------------
  
  observeEvent(input$start_assessment, {
    
    updateNavbarPage(
      session,
      "main_nav",
      selected = "Assess data"
    )
  })
  
  
  # ----------------------------------------------------------
  # INPUT INSTRUCTIONS
  # ----------------------------------------------------------
  
  output$input_instructions <- renderUI({
    
    if (input$input_method == "checklist") {
      
      tagList(
        
        h4("Checklist-level file"),
        
        p("One row per checklist."),
        
        p(
          paste(
            "Required columns: SAMPLING EVENT IDENTIFIER,",
            "OBSERVER ID, DURATION MINUTES, EFFORT DISTANCE KM,",
            "observed richness and observed abundance."
          )
        ),
        
        p(
          paste(
            "For CAS scoring, STATE CODE and OBSERVATION DATE",
            "are also required. Season can alternatively be",
            "provided directly."
          )
        ),
        
        p(
          paste(
            "CAS expected values and residual statistics are",
            "automatically retrieved from the reference CSV files."
          )
        )
      )
      
    } else if (input$input_method == "species") {
      
      tagList(
        
        h4("Species-level eBird EBD"),
        
        p("One row per species observation."),
        
        p(
          paste(
            "Required eBird columns: SAMPLING EVENT IDENTIFIER,",
            "OBSERVER ID, SCIENTIFIC NAME and OBSERVATION COUNT."
          )
        ),
        
        p(
          paste(
            "Required checklist metadata: DURATION MINUTES,",
            "EFFORT DISTANCE KM, STATE CODE and OBSERVATION DATE."
          )
        ),
        
        p(
          paste(
            "The application aggregates species records into",
            "checklist-level richness and abundance."
          )
        ),
        
        p(
          paste(
            "If any species count is recorded as X, total",
            "abundance is marked '1' and CAS is scored accordingly."
          )
        )
      )
      
    } else {
      
      tagList(
        
        h4("eBird checklist URLs"),
        
        p(
          paste(
            "URL retrieval requires official eBird API integration,",
            "which is not active in this prototype."
          )
        )
      )
    }
  })
  
  # ----------------------------------------------------------
  # PROCESS UPLOADED DATA OR RETRIEVE EBIRD CHECKLISTS
  # ----------------------------------------------------------
  
  observeEvent(input$run_assessment, {
    
    output$processing_status <- renderText(
      "Preparing assessment..."
    )
    
    tryCatch({
      
      api_errors <- character(0)
      
      # ======================================================
      # A. EBIRD URL INPUT
      # ======================================================
      
      if (identical(input$input_method, "url")) {
        
        # Validate API key
        if (
          is.null(input$ebird_api_key) ||
          !nzchar(trimws(input$ebird_api_key))
        ) {
          stop("Please enter your eBird API key.")
        }
        
        # Validate checklist URLs
        if (
          is.null(input$checklist_urls) ||
          !nzchar(trimws(input$checklist_urls))
        ) {
          stop("Please paste at least one eBird checklist URL.")
        }
        
        output$processing_status <- renderText(
          "Retrieving checklists from eBird..."
        )
        
        # Retrieve checklist data using the existing API function
        checklist_data <- fetch_ebird_checklists(
          url_text = input$checklist_urls,
          api_key = trimws(input$ebird_api_key),
          max_checklists = 50L
        )
        
        # Retain any retrieval warnings
        api_errors <- attr(
          checklist_data,
          "api_errors"
        ) %||% character(0)
        
        species_data <- attr(
          checklist_data,
          "api_species_records"
        )
        
        # Standard checklist processing
        checklist_data <- process_checklist_data(
          checklist_data
        )
        
      } else {
        
        # ====================================================
        # B. CSV / EBD INPUT
        # ====================================================
        
        if (is.null(input$csv_file)) {
          stop("Please upload a CSV, TSV or TXT file.")
        }
        
        output$processing_status <- renderText(
          "Reading uploaded data..."
        )
        
        raw_data <- read_ebird_file(
          input$csv_file$datapath,
          input$csv_file$name
        )
        
        # -----------------------------------------------
        # CHECKLIST-LEVEL INPUT
        # -----------------------------------------------
        
        if (identical(input$input_method, "checklist")) {
          
          checklist_data <- process_checklist_data(
            raw_data
          )
          
          species_data <- NULL
          
          # -----------------------------------------------
          # SPECIES-LEVEL INPUT
          # -----------------------------------------------
          
        } else if (identical(input$input_method, "species")) {
          
          processed_species <- process_species_data(
            raw_data
          )
          
          checklist_data <- processed_species$checklists
          species_data <- processed_species$species
          
        } else {
          
          stop("Invalid input method selected.")
          
        }
      }
      
      # ======================================================
      # C. STORE PROCESSED DATA
      # ======================================================
      
      # ======================================================
      # D. CAS IS ALREADY CALCULATED ONCE BY THE PROCESSOR
      # ======================================================

      results(checklist_data)
      species_records(species_data)
      
      assessed_n <- if ("CAS_scoring_status" %in% names(checklist_data)) {
        sum(checklist_data$CAS_scoring_status == "Scored", na.rm = TRUE)
      } else {
        0L
      }
      
      cas_unavailable <- if ("CAS_scoring_status" %in% names(checklist_data)) {
        any(grepl("^CAS unavailable:", checklist_data$CAS_scoring_status))
      } else {
        FALSE
      }
      
      source_text <- if (identical(input$input_method, "url")) {
        "Retrieved from eBird API and processed"
      } else {
        "Uploaded data processed"
      }
      
      status_text <- paste0(
        source_text, ": ",
        format(nrow(checklist_data), big.mark = ","),
        " checklists. ",
        if (cas_unavailable) {
          "CAS engine/reference scoring was unavailable; checklist rows and effort diagnostics were retained."
        } else {
          paste0("Fixed-reference CAS calculated for ",
                 format(assessed_n, big.mark = ","), " checklists.")
        }
      )
      
      # ======================================================
      # F. REPORT API RETRIEVAL WARNINGS
      # ======================================================
      
      if (length(api_errors) > 0) {
        
        status_text <- paste0(
          status_text,
          " Some checklists could not be retrieved: ",
          paste(api_errors, collapse = " | ")
        )
      }
      
      # ======================================================
      # G. DISPLAY SUCCESS AND OPEN RESULTS
      # ======================================================
      
      output$processing_status <- renderText(
        status_text
      )
      
      bslib::nav_select(
        id = "main_nav",
        selected = "Results",
        session = session
      )
      
    }, error = function(e) {
      
      error_text <- conditionMessage(e)
      
      output$processing_status <- renderText(
        paste("Assessment error:", error_text)
      )
      
      showNotification(
        paste("Assessment error:", error_text),
        type = "error",
        duration = NULL
      )
    })
  })
  
  # ----------------------------------------------------------
  # CLEAN AND FORMAT RESULTS
  # Used by BOTH the table and CSV download
  # ----------------------------------------------------------
  formatted_results <- reactive({
    
    req(results())
    
    data <- results()
    
    # --------------------------------------------------
    # REMOVE UNWANTED COLUMNS AND JOIN DUPLICATES
    # --------------------------------------------------
    
    remove_cols <- c(
      "state_code",
      "state_code.x",
      "state_code.y",
      "season",
      "season.x",
      "season.y",
      "abundance_incomplete",
      "effort_flag",
      "combined_effort_flag",
      "anomaly_flag",
      "either_anomaly_flag",
      "any_anomaly_flag",
      "CAS_scoring_status",
      "cas_scoring_status"
    )
    
    data <- data %>%
      dplyr::select(
        -dplyr::any_of(remove_cols)
      )
    
    # --------------------------------------------------
    # FORMAT START TIME TO 24-HOUR HH:MM
    # --------------------------------------------------
    
    if ("start_time" %in% names(data)) {
      
      data$start_time <- as.character(data$start_time)
      
      data$start_time <- vapply(data$start_time, function(x) {
        
        if (is.na(x) || trimws(x) == "") {
          return(NA_character_)
        }
        
        x <- trimws(x)
        
        formats <- c(
          "%H:%M:%S",
          "%H:%M",
          "%I:%M:%S %p",
          "%I:%M %p"
        )
        
        parsed <- NULL
        
        for (fmt in formats) {
          
          temp <- suppressWarnings(
            strptime(x, format = fmt, tz = "UTC")
          )
          
          if (!is.na(temp)) {
            parsed <- temp
            break
          }
        }
        
        if (is.null(parsed)) {
          return(x)
        }
        
        format(parsed, "%H:%M")
        
      }, character(1))
    }
    
    # --------------------------------------------------
    # ROUND NUMERIC COLUMNS TO THREE DECIMAL PLACES
    # --------------------------------------------------
    
    numeric_cols <- names(data)[
      vapply(data, is.numeric, logical(1))
    ]
    
    data[numeric_cols] <- lapply(
      data[numeric_cols],
      function(x) round(x, 3)
    )
    
    data
  })
  
  # ----------------------------------------------------------
  # RESULTS TABLE
  # ----------------------------------------------------------
  
  output$results_table <- renderDT({
    
    req(results())
    
    results_data <- results()
    
    numeric_columns <- names(results_data)[
      vapply(results_data, is.numeric, logical(1))
    ]
    
    datatable(
      
      results_data,
      
      rownames = FALSE,
      
      extensions = c("Buttons"),
      
      options = list(
        
        dom = "Bfrtip",
        
        scrollX = TRUE,
        
        pageLength = 25,
        
        lengthMenu = c(10, 25, 50, 100),
        
        autoWidth = TRUE,
        
        fixedHeader = FALSE,
        
        buttons = c("copy", "csv", "excel"),
        
        deferRender = TRUE
        
      )
      
    ) %>%
      
      formatRound(
        columns = numeric_columns,
        digits = 3
      )
    
  }, server = TRUE)
  
  # ----------------------------------------------------------
  # DOWNLOAD CLEANED RESULTS AS CSV
  # ----------------------------------------------------------
  
  output$download_results <- downloadHandler(
    
    filename = function() {
      
      paste0(
        "AviCRED_results_",
        Sys.Date(),
        ".csv"
      )
    },
    
    content = function(file) {
      
      data <- formatted_results()
      
      readr::write_csv(
        data,
        file,
        na = ""
      )
    }
  )
  
}

# ============================================================
# FINAL STRUCTURAL CHECKS
# ============================================================

# 1. Reference key uniqueness

stopifnot(
  !anyDuplicated(cas_reference$state_season)
)

# 2. Required reference columns

stopifnot(
  all(
    c(
      "state_season",
      "expected_richness",
      "expected_abundance",
      "richness_residual_median",
      "richness_residual_MAD",
      "abundance_residual_median",
      "abundance_residual_MAD"
    ) %in% names(cas_reference)
  )
)

# 3. Required functions

stopifnot(
  is.function(safe_numeric),
  is.function(standardise_names),
  is.function(validate_location_season),
  is.function(lookup_state_season_key),
  is.function(calculate_CAS),
  is.function(cas_fallback),
  is.function(process_checklist_data),
  is.function(process_species_data),
  is.function(read_ebird_file)
)

cat("\nAviCRED functions and CAS references loaded successfully.\n")

# ============================================================
# 9. RUN APPLICATION
# ============================================================

shinyApp(
  ui = ui,
  server = server
)
