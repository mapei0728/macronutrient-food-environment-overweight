# ============================================================
# Univariate Bayesian GLMM 
# ============================================================

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidybayes)
  library(car)
  library(posterior)
})

# -----------------------------
# 0) CONFIGURATION
# -----------------------------
DATA_PATH   <- "Physical/Data20260116.csv"
OUT_DIR     <- "results_univariate_bayes"
MODEL_DIR   <- file.path(OUT_DIR, "models_rds")

DIR.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DIR.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

SEED        <- 20260130
CHAINS      <- 4
ITER        <- 4000
WARMUP      <- 2000
CORES       <- 4
ADAPT_DELTA <- 0.95

RESPONSE_VAR  <- "OW"
RANDOM_EFFECT <- "SpatialZone"   # will be created from MinizoneID

SAVE_MODELS_RDS <- TRUE          # set FALSE if you do not want to save brmsfit objects

set.seed(SEED)

# -----------------------------
# 1) LOAD DATA
# -----------------------------
FD <- read.csv(DATA_PATH)

# -----------------------------
# 2) PREPROCESSING
# -----------------------------
FD <- FD %>%
  mutate(
    # random effect grouping
    SpatialZone = as.factor(MinizoneID),
    
    Gender        = as.factor(Gender),
    Ethnicity     = as.factor(Ethnicity),
    SmokingStatus = as.factor(SmokingStatus),
    HBP           = as.factor(HBP),
    HighEdu       = as.factor(HighEdu),
    
    # Collapse categories (adjust rules if needed)
    MaritalStatus = car::recode(as.factor(MaritalStatus), "c(3,4,5)=3"),
    WorkStatus    = car::recode(as.factor(WorkStatus),  "c(5,6,7)=5"),
    
    # AgeGroup: set reference level only if "65+" exists
    AgeGroup = {
      ag <- as.factor(AgeGroup)
      if ("65+" %in% levels(ag)) relevel(ag, ref = "65+") else ag
    },
    
    # Ensure OW is numeric 0/1
    OW = as.integer(as.character(OW))
  ) %>%
  mutate(
    # Example Green binary factor derived from GreenA_P
    Green = as.factor(dplyr::case_when(
      GreenA_P < 0.19 ~ 0,
      TRUE            ~ 1
    ))
  )

# Defensive check for outcome
if (!all(FD[[RESPONSE_VAR]] %in% c(0, 1, NA))) {
  stop("OW must be coded as 0/1 (or NA). Please check the original OW column.")
}

# -----------------------------
# 3) DEFINE PREDICTORS 
# -----------------------------
predictors <- c(
  "Gender","AgeGroup","Ethnicity","HighEdu","EduYr","MaritalStatus","WorkStatus",
  "SmokingStatus","HBP","TotalPAL",
  "AverageCaloriesPerMeal","AverageCarbohydratesPerMeal","AverageFatPerMeal",
  "AverageSaturatedFatPerMeal","AverageFibrePerMeal","AverageProteinPerMeal",
  "NutritionalValueScore","Sat11p","IntersectionDens","NLI","BusDens","MRTDens",
  "TaxiDens","Walkability","MrtWalk","HighW_P","ComA_P","ResA_P","GreenA_P","Green",
  "Fre_fast","Div_SSB","Diversity","Div_P","Div_SSB_P","Sat11p_Level",
  "Mini_D","ResaleP_SES","PD_2020","cluster_raw","Postal_N","FoodT_N","FoodO_N"
)

# Filter out predictors not present in the dataset
missing_vars <- setdiff(predictors, names(FD))
if (length(missing_vars) > 0) {
  message("Predictors NOT found in the dataset and will be skipped:\n- ",
          paste(missing_vars, collapse = ", "))
}
predictors <- intersect(predictors, names(FD))

if (length(predictors) == 0) stop("No valid predictors found in the dataset.")

# -----------------------------
# 4) MODEL FITTING FUNCTION
# -----------------------------
fit_univariate_bayes_glmm <- function(data,
                                      response_var,
                                      random_effect,
                                      predictors,
                                      do_scale = TRUE,
                                      chains = 4,
                                      iter = 4000,
                                      warmup = 2000,
                                      cores = 4,
                                      seed = 1,
                                      adapt_delta = 0.95,
                                      save_models_rds = TRUE,
                                      model_dir = NULL) {
  
  results <- vector("list", length(predictors))
  names(results) <- predictors
  
  for (var in predictors) {
    
    # Subset only needed columns
    df <- data[, c(response_var, random_effect, var), drop = FALSE]
    df <- df[complete.cases(df), , drop = FALSE]
    
    # Skip if no variation
    if (nrow(df) < 50) {
      message("Skipping ", var, " (too few complete cases: n=", nrow(df), ")")
      results[[var]] <- NULL
      next
    }
    if (length(unique(df[[var]])) < 2) {
      message("Skipping ", var, " (no variation in predictor after filtering)")
      results[[var]] <- NULL
      next
    }
    
    # Standardize numeric predictor
    if (do_scale && is.numeric(df[[var]])) {
      df[[var]] <- as.numeric(scale(df[[var]]))
    }
    
    # Build formula
    fml <- as.formula(paste0(response_var, " ~ ", var, " + (1|", random_effect, ")"))
    
    # Fit model
    fit <- brm(
      formula = fml,
      data    = df,
      family  = bernoulli(),
      chains  = chains,
      iter    = iter,
      warmup  = warmup,
      cores   = cores,
      seed    = seed,
      prior   = c(
        set_prior("normal(0, 2.5)", class = "b"),
        set_prior("student_t(3, 0, 5)", class = "Intercept"),
        set_prior("exponential(1)", class = "sd")
      ),
      control   = list(adapt_delta = adapt_delta),
      save_pars = save_pars(all = TRUE),
      refresh   = 0
    )
    
    results[[var]] <- fit
    message("✓ Finished model: ", var)
    
    # Optionally save model object
    if (save_models_rds && !is.null(model_dir)) {
      saveRDS(fit, file.path(model_dir, paste0("univ_", var, ".rds")))
    }
  }
  
  results
}

univ_models <- fit_univariate_bayes_glmm(
  data            = FD,
  response_var    = RESPONSE_VAR,
  random_effect   = RANDOM_EFFECT,
  predictors      = predictors,
  do_scale        = TRUE,
  chains          = CHAINS,
  iter            = ITER,
  warmup          = WARMUP,
  cores           = CORES,
  seed            = SEED,
  adapt_delta     = ADAPT_DELTA,
  save_models_rds = SAVE_MODELS_RDS,
  model_dir       = MODEL_DIR
)

# -----------------------------
# 5) EXTRACT RESULTS (OR + 95% CrI + posterior probability)
# -----------------------------
extract_or_summary <- function(model_list) {
  
  result_list <- list()
  failed_vars <- character(0)
  
  for (var in names(model_list)) {
    
    fit <- model_list[[var]]
    if (is.null(fit)) {
      failed_vars <- c(failed_vars, var)
      next
    }
    
    tryCatch({
      
      draws_df <- posterior::as_draws_df(fit)
      b_terms  <- grep("^b_", colnames(draws_df), value = TRUE)
      
      if (length(b_terms) == 0) {
        failed_vars <- c(failed_vars, var)
        next
      }
      
      for (term in b_terms) {
        beta <- draws_df[[term]]
        OR   <- exp(beta)
        
        prob_pos <- mean(beta > 0)
        prob_neg <- mean(beta < 0)
        
        direction <- ifelse(prob_pos >= 0.5, "Positive", "Negative")
        post_prob <- ifelse(direction == "Positive", prob_pos, prob_neg)
        pseudo_p  <- 2 * min(prob_pos, prob_neg)
        
        result_list[[paste0(var, "::", term)]] <- data.frame(
          Predictor        = var,
          Term             = term,
          OR_median        = round(median(OR), 3),
          OR_CrI_lower_95  = round(quantile(OR, 0.025), 3),
          OR_CrI_upper_95  = round(quantile(OR, 0.975), 3),
          PosteriorProbDir = round(post_prob, 4),
          Direction        = direction,
          PseudoP          = round(pseudo_p, 4),
          stringsAsFactors = FALSE
        )
      }
      
      message("✓ Extracted summary: ", var)
      
    }, error = function(e) {
      message("✗ Extraction failed for ", var, ": ", e$message)
      failed_vars <<- c(failed_vars, var)
    })
  }
  
  summary_table <- if (length(result_list) > 0) do.call(rbind, result_list) else data.frame()
  list(summary = summary_table, failed = failed_vars)
}

univ_results <- extract_or_summary(univ_models)

# Save table
OUT_CSV <- file.path(OUT_DIR, "univariate_bayes_or_summary.csv")
write.csv(univ_results$summary, OUT_CSV, row.names = FALSE)

# Save a simple log of failed predictors
if (length(univ_results$failed) > 0) {
  writeLines(univ_results$failed, con = file.path(OUT_DIR, "failed_predictors.txt"))
}

message("Done. Results saved to: ", OUT_CSV)
