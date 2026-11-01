# ==============================================================================
# SCIENTIFIC CODE FOR BEE MORTALITY ANALYSIS: PART 1
# Data Preprocessing, Multivariate DTW Clustering, and Cluster Statistics
# ==============================================================================

# ==============================================================================
# 1. LIBRARY LOADING & ENVIRONMENT SETUP
# ==============================================================================
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(zoo)
library(imputeTS)
library(dtwclust)
library(knitr)
library(kableExtra)

# Install and load ComplexHeatmap if needed via Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)


# Graphics and scientific output settings
theme_set(theme_minimal(base_size = 12))
options(scipen = 999)   # Prevent scientific notation

# ==============================================================================
# 2. DATA LOADING & INITIAL TYPE CONVERSION
# ==============================================================================
# Read the primary dataset
df <- read.csv("dataf.csv", header = TRUE, stringsAsFactors = FALSE, 
               na.strings = c("NA", "#DIV/0!", ""))

# Convert time column to Date format and sort data chronologically
df$time <- dmy(df$time)
df <- df %>% arrange(time)

# Create a working copy for cleaning
df_clean <- df

# Define continuous variables for imputation and outlier treatment
continuous_vars <- c("temp", "t_o", "weight_kg", "weight_d", "h", "p")
numeric_cols <- continuous_vars

# ==============================================================================
# 3. MISSING DATA IMPUTATION & OUTLIER TREATMENT
# ==============================================================================
for (col in continuous_vars) {
  if (col %in% names(df_clean)) {
    # Calculate 7-day moving average for slow-changing or baseline metrics
    ma7 <- rollapply(df_clean[[col]], width = 7, 
                     FUN = function(x) mean(x, na.rm = TRUE), 
                     fill = NA, align = "center", partial = TRUE)
    
    # Linear interpolation for general gaps
    interpolated <- na.approx(df_clean[[col]], na.rm = FALSE, rule = 2)
    
    # Apply specific imputation strategy based on variable type
    if (col == "weight_kg") {
      df_clean[[col]] <- ifelse(is.na(df_clean[[col]]), ma7, df_clean[[col]])
      df_clean[[col]] <- na.approx(df_clean[[col]], na.rm = FALSE, rule = 2)
    } else {
      df_clean[[col]] <- interpolated
    }
  }
}

# Outlier removal (beyond 3 standard deviations) and neighborhood smoothing
for (col in continuous_vars) {
  if (col %in% names(df_clean)) {
    mean_val <- mean(df_clean[[col]], na.rm = TRUE)
    sd_val <- sd(df_clean[[col]], na.rm = TRUE)
    outliers <- abs(df_clean[[col]] - mean_val) > 3 * sd_val
    
    if (sum(outliers, na.rm = TRUE) > 0) {
      df_clean[[col]][outliers] <- NA
      df_clean[[col]] <- na.approx(df_clean[[col]], na.rm = FALSE, rule = 2)
    }
  }
}

# Quality control summary report
cat("\n===== DATA QUALITY CONTROL & IMPUTATION REPORT =====\n")
for (col in numeric_cols) {
  if (col %in% names(df_clean)) {
    na_before <- sum(is.na(df[[col]]))
    na_after <- sum(is.na(df_clean[[col]]))
    cat(sprintf("Variable %s: NaNs before=%d, after=%d, imputed=%d\n", 
                col, na_before, na_after, na_before - na_after))
  }
}

# Save cleaned dataset
write.csv(df_clean, "dataf_processed_cleaned.csv", row.names = FALSE)

# ==============================================================================
# 4. FEATURE ENGINEERING & ANALYTICAL DATASET PREPARATION
# ==============================================================================
df_r <- df_clean
df_r$time <- as.Date(df_r$time)
df_r <- df_r %>% arrange(time)

# Convert binary variables to factors
binary_cols <- c("died", "treatment", "feeding", "swarming", "gueencell", "key")
for (bcol in binary_cols) {
  if (bcol %in% names(df_r)) df_r[[bcol]] <- as.factor(df_r[[bcol]])
}

# Extract temporal indicators
df_r <- df_r %>%
  mutate(
    year = year(time),
    month = month(time),
    quarter = quarter(time),
    week = week(time)
  )

# Build main analytical dataset
df_analysis <- df_r %>%
  mutate(
    died = as.numeric(died),
    temp.in = temp,       # Internal temperature
    temp.out = t_o,       # External temperature
    weight = weight_kg,   # Hive weight
    weight_d = weight_d,  # Hive weight delta
    humidity = h,         # Relative humidity
    pressure = p,         # Atmospheric pressure
    treatment = as.numeric(treatment),
    feeding = as.numeric(feeding),
    swarming = as.numeric(swarming),
    gueencell = as.numeric(gueencell),
    key = as.numeric(key),
    year_month = floor_date(time, "month"),
    year_quarter = floor_date(time, "quarter")
  )

df_analysis$date <- as.Date(df_analysis$time)
df_filtered <- df_analysis

# ==============================================================================
# 5. TIME-SERIES MATRIX CONSTRUCTION & STANDARDIZATION
# ==============================================================================
date_range <- seq(min(df_filtered$date), max(df_filtered$date), by = "day")

# Helper function to create regular time-series matrices per hive
create_time_matrix <- function(df_filtered, date_range, variable_name) {
  unique_hives <- unique(df_filtered$key)
  n_hives <- length(unique_hives)
  n_dates <- length(date_range)
  
  cat("Creating telemetry matrix for variable:", variable_name, "\n")
  
  matrix_result <- matrix(NA, nrow = n_hives, ncol = n_dates)
  rownames(matrix_result) <- as.character(unique_hives)
  colnames(matrix_result) <- as.character(date_range)
  
  date_range_numeric <- as.numeric(date_range)
  
  for (hive in unique_hives) {
    hive_data <- df_filtered %>% filter(key == hive) %>% arrange(date)
    hive_dates_numeric <- as.numeric(hive_data$date)
    
    for (i in 1:nrow(hive_data)) {
      date_diff <- abs(date_range_numeric - hive_dates_numeric[i])
      closest_idx <- which.min(date_diff)
      
      if (date_diff[closest_idx] <= 86400) { # 1-day tolerance threshold
        matrix_result[as.character(hive), closest_idx] <- hive_data[[variable_name]][i]
      }
    }
  }
  
  # Impute internal matrix missing values linearly
  matrix_filled <- t(apply(matrix_result, 1, function(row) {
    if (all(is.na(row))) return(rep(0, length(row)))
    imputeTS::na_interpolation(row, option = "linear")
  }))
  
  matrix_filled[is.na(matrix_filled)] <- 0
  return(matrix_filled)
}

# Generate matrices for all telemetry variables
weight_matrix   <- create_time_matrix(df_filtered, date_range, "weight")
weight_d_matrix <- create_time_matrix(df_filtered, date_range, "weight_d")
temp_in_matrix  <- create_time_matrix(df_filtered, date_range, "temp")
temp_out_matrix <- create_time_matrix(df_filtered, date_range, "t_o")
h_matrix        <- create_time_matrix(df_filtered, date_range, "h")
p_matrix        <- create_time_matrix(df_filtered, date_range, "p")

# Standardization function (Z-score scaling per row)
standardize_matrix <- function(mat) {
  t(scale(t(mat)))
}

weight_scaled   <- standardize_matrix(weight_matrix)
weight_d_scaled <- standardize_matrix(weight_d_matrix)
temp_in_scaled  <- standardize_matrix(temp_in_matrix)
temp_out_scaled <- standardize_matrix(temp_out_matrix)
h_scaled        <- standardize_matrix(h_matrix)
p_scaled        <- standardize_matrix(p_matrix)

# Assemble multivariate list for clustering
series_list_multivariate <- lapply(1:nrow(weight_scaled), function(i) {
  mat_i <- cbind(
    weight_scaled[i, ],
    weight_d_scaled[i, ],
    temp_in_scaled[i, ], 
    temp_out_scaled[i, ], 
    h_scaled[i, ], 
    p_scaled[i, ]
  )
  mat_i[is.na(mat_i)] <- 0
  return(mat_i)
})

# ==============================================================================
# 6. MULTIVARIATE DTW CLUSTERING & OPTIMAL K EVALUATION (ELBOW METHOD)
# ==============================================================================
set.seed(123)
n_clusters_test <- 2:8
withinss_multivariate <- numeric(length(n_clusters_test))

time_series_length <- ncol(weight_scaled)
window_size_10perc <- ceiling(0.10 * time_series_length)

cat("\n--- Running Multivariate DTW Clustering Evaluation (Elbow Method) ---\n")

for (k in seq_along(n_clusters_test)) {
  current_k <- n_clusters_test[k]
  cat("Testing cluster size k =", current_k, "...\n")
  
  cluster_result <- tsclust(
    series_list_multivariate, 
    type = "partitional",
    k = current_k, 
    distance = "dtw_basic", 
    centroid = "dba",
    seed = 123,
    trace = FALSE,
    args = tsclust_args(dist = list(window.size = window_size_10perc))
  )
  
  withinss_multivariate[k] <- sum(cluster_result@cldist)
}

# Plot and save Elbow Curve
elbow_data <- data.frame(k = n_clusters_test, withinss = withinss_multivariate)
elbow_plot_multivariate <- ggplot(elbow_data, aes(x = k, y = withinss)) +
  geom_line(color = "#1f77b4", size = 1.2) +
  geom_point(color = "#d62728", size = 3.5) +
  scale_x_continuous(breaks = n_clusters_test) +
  labs(
    title = "Optimal Cluster Evaluation via Elbow Method",
    subtitle = "Multivariate DTW Distance for Hive Telemetry Metrics",
    x = "Number of Clusters (k)",
    y = "Total Within-Cluster Sum of Distances"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

print(elbow_plot_multivariate)

# ==============================================================================
# 7. FINAL CLUSTERING EXECUTION (k = 4) & TABLE 1 GENERATION
# ==============================================================================
optimal_k_multivariate <- 4  # Set optimal cluster representation
cat("Selected optimal number of clusters:", optimal_k_multivariate, "\n")

set.seed(123)
clustering_result_multivariate <- tsclust(
  series_list_multivariate, 
  type = "partitional",
  k = optimal_k_multivariate, 
  distance = "dtw_basic", 
  centroid = "dba",
  seed = 123,
  trace = TRUE,
  control = partitional_control(iter.max = 30, nrep = 5)
)

# Extract cluster labels and merge with dataset
cluster_labels_multivariate <- clustering_result_multivariate[[1]]@cluster
names(cluster_labels_multivariate) <- rownames(weight_matrix)

df_filtered <- df_filtered %>%
  mutate(cluster_multivariate = cluster_labels_multivariate[as.character(key)])

# Compute cluster statistics and mortality distribution (Table 1 Metrics)
df_filtered$died_num <- as.numeric(as.factor(df_filtered$died))

cluster_stats <- df_filtered %>%
  group_by(cluster_multivariate) %>%
  summarise(
    keys = list(unique(key)),
    dead_keys = list(unique(key[died_num == 2])),
    weight_mean = mean(weight, na.rm = TRUE),
    weight_sd = sd(weight, na.rm = TRUE),
    weight_d_mean = mean(weight_d, na.rm = TRUE) * 1000, # convert kg/day to g/day
    weight_d_sd = sd(weight_d, na.rm = TRUE) * 1000, # convert kg/day to g/day
    temp_in_mean = mean(temp.in, na.rm = TRUE),
    temp_in_sd = sd(temp.in, na.rm = TRUE),
    temp_out_mean = mean(temp.out, na.rm = TRUE),
    temp_out_sd = sd(temp.out, na.rm = TRUE),
    h_mean = mean(h, na.rm = TRUE),
    h_sd = sd(h, na.rm = TRUE),
    p_mean = mean(p, na.rm = TRUE),
    p_sd = sd(p, na.rm = TRUE),
    total_hives = n_distinct(key),
    dead_hives = n_distinct(key[died_num == 2]),
    mortality_rate = (n_distinct(key[died_num == 2]) / n_distinct(key)) * 100,
    .groups = "drop"
  )

# Format for clean display and table rendering
cluster_stats_display <- cluster_stats %>%
  mutate(
    keys = sapply(keys, paste, collapse = ", "),
    dead_keys = sapply(dead_keys, function(x) if(length(x) == 0) "None" else paste(x, collapse = ", "))
  )

print(cluster_stats_display)

# Generate Markdown Table matching Table 1 in manuscript
kable(cluster_stats_display, 
      caption = "Descriptive Telemetry Metrics and Mortality Distribution Across DTW Clusters",
      format = "markdown") %>%
  kable_styling(full_width = TRUE)

# Export processed cluster statistics
cluster_stats_clean <- cluster_stats %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = ", "))))

write.csv(cluster_stats_clean, "cluster_stats.csv", row.names = FALSE)






# ==============================================================================
# SCIENTIFIC CODE FOR BEE MORTALITY ANALYSIS: PART 2
# Machine Learning Classification (Random Forest, XGBoost, SVM) 
# Evaluated With and Without Multivariate DTW Clusters
# ==============================================================================

# ==============================================================================
# 1. ESSENTIAL LIBRARY LOADING FOR MACHINE LEARNING
# ==============================================================================
library(randomForest)
library(xgboost)
library(e1071)
library(caret)
library(dplyr)

# Setup 5-fold cross-validation control with ROC and probabilities
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary, # Computes ROC, Sensitivity, and Specificity
  sampling = "down"                  # Handles class imbalance via down-sampling
)

# ==============================================================================
# SCENARIO A: WITH MULTIVARIATE DTW CLUSTER FEATURE
# ==============================================================================

# Prepare dataset for Scenario A
set.seed(123)
rf_data <- df_filtered %>%
  select(died, cluster_multivariate, weight, weight_d, temp.in, temp.out, h, p) %>%
  na.omit() %>%
  mutate(
    died = factor(ifelse(died == 2, "Died", "Alive"), levels = c("Died", "Alive")),
    cluster_multivariate = as.factor(cluster_multivariate)
  )

# ------------------------------------------------------------------------------
# A.1. Random Forest Classifier (with CV)
# ------------------------------------------------------------------------------
set.seed(123)
rf_cv <- train(
  died ~ ., 
  data = rf_data, 
  method = "rf", 
  metric = "ROC",
  trControl = train_control,
  ntree = 500,
  importance = TRUE
)

# Extract best ROC and SD for Random Forest
best_idx <- which.max(rf_cv$results$ROC)
rf_roc_cv <- rf_cv$results$ROC[best_idx]
rf_roc_sd <- rf_cv$results$ROCSD[best_idx]

# Calculate normalized percentage importance
rf_imp <- importance(rf_cv$finalModel)[, 1]
rf_imp_pct <- rf_imp / sum(rf_imp) * 100

# Consolidate cluster dummy variables into a single unified 'cluster_multivariate' importance metric
cluster_cols <- grep("^cluster_multivariate", names(rf_imp_pct), value = TRUE)
cluster_total <- sum(rf_imp_pct[cluster_cols])

rf_imp_pct_unified <- c(
  cluster_multivariate = cluster_total,
  rf_imp_pct[!names(rf_imp_pct) %in% cluster_cols]
)

# ------------------------------------------------------------------------------
# A.2. XGBoost Classifier (with CV & scale_pos_weight)
# ------------------------------------------------------------------------------
set.seed(123)
y_num <- ifelse(rf_data$died == "Died", 1, 0)
X_matrix <- data.matrix(rf_data %>% select(-died))
weight_ratio <- sum(y_num == 0) / sum(y_num == 1)

folds <- createFolds(y_num, k = 5, list = TRUE)
xgb_acc_cv <- numeric(5)

for(i in 1:5) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(1:nrow(X_matrix), test_idx)
  
  dtrain <- xgb.DMatrix(X_matrix[train_idx, ], label = y_num[train_idx])
  dtest <- xgb.DMatrix(X_matrix[test_idx, ], label = y_num[test_idx])
  
  xgb_model <- xgb.train(
    params = list(
      objective = "binary:logistic", 
      eval_metric = "logloss",
      scale_pos_weight = weight_ratio,
      eta = 0.1, 
      max_depth = 5
    ), 
    data = dtrain, 
    nrounds = 50, 
    verbose = 0
  )
  
  pred <- predict(xgb_model, dtest)
  pred_class <- ifelse(pred > 0.5, 1, 0)
  xgb_acc_cv[i] <- mean(pred_class == y_num[test_idx])
}

xgb_acc_mean <- mean(xgb_acc_cv)
xgb_acc_sd <- sd(xgb_acc_cv)

# Train full XGBoost model for feature importance
dtrain_full <- xgb.DMatrix(X_matrix, label = y_num)
xgb_full <- xgb.train(
  params = list(
    objective = "binary:logistic", 
    eval_metric = "logloss",
    scale_pos_weight = weight_ratio,
    eta = 0.1, 
    max_depth = 5
  ), 
  data = dtrain_full, 
  nrounds = 50, 
  verbose = 0
)

xgb_importance <- xgb.importance(model = xgb_full, feature_names = colnames(X_matrix))
xgb_imp <- data.frame(Overall = xgb_importance$Gain, row.names = xgb_importance$Feature)
xgb_imp_pct <- (xgb_imp$Overall / sum(xgb_imp$Overall)) * 100
names(xgb_imp_pct) <- rownames(xgb_imp)

# ------------------------------------------------------------------------------
# A.3. Support Vector Machines - SVM (with CV)
# ------------------------------------------------------------------------------
set.seed(123)
rf_data_features <- subset(rf_data, select = -c(died))

# 2. أنشئ الـ dummies باستخدام علامة النقطة (.) فقط لتعني "كل الأعمدة المتبقية"
dummies <- dummyVars(~ ., data = rf_data_features, fullRank = FALSE)

# 3. طبق التنبؤ الآن ولن يظهر أي خطأ
X_svm <- as.data.frame(predict(dummies, newdata = rf_data_features))
                       
y_svm <- rf_data$died

svm_cv <- train(
  x = X_svm,
  y = y_svm,
  method = "svmRadial",
  metric = "ROC",
  trControl = train_control,
  preProcess = c("center", "scale")
)

best_idx_svm <- which.max(svm_cv$results$ROC)
svm_roc_cv <- svm_cv$results$ROC[best_idx_svm]
svm_roc_sd <- svm_cv$results$ROCSD[best_idx_svm]

svm_imp <- varImp(svm_cv, useModel = FALSE)$importance
svm_imp_mean <- rowMeans(svm_imp)
svm_imp_pct <- (svm_imp_mean / sum(svm_imp_mean)) * 100
names(svm_imp_pct) <- rownames(svm_imp)

cluster_cols_svm <- grep("^cluster_multivariate", names(svm_imp_pct), value = TRUE)
if(length(cluster_cols_svm) > 0) {
  cluster_total_svm <- sum(svm_imp_pct[cluster_cols_svm])
  svm_imp_pct <- c(
    cluster_multivariate = cluster_total_svm,
    svm_imp_pct[!names(svm_imp_pct) %in% cluster_cols_svm]
  )
}
svm_imp_pct <- sort(svm_imp_pct, decreasing = TRUE)


# ==============================================================================
# SCENARIO B: WITHOUT MULTIVARIATE DTW CLUSTER FEATURE
# ==============================================================================

set.seed(123)
rf_data_o <- df_filtered %>%
  select(died, weight, weight_d, temp.in, temp.out, h, p) %>%
  na.omit() %>%
  mutate(
    died = factor(ifelse(died == 2, "Died", "Alive"), levels = c("Died", "Alive"))
  )

# ------------------------------------------------------------------------------
# B.1. Random Forest Classifier (Without Cluster)
# ------------------------------------------------------------------------------
set.seed(123)
rf_cv_o <- train(
  died ~ ., 
  data = rf_data_o, 
  method = "rf", 
  metric = "ROC",
  trControl = train_control,
  ntree = 500,
  importance = TRUE
)

best_idx_o <- which.max(rf_cv_o$results$ROC)
rf_roc_cv_o <- rf_cv_o$results$ROC[best_idx_o]
rf_roc_sd_o <- rf_cv_o$results$ROCSD[best_idx_o]

rf_imp_o <- importance(rf_cv_o$finalModel)[, 1]
rf_imp_pct_o <- rf_imp_o / sum(rf_imp_o) * 100

# ------------------------------------------------------------------------------
# B.2. XGBoost Classifier (Without Cluster)
# ------------------------------------------------------------------------------
set.seed(123)
y_num_o <- ifelse(rf_data_o$died == "Died", 1, 0)
X_matrix_o <- data.matrix(rf_data_o %>% select(-died))
weight_ratio_o <- sum(y_num_o == 0) / sum(y_num_o == 1)

folds_o <- createFolds(y_num_o, k = 5, list = TRUE)
xgb_acc_cv_o <- numeric(5)

for(i in 1:5) {
  test_idx_o <- folds_o[[i]]
  train_idx_o <- setdiff(1:nrow(X_matrix_o), test_idx_o)
  
  dtrain_o <- xgb.DMatrix(X_matrix_o[train_idx_o, ], label = y_num_o[train_idx_o])
  dtest_o <- xgb.DMatrix(X_matrix_o[test_idx_o, ], label = y_num_o[test_idx_o])
  
  xgb_model_o <- xgb.train(
    params = list(
      objective = "binary:logistic", 
      eval_metric = "logloss",
      scale_pos_weight = weight_ratio_o,
      eta = 0.1, 
      max_depth = 5
    ), 
    data = dtrain_o, 
    nrounds = 50, 
    verbose = 0
  )
  
  pred_o <- predict(xgb_model_o, dtest_o)
  pred_class_o <- ifelse(pred_o > 0.5, 1, 0)
  xgb_acc_cv_o[i] <- mean(pred_class_o == y_num_o[test_idx_o])
}

xgb_acc_mean_o <- mean(xgb_acc_cv_o)
xgb_acc_sd_o <- sd(xgb_acc_cv_o)

dtrain_full_o <- xgb.DMatrix(X_matrix_o, label = y_num_o)
xgb_full_o <- xgb.train(
  params = list(
    objective = "binary:logistic", 
    eval_metric = "logloss",
    scale_pos_weight = weight_ratio_o,
    eta = 0.1, 
    max_depth = 5
  ), 
  data = dtrain_full_o, 
  nrounds = 50, 
  verbose = 0
)

xgb_importance_o <- xgb.importance(model = xgb_full_o, feature_names = colnames(X_matrix_o))
xgb_imp_o <- data.frame(Overall = xgb_importance_o$Gain, row.names = xgb_importance_o$Feature)
xgb_imp_pct_o <- (xgb_imp_o$Overall / sum(xgb_imp_o$Overall)) * 100
names(xgb_imp_pct_o) <- rownames(xgb_imp_o)

# ------------------------------------------------------------------------------
# B.3. Support Vector Machines - SVM (Without Cluster)
# ------------------------------------------------------------------------------
set.seed(123)
rf_data_o$died <- as.factor(rf_data_o$died)
dummies_o <- dummyVars(~ . - died, data = rf_data_o)
X_svm_o <- as.data.frame(predict(dummies_o, newdata = rf_data_o))
y_svm_o <- rf_data_o$died

svm_cv_o <- train(
  x = X_svm_o,
  y = y_svm_o,
  method = "svmRadial",
  metric = "ROC",
  trControl = train_control,
  preProcess = c("center", "scale")
)

best_idx_svm_o <- which.max(svm_cv_o$results$ROC)
svm_roc_cv_o <- svm_cv_o$results$ROC[best_idx_svm_o]
svm_roc_sd_o <- svm_cv_o$results$ROCSD[best_idx_svm_o]

svm_imp_o <- varImp(svm_cv_o, useModel = FALSE)$importance
svm_imp_mean_o <- rowMeans(svm_imp_o)
svm_imp_pct_o <- (svm_imp_mean_o / sum(svm_imp_mean_o)) * 100
names(svm_imp_pct_o) <- rownames(svm_imp_o)
svm_imp_pct_o <- sort(svm_imp_pct_o, decreasing = TRUE)


# ==============================================================================
# 3. ALIGNMENT AND TABLE GENERATION (MATCHING TABLE 2)
# ==============================================================================

all_vars <- c("cluster_multivariate", "weight", "weight_d", "temp.in", "temp.out", "h", "p")
all_vars_o <- c("weight", "weight_d", "temp.in", "temp.out", "h", "p")

get_imp_aligned <- function(imp_vec, vars) {
  sapply(vars, function(v) {
    if (v %in% names(imp_vec)) unname(imp_vec[v]) else NA
  })
}

# Align importance vectors
rf_imp_aligned   <- get_imp_aligned(rf_imp_pct_unified, all_vars)
xgb_imp_aligned  <- get_imp_aligned(xgb_imp_pct, all_vars)
svm_imp_aligned  <- get_imp_aligned(svm_imp_pct, all_vars)

rf_imp_aligned_o  <- get_imp_aligned(rf_imp_pct_o, all_vars_o)
xgb_imp_aligned_o <- get_imp_aligned(xgb_imp_pct_o, all_vars_o)
svm_imp_aligned_o <- get_imp_aligned(svm_imp_pct_o, all_vars_o)

# ------------------------------------------------------------------------------
# COMPILING PERFORMANCE & ACCURACY TABLES
# ------------------------------------------------------------------------------

# Tables for Scenario A (With Cluster)
accuracy_table <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM"),
  ROC = round(c(rf_roc_cv, xgb_acc_mean, svm_roc_cv), 4),
  SD = round(c(rf_roc_sd, xgb_acc_sd, svm_roc_sd), 4)
)

importance_table <- data.frame(
  Variable = all_vars,
  RandomForest = round(rf_imp_aligned, 2),
  XGBoost = round(xgb_imp_aligned, 2),
  SVM = round(svm_imp_aligned, 2)
)

# Tables for Scenario B (Without Cluster)
accuracy_table_o <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM"),
  ROC = round(c(rf_roc_cv_o, xgb_acc_mean_o, svm_roc_cv_o), 4),
  SD = round(c(rf_roc_sd_o, xgb_acc_sd_o, svm_roc_sd_o), 4)
)

importance_table_o <- data.frame(
  Variable = all_vars_o,
  RandomForest = round(rf_imp_aligned_o, 2),
  XGBoost = round(xgb_imp_aligned_o, 2),
  SVM = round(svm_imp_aligned_o, 2)
)

# ------------------------------------------------------------------------------
# PRINTING RESULTS TO CONSOLE
# ------------------------------------------------------------------------------
cat("--- Performance Metrics WITH MULTIVARIATE (ROC ± SD) ---\n")
print(accuracy_table)
cat("--- Performance Metrics WITHOUT MULTIVARIATE (ROC ± SD) ---\n")
print(accuracy_table_o)
cat("\n--- Variable Importance WITH MULTIVARIATE (%) ---\n")
print(importance_table)
cat("\n--- Variable Importance WITHOUT MULTIVARIATE (%) ---\n")
print(importance_table_o)

# ------------------------------------------------------------------------------
# EXPORTING TABLES TO CSV
# ------------------------------------------------------------------------------
write.csv(accuracy_table, "table_performance_with_cluster.csv", row.names = FALSE)
write.csv(importance_table, "table_importance_with_cluster.csv", row.names = FALSE)
write.csv(accuracy_table_o, "table_performance_without_cluster.csv", row.names = FALSE)
write.csv(importance_table_o, "table_importance_without_cluster.csv", row.names = FALSE)


# ==============================================================================
# SCIENTIFIC CODE FOR BEE MORTALITY ANALYSIS: PART 3
# Model Interpretability (SHAP Analysis) & Statistical Risk Modeling (Firth's Regression)
# ==============================================================================

# ==============================================================================
# 1. SHAP ANALYSIS & COMBINED PUBLICATION PLOT (A & B)
# ==============================================================================
# ==========================================
# 0. LOAD REQUIRED LIBRARIES
# ==========================================
library(fastshap)
library(shapviz)
library(ggplot2)
library(patchwork)
library(dplyr)

# ==========================================
# 1. PREPARE PREDICTOR DATA
# ==========================================
# Exclude the target variable 'died' from the dataset to create X_rf
X_rf <- rf_data %>% select(-died)

# This protects against "subscript out of bounds" errors during SHAP calculation
p_fun_rf_safe <- function(object, newdata) {
  preds <- predict(object, newdata = newdata, type = "prob")
  
  # Ensure the probability of the "Died" class is extracted accurately
  if ("Died" %in% colnames(preds)) {
    return(preds[, "Died"])
  } else {
    return(preds[, 1])
  }
}

# ==========================================
# 3. CALCULATE SHAP VALUES
# ==========================================
set.seed(123)
shap_rf <- fastshap::explain(
  object = rf_cv,             # Direct caret model object
  X = X_rf,
  pred_wrapper = p_fun_rf_safe,
  nsim = 50                   # Number of Monte Carlo simulations
)

# ==========================================
# 4. CONVERT TO SHAPVIZ OBJECT
# ==========================================
shp_rf_obj <- shapviz(shap_rf, X = X_rf)

# ==========================================
# 5. DEFINE UNIFIED FEATURE LABELS DICTIONARY
# ==========================================
new_labels <- c(
  "weight"               = "Weight (kg)",
  "temp.in"              = "Internal Temp.",
  "cluster_multivariate" = "Cluster",
  "p"                    = "Pressure",
  "h"                    = "Humidity",
  "temp.out"             = "External Temp.",
  "weight_d"             = "Delta weight (kg/day)"
)

# ==========================================
# 6. PREPARE AND CREATE PLOT (A): BAR PLOT
# ==========================================
# Extract base bar plot data from shapviz
p_shap_bar_base <- sv_importance(shp_rf_obj, kind = "bar")
plot_data <- p_shap_bar_base$data

# Apply factor levels and new labels in reverse order for proper vertical alignment
plot_data$feature <- factor(
  plot_data$feature, 
  levels = rev(names(new_labels)), 
  labels = rev(new_labels)
)

# Build Plot A
p_a <- ggplot(plot_data, aes(x = value, y = reorder(feature, value))) +
  geom_col(fill = "#2B5C8F", width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", value)), hjust = -0.15, size = 3.5, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Mean |SHAP value|", y = NULL, tag = "(A)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", margin = ggplot2::margin(t = 8)),
    plot.tag = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  )

# ==========================================
# 7. CREATE PLOT (B): BEESWARM PLOT
# ==========================================
p_b <- sv_importance(shp_rf_obj, kind = "beeswarm") +
  scale_y_discrete(
    breaks = names(new_labels),
    labels = new_labels
  ) +
  labs(x = "SHAP value (impact on model output)", y = NULL, tag = "(B)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(), # Remove Y-axis labels to avoid duplication
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", margin = ggplot2::margin(t = 8)),
    plot.tag = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  )

# ==========================================
# 8. COMBINE AND SAVE FINAL PLOT
# ==========================================
# Combine both plots side by side using patchwork
p_combined <- p_a + p_b + plot_layout(widths = c(1, 1.4))

# Display the final combined plot
print(p_combined)

# Save high-resolution image for scientific publication (300 DPI)
ggsave(
  filename = "SHAP_Importance_Combined.png", 
  plot = p_combined, 
  width = 10, 
  height = 4.5, 
  dpi = 300
)

cat("\n========== COMBINED SHAP PLOT CREATED AND SAVED SUCCESSFULLY ==========\n")


# ==============================================================================
# SECTION 3: STATISTICAL RISK MODELING (FIRTH'S PENALIZED LOGISTIC REGRESSION)
# MULTICOLLINEARITY ANALYSIS & MODEL COMPARISON PIPELINE
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD REQUIRED LIBRARIES
# ------------------------------------------------------------------------------
library(logistf)     # Firth's penalized likelihood logistic regression
library(car)         # Variance Inflation Factor (VIF) and diagnostics
library(entropy)     # Information theory metrics
library(sandwich)   # Cluster-robust standard errors
library(lmtest)      # Hypothesis testing
library(dplyr)       # Data manipulation
library(Metrics)     # Performance metrics (MAE, RMSE)
library(caret)       # Confusion matrix & cross-validation metrics
library(pROC)        # ROC curve calculation and AUC assessment
library(ggplot2)     # Publication-grade visualizations
library(ggrepel)     # Advanced text and label repelling for plots

# ------------------------------------------------------------------------------
# 2. DATA PREPROCESSING & FACTOR ENCODING
# ------------------------------------------------------------------------------
# Binary encoding and factor level conversions
df_filtered <- df_filtered %>%
  mutate(
    # Recode binary indicators: 2 -> 1 (Event), 1 -> 0 (Control)
    feeding   = factor(ifelse(feeding == 2, 1, 0)),
    treatment = factor(ifelse(treatment == 2, 1, 0)),
    swarming  = factor(ifelse(swarming == 2, 1, 0)),
    died_binary = ifelse(died == 2, 1, 0),
    
    # Factor variables setup
    year                 = as.factor(year),
    quarter              = as.factor(quarter),
    cluster_multivariate = as.factor(cluster_multivariate)
  )

# Inspect unique category counts per categorical/factor variable
sapply(df_filtered, function(x) if (is.factor(x) || is.character(x)) length(unique(na.omit(x))))

# Check outcome distribution across multivariate clusters
table(Cluster = df_filtered$cluster_multivariate, Outcome = df_filtered$died_binary)

# Set Cluster 2 as the baseline reference category
df_filtered$cluster_multivariate <- relevel(df_filtered$cluster_multivariate, ref = "2")

# Continuous feature transformation: Z-score standardization for weight delta
summary(df_filtered$weight_d)
df_filtered$weight_d_scaled <- as.numeric(scale(df_filtered$weight_d))

# ------------------------------------------------------------------------------
# 3. HELPER FUNCTIONS FOR STATISTICAL REPORTING
# ------------------------------------------------------------------------------
#' Extract Summary Statistics from Firth's Logistic Model
#'
#' @param model A fitted 'logistf' model object.
#' @param label A character string defining the model specification name.
#' @return A formatted data frame containing coefficients, odds ratios, 95% CIs, and p-values.
extract_firth_summary <- function(model, label) {
  coefs <- coef(model)
  ci    <- confint(model)  # Profile penalization confidence intervals
  
  data.frame(
    Model      = label,
    Variable   = names(coefs),
    Estimate   = round(coefs, 4),
    Odds_Ratio = round(exp(coefs), 3),
    CI_Lower   = round(exp(ci[, 1]), 3),
    CI_Upper   = round(exp(ci[, 2]), 3),
    p_value    = round(model$prob, 4)
  )
}

# Define stable convergence control parameters for profile likelihood optimization
stable_control <- logistf.control(maxit = 1000, maxstep = 0.2, lconv = 0.0001)

# ------------------------------------------------------------------------------
# 4. MODEL SPECIFICATIONS & ESTIMATION
# ------------------------------------------------------------------------------

# --- Model A: Direct Physical Sensors ---
glm_model_A <- logistf(
  died_binary ~ weight + weight_d_scaled + temp.in + temp.out + h + p + 
    year + quarter + gueencell + feeding + treatment + swarming,
  data    = df_filtered,
  pl      = TRUE,
  control = stable_control
)

results_A <- extract_firth_summary(glm_model_A, "Model A (Direct Sensors)")
print(results_A)
write.csv(results_A, "Model_A_Direct_Sensors.csv", row.names = FALSE)


# --- Model B: Cluster Profiles (Multivariate Cluster Only) ---
glm_model_B <- logistf(
  died_binary ~ cluster_multivariate + 
    year + quarter + gueencell + feeding + treatment + swarming,
  data    = df_filtered,
  pl      = TRUE,
  control = stable_control
)

results_B <- extract_firth_summary(glm_model_B, "Model B (Cluster Profiles)")
print(results_B)
write.csv(results_B, "Model_B_Cluster_Profiles.csv", row.names = FALSE)


# --- Model C: Fully Integrated (Sensors + Cluster Profiles) ---
glm_model_C <- logistf(
  died_binary ~ weight + weight_d_scaled + temp.in + temp.out + h + p + 
    cluster_multivariate + 
    year + quarter + gueencell + feeding + treatment + swarming,
  data    = df_filtered,
  pl      = TRUE,
  control = stable_control
)

results_C <- extract_firth_summary(glm_model_C, "Model C (Fully Integrated)")
print(results_C)
write.csv(results_C, "Model_C_Fully_Integrated.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. PERFORMANCE METRICS EVALUATION PIPELINE
# ------------------------------------------------------------------------------

# Generate predicted probabilities for all candidate models
prob_A <- predict(glm_model_A, newdata = df_filtered, type = "response")
prob_B <- predict(glm_model_B, newdata = df_filtered, type = "response")
prob_C <- predict(glm_model_C, newdata = df_filtered, type = "response")

actual_class <- df_filtered$died_binary

#' Comprehensive Model Performance Evaluator
#'
#' @param prob Vector of predicted probabilities.
#' @param actual Vector of true class labels (0/1).
#' @param cutoff Classification decision threshold.
#' @param model_label Model descriptive tag.
#' @return Data frame of performance indicators.
evaluate_model_performance <- function(prob, actual, cutoff, model_label) {
  pred_class    <- factor(ifelse(prob >= cutoff, 1, 0), levels = c(0, 1))
  actual_factor <- factor(actual, levels = c(0, 1))
  
  cm <- caret::confusionMatrix(pred_class, actual_factor, positive = "1")
  
  roc_obj <- pROC::roc(actual, prob, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  data.frame(
    Model             = model_label,
    Optimal_Cutoff    = round(cutoff, 4),
    ROC_AUC           = round(auc_val, 4),
    Accuracy          = round(cm$overall["Accuracy"], 4),
    Balanced_Accuracy = round(cm$byClass["Balanced Accuracy"], 4),
    Kappa             = round(cm$overall["Kappa"], 4),
    Sensitivity       = round(cm$byClass["Sensitivity"], 4),
    Specificity       = round(cm$byClass["Specificity"], 4),
    Precision         = round(cm$byClass["Precision"], 4),
    F1_Score          = round(cm$byClass["F1"], 4)
  )
}

# Compute ROC objects
roc_A <- pROC::roc(actual_class, prob_A, quiet = TRUE)
roc_B <- pROC::roc(actual_class, prob_B, quiet = TRUE)
roc_C <- pROC::roc(actual_class, prob_C, quiet = TRUE)

# Calculate optimal decision threshold based on Youden's J statistic
best_cutoff_A <- coords(roc_A, "best", best.method = "youden")$threshold
best_cutoff_B <- coords(roc_B, "best", best.method = "youden")$threshold
best_cutoff_C <- coords(roc_C, "best", best.method = "youden")$threshold

cat("Optimal Decision Threshold (Model A):", best_cutoff_A, "\n")
cat("Optimal Decision Threshold (Model B):", best_cutoff_B, "\n")
cat("Optimal Decision Threshold (Model C):", best_cutoff_C, "\n")

# Compute comprehensive metrics table
metrics_A <- evaluate_model_performance(prob_A, actual_class, best_cutoff_A, "Model A (Direct Sensors)")
metrics_B <- evaluate_model_performance(prob_B, actual_class, best_cutoff_B, "Model B (Cluster Profiles)")
metrics_C <- evaluate_model_performance(prob_C, actual_class, best_cutoff_C, "Model C (Fully Integrated)")

# Consolidate and export summary metrics
performance_summary <- rbind(metrics_A, metrics_B, metrics_C)
print(performance_summary)
write.csv(performance_summary, "Performance_Summary_Optimized.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. PUBLICATION-GRADE ROC CURVE VISUALIZATION
# ------------------------------------------------------------------------------

# Consolidate ROC curves data
roc_all <- rbind(
  data.frame(Specificity = roc_A$specificities, Sensitivity = roc_A$sensitivities, Model = "Model A (Direct Sensors)"),
  data.frame(Specificity = roc_B$specificities, Sensitivity = roc_B$sensitivities, Model = "Model B (Cluster Profiles)"),
  data.frame(Specificity = roc_C$specificities, Sensitivity = roc_C$sensitivities, Model = "Model C (Fully Integrated)")
)

# Extract coordinates for optimal threshold points
best_A <- coords(roc_A, "best", best.method = "youden")
best_B <- coords(roc_B, "best", best.method = "youden")
best_C <- coords(roc_C, "best", best.method = "youden")

optimal_points <- data.frame(
  Specificity = c(best_A$specificity, best_B$specificity, best_C$specificity),
  Sensitivity = c(best_A$sensitivity, best_B$sensitivity, best_C$sensitivity),
  Cutoff      = c(best_cutoff_A, best_cutoff_B, best_cutoff_C),
  Model       = c("Model A (Direct Sensors)", 
                  "Model B (Cluster Profiles)", 
                  "Model C (Fully Integrated)")
)

# Plot ROC Curves via ggplot2
p_roc_optimized <- ggplot() +
  # 1. Base ROC Lines
  geom_line(data = roc_all, aes(x = Specificity, y = Sensitivity, color = Model), linewidth = 1.2) +
  
  # 2. Random Chance Line (Diagonal)
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  
  # 3. Optimal Cutoff Points
  geom_point(data = optimal_points, aes(x = Specificity, y = Sensitivity, color = Model), size = 3.5, shape = 19) +
  
  # 4. Repelled Threshold Labels
  geom_label_repel(
    data              = optimal_points, 
    aes(x = Specificity, y = Sensitivity, label = sprintf("Cutoff: %.3f", Cutoff), color = Model),
    box.padding       = 0.8,
    point.padding     = 0.6,
    min.segment.length = 0,
    force             = 5,
    max.overlaps      = Inf,
    show.legend       = FALSE,
    fontface          = "bold",
    size              = 3.8,
    fill              = "white",
    label.size        = 0.3,
    segment.color     = "grey50",
    segment.linewidth = 0.5
  ) +
  
  # 5. Axes Configuration
  scale_x_reverse(limits = c(1, 0), expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.08, 0.12))) +
  
  # 6. High-Contrast Publication Color Palette
  scale_color_manual(values = c(
    "Model A (Direct Sensors)"   = "#2563EB", # Royal Blue
    "Model B (Cluster Profiles)" = "#7C3AED", # Deep Purple
    "Model C (Fully Integrated)" = "#DB2777"  # Vivid Pink/Magenta
  )) +
  
  # 7. Labels & Clean Academic Theme
  labs(
    x     = "1 - Specificity (False Positive Rate)",
    y     = "Sensitivity (True Positive Rate)",
    color = "Model Specifications",
    title = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position      = c(0.75, 0.20),
    legend.justification = c(0.5, 0.5),
    legend.background    = element_rect(fill = "white", color = "black", linewidth = 0.4),
    legend.title         = element_text(face = "bold", size = 10),
    legend.text          = element_text(size = 9),
    axis.title           = element_text(face = "bold", size = 11),
    axis.text            = element_text(color = "black", size = 10),
    panel.grid.major     = element_line(color = "grey92", linewidth = 0.5),
    panel.grid.minor     = element_blank()
  )

# Display Plot
print(p_roc_optimized)

# Save high-resolution figures for submission (600 DPI)
ggsave("Figure_ROC_Comparison_Firth.png", plot = p_roc_optimized, width = 6.5, height = 5.5, dpi = 600)
ggsave("Figure_ROC_Comparison_Firth.tiff", plot = p_roc_optimized, width = 6.5, height = 5.5, dpi = 600, compression = "lzw")


# Clear environment and free memory
rm(list = ls())
gc()
