# install.packages(c("rms", "pROC", "rmda"))
library(rms)
library(pROC)
library(rmda)

output_dir <- "D:/Desktop/Rproject/shen protein/A2s_model_output/three_model_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
B_calibration <- 1000L
seed <- 123L

panel <- c("LEFTY2", 'LAMA2')#"THSD4", "COL2A1", "MFAP5"
clinical_vars <- c("current_smoking", "primary_site_subregion",
                   "extranodal_extension", "cN_01_23_group")
all_predictors <- c(clinical_vars, panel)

dat_work <- dat
recode_known <- function(x, dictionary, label) {
  z <- trimws(as.character(x))
  z[z == ""] <- NA_character_
  bad <- setdiff(unique(z[!is.na(z)]), names(dictionary))
  if (length(bad)) stop(label, " 存在未知编码：", paste(bad, collapse = ", "))
  unname(dictionary[z])
}

if ("drug_response_group" %in% names(dat_work)) {
  response <- tolower(trimws(as.character(dat_work$drug_response_group)))
  dat_work$ict_failure <- as.integer(recode_known(
    response, c(sensitive = "0", insensitive = "1"), "drug_response_group"))
} else if ("ict_failure" %in% names(dat_work)) {
  dat_work$ict_failure <- as.integer(recode_known(
    dat_work$ict_failure, c("0" = "0", "1" = "1"), "ict_failure"))
} else stop("dat 缺少 drug_response_group 或 ict_failure")

needed <- c("current_smoking", "primary_site_subregion", "extranodal_extension", panel)
missing_columns <- setdiff(needed, names(dat_work))
if (length(missing_columns)) stop("dat 缺少：", paste(missing_columns, collapse = ", "))

dat_work$current_smoking <- factor(recode_known(
  dat_work$current_smoking, c("0" = "No", "1" = "Current", No = "No", Current = "Current"),
  "current_smoking"), levels = c("No", "Current"))
dat_work$extranodal_extension <- factor(recode_known(
  dat_work$extranodal_extension, c("0" = "No", "1" = "Yes", No = "No", Yes = "Yes"),
  "extranodal_extension"), levels = c("No", "Yes"))
site_levels <- c("Pyriform sinus", "Posterior hypopharyngeal wall", "Postcricoid")
dat_work$primary_site_subregion <- factor(recode_known(
  dat_work$primary_site_subregion,
  c(setNames(site_levels, c("1", "2", "3")), setNames(site_levels, site_levels)),
  "primary_site_subregion"), levels = site_levels)

if ("cN_01_23_group" %in% names(dat_work)) {
  dat_work$cN_01_23_group <- factor(recode_known(
    dat_work$cN_01_23_group,
    c("0" = "N0-1", "1" = "N2-3", "N0-1" = "N0-1", "N2-3" = "N2-3"),
    "cN_01_23_group"), levels = c("N0-1", "N2-3"))
} else if ("cN" %in% names(dat_work)) {
  dat_work$cN_01_23_group <- factor(recode_known(
    dat_work$cN, c("0" = "N0-1", "1" = "N0-1", "2" = "N2-3", "3" = "N2-3"),
    "cN"), levels = c("N0-1", "N2-3"))
} else stop("dat 缺少 cN_01_23_group 或 cN")

if (!all(vapply(dat_work[panel], is.numeric, logical(1)))) {
  stop("蛋白表达必须为 numeric，请检查 str(dat[panel])；不要把连续表达转为 factor")
}
if (any(vapply(dat_work[panel], function(x) any(!is.finite(x) & !is.na(x)), logical(1)))) {
  stop("蛋白表达中存在 Inf/-Inf")
}

keep <- complete.cases(dat_work[, c("ict_failure", all_predictors)])
analysis_dat <- droplevels(dat_work[keep, c("ict_failure", all_predictors), drop = FALSE])
source_rows <- which(keep)
if (!setequal(unique(analysis_dat$ict_failure), c(0L, 1L))) stop("有效样本必须包含两类结局")
if (any(vapply(analysis_dat[all_predictors], function(x) length(unique(x)) < 2L, logical(1)))) {
  stop("某个预测变量在完整样本中没有变异")
}
refs <- c(current_smoking = "No", primary_site_subregion = "Pyriform sinus",
          extranodal_extension = "No", cN_01_23_group = "N0-1")
for (v in names(refs)) {
  if (!refs[v] %in% levels(analysis_dat[[v]])) stop(v, " 的参照组没有有效样本")
  analysis_dat[[v]] <- relevel(analysis_dat[[v]], ref = refs[v])
}
sample_summary <- data.frame(N_input = nrow(dat_work), N_complete = nrow(analysis_dat),
                             N_excluded = sum(!keep), Resistant = sum(analysis_dat$ict_failure == 1L),
                             Sensitive = sum(analysis_dat$ict_failure == 0L))
print(sample_summary)
write.csv(sample_summary, file.path(output_dir, "sample_summary.csv"), row.names = FALSE)

dd_comparison <- datadist(analysis_dat)
options(datadist = "dd_comparison")
fit_clin <- lrm(
  ict_failure ~ current_smoking + primary_site_subregion +
    extranodal_extension + cN_01_23_group,
  data = analysis_dat, x = TRUE, y = TRUE)
fit_protein <- lrm(
  ict_failure ~ LEFTY2 + LAMA2,
  data = analysis_dat, x = TRUE, y = TRUE)
fit <- lrm(
  ict_failure ~ current_smoking + primary_site_subregion +
    extranodal_extension + cN_01_23_group + LEFTY2 + LAMA2,
  data = analysis_dat, x = TRUE, y = TRUE)
models <- list(Clinical = fit_clin, Protein = fit_protein, Joint = fit)
if (any(vapply(models, function(m) isTRUE(m$fail), logical(1)))) stop("模型拟合失败，请检查输出")
if (any(vapply(models, function(m) any(!is.finite(coef(m))), logical(1)))) stop("存在不可估计的系数")
lapply(models, print)

fit_joint <- fit
model_dat_joint <- analysis_dat
fit_A <- fit
model_dat_A <- analysis_dat

or_table <- do.call(rbind, lapply(names(models), function(nm) {
  m <- models[[nm]]
  beta <- coef(m)[-1]
  se <- sqrt(diag(vcov(m)))[-1]
  data.frame(Model = nm, Coefficient = names(beta), Beta = unname(beta), SE = unname(se),
             OR = exp(beta), CI_low = exp(beta - qnorm(.975) * se),
             CI_high = exp(beta + qnorm(.975) * se),
             P_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE), row.names = NULL)
}))
write.csv(or_table, file.path(output_dir, "multivariable_OR.csv"), row.names = FALSE)
for (nm in names(models)) {
  write.csv(as.data.frame(unclass(anova(models[[nm]]))),
            file.path(output_dir, paste0(nm, "_overall_Wald.csv")))
}

y <- analysis_dat$ict_failure
lp_list <- lapply(models, function(m) as.numeric(predict(m, newdata = analysis_dat, type = "lp")))
prob_list <- lapply(lp_list, plogis)
roc_list <- lapply(prob_list, function(p) pROC::roc(y, p, levels = c(0, 1),
                                                    direction = "<", quiet = TRUE))
loglik <- vapply(lp_list, function(z) sum(y * z - pmax(z, 0) - log1p(exp(-abs(z)))), numeric(1))
n_parameters <- vapply(models, function(m) length(coef(m)), integer(1))
metrics <- do.call(rbind, lapply(names(models), function(nm) {
  ci <- as.numeric(pROC::ci.auc(roc_list[[nm]], method = "delong"))
  data.frame(Model = nm, N = length(y), Events = sum(y),
             N_slopes = n_parameters[nm] - 1L, AUC = as.numeric(pROC::auc(roc_list[[nm]])),
             AUC_low = ci[1], AUC_high = ci[3], Brier = mean((prob_list[[nm]] - y)^2),
             Log_loss = -loglik[nm] / length(y), AIC = -2 * loglik[nm] + 2 * n_parameters[nm],
             Evaluation = "Apparent_training", row.names = NULL)
}))
write.csv(metrics, file.path(output_dir, "apparent_model_metrics.csv"), row.names = FALSE)
print(metrics)

lrt <- do.call(rbind, lapply(c("Clinical", "Protein"), function(reduced) {
  statistic <- 2 * (loglik["Joint"] - loglik[reduced])
  df <- n_parameters["Joint"] - n_parameters[reduced]
  data.frame(Comparison = paste("Joint vs", reduced), LR_chisq = unname(statistic),
             df = unname(df), P_value = pchisq(statistic, df = df, lower.tail = FALSE),
             row.names = NULL)
}))
write.csv(lrt, file.path(output_dir, "nested_likelihood_ratio_tests.csv"), row.names = FALSE)
print(lrt)

pairs <- list(c("Joint", "Clinical"), c("Joint", "Protein"), c("Protein", "Clinical"))
auc_comparisons <- do.call(rbind, lapply(pairs, function(pair) {
  test <- pROC::roc.test(roc_list[[pair[1]]], roc_list[[pair[2]]],
                         method = "delong", paired = TRUE, alternative = "two.sided", conf.level = .95)
  delta <- as.numeric(pROC::auc(roc_list[[pair[1]]]) - pROC::auc(roc_list[[pair[2]]]))
  data.frame(Comparison = paste(pair, collapse = " vs "),
             Delta_AUC = delta,
             Delta_AUC_CI_low = as.numeric(test$conf.int[1]),
             Delta_AUC_CI_high = as.numeric(test$conf.int[2]),
             Exploratory_DeLong_P = test$p.value,
             Positive_increment_P_lt_0.05 = delta > 0 && test$p.value < .05,
             Evaluation = "Apparent_training")
}))
print(auc_comparisons)
write.csv(auc_comparisons, file.path(output_dir, "exploratory_AUC_comparisons.csv"), row.names = FALSE)
prediction_table <- data.frame(Source_row = source_rows, Source_rownames = rownames(analysis_dat),
                               ict_failure = y, as.data.frame(prob_list))
write.csv(prediction_table, file.path(output_dir, "three_model_predictions.csv"), row.names = FALSE)

colors <- c(Joint = "#E41A1C", Clinical = "#377EB8", Protein = "#4DAF4A")
display_names <- c(Joint = "Clinical + protein", Clinical = "Clinical", Protein = "Protein panel")
plot_order <- c("Joint", "Clinical", "Protein")
save_plot <- function(draw, name, width = 8, height = 7) {
  for (format in c("pdf", "png", "tiff")) {
    dest <- file.path(output_dir, paste0(name, ".", format))
    if (format == "pdf") pdf(dest, width = width, height = height, useDingbats = FALSE)
    if (format == "png") png(dest, width = width, height = height, units = "in", res = 600)
    if (format == "tiff") tiff(dest, width = width, height = height, units = "in", res = 600,
                               compression = "lzw")
    tryCatch(draw(), finally = dev.off())
  }
}

draw_ROC <- function() {
  op <- par(mar = c(5, 5, 3, 1), family = "sans"); on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC curves (apparent)")
  abline(0, 1, col = "grey65", lty = 2)
  for (nm in plot_order) lines(1 - roc_list[[nm]]$specificities, roc_list[[nm]]$sensitivities,
                               col = colors[nm], lwd = 2)
  auc_labels <- vapply(plot_order, function(nm) {
    r <- metrics[metrics$Model == nm, ]
    sprintf("%s: %.3f (%.3f-%.3f)", display_names[nm], r$AUC, r$AUC_low, r$AUC_high)
  }, character(1))
  format_p <- function(p) {
    if (is.na(p)) "P = NA" else if (p < .001) "P < 0.001" else sprintf("P = %.3f", p)
  }
  increment_labels <- unlist(lapply(c("Joint vs Clinical", "Joint vs Protein"), function(nm) {
    r <- auc_comparisons[auc_comparisons$Comparison == nm, ]
    c(sprintf("%s: dAUC = %.3f (%.3f to %.3f)", nm,
              r$Delta_AUC, r$Delta_AUC_CI_low, r$Delta_AUC_CI_high),
      paste("Paired DeLong", format_p(r$Exploratory_DeLong_P)))
  }))
  legend("bottomright", legend = c(auc_labels, "", increment_labels),
         col = c(unname(colors[plot_order]), rep(NA, 5)),
         lty = c(rep(1, 3), rep(NA, 5)), lwd = c(rep(2, 3), rep(NA, 5)),
         bty = "n", cex = .75, y.intersp = 1.1)
}
draw_ROC()
save_plot(draw_ROC, "ROC_three_models")

risk_dat <- data.frame(ict_failure = y, as.data.frame(prob_list))
dca_list <- setNames(lapply(plot_order, function(nm) {
  rmda::decision_curve(reformulate(nm, response = "ict_failure"), data = risk_dat,
                       fitted.risk = TRUE, policy = "opt-in", study.design = "cohort",
                       thresholds = seq(0, .99, .01), confidence.intervals = NA)
}), plot_order)
draw_DCA <- function() {
  op <- par(mar = c(5, 5, 3, 1), family = "sans"); on.exit(par(op))
  rmda::plot_decision_curve(dca_list, curve.names = unname(display_names[plot_order]),
                            col = c(unname(colors[plot_order]), "black", "grey65"), lty = c(1, 1, 1, 2, 1),
                            lwd = c(2.5, 2, 2, 1.2, 1.2), standardize = FALSE, confidence.intervals = FALSE,
                            cost.benefit.axis = FALSE, xlim = c(0, 1), ylim = c(-.05, mean(y) + .05),
                            xlab = "High Risk Threshold", ylab = "Net Benefit", legend.position = "topright",
                            main = "Decision curves (apparent)")
}
draw_DCA()
save_plot(draw_DCA, "DCA_three_models", width = 9)

saveRDS(list(models = models, data = analysis_dat, source_rows = source_rows,
             predictions = prediction_table, metrics = metrics, LRT = lrt, DCA = dca_list),
        file.path(output_dir, "three_logistic_models.rds"))

# cal_list <- setNames(lapply(names(models), function(nm) {
#   message("Bootstrap calibration: ", nm, " (B = ", B_calibration, ")")
#   set.seed(seed)
#   obj <- rms::calibrate(models[[nm]], method = "boot", B = B_calibration)
#   saveRDS(obj, file.path(output_dir, paste0("calibration_", nm, ".rds")))
#   obj
# }), names(models))
# calibration_metrics <- do.call(rbind, lapply(names(cal_list), function(nm) {
#   obj <- cal_list[[nm]]
#   p <- as.numeric(attr(obj, "predicted"))
#   good <- is.finite(obj[, "predy"]) & is.finite(obj[, "calibrated.corrected"])
#   pc <- approx(obj[good, "predy"], obj[good, "calibrated.corrected"],
#                xout = p, ties = mean)$y
#   data.frame(Model = nm, N = length(p), N_in_calibration_range = sum(is.finite(pc)),
#              Bootstrap_B = attr(obj, "B"), Corrected_MAE = mean(abs(p - pc), na.rm = TRUE))
# }))
# write.csv(calibration_metrics, file.path(output_dir, "bootstrap_calibration_metrics.csv"), row.names = FALSE)
# draw_calibration_comparison <- function() {
#   op <- par(mar = c(5, 5, 3, 1), family = "sans"); on.exit(par(op))
#   plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
#        xlab = "Predicted probability of resistance", ylab = "Calibrated probability of resistance",
#        main = paste0("Bootstrap-corrected calibration (B = ", B_calibration, ")"))
#   abline(0, 1, lty = 2, col = "grey50", lwd = 1.5)
#   for (nm in plot_order) {
#     obj <- cal_list[[nm]]
#     lines(obj[, "predy"], obj[, "calibrated.corrected"], col = colors[nm], lwd = 2)
#   }
#   legend("bottomright", legend = c(unname(display_names[plot_order]), "Ideal"),
#          col = c(unname(colors[plot_order]), "grey50"), lty = c(1, 1, 1, 2),
#          lwd = c(2, 2, 2, 1.5), bty = "n", cex = .9)
# }
# draw_calibration_comparison()
# save_plot(draw_calibration_comparison, "Calibration_three_models")
# print(calibration_metrics)
# 
# 
# 
# library(pROC)
# 
# comparisons <- c(
#   "Joint vs Clinical" = "Clinical",
#   "Joint vs Protein"  = "Protein"
# )
# 
# delong_results <- do.call(
#   rbind,
#   lapply(names(comparisons), function(label) {
#     
#     other <- comparisons[[label]]
#     
#     test <- roc.test(
#       roc_list$Joint,
#       roc_list[[other]],
#       method = "delong",
#       paired = TRUE,
#       alternative = "two.sided",
#       conf.level = 0.95
#     )
#     
#     delta_auc <- as.numeric(
#       auc(roc_list$Joint) - auc(roc_list[[other]])
#     )
#     
#     data.frame(
#       Comparison = label,
#       Delta_AUC = delta_auc,
#       CI_low = as.numeric(test$conf.int[1]),
#       CI_high = as.numeric(test$conf.int[2]),
#       P_value = test$p.value,
#       Positive_increment = delta_auc > 0 & test$p.value < 0.05
#     )
#   })
# )
# 
# print(delong_results)
# 
# write.csv(
#   delong_results,
#   file.path(output_dir, "DeLong_increment_results.csv"),
#   row.names = FALSE
# )





# ============================================================
# LEFTY2 + LAMA2
# ============================================================

library(rms)

# fit_protein：LEFTY2 + LAMA2

fit_target <- fit
model_name <- "Clinical_protein"
B_calibration <- 1000
set.seed(123)

cal_obj <- rms::calibrate(
  fit_target,
  method = "boot",
  B = B_calibration
)

saveRDS(
  cal_obj,
  file.path(
    output_dir,
    paste0("Calibration_", model_name, ".rds")
  )
)



cal_df <- as.data.frame(unclass(cal_obj))

print(colnames(cal_df))

required_cols <- c(
  "predy",
  "calibrated.orig",
  "calibrated.corrected"
)

stopifnot(
  all(required_cols %in% colnames(cal_df))
)

cal_df <- cal_df[
  is.finite(cal_df[["predy"]]) &
    is.finite(cal_df[["calibrated.orig"]]) &
    is.finite(cal_df[["calibrated.corrected"]]),
  ,
  drop = FALSE
]

cal_df <- cal_df[
  order(cal_df[["predy"]]),
  ,
  drop = FALSE
]

head(cal_df)


corrected_y <- approx(
  x = cal_df$predy,
  y = cal_df$calibrated.corrected,
  xout = pred_prob,
  ties = mean
)$y

mae <- mean(
  abs(pred_prob - corrected_y),
  na.rm = TRUE
)

draw_calibration <- function() {
  
  op <- par(
    mar = c(5.8, 4.8, 1.0, 0.8),
    mgp = c(2.7, 0.8, 0),
    family = "sans",
    las = 1,
    cex.axis = 1.05,
    cex.lab = 1.2
  )
  
  on.exit(par(op))
  
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xlab = "Predicted Probability",
    ylab = "Observed Probability",
    bty = "o"
  )
  
  abline(
    a = 0,
    b = 1,
    col = "black",
    lty = 2,
    lwd = 1.3
  )
  
  lines(
    cal_df$predy,
    cal_df$calibrated.orig,
    col = "red",
    lwd = 1.5
  )
  
  lines(
    cal_df$predy,
    cal_df$calibrated.corrected,
    col = "#2166D1",
    lwd = 1.5
  )
  
  rug(
    pred_prob,
    side = 3,
    ticksize = -0.015,
    col = adjustcolor("black", alpha.f = 0.35)
  )
  
  legend(
    "bottomright",
    legend = c(
      "Ideal",
      "Apparent",
      "Bias-corrected"
    ),
    col = c(
      "black",
      "red",
      "#2166D1"
    ),
    lty = c(2, 1, 1),
    lwd = c(1.3, 1.5, 1.5),
    bty = "n",
    cex = 0.95,
    inset = c(0.02, 0.04),
    y.intersp = 1.1
  )
  
  mtext(
    paste0(
      "B= ", B_calibration,
      " repetitions, boot"
    ),
    side = 1,
    line = 4.3,
    adj = 0,
    cex = 0.85
  )
  
  mtext(
    sprintf(
      "Mean absolute error=%.3f n=%d",
      mae,
      n
    ),
    side = 1,
    line = 4.3,
    adj = 1,
    cex = 0.85
  )
}


pred_prob <- as.numeric(
  predict(
    fit_target,
    newdata = analysis_dat,
    type = "fitted"
  )
)

pred_prob <- pred_prob[is.finite(pred_prob)]

n <- length(pred_prob)

corrected_y <- approx(
  x = cal_df[["predy"]],
  y = cal_df[["calibrated.corrected"]],
  xout = pred_prob,
  ties = mean
)$y

mae <- mean(
  abs(pred_prob - corrected_y),
  na.rm = TRUE
)

cat("N =", n, "\n")
cat("Mean absolute error =", round(mae, 3), "\n")

draw_calibration()
# draw_calibration()

tiff(
  filename = file.path(
    output_dir,
    paste0("Calibration_", model_name, ".tiff")
  ),
  width = 8,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

draw_calibration()
dev.off()

pdf(
  file = file.path(
    output_dir,
    paste0("Calibration_", model_name, ".pdf")
  ),
  width = 8,
  height = 6,
  useDingbats = FALSE
)

draw_calibration()
dev.off()

cat(
  "Model:", model_name,
  "\nN:", n,
  "\nBootstrap:", B_calibration,
  "\nMean absolute error:", round(mae, 3),
  "\n"
)














