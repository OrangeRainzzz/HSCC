library(dplyr)
library(tibble)
library(purrr)
library(glmnet)
library(pROC)
hpa_secretome <- read.delim(
  "./data/proteinatlas.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

head(hpa_secretome)
colnames(hpa_secretome)


secreted_protein <- hpa_secretome$Gene[grepl("secreted",hpa_secretome$`Protein class`)]

expr <- read.csv(
  "./data/expr_impute2.csv",
  row.names = 1,
  check.names = FALSE
)
# expr <- read.csv(
#   "./data/expr_input3.csv",
#   row.names = 1,
#   check.names = FALSE
# )

clin <- read.csv(
  "./data/clin.csv",
  row.names = 1,
  check.names = FALSE
)

expr <- expr[, clin$Sample]



## Outcome
clin$drug <- factor(
  clin$drug,
  levels = c("sensitive", "insensitive")
)

## Candidate proteins


up_genes <- read.csv("./differential analysis/input2/all_sample sense vs insense sig.csv")


# up_genes <- common_genes
# up_genes <- up_genes[up_genes$logFC > 1,]
up_genes <- up_genes[up_genes$logFC > 0.5,]
# up_genes <- deg_sig[deg_sig$logFC > 0.5,]

# up_genes <- up_genes[abs(up_genes$logFC) > 1,]
candidate_genes <- up_genes$gene[up_genes$gene %in% secreted_protein]

# candidate_genes <- up_genes$gene[abs(up_genes$logFC) > 1]
# candidate_genes <- up_genes$gene[up_genes$logFC > 1 ]

# candidate_genes <- common_genes
# candidate_genes <- c('SFRP4','COL2A1','SFRP2')
length(candidate_genes)
# up_genes
# candidate_genes <- up_genes$gene
## Training and validation samples
train_ids <- intersect(
  clin$Sample[clin$group == "re" |clin$group == "pro" ],
  colnames(expr)
)

test_ids <- intersect(
  clin$Sample[clin$group == "re" |clin$group == "pro" ],
  colnames(expr)
)

train_clin <- clin[
  match(train_ids, clin$Sample),
]

test_clin <- clin[
  match(test_ids, clin$Sample),
]

stopifnot(
  all(train_clin$Sample == train_ids),
  all(test_clin$Sample == test_ids)
)



library(tidyverse)

expr_train_raw <- expr_raw[
  candidate_genes,
  train_ids,
  drop = FALSE
]

sensitive_ids <- train_clin$Sample[
  train_clin$drug == "sensitive"
]

insensitive_ids <- train_clin$Sample[
  train_clin$drug == "insensitive"
]

missing_result <- tibble(
  gene = candidate_genes,
  
  missing_all = rowMeans(
    is.na(expr_train_raw)
  ),
  
  missing_sensitive = rowMeans(
    is.na(
      expr_train_raw[
        ,
        sensitive_ids,
        drop = FALSE
      ]
    )
  ),
  
  missing_insensitive = rowMeans(
    is.na(
      expr_train_raw[
        ,
        insensitive_ids,
        drop = FALSE
      ]
    )
  )
) %>%
  arrange(missing_all)

head(missing_result, 20)




keep_genes <- missing_result %>%
  filter(
    missing_all <= 0.25,
    # missing_sensitive <= 0.25,
    # missing_insensitive <= 0.25
  ) %>%
  pull(gene)


keep_genes <- candidate_genes

length(keep_genes)
keep_genes
# keep_genes <- common_genes

x_train <- t(
  expr[
    keep_genes,
    train_ids,
    drop = FALSE
  ]
)

x_test <- t(
  expr[
    keep_genes,
    test_ids,
    drop = FALSE
  ]
)



y_train <- as.integer(
  train_clin$drug == "insensitive"
)

test_clin <- clin[
  match(test_ids, clin$Sample),
]

y_test <- as.integer(
  test_clin$drug == "insensitive"
)





run_single_logistic <- function(gene) {
  
  dat <- data.frame(
    outcome = y_train,
    protein = x_train[, gene],
    Smoke = factor(train_clin$Smoke)
  )
  
  fit <- tryCatch(
    glm(
      # outcome ~ protein + Smoke,
      outcome ~ protein,
      data = dat,
      family = binomial()
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(
      data.frame(
        gene = gene,
        beta = NA,
        OR = NA,
        lower95 = NA,
        upper95 = NA,
        pvalue = NA
      )
    )
  }
  
  coef_table <- summary(fit)$coefficients
  
  if (!"protein" %in% rownames(coef_table)) {
    return(
      data.frame(
        gene = gene,
        beta = NA,
        OR = NA,
        lower95 = NA,
        upper95 = NA,
        pvalue = NA
      )
    )
  }
  
  beta <- coef_table[
    "protein",
    "Estimate"
  ]
  
  se <- coef_table[
    "protein",
    "Std. Error"
  ]
  
  data.frame(
    gene = gene,
    beta = beta,
    OR = exp(beta),
    lower95 = exp(beta - 1.96 * se),
    upper95 = exp(beta + 1.96 * se),
    pvalue = coef_table[
      "protein",
      "Pr(>|z|)"
    ]
  )
}

single_result <- do.call(
  rbind,
  lapply(
    keep_genes,
    run_single_logistic
  )
)

single_result$FDR <- p.adjust(
  single_result$pvalue,
  method = "BH"
)

single_result <- single_result[
  order(single_result$FDR),
]





library(glmnet)

set.seed(123)



single_result_valid <- single_result[single_result$pvalue < 0.05 ,]
single_result_valid <- single_result_valid[order(single_result_valid$pvalue,decreasing = F)  ,]


candidate_panel_genes <- single_result$gene

candidate_panel_genes <- single_result_valid$gene
# candidate_panel_genes <- single_result_valid$gene[1:5]

# candidate_panel_genes <- keep_genes
length(candidate_panel_genes)
candidate_panel_genes

library(glmnet)
library(dplyr)



x_lasso <- as.matrix(
  x_train[, candidate_panel_genes, drop = FALSE]
)

set.seed(123)

cv_lasso <- cv.glmnet(
  x = x_lasso,
  y = y_train,
  family = "binomial",
  alpha = 1,
  nfolds = 5,
  type.measure = "deviance",
  standardize = TRUE,
  keep = TRUE
)

selected_lambda <- cv_lasso$lambda.min

selected_lambda

coef_mat <- as.matrix(
  coef(
    cv_lasso,
    s = "lambda.1se"
  )
)

coef_table <- data.frame(
  protein = rownames(coef_mat),
  coefficient = as.numeric(coef_mat[, 1]),
  row.names = NULL
) %>%
  filter(protein != "(Intercept)") %>%
  mutate(
    status = ifelse(
      abs(coefficient) > 1e-10,
      "Non-zero",
      "Zero"
    ),
    direction = case_when(
      coefficient > 1e-10  ~ "Positive",
      coefficient < -1e-10 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) %>%
  arrange(desc(abs(coefficient)))

coef_table

final_panel <- coef_table %>%
  filter(status == "Non-zero")

panel_genes <- final_panel$protein

length(panel_genes)
panel_genes
final_panel


