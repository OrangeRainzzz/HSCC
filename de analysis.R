


####all ####

expr <- read.csv("./data/expr_select2.csv",row.names = 1)
clin <- read.csv("./data/clin.csv",row.names = 1)

#####smoke vs no#####

# colnames(expr)[1:5]
# clin$Sample[1:5]

# clin <- clin[clin$group == "re",]

expr <- expr[, clin$Sample]
expr <- as.matrix(expr)

# clin$Smoke <- factor(
#   clin$Smoke,
#   levels = c("No", "Current_smoke")
# )
# 
# clin$drug <- factor(
#   clin$drug,
#   levels = c("sensitive", "insensitive")
# )


group <- factor(clin$Smoke, levels = c("No", "Current_smoke"))

design <- model.matrix(
  ~ Smoke +drug,
  data = clin
)
colnames(design)

fit <- lmFit(expr, design)
fit2 <- eBayes(fit)

deg <- topTable(
  fit2,
  coef = "SmokeCurrent_smoke",
  number = Inf,
  adjust.method = "BH"
)

deg$gene <- rownames(deg)

head(deg)


deg$change <- "NS"
deg$change[deg$P.Value < 0.05 & deg$logFC > 1]  <- "Up"
deg$change[deg$P.Value < 0.05 & deg$logFC < -1] <- "Down"

table(deg$change)

# deg <- read.csv("./differential analysis/input2/impute/re_sample Smoke vs No.csv")
deg_sig <- deg[deg$P.Value < 0.05,]
deg_sig <- deg_sig[order(deg_sig$logFC, decreasing = TRUE),]
# deg_sig <- deg_sig[!is.na(deg_sig$X),]
deg_sig <- deg_sig[  abs(deg_sig$logFC)  > 0.5 ,]

write.csv(deg_sig,"./differential analysis/input2/all_sample Smoke vs No sig.csv")
write.csv(deg,"./differential analysis/input2/all_sample Smoke vs No.csv")




#####sense vs non ####


expr <- read.csv("./data/expr_select2.csv",row.names = 1)
clin <- read.csv("./data/clin.csv",row.names = 1)

# clin <- clin[clin$group == "re",]

expr <- expr[, clin$Sample]
expr <- as.matrix(expr)

clin$Smoke <- factor(
  clin$Smoke,
  levels = c("No", "Current_smoke")
)

clin$drug <- factor(
  clin$drug,
  levels = c("sensitive", "insensitive")
)


design <- model.matrix(
  ~ drug+Smoke ,
  data = clin
)
colnames(design)



fit <- lmFit(expr, design)
fit2 <- eBayes(fit)

deg <- topTable(
  fit2,
  coef = "druginsensitive",
  number = Inf,
  adjust.method = "BH"
)

deg$gene <- rownames(deg)

head(deg)


deg$change <- "NS"
deg$change[deg$P.Value < 0.05 & deg$logFC > 1]  <- "Up"
deg$change[deg$P.Value < 0.05 & deg$logFC < -1] <- "Down"

table(deg$change)




write.csv(deg,"./differential analysis/input2/all_sample sense vs insense.csv")
# deg <- read.csv("./differential analysis/input2/re_sample sense vs insense.csv")

deg_sig <- deg[deg$P.Value < 0.05,]
deg_sig <- deg_sig[order(deg_sig$logFC,decreasing = T),]
# deg_sig <- deg_sig[!is.na(deg_sig$gene),]
deg_sig <- deg_sig[  abs(deg_sig$logFC)  > 0.5 ,]
library(ggplot2)
write.csv(deg_sig,"./differential analysis/input2/all_sample sense vs insense sig.csv",row.names = F)



##### smoke insense vs non-smoke sense ####

expr <- read.csv("./data/expr_impute2.csv",row.names = 1)
clin <- read.csv("./data/clin.csv",row.names = 1)
clin <- clin[clin$group == "re",]

clin$group2 <- paste0(clin$Smoke,"_",clin$drug)
buffer <- clin[clin$group2 %in% c("No_sensitive","Current_smoke_insensitive"),]
table(buffer$group2)



expr <- expr[, buffer$Sample]
expr <- as.matrix(expr)


group <- factor(clin$group2, levels = c("Current_smoke_insensitive", "No_sensitive"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr, design)
cont.matrix <- makeContrasts(Current_smoke_insensitive - No_sensitive, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(fit2, number = Inf, adjust.method = "BH")
deg$gene <- rownames(deg)

head(deg)


deg$change <- "NS"
deg$change[deg$P.Value < 0.05 & deg$logFC > 1]  <- "Up"
deg$change[deg$P.Value < 0.05 & deg$logFC < -1] <- "Down"

table(deg$change)




write.csv(deg,"./differential analysis/input2/re_sample smoke_insense vs non_smoke_sense.csv")
# deg <- read.csv("./differential analysis/input2/re_sample sense vs insense.csv")

deg_sig <- deg[deg$P.Value < 0.05,]
deg_sig <- deg_sig[order(deg_sig$logFC,decreasing = T),]
# deg_sig <- deg_sig[!is.na(deg_sig$gene),]
deg_sig <- deg_sig[  abs(deg_sig$logFC)  > 0.5 ,]
library(ggplot2)
write.csv(deg_sig,"./differential analysis/input2/re_sample smoke_insense vs non_smoke_sense sig.csv",row.names = F)

##### smoke inse####


expr <- read.csv("./data/expr_select2.csv",row.names = 1)
clin <- read.csv("./data/clin.csv",row.names = 1)


# colnames(expr)[1:5]
# clin$Sample[1:5]

clin <- clin[clin$group == "re",]

expr <- expr[, clin$Sample]
expr <- as.matrix(expr)

clin$Smoke <- factor(
  clin$Smoke,
  levels = c("No", "Current_smoke")
)

clin$drug <- factor(
  clin$drug,
  levels = c("sensitive", "insensitive")
)



group <- factor(clin$Smoke, levels = c("No", "Current_smoke"))

design <- model.matrix(
  ~ drug * Smoke,
  data = clin
)
colnames(design) <- c(
  "Intercept",
  "Drug",
  "Smoke",
  "Interaction"
)

design
colnames(design)
#[1] "Intercept"   "Drug"        "Smoke"       "Interaction"
fit <- lmFit(expr, design)
fit2 <- eBayes(fit)

deg <- topTable(
  fit2,
  coef = "Interaction",
  number = Inf,
  adjust.method = "BH"
)

deg$gene <- rownames(deg)

head(deg)


deg$change <- "NS"
deg$change[deg$P.Value < 0.05 & deg$logFC > 1]  <- "Up"
deg$change[deg$P.Value < 0.05 & deg$logFC < -1] <- "Down"

table(deg$change)

# deg <- read.csv("./differential analysis/input2/impute/re_sample Smoke vs No.csv")
deg_sig <- deg[deg$P.Value < 0.05,]
deg_sig <- deg_sig[order(deg_sig$logFC, decreasing = TRUE),]
# deg_sig <- deg_sig[!is.na(deg_sig$X),]
deg_sig <- deg_sig[  abs(deg_sig$logFC)  > 0.5 ,]

write.csv(deg_sig,"./differential analysis/input2/re_sample interaction sig.csv")
write.csv(deg,"./differential analysis/input2/re_sample interaction.csv")












library(org.Hs.eg.db)
library(clusterProfiler)
files <- list.files("./differential analysis/input2/")
files <- files[grepl("sig",files)]
files_save <- gsub(".csv",'',files)


files_length <- length(files)

for (i in 1:files_length) {
  
  file_select <- files[i]
  files_save_select <- files_save[i]
  
  file_dir <- paste0("./differential analysis/input2/",file_select)
  
  buffer <- read.csv(file_dir)
  
  # buffer <- read.csv("./differential analysis/logfc1/re_sample sense vs insense sig.csv")
  
  # upgene1 <- buffer$gene
  upgene <- buffer$gene[buffer$logFC  > 1]
  dwgene <- buffer$gene[buffer$logFC  < -1]
  
  
  
  
  if (length(upgene) > 0) {
    
    genelist_enrich <- bitr(
      upgene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    if (!is.null(genelist_enrich) && nrow(genelist_enrich) > 0) {
      
      gobp <- enrichGO(
        gene = genelist_enrich$ENTREZID,
        universe = universe_df$ENTREZID,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        readable = TRUE
      )
      
      gobp_res <- gobp@result
      gobp_res <- gobp_res[gobp_res$pvalue < 0.05, ]
      
      write.csv(
        gobp_res,
        paste0(
          "./differential analysis/input2/enrich/",
          files_save_select,
          " up GOBP.csv"
        ),
        row.names = FALSE
      )
    }
    
  } else {
    message(files_save_select, ": 没有上调基因，跳过 up GOBP")
  }
  
  
  
  if (length(dwgene) > 0) {
    
    genelist_enrich <- bitr(
      dwgene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    if (!is.null(genelist_enrich) && nrow(genelist_enrich) > 0) {
      
      gobp <- enrichGO(
        gene = genelist_enrich$ENTREZID,
        universe = universe_df$ENTREZID,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        readable = TRUE
      )
      
      gobp_res <- gobp@result
      gobp_res <- gobp_res[gobp_res$pvalue < 0.05, ]
      
      write.csv(
        gobp_res,
        paste0(
          "./differential analysis/input2/enrich/",
          files_save_select,
          " dw GOBP.csv"
        ),
        row.names = FALSE
      )
    }
    
  } else {
    message(files_save_select, ": 没有下调基因，跳过 dw GOBP")
  }
}





file1 <- read.csv("./differential analysis/input2/enrich/re_sample sense vs insense sig up GOBP.csv")
file2 <- read.csv("./differential analysis/input2/enrich/re_sample smoke_insense vs non_smoke_sense sig up GOBP.csv")





