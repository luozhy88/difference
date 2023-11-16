

library(ggpubr)
load("../01_STAMP_Groups/input/Metat.rds")

###################################################### metacyc
data_df=read.csv("../01_STAMP_Groups/input/metacyc_data.spf",sep = "\t") %>% tibble::column_to_rownames("metacyc_name") %>% t() %>% data.frame()
data_df <- log2(data_df + 0.00001)

data_df_meta=merge(Metat,data_df,by=0) %>%  tibble::column_to_rownames("Row.names") %>% dplyr::select(-sampleID)
data_df_meta$Groups = factor(data_df_meta$Groups,levels = c("HC","DP")) #Groups Age Sex  BMI+吸烟

## linear regression (adjust sex age BMI+吸烟)
groups =data_df_meta$Groups %>% unique()
dir.create(paste0(groups, collapse = "_vs_"))
lm_fit <- list()
formula="~ Groups + Age + Sex +BMI+吸烟"
for (i in colnames(data_df_meta)[-c(1:5)]) {
  lm_fit[[i]] <- lm(
    formula = as.formula(paste0(i, formula)),
    data = data_df_meta
  )
}

lm_list <- lapply(lm_fit, model_parameters)
lm_dat <- do.call(rbind, lm_list)
lm_res <- subset(lm_dat, Parameter == "GroupsDP")
lm_res <- rownames_to_column(lm_res, var = "features")
lm_res$features <- gsub(".2$", "", lm_res$features)
lm_res$p.adj <- p.adjust(lm_res$p, method = "fdr")
openxlsx::write.xlsx(lm_res, "HC_vs_DP/immune_lm_sex_age_BMI_smoking_adjusted.xlsx")

## Plot

lm_res$coeff <- ifelse(lm_res$Coefficient > 0, "Increase", "Decresase")

data2 <- subset(lm_res, p < 0.01)

bar.plot <- ggbarplot(data2,
                      x = "features", y = "Coefficient",
                      fill = "coeff",
                      color = "white",
                      # palette = "",
                      sort.val = "asc",
                      sort.by.groups = FALSE,
                      x.text.angle = 90,
                      rotate = TRUE,
                      ylab = "Coefficient",
                      xlab = "",
                      title = paste0(groups, collapse = "_vs_"),
                      palette = "aaas",
                      ggtheme = theme_minimal(),
                      legend.title = ""
)

ggsave(plot = bar.plot,filename  = paste0(paste0(groups, collapse = "_vs_"), "/metacyc_barplot_p001_age_sex_BMI_smoking_adjusted.pdf"), width =10, height = 5)
