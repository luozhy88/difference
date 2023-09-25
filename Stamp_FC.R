library(dplyr)
library(tibble)
library(DescTools)
library(purrr)
library(ggplot2)
library(patchwork)
dir.create("output")
rm(list = ls())
# function

pathway_pdf=function(Meta_data=Meta_data,meta_col=meta_col,out.name=out.name,height=height,padj=NULL,pvalue=NULL,width=width){
  
  ####args
  # The rowname of Meta_data is sampleID;
  # the column of Meta_data is Features;
  # the first column of Meta_data is Groups;

  ###
  DF <-  Meta_data %>% dplyr::select(meta_col, everything())
  DF_features_col<-colnames(DF)[-c(1:length(meta_col))]
  
  ##
  table(DF$Group)
  
  DF$Group <- as.character(DF$Group)
  ### create all combinations
  # truncate_label <- function(x, n) {
  #   str_trunc(x, width = n)
  # }
  pair <- CombPairs(unique(DF$Group))
  
  for (j in 1:nrow(pair)) {
    # j=1
    print(j)
    DF2 <- subset(DF, Group == pair[j, 1] | Group == pair[j, 2])
    
    comparison <- paste0(pair[j, 1], "_", pair[j, 2])
    
    
    ## or ignore the subgroups
    DF.tss <- DF2
    
    # 计算每列的方差
    DF.tss_var <- DF.tss%>% 
      select_if(is.numeric)       %>%
      select(where(~ var(.x) > 1e-12))
    DF.tss=DF.tss %>% select(1, colnames(DF.tss_var))
    
    
    # DF.tss <- DF2 %>%janitor::adorn_percentages("row") %>% data.frame()
    DF.tss <- DF.tss[, colMeans(is.na(DF.tss)) <= 0.3]
    DF.tss_features_col =colnames(DF.tss)
    
    # Log2 Fold Change
    groupA_data <- subset(DF.tss, Groups == pair[j, 1])
    groupB_data <- subset(DF.tss, Groups == pair[j, 2])
    
    # Calculate the log2 fold change for each feature
    log2FC <- log2(colMeans(groupB_data[, -c(1, 2)]) / colMeans(groupA_data[, -c(1, 2)]))
    
    # Add log2 fold change column to the data frame
    lfc <- as.data.frame(log2FC)
    
    lfc <- tibble::rownames_to_column(lfc, var = "var")
    
    lfc$comparison_lfc <- paste0(pair[j, 2], "-", pair[j, 1])
    
    # Calculate mean differences with confidence intervals
    diff.real <- DF.tss %>%
      dplyr::select_if(is.numeric) %>%
      purrr::map_df(~ DescTools::MeanDiffCI(. ~ Groups, data = DF.tss), .id = "var")
    
    # Perform Wilcoxon rank-sum test and tidy the results
    diff.model <- DF.tss %>%
      dplyr::select_if(is.numeric) %>%
      purrr::map_df(~ broom::tidy(
        wilcox.test(
          . ~ Groups,
          data = DF.tss,
          conf.int = T,
          conf.level = .95
        )
      ), .id = "var")
    
    # Adjust p-values using Benjamini-Hochberg method
    diff.model$p.adj <- p.adjust(diff.model$p.value, "BH")
    
    # Merge mean difference and test results
    diff <- dplyr::left_join(diff.model, diff.real)
    
    diff <- dplyr::left_join(diff, lfc)
    diff$Log2fc=abs(diff$log2FC)
    
    if( is.null(pvalue)  ){
      print("Use:padj")
      diff.sig <- subset(diff, p.adj <= pvalue)
      out.name=paste0(out.name,"padj0.05")
      }else{
      print("Use:pvalue")
      diff.sig <- subset(diff, pvalue <= pvalue)
      out.name=paste0(out.name,"pvalue0.05")
      }
    write.csv(diff,glue::glue("output/", out.name, comparison, "_wilcox.csv"))
    diff.sig=diff.sig %>%dplyr::arrange(desc(log2FC))
    
    if (nrow(diff.sig) != 0) {
      # Add color column based on mean difference sign
      diff.sig$col <-
        ifelse(diff.sig$meandiff > 0, "01plus", "02minus")
      
      DF.tss2 <- reshape2::melt(  DF.tss[,  c("Groups", DF.tss_features_col) ]  )
      #  DF.tss2$Group <- factor(DF.tss2$Group, levels = c("young", str_replace(comparison, pattern = "young-vs-", replacement = "")))
      DF.tss2$variable <- as.character(DF.tss2$variable)
      ##
      DF.tss3 <- DF.tss2[DF.tss2$variable %in% diff.sig$var, ]
      
      level_order <- diff.sig %>%
        dplyr::arrange(log2FC) %>%
        dplyr::pull(var)
      ## base plot
      p1 <- ggplot(data = DF.tss3, aes(x = variable, y = value, fill = Groups))
      
      ### add alternated background
      for (i in 1:(nrow(diff.sig))) {
        p1 <- p1 + annotate(
          "rect",
          xmin = i - 0.5,
          xmax = i + 0.5,
          ymin = 0,
          ymax = Inf,
          fill = ifelse(i %% 2 == 0, "white", "gray95")
        ) + scale_x_discrete()
      }
      p1 <- p1 +
        geom_bar(
          color = "black",
          position = "dodge",
          stat = "summary",
          fun = "mean"
        ) +
        xlab("") +
        ylab("Proportion (%)") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(color = "black"),
          axis.ticks = element_line(color = "black", linewidth = 0.5)
        ) +
        scale_x_discrete(
          limits = level_order
          #   labels = function(x) {
          #     truncate_label(x, 40)
          #   }
        ) +
        coord_flip() +
        theme(legend.position = "top")
      
      # Set custom fill colors
      p1 <-
        p1 + scale_fill_manual(values = c("darkorange", "deepskyblue"))
      
      # Base plot for scatter plot
      p2 <- ggplot(diff.sig, aes(
        x = reorder(var, log2FC),
        y = log2FC,
        fill = col
      ))
      
      # Add alternated background to scatter plot
      for (i in 1:(nrow(diff.sig))) {
        p2 <- p2 + annotate(
          "rect",
          xmin = i - 0.5,
          xmax = i + 0.5,
          ymin = -Inf,
          ymax = Inf,
          fill = ifelse(i %% 2 == 0, "white", "gray95")
        ) + scale_x_discrete()
      }
      ## add the other components
      p2 <- p2 +
        geom_point(stat = "identity", size = 3, aes(col = col), position = position_dodge(width = 1)) +
        # geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci),
        #               width = .2, size = 0.5,
        #               position = "dodge"
        # ) +
        xlab("") +
        ylab("log2FC") +
        geom_hline(yintercept = 0, lty = 2) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(color = "black"),
          axis.ticks.x = element_line(color = "black", linewidth = 0.5)
        ) +
        coord_flip() +
        theme(legend.position = "none") +
        theme(axis.text.y = element_blank())
      p2 <- p2 + scale_color_manual(values = c("darkorange", "deepskyblue"))
      
      p1 + p2
      
      ###
      ##
      dir.create("output",recursive = T)
      # height=2
      # width=12
      ggsave(filename = glue::glue("output/", out.name, comparison, "_metagenomics.pdf"), width = width, height = height)
      
      #  }
    }
  }
}

############################################################################################################


##meta
load("../../../DP_TOFU_shotgun_HM32_20230925/input/human3/human3_DP_Meta_cazy.rData")
Meta_data=Meta_cazy


## Cazy
CLass="Cazy"
pathway_pdf(Meta_data=Meta_cazy,meta_col=c("Groups"),out.name= glue::glue("Human3_class_" ,CLass %>% make.names(),"_"),height=15,pvalue=0.05,width=15)

# Meta_data=Meta_cazy
# meta_col=c("Groups")
# out.name=glue::glue("Human3_class_" ,CLass %>% make.names(),"_")
# height=10
# pvalue=0.05
# width=15




