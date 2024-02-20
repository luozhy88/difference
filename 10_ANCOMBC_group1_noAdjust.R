library(microbiomeMarker)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
library(mia)
library(ANCOMBC)


# Create a directory named "ancombc"
dir.create("../output/differential/ancombc2_group1_noAdjust", recursive = TRUE)
path <- "../output/differential/ancombc2_group1_noAdjust/"


# Read the phyloseq object from a file and filter out taxa with zero counts
pseq <- readRDS("../../ampliseq_output_80subjects/phyloseq/dada2_phyloseq.rds")


sample_data(pseq)$age <- scale(sample_data(pseq)$age)
###############################################################
############### empty values for ranks replaced with "undefined"##
###############################################################
# Extract the tax_table from the phyloseq object
taxa <- as.data.frame(tax_table(pseq))

# Replace empty or NA values with "undefined"
taxa[is.na(taxa)] <- "undefined"
taxa[taxa == ""] <- "undefined"

# Convert the modified data frame back to a tax_table format
tax_table_modified <- tax_table(as.matrix(taxa))

# Assign the modified tax_table back to the phyloseq object
tax_table(pseq) <- tax_table_modified
###############################################################
###############################################################

pseq <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)
pseq <- phyloseq::prune_samples(sample_sums(pseq) > 0, pseq)


###### group1 ---> group

names(sample_data(pseq))[1] <- "group"

## GROUP_VAR <- paste0("sex+age+", "group")
GROUP_VAR <- paste0("group")

## here we refer to the main interest even though one has to add several covariates like above
MAIN <- trimws(GROUP_VAR, whitespace = ".*\\+")

### here I rename all names (unreliable at all) with ASVs. But we must know they are ASVs not Species
# ## better to make sure there are no improper names (empty) at Genus level
phyloseq::tax_table(pseq)[, "Species"] <- sprintf("%s%04d", "ASV", seq_along(taxa_names(pseq)))

############# split the phyloseq into a list based on the variable of interest
#######
agegut_split <- metagMisc::phyloseq_sep_variable(pseq, MAIN)
############ make pairwise comparison
pair <- combn(names(agegut_split), 2)


for (i in 1:(length(pair) / 2)) {
  ## new phyloseq
  pseq.sub <- phyloseq::merge_phyloseq(agegut_split[[pair[1, i]]], agegut_split[[pair[2, i]]])


  pseq.temp <- core(pseq.sub, detection = 1, prevalence = .1)

  ##
  for (my_rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    pseq.new <- phyloseq::tax_glom(pseq.temp, taxrank = my_rank)
    tax_table(pseq.new)[, my_rank] <- make.names(tax_table(pseq.new)[, my_rank])

    fil.names <- paste0(unique(microbiome::meta(pseq.new)[, MAIN])[1], "_vs_", unique(microbiome::meta(pseq.new)[, MAIN])[2])

    tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(pseq.new)
    sum_res_ancombc <- ANCOMBC::ancombc2(
      data = tse, assay_name = "counts", tax_level = my_rank,
      fix_formula = GROUP_VAR, rand_formula = NULL,
      p_adj_method = "fdr",
      prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
      group = MAIN, struc_zero = TRUE, neg_lb = TRUE,
      alpha = 0.05, n_cl = 2, verbose = TRUE,
      global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
      iter_control = list(
        tol = 1e-2, max_iter = 20,
        verbose = TRUE
      ),
      em_control = list(tol = 1e-5, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
      trend_control = list(
        contrast = list(
          matrix(c(1, 0, -1, 1),
            nrow = 2,
            byrow = TRUE
          ),
          matrix(c(-1, 0, 1, -1),
            nrow = 2,
            byrow = TRUE
          )
        ),
        node = list(2, 2),
        solver = "ECOS",
        B = 100
      )
    )

    tab_zero <- sum_res_ancombc$zero_ind

    sum_res <- sum_res_ancombc$res

    TwoGroups <- unique(microbiome::meta(pseq.new)[, MAIN])
    
    
    match <- colnames(sum_res)[grep("diff_group", colnames(sum_res))]

    REF <- gsub(pattern = "diff_group", replacement = "", match)
    
    notREF <- TwoGroups[TwoGroups != REF]

    colnames(sum_res) <- gsub("_group.*", "", colnames(sum_res))

    
    fil.names.new <- paste0(fil.names, "_***", REF, "_as_reference***")
    
    openxlsx::write.xlsx(tab_zero, file = paste0(path, fil.names.new, my_rank, "__ancombc_structural_zero.xlsx"))
    
    
    sum_res$group <- ifelse(sum_res$lfc < 0, "REF", "notREF")

    rownames(sum_res) <- sum_res$taxon
    ## we can take only the significantly differential taxa
    sig_ancombc <- subset(sum_res, q < 0.05)
    # Adding taxonomic labels
    taxa_info <- data.frame(tax_table(pseq.new))
    rownames(taxa_info) <- make.unique(taxa_info[, my_rank])

    sum_res <- merge(sum_res, taxa_info, by = 0)
    sig_ancombc <- merge(sig_ancombc, taxa_info, by = 0)

    ##
    ## save the output for each comparison, reference shown on the file name
    openxlsx::write.xlsx(sum_res, file = paste0(path, fil.names.new, my_rank, "__all_ancombc_cal.xlsx"))
    openxlsx::write.xlsx(sig_ancombc, file = paste0(path, fil.names.new, my_rank, "__sig_ancombc_cal.xlsx"))

    ## reformat the output for barplot
    sig_ancombc$q <-
      ifelse(sig_ancombc$q <= 0.05,
        "< 0.05",
        ifelse(sig_ancombc$q <= 0.1 & sig_ancombc$q > 0.05,
          "0.05 - 0.1",
          "0.1 - 0.2 "
        )
      )
    ##

    ### do not need to plot if there were no significantly differential taxa
    if (nrow(sig_ancombc) > 0) {
      ##
      sig_ancombc$taxa <- sig_ancombc$Row.names

      sig_ancombc$class <- ifelse(sig_ancombc$group == "REF", REF, notREF)

      ## one can also plot the results based on log fold change LFC but colored by padj range
      ggpubr::ggbarplot(sig_ancombc,
        x = "taxa", y = "lfc",
        fill = "class",
        color = "white",
        palette = c("#DD498D", "#00A4D3"),
        sort.val = "desc",
        sort.by.sex = FALSE,
        x.text.angle = 90,
        ylab = "Log2 Fold Change",
        xlab = "",
        rotate = TRUE,
        title = paste0("[", notREF, "]", "-", "[", REF, "]"),
        ggtheme = theme_minimal()
      )

      ## save the plot
      ##
      ggsave(filename = paste0(path, fil.names.new, my_rank, "__ancombc_cal_barplot.pdf"), width = 6, height = nrow(sig_ancombc) * 0.1 + 2)
    }
  }
}
