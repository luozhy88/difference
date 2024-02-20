   library(tidyverse)

   ## data_6 is a data frame with two columns named "sample" and "group", and all the
   ## other columns are features

   ## select top variance
    
    data_6 %>%
      pivot_longer(-c(sample, group)) %>%
      group_by(name) %>%
      summarise(variance = var(value)) %>%
      arrange(desc(variance)) %>%
      top_frac(0.1) %>% 
      pull(name) -> top_features
    
    data_6 <- data_6 %>%
      select(sample, group, all_of(top_features))




  # Perform Wilcoxon rank-sum test and tidy the results
    diff.model <- data_6 %>%
      dplyr::select_if(is.numeric) %>%
      purrr::map_df(~ broom::tidy(
        wilcox.test(
          . ~ group,
          data = data_6,
          conf.int = T,
          conf.level = .95
        )
      ), .id = "var")

    # Adjust p-values using Benjamini-Hochberg method
    diff.model$p.adj <- p.adjust(diff.model$p.value, "BH")
    
