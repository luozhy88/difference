# 生成分组数据
groups <- rep(c("G1", "G2"), each = 25)
gender <- sample(c("Male", "Female"), 50, replace = TRUE)
age <- sample(18:60, 50, replace = TRUE)
data_df <- data.frame(Groups = groups, Gender = gender, Age = age)


# tableone

library(tibble)
library(tidyverse)
library(tableone)
library(table1)
`%+%` <- function(a,b) {paste0(a,b)}
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  # print(x)
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}


f_2 <- as.formula("~ " %+%   paste(colnames(data_df)[-c(1)],collapse = " + ")  %+% " | " %+% "Groups"  )
f_2


gender_cha=""
w<-table1(f_2, data=data_df, overall=F, extra.col=list(`P-value`=pvalue),footnote=gender_cha)
w
dir.create("output")
htmltools::save_html(w, file ="output/"  %+%  "groups_tableone.html"  )
