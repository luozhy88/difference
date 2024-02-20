
############################lefse docker ###############################
library(microbiome)
library(dplyr)
######args:you need the column of sample and group in meta for phyloseq.
pseq <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)
pseq <- phyloseq::prune_samples(sample_sums(pseq) > 0, pseq)

Meta=meta(pseq) %>% dplyr::select(sample,group) %>% data.frame() %>% t()
data_df=pseq@otu_table  %>% data.frame()

colnames(Meta)== colnames(data_df)

lefse_input=rbind(Meta,data_df)
rownames(lefse_input)[2]="class"
write.table(lefse_input,"lefse_input.txt",row.names = T,col.names = F,quote = F,sep = "\t")

bash.txt_str=glue::glue("docker run -it -v ",getwd(),":/home/linuxbrew/.cache biobakery/lefse bash",
"
cd /home/linuxbrew/.cache
format_input.py lefse_input.txt  lefse_input.in -c 2  -u 1 -o 1000000
run_lefse.py lefse_input.in lefse_input.res
plot_res.py lefse_input.res lefse_input.png
")

writeLines(bash.txt_str,"lefse_docker_bash.sh")

