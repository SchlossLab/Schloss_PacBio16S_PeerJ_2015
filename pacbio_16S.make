#write_paper
sh get_references.bash
sh mothur_raw.bash
R CMD BATCH data_wrangling.R
sh mothur_process.bash
