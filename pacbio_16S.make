write_paper : get_references pipeline_dev data_wrangling

get_data : raw_data/* HMP_MOCK.fasta
	sh get_data.sh 

pipeline_dev : mothur_raw.bash get_data
	sh mothur_raw.bash

data_wrangling : data_wrangling.R mothur_raw.bash get_data
	R CMD BATCH data_wrangling.R

