write_paper : get_references pipeline_dev data_wrangling

pipeline_dev : mothur_raw.bash get_references
	sh mothur_raw.bash

data_wrangling : data_wrangling.R mothur_raw.bash get_references
	R CMD BATCH data_wrangling.R

get_references : get_references.bash
	sh get_references.bash
