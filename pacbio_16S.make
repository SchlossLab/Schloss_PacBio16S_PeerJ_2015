write_paper : got_references pipeline_dev data_wrangling

got_references : get_references.bash
	sh get_references.bash

pipeline_dev : mothur_raw.bash got_references
	sh mothur_raw.bash

data_wrangling : data_wrangling.R mothur_raw.bash got_references
	R CMD BATCH data_wrangling.R
