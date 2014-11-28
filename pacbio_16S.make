write_paper : get_references pipeline_dev data_wrangling

pipeline_dev : mothur_raw.bash get_references
	sh mothur_raw.bash

data_wrangling : data_wrangling.R mothur_raw.bash get_references
	R CMD BATCH data_wrangling.R

get_references : get_references.bash references/got.gg references/got.rdp references/got.silva references/HMP_MOCK.align

references/got.gg : get_references.bash references/gg_13_8_99.fasta references/gg_13_8_99.gg.tax
	sh get_references.bash

references/got.silva : get_references.bash references/silva.bacteria.fasta references/silva.bacteria.tax
	sh get_references.bash

references/got.rdp : get_references.bash references/trainset10_082014.pds.fasta references/trainset10_082014.pds.tax
	sh get_references.bash

references/HMP_MOCK.align : get_references.bash references/HMP_MOCK.fasta
	sh get_references.bash
