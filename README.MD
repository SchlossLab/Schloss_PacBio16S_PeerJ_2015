## The paper

This is the repository for the manuscript "Sequencing 16S rRNA gene fragments
using the PacBio SMRT DNA sequencing system" written by Patrick D. Schloss,
Sarah L. Westcott, Matthew L. Jenior, and Sarah K. Highlander.The raw data can
be obtained from the Sequence Read Archive at NCBI under accession SRP051686,
which are associated with BioProject PRJNA271568 and the
processed movies can be obtained using the commands in mothur_raw.bash. The
entire manuscript can be generated on a Mac OS X or Linux computer by running
the following commands:

```
git clone https://github.com/SchlossLab/PacBio_16S.git
cd PacBio_16S
sh write.paper
```

In running these commands we assume that you have a copy of R and mothur
installed in your PATH. You can then convert the resulting markdown file to a
Word docx file by doing:

```
sh compile_docx.sh PacBioManuscript.md
open PacBioManuscript.docx
```



## Background

We ran 8 SMRT cells and their data were located in the following folders that
were housed within a folder called `raw_data`:

```
B01_1_Cell2_PacBioRun109_PSchloss311_p15  
C01_1_Cell3_PacBioRun109_PSchloss315_p16  
D01_1_Cell4_PacBioRun108_PSchloss310_p13  
D01_1_Cell4_PacBioRun109_PSchloss316_p4  
E01_1_Cell5_PacBioRun108_PSchloss311_p15  
E01_1_Cell5_PacBioRun109_PSchloss317_p19  
F01_1_Cell6_PacBioRun108_PSchloss312_p35  
H01_1_Cell8_PacBioRun112_PSchloss319_p19
```

If you substitute the "\_p" for "\_v" you'll see what region each run
represented. There were two V15 runs and two V19 runs. These regions were
amplified with the following primer sets:

```
V19	AGRGTTTGATYMTGGCTCAG	GGYTACCTTGTTACGACTT  
V16	AGRGTTTGATYMTGGCTCAG	ACRACACGAGCTGACGAC  
V15	AGRGTTTGATYMTGGCTCAG	CCCGTCAATTCMTTTRAGT  
V13	AGRGTTTGATYMTGGCTCAG	ATTACCGCGGCTGCTGG  
V35	CCTACGGGAGGCAGCAG	CCCGTCAATTCMTTTRAGT  
V4	 GTGCCAGCMGCCGCGGTAA	GGACTACHVGGGTWTCTAAT
```

We also put sequencing indices at the 5' and 3' of the inserts to correspond to
the different samples:

```
mock1	ccaac	ccaac  
mock2	ggttg	ccaac  
mock3	ttggt	ccaac  
mouse1	aggtg	ccaac  
mouse2	cttac	ccaac  
mouse3	gaact	ccaac  
human1	tccga	ccaac  
human2	acagt	ccaac  
human3	cactg	ccaac  
soil1	aacca	ccaac  
soil2	tgtca	ccaac  
soil3	aaacc	ccaac
```


Putting it all together in an oligos file (pacbio.oligos), we get...

```
primer	AGRGTTTGATYMTGGCTCAG	GGYTACCTTGTTACGACTT	v19  
primer	AGRGTTTGATYMTGGCTCAG	ACRACACGAGCTGACGAC	v16  
primer	AGRGTTTGATYMTGGCTCAG	CCCGTCAATTCMTTTRAGT	v15  
primer	AGRGTTTGATYMTGGCTCAG	ATTACCGCGGCTGCTGG	v13  
primer	CCTACGGGAGGCAGCAG	CCCGTCAATTCMTTTRAGT	v35  
primer	GTGCCAGCMGCCGCGGTAA	GGACTACHVGGGTWTCTAAT	v4  
barcode	ccaac	ccaac	mock1  
barcode	ggttg	ccaac	mock2  
barcode	ttggt	ccaac	mock3  
barcode	aggtg	ccaac	mouse1  
barcode	cttac	ccaac	mouse2  
barcode	gaact	ccaac	mouse3  
barcode	tccga	ccaac	human1  
barcode	acagt	ccaac	human2  
barcode	cactg	ccaac	human3  
barcode	aacca	ccaac	soil1  
barcode	tgtca	ccaac	soil2  
barcode	aaacc	ccaac	soil3
```
