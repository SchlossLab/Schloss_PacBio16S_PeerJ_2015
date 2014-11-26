# Let's use bash to clean up the folders a bit since these file names are long
# and unwieldy. We'll start by creating a folder called  `ccs.fastq`. We'll
# transfer the fastqs from the `raw_data` folder into into region-specific
# folders:

mkdir ccs.fastqs
cd ccs.fastqs
mkdir v19 v16 v15 v13 v35 v4
for REGION in 19 16 15 13 35 4
do
    cp ../raw_data/*_p$REGION/Analysis_Results/*ccs.fastq v$REGION
done

ls */*fastq | cut -f 1 -d "/" | sort | uniq -c

# This last command allows us to see that every region has 3 \*.ccs.fastq files,
# except for the v15 (6 files) and v19 (4 files) regions. Let's go ahead and
# concatenate those fastq files and dump them into individual folders.

cd ../
mkdir v19 v16 v15 v13 v35 v4

for REGION in $(ls -d v*)
do
    > $REGION/$REGION.fastq
    for FILE in ccs.fastqs/$REGION/*.ccs.fastq
    do
        cat $FILE >> $REGION/$REGION.fastq
    done
done

# We'd also like to concatenate the subread fasta files so that we can see how
# many base pairs were sequenced for each fragment:

mkdir subreads.fasta
cd subreads.fasta

mkdir v19 v16 v15 v13 v35 v4
for REGION in 19 16 15 13 35 4
do
    cp ../raw_data/*_p$REGION/Analysis_Results/*subreads.fasta v$REGION > v$REGION/v$REGION.subreads.fasta

    for FILE in v$REGION/*.subreads.fasta
    do
        cat $FILE >> v$REGION/v$REGION.subreads.fasta
    done
done


# Now we'd like to generate fasta and quality score files for each ccs file.
# When we do this, it is critical that we set the pacbio parameter to T. This is
# because by default PacBio will assign a base call even if the quality score is
# zero:

cd ../

for REGION in v*
do
    mothur "#fastq.info(fastq=$REGION/$REGION.fastq, pacbio=T, outputdir=./$REGION/)"
done


# Using the primer and barcode information from above, we can make
# region-specific oligos files and put them into the regional folders:

for REGION in v*
do
    grep ""$REGION"" pacbio.oligos > $REGION/$REGION.oligos
    grep "barcode" pacbio.oligos >> $REGION/$REGION.oligos
done


# Now we're all set to run some mothur commands. Since each file is a mixture of
# our mock community and data from soil, human feces, and mouse feces, we need
# to split the fasta and qual files by barcode. Let's initially be generous and
# allow for 2 mismatches to each barcode and 4 mismatches to each primer. To
# keep things simple, we'll concatenate the three mock community fasta, quality
# score, and groups files. We modified the source code to output the total
# number of mismatches to the barcodes and primers on the header line for each
# sequence in the trim file.

for REGION in v*
do
    cd $REGION
    mothur "#trim.seqs(fasta=$REGION.fasta, qfile=$REGION.qual, oligos=$REGION.oligos, checkorient=T, pdiffs=6, bdiffs=4, allfiles=T, processors=8)"
    cat $REGION*mock?.$REGION.fasta > $REGION.mock.fasta
    cat $REGION*mock?.$REGION.qual > $REGION.mock.qual
    cat $REGION*mock?.$REGION.groups > $REGION.mock.groups
    grep ">" $REGION.mock.fasta > $REGION.mismatches
    cd ../
done


# We are now ready to calculate the error rate for the various regions using the
# mock community and the HMP_MOCK sequence data. We will use mothur to align the
# sequences to HMP_MOCK.align, determine the start and end positions of the
# alignment, and calculate the error rate.

for REGION in v*
do
    mothur "#align.seqs(fasta=$REGION/$REGION.mock.fasta, reference=HMP_MOCK.align, processors=8, outputdir=./$REGION/);
        filter.seqs(fasta=$REGION/$REGION.mock.align-HMP_MOCK.align, vertical=T);
        summary.seqs();
        seq.error(fasta=$REGION/$REGION.mock.filter.fasta, reference=$REGION/HMP_MOCK.filter.fasta, report=$REGION/$REGION.mock.align.report, qfile=$REGION/$REGION.mock.qual, processors=8);"
done
