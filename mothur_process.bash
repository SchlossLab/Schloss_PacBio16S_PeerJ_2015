# To recap, the pipeline consists of removing any fragment that has more than
# one mismatch to the barcodes or primers, contains an ambiguous base call, has
# an average quality score below 60, has a homopolymer longer than 8 nt, and
# does not align correctly to the region of interest. We would now like to run
# that protocol from the beginning on all of our samples. Let's create some new
# folders and make links to the original fasta files:

mkdir analysis
cd analysis
mkdir v13  v15  v16  v19  v35  v4
cd ../pipeline_dev

for REGION in v*
do
    cp -l $REGION/$REGION.fasta $REGION/$REGION.qual $REGION/$REGION.oligos ../analysis/$REGION
done
cd ../analysis


# Let's run this using trim.seqs in mothur on each of the files applying the
# criteria we developed earlier. While we're at it, we'll go ahead and unique,
# align the sequences, and summarize the alignments.

for REGION in v*
do
    cd $REGION
    mothur "#trim.seqs(fasta=$REGION.fasta, qfile=$REGION.qual, oligos=$REGION.oligos, checkorient=T, pdiffs=1, bdiffs=1, tdiffs=1, maxambig=0, maxhomop=8, qaverage=60, processors=8);
            unique.seqs(fasta=current);
            align.seqs(fasta=current, reference=../../references/silva.bacteria.align);
            summary.seqs(name=current)"
    cd ../
done

# Looking at the output from the summary.seqs commands we come up with the
# following positions within the alignment for our start and end values:
#
# | Region | Start | End   |
# |-------------------------
# | v13    | 1046  | 13125 |
# | v15    | 1046  | 27659 |
# | v16    | 1046  | 34113 |
# | v19    | 1046  | 43116 |
# | v35    | 6428  | 27659 |
# | v4     | 13862 | 23444 |

# Now we want to run screen.seqs to remove sequences that do not start at or
# before the start position or end at or after end.

mothur "#screen.seqs(fasta=v13.trim.unique.align, name=v13.trim.names, group=v13.groups, start=1046, end=13125, processors=8, inputdir=./v13)"
mothur "#screen.seqs(fasta=v15.trim.unique.align, name=v15.trim.names, group=v15.groups, start=1046, end=27659, processors=8, inputdir=./v15)"
mothur "#screen.seqs(fasta=v16.trim.unique.align, name=v16.trim.names, group=v16.groups, start=1046, end=34113, processors=8, inputdir=./v16)"
mothur "#screen.seqs(fasta=v19.trim.unique.align, name=v19.trim.names, group=v19.groups, start=1046, end=43116, processors=8, inputdir=./v19)"
mothur "#screen.seqs(fasta=v35.trim.unique.align, name=v35.trim.names, group=v35.groups, start=6428, end=27659, processors=8, inputdir=./v35)"
mothur "#screen.seqs(fasta=v4.trim.unique.align, name=v4.trim.names, group=v4.groups, start=13862, end=23444, processors=8, inputdir=./v4)"


# After we run screen.seqs, we'll go ahead and complete the pipeline by using
# filter.seqs to remove any columns that only contain gaps (vertical=T) or that
# contain missing data (trump=.) and we'll run unique.seqs and summary.seqs to
# get the sequence lengths:

for REGION in v*
do
    mothur "#set.dir(input=./$REGION, output=./$REGION);
            filter.seqs(fasta=$REGION.trim.unique.good.align-../references/HMP_MOCK.align, vertical=T, trump=., processors=8);
            unique.seqs(fasta=$REGION.trim.unique.good.filter.fasta, name=$REGION.trim.good.names);
            summary.seqs(name=current)"
done


# Next we'll want to run pre.cluster using the diffs parameter equal to 1
# difference per 100 bp. In the previous step we got those lengths:
#
# | Region | Length | diffs |
# |--------|--------|-------|
# | v13    | 489    |  4    |
# | v15    | 879    |  8    |
# | v16    | 1028   | 10    |
# | v19    | 1458   | 14    |
# | v35    | 545    |  5    |
# | v4     | 253    |  2    |

# Once we run pre.cluster, we'll check for chimeras, remove them, cluster the
# sequences and then make a shared file from them

mothur "#pre.cluster(fasta=v13/v13.trim.unique.good.filter.unique.fasta, name=v13/v13.trim.unique.good.filter.names, group=v13/v13.good.groups, diffs=4, processors=8);"
mothur "#pre.cluster(fasta=v15/v15.trim.unique.good.filter.unique.fasta, name=v15/v15.trim.unique.good.filter.names, group=v15/v15.good.groups, diffs=10, processors=8);"
mothur "#pre.cluster(fasta=v16/v16.trim.unique.good.filter.unique.fasta, name=v16/v16.trim.unique.good.filter.names, group=v16/v16.good.groups, diffs=10, processors=8);"
mothur "#pre.cluster(fasta=v19/v19.trim.unique.good.filter.unique.fasta, name=v19/v19.trim.unique.good.filter.names, group=v19/v19.good.groups, diffs=14, processors=8);"
mothur "#pre.cluster(fasta=v35/v35.trim.unique.good.filter.unique.fasta, name=v35/v35.trim.unique.good.filter.names, group=v35/v35.good.groups, diffs=5, processors=8);"
mothur "#pre.cluster(fasta=v4/v4.trim.unique.good.filter.unique.fasta, name=v4/v4.trim.unique.good.filter.names, group=v4/v4.good.groups, diffs=2, processors=8);"

for REGION in v*
do
    mothur "#set.dir(input=./$REGION, output=./$REGION);
      chimera.uchime(fasta=$REGION.trim.unique.good.filter.unique.precluster.fasta, name=$REGION.trim.unique.good.filter.unique.precluster.names, group=$REGION.good.groups, dereplicate=T, processors=8);
      remove.seqs(fasta=current, group=current, name=current, accnos=current, dups=F);
      dist.seqs(fasta=current, cutoff=0.15, processors=8);
      cluster();
      make.shared(list=current, label=0.03, group=current)"
done



# Let's merge the three replicates for each library and recreate a shared file
# based on the merged group files and get the number of sequences per library

for REGION in v*
do
  cut -f 1 -d "." $REGION/$REGION.good.pick.groups | sed "s/.$//" > $REGION/$REGION.good.pick.merge.groups
  cp $REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.an.list $REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.an.merge.list
  mothur "#make.shared(list=$REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.an.merge.list, group=$REGION/$REGION.good.pick.merge.groups, label=0.03);
        summary.single(shared=$REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.an.merge.shared, calc=nseqs)"
done


# Let's get the error rates for our mock communities from before and after
# running the pre.cluster steps.

for REGION in v*
do
  for REP in 1 2 3
  do
    mothur "#set.dir(input=./$REGION, output=./$REGION);
        get.groups(fasta=$REGION.trim.unique.good.filter.unique.fasta, group=$REGION.good.pick.groups, name=$REGION.trim.unique.good.filter.names, groups=mock$REP.$REGION);
        system(mv $REGION/$REGION.trim.unique.good.filter.pick.names $REGION/$REGION.mock$REP.unique.names);
        system(mv $REGION/$REGION.trim.unique.good.filter.unique.pick.fasta $REGION/$REGION.mock$REP.unique.fasta);
        seq.error(fasta=$REGION.mock$REP.unique.fasta, name=$REGION.mock$REP.unique.names, reference=HMP_MOCK.filter.fasta, aligned=T, processors=8)"

    mothur "#set.dir(input=./$REGION, output=./$REGION);
        get.groups(fasta=$REGION.trim.unique.good.filter.unique.precluster.pick.fasta, group=$REGION.good.pick.groups, name=$REGION.trim.unique.good.filter.unique.precluster.pick.names, groups=mock$REP.$REGION);
        system(mv $REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.pick.names $REGION/$REGION.mock$REP.precluster.names);
        system(mv $REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.pick.fasta $REGION/$REGION.mock$REP.precluster.fasta);
        seq.error(fasta=$REGION.mock$REP.precluster.fasta, name=$REGION.mock$REP.precluster.names, reference=HMP_MOCK.filter.fasta, aligned=T, processors=8)"
  done
done


# Let's remove all of the chimeras we didn't detect from the mock communities

for REGION in v*
do
  grep "2$" $REGION/$REGION.mock?.precluster*error.summary | cut -f 2 -d ":" | cut -f 1 > $REGION/$REGION.extra.chimeras
  cat $REGION/$REGION.mock?.precluster.fasta > $REGION/$REGION.mock.precluster.perfect.fasta
  cat $REGION/$REGION.mock?.precluster.names > $REGION/$REGION.mock.precluster.perfect.names
  mothur "#remove.seqs(fasta=$REGION/$REGION.mock.precluster.perfect.fasta, name=$REGION/$REGION.mock.precluster.perfect.names, accnos=$REGION/$REGION.extra.chimeras);
          dist.seqs(cutoff=0.15);
          cluster();
          summary.single(calc=nseqs, label=0.03)"
done



# Let's count how many times each sequence shows up in each of the
# merged groups

for REGION in v*
do
  mothur "#count.seqs(group=$REGION/$REGION.good.pick.merge.groups, name=$REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.names)"
done


# Let's find the smallest number of sequences in each region
MIN=$(cat v*/*mock.precluster.perfect.pick.an.summary | cut -f 2 | grep "\\." | sort -n | head -n 1 | cut -f 1 -d ".")
SMALL_LIB=$(cat */*.trim.unique.good.filter.unique.precluster.pick.an.merge.groups.summary | sort -k 3 -n | head -n 7 | tail -n 1 | cut -f 3 | cut -f 1 -d ".")

MIN=$(echo $(($MIN<$SMALL_LIB?$MIN:$SMALL_LIB)))



# Let's rarefy everything to $MIN reads per sample since this was the size of
# the smallest v19 library

for REGION in v*
do
  # the observed data...
  mothur "#summary.single(shared=$REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.an.merge.shared, subsample=$MIN, calc=sobs-coverage-nseqs, label=0.03, iters=1000)"

  # the perfect chimera removal data
  mothur "#summary.single(list=$REGION/$REGION.mock.precluster.perfect.pick.an.list, subsample=$MIN, calc=sobs-coverage-nseqs, label=0.03, iters=1000)"
done





# Let's use the reference sequences we sampled to see what the "correct" number
# of OTUs would have been

for REGION in v*
do
  grep "1$" $REGION/$REGION.mock*error.summary | cut -f 2 | sort | uniq > $REGION/mock.ref.accnos
  mothur "#get.seqs(fasta=$REGION/HMP_MOCK.filter.fasta, accnos=$REGION/mock.ref.accnos);
          dist.seqs(cutoff=0.15, output=lt); cluster(phylip=current); summary.single(label=0.03, calc=sobs)"
done


# Now we'd like to go ahead and classify all of our sequences from the different
# libraries and the mock community reference using the RDP, greengenes, and
# SILVA training sets.

for REGION in v*
do
  cat $REGION/$REGION.mock?.precluster.fasta > $REGION/$REGION.mock.precluster.fasta
  mothur "#classify.seqs(fasta=$REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.fasta-$REGION/$REGION.mock.precluster.fasta-$REGION/HMP_MOCK.filter.fasta, reference=../references/trainset10_082014.pds.fasta, taxonomy=../references/trainset10_082014.pds.tax, processors=8);
  classify.seqs(fasta=$REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.fasta-$REGION/$REGION.mock.precluster.fasta-$REGION/HMP_MOCK.filter.fasta, reference=../references/gg_13_8_99.fasta, taxonomy=../references/gg_13_8_99.gg.tax, processors=8);
  classify.seqs(fasta=$REGION/$REGION.trim.unique.good.filter.unique.precluster.pick.fasta-$REGION/$REGION.mock.precluster.fasta-$REGION/HMP_MOCK.filter.fasta, reference=../references/silva.bacteria.fasta, taxonomy=../references/silva.bacteria.tax, processors=8);"
done






# It is commonly said that the error profile in PacBio-generated sequence data
# are random. As we saw above, the substitution preference was random. However,
# if the errors are truly random, then we would expect to have a bunch of
# singleton error sequences, not sequences with errors that show up multiple
# times. Let's see what type of error profile we have. Because they have the
# most reads, let's consider the v4, v35, and v15 mock community samples from
# before the pre.cluster step:

for REGION in v*
do
  mothur "#get.groups(fasta=$REGION/$REGION.trim.unique.good.filter.unique.fasta, name=$REGION/$REGION.trim.unique.good.filter.names, group=$REGION/$REGION.good.groups, groups=mock1.$REGION-mock2.$REGION-mock3.$REGION);
            seq.error(fasta=current, name=current, reference=$REGION/HMP_MOCK.filter.fasta, processors=8)"
done


#Back to the parent directory!
cd ../
