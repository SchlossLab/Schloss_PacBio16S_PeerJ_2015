# This command will pull down a copy of the HMP_MOCK.fasta reference sequence
# file and the necessary reference databases. Yes, I realize that this is a hack
# to avoid figuring out how to use make correctly.

mkdir -p references
cd references


# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v119 and described at http://blog.mothur.org/2014/08/08/SILVA-v119-reference-files/.
# We will use the full-length version of the database, which contains 137,879
# bacterial sequences. This also contains the reference taxonomy. We will limit
# the databases to only include bacterial sequences.


if [ ! -e got.silva ] || [ got.silva -ot silva.bacteria.align ] || [ got.silva -ot silva.bacteria.fasta ] || [ got.silva -ot silva.bacteria.tax ]; then
  wget -N http://www.mothur.org/w/images/2/27/Silva.nr_v119.tgz
  tar xvzf Silva.nr_v119.tgz
  mothur "#get.lineage(fasta=silva.nr_v119.align, taxonomy=silva.nr_v119.tax, taxon=Bacteria)"
  mothur "#degap.seqs(fasta=silva.nr_v119.pick.align)"
  mv silva.nr_v119.pick.align silva.bacteria.align
  mv silva.nr_v119.pick.tax silva.bacteria.tax
  mv silva.nr_v119.pick.ng.fasta silva.bacteria.fasta
  touch got.silva
  #silva.bacteria.align
  #silva.bacteria.fasta
  #silva.bacteria.tax
fi

# We also want the greengenes reference taxonomy. This version is from the
# greengenes v13_8_99 and is described at
# http://blog.mothur.org/2014/08/12/greengenes-v13_8_99-reference-files/

if [ ! -e got.gg ] || [ got.gg -ot gg_13_8_99.fasta ] || [ got.gg -ot gg_13_8_99.gg.tax ]; then
  wget -N http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
  tar xvzf Gg_13_8_99.taxonomy.tgz
  touch got.gg
  #gg_13_8_99.fasta
  #gg_13_8_99.gg.tax
fi

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are
# described at http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

if [ ! -e got.rdp ] || [ got.rdp -ot trainset10_082014.pds.fasta ] || [ got.rdp -ot trainset10_082014.pds.tax ]; then
  wget -N http://www.mothur.org/w/images/2/24/Trainset10_082014.pds.tgz
  tar xvzf Trainset10_082014.pds.tgz
  mv trainset10_082014.pds/* ./
  rmdir trainset10_082014.pds
  touch got.rdp
  #trainset10_082014.pds.fasta
  #trainset10_082014.pds.tax
fi

# Finally, we want to align the mock community reference sequences to our newly
# created silva.bacteria.fasta file...

if [ ! -e HMP_MOCK.align ] || [ HMP_MOCK.align -ot HMP_MOCK.fasta ]; then
  mothur "#align.seqs(fasta=HMP_MOCK.fasta, reference=silva.bacteria.align)"
  #HMP_MOCK.align
fi

# Clean up and move out...

rm README*

cd ../
