#!/bin/bash
TXTDIR=/hdfs/local/stpol/joosep/susy/run4/
ODIR=/home/joosep/susy/run4/
rm -Rf $ODIR
mkdir $ODIR
rsync --rsh=ssh --partial -u jpata@lxplus.cern.ch:~/work/public/*.bz2 $TXTDIR
python make_table.py "$TXTDIR/*.bz2" $ODIR 16
python merge_tables.py "$ODIR/*.h5"
mv out.h5 nmssm4.h5
