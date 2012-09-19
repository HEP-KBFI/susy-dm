a=$(echo $RANDOM | md5sum)
#TAG="${a:0:8}"
TAG=$LSB_JOBID
WORKDIR="/afs/cern.ch/user/j/jpata/work/public/gridjob/"
NMSSMDIR=$WORKDIR/NMSSMTools_3.2.1_$TAG
BASEDIR=/afs/cern.ch/user/j/jpata/private/susy-dm/NMSSMTools_3.2.1
OUTDIR=/afs/cern.ch/user/j/jpata/work/public/

cp -R $BASEDIR $NMSSMDIR
cd $NMSSMDIR
make clean
make init
make
cat inp.template | sed "s/SEEDVAL/$RANDOM/" > randinp_$TAG.dat 
cd main
./nmspec_rand < ../randinp_$TAG.dat
cat out.txt | python binaryOutput.py
rm out.txt
bzip2 -c9 test.npy > randout_$TAG.npy.bz2
rm test.npy
rsync randout_$TAG.npy.bz2 $OUTDIR/
cd ..
rm -Rf $NMSSMDIR
