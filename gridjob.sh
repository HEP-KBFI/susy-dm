a=$(echo $RANDOM | md5sum)
#TAG="${a:0:8}"
TAG=$LSB_JOBID
WORKDIR="/afs/cern.ch/user/j/jpata/work/public/gridjob/"
NMSSMDIR=$WORKDIR/NMSSMTools_3.2.1_$TAG
BASEDIR=/afs/cern.ch/user/j/jpata/private/susy-dm/NMSSMTools_3.2.1
OUTDIR=/afs/cern.ch/user/j/jpata/work/public/

cp -R $BASEDIR $NMSSMDIR
cd $NMSSMDIR
make
cat randinp_NMSSM2.dat.template | sed "s/SEEDVAL/$RANDOM/" > randinp_$TAG.dat 
cd main
./nmhdecay_rand < ../randinp_$TAG.dat | grep "[\d\s\+-E\.]+" | bzcat -z9 > ../randout_$TAG.dat.bz2
cd ..
rsync randout_$TAG.dat.bz2 $OUTDIR/
rm -Rf $NMSSMDIR
