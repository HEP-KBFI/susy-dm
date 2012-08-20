a=$(echo $RANDOM | md5sum)
TAG="${a:0:8}"
WORKDIR="/tmp/jpata"
NMSSMDIR=$WORKDIR/NMSSMTools_3.2.1_$TAG
BASEDIR=/afs/cern.ch/user/j/jpata/private/susy-dm/NMSSMTools_3.2.1
OUTDIR=/afs/cern.ch/user/j/jpata/work/public/

cp -R $BASEDIR $NMSSMDIR
cd $NMSSMDIR
make
cat randinp_NMSSM2.dat.template | sed "s/SEEDVAL/$RANDOM/" > randinp_$TAG.dat 
./run randinp_$TAG.dat
ls -al
grep "[\s\d\.E\+-]+\n*" randout_$TAG.dat | bzip -c9 > randout_$TAG.dat.bz2
du -sh *.bz2
rsync randout_$TAG.dat.bz2 $OUTDIR/
rm -Rf $NMSSMDIR
