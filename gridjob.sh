a=$(echo $RANDOM | md5sum)
TAG="${a:0:8}"
WORKDIR="/tmp/jpata"
NMSSMDIR=$WORKDIR/NMSSMTools_3.2.1_$TAG
BASEDIR=/afs/cern.ch/user/j/jpata/private/NMSSMTools_3.2.1
OUTDIR=/afs/cern.ch/user/j/jpata/work/public/

cp -R $BASEDIR $NMSSMDIR
cd $NMSSMDIR
make
cat randinp.dat.template | sed "s/SEEDVAL/$RANDOM/" > randinp_$TAG.dat 
./run randinp_$TAG.dat
bzip2 -c -9 randout_$TAG.dat > randout_$TAG.dat.bz2
cp randout_$TAG.dat.bz2 $OUTDIR/
rm -Rf $NMSSMDIR
