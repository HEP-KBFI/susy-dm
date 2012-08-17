a=$(echo $RANDOM | md5sum)
TAG="${a:0:8}"
WORKDIR="/tmp/jpata"
NMSSMDIR=$WORKDIR/NMSSMTools_3.2.1_$TAG

cp -R /afs/cern.ch/user/j/jpata/private/NMSSMTools_3.2.1 $NMSSMDIR
cd $NMSSMDIR
make
cat randinp.dat.template | sed "s/SEEDVAL/$RANDOM/" > randinp_$TAG.dat 
./run randinp_$TAG.dat
cp rand*$TAG.dat /afs/cern.ch/user/j/jpata/work/public/
rm -Rf $NMSSMDIR
