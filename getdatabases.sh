#!/bin/bash
# The script downloads, extracts converts databases used in the paper
# such that the Matlab scripts can find the files.

SQAMARCHIVE=SQAM_FLAC.zip
MOCHAARCHIVE1=fsew0_v1.1.tar.gz
MOCHAARCHIVE2=maps0.tar.gz


SQAMLINK=https://tech.ebu.ch/files/live/sites/tech/files/shared/testmaterial/$SQAMARCHIVE
MOCHALINK1=http://data.cstr.ed.ac.uk/mocha/$MOCHAARCHIVE1
MOCHALINK2=http://data.cstr.ed.ac.uk/mocha/$MOCHAARCHIVE2

CURRDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

wget -N $SQAMLINK
wget -N $MOCHALINK1
wget -N $MOCHALINK2


mkdir -p $CURRDIR/Databases/SQAM
unzip -o $SQAMARCHIVE -d $CURRDIR/Databases/SQAM
cd $CURRDIR/Databases/SQAM
ls *.flac | xargs flac -df
find . -type f -not -name '*.wav' -delete

cd $CURRDIR
mkdir -p $CURRDIR/Databases/maps0
tar -C $CURRDIR/Databases/maps0 -zxvf $MOCHAARCHIVE2
find $CURRDIR/Databases/maps0 -type f -not -name '*.wav' -delete

mkdir -p $CURRDIR/Databases/fsew0
tar -C $CURRDIR/Databases/fsew0 -zxvf $MOCHAARCHIVE1
find $CURRDIR/Databases/fsew0 -type f -not -name '*.wav' -delete

