#!/usr/bin/env bash

databases=${databases:-/media/bigdrive1/Databases/dbCAN}

echo "Downloading and formatting databases in: ${databases}"

mkdir $databases && cd $databases
wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/readme.txt
wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/CAZyDB.07312020.fa
diamond makedb --in CAZyDB.07312020.fa -d CAZy
wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/dbCAN-HMMdb-V9.txt
mv dbCAN-HMMdb-V9.txt dbCAN.txt
hmmpress dbCAN.txt
wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/tcdb.fa
diamond makedb --in tcdb.fa -d tcdb
wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/tf-1.hmm
hmmpress tf-1.hmm
wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/tf-2.hmm
hmmpress tf-2.hmm
wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/stp.hmm
hmmpress stp.hmm
