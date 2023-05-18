#!/bin/bash

# second step
# convert txt file to python file that can be used

# export infile="MinBias_Run2_13TeV_100k_Neutron_XS_1E4s.txt"
# export outfile="MinBias_Run2_13TeV_100k_Neutron_XS_1E4s.py"

# 13X 2026
# export infile="MinBias_Phase2_14TeV_TuneCP5_100k_Neutron_XS_2026D99_1E4s.txt"
# export outfile="MinBias_Phase2_14TeV_TuneCP5_100k_Neutron_XS_2026D99_1E4s.py"

# 13X 2023
export infile="MinBias_Phase1_14TeV_TuneCP5_100k_Neutron_XS_2023DB_1E4s.txt"
export outfile="MinBias_Phase1_14TeV_TuneCP5_100k_Neutron_XS_2023DB_1E4s.py"

echo "import FWCore.ParameterSet.Config as cms" >> $outfile
echo "maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )" >> $outfile
echo "readFiles = cms.untracked.vstring()" >> $outfile
echo "secFiles = cms.untracked.vstring()" >> $outfile
echo "source = cms.Source (\"PoolSource\",fileNames = readFiles, secondaryFileNames = secFiles)" >> $outfile
echo "" >> $outfile
echo "readFiles.extend( [" >> $outfile

export infilelist=`cat $infile`
for i in $infilelist;
do echo "'${i}'," >> $outfile;
done

echo "]);" >> $outfile
echo "" >> $outfile
echo "secFiles.extend( [] )" >> $outfile
