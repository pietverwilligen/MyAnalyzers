#!/bin/bash
set -x
set +x 
export FILE1="$1.txt"
export FILE2="$1.py"

touch ${FILE2}
echo "created ${FILE2}"

echo "import FWCore.ParameterSet.Config as cms" >> ${FILE2} 
echo "maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )" >> ${FILE2}
echo "readFiles = cms.untracked.vstring()" >> ${FILE2}
echo "secFiles = cms.untracked.vstring()" >> ${FILE2}
echo "source = cms.Source (\"PoolSource\",fileNames = readFiles, secondaryFileNames = secFiles)" >> ${FILE2}
echo "readFiles.extend( [" >> ${FILE2}

for i in `cat ${FILE1}`; 
do
    echo "'${i}'," >> ${FILE2};
done

echo "] );" >> ${FILE2}
echo "secFiles.extend( [" >> ${FILE2}
echo "] );" >> ${FILE2}
