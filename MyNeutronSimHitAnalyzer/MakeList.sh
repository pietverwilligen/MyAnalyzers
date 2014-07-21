#!/bin/bash

for i in `ls -lrth /lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/ | awk '{print $9}'`; do echo "/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/${i}" >> MinBias_13TeV_2023Muon_1E5Events_filtered.txt; done
