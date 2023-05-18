#!/bin/bash

# for i in `ls -lrth /lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/ | awk '{print $9}'`; do echo "/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/${i}" >> MinBias_13TeV_2023Muon_1E5Events_filtered.txt; done

### alternatively using the dasgoclient
### 100k events #######################
# dasgoclient --query 'file dataset=/MinBias_14TeV_GEN_SIM_XS_2026D49_100k_11X_ProdCutsDefault_v1/piet-crab_MinBias_14TeV_100k_11X_ProdCutsDefault_v1-7df21a947851a340cb0cd905c371c7a6/USER instance=prod/phys03' &> MinBias_14TeV_GEN_SIM_XS_2026D49_100k_11X_ProdCutsDefault_v1_allfiles.txt

### 1M events #########################
# dasgoclient --query 'file dataset=/MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1/piet-crab_MinBias_14TeV_1M_11X_ProdCutsDefault_v1_pt1-7df21a947851a340cb0cd905c371c7a6/USER instance=prod/phys03' &> MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p1_allfiles.txt
# source ConvertToPython.sh MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p1_allfiles

# dasgoclient --query 'file dataset=/MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1/piet-crab_MinBias_14TeV_1M_11X_ProdCutsDefault_v1_pt2-7df21a947851a340cb0cd905c371c7a6/USER instance=prod/phys03' &> MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p2_allfiles.txt
# source ConvertToPython.sh MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p2_allfiles

# dasgoclient --query 'file dataset=/MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1/piet-crab_MinBias_14TeV_1M_11X_ProdCutsDefault_v1_pt3-7df21a947851a340cb0cd905c371c7a6/USER instance=prod/phys03' &> MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p3_allfiles.txt
# source ConvertToPython.sh MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p3_allfiles

# dasgoclient --query 'file dataset=/MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1/piet-crab_MinBias_14TeV_1M_11X_ProdCutsDefault_v1_pt4-7df21a947851a340cb0cd905c371c7a6/USER instance=prod/phys03' &> MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p4_allfiles.txt
# source ConvertToPython.sh MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p4_allfiles

# dasgoclient --query 'file dataset=/MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1/piet-crab_MinBias_14TeV_1M_11X_ProdCutsDefault_v1_pt5-7df21a947851a340cb0cd905c371c7a6/USER instance=prod/phys03' &> MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p5_allfiles.txt
# source ConvertToPython.sh MinBias_14TeV_GEN_SIM_XS_2026D49_1M_11X_ProdCutsDefault_v1_p5_allfiles


### No Neutron Background === 100k events ###
# dasgoclient --query 'file dataset=/MinBias_14TeV_GEN_SIM_NoBkg_2026D49_100k_11X_ProdCutsDefault/piet-crab_MinBias_14TeV_100k_11X_ProdCutsDefault_v2-401867b49e0939e637da3bf19603d355/USER instance=prod/phys03' &> MinBias_14TeV_GEN_SIM_NoBkg_2026D49_100k_11X_ProdCutsDefault_v2_allfiles.txt
# source ConvertToPython.sh MinBias_14TeV_GEN_SIM_NoBkg_2026D49_100k_11X_ProdCutsDefault_v2_allfiles


### Neutron Background in Run-3 Geometry (2021) ###
# dasgoclient --query 'file dataset=/MinBias_13TeV_GEN_SIM_XS_2021_100k_11X_ProdCutsDefault/piet-crab_MinBias_13TeV_100k_11X_2021_ProdCutsDefault-57e6adbb11f628ac8185af7dc72e0ac6/USER instance =prod/phys03' &> MinBias_13TeV_GEN_SIM_2021_100k_11X_ProdCutsDefault_v1_allfiles.txt
# source ConvertToPython.sh MinBias_13TeV_GEN_SIM_2021_100k_11X_ProdCutsDefault_v1_allfiles

### ZMM in Run-3 Geometry (2021) ###
dasgoclient --query 'file dataset=/ZMM_13TeV_GEN_SIM_XS_2021_100k_11X_ProdCutsDefault/piet-crab_ZMM_14TeV_100k_11X_2021_ProdCutsDefault-bf470c918fd3e703c3904d56fca9a467/USER instance =prod/phys03' &> ZMM_14TeV_GEN_SIM_2021_100k_11X_ProdCutsDefault_v1_allfiles.txt
source ConvertToPython.sh ZMM_14TeV_GEN_SIM_2021_100k_11X_ProdCutsDefault_v1_allfiles
