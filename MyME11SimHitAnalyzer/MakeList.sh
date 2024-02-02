#!/bin/bash

### 2016 Filtered ###
#####################
# for i in `ls -lrth /lustre/cms/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/ | awk '{print $9}'`; do echo "/store/user/piet/NeutronBackground/MinBias_13TeV_2023Muon_1E5Events_filtered/${i}" >> MinBias_13TeV_2023Muon_1E5Events_filtered.txt; done

### 2019 ############
#####################
# for i in `ls -lrth /lustre/cms/store/user/piet/HGCAL-geom-study/MinBias_13TeV_GEN_SIM_HGCAL_TDRgeom_100k_2023D17_Neutron_XS_100ms_v1/crab_MinBias_13TeV_GEN-SIM_HGCAL-TDRgeom_100k_1040_2023D17_Neutron_XS_100ms_v2/190206_130559/0000/ | awk '{print $9}'`; do echo "/store/user/piet/HGCAL-geom-study/MinBias_13TeV_GEN_SIM_HGCAL_TDRgeom_100k_2023D17_Neutron_XS_100ms_v1/crab_MinBias_13TeV_GEN-SIM_HGCAL-TDRgeom_100k_1040_2023D17_Neutron_XS_100ms_v2/190206_130559/0000/${i}" >> MinBias_13TeV_HGCAL_TDRgeom_100k_2023D17_Neutron_XS_100ms.txt; done

# for i in `ls -lrth /lustre/cms/store/user/piet/HGCAL-geom-study/MinBias_13TeV_GEN_SIM_HGCAL_TDRgeom_100k_2023D28_Neutron_XS_100ms_v1/crab_MinBias_13TeV_GEN-SIM_HGCAL-TDRgeom_100k_1040_2023D28_Neutron_XS_100ms_v2/190206_130550/0000/ | awk '{print $9}'`; do echo "/store/user/piet/HGCAL-geom-study/MinBias_13TeV_GEN_SIM_HGCAL_TDRgeom_100k_2023D28_Neutron_XS_100ms_v1/crab_MinBias_13TeV_GEN-SIM_HGCAL-TDRgeom_100k_1040_2023D28_Neutron_XS_100ms_v2/190206_130550/0000/${i}" >> MinBias_13TeV_HGCAL_TDRgeom_100k_2023D28_Neutron_XS_100ms.txt; done

### 2022 #############
######################
# for i in `ls -lrth /lustre/cms/store/user/piet/NeutronBackground/MinBias_Run2_13TeV_GEN_SIM_XS_Run2_2018_100k_1E4s_12X_ProdCutsDefault_v1/crab_MinBias_Run2_13TeV_100k_1E4s_XS_12X_ProdCutsDefault_v1/220322_101204/0000/ | grep root | awk '{print $9}'`; do echo "/store/user/piet/NeutronBackground/MinBias_Run2_13TeV_GEN_SIM_XS_Run2_2018_100k_1E4s_12X_ProdCutsDefault_v1/crab_MinBias_Run2_13TeV_100k_1E4s_XS_12X_ProdCutsDefault_v1/220322_101204/0000/${i}" >> MinBias_Run2_13TeV_100k_Neutron_XS_1E4s.txt; done


# - - -  13X - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### 2026 ### 13X #####
######################
# for i in `ls -lrth /lustre/cms/store/user/piet/NeutronBackground/MinBias_Phase2_14TeV_GEN_SIM_XS_2026D99mod_100k_1E4s_13X_v1/crab_MinBias_Phase2_14TeV_100k_1E4s_XS_13X_v1/230504_162117/0000/ | grep root | awk '{print $9}'`; do echo "/store/user/piet/NeutronBackground/MinBias_Phase2_14TeV_GEN_SIM_XS_2026D99mod_100k_1E4s_13X_v1/crab_MinBias_Phase2_14TeV_100k_1E4s_XS_13X_v1/230504_162117/0000/${i}" >> MinBias_Phase2_14TeV_TuneCP5_100k_Neutron_XS_2026D99_1E4s.txt; done

### 2022 ### 13X #####
######################
for i in `ls -lrth /lustre/cms/store/user/piet/NeutronBackground/MinBias_Phase1_14TeV_GEN_SIM_XS_2023DB_100k_1E4s_13X_v1/crab_MinBias_Phase1_14TeV_100k_1E4s_XS_13X_v1/230512_112600/0000/ | grep root | awk '{print $9}'`; do echo "/store/user/piet/NeutronBackground/MinBias_Phase1_14TeV_GEN_SIM_XS_2023DB_100k_1E4s_13X_v1/crab_MinBias_Phase1_14TeV_100k_1E4s_XS_13X_v1/230512_112600/0000/${i}" >> MinBias_Phase1_14TeV_TuneCP5_100k_Neutron_XS_2023DB_1E4s.txt; done

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
