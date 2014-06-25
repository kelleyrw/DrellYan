#!/bin/bash

git clone https://github.com/kelleyrw/AnalysisTools.git 
# git checkout -t dy_V00-00-01 
git clone https://github.com/cmstas/Dictionaries.git CMS2/Dictionaries
git clone https://github.com/cmstas/CORE.git CMS2/NtupleMacrosCore
cd CMS2/NtupleMacrosCore
git pull
git checkout scram_compliant
./setup/setupCoreForSCRAM.sh setup/cms2_ntuple_postprocess_05.03.23.root 
cd $CMSSW_BASE/src
scram b -j20
