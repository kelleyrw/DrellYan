import FWCore.ParameterSet.Config as cms
import sys
import os
import copy
from sys import platform as platform

## add the configuration path (THIS SHOULD NOT CHANGE)
sys.path.append(os.getenv("CMSSW_BASE") + "/src/Analysis/DrellYan/psets")
from dy_samples_cfg import process as orig
process = copy.deepcopy(orig)

## ------------------------------------------------------------- #
## Parameters for the samples 
## ------------------------------------------------------------- #

data_path = "/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24"
process.dy_samples[0 ].ntuple_path = cms.string(data_path+ "/DoubleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root,"+data_path+"/DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root") 
