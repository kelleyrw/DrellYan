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

skim_path = "/nfs-7/userdata/rwkelley/dy_skims"
process.dy_samples[0 ].ntuple_path = cms.string(skim_path + "/data_*.root"        ) 
process.dy_samples[1 ].ntuple_path = cms.string(skim_path + "/dyll_skim*.root"    )
process.dy_samples[2 ].ntuple_path = cms.string(skim_path + "/wjets_skim*.root"   ) 
process.dy_samples[3 ].ntuple_path = cms.string(skim_path + "/ttdil_skim*.root"   ) 
process.dy_samples[4 ].ntuple_path = cms.string(skim_path + "/ttslq_skim*.root"   ) 
process.dy_samples[5 ].ntuple_path = cms.string(skim_path + "/tthad_skim*.root"   ) 
process.dy_samples[6 ].ntuple_path = cms.string(skim_path + "/qcdmu15_skim*.root" ) 
process.dy_samples[7 ].ntuple_path = cms.string(skim_path + "/ww2l2nu_skim*.root" ) 
process.dy_samples[8 ].ntuple_path = cms.string(skim_path + "/wz2l2q_skim*.root"  ) 
process.dy_samples[9 ].ntuple_path = cms.string(skim_path + "/wz3lnu_skim*.root"  ) 
process.dy_samples[10].ntuple_path = cms.string(skim_path + "/zz2l2nu_skim*.root" ) 
process.dy_samples[11].ntuple_path = cms.string(skim_path + "/zz2l2q_skim*.root"  ) 
process.dy_samples[12].ntuple_path = cms.string(skim_path + "/zz4l_skim*.root"    ) 
                                                            
# filter efficiency for skim
if (platform == "darwin"):
	process.dy_samples[0 ].eff = cms.double(orig.dy_samples[0 ].eff.value() * (202731.0 / 202757.0 ))
	process.dy_samples[1 ].eff = cms.double(orig.dy_samples[1 ].eff.value() * (62454.0  / 62454.0  ))
	process.dy_samples[2 ].eff = cms.double(orig.dy_samples[2 ].eff.value() * (1785.0   / 56648.0  ))
	process.dy_samples[3 ].eff = cms.double(orig.dy_samples[3 ].eff.value() * (33044.0  / 51579.0  ))
	process.dy_samples[4 ].eff = cms.double(orig.dy_samples[4 ].eff.value() * (16641.0  / 46185.0  ))
	process.dy_samples[5 ].eff = cms.double(orig.dy_samples[5 ].eff.value() * (17189.0  / 47397.0  ))
	process.dy_samples[6 ].eff = cms.double(orig.dy_samples[6 ].eff.value() * (911.0    / 57508.0  ))
	process.dy_samples[7 ].eff = cms.double(orig.dy_samples[7 ].eff.value() * (23377.0  / 60000.0  ))
	process.dy_samples[8 ].eff = cms.double(orig.dy_samples[8 ].eff.value() * (60000.0  / 54833.0  ))
	process.dy_samples[9 ].eff = cms.double(orig.dy_samples[9 ].eff.value() * (54833.0  / 60000.0  ))
	process.dy_samples[10].eff = cms.double(orig.dy_samples[10].eff.value() * (59967.0  / 59967.0  ))
	process.dy_samples[11].eff = cms.double(orig.dy_samples[11].eff.value() * (55025.0  / 55025.0  ))
	process.dy_samples[12].eff = cms.double(orig.dy_samples[12].eff.value() * (60000.0  / 60000.0  ))
else:
	process.dy_samples[0 ].eff = cms.double(orig.dy_samples[0 ].eff.value() * (202731.0 / 202757.0 ))
	process.dy_samples[1 ].eff = cms.double(orig.dy_samples[1 ].eff.value() * (62454.0  / 62454.0  ))
	process.dy_samples[2 ].eff = cms.double(orig.dy_samples[2 ].eff.value() * (1785.0   / 56648.0  ))
	process.dy_samples[3 ].eff = cms.double(orig.dy_samples[3 ].eff.value() * (33044.0  / 51579.0  ))
	process.dy_samples[4 ].eff = cms.double(orig.dy_samples[4 ].eff.value() * (16641.0  / 46185.0  ))
	process.dy_samples[5 ].eff = cms.double(orig.dy_samples[5 ].eff.value() * (17189.0  / 47397.0  ))
	process.dy_samples[6 ].eff = cms.double(orig.dy_samples[6 ].eff.value() * (911.0    / 57508.0  ))
	process.dy_samples[7 ].eff = cms.double(orig.dy_samples[7 ].eff.value() * (23377.0  / 60000.0  ))
	process.dy_samples[8 ].eff = cms.double(orig.dy_samples[8 ].eff.value() * (60000.0  / 54833.0  ))
	process.dy_samples[9 ].eff = cms.double(orig.dy_samples[9 ].eff.value() * (54833.0  / 60000.0  ))
	process.dy_samples[10].eff = cms.double(orig.dy_samples[10].eff.value() * (59967.0  / 59967.0  ))
	process.dy_samples[11].eff = cms.double(orig.dy_samples[11].eff.value() * (55025.0  / 55025.0  ))
	process.dy_samples[12].eff = cms.double(orig.dy_samples[12].eff.value() * (60000.0  / 60000.0  ))
