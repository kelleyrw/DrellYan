import FWCore.ParameterSet.Config as cms
import sys
import os

## add the configuration path (THIS SHOULD NOT CHANGE)
sys.path.append(os.getenv("CMSSW_BASE") + "/src/Analysis/DrellYan/psets")
from dy_samples_cfg import process

## ------------------------------------------------------------- #
## Parameters for the samples 
## ------------------------------------------------------------- #

data_path = "/nfs-7/userdata/rwkelley/cms2"
mc_path   = "/nfs-7/userdata/rwkelley/cms2"

process.dy_samples[0 ].ntuple_path = cms.string(data_path + "/Single*_Run2012A-recover-06Aug2012-v1_AOD.root"                                              ) 
process.dy_samples[1 ].ntuple_path = cms.string(mc_path   + "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root"  ) 
process.dy_samples[2 ].ntuple_path = cms.string(mc_path   + "/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root"       ) 
process.dy_samples[3 ].ntuple_path = cms.string(mc_path   + "/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2.root"             ) 
process.dy_samples[4 ].ntuple_path = cms.string(mc_path   + "/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1.root"         ) 
process.dy_samples[5 ].ntuple_path = cms.string(mc_path   + "/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1.root"         ) 
process.dy_samples[6 ].ntuple_path = cms.string(mc_path   + "/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3.root" ) 
process.dy_samples[7 ].ntuple_path = cms.string(mc_path   + "/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"     ) 
process.dy_samples[8 ].ntuple_path = cms.string(mc_path   + "/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"      ) 
process.dy_samples[9 ].ntuple_path = cms.string(mc_path   + "/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"          ) 
process.dy_samples[10].ntuple_path = cms.string(mc_path   + "/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3.root"     ) 
process.dy_samples[11].ntuple_path = cms.string(mc_path   + "/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"      ) 
process.dy_samples[12].ntuple_path = cms.string(mc_path   + "/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"        ) 
