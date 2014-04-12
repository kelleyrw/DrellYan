
import FWCore.ParameterSet.Config as cms

## process to parse (THIS SHOULD NOT CHANGE)
process = cms.PSet()

## ------------------------------------------------------------- #
## Parameters for the cms2tools_keep_branches utility 
## ------------------------------------------------------------- #

data_path = "/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24"
mc_path   = "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC"

process.dy_samples = cms.VPSet(
    cms.PSet(
    	name        = cms.string("data"),
    	title       = cms.string("data"),
    	latex       = cms.string("data"),
    	ntuple_path = cms.string(data_path+ "/Single*_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root"),
    	color       = cms.int32 (1), # kBlack
    	eff         = cms.double(1.0)
    ),
    cms.PSet(
    	name        = cms.string("dyll"),
    	title       = cms.string("Z/#gamma #rightarrow l^{+}l^{-}"),
    	latex       = cms.string("Z/#gamma \\rightarrow \\ell \\ell$"),
    	ntuple_path = cms.string(mc_path + "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (798), # kOrange-2
    	eff         = cms.double(27137253.0/30459500.0)
    ),
    cms.PSet(
    	name        = cms.string("wjets"), 
    	title       = cms.string("W+jets #rightarrow l#nu"), 
    	latex       = cms.string("$W+jets \\rightarrow \\ell \\nu$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-28/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (921), #kGray+1,
    	eff         = cms.double(14890630.0/18393090.0),
    ),
    cms.PSet(
    	name        = cms.string("ttdil"), 
    	title       = cms.string("t#bar{t} #rightarrow llX"), 
    	latex       = cms.string("$t\\overline{t} \\rightarrow \\ell \\ell X$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-24/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (807), #kOrange+7,
    	eff         = cms.double(11902045.0/12119013.0),
    ),
    cms.PSet(
    	name        = cms.string("ttslq"), 
    	title       = cms.string("t#bar{t} #rightarrow l(q #rightarrow l)X"), 
    	latex       = cms.string("$t\\overline{t} \\rightarrow \\ell (q \\rightarrow \\ell) X$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/V05-03-24/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (806), #kOrange+6,
    	eff         = cms.double(23453443.0/25384818.0),
    ),
    cms.PSet(
    	name        = cms.string("tthad"), 
    	title       = cms.string("t#bar{t} #rightarrow hadrons"), 
    	latex       = cms.string("$t\\overline{t} \\rightarrow hadrons$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/V05-03-24/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (803), #kOrange+3,
    	eff         = cms.double(21211839.0/31223821.0),
    ),
    cms.PSet(
    	name        = cms.string("qcdmu15"), 
    	title       = cms.string("QCD (#mu15 enriched)"), 
    	latex       = cms.string("QCD ($\\mu$15 enriched)"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-18_slim/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (616), #kMagenta,
    	eff         = cms.double(19770457.0/20614602.0),
    ),
    cms.PSet(
    	name        = cms.string("ww2l2nu"), 
    	title       = cms.string("WW #rightarrow 2l + 2#nu"), 
    	latex       = cms.string("$WW \\rightarrow 2\\ell + 2\\nu$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (418), #kGreen+2,
    	eff         = cms.double(1933235.0/1933235.0),
    ),
    cms.PSet(
    	name        = cms.string("wz2l2q"), 
    	title       = cms.string("WZ #rightarrow 2l + 2q"), 
    	latex       = cms.string("$WZ \\rightarrow 2\\ell + 2q$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (634), #kRed+2,
    	eff         = cms.double(2937874.0/3215990.0),
    ),
    cms.PSet(
    	name        = cms.string("wz3lnu"), 
    	title       = cms.string("WZ #rightarrow 3l + #nu"), 
    	latex       = cms.string("$WZ \\rightarrow 3\\ell + \\nu$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (625), #kRed-7,
    	eff         = cms.double(2937874.0/3215990.0),
    ),
    cms.PSet(
    	name        = cms.string("zz2l2nu"), 
    	title       = cms.string("ZZ #rightarrow 2l + 2#nu"), 
    	latex       = cms.string("$ZZ \\rightarrow 2\\ell + 2\\nu$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-23/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (600), #kBlue,
    	eff         = cms.double(857982.0/954911.0),
    ),
    cms.PSet(
    	name        = cms.string("zz2l2q"), 
    	title       = cms.string("ZZ #rightarrow 2l + 2q"), 
    	latex       = cms.string("$ZZ \\rightarrow 2\\ell + 2q$"), 
    	ntuple_path = cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root"),
    	color       = cms.int32 (602), #kBlue+2,
    	eff         = cms.double(1777571.0/1936727.0),
    ),
    cms.PSet(
    	name        =cms.string("zz4l"), 
    	title       =cms.string("ZZ #rightarrow 4l"), 
    	latex       =cms.string("$ZZ \\rightarrow 4\\ell$"), 
    	ntuple_path =cms.string("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root"),
    	color       =cms.int32 (595), #kBlue-5,
    	eff         =cms.double(4807893.0/4807893.0),
	),
)
