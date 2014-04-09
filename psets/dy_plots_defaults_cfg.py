import FWCore.ParameterSet.Config as cms

## process to parse (THIS SHOULD NOT CHANGE)
process = cms.PSet()

## ------------------------------------------------------------- #
## Parameters for the cms2tools_keep_branches utility 
## ------------------------------------------------------------- #

process.dy_plots = cms.PSet(

	## max number of events to run on (-1 means all)
	nevts = cms.untracked.int64(-1),

	## luminosity in fb^-1 
	lumi = cms.untracked.double(0.082),

	## label name -- to give outputs unique dir
	label = cms.untracked.string("v1"),

	## sample name (from dy::Sample)
# 	sample = cms.untracked.string("dyll"),

	## tree name
	tree_name = cms.untracked.string("Events"),
	
	## path to the analysis (can be csv)
# 	input = cms.untracked.string(
# 		"/nfs-7/userdata/rwkelley/cms2/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root"
# 	),

	## output file name 
# 	output = cms.untracked.string("plots/dy_plots.root"),

	## run_list
	run_list = cms.untracked.string("json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_cms2.txt"),

	## verbosity
	verbose = cms.untracked.bool(False)
)
