#!/usr/bin/env python

# ------------------#
# imports
# ------------------#

from optparse import OptionParser
import os
import sys

# ------------------#
# parse inputs
# ------------------#

default_pset     = "psets/dy_plots_defaults_cfg.py"
default_run_list = "json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_cms2.txt"
default_lumi     = 0.082 # fb^-1

# parameter options
parser = OptionParser()
parser.add_option("--pset"     , dest="pset"      , default=default_pset     , help="path to the razor cards"                             )
parser.add_option("--nevts"    , dest="nevts"     , default=-1               , help="REQUIRED: python configuration file"                 )
parser.add_option("--label"    , dest="label"     , default="test"           , help="unique output label to keep differnet jobs straight" )
parser.add_option("--run_list" , dest="run_list"  , default=default_run_list , help="good run list (empty == none)"                       )
parser.add_option("--lumi"     , dest="lumi"      , default=default_lumi     , help="luminosity (default 0.082 fb^-1)"                    )

# boolean options
parser.add_option("--test"       , action="store_true"  , dest="test"        , default=False , help="test script -- print commands but do nothing")
parser.add_option("--verbose"    , action="store_true"  , dest="verbose"     , default=False , help="verbose print out"                           )
# parser.add_option("--no_hist"    , action="store_true"  , dest="no_hist"     , default=False , help="do not create histograms, do everything else"   )

(options, args) = parser.parse_args()

# check for validity
def CheckOptions():
	# pset
	if (not options.pset or not os.path.isfile(options.pset)):
		raise Exception("required pset is missing")

	return

# ------------------#
# available samples 
# ------------------#

samples = [
	"data",
	"dyll",
	"wjets",
	"ttdil",
	"ttslq",
	"tthad",
	"qcdmu15",
	"ww2l2nu",
	"wz2l2q",
	"wz3lnu",
	"zz2l2nu",
	"zz2l2q",
	"zz4l",
]

# mac file names (hack -- FIXME)
from sys import platform as platform
mac_sample_path = "/nfs-7/userdata/rwkelley/cms2" 
mac_sample_files = {
	"data"    : "%s/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD.root,%s/SingleMu_Run2012A-recover-06Aug2012-v1_AOD.root" % (mac_sample_path, mac_sample_path), 
	"dyll"    : "%s/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root"               % (mac_sample_path), 
	"wjets"   : "%s/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root"                    % (mac_sample_path), 
	"ttdil"   : "%s/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2.root"                          % (mac_sample_path), 
	"ttslq"   : "%s/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1.root"                      % (mac_sample_path), 
	"tthad"   : "%s/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1.root"                      % (mac_sample_path), 
	"qcdmu15" : "%s/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3ple_1.root"         % (mac_sample_path), 
	"ww2l2nu" : "%s/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"                  % (mac_sample_path), 
	"wz2l2q"  : "%s/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"                       % (mac_sample_path), 
	"wz3lnu"  : "%s/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"                   % (mac_sample_path), 
	"zz2l2nu" : "%s/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3.root"                  % (mac_sample_path), 
	"zz2l2q"  : "%s/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3.root"                  % (mac_sample_path), 
	"zz4l"    : "%s/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"                     % (mac_sample_path), 
}

# ---------------------------------------------------------------------------------- #
# make the histograms for a particular sample and signal region
# ---------------------------------------------------------------------------------- #

def MakeHists(sample):

	# start the command
	cmd = "dy_plots --pset %s" % options.pset

	# luminosity
	cmd += " --lumi %1.3f" % float(options.lumi)
	
	# the input sample
	cmd += " --sample %s" % sample

	# output label
	cmd += " --label %s" % options.label
	
	# number of events
	cmd += " --nevts %s" % int(options.nevts)

	# sample input file (if mac)
	if (platform == "darwin"):
		cmd += " --input \"%s\"" % mac_sample_files[sample]

	# verbose	
	if (options.verbose):
		cmd += " --verbose 1"

	# logname
	log_dir_name  = "logs/%s" % options.label
	log_file_name = "%s/%s.log" % (log_dir_name, sample)
	cmd += " >& %s &" % log_file_name
	if (not options.test and not os.path.exists(log_dir_name)):
		os.makedirs(log_dir_name)		

	print "[dy_plots] making plots for %s..." % (sample)
	print cmd
	if (not options.test):
		os.system(cmd)

	return	

# ------------------#
# "main program" 
# ------------------#

def main():
	try:
		# check the options
		CheckOptions()

		# samples to run on
		for sample in samples:	
			MakeHists(sample)

	except Exception, e:
		print "[dy_create_plots] ERROR:", e
		return 1
	except:
		print "[dy_create_plots] ERROR:", sys.exc_info()[0]
		return 1

# do it
if __name__ == '__main__':
	main()

