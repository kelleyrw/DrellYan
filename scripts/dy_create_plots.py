#!/usr/bin/env python

# ------------------#
# imports
# ------------------#

from optparse import OptionParser
import os
import sys
from sys import platform as platform

# ------------------#
# parse inputs
# ------------------#

default_run_list = "json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_cms2.txt"
default_lumi     = 0.082 # fb^-1
default_pset     = "psets/dy_samples_cfg.py"

# parameter options
parser = OptionParser()
parser.add_option("--nevts"      , dest="nevts"      , default=-1               , help="REQUIRED: python configuration file"                 )
parser.add_option("--label"      , dest="label"      , default="test"           , help="unique output label to keep differnet jobs straight" )
parser.add_option("--run_list"   , dest="run_list"   , default=default_run_list , help="good run list (empty == none)"                       )
parser.add_option("--lumi"       , dest="lumi"       , default=default_lumi     , help="luminosity (default 0.082 fb^-1)"                    )
parser.add_option("--sample_pset", dest="sample_pset", default=default_pset     , help="sample pset to use"                                  )

# boolean options
parser.add_option("--test"       , action="store_true"  , dest="test"        , default=False , help="test script -- print commands but do nothing")
parser.add_option("--verbose"    , action="store_true"  , dest="verbose"     , default=False , help="verbose print out"                           )

(options, args) = parser.parse_args()

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

# ---------------------------------------------------------------------------------- #
# make the histograms for a particular sample and signal region
# ---------------------------------------------------------------------------------- #

def MakeHists(sample):

	# start the command
	cmd = "dy_plots "

	# luminosity
	cmd += " --lumi %1.3f" % float(options.lumi)
	
	# the input sample
	cmd += " --sample %s" % sample

	# output label
	cmd += " --label %s" % options.label
	
	# number of events
	cmd += " --nevts %s" % int(options.nevts)

	# sample pset 
	cmd += " --sample_pset %s" % options.sample_pset

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

