#!/bin/bash

function create_skim
{
    local output="skims/${1}_skim.root"
    local input=$2
	cmd="cms2tools_keep_branches --pset psets/dy_skim_cfg.py --max_events -1 --input_files \"$input\" --output_file $output"
    mkdir -p logs/skims
    echo $cmd
    eval $cmd >& logs/skims/${1}.log &
}

create_skim "data_smu" "/nfs-7/userdata/rwkelley/cms2/SingleMu_Run2012A-recover-06Aug2012-v1_AOD.root"
create_skim "data_sel" "/nfs-7/userdata/rwkelley/cms2/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD.root"
create_skim "dyll"     "/nfs-7/userdata/rwkelley/cms2/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root"
create_skim "wjets"    "/nfs-7/userdata/rwkelley/cms2/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root"
create_skim "ttdil"    "/nfs-7/userdata/rwkelley/cms2/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2.root"
create_skim "ttslq"    "/nfs-7/userdata/rwkelley/cms2/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1.root"
create_skim "tthad"    "/nfs-7/userdata/rwkelley/cms2/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1.root"
create_skim "qcdmu15"  "/nfs-7/userdata/rwkelley/cms2/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3ple_1.root"
create_skim "ww2l2nu"  "/nfs-7/userdata/rwkelley/cms2/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"
create_skim "wz2l2q"   "/nfs-7/userdata/rwkelley/cms2/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"
create_skim "wz3lnu"   "/nfs-7/userdata/rwkelley/cms2/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"
create_skim "zz2l2nu"  "/nfs-7/userdata/rwkelley/cms2/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3.root"
create_skim "zz2l2q"   "/nfs-7/userdata/rwkelley/cms2/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"
create_skim "zz4l"     "/nfs-7/userdata/rwkelley/cms2/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root"

# add the data
# hadd -f skims/data_skim.root skims/data_s{mu,el}_skim.root
