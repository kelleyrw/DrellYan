import FWCore.ParameterSet.Config as cms

## process to parse (THIS SHOULD NOT CHANGE)
process = cms.PSet()

## ------------------------------------------------------------- #
## Parameters for the cms2tools_keep_branches utility 
## ------------------------------------------------------------- #

process.cms2tools_drop_branches = cms.PSet(

	## max number of events to run on (-1 means all)
	max_events = cms.int64(-1),

	## tree name
	tree_name = cms.string("Events"),
	
	## path to the analysis
	input_files = cms.vstring(
		"/nfs-7/userdata/rwkelley/cms2/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
	),

	## output file name 
	output_file = cms.string("skims/dy_skim.root"),

	## selection (same as TTree::Draw)
	selection = cms.string(
		"(!evt_isRealData && Sum$(genps_status==3 && genps_id==11)>=1 && Sum$(genps_status==3 && genps_id==-11)>=1) || "
		"(!evt_isRealData && Sum$(genps_status==3 && genps_id==13)>=1 && Sum$(genps_status==3 && genps_id==-13)>=1) || "
		"(!evt_isRealData && Sum$(genps_status==3 && genps_id==15)>=1 && Sum$(genps_status==3 && genps_id==-15)>=1) || "
		"!hyp_type@.empty()"
	),

	## aliases to branches to keep (overides anything in drop_alias_names)
	## have to use regegular expression
	## see: http://www.cplusplus.com/reference/regex/ECMAScript/
	keep_alias_names = cms.vstring(
		"evt_xsec_excl",
# 		"evt_CMS2tag",
# 		"evt_dataset",
		"evt_event",
		"evt_filt_eff",
		"evt_isRealData",
		"evt_kt6calo_rho",
		"evt_lumiBlock",
		"evt_met",
		"evt_metPhi",
		"evt_nEvts",
		"evt_pfmet",
		"evt_pfmetPhi",
		"evt_kt6pf_foregiso_rho",
		"evt_run",
		"evt_scale1fb",
		"gsftrks_charge",
		"gsftrks_d0",
		"gsftrks_d0Err",
		"gsftrks_d0phiCov",
		"gsftrks_etaErr",
		"gsftrks_p4",
		"gsftrks_phiErr",
		"gsftrks_z0",
		"gsftrks_z0Err",
		"genps_id",
		"genps_id_mother",
		"genps_lepdaughter_id",
		"genps_lepdaughter_idx",
		"genps_lepdaughter_p4",
		"genps_p4",
		"genps_signalProcessID",
		"genps_status",
		"genps_weight",
		"hyp_(.*)",
		"mus_charge",
		"mus_chi2",
		"mus_d0",
		"mus_d0Err",
		"mus_gfit(.*)",
		"mus_ip3d",
		"mus_ip3derr",
		"mus_isoR04(.*)",
		"mus_mc(.*)",
		"mus_ndof",
		"mus_numberOfMatchedStations",
		"mus_p4",
		"mus_trk_charge",
		"mus_trk_p4",
		"mus_trkidx",
		"mus_type",
		"mus_z0",
		"mus_z0Err",
		"els_charge",
		"els_chi2",
		"els_dEtaIn",
		"els_dEtaOut",
		"els_dPhiIn",
		"els_dPhiInPhiOut",
		"els_dPhiOut",
		"els_e3x3",
		"els_e5x5",
		"els_eOverPIn",
		"els_eOverPOut",
		"els_ecalEnergy",
		"els_eSC",
		"els_etaSC",
		"els_exp_innerlayers",
		"els_fbrem",
		"els_fiduciality",
		"els_gsftrkidx",
		"els_hOverE",
		"els_inner_position",
		"els_ip3d",
		"els_ip3derr",
		"els_iso03_pf2012ext_(.*)",
		"els_mc(.*)",
		"els_ndof",
		"els_p4",
		"els_sigmaIEtaIEta",
		"els_trk_charge",
		"els_trk_p4",
		"els_trkidx",
		"els_z0",
		"els_z0Err",
		"trks_(.*)",
		"vtxs_(.*)",
		"hlt_(.*)",
		"convs_(.*)",
	),
)
