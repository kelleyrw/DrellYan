#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "CTable.h"
#include "TCanvas.h"
#include "TEventList.h"
#include "TCut.h"
#include <limits>
#include <string>
#include <iostream>

static const double s_SDEfficiency[] =
{
     /*data    =*/ (1.0                  ),
     /*dyll    =*/ (27137253.0/30459500.0),
     /*wjets   =*/ (14890630.0/18393090.0),
     /*ttdil   =*/ (11902045.0/12119013.0),
     /*ttslq   =*/ (23453443.0/25384818.0),
     /*tthad   =*/ (21211839.0/31223821.0),
     /*qcd     =*/ (19770457.0/20614602.0),
     /*ww2l2q  =*/ (1933235.0 /1933235.0 ),
     /*wz2l2q  =*/ (2937874.0 /3215990.0 ),
     /*wz3l    =*/ (2017979.0 /2017979.0 ),
     /*zz2l2q  =*/ (1777571.0 /1936727.0 ),
     /*zz2l2nu =*/ (857982.0  /954911.0  ),
     /*zz4l    =*/ (4807893.0 /4807893.0 )
};
static const size_t s_SamplesCount = sizeof(s_SDEfficiency)/sizeof(s_SDEfficiency[0]);

void FillHists(rt::TH1Container& hc, TChain& chain, const std::string& sample_name, long num_events = std::numeric_limits<long>::max())
{
    // number of events
    num_events = (chain.GetEntries() < num_events ? chain.GetEntries() : num_events);
    std::cout << "Beginning " << sample_name << " on " << num_events << " of " << chain.GetEntries() << " events..." << std::endl; 

     // scale the # of events
    const float nevts_aod   = chain.GetMaximum("evt_nEvts"); // # files in AOD
    const float nevts_file  = chain.GetEntries();            // # files in this sample
    const float nevts_scale = nevts_aod/nevts_file;          // scale factor to "fix" evt_scale1fb since we are using a subset of the data 

    // lumi
    static const float lumi = 0.082; //fb^-1

    // efficiency scale
    float eff_scale = 1.0;
    if (sample_name == "data"   ) {eff_scale = s_SDEfficiency[ 0];}
    if (sample_name == "dyll"   ) {eff_scale = s_SDEfficiency[ 1];}
    if (sample_name == "wjets"  ) {eff_scale = s_SDEfficiency[ 2];}
    if (sample_name == "ttdil"  ) {eff_scale = s_SDEfficiency[ 3];}
    if (sample_name == "ttslq"  ) {eff_scale = s_SDEfficiency[ 4];}
    if (sample_name == "tthad"  ) {eff_scale = s_SDEfficiency[ 5];}
    if (sample_name == "qcd"    ) {eff_scale = s_SDEfficiency[ 6];}
    if (sample_name == "ww"     ) {eff_scale = s_SDEfficiency[ 7];}
    if (sample_name == "wz2l2q" ) {eff_scale = s_SDEfficiency[ 8];}
    if (sample_name == "wz3l"   ) {eff_scale = s_SDEfficiency[ 9];}
    if (sample_name == "zz2l2q" ) {eff_scale = s_SDEfficiency[10];}
    if (sample_name == "zz2l2nu") {eff_scale = s_SDEfficiency[11];}
    if (sample_name == "zz4l"   ) {eff_scale = s_SDEfficiency[12];}

    // overall scale
    const float scale = lumi * eff_scale * nevts_scale; 
    std::cout << "scale factor = " << Form("%f * %f * %f = %f", lumi, eff_scale, nevts_scale, scale) << std::endl;

    // book hists
    hc.Add(new TH1F(Form("h_count_mm_%s", sample_name.c_str()), Form("Event ee count (%s)"                               , sample_name.c_str()), 3 , -0.5, 2.5));
    hc.Add(new TH1F(Form("h_count_ee_%s", sample_name.c_str()), Form("Event #mu#mu count (%s)"                           , sample_name.c_str()), 3 , -0.5, 2.5));
    hc.Add(new TH1F(Form("h_mee_%s"     , sample_name.c_str()), Form("Dimuon Electron (%s);m_{ee}(GeV);# Events Expected", sample_name.c_str()), 60,   60, 120));
    hc.Add(new TH1F(Form("h_mmm_%s"     , sample_name.c_str()), Form("Dimuon Mass (%s);m_{#mu#mu}(GeV);# Events Expected", sample_name.c_str()), 60,   60, 120));
    hc.Sumw2();

    // alias for dilepton mass
    chain.SetAlias("mass_mm_px", "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.Px() : 0.0)");
    chain.SetAlias("mass_mm_py", "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.Py() : 0.0)");
    chain.SetAlias("mass_mm_pz", "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.Pz() : 0.0)");
    chain.SetAlias("mass_mm_e" , "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.E()  : 0.0)");
    chain.SetAlias("mass_mm"   , "sqrt(mass_mm_e*mass_mm_e - mass_mm_px*mass_mm_px - mass_mm_py*mass_mm_py - mass_mm_pz*mass_mm_pz)");

    chain.SetAlias("mass_ee_px", "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.Px() : 0.0)");
    chain.SetAlias("mass_ee_py", "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.Py() : 0.0)");
    chain.SetAlias("mass_ee_pz", "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.Pz() : 0.0)");
    chain.SetAlias("mass_ee_e" , "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.E()  : 0.0)");
    chain.SetAlias("mass_ee"   , "sqrt(mass_ee_e*mass_ee_e - mass_ee_px*mass_ee_px - mass_ee_py*mass_ee_py - mass_ee_pz*mass_ee_pz)");

    // fill hist
    hc.SetDirectory(gDirectory);
    const TCut selection_ee = Form("%1.4f*evt_scale1fb*(Sum$(genps_status==3 && genps_id==11)>=1 && Sum$(genps_status==3 && genps_id==-11)>=1)", scale);
    chain.Draw(">>event_list_ee", selection_ee, "goff", num_events);
    TEventList * const event_list_ee = dynamic_cast<TEventList*>(gDirectory->Get("event_list_ee"));
    chain.SetEventList(event_list_ee);
    std::cout << "filling ee hists..." << std::endl;
    chain.Draw(Form("1>>h_count_ee_%s" , sample_name.c_str()), selection_ee, "goff", num_events);
    chain.Draw(Form("mass_ee>>h_mee_%s", sample_name.c_str()), selection_ee, "goff", num_events);

    const TCut selection_mm = Form("%1.4f*evt_scale1fb*(Sum$(genps_status==3 && genps_id==13)>=1 && Sum$(genps_status==3 && genps_id==-13)>=1)", scale);
    chain.SetEventList(NULL);
    chain.Draw(">>event_list_mm", selection_mm, "goff", num_events);
    TEventList * const event_list_mm = dynamic_cast<TEventList*>(gDirectory->Get("event_list_mm"));
    chain.SetEventList(event_list_mm);
    std::cout << "filling mm hists..." << std::endl;
    chain.Draw(Form("1>>h_count_mm_%s" , sample_name.c_str()), selection_mm, "goff", num_events);
    chain.Draw(Form("mass_mm>>h_mmm_%s", sample_name.c_str()), selection_mm, "goff", num_events);
    hc.SetDirectory(NULL);

    std::cout << "Complete " << sample_name  << ":  ";
    std::cout << "mm count = " << rt::Integral(hc["h_count_mm_"+sample_name]) << " (" << hc["h_count_mm_"+sample_name]->GetEntries() << ") : "; 
    std::cout << "ee count = " << rt::Integral(hc["h_count_ee_"+sample_name]) << " (" << hc["h_count_ee_"+sample_name]->GetEntries() << ")\n" << std::endl; 
   
    // done
    return;
}

void create_dilep_count(const long num_events = std::numeric_limits<long>::max())
{
    // get data
    TChain ch_dyll("Events");
    ch_dyll.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root");

    TChain ch_data("Events");
    ch_data.Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root");
    ch_data.Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root");

    TChain ch_wjets("Events");
    ch_wjets.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-28/merged_ntuple_1[0-9].root");

    TChain ch_ttdil("Events");
    ch_ttdil.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-24/merged_ntuple_1[0-9].root");

    TChain ch_ttslq("Events");
    ch_ttslq.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/V05-03-24/merged_ntuple_1[0-9].root");

    TChain ch_tthad("Events");
    ch_tthad.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/V05-03-24/merged_ntuple_1[0-9].root");

    TChain ch_qcd("Events");
    ch_qcd.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-18_slim/merged_ntuple_1[0-9].root");

    TChain ch_ww("Events");
    ch_ww.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root");

    TChain ch_wz2l2q("Events");
    ch_wz2l2q.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root");

    TChain ch_wz3l("Events");
    ch_wz3l.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root");

    TChain ch_zz2l2q("Events");
    ch_zz2l2q.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-23/merged_ntuple_1[0-9].root");

    TChain ch_zz2l2nu("Events");
    ch_zz2l2nu.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root");

    TChain ch_zz4l("Events");
    ch_zz4l.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root");

    // get hists
    rt::TH1Container hc;
    FillHists(hc, ch_dyll    , "dyll"    , num_events);
    FillHists(hc, ch_wjets   , "wjets"   , num_events);
    FillHists(hc, ch_ttdil   , "ttdil"   , num_events);
    FillHists(hc, ch_ttslq   , "ttslq"   , num_events);
    FillHists(hc, ch_ttdil   , "tthad"   , num_events);
    FillHists(hc, ch_qcd     , "qcd"     , num_events);
    FillHists(hc, ch_ww      , "ww"      , num_events);
    FillHists(hc, ch_wz2l2q  , "wz2l2q"  , num_events);
    FillHists(hc, ch_wz3l    , "wz3l"    , num_events);
    FillHists(hc, ch_zz2l2q  , "zz2l2q"  , num_events);
    FillHists(hc, ch_zz2l2nu , "zz2l2nu" , num_events);
    FillHists(hc, ch_zz4l    , "zz4l"    , num_events);

    // write output
    hc.Write("plots/dilep_counts.root");
}

void print_dilep_count()
{
    rt::TH1Container hc("plots/dilep_counts.root");

    const float wz_mm = rt::Integral(hc["h_count_mm_wz2l2q" ]) + rt::Integral(hc["h_count_mm_wz3l"]);
    const float wz_ee = rt::Integral(hc["h_count_ee_wz2l2q" ]) + rt::Integral(hc["h_count_ee_wz3l"]);
    const float zz_mm = rt::Integral(hc["h_count_mm_zz2l2q" ]) + rt::Integral(hc["h_count_mm_zz2l2nu"]) + rt::Integral(hc["h_count_mm_zz4l"]);
    const float zz_ee = rt::Integral(hc["h_count_ee_zz2l2q" ]) + rt::Integral(hc["h_count_ee_zz2l2nu"]) + rt::Integral(hc["h_count_ee_zz4l"]);
    
    CTable t1;
    t1.useTitle();
    t1.setTitle("Generator level number of \\ell\\ell events (no selections)");
    t1.setTable() (                                                 "mm ",                                   "ee")
                  ( "DY --> mm"   , rt::Integral(hc["h_count_mm_dyll"   ]), rt::Integral(hc["h_count_ee_dyll"   ]))
                  ( "W --> mnu"   , rt::Integral(hc["h_count_mm_wjets"  ]), rt::Integral(hc["h_count_ee_wjets"  ]))
                  ( "tt"          , rt::Integral(hc["h_count_mm_ttdil"  ]), rt::Integral(hc["h_count_ee_ttdil"  ]))
                  ( "QCD"         , rt::Integral(hc["h_count_mm_qcd"    ]), rt::Integral(hc["h_count_ee_qcd"    ]))
                  ( "WW"          , rt::Integral(hc["h_count_mm_ww"     ]), rt::Integral(hc["h_count_ee_ww"     ]))
                  ( "WZ"          , wz_mm                                 , wz_ee                                 )
                  ( "WZ --> 2l2q" , rt::Integral(hc["h_count_mm_wz2l2q" ]), rt::Integral(hc["h_count_ee_wz2l2q" ]))
                  ( "WZ --> 3lnu" , rt::Integral(hc["h_count_mm_wz3l"   ]), rt::Integral(hc["h_count_ee_wz3l"   ]))
                  ( "ZZ"          , zz_mm                                 , zz_ee                                 ) 
                  ( "ZZ --> 2l2q" , rt::Integral(hc["h_count_mm_zz2l2q" ]), rt::Integral(hc["h_count_ee_zz2l2q" ]))
                  ( "ZZ --> 2lnu" , rt::Integral(hc["h_count_mm_zz2l2nu"]), rt::Integral(hc["h_count_ee_zz2l2nu"]))
                  ( "ZZ --> 4l"   , rt::Integral(hc["h_count_mm_zz4l"   ]), rt::Integral(hc["h_count_ee_zz4l"   ]))
                  ;
    std::cout << t1 << std::endl;
}

void print_dilep_plots(const std::string& suffix = "png")
{
    rt::TH1Container hc("plots/dilep_counts.root");

    // formatting
    rt::TH1Overlay::legend_height_per_entry_default = 0.050;
    rt::TH1Overlay::legend_text_size_default        = 0.025;
    const float ymin = 1e0;
    const float ymax = 3e4;
    const std::string option = "sb::off dt::stack lg::top_left";

    rt::TH1Overlay p1("Generator Level Dimuon distributions - #intL = 82 pb^{-1};m_{#mu^{+}#mu^{-}} (GeV/c^{2});Events/1.0 GeV/c^{2}", option);
    p1.Add(hc["h_mmm_dyll"   ],  "DY #rightarrow #mu#mu"        , kOrange-2);
    p1.Add(hc["h_mmm_wz2l2q" ],  "WZ #rightarrow 2l2q"          , kRed+2   );
    p1.Add(hc["h_mmm_wz3l"   ],  "WZ #rightarrow 3l#nu"         , kRed-7   );
    p1.Add(hc["h_mmm_zz2l2q" ],  "ZZ #rightarrow 2l2q"          , kBlue    );
    p1.Add(hc["h_mmm_zz2l2nu"],  "ZZ #rightarrow 2l#nu"         , kBlue+2  );
    p1.Add(hc["h_mmm_zz4l"   ],  "ZZ #rightarrow 4l"            , kBlue-5  );
    p1.Add(hc["h_mmm_ttdil"  ],  "t#bar{t} #rightarrow 2l2#nu2b", kOrange+7);
    p1.Add(hc["h_mmm_ww"     ],  "WW #rightarrow 2l2#nu"        , kGreen+2 );
    p1.Add(hc["h_mmm_wjets"  ],  "W #rightarrow l#nu"           , kGray+1  );
    p1.Add(hc["h_mmm_qcd"    ],  "QCD"                          , kMagenta );
    p1.SetYAxisRange(ymin, ymax);

    rt::TH1Overlay p2("Generator Level Dielectron distributions - #intL = 82 pb^{-1};m_{e^{+}e^{-}} (GeV/c^{2});Events/1.0 GeV/c^{2}", option);
    p2.Add(hc["h_mee_dyll"   ],  "DY #rightarrow ee"            , kOrange-2);
    p2.Add(hc["h_mee_wz2l2q" ],  "WZ #rightarrow 2l2q"          , kRed+2   );
    p2.Add(hc["h_mee_wz3l"   ],  "WZ #rightarrow 3l#nu"         , kRed-7   );
    p2.Add(hc["h_mee_zz2l2q" ],  "ZZ #rightarrow 2l2q"          , kBlue    );
    p2.Add(hc["h_mee_zz2l2nu"],  "ZZ #rightarrow 2l#nu"         , kBlue+2  );
    p2.Add(hc["h_mee_zz4l"   ],  "ZZ #rightarrow 4l"            , kBlue-5  );
    p2.Add(hc["h_mee_ttdil"  ],  "t#bar{t} #rightarrow 2l2#nu2b", kOrange+7);
    p2.Add(hc["h_mee_ww"     ],  "WW #rightarrow 2l2#nu"        , kGreen+2 );
    p2.Add(hc["h_mee_wjets"  ],  "W #rightarrow l#nu"           , kGray+1  );
    p2.Add(hc["h_mee_qcd"    ],  "QCD"                          , kMagenta );
    p2.SetYAxisRange(ymin, ymax);

    // output
    rt::Print(p1, "plots/p_gen_mmm", suffix);
    rt::Print(p2, "plots/p_gen_mee", suffix);

    // log output
    p1.SetLogy(true);
    rt::Print(p1, "plots/p_gen_mmm_log", suffix);
    p2.SetLogy(true);
    rt::Print(p2, "plots/p_gen_mee_log", suffix);
}

void quick_dilep_count(const long num_events = std::numeric_limits<long>::max())
{
    create_dilep_count(num_events);
    print_dilep_count();
    print_dilep_plots();
}
