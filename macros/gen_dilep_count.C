#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "Analysis/DrellYan/interface/Sample.h"
#include "CTable.h"
#include "TCanvas.h"
#include "TEventList.h"
#include "TCut.h"
#include <limits>
#include <string>
#include <iostream>

void FillHists(rt::TH1Container& hc, const dy::Sample::value_type sample, long num_events = std::numeric_limits<long>::max())
{
    // number of events
    const dy::Sample::Info sample_info = dy::GetSampleInfo(sample);
    TChain* const chain = dy::GetSampleTChain(sample_info.sample);
    num_events = (chain->GetEntries() < num_events ? chain->GetEntries() : num_events);
    std::cout << "Beginning " << sample_info.name << " on " << num_events << " of " << chain->GetEntries() << " events..." << std::endl; 

     // scale the # of events
    const float nevts_full  = chain->GetMaximum("evt_nEvts");
    const float nevts_file  = chain->GetEntries(); 
    const float nevts_scale = nevts_full/nevts_file;

    // lumi
    static const float lumi = 0.082; //fb^-1

    // efficiency scale
    float eff_scale = sample_info.filter_eff;

    // overall scale
    const float scale = lumi * eff_scale * nevts_scale; 
    std::cout << "scale factor = " << Form("%f * %f * %f = %f", lumi, eff_scale, nevts_scale, scale) << std::endl;

    // book hists
    hc.Add(new TH1F(Form("h_count_mm_%s", sample_info.name.c_str()), Form("Event ee count (%s)"                               , sample_info.title.c_str()), 3 , -0.5, 2.5));
    hc.Add(new TH1F(Form("h_count_ee_%s", sample_info.name.c_str()), Form("Event #mu#mu count (%s)"                           , sample_info.title.c_str()), 3 , -0.5, 2.5));
    hc.Add(new TH1F(Form("h_mee_%s"     , sample_info.name.c_str()), Form("Dimuon Electron (%s);m_{ee}(GeV);# Events Expected", sample_info.title.c_str()), 60,   60, 120));
    hc.Add(new TH1F(Form("h_mmm_%s"     , sample_info.name.c_str()), Form("Dimuon Mass (%s);m_{#mu#mu}(GeV);# Events Expected", sample_info.title.c_str()), 60,   60, 120));
    hc.Sumw2();

    // alias for dilepton mass
    chain->SetAlias("mass_mm_px", "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.Px() : 0.0)");
    chain->SetAlias("mass_mm_py", "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.Py() : 0.0)");
    chain->SetAlias("mass_mm_pz", "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.Pz() : 0.0)");
    chain->SetAlias("mass_mm_e" , "Sum$(abs(genps_id)==13 && genps_status==3 ? genps_p4.E()  : 0.0)");
    chain->SetAlias("mass_mm"   , "sqrt(mass_mm_e*mass_mm_e - mass_mm_px*mass_mm_px - mass_mm_py*mass_mm_py - mass_mm_pz*mass_mm_pz)");

    chain->SetAlias("mass_ee_px", "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.Px() : 0.0)");
    chain->SetAlias("mass_ee_py", "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.Py() : 0.0)");
    chain->SetAlias("mass_ee_pz", "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.Pz() : 0.0)");
    chain->SetAlias("mass_ee_e" , "Sum$(abs(genps_id)==11 && genps_status==3 ? genps_p4.E()  : 0.0)");
    chain->SetAlias("mass_ee"   , "sqrt(mass_ee_e*mass_ee_e - mass_ee_px*mass_ee_px - mass_ee_py*mass_ee_py - mass_ee_pz*mass_ee_pz)");

    // fill hist
    hc.SetDirectory(gDirectory);
    const TCut selection_ee = Form("%1.4f*evt_scale1fb*(Sum$(genps_status==3 && genps_id==11)>=1 && Sum$(genps_status==3 && genps_id==-11)>=1)", scale);
    chain->Draw(">>event_list_ee", selection_ee, "goff", num_events);
    TEventList * const event_list_ee = dynamic_cast<TEventList*>(gDirectory->Get("event_list_ee"));
    chain->SetEventList(event_list_ee);
    std::cout << "filling ee hists..." << std::endl;
    chain->Draw(Form("1>>h_count_ee_%s" , sample_info.name.c_str()), selection_ee, "goff", num_events);
    chain->Draw(Form("mass_ee>>h_mee_%s", sample_info.name.c_str()), selection_ee, "goff", num_events);

    const TCut selection_mm = Form("%1.4f*evt_scale1fb*(Sum$(genps_status==3 && genps_id==13)>=1 && Sum$(genps_status==3 && genps_id==-13)>=1)", scale);
    chain->SetEventList(NULL);
    chain->Draw(">>event_list_mm", selection_mm, "goff", num_events);
    TEventList * const event_list_mm = dynamic_cast<TEventList*>(gDirectory->Get("event_list_mm"));
    chain->SetEventList(event_list_mm);
    std::cout << "filling mm hists..." << std::endl;
    chain->Draw(Form("1>>h_count_mm_%s" , sample_info.name.c_str()), selection_mm, "goff", num_events);
    chain->Draw(Form("mass_mm>>h_mmm_%s", sample_info.name.c_str()), selection_mm, "goff", num_events);
    hc.SetDirectory(NULL);

    std::cout << "Complete " << sample_info.name  << ":  ";
    std::cout << "mm count = " << rt::Integral(hc["h_count_mm_"+sample_info.name]) << " (" << hc["h_count_mm_"+sample_info.name]->GetEntries() << ") : "; 
    std::cout << "ee count = " << rt::Integral(hc["h_count_ee_"+sample_info.name]) << " (" << hc["h_count_ee_"+sample_info.name]->GetEntries() << ")\n" << std::endl; 
   
    // done
    return;
}

void create_dilep_count(const long num_events = std::numeric_limits<long>::max())
{
    // get hists
    rt::TH1Container hc;
    FillHists(hc, dy::Sample::dyll    , num_events);
    FillHists(hc, dy::Sample::wjets   , num_events);
    FillHists(hc, dy::Sample::ttdil   , num_events);
    FillHists(hc, dy::Sample::ttslq   , num_events);
    FillHists(hc, dy::Sample::tthad   , num_events);
    FillHists(hc, dy::Sample::qcdmu15 , num_events);
    FillHists(hc, dy::Sample::ww2l2nu , num_events);
    FillHists(hc, dy::Sample::wz2l2q  , num_events);
    FillHists(hc, dy::Sample::wz3lnu  , num_events);
    FillHists(hc, dy::Sample::zz2l2q  , num_events);
    FillHists(hc, dy::Sample::zz2l2nu , num_events);
    FillHists(hc, dy::Sample::zz4l    , num_events);

    // write output
    hc.Write("plots/gen/dilep_counts.root");
}

void print_dilep_count()
{
    rt::TH1Container hc("plots/gen/dilep_counts.root");

    const float wz_mm = rt::Integral(hc["h_count_mm_wz2l2q" ]) + rt::Integral(hc["h_count_mm_wz3lnu"]);
    const float wz_ee = rt::Integral(hc["h_count_ee_wz2l2q" ]) + rt::Integral(hc["h_count_ee_wz3lnu"]);
    const float zz_mm = rt::Integral(hc["h_count_mm_zz2l2q" ]) + rt::Integral(hc["h_count_mm_zz2l2nu"]) + rt::Integral(hc["h_count_mm_zz4l"]);
    const float zz_ee = rt::Integral(hc["h_count_ee_zz2l2q" ]) + rt::Integral(hc["h_count_ee_zz2l2nu"]) + rt::Integral(hc["h_count_ee_zz4l"]);
    
    CTable t1;
    t1.useTitle();
    t1.setTitle("Generator level number of \\ell\\ell events (no selections)");
    t1.setTable() (                                                 "mm ",                                   "ee")
                  ( "DY --> mm"   , rt::Integral(hc["h_count_mm_dyll"   ]), rt::Integral(hc["h_count_ee_dyll"   ]))
                  ( "W --> mnu"   , rt::Integral(hc["h_count_mm_wjets"  ]), rt::Integral(hc["h_count_ee_wjets"  ]))
                  ( "tt"          , rt::Integral(hc["h_count_mm_ttdil"  ]), rt::Integral(hc["h_count_ee_ttdil"  ]))
                  ( "QCD"         , rt::Integral(hc["h_count_mm_qcdmu15"]), rt::Integral(hc["h_count_ee_qcdmu15"]))
                  ( "WW"          , rt::Integral(hc["h_count_mm_ww2l2nu"]), rt::Integral(hc["h_count_ee_ww2l2nu"]))
                  ( "WZ"          , wz_mm                                 , wz_ee                                 )
                  ( "WZ --> 2l2q" , rt::Integral(hc["h_count_mm_wz2l2q" ]), rt::Integral(hc["h_count_ee_wz2l2q" ]))
                  ( "WZ --> 3lnu" , rt::Integral(hc["h_count_mm_wz3lnu" ]), rt::Integral(hc["h_count_ee_wz3lnu" ]))
                  ( "ZZ"          , zz_mm                                 , zz_ee                                 ) 
                  ( "ZZ --> 2l2q" , rt::Integral(hc["h_count_mm_zz2l2q" ]), rt::Integral(hc["h_count_ee_zz2l2q" ]))
                  ( "ZZ --> 2lnu" , rt::Integral(hc["h_count_mm_zz2l2nu"]), rt::Integral(hc["h_count_ee_zz2l2nu"]))
                  ( "ZZ --> 4l"   , rt::Integral(hc["h_count_mm_zz4l"   ]), rt::Integral(hc["h_count_ee_zz4l"   ]))
                  ;
    std::cout << t1 << std::endl;
}

void print_dilep_plots(const std::string& suffix = "png")
{
    rt::TH1Container hc("plots/gen/dilep_counts.root");

    // formatting
    rt::TH1Overlay::legend_height_per_entry_default = 0.050;
    rt::TH1Overlay::legend_text_size_default        = 0.025;
    const float ymin = 1e0;
    const float ymax = 3e4;
    const std::string option = "sb::off dt::stack lg::top_left";

    rt::TH1Overlay p1("Generator Level Dimuon distributions - #intL = 82 pb^{-1};m_{#mu^{+}#mu^{-}} (GeV/c^{2});Events/1.0 GeV/c^{2}", option);
    p1.Add(hc["h_mmm_dyll"   ],  "DY #rightarrow #mu#mu"        , kOrange-2);
    p1.Add(hc["h_mmm_wz2l2q" ],  "WZ #rightarrow 2l2q"          , kRed+2   );
    p1.Add(hc["h_mmm_wz3lnu" ],  "WZ #rightarrow 3l#nu"         , kRed-7   );
    p1.Add(hc["h_mmm_zz2l2q" ],  "ZZ #rightarrow 2l2q"          , kBlue    );
    p1.Add(hc["h_mmm_zz2l2nu"],  "ZZ #rightarrow 2l#nu"         , kBlue+2  );
    p1.Add(hc["h_mmm_zz4l"   ],  "ZZ #rightarrow 4l"            , kBlue-5  );
    p1.Add(hc["h_mmm_ttdil"  ],  "t#bar{t} #rightarrow 2l2#nu2b", kOrange+7);
    p1.Add(hc["h_mmm_ww2l2nu"],  "WW #rightarrow 2l2#nu"        , kGreen+2 );
    p1.Add(hc["h_mmm_wjets"  ],  "W #rightarrow l#nu"           , kGray+1  );
    p1.Add(hc["h_mmm_qcdmu15"],  "QCD"                          , kMagenta );
    p1.SetYAxisRange(ymin, ymax);

    rt::TH1Overlay p2("Generator Level Dielectron distributions - #intL = 82 pb^{-1};m_{e^{+}e^{-}} (GeV/c^{2});Events/1.0 GeV/c^{2}", option);
    p2.Add(hc["h_mee_dyll"   ],  "DY #rightarrow ee"            , kOrange-2);
    p2.Add(hc["h_mee_wz2l2q" ],  "WZ #rightarrow 2l2q"          , kRed+2   );
    p2.Add(hc["h_mee_wz3lnu" ],  "WZ #rightarrow 3l#nu"         , kRed-7   );
    p2.Add(hc["h_mee_zz2l2q" ],  "ZZ #rightarrow 2l2q"          , kBlue    );
    p2.Add(hc["h_mee_zz2l2nu"],  "ZZ #rightarrow 2l#nu"         , kBlue+2  );
    p2.Add(hc["h_mee_zz4l"   ],  "ZZ #rightarrow 4l"            , kBlue-5  );
    p2.Add(hc["h_mee_ttdil"  ],  "t#bar{t} #rightarrow 2l2#nu2b", kOrange+7);
    p2.Add(hc["h_mee_ww2l2nu"],  "WW #rightarrow 2l2#nu"        , kGreen+2 );
    p2.Add(hc["h_mee_wjets"  ],  "W #rightarrow l#nu"           , kGray+1  );
    p2.Add(hc["h_mee_qcdmu15"],  "QCD"                          , kMagenta );
    p2.SetYAxisRange(ymin, ymax);

    // output
    rt::Print(p1, "plots/gen/p_gen_mmm", suffix);
    rt::Print(p2, "plots/gen/p_gen_mee", suffix);

    // log output
    p1.SetLogy(true);
    rt::Print(p1, "plots/gen/p_gen_mmm_log", suffix);
    p2.SetLogy(true);
    rt::Print(p2, "plots/gen/p_gen_mee_log", suffix);
}

void gen_dilep_count(const long num_events = std::numeric_limits<long>::max())
{
/*     create_dilep_count(num_events); */
    print_dilep_count();
    print_dilep_plots("eps");
}
