#include <fstream>
#include <sstream>
#include <iostream>
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "TChain.h"
#include "TCut.h"
#include "TDirectory.h"

void SetYieldAxisLabel(TH1* const hist)
{
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetBinLabel(1, ""      );
    hist->GetXaxis()->SetBinLabel(2, "ll"    );
    hist->GetXaxis()->SetBinLabel(3, "#mu#mu");
    hist->GetXaxis()->SetBinLabel(4, "ee"    );
}

void BookHists(rt::TH1Container& hc)
{
    // gen level plots
    hc.Add(new TH1D("h_gen_yield"          , "Yield count of gen level l^{+}l^{-}"                              ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_notau_yield"    , "Yield count of gen level l^{+}l^{-} (no #taus)"                   ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_tau_yield"      , "Yield count of gen level l^{+}l^{-} (#taus #rightarrow e/#mu)"    ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_raw_yield"      , "Yield raw count of gen level l^{+}l^{-} (#taus #rightarrow e/#mu)",   4, -1,  3));
    hc.Add(new TH1D("h_gen_raw_notau_yield", "Yield raw count of gen level l^{+}l^{-} (no #taus)"               ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_raw_tau_yield"  , "Yield raw count of gen level l^{+}l^{-} (#taus #rightarrow e/#mu)",   4, -1,  3));

    hc.Add(new TH1D("h_gen_mee", "Generator level dielectron mass;m_{ee} (GeV)" , 150, 0, 150));
    hc.Add(new TH1D("h_gen_mmm", "Generator level dilmuon mass;m_{#mu#mu} (GeV)", 150, 0, 150));
    hc.Add(new TH1D("h_gen_mll", "Generator level dilepton mass;m_{ll} (GeV)"   , 150, 0, 150));

    // acceptance plots
    hc.Add(new TH1D("h_acc_gen_den", "Acceptence denominator;Channel;Event Count"        , 4, -1, 3));
    hc.Add(new TH1D("h_acc_gen_num", "Acceptence generator numerator;Channel;Event Count", 4, -1, 3));
    hc.Add(new TH1D("h_acc_rec_num", "Acceptence numerator;Channel;Event Count"          , 4, -1, 3));

    // reco level plots
    hc.Add(new TH1D("h_reco_ossf_yield" , "OSSF yield count of reco level l^{+}l^{-}"  ,  4, -1,  3));
    hc.Add(new TH1D("h_reco_mwin_yield" , "mwin yield count of reco level l^{+}l^{-}"  ,  4, -1,  3));
    hc.Add(new TH1D("h_reco_svtx_yield" , "svtx yield count of reco level l^{+}l^{-}"  ,  4, -1,  3));
    hc.Add(new TH1D("h_reco_trig_yield" , "trig yield count of reco level l^{+}l^{-}"  ,  4, -1,  3));
    hc.Add(new TH1D("h_reco_idiso_yield", "Id/ISO yield count of reco level l^{+}l^{-}",  4, -1,  3));
    hc.Add(new TH1D("h_reco_full_yield" , "Full yield count of reco level l^{+}l^{-}"  ,  4, -1,  3));

    hc.Add(new TH1D("h_reco_full_mee"  , "Full dielectron mass;m_{ee} (GeV)" , 80, 40, 120));
    hc.Add(new TH1D("h_reco_full_mmm"  , "Full dilmuon mass;m_{#mu#mu} (GeV)", 80, 40, 120));
    hc.Add(new TH1D("h_reco_full_mll"  , "Full dilepton mass;m_{ll} (GeV)"   , 80, 40, 120));

    hc.Add(new TH1D("h_reco_ossf_mee"  , "OSSF reco dielectron mass;m_{ee} (GeV)" , 80, 40, 120));
    hc.Add(new TH1D("h_reco_ossf_mmm"  , "OSSF reco dilmuon mass;m_{#mu#mu} (GeV)", 80, 40, 120));
    hc.Add(new TH1D("h_reco_ossf_mll"  , "OSSF reco dilepton mass;m_{ll} (GeV)"   , 80, 40, 120));

    // change axis labels
    SetYieldAxisLabel(hc["h_gen_yield"          ]);
    SetYieldAxisLabel(hc["h_gen_notau_yield"    ]);
    SetYieldAxisLabel(hc["h_gen_tau_yield"      ]);
    SetYieldAxisLabel(hc["h_gen_raw_yield"      ]);
    SetYieldAxisLabel(hc["h_gen_raw_notau_yield"]);
    SetYieldAxisLabel(hc["h_gen_raw_tau_yield"  ]);
    SetYieldAxisLabel(hc["h_reco_ossf_yield"    ]);
    SetYieldAxisLabel(hc["h_reco_mwin_yield"    ]);
    SetYieldAxisLabel(hc["h_reco_svtx_yield"    ]);
    SetYieldAxisLabel(hc["h_reco_trig_yield"    ]);
    SetYieldAxisLabel(hc["h_reco_idiso_yield"   ]);
    SetYieldAxisLabel(hc["h_reco_full_yield"    ]);
    SetYieldAxisLabel(hc["h_acc_gen_den"        ]);
    SetYieldAxisLabel(hc["h_acc_gen_num"        ]);
    SetYieldAxisLabel(hc["h_acc_rec_num"        ]);

    // sumw2()
    hc.Sumw2();
}

TCut ApplyScale(const TCut& cut, const bool do_sf = false)
{
    const double sf   = (do_sf ? 0.95*0.95 : 1.0);
    const TCut result = Form("lumi*scale1fb*%f*(%s)", sf, cut.GetTitle());
    return result;
}

void CreatePlots
(
    const dy::Sample::Info& sample_info,
    const std::string& label,
    const bool apply_sf = false
)
{
    // baby
    // ----------------------------- // 

    TChain chain("tree");
    chain.Add(Form("babies/%s_baby.root", sample_info.name.c_str()));

    // hists
    // ----------------------------- // 

    rt::TH1Container hc;
    BookHists(hc);

    // selections
    // ----------------------------- // 

    const TCut is_gen_ee            = "is_gen_ee";
    const TCut is_gen_mm            = "is_gen_mm";
    const TCut is_gen_tt            = "is_gen_tt";
    const TCut is_gen_ll            = is_gen_ee || is_gen_mm; 
    const TCut is_gen_ee_includetau = "is_gen_ee_includetau";
    const TCut is_gen_mm_includetau = "is_gen_mm_includetau";
    const TCut is_gen_ll_includetau = is_gen_ee_includetau || is_gen_mm_includetau;
    const TCut is_gen_ee_onlytau     = "is_gen_tt && is_gen_ee_includetau";
    const TCut is_gen_mm_onlytau     = "is_gen_tt && is_gen_mm_includetau";
    const TCut is_gen_ll_onlytau     = is_gen_ee_onlytau || is_gen_mm_onlytau;
    const TCut is_ee                 = "is_ee";
    const TCut is_mm                 = "is_mm";
    const TCut is_ll                 = "is_ll";
    const TCut gen_mwin              = "60 < gen_p4.mass() && gen_p4.mass() < 120";
    const TCut is_gen_z              = "is_gen_z";
    const TCut gen_acc_den           = "is_gen_acc_den";
    const TCut gen_acc_num           = "is_gen_acc_num";
    const TCut reco_acc_num          = gen_acc_den && ((is_gen_ee && is_ee) || (is_gen_mm && is_mm)) && "passes_full";

    // decide whether to apply the scale factor
    const float do_sf = (apply_sf && sample_info.sample != dy::Sample::data);

    // Fill hists
    // ----------------------------- // 

    hc.SetDirectory(gDirectory);

    // raw gen yields
    chain.Draw("0>> h_gen_yield", ApplyScale(is_gen_ll_includetau), "goff");
    chain.Draw("1>>+h_gen_yield", ApplyScale(is_gen_mm_includetau), "goff");
    chain.Draw("2>>+h_gen_yield", ApplyScale(is_gen_ee_includetau), "goff");

    chain.Draw("0>> h_gen_raw_yield", is_gen_ll_includetau, "goff");
    chain.Draw("1>>+h_gen_raw_yield", is_gen_mm_includetau, "goff");
    chain.Draw("2>>+h_gen_raw_yield", is_gen_ee_includetau, "goff");

    chain.Draw("gen_p4.mass()>>h_gen_mll", ApplyScale(is_gen_ll), "goff");
    chain.Draw("gen_p4.mass()>>h_gen_mmm", ApplyScale(is_gen_mm), "goff");
    chain.Draw("gen_p4.mass()>>h_gen_mee", ApplyScale(is_gen_ee), "goff");

    // acceptance
    chain.Draw("0>> h_acc_gen_den", ApplyScale(gen_acc_den && is_gen_ll), "goff");
    chain.Draw("1>>+h_acc_gen_den", ApplyScale(gen_acc_den && is_gen_mm), "goff");
    chain.Draw("2>>+h_acc_gen_den", ApplyScale(gen_acc_den && is_gen_ee), "goff");

    chain.Draw("0>> h_acc_gen_num", ApplyScale(gen_acc_num && is_gen_ll), "goff");
    chain.Draw("1>>+h_acc_gen_num", ApplyScale(gen_acc_num && is_gen_mm), "goff");
    chain.Draw("2>>+h_acc_gen_num", ApplyScale(gen_acc_num && is_gen_ee), "goff");

    chain.Draw("0>> h_acc_rec_num", ApplyScale(reco_acc_num && is_gen_ll), "goff");
    chain.Draw("1>>+h_acc_rec_num", ApplyScale(reco_acc_num && is_gen_mm), "goff");
    chain.Draw("2>>+h_acc_rec_num", ApplyScale(reco_acc_num && is_gen_ee), "goff");

    // gen yields not including taus
    chain.Draw("0>> h_gen_notau_yield", ApplyScale(is_gen_ll), "goff");
    chain.Draw("1>>+h_gen_notau_yield", ApplyScale(is_gen_mm), "goff");
    chain.Draw("2>>+h_gen_notau_yield", ApplyScale(is_gen_ee), "goff");
    chain.Draw("0>> h_gen_raw_notau_yield", is_gen_ll, "goff");
    chain.Draw("1>>+h_gen_raw_notau_yield", is_gen_mm, "goff");
    chain.Draw("2>>+h_gen_raw_notau_yield", is_gen_ee, "goff");

    // gen yields only including taus
    chain.Draw("0>>h_gen_tau_yield" , ApplyScale(is_gen_ll_onlytau), "goff");
    chain.Draw("1>>+h_gen_tau_yield", ApplyScale(is_gen_mm_onlytau), "goff");
    chain.Draw("2>>+h_gen_tau_yield", ApplyScale(is_gen_ee_onlytau), "goff");
    chain.Draw("0>> h_gen_raw_tau_yield", is_gen_ll_onlytau, "goff");
    chain.Draw("1>>+h_gen_raw_tau_yield", is_gen_mm_onlytau, "goff");
    chain.Draw("2>>+h_gen_raw_tau_yield", is_gen_ee_onlytau, "goff");

    // raw full yields
    chain.Draw("0>> h_reco_ossf_yield" , ApplyScale(is_ll && "passes_ossf", do_sf), "goff");
    chain.Draw("1>>+h_reco_ossf_yield" , ApplyScale(is_mm && "passes_ossf", do_sf), "goff");
    chain.Draw("2>>+h_reco_ossf_yield" , ApplyScale(is_ee && "passes_ossf", do_sf), "goff");
                                                                       
    chain.Draw("0>> h_reco_mwin_yield" , ApplyScale(is_ll && "passes_mwin", do_sf), "goff");
    chain.Draw("1>>+h_reco_mwin_yield" , ApplyScale(is_mm && "passes_mwin", do_sf), "goff");
    chain.Draw("2>>+h_reco_mwin_yield" , ApplyScale(is_ee && "passes_mwin", do_sf), "goff");
                                                                       
    chain.Draw("0>> h_reco_svtx_yield" , ApplyScale(is_ll && "passes_svtx", do_sf), "goff");
    chain.Draw("1>>+h_reco_svtx_yield" , ApplyScale(is_mm && "passes_svtx", do_sf), "goff");
    chain.Draw("2>>+h_reco_svtx_yield" , ApplyScale(is_ee && "passes_svtx", do_sf), "goff");
                                                                       
    chain.Draw("0>> h_reco_trig_yield" , ApplyScale(is_ll && "passes_trig", do_sf), "goff");
    chain.Draw("1>>+h_reco_trig_yield" , ApplyScale(is_mm && "passes_trig", do_sf), "goff");
    chain.Draw("2>>+h_reco_trig_yield" , ApplyScale(is_ee && "passes_trig", do_sf), "goff");
                                                                        
    chain.Draw("0>> h_reco_idiso_yield", ApplyScale(is_ll && "passes_idiso", do_sf), "goff");
    chain.Draw("1>>+h_reco_idiso_yield", ApplyScale(is_mm && "passes_idiso", do_sf), "goff");
    chain.Draw("2>>+h_reco_idiso_yield", ApplyScale(is_ee && "passes_idiso", do_sf), "goff");

    chain.Draw("0>> h_reco_full_yield" , ApplyScale(is_ll && "passes_full", do_sf), "goff");
    chain.Draw("1>>+h_reco_full_yield" , ApplyScale(is_mm && "passes_full", do_sf), "goff");
    chain.Draw("2>>+h_reco_full_yield" , ApplyScale(is_ee && "passes_full", do_sf), "goff");

    chain.Draw("hyp_p4.mass()>>h_reco_full_mee", ApplyScale(is_ll && "passes_full", do_sf), "goff");
    chain.Draw("hyp_p4.mass()>>h_reco_full_mmm", ApplyScale(is_mm && "passes_full", do_sf), "goff");
    chain.Draw("hyp_p4.mass()>>h_reco_full_mll", ApplyScale(is_ee && "passes_full", do_sf), "goff");

    chain.Draw("hyp_p4.mass()>>h_reco_ossf_mee", ApplyScale(is_ll && "passes_ossf", do_sf), "goff");
    chain.Draw("hyp_p4.mass()>>h_reco_ossf_mmm", ApplyScale(is_mm && "passes_ossf", do_sf), "goff");
    chain.Draw("hyp_p4.mass()>>h_reco_ossf_mll", ApplyScale(is_ee && "passes_ossf", do_sf), "goff");

    hc.SetDirectory(NULL);

    // post processing
    // ----------------------------- // 

    hc.Add(rt::DivideHists(hc["h_acc_rec_num"], hc["h_acc_gen_den"], "h_acc_rec", "Reco Acceptence;Channel;Event Count"));
    hc.Add(rt::DivideHists(hc["h_acc_gen_num"], hc["h_acc_gen_den"], "h_acc_gen", "Gen Acceptence;Channel;Event Count" ));
    SetYieldAxisLabel(hc["h_acc_rec"]);
    SetYieldAxisLabel(hc["h_acc_gen"]);

    // write plots
    // ----------------------------- // 

    hc.Write(Form("plots/%s/%s_plots.root", label.c_str(), sample_info.name.c_str()));
}

void CreatePlotsFromBabies(const std::string& label, const bool apply_sf = false)
{
    const auto& sample_infos = dy::Sample::GetInfos();
    for (const auto& sample_info : sample_infos)
    {
        std::cout << "[CreatePlotsFromBabies] create histograms for " << sample_info.name << "\n";
        CreatePlots(sample_info, label, apply_sf);
    }
}
