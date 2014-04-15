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
    hc.Add(new TH1D("h_gen_yield" , "Yield count of gen level l^{+}l^{-}"          ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_mee"   , "Generator level dielectron mass;m_{ee} (GeV)" , 150, 0, 150));
    hc.Add(new TH1D("h_gen_mmm"   , "Generator level dilmuon mass;m_{#mu#mu} (GeV)", 150, 0, 150));
    hc.Add(new TH1D("h_gen_mll"   , "Generator level dilepton mass;m_{ll} (GeV)"   , 150, 0, 150));

    hc.Add(new TH1D("h_gen_notau_yield", "Yield count of gen level l^{+}l^{-} (no #taus)"          ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_notau_mee"  , "Generator level dielectron mass (no #taus);m_{ee} (GeV)" , 150, 0, 150));
    hc.Add(new TH1D("h_gen_notau_mmm"  , "Generator level dilmuon mass (no #taus);m_{#mu#mu} (GeV)", 150, 0, 150));
    hc.Add(new TH1D("h_gen_notau_mll"  , "Generator level dilepton mass (no #taus);m_{ll} (GeV)"   , 150, 0, 150));

    hc.Add(new TH1D("h_gen_tau_yield", "Yield count of gen level l^{+}l^{-} (#taus #rightarrow e/#mu)"          ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_tau_mee"  , "Generator level dielectron mass (#taus #rightarrow e/#mu);m_{ee} (GeV)" , 150, 0, 150));
    hc.Add(new TH1D("h_gen_tau_mmm"  , "Generator level dilmuon mass (#taus #rightarrow e/#mu);m_{#mu#mu} (GeV)", 150, 0, 150));
    hc.Add(new TH1D("h_gen_tau_mll"  , "Generator level dilepton mass (#taus #rightarrow e/#mu);m_{ll} (GeV)"   , 150, 0, 150));

    // acceptance plots
    hc.Add(new TH1D("h_acc_gen_den", "Acceptence denominator;Channel;Event Count"        , 4, -1, 3));
    hc.Add(new TH1D("h_acc_gen_num", "Acceptence generator numerator;Channel;Event Count", 4, -1, 3));
    hc.Add(new TH1D("h_acc_rec_num", "Acceptence numerator;Channel;Event Count"          , 4, -1, 3));

/*     hc.Add(new TH1D("h_reco_ossf_yield", "OSSF yield count of reco level l^{+}l^{-}",   4, -1,  3)); */
/*     hc.Add(new TH1D("h_reco_ossf_mee"  , "OSSF dielectron mass;m_{ee} (GeV)"        , 150, 0, 150)); */
/*     hc.Add(new TH1D("h_reco_ossf_mmm"  , "OSSF dilmuon mass;m_{#mu#mu} (GeV)"       , 150, 0, 150)); */
/*     hc.Add(new TH1D("h_reco_ossf_mll"  , "OSSF dilepton mass;m_{ll} (GeV)"          , 150, 0, 150)); */

/*     hc.Add(new TH1D("h_reco_full_yield", "Full yield count of reco level l^{+}l^{-}",   4, -1,  3)); */
/*     hc.Add(new TH1D("h_reco_full_mee"  , "Full dielectron mass;m_{ee} (GeV)"        , 150, 0, 150)); */
/*     hc.Add(new TH1D("h_reco_full_mmm"  , "Full dilmuon mass;m_{#mu#mu} (GeV)"       , 150, 0, 150)); */
/*     hc.Add(new TH1D("h_reco_full_mll"  , "Full dilepton mass;m_{ll} (GeV)"          , 150, 0, 150)); */
    // reco level plots
    /*     hc.Add(new TH1D("h_reco_cutflow_ee", "e^{+}e^{-}  Cut Flow;Selection;Event Count"    , Selection::static_size, 0, Selection::static_size)); */
    /*     hc.Add(new TH1D("h_reco_cutflow_mm", "#mu^{+}#mu^{-}  Cut Flow;Selection;Event Count", Selection::static_size, 0, Selection::static_size)); */
    /*     hc.Add(new TH1D("h_reco_cutflow_ll", "l^{+}l^{-} Cut Flow;Selection;Event Count"     , Selection::static_size, 0, Selection::static_size)); */
    /*     SetCutflowAxisLabel(hc["h_reco_cutflow_ee"]); */
    /*     SetCutflowAxisLabel(hc["h_reco_cutflow_mm"]); */
    /*     SetCutflowAxisLabel(hc["h_reco_cutflow_ll"]); */
    /*     for (const auto& sel : Selections()) */
    /*     { */
    /*         hc.Add(new TH1D(Form("h_reco_%s_yield", sel.name.c_str()), Form("%s yield count of reco level l^{+}l^{-}", sel.title.c_str()),   4, -1,  3)); */
    /*         hc.Add(new TH1D(Form("h_reco_%s_mee"  , sel.name.c_str()), Form("%s dielectron mass;m_{ee} (GeV)"        , sel.title.c_str()), 150, 0, 150)); */
    /*         hc.Add(new TH1D(Form("h_reco_%s_mmm"  , sel.name.c_str()), Form("%s dilmuon mass;m_{#mu#mu} (GeV)"       , sel.title.c_str()), 150, 0, 150)); */
    /*         hc.Add(new TH1D(Form("h_reco_%s_mll"  , sel.name.c_str()), Form("%s dilepton mass;m_{ll} (GeV)"          , sel.title.c_str()), 150, 0, 150)); */
    /*         SetYieldAxisLabel(hc["h_reco_"+sel.name+"_yield"]); */
    /*     } */

    // change axis labels
    SetYieldAxisLabel(hc["h_gen_yield"      ]);
    SetYieldAxisLabel(hc["h_gen_notau_yield"]);
    SetYieldAxisLabel(hc["h_gen_tau_yield"  ]);
    SetYieldAxisLabel(hc["h_acc_gen_den"    ]);
    SetYieldAxisLabel(hc["h_acc_gen_num"    ]);
    SetYieldAxisLabel(hc["h_acc_rec_num"    ]);

    // sumw2()
    hc.Sumw2();
}

TCut WrapCut(const TCut& cut)
{
    const TCut result = Form("lumi*scale1fb*(%s)", cut.GetTitle());
    return result;
}

void CreatePlots
(
    const dy::Sample::Info& sample_info,
    const std::string& label
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
    const TCut is_gen_ll            = is_gen_ee || is_gen_mm; 
    const TCut is_gen_ee_includetau = "is_gen_ee_includetau";
    const TCut is_gen_mm_includetau = "is_gen_mm_includetau";
    const TCut is_gen_ll_includetau = is_gen_ee_includetau || is_gen_mm_includetau;

    // Fill hists
    // ----------------------------- // 

    hc.SetDirectory(gDirectory);

    // raw gen yields
    chain.Draw("0>>h_gen_yield" , WrapCut(is_gen_ll), "goff");
    chain.Draw("1>>+h_gen_yield", WrapCut(is_gen_ee), "goff");
    chain.Draw("2>>+h_gen_yield", WrapCut(is_gen_mm), "goff");

    chain.Draw("gen_p4.mass()>>h_gen_mee", WrapCut(is_gen_ee), "goff");
    chain.Draw("gen_p4.mass()>>h_gen_mmm", WrapCut(is_gen_mm), "goff");
    chain.Draw("gen_p4.mass()>>h_gen_mll", WrapCut(is_gen_ll), "goff");

    // gen yields including taus
    chain.Draw("0>>h_gen_tau_yield" , WrapCut(is_gen_ll_includetau), "goff");
    chain.Draw("1>>+h_gen_tau_yield", WrapCut(is_gen_ee_includetau), "goff");
    chain.Draw("2>>+h_gen_tau_yield", WrapCut(is_gen_mm_includetau), "goff");

    chain.Draw("gen_p4.mass()>>h_gen_tau_mee", WrapCut(is_gen_ee_includetau), "goff");
    chain.Draw("gen_p4.mass()>>h_gen_tau_mmm", WrapCut(is_gen_mm_includetau), "goff");
    chain.Draw("gen_p4.mass()>>h_gen_tau_mll", WrapCut(is_gen_ll_includetau), "goff");
    hc.SetDirectory(NULL);

    // write plots
    // ----------------------------- // 

    hc.Write(Form("plots/%s/%s.root", label.c_str(), sample_info.name.c_str()));
}

void CreatePlotsFromBabies(const std::string& label)
{
    const auto& sample_infos = dy::Sample::GetInfos();
    for (const auto& sample_info : sample_infos)
    {
        std::cout << "[CreatePlotsFromBabies] create histograms for " << sample_info.name << "\n";
        CreatePlots(sample_info, label);
    }
}
