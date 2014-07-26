#include <vector>
#include <string>
#include <map>
#include "THStack.h"
#include "TStyle.h" 
#include "TLegend.h" 
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"

rt::TH1Container GetSampleHists(dy::Sample::Info sample_info, const std::string& label)
{
    rt::TH1Container hc(Form("plots/%s/%s_plots.root", label.c_str(), sample_info.name.c_str()));
    return hc;
}

void test_stack()
{
    // get the histograms
    std::map<dy::Sample::Info, rt::TH1Container> sample_hist_map;
    for (const auto& s : dy::Sample::GetInfos())
    {
        sample_hist_map[s] = GetSampleHists(s, "full_babies");
        sample_hist_map[s].SetLineColor(kBlack);
        sample_hist_map[s].SetFillColor(s.color);
        sample_hist_map[s]["h_gen_mee"]->SetMinimum(1e-1);
        sample_hist_map[s]["h_gen_mee"]->SetMaximum(1e5);
    }

    // build the statck
    THStack* const h_stack = new THStack;
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::dytt   )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::wjets  )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::ttdil  )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::ttslq  )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::tthad  )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::qcdmu15)]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::ww2l2nu)]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::wz2l2q )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::wz3lnu )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::zz2l2nu)]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::zz2l2q )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::zz4l   )]["h_gen_mee"]->Clone()), "hist");
    h_stack->Add(dynamic_cast<TH1*>(sample_hist_map[dy::GetSampleInfo(dy::Sample::dyll   )]["h_gen_mee"]->Clone()), "hist");
    h_stack->SetMinimum(1e-1);
    h_stack->SetMaximum(1e5);
    h_stack->SetTitle(";m_{ee} (GeV);Events/1.0 GeV");
    
    // legend
    TLegend* const legend = new TLegend(0.15, 0.4, 0.4, 0.9);
    legend->SetFillColor(0);  // 0 makes it the background clear on the pad
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    for (const auto& s : dy::Sample::GetInfos())
    {
        if (s.sample == dy::Sample::data) continue;
        legend->AddEntry(sample_hist_map[s]["h_gen_mee"]->Clone(), s.title.c_str(), "f"); 
    }

    TCanvas* const canvas = new TCanvas("canvas", "cavnas", 600, 600);
    canvas->SetLogy();
    h_stack->Draw();
    legend->Draw();
}


