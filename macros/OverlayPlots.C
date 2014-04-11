#include <vector>
#include <string>
#include <map>
#include "TStyle.h" 
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"

rt::TH1Overlay CreateOverlay
(
    const std::map<dy::Sample::Info, rt::TH1Container>& sample_hist_map,
    const std::string& hist_stem, 
    const std::string& title, 
    std::string option = "sb::off dt::stack lg::top_right",
    const double ymin = 1.0,
    const double ymax = -1.0,
    const int ncol = 1
)
{
    // colors
    const Style_t data_marker = 20;
    const Width_t data_width  = 2;
    const float marker_size   = 0.9;
/*     const Style_t shade_style = 3335; */
/*     const string unc_legend   = "Total Uncertainty"; */
    const std::string hist_name = "h_" + hist_stem;

    if (lt::string_contains(hist_stem, "yield"))
    {
        option.append(" B");
        const float width = 0.6;
        const float offset = 0.2;
        const float label_size = 0.07;
        for (const auto& hc : sample_hist_map)
        {
            TH1* const hist = hc.second[hist_name];
            hist->SetBarWidth(width);
            hist->SetBarOffset(offset);
            hist->SetLabelSize(label_size);
            hist->GetXaxis()->SetLabelSize(0.06);
            hist->GetXaxis()->SetBinLabel(1, ""      );
            hist->GetXaxis()->SetBinLabel(2, "ll"    );
            hist->GetXaxis()->SetBinLabel(3, "#mu#mu");
            hist->GetXaxis()->SetBinLabel(4, "ee"    );
        }
    }

    // uncertainy histogram
    TH1* const h_pred = static_cast<TH1*>(sample_hist_map.at(dy::GetSampleInfo(dy::Sample::data))[hist_name]->Clone());
    h_pred->Reset();
    h_pred->SetLineColor(kWhite);
    h_pred->SetFillColor(kBlack);
    h_pred->SetFillStyle(3335);
    h_pred->SetDrawOption("hist E5");
    for (const auto& hc : sample_hist_map)
    {
        if (hc.first.sample == dy::Sample::data) {continue;}
        TH1* const hist = hc.second[hist_name];
        h_pred->Add(hist);
    }

    // create the overlay
    rt::TH1Overlay::legend_height_per_entry_default = 0.040;
    rt::TH1Overlay::legend_text_size_default        = 0.025;
    rt::TH1Overlay::legend_ncol_default             = ncol;
    rt::TH1Overlay p(title, option);
    for (const auto& hc : sample_hist_map)
    {
        if (hc.first.sample == dy::Sample::data && lt::string_contains(hist_stem, "gen")) {continue;}

        TH1* const hist          = hc.second[hist_name];
        const std::string legend = hc.first.title;
        const Color_t color      = hc.first.color;
        if (hc.first.sample == dy::Sample::data)
        {
            hist->SetMarkerSize(marker_size);
            p.Add(hist, /*no_stack=*/true, legend, color, data_width, data_marker);
        }
        else
        {
            p.Add(hist, legend, color);
        }
    }
    if (ymin < ymax)
    {
        p.SetYAxisRange(ymin, ymax);
    }
    //p.Add(h_pred, /*no_stack=*/true, unc_legend, 1, 2, 1, shade_style);
    return p;
}

rt::TH1Container GetSampleHists(dy::Sample::Info sample_info, const std::string& label)
{
    rt::TH1Container hc(Form("plots/%s/%s_plots.root", label.c_str(), sample_info.name.c_str()));
    return hc;
}

// Overlay the kinematic plots
void OverlayPlots
(
    const std::string& label, 
    const std::string& suffix = "png"
)
{
    // get the histograms
    dy::SampleMap sm = dy::GetSampleMap();
    std::map<dy::Sample::Info, rt::TH1Container> sample_hist_map;
    for (const auto& s : sm)
    {
        sample_hist_map[s.second] = GetSampleHists(s.second, label);
    }
    // set style
    rt::SetTDRStyle();
    gStyle->SetHatchesSpacing(1.00);

    // title
    const double lumi = 0.082;
    const std::string title = Form("CMS2 Drell-Yan Exercise, #sqrt{s} = 8 TeV, L_{int} = %1.3f fb^{-1}", lumi);

    // overlays
    map<string, rt::TH1Overlay> p;
    rt::TH1Overlay::profile_marker_size_default = 10.0;
    p["p_gen_yield"         ] = CreateOverlay(sample_hist_map, "gen_yield"       , Form("%s;channel;Generator Events", title.c_str()), "sb::off dt::stack lg::top_left logy", 10, 3e6);
    p["p_reco_yield"        ] = CreateOverlay(sample_hist_map, "reco_yield"      , Form("%s;channel;Selected Events" , title.c_str()), "sb::off dt::stack lg::top_left logy", 10, 3e6);
    p["p_reco_nosel_yield"  ] = CreateOverlay(sample_hist_map, "reco_nosel_yield", Form("%s;channel;Selected Events" , title.c_str()), "sb::off dt::stack lg::top_left logy", 10, 3e6);
    p["p_gen_mmm"           ] = CreateOverlay(sample_hist_map, "gen_mmm"         , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_gen_mmm_log"       ] = CreateOverlay(sample_hist_map, "gen_mmm"         , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left logy", 5, 3e5, 1);
    p["p_reco_mmm"          ] = CreateOverlay(sample_hist_map, "reco_mmm"        , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_reco_mmm_log"      ] = CreateOverlay(sample_hist_map, "reco_mmm"        , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left logy", 5, 3e5, 1);
    p["p_reco_nosel_mmm"    ] = CreateOverlay(sample_hist_map, "reco_nosel_mmm"  , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_reco_nosel_mmm_log"] = CreateOverlay(sample_hist_map, "reco_nosel_mmm"  , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left logy", 5, 3e6, 1);
    p["p_gen_mee"           ] = CreateOverlay(sample_hist_map, "gen_mee"         , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_gen_mee_log"       ] = CreateOverlay(sample_hist_map, "gen_mee"         , Form("%s;m_{#mu#mu} (GeV);Events" , title.c_str()), "sb::off dt::stack lg::top_left logy", 5, 3e5, 1);
    p["p_reco_mee"          ] = CreateOverlay(sample_hist_map, "reco_mee"        , Form("%s;m_{ee} (GeV);Events"     , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_reco_mee_log"      ] = CreateOverlay(sample_hist_map, "reco_mee"        , Form("%s;m_{ee} (GeV);Events"     , title.c_str()), "sb::off dt::stack lg::top_left logy", 5, 3e5, 1);
    p["p_reco_nosel_mee"    ] = CreateOverlay(sample_hist_map, "reco_nosel_mee"  , Form("%s;m_{ee} (GeV);Events"     , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_reco_nosel_mee_log"] = CreateOverlay(sample_hist_map, "reco_nosel_mee"  , Form("%s;m_{ee} (GeV);Events"     , title.c_str()), "sb::off dt::stack lg::top_left logy", 10, 3e6, 1);

     // write
    const string plots_path = Form("plots/%s/overlays", label.c_str());
    rt::Print(p, plots_path, suffix);

    // print yield explicitly
    // this is a kludge to the the x error bars the right size for the yeild plot
    gStyle->SetErrorX(0.3);
    rt::Print(p["p_reco_yield"      ], Form("%s/p_reco_yield"      , plots_path.c_str()), suffix);
    rt::Print(p["p_reco_nosel_yield"], Form("%s/p_reco_nosel_yield", plots_path.c_str()), suffix);
    gStyle->SetErrorX();
}
