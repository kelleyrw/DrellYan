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
    std::string option = "sb::off dt::stack lg::top_right"
)
{
    // colors
    const Style_t data_marker = 20;
    const Width_t data_width  = 2;
    const float marker_size   = 0.9;
/*     const Style_t shade_style = 3335; */
/*     const string unc_legend   = "Total Uncertainty"; */
    const std::string hist_name = "h_" + hist_stem;

    if (hist_stem == "reco_yield" or hist_stem == "gen_yield")
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
            hist->GetXaxis()->SetRangeUser(-2, 4);
        }
    }

    // create the overlay
    rt::TH1Overlay::legend_height_per_entry_default = 0.040;
    rt::TH1Overlay::legend_text_size_default        = 0.025;
    rt::TH1Overlay::legend_ncol_default             = 1;
    rt::TH1Overlay p(title, option);
    for (const auto& hc : sample_hist_map)
    {
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
    //p.SetYAxisRange(0.0, max);
    if (p.GetLogy())
    {
        p.SetYAxisRange(1e-1, 1e5);
    }
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
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::data   )] = GetSampleHists(dy::Sample::data    , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::dyll   )] = GetSampleHists(dy::Sample::dyll    , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::wjets  )] = GetSampleHists(dy::Sample::wjets   , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::ttdil  )] = GetSampleHists(dy::Sample::ttdil   , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::ttslq  )] = GetSampleHists(dy::Sample::ttslq   , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::tthad  )] = GetSampleHists(dy::Sample::tthad   , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::qcdmu15)] = GetSampleHists(dy::Sample::qcdmu15 , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::ww2l2nu)] = GetSampleHists(dy::Sample::ww2l2nu , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::wz2l2q )] = GetSampleHists(dy::Sample::wz2l2q  , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::wz3lnu )] = GetSampleHists(dy::Sample::wz3lnu  , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::zz2l2nu)] = GetSampleHists(dy::Sample::zz2l2nu , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::zz2l2q )] = GetSampleHists(dy::Sample::zz2l2q  , label); */
/*     sample_hist_map[dy::GetSampleInfo(dy::Sample::zz4l   )] = GetSampleHists(dy::Sample::zz4l    , label); */
    
    // set style
    rt::SetTDRStyle();
    gStyle->SetHatchesSpacing(1.00);

    // title
    const double lumi = 0.082;
    std::string title = Form("CMS2 Drell-Yan Exercise, #sqrt{s} = 8 TeV, L_{int} = %1.3f fb^{-1}", lumi);

    // overlays
    map<string, rt::TH1Overlay> p;
    rt::TH1Overlay::profile_marker_size_default = 10.0;
    p["p_gen_yield"   ] = CreateOverlay(sample_hist_map, "gen_yield"  , Form("%s;channel;Generator Events"                , title.c_str()), "sb::off dt::stack lg::top_left" );
    p["p_reco_yield"  ] = CreateOverlay(sample_hist_map, "reco_yield" , Form("%s;channel;Selected Events"                 , title.c_str()), "sb::off dt::stack lg::top_left" );
    p["p_reco_mmm"    ] = CreateOverlay(sample_hist_map, "reco_mmm"   , Form("%s;m_{#mu#mu} (GeV);Events"             , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_reco_mmm_log"] = CreateOverlay(sample_hist_map, "reco_mmm"   , Form("%s;m_{#mu#mu} (GeV);Events"             , title.c_str()), "sb::off dt::stack lg::top_left logy");
    p["p_reco_nosel_mmm"    ] = CreateOverlay(sample_hist_map, "reco_nosel_mmm"   , Form("%s;m_{#mu#mu} (GeV);Events"             , title.c_str()), "sb::off dt::stack lg::top_left");
    p["p_reco_nosel_mmm_log"] = CreateOverlay(sample_hist_map, "reco_nosel_mmm"   , Form("%s;m_{#mu#mu} (GeV);Events"             , title.c_str()), "sb::off dt::stack lg::top_left logy");
/*     p["p_pt1"         ] = CreateOverlay(sample_hist_map, "pt1"        , Form("%s;p^{lep1}_{T} (GeV);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_pt2"         ] = CreateOverlay(sample_hist_map, "pt2"        , Form("%s;p^{lep2}_{T} (GeV);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_met"         ] = CreateOverlay(sample_hist_map, "met"        , Form("%s;E^{miss}_{T} (GeV);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_ht"          ] = CreateOverlay(sample_hist_map, "ht"         , Form("%s;H_{T} (GeV);Events"              , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_mt"          ] = CreateOverlay(sample_hist_map, "mt"         , Form("%s;m_{T} (GeV);Events"              , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_nbtags"      ] = CreateOverlay(sample_hist_map, "nbtags"     , Form("%s;number of b-tagged jets;Events"  , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_njets"       ] = CreateOverlay(sample_hist_map, "njets"      , Form("%s;number of jets;Events"           , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_lepdphi"     ] = CreateOverlay(sample_hist_map, "lepdphi"    , Form("%s;#Delta#Phi(lep1, lep2);Events"   , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_lepdeta"     ] = CreateOverlay(sample_hist_map, "lepdeta"    , Form("%s;#Delta#eta(lep1, lep2);Events"   , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_lepdr"       ] = CreateOverlay(sample_hist_map, "lepdr"      , Form("%s;#DeltaR(lep1, lep2);Events"      , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_drlepb"      ] = CreateOverlay(sample_hist_map, "drlepb"     , Form("%s;#DeltaR(lep, btag);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_btagdr"      ] = CreateOverlay(sample_hist_map, "btagdr"     , Form("%s;#DeltaR(btag1, btag2);Events"    , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_drjetb"      ] = CreateOverlay(sample_hist_map, "drjetb"     , Form("%s;#DeltaR(btag, jet);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_ptjetlep"    ] = CreateOverlay(sample_hist_map, "ptjetlep"   , Form("%s;jet p_{T} / lep p_{T} - 1;Events", title.c_str()), "sb::off dt::stack lg::right"); */
/*     p["p_drlep3rdlep" ] = CreateOverlay(sample_hist_map, "drlep3rdlep", Form("%s;#DeltaR(lep, 3rd lep);Events"    , title.c_str()), "sb::off dt::stack lg::right"); */

    // overlay individual channels
/*     for (size_t i = 1; i != at::DileptonHypType::static_size; i++) */
/*     { */
/*         at::DileptonHypType::value_type hyp_type = static_cast<at::DileptonHypType::value_type>(i); */
/*  */
/*         // name and title suffixes */
/*         string hn = Form("_%s" ,  GetDileptonHypTypeName(hyp_type).c_str()); */
/*         //string ht = Form(" (%s)",  GetDileptonHypTypeTitle(hyp_type).c_str()); */
/*  */
/*         p["p_dilep_mass"+hn] = CreateOverlay(hc_data, hc_mc, "dilep_mass"+hn , Form("%s;m_{ll} (GeV);Events"             , title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_pt1"       +hn] = CreateOverlay(hc_data, hc_mc, "pt1"       +hn , Form("%s;p^{lep1}_{T} (GeV);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_pt2"       +hn] = CreateOverlay(hc_data, hc_mc, "pt2"       +hn , Form("%s;p^{lep2}_{T} (GeV);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_met"       +hn] = CreateOverlay(hc_data, hc_mc, "met"       +hn , Form("%s;E^{miss}_{T} (GeV);Events"       , title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_ht"        +hn] = CreateOverlay(hc_data, hc_mc, "ht"        +hn , Form("%s;H_{T} (GeV);Events"              , title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_mt"        +hn] = CreateOverlay(hc_data, hc_mc, "mt"        +hn , Form("%s;m_{T} (GeV);Events"              , title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_nbtags"    +hn] = CreateOverlay(hc_data, hc_mc, "nbtags"    +hn , Form("%s;number of b-tagged jets;Events"  , title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_njets"     +hn] = CreateOverlay(hc_data, hc_mc, "njets"     +hn , Form("%s;number of jets;Events"           , title.c_str()), "sb::off dt::stack lg::right"); */
/*  */
/*         p["p_dilep_mass_nj0"+hn] = CreateOverlay(hc_data, hc_mc, "dilep_mass_nj0"+hn , Form("%s;m_{ll} (GeV);Events", title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_dilep_mass_nj1"+hn] = CreateOverlay(hc_data, hc_mc, "dilep_mass_nj1"+hn , Form("%s;m_{ll} (GeV);Events", title.c_str()), "sb::off dt::stack lg::right"); */
/*         p["p_dilep_mass_nj2"+hn] = CreateOverlay(hc_data, hc_mc, "dilep_mass_nj2"+hn , Form("%s;m_{ll} (GeV);Events", title.c_str()), "sb::off dt::stack lg::right"); */
/*     } */

     // write
    const string plots_path = Form("plots/%s/overlays", label.c_str());
    rt::Print(p, plots_path, suffix);

    // print yield explicitly
    // this is a kludge to the the x error bars the right size for the yeild plot
/*     gStyle->SetErrorX(0.3); */
/*     rt::Print(p["p_yield"], Form("%s/kin/p_yield", plots_path.c_str()), suffix); */
/*     gStyle->SetErrorX(); */
}
