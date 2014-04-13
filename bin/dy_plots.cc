// C++
#include <iostream>
#include <vector>
#include <functional>

// ROOT
#include "TChain.h"
#include "Math/VectorUtil.h"

// CMS2
#include "CMS2/NtupleMacrosHeader/interface/CMS2.h"

// CORE
#include "CMS2/NtupleMacrosCore/interface/mcSelections.h"
#include "CMS2/NtupleMacrosCore/interface/eventSelections.h"

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

// tools 
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "Analysis/DrellYan/interface/dySelections.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "AnalysisTools/CMS2Tools/interface/DileptonHypType.h"
#include "AnalysisTools/CMS2Tools/interface/GenHypType.h"

// -------------------------------------------------//
// Simple class to hold your analysis data
// -------------------------------------------------//

class DrellYanLooper
{
    public:
        // construct:
        DrellYanLooper() = delete;
        DrellYanLooper
        (
            const dy::Sample::Info sample_info, 
            const std::string& output_filename, 
            const double lumi,
            const long num_events,
            const bool verbose
        );

        // destroy:
        ~DrellYanLooper();

        // methods:
        void SetFilename(const std::string& output_filename);

        // basic methods:
        void BeginJob();
        void Analyze(const long event);
        void EndJob();

    private:

        // members:
        dy::Sample::Info m_sample_info; 
        std::string m_output_filename;
        std::string m_runlist_filename;
        double m_lumi;
        long m_num_events;
        bool m_verbose;
        rt::TH1Container hc;
};

// construct:
DrellYanLooper::DrellYanLooper
(
    const dy::Sample::Info sample_info, 
    const std::string& output_filename, 
    const double lumi,
    const long num_events,
    const bool verbose
)
    : m_sample_info(sample_info)
    , m_output_filename(output_filename)
    , m_lumi(lumi)
    , m_num_events(num_events)
    , m_verbose(verbose)
{
}

// destroy:
DrellYanLooper::~DrellYanLooper()
{
}

// ------------------------------------ //
// Stuff to do before job starts
// ------------------------------------ //

struct Selection
{
    enum value_type
    {
        ossf,
        mwin,
        svtx,
        trig,
        idiso,
        full,
        static_size
    };

    struct Info
    {
        std::string name;
        std::string title;
    };
};

// typedef const Selection::Info (&SelectionArrayRef)[Selection::static_size];
// SelectionArrayRef Selections()
std::vector<Selection::Info> Selections()
{
//     static const Selection::Info s_Selections[Selection::static_size] = 
    static const std::vector<Selection::Info> s_Selections = 
    {
        {"ossf" , "OSSF"                         },
        {"mwin" , "Mass Window"                  },
        {"svtx" , "Same Vertex"                  },
        {"trig" , "Passes Trigger"               },
        {"idiso", "Passes ID/Isolation"          },
        {"full" , "Full Selection"               }
    };
    assert(s_Selections.size()==Selection::static_size);
    return s_Selections;
}

void SetYieldAxisLabel(TH1* const hist)
{
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetBinLabel(1, ""      );
    hist->GetXaxis()->SetBinLabel(2, "ll"    );
    hist->GetXaxis()->SetBinLabel(3, "#mu#mu");
    hist->GetXaxis()->SetBinLabel(4, "ee"    );
}

void SetCutflowAxisLabel(TH1* const hist)
{
    hist->GetXaxis()->SetLabelSize(0.05);
    int bin = 1;
    for (const auto& sel : Selections())
    {
        hist->GetXaxis()->SetBinLabel(bin++, sel.name.c_str());
    }
}

void FillRecoHists(rt::TH1Container& hc, const int hyp_idx, const Selection::value_type sel, const double event_scale)
{
    // convenience variables
    const float dilep_mass = tas::hyp_p4().at(hyp_idx).mass();
    const int flavor_type  = tas::hyp_type().at(hyp_idx);
    const bool is_ee       = (flavor_type == 3);
    const bool is_mm       = (flavor_type == 0);

    // fill
    const std::string hist_prefix = "h_reco_" + Selections()[sel].name + "_";
    if (is_mm)
    {
        rt::Fill1D(hc[hist_prefix+"yield"], 1.0       , event_scale);
        rt::Fill1D(hc[hist_prefix+"mmm"  ], dilep_mass, event_scale);
        rt::Fill1D(hc["h_reco_cutflow_mm"], sel       , event_scale);
    }
    if (is_ee)
    {
        rt::Fill1D(hc[hist_prefix+"yield"], 2.0       , event_scale);
        rt::Fill1D(hc[hist_prefix+"mee"  ], dilep_mass, event_scale);
        rt::Fill1D(hc["h_reco_cutflow_ee"], sel       , event_scale);
    }
    rt::Fill1D(hc[hist_prefix+"yield"], 0.0, event_scale);
    rt::Fill1D(hc["h_reco_cutflow_ll"], sel, event_scale);
}

void DrellYanLooper::BeginJob()
{
    // gen level plots
    hc.Add(new TH1D("h_gen_yield" , "Yield count of gen level l^{+}l^{-}"          ,   4, -1,  3));
    hc.Add(new TH1D("h_gen_mee"   , "Generator level dielectron mass;m_{ee} (GeV)" , 150, 0, 150));
    hc.Add(new TH1D("h_gen_mmm"   , "Generator level dilmuon mass;m_{#mu#mu} (GeV)", 150, 0, 150));
    hc.Add(new TH1D("h_gen_mll"   , "Generator level dilepton mass;m_{ll} (GeV)"   , 150, 0, 150));

    // acceptance plots
    hc.Add(new TH1D("h_acc_gen_den", "Acceptence denominator;Channel;Event Count"        , 4, -1, 3));
    hc.Add(new TH1D("h_acc_gen_num", "Acceptence generator numerator;Channel;Event Count", 4, -1, 3));
    hc.Add(new TH1D("h_acc_rec_num", "Acceptence numerator;Channel;Event Count"          , 4, -1, 3));

    // reco level plots
    hc.Add(new TH1D("h_reco_cutflow_ee", "e^{+}e^{-}  Cut Flow;Selection;Event Count"    , Selection::static_size, 0, Selection::static_size));
    hc.Add(new TH1D("h_reco_cutflow_mm", "#mu^{+}#mu^{-}  Cut Flow;Selection;Event Count", Selection::static_size, 0, Selection::static_size));
    hc.Add(new TH1D("h_reco_cutflow_ll", "l^{+}l^{-} Cut Flow;Selection;Event Count"     , Selection::static_size, 0, Selection::static_size));
    SetCutflowAxisLabel(hc["h_reco_cutflow_ee"]);
    SetCutflowAxisLabel(hc["h_reco_cutflow_mm"]);
    SetCutflowAxisLabel(hc["h_reco_cutflow_ll"]);
    for (const auto& sel : Selections())
    {
        hc.Add(new TH1D(Form("h_reco_%s_yield", sel.name.c_str()), Form("%s yield count of reco level l^{+}l^{-}", sel.title.c_str()),   4, -1,  3));
        hc.Add(new TH1D(Form("h_reco_%s_mee"  , sel.name.c_str()), Form("%s dielectron mass;m_{ee} (GeV)"        , sel.title.c_str()), 150, 0, 150));
        hc.Add(new TH1D(Form("h_reco_%s_mmm"  , sel.name.c_str()), Form("%s dilmuon mass;m_{#mu#mu} (GeV)"       , sel.title.c_str()), 150, 0, 150));
        hc.Add(new TH1D(Form("h_reco_%s_mll"  , sel.name.c_str()), Form("%s dilepton mass;m_{ll} (GeV)"          , sel.title.c_str()), 150, 0, 150));
        SetYieldAxisLabel(hc["h_reco_"+sel.name+"_yield"]);
    }

    // change axis labels
    SetYieldAxisLabel(hc["h_gen_yield"  ]);
    SetYieldAxisLabel(hc["h_acc_gen_den"]);
    SetYieldAxisLabel(hc["h_acc_gen_num"]);
    SetYieldAxisLabel(hc["h_acc_rec_num"]);

    // sumw2()
    hc.Sumw2();
}

// ------------------------------------ //
// Stuff to do on each event 
// ------------------------------------ //

void DrellYanLooper::Analyze(const long event)
{
    if (m_verbose)
    {
        std::cout << "[DrellYanLooper] Running on run, ls, event: " 
            << tas::evt_run()       << ", "
            << tas::evt_lumiBlock() << ", "
            << tas::evt_event()     << ", "
            << std::endl;
    }

    // Event Cleaning
    // ---------------------- // 

    // require at least 3 tracks in the event
    if (tas::trks_trk_p4().size() < 3)
    {
        if (m_verbose) {std::cout << "fails # trks >= 3 requirement" << std::endl;}
        return;
    }

    // require standard cleaning 
    if (!cleaning_standardNovember2011()) 
    {
        if (m_verbose) {std::cout << "fails November2011 cleaning requirement" << std::endl;}
        return;
    }

    // event scale factors 
    // ---------------------- // 

    double event_scale = 1.0;
    if (!tas::evt_isRealData())
    {
        const double nevts_cms2  = m_sample_info.filter_eff * tas::evt_nEvts(); // number of events run in CMSSW job to make ntuple
        const double nevts_scale = nevts_cms2/m_num_events;                     // scale up the weight to account fo lower stats
        const double scale1fb    = tas::evt_scale1fb();                         // scale1fb stored in event
        event_scale              = m_lumi * scale1fb * nevts_scale;
    }
    if (m_verbose) {std::cout << "event_scale = " << event_scale << "\n";}

    // generator level plots
    // ---------------------- // 

    bool is_gen_ee      = false;
    bool is_gen_mm      = false;
    bool passes_acc_den = false;

    if (!tas::evt_isRealData())
    {
        const std::vector<at::GenHyp> gen_hyps = at::GetGenHyps(/*min_pt=*/0.0, /*max_eta=*/1000.0);
        const std::vector<at::GenHyp> gen_hyps_clean = lt::filter_container(gen_hyps,
            [](const at::GenHyp& h)
            {
                return (h.IsOS() and (h.IsEE_IncludeTaus() or h.IsMuMu_IncludeTaus()));
            }
        );

        if (!gen_hyps_clean.empty())
        {
            // observables:
            passes_acc_den            = true; 
            const at::GenHyp& gen_hyp = gen_hyps_clean.front();
            const double gen_mass     = gen_hyp.P4().mass();

            // fill hists
            if (gen_hyp.IsMuMu_IncludeTaus()) {hc["h_gen_yield"]->Fill(1.0, event_scale);}
            if (gen_hyp.IsEE_IncludeTaus()  ) {hc["h_gen_yield"]->Fill(2.0, event_scale);}
            hc["h_gen_yield"]->Fill(0.0, event_scale);

            // kinematics
            if (gen_hyp.IsMuMu_IncludeTaus()) {rt::Fill1D(hc["h_gen_mee"], gen_mass, event_scale);}
            if (gen_hyp.IsEE_IncludeTaus()  ) {rt::Fill1D(hc["h_gen_mmm"], gen_mass, event_scale);}
            rt::Fill1D(hc["h_gen_mll"], gen_mass, event_scale);
        }

        // acceptence denominator
        const std::vector<at::GenHyp> gen_hyps_acc_den = lt::filter_container(gen_hyps_clean,
            [](const at::GenHyp& h)
            {
                return (h.IsFromZ() and (60 < h.P4().mass() && h.P4().mass() < 120));
            }
        );
        if (!gen_hyps_acc_den.empty())
        {
            // observables:
            passes_acc_den            = true; 
            const at::GenHyp& gen_hyp = gen_hyps_acc_den.front();
            is_gen_mm                 = gen_hyp.IsMuMu_IncludeTaus();
            is_gen_ee                 = gen_hyp.IsEE_IncludeTaus();

            // fill hists 
            if (is_gen_mm) {hc["h_acc_gen_den"]->Fill(1.0, event_scale);}
            if (is_gen_ee) {hc["h_acc_gen_den"]->Fill(2.0, event_scale);}
            hc["h_acc_gen_den"]->Fill(0.0, event_scale);
        }

        // acceptence numerator
        const std::vector<at::GenHyp> gen_hyps_acc_num = lt::filter_container(gen_hyps_acc_den,
            [](const at::GenHyp& h)
            {
                return ((h.IsEE() or h.IsMuMu()) and h.IsAccepted(/*min_pt=*/25.0, /*max_eta=*/2.5));
            }
        );

        if (passes_acc_den and !gen_hyps_acc_num.empty())
        {
            // observables:
            const at::GenHyp& gen_hyp = gen_hyps_acc_num.front();

            // fill hists
            if (gen_hyp.IsMuMu_IncludeTaus()) {hc["h_acc_gen_num"]->Fill(1.0, event_scale);}
            if (gen_hyp.IsEE_IncludeTaus()  ) {hc["h_acc_gen_num"]->Fill(2.0, event_scale);}
            hc["h_acc_gen_num"]->Fill(0.0, event_scale);
        }
    }

    // reco level plots
    // ---------------------- // 

    // loop over hypotheses
    int best_hyp = -1;
    for (size_t hyp_idx = 0; hyp_idx < tas::hyp_type().size(); ++hyp_idx)
    {                
        // convenience variables
        const float dilep_mass = tas::hyp_p4().at(hyp_idx).mass();
        const int flavor_type  = tas::hyp_type().at(hyp_idx);
        const bool is_ee       = (flavor_type == 3);
        const bool is_mm       = (flavor_type == 0);

        // apply selections

        // OSSF
        if (tas::hyp_lt_charge().at(hyp_idx) == tas::hyp_ll_charge().at(hyp_idx))
        {
            if (m_verbose) {std::cout << "not OS" << std::endl;}
            continue;
        }
        if (not(is_ee or is_mm))
        {
            if (m_verbose) {std::cout << "not SF" << std::endl;}
            continue;
        }
        FillRecoHists(hc, hyp_idx, Selection::ossf, event_scale);
    
        // 60 < m_ll << 120 GeV
        if (not (60 < dilep_mass && dilep_mass < 120.0))        
        {
            if (m_verbose) {std::cout << "not SF" << std::endl;}
            continue;
        }
        FillRecoHists(hc, hyp_idx, Selection::mwin, event_scale);
    
        // both leptons from first vertex
        if (not hypsFromFirstGoodVertex(hyp_idx))
        {
            if (m_verbose) {std::cout << "did not pass same vertex" << std::endl;}
            continue;
        }
        FillRecoHists(hc, hyp_idx, Selection::svtx, event_scale);
    
        // trigger
        if (not dy::passesTrigger(tas::hyp_type().at(hyp_idx)))
        {
            if (m_verbose) {std::cout << "did not pass trigger" << std::endl;}
            continue;
        }
        FillRecoHists(hc, hyp_idx, Selection::trig, event_scale);
    
        // l1 and l2 pass selection
        if (not dy::isSelectedHypothesis(hyp_idx))
        {
            if (m_verbose) {std::cout << "not selected" << std::endl;}
            continue;
        }
        FillRecoHists(hc, hyp_idx, Selection::idiso, event_scale);
    
        // if we're here, then good event :)
        best_hyp = dy::ChooseBetterHypothesis(best_hyp, hyp_idx);

    } // end loop over hypothesis

    // only continue if hyp has been selected
    // all: 0, mm: 1, em: 2, ee: 3
    const int hyp_idx = best_hyp;
    if (hyp_idx < 0)
    {
        if (m_verbose) {std::cout << "no good hypthesis chosen" << std::endl;}
        return;
    }
    else
    {
        if (m_verbose) {std::cout << "hypthesis chosen: " << hyp_idx << " of " << tas::hyp_p4().size() << std::endl;}
    }
    FillRecoHists(hc, hyp_idx, Selection::full, event_scale);

    // acceptance
    const at::DileptonHypType::value_type reco_flavor_type = at::hyp_typeToHypType(tas::hyp_type().at(hyp_idx));
    const bool is_ee          = (reco_flavor_type == at::DileptonHypType::EE);
    const bool is_mm          = (reco_flavor_type == at::DileptonHypType::MUMU);
    const bool passes_acc_num = passes_acc_den and ((is_gen_ee and is_ee) or (is_gen_mm and is_mm));
    if (passes_acc_num)
    {
        if (is_mm) {hc["h_acc_rec_num"]->Fill(1.0, event_scale);}
        if (is_ee) {hc["h_acc_rec_num"]->Fill(2.0, event_scale);}
        hc["h_acc_rec_num"]->Fill(0.0             , event_scale);
    }

    // done with event
    // ---------------------- // 

    if (m_verbose)
    {
        std::cout << "[DrellYanLooper] Finished with run, ls, event: " 
            << tas::evt_run()       << ", "
            << tas::evt_lumiBlock() << ", "
            << tas::evt_event()     << ", "
            << "\n" << std::endl;
    }
}

// ------------------------------------ //
// Stuff to do after job finishes
// ------------------------------------ //
void DrellYanLooper::EndJob()
{
    // calculate accecptance
    hc.Add(rt::DivideHists(hc["h_acc_rec_num"], hc["h_acc_gen_den"], "h_acc_rec", "Reco Acceptence;Channel;Event Count"));
    hc.Add(rt::DivideHists(hc["h_acc_gen_num"], hc["h_acc_gen_den"], "h_acc_gen", "Gen Acceptence;Channel;Event Count" ));
    SetYieldAxisLabel(hc["h_acc_rec"]);
    SetYieldAxisLabel(hc["h_acc_gen"]);

    // acceptance 
    if (m_sample_info.sample != dy::Sample::data)
    {
        const dy::Yield y_acc_gen_den = dy::GetYieldFromHist(*hc["h_acc_gen_den"]);
        const dy::Yield y_acc_gen_num = dy::GetYieldFromHist(*hc["h_acc_gen_num"]);
        const dy::Yield y_acc_rec_num = dy::GetYieldFromHist(*hc["h_acc_rec_num"]);
        const dy::Yield y_acc_gen     = dy::GetYieldFromHist(*hc["h_acc_gen"    ]);
        const dy::Yield y_acc_rec     = dy::GetYieldFromHist(*hc["h_acc_rec"    ]);
        CTable t_gen_yields;
        t_gen_yields.useTitle();
        t_gen_yields.setTitle("Acceptance for Drell-Yan Exercise");
        t_gen_yields.setTable()
        (                                                     "ee",                        "mm",                       "ll")
        ("Acceptance Denominator"     , y_acc_gen_den.ee.pm("4.0"),  y_acc_gen_den.ee.pm("4.0"), y_acc_gen_den.ll.pm("4.0"))
        ("Acceptance Numerator (gen)" , y_acc_gen_num.ee.pm("4.0"),  y_acc_gen_num.ee.pm("4.0"), y_acc_gen_num.ll.pm("4.0"))
        ("Acceptance Numerator (reco)", y_acc_rec_num.ee.pm("4.0"),  y_acc_rec_num.ee.pm("4.0"), y_acc_rec_num.ll.pm("4.0"))
        ("Gen Acceptance"             ,     y_acc_gen.ee.pm("4.3"),      y_acc_gen.ee.pm("4.3"),     y_acc_gen.ll.pm("4.3"))
        ("Reco Acceptance"            ,     y_acc_rec.ee.pm("4.3"),      y_acc_rec.ee.pm("4.3"),     y_acc_rec.ll.pm("4.3"))
        ;
        std::cout << t_gen_yields << std::endl;
    }

    // yields 
    CTable t_reco_yields;
    t_reco_yields.useTitle();
    t_reco_yields.setTitle("yields for Drell-Yan Exercise");
    t_reco_yields.setTable() ("ee", "mm", "ll");
    for (size_t i = 0; i < Selections().size(); ++i)
    {
        const auto s = Selections()[i];
        const dy::Yield yield = dy::GetYieldFromHist(*hc["h_reco_"+s.name+"_yield"]);
        t_reco_yields.setCell(yield.ee.pm("4.1"), i, 0);
        t_reco_yields.setCell(yield.mm.pm("4.1"), i, 1);
        t_reco_yields.setCell(yield.ll.pm("4.1"), i, 2);
        t_reco_yields.setRowLabel(s.title, i);
    }
    std::cout << t_reco_yields << std::endl;

    std::cout << "[DrellYanLooper] Saving hists to output file: " << m_output_filename << std::endl;
    hc.Write(m_output_filename);
    if (m_verbose) hc.List();
    return;
}

void DrellYanLooper::SetFilename(const std::string& filename)
{
    m_output_filename = filename;
}
   
// -------------------------------------------------//
// main program
// -------------------------------------------------//

#include "AnalysisTools/CMS2Tools/interface/ScanChain.h"
#include "AnalysisTools/CMS2Tools/interface/CMS2Wrapper.h"
#include "AnalysisTools/CMS2Tools/interface/LoadFWLite.h"
#include <algorithm> 
#include "boost/program_options.hpp"

int main(int argc, char **argv)
try
{
    // parse the inputs
    // -------------------------------------------------------------------------------------------------//

    // defaults 
    long long number_of_events = -1;
    std::string sample_name    = "";
    std::string sample_pset    = "psets/dy_samples_cfg.py";
    std::string input_file     = "";
    std::string output_file    = "";
    std::string label          = "test";
    std::string run_list       = "json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_cms2.txt";
    double lumi                = 0.082;
    bool verbose               = false;
    int event                  = -1;
    int run                    = -1;
    int ls                     = -1;
    
    // parse arguments
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help"       , "print this menu")
        ("sample"     , po::value<std::string>(&sample_name)->required(), "REQUIRED: Sample name to run on (from Sample.h)"        )
        ("sample_pset", po::value<std::string>(&sample_pset)            , "pset to change the value of Sample::Info(from Sample.h)")
        ("input"      , po::value<std::string>(&input_file)             , "input ROOT file (or csv)"                               )
        ("output"     , po::value<std::string>(&output_file)            , "output ROOT file"                                       )
        ("nevts"      , po::value<long long>(&number_of_events)         , "maximum number of events to skim"                       )
        ("label"      , po::value<std::string>(&label)                  , "unique output label to keep differnet jobs straight"    )
        ("run_list"   , po::value<std::string>(&run_list)               , "good run list (empty == none)"                          )
        ("lumi"       , po::value<double>(&lumi)                        , "luminosity (default -1)"                                )
        ("verbose"    , po::value<bool>(&verbose)                       , "verbosity toggle"                                       )
        ("event"      , po::value<int>(&event)                          , "specific event to run on (-1 == all events)"            )
        ("run"        , po::value<int>(&run)                            , "specific run to run on (-1 == all events)"              )
        ("ls"         , po::value<int>(&ls)                             , "specific lumi section to run on (-1 == all events)"     )
        ;
    try
    {
        // first parse command line options
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) 
        {
            std::cout << desc << "\n";
            return 1;
        }
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << "\nexiting" << std::endl;
        std::cout << desc << "\n";
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown error!" << "\n";
        return 1;
    }

    // print the inputs  
    std::cout << "[dy_plots] inputs:\n";
    std::cout << "sample           = " << sample_name      << "\n";
    std::cout << "sample_pset      = " << sample_pset      << "\n";
    std::cout << "input            = " << input_file       << "\n";
    std::cout << "output           = " << output_file      << "\n";
    std::cout << "label            = " << label            << "\n";
    std::cout << "nevts            = " << number_of_events << "\n";
    std::cout << "run_list         = " << run_list         << "\n";
    std::cout << "lumi             = " << lumi             << "\n";
    std::cout << "verbose          = " << verbose          << "\n";
    std::cout << "event            = " << event            << "\n";
    std::cout << "run              = " << run              << "\n";
    std::cout << "ls               = " << ls               << "\n";
    std::cout << std::endl;

    // check inputs 
    // -------------------------------------------------------------------------------------------------//

    // sample info
    if (sample_pset.empty() || !lt::file_exists(sample_pset))
    {
        throw std::invalid_argument(Form("[dy_plots] pset %s does not exist.\n", sample_pset.c_str()));
    }
    dy::Sample::SetPsetPath(sample_pset);
    dy::Sample::Info sample_info = dy::GetSampleInfo(sample_name);

    // get the chain
    at::LoadFWLite();
    std::vector<std::string> chain_files = lt::string_split(input_file.empty() ? sample_info.ntuple_path : input_file);
    TChain * const chain = rt::CreateTChain("Events", chain_files);

    // check the output
    if (output_file.empty())
    {
        output_file = Form("plots/%s/%s_plots.root", label.c_str(), sample_info.name.c_str());
    }

    // run the looper
    // -------------------------------------------------------------------------------------------------//

    std::cout << "[dy_plots] TChain set to run on:\n";
    rt::PrintFilesFromTChain(chain);

    // create looper
    const long long num_events_to_run = (number_of_events < 0 ? chain->GetEntries() : std::min(number_of_events, chain->GetEntries()));
    DrellYanLooper looper
    (
        sample_info,
        output_file, 
        lumi,
        num_events_to_run,
        verbose
    );
    std::cout << "[dy_plots] running drell-yan plotting looper...\n";

    // simple style
    rt::SetStyle();

    // scan the chain
    at::ScanChain<CMS2>
    (
        chain, 
        looper,
        cms2,
        num_events_to_run,
        run_list,
        /*fast=*/true,
        verbose,
        run,
        ls,
        event
    ); 

    // cleanup
    delete chain;

    return 0;
}
catch (std::exception& e)
{
    std::cerr << "[dy_plots] Error: failed..." << std::endl;
    std::cerr << e.what() << std::endl;
    return 1;
}

