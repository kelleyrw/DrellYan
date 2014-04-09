// C++
#include <iostream>
#include <vector>

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
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// tools 
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "Analysis/DrellYan/interface/dySelections.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "AnalysisTools/CMS2Tools/interface/DileptonChargeType.h"
#include "AnalysisTools/CMS2Tools/interface/DileptonHypType.h"

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

void SetYieldAxisLabel(TH1* const hist)
{
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetBinLabel(1, "ll"    );
    hist->GetXaxis()->SetBinLabel(2, "#mu#mu");
    hist->GetXaxis()->SetBinLabel(3, "e#mu"  );
    hist->GetXaxis()->SetBinLabel(4, "ee"    );
}

void DrellYanLooper::BeginJob()
{
    // gen level plots
    hc.Add(new TH1D("h_gen_yield" , "Yield count of gen  level l^{+}l^{-}"         ,   4, 0,   4));
    hc.Add(new TH1D("h_gen_mee"   , "Generator level dielectron mass;m_{ee} (GeV)" , 150, 0, 150));
    hc.Add(new TH1D("h_gen_mmm"   , "Generator level dilmuon mass;m_{#mu#mu} (GeV)", 150, 0, 150));
    hc.Add(new TH1D("h_gen_mll"   , "Generator level dilepton mass;m_{ll} (GeV)"   , 150, 0, 150));

    // reco level plots
    hc.Add(new TH1D("h_reco_yield", "Yield count of reco level l^{+}l^{-}",   4,  0,   4));
    hc.Add(new TH1D("h_reco_mee"  , "Final dielectron mass;m_{ee} (GeV)"  ,  60, 60, 120));
    hc.Add(new TH1D("h_reco_mmm"  , "Final dilmuon mass;m_{#mu#mu} (GeV)" ,  60, 60, 120));
    hc.Add(new TH1D("h_reco_mll"  , "Final dilepton mass;m_{ll} (GeV)"    ,  60, 60, 120));

    // change axis labels
    SetYieldAxisLabel(hc["h_gen_yield" ]);
    SetYieldAxisLabel(hc["h_reco_yield"]);

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

    // event scale factors 
    // ---------------------- // 

    double event_scale = 1.0;
    if (!tas::evt_isRealData())
    {
        const double nevts_aod     = tas::evt_nEvts();          // number of events run in CMSSW job to make ntuple
        const double nevts_file    = m_num_events;              // number of events in the current job
        const double nevts_scale   = nevts_aod/nevts_file;      // scale up the weight to account fo lower stats
        const double  scale1fb     = tas::evt_scale1fb();       // scale1fb stored in event
        const double sd_filter_eff = m_sample_info.filter_eff;  // number of events after the "EDFilter"
        event_scale                = m_lumi * scale1fb * nevts_scale * sd_filter_eff;
    }
    if (m_verbose) {std::cout << "event_scale = " << event_scale << "\n";}

    // generator level plots
    // ---------------------- // 

    if (!tas::evt_isRealData())
    {
        const int gen_flavor_type = dy::GenDileptonType(); 
        if (m_verbose) {std::cout << "gen_flavor_type = " << gen_flavor_type << "\n";}

        // flavor
        if (gen_flavor_type > -1)
        {
            hc["h_gen_yield"]->Fill(gen_flavor_type, event_scale);
            hc["h_gen_yield"]->Fill(0.0            , event_scale);
        }

        // kinematics
        //         if (gen_flavor_type > -1)
        //         {
        //             hc["h_gen_"]->Fill(gen_flavor_type, event_scale);
        //             hc["h_gen_"]->Fill(0.0            , event_scale);
        //         }
    }

    // reco level plots
    // ---------------------- // 

    // loop over hypotheses
    int best_hyp = -1;
    for (size_t hyp_idx = 0; hyp_idx != tas::hyp_type().size(); hyp_idx++)
    {                
        // convenience variables
//         const int lt_id                                   = tas::hyp_lt_id().at(hyp_idx);
//         const int ll_id                                   = tas::hyp_ll_id().at(hyp_idx);
//         const int lt_idx                                  = tas::hyp_lt_index().at(hyp_idx);
//         const int ll_idx                                  = tas::hyp_ll_index().at(hyp_idx);
//         const float dilep_mass                            = tas::hyp_p4().at(hyp_idx).mass();
        const at::DileptonHypType::value_type flavor_type = at::hyp_typeToHypType(tas::hyp_type().at(hyp_idx));

        // apply selections
//         if (tas::evt_isRealData() && !dy::passesTrigger(tas::hyp_type().at(hyp_idx)))                {continue;}
        if ((tas::hyp_lt_charge().at(hyp_idx) * tas::hyp_ll_charge().at(hyp_idx)) > 0)               {continue;}
        if (not(flavor_type == at::DileptonHypType::EE or flavor_type == at::DileptonHypType::MUMU)) {continue;}
//         if (not (60 < dilep_mass && dilep_mass < 120.0))                                             {continue;}
//         if (not dy::isSelectedLepton(lt_id, lt_idx))                                                 {continue;}
//         if (not dy::isSelectedLepton(ll_id, ll_idx))                                                 {continue;}
//         if (not hypsFromFirstGoodVertex(hyp_idx))                                                    {continue;}

        // if we're here, then good event :)
        best_hyp = dy::ChooseBetterHypothesis(best_hyp, hyp_idx);

    } // end looper over hypothesis

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

    // observables
    const at::DileptonHypType::value_type reco_flavor_type = at::hyp_typeToHypType(tas::hyp_type().at(hyp_idx));
//     const LorentzVector& lt_p4  = tas::hyp_lt_p4().at(hyp_idx);
//     const LorentzVector& ll_p4  = tas::hyp_ll_p4().at(hyp_idx);
    const LorentzVector& hyp_p4 = tas::hyp_p4().at(hyp_idx);
//     const LorentzVector& l1_p4  = (lt_p4.pt() > ll_p4.pt() ? tas::hyp_lt_p4().at(hyp_idx)    : tas::hyp_ll_p4().at(hyp_idx)   );
//     const LorentzVector& l2_p4  = (lt_p4.pt() > ll_p4.pt() ? tas::hyp_ll_p4().at(hyp_idx)    : tas::hyp_lt_p4().at(hyp_idx)   );
//     const int l1_idx            = (lt_p4.pt() > ll_p4.pt() ? tas::hyp_lt_index().at(hyp_idx) : tas::hyp_ll_index().at(hyp_idx));
//     const int l2_idx            = (lt_p4.pt() > ll_p4.pt() ? tas::hyp_ll_index().at(hyp_idx) : tas::hyp_lt_index().at(hyp_idx));
//     const int l1_id             = (lt_p4.pt() > ll_p4.pt() ? tas::hyp_lt_id().at(hyp_idx)    : tas::hyp_ll_id().at(hyp_idx)   );
//     const int l2_id             = (lt_p4.pt() > ll_p4.pt() ? tas::hyp_ll_id().at(hyp_idx)    : tas::hyp_lt_id().at(hyp_idx)   );

    // flavor bools
    const bool is_ee = (reco_flavor_type == at::DileptonHypType::EE);
    const bool is_mm = (reco_flavor_type == at::DileptonHypType::MUMU);

    // fill hist
    hc["h_reco_yield"]->Fill(reco_flavor_type, event_scale);
    hc["h_reco_yield"]->Fill(0.0             , event_scale);

    if (is_mm) {rt::Fill1D(hc["h_reco_mmm"], hyp_p4.mass(), event_scale);}
    if (is_ee) {rt::Fill1D(hc["h_reco_mee"], hyp_p4.mass(), event_scale);}
    rt::Fill1D(hc["h_reco_mll"], hyp_p4.mass(), event_scale);

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
    // output counts
    dy::Yield gen_ = dy::GetYieldFromHist(*hc["h_gen_yield" ]);
    std::cout << dy::GetYieldString(gen_, "Gen level Yields") << std::endl;

    dy::Yield reco_ = dy::GetYieldFromHist(*hc["h_reco_yield" ]);
    std::cout << dy::GetYieldString(reco_, "Reco level Yields") << std::endl;

    std::cout << "[DrellYanLooper] Saving hists to output file: " << m_output_filename << std::endl;
    hc.Write(m_output_filename);
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

    // get the inputs 
    long long number_of_events = -1; 
    std::string sample_name = "";
    std::string pset_filename = "";
    std::string input_file = ""; 
    std::string output_file = ""; 
    std::string label = "test"; 
    std::string run_list = ""; 
    double lumi = 1.0;
    bool verbose = false;
    int event = -1;
    int run = -1;
    int ls = -1;
    
    // parse arguments
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help"     , "print this menu")
        ("pset"     , po::value<std::string>(&pset_filename)->required() , "REQUIRED: python configuration file"                )
        ("nevts"    , po::value<long long>(&number_of_events)            , "maximum number of events to skim"                   )
        ("sample"   , po::value<std::string>(&sample_name)               , "Sample name to run on (from Sample.h)"              )
        ("input"    , po::value<std::string>(&input_file)                , "input ROOT file (or csv)"                           )
        ("output"   , po::value<std::string>(&output_file)               , "output ROOT file"                                   )
        ("label"    , po::value<std::string>(&label)                     , "unique output label to keep differnet jobs straight")
        ("run_list" , po::value<std::string>(&run_list)                  , "good run list (empty == none)"                      )
        ("lumi"     , po::value<double>(&lumi)                           , "luminosity (default -1)"                            )
        ("verbose"  , po::value<bool>(&verbose)                          , "verbosity toggle"                                   )
        ("event"    , po::value<int>(&event)                             , "specific event to run on (-1 == all events)"        )
        ("run"      , po::value<int>(&run)                               , "specific run to run on (-1 == all events)"          )
        ("ls"       , po::value<int>(&ls)                                , "specific lumi section to run on (-1 == all events)" )
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

        // get the inputs from pset

        // check that pset contains "process" 
        if (!edm::readPSetsFrom(pset_filename)->existsAs<edm::ParameterSet>("process"))
        {
            throw std::invalid_argument(Form("[dy_plots] Error: ParametersSet 'process' is missing in your configuration file"));
        }

        // get the python configuration
        const edm::ParameterSet& process = edm::readPSetsFrom(pset_filename)->getParameter<edm::ParameterSet>("process");
        const edm::ParameterSet& cfg     = process.getParameter<edm::ParameterSet>("dy_plots");

        number_of_events = cfg.getUntrackedParameter<long long>   ( "nevts"    , number_of_events);
        sample_name      = cfg.getUntrackedParameter<std::string> ( "sample"   , sample_name     );
        input_file       = cfg.getUntrackedParameter<std::string> ( "input"    , input_file      );
        output_file      = cfg.getUntrackedParameter<std::string> ( "output"   , output_file     );
        label            = cfg.getUntrackedParameter<std::string> ( "label"    , label           );
        run_list         = cfg.getUntrackedParameter<std::string> ( "run_list" , run_list        );
        lumi             = cfg.getUntrackedParameter<double>      ( "lumi"     , lumi            );
        verbose          = cfg.getUntrackedParameter<bool>        ( "verbose"  , verbose         );
        event            = cfg.getUntrackedParameter<int>         ( "event"    , event           );
        run              = cfg.getUntrackedParameter<int>         ( "run"      , run             );
        ls               = cfg.getUntrackedParameter<int>         ( "ls"       , ls              );

        // now parse command line again to see if there are overrides
        po::store(po::parse_command_line(argc, argv, desc), vm);
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

