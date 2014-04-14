// C++
#include <iostream>
#include <vector>
#include <functional>

// ROOT
#include "TChain.h"
#include "Math/LorentzVector.h"
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

// typedefs
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<LorentzVector> vecLorentzVector;
typedef std::vector<float> vecd;
typedef std::vector<int> veci;

// -------------------------------------------------//
// Class to hold the ntuple information 
// -------------------------------------------------//

class DrellYanInfo
{
public:

    // constructors and destructor
    DrellYanInfo();

    void Reset();
    void SetBranches(TTree& tree);
//     void SetAliases(TTree& tree) const;

    // event level info
    int run;
    int ls;
    int evt;
    int sample;
    std::string dataset;
    std::string filename;
    bool is_real_data;
    float scale1fb;
    float scale1fb_cms2;
    float xsec;
    float nevts_aod;
    float nevts_cms2;
    float nevts_file;
    float kfactor;

    // some reco event level info
    int nvtxs;
    float pfmet;
    float pfmet_phi;
    int pu_nvtxs;
    float pu_ntrueint;

    // trigger
    bool trig_dmu;
    bool trig_del;
    bool trig_smu;
    bool trig_sel;
    bool trig;

    // Gen hypothesis specific info 
    int gen_hyp_type;
    bool is_gen_ee;
    bool is_gen_mm;
    bool is_gen_tt;
    bool is_gen_ee_includetau;
    bool is_gen_mm_includetau;
    bool is_gen_fromz;
    bool is_gen_acc_den;
    bool is_gen_acc_num;
    LorentzVector gen_p4;
    LorentzVector gen_lep1_p4;
    int gen_lep1_id;
    int gen_lep1_charge;
    LorentzVector gen_lep2_p4;
    int gen_lep2_id;
    int gen_lep2_charge;
};

DrellYanInfo::DrellYanInfo()
    : run                  (-999999   ) 
    , ls                   (-999999   ) 
    , evt                  (-999999   ) 
    , sample               (-999999   ) 
    , dataset              (""        ) 
    , filename             (""        ) 
    , is_real_data         (false     ) 
    , scale1fb             (-999999   ) 
    , scale1fb_cms2        (-999999   ) 
    , xsec                 (-999999   ) 
    , nevts_aod            (-999999   ) 
    , nevts_cms2           (-999999   ) 
    , nevts_file           (-999999   ) 
    , kfactor              (-999999   ) 
    , nvtxs                (-999999   ) 
    , pfmet                (-999999   ) 
    , pfmet_phi            (-999999   ) 
    , pu_nvtxs             (-999999   ) 
    , pu_ntrueint          (-999999   ) 
    , trig_dmu             (-999999   ) 
    , trig_del             (-999999   ) 
    , trig_smu             (-999999   ) 
    , trig_sel             (-999999   ) 
    , trig                 (-999999   ) 
    , gen_hyp_type         (-999999   ) 
    , is_gen_ee            (-999999   ) 
    , is_gen_mm            (-999999   ) 
    , is_gen_tt            (-999999   ) 
    , is_gen_ee_includetau (-999999   ) 
    , is_gen_mm_includetau (-999999   ) 
    , is_gen_fromz         (-999999   ) 
    , is_gen_acc_den       (-999999   ) 
    , is_gen_acc_num       (-999999   ) 
    , gen_p4               (0, 0, 0, 0) 
    , gen_lep1_p4          (0, 0, 0, 0) 
    , gen_lep1_id          (-999999   ) 
    , gen_lep1_charge      (-999999   ) 
    , gen_lep2_p4          (0, 0, 0, 0) 
    , gen_lep2_id          (-999999   ) 
    , gen_lep2_charge      (-999999   ) 
{
}

void DrellYanInfo::Reset()
{
    run                  = -999999; 
    ls                   = -999999; 
    evt                  = -999999; 
    sample               = -999999; 
    dataset              = ""; 
    filename             = ""; 
    is_real_data         = false; 
    scale1fb             = -999999; 
    scale1fb_cms2        = -999999; 
    xsec                 = -999999; 
    nevts_aod            = -999999; 
    nevts_cms2           = -999999; 
    nevts_file           = -999999; 
    kfactor              = -999999; 
    nvtxs                = -999999; 
    pfmet                = -999999; 
    pfmet_phi            = -999999; 
    pu_nvtxs             = -999999; 
    pu_ntrueint          = -999999; 
    trig_dmu             = -999999; 
    trig_del             = -999999; 
    trig_smu             = -999999; 
    trig_sel             = -999999; 
    trig                 = -999999; 
    gen_hyp_type         = -999999; 
    is_gen_ee            = -999999; 
    is_gen_mm            = -999999; 
    is_gen_tt            = -999999; 
    is_gen_ee_includetau = -999999; 
    is_gen_mm_includetau = -999999; 
    is_gen_fromz         = -999999; 
    is_gen_acc_den       = -999999; 
    is_gen_acc_num       = -999999; 
    gen_p4               = LorentzVector(0, 0, 0, 0); 
    gen_lep1_p4          = LorentzVector(0, 0, 0, 0); 
    gen_lep1_id          = -999999; 
    gen_lep1_charge      = -999999; 
    gen_lep2_p4          = LorentzVector(0, 0, 0, 0); 
    gen_lep2_id          = -999999; 
    gen_lep2_charge      = -999999; 
}
    
void DrellYanInfo::SetBranches(TTree& tree)
{
    tree.Branch("run"                  , &run                  , "run/I");
    tree.Branch("ls"                   , &ls                   , "ls/I" );
    tree.Branch("evt"                  , &evt                  , "evt/I");
    tree.Branch("sample"               , &sample               );
    tree.Branch("dataset"              , &dataset              );
    tree.Branch("filename"             , &filename             );
    tree.Branch("is_real_data"         , &is_real_data         );
    tree.Branch("nvtxs"                , &nvtxs                );
    tree.Branch("pfmet"                , &pfmet                );
    tree.Branch("pfmet_phi"            , &pfmet_phi            );
    tree.Branch("pu_nvtxs"             , &pu_nvtxs             );
    tree.Branch("pu_ntrueint"          , &pu_ntrueint          );
    tree.Branch("scale1fb"             , &scale1fb             );
    tree.Branch("scale1fb_cms2"        , &scale1fb_cms2        );
    tree.Branch("xsec"                 , &xsec                 );
    tree.Branch("nevts_aod"            , &nevts_aod            );
    tree.Branch("nevts_cms2"           , &nevts_cms2           );
    tree.Branch("nevts_file"           , &nevts_file           );
    tree.Branch("kfactor"              , &kfactor              );
    tree.Branch("trig_dmu"             , &trig_dmu             );
    tree.Branch("trig_del"             , &trig_del             );
    tree.Branch("trig_smu"             , &trig_smu             );
    tree.Branch("trig_sel"             , &trig_sel             );
    tree.Branch("trig"                 , &trig                 );
    tree.Branch("gen_hyp_type"         , &gen_hyp_type         );
    tree.Branch("is_gen_ee"            , &is_gen_ee            );
    tree.Branch("is_gen_mm"            , &is_gen_mm            );
    tree.Branch("is_gen_tt"            , &is_gen_tt            );
    tree.Branch("is_gen_ee_includetau" , &is_gen_ee_includetau );
    tree.Branch("is_gen_mm_includetau" , &is_gen_mm_includetau );
    tree.Branch("is_gen_fromz"         , &is_gen_fromz         );
    tree.Branch("is_gen_acc_den"       , &is_gen_acc_den       );
    tree.Branch("is_gen_acc_num"       , &is_gen_acc_num       );
    tree.Branch("gen_lep1_id"          , &gen_lep1_id          );
    tree.Branch("gen_lep1_charge"      , &gen_lep1_charge      );
    tree.Branch("gen_lep2_id"          , &gen_lep2_id          );
    tree.Branch("gen_lep2_charge"      , &gen_lep2_charge      );

    tree.Branch("gen_p4"     , "LorentzVector", &gen_p4     );
    tree.Branch("gen_lep1_p4", "LorentzVector", &gen_lep1_p4);
    tree.Branch("gen_lep2_p4", "LorentzVector", &gen_lep2_p4);
}

std::ostream& operator<< (std::ostream& out, const DrellYanInfo& info)
{
    out << "run                  = " << info.run                  << std::endl;
    out << "ls                   = " << info.ls                   << std::endl;
    out << "evt                  = " << info.evt                  << std::endl;
    out << "sample               = " << info.sample               << std::endl;
    out << "dataset              = " << info.dataset              << std::endl;
    out << "filename             = " << info.filename             << std::endl;
    out << "is_real_data         = " << info.is_real_data         << std::endl;
    out << "scale1fb             = " << info.scale1fb             << std::endl;
    out << "scale1fb_cms2        = " << info.scale1fb_cms2        << std::endl;
    out << "xsec                 = " << info.xsec                 << std::endl;
    out << "nevts_aod            = " << info.nevts_aod            << std::endl;
    out << "nevts_cms2           = " << info.nevts_cms2           << std::endl;
    out << "nevts_file           = " << info.nevts_file           << std::endl;
    out << "kfactor              = " << info.kfactor              << std::endl;
    out << "nvtxs                = " << info.nvtxs                << std::endl;
    out << "pfmet                = " << info.pfmet                << std::endl;
    out << "pfmet_phi            = " << info.pfmet_phi            << std::endl;
    out << "pu_nvtxs             = " << info.pu_nvtxs             << std::endl;
    out << "pu_ntrueint          = " << info.pu_ntrueint          << std::endl;
    out << "trig_dmu             = " << info.trig_dmu             << std::endl;
    out << "trig_del             = " << info.trig_del             << std::endl;
    out << "trig_smu             = " << info.trig_smu             << std::endl;
    out << "trig_sel             = " << info.trig_sel             << std::endl;
    out << "trig                 = " << info.trig                 << std::endl;
    out << "gen_hyp_type         = " << info.gen_hyp_type         << std::endl;
    out << "is_gen_ee            = " << info.is_gen_ee            << std::endl;
    out << "is_gen_mm            = " << info.is_gen_mm            << std::endl;
    out << "is_gen_tt            = " << info.is_gen_tt            << std::endl;
    out << "is_gen_ee_includetau = " << info.is_gen_ee_includetau << std::endl;
    out << "is_gen_mm_includetau = " << info.is_gen_mm_includetau << std::endl;
    out << "is_gen_fromz         = " << info.is_gen_fromz         << std::endl;
    out << "is_gen_acc_den       = " << info.is_gen_acc_den       << std::endl;
    out << "is_gen_acc_num       = " << info.is_gen_acc_num       << std::endl;
    out << "gen_p4.mass()        = " << info.gen_p4.mass()        << std::endl;
    out << "gen_lep1_p4.pt()     = " << info.gen_lep1_p4.pt()     << std::endl;
    out << "gen_lep1_id          = " << info.gen_lep1_id          << std::endl;
    out << "gen_lep1_charge      = " << info.gen_lep1_charge      << std::endl;
    out << "gen_lep2_p4.pt()     = " << info.gen_lep2_p4.pt()     << std::endl;
    out << "gen_lep2_id          = " << info.gen_lep2_id          << std::endl;
    out << "gen_lep2_charge      = " << info.gen_lep2_charge      << std::endl;
    return out;
}

// -------------------------------------------------//
// Simple class to hold your fill the DrellYanInfo 
// -------------------------------------------------//

class DrellYanNtupleMaker
{
    public:
        // construct:
        DrellYanNtupleMaker() = delete;
        DrellYanNtupleMaker
        (
            const dy::Sample::Info sample_info, 
            const std::string& output_filename, 
            const double lumi,
            const long num_events,
            const float min_pt,
            const float max_eta,
            const bool verbose
        );

        // destroy:
        ~DrellYanNtupleMaker();

        // basic methods:
        void BeginJob();
        void Analyze(const long event, const std::string& current_file);
        void EndJob();

    private:

        // members:
        dy::Sample::Info m_sample_info; 
        std::string m_output_filename;
        double m_lumi;
        long m_num_events;
        float m_min_pt;
        float m_max_eta;
        bool m_verbose;
        DrellYanInfo m_info;
        TFile& m_file;
        TTree& m_tree;
};

// construct:
DrellYanNtupleMaker::DrellYanNtupleMaker
(
    const dy::Sample::Info sample_info, 
    const std::string& output_filename, 
    const double lumi,
    const long num_events,
    const float min_pt,
    const float max_eta,
    const bool verbose
)
    : m_sample_info(sample_info)
    , m_output_filename(output_filename)
    , m_lumi(lumi)
    , m_num_events(num_events)
    , m_min_pt(min_pt)
    , m_max_eta(max_eta)
    , m_verbose(verbose)
    , m_info()
    , m_file(*TFile::Open(output_filename.c_str(), "RECREATE"))
    , m_tree(*new TTree("tree", "DY Exercise TTree"))
{
//     std::cout << gDirectory->GetName() << std::endl;
    // setup TTree
//     m_tree.SetDirectory(&m_file);
}

// destroy:
DrellYanNtupleMaker::~DrellYanNtupleMaker()
{
}

// ------------------------------------ //
// Stuff to do before job starts
// ------------------------------------ //

void DrellYanNtupleMaker::BeginJob()
{
    m_info.SetBranches(m_tree);
}

// ------------------------------------ //
// Stuff to do on each event 
// ------------------------------------ //

void DrellYanNtupleMaker::Analyze(const long event, const std::string& current_filename)
{
    if (m_verbose)
    {
        std::cout << "[DrellYanNtupleMaker] Running on run, ls, event: " 
            << tas::evt_run()       << ", "
            << tas::evt_lumiBlock() << ", "
            << tas::evt_event()     << ", "
            << std::endl;
    }

    // reset the TTree variables
    // ---------------------- // 
    m_info.Reset();

    // event information 
    // ---------------------- // 

    m_info.run          = tas::evt_run();
    m_info.ls           = tas::evt_lumiBlock();
    m_info.evt          = tas::evt_event();
    m_info.sample       = m_sample_info.sample; 
    m_info.dataset      = tas::evt_dataset().front().Data();
    m_info.filename     = current_filename; 
    m_info.is_real_data = tas::evt_isRealData(); 
    if (!tas::evt_isRealData())
    {
        m_info.nevts_aod     = tas::evt_nEvts();
        m_info.nevts_cms2    = m_sample_info.filter_eff * tas::evt_nEvts();
        m_info.nevts_file    = m_num_events;
        m_info.scale1fb_cms2 = tas::evt_scale1fb();
        m_info.scale1fb      = m_info.scale1fb_cms2 * m_info.nevts_cms2/m_info.nevts_file;
        m_info.xsec          = tas::evt_xsec_excl();
    }
    if (m_verbose)
    {
        std::cout << "\n[DrellYanNtupleMaker] Running on run, ls, event: " 
            << tas::evt_run()       << ", "
            << tas::evt_lumiBlock() << ", "
            << tas::evt_event()     << std::endl;
    }

    // gen information 
    // ---------------------- // 

    if (!tas::evt_isRealData())
    {
        // gen gen hypotheses
        const std::vector<at::GenHyp> gen_hyps       = at::GetGenHyps(/*min_pt=*/0.0, /*max_eta=*/1000.0);
        const std::vector<at::GenHyp> gen_hyps_clean = lt::filter_container(gen_hyps,
            [](const at::GenHyp& h)
            {
                return (h.IsOS() and (h.IsEE_IncludeTaus() or h.IsMuMu_IncludeTaus()));
            }
        );

        if (m_verbose) {std::cout << "number of gen hyps = " << gen_hyps_clean.size() << std::endl;}
        if (!gen_hyps_clean.empty())
        {
            const at::GenHyp& gen_hyp = gen_hyps_clean.front();

            // gen info:
            m_info.gen_hyp_type         = static_cast<int>(gen_hyp.Type());
            m_info.is_gen_ee            = gen_hyp.IsEE();
            m_info.is_gen_mm            = gen_hyp.IsMuMu();
            m_info.is_gen_tt            = gen_hyp.IsTauTau();
            m_info.is_gen_ee_includetau = gen_hyp.IsEE_IncludeTaus();
            m_info.is_gen_mm_includetau = gen_hyp.IsMuMu_IncludeTaus();
            m_info.is_gen_fromz         = gen_hyp.IsFromZ();
            m_info.is_gen_acc_den       = (gen_hyp.IsFromZ() and not gen_hyp.IsTauTau() and (60 < gen_hyp.P4().mass() && gen_hyp.P4().mass() < 120));
            m_info.is_gen_acc_num       = (m_info.is_gen_acc_den and gen_hyp.IsAccepted(m_min_pt, m_max_eta));
            m_info.gen_p4               = gen_hyp.P4();
            m_info.gen_lep1_p4          = gen_hyp.Lep1().p4;
            m_info.gen_lep1_id          = gen_hyp.Lep1().id;
            m_info.gen_lep1_charge      = gen_hyp.Lep1().charge;
            m_info.gen_lep2_p4          = gen_hyp.Lep2().p4;
            m_info.gen_lep2_id          = gen_hyp.Lep2().id;
            m_info.gen_lep2_charge      = gen_hyp.Lep2().charge;
        }
    }

    // fill the tree
    // ---------------------- // 
    if (m_verbose) {std::cout << m_info << std::endl;}

    m_tree.Fill();

    // done
    return;
}

// ------------------------------------ //
// Stuff to do after job finishes
// ------------------------------------ //

void DrellYanNtupleMaker::EndJob()
{
    std::cout << "[DrellYanNtupleMaker] Saving TTree to output file: " << m_output_filename << std::endl;
    m_file.cd(); 
    m_tree.Write();
    m_file.Close(); 
    return;
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
    std::cout << "[dy_create_ntuple] inputs:\n";
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
        throw std::invalid_argument(Form("[dy_create_ntuple] pset %s does not exist.\n", sample_pset.c_str()));
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
        output_file = Form("babies/%s/%s_baby.root", label.c_str(), sample_info.name.c_str());
    }
    lt::mkdir(lt::dirname(output_file), /*force=*/true);

    // run the looper
    // -------------------------------------------------------------------------------------------------//

    std::cout << "[dy_create_ntuple] TChain set to run on:\n";
    rt::PrintFilesFromTChain(chain);

    // create looper
    const long long num_events_to_run = (number_of_events < 0 ? chain->GetEntries() : std::min(number_of_events, chain->GetEntries()));
    DrellYanNtupleMaker looper
    (
        sample_info,
        output_file, 
        lumi,
        num_events_to_run,
        /*min_pt=*/25.0,
        /*max_eta=*/2.5,
        verbose
    );
    std::cout << "[dy_create_ntuple] running drell-yan plotting looper...\n";

    // scan the chain
    at::ScanChainWithFilename<CMS2>
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
    std::cerr << "[dy_create_ntuple] Error: failed..." << std::endl;
    std::cerr << e.what() << std::endl;
    return 1;
}
