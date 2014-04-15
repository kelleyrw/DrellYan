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
    float lumi;
    float xsec;
    float nevts_aod;
    float nevts_cms2;
    float nevts_file;
    float kfactor;
    float filt_eff;

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

    // selections
    bool passes_ossf;
    bool passes_mwin;
    bool passes_svtx;
    bool passes_trig;
    bool passes_idiso;
    bool passes_full;

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

    // reco hyp variables
    int hyp_type;
    LorentzVector hyp_p4;
    bool is_ee;
    bool is_mm;
    bool is_ll;

    // reco lepton variables
    LorentzVector lep1_p4;
    int lep1_id;
    int lep1_charge;
    double lep1_iso;
    double lep1_d0;
    double lep1_dz;
    bool lep1_passes_id;
    bool lep1_passes_iso;
    bool lep1_passes_sel;

    LorentzVector lep2_p4;
    int lep2_id;
    int lep2_charge;
    double lep2_iso;
    double lep2_d0;
    double lep2_dz;
    bool lep2_passes_id;
    bool lep2_passes_iso;
    bool lep2_passes_sel;
};

DrellYanInfo::DrellYanInfo()
    : run                  ( -999999    ) 
    , ls                   ( -999999    ) 
    , evt                  ( -999999    ) 
    , sample               ( -999999    ) 
    , dataset              ( ""         ) 
    , filename             ( ""         ) 
    , is_real_data         ( false      ) 
    , scale1fb             ( 1.0        ) 
    , scale1fb_cms2        ( 1.0        ) 
    , lumi                 ( 1.0        ) 
    , xsec                 ( -999999    ) 
    , nevts_aod            ( -999999    ) 
    , nevts_cms2           ( -999999    ) 
    , nevts_file           ( -999999    ) 
    , kfactor              ( -999999    ) 
    , filt_eff             ( -999999    ) 
    , gen_hyp_type         ( false      ) 
    , is_gen_ee            ( false      ) 
    , is_gen_mm            ( false      ) 
    , is_gen_tt            ( false      ) 
    , is_gen_ee_includetau ( false      ) 
    , is_gen_mm_includetau ( false      ) 
    , is_gen_fromz         ( false      ) 
    , is_gen_acc_den       ( false      ) 
    , is_gen_acc_num       ( false      ) 
    , gen_p4               ( 0, 0, 0, 0 ) 
    , gen_lep1_p4          ( 0, 0, 0, 0 ) 
    , gen_lep1_id          ( -999999    ) 
    , gen_lep1_charge      ( -999999    ) 
    , gen_lep2_p4          ( 0, 0, 0, 0 ) 
    , gen_lep2_id          ( -999999    ) 
    , gen_lep2_charge      ( -999999    ) 
    , passes_ossf          ( false      )
    , passes_mwin          ( false      )
    , passes_svtx          ( false      )
    , passes_trig          ( false      )
    , passes_idiso         ( false      )
    , passes_full          ( false      )
    , nvtxs                ( -999999    ) 
    , pfmet                ( -999999    ) 
    , pfmet_phi            ( -999999    ) 
    , pu_nvtxs             ( -999999    ) 
    , pu_ntrueint          ( -999999    ) 
    , trig_dmu             ( false      ) 
    , trig_del             ( false      ) 
    , trig_smu             ( false      ) 
    , trig_sel             ( false      ) 
    , trig                 ( false      ) 
    , hyp_type             ( -99999     ) 
    , hyp_p4               ( 0, 0, 0, 0 ) 
    , is_ee                ( false      ) 
    , is_mm                ( false      ) 
    , is_ll                ( false      ) 
    , lep1_p4              ( 0, 0, 0, 0 ) 
    , lep1_id              ( -999999    ) 
    , lep1_charge          ( -999999    ) 
    , lep1_iso             ( -999999    ) 
    , lep1_d0              ( -999999    ) 
    , lep1_dz              ( -999999    ) 
    , lep1_passes_id       ( false      ) 
    , lep1_passes_iso      ( false      ) 
    , lep1_passes_sel      ( false      ) 
    , lep2_p4              ( 0, 0, 0, 0 ) 
    , lep2_id              ( -999999    ) 
    , lep2_charge          ( -999999    ) 
    , lep2_iso             ( -999999    ) 
    , lep2_d0              ( -999999    ) 
    , lep2_dz              ( -999999    ) 
    , lep2_passes_id       ( false      ) 
    , lep2_passes_iso      ( false      ) 
    , lep2_passes_sel      ( false      ) 
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
    scale1fb             = 1.0; 
    scale1fb_cms2        = 1.0; 
    lumi                 = 1.0; 
    xsec                 = -999999; 
    nevts_aod            = -999999; 
    nevts_cms2           = -999999; 
    nevts_file           = -999999; 
    kfactor              = -999999; 
    filt_eff             = -999999; 
    gen_hyp_type         = false; 
    is_gen_ee            = false; 
    is_gen_mm            = false; 
    is_gen_tt            = false; 
    is_gen_ee_includetau = false; 
    is_gen_mm_includetau = false; 
    is_gen_fromz         = false; 
    is_gen_acc_den       = false; 
    is_gen_acc_num       = false; 
    gen_p4               = LorentzVector(0, 0, 0, 0); 
    gen_lep1_p4          = LorentzVector(0, 0, 0, 0); 
    gen_lep1_id          = -999999; 
    gen_lep1_charge      = -999999; 
    gen_lep2_p4          = LorentzVector(0, 0, 0, 0); 
    gen_lep2_id          = -999999; 
    gen_lep2_charge      = -999999; 
    passes_ossf          = false;
    passes_mwin          = false;
    passes_svtx          = false;
    passes_trig          = false;
    passes_idiso         = false;
    passes_full          = false;
    nvtxs                = -999999; 
    pfmet                = -999999; 
    pfmet_phi            = -999999; 
    pu_nvtxs             = -999999; 
    pu_ntrueint          = -999999; 
    trig_dmu             = false; 
    trig_del             = false; 
    trig_smu             = false; 
    trig_sel             = false; 
    trig                 = false; 
    hyp_type             = -99999;
    hyp_p4               = LorentzVector(0, 0, 0, 0 );
    is_ee                = false;
    is_mm                = false; 
    is_ll                = false; 
    lep1_p4              = LorentzVector(0, 0, 0, 0); 
    lep1_id              = -999999;
    lep1_charge          = -999999;
    lep1_iso             = -999999;
    lep1_d0              = -999999;
    lep1_dz              = -999999;
    lep1_passes_id       = false; 
    lep1_passes_iso      = false; 
    lep1_passes_sel      = false; 
    lep2_p4              = LorentzVector(0, 0, 0, 0); 
    lep2_id              = -999999;
    lep2_charge          = -999999;
    lep2_iso             = -999999;
    lep2_d0              = -999999;
    lep2_dz              = -999999;
    lep2_passes_id       = false;
    lep2_passes_iso      = false;
    lep2_passes_sel      = false;
}
    
void DrellYanInfo::SetBranches(TTree& tree)
{
    tree.Branch("run"                  , &run                  );
    tree.Branch("ls"                   , &ls                   );
    tree.Branch("evt"                  , &evt                  );
    tree.Branch("sample"               , &sample               );
    tree.Branch("dataset"              , &dataset              );
    tree.Branch("filename"             , &filename             );
    tree.Branch("is_real_data"         , &is_real_data         );
    tree.Branch("scale1fb"             , &scale1fb             );
    tree.Branch("scale1fb_cms2"        , &scale1fb_cms2        );
    tree.Branch("lumi"                 , &lumi                 );
    tree.Branch("xsec"                 , &xsec                 );
    tree.Branch("nevts_aod"            , &nevts_aod            );
    tree.Branch("nevts_cms2"           , &nevts_cms2           );
    tree.Branch("nevts_file"           , &nevts_file           );
    tree.Branch("filt_eff"             , &filt_eff             );
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
    tree.Branch("passes_ossf"          , &passes_ossf          );
    tree.Branch("passes_mwin"          , &passes_mwin          );
    tree.Branch("passes_svtx"          , &passes_svtx          );
    tree.Branch("passes_trig"          , &passes_trig          );
    tree.Branch("passes_idiso"         , &passes_idiso         );
    tree.Branch("passes_full"          , &passes_full          );
    tree.Branch("nvtxs"                , &nvtxs                );
    tree.Branch("pfmet"                , &pfmet                );
    tree.Branch("pfmet_phi"            , &pfmet_phi            );
    tree.Branch("pu_nvtxs"             , &pu_nvtxs             );
    tree.Branch("pu_ntrueint"          , &pu_ntrueint          );
    tree.Branch("trig_dmu"             , &trig_dmu             );
    tree.Branch("trig_del"             , &trig_del             );
    tree.Branch("trig_smu"             , &trig_smu             );
    tree.Branch("trig_sel"             , &trig_sel             );
    tree.Branch("trig"                 , &trig                 );
    tree.Branch("hyp_type"             , &hyp_type             );
    tree.Branch("is_ee"                , &is_ee                );
    tree.Branch("is_mm"                , &is_mm                );
    tree.Branch("is_ll"                , &is_ll                );
    tree.Branch("lep1_id"              , &lep1_id              );
    tree.Branch("lep1_charge"          , &lep1_charge          );
    tree.Branch("lep1_iso"             , &lep1_iso             );
    tree.Branch("lep1_d0"              , &lep1_d0              );
    tree.Branch("lep1_dz"              , &lep1_dz              );
    tree.Branch("lep1_passes_id"       , &lep1_passes_id       );
    tree.Branch("lep1_passes_iso"      , &lep1_passes_iso      );
    tree.Branch("lep1_passes_sel"      , &lep1_passes_sel      );
    tree.Branch("lep2_id"              , &lep2_id              );
    tree.Branch("lep2_charge"          , &lep2_charge          );
    tree.Branch("lep2_iso"             , &lep2_iso             );
    tree.Branch("lep2_d0"              , &lep2_d0              );
    tree.Branch("lep2_dz"              , &lep2_dz              );
    tree.Branch("lep2_passes_id"       , &lep2_passes_id       );
    tree.Branch("lep2_passes_iso"      , &lep2_passes_iso      );
    tree.Branch("lep2_passes_sel"      , &lep2_passes_sel      );

    tree.Branch("hyp_p4"      , "LorentzVector" , &hyp_p4     );
    tree.Branch("lep1_p4"     , "LorentzVector" , &lep1_p4    );
    tree.Branch("lep2_p4"     , "LorentzVector" , &lep2_p4    );
    tree.Branch("gen_p4"      , "LorentzVector" , &gen_p4     );
    tree.Branch("gen_lep1_p4" , "LorentzVector" , &gen_lep1_p4);
    tree.Branch("gen_lep2_p4" , "LorentzVector" , &gen_lep2_p4);
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
    out << "lumi                 = " << info.lumi                 << std::endl;
    out << "xsec                 = " << info.xsec                 << std::endl;
    out << "nevts_aod            = " << info.nevts_aod            << std::endl;
    out << "nevts_cms2           = " << info.nevts_cms2           << std::endl;
    out << "nevts_file           = " << info.nevts_file           << std::endl;
    out << "kfactor              = " << info.kfactor              << std::endl;
    out << "filt_eff             = " << info.filt_eff             << std::endl;
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
    out << "passes_ossf          = " << info.passes_ossf          << std::endl;
    out << "passes_mwin          = " << info.passes_mwin          << std::endl;
    out << "passes_svtx          = " << info.passes_svtx          << std::endl;
    out << "passes_trig          = " << info.passes_trig          << std::endl;
    out << "passes_idiso         = " << info.passes_idiso         << std::endl;
    out << "passes_full          = " << info.passes_full          << std::endl;
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
    out << "hyp_type             = " << info.hyp_type             << std::endl;
    out << "hyp_p4.mass()        = " << info.hyp_p4.mass()        << std::endl;
    out << "is_ee                = " << info.is_ee                << std::endl;
    out << "is_mm                = " << info.is_mm                << std::endl;
    out << "is_ll                = " << info.is_ll                << std::endl;
    out << "lep1_p4.pt()         = " << info.lep1_p4.pt()         << std::endl;
    out << "lep1_id              = " << info.lep1_id              << std::endl;
    out << "lep1_charge          = " << info.lep1_charge          << std::endl;
    out << "lep1_iso             = " << info.lep1_iso             << std::endl;
    out << "lep1_d0              = " << info.lep1_d0              << std::endl;
    out << "lep1_dz              = " << info.lep1_dz              << std::endl;
    out << "lep1_passes_id       = " << info.lep1_passes_id       << std::endl;
    out << "lep1_passes_iso      = " << info.lep1_passes_iso      << std::endl;
    out << "lep1_passes_sel      = " << info.lep1_passes_sel      << std::endl;
    out << "lep2_p4.pt()         = " << info.lep2_p4.pt()         << std::endl;
    out << "lep2_id              = " << info.lep2_id              << std::endl;
    out << "lep2_charge          = " << info.lep2_charge          << std::endl;
    out << "lep2_iso             = " << info.lep2_iso             << std::endl;
    out << "lep2_d0              = " << info.lep2_d0              << std::endl;
    out << "lep2_dz              = " << info.lep2_dz              << std::endl;
    out << "lep2_passes_id       = " << info.lep2_passes_id       << std::endl;
    out << "lep2_passes_iso      = " << info.lep2_passes_iso      << std::endl;
    out << "lep2_passes_sel      = " << info.lep2_passes_sel      << std::endl;
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
}

// destroy:
DrellYanNtupleMaker::~DrellYanNtupleMaker()
{
}

// ------------------------------------ //
// helper functions and classes 
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
        bool passes;
        int hyp_idx;
    };
};

std::ostream& operator<<(std::ostream& out, Selection::Info si)
{
    out << "{"  << si.name
        << ", " << si.title
        << ", " << si.passes
        << ", " << si.hyp_idx
        << "}";
    return out;
}

void UpdateSelection(Selection::Info& sel, const int hyp_idx)
{
    sel.passes  = (hyp_idx >= 0); 
    sel.hyp_idx = dy::ChooseBetterHypothesis(hyp_idx, sel.hyp_idx); 
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
        std::cout << "\n[DrellYanNtupleMaker] Running on run, ls, event: " 
            << tas::evt_run()       << ", "
            << tas::evt_lumiBlock() << ", "
            << tas::evt_event()     << std::endl;
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
        m_info.lumi          = m_lumi; 
        m_info.xsec          = tas::evt_xsec_excl();
        m_info.kfactor       = tas::evt_kfactor();
        m_info.filt_eff      = tas::evt_filt_eff();
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

            // for dytt, only require tau tau events
            if (m_sample_info.sample == dy::Sample::dytt and not gen_hyp.IsTauTau())
            {
                if (m_verbose) {std::cout << "fails dy --> tau tau requirement" << std::endl;}
                return;
            }

            // gen info:
            if (tas::puInfo_nPUvertices().size()==3)
            {
                m_info.pu_nvtxs    = tas::puInfo_nPUvertices().at(1);
                m_info.pu_ntrueint = tas::puInfo_trueNumInteractions().at(1);
            }
            else
            {   
                m_info.pu_nvtxs    = tas::puInfo_nPUvertices().at(0);
                m_info.pu_ntrueint = tas::puInfo_trueNumInteractions().at(0);
            }
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

    // reco information 
    // ---------------------- // 

    // selections
    std::vector<Selection::Info> selections;
    selections.push_back(Selection::Info{"ossf"  , "OSSF"                , false , -1});
    selections.push_back(Selection::Info{"mwin"  , "Mass Window"         , false , -1});
    selections.push_back(Selection::Info{"svtx"  , "Same Vertex"         , false , -1});
    selections.push_back(Selection::Info{"trig"  , "Passes Trigger"      , false , -1});
    selections.push_back(Selection::Info{"idiso" , "Passes ID/Isolation" , false , -1});
    selections.push_back(Selection::Info{"full"  , "Full Selection"      , false , -1});
    assert(selections.size()==Selection::static_size);

    // loop over hypotheses
    if (m_verbose) {std::cout << "looping over " << tas::hyp_type().size() << " hypotheses" << std::endl;}

    for (size_t hyp_idx = 0; hyp_idx < tas::hyp_type().size(); ++hyp_idx)
    {                
        if (m_verbose)
        {
            std::cout << "hyp " << hyp_idx << " of " << tas::hyp_type().size() << std::endl;
        }

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
        else if (not(is_ee or is_mm))
        {
            if (m_verbose) {std::cout << "not SF" << std::endl;}
            continue;
        }
        UpdateSelection(selections[Selection::ossf], hyp_idx);
    
        // 60 < m_ll << 120 GeV
        if (not (60 < dilep_mass && dilep_mass < 120.0))        
        {
            if (m_verbose) {std::cout << "not SF" << std::endl;}
            continue;
        }
        UpdateSelection(selections[Selection::mwin], hyp_idx);
    
        // both leptons from first vertex
        if (not hypsFromFirstGoodVertex(hyp_idx))
        {
            if (m_verbose) {std::cout << "did not pass same vertex" << std::endl;}
            continue;
        }
        UpdateSelection(selections[Selection::svtx], hyp_idx);

        // trigger
        if (not dy::passesTrigger(flavor_type))
        {
            if (m_verbose) {std::cout << "did not pass trigger" << std::endl;}
            continue;
        }
        UpdateSelection(selections[Selection::trig], hyp_idx);

        // l1 and l2 pass selection
        if (not dy::isSelectedHypothesis(hyp_idx))
        {
            if (m_verbose) {std::cout << "not selected" << std::endl;}
            continue;
        }
        UpdateSelection(selections[Selection::idiso], hyp_idx);
    
        // if here, fully selected
        UpdateSelection(selections[Selection::full], hyp_idx);

    } // end loop over hypothesis

    if (m_verbose)
    {
        for (const auto& selection : selections)
        {
            std::cout << "selection " << selection.name << " = " 
                << selection.passes << "\t" 
                << selection.hyp_idx << "\n";
        }
    }


    // only continue if a hyp has been selected
    // the order of the selction matters
    // all: 0, mm: 1, em: 2, ee: 3
    int best_hyp = -1;
    for (const auto& s : selections)
    {
        best_hyp = (s.passes ? s.hyp_idx : best_hyp);
    }

    // require at least 3 tracks in the event
    const bool clean_tracks   = (tas::trks_trk_p4().size() >= 3);
    if (not clean_tracks)
    {
        if (m_verbose) {std::cout << "fails # trks >= 3 requirement" << std::endl;}
    }

    // require standard cleaning 
    const bool clean_standard = cleaning_standardNovember2011(); 
    if (not clean_standard)
    {
        if (m_verbose) {std::cout << "fails November2011 cleaning requirement" << std::endl;}
    }

    const int hyp_idx = best_hyp;
    if (not clean_tracks or not clean_standard or hyp_idx < 0)
    {
        if (m_verbose) {std::cout << "failed cleaning or no good hypthesis chosen" << std::endl;}
    }
    else
    {
        if (m_verbose) {std::cout << "best hypthesis chosen = " << hyp_idx << std::endl;}

        // selections:
        m_info.passes_ossf  = selections[Selection::ossf ].passes;
        m_info.passes_mwin  = selections[Selection::mwin ].passes;
        m_info.passes_svtx  = selections[Selection::svtx ].passes;
        m_info.passes_trig  = selections[Selection::trig ].passes;
        m_info.passes_idiso = selections[Selection::idiso].passes;
        m_info.passes_full  = selections[Selection::full ].passes;

        // reco hyp variables
        m_info.pfmet     = tas::evt_pfmet_type1cor();
        m_info.pfmet_phi = tas::evt_pfmetPhi_type1cor();
        m_info.nvtxs     = numberOfGoodVertices();

        // reco hyp variables
        m_info.hyp_type = tas::hyp_type().at(hyp_idx);
        m_info.hyp_p4   = tas::hyp_p4().at(hyp_idx);
        m_info.is_ee    = (m_info.hyp_type==3);
        m_info.is_mm    = (m_info.hyp_type==0);
        m_info.is_ll    = (m_info.is_ee or m_info.is_mm);

        // trigger
        m_info.trig_dmu = (m_info.is_mm && dy::passesTriggerDoubleLep(m_info.hyp_type));
        m_info.trig_del = (m_info.is_ee && dy::passesTriggerDoubleLep(m_info.hyp_type));
        m_info.trig_smu = (m_info.is_mm && dy::passesTriggerSingleLep(m_info.hyp_type));
        m_info.trig_sel = (m_info.is_ee && dy::passesTriggerSingleLep(m_info.hyp_type));
        m_info.trig     = dy::passesTrigger(m_info.hyp_type);

        // reco lepton variables
        // NOTE: set lepton info (lep1 is higher pT lepton, lep2 is lower pT lepton)
        LorentzVector lep1_p4;
        LorentzVector lep2_p4;
        int lep1_id;
        int lep1_idx;
        int lep2_id;
        int lep2_idx;
        if (tas::hyp_lt_p4().at(hyp_idx).pt() > tas::hyp_ll_p4().at(hyp_idx).pt())
        {
            lep1_p4  = cms2.hyp_lt_p4().at(hyp_idx);
            lep1_id  = cms2.hyp_lt_id().at(hyp_idx);
            lep1_idx = cms2.hyp_lt_index().at(hyp_idx); 
            lep2_p4  = cms2.hyp_ll_p4().at(hyp_idx);
            lep2_id  = cms2.hyp_ll_id().at(hyp_idx);    
            lep2_idx = cms2.hyp_ll_index().at(hyp_idx); 
        }
        else
        {
            lep1_p4  = cms2.hyp_ll_p4().at(hyp_idx);
            lep1_id  = cms2.hyp_ll_id().at(hyp_idx);
            lep1_idx = cms2.hyp_ll_index().at(hyp_idx); 
            lep2_p4  = cms2.hyp_lt_p4().at(hyp_idx);
            lep2_id  = cms2.hyp_lt_id().at(hyp_idx);    
            lep2_idx = cms2.hyp_lt_index().at(hyp_idx); 
        }
        m_info.lep1_p4         = lep1_p4; 
        m_info.lep1_id         = lep1_id;
        m_info.lep1_charge     = -1*lep1_id/abs(lep1_id);
        m_info.lep1_iso        = dy::leptonIsolation(lep1_id, lep1_idx);
        m_info.lep1_d0         = dy::leptonD0(lep1_id, lep1_idx);
        m_info.lep1_dz         = dy::leptonDz(lep1_id, lep1_idx);
        m_info.lep1_passes_id  = dy::isGoodLepton(lep1_id, lep1_idx);
        m_info.lep1_passes_iso = dy::isIsolatedLepton(lep1_id, lep1_idx);
        m_info.lep1_passes_sel = m_info.lep1_passes_id && m_info.lep1_passes_iso;

        m_info.lep2_p4         = lep2_p4; 
        m_info.lep2_id         = lep2_id;
        m_info.lep2_charge     = -1*lep2_id/abs(lep2_id);
        m_info.lep2_iso        = dy::leptonIsolation(lep2_id, lep2_idx);
        m_info.lep2_d0         = dy::leptonD0(lep2_id, lep2_idx);
        m_info.lep2_dz         = dy::leptonDz(lep2_id, lep2_idx);
        m_info.lep2_passes_id  = dy::isGoodLepton(lep2_id, lep2_idx);
        m_info.lep2_passes_iso = dy::isIsolatedLepton(lep2_id, lep2_idx);
        m_info.lep2_passes_sel = m_info.lep2_passes_id && m_info.lep2_passes_iso;
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
        output_file = Form("babies/%s_baby.root", sample_info.name.c_str());
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
