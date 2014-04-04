// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TStyle.h" 
#include "TDirectory.h"
#include "TFile.h"
#include "TSystem.h" 
#include "TROOT.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "TDatabasePDG.h"

// CMS2
#include "CMS2/NtupleMacrosHeader/interface/CMS2.h"

// CORE
#include "CMS2/NtupleMacrosCore/interface/mcSelections.h"

// tools 
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"

// -------------------------------------------------//
// Simple class to hold your analysis data
// -------------------------------------------------//

class DrellYanLooper
{
    public:
        // construct:
        DrellYanLooper();
        DrellYanLooper(const std::string& filename, const long num_events);

        // destroy:
        ~DrellYanLooper();

        // methods:
        void SetFilename(const std::string& filename);

        // basic methods:
        void BeginJob();
        void Analyze(const long event);
        void EndJob();

    private:

        // members:
        std::string m_output_filename;
        std::string m_runlist_filename;
        long m_num_events;
        rt::TH1Container hc;
};

// construct:
DrellYanLooper::DrellYanLooper()
    : m_output_filename("plots/dy_scale1fb_example.root")
    , m_num_events(0)
{
}

DrellYanLooper::DrellYanLooper(const std::string& filename, const long num_events)
    : m_output_filename(filename)
    , m_num_events(num_events)
{
}

// destroy:
DrellYanLooper::~DrellYanLooper()
{
}

// ------------------------------------ //
// Stuff to do before job starts
// ------------------------------------ //
void DrellYanLooper::BeginJob()
{
    hc.Add(new TH1D("h_ee_yield_full", "Yield count of gen level e^{+}e^{-} (full sample)"  , 3, -0.5, 2.5));
    hc.Add(new TH1D("h_ee_yield_sub" , "Yield count of gen level e^{+}e^{-} (subset sample)", 3, -0.5, 2.5));
    gDirectory->ls();
    return;
}

// ------------------------------------ //
// Stuff to do on each event 
// ------------------------------------ //
void DrellYanLooper::Analyze(const long event)
{
    // select e+e- events (NOTE: no test for Z)
    bool contains_eplus  = false;
    bool contains_eminus = false;
    bool is_ee           = false;
    for (size_t idx = 0; idx < tas::genps_id().size() && !is_ee; ++idx)
    {
        if (tas::genps_status().at(idx) != 3) {continue;}
        if (tas::genps_id().at(idx) ==  11  ) {contains_eminus = true;}
        if (tas::genps_id().at(idx) == -11  ) {contains_eplus  = true;}
        is_ee = (contains_eplus && contains_eminus);
    }
    if (!is_ee) {return;}

    // --------------------------------------------------------- //
    // scale factor to scale assuming running on full MC sample 
    // --------------------------------------------------------- //

    const double lumi = 0.082; //fb^-1
    const double scale_full = lumi * tas::evt_scale1fb();
    hc["h_ee_yield_full"]->Fill(1.0, scale_full);

    // --------------------------------------------------------- //
    // scale1fb correction since we used a subset of the events
    // --------------------------------------------------------- //

    // information from twiki: http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/Summer12MonteCarlo53X_Slim_Winter13#Drell_Yan
//     const double nevts_aod   = tas::evt_nEvts();       // number of events run in CMSSW job to make ntuple
    const double nevts_cms2  = 27137253;               // number of events after the CMS2 level filter
    const double nevts_file  = m_num_events;           // number of events in the current job
    const double nevts_scale = nevts_cms2/nevts_file;  // scale up the weight to account fo lower stats

    const float scale_subset = scale_full * nevts_scale;
    hc["h_ee_yield_sub"]->Fill(1.0, scale_subset);
}

// ------------------------------------ //
// Stuff to do after job finishes
// ------------------------------------ //
void DrellYanLooper::EndJob()
{
    // output counts
    std::cout << "raw  count   = " << hc["h_ee_yield_full"]->GetEntries() << std::endl;
    std::cout << "full count   = " << hc["h_ee_yield_full"]->Integral()   << std::endl;
    std::cout << "subset count = " << hc["h_ee_yield_sub" ]->Integral()   << std::endl;
    std::cout << std::endl;

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

int main()
{
    const long long number_of_events = -1;
    const std::string good_run_list = "json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_cms2.txt";
    const bool verbose = false;

    // get the chain
    at::LoadFWLite();
    TChain chain("Events");
    chain.Add("/nfs-7/userdata/rwkelley/cms2/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
//     chain.Add("skims/dy_skim.root");
    const long long num_events_to_run = (number_of_events < 0 ? chain.GetEntries() : std::min(number_of_events, chain.GetEntries()));

    DrellYanLooper looper("plots/dy_scale1fb_plots.root", num_events_to_run);
    std::cout << "running drell-yan looper..." << std::endl;

    // simple style
    rt::SetStyle();

    // scan the chain
    at::ScanChain<CMS2>
    (
        &chain, 
        looper,
        cms2,
        num_events_to_run,
        good_run_list,
        /*fast=*/true,
        verbose
    ); 

    return 0;
}

