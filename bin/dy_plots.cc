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

// tools 
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
    hc.Add(new TH1D("h_reco_yield", "Yield count of reco level e^{+}e^{-} (full sample)", 3, 0, 3));
    return;
}

// ------------------------------------ //
// Stuff to do on each event 
// ------------------------------------ //
void DrellYanLooper::Analyze(const long event)
{
    // loop over hypotheses
    bool selected = false;
    for (size_t hyp_idx = 0; hyp_idx != tas::hyp_type().size(); hyp_idx++)
    {                
        // convenience variables
        const int lt_id                                   = tas::hyp_lt_id().at(hyp_idx);
        const int ll_id                                   = tas::hyp_ll_id().at(hyp_idx);
        const int lt_idx                                  = tas::hyp_lt_index().at(hyp_idx);
        const int ll_idx                                  = tas::hyp_ll_index().at(hyp_idx);
        const float dilep_mass                            = tas::hyp_p4().at(hyp_idx).mass();
        const at::DileptonHypType::value_type flavor_type = at::hyp_typeToHypType(tas::hyp_type().at(hyp_idx));

        // apply selections
        if (tas::evt_isRealData() && !dy::passesTrigger(tas::hyp_type().at(hyp_idx)))                {continue;}
        if ((tas::hyp_lt_charge().at(hyp_idx) * tas::hyp_ll_charge().at(hyp_idx)) > 0)               {continue;}
        if (not(flavor_type == at::DileptonHypType::EE or flavor_type == at::DileptonHypType::MUMU)) {continue;}
        if (not (60 < dilep_mass && dilep_mass < 120.0))                                             {continue;}
        if (not dy::isSelectedLepton(lt_id, lt_idx))                                                 {continue;}
        if (not dy::isSelectedLepton(ll_id, ll_idx))                                                 {continue;}
        if (not hypsFromFirstGoodVertex(hyp_idx))                                                    {continue;}

        // if we're here, then good event :)
        selected = true; 
    } // end looper over hypothesis

    if (not selected) {return;}

    const double nevts_full     = tas::evt_nEvts();       // number of events run in CMSSW job to make ntuple
    const double nevts_file     = m_num_events;           // number of events in the current job
    const double nevts_scale    = nevts_full/nevts_file;  // scale up the weight to account fo lower stats
    const double  scale1fb      = tas::evt_scale1fb();
    const double  lumi          = 0.082;
    const double nevts_sdfilter = 27137253;            // number of events after the "SDFilter"
    const double sd_filter_eff  = nevts_sdfilter/nevts_full;
    const float scale_subset    = lumi * scale1fb * nevts_scale * sd_filter_eff;

    // fill hist
    hc["h_reco_yield"]->Fill(2, scale_subset);
}

// ------------------------------------ //
// Stuff to do after job finishes
// ------------------------------------ //
void DrellYanLooper::EndJob()
{
    // output counts
    std::cout << "ll yield = " << rt::Integral(hc["h_reco_yield"]) << std::endl;
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
//     chain.Add("/nfs-7/userdata/rwkelley/cms2/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
    chain.Add("skims/dy_skim.root");
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

