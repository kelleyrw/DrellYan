#include "TChain.h"
#include "TFile.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TH1.h"
#include <iostream>
#include <string>
#include <cassert>

void dy_scale1fb_example(const double lumi = 0.082/*fb^-1*/)
{
    TChain chain("Events");
    chain.Add("/nfs-7/userdata/rwkelley/cms2/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root");

    // -------------------------------------- //
    // straight up
    // -------------------------------------- //

    // get the raw number of entries for true e^+e^- events
    const TCut selection   = "Sum$(genps_status==3 && genps_id==11)>=1 && Sum$(genps_status==3 && genps_id==-11)>=1";
    const double raw_count = chain.GetEntries(selection);

    // information from twiki: http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/Summer12MonteCarlo53X_Slim_Winter13#Drell_Yan
    const double xsec       = 3532.81; // pb
    const double kfactor    = 1.0;     // kinematically dependent factor to scale to a higher order calc -- not really used as much
    const double filt_eff   = 1.0;     // Filter efficiency for the generator level (done during MC generation at the CMSSW level)
    const double nevts_file = chain.GetEntries();

    // scale1fb = the factor the weight each event to get scale # events / fb^-1
    const double scale1fb = (xsec * kfactor * 1000.0f * filt_eff)/nevts_file;

    // output
    std::cout << "raw_count    = " << raw_count << "\n";
    std::cout << "scale1fb     = " << scale1fb << "\n";
    std::cout << "scaled count = lumi * scale1fb * raw_count = " << lumi*scale1fb*raw_count << "\n\n";

    // -------------------------------------- //
    // slightly wrong -- we did not account
    // for the events filtered by the CMSSW job's 
    // EDFilter
    // -------------------------------------- //

    // information from twiki: http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/Summer12MonteCarlo53X_Slim_Winter13#Drell_Yan
    const double nevts_aod       = 30459503; 
    const double nevts_cms2      = 27137253;
    const double cms2_filter_eff = nevts_cms2/nevts_aod;
    std::cout << "properly scale count = lumi * scale1fb * raw_count * cms2_filter_eff = " << lumi*scale1fb*raw_count*cms2_filter_eff << "\n\n";

    // -------------------------------------- //
    // Using info from the Tree 
    // -------------------------------------- //

    // scale1fb correction since we used a subset of the events
    const double scale = nevts_cms2/nevts_file; 

    // get the scaled entries
    const TCut scaled_selection  = Form("%f*%f*evt_scale1fb*(Sum$(genps_status==3 && genps_id==11)>=1 && Sum$(genps_status==3 && genps_id==-11)>=1)", lumi, scale);
    TH1D* h_scale1fb_count const = new TH1D("h_scale1fb_count", "h_scale1fb_count", 3, -0.5, 2.5);
    chain.Draw("1>>h_scale1fb_count", scaled_selection, "goff");

    std::cout << "hist entries    = " << h_scale1fb_count->GetEntries()   << "\n";
    std::cout << "scale1fb        = " << chain.GetMaximum("evt_scale1fb") << "\n";
    std::cout << "applying weight = " << scaled_selection.GetTitle()      << "\n";
    std::cout << "scale1fb_count  = " << h_scale1fb_count->Integral()     << "\n";
}

