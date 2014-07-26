// C++
#include <iostream>
#include <vector>
#include <functional>

// ROOT
#include "TChain.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/PtEtaPhiE4D.h"

// CMS2
#include "CMS2/NtupleMacrosHeader/interface/CMS2.h"

// CORE
#include "CMS2/NtupleMacrosCore/interface/mcSelections.h"
#include "CMS2/NtupleMacrosCore/interface/eventSelections.h"
#include "CMS2/NtupleMacrosCore/interface/trackSelections.h"
#include "CMS2/NtupleMacrosCore/interface/electronSelections.h"
#include "CMS2/NtupleMacrosCore/interface/muonSelections.h"
#include "CMS2/NtupleMacrosCore/interface/mcSelections.h"

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

// -------------------------------------------------//
// Class to hold the ntuple information 
// -------------------------------------------------//

class SingleLeptonInfo
{
    public:
        SingleLeptonInfo();

        void Reset();
        void SetBranches(TTree& tree);
        void SetAliases(TTree& tree) const;
        void FillCommon(const int id, const int idx);

    public:    

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

        // lepton reco info
        LorentzVector p4;
        bool passes_id;
        bool passes_iso;
        bool is_num;
        bool is_den;
        bool is_fo;
        bool is_mu;
        bool is_el;
        bool exists;
        int is_fromw;
        int charge;
        int pdgid;
        int type;
        float d0;
        float d0err;
        float dz;
        float dzerr;
        float ip3d;
        float ip3derr;
        float mt;
        float corpfiso;
        float pfiso;
        float chiso;
        float emiso;
        float nhiso;
        float corpfiso04;
        float pfiso04;
        float chiso04;
        float emiso04;
        float nhiso04;
        float cordetiso;
        float detiso;
        float trkiso;
        float ecaliso;
        float hcaliso;
        float cordetiso04;
        float detiso04;
        float trkiso04;
        float ecaliso04;
        float hcaliso04;
        float effarea;
        float effarea04;
        float dbeta;
        float dbeta04;
        float sf_lepeff;
        float sf_trig;

        // lepton gen info
        LorentzVector mcp4;
        LorentzVector mc3p4;
        LorentzVector mc_momp4;
        int mcid;
        int mc3id;
        int momid;
        int mc3_momid;

        // electron specific info
        LorentzVector gsf_p4;
        LorentzVector ctf_p4;
        LorentzVector sc_p4;
        bool q3agree;
        bool is_conv;
        int qsc;
        int qctf;
        int qgsf;
        int nmhits;
        int eleid_veto;
        int eleid_loose;
        int eleid_medium;
        int eleid_tight;
        bool is_eleid_veto;
        bool is_eleid_loose;
        bool is_eleid_medium;
        bool is_eleid_tight;
        float dphiin;
        float detain;
        float sieie;
        float hoe;
        float ooemoop;
        float conv_dist;
        float conv_dcot;

        // muon specific info
        LorentzVector gfit_p4;
        bool is_global;
        bool is_tracker;
        bool is_stamu;
        bool is_pfmu;
        bool is_loosemu;
        bool is_tightmu;
        int npixelhits;
        int nsihits;    
        int nsilayers;
        int nstahits;
        int nstations;
        float chi2;
        float ndof;
        float pterr;
        float ecal_vetodep;
        float hcal_vetodep;    

};

SingleLeptonInfo::SingleLeptonInfo ()
    : run             ( -999999    ) 
    , ls              ( -999999    ) 
    , evt             ( -999999    ) 
    , sample          ( -999999    ) 
    , dataset         ( ""         ) 
    , filename        ( ""         ) 
    , is_real_data    ( false      ) 
    , scale1fb        ( 1.0        ) 
    , scale1fb_cms2   ( 1.0        ) 
    , lumi            ( 1.0        ) 
    , xsec            ( -999999    ) 
    , nevts_aod       ( -999999    ) 
    , nevts_cms2      ( -999999    ) 
    , nevts_file      ( -999999    ) 
    , kfactor         ( -999999    ) 
    , filt_eff        ( -999999    ) 
    , p4              ( 0, 0, 0, 0 )
    , passes_id       ( false      )
    , passes_iso      ( false      )
    , is_num          ( false      )
    , is_den          ( false      )
    , is_fo           ( false      )
    , is_mu           ( false      )
    , is_el           ( false      )
    , exists          ( false      )
    , is_fromw        ( -999999    )
    , charge          ( -999999    )
    , pdgid           ( -999999    )
    , type            ( -999999    )
    , d0              ( -999999.   )
    , d0err           ( -999999.   )
    , dz              ( -999999.   )
    , dzerr           ( -999999.   )
    , mt              ( -999999.   )
    , corpfiso        ( -999999.   )
    , pfiso           ( -999999.   )
    , chiso           ( -999999.   )
    , emiso           ( -999999.   )
    , nhiso           ( -999999.   )
    , corpfiso04      ( -999999.   )
    , pfiso04         ( -999999.   )
    , chiso04         ( -999999.   )
    , emiso04         ( -999999.   )
    , nhiso04         ( -999999.   )
    , cordetiso       ( -999999.   )
    , detiso          ( -999999.   )
    , trkiso          ( -999999.   )
    , ecaliso         ( -999999.   )
    , hcaliso         ( -999999.   )
    , cordetiso04     ( -999999.   )
    , detiso04        ( -999999.   )
    , trkiso04        ( -999999.   )
    , ecaliso04       ( -999999.   )
    , hcaliso04       ( -999999.   )
    , effarea         ( -999999.   )
    , effarea04       ( -999999.   )
    , dbeta           ( -999999.   )
    , dbeta04         ( -999999.   )
    , sf_lepeff       ( -999999.   )
    , sf_trig         ( -999999.   )
    , mcp4            ( 0, 0, 0, 0 )
    , mc3p4           ( 0, 0, 0, 0 )
    , mc_momp4        ( 0, 0, 0, 0 )
    , mcid            ( -999999    )
    , mc3id           ( -999999    )
    , momid           ( -999999    )
    , mc3_momid       ( -999999    )
    , gsf_p4          ( 0, 0, 0, 0 )
    , ctf_p4          ( 0, 0, 0, 0 )
    , sc_p4           ( 0, 0, 0, 0 )
    , q3agree         ( false      )
    , is_conv         ( false      )
    , qsc             ( -999999    )
    , qctf            ( -999999    )
    , qgsf            ( -999999    )
    , nmhits          ( -999999    )
    , eleid_veto      ( -999999    )
    , eleid_loose     ( -999999    )
    , eleid_medium    ( -999999    )
    , eleid_tight     ( -999999    )
    , is_eleid_veto   ( false      )
    , is_eleid_loose  ( false      )
    , is_eleid_medium ( false      )
    , is_eleid_tight  ( false      )
    , dphiin          ( -999999.   )
    , detain          ( -999999.   )
    , sieie           ( -999999.   )
    , hoe             ( -999999.   )
    , ooemoop         ( -999999.   )
    , conv_dist       ( -999999.   )
    , conv_dcot       ( -999999.   )
    , gfit_p4         ( 0, 0, 0, 0 )
    , is_global       ( false      )
    , is_tracker      ( false      )
    , is_stamu        ( false      )
    , is_pfmu         ( false      )
    , is_loosemu      ( false      )
    , is_tightmu      ( false      )
    , npixelhits      ( -999999    )
    , nsihits         ( -999999    )
    , nsilayers       ( -999999    )
    , nstahits        ( -999999    )
    , nstations       ( -999999    )
    , chi2            ( -999999.   )
    , ndof            ( -999999.   )
    , pterr           ( -999999.   )
    , ecal_vetodep    ( -999999.   )
    , hcal_vetodep    ( -999999.   )
{
}

//********************************************************//
//
// NOTE: This function only fills information that is
//       expected to be common to all analyses.
//
// The following information will need to be filled
// by the user for his specific analysis:
//
//      bool passes_id;
//      bool passes_iso;
//      bool is_num;
//      bool is_den;
//      bool is_fo;
//      float corpfiso, cordetiso, corpfiso04, cordetiso04;
//      float effarea, effarea04;
//      float sf_lepeff; (t&p scale factor)
//      float sf_trig; (trigger scale factor)
//      bool is_conv;

//      
//********************************************************//
void SingleLeptonInfo::Reset()
{
    run             = -999999; 
    ls              = -999999; 
    evt             = -999999; 
    sample          = -999999; 
    dataset         = ""; 
    filename        = ""; 
    is_real_data    = false; 
    scale1fb        = 1.0; 
    scale1fb_cms2   = 1.0; 
    lumi            = 1.0; 
    xsec            = -999999; 
    nevts_aod       = -999999; 
    nevts_cms2      = -999999; 
    nevts_file      = -999999; 
    kfactor         = -999999; 
    filt_eff        = -999999; 
    passes_id       = false;
    passes_iso      = false;
    is_num          = false;
    is_den          = false;
    is_fo           = false;
    is_mu           = false;
    is_el           = false;
    exists          = false;
    is_fromw        = -999999;
    charge          = -999999;
    pdgid           = -999999;
    type            = -999999;
    d0              = -999999.;
    d0err           = -999999.;
    dz              = -999999.;
    dzerr           = -999999.;
    mt              = -999999.;
    corpfiso        = -999999.;
    pfiso           = -999999.;
    chiso           = -999999.;
    emiso           = -999999.;
    nhiso           = -999999.;
    corpfiso04      = -999999.;
    pfiso04         = -999999.;
    chiso04         = -999999.;
    emiso04         = -999999.;
    nhiso04         = -999999.;
    cordetiso       = -999999.;
    detiso          = -999999.;
    trkiso          = -999999.;
    ecaliso         = -999999.;
    hcaliso         = -999999.;
    cordetiso04     = -999999.;
    detiso04        = -999999.;
    trkiso04        = -999999.;
    ecaliso04       = -999999.;
    hcaliso04       = -999999.;
    effarea         = -999999.;
    effarea04       = -999999.;
    dbeta           = -999999.;
    dbeta04         = -999999.;
    sf_lepeff       = -999999.;
    sf_trig         = -999999.;
    mcid            = -999999;
    mc3id           = -999999;
    momid           = -999999;
    mc3_momid       = -999999;
    q3agree         = false;
    is_conv         = false;
    qsc             = -999999;
    qctf            = -999999;
    qgsf            = -999999;
    nmhits          = -999999;
    eleid_veto      = -999999;
    eleid_loose     = -999999;
    eleid_medium    = -999999;
    eleid_tight     = -999999;
    is_eleid_veto   = false;
    is_eleid_loose  = false;
    is_eleid_medium = false;
    is_eleid_tight  = false;
    dphiin          = -999999.;
    detain          = -999999.;
    sieie           = -999999.;
    hoe             = -999999.;
    ooemoop         = -999999.;
    conv_dist       = -999999.;
    conv_dcot       = -999999.;
    is_global       = false;
    is_tracker      = false;
    is_stamu        = false;
    is_pfmu         = false;
    is_loosemu      = false;
    is_tightmu      = false;
    npixelhits      = -999999;
    nsihits         = -999999;
    nsilayers       = -999999;
    nstahits        = -999999;
    nstations       = -999999;
    chi2            = -999999.;
    ndof            = -999999.;
    pterr           = -999999.;
    ecal_vetodep    = -999999.;
    hcal_vetodep    = -999999.;

    p4           = LorentzVector(0, 0, 0, 0);
    mcp4         = LorentzVector(0, 0, 0, 0);
    mc3p4        = LorentzVector(0, 0, 0, 0);
    mc_momp4     = LorentzVector(0, 0, 0, 0);
    gsf_p4       = LorentzVector(0, 0, 0, 0);
    ctf_p4       = LorentzVector(0, 0, 0, 0);
    sc_p4        = LorentzVector(0, 0, 0, 0);
    gfit_p4      = LorentzVector(0, 0, 0, 0);
}

void SingleLeptonInfo::FillCommon(const int id, const int idx)
{
    if (idx < 0) return;

    // electron information 
    // ---------------------- // 

    const int vtxidx = firstGoodVertex();
    is_fromw = not tas::evt_isRealData() ? leptonIsFromW(idx, id, true) : -999999;

    if (abs(id) == 11)
    {
        exists = true;
        is_el  = true;
        p4     = tas::els_p4().at(idx);
        charge = tas::els_charge().at(idx);
        pdgid  = charge * -11;
        type   = tas::els_type().at(idx);

        const int gsfidx = tas::els_gsftrkidx().at(idx);
        if (gsfidx >= 0 && vtxidx >= 0)
        {
            std::pair<float, float> cord0 = gsftrks_d0_pv(gsfidx, vtxidx);
            std::pair<float, float> cordz = gsftrks_dz_pv(gsfidx, vtxidx);
            d0 = cord0.first;
            dz = cordz.first;
            d0err = cord0.second;
            dzerr = cordz.second;
        }
        ip3d    = tas::els_ip3d().at(idx);;
        ip3derr = tas::els_ip3derr().at(idx);;

        if (!tas::evt_isRealData())
        {
            mcp4       = tas::els_mc_p4().at(idx);
            mc_momp4   = tas::els_mc_motherp4().at(idx);
            mcid       = tas::els_mc_id().at(idx);
            momid      = tas::els_mc_motherid().at(idx);
            mc3id      = tas::els_mc3_id().at(idx);
            mc3_momid  = tas::els_mc3_motherid().at(idx);
            const int mc3idx = tas::els_mc3idx().at(idx);
            if (mc3idx >= 0)
            {
                mc3p4 = tas::genps_p4().at(mc3idx);           
            }
        }

        gsf_p4 = tas::els_trk_p4().at(idx);

        const int ctfidx = tas::els_trkidx().at(idx);
        if (ctfidx >= 0)
        {
            ctf_p4 = tas::trks_trk_p4().at(ctfidx);
            qctf   = tas::trks_charge().at(ctfidx);
        }

        const float tmp_eta = tas::els_etaSC().at(idx);
        const float tmp_phi = tas::els_phiSC().at(idx);
        const float tmp_e   = tas::els_eSC().at(idx);
        const float tmp_pt  = tmp_e / cosh(tmp_eta);
        const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > tmp_scp4(tmp_pt, tmp_eta, tmp_phi, tmp_e);
        sc_p4 = tmp_scp4;

        q3agree         = !isChargeFlip3agree(idx);
        qgsf            = tas::els_trk_charge().at(idx);
        qsc             = tas::els_sccharge().at(idx);
        nmhits          = tas::els_exp_innerlayers().at(idx);
        eleid_veto      = tas::els_id2012_veto().at(idx);
        eleid_loose     = tas::els_id2012_loose().at(idx);
        eleid_medium    = tas::els_id2012_medium().at(idx);
        eleid_tight     = tas::els_id2012_tight().at(idx);
        is_eleid_veto   = electronId_WP2012_v3(idx, VETO);
        is_eleid_loose  = electronId_WP2012_v3(idx, LOOSE);
        is_eleid_medium = electronId_WP2012_v3(idx, MEDIUM);
        is_eleid_tight  = electronId_WP2012_v3(idx, TIGHT);
        dphiin          = tas::els_dPhiIn().at(idx);
        detain          = tas::els_dEtaIn().at(idx);
        hoe             = tas::els_hOverE().at(idx);
        sieie           = tas::els_sigmaIEtaIEta().at(idx);
        ooemoop         = fabs( (1.0/tas::els_ecalEnergy().at(idx)) - (tas::els_eOverPIn().at(idx)/tas::els_ecalEnergy().at(idx)) );
        conv_dist       = tas::els_conv_dist().at(idx);
        conv_dcot       = tas::els_conv_dcot().at(idx);

        trkiso  = tas::els_tkIso().at(idx);
        ecaliso = tas::els_ecalIso().at(idx);
        hcaliso = tas::els_hcalIso().at(idx);
        detiso  = electronIsolation_rel_v1(idx, true);

        trkiso04  = tas::els_tkIso04().at(idx);
        ecaliso04 = tas::els_ecalIso04().at(idx);
        hcaliso04 = tas::els_hcalIso04().at(idx);
        detiso04  = trkiso04 + hcaliso04 + ecaliso04;
        if (fabs(sc_p4.eta()) >= 1.479) detiso += tas::els_ecalIso04().at(idx);
        else detiso += std::max(tas::els_ecalIso04().at(idx) - 1., 0.0);
        detiso /= p4.pt();

        if (vtxidx >= 0)
        {
            electronIsoValuePF2012(chiso, emiso, nhiso, 0.3, idx, vtxidx);
            pfiso = (chiso + emiso + nhiso) / p4.pt();

            electronIsoValuePF2012(chiso04, emiso04, nhiso04, 0.4, idx, vtxidx);
            pfiso04 = (chiso04 + emiso04 + nhiso04) / p4.pt();
        }

    } // end electron block

    // muons information 
    // ---------------------- // 

    if (abs(id) == 13)
    {
        exists = true;
        is_mu = true;
        p4 = tas::mus_p4().at(idx);
        charge = tas::mus_charge().at(idx);
        pdgid = charge * -13;
        type = tas::mus_type().at(idx);

        const int trkidx = tas::mus_trkidx().at(idx);
        if (trkidx >= 0 && vtxidx >= 0)
        {
            std::pair<float, float> cord0 = trks_d0_pv(trkidx, vtxidx);
            std::pair<float, float> cordz = trks_dz_pv(trkidx, vtxidx);
            d0 = cord0.first;
            dz = cordz.first;
            d0err = cord0.second;
            dzerr = cordz.second;
        }
        ip3d    = tas::mus_ip3d().at(idx);;
        ip3derr = tas::mus_ip3derr().at(idx);;

        if (!tas::evt_isRealData())
        {
            mcp4      = tas::mus_mc_p4().at(idx);
            mc_momp4  = tas::mus_mc_motherp4().at(idx);
            mcid      = tas::mus_mc_id().at(idx);
            momid     = tas::mus_mc_motherid().at(idx);
            mc3id     = tas::mus_mc3_id().at(idx);
            mc3_momid = tas::mus_mc3_motherid().at(idx);

            const int mc3idx = tas::mus_mc3idx().at(idx);
            if (mc3idx >= 0)
            {
                mc3p4 = tas::genps_p4().at(mc3idx);
            }
        }

        gfit_p4    = tas::mus_gfit_p4().at(idx);
        is_global  = ((tas::mus_type().at(idx) & (1<<1)) != 0);
        is_tracker = ((tas::mus_type().at(idx) & (1<<2)) != 0);
        is_stamu   = ((tas::mus_type().at(idx) & (1<<3)) != 0);
        is_pfmu    = ((tas::mus_type().at(idx) & (1<<5)) != 0);

        nsihits      = tas::mus_validHits().at(idx);
        nstahits     = tas::mus_gfit_validSTAHits().at(idx);
        nstations    = tas::mus_numberOfMatchedStations().at(idx);
        chi2         = tas::mus_gfit_chi2().at(idx);
        ndof         = tas::mus_gfit_ndof().at(idx);
        pterr        = tas::mus_ptErr().at(idx);
        ecal_vetodep = tas::mus_iso_ecalvetoDep().at(idx);
        hcal_vetodep = tas::mus_iso_hcalvetoDep().at(idx);

        int ctfidx = tas::mus_trkidx().at(idx);
        if (ctfidx >= 0)
        {
            npixelhits = tas::trks_valid_pixelhits().at(ctfidx);
            nsilayers  = tas::trks_nlayers().at(ctfidx);            
        }

        if (vtxidx >= 0)
        {
            is_loosemu = passes_muid_wp2012(idx, mu2012_tightness::LOOSE);
            is_tightmu = passes_muid_wp2012(idx, mu2012_tightness::TIGHT);
        }

        trkiso  = tas::mus_iso03_sumPt().at(idx);
        ecaliso = tas::mus_iso03_emEt().at(idx);
        hcaliso = tas::mus_iso03_hadEt().at(idx);
        detiso  = muonIsoValue(idx, false);

        chiso = tas::mus_isoR03_pf_ChargedHadronPt().at(idx);
        emiso = tas::mus_isoR03_pf_PhotonEt().at(idx);
        nhiso = tas::mus_isoR03_pf_NeutralHadronEt().at(idx);
        pfiso = (chiso + emiso + nhiso) / p4.pt();
        dbeta = tas::mus_isoR03_pf_PUPt().at(idx);

        chiso04 = tas::mus_isoR04_pf_ChargedHadronPt().at(idx);
        emiso04 = tas::mus_isoR04_pf_PhotonEt().at(idx);
        nhiso04 = tas::mus_isoR04_pf_NeutralHadronEt().at(idx);
        pfiso04 = (chiso04 + emiso04 + nhiso04) / p4.pt();        
        dbeta04 = tas::mus_isoR04_pf_PUPt().at(idx);

    } // end muon block
}

void SingleLeptonInfo::SetBranches(TTree &tree)
{
    tree.Branch("run"             , &run             );
    tree.Branch("ls"              , &ls              );
    tree.Branch("evt"             , &evt             );
    tree.Branch("sample"          , &sample          );
    tree.Branch("dataset"         , &dataset         );
    tree.Branch("filename"        , &filename        );
    tree.Branch("is_real_data"    , &is_real_data    );
    tree.Branch("scale1fb"        , &scale1fb        );
    tree.Branch("scale1fb_cms2"   , &scale1fb_cms2   );
    tree.Branch("lumi"            , &lumi            );
    tree.Branch("xsec"            , &xsec            );
    tree.Branch("nevts_aod"       , &nevts_aod       );
    tree.Branch("nevts_cms2"      , &nevts_cms2      );
    tree.Branch("nevts_file"      , &nevts_file      );
    tree.Branch("filt_eff"        , &filt_eff        );
    tree.Branch("passes_id"       , &passes_id       ); 
    tree.Branch("passes_iso"      , &passes_iso      ); 
    tree.Branch("is_num"          , &is_num          ); 
    tree.Branch("is_den"          , &is_den          ); 
    tree.Branch("is_fo"           , &is_fo           ); 
    tree.Branch("is_mu"           , &is_mu           ); 
    tree.Branch("is_el"           , &is_el           ); 
    tree.Branch("exists"          , &exists          ); 
    tree.Branch("is_fromw"        , &is_fromw        ); 
    tree.Branch("charge"          , &charge          ); 
    tree.Branch("pdgid"           , &pdgid           ); 
    tree.Branch("type"            , &type            ); 
    tree.Branch("d0"              , &d0              ); 
    tree.Branch("d0err"           , &d0err           ); 
    tree.Branch("dz"              , &dz              ); 
    tree.Branch("ip3d"            , &ip3d            ); 
    tree.Branch("ip3derr"         , &ip3derr         ); 
    tree.Branch("dzerr"           , &dzerr           ); 
    tree.Branch("mt"              , &mt              ); 
    tree.Branch("corpfiso"        , &corpfiso        ); 
    tree.Branch("pfiso"           , &pfiso           ); 
    tree.Branch("chiso"           , &chiso           ); 
    tree.Branch("emiso"           , &emiso           ); 
    tree.Branch("nhiso"           , &nhiso           ); 
    tree.Branch("corpfiso04"      , &corpfiso04      ); 
    tree.Branch("pfiso04"         , &pfiso04         ); 
    tree.Branch("chiso04"         , &chiso04         ); 
    tree.Branch("emiso04"         , &emiso04         ); 
    tree.Branch("nhiso04"         , &nhiso04         ); 
    tree.Branch("cordetiso"       , &cordetiso       ); 
    tree.Branch("detiso"          , &detiso          ); 
    tree.Branch("trkiso"          , &trkiso          ); 
    tree.Branch("ecaliso"         , &ecaliso         ); 
    tree.Branch("hcaliso"         , &hcaliso         ); 
    tree.Branch("cordetiso04"     , &cordetiso04     ); 
    tree.Branch("detiso04"        , &detiso04        ); 
    tree.Branch("trkiso04"        , &trkiso04        ); 
    tree.Branch("ecaliso04"       , &ecaliso04       ); 
    tree.Branch("hcaliso04"       , &hcaliso04       ); 
    tree.Branch("effarea"         , &effarea         ); 
    tree.Branch("effarea04"       , &effarea04       ); 
    tree.Branch("dbeta"           , &dbeta           ); 
    tree.Branch("dbeta04"         , &dbeta04         ); 
    tree.Branch("sf_lepeff"       , &sf_lepeff       ); 
    tree.Branch("sf_trig"         , &sf_trig         ); 
    tree.Branch("mcid"            , &mcid            ); 
    tree.Branch("mc3id"           , &mc3id           ); 
    tree.Branch("momid"           , &momid           ); 
    tree.Branch("mc3_momid"       , &mc3_momid       ); 
    tree.Branch("q3agree"         , &q3agree         ); 
    tree.Branch("is_conv"         , &is_conv         ); 
    tree.Branch("qsc"             , &qsc             ); 
    tree.Branch("qctf"            , &qctf            ); 
    tree.Branch("qgsf"            , &qgsf            ); 
    tree.Branch("nmhits"          , &nmhits          ); 
    tree.Branch("eleid_veto"      , &eleid_veto      ); 
    tree.Branch("eleid_loose"     , &eleid_loose     ); 
    tree.Branch("eleid_medium"    , &eleid_medium    ); 
    tree.Branch("eleid_tight"     , &eleid_tight     ); 
    tree.Branch("is_eleid_veto"   , &is_eleid_veto   ); 
    tree.Branch("is_eleid_loose"  , &is_eleid_loose  ); 
    tree.Branch("is_eleid_medium" , &is_eleid_medium ); 
    tree.Branch("is_eleid_tight"  , &is_eleid_tight  ); 
    tree.Branch("dphiin"          , &dphiin          ); 
    tree.Branch("detain"          , &detain          ); 
    tree.Branch("sieie"           , &sieie           ); 
    tree.Branch("hoe"             , &hoe             ); 
    tree.Branch("ooemoop"         , &ooemoop         ); 
    tree.Branch("conv_dist"       , &conv_dist       ); 
    tree.Branch("conv_dcot"       , &conv_dcot       ); 
    tree.Branch("is_global"       , &is_global       ); 
    tree.Branch("is_tracker"      , &is_tracker      ); 
    tree.Branch("is_stamu"        , &is_stamu        ); 
    tree.Branch("is_pfmu"         , &is_pfmu         ); 
    tree.Branch("is_loosemu"      , &is_loosemu      ); 
    tree.Branch("is_tightmu"      , &is_tightmu      ); 
    tree.Branch("npixelhits"      , &npixelhits      ); 
    tree.Branch("nsihits"         , &nsihits         ); 
    tree.Branch("nsilayers"       , &nsilayers       ); 
    tree.Branch("nstahits"        , &nstahits        ); 
    tree.Branch("nstations"       , &nstations       ); 
    tree.Branch("chi2"            , &chi2            ); 
    tree.Branch("ndof"            , &ndof            ); 
    tree.Branch("pterr"           , &pterr           ); 
    tree.Branch("ecal_vetodep"    , &ecal_vetodep    ); 
    tree.Branch("hcal_vetodep"    , &hcal_vetodep    ); 

    tree.Branch("p4"      , "LorentzVector" , &p4      );
    tree.Branch("mcp4"    , "LorentzVector" , &mcp4    );
    tree.Branch("mc3p4"   , "LorentzVector" , &mc3p4   );
    tree.Branch("mc_momp4", "LorentzVector" , &mc_momp4);
    tree.Branch("gsf_p4"  , "LorentzVector" , &gsf_p4  );
    tree.Branch("ctf_p4"  , "LorentzVector" , &ctf_p4  );
    tree.Branch("sc_p4"   , "LorentzVector" , &sc_p4   );
    tree.Branch("gfit_p4" , "LorentzVector" , &gfit_p4 );
}

std::ostream& operator<< (std::ostream& out, const SingleLeptonInfo& info)
{
    out << "run             = " << info.run             << std::endl;
    out << "ls              = " << info.ls              << std::endl;
    out << "evt             = " << info.evt             << std::endl;
    out << "sample          = " << info.sample          << std::endl;
    out << "dataset         = " << info.dataset         << std::endl;
    out << "filename        = " << info.filename        << std::endl;
    out << "is_real_data    = " << info.is_real_data    << std::endl;
    out << "scale1fb        = " << info.scale1fb        << std::endl;
    out << "scale1fb_cms2   = " << info.scale1fb_cms2   << std::endl;
    out << "lumi            = " << info.lumi            << std::endl;
    out << "xsec            = " << info.xsec            << std::endl;
    out << "nevts_aod       = " << info.nevts_aod       << std::endl;
    out << "nevts_cms2      = " << info.nevts_cms2      << std::endl;
    out << "nevts_file      = " << info.nevts_file      << std::endl;
    out << "kfactor         = " << info.kfactor         << std::endl;
    out << "filt_eff        = " << info.filt_eff        << std::endl;
    out << "p4.mass()       = " << info.p4.mass()       << std::endl;
    out << "passes_id       = " << info.passes_id       << std::endl;
    out << "passes_iso      = " << info.passes_iso      << std::endl;
    out << "is_num          = " << info.is_num          << std::endl;
    out << "is_den          = " << info.is_den          << std::endl;
    out << "is_fo           = " << info.is_fo           << std::endl;
    out << "is_mu           = " << info.is_mu           << std::endl;
    out << "is_el           = " << info.is_el           << std::endl;
    out << "exists          = " << info.exists          << std::endl;
    out << "is_fromw        = " << info.is_fromw        << std::endl;
    out << "charge          = " << info.charge          << std::endl;
    out << "pdgid           = " << info.pdgid           << std::endl;
    out << "type            = " << info.type            << std::endl;
    out << "d0              = " << info.d0              << std::endl;
    out << "d0err           = " << info.d0err           << std::endl;
    out << "dz              = " << info.dz              << std::endl;
    out << "dzerr           = " << info.dzerr           << std::endl;
    out << "mt              = " << info.mt              << std::endl;
    out << "corpfiso        = " << info.corpfiso        << std::endl;
    out << "pfiso           = " << info.pfiso           << std::endl;
    out << "chiso           = " << info.chiso           << std::endl;
    out << "emiso           = " << info.emiso           << std::endl;
    out << "nhiso           = " << info.nhiso           << std::endl;
    out << "corpfiso04      = " << info.corpfiso04      << std::endl;
    out << "pfiso04         = " << info.pfiso04         << std::endl;
    out << "chiso04         = " << info.chiso04         << std::endl;
    out << "emiso04         = " << info.emiso04         << std::endl;
    out << "nhiso04         = " << info.nhiso04         << std::endl;
    out << "cordetiso       = " << info.cordetiso       << std::endl;
    out << "detiso          = " << info.detiso          << std::endl;
    out << "trkiso          = " << info.trkiso          << std::endl;
    out << "ecaliso         = " << info.ecaliso         << std::endl;
    out << "hcaliso         = " << info.hcaliso         << std::endl;
    out << "cordetiso04     = " << info.cordetiso04     << std::endl;
    out << "detiso04        = " << info.detiso04        << std::endl;
    out << "trkiso04        = " << info.trkiso04        << std::endl;
    out << "ecaliso04       = " << info.ecaliso04       << std::endl;
    out << "hcaliso04       = " << info.hcaliso04       << std::endl;
    out << "effarea         = " << info.effarea         << std::endl;
    out << "effarea04       = " << info.effarea04       << std::endl;
    out << "dbeta           = " << info.dbeta           << std::endl;
    out << "dbeta04         = " << info.dbeta04         << std::endl;
    out << "sf_lepeff       = " << info.sf_lepeff       << std::endl;
    out << "sf_trig         = " << info.sf_trig         << std::endl;
    out << "mcp4.mass()     = " << info.mcp4.mass()     << std::endl;
    out << "mc3p4.mass()    = " << info.mc3p4.mass()    << std::endl;
    out << "mc_momp4.mass() = " << info.mc_momp4.mass() << std::endl;
    out << "mcid            = " << info.mcid            << std::endl;
    out << "mc3id           = " << info.mc3id           << std::endl;
    out << "momid           = " << info.momid           << std::endl;
    out << "mc3_momid       = " << info.mc3_momid       << std::endl;
    out << "gsf_p4.mass()   = " << info.gsf_p4.mass()   << std::endl;
    out << "ctf_p4.mass()   = " << info.ctf_p4.mass()   << std::endl;
    out << "sc_p4.mass()    = " << info.sc_p4.mass()    << std::endl;
    out << "q3agree         = " << info.q3agree         << std::endl;
    out << "is_conv         = " << info.is_conv         << std::endl;
    out << "qsc             = " << info.qsc             << std::endl;
    out << "qctf            = " << info.qctf            << std::endl;
    out << "qgsf            = " << info.qgsf            << std::endl;
    out << "nmhits          = " << info.nmhits          << std::endl;
    out << "eleid_veto      = " << info.eleid_veto      << std::endl;
    out << "eleid_loose     = " << info.eleid_loose     << std::endl;
    out << "eleid_medium    = " << info.eleid_medium    << std::endl;
    out << "eleid_tight     = " << info.eleid_tight     << std::endl;
    out << "is_eleid_veto   = " << info.is_eleid_veto   << std::endl;
    out << "is_eleid_loose  = " << info.is_eleid_loose  << std::endl;
    out << "is_eleid_medium = " << info.is_eleid_medium << std::endl;
    out << "is_eleid_tight  = " << info.is_eleid_tight  << std::endl;
    out << "dphiin          = " << info.dphiin          << std::endl;
    out << "detain          = " << info.detain          << std::endl;
    out << "sieie           = " << info.sieie           << std::endl;
    out << "hoe             = " << info.hoe             << std::endl;
    out << "ooemoop         = " << info.ooemoop         << std::endl;
    out << "conv_dist       = " << info.conv_dist       << std::endl;
    out << "conv_dcot       = " << info.conv_dcot       << std::endl;
    out << "gfit_p4.mass()  = " << info.gfit_p4.mass()  << std::endl;
    out << "is_global       = " << info.is_global       << std::endl;
    out << "is_tracker      = " << info.is_tracker      << std::endl;
    out << "is_stamu        = " << info.is_stamu        << std::endl;
    out << "is_pfmu         = " << info.is_pfmu         << std::endl;
    out << "is_loosemu      = " << info.is_loosemu      << std::endl;
    out << "is_tightmu      = " << info.is_tightmu      << std::endl;
    out << "npixelhits      = " << info.npixelhits      << std::endl;
    out << "nsihits         = " << info.nsihits         << std::endl;
    out << "nsilayers       = " << info.nsilayers       << std::endl;
    out << "nstahits        = " << info.nstahits        << std::endl;
    out << "nstations       = " << info.nstations       << std::endl;
    out << "chi2            = " << info.chi2            << std::endl;
    out << "ndof            = " << info.ndof            << std::endl;
    out << "pterr           = " << info.pterr           << std::endl;
    out << "ecal_vetodep    = " << info.ecal_vetodep    << std::endl;
    out << "hcal_vetodep    = " << info.hcal_vetodep    << std::endl;
    return out;
}

// -------------------------------------------------//
// Simple class to hold your fill the SingleLeptonInfo 
// -------------------------------------------------//

class SingleLeptonNtupleMaker
{
    public:
        // construct:
        SingleLeptonNtupleMaker() = delete;
        SingleLeptonNtupleMaker
        (
            const dy::Sample::Info sample_info, 
            const std::string& output_filename, 
            const double lumi,
            const long num_events,
            const bool verbose
        );

        // destroy:
        ~SingleLeptonNtupleMaker();

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
        bool m_verbose;
        SingleLeptonInfo m_info;
        TFile& m_file;
        TTree& m_tree;
};

// construct:
SingleLeptonNtupleMaker::SingleLeptonNtupleMaker
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
    , m_info()
    , m_file(*TFile::Open(output_filename.c_str(), "RECREATE"))
    , m_tree(*new TTree("tree", "DY Exercise TTree"))
{
}

// destroy:
SingleLeptonNtupleMaker::~SingleLeptonNtupleMaker()
{
}

// ------------------------------------ //
// helper functions and classes 
// ------------------------------------ //

// ------------------------------------ //
// Stuff to do before job starts
// ------------------------------------ //

void SingleLeptonNtupleMaker::BeginJob()
{
    m_info.SetBranches(m_tree);
}

// ------------------------------------ //
// Stuff to do on each event 
// ------------------------------------ //

void SingleLeptonNtupleMaker::Analyze(const long event, const std::string& current_filename)
{
    if (m_verbose)
    {
        std::cout << "\n[SingleLeptonNtupleMaker] Running on run, ls, event: " 
            << tas::evt_run()       << ", "
            << tas::evt_lumiBlock() << ", "
            << tas::evt_event()     << std::endl;
    }

    // electrons
    // ------------------------------------ //

    for (size_t el_idx = 0; el_idx < tas::els_p4().size(); ++el_idx)
    {
        // reset the TTree variables
        m_info.Reset();

        // event information 
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

        // electron information
        m_info.FillCommon(11, el_idx);

        // fill the tree
        if (m_verbose) {std::cout << m_info << std::endl;}
        m_tree.Fill();

    } // end electrons

    // muons
    // ------------------------------------ //

    for (size_t mu_idx = 0; mu_idx < tas::mus_p4().size(); ++mu_idx)
    {
        // reset the TTree variables
        m_info.Reset();

        // event information 
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

        // electron information
        m_info.FillCommon(13, mu_idx);

        // fill the tree
        if (m_verbose) {std::cout << m_info << std::endl;}
        m_tree.Fill();

    } // end muons

    // done
    return;
}

// ------------------------------------ //
// Stuff to do after job finishes
// ------------------------------------ //

void SingleLeptonNtupleMaker::EndJob()
{
    std::cout << "[SingleLeptonNtupleMaker] Saving TTree to output file: " << m_output_filename << std::endl;
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
        output_file = Form("babies/singlelep/%s_baby.root", sample_info.name.c_str());
    }
    lt::mkdir(lt::dirname(output_file), /*force=*/true);

    // run the looper
    // -------------------------------------------------------------------------------------------------//

    std::cout << "[dy_create_ntuple] TChain set to run on:\n";
    rt::PrintFilesFromTChain(chain);

    // create looper
    const long long num_events_to_run = (number_of_events < 0 ? chain->GetEntries() : std::min(number_of_events, chain->GetEntries()));
    SingleLeptonNtupleMaker looper
    (
        sample_info,
        output_file, 
        lumi,
        num_events_to_run,
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
