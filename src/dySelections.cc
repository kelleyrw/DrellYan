#include "Analysis/DrellYan/interface/dySelections.h"

// CMS2 Includes
#include "CMS2/NtupleMacrosHeader/interface/CMS2.h"
#include "CMS2/NtupleMacrosCore/interface/electronSelections.h"
#include "CMS2/NtupleMacrosCore/interface/electronSelectionsParameters.h"
#include "CMS2/NtupleMacrosCore/interface/muonSelections.h"
#include "CMS2/NtupleMacrosCore/interface/trackSelections.h"
#include "CMS2/NtupleMacrosCore/interface/eventSelections.h"
#include "CMS2/NtupleMacrosCore/interface/triggerUtils.h"
#include "CMS2/NtupleMacrosCore/interface/susySelections.h"
#include "CMS2/NtupleMacrosCore/interface/MITConversionUtilities.h"

/////////////////////////////////////////////////////////////////
///                                                           ///
///                                                           ///
///                                                           ///
///          Drell-Yan Selections                             ///
///                                                           ///
///                                                           ///
///                                                           ///
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 lepton impact parameters 
// uses CTF track for muons and GSF tracks for elections
// if no matching track found, return bogus value of -999999
// calc w.r.t first good vertex
////////////////////////////////////////////////////////////////////////////////////////////     
float dy::leptonD0(const int id, const int idx)
{
    const int vtxidx = firstGoodVertex();
    if (vtxidx < 0)
    {
        std::cout << "[dy::leptonD0] WARNING - first good vertex index < 0.  Returning bogus value 999999" << std::endl;
        return 999999.0;
    }
    if (abs(id)==13)
    {
        const int trkidx = tas::mus_trkidx().at(idx);
        if (trkidx >= 0)
        {
            return trks_d0_pv(trkidx, vtxidx).first;
        }
    }
    else if (abs(id)==11)
    {
        const int gsfidx = tas::els_gsftrkidx().at(idx);
        if (gsfidx >= 0) 
        {
            return gsftrks_d0_pv(gsfidx, vtxidx).first;
        }
    }

    // return bogus for non electon/muon
    return -999999.0;
}

float dy::leptonDz(const int id, const int idx)
{
    const int vtxidx = firstGoodVertex();
    if (vtxidx < 0)
    {
        std::cout << "[dy::leptonDz] WARNING - first good vertex index < 0.  Returning bogus value 999999" << std::endl;
        return 999999.0;
    }
    if (abs(id)==13)
    {
        const int trkidx = tas::mus_trkidx().at(idx);
        if (trkidx >= 0)
        {
            return trks_dz_pv(trkidx, vtxidx).first;
        }
    }
    else if (abs(id)==11)
    {
        const int gsfidx = tas::els_gsftrkidx().at(idx);
        if (gsfidx >= 0)
        {
            return gsftrks_dz_pv(gsfidx, vtxidx).first;
        }
    }

    // return bogus for non electon/muon
    return -999999.0;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// good lepton (passes ID)
////////////////////////////////////////////////////////////////////////////////////////////     
bool dy::isGoodLepton(const int id, const int idx)
{
    using namespace tas;

    // selections from AN-12-067
    // http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2012_067_v6.pdf

    // electrons
    if (abs(id) == 11)
    {
        const float d0                = leptonD0(id, idx);
        const float dz                = leptonDz(id, idx);
        const float aeta_sc           = fabs(els_p4().at(idx).eta());
        const bool vtx_fit_conversion = isMITConversion(idx, /*nWrongHitsMax*/0, /*probMin*/1e-6, /*dlMin*/2.0, /*matchCTF*/true, /*requireArbitratedMerged*/false);

        if (els_p4().at(idx).Et() < 25/*GeV*/)    {return false;} // E_T > 2.5 GeV
        if (1.4442 < aeta_sc and aeta_sc < 1.566) {return false;} // reject crack electons defined by 1.4442 < |eta_{sc}| < 1.566
        if (aeta_sc > 2.5)                        {return false;} // |eta_{sc}| < 2.5
        if (fabs(d0) > 0.02/*cm*/)                {return false;} // Its tracker track has transverse impact parameter dxy < 200 µm w.r.t. the primary vertex
        if (fabs(dz) > 0.1/*cm*/)                 {return false;} // The longitudinal distance of the tracker track wrt. the primary vertex is dz < 5 mm
        if (els_exp_innerlayers().at(idx) > 1)    {return false;} // # missing hits on the track < 1
        if (vtx_fit_conversion)                   {return false;} // reject convertions

        // |1/E - 1/p|
        const float ooemoop = fabs((1.0/els_ecalEnergy().at(idx)) - (els_eOverPIn().at(idx)/els_ecalEnergy().at(idx)));

        // shape variables
        if (aeta_sc < 1.4442) // barrel
        {
            if (els_sigmaIEtaIEta().at(idx) > 0.01 ) {return false;}
            if (fabs(els_dPhiIn().at(idx))  > 0.06 ) {return false;}
            if (fabs(els_dEtaIn().at(idx))  > 0.004) {return false;}
            if (ooemoop                     > 0.05 ) {return false;}
            if (els_hOverE().at(idx)        > 0.12 ) {return false;}
        }
        else // endcap
        {
            if (els_sigmaIEtaIEta().at(idx) > 0.03 ) {return false;}
            if (fabs(els_dPhiIn().at(idx))  > 0.03 ) {return false;}
            if (fabs(els_dEtaIn().at(idx))  > 0.007) {return false;}
            if (ooemoop                     > 0.05 ) {return false;}
            if (els_hOverE().at(idx)        > 0.10 ) {return false;}
        }
        
        // if we're here, then it passes
        return true;
    }

    // muons
    if (abs(id) == 13)
    {
        const bool is_global  = ((mus_type().at(idx) & (1<<1)) != 0);
        const bool is_pfmu    = ((mus_type().at(idx) & (1<<5)) != 0);
        const int ctfidx      = mus_trkidx().at(idx);
        const float d0        = leptonD0(id, idx);
        const float dz        = leptonDz(id, idx);
        const float chi2ndof  = mus_gfit_chi2().at(idx)/mus_gfit_ndof().at(idx);

        if (mus_p4().at(idx).pt() < 25/*GeV*/)            {return false;} // p_T > 2.5 GeV
        if (fabs(mus_p4().at(idx).eta()) > 2.1)           {return false;} // |eta| < 2.1
        if (not is_global)                                {return false;} // The candidate is reconstructed as a Global Muon
        if (not is_pfmu)                                  {return false;} // Particle-Flow muon id 
        if (chi2ndof >= 10)                               {return false;} // χ2/ndof of the global-muon track fit < 10
        if (mus_numberOfMatchedStations().at(idx) <= 1)   {return false;} // Muon segments in at least two muon stations
        if (mus_gfit_validSTAHits().at(idx) <= 0)         {return false;} // At least one muon chamber hit included in the global-muon track fit 
        if (fabs(d0) > 0.02/*cm*/)                        {return false;} // Its tracker track has transverse impact parameter dxy < 200 µm w.r.t. the primary vertex
        if (fabs(dz) > 0.5/*cm*/)                         {return false;} // The longitudinal distance of the tracker track wrt. the primary vertex is dz < 5 mm
        if (trks_valid_pixelhits().at(ctfidx) <= 0)       {return false;} // Number of pixel hits > 0
        if (trks_nlayers().at(ctfidx) <= 5)               {return false;} // Cut on number of tracker layers with hits > 5

        // if we're here, then it passes
        return true;
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 isolated lepton
////////////////////////////////////////////////////////////////////////////////////////////     
bool dy::isIsolatedLepton(const int id, const int idx)
{
    const float iso = dy::leptonIsolation(id, idx);

    // selections from AN-12-067
    // http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2012_067_v6.pdf

    // electrons
    if (abs(id) == 11)
    {
        return (iso < 0.15);
    }

    // muons
    if (abs(id) == 13)
    {
        return (iso < 0.12);
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 lepton isolation value
////////////////////////////////////////////////////////////////////////////////////////////     
double dy::leptonIsolation(const int id, const int idx)
{
    // electrons
    if (abs(id) == 11)
    {
        return dy::electronIsolationPF2012(idx);
    }

    // muons
    if (abs(id) == 13)
    {
        return dy::muonIsoValuePF2012(idx);
    }

    return -999999.0;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 effective area 
////////////////////////////////////////////////////////////////////////////////////////////     

float dy::EffectiveArea03(const int id, const int idx)
{
    if (abs(id)!=11)
        return -999990.0;

    // use SC eta
    float eta = fabs(tas::els_etaSC().at(idx));

    // get effective area from electronSelections.h
    //return fastJetEffArea03_v1(eta);  // used for HCP and ICHEP
    return fastJetEffArea03_v2(eta);    // 2013
}

float dy::EffectiveArea04(int id, int idx)
{
    if (abs(id)!=11)
        return -999990.0;

    // use SC eta
    float eta = fabs(tas::els_etaSC().at(idx));

    // get effective area from electronSelections.h
    //return fastJetEffArea04_v1(eta);  // used for HCP and ICHEP
    return fastJetEffArea04_v2(eta);    // 2013
}
///////////////////////////////////////////////////////////////////////////////////////////
// calculate PF-based isolation for electrons with rho*Aeff correction
// using cone size 03
///////////////////////////////////////////////////////////////////////////////////////////
float dy::electronIsolationPF2012(const int idx)
{
    return dy::electronIsolationPF2012_cone03(idx);
}

float dy::electronIsolationPF2012_cone03(const int idx)
{
    // electron pT
    const float pt = tas::els_p4().at(idx).pt();

    // get effective area
    const float AEff = EffectiveArea03(11, idx);

    // pf iso
    const float pfiso_ch = tas::els_iso03_pf2012ext_ch().at(idx);
    const float pfiso_em = tas::els_iso03_pf2012ext_em().at(idx);
    const float pfiso_nh = tas::els_iso03_pf2012ext_nh().at(idx);

    // rho
    const float rhoPrime = std::max(tas::evt_kt6pf_foregiso_rho(), 0.0f);
    const float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, 0.0f);
    const float pfiso = (pfiso_ch + pfiso_n) / pt;

    return pfiso;
}

float dy::electronIsolationPF2012_cone04(const int idx)
{
    // electron pT
    const float pt = tas::els_p4().at(idx).pt();

    // get effective area
    const float AEff = EffectiveArea04(11, idx);

    // pf iso
    const float pfiso_ch = tas::els_iso04_pf2012ext_ch().at(idx);
    const float pfiso_em = tas::els_iso04_pf2012ext_em().at(idx);
    const float pfiso_nh = tas::els_iso04_pf2012ext_nh().at(idx);

    // rho
    const float rhoPrime = std::max(tas::evt_kt6pf_foregiso_rho(), 0.0f);
    const float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, 0.0f);
    const float pfiso = (pfiso_ch + pfiso_n) / pt;

    return pfiso;
}

///////////////////////////////////////////////////////////////////////////////////////////
// calculate PF-based isolation for muon with Delta-Beta correction using cone size 04
///////////////////////////////////////////////////////////////////////////////////////////
float dy::muonIsoValuePF2012(const unsigned int imu)
{
    const float chiso     = tas::mus_isoR04_pf_ChargedHadronPt().at(imu);
    const float nhiso     = tas::mus_isoR04_pf_NeutralHadronEt().at(imu);
    const float emiso     = tas::mus_isoR04_pf_PhotonEt().at(imu);
    const float deltaBeta = tas::mus_isoR04_pf_PUPt().at(imu);
    const float pt        = tas::mus_p4().at(imu).pt();
    const float absiso    = chiso + max(0.0, nhiso + emiso - 0.5 * deltaBeta);
    return (absiso / pt);
}

///////////////////////////////////////////////////////////////////////////////////////////
// passes dilepton trigger
///////////////////////////////////////////////////////////////////////////////////////////

bool dy::passesTrigger(const int hyp_type)
{
    //----------------------------------------
    // no trigger requirements applied to MC
    //----------------------------------------

    if (!tas::evt_isRealData())
    {
        return true; 
    }

    //---------------------------------
    // triggers for dilepton datasets
    //---------------------------------

    // mm
    switch(hyp_type)
    {
        /*mu mu*/case 0: return passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v"); break;
        /*e mu*/ case 1: return false; break;
        /*e mu*/ case 2: return false; break;
        /*e e*/  case 3: return passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"); break;
        default: return false;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 selected lepton (passes ID and isolation)
////////////////////////////////////////////////////////////////////////////////////////////     
bool dy::isSelectedLepton(const int id, const int idx)
{
    return (dy::isGoodLepton(id, idx) && dy::isIsolatedLepton(id, idx));
}

////////////////////////////////////////////////////////////////////////////////////////////     
// 2012 selected hypothesis (passes ID and isolation)
////////////////////////////////////////////////////////////////////////////////////////////     
bool dy::isSelectedHypothesis(const int idx)
{
    if (!dy::isSelectedLepton(tas::hyp_lt_id().at(idx), tas::hyp_lt_index().at(idx)))
    {
        return false;
    }
    if (!dy::isSelectedLepton(tas::hyp_ll_id().at(idx), tas::hyp_ll_index().at(idx)))
    {
        return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// disambiguate between two passing hypotheses 
////////////////////////////////////////////////////////////////////////////////////////////     
int dy::ChooseBetterHypothesis(const int hyp1_idx, const int hyp2_idx)
{
    // both are dummy (< 0)
    if (hyp1_idx < 0 && hyp2_idx < 0)
    {
        return -99999;
    }

    // hyp1 dummy 
    if (hyp1_idx < 0)
    {
        return hyp2_idx;
    }

    // hyp2 dummy 
    if (hyp2_idx < 0)
    {
        return hyp1_idx;
    }

    // choose mumu over ee, reject emu
    // NOTE: MUMU 0, EMU 1, 2, EE 3
    const int hyp1_flavor_type = tas::hyp_type().at(hyp1_idx);
    const int hyp2_flavor_type = tas::hyp_type().at(hyp2_idx);
    
    const bool hyp1_is_em = (hyp1_flavor_type==1 or hyp1_flavor_type==2);
    const bool hyp2_is_em = (hyp2_flavor_type==1 or hyp2_flavor_type==2);
    const bool hyp1_is_ee = (hyp1_flavor_type==3);
    const bool hyp2_is_ee = (hyp2_flavor_type==3);
    const bool hyp1_is_mm = (hyp1_flavor_type==0);
    const bool hyp2_is_mm = (hyp2_flavor_type==0);

    if (hyp1_is_em and not hyp2_is_em)
    {
        return hyp2_idx;
    }
    if (not hyp1_is_em and hyp2_is_em)
    {
        return hyp1_idx;
    }
    if (hyp1_is_mm and hyp2_is_ee)
    {
        return hyp1_idx;
    }
    if (hyp1_is_ee and hyp2_is_mm)
    {
        return hyp2_idx;
    }

    // choose one closer to pdg value of mass of Z-boson
    static const float Mz = 91.1876;
    const float dm1 = fabs(tas::hyp_p4().at(hyp1_idx).mass() - Mz);
    const float dm2 = fabs(tas::hyp_p4().at(hyp2_idx).mass() - Mz);
    if (dm1 < dm2)
    {
        return hyp1_idx;
    }
    else
    {
        return hyp2_idx;
    }

    // if we're here, give bogus value
    return -99999;
}

////////////////////////////////////////////////////////////////////////////////////////////     
// Gen hypotheses 
// 1 mumu; 2 emu; 3 ee; -1 other;
////////////////////////////////////////////////////////////////////////////////////////////     
int dy::GenDileptonType()
{
    unsigned int nmus  = 0;
    unsigned int nels  = 0;
    for (size_t genps_idx = 0; genps_idx < tas::genps_id().size(); ++genps_idx)
    {
        if (tas::genps_status().at(genps_idx)  != 3 ) {continue;}
        if (abs(tas::genps_id().at(genps_idx)) == 11) {nels++;}
        if (abs(tas::genps_id().at(genps_idx)) == 13) {nmus++;}
        if (abs(tas::genps_id().at(genps_idx)) == 15)
        {
            for(size_t d_idx = 0; d_idx < tas::genps_lepdaughter_id().at(genps_idx).size(); ++d_idx)
            {
                const int daughter = abs(tas::genps_lepdaughter_id().at(genps_idx).at(d_idx));
                if (daughter == 11) {nels++;}
                if (daughter == 13) {nmus++;}
            }
        }
    }
    if ((nels + nmus) != 2    ) {return -1;}
    if (nmus == 2             ) {return 0;}
    if (nels == 1 && nmus == 1) {return 2;} 
    if (nels == 2             ) {return 3;}
    return -1;
}
