#include "Analysis/DrellYan/interface/Sample.h"
#include "AnalysisTools/LanguageTools/interface/StringTools.h"
#include <stdexcept>
#include <vector>
#include <string>
#include "TChain.h"
#include "TColor.h"

namespace dy
{
    // Array of SampleInfo's with the relevant metadata
    static const Sample::Info s_SampleInfos[Sample::static_size] = 
    {
        {
            // name
            "data", 
            // title
            "Run2012A-recover-06Aug2012-v1",
            // latex
            "Run2012A-recover-06Aug2012-v1",
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root"
            ",/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root",
            // color
            kBlack,
            // filter efficiency
            1.0,
            // Sample
            Sample::data
        },
        {
            // name
            "dyll", 
            // title
            "Z/#gamma #rightarrow l^{+}l^{-}",
            // latex
            "$DY \\rightarrow \\ell \\ell$",
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root",
            // color
            kOrange-2,
            // filter efficiency
            (27137253.0/30459500.0),
            // Sample
            Sample::dyll
        },
        {
            // name
            "wjets", 
            // title
            "W+jets #rightarrow l#nu", 
            // latex
            "$W+jets \\rightarrow \\ell \\nu$", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-28/merged_ntuple_1[0-9].root",
            // color
            kGray+1,
            // filter efficiency
            (14890630.0/18393090.0),
            // Sample
            Sample::wjets
        },
        {
            // name
            "ttdil", 
            // title
            "t#bar{t} #rightarrow llX", 
            // latex
            "$t\\overline{t} \\rightarrow \\ell \\ell X$", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-24/merged_ntuple_1[0-9].root",
            // color
            kOrange+7,
            // filter efficiency
            (11902045.0/12119013.0),
            // Sample
            Sample::ttdil
        },
        {
             // name
            "ttslq", 
            // title
            "t#bar{t} #rightarrow l(q #rightarrow l)X", 
            // latex
            "$t\\overline{t} \\rightarrow \\ell (q \\rightarrow \\ell) X$", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/V05-03-24/merged_ntuple_1[0-9].root",
            // color
            kOrange+6,
            // filter efficiency
            (23453443.0/25384818.0),
            // Sample
            Sample::ttslq
        },
        {
            // name
            "tthad", 
            // title
            "t#bar{t} #rightarrow hadrons", 
            // latex
            "$t\\overline{t} \\rightarrow hadrons$", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/V05-03-24/merged_ntuple_1[0-9].root",
            // color
            kOrange+3,
            // filter efficiency
            (21211839.0/31223821.0),
            // Sample
            Sample::tthad
        },
        {
            // name
            "qcdmu15", 
            // title
            "QCD (#mu15 enriched)", 
            // latex
            "QCD (\\mu15 enriched)", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-18_slim/merged_ntuple_1[0-9].root",
            // color
            kMagenta,
            // filter efficiency
            (19770457.0/20614602.0),
            // Sample
            Sample::qcdmu15
        },
        {
            // name
            "ww2l2nu", 
            // title
            "WW #rightarrow 2l + 2#nu", 
            // latex
            "WW \\rightarrow 2\\ell + 2\\nu", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root",
            // color
            kGreen+2,
            // filter efficiency
            (1933235.0/1933235.0),
            // Sample
            Sample::ww2l2nu
        },
        {
            // name
            "wz2l2q", 
            // title
            "WZ #rightarrow 2l + 2q", 
            // latex
            "WZ \\rightarrow 2\\ell + 2q", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root",
            // color
            kRed+2,
            // filter efficiency
            (2937874.0/3215990.0),
            // Sample
            Sample::wz2l2q
        },
        {
            // name
            "wz3lnu", 
            // title
            "WZ #rightarrow 3l + #nu", 
            // latex
            "WZ \\rightarrow 3\\ell + \\nu", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root",
            // color
            kRed-7,
            // filter efficiency
            (2937874.0/3215990.0),
            // Sample
            Sample::wz3lnu
        },
        {
            // name
            "zz2l2nu", 
            // title
            "ZZ #rightarrow 2l + 2#nu", 
            // latex
            "ZZ \\rightarrow 2\\ell + 2\\nu", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-23/merged_ntuple_1[0-9].root",
            // color
            kBlue,
            // filter efficiency
            (857982.0/954911.0),
            // Sample
            Sample::zz2l2nu
        },
        {
            // name
            "zz2l2q", 
            // title
            "ZZ #rightarrow 2l + 2q", 
            // latex
            "ZZ \\rightarrow 2\\ell + 2q", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root",
            // color
            kBlue+2,
            // filter efficiency
            (1777571.0/1936727.0),
            // Sample
            Sample::zz2l2q
        },
        {
            // name
            "zz4l", 
            // title
            "ZZ #rightarrow 4l", 
            // latex
            "ZZ \\rightarrow 4\\ell", 
            // ntuple_path
            "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root",
            // color
            kBlue-5,
            // filter efficiency
            (4807893.0/4807893.0),
            // Sample
            Sample::zz4l
        },
    };

    // Get the Sample from a string
    Sample::value_type GetSampleFromName(const std::string& sample_name)
    {
        Sample::value_type sample = Sample::static_size;
        for (int i = 0; i < Sample::static_size; ++i)
        {
            if (std::string(s_SampleInfos[i].name) == sample_name)
            {
                sample = static_cast<Sample::value_type>(i);
                break;
            }
        }

        // throw if not found
        if (sample == Sample::static_size)
        {
            throw std::domain_error(Form("[dy::GetSampleInfo] Error: sample %s not found!", sample_name.c_str()));
        }

        return sample; 
    }

    // Get the Sample from a number
    Sample::value_type GetSampleFromNumber(const int sample_num)
    {
        // throw if not found
        if (sample_num < 0 or sample_num >= Sample::static_size)
        {
            throw std::domain_error("[dy::GetSampleInfo] Error: sample number out of bounds!");
        }
        return static_cast<Sample::value_type>(sample_num);
    }


    // test if a string is a sample name 
    bool IsSample(const std::string& sample_name)
    {
        for (int i = 0; i < Sample::static_size; ++i)
        {
            if (std::string(s_SampleInfos[i].name) == sample_name)
            {
                return true; 
            }
        }
        return false;
    }

    // wrapper function to get the SampleInfo
    Sample::Info GetSampleInfo(const Sample::value_type& sample)
    {
        return s_SampleInfos[sample]; 
    }

    Sample::Info GetSampleInfo(const std::string& sample_name)
    {
        Sample::value_type sample = GetSampleFromName(sample_name); 
        return GetSampleInfo(sample); 
    }

    Sample::Info GetSampleInfo(const int sample_num)
    {
        Sample::value_type sample = GetSampleFromNumber(sample_num); 
        return GetSampleInfo(sample); 
    }

    // get the chain from the Sample
    TChain* GetSampleTChain(const Sample::value_type& sample)
    {
        const Sample::Info info = GetSampleInfo(sample);

        // build the list of files
        std::vector<std::string> vpath = lt::string_split(info.ntuple_path, ",");

        // build the chain
        TChain* chain = new TChain("Events");
        for (const auto& file : vpath)
        {
            chain->Add(file.c_str());
        }
        return chain;
    }

} // namespace dy
 
