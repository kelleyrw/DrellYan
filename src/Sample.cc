#include "Analysis/DrellYan/interface/Sample.h"
#include "AnalysisTools/LanguageTools/interface/StringTools.h"
#include "AnalysisTools/LanguageTools/interface/OSTools.h"
#include <stdexcept>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include "TChain.h"
#include "TColor.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

namespace dy
{
    // PSET that holds the sample metadata
    static std::string sample_pset_path = "psets/dy_samples_cfg.py";

    // create a sample Info from the pset
    Sample::Info CreateSampleInfo(const Sample::value_type sample)
    {
        assert(!sample_pset_path.empty());
        assert(lt::file_exists(sample_pset_path));
        const edm::ParameterSet& process = edm::readPSetsFrom(sample_pset_path)->getParameter<edm::ParameterSet>("process");
        const auto& pset_vec             = process.getParameter<std::vector<edm::ParameterSet> >("dy_samples");
        assert(pset_vec.size() == Sample::static_size);
        const auto& pset = pset_vec[sample];
        Sample::Info info
        {
            pset.getParameter<std::string>("name"),
            pset.getParameter<std::string>("title"),
            pset.getParameter<std::string>("title"),
            pset.getParameter<std::string>("ntuple_path"),
            static_cast<Color_t>(pset.getParameter<int>("color")),
            pset.getParameter<double>("eff"),
            sample
        };
        return info;
    }

    // Array of SampleInfo's with the relevant metadata
    static Sample::Info s_SampleInfos[Sample::static_size]
    {
        CreateSampleInfo(Sample::data    ), 
        CreateSampleInfo(Sample::dyll    ), 
        CreateSampleInfo(Sample::wjets   ), 
        CreateSampleInfo(Sample::ttdil   ), 
        CreateSampleInfo(Sample::ttslq   ), 
        CreateSampleInfo(Sample::tthad   ), 
        CreateSampleInfo(Sample::qcdmu15 ), 
        CreateSampleInfo(Sample::ww2l2nu ), 
        CreateSampleInfo(Sample::wz2l2q  ), 
        CreateSampleInfo(Sample::wz3lnu  ), 
        CreateSampleInfo(Sample::zz2l2nu ), 
        CreateSampleInfo(Sample::zz2l2q  ), 
        CreateSampleInfo(Sample::zz4l    ) 
    };

    // rest the pset 
    /*static*/ void Sample::SetPsetPath(const std::string& pset_path)
    {
        sample_pset_path = pset_path;
        for (int i = 0; i < Sample::static_size; ++i)
        {
            s_SampleInfos[i] = CreateSampleInfo(static_cast<Sample::value_type>(i));
        }
    }

    // get the current pset
    /*static*/ const std::string& Sample::GetPsetPath()
    {
        return sample_pset_path;
    }

    // operators:
    bool operator < (const Sample::Info& s1, const Sample::Info& s2)
    {
        return (s1.sample < s2.sample);
    }

    // Get the Sample from a string
    Sample::value_type GetSampleFromName(const std::string& sample_name)
    {
        for (const auto& sample_info : s_SampleInfos)
        {
            if (sample_info.name == sample_name)
            {
                return sample_info.sample; 
            }
        }

        // throw if not found
        throw std::domain_error(Form("[dy::GetSampleInfo] Error: sample %s not found!", sample_name.c_str()));
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
        for (const auto& sample_info : s_SampleInfos)
        {
            if (sample_info.name == sample_name)
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

    // get map of Sample::Info's
    std::map<Sample::value_type, Sample::Info> GetSampleMap()
    {
        std::map<Sample::value_type, Sample::Info> result; 
        for (int sample_num = 0; sample_num < Sample::static_size; ++sample_num)
        {
            const Sample::Info sample_info = GetSampleInfo(sample_num);
            result[sample_info.sample] = sample_info;
        }
        return result;
    }

} // namespace dy
 
