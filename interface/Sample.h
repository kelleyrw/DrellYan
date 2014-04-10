#ifndef DY_SAMPLE_H
#define DY_SAMPLE_H

#include <string>
#include "TChain.h"
#include "TColor.h"

namespace dy
{
    // simple Sample class
    struct Sample
    {
        // list of available samples
        enum value_type
        {
            data,
            dyll,
            wjets,
            ttdil,
            ttslq,
            tthad,
            qcdmu15,
            ww2l2nu,
            wz2l2q,
            wz3lnu,
            zz2l2nu,
            zz2l2q,
            zz4l,

            static_size
        };

        // sample information
        struct Info
        {
            std::string name;        // short name
            std::string title;       // ROOT TLatex title
            std::string latex;       // real latex title
            std::string ntuple_path; // logical name for path
            Color_t color;           // color for plots
            double filter_eff;       // SD filter efficiency
            value_type sample;       // redundant process enum
        };
    };

    // operator
    bool operator < (const Sample::Info& s1, const Sample::Info& s2);

    // Get the Sample from a string/number
    Sample::value_type GetSampleFromName(const std::string& sample_name);
    Sample::value_type GetSampleFromNumber(const int sample_num);

    // test if a string is on of the samples
    bool IsSample(const std::string& sample_name);

    // wrapper function to get the Sample::Info
    Sample::Info GetSampleInfo(const Sample::value_type& sample);
    Sample::Info GetSampleInfo(const std::string& sample_name);
    Sample::Info GetSampleInfo(const int sample_num);

    // get the chain from the Sample
    TChain* GetSampleTChain(const Sample::value_type& sample); 

    // get map of Sample::Info's
    typedef std::map<Sample::value_type, Sample::Info> SampleMap;
    SampleMap GetSampleMap();

} // namespace dy

#endif //DY_SAMPLE_H
