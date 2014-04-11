#ifndef DY_YIELD_H
#define DY_YIELD_H

#include <string> 
#include <map> 
#include "TH1.h"
#include "Analysis/DrellYan/interface/Sample.h"

namespace dy
{
    // simple struct to hold the yield
    struct Yield
    {
        struct value_t
        {
            // members:
            double value;
            double error;
        
            // methods:
            std::string pm(const std::string& fmt = "4.1") const;  // value +/- error
            std::string pp(const std::string& fmt = "4.1") const;  // value +/- error (perc error %)
            double frac_error() const;
        };

        value_t ee;
        value_t mm;
        value_t ll;

        // operators:
        Yield& operator+=(const Yield& y);
        Yield& operator-=(const Yield& y);
    };

    // non-members
    Yield operator-(const Yield& y1, const Yield& y2);
    Yield operator+(const Yield& y1, const Yield& y2);
    Yield operator*(const double scale, const Yield& y);

    // get yields from a yield histogram
    // NOTE: expect 4 bins with bin 0 --> ll, 1 --> ee, 3 --> mm
    // (bin 2 is em which is unused)
    Yield GetYieldFromHist(TH1& hist);
    Yield GetYieldFromLabel(const Sample::value_type, const std::string& label, const std::string& hist_name = "h_reco_yield");

    // get total background prediction
    Yield GetBackgroundPred(const std::string& label, const std::string& hist_name = "h_reco_yield");

    // human readable yield table
    std::string GetYieldString(const Yield& yield, const std::string& title = "", const std::string& fmt = "4.1");

    // get yield map for all samples 
    typedef std::map<dy::Sample::value_type, dy::Yield> YieldMap;
    YieldMap GetYieldMap(const std::string& label, const std::string& hist_name = "h_reco_yield");

} // namepsace dy

#endif // DY_YIELD_H
