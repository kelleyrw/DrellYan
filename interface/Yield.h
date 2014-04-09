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
    Yield GetRecoYieldFromLabel(const Sample::value_type, const std::string& label);
    Yield GetGenYieldFromLabel (const Sample::value_type, const std::string& label);

    // human readable yield table
    std::string GetYieldString(const Yield& yield, const std::string& title = "");

    // get yield map for all samples 
    typedef std::map<dy::Sample::value_type, dy::Yield> YieldMap;
    YieldMap GetRecoYieldMap(const std::string& label);
    YieldMap GetGenYieldMap (const std::string& label);

} // namepsace dy

#endif // DY_YIELD_H
