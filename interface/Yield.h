#ifndef DY_YIELD_H
#define DY_YIELD_H

#include <string> 
#include "TH1.h"

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
    };

    // get yields from a yield histogram
    // NOTE: expect 4 bins with bin 0 --> ll, 1 --> ee, 3 --> mm
    // (bin 2 is em which is unused)
    Yield GetYieldFromHist(TH1& hist);

    // human readable yield table
    std::string GetYieldString(const Yield& yield, const std::string& title = "");

} // namepsace dy

#endif // DY_YIELD_H
