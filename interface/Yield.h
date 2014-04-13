#ifndef DY_YIELD_H
#define DY_YIELD_H

#include <string> 
#include <map> 
#include <vector> 
#include "TH1.h"
#include "Analysis/DrellYan/interface/Sample.h"

namespace dy
{
    // simple struct to hold the yield
    class Yield
    {
        public:
            class value_t
            {
                public:
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
    Yield operator*(const Yield& y1, const Yield& y2);
    Yield operator/(const Yield& y1, const Yield& y2);

    Yield operator*(const double scale, const Yield& y);
    Yield operator*(const Yield& y, const double scale);

    // get yields from a yield histogram
    // NOTE: expects 3 bins with bin 0 --> ll, 1 --> mm, 2 --> ee
    Yield GetYieldFromHist(TH1& hist);
    Yield GetYieldFromLabel(const Sample::value_type, const std::string& label, const std::string& hist_name = "h_reco_full_yield");

    // get total background prediction
    Yield GetBackgroundPred(const std::string& label, const std::string& hist_name = "h_reco_full_yield");

    // human readable yield table
    std::string GetYieldString(const Yield& yield, const std::string& title = "", const std::string& fmt = "4.1");

    // yield array for all samples 
    typedef std::vector<dy::Yield> YieldVector;
    YieldVector GetYieldVector(const std::string& label, const std::string& hist_name = "h_reco_full_yield");

} // namepsace dy

#endif // DY_YIELD_H
