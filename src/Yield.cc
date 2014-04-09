#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include <sstream>

namespace dy 
{
    // simple struct to print +/-
    std::string Yield::value_t::pm(const std::string& fmt) const
    {
        return lt::pm(value, error, fmt);
    }

    // simple struct to print +/-
    std::string Yield::value_t::pp(const std::string& fmt) const
    {
        return Form("%s (%1.0f%%)", lt::pm(value, error, fmt).c_str(), 100.0*frac_error());
    }

    // calc the percent error 
    double Yield::value_t::frac_error() const
    {
        return error/value;
    }

    // get yields from a yield histogram
    // NOTE: expect 4 bins with bin 0 --> ll, 1 --> ee, 3 --> mm
    // (bin 2 is em which is unused)
    Yield GetYieldFromHist(TH1& hist)
    {
        auto ll = rt::IntegralAndError(&hist, 0.0, 1.0);
        auto mm = rt::IntegralAndError(&hist, 1.0, 2.0);
        auto ee = rt::IntegralAndError(&hist, 3.0, 4.0);
        Yield result = 
        {
            /*ee = */{ee.first, ee.second},
            /*mm = */{mm.first, mm.second},
            /*ll = */{ll.first, ll.second},
        };
        return result;
    }

    // human readable yield table
    std::string GetYieldString(const Yield& yield, const std::string& title)
    {
        CTable table;
        table.setTitle(title.empty() ? "Drell-Yan Analysis yield" : title);
        table.useTitle();
        table.setTable() (                  "ee",          "mm",          "ll")
                         ("yield", yield.ee.pm(), yield.mm.pm(), yield.ll.pm())
                         ;
        std::ostringstream os;
        os << table;
        return os.str();
    }

} // namespace dy 
