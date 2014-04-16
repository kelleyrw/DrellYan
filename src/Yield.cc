#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include <sstream>
#include <cassert>

namespace dy 
{
    // methods:
    Yield& Yield::operator+=(const Yield& y)
    {
        Yield temp_yield = (*this) + y;
        *this = temp_yield;
        return *this; 
    }

    Yield& Yield::operator-=(const Yield& y)
    {
        Yield temp_yield = (*this) - y;
        *this = temp_yield;
        return *this; 
    }

    // non-member methods
    Yield operator*(const double scale, const Yield& y)
    {
        Yield result =
        {
            {scale * y.ee.value, scale * y.ee.error}, 
            {scale * y.mm.value, scale * y.mm.error}, 
            {scale * y.ll.value, scale * y.ll.error}
        };
        return result; 
    }
    
    Yield operator*(const Yield& y, const double scale)
    {
        return scale * y;
    }
    
    Yield operator-(const Yield& y1, const Yield& y2)
    {
        Yield result =
        {
            {y1.ee.value - y2.ee.value, sqrt(y1.ee.error*y1.ee.error + y2.ee.error*y2.ee.error)}, 
            {y1.mm.value - y2.mm.value, sqrt(y1.mm.error*y1.mm.error + y2.mm.error*y2.mm.error)}, 
            {y1.ll.value - y2.ll.value, sqrt(y1.ll.error*y1.ll.error + y2.ll.error*y2.ll.error)}
        };
        return result; 
    }
    
    Yield operator+(const Yield& y1, const Yield& y2)
    {
        Yield result =
        {
            {y1.ee.value + y2.ee.value, sqrt(y1.ee.error*y1.ee.error + y2.ee.error*y2.ee.error)}, 
            {y1.mm.value + y2.mm.value, sqrt(y1.mm.error*y1.mm.error + y2.mm.error*y2.mm.error)}, 
            {y1.ll.value + y2.ll.value, sqrt(y1.ll.error*y1.ll.error + y2.ll.error*y2.ll.error)}
        };
        return result; 
    }

    Yield operator*(const Yield& y1, const Yield& y2)
    {
        Yield result =
        {
            {y1.ee.value * y2.ee.value, sqrt((y1.ee.error*y1.ee.error)/(y1.ee.value*y2.ee.value) + (y2.ee.error*y2.ee.error)/(y1.ee.value*y2.ee.value))}, 
            {y1.mm.value * y2.mm.value, sqrt((y1.mm.error*y1.mm.error)/(y1.mm.value*y2.mm.value) + (y2.mm.error*y2.mm.error)/(y1.mm.value*y2.mm.value))}, 
            {y1.ll.value * y2.ll.value, sqrt((y1.ll.error*y1.ll.error)/(y1.ll.value*y2.ll.value) + (y2.ll.error*y2.ll.error)/(y1.ll.value*y2.ll.value))}
        };
        return result; 
    }

    Yield operator/(const Yield& y1, const Yield& y2)
    {
        Yield::value_t ee =
        {
            lt::is_zero(y2.ee.value) ? 0.0 : y1.ee.value / y2.ee.value, 
            lt::is_zero(y2.ee.value) ? 0.0 : sqrt((y1.ee.error*y1.ee.error)/(y1.ee.value * y2.ee.value) + (y2.ee.error*y2.ee.error)/(y1.ee.value * y2.ee.value))
        }; 
        Yield::value_t mm =
        {
            lt::is_zero(y2.mm.value) ? 0.0 : y1.mm.value / y2.mm.value, 
            lt::is_zero(y2.mm.value) ? 0.0 : sqrt((y1.mm.error*y1.mm.error)/(y1.mm.value * y2.mm.value) + (y2.mm.error*y2.mm.error)/(y1.mm.value * y2.mm.value))
        }; 
        Yield::value_t ll =
        {
            lt::is_zero(y2.ll.value) ? 0.0 : y1.ll.value / y2.ll.value, 
            lt::is_zero(y2.ll.value) ? 0.0 : sqrt((y1.ll.error*y1.ll.error)/(y1.ll.value * y2.ll.value) + (y2.ll.error*y2.ll.error)/(y1.ll.value * y2.ll.value))
        }; 
        Yield result = {ee, mm, ll};
        return result; 
    }

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
        assert(hist.GetNbinsX() == 4);
        auto ll = rt::IntegralAndError(&hist, 0.0, 1.0);
        auto mm = rt::IntegralAndError(&hist, 1.0, 2.0);
        auto ee = rt::IntegralAndError(&hist, 2.0, 3.0);
        Yield result = 
        {
            /*ee = */{ee.first, ee.second},
            /*mm = */{mm.first, mm.second},
            /*ll = */{ll.first, ll.second},
        };
        return result;
    }

    Yield GetYieldFromLabel(const Sample::value_type sample, const std::string& label, const std::string& hist_name)
    { 
        if (label.empty())
        {
            throw std::runtime_error("[dy::GetYieldFromFile] Error: label is empty");
        }
        if (!lt::file_exists("plots/"+label))
        {
            throw std::runtime_error(Form("[dy::GetYieldFromFile] Error: label %s does not exist", label.c_str()));
        }
        const Sample::Info sample_info = GetSampleInfo(sample);
        const std::string file_name = Form("plots/%s/%s_plots.root", label.c_str(), sample_info.name.c_str());
        TH1 * const yield_hist = rt::GetHistFromRootFile<TH1>(file_name, hist_name);
        Yield result = GetYieldFromHist(*yield_hist);
        return result;
    }

    Yield GetBackgroundPred(const std::string& label, const std::string& hist_name)
    {
        const YieldVector yields = GetYieldVector(label, hist_name);
        Yield pred{{0, 0},{0,0},{0.0}};
        for (size_t i = 0; i < yields.size(); ++i)
        {
            const Sample::value_type sample = GetSampleFromNumber(i);
            if (sample != Sample::data && sample != Sample::dyll)
            {
                pred += yields[sample];
            }
        }
        return pred;
    }

    // human readable yield table
    std::string GetYieldString(const Yield& yield, const std::string& title, const std::string& fmt)
    {
        CTable table;
        table.setTitle(title.empty() ? "Drell-Yan Analysis yield" : title);
        table.useTitle();
        table.setTable() (                     "ee",             "mm",             "ll")
                         ("yield", yield.ee.pm(fmt), yield.mm.pm(fmt), yield.ll.pm(fmt))
                         ;
        std::ostringstream os;
        os << table;
        return os.str();
    }

    // get yield map for all samples 
    YieldVector GetYieldVector(const std::string& label, const std::string& hist_name)
    {
        std::vector<Yield> result; 
        result.reserve(Sample::static_size);
        for (const auto& si : Sample::GetInfos())
        {
            const Yield sample_yield = GetYieldFromLabel(si.sample, label, hist_name);
            result.push_back(sample_yield);
        }
        return result;
    }

} // namespace dy 
