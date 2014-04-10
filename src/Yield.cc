#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include <sstream>

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
        auto ee = rt::IntegralAndError(&hist, 2.0, 3.0);
        Yield result = 
        {
            /*ee = */{ee.first, ee.second},
            /*mm = */{mm.first, mm.second},
            /*ll = */{ll.first, ll.second},
        };
        return result;
    }

    Yield GetRecoYieldFromLabel(const Sample::value_type sample, const std::string& label)
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
        TH1 * const yield_hist = rt::GetHistFromRootFile<TH1>(file_name, "h_reco_yield");
        Yield result = GetYieldFromHist(*yield_hist);
        return result;
    }

    Yield GetGenYieldFromLabel(const Sample::value_type sample, const std::string& label)
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
        TH1 * const yield_hist = rt::GetHistFromRootFile<TH1>(file_name, "h_gen_yield");
        Yield result = GetYieldFromHist(*yield_hist);
        return result;
    }

    Yield GetBackgroundPred(const std::string& label)
    {
        const dy::YieldMap ym = dy::GetRecoYieldMap(label);
        dy::Yield pred;
        for (const auto& s : ym)
        {
            if (s.first != dy::Sample::data && s.first != dy::Sample::dyll)
            {
                pred += s.second;
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
    std::map<Sample::value_type, Yield> GetRecoYieldMap(const std::string& label)
    {
        std::map<Sample::value_type, Yield> result; 
        for (int sample_num = 0; sample_num < Sample::static_size; ++sample_num)
        {
            const Sample::Info sample_info = GetSampleInfo(sample_num);
            const Yield sample_yield = GetRecoYieldFromLabel(sample_info.sample, label);
            result[sample_info.sample] = sample_yield;
        }
        return result;
    }

    std::map<Sample::value_type, Yield> GetGenYieldMap(const std::string& label)
    {
        std::map<Sample::value_type, Yield> result; 
        for (int sample_num = 0; sample_num < Sample::static_size; ++sample_num)
        {
            const Sample::Info sample_info = GetSampleInfo(sample_num);
            const Yield sample_yield = GetGenYieldFromLabel(sample_info.sample, label);
            result[sample_info.sample] = sample_yield;
        }
        return result;
    }

} // namespace dy 
