#include <fstream>
#include <sstream>
#include <iostream>
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"

std::string GetLatex(const std::string& title, const dy::Yield& yield, const bool data = false)
{

    std::string result;
    if (data)
    {
        result = Form
        (
            "%35s & %5.0f & %5.0f & %5.0f", 
            title.c_str(),
            yield.ee.value,
            yield.mm.value,
            yield.ll.value
        );
    }
    else
    {
        result = Form
        (
            "%35s & %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f", 
            title.c_str(),
            yield.ee.value,
            yield.ee.error,
            yield.mm.value,
            yield.mm.error,
            yield.ll.value,
            yield.ll.error
        );
    }
    return result;
}

std::string GetLatex
(
    const dy::Sample::value_type sample,
    const dy::YieldMap&  ym,
    const bool data = false
)
{
    return GetLatex(GetSampleInfo(sample).latex, ym.at(sample), data);
}

// print the yields
void CreateYieldTable 
(
    const std::string& label, 
    const std::string& hist_name = "h_reco_yield",
    const std::string& output_file = "", 
    bool print_latex = false
)
{
    // map of samples and yields
    dy::YieldMap  ym    = dy::GetYieldMap(label, hist_name);
    dy::Yield bkgd_pred = dy::GetBackgroundPred(label, hist_name);
    dy::Yield dy_pred   = ym[dy::Sample::dyll] + bkgd_pred;

    std::string table;
    if (print_latex)
    {
        string latex("\\begin{table}[ht!]\n"                            );
        latex.append("\\begin{center}\n"                                );
        latex.append("\\begin{tabular}{l|cccc} \\hline\\hline\n"        );
        latex.append("source & $ee$ & $\\mu\\mu$ & $\\ell\\ell $ \\\\\n");
        latex.append("\\hline\n"                                                            );
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::dyll       , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::wjets      , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::ttdil      , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::ttslq      , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::tthad      , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::qcdmu15    , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::ww2l2nu    , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::wz2l2q     , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::wz3lnu     , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::zz2l2nu    , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::zz2l2q     , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::zz4l       , ym       ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("Background Prediction", bkgd_pred).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("MC Prediction"        , dy_pred  ).c_str()));
        latex.append("\\hline\\hline\n");
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::data, ym, /*data=*/true).c_str()));
        latex.append("\\hline\\hline\n"                      );
        latex.append("\\end{tabular}\n"                      );
        latex.append("\\caption{Drell-Yan Exercise Yields}\n");
        latex.append("\\end{center}\n"                       );
        latex.append("\\end{table}"                          );

        // print it
        table = latex; 
    }
    else
    {
        // make the table
        CTable t_yields;
        t_yields.useTitle();
        t_yields.setTitle("yields for Drell-Yan Exercise");
        t_yields.setTable()
        (                                                                                "ee",                                   "mm",                                  "ll")
        (dy::GetSampleInfo(dy::Sample::dyll    ).name, ym[dy::Sample::dyll    ].ee.pm("4.1") , ym[dy::Sample::dyll    ].mm.pm("4.1") , ym[dy::Sample::dyll    ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::wjets   ).name, ym[dy::Sample::wjets   ].ee.pm("4.1") , ym[dy::Sample::wjets   ].mm.pm("4.1") , ym[dy::Sample::wjets   ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::ttdil   ).name, ym[dy::Sample::ttdil   ].ee.pm("4.1") , ym[dy::Sample::ttdil   ].mm.pm("4.1") , ym[dy::Sample::ttdil   ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::ttslq   ).name, ym[dy::Sample::ttslq   ].ee.pm("4.1") , ym[dy::Sample::ttslq   ].mm.pm("4.1") , ym[dy::Sample::ttslq   ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::tthad   ).name, ym[dy::Sample::tthad   ].ee.pm("4.1") , ym[dy::Sample::tthad   ].mm.pm("4.1") , ym[dy::Sample::tthad   ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::qcdmu15 ).name, ym[dy::Sample::qcdmu15 ].ee.pm("4.1") , ym[dy::Sample::qcdmu15 ].mm.pm("4.1") , ym[dy::Sample::qcdmu15 ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::ww2l2nu ).name, ym[dy::Sample::ww2l2nu ].ee.pm("4.1") , ym[dy::Sample::ww2l2nu ].mm.pm("4.1") , ym[dy::Sample::ww2l2nu ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::wz2l2q  ).name, ym[dy::Sample::wz2l2q  ].ee.pm("4.1") , ym[dy::Sample::wz2l2q  ].mm.pm("4.1") , ym[dy::Sample::wz2l2q  ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::wz3lnu  ).name, ym[dy::Sample::wz3lnu  ].ee.pm("4.1") , ym[dy::Sample::wz3lnu  ].mm.pm("4.1") , ym[dy::Sample::wz3lnu  ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::zz2l2nu ).name, ym[dy::Sample::zz2l2nu ].ee.pm("4.1") , ym[dy::Sample::zz2l2nu ].mm.pm("4.1") , ym[dy::Sample::zz2l2nu ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::zz2l2q  ).name, ym[dy::Sample::zz2l2q  ].ee.pm("4.1") , ym[dy::Sample::zz2l2q  ].mm.pm("4.1") , ym[dy::Sample::zz2l2q  ].ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::zz4l    ).name, ym[dy::Sample::zz4l    ].ee.pm("4.1") , ym[dy::Sample::zz4l    ].mm.pm("4.1") , ym[dy::Sample::zz4l    ].ll.pm("4.1"))
        ("Background Pred"                           ,                bkgd_pred.ee.pm("4.1") ,                bkgd_pred.mm.pm("4.1") ,                bkgd_pred.ll.pm("4.1"))
        ("MC Pred"                                   ,                  dy_pred.ee.pm("4.1") ,                  dy_pred.mm.pm("4.1") ,                  dy_pred.ll.pm("4.1"))
        (dy::GetSampleInfo(dy::Sample::data    ).name, ym[dy::Sample::data    ].ee.pm("4.0") , ym[dy::Sample::data    ].mm.pm("4.0") , ym[dy::Sample::data    ].ll.pm("4.0"))
        ;

        // print it
        std::ostringstream os;
        os << t_yields;
        table = os.str();
    }

    // output
    if (output_file.empty())
    {
        std::cout << table << std::endl;
    }
    else
    {
        std::ofstream fout(output_file);
        fout << table << std::endl;
    }
}
