#include <fstream>
#include <sstream>
#include <iostream>
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"

std::string GetLatex(const std::string& title, const dy::Yield& yield)
{
    const std::string result = Form
    (
        "%35s & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f", 
        title.c_str(),
        yield.ee.value,
        yield.ee.error,
        yield.mm.value,
        yield.mm.error,
        yield.ll.value,
        yield.ll.error
    );
    return result;
}

std::string GetLatex
(
    const dy::Sample::value_type sample,
    const dy::YieldMap&  ym,
    const dy::SampleMap& sm
)
{
    return GetLatex(sm.at(sample).latex, ym.at(sample));
}

// print the yields
void CreateYieldTable 
(
    const std::string& label, 
    const std::string& output_file = "", 
    bool print_latex = false
)
{
    // map of samples and yields
    dy::YieldMap  ym    = dy::GetRecoYieldMap(label);
    dy::SampleMap sm    = dy::GetSampleMap();
    dy::Yield bkgd_pred = dy::GetBackgroundPred(label);
    dy::Yield dy_pred   = ym[dy::Sample::dyll] - bkgd_pred;

    std::string table;
    if (print_latex)
    {
        string latex("\\begin{table}[ht!]\n"                            );
        latex.append("\\begin{center}\n"                                );
        latex.append("\\begin{tabular}{l|cccc} \\hline\\hline\n"        );
        latex.append("source & $ee$ & $\\mu\\mu$ & $\\ell\\ell $ \\\\\n");
        latex.append("\\hline\n"                                                            );
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::dyll       , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::wjets      , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::ttdil      , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::ttslq      , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::tthad      , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::qcdmu15    , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::ww2l2nu    , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::wz2l2q     , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::wz3lnu     , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::zz2l2nu    , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::zz2l2q     , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::zz4l       , ym, sm   ).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("Background Prediction", bkgd_pred).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("MC Prediction"        , dy_pred  ).c_str()));
        latex.append("\\hline\\hline\n");
        latex.append(Form("%s \\\\\n", GetLatex(dy::Sample::data, ym, sm).c_str()));
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
        (                                                            "ee",                              "mm",                             "ll")
        (sm[dy::Sample::dyll    ].name, ym[dy::Sample::dyll    ].ee.pm() , ym[dy::Sample::dyll    ].mm.pm() , ym[dy::Sample::dyll    ].ll.pm())
        (sm[dy::Sample::wjets   ].name, ym[dy::Sample::wjets   ].ee.pm() , ym[dy::Sample::wjets   ].mm.pm() , ym[dy::Sample::wjets   ].ll.pm())
        (sm[dy::Sample::ttdil   ].name, ym[dy::Sample::ttdil   ].ee.pm() , ym[dy::Sample::ttdil   ].mm.pm() , ym[dy::Sample::ttdil   ].ll.pm())
        (sm[dy::Sample::ttslq   ].name, ym[dy::Sample::ttslq   ].ee.pm() , ym[dy::Sample::ttslq   ].mm.pm() , ym[dy::Sample::ttslq   ].ll.pm())
        (sm[dy::Sample::tthad   ].name, ym[dy::Sample::tthad   ].ee.pm() , ym[dy::Sample::tthad   ].mm.pm() , ym[dy::Sample::tthad   ].ll.pm())
        (sm[dy::Sample::qcdmu15 ].name, ym[dy::Sample::qcdmu15 ].ee.pm() , ym[dy::Sample::qcdmu15 ].mm.pm() , ym[dy::Sample::qcdmu15 ].ll.pm())
        (sm[dy::Sample::ww2l2nu ].name, ym[dy::Sample::ww2l2nu ].ee.pm() , ym[dy::Sample::ww2l2nu ].mm.pm() , ym[dy::Sample::ww2l2nu ].ll.pm())
        (sm[dy::Sample::wz2l2q  ].name, ym[dy::Sample::wz2l2q  ].ee.pm() , ym[dy::Sample::wz2l2q  ].mm.pm() , ym[dy::Sample::wz2l2q  ].ll.pm())
        (sm[dy::Sample::wz3lnu  ].name, ym[dy::Sample::wz3lnu  ].ee.pm() , ym[dy::Sample::wz3lnu  ].mm.pm() , ym[dy::Sample::wz3lnu  ].ll.pm())
        (sm[dy::Sample::zz2l2nu ].name, ym[dy::Sample::zz2l2nu ].ee.pm() , ym[dy::Sample::zz2l2nu ].mm.pm() , ym[dy::Sample::zz2l2nu ].ll.pm())
        (sm[dy::Sample::zz2l2q  ].name, ym[dy::Sample::zz2l2q  ].ee.pm() , ym[dy::Sample::zz2l2q  ].mm.pm() , ym[dy::Sample::zz2l2q  ].ll.pm())
        (sm[dy::Sample::zz4l    ].name, ym[dy::Sample::zz4l    ].ee.pm() , ym[dy::Sample::zz4l    ].mm.pm() , ym[dy::Sample::zz4l    ].ll.pm())
        ("Background Pred"            ,                bkgd_pred.ee.pm() ,                bkgd_pred.mm.pm() ,                bkgd_pred.ll.pm())
        ("MC Pred"                    ,                  dy_pred.ee.pm() ,                  dy_pred.mm.pm() ,                  dy_pred.ll.pm())
        (sm[dy::Sample::data    ].name, ym[dy::Sample::data    ].ee.pm() , ym[dy::Sample::data    ].mm.pm() , ym[dy::Sample::data    ].ll.pm())
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
