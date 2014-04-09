#include <fstream>
#include <sstream>
#include <iostream>
#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"

// print the yields
void CreateYieldTable 
(
    const std::string& label, 
    bool print_latex = false
)
{
    // map of samples and yields
    dy::YieldMap  ym = dy::GetRecoYieldMap(label);
    dy::SampleMap sm = dy::GetSampleMap();

    // total prediction
    dy::Yield pred;
    for (const auto& s : ym)
    {
        if (s.first != dy::Sample::data) pred += s.second;
    }

    std::string table;
    if (print_latex)
    {
        string latex;
/*         // before table */
/*         string latex("\\begin{table}[ht!]\n"                                       );  */
/*         latex.append("\\begin{center}\n"                                           );  */
/*         latex.append("\\begin{tabular}{l|cccc} \\hline\\hline\n"                   );  */
/*         latex.append(do_caption ? Form("\\caption{%s}\n", sr_info.latex.c_str()) : ""); */
/*         latex.append("source & $ee$ & $\\mu\\mu$ & $\\ell\\ell $ \\\\\n" );  */
/*         latex.append("\\hline\n"                                                   );  */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append("\\hline\\hline\n"); */
/*         latex.append(Form("%s \\\\\n", ym[dy::Sample::data].GetLatexLine().c_str())); */
/*         latex.append("\\hline\\hline\n" );  */
/*         latex.append("\\end{tabular}\n" );  */
/*         latex.append("\\end{center}\n"  );  */
/*         latex.append("\\end{table}"     );  */
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
        ("Total Pred"                 ,                     pred.ee.pm() ,                     pred.mm.pm() ,                     pred.ll.pm())
        (sm[dy::Sample::data    ].name, ym[dy::Sample::data    ].ee.pm() , ym[dy::Sample::data    ].mm.pm() , ym[dy::Sample::data    ].ll.pm())
        ;

        // print it
        std::ostringstream os;
        os << t_yields;
        table = os.str();
    }

    std::cout << table << std::endl;
}
