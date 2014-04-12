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
            "%35s & %5.0f & %5.0f", 
            title.c_str(),
            yield.ee.value,
            yield.mm.value
        );
    }
    else
    {
        result = Form
        (
            "%35s & %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f", 
            title.c_str(),
            yield.ee.value,
            yield.ee.error,
            yield.mm.value,
            yield.mm.error
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
void CreateCrossSectionTable 
(
    const std::string& label, 
    const std::string& hist_name = "h_reco_yield",
    const std::string& output_file = "", 
    bool print_latex = false
)
{
    const double lumi = 0.082; // fb^-1

    // map of samples and yields
    dy::YieldMap  ym    = dy::GetYieldMap(label, hist_name);
    dy::Yield y_data    = ym[dy::Sample::data];
    dy::Yield y_bk_pred = dy::GetBackgroundPred(label, hist_name);
    dy::Yield y_dy_pred = ym[dy::Sample::dyll];
    dy::Yield y_mc_pred = y_dy_pred + y_bk_pred;
    dy::Yield y_acc     = dy::GetYieldFromLabel(dy::Sample::dyll, label, "h_acc");
    dy::Yield y_nsig    = y_data - y_bk_pred; 

    // xsec = (Nobs - Nbkg)/(lumi * acc)
    dy::Yield xsec = y_nsig/(lumi * y_acc); 
    xsec = xsec * (1e-6); // fb --> nb


    std::string result;
    if (print_latex)
    {
        string latex("\\begin{table}[ht!]\n"                    ); 
        latex.append("\\begin{center}\n"                        ); 
        latex.append("\\begin{tabular}{l|ccc} \\hline\\hline\n" ); 
        latex.append("source & $ee$ & $\\mu\\mu$$ \\\\\n"       ); 
        latex.append("\\hline\n"                                ); 
        latex.append(Form("%s \\\\\n", GetLatex("N_obs"     ,    y_data).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("N_mc"      , y_mc_pred).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("N_dyll"    , y_dy_pred).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("N_bkgd"    , y_bk_pred).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("N_sig"     ,    y_nsig).c_str()));
        latex.append(Form("%s \\\\\n", GetLatex("Acc"       ,     y_acc).c_str()));
        latex.append("\\hline\\hline\n");
        latex.append(Form("%s \\\\\n", GetLatex("Sigma (nb)",      xsec).c_str()));
        latex.append("\\hline\\hline\n"                              ); 
        latex.append("\\end{tabular}\n"                              ); 
        latex.append("\\caption{Drell-Yan Exercise Cross-Section}\n" ); 
        latex.append("\\end{center}\n"                               ); 
        latex.append("\\end{table}"                                  ); 
        result = latex;
    }
    else
    {
        // make the table
        CTable t_yields;
        t_yields.useTitle();
        t_yields.setTitle("Cross-Sectoin for Drell-Yan Exercise");
        t_yields.setTable()
        (                                   "ee",                    "mm") 
        ("N_obs"       ,    y_data.ee.pm("4.0") ,     y_data.mm.pm("4.0")) 
        ("N_mc"        , y_mc_pred.ee.pm("4.1") ,  y_mc_pred.mm.pm("4.1")) 
        ("N_dyll"      , y_dy_pred.ee.pm("4.1") ,  y_dy_pred.mm.pm("4.1")) 
        ("N_bkgd"      , y_bk_pred.ee.pm("4.1") ,  y_bk_pred.mm.pm("4.1")) 
        ("N_sig"       ,    y_nsig.ee.pm("4.1") ,     y_nsig.mm.pm("4.1")) 
        ("Acc"         ,     y_acc.ee.pm("4.3") ,      y_acc.mm.pm("4.3")) 
        ("Sigma (nb)"  ,      xsec.ee.pm("4.3") ,       xsec.mm.pm("4.3")) 
        ;

        // print it
        std::ostringstream os;
        os << t_yields;
        result = os.str();
    }

    // output
    if (output_file.empty())
    {
        std::cout << result << std::endl;
    }
    else
    {
        std::ofstream fout(output_file);
        fout << result << std::endl;
    }
}
