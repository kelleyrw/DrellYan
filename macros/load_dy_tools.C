{
    // load the AnalysisTools for interactive use 
    const TString cmssw_base = gSystem->Getenv("CMSSW_BASE");
    const bool cmssw = cmssw_base.Length() > 0;
    load_cmssw();
    gROOT->ProcessLine(".x $CMSSW_BASE/src/AnalysisTools/utils/load_analysis_tools.C");
    if (cmssw)
    {
        std::cout << "Loading DY tools" << std::endl;
        const TString scram_arch = gSystem->Getenv("SCRAM_ARCH");
        const TString cmssw_lib  = Form("%s/lib/%s", cmssw_base.Data(), scram_arch.Data());
        const TString dy_suffix  = (scram_arch.Contains("osx") ? "dylib" : "so");
        gSystem->AddIncludePath(Form("-I\"%s/src\"", cmssw_base.Data()));
        gSystem->AddIncludePath(Form("-I\"%s/src/Analysis/DrellYan/interface\"", cmssw_base.Data()));
        gSystem->Load(Form("%s/libAnalysisDrellYan.%s", cmssw_lib.Data(), dy_suffix.Data()));
    }
}
