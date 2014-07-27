{
    // data
    TChain chain("Events");
    chain.Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_1[0-9].root");
    const double nevts_cms2 = 27137253;
    const double nevts_file = chain.GetEntries(); 
    const double lumi       = 0.082;

    // selections
    TCut raw_ee    = "Sum$(genps_status==3 && genps_id==11)>=1 && Sum$(genps_status==3 && genps_id==-11)>=1";
    TCut scaled_ee = Form("%f*evt_scale1fb*(Sum$(genps_status==3 && genps_id==11)>=1 && Sum$(genps_status==3 && genps_id==-11)>=1)", lumi*nevts_cms2/nevts_file);
    cout << "raw_ee = " << raw_ee.GetTitle() << endl;
    cout << "scaled_ee = " << scaled_ee.GetTitle() << endl;

    // raw yields
    const double y_raw_ee = chain.GetEntries(raw_ee);
    cout << "raw ee yield = " << y_raw_ee << endl;

    // scaled yields
    TH1D* h_gen_yield = new TH1D("h_gen_yield", "h_gen_yield;m_{ee}(GeV);# events", 3, 0, 3);
    chain.Draw("1>>h_gen_yield", scaled_ee, "goff");
    h_gen_yield.Draw();
    cout << "scaled ee yield = " << h_gen_yield->Integral(0, -1) << endl;

    // another way to get scaled yields
    const double evt_scale1fb = chain.GetMaximum("evt_scale1fb");
    const double event_scale = lumi * evt_scale1fb * nevts_cms2/nevts_file;
    cout << "scaled ee yield (v2) = " << y_raw_ee * event_scale << endl;
}
