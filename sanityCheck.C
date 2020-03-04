{
    TFile* f2 = TFile::Open("root://cmseos.fnal.gov//store/user/chpapage/DeadCellsSamples/TrainingSamples/out_E0to120Eta1p7_df01_0.root");
    TFile* f1 = TFile::Open("root://cmseos.fnal.gov//store/user/chpapage/SingleGamma_E0to120Eta1p7/SingleGamma_E0to120Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/200106_015819/0000/ntuples_pt5_1.root");

    TTree* t1 = dynamic_cast< TTree* >(f1->Get("hgcalTupleTree/tree"));
    TTree* t2 = dynamic_cast< TTree* >(f2->Get("t1"));

    float rechit, waferU, waferV, cellU, cellV, event, layer;
    float n1, n2, n3, n4, n5, n6;

    t2->SetBranchAddress("MLndown",&rechit);
    t2->SetBranchAddress("MLwaferU",&waferU);
    t2->SetBranchAddress("MLwaferV",&waferV);
    t2->SetBranchAddress("MLcellU" ,&cellU);
    t2->SetBranchAddress("MLcellV" ,&cellV);
    t2->SetBranchAddress("MLevent",&event);
    t2->SetBranchAddress("MLlayer",&layer);
    t2->SetBranchAddress("MLdn1",&n1);
    t2->SetBranchAddress("MLdn2",&n2);
    t2->SetBranchAddress("MLdn3",&n3);
    t2->SetBranchAddress("MLdn4",&n4);
    t2->SetBranchAddress("MLdn5",&n5);
    t2->SetBranchAddress("MLdn6",&n6);

    int n = t2->GetEntries();

    TCanvas* c = new TCanvas("c","c",1200,600);
    c->Divide(2,1);
    c->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.15);
    TH2Poly* h_honeycomb = new TH2Poly();
    h_honeycomb->Honeycomb(0,0,1,5,5);
    h_honeycomb->SetTitle("Recovered RecHits;Energy[GeV]");
    h_honeycomb->GetXaxis()->SetLabelOffset(999);
    h_honeycomb->GetXaxis()->SetLabelSize(0);
    h_honeycomb->GetXaxis()->SetTickLength(0);
    h_honeycomb->GetYaxis()->SetLabelOffset(999);
    h_honeycomb->GetYaxis()->SetLabelSize(0);
    h_honeycomb->GetYaxis()->SetTickLength(0);

    for (int i = 436; i <= n; ++i) {
        t2->GetEntry(i);
        h_honeycomb->Fill(3.5,2.5,n1);
        h_honeycomb->Fill(2.5,4,n3);
        h_honeycomb->Fill(3.5,5.5,n5);
        h_honeycomb->Fill(5,5.5,n6);
        h_honeycomb->Fill(6,4,n4);
        h_honeycomb->Fill(5,2.5,n2);
        h_honeycomb->Fill(5,4,rechit);
        break;
    }
    h_honeycomb->Draw("colz 0 text");

    cout << event << endl;
    cout << layer << endl;
    cout << waferU << endl;
    cout << waferV << endl;
    cout << cellU << endl;
    cout << cellV << endl;
    //std::ostringstream s1;
    //s1 << "HGCRecHitWaferU == " << waferU;

    TCut cut0 = "HGCRecHitIndex == 0";
    //char* tab1 = new char [s1.length()+1];
    //strcpy(tab1, s1.c_str());
    TCut cut1 = "HGCRecHitWaferU == -5";//tab1;
    TCut cut2 = "HGCRecHitWaferV == -9";
    TCut cut3 = "HGCRecHitLayer == 7";
    TCut cut4 = "event == 92";
    c->cd(2);
    gPad->SetRightMargin(0.17);
    t1->Draw(
        "HGCRecHitCellU:HGCRecHitCellV",
        ("HGCRecHitEnergy")*(cut0&&cut1&&cut2&&cut3&&cut4),
        "colz 0 text"
    );
}
