void rechitSum() {
    TCanvas* c = new TCanvas("c","c",1);
    TFile* f = TFile::Open("SingleGammaPt100Eta1p7.root");
    TTree* t = (TTree*)f->Get("hgcalTupleTree/tree");

    TFile* file_in = TFile::Open("/Users/chrispap/Documents/PFCal/PFCalEE/analysis/rootFiles/et100_eta1.7_test_deadfrac050.root");
    TH1F* AAA = static_cast<TH1F*>(file_in->Get("h_rechitsum")->Clone());
    AAA->Rebin(3);
    AAA->GetXaxis()->SetRangeUser(200,350);
    AAA->SetNormFactor(1);
    //AAA->SetMaximum(1);

    std::vector<float> *rechit = 0;
    std::vector<float> *rechiteta = 0;
    std::vector<float> *rechitphi = 0;
    std::vector<float> *rechitz = 0;
    std::vector<float> *rechitx = 0;
    std::vector<float> *rechity = 0;
    std::vector<int> *rechitlayer = 0;
    std::vector<float> *etagen = 0;
    std::vector<float> *phigen = 0;
    std::vector<float> *ptgen = 0;
    std::vector<float> *pgen = 0;

    t->SetBranchAddress("HGCRecHitEnergy",&rechit);
    t->SetBranchAddress("HGCRecHitEta",&rechiteta);
    t->SetBranchAddress("HGCRecHitPhi",&rechitphi);
    t->SetBranchAddress("HGCRecHitLayer",&rechitlayer);
    t->SetBranchAddress("HGCRecHitPosz",&rechitz);
    t->SetBranchAddress("HGCRecHitPosx",&rechitx);
    t->SetBranchAddress("HGCRecHitPosy",&rechity);
    t->SetBranchAddress("GenParEta",&etagen);
    t->SetBranchAddress("GenParPhi",&phigen);
    t->SetBranchAddress("GenParP",&pgen);

    TRandom3 r(0);

    TH1F* h1 = new TH1F("h1","sum of HGCRecHitEnergy (GenParEta > 1.5 && #Delta R<0.2 && layer < 35);rechitsum [GeV]",50,200,350);
    TH1F* h2 = new TH1F("h2","sum of HGCRecHitEnergy (GenParEta > 1.5 && #Delta R<0.2 && layer < 35) dead Si;rechitsum [GeV]",50,200,350);
    TH1F* h3 = new TH1F("h3","GenParP (GenParEta > 1.5 && #Delta R<0.2 && layer < 35);rechitsum [GeV]",50,150,400);
    //h2->SetMaximum(1);

    for (long i = 1; i <= t->GetEntries(); ++i) {
        t->GetEntry(i-1);
        float rechitsum = 0;
        float rechitsumDeadSi = 0;
        for (int j = 0; j<rechit->size();j++){
            if(rechitz->at(j)>0) {
                float deltaR = sqrt(pow(etagen->at(0)-rechiteta->at(j),2)+pow(phigen->at(0)-rechitphi->at(j),2));
                float thetagen = 2.*atan(exp(-etagen->at(0)));
                float rgen = rechitz->at(j)*tan(thetagen);
                float xgen = rgen*cos(phigen->at(0));
                float ygen = rgen*sin(phigen->at(0));
                float dR1 = fabs(sqrt((xgen-rechitx->at(j))*(xgen-rechitx->at(j))+(ygen-rechity->at(j))*(ygen-rechity->at(j))));
                if(deltaR < 0.3 && dR1 < 5.3) {
                    rechitsum+=rechit->at(j);
                    rechitsumDeadSi+=rechit->at(j);
                    if (r.Rndm() < 0.05) rechitsumDeadSi-=rechit->at(j);
                }
            }
        }
        if (rechitsum > 1) {
            h1->Fill(rechitsum);
            h2->Fill(rechitsumDeadSi);
        }
        h3->Fill(pgen->at(0));
    }
    h2->SetNormFactor(1);
    h1->SetNormFactor(1);
    THStack* hs = new THStack("hs","Energy;E [GeV]");
    hs->Add(AAA);
    hs->Add(h2);
    h2->SetLineColor(2);
    //h2->Draw();
    //c->SetLogy();

    TF1* f1 = new TF1("f1","gaus");
    TF1* f2 = new TF1("f2","gaus");
    f1->SetNpx(200);
    f2->SetNpx(200);
    f1->SetLineColor(4);

    //AAA->Fit("f1","");
//    h2->Fit("f2","");

    //AAA->Fit("f1","","",f1->GetParameter(1)-1*f1->GetParameter(2),f1->GetParameter(1)+3*f1->GetParameter(2));
    //h2->Fit("f2","","",f2->GetParameter(1)-1*f2->GetParameter(2),f2->GetParameter(1)+3*f2->GetParameter(2));


    //hs->Draw("nostack");
    //hs->GetXaxis()->SetRangeUser(200,350);
    //hs->SetMaximum(1);
    h1->SetLineColor(2);
    AAA->Draw();
    h1->Draw("same");
    c->Update();
}
