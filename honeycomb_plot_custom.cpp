{
    // Import tree & variables
    TCanvas* c = new TCanvas("c","c",1100,400);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.04);
    TFile* f = TFile::Open("data/training_sample_full_8samples.root");
    TTree* t = dynamic_cast< TTree* >(f->Get("t1"));
    float n1, n2, n3, n4, n5, n6, dead, down, up, phi;
    float un1, un2, un3, un4, un5, un6;
    float dn1, dn2, dn3, dn4, dn5, dn6, layer;
    t->SetBranchAddress("MLlayer",&layer);
    t->SetBranchAddress("MLphi",&phi);
    t->SetBranchAddress("MLn1",&n1);
    t->SetBranchAddress("MLn2",&n2);
    t->SetBranchAddress("MLn3",&n3);
    t->SetBranchAddress("MLn4",&n4);
    t->SetBranchAddress("MLn5",&n5);
    t->SetBranchAddress("MLn6",&n6);
    t->SetBranchAddress("MLdead",&dead);
    t->SetBranchAddress("MLndown",&down);
    t->SetBranchAddress("MLnup",&up);
    t->SetBranchAddress("MLun1",&un1);
    t->SetBranchAddress("MLun2",&un2);
    t->SetBranchAddress("MLun3",&un3);
    t->SetBranchAddress("MLun4",&un4);
    t->SetBranchAddress("MLun5",&un5);
    t->SetBranchAddress("MLun6",&un6);
    t->SetBranchAddress("MLdn1",&dn1);
    t->SetBranchAddress("MLdn2",&dn2);
    t->SetBranchAddress("MLdn3",&dn3);
    t->SetBranchAddress("MLdn4",&dn4);
    t->SetBranchAddress("MLdn5",&dn5);
    t->SetBranchAddress("MLdn6",&dn6);

    // Make histo
    TH2Poly* h_honeycomb = new TH2Poly();
    h_honeycomb->Honeycomb(0,0,1,13,5);
    h_honeycomb->SetTitle("Rechits for high deviation case (>20 GeV)");
    h_honeycomb->GetXaxis()->SetLabelOffset(999);
    h_honeycomb->GetXaxis()->SetLabelSize(0);
    h_honeycomb->GetXaxis()->SetTickLength(0);
    h_honeycomb->GetYaxis()->SetLabelOffset(999);
    h_honeycomb->GetYaxis()->SetLabelSize(0);
    h_honeycomb->GetYaxis()->SetTickLength(0);

    float a = 0;
    int iC = 3;
    for (long i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if(dead-up/2-down/2>20) {
            if(iC){
                iC--;
                continue;
            }
            h_honeycomb->Fill(3.5+a,2.5,dn1);
            h_honeycomb->Fill(2.5+a,4,dn3);
            h_honeycomb->Fill(3.5+a,5.5,dn5);
            h_honeycomb->Fill(5+a,5.5,dn6);
            h_honeycomb->Fill(6+a,4,dn4);
            h_honeycomb->Fill(5+a,2.5,dn2);
            h_honeycomb->Fill(5+a,4,down);
            a = 7;
            h_honeycomb->Fill(3.5+a,2.5,n1);
            h_honeycomb->Fill(2.5+a,4,n3);
            h_honeycomb->Fill(3.5+a,5.5,n5);
            h_honeycomb->Fill(5+a,5.5,n6);
            h_honeycomb->Fill(6+a,4,n4);
            h_honeycomb->Fill(5+a,2.5,n2);
            h_honeycomb->Fill(5+a,4,dead);
            a = 14;
            h_honeycomb->Fill(3.5+a,2.5,un1);
            h_honeycomb->Fill(2.5+a,4,un3);
            h_honeycomb->Fill(3.5+a,5.5,un5);
            h_honeycomb->Fill(5+a,5.5,un6);
            h_honeycomb->Fill(6+a,4,un4);
            h_honeycomb->Fill(5+a,2.5,un2);
            h_honeycomb->Fill(5+a,4,up);
            break;
        }
    }
    h_honeycomb->GetZaxis()->SetTitle("Rechit [GeV]");
    h_honeycomb->GetZaxis()->SetTitleOffset(0.54);
    h_honeycomb->Draw("colz 0");

    // Add necessary text
    TLatex ltx;
    ltx.DrawLatex(3,7,"Previous layer");
    ltx.DrawLatex(10,7,"Dead cell layer");
    ltx.DrawLatex(17.5,7,"Next layer");

    // Add separating lines
    TLine l1(7.5,0,7.5,8);
    l1.Draw();
    TLine l2(15,0,15,8);
    l2.Draw();
    TMarker mr1(4.35,4,3);
    mr1.Draw();
    TMarker mr2(11.3,4,3);
    mr2.Draw();
    TMarker mr3(18.23,4,3);
    mr3.Draw();

    /* Add phi angle lines (deprecated)
    float ax = 0.866;
    float ay = 1.25;
    cout << "\nphi = " << phi*180/3.14 << endl;
    float f1 = 4.09;
    TLine lines1(5*ax-1.1*f1*sin(phi),4-1.1*f1*cos(phi),5*ax+0.8*f1*sin(phi),4+0.8*f1*cos(phi));
    lines1.SetLineColor(kRed);
    lines1.Draw();
    TLine lines2(5*ax-0.97*f1*sin(phi)+7,4-0.97*f1*cos(phi),5*ax+0.93*f1*sin(phi)+7,4+0.93*f1*cos(phi));
    lines2.SetLineColor(kRed);
    lines2.Draw();
    TLine lines3(5*ax-0.82*f1*sin(phi)+13.9,4-0.82*f1*cos(phi),(5*ax+1.08*f1*sin(phi))+13.9,4+1.08*f1*cos(phi));
    lines3.SetLineColor(kRed);
    lines3.Draw();
    char angletxt[15];
    sprintf(angletxt,"#phi = %d^{o}",(int)(round(phi*180/TMath::Pi())));
    TLatex anltx;
    anltx.SetTextAngle(15);
    anltx.DrawLatex(6,4.6,angletxt);
    anltx.DrawLatex(6.3+7,4.7,angletxt);
    anltx.DrawLatex(6.3+13.9,4.7,angletxt);*/
    }
    h1->Draw("colz");
}
