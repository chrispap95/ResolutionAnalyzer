/****************************************************
**  Channelog:
**      Moved to CMSSW output.
****************************************************/

#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<map>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"
#ifdef _DEBUG
#include "debug_new.h"
#endif
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

using boost::lexical_cast;
namespace po=boost::program_options;

double DeltaR(double eta1,double phi1,double eta2,double phi2){
    double dr=99999.;
    double deta=fabs(eta1-eta2);
    double dphi=fabs(phi1-phi2);
    if(dphi>TMath::Pi()) dphi=2.*TMath::Pi()-dphi;
    dr=sqrt(deta*deta+dphi*dphi);
    return dr;
}

int main(int argc, char** argv){
    /**********************************
    ** initialize some variables
    **********************************/

    //Input output and config options
    std::string cfg;
    unsigned pNevts;
    std::string outFilePath;
    std::string filePath;
    std::string digifilePath;
    unsigned nRuns;
    std::string recoFileName;
    std::string MLFilePath;
    unsigned debug;
    double deadfrac;
    bool adjacent, MLsample;
    po::options_description preconfig("Configuration");
    preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
    po::notify(vm);
    po::options_description config("Configuration");
    config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",        po::value<unsigned>(&pNevts)->default_value(0))
    ("outFilePath,o",   po::value<std::string>(&outFilePath)->required())
    ("filePath,i",      po::value<std::string>(&filePath)->required())
    ("recoFileName,r",  po::value<std::string>(&recoFileName)->required())
    ("nRuns",           po::value<unsigned>(&nRuns)->default_value(0))
    ("debug,d",         po::value<unsigned>(&debug)->default_value(0))
    ("deadfrac",        po::value<double>(&deadfrac)->default_value(0))
    //Restrict number of adjacent dead cells
    ("adjacent",        po::value<bool>(&adjacent)->default_value(0))
    //Generate ML study training sample
    ("MLsample",        po::value<bool>(&MLsample)->default_value(1))
    //File to export data for ML
    ("MLFilePath",      po::value<std::string>(&MLFilePath)->default_value("training_sample.root"))

    ;
    po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
    po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
    po::notify(vm);

    std::cout << " -- Input parameters: " << std::endl
    << " -- Input file path: " << filePath << std::endl
    << " -- Output file path: " << outFilePath << std::endl
    << std::endl
    << " -- Processing ";
    if (pNevts == 0) std::cout << "all events." << std::endl;
    else std::cout << pNevts << " events per run." << std::endl;

    TRandom3 lRndm(0);
    std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

    /**********************************
    ** Input
    **********************************/
    // Get Reco files
    std::ostringstream inputrec;
    if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
    else
    inputrec << digifilePath << "/" << recoFileName;

    TChain *lRecTree = 0;

    lRecTree = new TChain("hgcalTupleTree/tree");

    if (nRuns == 0){
        lRecTree->AddFile(inputrec.str().c_str());
    }
    else {
        for (unsigned i(1);i<=nRuns;++i){
            std::ostringstream lstrrec;
            lstrrec << inputrec.str() << "_" << i << ".root";
            lRecTree->AddFile(lstrrec.str().c_str());
        }
    }

    if (!lRecTree){
        std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
        return 1;
    }

    /**********************************
    ** Output
    **********************************/
    // Define output file and the histograms contained
    TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");

    if (!outputFile) {
        std::cout << " -- Error, output file " << outFilePath
        << " cannot be opened. Please create output directory. Exiting..." << std::endl;
        return 1;
    }
    else {
        std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
    }
    outputFile->cd();

    TH1F* h_rechitsum = new TH1F("h_rechitsum","Rechitsum silicon;E[GeV]",100,0,20.);
    TH1F* h_rechitsumdead_Si = new TH1F("h_rechitsumdead_Si","Rechitsum dead silicon;E[GeV]",100,0,20.);
    TH1F* h_rechitsumave = new TH1F("h_rechitsumave","Sum energy average method;E[GeV]",100,0,20.);

    /**********************************
    ** ML Study output section
    **     - MLlayer is dead cell layer
    **     - MLeta is gen eta
    **     - MLphi is gen phi
    **     - MLni is ith dead cell neighbor
    **     - MLuni is 6 neighbors at layer+1
    **     - MLdni is 6 neighbors at layer-1
    **     - MLdead is dead cell rechit
    ** We also need to create a ?set? container to store the values before writing to TTree
    **********************************/
    TFile* fout = new TFile(MLFilePath.c_str(),"RECREATE");
    float MLlayer,MLcellid, MLeta, MLphi, MLdead, MLnup, MLndown, MLevent;
    float MLn1, MLn2, MLn3, MLn4, MLn5, MLn6;
    float MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6;
    float MLun1, MLun2, MLun3, MLun4, MLun5, MLun6;
    float MLrechitsum;
    TTree* t1 = new TTree("t1","sample");
    t1->Branch("MLlayer",&MLlayer,"MLlayer/F");
    t1->Branch("MLcellid",&MLcellid,"MLcellid/F");
    t1->Branch("MLeta",&MLeta,"MLeta/F");
    t1->Branch("MLphi",&MLphi,"MLphi/F");
    t1->Branch("MLn1",&MLn1,"MLn1/F");
    t1->Branch("MLn2",&MLn2,"MLn2/F");
    t1->Branch("MLn3",&MLn3,"MLn3/F");
    t1->Branch("MLn4",&MLn4,"MLn4/F");
    t1->Branch("MLn5",&MLn5,"MLn5/F");
    t1->Branch("MLn6",&MLn6,"MLn6/F");
    t1->Branch("MLdead",&MLdead,"MLdead/F");
    t1->Branch("MLnup",&MLnup,"MLnup/F");
    t1->Branch("MLndown",&MLndown,"MLndown/F");
    t1->Branch("MLun1",&MLun1,"MLun1/F");
    t1->Branch("MLun2",&MLun2,"MLun2/F");
    t1->Branch("MLun3",&MLun3,"MLun3/F");
    t1->Branch("MLun4",&MLun4,"MLun4/F");
    t1->Branch("MLun5",&MLun5,"MLun5/F");
    t1->Branch("MLun6",&MLun6,"MLun6/F");
    t1->Branch("MLdn1",&MLdn1,"MLdn1/F");
    t1->Branch("MLdn2",&MLdn2,"MLdn2/F");
    t1->Branch("MLdn3",&MLdn3,"MLdn3/F");
    t1->Branch("MLdn4",&MLdn4,"MLdn4/F");
    t1->Branch("MLdn5",&MLdn5,"MLdn5/F");
    t1->Branch("MLdn6",&MLdn6,"MLdn6/F");
    t1->Branch("MLevent",&MLevent,"MLevent/F");
    t1->Branch("MLrechitsum",&MLrechitsum,"MLrechitsum/F");

    /*
    ** Define a vector of the array:
    ** {dead cell(dc) layer, dc id, dc eta, dc phi, MLn1, MLn2, MLn3, MLn4, MLn5, MLn6, dc rechit}
    */
    std::vector<std::array<float, 26>> MLvectorev;

    /**********************************
    ** for missing channel study
    **********************************/
    // SILICON
    std::set<std::tuple<unsigned, int, int, unsigned, unsigned>> deadlistsi;

    // Kill cells and calculate statistics on adjacent dead cells
    unsigned N_try_success = 0; // Number of killed cells
    unsigned N_try_all = 0; // Number of trials to kill cells
    /*
    float N_cluster2 = 0; // Number of dead cells clusters (n_dead = 2)
    float N_clusters = 0; // Number of dead cells clusters (n_dead > 2)
    */

    TRandom3 r(0);
    for(unsigned lr = 1; lr <= 28; ++lr) {
        for(int waferU = -12; waferU <= 12; ++waferU) {
            for(int waferV = -12; waferV <= 12; ++waferV) {
                for(unsigned cellU = 0; cellU <= 16; ++cellU) {
                    for(unsigned cellV = 0; cellV <=16; ++cellV){
                        //if((cellU > cellV+9) || (cellV > cellU+8)) {
                            N_try_all++;
                            if(r.Rndm() < deadfrac){
                                N_try_success++;
                                deadlistsi.insert(std::make_tuple(
                                    lr,
                                    waferU,
                                    waferV,
                                    cellU,
                                    cellV
                                ));
                                std::array<float, 26> temp_vector;
                                for(unsigned k(0); k < 26; ++k) temp_vector[k] = 0;
                                temp_vector[0] = (float)(*itr).first; //layer
                                temp_vector[1] = (float)(*itr).second; //dead cell's id
                                MLvectorev.push_back(temp_vector);
                            }
                        //}
                    }
                }
            }
        }
    }

    std::cout << "List of dead Si cells was created successfully. \n"
    << "Killed " << N_try_success << " cells using " << N_try_all << " trials.\n"
    << std::endl;

    /* Old code
    for(unsigned i(0);i<nsidead;i++) {
        N_try_success++;
        ld_si=lRndm.Integer(nsilayer);
        range_si=simaxid[ld_si]-siminid[ld_si];
        cd_si=siminid[ld_si]+(lRndm.Integer(range_si));
        // Enforce that any dead cell has no more than one adjacent dead cell
        unsigned adj_ok = 0;
        //Need to switch this off for the moment
        if (deadlistsi.find(std::make_pair(ld_si, cd_si)) != deadlistsi.end()) {
            --i;
            continue;
        } //Insert end of comment here
        if (adjacent) {
            if (deadlistsi.find(std::make_pair(ld_si, cd_si-497)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si-496)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si-1)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si+1)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si+496)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si+497)) != deadlistsi.end()) adj_ok++;
            if (adj_ok > 1) {
                --i;
                N_cluster2++;
                //std::cout << "N_try_success = " << N_try_success
                << ": Found two or more adjascent dead cells for " << cd_si
                << " at layer " << ld_si << ". Aborting cell killing..."
                << std::endl;//Insert end of comment here
                continue;
            }
            if (adj_ok < 2) {
                deadlistsi.insert(std::make_pair(ld_si,cd_si));
            }
            if (adj_ok == 1) N_clusters++;
        }
        else {
            deadlistsi.insert(std::make_pair(ld_si,cd_si));
        }
    }

    // Print statistics on adjacent dead cells
    if (adjacent) {
        std::cout << std::string(120,'-') << std::endl
        << std::string(49,'-') << " Dead cells statistics "
        << std::string(48,'-') << std::endl
        << std::string(120,'-') << std::endl
        << "Number of dead cells clusters: " << N_clusters << std::endl
        << "Fraction of dead cluster cells: "
        << N_clusters*2./13983553. << std::endl
        << "Fraction of dead cells having a dead neighbor: "
        << N_clusters*2./nsidead << std::endl
        << "Dead fraction: " << deadfrac << std::endl
        << "Times the code tried to create clusters with more than 2 dead cells: "
        << N_cluster2 << std::endl
        << std::string(120,'-') << std::endl;
    }
    */

    // Define average energy in layers plus and minus 1
    std::set<std::tuple<unsigned, int, int, unsigned, unsigned>> adj_to_dead;
    // Define average energy in layer in cells plus and minus 1
    //std::vector<std::pair<unsigned, unsigned>> adj_to_dead_inlay;
    for(auto itr=deadlistsi.begin();itr!=deadlistsi.end();itr++ ) {
        adj_to_dead.insert({
            std::get<0>(*itr)-1,
            std::get<1>(*itr),
            std::get<2>(*itr),
            std::get<3>(*itr),
            std::get<4>(*itr)
        });
        adj_to_dead.insert({
            std::get<0>(*itr)+1,
            std::get<1>(*itr),
            std::get<2>(*itr),
            std::get<3>(*itr),
            std::get<4>(*itr)
        });

        /*
        adj_to_dead_inlay.push_back({(*itr).first, (*itr).second-497});
        adj_to_dead_inlay.push_back({(*itr).first, (*itr).second-496});
        adj_to_dead_inlay.push_back({(*itr).first, (*itr).second-1});
        adj_to_dead_inlay.push_back({(*itr).first, (*itr).second+1});
        adj_to_dead_inlay.push_back({(*itr).first, (*itr).second+496});
        adj_to_dead_inlay.push_back({(*itr).first, (*itr).second+497});
        */
    }

    /**********************************
    **  start event loop
    **********************************/
    const unsigned nEvts = (
        (pNevts > lRecTree->GetEntries() || pNevts==0)
        ? static_cast<unsigned>(lRecTree->GetEntries())
        : pNevts
    );
    std::cout << " -- Processing " << nEvts << " events out of "
    << lRecTree->GetEntries() << std::endl;

    //loop on events
    std::vector<float   > *rechitEnergy = 0;
    std::vector<float   > *rechitEta    = 0;
    std::vector<float   > *rechitPhi    = 0;
    std::vector<float   > *rechitPosx   = 0;
    std::vector<float   > *rechitPosy   = 0;
    std::vector<float   > *rechitPosz   = 0;
    std::vector<unsigned> *rechitLayer  = 0;
    std::vector<unsigned> *rechitIndex  = 0;
    std::vector<int     > *rechitWaferU = 0;
    std::vector<int     > *rechitWaferV = 0;
    std::vector<unsigned> *rechitCellU  = 0;
    std::vector<unsigned> *rechitCellV  = 0;
    std::vector<float   > *genEta       = 0;
    std::vector<float   > *genPhi       = 0;

    lRecTree->SetBranchAddress("HGCRecHitEnergy" ,&rechitEnergy);
    lRecTree->SetBranchAddress("HGCRecHitEta"    ,&rechitEta);
    lRecTree->SetBranchAddress("HGCRecHitPhi"    ,&rechitPhi);
    lRecTree->SetBranchAddress("HGCRecHitPosx"   ,&rechitPosx);
    lRecTree->SetBranchAddress("HGCRecHitPosy"   ,&rechitPosy);
    lRecTree->SetBranchAddress("HGCRecHitPosz"   ,&rechitPosz);
    lRecTree->SetBranchAddress("HGCRecHitLayer"  ,&rechitLayer);
    lRecTree->SetBranchAddress("HGCRecHitIndex"  ,&rechitIndex);
    lRecTree->SetBranchAddress("HGCRecHitWaferU" ,&rechitWaferU);
    lRecTree->SetBranchAddress("HGCRecHitWaferV" ,&rechitWaferV);
    lRecTree->SetBranchAddress("HGCRecHitCellU"  ,&rechitCellU);
    lRecTree->SetBranchAddress("HGCRecHitCellV"  ,&rechitCellV);
    lRecTree->SetBranchAddress("GenParEta"       ,&genEta);
    lRecTree->SetBranchAddress("GenParPhi"       ,&genPhi);

    unsigned ievtRec = 0;

    // Loop over entries (events)
    for (unsigned ievt(0); ievt<nEvts; ++ievt){
        if (ievtRec>=lRecTree->GetEntries()) continue;
        Long64_t local_entry = lRecTree->LoadTree(ievt);

        if (debug) std::cout << std::endl<<std::endl<<"... Processing entry: " << ievt << std::endl;
        else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

        if (local_entry < 0) continue;
        if (local_entry == 0) {
            lRecTree->SetBranchAddress("HGCRecHitEnergy" ,&rechitEnergy);
            lRecTree->SetBranchAddress("HGCRecHitEta"    ,&rechitEta);
            lRecTree->SetBranchAddress("HGCRecHitPhi"    ,&rechitPhi);
            lRecTree->SetBranchAddress("HGCRecHitPosx"   ,&rechitPosx);
            lRecTree->SetBranchAddress("HGCRecHitPosy"   ,&rechitPosy);
            lRecTree->SetBranchAddress("HGCRecHitPosz"   ,&rechitPosz);
            lRecTree->SetBranchAddress("HGCRecHitLayer"  ,&rechitLayer);
            lRecTree->SetBranchAddress("HGCRecHitIndex"  ,&rechitIndex);
            lRecTree->SetBranchAddress("HGCRecHitWaferU" ,&rechitWaferU);
            lRecTree->SetBranchAddress("HGCRecHitWaferV" ,&rechitWaferV);
            lRecTree->SetBranchAddress("HGCRecHitCellU"  ,&rechitCellU);
            lRecTree->SetBranchAddress("HGCRecHitCellV"  ,&rechitCellV);
            lRecTree->SetBranchAddress("GenParEta"       ,&genEta);
            lRecTree->SetBranchAddress("GenParPhi"       ,&genPhi);
        }

        lRecTree->GetEntry(ievtRec);

        double etagen   = 99999.;
        double phigen   = 99999.;
        //double thetagen = -1.;
        if((*genEta).size()>0) {
            etagen   = (*genEta)[0];
            phigen   = (*genPhi)[0];
            //thetagen = 2*TMath::ATan(exp(-etagen));
        }

        if (debug) std::cout << " - Event contains " << (*rechitEnergy).size()
        << " rechits." << std::endl;
        double coneSize = 0.3;
        double rechitsum = 0;
        double rechitsumdead_Si = 0;
        double rechitsumlaypn = 0;

        // Loop over hits of event
        for (unsigned iH(0); iH<(*rechitEnergy).size(); ++iH){
            unsigned layer   = (*rechitLayer)[iH];
            //double   xh      = (*rechitPosx)[iH];
            //double   yh      = (*rechitPosy)[iH];
            double   zh      = (*rechitPosz)[iH];
            double   lenergy = (*rechitEnergy)[iH];
            double   leta    = (*rechitEta)[iH];
            double   lphi    = (*rechitPhi)[iH];
            double   dR      = DeltaR(etagen,phigen,leta,lphi);
            //double   rgen    = zh*tan(thetagen);
            //double   xgen    = rgen*cos(phigen);
            //double   ygen    = rgen*sin(phigen);
            //double   dR1     = fabs(sqrt((xgen-xh)*(xgen-xh)+(ygen-yh)*(ygen-yh)));

            int waferU  = (*rechitWaferU)[iH];
            int waferV  = (*rechitWaferV)[iH];
            unsigned cellU   = (*rechitCellU)[iH];
            unsigned cellV   = (*rechitCellV)[iH];
            unsigned index   = (*rechitIndex)[iH];

            /* Select hits that are:
            **     - in CE-E
            **     - within DeltaR < 0.3 wrt gen particle
            **     - in positive endcap
            */
            if(!index && zh > 0 && dR < coneSize) {
                rechitsum += lenergy;
                std::tuple<unsigned, int, int, unsigned, unsigned> tempsi(layer,waferU,waferV,cellU,cellV);
                std::set<std::tuple<unsigned, int, int, unsigned, unsigned>>::iterator ibc=deadlistsi.find(tempsi);

                // Calculate energy without dead Si cells
                if(ibc == deadlistsi.end()) {
                    rechitsumdead_Si += lenergy;
                }

                for(auto itr=deadlistsi.begin();itr!=deadlistsi.end();itr++ ) {
                    // Perform Simple average method
                    bool simpleAverage1 = (
                        layer  == std::get<0>(*itr)+1 &&
                        waferU == std::get<1>(*itr)   &&
                        waferV == std::get<2>(*itr)   &&
                        cellU  == std::get<3>(*itr)   &&
                        cellV  == std::get<4>(*itr)
                    );
                    bool simpleAverage2 = (
                        layer  == std::get<0>(*itr)-1 &&
                        waferU == std::get<1>(*itr)   &&
                        waferV == std::get<2>(*itr)   &&
                        cellU  == std::get<3>(*itr)   &&
                        cellV  == std::get<4>(*itr)
                    );
                    if(simpleAverage1 || simpleAverage2){
                        rechitsumlaypn += lenergy/2;
                    }
                }
            }
        }
        // Fill histograms
        double rechitsumave=rechitsumlaypn+rechitsumdead_Si;
        h_rechitsumave->Fill(rechitsumave);
        h_rechitsum->Fill(rechitsum);
        h_rechitsumdead_Si->Fill(rechitsumdead_Si);

        ievtRec++;
    }

    if(debug) std::cout << "Writing files ..." << std::endl;
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
