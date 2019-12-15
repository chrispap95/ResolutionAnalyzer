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

std::vector<std::tuple<int, int, int, int, int>> getNeighbors(
    std::tuple<int, int, int, int, int> deadCell)
{
    std::vector<std::tuple<int, int, int, int, int>> neighbors;
    // Find same-layer neighboring cells
    // cell ( 0,-1) wrt given
    std::tuple<int, int, int, int, int> n1(deadCell);
    std::get<4>(n1) -= 1;
    // cell (-1,-1) wrt given
    std::tuple<int, int, int, int, int> n2(deadCell);
    std::get<3>(n2) -= 1;
    std::get<4>(n2) -= 1;
    // cell (-1, 0) wrt given
    std::tuple<int, int, int, int, int> n3(deadCell);
    std::get<3>(n3) -= 1;
    // cell ( 0,+1) wrt given
    std::tuple<int, int, int, int, int> n4(deadCell);
    std::get<4>(n4) += 1;
    // cell (+1, 0) wrt given
    std::tuple<int, int, int, int, int> n5(deadCell);
    std::get<3>(n5) += 1;
    // cell (+1,+1) wrt given
    std::tuple<int, int, int, int, int> n6(deadCell);
    std::get<3>(n6) += 1;
    std::get<4>(n6) += 1;

    // Check boundary conditions and make transitions between wafers when on the edge
    // For n1
    if (std::get<3>(n1) > -1 && std::get<3>(n1) < 8 && std::get<4>(n1) == -1){
        std::get<1>(n1) += 1;
        std::get<3>(n1) += 8;
        std::get<4>(n1) = 15;
    }else if (std::get<3>(n1)-std::get<4>(n1) == 9){
        std::get<1>(n1) += 1;
        std::get<2>(n1) += 1;
        std::get<3>(n1) -= 8;
        std::get<4>(n1) += 8;
    }
    // For n2
    if (std::get<3>(n2) > -1 && std::get<3>(n2) < 8 && std::get<4>(n2) == -1){
        std::get<1>(n2) += 1;
        std::get<3>(n2) += 8;
        std::get<4>(n2) = 15;
    }else if (std::get<4>(n2) > -1 && std::get<4>(n2) < 8 && std::get<3>(n2) == -1){
        std::get<2>(n2) -= 1;
        std::get<3>(n2) = 15;
        std::get<4>(n2) += 8;
    }
    // For n3
    if (std::get<4>(n3) > -1 && std::get<4>(n3) < 8 && std::get<3>(n3) == -1){
        std::get<2>(n3) -= 1;
        std::get<3>(n3) = 15;
        std::get<4>(n3) += 8;
    }else if (std::get<4>(n3)-std::get<3>(n3) == 8){
        std::get<1>(n3) -= 1;
        std::get<2>(n3) -= 1;
        std::get<3>(n3) += 8;
        std::get<4>(n3) -= 8;
    }
    // For n4
    if (std::get<3>(n4) > 7 && std::get<3>(n4) < 16 && std::get<4>(n4) == 16){
        std::get<1>(n4) -= 1;
        std::get<3>(n4) -= 8;
        std::get<4>(n4) = 0;
    }else if (std::get<4>(n4)-std::get<3>(n4) == 8){
        std::get<1>(n4) -= 1;
        std::get<2>(n4) -= 1;
        std::get<3>(n4) += 8;
        std::get<4>(n4) -= 8;
    }
    // For n5
    if (std::get<4>(n5) > 7 && std::get<4>(n5) < 16 && std::get<3>(n5) == 16){
        std::get<2>(n5) += 1;
        std::get<3>(n5) = 0;
        std::get<4>(n5) -= 8;
    }else if (std::get<3>(n5) > 7 && std::get<3>(n5) < 16 && std::get<4>(n5) == 16){
        std::get<1>(n5) -= 1;
        std::get<3>(n5) -= 8;
        std::get<4>(n5) = 0;
    }
    // For n6
    if (std::get<4>(n6) > 7 && std::get<4>(n6) < 16 && std::get<3>(n6) == 16){
        std::get<2>(n6) += 1;
        std::get<3>(n6) = 0;
        std::get<4>(n6) -= 8;
    }else if (std::get<3>(n6)-std::get<4>(n6) == 9){
        std::get<1>(n6) += 1;
        std::get<2>(n6) += 1;
        std::get<3>(n6) -= 8;
        std::get<4>(n6) += 8;
    }

    neighbors.push_back(n1);
    neighbors.push_back(n2);
    neighbors.push_back(n3);
    neighbors.push_back(n4);
    neighbors.push_back(n5);
    neighbors.push_back(n6);
    return neighbors;
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
    float MLlayer, MLeta, MLphi, MLdead, MLnup, MLndown, MLevent;
    float MLwaferU, MLwaferV, MLcellU, MLcellV;
    float MLn1, MLn2, MLn3, MLn4, MLn5, MLn6;
    float MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6;
    float MLun1, MLun2, MLun3, MLun4, MLun5, MLun6;
    float MLrechitsum;
    TTree* t1 = new TTree("t1","sample");
    t1->Branch("MLlayer"    ,&MLlayer    ,"MLlayer/F"    );
    t1->Branch("MLwaferU"   ,&MLwaferU   ,"MLwaferU/F"   );
    t1->Branch("MLwaferV"   ,&MLwaferV   ,"MLwaferV/F"   );
    t1->Branch("MLcellU"    ,&MLcellU    ,"MLcellU/F"    );
    t1->Branch("MLcellV"    ,&MLcellV    ,"MLcellV/F"    );
    t1->Branch("MLeta"      ,&MLeta      ,"MLeta/F"      );
    t1->Branch("MLphi"      ,&MLphi      ,"MLphi/F"      );
    t1->Branch("MLn1"       ,&MLn1       ,"MLn1/F"       );
    t1->Branch("MLn2"       ,&MLn2       ,"MLn2/F"       );
    t1->Branch("MLn3"       ,&MLn3       ,"MLn3/F"       );
    t1->Branch("MLn4"       ,&MLn4       ,"MLn4/F"       );
    t1->Branch("MLn5"       ,&MLn5       ,"MLn5/F"       );
    t1->Branch("MLn6"       ,&MLn6       ,"MLn6/F"       );
    t1->Branch("MLdead"     ,&MLdead     ,"MLdead/F"     );
    t1->Branch("MLnup"      ,&MLnup      ,"MLnup/F"      );
    t1->Branch("MLndown"    ,&MLndown    ,"MLndown/F"    );
    t1->Branch("MLun1"      ,&MLun1      ,"MLun1/F"      );
    t1->Branch("MLun2"      ,&MLun2      ,"MLun2/F"      );
    t1->Branch("MLun3"      ,&MLun3      ,"MLun3/F"      );
    t1->Branch("MLun4"      ,&MLun4      ,"MLun4/F"      );
    t1->Branch("MLun5"      ,&MLun5      ,"MLun5/F"      );
    t1->Branch("MLun6"      ,&MLun6      ,"MLun6/F"      );
    t1->Branch("MLdn1"      ,&MLdn1      ,"MLdn1/F"      );
    t1->Branch("MLdn2"      ,&MLdn2      ,"MLdn2/F"      );
    t1->Branch("MLdn3"      ,&MLdn3      ,"MLdn3/F"      );
    t1->Branch("MLdn4"      ,&MLdn4      ,"MLdn4/F"      );
    t1->Branch("MLdn5"      ,&MLdn5      ,"MLdn5/F"      );
    t1->Branch("MLdn6"      ,&MLdn6      ,"MLdn6/F"      );
    t1->Branch("MLevent"    ,&MLevent    ,"MLevent/F"    );
    t1->Branch("MLrechitsum",&MLrechitsum,"MLrechitsum/F");

    /*
    ** Define a vector of the array:
    ** {dead cell:
    **      layer, waferU, waferV, cellU, cellV,
    **      eta, phi,
    **      MLn1, MLn2, MLn3, MLn4, MLn5, MLn6,
    **      rechit,
    **      MLup, MLdown,
    **      MLun1, MLun2, MLun3, MLun4, MLun5, MLun6,
    **      MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6,
    **      MLevent,
    **      MLrechitsum
    ** }
    */
    std::vector<std::array<float, 30>> MLvectorev;

    /**********************************
    ** for missing channel study
    **********************************/
    // SILICON
    std::set<std::tuple<int, int, int, int, int>> deadlistsi;

    // Define average energy in layers plus and minus 1
    std::set<std::tuple<int, int, int, int, int, int>> adj_to_dead;
    std::set<std::tuple<int, int, int, int, int, int>> adj_to_dead_inlay;

    // Kill cells and calculate statistics on adjacent dead cells
    unsigned N_try_success = 0; // Number of killed cells
    unsigned N_try_all = 0; // Number of trials to kill cells
    /*
    float N_cluster2 = 0; // Number of dead cells clusters (n_dead = 2)
    float N_clusters = 0; // Number of dead cells clusters (n_dead > 2)
    */

    TRandom3 r(0);
    for(int lr = 1; lr <= 28; ++lr) {
        for(int waferU = -12; waferU <= 12; ++waferU) {
            for(int waferV = -12; waferV <= 12; ++waferV) {
                for(int cellU = 0; cellU <= 16; ++cellU) {
                    for(int cellV = 0; cellV <=16; ++cellV){
                        N_try_all++;
                        if(r.Rndm() < deadfrac){
                            N_try_success++;
                            std::tuple<int,int,int,int,int> deadCell(
                                lr,
                                waferU,
                                waferV,
                                cellU,
                                cellV
                            );
                            deadlistsi.insert(deadCell);

                            adj_to_dead.insert({
                                0, //corresponds to cell bellow
                                std::get<0>(deadCell)-1,
                                std::get<1>(deadCell),
                                std::get<2>(deadCell),
                                std::get<3>(deadCell),
                                std::get<4>(deadCell)
                            });
                            adj_to_dead.insert({
                                1, //corresponds to cell above
                                std::get<0>(deadCell)+1,
                                std::get<1>(deadCell),
                                std::get<2>(deadCell),
                                std::get<3>(deadCell),
                                std::get<4>(deadCell)
                            });

                            std::vector<std::tuple<int,int,int,int,int>> inLayerNeighbors;
                            inLayerNeighbors = getNeighbors(deadCell);
                            int iN = 0;
                            for(auto itr = inLayerNeighbors.begin(); itr!=inLayerNeighbors.end(); ++itr){
                                adj_to_dead_inlay.insert({
                                    iN,
                                    std::get<0>(*itr),
                                    std::get<1>(*itr),
                                    std::get<2>(*itr),
                                    std::get<3>(*itr),
                                    std::get<4>(*itr)
                                });
                            }

                            std::array<float, 30> temp_vector;
                            for(unsigned k(0); k < 30; ++k) temp_vector[k] = 0;
                            temp_vector[0] = (float)lr; //layer
                            temp_vector[1] = (float)waferU; //dead cell's waferU
                            temp_vector[2] = (float)waferV; //dead cell's waferV
                            temp_vector[3] = (float)cellU;  //dead cell's cellU
                            temp_vector[4] = (float)cellV;  //dead cell's cellV
                            MLvectorev.push_back(temp_vector);
                        }
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
    std::vector<int     > *rechitLayer  = 0;
    std::vector<int     > *rechitIndex  = 0;
    std::vector<int     > *rechitWaferU = 0;
    std::vector<int     > *rechitWaferV = 0;
    std::vector<int     > *rechitCellU  = 0;
    std::vector<int     > *rechitCellV  = 0;
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
        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
            for(unsigned k(5); k < 30; ++k) (*itr)[k] = 0;
        }
        if (ievtRec>=lRecTree->GetEntries()) continue;
        Long64_t local_entry = lRecTree->LoadTree(ievt);

        if (debug) std::cout << std::endl<<std::endl << "... Processing entry: " << ievt << std::endl;
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
            int layer   = (*rechitLayer)[iH];
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
            int cellU   = (*rechitCellU)[iH];
            int cellV   = (*rechitCellV)[iH];
            int index   = (*rechitIndex)[iH];

            /* Select hits that are:
            **     - in CE-E
            **     - within DeltaR < 0.3 wrt gen particle
            **     - in positive endcap
            */
            if(!index && zh > 0 && dR < coneSize) {
                rechitsum += lenergy;
                std::tuple<int, int, int, int, int> tempsi(layer,waferU,waferV,cellU,cellV);
                std::set<std::tuple<int, int, int, int, int>>::iterator ibc=deadlistsi.find(tempsi);

                // Calculate energy without dead Si cells
                if(ibc == deadlistsi.end()) {
                    rechitsumdead_Si += lenergy;
                    MLrechitsum += lenergy;
                }else {
                    // Do stuff with dead cells
                    /* ML code
                    ** Input dead cells eta, phi and rechits
                    */
                    for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                        if( (*itr)[0] == layer &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[5] = leta;
                            (*itr)[6] = lphi;
                            (*itr)[13] = lenergy;
                            (*itr)[29] = (float)ievt;
                        }
                    }
                }

                /* Perform Simple Average
                ** First, check if the cell is in a neighbors list
                */
                std::tuple<int, int, int, int, int, int> tempsiU(
                    1,layer,waferU,waferV,cellU,cellV
                );
                std::tuple<int, int, int, int, int, int> tempsiD(
                    0,layer,waferU,waferV,cellU,cellV
                );
                std::set<std::tuple<int, int, int, int, int, int>>::iterator itrU=adj_to_dead.find(tempsiU);
                std::set<std::tuple<int, int, int, int, int, int>>::iterator itrD=adj_to_dead.find(tempsiD);

                if(itrU!=adj_to_dead.end()) {
                    rechitsumlaypn += lenergy/2;
                    for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                        if( (*itr)[0] == layer-1 &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[14] = lenergy;
                        }
                    }
                }
                if(itrD!=adj_to_dead.end()) {
                    rechitsumlaypn += lenergy/2;
                    for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                        if( (*itr)[0] == layer+1 &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[15] = lenergy;
                        }
                    }
                }

                // Get rechits of dead cells' neighbors
                for(int n = 0; n < 6; ++n){
                    // Same layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiNn(
                        n,layer,waferU,waferV,cellU,cellV
                    );
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrNn=adj_to_dead_inlay.find(tempsiNn);
                    if(itrNn!=adj_to_dead_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int>> sameLayerNeighbors;
                        sameLayerNeighbors = getNeighbors(tempsi);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrNn)+3)%6;
                        std::tuple<int, int, int, int> deadCell;
                        std::get<0>(deadCell) = std::get<1>(sameLayerNeighbors[nn]);
                        std::get<1>(deadCell) = std::get<2>(sameLayerNeighbors[nn]);
                        std::get<2>(deadCell) = std::get<3>(sameLayerNeighbors[nn]);
                        std::get<3>(deadCell) = std::get<4>(sameLayerNeighbors[nn]);
                        std::cout << "Debug: " << n << ", " << nn << ", " << lenergy << std::endl;
                        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                            if( (*itr)[0] == layer &&
                                (*itr)[1] == std::get<0>(deadCell) && (*itr)[2] == std::get<1>(deadCell) &&
                                (*itr)[3] == std::get<2>(deadCell) && (*itr)[4] == std::get<3>(deadCell)
                            ){
                                (*itr)[n+7] = lenergy;
                            }
                        }
                    }
                    // Next layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiUNn(
                        n,layer,waferU,waferV,cellU,cellV
                    );
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrUNn=adj_to_dead_inlay.find(tempsiUNn);
                    if(itrUNn!=adj_to_dead_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int>> nextLayerNeighbors;
                        nextLayerNeighbors = getNeighbors(tempsi);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrNn)+3)%6;
                        std::tuple<int, int, int, int> deadCell;
                        std::get<0>(deadCell) = std::get<1>(nextLayerNeighbors[nn]);
                        std::get<1>(deadCell) = std::get<2>(nextLayerNeighbors[nn]);
                        std::get<2>(deadCell) = std::get<3>(nextLayerNeighbors[nn]);
                        std::get<3>(deadCell) = std::get<4>(nextLayerNeighbors[nn]);
                        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                            if( (*itr)[0] == layer-1 &&
                            (*itr)[1] == std::get<0>(deadCell) && (*itr)[2] == std::get<1>(deadCell) &&
                            (*itr)[3] == std::get<2>(deadCell) && (*itr)[4] == std::get<3>(deadCell)
                            ){
                                (*itr)[n+16] = lenergy;
                            }
                        }
                    }
                    // Previous layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiDNn(
                        n,layer,waferU,waferV,cellU,cellV
                    );
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrDNn=adj_to_dead_inlay.find(tempsiDNn);
                    if(itrDNn!=adj_to_dead_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int>> prevLayerNeighbors;
                        prevLayerNeighbors = getNeighbors(tempsi);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrNn)+3)%6;
                        std::tuple<int, int, int, int> deadCell;
                        std::get<0>(deadCell) = std::get<1>(prevLayerNeighbors[nn]);
                        std::get<1>(deadCell) = std::get<2>(prevLayerNeighbors[nn]);
                        std::get<2>(deadCell) = std::get<3>(prevLayerNeighbors[nn]);
                        std::get<3>(deadCell) = std::get<4>(prevLayerNeighbors[nn]);
                        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                            if( (*itr)[0] == layer+1 &&
                            (*itr)[1] == std::get<0>(deadCell) && (*itr)[2] == std::get<1>(deadCell) &&
                            (*itr)[3] == std::get<2>(deadCell) && (*itr)[4] == std::get<3>(deadCell)
                            ){
                                (*itr)[n+22] = lenergy;
                            }
                        }
                    }
                }
            }
        }
        // Fill histograms
        double rechitsumave=rechitsumlaypn+rechitsumdead_Si;
        h_rechitsumave->Fill(rechitsumave);
        h_rechitsum->Fill(rechitsum);
        h_rechitsumdead_Si->Fill(rechitsumdead_Si);

        //Export the ML dataset values to the TTree
        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); ++itr) {
            if ((*itr)[5] > 0) {
                /* This condition is necessary to ensure the cell was within
                ** the cone.
                */
                MLlayer  = (*itr)[0];
                MLwaferU = (*itr)[1];
                MLwaferV = (*itr)[2];
                MLcellU  = (*itr)[3];
                MLcellV  = (*itr)[4];
                MLeta    = (*itr)[5];
                MLphi    = (*itr)[6];
                MLn1     = (*itr)[7];
                MLn2     = (*itr)[8];
                MLn3     = (*itr)[9];
                MLn4     = (*itr)[10];
                MLn5     = (*itr)[11];
                MLn6     = (*itr)[12];
                MLdead   = (*itr)[13];
                MLnup    = (*itr)[14];
                MLndown  = (*itr)[15];
                MLun1    = (*itr)[16];
                MLun2    = (*itr)[17];
                MLun3    = (*itr)[18];
                MLun4    = (*itr)[19];
                MLun5    = (*itr)[20];
                MLun6    = (*itr)[21];
                MLdn1    = (*itr)[22];
                MLdn2    = (*itr)[23];
                MLdn3    = (*itr)[24];
                MLdn4    = (*itr)[26];
                MLdn5    = (*itr)[27];
                MLdn6    = (*itr)[28];
                MLevent  = (*itr)[29];
                MLrechitsum = rechitsumdead_Si;
                t1->Fill();
            }
        }
        ievtRec++;
    }

    if(debug) std::cout << "Writing files ..." << std::endl;
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
