/****************************************************
**  Channelog:
**      Moved to CMSSW output.
**      Editor: Christos Papageorgakis
**
**      This code uses ntuples to produce rechitsums
**      for a given dead Si cell fraction.
**      The output also contains an ntuple of all the
**      dead cells that can be used to train a DNN to
**      estimate the lost energy.
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

/* Function that gives a vector of tuples that describe the neighbors
** of the tuple that is given as the input.
*/
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
    std::get<4>(n5) += 1;
    // cell (+1,+1) wrt given
    std::tuple<int, int, int, int, int> n6(deadCell);
    std::get<3>(n6) += 1;

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
    //File to export data for ML **********Attention: Obsolete!
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
    float MLlayer, MLeta, MLphi, MLdead, MLevent;
    float MLwaferU, MLwaferV, MLcellU, MLcellV;
    float MLrechitsum;
    TTree* t1 = new TTree("t1","sample");
    t1->Branch("MLlayer"    ,&MLlayer    ,"MLlayer/F"    );
    t1->Branch("MLwaferU"   ,&MLwaferU   ,"MLwaferU/F"   );
    t1->Branch("MLwaferV"   ,&MLwaferV   ,"MLwaferV/F"   );
    t1->Branch("MLcellU"    ,&MLcellU    ,"MLcellU/F"    );
    t1->Branch("MLcellV"    ,&MLcellV    ,"MLcellV/F"    );
    t1->Branch("MLeta"      ,&MLeta      ,"MLeta/F"      );
    t1->Branch("MLphi"      ,&MLphi      ,"MLphi/F"      );
    t1->Branch("MLdead"     ,&MLdead     ,"MLdead/F"     );
    t1->Branch("MLevent"    ,&MLevent    ,"MLevent/F"    );
    t1->Branch("MLrechitsum",&MLrechitsum,"MLrechitsum/F");

    /*
    ** Define a vector of the array:
    ** {dead cell:
    **      layer, waferU, waferV, cellU, cellV,
    **      eta, phi,
    **      rechit,
    **      MLevent,
    **      MLrechitsum
    ** }
    */
    std::vector<std::array<float, 10>> MLvectorev;

    /**********************************
    ** for missing channel study
    **********************************/
    // SILICON
    std::set<std::tuple<int, int, int, int, int>> deadlistsi;

    // Kill cells and calculate statistics on adjacent dead cells
    unsigned N_try_success = 0; // Number of killed cells
    unsigned N_try_all = 0; // Number of trials to kill cells

    /* Loops over all possible cells and kills them with a probability
    ** given by the dead fraction.
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

                            std::array<float, 10> temp_vector;
                            for(unsigned k(0); k < 10; ++k) temp_vector[k] = 0;
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
            for(unsigned k(5); k < 10; ++k) (*itr)[k] = 0;
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
        if((*genEta).size()>0) {
            etagen   = (*genEta)[0];
            phigen   = (*genPhi)[0];
        }

        if (debug) std::cout << " - Event contains " << (*rechitEnergy).size()
        << " rechits." << std::endl;
        double coneSize = 0.3;

        // Loop over hits of event
        for (unsigned iH(0); iH<(*rechitEnergy).size(); ++iH){
            int layer   = (*rechitLayer)[iH];
            double   zh      = (*rechitPosz)[iH];
            double   lenergy = (*rechitEnergy)[iH];
            double   leta    = (*rechitEta)[iH];
            double   lphi    = (*rechitPhi)[iH];
            double   dR      = DeltaR(etagen,phigen,leta,lphi);

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
                std::tuple<int, int, int, int, int> tempsi(layer,waferU,waferV,cellU,cellV);
                std::set<std::tuple<int, int, int, int, int>>::iterator ibc=deadlistsi.find(tempsi);

                // Calculate energy without dead Si cells
                if(ibc == deadlistsi.end()) {
                    for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); itr++) {
                        (*itr)[9] += lenergy;
                    }
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
                            (*itr)[7] = lenergy;
                            (*itr)[8] = (float)ievt;
                        }
                    }
                }
            }
        }

        //Export the ML dataset values to the TTree
        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); ++itr) {
            bool check = 1;
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
                MLdead   = (*itr)[7];
                MLevent  = (*itr)[8];
                MLrechitsum = (*itr)[9];
                t1->Fill();
                if (check) {
                    check = 0;
                    std::cout << MLrechitsum << ", " << MLdead << std::endl;
                }
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
