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
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
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

/* Returns a vector of tuples that describe the neighbors
** of the tuple that is given as the input.
**     - offset: is 0 for low density areas and 4 for high density areas
*/
std::vector<std::tuple<int, int, int, int, int>> getNeighborsSi(
    std::tuple<int, int, int, int, int> deadCell, bool isDense)
{
    int offset = 0;
    if (isDense) offset = 4;
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
    if (std::get<3>(n1) > -1 && std::get<3>(n1) < (8 + offset) && std::get<4>(n1) == -1){
        std::get<1>(n1) += 1;
        std::get<3>(n1) += 8 + offset;
        std::get<4>(n1) = 15 + 2*offset;
    }else if (std::get<3>(n1)-std::get<4>(n1) == 9 + offset){
        std::get<1>(n1) += 1;
        std::get<2>(n1) += 1;
        std::get<3>(n1) -= 8 + offset;
        std::get<4>(n1) += 8 + offset;
    }
    // For n2
    if (std::get<3>(n2) > -1 && std::get<3>(n2) < (8 + offset) && std::get<4>(n2) == -1){
        std::get<1>(n2) += 1;
        std::get<3>(n2) += 8 + offset;
        std::get<4>(n2) = 15 + 2*offset;
    }else if (std::get<4>(n2) > -1 && std::get<4>(n2) < (8 + offset) && std::get<3>(n2) == -1){
        std::get<2>(n2) -= 1;
        std::get<3>(n2) = 15 + 2*offset;
        std::get<4>(n2) += 8 + offset;
    }
    // For n3
    if (std::get<4>(n3) > -1 && std::get<4>(n3) < (8 + offset) && std::get<3>(n3) == -1){
        std::get<2>(n3) -= 1;
        std::get<3>(n3) = 15 + 2*offset;
        std::get<4>(n3) += 8 + offset;
    }else if (std::get<4>(n3)-std::get<3>(n3) == 8 + offset){
        std::get<1>(n3) -= 1;
        std::get<2>(n3) -= 1;
        std::get<3>(n3) += 8 + offset;
        std::get<4>(n3) -= 8 + offset;
    }
    // For n4
    if (std::get<3>(n4) > (7 + offset) && std::get<3>(n4) < (16 + 2*offset) && std::get<4>(n4) == (16 + 2*offset)){
        std::get<1>(n4) -= 1;
        std::get<3>(n4) -= 8 + offset;
        std::get<4>(n4) = 0;
    }else if (std::get<4>(n4)-std::get<3>(n4) == (8 + offset)){
        std::get<1>(n4) -= 1;
        std::get<2>(n4) -= 1;
        std::get<3>(n4) += 8 + offset;
        std::get<4>(n4) -= 8 + offset;
    }
    // For n5
    if (std::get<4>(n5) > (7 + offset) && std::get<4>(n5) < (16 + 2*offset) && std::get<3>(n5) == (16 + 2*offset)){
        std::get<2>(n5) += 1;
        std::get<3>(n5) = 0;
        std::get<4>(n5) -= 8 + offset;
    }else if (std::get<3>(n5) > (7 + offset) && std::get<3>(n5) < (16 + 2*offset) && std::get<4>(n5) == (16 + 2*offset)){
        std::get<1>(n5) -= 1;
        std::get<3>(n5) -= 8 + offset;
        std::get<4>(n5) = 0;
    }
    // For n6
    if (std::get<4>(n6) > (7 + offset) && std::get<4>(n6) < (16 + 2*offset) && std::get<3>(n6) == (16 + 2*offset)){
        std::get<2>(n6) += 1;
        std::get<3>(n6) = 0;
        std::get<4>(n6) -= 8 + offset;
    }else if (std::get<3>(n6)-std::get<4>(n6) == (9 + offset)){
        std::get<1>(n6) += 1;
        std::get<2>(n6) += 1;
        std::get<3>(n6) -= 8 + offset;
        std::get<4>(n6) += 8 + offset;
    }

    neighbors.push_back(n1);
    neighbors.push_back(n2);
    neighbors.push_back(n3);
    neighbors.push_back(n4);
    neighbors.push_back(n5);
    neighbors.push_back(n6);
    return neighbors;
}

std::vector<std::tuple<int, int, int>> getNeighborsScint(
    std::tuple<int, int, int> deadChannel)
{
  std::vector<std::tuple<int, int, int>> neighbors;
  // Find same-layer neighboring channel
  std::tuple<int, int, int> n1(deadChannel);
  std::get<2>(n1) -= 1;
  std::tuple<int, int, int> n2(deadChannel);
  std::get<1>(n2) -= 1;
  std::get<2>(n2) -= 1;
  std::tuple<int, int, int> n3(deadChannel);
  std::get<1>(n3) -= 1;
  std::tuple<int, int, int> n4(deadChannel);
  std::get<1>(n4) -= 1;
  std::get<2>(n4) += 1;
  std::tuple<int, int, int> n5(deadChannel);
  std::get<2>(n5) += 1;
  std::tuple<int, int, int> n6(deadChannel);
  std::get<1>(n6) += 1;
  std::get<2>(n6) += 1;
  std::tuple<int, int, int> n7(deadChannel);
  std::get<1>(n7) += 1;
  std::tuple<int, int, int> n8(deadChannel);
  std::get<1>(n8) += 1;
  std::get<2>(n8) -= 1;

  if (std::get<2>(n1) == 289) std::get<2>(n1) = 1;
  if (std::get<2>(n1) == 0) std::get<2>(n1) = 288;
  if (std::get<2>(n2) == 289) std::get<2>(n2) = 1;
  if (std::get<2>(n2) == 0) std::get<2>(n2) = 288;
  if (std::get<2>(n3) == 289) std::get<2>(n3) = 1;
  if (std::get<2>(n3) == 0) std::get<2>(n3) = 288;
  if (std::get<2>(n4) == 289) std::get<2>(n4) = 1;
  if (std::get<2>(n4) == 0) std::get<2>(n4) = 288;
  if (std::get<2>(n5) == 289) std::get<2>(n5) = 1;
  if (std::get<2>(n5) == 0) std::get<2>(n5) = 288;
  if (std::get<2>(n6) == 289) std::get<2>(n6) = 1;
  if (std::get<2>(n6) == 0) std::get<2>(n6) = 288;
  if (std::get<2>(n7) == 289) std::get<2>(n7) = 1;
  if (std::get<2>(n7) == 0) std::get<2>(n7) = 288;
  if (std::get<2>(n8) == 289) std::get<2>(n8) = 1;
  if (std::get<2>(n8) == 0) std::get<2>(n8) = 288;

  neighbors.push_back(n1);
  neighbors.push_back(n2);
  neighbors.push_back(n3);
  neighbors.push_back(n4);
  neighbors.push_back(n5);
  neighbors.push_back(n6);
  neighbors.push_back(n7);
  neighbors.push_back(n8);

  return neighbors;
}

int main(int argc, char** argv){
    //Input output and config options
    std::string cfg;
    unsigned pNevts;
    std::string outFilePath;
    std::string filePath;
    std::string digifilePath;
    unsigned nRuns,firstRun;
    std::string recoFileName;
    unsigned debug;
    double deadfrac;
    bool adjacent;
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
    ("firstRun",        po::value<unsigned>(&firstRun)->default_value(1))
    ("debug,d",         po::value<unsigned>(&debug)->default_value(0))
    ("deadfrac",        po::value<double>(&deadfrac)->default_value(0))
    //Restrict number of adjacent dead cells
    ("adjacent",        po::value<bool>(&adjacent)->default_value(0))
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

    lRecTree = new TChain("ana/hgc");

    if (nRuns == 0){
        lRecTree->AddFile(inputrec.str().c_str());
    }
    else {
        for (unsigned i = firstRun;i<nRuns+firstRun;++i){
            std::ostringstream lstrrec;
            lstrrec << inputrec.str() << "_" << i << ".root";
            lRecTree->AddFile(lstrrec.str().c_str());
            std::cout << "Opening file: " << lstrrec.str() << std::endl;
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
    float MLlayer, MLeta, MLphi, MLdead, MLnup, MLndown, MLevent;
    float MLwaferU, MLwaferV, MLcellU, MLcellV, MLieta, MLiphi;
    float MLn1, MLn2, MLn3, MLn4, MLn5, MLn6, MLn7, MLn8;
    float MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6, MLdn7, MLdn8;
    float MLun1, MLun2, MLun3, MLun4, MLun5, MLun6, MLun7, MLun8;
    float MLrechitsum, MLthickness;
    TTree* t1 = new TTree("t1","sample");
    t1->Branch("layer"    ,&MLlayer    ,"layer/F"    );
    t1->Branch("waferU"   ,&MLwaferU   ,"waferU/F"   );
    t1->Branch("waferV"   ,&MLwaferV   ,"waferV/F"   );
    t1->Branch("cellU"    ,&MLcellU    ,"cellU/F"    );
    t1->Branch("cellV"    ,&MLcellV    ,"cellV/F"    );
    t1->Branch("eta"      ,&MLeta      ,"eta/F"      );
    t1->Branch("phi"      ,&MLphi      ,"phi/F"      );
    t1->Branch("n1"       ,&MLn1       ,"n1/F"       );
    t1->Branch("n2"       ,&MLn2       ,"n2/F"       );
    t1->Branch("n3"       ,&MLn3       ,"n3/F"       );
    t1->Branch("n4"       ,&MLn4       ,"n4/F"       );
    t1->Branch("n5"       ,&MLn5       ,"n5/F"       );
    t1->Branch("n6"       ,&MLn6       ,"n6/F"       );
    t1->Branch("dead"     ,&MLdead     ,"dead/F"     );
    t1->Branch("nup"      ,&MLnup      ,"nup/F"      );
    t1->Branch("ndown"    ,&MLndown    ,"ndown/F"    );
    t1->Branch("un1"      ,&MLun1      ,"un1/F"      );
    t1->Branch("un2"      ,&MLun2      ,"un2/F"      );
    t1->Branch("un3"      ,&MLun3      ,"un3/F"      );
    t1->Branch("un4"      ,&MLun4      ,"un4/F"      );
    t1->Branch("un5"      ,&MLun5      ,"un5/F"      );
    t1->Branch("un6"      ,&MLun6      ,"un6/F"      );
    t1->Branch("dn1"      ,&MLdn1      ,"dn1/F"      );
    t1->Branch("dn2"      ,&MLdn2      ,"dn2/F"      );
    t1->Branch("dn3"      ,&MLdn3      ,"dn3/F"      );
    t1->Branch("dn4"      ,&MLdn4      ,"dn4/F"      );
    t1->Branch("dn5"      ,&MLdn5      ,"dn5/F"      );
    t1->Branch("dn6"      ,&MLdn6      ,"dn6/F"      );
    t1->Branch("event"    ,&MLevent    ,"event/F"    );
    t1->Branch("rechitsum",&MLrechitsum,"rechitsum/F");
    t1->Branch("thickness",&MLthickness,"thickness/F");

    TTree* t2 = new TTree("t2","sample");
    t2->Branch("layer"    ,&MLlayer    ,"layer/F"    );
    t2->Branch("ieta"     ,&MLieta     ,"ieta/F"     );
    t2->Branch("iphi"     ,&MLiphi     ,"iphi/F"     );
    t2->Branch("eta"      ,&MLeta      ,"eta/F"      );
    t2->Branch("phi"      ,&MLphi      ,"phi/F"      );
    t2->Branch("n1"       ,&MLn1       ,"n1/F"       );
    t2->Branch("n2"       ,&MLn2       ,"n2/F"       );
    t2->Branch("n3"       ,&MLn3       ,"n3/F"       );
    t2->Branch("n4"       ,&MLn4       ,"n4/F"       );
    t2->Branch("n5"       ,&MLn5       ,"n5/F"       );
    t2->Branch("n6"       ,&MLn6       ,"n6/F"       );
    t2->Branch("n7"       ,&MLn7       ,"n7/F"       );
    t2->Branch("n8"       ,&MLn8       ,"n8/F"       );
    t2->Branch("dead"     ,&MLdead     ,"dead/F"     );
    t2->Branch("nup"      ,&MLnup      ,"nup/F"      );
    t2->Branch("ndown"    ,&MLndown    ,"ndown/F"    );
    t2->Branch("un1"      ,&MLun1      ,"un1/F"      );
    t2->Branch("un2"      ,&MLun2      ,"un2/F"      );
    t2->Branch("un3"      ,&MLun3      ,"un3/F"      );
    t2->Branch("un4"      ,&MLun4      ,"un4/F"      );
    t2->Branch("un5"      ,&MLun5      ,"un5/F"      );
    t2->Branch("un6"      ,&MLun6      ,"un6/F"      );
    t2->Branch("un7"      ,&MLun7      ,"un7/F"      );
    t2->Branch("un8"      ,&MLun8      ,"un8/F"      );
    t2->Branch("dn1"      ,&MLdn1      ,"dn1/F"      );
    t2->Branch("dn2"      ,&MLdn2      ,"dn2/F"      );
    t2->Branch("dn3"      ,&MLdn3      ,"dn3/F"      );
    t2->Branch("dn4"      ,&MLdn4      ,"dn4/F"      );
    t2->Branch("dn5"      ,&MLdn5      ,"dn5/F"      );
    t2->Branch("dn6"      ,&MLdn6      ,"dn6/F"      );
    t2->Branch("dn7"      ,&MLdn7      ,"dn7/F"      );
    t2->Branch("dn8"      ,&MLdn8      ,"dn8/F"      );
    t2->Branch("event"    ,&MLevent    ,"event/F"    );
    t2->Branch("rechitsum",&MLrechitsum,"rechitsum/F");

    /*
    ** Define a vector of the array for Si:
    ** {dead cell:
    **      layer, waferU, waferV, cellU, cellV,
    **      eta, phi,
    **      MLn1, MLn2, MLn3, MLn4, MLn5, MLn6,
    **      rechit,
    **      MLup, MLdown,
    **      MLun1, MLun2, MLun3, MLun4, MLun5, MLun6,
    **      MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6,
    **      MLthickness
    ** }
    ** and for scintillator:
    ** {dead channel:
    **      layer, ieta, iphim,
    **      eta, phi
    **      MLn1, MLn2, MLn3, MLn4, MLn5, MLn6, MLn7, MLn8,
    **      rechit,
    **      MLup, MLdown,
    **      MLun1, MLun2, MLun3, MLun4, MLun5, MLun6, MLun7, MLun8,
    **      MLdn1, MLdn2, MLdn3, MLdn4, MLdn5, MLdn6, MLdn7, MLdn8
    ** }
    */
    std::vector<std::array<float, 29>> MLvectorevSi;
    std::vector<std::array<float, 32>> MLvectorevScint;

    /**********************************
    ** for missing channel study
    **********************************/
    std::set<std::tuple<int, int, int, int, int>> deadlistSi;
    std::set<std::tuple<int, int, int>> deadlistScint;

    // Define average energy in layers plus and minus 1
    std::set<std::tuple<int, int, int, int, int, int>> adj_to_dead_Si;
    std::set<std::tuple<int, int, int, int, int, int>> adj_to_dead_Si_inlay;
    std::set<std::tuple<int, int, int, int>> adj_to_dead_Scint;
    std::set<std::tuple<int, int, int, int>> adj_to_dead_Scint_inlay;

    // Kill cells and calculate statistics on adjacent dead cells
    unsigned N_try_success_Si = 0; // Number of killed cells
    unsigned N_try_all_Si = 0; // Number of trials to kill cells
    unsigned N_try_success_Scint = 0; // Number of killed cells
    unsigned N_try_all_Scint = 0; // Number of trials to kill cells
    /*
    float N_cluster2 = 0; // Number of dead cells clusters (n_dead = 2)
    float N_clusters = 0; // Number of dead cells clusters (n_dead > 2)
    */

    /* Loops over all possible cells and kills them with a probability
    ** given by the dead fraction.
    ** offset: if in dense area then it is 4
    */
    TRandom3 r(0);
    int offset = 0;
    float r_denseLimit = 10.5;
    for(int lr = 1; lr < 51; ++lr) {
        for(int waferU = -11; waferU < 12; ++waferU) {
            for(int waferV = -11; waferV < 12; ++waferV) {
                int waferX = -2*waferU+waferV;
                int waferY = 2*waferV;
                if (lr > 28) r_denseLimit = 16.3-0.2*lr;
                if(sqrt(pow(waferX, 2)+pow(waferY, 2)) < r_denseLimit) {
                    if((abs(waferU-waferV) == 1) || (abs(waferU-waferV) == 4) ||
                       (abs(waferU-waferV) == 5) || (lr > 40)
                    ){
                        offset = 0;
                    } else {
                        offset = 4;
                    }
                } else {
                    offset = 0;
                }
                for(int cellU = 0; cellU < 16+2*offset; ++cellU) {
                    for(int cellV = 0; cellV < 16+2*offset; ++cellV) {
                        N_try_all_Si++;
                        if(r.Rndm() < deadfrac){
                            N_try_success_Si++;
                            std::tuple<int,int,int,int,int> deadCell(
                                lr,
                                waferU,
                                waferV,
                                cellU,
                                cellV
                            );
                            deadlistSi.insert(deadCell);

                            adj_to_dead_Si.insert({
                                0, //corresponds to cell bellow
                                std::get<0>(deadCell)-1,
                                std::get<1>(deadCell),
                                std::get<2>(deadCell),
                                std::get<3>(deadCell),
                                std::get<4>(deadCell)
                            });
                            adj_to_dead_Si.insert({
                                1, //corresponds to cell above
                                std::get<0>(deadCell)+1,
                                std::get<1>(deadCell),
                                std::get<2>(deadCell),
                                std::get<3>(deadCell),
                                std::get<4>(deadCell)
                            });

                            std::vector<std::tuple<int,int,int,int,int>> inLayerNeighbors;
                            bool isDense = (offset == 4) ? 1 : 0;
                            inLayerNeighbors = getNeighborsSi(deadCell, isDense);
                            int iN = 0;
                            for(auto itr = inLayerNeighbors.begin(); itr!=inLayerNeighbors.end(); ++itr){
                                adj_to_dead_Si_inlay.insert({
                                    iN,
                                    std::get<0>(*itr),
                                    std::get<1>(*itr),
                                    std::get<2>(*itr),
                                    std::get<3>(*itr),
                                    std::get<4>(*itr)
                                });
                                iN++;
                            }

                            std::array<float, 29> temp_vector;
                            for(unsigned k(0); k < 29; ++k) temp_vector[k] = 0;
                            temp_vector[0] = (float)lr; //layer
                            temp_vector[1] = (float)waferU; //dead cell's waferU
                            temp_vector[2] = (float)waferV; //dead cell's waferV
                            temp_vector[3] = (float)cellU;  //dead cell's cellU
                            temp_vector[4] = (float)cellV;  //dead cell's cellV
                            MLvectorevSi.push_back(temp_vector);
                        }
                    }
                }
            }
        }
    }

    for(int lr = 37; lr < 51; ++lr) {
        for(int ie = 10; ie < 40; ++ie) {
            for(int ip = 1; ip < 289; ++ip) {
                N_try_all_Scint++;
                if(r.Rndm() < deadfrac) {
                    N_try_success_Scint++;
                    std::tuple<int,int,int> deadChannel(
                        lr,
                        ie,
                        ip
                    );
                    deadlistScint.insert(deadChannel);

                    adj_to_dead_Scint.insert({
                        0, //corresponds to channel bellow
                        std::get<0>(deadChannel)-1,
                        std::get<1>(deadChannel),
                        std::get<2>(deadChannel),
                    });
                    adj_to_dead_Scint.insert({
                        1, //corresponds to channel above
                        std::get<0>(deadChannel)+1,
                        std::get<1>(deadChannel),
                        std::get<2>(deadChannel),
                    });

                  std::vector<std::tuple<int,int,int>> inLayerNeighbors;
                  inLayerNeighbors = getNeighborsScint(deadChannel);
                  int iN = 0;
                  for(auto itr = inLayerNeighbors.begin(); itr!=inLayerNeighbors.end(); ++itr){
                      adj_to_dead_Scint_inlay.insert({
                          iN,
                          std::get<0>(*itr),
                          std::get<1>(*itr),
                          std::get<2>(*itr),
                      });
                      iN++;
                  }

                  std::array<float, 32> temp_vector;
                  for(unsigned k(0); k < 32; ++k) temp_vector[k] = 0;
                  temp_vector[0] = (float)lr; //layer
                  temp_vector[1] = (float)ie; //dead channel's ieta
                  temp_vector[2] = (float)ip; //dead cell's iphi
                  MLvectorevScint.push_back(temp_vector);
                }
            }
        }
    }

    /* This extra vector makes sure the information is passed even if there are
    ** no available dead rechits.
    */
    std::array<float, 29> buffer_vector_Si;
    for(unsigned k(0); k < 29; ++k) buffer_vector_Si[k] = -1;
    MLvectorevSi.push_back(buffer_vector_Si);

    std::array<float, 32> buffer_vector_Scint;
    for(unsigned k(0); k < 32; ++k) buffer_vector_Scint[k] = -1;
    MLvectorevScint.push_back(buffer_vector_Scint);

    std::cout << "List of dead Si cells and scintillator channels was created successfully. \n"
    << "Killed " << N_try_success_Si << " Si cells using " << N_try_all_Si << " trials\n"
    << "and " << N_try_success_Scint << " scintillator channels using " << N_try_all_Scint << " trials.\n"
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
    ULong64_t event = 0;
    std::vector<float> *rechitEnergy = 0;
    std::vector<float> *rechitEta = 0;
    std::vector<float> *rechitPhi = 0;
    std::vector<float> *rechitPosx = 0;
    std::vector<float> *rechitPosy = 0;
    std::vector<float> *rechitPosz = 0;
    std::vector<int> *rechitLayer = 0;
    std::vector<int> *rechitWaferU = 0;
    std::vector<int> *rechitWaferV = 0;
    std::vector<int> *rechitCellU = 0;
    std::vector<int> *rechitCellV = 0;
    std::vector<int> *rechitThickness = 0;
    std::vector<float> *genEta = 0;
    std::vector<float> *genPhi = 0;

    lRecTree->SetBranchAddress("event" ,&event);
    lRecTree->SetBranchAddress("rechit_energy" ,&rechitEnergy);
    lRecTree->SetBranchAddress("rechit_eta" ,&rechitEta);
    lRecTree->SetBranchAddress("rechit_phi" ,&rechitPhi);
    lRecTree->SetBranchAddress("rechit_x" ,&rechitPosx);
    lRecTree->SetBranchAddress("rechit_y" ,&rechitPosy);
    lRecTree->SetBranchAddress("rechit_z" ,&rechitPosz);
    lRecTree->SetBranchAddress("rechit_layer" ,&rechitLayer);
    lRecTree->SetBranchAddress("rechit_wafer_u" ,&rechitWaferU);
    lRecTree->SetBranchAddress("rechit_wafer_v" ,&rechitWaferV);
    lRecTree->SetBranchAddress("rechit_cell_u" ,&rechitCellU);
    lRecTree->SetBranchAddress("rechit_cell_v" ,&rechitCellV);
    lRecTree->SetBranchAddress("rechit_thickness" ,&rechitThickness);
    lRecTree->SetBranchAddress("gen_eta" ,&genEta);
    lRecTree->SetBranchAddress("gen_phi" ,&genPhi);

    unsigned ievtRec = 0;

    // Loop over entries (events)
    for (unsigned ievt(0); ievt<nEvts; ++ievt){
        // Flush MLvectorevSi and MLvectorevScint contents while keeping the list of dead cells intact
        for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); itr++) {
            for(unsigned k(5); k < 29; ++k) (*itr)[k] = 0;
        }
        for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); itr++) {
            for(unsigned k(3); k < 32; ++k) (*itr)[k] = 0;
        }

        if (ievtRec>=lRecTree->GetEntries()) continue;
        Long64_t local_entry = lRecTree->LoadTree(ievt);

        if (debug) std::cout << std::endl<<std::endl << "... Processing entry: " << ievt << std::endl;
        else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

        if (local_entry < 0) continue;
        if (local_entry == 0) {
            lRecTree->SetBranchAddress("event" ,&event);
            lRecTree->SetBranchAddress("rechit_energy" ,&rechitEnergy);
            lRecTree->SetBranchAddress("rechit_eta" ,&rechitEta);
            lRecTree->SetBranchAddress("rechit_phi" ,&rechitPhi);
            lRecTree->SetBranchAddress("rechit_x" ,&rechitPosx);
            lRecTree->SetBranchAddress("rechit_y" ,&rechitPosy);
            lRecTree->SetBranchAddress("rechit_z" ,&rechitPosz);
            lRecTree->SetBranchAddress("rechit_layer" ,&rechitLayer);
            lRecTree->SetBranchAddress("rechit_wafer_u" ,&rechitWaferU);
            lRecTree->SetBranchAddress("rechit_wafer_v" ,&rechitWaferV);
            lRecTree->SetBranchAddress("rechit_cell_u" ,&rechitCellU);
            lRecTree->SetBranchAddress("rechit_cell_v" ,&rechitCellV);
            lRecTree->SetBranchAddress("rechit_thickness" ,&rechitThickness);
            lRecTree->SetBranchAddress("gen_eta" ,&genEta);
            lRecTree->SetBranchAddress("gen_phi" ,&genPhi);
        }

        lRecTree->GetEntry(ievtRec);

        double etagen = 99999.;
        double phigen = 99999.;
        if((*genEta).size() > 0) {
            etagen = (*genEta)[0];
            phigen = (*genPhi)[0];
        }

        double coneSize = 0.15;
        MLrechitsum = 0;

        // Loop over hits of event
        for (unsigned iH(0); iH < (*rechitEnergy).size(); ++iH){
            int layer = (*rechitLayer)[iH];
            double zh = (*rechitPosz)[iH];
            double lenergy = (*rechitEnergy)[iH];
            double leta = (*rechitEta)[iH];
            double lphi = (*rechitPhi)[iH];
            double dR = DeltaR(etagen,phigen,leta,lphi);

            int waferU = (*rechitWaferU)[iH];
            int waferV = (*rechitWaferV)[iH];
            int cellU = (*rechitCellU)[iH];
            int cellV = (*rechitCellV)[iH];
            int thickness = (*rechitThickness)[iH];
            if (ievt == 1) std::cout << thickness << std::endl;
            bool isDense = (thickness == 120) ? 1 : 0;
            bool isScint = (thickness > 400) ? 1 : 0;
            int ieta = (isScint) ? waferU : std::numeric_limits<int>::max();
            int iphi = (isScint) ? waferV : std::numeric_limits<int>::max();

            /* Select hits that are:
            **     - in Si cells
            **     - within DeltaR < 0.15 wrt gen particle
            **     - in positive endcap
            */
            if(!isScint && zh > 0 && dR < coneSize) {
                std::tuple<int, int, int, int, int> tempsi(layer,waferU,waferV,cellU,cellV);
                std::set<std::tuple<int, int, int, int, int>>::iterator ibc=deadlistSi.find(tempsi);
                bool isDead = false;

                // Calculate energy without dead Si cells
                if(ibc == deadlistSi.end()) {
                    MLrechitsum += lenergy;
                }else {
                    // Do stuff with dead cells
                    /* ML code
                    ** Input dead cells eta, phi and rechits
                    */
                    isDead = true;
                    for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); itr++) {
                        if( (*itr)[0] == layer &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[5] = leta;
                            (*itr)[6] = lphi;
                            (*itr)[13] = lenergy;
                            (*itr)[28] = thickness;
                        }
                    }
                }

                /* Get rechits for adjacent cells in neighboring layers.
                ** If cells are also dead then assign a rechit value of -100
                */
                std::tuple<int, int, int, int, int, int> tempsiU(
                    1,layer,waferU,waferV,cellU,cellV
                );
                std::tuple<int, int, int, int, int, int> tempsiD(
                    0,layer,waferU,waferV,cellU,cellV
                );
                std::set<std::tuple<int, int, int, int, int, int>>::iterator itrU=adj_to_dead_Si.find(tempsiU);
                std::set<std::tuple<int, int, int, int, int, int>>::iterator itrD=adj_to_dead_Si.find(tempsiD);
                if(itrU!=adj_to_dead_Si.end()) {
                    for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); itr++) {
                        if( (*itr)[0] == layer-1 &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[14] = lenergy;
                            if(isDead) (*itr)[14] = -100;
                        }
                    }
                }
                if(itrD!=adj_to_dead_Si.end()) {
                    for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); itr++) {
                        if( (*itr)[0] == layer+1 &&
                            (*itr)[1] == waferU && (*itr)[2] == waferV &&
                            (*itr)[3] == cellU  && (*itr)[4] == cellV
                        ){
                            (*itr)[15] = lenergy;
                            if(isDead) (*itr)[15] = -100;
                        }
                    }
                }

                /* Get rechits of same layer neighbors for the dead cell
                ** and its adjacent cell neighbors.
                */
                for(int n = 0; n < 6; ++n){
                    // Same layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiNn(
                        n,layer,waferU,waferV,cellU,cellV
                    );
                    //Check if cell is an nth neighbor of some dead cell
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrNn=adj_to_dead_Si_inlay.find(tempsiNn);
                    if(itrNn!=adj_to_dead_Si_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int>> sameLayerNeighbors;
                        sameLayerNeighbors = getNeighborsSi(tempsi, isDense);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrNn)+3)%6;
                        std::tuple<int, int, int, int> deadCell;
                        std::get<0>(deadCell) = std::get<1>(sameLayerNeighbors[nn]);
                        std::get<1>(deadCell) = std::get<2>(sameLayerNeighbors[nn]);
                        std::get<2>(deadCell) = std::get<3>(sameLayerNeighbors[nn]);
                        std::get<3>(deadCell) = std::get<4>(sameLayerNeighbors[nn]);
                        for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); itr++) {
                            if( (*itr)[0] == layer &&
                                (*itr)[1] == std::get<0>(deadCell) && (*itr)[2] == std::get<1>(deadCell) &&
                                (*itr)[3] == std::get<2>(deadCell) && (*itr)[4] == std::get<3>(deadCell)
                            ){
                                (*itr)[n+7] = lenergy;
                                if(isDead) (*itr)[n+7] = -100;
                            }
                        }
                    }

                    // Next layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiUNn(
                        n,layer-1,waferU,waferV,cellU,cellV
                    );
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrUNn=adj_to_dead_Si_inlay.find(tempsiUNn);
                    if(itrUNn!=adj_to_dead_Si_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int>> nextLayerNeighbors;
                        nextLayerNeighbors = getNeighborsSi(tempsi, isDense);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrUNn)+3)%6;
                        std::tuple<int, int, int, int> deadCell;
                        std::get<0>(deadCell) = std::get<1>(nextLayerNeighbors[nn]);
                        std::get<1>(deadCell) = std::get<2>(nextLayerNeighbors[nn]);
                        std::get<2>(deadCell) = std::get<3>(nextLayerNeighbors[nn]);
                        std::get<3>(deadCell) = std::get<4>(nextLayerNeighbors[nn]);
                        for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); itr++) {
                            if( (*itr)[0] == layer-1 &&
                            (*itr)[1] == std::get<0>(deadCell) && (*itr)[2] == std::get<1>(deadCell) &&
                            (*itr)[3] == std::get<2>(deadCell) && (*itr)[4] == std::get<3>(deadCell)
                            ){
                                (*itr)[n+16] = lenergy;
                                if(isDead) (*itr)[n+16] = -100;
                            }
                        }
                    }

                    // Previous layer neighbors
                    std::tuple<int, int, int, int, int, int> tempsiDNn(
                        n,layer+1,waferU,waferV,cellU,cellV
                    );
                    std::set<std::tuple<int, int, int, int, int, int>>::iterator itrDNn=adj_to_dead_Si_inlay.find(tempsiDNn);
                    if(itrDNn!=adj_to_dead_Si_inlay.end()) {
                        std::vector<std::tuple<int,int,int,int,int>> prevLayerNeighbors;
                        prevLayerNeighbors = getNeighborsSi(tempsi, isDense);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrDNn)+3)%6;
                        std::tuple<int, int, int, int> deadCell;
                        std::get<0>(deadCell) = std::get<1>(prevLayerNeighbors[nn]);
                        std::get<1>(deadCell) = std::get<2>(prevLayerNeighbors[nn]);
                        std::get<2>(deadCell) = std::get<3>(prevLayerNeighbors[nn]);
                        std::get<3>(deadCell) = std::get<4>(prevLayerNeighbors[nn]);
                        for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); itr++) {
                            if( (*itr)[0] == layer+1 &&
                            (*itr)[1] == std::get<0>(deadCell) && (*itr)[2] == std::get<1>(deadCell) &&
                            (*itr)[3] == std::get<2>(deadCell) && (*itr)[4] == std::get<3>(deadCell)
                            ){
                                (*itr)[n+22] = lenergy;
                                if(isDead) (*itr)[n+22] = -100;
                            }
                        }
                    }
                }
            }

            /* Select hits that are:
            **     - in Scint channels
            **     - within DeltaR < 0.15 wrt gen particle
            **     - in positive endcap
            */
            //if (ievt==1) std::cout << "Check 1: isScint = " << isScint << ", zh = " << zh << ", dR = " << dR << std::endl;
            if(isScint && zh > 0 && dR < coneSize) {
                std::cout << "Check 2" << std::endl;
                std::tuple<int, int, int> tempscint(layer,ieta,iphi);
                std::cout << "Check 3" << std::endl;
                std::set<std::tuple<int, int, int>>::iterator ibc=deadlistScint.find(tempscint);
                bool isDead = false;

                std::cout << "Check 4: ieta = " << ieta << ", iphi = " << iphi << std::endl;
                // Calculate energy without dead Scint channels
                if(ibc == deadlistScint.end()) {
                    MLrechitsum += lenergy;
                }else {
                    std::cout << "Check 2" << std::endl;
                    // Do stuff with dead channels
                    /* ML code
                    ** Input dead channels eta, phi and rechits
                    */
                    isDead = true;
                    for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); itr++) {
                        if( (*itr)[0] == layer && (*itr)[1] == ieta && (*itr)[2] == iphi ){
                            (*itr)[3] = leta;
                            (*itr)[4] = lphi;
                            (*itr)[13] = lenergy;
                        }
                    }
                }

                /* Get rechits for adjacent channels in neighboring layers.
                ** If channels are also dead then assign a rechit value of -100
                */
                std::tuple<int, int, int, int> tempscintU(1,layer,ieta,iphi);
                std::tuple<int, int, int, int> tempscintD(0,layer,ieta,iphi);
                std::set<std::tuple<int, int, int, int>>::iterator itrU=adj_to_dead_Scint.find(tempscintU);
                std::set<std::tuple<int, int, int, int>>::iterator itrD=adj_to_dead_Scint.find(tempscintD);
                if(itrU!=adj_to_dead_Scint.end()) {
                    for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); itr++) {
                        if( (*itr)[0] == layer-1 && (*itr)[1] == ieta && (*itr)[2] == iphi ){
                            (*itr)[14] = lenergy;
                            if(isDead) (*itr)[14] = -100;
                        }
                    }
                }
                if(itrD!=adj_to_dead_Scint.end()) {
                    for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); itr++) {
                        if( (*itr)[0] == layer+1 && (*itr)[1] == ieta && (*itr)[2] == iphi ){
                            (*itr)[15] = lenergy;
                            if(isDead) (*itr)[15] = -100;
                        }
                    }
                }

                /* Get rechits of same layer neighbors for the dead channel
                ** and its adjacent channel neighbors.
                */
                for(int n = 0; n < 8; ++n){
                    // Same layer neighbors
                    std::tuple<int, int, int, int> tempscintNn(n,layer,ieta,iphi);
                    //Check if channel is an nth neighbor of some dead channel
                    std::set<std::tuple<int, int, int, int>>::iterator itrNn=adj_to_dead_Scint_inlay.find(tempscintNn);
                    if(itrNn!=adj_to_dead_Scint_inlay.end()) {
                        std::vector<std::tuple<int,int,int>> sameLayerNeighbors;
                        sameLayerNeighbors = getNeighborsScint(tempscint);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrNn)+4)%8;
                        std::tuple<int, int> deadchannel;
                        std::get<0>(deadchannel) = std::get<1>(sameLayerNeighbors[nn]);
                        std::get<1>(deadchannel) = std::get<2>(sameLayerNeighbors[nn]);
                        for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); itr++) {
                            if( (*itr)[0] == layer && (*itr)[1] == std::get<0>(deadchannel) &&
                                (*itr)[2] == std::get<1>(deadchannel)
                            ){
                                (*itr)[n+5] = lenergy;
                                if(isDead) (*itr)[n+5] = -100;
                            }
                        }
                    }

                    // Next layer neighbors
                    std::tuple<int, int, int, int> tempscintUNn(n,layer-1,ieta,iphi);
                    std::set<std::tuple<int, int, int, int>>::iterator itrUNn=adj_to_dead_Scint_inlay.find(tempscintUNn);
                    if(itrUNn!=adj_to_dead_Scint_inlay.end()) {
                        std::vector<std::tuple<int,int,int>> nextLayerNeighbors;
                        nextLayerNeighbors = getNeighborsScint(tempscint);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrUNn)+4)%8;
                        std::tuple<int, int> deadchannel;
                        std::get<0>(deadchannel) = std::get<1>(nextLayerNeighbors[nn]);
                        std::get<1>(deadchannel) = std::get<2>(nextLayerNeighbors[nn]);
                        for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); itr++) {
                            if( (*itr)[0] == layer-1 && (*itr)[1] == std::get<0>(deadchannel) &&
                                (*itr)[2] == std::get<1>(deadchannel)
                            ){
                                (*itr)[n+16] = lenergy;
                                if(isDead) (*itr)[n+16] = -100;
                            }
                        }
                    }

                    // Previous layer neighbors
                    std::tuple<int, int, int, int> tempscintDNn(n,layer+1,ieta,iphi);
                    std::set<std::tuple<int, int, int, int>>::iterator itrDNn=adj_to_dead_Scint_inlay.find(tempscintDNn);
                    if(itrDNn!=adj_to_dead_Scint_inlay.end()) {
                        std::vector<std::tuple<int,int,int>> prevLayerNeighbors;
                        prevLayerNeighbors = getNeighborsScint(tempscint);
                        // Get neighbor number
                        int nn = (std::get<0>(*itrDNn)+4)%8;
                        std::tuple<int, int> deadchannel;
                        std::get<0>(deadchannel) = std::get<1>(prevLayerNeighbors[nn]);
                        std::get<1>(deadchannel) = std::get<2>(prevLayerNeighbors[nn]);
                        for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); itr++) {
                            if( (*itr)[0] == layer+1 && (*itr)[1] == std::get<0>(deadchannel) &&
                                (*itr)[2] == std::get<1>(deadchannel)
                            ){
                                (*itr)[n+24] = lenergy;
                                if(isDead) (*itr)[n+24] = -100;
                            }
                        }
                    }
                }
            }
        }

        //Export the ML dataset values to the TTree
        for(auto itr = MLvectorevSi.begin(); itr != MLvectorevSi.end(); ++itr) {
            if ((*itr)[5] > 0 || (*itr)[0]==-1) {
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
                MLdn4    = (*itr)[25];
                MLdn5    = (*itr)[26];
                MLdn6    = (*itr)[27];
                MLthickness = (*itr)[28];
                MLevent  = (float)(event);
                t1->Fill();
            }
        }
        for(auto itr = MLvectorevScint.begin(); itr != MLvectorevScint.end(); ++itr) {
            if ((*itr)[3] > 0 || (*itr)[0]==-1) {
                /* This condition is necessary to ensure the cell was within
                ** the cone.
                */
                MLlayer  = (*itr)[0];
                MLieta   = (*itr)[1];
                MLiphi   = (*itr)[2];
                MLeta    = (*itr)[3];
                MLphi    = (*itr)[4];
                MLn1     = (*itr)[5];
                MLn2     = (*itr)[6];
                MLn3     = (*itr)[7];
                MLn4     = (*itr)[8];
                MLn5     = (*itr)[9];
                MLn6     = (*itr)[10];
                MLn7     = (*itr)[11];
                MLn8     = (*itr)[12];
                MLdead   = (*itr)[13];
                MLnup    = (*itr)[14];
                MLndown  = (*itr)[15];
                MLun1    = (*itr)[16];
                MLun2    = (*itr)[17];
                MLun3    = (*itr)[18];
                MLun4    = (*itr)[19];
                MLun5    = (*itr)[20];
                MLun6    = (*itr)[21];
                MLun7    = (*itr)[22];
                MLun8    = (*itr)[23];
                MLdn1    = (*itr)[24];
                MLdn2    = (*itr)[25];
                MLdn3    = (*itr)[26];
                MLdn4    = (*itr)[27];
                MLdn5    = (*itr)[28];
                MLdn6    = (*itr)[29];
                MLdn7    = (*itr)[30];
                MLdn8    = (*itr)[31];
                MLevent  = (float)(event);
                t2->Fill();
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
