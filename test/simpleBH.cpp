/****************************************************
**      This code uses ntuples to produce rechitsums
**      for a given dead Si cell fraction.
**      The output also contains an ntuple of all the
**      dead cells that can be used to train a DNN to
**      estimate the lost energy.
**
**      Edit: This is an version for the DAQ study.
**            Module are switched off instead of
**            individual cells.
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
    **     - MLdead is dead cell rechit
    ** We also need to create a set container to store the values before writing to TTree
    **********************************/
    float MLlayer, MLeta, MLphi, MLdead, MLevent;
    float MLwaferU, MLwaferV, MLcellU, MLcellV;
    float MLieta, MLiphi, MLrechitsum, MLthickness;
    TTree* t1 = new TTree("t1","sample");
    t1->Branch("layer"    ,&MLlayer    ,"layer/F"    );
    t1->Branch("waferU"   ,&MLwaferU   ,"waferU/F"   );
    t1->Branch("waferV"   ,&MLwaferV   ,"waferV/F"   );
    t1->Branch("cellU"    ,&MLcellU    ,"cellU/F"    );
    t1->Branch("cellV"    ,&MLcellV    ,"cellV/F"    );
    t1->Branch("ieta"     ,&MLieta     ,"ieta/F"     );
    t1->Branch("iphi"     ,&MLiphi     ,"iphi/F"     );
    t1->Branch("eta"      ,&MLeta      ,"eta/F"      );
    t1->Branch("phi"      ,&MLphi      ,"phi/F"      );
    t1->Branch("dead"     ,&MLdead     ,"dead/F"     );
    t1->Branch("event"    ,&MLevent    ,"event/F"    );
    t1->Branch("rechitsum",&MLrechitsum,"rechitsum/F");
    t1->Branch("thickness",&MLthickness,"thickness/F");

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
        /*
        ** Define a vector of the array:
        ** {dead cell:
        **      layer, waferU, waferV, cellU, cellV,
        **      ieta, iphim, eta, phi, rechit, thickness
        ** }
        ** Include a buffer vector that makes sure the information is passed even if there are
        ** no available dead rechits.
        */
        std::vector<std::array<float, 11>> MLvectorev;
        std::array<float, 11> buffer_vector;
        for(unsigned k(0); k < 11; ++k) buffer_vector[k] = -1;
        MLvectorev.push_back(buffer_vector);

        if (ievtRec>=lRecTree->GetEntries()) continue;
        Long64_t local_entry = lRecTree->LoadTree(ievt);

        if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

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
            bool isDense = (thickness == 120) ? 1 : 0;
            bool isScint = (thickness != 120 && thickness != 200 && thickness != 300) ? 1 : 0;
            int ieta = (isScint) ? waferU : std::numeric_limits<int>::max();
            int iphi = (isScint) ? waferV : std::numeric_limits<int>::max();

            /* Select hits that are:
            **     - in Si cells
            **     - within DeltaR < 0.15 wrt gen particle
            **     - in positive endcap
            */
            if(!isScint && zh > 0 && dR < coneSize) {
                // Calculate energy without dead Si cells
                if(/* Some condition that sets the cell on or off */) {
                    MLrechitsum += lenergy;
                }else {
                    std::array<float, 11> temp_vector;
                    for(unsigned k(0); k < 11; ++k) temp_vector[k] = -1;
                    temp_vector[0] = (float)lr; //layer
                    temp_vector[1] = (float)waferU;  //dead cell's waferU
                    temp_vector[2] = (float)waferV;  //dead cell's waferV
                    temp_vector[3] = (float)cellU;  //dead cell's cellU
                    temp_vector[4] = (float)cellV;  //dead cell's cellV
                    temp_vector[7] = leta;  //dead cell's eta
                    temp_vector[8] = lphi;  //dead cell's phi
                    temp_vector[9] = lenergy;  //dead cell's energy
                    temp_vector[10] = thickness;  //dead cell's thickness
                    MLvectorev.push_back(temp_vector);
                }
            }

            /* Select hits that are:
            **     - in Scint channels
            **     - within DeltaR < 0.15 wrt gen particle
            **     - in positive endcap
            */
            if(isScint && zh > 0 && dR < coneSize) {
                // Calculate energy without dead Scint channels
                if(/* Some condition that sets the cell on or off */) {
                    MLrechitsum += lenergy;
                }else {
                    std::array<float, 11> temp_vector;
                    for(unsigned k(0); k < 11; ++k) temp_vector[k] = -1;
                    temp_vector[0] = (float)lr; //layer
                    temp_vector[5] = (float)ieta;  //dead cell's ieta
                    temp_vector[6] = (float)iphi;  //dead cell's iphi
                    temp_vector[7] = leta;  //dead cell's eta
                    temp_vector[8] = lphi;  //dead cell's phi
                    temp_vector[9] = lenergy;  //dead cell's energy
                    MLvectorev.push_back(temp_vector);
                }
            }
        }

        //Export the ML dataset values to the TTree
        for(auto itr = MLvectorev.begin(); itr != MLvectorev.end(); ++itr) {
            if ((*itr)[7] > 0 || (*itr)[0]==-1) {
                /* This condition is necessary to ensure the cell was within
                ** the cone.
                */
                MLlayer  = (*itr)[0];
                MLwaferU = (*itr)[1];
                MLwaferV = (*itr)[2];
                MLcellU  = (*itr)[3];
                MLcellV  = (*itr)[4];
                MLieta   = (*itr)[5];
                MLiphi   = (*itr)[6];
                MLeta    = (*itr)[7];
                MLphi    = (*itr)[8];
                MLdead   = (*itr)[9];
                MLthickness = (*itr)[10];
                MLevent  = (float)(event);
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
