#include "TFile.h"
#include "TTree.h"
#include <iostream>

void print_hits()
{
    TFile* f = new TFile("trkNtup.root");
    TTree* t = (TTree*) ((TDirectoryFile*) f->Get("trackingNtuple"))->Get("tree");

    std::vector<float>* ph2_x = 0;
    std::vector<float>* ph2_y = 0;
    std::vector<float>* ph2_z = 0;
    std::vector<int>* ph2_detId = 0;
    std::vector<int>* ph2_moduleType = 0;

    t->SetBranchAddress("ph2_x", &ph2_x);
    t->SetBranchAddress("ph2_y", &ph2_y);
    t->SetBranchAddress("ph2_z", &ph2_z);
    t->SetBranchAddress("ph2_detId", &ph2_detId);
    t->SetBranchAddress("ph2_moduleType", &ph2_moduleType);

    for (int i = 0; i < t->GetEntries(); ++i)
    {
        t->GetEntry(i);
        for (int ihit = 0; ihit < ph2_x->size(); ++ihit)
        {
            float x = ph2_x->at(ihit);
            float y = ph2_y->at(ihit);
            float z = ph2_z->at(ihit);
            int detId = ph2_detId->at(ihit);
            int moduleType = ph2_moduleType->at(ihit);
            std::cout <<  " x: " << x <<  " y: " << y <<  " z: " << z <<  " detId: " << detId << " moduleType: " << moduleType <<  std::endl;
        }
    }
}
