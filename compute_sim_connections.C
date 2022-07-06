#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "ModuleDetIdParser.cxx"
#include <iostream>

void compute_sim_connections(TString fileDir="/data2/segmentlinking/CMSSW_12_2_0_pre2/trackingNtuple_10mu_pt_0p5_50.root")
{
  bool debug = false;

  TFile* f = new TFile(fileDir,"read");
  TTree* t = (TTree*) ((TDirectoryFile*) f->Get("trackingNtuple"))->Get("tree");

  std::unordered_map<unsigned int,std::set<unsigned int>> mod_map;

  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");

  std::vector<std::vector<int>>* sim_simHitIdx = 0;
  std::vector<unsigned int>* simhit_detId = 0;
  std::vector<unsigned short>* simhit_layer = 0;
  std::vector<unsigned short>* simhit_ring = 0;
  std::vector<unsigned short>* simhit_subdet = 0;

  t->SetBranchAddress("sim_simHitIdx", &sim_simHitIdx);
  t->SetBranchAddress("simhit_detId", &simhit_detId);
  t->SetBranchAddress("simhit_layer", &simhit_layer);
  t->SetBranchAddress("simhit_ring", &simhit_ring);
  t->SetBranchAddress("simhit_subdet", &simhit_subdet);

  int maxEvents = debug ? 1 : t->GetEntries();
  for (int i = 0; i < maxEvents; ++i)
  {
    t->GetEntry(i);
    for (unsigned int isim = 0; isim < sim_simHitIdx->size(); ++isim)
    {
      std::vector<int> simhitIdx = sim_simHitIdx->at(isim);
      if (simhitIdx.size()<1) continue;
      if (debug) std::cout << "New track!\n";

      int  prev_simhitIdx = simhitIdx[0];
      unsigned int  prev_detId = simhit_detId->at(prev_simhitIdx);
      unsigned short  prev_layer = simhit_layer->at(prev_simhitIdx);
      unsigned short  prev_ring = simhit_ring->at(prev_simhitIdx);
      unsigned short  prev_subdet = simhit_subdet->at(prev_simhitIdx);
      unsigned short  prev_isLower = (prev_subdet==5 || prev_subdet==4) ? SDL::Module(prev_detId).isLower() : 0; // Defined only for non-pixel modules

      for (unsigned int isimhit = 1; isimhit < simhitIdx.size(); ++isimhit)
      {
        if (debug && (prev_subdet == 5 || prev_subdet == 4) && prev_isLower != 1) std::cout << prev_detId << " not lower!\n";

        if (!( ( prev_subdet == 5 && prev_isLower == 1 && prev_layer != 6 ) || // Verify that the module is a non-pixel one,
               ( prev_subdet == 4 && prev_isLower == 1 && prev_layer != 5 &&   // is a lower one and 
                 !(prev_ring == 15 && prev_layer == 1) &&                      // is not on the outer side of the tracker.
                 !(prev_ring == 15 && prev_layer == 2) &&                      // If it is not, update to the next hit and move on.
                 !(prev_ring == 12 && prev_layer == 3) &&
                 !(prev_ring == 12 && prev_layer == 4)
               )
             )
           )
        {
          prev_simhitIdx = simhitIdx[isimhit];
          prev_detId = simhit_detId->at(prev_simhitIdx);
          prev_layer = simhit_layer->at(prev_simhitIdx);
          prev_ring = simhit_ring->at(prev_simhitIdx);
          prev_subdet = simhit_subdet->at(prev_simhitIdx);
          prev_isLower = (prev_subdet==5 || prev_subdet==4) ? SDL::Module(prev_detId).isLower() : 0; // Defined only for non-pixel modules
          continue;
        }
        if(debug) std::cout << prev_detId << " : " << prev_isLower << "\n";

        unsigned int cur_detId = simhit_detId->at(simhitIdx[isimhit]);
        if(debug) std::cout << cur_detId << " : " << SDL::Module(cur_detId).isLower() << "\n";

        if (prev_detId == cur_detId || SDL::Module(prev_detId).partnerDetId() == cur_detId ) continue; // Do not consider connections with itself or its partner module. Do not update to the next hit in this case - we want to see the non-tautological connection of this module.

        if(simhit_subdet->at(simhitIdx[isimhit])!=5 && simhit_subdet->at(simhitIdx[isimhit])!=4) // Loopers can create connections from outer to the inner tracker.
        {                                                                                        // We ignore these connections - we update to the next hit and move on.
          if(debug) std::cout << "Connection to pixel modules --> Skipping...\n";
          prev_simhitIdx = simhitIdx[isimhit];
          prev_detId = simhit_detId->at(prev_simhitIdx);
          prev_layer = simhit_layer->at(prev_simhitIdx);
          prev_ring = simhit_ring->at(prev_simhitIdx);
          prev_subdet = simhit_subdet->at(prev_simhitIdx);
          prev_isLower = (prev_subdet==5 || prev_subdet==4) ? SDL::Module(prev_detId).isLower() : 0; // Defined only for non-pixel modules
          continue;
        }
        if(!SDL::Module(cur_detId).isLower()) // Consider connections only with lower modules - if the hit happened to be in the upper module, register the connection to the lower module instead.
        {
          cur_detId = SDL::Module(cur_detId).partnerDetId();
          if(debug) std::cout << "Switch to lower module --> " << cur_detId << " : " << SDL::Module(cur_detId).isLower() << "\n";
        }

        if (mod_map.find(prev_detId) != mod_map.end()) mod_map[prev_detId].insert(cur_detId);
        else
        {
           std::set<unsigned int> set_temp{cur_detId};
           mod_map[prev_detId] = set_temp;
        }
        prev_simhitIdx = simhitIdx[isimhit];
        prev_detId = simhit_detId->at(prev_simhitIdx);
        prev_layer = simhit_layer->at(prev_simhitIdx);
        prev_ring = simhit_ring->at(prev_simhitIdx);
        prev_subdet = simhit_subdet->at(prev_simhitIdx);
        prev_isLower = (prev_subdet==5 || prev_subdet==4) ? SDL::Module(prev_detId).isLower() : 0; // Defined only for non-pixel modules
        if(debug) std::cout << "\n";
      }
    }
  }

  std::cout << "\n";
  for (auto kv : mod_map)
  {
    std::cout << kv.first << " " << kv.second.size();
    for (auto v : kv.second) std::cout << " " << v;
    std::cout << "\n";
  }
}
