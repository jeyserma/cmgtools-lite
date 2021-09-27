#ifndef FUNCTIONS_LOWPU_H
#define FUNCTIONS_LOWPU_H

#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TVector2.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib> //as stdlib.h      
#include <cstdio>
#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include "defines.h"
#include <limits>
#include <map>

#include "RoccoR.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;

// lepton SF
/*
TFile *_leptonSF = NULL;
TH2 * lepSF_mu_SIT = NULL;
TH2 * lepSF_mu_STA = NULL;
TH2 * lepSF_mu_HLT_pos = NULL;
TH2 * lepSF_mu_HLT_neg = NULL;
*/

TFile *_leptonSF = new TFile("lowPU/leptonSF/lepton_SF.root", "READ");;
TH2 * lepSF_mu_SIT = (TH2F*)(_leptonSF->Get("mu_SIT"));
TH2 * lepSF_mu_STA = (TH2F*)(_leptonSF->Get("mu_STA"));
TH2 * lepSF_mu_HLT_pos = (TH2F*)(_leptonSF->Get("mu_HLT_pos"));
TH2 * lepSF_mu_HLT_neg = (TH2F*)(_leptonSF->Get("mu_HLT_neg"));


TH2 * lepSF_mu_HLT_pos_DATA = (TH2F*)(_leptonSF->Get("mu_HLT_pos_DATA"));
TH2 * lepSF_mu_HLT_neg_DATA = (TH2F*)(_leptonSF->Get("mu_HLT_neg_DATA"));
TH2 * lepSF_mu_HLT_pos_MC = (TH2F*)(_leptonSF->Get("mu_HLT_pos_MC"));
TH2 * lepSF_mu_HLT_neg_MC = (TH2F*)(_leptonSF->Get("mu_HLT_neg_MC"));

RoccoR * rochester = new RoccoR("lowPU/RoccoR/RoccoR2017.txt"); // https://gitlab.cern.ch/akhukhun/roccor

Vec_b testTrigger(Vec_f eta, Vec_f phi, Vec_f TrigObj_eta, Vec_f TrigObj_phi) {
    
   Vec_b res(eta.size(),false); // initialize to 0

   return res;

}


double applyMuonSF(Vec_f pt, Vec_f eta, Vec_i q, Vec_b trg) {


    /*
    if(!_leptonSF) {
        
        _leptonSF = new TFile("lowPU/leptonSF/lepton_SF.root", "READ");
        lepSF_mu_SIT = (TH2F*)(_leptonSF->Get("mu_SIT"));
        lepSF_mu_STA = (TH2F*)(_leptonSF->Get("mu_STA"));
        lepSF_mu_HLT_pos = (TH2F*)(_leptonSF->Get("mu_HLT_pos"));
        lepSF_mu_HLT_neg = (TH2F*)(_leptonSF->Get("mu_HLT_neg"));
    }

    */

    unsigned int size = pt.size();
    double corr = 1.0;
  
    
    bool trg_applied = false;
    for (unsigned int i = 0; i < size; ++i) {
        
        double eta_ = eta.at(i);
        double pt_ = pt.at(i);
        int q_ = q.at(i);

        int etabin_SIT = std::max(1, std::min(lepSF_mu_SIT->GetNbinsX(), lepSF_mu_SIT->GetXaxis()->FindBin(eta_)));
        int ptbin_SIT = std::max(1, std::min(lepSF_mu_SIT->GetNbinsY(), lepSF_mu_SIT->GetYaxis()->FindBin(pt_)));
        
        int etabin_STA = std::max(1, std::min(lepSF_mu_STA->GetNbinsX(), lepSF_mu_STA->GetXaxis()->FindBin(eta_)));
        int ptbin_STA = std::max(1, std::min(lepSF_mu_STA->GetNbinsY(), lepSF_mu_STA->GetYaxis()->FindBin(pt_)));
        
        int etabin_HLT = std::max(1, std::min(lepSF_mu_HLT_pos->GetNbinsX(), lepSF_mu_HLT_pos->GetXaxis()->FindBin(eta_)));
        int ptbin_HLT = std::max(1, std::min(lepSF_mu_HLT_pos->GetNbinsY(), lepSF_mu_HLT_pos->GetYaxis()->FindBin(pt_)));
        
        
        //cout << "etabin_SIT=" << etabin_SIT << " ptbin_SIT=" << ptbin_SIT << " eta=" << eta_ << " pt=" << pt_ << endl;
        
        corr *= lepSF_mu_SIT->GetBinContent(etabin_SIT, ptbin_SIT);
        corr *= lepSF_mu_STA->GetBinContent(etabin_STA, ptbin_STA);
        
        if(trg.at(i) and not trg_applied) {
            if(q_ > 0) corr *= lepSF_mu_HLT_pos->GetBinContent(etabin_HLT, ptbin_HLT);
            else corr *= lepSF_mu_HLT_neg->GetBinContent(etabin_HLT, ptbin_HLT);
            trg_applied = true;
        }
        
        /*
        double trgSF_DATA = 1.0;
        double trgSF_MC = 1.0;
        if(q_ > 0) {
            trgSF_DATA *= 1. - lepSF_mu_HLT_pos_DATA->GetBinContent(etabin_HLT, ptbin_HLT);
            trgSF_MC *= 1. - lepSF_mu_HLT_pos_MC->GetBinContent(etabin_HLT, ptbin_HLT);
        }
        else {
                
            trgSF_DATA *= 1. - lepSF_mu_HLT_neg_DATA->GetBinContent(etabin_HLT, ptbin_HLT);
            trgSF_MC *= 1. - lepSF_mu_HLT_neg_MC->GetBinContent(etabin_HLT, ptbin_HLT);
        }
        corr *= (1. - trgSF_DATA) / (1. - trgSF_MC);
        */      
        
        
    }
    
    //cout << "corr=" << corr << endl;
    return corr;

}


Vec_f applyRochesterMC(Vec_f pt, Vec_f eta, Vec_f phi, Vec_f ch, Vec_i gen_idx, Vec_f gen_pt, Vec_i nTrackerLayers) {
    
    unsigned int size = pt.size();
    Vec_f res(size, 0.0);
  
    // https://gitlab.cern.ch/akhukhun/roccor
    for (unsigned int i = 0; i < size; ++i) {
            
        if(gen_idx.at(i) >= 0) res[i] = pt[i] * rochester->kSpreadMC(ch[i], pt[i], eta[i], phi[i], gen_pt.at(gen_idx.at(i)), 0, 0);
        else res[i] = pt[i] * rochester->kSmearMC(ch[i], pt[i], eta[i], phi[i], nTrackerLayers.at(i), gRandom->Rndm(), 0, 0);
    }        

    return res;

}


Vec_f applyRochesterData(Vec_f pt, Vec_f eta, Vec_f phi, Vec_f ch) {
    
    unsigned int size = pt.size();
    Vec_f res(size, 0.0);
  
    // https://gitlab.cern.ch/akhukhun/roccor
    for (unsigned int i = 0; i < size; ++i) {
            
        res[i] = pt[i] * rochester->kScaleDT(ch[i], pt[i], eta[i], phi[i]);
    }
    return res;

}


Vec_f sortPt(Vec_f pt) {

    unsigned int size = pt.size();
    Vec_f res(size, 0.0); // 2 elements initialized to 0
  
    for (unsigned int i = 0; i < size; ++i)  res[i] = pt[i];
    std::sort(res.begin(), res.end());
    
    return res;
}


double leadingPt(Vec_f pt) {

    unsigned int size = pt.size();
    if(size == 0) return -1;
    Vec_f res(size, 0.0);
  
    for (unsigned int i = 0; i < size; ++i)  res[i] = pt[i];
    std::sort(res.begin(), res.end());
    
    return res.at(0);
}

double subLeadingPt(Vec_f pt) {

    unsigned int size = pt.size();
    if(size <= 1) return -1;
    Vec_f res(size, 0.0);
  
    for (unsigned int i = 0; i < size; ++i)  res[i] = pt[i];
    std::sort(res.begin(), res.end());
    
    return res.at(1);
}

Vec_b goodMuonTriggerCandidateLowPU(const Vec_i& TrigObj_id, const Vec_f& TrigObj_pt, const Vec_f& TrigObj_l1pt, const Vec_f& TrigObj_l2pt, const Vec_i& TrigObj_filterBits) {

   Vec_b res(TrigObj_id.size(),false); // initialize to 0   
   for (unsigned int i = 0; i < res.size(); ++i) {
       if (TrigObj_id[i]  != 13 ) continue;
       if (TrigObj_pt[i]   < 18.) continue;
       if (TrigObj_l1pt[i] < 16.) continue;
       //if (! (( TrigObj_filterBits[i] & 8) || (TrigObj_l2pt[i] > 10. && (TrigObj_filterBits[i] & 2) )) ) continue;
       res[i] = true;
   }
   // res will be goodTrigObjs in RDF
   // e.g. RDF::Define("goodTrigObjs","goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
   return res;
}


int test(Vec_f pt1) {

    int res = -1;
    if(pt1.at(0) > pt1.at(1)) res = 1;
    else res = 0;
    if(res == 0) std::cout << pt1.at(0) << " " << pt1.at(1) << " " << res << std::endl;
    
    return 0;
}


float transversemomentumBoson(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).Pt();
}






float lowPU_pu_reweight_sto[100] = {1, 0.041879317950684, 0.2906403692601066, 2.926955277913856, 0.604837179882255, 0.13291723160425126, 0.06483328827710513, 0.03399077851969202, 0.01895726346472344, 0.010933516524116909, 0.006426849386948908, 0.0033939880542484023, 0.00174393639715267, 0.0004962683949237972, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
float lowPU_pu_reweight(int nTrueInt) { if (nTrueInt<100) return lowPU_pu_reweight_sto[nTrueInt]; else return 0; }

#endif
