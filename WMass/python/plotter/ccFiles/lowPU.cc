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

TFile *_leptonSF = new TFile("lowPU/leptonSF/lepton_SF.root", "READ");
TH2 * lepSF_mu_SIT = (TH2F*)(_leptonSF->Get("mu_SIT"));
TH2 * lepSF_mu_STA = (TH2F*)(_leptonSF->Get("mu_STA"));
TH2 * lepSF_mu_HLT_pos = (TH2F*)(_leptonSF->Get("mu_HLT_pos"));
TH2 * lepSF_mu_HLT_neg = (TH2F*)(_leptonSF->Get("mu_HLT_neg"));

TH2 * lepSF_el_GSFSI = (TH2F*)(_leptonSF->Get("el_GSFSI"));
TH2 * lepSF_el_HLT_pos = (TH2F*)(_leptonSF->Get("el_HLT_pos"));
TH2 * lepSF_el_HLT_neg = (TH2F*)(_leptonSF->Get("el_HLT_neg"));


TH2 * lepSF_mu_HLT_pos_DATA = (TH2F*)(_leptonSF->Get("mu_HLT_pos_DATA"));
TH2 * lepSF_mu_HLT_neg_DATA = (TH2F*)(_leptonSF->Get("mu_HLT_neg_DATA"));
TH2 * lepSF_mu_HLT_pos_MC = (TH2F*)(_leptonSF->Get("mu_HLT_pos_MC"));
TH2 * lepSF_mu_HLT_neg_MC = (TH2F*)(_leptonSF->Get("mu_HLT_neg_MC"));

TH2 * lepSF_el_HLT_pos_DATA = (TH2F*)(_leptonSF->Get("el_HLT_pos_DATA"));
TH2 * lepSF_el_HLT_neg_DATA = (TH2F*)(_leptonSF->Get("el_HLT_neg_DATA"));
TH2 * lepSF_el_HLT_pos_MC = (TH2F*)(_leptonSF->Get("el_HLT_pos_MC"));
TH2 * lepSF_el_HLT_neg_MC = (TH2F*)(_leptonSF->Get("el_HLT_neg_MC"));

RoccoR * rochester = new RoccoR("lowPU/RoccoR/RoccoR2017.txt"); // https://gitlab.cern.ch/akhukhun/roccor


// electron SF
TFile *_electronSF = new TFile("lowPU/leptonSF/egammaEffi.txt_EGM2D_Medium_UL17.root", "READ");
TH2 * lepSF_el_ID = (TH2F*)(_electronSF->Get("EGamma_SF2D"));

// pre-fire
TFile *_preFire = new TFile("lowPU/preFire/All2017Gand2017HPrefiringMaps.root", "READ");
TH2 * preFire_jet = (TH2F*)(_preFire->Get("L1prefiring_jetpt_2017H"));
TH2 * preFire_photon = (TH2F*)(_preFire->Get("L1prefiring_photonpt_2017H"));

Vec_b testTrigger(Vec_f eta, Vec_f phi, Vec_f TrigObj_eta, Vec_f TrigObj_phi) {
    
   Vec_b res(eta.size(),false); // initialize to 0

   return res;

}

// nano scalefactors
//auto csetEl = correction::CorrectionSet::from_file("electron.json");


double prefireCorr_getPrefireProbability(int type, double eta, double pt, double maxpt) {
    
    double pref_prob = 0.0;
    double pt_ = std::min(pt, maxpt-0.01);
    
 
    if(type == 0) { // jet map
        int etabin = std::max(1, std::min(preFire_jet->GetNbinsX(), preFire_jet->GetXaxis()->FindBin(eta)));
        int ptbin = std::max(1, std::min(preFire_jet->GetNbinsY(), preFire_jet->GetYaxis()->FindBin(pt_)));
        pref_prob = preFire_jet->GetBinContent(etabin, ptbin);
    }
    else { // photon map

        int etabin = std::max(1, std::min(preFire_photon->GetNbinsX(), preFire_photon->GetXaxis()->FindBin(eta)));
        int ptbin = std::max(1, std::min(preFire_photon->GetNbinsY(), preFire_photon->GetYaxis()->FindBin(pt_)));
        pref_prob = preFire_photon->GetBinContent(etabin, ptbin);
    }

    
    return pref_prob;
}

double prefireCorr_eg(int jid, Vec_f photon_pt, Vec_f photon_eta, Vec_i photon_jetIdx, Vec_i photon_electronIdx, Vec_f electron_pt, Vec_f electron_eta, Vec_i electron_jetIdx, Vec_i electron_photonIdx) {

    double phopf = 1.0;
    std::vector<int> PhotonInJet;
    
    double PhotonMinPt = 20;
    double PhotonMaxPt = 500;
    double PhotonMinEta = 2.0;
    double PhotonMaxEta = 3.0;

    for(unsigned pid = 0; pid < photon_pt.size(); ++pid) {
        
        if(photon_jetIdx.at(pid) == jid) { // photon in jet
            
            if(photon_pt.at(pid) >= PhotonMinPt and std::abs(photon_eta.at(pid)) <= PhotonMaxEta and std::abs(photon_eta.at(pid)) >= PhotonMinEta) {

                double phopf_temp = 1.0 - prefireCorr_getPrefireProbability(1, photon_eta.at(pid), photon_pt.at(pid), PhotonMaxPt);
                double elepf_temp = 1.0;
                if(photon_electronIdx.at(pid) > -1) { // what if the electron corresponding to the photon would return a different value?

                    if(electron_pt.at(photon_electronIdx.at(pid)) >= PhotonMinPt and std::abs(electron_eta.at(photon_electronIdx.at(pid))) <= PhotonMaxEta and std::abs(electron_eta.at(photon_electronIdx.at(pid))) >= PhotonMinEta) {
                        
                        elepf_temp = 1.0 - prefireCorr_getPrefireProbability(1, electron_eta.at(photon_electronIdx.at(pid)), electron_pt.at(photon_electronIdx.at(pid)), PhotonMaxPt);
                    }
                }
                
                phopf *= std::min(phopf_temp, elepf_temp); // the highest prefire-probablity between the photon and corresponding electron is chosen
                PhotonInJet.push_back(pid);    
            }
        }
    }
        
    // loop over electrons
    for(unsigned eid = 0; eid < electron_pt.size(); ++eid) {

        // electron should not be in PhotonInJet collection
        if(electron_jetIdx.at(eid) == jid and std::find(PhotonInJet.begin(), PhotonInJet.end(), electron_photonIdx.at(eid)) == PhotonInJet.end()) {

            if(electron_pt.at(eid) >= PhotonMinPt and std::abs(electron_eta.at(eid)) <= PhotonMaxEta and std::abs(electron_eta.at(eid)) >= PhotonMinEta) {
                    
                phopf *= 1.0 - prefireCorr_getPrefireProbability(1, electron_eta.at(eid), electron_pt.at(eid), PhotonMaxPt);
            }
        }
    }

    return phopf;
}


double prefireCorr(Vec_f jet_pt, Vec_f jet_eta, Vec_i jet_muef, Vec_f photon_pt, Vec_f photon_eta, Vec_i photon_jetIdx, Vec_i photon_electronIdx, Vec_f electron_pt, Vec_f electron_eta, Vec_i electron_jetIdx, Vec_i electron_photonIdx) {
    
    // jet_muef = muon energy fraction
    
    double prefw = 1.0;
    
    double JetMinPt = 20;
    double JetMaxPt = 500;
    double JetMinEta = 2.0;
    double JetMaxEta = 3.0;
    
    // loop over jets
    for(unsigned jid = 0; jid < jet_pt.size(); ++jid) {

        double jetpf = 1.0;

        double jetpt = jet_pt.at(jid);
        double jeteta = jet_eta.at(jid);
        double jetmuef = jet_muef.at(jid);

        
        if(jetpt >= JetMinPt and std::abs(jeteta) <= JetMaxEta and std::abs(jeteta) >= JetMinEta and jetmuef < 0.5) {
            
            jetpf *= 1.0 - prefireCorr_getPrefireProbability(0, jeteta, jetpt, JetMaxPt);
        }
        
        double phopf = prefireCorr_eg(jid, photon_pt, photon_eta, photon_jetIdx, photon_electronIdx, electron_pt, electron_eta, electron_jetIdx, electron_photonIdx);
        prefw *= std::min(jetpf, phopf); // the highest prefire-probablity between the jet and the lower-pt photon(s)/elecron(s) from the jet is chosen
    }
    
    // loop over photons/electrons not associated to jets
    prefw *= prefireCorr_eg(-1, photon_pt, photon_eta, photon_jetIdx, photon_electronIdx, electron_pt, electron_eta, electron_jetIdx, electron_photonIdx);
    
    return prefw;
}



double deltaR_(double eta1, double phi1, double eta2, double phi2) {
    
    const double pi = 3.14159265358979;
    double dphi = fabs(phi1-phi2);
    while (dphi>pi) dphi = fabs(dphi - 2.0*pi);
    double deta = eta1-eta2;
  
    return std::sqrt(dphi*dphi + deta*deta);
}




double prefireCorr_MITEWK(Vec_f jet_pt, Vec_f jet_eta, Vec_f jet_phi, Vec_i jet_muef, Vec_f photon_pt, Vec_f photon_eta, Vec_f photon_phi) {

    // Mit-EWK implementation
    // https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/Utils/PrefiringEfficiency.cc

    double prefw = 1.0;
 
    double JetMinPt = 20;
    double JetMaxPt = 500;
    double JetMinEta = 2.0;
    double JetMaxEta = 3.0; 
    double PhotonMinPt = 20;
    double PhotonMaxPt = 500;
    double PhotonMinEta = 2.0;
    double PhotonMaxEta = 3.0;

    unsigned njets = jet_pt.size();

    if(njets == 0) { // photons only
     
        for(unsigned pid = 0; pid < photon_pt.size(); ++pid) {
            
            if(photon_pt.at(pid) >= PhotonMinPt and std::abs(photon_eta.at(pid)) <= PhotonMaxEta and std::abs(photon_eta.at(pid)) >= PhotonMinEta) {
     
                double phopf = prefireCorr_getPrefireProbability(1, photon_eta.at(pid), photon_pt.at(pid), PhotonMaxPt);
                prefw *= 1.0 - phopf;
            }
        }
        
        return prefw;
    }

    for(unsigned jid = 0; jid < jet_pt.size(); ++jid) {

        double jetpf = 0.;
        double phopf = 0.;
        double effObj = 0.;

        double jetpt = jet_pt.at(jid);
        double jeteta = jet_eta.at(jid);
        double jetphi = jet_phi.at(jid);
        double jetmuef = jet_muef.at(jid);
        
        if(jetpt >= JetMinPt and std::abs(jeteta) <= JetMaxEta and std::abs(jeteta) >= JetMinEta and jetmuef < 0.5) {
            
            jetpf = prefireCorr_getPrefireProbability(0, jeteta, jetpt, JetMaxPt);
            effObj = jetpf;
        }
        
        for(unsigned pid = 0; pid < photon_pt.size(); ++pid) {
        
            if(photon_pt.at(pid) >= PhotonMinPt and std::abs(photon_eta.at(pid)) <= PhotonMaxEta and std::abs(photon_eta.at(pid)) >= PhotonMinEta) {
                
                // check if photon in jet
                if(deltaR_(jeteta, jetphi, photon_eta.at(pid), photon_phi.at(pid)) < 0.4)  {
                    
                    phopf = prefireCorr_getPrefireProbability(1, photon_eta.at(pid), photon_pt.at(pid), PhotonMaxPt);
                    if(phopf > jetpf) effObj = phopf; // select maximum
                }
            }
        }
        prefw *= 1.0 - effObj;
    }

    return prefw;
   
}
    


double prefireCorr_PAT(Vec_f jet_pt, Vec_f jet_eta, Vec_f jet_phi, Vec_i jet_muef, Vec_f photon_pt, Vec_f photon_eta, Vec_f photon_phi, Vec_f tag_pt, Vec_f tag_eta, Vec_f tag_phi) {

    // PAT implementation
    // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/plugins/L1PrefiringWeightProducer.cc
    
    // tag_[pt,eta,phi] the object to be filtered from the jets

    
    double nonPrefiringProba = 1.0; // final one
    double nonPrefiringProbaECAL = 1.0;
    double nonPrefiringProbaMuon = 1.0;
 
    double JetMinPt = 20;
    double JetMaxPt = 500;
    double JetMinEta = 2.0;
    double JetMaxEta = 3.0; 
    double PhotonMinPt = 20;
    double PhotonMaxPt = 500;
    double PhotonMinEta = 2.0;
    double PhotonMaxEta = 3.0;

    unsigned njets = jet_pt.size();


    // loop over photons
    for(unsigned pid = 0; pid < photon_pt.size(); ++pid) {
        
        double p_pt = photon_pt.at(pid);
        double p_eta = photon_eta.at(pid);
        
        if(p_pt >= PhotonMinPt and std::abs(p_eta) <= PhotonMaxEta and std::abs(p_eta) >= PhotonMinEta) {
            
            double prefiringprob_gam = prefireCorr_getPrefireProbability(1, photon_eta.at(pid), photon_pt.at(pid), PhotonMaxPt);
            nonPrefiringProbaECAL *= (1.0 - prefiringprob_gam);
        }
    }

    // loop over jets
    for(unsigned jid = 0; jid < jet_pt.size(); ++jid) {

        double j_pt = jet_pt.at(jid);
        double j_eta = jet_eta.at(jid);
        double j_phi = jet_phi.at(jid);
        double j_muef = jet_muef.at(jid);
        
        // check if jet is in the tag collection
        bool skipJet = false;
        for(unsigned tid = 0; tid < tag_pt.size(); ++tid) {
            
            double t_eta = tag_eta.at(tid);
            double t_phi = tag_phi.at(tid);
            if(deltaR_(j_eta, j_phi, t_eta, t_phi) < 0.16) {
                skipJet = true;
                break;
            }
        }
        if(skipJet) continue;
        
        
        
        double nonprefiringprobfromoverlappingphotons = 1.;
        bool foundOverlappingPhotons = false;
        
        if(j_pt < JetMinPt or std::abs(j_eta) > JetMaxEta or std::abs(j_eta) < JetMinEta or j_muef > 0.5) continue;
            
        // loop over photons to remove overlap
        for(unsigned pid = 0; pid < photon_pt.size(); ++pid) {
            
            double p_pt = photon_pt.at(pid);
            double p_eta = photon_eta.at(pid);
            double p_phi = photon_phi.at(pid);
                
            if(p_pt < PhotonMinPt and std::abs(p_eta) > PhotonMaxEta and std::abs(p_eta) < PhotonMinEta) continue;
                    
            // check if photon in jet
            if(deltaR_(j_eta, j_phi, p_eta, p_phi) > 0.16) continue; // overlap criteria
            
            double prefiringprob_gam = prefireCorr_getPrefireProbability(1, p_eta, p_pt, PhotonMaxPt);
            nonprefiringprobfromoverlappingphotons *= (1.0 - prefiringprob_gam);
            foundOverlappingPhotons = true;

        }
            
        double nonprefiringprobfromoverlappingjet = 1.0 - prefireCorr_getPrefireProbability(0, j_eta, j_pt, JetMaxPt);

        
        if(!foundOverlappingPhotons) {
            
            nonPrefiringProbaECAL *= nonprefiringprobfromoverlappingjet;
        }
        else if(nonprefiringprobfromoverlappingphotons > nonprefiringprobfromoverlappingjet) {
          
            // if overlapping photons have a non prefiring rate larger than the jet, then replace these weights by the jet one
            // i.e. select the maximum prefiring probability (or minimum non-prefireing probability)
            if (nonprefiringprobfromoverlappingphotons > 0.) {
                
                nonPrefiringProbaECAL *= nonprefiringprobfromoverlappingjet / nonprefiringprobfromoverlappingphotons;
            } 
            else {
            
                nonPrefiringProbaECAL = 0.;
            }
        }
        else {
            
            // if overlapping photons have a non prefiring rate smaller than the jet, don't consider the jet in the event weight, and do nothing.
        }
    }
    
    nonPrefiringProba = nonPrefiringProbaECAL;

    return nonPrefiringProba;
}
    





double applyElectronSF(Vec_f pt, Vec_f eta, Vec_i q, Vec_b trg) {

    unsigned int size = pt.size();
    double corr = 1.0;
    
    double trgSF_DATA = 1.0;
    double trgSF_MC = 1.0;
  
    
    bool trg_applied = false;
    for (unsigned int i = 0; i < size; ++i) {
        
        double eta_ = eta.at(i);
        double pt_ = pt.at(i);
        int q_ = q.at(i);

        int etabin_ID = std::max(1, std::min(lepSF_el_ID->GetNbinsX(), lepSF_el_ID->GetXaxis()->FindBin(eta_)));
        int ptbin_ID = std::max(1, std::min(lepSF_el_ID->GetNbinsY(), lepSF_el_ID->GetYaxis()->FindBin(pt_)));
        
        int etabin_HLT = std::max(1, std::min(lepSF_el_HLT_pos_DATA->GetNbinsX(), lepSF_el_HLT_pos_DATA->GetXaxis()->FindBin(eta_)));
        int ptbin_HLT = std::max(1, std::min(lepSF_el_HLT_pos_DATA->GetNbinsY(), lepSF_el_HLT_pos_DATA->GetYaxis()->FindBin(pt_)));
                
        corr *= lepSF_el_ID->GetBinContent(etabin_ID, ptbin_ID);
        
        if(q_ > 0) {
            trgSF_DATA *= 1. - lepSF_el_HLT_pos_DATA->GetBinContent(etabin_HLT, ptbin_HLT);
            trgSF_MC *= 1. - lepSF_el_HLT_pos_MC->GetBinContent(etabin_HLT, ptbin_HLT);
        }
        else {
                
            trgSF_DATA *= 1. - lepSF_el_HLT_neg_DATA->GetBinContent(etabin_HLT, ptbin_HLT);
            trgSF_MC *= 1. - lepSF_el_HLT_neg_MC->GetBinContent(etabin_HLT, ptbin_HLT);
        }
        
        /*
        if(trg.at(i) and not trg_applied) {
            if(q_ > 0) corr *= lepSF_el_HLT_pos->GetBinContent(etabin_HLT, ptbin_HLT);
            else corr *= lepSF_el_HLT_neg->GetBinContent(etabin_HLT, ptbin_HLT);
            trg_applied = true;
        }
        */
 
        
        
    }
    corr *= (1. - trgSF_DATA) / (1. - trgSF_MC);
    
    //cout << "corr=" << corr << endl;
    return corr;

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
  
  
    double trgSF_DATA = 1.0;
    double trgSF_MC = 1.0;
    
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
        
        /*
        if(trg.at(i) and not trg_applied) {
            if(q_ > 0) corr *= lepSF_mu_HLT_pos->GetBinContent(etabin_HLT, ptbin_HLT);
            else corr *= lepSF_mu_HLT_neg->GetBinContent(etabin_HLT, ptbin_HLT);
            trg_applied = true;
        }
        */
        
       
        
        if(q_ > 0) {
            trgSF_DATA *= 1. - lepSF_mu_HLT_pos_DATA->GetBinContent(etabin_HLT, ptbin_HLT);
            trgSF_MC *= 1. - lepSF_mu_HLT_pos_MC->GetBinContent(etabin_HLT, ptbin_HLT);
        }
        else {
                
            trgSF_DATA *= 1. - lepSF_mu_HLT_neg_DATA->GetBinContent(etabin_HLT, ptbin_HLT);
            trgSF_MC *= 1. - lepSF_mu_HLT_neg_MC->GetBinContent(etabin_HLT, ptbin_HLT);
        }
        
            
        
        
    }
    corr *= (1. - trgSF_DATA) / (1. - trgSF_MC);
    
    //cout << "corr=" << corr << endl;
    return corr;

}


Vec_f applyRochesterMC(Vec_f pt, Vec_f eta, Vec_f phi, Vec_f ch, Vec_i gen_idx, Vec_f gen_pt, Vec_i nTrackerLayers) {
    
    unsigned int size = pt.size();
    Vec_f res(size, 0.0);
  
    // https://gitlab.cern.ch/akhukhun/roccor
    for (unsigned int i = 0; i < size; ++i) {
            
        if(gen_idx.at(i) >= 0) res[i] = 1.0000*pt[i] * rochester->kSpreadMC(ch[i], pt[i], eta[i], phi[i], gen_pt.at(gen_idx.at(i)), 0, 0);
        else res[i] = 1.0000*pt[i] * rochester->kSmearMC(ch[i], pt[i], eta[i], phi[i], nTrackerLayers.at(i), gRandom->Rndm(), 0, 0);
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

Vec_b goodElectronTriggerCandidateLowPU(const Vec_i& TrigObj_id, const Vec_f& TrigObj_pt, const Vec_f& TrigObj_l1pt, const Vec_f& TrigObj_l2pt, const Vec_i& TrigObj_filterBits) {

   Vec_b res(TrigObj_id.size(),false); // initialize to 0   
   for (unsigned int i = 0; i < res.size(); ++i) {
       if (TrigObj_id[i]  != 11 ) continue;
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


/*
Vec_b passEleMediumID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho) { 

    // from https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/Utils/LeptonIDCuts.hh
    // Tight Electron ID from 2017

    const Double_t ECAL_GAP_LOW  = 1.4442;
    const Double_t ECAL_GAP_HIGH = 1.566;

    // keeping these for now, can always remove at later step
    // if((fabs(electron->eta)>ECAL_GAP_LOW) && (fabs(electron->eta)<ECAL_GAP_HIGH)) return kFALSE;

    if(electron->isConv)            return kFALSE; 
    
    Electron_convVeto
    Electron_eta
    Electron_pt
    Electron_lostHits = missingHits
    Electron_pfIso03_all 
    
    Electron_dxy
    Electron_dz

    Double_t ea = getEffAreaEl(tag.Eta());
    Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

    if(fabs(electron->eta)<=ECAL_GAP_LOW) {
    // // // if(iso >= 0.0354*(tag.Pt()))                                      return kFALSE; // regular ISO
    if(iso/(tag.Pt()) >= 0.0478+0.506/(tag.Pt()))                                      return kFALSE; // original
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie >= 0.0106)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.0547)                              return kFALSE; // original
    if(fabs(electron->dEtaIn) >= 0.0032)                              return kFALSE; // original
    if(electron->hovere >= 0.046+1.16/electron->ecalEnergy+0.0324*rho/electron->ecalEnergy)          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.184*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.05)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.10)                                  return kFALSE;
  } else {
    // // // if(iso >= 0.0646*(tag.Pt()))                                      return kFALSE; // regular ISO
    if(iso/(tag.Pt()) >= 0.0658+0.963/(tag.Pt()))                                      return kFALSE; // original
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie            >= 0.0387)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.0394)                           return kFALSE; // original
    if(fabs(electron->dEtaIn)     >= 0.00632)                         return kFALSE; // original
    if(electron->hovere           >= 0.0275+2.52/electron->ecalEnergy+0.183*rho/electron->ecalEnergy)         return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0721*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.10)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.20)                                 return kFALSE;
  }

  return kTRUE;

}
*/


Vec_b selectElectronEta(Vec_f eta, Vec_f deltaEtaSC) {
    
    Vec_b res(eta.size(), false);
    for(unsigned int i = 0; i < res.size(); ++i) {
        
        double e = std::abs(eta.at(i) + deltaEtaSC.at(i));
        if(e > 2.4 or (e < 1.566 and e > 1.4442)) res[i] = false;
        else res[i] = true;
    }
  
    return res;
}


float lowPU_pu_reweight_sto[100] = {1, 0.041879317950684, 0.2906403692601066, 2.926955277913856, 0.604837179882255, 0.13291723160425126, 0.06483328827710513, 0.03399077851969202, 0.01895726346472344, 0.010933516524116909, 0.006426849386948908, 0.0033939880542484023, 0.00174393639715267, 0.0004962683949237972, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
float lowPU_pu_reweight(int nTrueInt) { if (nTrueInt<100) return lowPU_pu_reweight_sto[nTrueInt]; else return 0; }

#endif
