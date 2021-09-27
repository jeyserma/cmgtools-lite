
import sys, glob, os, copy
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)




if __name__ == "__main__":

    outDir = "/eos/user/j/jaeyserm/www/wmass/LeptonSF"

    toPlot = ["mu_SIT", "mu_STA", "mu_HLT_pos", "mu_HLT_neg"]    
    titles = ["Muon Selection+ID+Isolation", "Muon Standalone", "Muon HLT q=#plus1", "Muon HLT q=#minus1"]
    
    coll = []
    
    for i,p in enumerate(toPlot):
    
        fIn = ROOT.TFile("%s_DATA.root" % p)
        h_data = copy.deepcopy(fIn.Get("hEffEtaPt"))
        fIn.Close()
        
        
     
        fIn = ROOT.TFile("%s_MC.root" % p)
        h_mc = copy.deepcopy(fIn.Get("hEffEtaPt"))
        fIn.Close()
        
        print(h_data.GetNbinsX(), h_data.GetNbinsY())
        print(h_mc.GetNbinsX(), h_mc.GetNbinsY())
        
        '''
        if "HLT" in toPlot:
        
            for a in range(0, h_data.GetNbinsX())
                for b in range(0, h_data.GetNbinsY()):
                    
                    eff_data = 1. - h_data.GetBinContent(a, b)
                    eff_mc = 1. - h_mc.GetBinContent(a, b)
                    h_data.SetBinContent(a, b)
            
        else: h_data.Divide(h_mc)
        '''
        h_data.Divide(h_mc)
        h_data.SetName(p)
        h_data.SetTitle(titles[i])
        coll.append(h_data)
        
        c = ROOT.TCanvas("c", "c", 1000, 1000)
        c.SetTopMargin(0.07)
        c.SetRightMargin(0.2)
        c.SetLeftMargin(0.15)
        c.SetBottomMargin(0.14)
        c.SetLogy()
        h_data.GetYaxis().SetRangeUser(25, 1e4)
        h_data.GetXaxis().SetTitle("Muon #eta")
        h_data.GetYaxis().SetTitle("Muon p_{T} (GeV)")
        h_data.Draw("COLZ")
        c.SaveAs("%s/%s.png" % (outDir, p))
        c.SaveAs("%s/%s.pdf" % (outDir, p))
        
        del c
       
    for i,p in enumerate(["mu_HLT_pos", "mu_HLT_neg"] ):
    
        fIn = ROOT.TFile("%s_DATA.root" % p)
        h_data = copy.deepcopy(fIn.Get("hEffEtaPt"))
        fIn.Close()
        h_data.SetName(p+"_DATA")
        
     
        fIn = ROOT.TFile("%s_MC.root" % p)
        h_mc = copy.deepcopy(fIn.Get("hEffEtaPt"))
        fIn.Close()
        h_mc.SetName(p+"_MC")
        
        
        coll.append(h_data)
        coll.append(h_mc)
       
    fOut = ROOT.TFile("lepton_SF.root", "RECREATE")
    for c in coll: c.Write()
    fOut.Close()