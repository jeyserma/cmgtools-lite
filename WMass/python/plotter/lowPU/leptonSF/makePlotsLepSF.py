
import sys, glob, os, copy
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)




if __name__ == "__main__":

    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/LeptonSF_MiNNLO"

    toPlot = ["mu_SF2D_nominal_idip_lowPU_both", "mu_SF2D_nominal_iso_lowPU_both", "mu_SF2D_nominal_trigger_lowPU_both", "el_SF2D_nominal_trigger_lowPU_both", "el_SF2D_nominal_idiso_lowPU_both"]    
    titles = ["Muon Selection+ID+Isolation", "Muon Standalone", "Muon HLT q=#plus1", "Muon HLT q=#minus1", "Electron GSF+ID+Isolation", "Electron HLT q=#plus1", "Electron HLT q=#minus1"]
    
    

    fIn = ROOT.TFile("2021-10-22_allSFs_nodz_dxybs_lowPU.root")
    

    for i,p in enumerate(toPlot):
    
        h = fIn.Get(p)
        
        zMin, zMax = 1e9, -1e9
        
        for ix in range(0, h.GetNbinsX()+1):
            for iy in range(0, h.GetNbinsY()+1):
            
                binc = h.GetBinContent(ix, iy)
                if binc < 0.7 or binc > 1.3: continue
                #print(binc)
                
                if binc < zMin: zMin = binc
                if binc > zMax: zMax = binc
   
        
        print(zMin, zMax)
        
        c = ROOT.TCanvas("c", "c", 1000, 1000)
        c.SetTopMargin(0.07)
        c.SetRightMargin(0.2)
        c.SetLeftMargin(0.16)
        c.SetBottomMargin(0.15)
        
        h.GetXaxis().SetRangeUser(-2.4, 2.4)
        h.GetYaxis().SetRangeUser(25, 65)
        h.GetZaxis().SetRangeUser(0.99*zMin, 1.01*zMax)
        h.GetXaxis().SetTitle("Muon #eta")
        h.GetYaxis().SetTitle("Muon p_{T} (GeV)")
        if "el_" in p:
            h.GetXaxis().SetTitle("Electron #eta")
            h.GetYaxis().SetTitle("Electron p_{T} (GeV)")
        h.Draw("COLZ")
        c.SaveAs("%s/%s.png" % (outDir, p))
        c.SaveAs("%s/%s.pdf" % (outDir, p))
        
        del c
       
   