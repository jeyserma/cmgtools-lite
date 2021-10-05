
import sys, glob, os, copy
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)



if __name__ == "__main__":

    # data
    fIn = ROOT.TFile("lowPU/lumiPileup/pilup_data.root")
    fIn.ls()
    h_data = copy.deepcopy(fIn.Get("pileup"))
    h_data.SetName("pu_data")
    fIn.Close()

    
    # mc
    fIn = ROOT.TFile("lowPU/lumiPileup/pilup_mc.root")
    fIn.ls()
    h_mc = copy.deepcopy(fIn.Get("h_pu"))
    h_mc.SetName("pu_mc")
    fIn.Close()
    
    h_data.Scale(1./h_data.Integral())
    h_mc.Scale(1./h_mc.Integral())
    
    print(h_data.GetMean())
    print(h_mc.GetMean())
    
    print(h_data.GetNbinsX())
    print(h_mc.GetNbinsX())
    
    weights = []
    #for i in range(1, h_data.GetNbinsX()): 
    for i in range(1, h_data.GetNbinsX()+1): 
    
        # mc: i=1 correspond to PU=0
        # data:
    
        #print(i-1, h_data.GetBinContent(i-1), h_mc.GetBinContent(i))
        print(i, h_data.GetBinContent(i), h_mc.GetBinContent(i))
        if(h_mc.GetBinContent(i) > 0): weights.append(h_data.GetBinContent(i-1)/h_mc.GetBinContent(i))
        else: weights.append(1)
        
    print(weights)
    
    
    sys.exit()
    

    tf = ROOT.TFile.Open(datafile)
    datahist = tf.Get(dataname)
    datahist.SetDirectory(0)
    tf.Close()

    mchist = None
    tf = ROOT.TFile.Open(mcfile)
    for i,mcn in enumerate(mcname.split(',')):
        tmp = tf.Get(mcn)
        if i == 0:
            mchist = copy.deepcopy(tmp.Clone("mchist"))
            mchist.SetDirectory(0)
        else:
            mchist.Add(tmp)
    tf.Close()

    if not datahist or not mchist:
        print("Error getting histograms")
        quit()

    if datahist.GetNbinsX() != mchist.GetNbinsX():
        print("Warning: histograms have different number of bins")
        quit()
    nbins = datahist.GetNbinsX()
    if datahist.GetXaxis().GetBinLowEdge(1) != mchist.GetXaxis().GetBinLowEdge(1):
        print("Warning: histograms have different low edge")
        quit()
    if datahist.GetXaxis().GetBinLowEdge(1+nbins) != mchist.GetXaxis().GetBinLowEdge(1+nbins):
        print("Warning: histograms have different up edge")
        quit()

    #if not ROOT.TH1.CheckConsistency(datahist,mchist):
    #    print("Warning: histograms not compatible"
    #    quit()

    integralDataFull = datahist.Integral(0,1+datahist.GetNbinsX()) 
    integralMCFull = mchist.Integral(0,1+mchist.GetNbinsX()) 
    integralDataRange = datahist.Integral(0,normalizeIntegralUpToBin) 
    integralMCRange = mchist.Integral(0,normalizeIntegralUpToBin) 
    
    mchist.Scale(datahist.Integral(0,normalizeIntegralUpToBin)/mchist.Integral(0,normalizeIntegralUpToBin))
    integralMCafterNormData = mchist.Integral(0,1+mchist.GetNbinsX())
    ratio = datahist.Clone("puWeight_2016_UL")
    oldratio = datahist.Clone("puWeight_2016_ReReco")
    ratio.Reset("ICESM")
    oldratio.Reset("ICESM")

    for i in range(1,1+nbins):
        puratio = 1.0
        if mchist.GetBinContent(i) == 0:
            ratio.SetBinContent(i,puratio)
        else:
            puratio = datahist.GetBinContent(i)/mchist.GetBinContent(i)
            ratio.SetBinContent(i,puratio)
        newWeights.append(puratio)
        oldratio.SetBinContent(i,oldWeights[i-1])