
import sys,array,math,os,copy
import datasets as ds
import plotter

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def Rebin(h, newbins):

    ''' Rebin with un-even bins '''
    mybins = array.array('d', newbins)
    h1 = h.Rebin(len(mybins)-1, h.GetName(), mybins)
    return h1



def parseHists(hName, leg, rebin=1):

    print("Probe %s" % hName)
    # get data
    h_data = None
    for i,d in enumerate(data):
		
        print("Get %s" % ("%s_%s" % (hName, d)))
        h = fIn.Get("%s_%s" % (hName, d))
        if isinstance(rebin, int): h.Rebin(rebin)
        else: h = Rebin(h, rebin)
        if h_data == None: h_data = h
        else: h_data.Add(h)
        
    h_data.SetName("data")
    leg.AddEntry(h_data, labels_data, "PE")
    
    # get Monte Carlo
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg = h_data.Clone("bkg_nominal") # total MC background
    h_bkg.Reset("ACE")
    for i,bkg in enumerate(bkgs):
		
        hist = None
        for proc in bkg:
            
            print("Get %s" % ("%s_%s" % (hName, proc)))
            h = fIn.Get("%s_%s" % (hName, proc))
            if isinstance(rebin, int): h.Rebin(rebin)
            else: h = Rebin(h, rebin)
            h.Scale(getattr(ds, proc)['xsec']*lumi) # /getattr(ds, proc)['sumw']
            
            if hist == None: hist = h
            else: hist.Add(h)
            
            print(bkg, h.Integral())
            

        hist.SetName("bkg%d" % i)
        hist.SetFillColor(colors_bkgs[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
		
        leg.AddEntry(hist, labels_bkgs[i], "F")
        st.Add(hist)
        h_bkg.Add(hist)
       


    
    # style sum of backgrounds
    h_bkg.SetLineColor(ROOT.kBlack)
    h_bkg.SetFillColor(0)
    h_bkg.SetLineWidth(2)

    # style data
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)
    
    # style background error
    h_bkg_err = h_bkg.Clone("bkg_err")
    h_bkg_err.SetFillColor(ROOT.kBlack)
    h_bkg_err.SetMarkerSize(0)
    h_bkg_err.SetLineWidth(0)
    h_bkg_err.SetFillStyle(3004)
    leg.AddEntry(h_bkg_err, "Stat. Unc.", "F")
    
    
    
    # ratios
    
    # compute MC-only errors
    hRatio = h_data.Clone("hRatio")
    for i in range(hRatio .GetNbinsX()+1): hRatio .SetBinError(i, 0) # set data errors to zero
    hRatio.Divide(h_bkg)
    

    hRatio.SetMarkerStyle(20)
    hRatio.SetMarkerSize(0.7)
    hRatio.SetMarkerColor(ROOT.kBlack)
    hRatio.SetLineColor(ROOT.kBlack)
    hRatio.SetLineWidth(1)
    hRatio.SetFillColor(ROOT.kBlack)
    hRatio.SetFillStyle(3004)    
    
    h_bkg_ratio = h_bkg.Clone("h_bkg_ratio") # nominal point, need to remove stat. error
    h_bkg_err_ratio = h_bkg_err.Clone("h_bkg_err")

    for i in range(h_bkg_ratio.GetNbinsX()+1): 
    
        e, e0, n = h_bkg.GetBinError(i), h_bkg_err.GetBinError(i), h_bkg.GetBinContent(i)
        e0 = math.sqrt(n) # shaded area represents the Poisson error on data
        h_bkg_err_ratio.SetBinContent(i,  1 if n > 0 else 0)
        h_bkg_err_ratio.SetBinError(i, e0/n if n > 0 else 0)
    
   
    #for i in range(h_bkg_ratio.GetNbinsX()+1): h_bkg_ratio.SetBinError(i, 0)
    #h_bkg_err_ratio.Divide(h_bkg_ratio)
    h_bkg_err_ratio.SetMarkerColor(ROOT.kBlack)
    h_bkg_err_ratio.SetLineWidth(2)
    h_bkg_err_ratio.SetFillColor(ROOT.kBlack)
    h_bkg_err_ratio.SetFillStyle(3004)

    

    return h_data, st, h_bkg, h_bkg_err, hRatio, h_bkg_err_ratio



def singlePlot(hName, fOut, xMin, xMax, yMin, yMax, xLabel, yLabel, logY=True, rebin=1, legPos=[]):
    
    leg = ROOT.TLegend(.20, 0.85-(len(bkgs)+2)*0.055, .5, .85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.040)
    
    
    h_data, st, h_bkg, h_bkg_err, hRatio, h_bkg_err_ratio = parseHists(hName, leg, rebin)

    cfg['logy'] = logY
    cfg['xmin'] = xMin
    cfg['xmax'] = xMax
    cfg['ymin'] = yMin
    cfg['ymax'] = yMax
    
    cfg['xtitle'] = xLabel
    cfg['ytitle'] = yLabel
    

    plotter.cfg = cfg
    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio()
        
    ## top panel
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    dummyT.Draw("HIST")
        
    st.Draw("HIST SAME")
    h_bkg_err.Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    h_data.Draw("PE SAME")
    leg.Draw("SAME")
    
    plotter.auxRatio()  
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()  
    


    ## bottom panel
    canvas.cd()
    padB.Draw()
    padB.SetFillStyle(0)
    padB.cd()
    dummyB.Draw("HIST")
    dummyL.Draw("SAME")
     
    hRatio.Draw("P SAME") # E2 SAME
    h_bkg_err_ratio.Draw("E2 SAME")

    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("%s/%s.png" % (outDir, fOut))
    canvas.SaveAs("%s/%s.pdf" % (outDir, fOut))
    canvas.Close()

	
if __name__ == "__main__":

    lumi = 0.199269742 # /fb
    tag = "test"
    
    # default cfg
    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 60,
        'xmax'              : 120,
        'ymin'              : 1e1,
        'ymax'              : 1e5, # 3e6
            
        'xtitle'            : "m(#mu,#mu) (GeV)",
        'ytitle'            : "Events",
            
        'topRight'          : "199 pb^{#minus1} (13 TeV)", 
        'topLeft'           : "#bf{CMS} #scale[0.7]{#it{Preliminary}}",

        'ratiofraction'     : 0.3,
        'ytitleR'           : "Data/MC",
        'yminR'             : 0.85,
        'ymaxR'             : 1.15,
    }


    
    
    data = ["data"]
    bkgs = [["DY_MiNNLO"]]
    
    labels_bkgs = ["DY (MiNNLO+N3LL)"]
    labels_data = "Data"

    colors_bkgs = [ROOT.TColor.GetColor(248, 206, 104)]
    
    
    if True:
    
        labels_bkgs = ["DY #rightarrow #mu^{#plus}#mu^{#minus} (MiNNLO+N3LL)"]
    
        outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Zmumu_lepSF_preFire_globalMuon/"
        if not os.path.exists(outDir): os.makedirs(outDir)
        fIn = ROOT.TFile("/eos/user/j/jaeyserm/www/wmass/test/Zmumu.root")
        
        
        singlePlot("m_mumu", "m_mumu", 60, 120, 1e0, 1e6, "m(#mu^{#plus}, #mu^{#minus}) (GeV)", "Events")
        singlePlot("m_mumu", "m_mumu_nolog", 60, 120, 0, 1.5e4, "m(#mu^{#plus}, #mu^{#minus}) (GeV)", "Events", logY=False)
     
        singlePlot("pt_mumu", "pt_mumu", 0, 100, 1e1, 1e6, "p_{T}(#mu^{#plus}, #mu^{#minus}) (GeV)", "Events", rebin=2)
        singlePlot("eta_mumu", "eta_mumu", -2.4, 2.4, 1e1, 1e6, "Rapidity (#mu^{#plus}, #mu^{#minus})", "Events", rebin=2)
        

        singlePlot("MuonPlus_eta", "Muon_plus_eta", -2.4, 2.4, 1e1, 1e6, "Triggered #mu^{#plus} #eta", "Events", rebin=1000)
        singlePlot("MuonMinus_eta", "Muon_minus_eta", -2.4, 2.4, 1e1, 1e6, "Non-triggered #mu^{#minus} #eta", "Events", rebin=1000)
        
        singlePlot("MuonPlus_abseta", "Muon_plus_abseta", 0, 2.4, 1e1, 1e6, "Triggered #mu^{#plus} #eta", "Events", rebin=1000)
        singlePlot("MuonMinus_abseta", "Muon_minus_abseta", 0, 2.4, 1e1, 1e6, "Non-triggered #mu^{#minus} #eta", "Events", rebin=1000)
        
        bins_sf = [-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0.,0.3,0.9,1.2,1.6,2.1,2.4]
        singlePlot("MuonPlus_eta", "MuonPlus_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Triggered #mu^{#plus} #eta", "Events", rebin=bins_sf)
        singlePlot("MuonMinus_eta", "MuonMinus_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Non-triggered #mu^{#minus} #eta", "Events", rebin=bins_sf)
        
        singlePlot("muon_eta", "Muon_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Electrons (#mu^{#plus} + #mu^{#minus}) #eta", "Events", rebin=bins_sf)

        
        singlePlot("MuonPlus_pt", "Muon_plus_pt", 0, 100, 1e1, 1e5, "Triggered #mu^{#plus} p_{T} (GeV)", "Events")
        singlePlot("MuonMinus_pt", "Muon_minus_pt", 0, 100, 1e1, 1e5, "Non-triggered #mu^{#minus} p_{T} (GeV)", "Events")



        
        sys.exit()
        
        
        
        
        
       
        singlePlot("m_mumu", "m_mumu", 60, 120, 1e0, 1e6, "m(#mu,#mu) (GeV)", "Events")
        singlePlot("m_mumu", "m_mumu_nolog", 60, 120, 0, 1e4, "m(#mu,#mu) (GeV)", "Events", logY=False)
     
        singlePlot("pt_mumu", "pt_mumu", 0, 100, 1e1, 1e6, "p_{T}(#mu,#mu) (GeV)", "Events", rebin=2)
        singlePlot("eta_mumu", "eta_mumu", -2.4, 2.4, 1e1, 1e6, "#eta (#mu,#mu)", "Events", rebin=2)
        
        #singlePlot("leading_muon_pt", "leading_muon_pt", 0, 100, 1e1, 1e5, "Leading muon p_{T} (GeV)", "Events", rebin=rebin)
        #singlePlot("subleading_muon_pt", "subleading_muon_pt", 0, 100, 1e1, 1e5, "Subleading muon p_{T} (GeV)", "Events", rebin=rebin)
        
        #singlePlot("leading_muon_eta", "leading_muon_eta", -2.4, 2.4, 1e1, 1e6, "Leading muon #eta", "Events", rebin=rebin)
        #singlePlot("subleading_muon_eta", "subleading_muon_eta", -2.4, 2.4, 1e1, 1e6, "Subleading muon #eta", "Events", rebin=rebin)
        #singlePlot("leading_muon_abs_eta", "leading_muon_abs_eta", 0, 2.4, 1e1, 1e6, "Leading muon |#eta|", "Events", rebin=rebin)
        #singlePlot("subleading_muon_abs_eta", "subleading_muon_abs_eta", 0, 2.4, 1e1, 1e6, "Subleading muon |#eta|", "Events", rebin=rebin)
        
        singlePlot("MuonPlus_eta", "muon_plus_eta", -2.4, 2.4, 1e1, 1e6, "Triggered muon (+) #eta", "Events", rebin=10)
        singlePlot("MuonMinus_eta", "muon_minus_eta", -2.4, 2.4, 1e1, 1e6, "Non-triggered muon (-) #eta", "Events", rebin=10)
        
        singlePlot("MuonPlus_abseta", "muon_plus_abseta", 0, 2.4, 1e1, 1e6, "Triggered muon (+) #eta", "Events", rebin=10)
        singlePlot("MuonMinus_abseta", "muon_minus_abseta", 0, 2.4, 1e1, 1e6, "Non-triggered muon (-) #eta", "Events", rebin=10)
        
        
        singlePlot("MuonPlus_eta", "MuonPlus_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Triggered muon (+) #eta", "Events", rebin=bins_sf)
        singlePlot("MuonMinus_eta", "MuonMinus_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Non-triggered muon (-) #eta", "Events", rebin=bins_sf)
        
        singlePlot("MuonPlus_pt", "muon_plus_pt", 0, 100, 1e1, 1e5, "Triggered muon (+) p_{T} (GeV)", "Events")
        singlePlot("MuonMinus_pt", "muon_minus_pt", 0, 100, 1e1, 1e5, "Non-triggered muon (-) p_{T} (GeV)", "Events")
        
    if False:
        labels_bkgs = ["DY #rightarrow e^{#plus}e^{#minus} (MiNNLO+N3LL)"]
        bkgs = [["DY_MiNNLO_ee"]]
    
        outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Zee_testSF/"
        if not os.path.exists(outDir): os.makedirs(outDir)
        fIn = ROOT.TFile("/eos/user/j/jaeyserm/www/wmass/test/Zee.root")
        
       
        singlePlot("m_ee", "m_ee", 60, 120, 1e0, 1e6, "m(e^{#plus}, e^{#minus}) (GeV)", "Events")
        singlePlot("m_ee", "m_ee_nolog", 60, 120, 0, 1.5e4, "m(e^{#plus}, e^{#minus}) (GeV)", "Events", logY=False)
     
        singlePlot("pt_ee", "pt_ee", 0, 100, 1e1, 1e6, "p_{T}(e^{#plus}, e^{#minus}) (GeV)", "Events", rebin=2)
        singlePlot("eta_ee", "eta_ee", -2.4, 2.4, 1e1, 1e6, "Rapidity (e^{#plus}, e^{#minus})", "Events", rebin=2)
        
        #singlePlot("leading_Electron_pt", "leading_Electron_pt", 0, 100, 1e1, 1e5, "Leading Electron p_{T} (GeV)", "Events", rebin=rebin)
        #singlePlot("subleading_Electron_pt", "subleading_Electron_pt", 0, 100, 1e1, 1e5, "Subleading Electron p_{T} (GeV)", "Events", rebin=rebin)
        
        #singlePlot("leading_Electron_eta", "leading_Electron_eta", -2.4, 2.4, 1e1, 1e6, "Leading Electron #eta", "Events", rebin=rebin)
        #singlePlot("subleading_Electron_eta", "subleading_Electron_eta", -2.4, 2.4, 1e1, 1e6, "Subleading Electron #eta", "Events", rebin=rebin)
        #singlePlot("leading_Electron_abs_eta", "leading_Electron_abs_eta", 0, 2.4, 1e1, 1e6, "Leading Electron |#eta|", "Events", rebin=rebin)
        #singlePlot("subleading_Electron_abs_eta", "subleading_Electron_abs_eta", 0, 2.4, 1e1, 1e6, "Subleading Electron |#eta|", "Events", rebin=rebin)
        
        singlePlot("ElectronPlus_eta", "Electron_plus_eta", -2.4, 2.4, 1e1, 1e6, "Triggered e^{#plus} #eta", "Events", rebin=1000)
        singlePlot("ElectronMinus_eta", "Electron_minus_eta", -2.4, 2.4, 1e1, 1e6, "Non-triggered e^{#minus} #eta", "Events", rebin=1000)
        
        singlePlot("ElectronPlus_abseta", "Electron_plus_abseta", 0, 2.4, 1e1, 1e6, "Triggered e^{#plus} #eta", "Events", rebin=1000)
        singlePlot("ElectronMinus_abseta", "Electron_minus_abseta", 0, 2.4, 1e1, 1e6, "Non-triggered e^{#minus} #eta", "Events", rebin=1000)
        
        bins_sf = [-2.4, -2.25, -2.0, -1.566, -1.4442, -1.0, -0.5, 0, 0.5, 1.0, 1.4442, 1.566, 2.0, 2.25, 2.4]
        singlePlot("ElectronPlus_eta", "ElectronPlus_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Triggered e^{#plus} #eta", "Events", rebin=bins_sf)
        singlePlot("ElectronMinus_eta", "ElectronMinus_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Non-triggered e^{#minus} #eta", "Events", rebin=bins_sf)
        
        singlePlot("electron_eta", "Electron_eta_SFBinning", -2.4, 2.4, 1e1, 1e6, "Electrons (e^{#plus} + e^{#minus}) #eta", "Events", rebin=bins_sf)

        
        singlePlot("ElectronPlus_pt", "Electron_plus_pt", 0, 100, 1e1, 1e5, "Triggered e^{#plus} p_{T} (GeV)", "Events")
        singlePlot("ElectronMinus_pt", "Electron_minus_pt", 0, 100, 1e1, 1e5, "Non-triggered e^{#minus} p_{T} (GeV)", "Events")
        