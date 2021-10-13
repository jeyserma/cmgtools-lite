
import sys,array,math,os,copy
import datasets as ds
import plotter

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

#sys.path.insert(0, '/afs/cern.ch/work/j/jaeyserm/pythonlibs')




'''
muon_pt             : Muon_pt[goodMuonsCharge][0]       : 100,0,100; XTitle="Muon p_{T} (GeV)", Legend='TR', IncludeOverflows=True
muon_eta            : Muon_eta[goodMuonsCharge][0]      : 25,-2.4,2.4; XTitle="Muon #eta (GeV)", Legend='TR', IncludeOverflows=True
muon_phi            : Muon_phi[goodMuonsCharge][0]      : 25,-4,4; XTitle="Muon #phi (GeV)", Legend='TR', IncludeOverflows=True
m_mumu              : m_mumu                            : 60,60,120; XTitle="m(#mu #mu) (GeV)", Legend='TR', IncludeOverflows=True
pt_mumu             : pt_mumu                           : 100,0,100; XTitle="p_{T}(#mu #mu) (GeV)", Legend='TR', IncludeOverflows=True
    
muon_leading_pt     : leadingMuonPt                     : 100,0,100; XTitle="Leading muon p_{T} (GeV)", Legend='TR', IncludeOverflows=True
muon_subleading_pt  : subLeadingMuonPt                  : 100,0,100; XTitle="Subleading muon p_{T} (GeV)", Legend='TR', IncludeOverflows=True


'''





def getHist(proc, sel, hName, rebin):


    fIn = ROOT.TFile(dc.datasets[proc]['file'].format(sel=sel))
    h = copy.deepcopy(fIn.Get(hName))
    h.Rebin(rebin)
    h.Scale(lumi*dc.datasets[proc]['xsec']*1e6/dc.datasets[proc]['nevents'])
    fIn.Close()
    return h



def recoil_log(sel):

    fOut = "leptonic_recoil_m"
    xMin, xMax = 0, 150
    rebin = 500

    hName = "leptonic_recoil_m"

    st = ROOT.THStack()
    st.SetName("stack")
        
    leg = ROOT.TLegend(.4, 0.97-(len(bkgs)+2)*0.055, .7, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    


    h_sig = None
    for sig in sigs:
    
        h = getHist(sig, sel, hName, rebin)
        if h_sig == None: h_sig = h
        else: h_sig.Add(h)
		
    h_sig.SetLineColor(ROOT.TColor.GetColor("#BF2229"))
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    leg.AddEntry(h_sig, sigLegend, "L") #  (10#times)

    
    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg = None
    for i,bkg in enumerate(bkgs):
		
        hist = None
        for x in bgks_cfg[bkg]:
            
            h = getHist(x, sel, hName, rebin)
		
            if hist == None: hist = h
            else: hist.Add(h)
		
        hist.SetName(bkg)
        hist.SetFillColor(bkgs_colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
		
        leg.AddEntry(hist, bkgs_legends[i], "F")
        st.Add(hist)
        if h_bkg == None:
            h_bkg = copy.deepcopy(hist)
            h_bkg.SetName("h_bkg")
        else: h_bkg.Add(hist)
        

    
    yMax = math.ceil(h_bkg.GetMaximum()*100)/10.
    

    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : 1e2,
        'ymax'              : 1e5,
            
        'xtitle'            : "m_{rec} (GeV)",
        'ytitle'            : "Events / 0.5 GeV",
            
        'topRight'          : "#sqrt{s} = 240 GeV, 5 ab^{#minus1}",
        'topLeft'           : "#bf{FCCee} #scale[0.7]{#it{Simulation}}",

    }
        

    plotter.cfg = cfg
    canvas = plotter.canvas()
        
    dummy = plotter.dummy()
    dummy.Draw("HIST")
        
    st.Draw("HIST SAME")
    
    
    '''
    hTot_err = hTot.Clone("hTot_err")
    hTot_err.SetFillColor(ROOT.kBlack)
    hTot_err.SetMarkerColor(ROOT.kBlack)
    hTot_err.SetFillStyle(3004)
    leg.AddEntry(hTot_err, "Stat. Unc.", "F")
    '''
    
    h_bkg.SetLineColor(ROOT.kBlack)
    h_bkg.SetFillColor(0)
    h_bkg.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    
    h_sig.Draw("HIST SAME")
    
    leg.Draw("SAME")
        
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/recoil_%s.png" % sel)
    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/recoil_%s.pdf" % sel)
    canvas.Close()



def recoil_lin(sel):

    fOut = "leptonic_recoil_m"
    xMin, xMax = 0, 150
    rebin = 500

    hName = "leptonic_recoil_m"

    st = ROOT.THStack()
    st.SetName("stack")
        
    leg = ROOT.TLegend(.4, 0.97-(len(bkgs)+2)*0.055, .7, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    


    h_sig = None
    for sig in sigs:
    
        h = getHist(sig, sel, hName, rebin)
        if h_sig == None: h_sig = h
        else: h_sig.Add(h)
		
    h_sig.SetLineColor(ROOT.TColor.GetColor("#BF2229"))
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    leg.AddEntry(h_sig, sigLegend, "L") #  (10#times)

    
    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg = None
    for i,bkg in enumerate(bkgs):
		
        hist = None
        for x in bgks_cfg[bkg]:
            
            h = getHist(x, sel, hName, rebin)
		
            if hist == None: hist = h
            else: hist.Add(h)
		
        hist.SetName(bkg)
        hist.SetFillColor(bkgs_colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
		
        leg.AddEntry(hist, bkgs_legends[i], "F")
        st.Add(hist)
        if h_bkg == None:
            h_bkg = copy.deepcopy(hist)
            h_bkg.SetName("h_bkg")
        else: h_bkg.Add(hist)
        

    
    yMax = math.ceil(h_bkg.GetMaximum()*100)/10.
    

    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : 0,
        'ymax'              : 1.5e4,
            
        'xtitle'            : "m_{rec} (GeV)",
        'ytitle'            : "Events / 0.5 GeV",
            
        'topRight'          : "#sqrt{s} = 240 GeV, 5 ab^{#minus1}",
        'topLeft'           : "#bf{FCCee} #scale[0.7]{#it{Simulation}}",

    }
        

    plotter.cfg = cfg
    canvas = plotter.canvas()
        
    dummy = plotter.dummy()
    dummy.Draw("HIST")
        
    st.Draw("HIST SAME")
    
    
    '''
    hTot_err = hTot.Clone("hTot_err")
    hTot_err.SetFillColor(ROOT.kBlack)
    hTot_err.SetMarkerColor(ROOT.kBlack)
    hTot_err.SetFillStyle(3004)
    leg.AddEntry(hTot_err, "Stat. Unc.", "F")
    '''
    
    h_bkg.SetLineColor(ROOT.kBlack)
    h_bkg.SetFillColor(0)
    h_bkg.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    
    h_sig.Draw("HIST SAME")
    
    leg.Draw("SAME")
        
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/recoil_lin_%s.png" % sel)
    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/recoil_lin_%s.pdf" % sel)
    canvas.Close()



def recoil_lin_zoom(sel):

    fOut = "leptonic_recoil_m"
    xMin, xMax = 120, 140
    rebin = 100

    hName = "leptonic_recoil_m"

    st = ROOT.THStack()
    st.SetName("stack")
        
    leg = ROOT.TLegend(.4, 0.97-(len(bkgs)+2)*0.055, .7, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    


    h_sig = None
    for sig in sigs:
    
        h = getHist(sig, sel, hName, rebin)
        if h_sig == None: h_sig = h
        else: h_sig.Add(h)
		
    h_sig.SetLineColor(ROOT.TColor.GetColor("#BF2229"))
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    leg.AddEntry(h_sig, sigLegend, "L") #  (10#times)

    
    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg = None
    for i,bkg in enumerate(bkgs):
		
        hist = None
        for x in bgks_cfg[bkg]:
            
            h = getHist(x, sel, hName, rebin)
		
            if hist == None: hist = h
            else: hist.Add(h)
		
        hist.SetName(bkg)
        hist.SetFillColor(bkgs_colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
		
        leg.AddEntry(hist, bkgs_legends[i], "F")
        st.Add(hist)
        if h_bkg == None:
            h_bkg = copy.deepcopy(hist)
            h_bkg.SetName("h_bkg")
        else: h_bkg.Add(hist)
        

    
    yMax = math.ceil(h_bkg.GetMaximum()*100)/10.
    

    cfg = {

        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : 0,
        'ymax'              : 1000,
            
        'xtitle'            : "m_{rec} (GeV)",
        'ytitle'            : "Events / 0.1 GeV",
            
        'topRight'          : "#sqrt{s} = 240 GeV, 5 ab^{#minus1}",
        'topLeft'           : "#bf{FCCee} #scale[0.7]{#it{Simulation}}",

    }
        

    plotter.cfg = cfg
    canvas = plotter.canvas()
        
    dummy = plotter.dummy()
    dummy.Draw("HIST")
        
    st.Draw("HIST SAME")
    
    
    '''
    hTot_err = hTot.Clone("hTot_err")
    hTot_err.SetFillColor(ROOT.kBlack)
    hTot_err.SetMarkerColor(ROOT.kBlack)
    hTot_err.SetFillStyle(3004)
    leg.AddEntry(hTot_err, "Stat. Unc.", "F")
    '''
    
    h_bkg.SetLineColor(ROOT.kBlack)
    h_bkg.SetFillColor(0)
    h_bkg.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    
    h_sig.Draw("HIST SAME")
    
    leg.Draw("SAME")
        
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/recoil_lin_zoom_%s.png" % sel)
    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/recoil_lin_zoom_%s.pdf" % sel)
    canvas.Close()


def mll(sel):

    fOut = "zed_leptonic_m"
    xMin, xMax = 0, 250
    rebin = 1500 # 3000 = /1GeV, 300 = 0.1 GeV

    hName = "zed_leptonic_m"

    st = ROOT.THStack()
    st.SetName("stack")
        
    leg = ROOT.TLegend(.47, 0.97-(len(bkgs)+2)*0.055, 0.47+0.3, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    


    h_sig = None
    for sig in sigs:
    
        h = getHist(sig, sel, hName, rebin)
        if h_sig == None: h_sig = h
        else: h_sig.Add(h)
		
    h_sig.SetLineColor(ROOT.TColor.GetColor("#BF2229"))
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    #h_sig.Scale(10)
    leg.AddEntry(h_sig, sigLegend, "L") #  (10#times)

    
    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg = None
    for i,bkg in enumerate(bkgs):
		
        hist = None
        for x in bgks_cfg[bkg]:
            
            h = getHist(x, sel, hName, rebin)
		
            if hist == None: hist = h
            else: hist.Add(h)
		
        hist.SetName(bkg)
        hist.SetFillColor(bkgs_colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
		
        leg.AddEntry(hist, bkgs_legends[i], "F")
        st.Add(hist)
        if h_bkg == None:
            h_bkg = copy.deepcopy(hist)
            h_bkg.SetName("h_bkg")
        else: h_bkg.Add(hist)
        

    
    yMax = math.ceil(h_bkg.GetMaximum()*100)/10.
    

    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : 1e3,
        'ymax'              : 1e7, # 3e6
            
        'xtitle'            : "m_{#mu^{+},#mu^{#minus}} (GeV)",
        'ytitle'            : "Events / 0.5 GeV",
            
        'topRight'          : "#sqrt{s} = 240 GeV, 5 ab^{#minus1}", 
        'topLeft'           : "#bf{FCCee} #scale[0.7]{#it{Simulation}}",

    }
        

    plotter.cfg = cfg
    canvas = plotter.canvas()
        
    dummy = plotter.dummy()
    dummy.Draw("HIST")
        
    st.Draw("HIST SAME")
    
    ROOT.TGaxis.SetExponentOffset(-0.07,0.015)
    
    
    '''
    hTot_err = hTot.Clone("hTot_err")
    hTot_err.SetFillColor(ROOT.kBlack)
    hTot_err.SetMarkerColor(ROOT.kBlack)
    hTot_err.SetFillStyle(3004)
    leg.AddEntry(hTot_err, "Stat. Unc.", "F")
    '''
    
    h_bkg.SetLineColor(ROOT.kBlack)
    h_bkg.SetFillColor(0)
    h_bkg.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    
    h_sig.Draw("HIST SAME")
    
    leg.Draw("SAME")
        
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/mll_%s.png" % sel)
    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/mll_%s.pdf" % sel)
    canvas.Close()





def ptll(sel):

    fOut = "zed_leptonic_pt"
    xMin, xMax = 0, 120
    rebin = 1000

    hName = "zed_leptonic_pt"

    st = ROOT.THStack()
    st.SetName("stack")
        
    leg = ROOT.TLegend(.4, 0.97-(len(bkgs)+2)*0.055, .7, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    


    h_sig = None
    for sig in sigs:
    
        h = getHist(sig, sel, hName, rebin)
        if h_sig == None: h_sig = h
        else: h_sig.Add(h)
		
    h_sig.SetLineColor(ROOT.TColor.GetColor("#BF2229"))
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    #h_sig.Scale(10)
    leg.AddEntry(h_sig, sigLegend, "L") #  (10#times)

    
    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg = None
    for i,bkg in enumerate(bkgs):
		
        hist = None
        for x in bgks_cfg[bkg]:
            
            h = getHist(x, sel, hName, rebin)
		
            if hist == None: hist = h
            else: hist.Add(h)
		
        hist.SetName(bkg)
        hist.SetFillColor(bkgs_colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
		
        leg.AddEntry(hist, bkgs_legends[i], "F")
        st.Add(hist)
        if h_bkg == None:
            h_bkg = copy.deepcopy(hist)
            h_bkg.SetName("h_bkg")
        else: h_bkg.Add(hist)
        


    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : 1e2,
        'ymax'              : 1e7, # 3e6
            
        'xtitle'            : "p_{T}^{#mu^{+},#mu^{#minus}} (GeV)",
        'ytitle'            : "Events / 0.5 GeV",
            
        'topRight'          : "#sqrt{s} = 240 GeV, 5 ab^{#minus1}", 
        'topLeft'           : "#bf{FCCee} #scale[0.7]{#it{Simulation}}",

    }
        

    plotter.cfg = cfg
    canvas = plotter.canvas()
        
    dummy = plotter.dummy()
    dummy.Draw("HIST")
        
    st.Draw("HIST SAME")
    
    ROOT.TGaxis.SetExponentOffset(-0.07,0.015)
    
    
    '''
    hTot_err = hTot.Clone("hTot_err")
    hTot_err.SetFillColor(ROOT.kBlack)
    hTot_err.SetMarkerColor(ROOT.kBlack)
    hTot_err.SetFillStyle(3004)
    leg.AddEntry(hTot_err, "Stat. Unc.", "F")
    '''
    
    h_bkg.SetLineColor(ROOT.kBlack)
    h_bkg.SetFillColor(0)
    h_bkg.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg.Draw("HIST SAME")
    
    h_sig.Draw("HIST SAME")
    
    leg.Draw("SAME")
        
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/ptll_%s.png" % sel)
    canvas.SaveAs("/afs/cern.ch/user/j/jaeyserm/www/FCCee/StatAnalysis/FinalPlots/ptll_%s.pdf" % sel)
    canvas.Close()



def Rebin(h, newbins):

    ''' Rebin with un-even bins '''
    mybins = array.array('d', newbins)
    h1 = h.Rebin(len(mybins)-1, h.GetName(), mybins)
    return h1



def parseHists(hName, leg, rebin=1):

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

def mZ_mumu():

    fOut = "mZ_mumu"
    rebin = 1
    hName = "m_mumu"

        
    leg = ROOT.TLegend(.20, 0.85-(len(bkgs)+2)*0.055, .5, .85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.040)
    
    
    h_data, st, h_bkg, h_bkg_err, hRatio, h_bkg_err_ratio = parseHists(hName, leg)

    cfg['logy'] = True
    cfg['xmin'] = 60
    cfg['xmax'] = 120
    cfg['ymin'] = 1e1
    cfg['ymax'] = 1e5
    
    cfg['xtitle'] = "m(#mu,#mu) (GeV)"
    cfg['ytitle'] = "Events"
    

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

    canvas.SaveAs("%s/%s.png" % (outDir, hName))
    canvas.SaveAs("%s/%s.pdf" % (outDir, hName))
    canvas.Close()




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
    #bkgs = [["TT2l", "TT1l", "TT0l"], ["WJets0J", "WJets1J", "WJets2J", "WW", "WZ", "ZZ"], ["DY"]]
    bkgs = [["TT2l", "TT1l", "TT0l"], ["WJets0J", "WJets1J", "WJets2J", "WW", "WZ", "ZZ"], ["DY_MiNNLO"]]
    bkgs = [["TT2l", "TT1l", "TT0l"], ["WJets0J", "WJets1J", "WJets2J", "WW", "WZ", "ZZ"], ["DY_MiNNLO"]]
    #bkgs = [["DY_MiNNLO"]]
    
    labels_bkgs = ["TT", "EWK", "DY (Mg5_amc)"]
    labels_bkgs = ["TT", "EWK", "DY (MiNNLO+N3LL)"]
    #labels_bkgs = ["DY"]
    labels_data = "Data"

    colors_bkgs = [ROOT.kAzure+2]
    colors_bkgs = [ROOT.TColor.GetColor(222, 90, 106), ROOT.TColor.GetColor(100, 192, 232), ROOT.TColor.GetColor(248, 206, 104), ROOT.TColor.GetColor(155, 152, 204)]
    #colors_bkgs = [ROOT.TColor.GetColor(248, 206, 104)]
    

    # muon

    #rebin = [60, 63, 66, 69, 72, 75, 78, 80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,104,106,108,110,112,114,116,118,120]


    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Zmumu/"
    fIn = ROOT.TFile("/eos/user/j/jaeyserm/www/wmass/test/Zmumu.root")
    

    rebin=1
    singlePlot("m_mumu", "m_mumu", 60, 120, 1e0, 1e6, "m(#mu,#mu) (GeV)", "Events", rebin=rebin)
    #singlePlot("m_mumu", "m_mumu_nolog", 60, 120, 0, 2e4, "m(#mu,#mu) (GeV)", "Events", rebin=rebin, logY=False)
    sys.exit()
    rebin=2
    singlePlot("pt_mumu", "pt_mumu", 0, 100, 1e1, 1e6, "p_{T}(#mu,#mu) (GeV)", "Events", rebin=rebin)
    
    rebin=1
    singlePlot("leading_muon_pt", "leading_muon_pt", 0, 100, 1e1, 1e5, "Leading muon p_{T} (GeV)", "Events", rebin=rebin)
    singlePlot("subleading_muon_pt", "subleading_muon_pt", 0, 100, 1e1, 1e5, "Subleading muon p_{T} (GeV)", "Events", rebin=rebin)
    
    singlePlot("leading_muon_eta", "leading_muon_eta", -2.4, 2.4, 1e1, 1e6, "Leading muon #eta", "Events", rebin=rebin)
    singlePlot("subleading_muon_eta", "subleading_muon_eta", -2.4, 2.4, 1e1, 1e6, "Subleading muon #eta", "Events", rebin=rebin)
    singlePlot("leading_muon_abs_eta", "leading_muon_abs_eta", 0, 2.4, 1e1, 1e6, "Leading muon |#eta|", "Events", rebin=rebin)
    singlePlot("subleading_muon_abs_eta", "subleading_muon_abs_eta", 0, 2.4, 1e1, 1e6, "Subleading muon |#eta|", "Events", rebin=rebin)
    
    singlePlot("leading_muon_phi", "leading_muon_phi", -3.5, 3.5, 1e1, 1e6, "Leading muon #phi", "Events", rebin=rebin)
    singlePlot("subleading_muon_phi", "subleading_muon_phi", -3.5, 3.5, 1e1, 1e6, "Subleading muon #phi", "Events", rebin=rebin)
    
    #singlePlot("PV_npvs", "PV_npvs", 0, 10, 1e1, 1e5, "Number of primary vertices", "Events", rebin=rebin)
    
    fIn.Close()
    sys.exit()
    
    ## electron
    data = ["data"]
    bkgs = [["DY"]]
    labels_bkgs = ["DY"]
    colors_bkgs = [ROOT.TColor.GetColor(248, 206, 104)]
    labels_data = "Data"
    outDir = "/eos/user/j/jaeyserm/www/wmass/lowPU/Zee/"
    fIn = ROOT.TFile("/eos/user/j/jaeyserm/www/wmass/test/Zee.root")
    
    
    rebin=1
    singlePlot("m_ee", "m_ee", 60, 120, 1e0, 1e6, "m(e,e) (GeV)", "Events", rebin=rebin)
    singlePlot("pt_ee", "pt_ee", 0, 100, 1e1, 1e5, "p_{T}(e,e) (GeV)", "Events", rebin=rebin)
    
    singlePlot("leading_electron_pt", "leading_electron_pt", 0, 100, 1e1, 1e5, "Leading electron p_{T} (GeV)", "Events", rebin=rebin)
    singlePlot("subleading_electron_pt", "subleading_electron_pt", 0, 100, 1e1, 1e5, "Subleading electron p_{T} (GeV)", "Events", rebin=rebin)
    
    singlePlot("leading_electron_eta", "leading_electron_eta", -2.4, 2.4, 1e1, 1e6, "Leading electron #eta", "Events", rebin=rebin)
    singlePlot("subleading_electron_eta", "subleading_electron_eta", -2.4, 2.4, 1e1, 1e6, "Subleading electron #eta", "Events", rebin=rebin)
    singlePlot("leading_electron_abs_eta", "leading_electron_abs_eta", 0, 2.4, 1e1, 1e6, "Leading electron |#eta|", "Events", rebin=rebin)
    singlePlot("leading_electron_abs_eta", "leading_electron_abs_eta_lin", 0, 2.4, 0, 6e3, "Leading electron |#eta|", "Events", rebin=rebin, logY=False)
    singlePlot("subleading_electron_abs_eta", "subleading_electron_abs_eta", 0, 2.4, 1e1, 1e6, "Subleading electron |#eta|", "Events", rebin=rebin)
    
    singlePlot("leading_electron_phi", "leading_electron_phi", -3.5, 3.5, 1e1, 1e6, "Leading electron #phi", "Events", rebin=rebin)
    singlePlot("subleading_electron_phi", "subleading_electron_phi", -3.5, 3.5, 1e1, 1e6, "Subleading electron #phi", "Events", rebin=rebin)
    
    
    #singlePlot("PV_npvs", "PV_npvs", 0, 10, 1e1, 1e5, "Number of primary vertices", "Events", rebin=rebin)
    
    fIn.Close()
    
    sys.exit()

    #sel = "sel5_0" #["sel5_0", "sel5_1", "sel5_2", "sel5_3", "sel5_4", "sel5_5"]
    labels = ["All events", "#mu^{+}#mu^{#minus} pair", "86 < m_{#mu^{+}#mu^{#minus}} < 96", "20 < p_{T}^{#mu^{+}#mu^{#minus}} < 70", "|cos#theta_{missing}| < 0.98", "120 < m_{rec} < 140"]
    
    sigs = ["wzp6_ee_mumuH_ecm240", "wzp6_ee_tautauH_ecm240", "wzp6_ee_eeH_ecm240", "wzp6_ee_nunuH_ecm240", "wzp6_ee_qqH_ecm240"]
    sigLegend = "ZH" # "ZH #rightarrow #mu^{+}#mu^{#minus}H"
	
    bkgs = ["WW", "ZZ", "Zg", "RARE"] # , 
    bkgs_legends = ["W^{+}(#nu#mu^{+})W^{#minus}(#bar{#nu}#mu^{#minus})", "ZZ", "Z/#gamma^{*} #rightarrow l#bar{l}", "Rare (e(e)Z, #gamma#gamma#rightarrow#mu#mu,#tau#tau)"]
    bkgs_colors = [ROOT.kGreen+2, ROOT.kOrange+1, ROOT.kMagenta]
    bkgs_colors = [ROOT.TColor.GetColor(248, 206, 104), ROOT.TColor.GetColor(222, 90, 106), ROOT.TColor.GetColor(100, 192, 232), ROOT.TColor.GetColor(155, 152, 204)] # from https://github.com/cms-opendata-analyses/HiggsTauTauNanoAODOutreachAnalysis/blob/master/plot.py
    bgks_cfg = { # this is the order of the plot
        "WW"	: ["p8_ee_WW_mumu_ecm240"],
        "ZZ"	: ["p8_ee_ZZ_ecm240"],
        "Zg"    : ["p8_ee_Zll_ecm240"],
        "RARE"	: ["wzp6_egamma_eZ_Zmumu_ecm240", "wzp6_gammae_eZ_Zmumu_ecm240", "wzp6_gaga_mumu_60_ecm240", "wzp6_gaga_tautau_60_ecm240"],
    }

    mll("sel5_0")
    ptll("sel5_2")
    costhetamissing("sel5_3")
    recoil_log("sel5_4")
    recoil_lin("sel5_4")
    recoil_lin_zoom("sel5_4")
