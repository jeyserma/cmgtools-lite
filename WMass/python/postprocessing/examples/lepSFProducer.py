import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from CMGTools.WMass.postprocessing.framework.datamodel import Collection 
from CMGTools.WMass.postprocessing.framework.eventloop import Module

class lepSFProducer(Module):
    def __init__(self, muonSelectionTag, electronSelectionTag):
        if muonSelectionTag=="TightWP_2016":
            mu_f=["Mu_Trg.root","Mu_ID.root","Mu_Iso.root"]
            mu_h = ["IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio",
                    "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio",
                    "TightISO_TightID_pt_eta/pt_abseta_ratio"]
        if electronSelectionTag=="GPMVA90_2016":
            el_f = ["EGM2D_eleGSF.root","EGM2D_eleMVA90.root"]
        elif electronSelectionTag=="CutBasedTight_2016":
            el_f = ["EGM2D_eleGSF.root","EGM2D_eleCutBasedTightWP.root"]
        elif electronSelectionTag=="CutBasedMedium_2016":
            el_f = ["EGM2D_eleGSF.root","EGM2D_eleCutBasedMediumWP.root"]
        else:
            print "Not foreseen WP: ",electronSelectionTag
            el_f = []
        el_h = ["EGamma_SF2D", "EGamma_SF2D"]
        mu_f = ["%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in mu_f]
        el_f = ["%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in el_f]

        self.mu_f = ROOT.std.vector(str)(len(mu_f))
        self.mu_h = ROOT.std.vector(str)(len(mu_f))
        for i in range(len(mu_f)): self.mu_f[i] = mu_f[i]; self.mu_h[i] = mu_h[i];
        self.el_f = ROOT.std.vector(str)(len(el_f))
        self.el_h = ROOT.std.vector(str)(len(el_f))
        for i in range(len(el_f)): self.el_f[i] = el_f[i]; self.el_h[i] = el_h[i];

        if "/LeptonEfficiencyCorrector_cc.so" not in ROOT.gSystem.GetLibraries():
            print "Load C++ Worker"
            ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/LeptonEfficiencyCorrector.cc+" % os.environ['CMSSW_BASE'])
    def beginJob(self):
        self._worker_mu = ROOT.LeptonEfficiencyCorrector(self.mu_f,self.mu_h)
        self._worker_el = ROOT.LeptonEfficiencyCorrector(self.el_f,self.el_h)
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LepGood_effSF", "F", lenVar="nLepGood")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        leps = Collection(event, "LepGood")
        sf = []
        for l in leps:
            if event.isData:
                sf.append(1.)
            else:
                worker = self._worker_el if abs(l.pdgId)==11 else self._worker_mu
                sf.append(worker.getSF(l.pdgId,l.pt,l.eta))
        self.out.fillBranch("LepGood_effSF", sf)
        return True

class lepTrgSFProducer(Module):
    def __init__(self,dimensions=2,prefer1Dvar="eta",maxrun=278808):
        # from  https://indico.cern.ch/event/570616/contributions/2354285/attachments/1365274/2067939/tnpBuda_nov03_v1.pdf
        # should make the weighted average in the data chunk considered.
        # here take the last measurement, since similar
        self.versions = {(273158,274442): "v1",
                         (274954,275066): "v2",
                         (275067,275311): "v3",
                         (275319,276834): "v4",
                         (276870,278240): "v5",
                         (278273,280385): "v6",
                         (281639,283059): "v7"}
        self.maxrun = maxrun
        self.dim = dimensions
        self.var = prefer1Dvar
        if "/WeightCalculatorFromHistogram_cc.so" not in ROOT.gSystem.GetLibraries():
            print "Load C++ Worker"
            ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/WeightCalculatorFromHistogram.cc+" % os.environ['CMSSW_BASE'])
    def beginJob(self):
        ver=None
        for k,v in self.versions.iteritems(): 
            if k[0] < self.maxrun < k[1]: ver = v
        if ver:
            f = "%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/el_trg/%s/%s/passHLT/eff1D.root" % (os.environ['CMSSW_BASE'],ver,self.var)
            h = "s1c_eff"
            if self.dim==2:
                f2D = "%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/el_trg/%s/sf/passHLT/eff2D.root" % (os.environ['CMSSW_BASE'],ver)
                if os.path.isfile(f2D):
                    f = f2D
                    h = "s2c_eff"
                else: self.dim = 1
        else: raise Exception('No suitable version of trigger scale factors found!')
        print "Reading trigger scale factors in %d dimensions from file %s..." % (self.dim,f)
        tf = ROOT.TFile.Open(f)
        th = tf.Get(h).Clone("sf_%s" % v)
        th.SetDirectory(None)
        self._worker = ROOT.WeightCalculatorFromHistogram(th)
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LepGood_trgSF", "F", lenVar="nLepGood")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        leps = Collection(event, "LepGood")
        sf = []
        for l in leps:
            if event.isData or abs(l.pdgId)!=11:
                sf.append(1.)
            else:
                if self.dim == 2:
                    wgt = self._worker.getWeight(l.pt,l.eta)
                    sf.append(wgt if wgt>0 else 1.)
                else:
                    sf.append(self._worker.getWeight(getattr(l,self.var)))
        self.out.fillBranch("LepGood_trgSF", sf)
        return True

####################################################################
####################################################################

class scaleFactorManager:
    def __init__(self,filename,path,hname):
        self.hname = hname
        self.fname = path + ("" if path.endswith("/") else "/") + filename
        self.hist = 0
        self.epsilon = 0.0001 # used to get correct bin number given value on axis (just in case we are picking a value on a bin edge)
        self.hasLoadedHisto = False

    def printFile():
        tf = ROOT.TFile.Open(self.fname)
        print "-"*20
        tf.Print()
        tf.Close()
        print "-"*20

    def loadHist(self,newfname=None,newhname=None):
        # can also change file and/or histogram name with this method
        if newfname:
            self.fname = newfname
        if newhname:
            self.hname = newhname
        tf = ROOT.TFile.Open(self.fname)
        self.hist = tf.Get(self.hname)
        self.hist.SetDirectory(0)
        tf.Close()
        if not self.hist:
            print "*"*20
            print "WARNING: could not load histogram {n} from file {f}. Will return 1".format(n=self.hname,f=self.fname)
            print "*"*20
            self.printFile()
        else:
            print "Histogram {n} successfully loaded from file {f}  ;-)".format(n=self.hname,f=self.fname)
            self.hasLoadedHisto = True

    def getSF(self,pt,eta):
        if not self.hist:
            #print "Warning in scaleFactorManager.getSF(): histogram not found, I will try to load it now. If not successfull, I will return 1."
            # just in case it was not loaded explicitly
            if not self.hasLoadedHisto:
                self.loadHist()
                return self.getSF(pt,eta)
            else: 
                return 1.
        else:
            ieta = self.hist.GetXaxis().FindFixBin(eta+self.epsilon)
            ipt =  self.hist.GetYaxis().FindFixBin(pt +self.epsilon)
            # protect against underflow and overflow
            return self.hist.GetBinContent(min(max(1,ieta),self.hist.GetNbinsX()),min(max(1,ipt),self.hist.GetNbinsY()))

    def getSF_err(self,pt,eta):
        if not self.hist:
            #print "Warning in scaleFactorManager.getSF_err(): histogram not found, I will try to load it now. If not successfull, I will return 1."
            # just in case it was not loaded explicitly
            if not self.hasLoadedHisto:
                self.loadHist()
                return self.getSF(pt,eta)
            else: 
                return 1.
        else:
            ieta = self.hist.GetXaxis().FindFixBin(eta+self.epsilon)
            ipt =  self.hist.GetYaxis().FindFixBin(pt +self.epsilon)
            return self.hist.GetBinError(min(max(1,ieta),self.hist.GetNbinsX()),min(max(1,ipt),self.hist.GetNbinsY()))
        
#-----------------------------------------

class lep2016SFProducer(Module):
    def __init__(self):

        self.mu_f = {"trigger"       :"smoothEfficiency_muons_trigger.root", 
                     "identification":"smoothEfficiency_muons_full2016_ID.root", 
                     "isolation"     :"smoothEfficiency_muons_full2016_ISO.root"
                     }
        self.el_f = {"trigger"       :"smoothEfficiency_electrons_trigger.root"
                     }
        self.filePath = "%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/" % os.environ['CMSSW_BASE']

    def beginJob(self):
        # create muon scale factor manager: pass file name and location, and then the name of histogram to read
        # for muon ID, might also want to use "scaleFactorOriginal" which has the unsmoothed version of the scale factors
        self.trgSF_manager_mu = scaleFactorManager(self.mu_f["trigger"],       self.filePath,"scaleFactor")
        self.idSF_manager_mu  = scaleFactorManager(self.mu_f["identification"],self.filePath,"scaleFactor")  
        self.isoSF_manager_mu = scaleFactorManager(self.mu_f["isolation"],     self.filePath,"scaleFactor")

        # create electron scale factor manager        
        self.trgSF_manager_el = scaleFactorManager(self.el_f["trigger"],       self.filePath,"scaleFactor")

        # load histograms
        self.trgSF_manager_mu.loadHist()
        self.idSF_manager_mu.loadHist()
        self.isoSF_manager_mu.loadHist()
        self.trgSF_manager_el.loadHist()

    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LepGood_IDeffSF",      "F", lenVar="nLepGood")
        self.out.branch("LepGood_ISOeffSF",     "F", lenVar="nLepGood")
        self.out.branch("LepGood_trgSF",        "F", lenVar="nLepGood")
        self.out.branch("LepGood_IDeffSF_err",  "F", lenVar="nLepGood")
        self.out.branch("LepGood_ISOeffSF_err", "F", lenVar="nLepGood")
        self.out.branch("LepGood_trgSF_err",    "F", lenVar="nLepGood")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        leps = Collection(event, "LepGood")
        sf_id  = []
        sf_iso = []
        sf_trg = []
        sf_id_err  = []
        sf_iso_err = []
        sf_trg_err = []
        #print ">>>>>>> check"
        #i = 0
        for l in leps:
            #print "Inside loop: %d" % i
            #i += 1
            if event.isData:
                sf_trg.append(1.)
                sf_id.append(1.)
                sf_iso.append(1.)
                sf_trg_err.append(0.)
                sf_id_err.append(0.)
                sf_iso_err.append(0.)
            else:
                if abs(l.pdgId)==11:
                    #print "Electron found"
                    sf_trg.append(float(self.trgSF_manager_el.getSF(l.pt,l.eta)))                    
                    sf_id.append(1.)
                    sf_iso.append(1.)
                    sf_trg_err.append(float(self.trgSF_manager_el.getSF_err(l.pt,l.eta)))                    
                    sf_id_err.append(1.)
                    sf_iso_err.append(1.)
                else:
                    # print "Muon found"
                    # print "SF: trg - id - iso --> %.3f - %.3f - %.3f " % (float(self.trgSF_manager_mu.getSF(l.pt,l.eta)),
                    #                                                       float(self.idSF_manager_mu.getSF(l.pt,l.eta)),
                    #                                                       float(self.isoSF_manager_mu.getSF(l.pt,l.eta))
                    #                                                       )                
                    sf_trg.append(    float(self.trgSF_manager_mu.getSF(l.pt,l.eta)))
                    sf_id.append(     float(self.idSF_manager_mu.getSF(l.pt,l.eta)))
                    sf_iso.append(    float(self.isoSF_manager_mu.getSF(l.pt,l.eta)))
                    sf_trg_err.append(float(self.trgSF_manager_mu.getSF_err(l.pt,l.eta)))
                    sf_id_err.append( float(self.idSF_manager_mu.getSF_err(l.pt,l.eta)))
                    sf_iso_err.append(float(self.isoSF_manager_mu.getSF_err(l.pt,l.eta)))
        self.out.fillBranch("LepGood_IDeffSF",  sf_id)
        self.out.fillBranch("LepGood_ISOeffSF", sf_iso)
        self.out.fillBranch("LepGood_trgSF",    sf_trg)
        self.out.fillBranch("LepGood_IDeffSF_err",  sf_id_err)
        self.out.fillBranch("LepGood_ISOeffSF_err", sf_iso_err)
        self.out.fillBranch("LepGood_trgSF_err",    sf_trg_err)
        #print ">>>>>>> check after fill branch"        
        return True

#########################################################


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

lepSF = lambda : lepSFProducer( "TightWP_2016", "CutBasedTight_2016")
trgSF = lambda : lepTrgSFProducer()
lep2016SF = lambda : lep2016SFProducer()
