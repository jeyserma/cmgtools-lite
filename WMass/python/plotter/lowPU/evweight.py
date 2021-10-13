
import sys, glob, os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def evweight(iFiles):
    
    

    # load all the input files to the tree
    tree = ROOT.TChain("Events")
    for file in iFiles: tree.Add(file)
    
    h_evweight = ROOT.TH1D("h_evweight", "h_evweight", 1, 0, 2)
    h_pu = ROOT.TH1D("h_pu", "h_pu", 100, 0, 100)
    #tree.Draw("1 >> h_evweight", "genWeight", "goff")
    tree.Draw("1 >> h_evweight", "genWeight", "goff")
    
    tree.Draw("Pileup_nPU >> h_pu", "Pileup_nPU", "goff")
    
    print(h_evweight.Integral(), h_evweight.GetBinContent(1))
    
    
    fOut = ROOT.TFile("lowPU/lumiPileup/pilup_mc.root", "RECREATE")
    h_pu.Write()
    
    print(h_pu.GetMean())
    #h_evweight.Write()
    fOut.Close()
    
    sys.exit()
    # loop over all the events in the inputfiles
    nEvents = tree.GetEntries()
    evw_tot = 0.



    #weightsum = ROOT.TH1D("weightsum", "Sum of mcWeights", 1, 0, 2)
    #tree.Draw("1 >> weightsum", "genWeight")


    #print("%s" % iFiles[0].split("/")[6])
    #print("evts = %d" % nEvents)
    #print("sumw = %f" % weightsum.GetBinContent(1))
    #print("xsec = %f" % xsecsum.GetBinContent(1))


      
    for i in range(0, nEvents):

        tree.GetEntry(i) # load event i in the memory
        
        evw = tree.genWeight
        evw_tot += evw
        
        #if i%100 == 0: print(evw)
       
    print()
        
    print("Events:", nEvents)  
    print("Event weight:", evw_tot)
    print(iFiles)

# From Mit-HEP
def findEOS(name, mount=""):

    #return FindEOSOld(name)
    new = []
    for root, directories, filenames in os.walk(name):
        for f in filenames:

            filePath = os.path.join(os.path.abspath(root), f)
            if "failed/" in filePath: continue
            if "log/" in filePath: continue
            new.append(filePath)

    return new


if __name__ == "__main__":


    files = findEOS("/eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/")
    print(len(files))
    evweight(files)


    '''
    data_sm:            .*;     SubPath="SingleMuon/",  FillColor=ROOT.kBlue, Label="Data (single Muon)"
    data_eg:            .*;     SubPath="HighEGJet/",   FillColor=ROOT.kBlue, Label="Data (EG)"

    TT:         .* :    87.31;      FillColor=ROOT.kGreen+2,    Label="t quark", SubPath="TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/"
    TT:         .* :    380.10;     FillColor=ROOT.kGreen+2,    Label="t quark", SubPath="TTToHadronic_TuneCP5_13TeV-powheg-pythia8/"
    TT:         .* :    364.35;     FillColor=ROOT.kGreen+2,    Label="t quark", SubPath="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/" 

    DY:         .* :    6077.22;    FillColor=ROOT.kAzure+2,    Label="DY", SubPath="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/"                   
                         
    WJets:      .* :    49155   ;   FillColor=ROOT.kRed+1,      Label="WJets",  SubPath="WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/"
    WJets:      .* :    8047.00 ;   FillColor=ROOT.kRed+1,      Label="WJets",  SubPath="WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/"
    WJets:      .* :    3160.00 ;   FillColor=ROOT.kRed+1,      Label="WJets",  SubPath="WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/"

    Diboson:    .* :    12.178;     FillColor=ROOT.kViolet,     Label="Diboson", SubPath="WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/"
    Diboson:    .* :    4.42965;    FillColor=ROOT.kViolet,     Label="Diboson", SubPath="WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/"
    Diboson:    .* :    0.564;      FillColor=ROOT.kViolet,     Label="Diboson", SubPath="ZZ_TuneCP5_13TeV-pythia8/"  

    '''