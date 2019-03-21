import ROOT, os, sys, re, array, math
import time

ROOT.gROOT.SetBatch(True)

dryrun = 0
doMuons = 0

if doMuons:

    inputFolder = "fr_10_03_2019_eta_pt_finer_mT40_35p9fb_signedEta_subtrAllMC_L1EGprefire_jetPt30"
    inputFullPath = "www/wmass/13TeV/fake-rate/test_mu/testFRv8/" + inputFolder + "/mu/comb/"
    outputFolder = "www/wmass/13TeV/fake-rate/muon/FR_graphs_tests/fr_fitPol1_xMaxFit55_withLepSF"
    outfileTag = outputFolder.split('/')[-2]
    histPrefix = "fakeRateNumerator_mu_vs_etal1mu_pt_finer"
    isMuon = "true"
    showMergedEWK = "true"
    saveToFile = "false"
    noDrawQCD = "false"
    etaBinBoundariesList = "-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4"
    #etaBinBoundariesList = "-2.4,-2.2,-2.05,-1.9,-1.75,-1.6,-1.45,-1.3,-1.15,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.15,1.3,1.45,1.6,1.75,1.9,2.05,2.2,2.4"

else:

    inputFolder = "fr_12_03_2019_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_L1EGprefire_jetPt30_Zveto_newSkim_noDensity"
    inputFullPath = "www/wmass/13TeV/fake-rate/test/testFRv8/" + inputFolder + "/el/comb/"
    outputFolder = "www/wmass/13TeV/fake-rate/electron/FR_graphs_tests/fr_fitPol2/"
    outfileTag = outputFolder.split('/')[-2]
    histPrefix = "fakeRateNumerator_el_vs_etal1_pt_granular"
    isMuon = "false"
    showMergedEWK = "true"
    saveToFile = "false"
    noDrawQCD = "true"
    etaBinBoundariesList = "-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.566,-1.4442,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4442,1.566,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5"


# work in progress 

options = '"{i}","{o}","{f}",'.format(i=inputFullPath, o=outputFolder, f=outfileTag)
options += '"{h}",{ism}, {smewk}, {save}, "{etabin}", {qcd}'.format(h=histPrefix, ism=isMuon, smewk=showMergedEWK, 
                                                                    save=saveToFile, etabin=etaBinBoundariesList, qcd=noDrawQCD)

cmd = "root -l -b -q 'makeFakeRateGraphPlotsAndSmoothing.C++({opt})' ".format(opt=options)
print "-"*30
print cmd
print "-"*30

if not dryrun:
    os.system(cmd)
