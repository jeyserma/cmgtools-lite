#!/usr/bin/env python3

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT
import numpy as np
import argparse

# processes to group
def getProcessGroup(proc, era):
    # check whether Wmunu is split by charge or includes both samples, same for Wtaunu
    if any(x in proc for x in ["Wmunu", "Wtaunu"]) and all(c not in proc for c in ["plus", "minus"]):
        if era == "all":
            ret = f"--pg '{proc} := {proc}_plus_preVFP,{proc}_plus_postVFP,{proc}_minus_preVFP,{proc}_minus_postVFP'"
        else:
            ret = f"--pg '{proc} := {proc}_plus_{era},{proc}_minus_{era}'"
    else:
        if era == "all":
            ret = f"--pg '{proc} := {proc}_preVFP,{proc}_postVFP'"
        else:
            ret = f"--pg '{proc} := {proc}_{era}'"
    return ret


def main(args):

    era = args.era
    postfix = args.postfix
    if len(args.postfix) and not postfix.startswith('_'):
        postfix = "_" + postfix
    
    ######################################################################
    ## Configurations below
    ######################################################################
    lumi = 16.8 if era == "postVFP" else 19.5 if era == "preVFP" else 36.3

    # cfg files
    cfgFolder = args.cfgFolder
    mca    = cfgFolder + "mca-wmass.txt"
    cut    = cfgFolder + "test/cuts_fakerate.txt"
    plot   = cfgFolder + args.plotFile
    define = cfgFolder + "test/rdfDefine_fakerate.txt"

    # additional cuts
    addcut = ""
    if args.charge != "all":
        chsign = ">" if args.charge == "plus" else "<"
        addcut = f"-A onemuon charge{args.charge} 'Muon_charge[goodMuons][0] {chsign} 0'"
    #addcut = "-A awayjet pfRelIso04 'Muon_pfRelIso04_all[goodMuons][0] < 0.15' "

    # additional defines and aliases
    otherDefines = " --rdf-alias 'goodMuonsCharge: goodMuons:.*'"

    # input samples
    samples = "/data/shared/originalNANO/"

    # output
    plotdir = f"plots/testNanoAOD/WmassPlots/fakeRateRegion_{era}_{args.charge}{postfix}/"

    # histograms to make (use .* to activate all those in plot file)
    hists = args.variables
    #hists = "muon_pt_eta_isoMtRegions"
    
    # processes and related options (also for ratio plots when customizing names)
    processes = "data,Zmumu,Ztautau,Wmunu,Wtaunu,Top,Diboson"
    procGroups = " ".join([getProcessGroup(p,era) for p in processes.split(',')])
    procOptions = f"-p '{processes}' {procGroups}"

    # ratio (settings, while customization of processes should go in procOptions)
    ratio = "--showRatio --maxRatioRange 0.5 1.5 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel 'Data/MC'"

    # event weight (global one for MC)
    weight = f"puw_2016UL_era(Pileup_nTrueInt,eraVFP)*_get_fullMuonSF(Muon_pt[goodMuonsCharge][0],Muon_eta[goodMuonsCharge][0],Muon_charge[goodMuonsCharge][0],-1,-1,eraVFP,Muon_pfRelIso04_all[goodMuons][0]<0.15)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)"

    # gen weight customization
    genweight = "--max-genWeight-procs 'W|Z' '50118.72' --clip-genWeight-toMax"

    # whatever with legend
    legOptions = "--legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92"

    # general options not in a specific group
    general = "-f  -j 8  --nanoaod-tree -v 3"

    ######################################################################
    ## Finally the command
    ######################################################################
    
    command = f"python mcPlots.py -l {lumi} {mca} {cut} {plot} --noCms -P {samples} --sP '{hists}'   -W '{weight}' --pdir {plotdir} {procOptions} {legOptions}  {genweight}  --rdf-define-file {define} {otherDefines} {addcut} {ratio} {general} {args.options}"

    if args.dryRun:
        print(command)
    else:
        os.system(command)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dry-run', dest="dryRun", action='store_true', help='Print, but do not execute')
    parser.add_argument('-s', '--run-systs',   dest="systs",  action='store_true', help='Special config for analysis, using all systematics on other processes, otherwise make plots in all the 4 regions separately')
    parser.add_argument('-e', '--era',     type=str, default="all", choices=["all","preVFP","postVFP"], help='Era')
    parser.add_argument('-c', '--charge',  type=str, default="all", choices=["all","plus","minus"], help='Era')
    parser.add_argument('-p', '--postfix', type=str, default="", help='Postfix to output folder name')
    parser.add_argument(      '--options', type=str, default="", help='Other options to pass to command, if not already present')
    parser.add_argument(      '--variables', type=str, default="muon_pt_eta,muon_pt,muon_eta_fine,mt_MET,MET_pt,nJetClean", help='Histograms to make')
    parser.add_argument(      '--cfg-folder', dest="cfgFolder", type=str, default="w-mass-13TeV/testingNano/cfg/", help='Folder where cfg files are taken. Can leave this default')
    parser.add_argument(      '--plot-file', dest="plotFile", type=str, default="plots_fakerate.txt", help='File with histogram definition (inside cfgFolder)')
    args = parser.parse_args()
    
    if args.systs:
        args.postfix += "_systTH3"
        regionCut = "transverseMass >= 40.0 || Sum(goodCleanJets)>=1"
        args.options += f" -A trigMatch regionCut '{regionCut}' "
        print("")
        main(args)
        print('-'*30)
        print("")
    else:
        mT_expr = "transverseMass" #"mt_2(Muon_pt[goodMuons][0],Muon_phi[goodMuons][0],MET_pt,MET_phi)"
        iso_expr = "Muon_pfRelIso04_all[goodMuons][0]"
        regionIsoMt_cuts = {"lowIso_lowMt"   : f"{iso_expr} < 0.15 && {mT_expr} < 40 && Sum(goodCleanJets)>=1",
                            "highIso_lowMt"  : f"{iso_expr} > 0.15 && {mT_expr} < 40 && Sum(goodCleanJets)>=1",
                            "lowIso_highMt"   : f"{iso_expr} < 0.15 && {mT_expr} > 40",
                            "highIso_highMt"  : f"{iso_expr} > 0.15 && {mT_expr} > 40"
        }
        postfix_backup = args.postfix
        options_backup = args.options
        for k in regionIsoMt_cuts.keys():
            args.postfix = f"{postfix_backup}_{k}" 
            args.options = options_backup + f" -A trigMatch regionCut '{regionIsoMt_cuts[k]}' "
            print("")
            main(args)
            print('-'*30)
            print("")

