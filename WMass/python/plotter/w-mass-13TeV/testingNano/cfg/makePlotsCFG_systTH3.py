#!/usr/bin/env python3

import re
import os, os.path
import logging
import argparse
import shutil

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

def wptBinsScales(i):
    # 5% quantiles (to be redone on the new MC)
    wptbins = [0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 13000.0]
    if len(wptbins)<2*i:
        print('You are asking too much from the wpt binning for decorrelation of scales')
    ptlo = wptbins[2*(i-1)]
    pthi = wptbins[2*i]
    return [ptlo, pthi]

# TODO: Use this for other systematics (it needs to be generalized a bit)
def write3DHist(label, pt_expr, eta_expr, nsyst, etapt_binning, xylabels, weight_axis, regex, outfile=None, systBinStart=-0.5, indexStart=0, addWeight=None, replaceWeight=None, nBinsZaxis=None):

    if nBinsZaxis == None:
        nBinsZaxis = nsyst
    syst_binning = "%d,%.1f,%.1f" % (nBinsZaxis, systBinStart, nBinsZaxis+systBinStart)
    expr_string = f"indices({nsyst},{indexStart})\:scalarToRVec({pt_expr},{nsyst})\:scalarToRVec({eta_expr},{nsyst})"
    weight_items = []
    if addWeight:
        weight_items.append(f"AddWeight='{addWeight}'")
    if replaceWeight:
        weight_items.append(f"ReplaceWeight='{replaceWeight}'")
    weight_str = ", ".join(weight_items) if len(weight_items) else ""

    line = f"{label}_: {expr_string} : {etapt_binning},{syst_binning};" \
        f" {xylabels}, ZTitle='{weight_axis}', {weight_str}, ProcessRegexp='{regex}'\n"
    print(line)
    if outfile:
        outfile.write(line+'\n')

parser = argparse.ArgumentParser()
#parser.add_argument('-n', '--name', dest='baseHistName', default='muon_eta_pt', type=str, help='Base name for histograms')
parser.add_argument('-x', '--xAxisName', default='Muon #eta', type=str, help='x axis name')
parser.add_argument('-y', '--yAxisName', default='Muon p_{T} (GeV)', type=str, help='y axis name')
parser.add_argument('-b', '--bins', dest="etaptBins", default='48,-2.4,2.4,29,26,55', type=str, help='Bins for eta-pt, passed as to TH2 (only supports uniform binning for now)')
parser.add_argument('--ptVar', default='Muon_pt[goodMuonsCharge][0]', type=str, help='Expression for variable on pt axis')
parser.add_argument('--etaVar', default='Muon_eta[goodMuonsCharge][0]', type=str, help='Expression for variable on eta axis')
parser.add_argument('-a', '--analysis', choices=["wlike","wmass"], default="wmass", help='Analysis type (some settings are customized accordingly)')
parser.add_argument('-o', '--output', dest="outputFile", default='', type=str, help='Output file to store lines (they are also printed on stdout anyway)')
args = parser.parse_args()

###################################
# SOME BASIC CONFIGS #
######################
#baseHistName = args.baseHistName
axisNames = "XTitle='{x}', YTitle='{y}'".format(x=args.xAxisName,y=args.yAxisName)
etaptBins = args.etaptBins
isWlike = args.analysis == "wlike"
####################################

printToFile = False
outf = None
if args.outputFile != "":
    outf = open(args.outputFile,"w")
    printToFile = True
    
#nominal
#line = "{n}: {y}\:Muon_eta[goodMuonsCharge][0]: {b}; {axis} \n".format(n=baseHistName,y=args.ptVar,b=etaptBins,axis=axisNames)
line = f"nominal_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames} \n"

print(line)
if printToFile: outf.write(line+'\n')

# pdf + alphaS
write3DHist(label = "pdf",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 103, # for PDFs, len(LHEPdfWeight) == 103 because it has nominal + 102 weights (100 pdf + 2 alphaS)
            xylabels = axisNames,
            weight_axis = "PDF/alpha_{S} index",
            etapt_binning = etaptBins,
            regex = "W.*|Z.*",
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 0,
            addWeight = "LHEPdfWeight",
            nBinsZaxis = 102 # we want 102,0.5,102.5 (could have been 103,-0.5,102.5, but then histogram bin 1 would not be pdf1 but nominal, histgram bin 2 would not be pdf2 but pdf1 and so on)
)


# qcd scales (not Vpt binned, that one comes just afterward)
write3DHist(label = "qcdScale",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 18,
            xylabels = axisNames,
            weight_axis = "QCD scale index",
            etapt_binning = etaptBins,
            regex = "W.*" if isWlike else "Z.*",
            outfile = outf,
            systBinStart = -0.5,
            indexStart = 0,
            addWeight = "LHEScaleWeight" 
)

# qcd scales (Vpt binned)
NVTPBINS = 10
process_regexpr =  "Z.*" if isWlike else "W.*" # opposite with respect to unbinned QCD scales
for ipt in range(1,1+NVTPBINS):
    ptcut = wptBinsScales(ipt)
    write3DHist(label = "qcdScaleVptBin%d" % ipt,
                pt_expr = args.ptVar,
                eta_expr = args.etaVar,
                nsyst = 18,
                xylabels = axisNames,
                weight_axis = "QCD scale index",
                etapt_binning = etaptBins,
                regex = "Z.*" if isWlike else "W.*", # opposite with respect to unbinned QCD scales 
                outfile = outf,
                systBinStart = -0.5,
                indexStart = 0,
                addWeight = f"qcdScaleWeight_VptBinned(LHEScaleWeight\,ptVgen\,{ptcut[0]}\,{ptcut[1]})" 
    )

# eff. stat. nuisances, one nuisance per TnP bin, treated as uncorrelated
# function to use is _get_fullMuonSFvariation, which replace _get_fullMuonSF in the nominal weight, using ReplaceWeight
write3DHist(label = "effStatTnP",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 576, # remember to edit weight below if changing this 
            xylabels = axisNames,
            weight_axis = "Eff. stat. nuisance index",
            etapt_binning = etaptBins,
            regex = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            replaceWeight = f"_get_fullMuonSF(->_get_fullMuonSFvariation(576\," 
)


# prefiring uncertainty, uncorrelated for each eta bin. Here the function is _get_MuonPrefiringSFvariation, which replace _get_MuonPrefiringSF in the nominal weight, using ReplaceWeight
write3DHist(label = "muonL1Prefire",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 16, # remember to edit weight below if changing this 
            xylabels = axisNames,
            weight_axis = "Muon L1 prefiring nuisance index",
            etapt_binning = etaptBins,
            regex = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            replaceWeight = f"_get_MuonPrefiringSF(Muon_eta\,Muon_pt\,Muon_looseId\,eraVFP)->_get_MuonPrefiringSFvariation(16\,Muon_eta\,Muon_pt\,Muon_looseId\,eraVFP)" 
)


# mass weights from LHEReweightingWeightCorrectMass
write3DHist(label = "massWeight",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst=23 if isWlike else 21,
            xylabels = axisNames,
            weight_axis = "Mass weight index",
            etapt_binning = etaptBins,
            regex = "Z.*" if isWlike else "W.*",
            outfile = outf,
            systBinStart = -0.5,
            indexStart = 0,
            addWeight = "LHEReweightingWeightCorrectMass"
)


for chan in ["Wm", "Wp", "Z"]:
    process_regexpr = "Z.*" 
    if "W" in chan:
        process_regexpr = ".*".join(["W", "minus" if "m" in chan else "plus", ])
    START = -0.5
    NWEIGHTS = 45
    write3DHist(label ="scetlibWeights", # no trailing '_', it is added inside directly (and mcPlots will add another one) 
                pt_expr =args.ptVar,
                eta_expr = args.etaVar,
                nsyst = NWEIGHTS,
                xylabels = axisNames,
                weight_axis = "SCETlib variation index",
                etapt_binning = etaptBins,
                regex = process_regexpr,
                outfile = outf,
                systBinStart = START,
                indexStart = 0,
                addWeight = f"scetlibWeights{chan}"
    )

print('-'*30)
print("SUMMARY")
print('-'*30)
print("Analysis:       %s" % args.analysis)
#print("Base hist name: %s" % baseHistName)
print("Binning eta-pt: %s" % etaptBins)
print('-'*30)
if printToFile:
    outf.close()
    print("Output saved in file %s" % args.outputFile)
    print()
