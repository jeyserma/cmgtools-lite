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

parser = argparse.ArgumentParser()
#parser.add_argument('-n', '--name', dest='baseHistName', default='muon_eta_pt', type=str, help='Base name for histograms')
parser.add_argument('-x', '--xAxisName', default='Muon #eta', type=str, help='x axis name')
parser.add_argument('-y', '--yAxisName', default='Muon p_{T} (GeV)', type=str, help='y axis name')
parser.add_argument('-b', '--bins', dest="etaptBins", default='48,-2.4,2.4,29,26,55', type=str, help='Bins for eta-pt, passed as to TH2 (only supports uniform binning for now)')
parser.add_argument('--ptVar', default='Muon_pt[goodMuonsCharge][0]', type=str, help='Expression for variable on pt axis')
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
if args.outputFile != "":
    outf = open(args.outputFile,"w")
    printToFile = True
    
#nominal
#line = "{n}: {y}\:Muon_eta[goodMuonsCharge][0]: {b}; {axis} \n".format(n=baseHistName,y=args.ptVar,b=etaptBins,axis=axisNames)
line = f"nominal_: {args.ptVar}\:Muon_eta[goodMuonsCharge][0]: {etaptBins}; {axisNames} \n"

print(line)
if printToFile: outf.write(line+'\n')

# pdf
syst_key = "pdf"
syst_expr = f"indices(LHEPdfWeight)\:scalarToRVec({args.ptVar},LHEPdfWeight)\:scalarToRVec(Muon_eta[goodMuonsCharge][0],LHEPdfWeight)"
syst_binning = "102,0.5,102.5"
syst_axisname = "ZTitle='PDF index'"
syst_weight  = "LHEPdfWeight"
process_regexpr = "W.*|Z.*"

#line = "{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)
line = "{sk}_: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)

print(line)
if printToFile: outf.write(line+'\n')


# qcd scales (not Vpt binned, that one comes just afterward)
syst_key = "qcdScale"
syst_expr = f"indices(LHEScaleWeight)\:scalarToRVec({args.ptVar},LHEScaleWeight)\:scalarToRVec(Muon_eta[goodMuonsCharge][0],LHEScaleWeight)"
syst_binning = "18,-0.5,17.5"
syst_axisname = "ZTitle='QCD scale index'"
syst_weight  = "LHEScaleWeight"
process_regexpr = "W.*" if isWlike else "Z.*"

#line = "{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)
line = "{sk}_: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)

print(line)
if printToFile: outf.write(line+'\n')

# qcd scales (Vpt binned)
NVTPBINS = 10
process_regexpr =  "Z.*" if isWlike else "W.*" # opposite with respect to unbinned QCD scales
for ipt in range(1,1+NVTPBINS):
    syst_key = "qcdScaleVptBin%d" % ipt
    ptcut = wptBinsScales(ipt)
    syst_weight  = "qcdScaleWeight_VptBinned(LHEScaleWeight\,ptVgen\,{ptlo}\,{pthi})".format(ptlo=ptcut[0],pthi=ptcut[1])

    #line = "{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)
    line = "{sk}_: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)
    print(line)
    if printToFile: outf.write(line+'\n')

# eff. stat. nuisances, one nuisance per TnP bin, treated as uncorrelated
# function to use is _get_fullMuonSFvariation, which replace _get_fullMuonSF in the nominal weight, using ReplaceWeight
NTNPBINS = 576 # 48 eta * 12 pt, hardcoded for now. SF have 13 pt bins, but the last one is between 55 and 65, so outside the analysis acceptance. Lowest pt is 26, in case the acceptance starts at 30, it is better to let the bin start from 1, and then one can just drop some nuisances afterwards
START = 0.5
syst_key = "effStatTnP"
syst_expr = f"indices({NTNPBINS}, 1)\:scalarToRVec({args.ptVar},{NTNPBINS})\:scalarToRVec(Muon_eta[goodMuonsCharge][0],{NTNPBINS})"
syst_binning = "%d,%.1f,%.1f" % (NTNPBINS, START, NTNPBINS+START)
syst_axisname = "ZTitle='Eff. stat. nuisance index'"
ptOther  = "Muon_pt[goodMuonsOther][0]"  if isWlike else "-1"
etaOther = "Muon_eta[goodMuonsOther][0]" if isWlike else "-1"
syst_weight  = f"_get_fullMuonSF(Muon_pt[goodMuonsCharge][0]\,Muon_eta[goodMuonsCharge][0]\,Muon_charge[goodMuonsCharge][0]\,{ptOther}\,{etaOther}\,eraVFP)->_get_fullMuonSFvariation({NTNPBINS}\,Muon_pt[goodMuonsCharge][0]\,Muon_eta[goodMuonsCharge][0]\,Muon_charge[goodMuonsCharge][0]\,{ptOther}\,{etaOther}\,eraVFP)"
process_regexpr = "W.*|Z.*|Top|Dibosons"

#line = "{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)
line = "{sk}_: {se}: {b},{sb}; {axis}, {sa}, ReplaceWeight='{sw}', ProcessRegexp='{prg}' \n".format(sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)

print(line)
if printToFile: outf.write(line+'\n')

# prefiring uncertainty, uncorrelated for each eta bin. Here the function is _get_MuonPrefiringSFvariation, which replace _get_MuonPrefiringSF in the nominal weight, using ReplaceWeight
NETABINS = 16
START = 0.5
syst_key = "muonL1Prefire"
syst_expr = f"indices({NETABINS}, 1)\:scalarToRVec({args.ptVar},{NETABINS})\:scalarToRVec(Muon_eta[goodMuonsCharge][0],{NETABINS})"
syst_binning = "%d,%.1f,%.1f" % (NETABINS, START, NETABINS+START)
syst_axisname = "ZTitle='Muon L1 prefiring nuisance index'"
syst_weight  = f"_get_MuonPrefiringSF(Muon_eta\,Muon_pt\,Muon_looseId\,eraVFP)->_get_MuonPrefiringSFvariation({NETABINS}\,Muon_eta\,Muon_pt\,Muon_looseId\,eraVFP)"
process_regexpr = "W.*|Z.*|Top|Dibosons"

#line = "{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)
line = "{sk}_: {se}: {b},{sb}; {axis}, {sa}, ReplaceWeight='{sw}', ProcessRegexp='{prg}' \n".format(sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)

print(line)
if printToFile: outf.write(line+'\n')

# mass weights from LHEReweightingWeight (might change at some point)
NWEIGHTS = 44 if isWlike else 33
START = -0.5
syst_key = "massWeight"
syst_expr = f"indices(LHEReweightingWeight)\:scalarToRVec({args.ptVar},LHEReweightingWeight)\:scalarToRVec(Muon_eta[goodMuonsCharge][0],LHEReweightingWeight)"
syst_binning = "%d,%.1f,%.1f" % (NWEIGHTS, START, NWEIGHTS+START)
syst_axisname = "ZTitle='Mass weight index'"
syst_weight  = "LHEReweightingWeight"
process_regexpr = "Z.*" if isWlike else "W.*"

#line = "{n}_{sk}: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(n=baseHistName, sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)
line = "{sk}_: {se}: {b},{sb}; {axis}, {sa}, AddWeight='{sw}', ProcessRegexp='{prg}' \n".format(sk=syst_key, se=syst_expr, b=etaptBins, sb=syst_binning, axis=axisNames, sa=syst_axisname, sw=syst_weight, prg=process_regexpr)

print(line)
if printToFile: outf.write(line+'\n')


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
