#!/bin/env python

# script to make ratio of two TH2 histograms
# the two objects can have different bin ranges, but in order for the ratio to make sense, 
# the binning of the smaller TH2 should be a subset of the other histogram's binning
# In any case, the ratio is made looping on the first histogram passed (which is cloned to get a new one with the same binning) and getting the bin center
# for each bin, then that value is used to look for the bin of the other histogram
# therefore, in case the two TH2 have different binning, make sure that the first one has the finer granularity

###################
# >>>>>>>>> WARNING <<<<<<<<<<<< 
###################
# Apparently it works in release CMSSW_10_2_0_pre4
# while in CMSSW_8_0_25 it produces a segmentatio fault from the setTDRStyle()
# Maybe it is just due to the root version

################################
# Exaples
################################

# python w-helicity-13TeV/makeRatioTH2.py ~/www/wmass/13TeV/scaleFactors/electron/fullID_extPt/smoothEfficiency_electrons_fullID.root scaleFactor ~/www/wmass/13TeV/scaleFactors/electron/fullID_noErfPlusLine/smoothEfficiency_electrons_fullID.root scaleFactor -o ~/www/wmass/13TeV/scaleFactors/ratio/electron/fullID_extPt__over__fullID/ -f ratio.root -n ratio2D -t "full ID scale factor ratio" -z ratio --ratioRange 0.9 1.1

# python w-helicity-13TeV/makeRatioTH2.py ../../data/fakerate/fakeRateSmoothed_el_mT40_35p9fb_signedEta_pt65_fullWMC_newTrigSF_fitpol2.root fr_pt_eta_ewk ~/www/wmass/13TeV/fake-rate/electron/FR_graphs/fakeRate_eta_pt_granular_mT40_35p9fb_signedEta_pt65_subtrAllMC_newTrigSF_fitpol2/fakeRateSmoothed_el_mT40_35p9fb_signedEta_pt65_subtrAlllMC_newTrigSF_fitpol2.root fr_pt_eta_ewk -o ~/www/wmass/13TeV/fake-rate/electron/FR_graphs/fakeRate_eta_pt_granular_mT40_35p9fb_signedEta_pt65_subtrAllMC_newTrigSF_fitpol2/ -f ratio_PR_WZ_allMC.root -n ratioPR_WZ_allMC -t "PR W,Z MC / all MC" -z "PR ratio" -r 0.98 1.02

# python w-helicity-13TeV/makeRatioTH2.py ../../data/fakerate/fakeRateSmoothed_el_mT40_35p9fb_signedEta_pt65_subtrAlllMC_newTrigSF_fitpol2.root fr_pt_eta_ewk ~/www/wmass/13TeV/fake-rate/electron/FR_graphs/fakeRate_eta_pt_granular_mT40_35p9fb_signedEta_pt65_subtrAllMC_newTrigSF_fitpol2_testTrigSF/fakeRateSmoothed_el_mT40_35p9fb_signedEta_pt65_subtrAlllMC_newTrigSF_fitpol2_testTrigSF.root fr_pt_eta_ewk -o ~/www/wmass/13TeV/fake-rate/electron/FR_graphs/fakeRate_eta_pt_granular_mT40_35p9fb_signedEta_pt65_subtrAllMC_newTrigSF_fitpol2/ -f ratioPR_nomi_test.root -n FILE -t "PR nomi / test" -z "PR ratio" -r 0.99 1.01

# python w-helicity-13TeV/makeRatioTH2.py /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/scaleFactors/electron/trigger_extPt_30_45/smoothEfficiency_electrons_trigger.root scaleFactor /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/scaleFactors/electron/trigger_extPt/smoothEfficiency_electrons_trigger.root scaleFactor -o /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/scaleFactors/ratio/electron/trigger_ptUpTo45__Over__ptUpTo55/ -f trgSF_pt45_pt55.root -n FILE -t "p_{T} < 45 / p_{T} < 55" -z "trigger scale factor ratio" -r 0.95 1.05

# FR and PR ratio

# python w-helicity-13TeV/makeRatioTH2.py /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/fake-rate/electron/FR_graphs_tests/fr_16_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32/fakeRateSmoothed_el_fr_16_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32.root fr_pt_eta_data /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/fake-rate/electron/FR_graphs_tests/fr_18_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32_jetPt45/fakeRateSmoothed_el_fr_18_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32_jetPt45.root fr_pt_eta_data -o /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/fake-rate/electron/ratio_FR_PR/subtrAllMC_allNewSF_fitPol1_minPtData32__jetPt30_Over_jetPt45/ -f ratio_FR -n FILE -t "FR jet p_{T} > 30 / jet p_{T} > 45" -z "Fake rate ratio" -r 0.9 1.1

# python w-helicity-13TeV/makeRatioTH2.py /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/fake-rate/electron/FR_graphs_tests/fr_16_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32/fakeRateSmoothed_el_fr_16_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32.root fr_pt_eta_ewk /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/fake-rate/electron/FR_graphs_tests/fr_18_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32_jetPt45/fakeRateSmoothed_el_fr_18_09_2018_eta_pt_granular_mT40_35p9fb_signedEta_subtrAllMC_allNewSF_fitPol1_minPtData32_jetPt45.root fr_pt_eta_ewk -o /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/fake-rate/electron/ratio_FR_PR/subtrAllMC_allNewSF_fitPol1_minPtData32__jetPt30_Over_jetPt45/ -f ratio_PR -n FILE -t "PR jet p_{T} > 30 / jet p_{T} > 45" -z "Prompt rate ratio" -r 0.98 1.04

# muon FR (need to build the FR versu pt and eta
# python w-helicity-13TeV/makeRatioTH2.py ../../data/fakerate/frAndPr_fit_mu_2018-09-13_finerETA.root fakerates_smoothed_data_interpolated ../../data/fakerate/frAndPr_fit_mu_2018-09-19_jetPt45_finerETA.root fakerates_smoothed_data_interpolated_awayJetPt45 -o /afs/cern.ch/user/m/mciprian/www/wmass/13TeV/fake-rate/muon/ratio_FR_PR/fromMarc_jetPt30_Over_jetPt45/ -f ratio_FR -n FILE -t "FR jet p_{T} > 30 / jet p_{T} > 45" -z "Fake rate ratio" -r 0.9 1.2   --buildFakeRate -x "muon p_{T} [GeV]" -y "muon #eta" --h1Dbinning "75,0.9,1.2"

# python w-helicity-13TeV/makeRatioTH2.py plots/distribution/muonPlots/SKIMS_muons_latest/ptVsEta_W_metJecUp_plus/test_plots.root ptl1__etal1_W plots/distribution/muonPlots/SKIMS_muons_latest/ptVsEta_W_metNominal_plus/test_plots.root ptl1__etal1_W -o plots/distribution/muonPlots/SKIMS_muons_latest/ratio_metJecUpMT/ -f ratio_MetJecUp_nominal_W -n FILE -z "yields ratio for W^{+}" --ratioRange 0.97 1.005 -t "E_{T}^{miss}(JEC Up) / nominal"

################################
################################


import ROOT, os, sys, re, array, math
import time

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)
        
if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] file1 hist1 file2 hist2')
    parser.add_option('-o','--outdir',      dest='outdir',      default='', type='string', help='output directory to save things')
    parser.add_option('-f','--outfilename', dest='outfilename', default='', type='string', help='Name of output file to save results')
    parser.add_option('-n','--outhistname', dest='outhistname', default='', type='string', help='Name of output histogram saved in output file. If FILE is used, take same name as the output file, removing the extension')
    parser.add_option('-x','--xAxisTitle',  dest='xAxisTitle',  default='', type='string', help='X axis title. If not given, use the one from hist1')
    parser.add_option('-y','--yAxisTitle',  dest='yAxisTitle',  default='', type='string', help='Y axis title. If not given, use the one from hist1')
    parser.add_option('-z','--zAxisTitle',  dest='zAxisTitle',  default='', type='string', help='Z axis title. If not given, use the one from hist1')
    parser.add_option('-t','--histTitle',   dest='histTitle',   default='', type='string', help='Title to assign to output histogram. It is used as a label for the canvas')
    parser.add_option('-r','--ratioRange',  dest='ratioRange',  default=(0, 2),type="float", nargs=2, help="Min and max for the ratio in the plot")
    parser.add_option(     '--h1Dbinning',  dest='h1Dbinning',  default='50,0.9,1.1', type='string', help='Comma separated list of 3 numbers: nbins,min,max')
    parser.add_option('-v','--valBadRatio', dest='valBadRatio', default='0', type='float', help='Value to be used in case of bad ratio (division by 0). The 1D histogram is not filled in case of bad ratio')
    parser.add_option(     '--buildFakeRate', dest="buildFakeRate", action="store_true", default=False, help="The input histograms have the parameters of the linear fits to fake-rate or prompt-rate versus eta: build the histogram with FR (PR) vs pt and eta (obsolete, no longer using pol1 to interpolate)")
    parser.add_option(     '--xRange'     , dest='xRange', default=(0,-1), type='float', nargs=2, help='Select range for X axis to plot. Also, bins outside this range are not considered in the 1D histogram. If min > max, the option is neglected')
    parser.add_option(     '--yRange'     , dest='yRange', default=(0,-1), type='float', nargs=2, help='Select range for Y axis to plot. Also, bins outside this range are not considered in the 1D histogram. If min > max, the option is neglected')
    parser.add_option('-e', '--divide-error', dest="divideError", action="store_true", default=False, help="Make ratio of uncertainties (the output histogram will have no error assigned to it)")
    parser.add_option('-E',  '--divide-relative-error', dest="divideRelativeError", action="store_true", default=False, help="Make ratio of relative uncertainties (the output histogram will have no error assigned to it)")
    parser.add_option(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use a negative number to select a built-in one, otherwise the default is 55 (kRainbow)')
    parser.add_option('-a', '--make-asymmetry', dest="makeAsymmetry", action="store_true", default=False, help="Make ratio of difference over the sum. For this to make sense, the binnign of the two inputs must be consistent")
    parser.add_option('-p', '--make-pulls', dest="makePulls", action="store_true", default=False, help="Make pulls of input histograms, i.e. (h1-h2)/error, where error is taken as the quadrature sum of the errors of the input")
    parser.add_option(       '--pull-error-ScaleFactor', dest='pullErrorScaleFactor', default='1.', type='float', help='Inflate the error by this factor when making the pulls (because it is assumed the inputs are uncorrelated, so the error might need a correction)')
    parser.add_option(      '--roll1Dto2D', dest="roll1Dto2D", action="store_true", default=False, help="Input histograms are 1D distributions to be unrolled into 2D. Need binning from option --binning-file-to-roll")
    parser.add_option(      '--binning-file-to-roll', dest="binFileToRoll", default="", help="File with binning to roll 1D into 2D (the reco binning is used)")
    (options, args) = parser.parse_args()

    if len(sys.argv) < 4:
        parser.print_usage()
        quit()

    cmssw_version = os.environ['CMSSW_VERSION']
    if int(cmssw_version.split('_')[1]) < 10:
        print ""
        print ">>>> WARNING:" 
        print "current version of CMSSW ({version}) might lead to segmentation fault when running this script!".format(version=cmssw_version) 
        print "Or, the Z axis range might not be set properly"
        print "If you experience any of these issues, try running the script from a release CMSSW_X_Y_Z with X > 8 (10_2 suggested)"
        print "Versions with Y > 2 seems to work, but lead to some warnings and print garbage on stdout"
        print "Continuing in 3 seconds ..."
        print ""
        time.sleep(3)

    f1 = args[0]
    h1 = args[1]
    f2 = args[2]
    h2 = args[3]

    ROOT.TH1.SetDefaultSumw2()

    print ""

    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()
    if not options.outfilename:
        print "Error: you should specify an output file name using option -f <name>. Exit"
        quit()
    if not options.outhistname:
        print "Error: you should specify an output histogram name using option -n <name>. "
        print "If FILE is used, take same name as the output file, removing the extension"
        print "Exit"
        quit()

    if options.makeAsymmetry:
        if options.divideRelativeError or options.divideError:
            print "Error: option -a incompatible with options -e and -E. Exit"
            quit()
    if options.buildFakeRate: options.makePulls = False
    if options.makePulls:
        if options.divideRelativeError or options.divideError:
            print "Error: option -p incompatible with options -e and -E. Exit"
            quit()
    if options.makeAsymmetry and options.makePulls:
            print "Error: option -a incompatible with options -p. Exit"
            quit()


    if options.outhistname == "FILE":
        options.outhistname = options.outfilename.split('.')[0]

    hratio = 0

    # file 1
    tf = ROOT.TFile.Open(f1)        
    hist1 =   tf.Get(h1)
    if (hist1 == 0):
        print "Error: could not retrieve %s from input file %s. Exit" % (h1,f1)
        quit()
    else:
        hist1.SetDirectory(0)
    tf.Close()

    # file2
    tf = ROOT.TFile.Open(f2)        
    hist2 =   tf.Get(h2)
    if (hist2 == 0):
        print "Error: could not retrieve %s from input file %s. Exit" % (h2,f2)
        quit()
    else:
        hist2.SetDirectory(0)
    tf.Close()

    hinput1 = hist1
    hinput2 = hist2

    if options.buildFakeRate:
        # in this case the input TH2 have eta on x axis and offset/slope on the other one
        # this is needed only for muons, for electrons I save directly the smoothed FR (PR)
        neta = hist1.GetNbinsX()
        etabins = [hist1.GetXaxis().GetBinLowEdge(i) for i in range(1,2+neta)]
        #etamin = hist1.GetXaxis().GetBinLowEdge(1)
        #etamax = hist1.GetXaxis().GetBinLowEdge(1+neta)
        #hFR1 = ROOT.TH2D(hist1.GetName()+"_FRorPR","",195,26,65,neta,etamin,etamax)
        #hFR2 = ROOT.TH2D(hist2.GetName()+"_FRorPR","",195,26,65,neta,etamin,etamax)
        hFR1 = ROOT.TH2D(hist1.GetName()+"_FRorPR","",195,26,65,neta,array('d',etabins))
        hFR2 = ROOT.TH2D(hist2.GetName()+"_FRorPR","",195,26,65,neta,array('d',etabins))
        for ix in range (1,1+hFR1.GetNbinsX()):
            for iy in range (1,1+hFR1.GetNbinsY()):
                fr1 = hist1.GetBinContent(iy,1) + hist1.GetBinContent(iy,2) * hFR1.GetXaxis().GetBinCenter(ix)
                hFR1.SetBinContent(ix,iy,fr1)
                fr2 = hist2.GetBinContent(iy,1) + hist2.GetBinContent(iy,2) * hFR2.GetXaxis().GetBinCenter(ix)
                hFR2.SetBinContent(ix,iy,fr2)
        hinput1 = hFR1
        hinput2 = hFR2

    if options.roll1Dto2D:
      # TO DO   
        if options.binFileToRoll:
            etaPtBinningVec = getDiffXsecBinning(options.binFileToRoll, "reco")
            recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
            #following array is used to call function dressed2D()         
            binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
            hinput1 = dressed2D(hist1,binning,hist1.GetName()+"_rolled2D",hist1.GetName()+"_rolled2D")
            hinput2 = dressed2D(hist2,binning,hist2.GetName()+"_rolled2D",hist2.GetName()+"_rolled2D")
        else:
            print "Error: you need to pass the binning to roll 1D histograms into 2D"
            quit()

    xMin = options.xRange[0]
    xMax = options.xRange[1]
    yMin = options.yRange[0]
    yMax = options.yRange[1]

    hratio = hinput1.Clone(options.outhistname)
    hratio.SetTitle(options.histTitle)

    hsum = None
    if options.makeAsymmetry:
        hsum = hinput1.Clone("diff")
        hsum.Add(hinput2)
        hratio.Add(hinput2, -1.)
        #hratio.Divide(hsum) # do this below
        

    nbins,minx,maxx = options.h1Dbinning.split(',')
    hratioDistr = ROOT.TH1D(options.outhistname+"_1D","Distribution of ratio values",int(nbins),float(minx),float(maxx))
    #print "minx, maxx = %s, %s" % (str(minx),str(maxx))    
    nout = 0

    for ix in range(1,1+hratio.GetNbinsX()):
        for iy in range(1,1+hratio.GetNbinsY()):
            xval = hratio.GetXaxis().GetBinCenter(ix) 
            yval = hratio.GetYaxis().GetBinCenter(iy) 
            hist2xbin = hinput2.GetXaxis().FindFixBin(xval)
            hist2ybin = hinput2.GetYaxis().FindFixBin(yval)
            if xMin < xMax:
                if xval < xMin or xval > xMax: continue
            if yMin < yMax:
                if yval < yMin or yval > yMax: continue
            denval = 0
            if options.makeAsymmetry:
                denval = hsum.GetBinContent(hist2xbin, hist2ybin)
            else:
                denval = hinput2.GetBinError(hist2xbin, hist2ybin) if options.divideError else hinput2.GetBinContent(hist2xbin, hist2ybin)
            if options.divideRelativeError:
                if hinput2.GetBinContent(hist2xbin, hist2ybin) != 0:
                    denval = hinput2.GetBinError(hist2xbin, hist2ybin) / hinput2.GetBinContent(hist2xbin, hist2ybin)
                else:
                    denval = 0
            if denval != 0:                
                numval = hratio.GetBinError(ix,iy) if options.divideError else hratio.GetBinContent(ix,iy)
                if options.divideRelativeError:
                    if hratio.GetBinContent(ix,iy) != 0:
                        numval = hratio.GetBinError(ix,iy) / hratio.GetBinContent(ix,iy)
                    else:
                        numval = 0
                ratio = numval / denval
                hratioDistr.Fill(ratio)
                hratio.SetBinContent(ix,iy,ratio)
                if ratio < float(minx) or ratio > float(maxx): nout += 1
                #profX.Fill()
            else: 
                print "Warning: found division by 0 in one bin: setting ratio to " + str(options.valBadRatio)
                hratio.SetBinContent(ix,iy,options.valBadRatio)

    print "nout = " + str(nout)

    if options.xAxisTitle: hratio.GetXaxis().SetTitle(options.xAxisTitle)    
    if options.yAxisTitle: hratio.GetYaxis().SetTitle(options.yAxisTitle)    
    if options.zAxisTitle: hratio.GetZaxis().SetTitle(options.zAxisTitle)    
    xAxisTitle = hratio.GetXaxis().GetTitle()
    yAxisTitle = hratio.GetYaxis().GetTitle()
    zAxisTitle = hratio.GetZaxis().GetTitle()

    if xMax > xMin and not "::" in xAxisTitle: 
        xAxisTitle = xAxisTitle + "::" + str(options.xRange[0]) + "," + str(options.xRange[1])
    if yMax > yMin and not "::" in yAxisTitle: 
        yAxisTitle = yAxisTitle + "::" + str(options.yRange[0]) + "," + str(options.yRange[1])

    # print "xAxisTitle = " + xAxisTitle
    # print "yAxisTitle = " + yAxisTitle
    # print "zAxisTitle = " + zAxisTitle

    adjustSettings_CMS_lumi()

    canvas2D = ROOT.TCanvas("canvas2D","",700,700)
    
    # the axis name can be used to set the range if it is in the format "name::min,maz"
    # if this is not already the case, use the selected range from the input option
    if not "::" in zAxisTitle:  
        zAxisTitle = zAxisTitle + "::" + str(options.ratioRange[0]) + "," + str(options.ratioRange[1])
    drawCorrelationPlot(hratio,xAxisTitle,yAxisTitle,zAxisTitle,
                        options.outhistname,"ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D)
    
    canvas = ROOT.TCanvas("canvas","",800,700)
    drawTH1(hratioDistr, 
            hratio.GetZaxis().GetTitle() if options.zAxisTitle else "ratio",
            "number of events",
            outname,
            "ratioDistribution",
            options.outhistname,
            passCanvas=canvas
            )

    # making distribution of pulls
    if options.makePulls:
        hpull = ROOT.TH1D("pulls","Distribution of pulls",100,-5,5)
        for ix in range(1,1+hinput1.GetNbinsX()):
            for iy in range(1,1+hinput1.GetNbinsY()):
                err = math.sqrt(pow(hinput1.GetBinError(ix,iy),2) + pow(hinput2.GetBinError(ix,iy),2))
                err *= options.pullErrorScaleFactor
                pull = hinput1.GetBinContent(ix,iy) - hinput2.GetBinContent(ix,iy)
                hpull.Fill(pull/err)

        drawTH1(hpull, 
                "pulls",
                "number of events",
                outname,
                "pullsDistribution",
                options.outhistname,
                passCanvas=canvas
                )
 
    ###########################
    # Now save things
    ###########################
    if not options.outfilename.endswith(".root"):
        options.outfilename = options.outfilename + ".root"
    tf = ROOT.TFile.Open(outname+options.outfilename,'recreate')
    hratio.Write(options.outhistname)
    tf.Close()
    print ""
    print "Created file %s" % (outname+options.outfilename)
    print ""

                               
         
