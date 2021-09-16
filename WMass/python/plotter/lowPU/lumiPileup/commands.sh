#!/bin/bash


######### lumi #########

# install brilcalc
# https://github.com/piedraj/instructions/blob/master/BRILCALC.md


#brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_299996-303819_13TeV_PromptReco_Collisions17_JSON_LowPU.txt

#brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU.txt



######### pileup #########



pileupCalc.py -i ~/bril/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 10 --numPileupBins 10 pileup.root


pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU.txt   --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 10 --numPileupBins 10 pileup_full.root