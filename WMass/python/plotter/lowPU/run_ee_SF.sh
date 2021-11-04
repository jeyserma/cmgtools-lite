#!/bin/bash

TAG="test"


## data

#python mcPlots.py -f lowPU/cfg_SF/mca.txt lowPU/cfg_SF/zee_data_cuts.txt lowPU/cfg_SF/zee_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/  -p "data" --pg "data := data_eg" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg_SF/zee_data_defines.txt -v 3 --json lowPU/lumiPileup/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt  --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//data_Zee.root --lumi-weight 1.0

python mcPlots.py -f lowPU/cfg_SF/mca.txt lowPU/cfg_SF/zee_mc_cuts.txt lowPU/cfg_SF/zee_plots.txt -P /data/shared/lowPU/Nano_0302/ -p "DY_MiNNLO_ee" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg_SF/zee_mc_defines.txt -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//DY_Zee.root --lumi-weight 1.0 -W "lepSF*prefireSFPAT"

# -W "lepSF"  *

#python mcPlots.py -f lowPU/cfg_SF/mca.txt lowPU/cfg_SF/zee_mc_cuts.txt lowPU/cfg_SF/zee_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ -p "TT2l,TT1l,TT0l" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg_SF/zee_mc_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//TT_Zee.root --lumi-weight 1.0 -W "lepSF"

#python mcPlots.py -f lowPU/cfg_SF/mca.txt lowPU/cfg_SF/zee_mc_cuts.txt lowPU/cfg_SF/zee_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ -p "WJets0J,WJets1J,WJets2J" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg_SF/zee_mc_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//WJets_Zee.root --lumi-weight 1.0 -W "lepSF"

#python mcPlots.py -f lowPU/cfg_SF/mca.txt lowPU/cfg_SF/zee_mc_cuts.txt lowPU/cfg_SF/zee_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ -p "WW,WZ,ZZ" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg_SF/zee_mc_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//VV_Zee.root --lumi-weight 1.0 -W "lepSF"

#-W "lepSF"

#-W "applyMuonSF(Muon_pt_corr[goodMuons], Muon_eta[goodMuons], Muon_charge[goodMuons])"
#-W "lowPU_pu_reweight(Pileup_nPU)"

#applyMuonSF(Vec_f pt, Vec_f eta, Vec_i q) 

#--unweight

rm /eos/user/j/jaeyserm/www/wmass/test/Zee.root
hadd /eos/user/j/jaeyserm/www/wmass/test/Zee.root  /eos/user/j/jaeyserm/www/wmass/test/DY_Zee.root /eos/user/j/jaeyserm/www/wmass/test/data_Zee.root
#hadd /eos/user/j/jaeyserm/www/wmass/test/Zee.root  /eos/user/j/jaeyserm/www/wmass/test/DY_Zee.root /eos/user/j/jaeyserm/www/wmass/test/TT_Zee.root /eos/user/j/jaeyserm/www/wmass/test/VV_Zee.root /eos/user/j/jaeyserm/www/wmass/test/WJets_Zee.root /eos/user/j/jaeyserm/www/wmass/test/data_Zee.root

## MC


#python mcPlots.py -f -l 0.199269742 --lumi-weight 0.199269742 lowPU/cfg_SF/mca.txt lowPU/cfg_SF/cuts_Zee.txt lowPU/cfg_SF/plots_Zee.txt --noCms -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ --sP "m_ee" --pdir /eos/user/j/jaeyserm/www/wmass/test/ -p "data,DY" --pg "data := data_sm"  --legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg_SF/rdfDefine.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --json lowPU/lumiPileup/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt --showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --ratioNums 'DY' --ratioDen "data" 

# muon_pt,m_ee,muon_eta,muon_phi,muon_leading_pt,muon_subleading_pt,pt_ee


#--showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --plotmode nostack --ratioNums 'Zee_postVFP' --ratioDen "data"


#python mcPlots.py -f -l 16.8 w-mass-13TeV/testingNano/cfg_SF/mca-wlike.txt w-mass-13TeV/testingNano/cfg_SF/test/cuts_wlike.txt w-mass-13TeV/testingNano/cfg_SF/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP "muon_pt,muon_eta_fine,muon_pt_other,muon_eta_other_fine,zmass,zpt,zy,mt_wlike_MET,MET_pt,met_wlike_MET"   -W "puw_2016UL_era_OLD(Pileup_nTrueInt,eraVFP)*_get_fullMuonSF(Muon_pt[goodMuonsCharge][0],Muon_eta[goodMuonsCharge][0],Muon_charge[goodMuonsCharge][0],Muon_pt[goodMuonsOther][0],Muon_eta[goodMuonsOther][0],eraVFP)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)" --pdir plots/testNanoAOD/testZ/allEvents_trigPlus_dataVSmc_noMT_postVFP/ -p "data,Zee_postVFP" --pg "data := data_postVFP" --legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree --rdf-define-file w-mass-13TeV/testingNano/cfg_SF/test/rdfDefine_wlike.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --plotmode nostack --ratioNums 'Zee_postVFP' --ratioDen "data"


# 
#puw_2016UL_era_OLD(Pileup_nTrueInt,eraVFP)
#_get_fullMuonSF(
#    Muon_pt[goodMuonsCharge][0],
#    Muon_eta[goodMuonsCharge][0],
#    Muon_charge[goodMuonsCharge][0],
#    Muon_pt[goodMuonsOther][0],
#    Muon_eta[goodMuonsOther][0],
#    raVFP)
#_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)