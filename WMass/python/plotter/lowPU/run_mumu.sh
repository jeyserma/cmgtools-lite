#!/bin/bash

TAG="test"


## data

#python mcPlots.py -f lowPU/cfg/mca.txt lowPU/cfg/zmumu_data_cuts.txt lowPU/cfg/zmumu_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/  -p "data" --pg "data := data_sm" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/zmumu_data_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --json lowPU/lumiPileup/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt  --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//data_Zmumu.root --lumi-weight 1.0

python mcPlots.py -f lowPU/cfg/mca.txt lowPU/cfg/zmumu_mc_cuts.txt lowPU/cfg/zmumu_plots.txt -P /data/shared/lowPU/Nano_0302/ -p "DY_MiNNLO" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/zmumu_mc_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//DY_Zmumu.root --lumi-weight 1.0 -W "lepSFNew*prefireSFPAT*scetlibWeightsZscalar"

# DY,

#python mcPlots.py -f lowPU/cfg/mca.txt lowPU/cfg/zmumu_mc_cuts.txt lowPU/cfg/zmumu_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ -p "TT2l,TT1l,TT0l" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/zmumu_mc_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//TT_Zmumu.root --lumi-weight 1.0 -W "lepSF*prefireSFPAT"

#python mcPlots.py -f lowPU/cfg/mca.txt lowPU/cfg/zmumu_mc_cuts.txt lowPU/cfg/zmumu_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ -p "WJets0J,WJets1J,WJets2J" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/zmumu_mc_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//WJets_Zmumu.root --lumi-weight 1.0 -W "lepSF*prefireSFPAT"

#python mcPlots.py -f lowPU/cfg/mca.txt lowPU/cfg/zmumu_mc_cuts.txt lowPU/cfg/zmumu_plots.txt -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ -p "WW,WZ,ZZ" -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/zmumu_mc_defines.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --skipPlot --out /eos/user/j/jaeyserm/www/wmass/test//VV_Zmumu.root --lumi-weight 1.0 -W "lepSF*prefireSFPAT"

#-W "lepSF"

#-W "applyMuonSF(Muon_pt_corr[goodMuons], Muon_eta[goodMuons], Muon_charge[goodMuons])"
#-W "lowPU_pu_reweight(Pileup_nPU)"

#applyMuonSF(Vec_f pt, Vec_f eta, Vec_i q) 

#--unweight

rm /eos/user/j/jaeyserm/www/wmass/test/Zmumu.root
#hadd /eos/user/j/jaeyserm/www/wmass/test/Zmumu.root  /eos/user/j/jaeyserm/www/wmass/test/DY_Zmumu.root /eos/user/j/jaeyserm/www/wmass/test/data_Zmumu.root

hadd /eos/user/j/jaeyserm/www/wmass/test/Zmumu.root  /eos/user/j/jaeyserm/www/wmass/test/DY_Zmumu.root /eos/user/j/jaeyserm/www/wmass/test/TT_Zmumu.root /eos/user/j/jaeyserm/www/wmass/test/VV_Zmumu.root /eos/user/j/jaeyserm/www/wmass/test/WJets_Zmumu.root /eos/user/j/jaeyserm/www/wmass/test/data_Zmumu.root

## MC


#python mcPlots.py -f -l 0.199269742 --lumi-weight 0.199269742 lowPU/cfg/mca.txt lowPU/cfg/cuts_Zmumu.txt lowPU/cfg/plots_Zmumu.txt --noCms -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ --sP "m_mumu" --pdir /eos/user/j/jaeyserm/www/wmass/test/ -p "data,DY" --pg "data := data_sm"  --legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/rdfDefine.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --json lowPU/lumiPileup/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt --showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --ratioNums 'DY' --ratioDen "data" 

# muon_pt,m_mumu,muon_eta,muon_phi,muon_leading_pt,muon_subleading_pt,pt_mumu


#--showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --plotmode nostack --ratioNums 'Zmumu_postVFP' --ratioDen "data"


#python mcPlots.py -f -l 16.8 w-mass-13TeV/testingNano/cfg/mca-wlike.txt w-mass-13TeV/testingNano/cfg/test/cuts_wlike.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP "muon_pt,muon_eta_fine,muon_pt_other,muon_eta_other_fine,zmass,zpt,zy,mt_wlike_MET,MET_pt,met_wlike_MET"   -W "puw_2016UL_era_OLD(Pileup_nTrueInt,eraVFP)*_get_fullMuonSF(Muon_pt[goodMuonsCharge][0],Muon_eta[goodMuonsCharge][0],Muon_charge[goodMuonsCharge][0],Muon_pt[goodMuonsOther][0],Muon_eta[goodMuonsOther][0],eraVFP)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)" --pdir plots/testNanoAOD/testZ/allEvents_trigPlus_dataVSmc_noMT_postVFP/ -p "data,Zmumu_postVFP" --pg "data := data_postVFP" --legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_wlike.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --plotmode nostack --ratioNums 'Zmumu_postVFP' --ratioDen "data"


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