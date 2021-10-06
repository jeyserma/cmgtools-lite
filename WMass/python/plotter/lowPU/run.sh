#!/bin/bash

eosdir=/eos/user/${USER::1}/${USER}/www/wmass/MiNNLOCorr
if [ ! -d $eosdir ]; then
    mkdir -p $eosdir
fi

runData=1
if [ $# -gt 0 ]; then
    runData=$1
fi

if [ $runData -gt 0 ]; then
    python mcPlots.py -f -l 0.199269742 --lumi-weight 0.199269742 lowPU/cfg/mca.txt lowPU/cfg/cuts_Zmumu.txt lowPU/cfg/plots_Zmumu.txt --noCms -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ --sP "muon_pt,m_mumu,muon_eta,muon_phi,muon_leading_pt,muon_subleading_pt,pt_mumu,y_mumu" --pdir $eosdir -p "data,DY" --pg "data := data_sm"  --legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/rdfDefine.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --json lowPU/lumiPileup/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt --showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --ratioNums 'DY' --ratioDen "data" 
else
    python mcPlots.py -f -l 0.199269742 --lumi-weight 0.199269742 lowPU/cfg/mca.txt lowPU/cfg/cuts_Zmumu.txt lowPU/cfg/plots_Zmumu.txt --noCms -P /eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/ --sP "muon_pt,m_mumu,muon_eta,muon_phi,muon_leading_pt,muon_subleading_pt,pt_mumu,y_mumu" --pdir $eosdir -p "DY,DY_N3LLCorr" --legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree --rdf-define-file lowPU/cfg/rdfDefine.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --json lowPU/lumiPileup/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU_suppressedHighPULS.txt -X mtl1pf40 --showRatio --maxRatioRange 0.9 1.3 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "Corr/Uncorr" -X mtl1pf40  --ratioNums 'DY_N3LLCorr' --ratioDen "DY" --plotmode nostack
fi




#--showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --plotmode nostack --ratioNums 'Zmumu_postVFP' --ratioDen "data"


#python mcPlots.py -f -l 16.8 w-mass-13TeV/testingNano/cfg/mca-wlike.txt w-mass-13TeV/testingNano/cfg/test/cuts_wlike.txt w-mass-13TeV/testingNano/cfg/plots_test.txt --noCms -P /data/shared/originalNANO/ --sP "muon_pt,muon_eta_fine,muon_pt_other,muon_eta_other_fine,zmass,zpt,zy,mt_wlike_MET,MET_pt,met_wlike_MET"   -W "puw_2016UL_era_OLD(Pileup_nTrueInt,eraVFP)*_get_fullMuonSF(Muon_pt[goodMuonsCharge][0],Muon_eta[goodMuonsCharge][0],Muon_charge[goodMuonsCharge][0],Muon_pt[goodMuonsOther][0],Muon_eta[goodMuonsOther][0],eraVFP)*_get_MuonPrefiringSF(Muon_eta,Muon_pt,Muon_looseId,eraVFP)" --pdir plots/testNanoAOD/testZ/allEvents_trigPlus_dataVSmc_noMT_postVFP/ -p "data,Zmumu_postVFP" --pg "data := data_postVFP" --legendFontSize 0.042 --allProcInLegend --n-column-legend 2 --setLegendCoordinates 0.2,0.76,0.9,0.92  -j 8  --nanoaod-tree --rdf-define-file w-mass-13TeV/testingNano/cfg/test/rdfDefine_wlike.txt --rdf-alias "goodMuonsCharge: goodMuonsPlus:.*" --rdf-alias "goodMuonsOther: goodMuonsMinus:.*" -v 3 --showRatio --maxRatioRange 0.9 1.1 --fixRatioRange  --onlyStatErrorOnRatio --noLegendRatioPlot --ratioYLabel "MC/Data" -X mtl1pf40  --plotmode nostack --ratioNums 'Zmumu_postVFP' --ratioDen "data"


