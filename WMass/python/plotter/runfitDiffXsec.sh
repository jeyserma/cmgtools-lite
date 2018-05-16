#!/bin/bash                                                                                                                                                                  

wmassDir="/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/"
plotterDir="${wmassDir}python/plotter/"
cardDir="${plotterDir}testDatacardWriter/"
combineDir="${plotterDir}combineDir"
#wwwDir="/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/"

nVerbose=0
nTotBins=570
justPrint="n"
log="p"  # p, t, s for print on stdout (default), throw away (to /dev/null) or save to a file
saveOnlyhiggsCombineFile="y"

dirToRunCombineFrom="${plotterDir}diffXsecFit_testScalLikelihood_freezeShapeNuis/"
#combineOptions=" --saveFitResult "
combineOptions=" --algo grid  --points 20 --keepFailures --setParameterRanges \"r=0.9,1.1\" --cminDefaultMinimizerType GSLMultiMin --cminDefaultMinimizerAlgo BFGS2  --freezeNuisanceGroups pdfs,scales,alphaS,wpt,frshape "
#selectedbins=`seq 1 $nTotBins`
selectedbins=(20 35 50 65 80 95 125 140 155 170 185 200 215 230 245 260 275 290 305 320 335 350 365 380 395 410 425 440 510 520)
#charges=("plus" "minus")
charges=("plus")

cd ${plotterDir}
eval `scramv1 runtime -sh`
echo "${CMSSW_VERSION} -->  ${PWD}"

cd ${combineDir}
eval `scramv1 runtime -sh`
echo "${CMSSW_VERSION} -->  ${PWD}"

mkdir -p ${dirToRunCombineFrom}
cd ${dirToRunCombineFrom}
echo "${CMSSW_VERSION} -->  ${PWD}"

if [[ "${saveOnlyhiggsCombineFile}" != "y" ]]; then
    combineOptions="${combineOptions=} --saveFitResult"
fi        

for charge in "${charges[@]}"
do

    echo "#########################################"

    for nbin in "${selectedbins[@]}"
    do
	
	echo "Charge: ${charge}, bin = ${nbin}"

	outTagName="_W${charge}_bin${nbin}"

	outlog=""
	if [[ "${log}" == "t" ]]; then
	    outlog=" &> /dev/null "
	elif [[ "${log}" == "s" ]]; then
	    outlog=" &> ${dirToRunCombineFrom}combine${outTagName}.log "
	fi

	runCombine="combine ${cardDir}Wel_${charge}_shapes_addInclW_card_bin${nbin}.txt -M MultiDimFit -t -1 --expectSignal=1 ${combineOptions} -n ${outTagName} -v ${nVerbose} ${outlog}"

	fitresName="multidimfit${outTagName}.root"

	echo "${runCombine}"
	if [[ "${justPrint}" != "y" ]]; then
            echo "${runCombine}" | bash

	    if [[ "${saveOnlyhiggsCombineFile}" != "y" ]]; then
		if [ -e "${dirToRunCombineFrom}${fitresName}" ]; then
		    echo "Fit for charge = ${charge} and bin = ${nbin} done!"
		else
                    echo "======================================================"
                    echo "======================================================"
                    echo "===  WARNING! Something wrong with the fit for charge = ${charge} and bin = ${nbin}."
                    echo "===  File ${fitresName} was not created! Skipping this configuration."
                    echo "======================================================"
                    echo "======================================================"
                    continue
		fi
	    else
		echo "Fit for charge = ${charge} and bin = ${nbin} done!"
	    fi

	fi
	
	echo "#########################################"

    done

done