
import ROOT


fIn = ROOT.TFile("/eos/cms/store/cmst3/group/wmass/LowPU/Nano_0302/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NanoAOD_0302/210303_060423/0000/NanoAOD_MC_1.root")


fIn.ls()

tree = fIn.Get("Events")

tree.Print()
