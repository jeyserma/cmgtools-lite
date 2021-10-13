
import ROOT

DY = {

    'tag':      "DY",
    'name':     "",
    'xsec':     6077.22,
    'color':    ROOT.kAzure+2,
    'sumw':     11828591.0, #311461009200.830078,

}

DY_MiNNLO = {

    'tag':      "DY_MiNNLO",
    'name':     "",
    'xsec':     6077.22/3.0,
    'color':    ROOT.kAzure+2,
    'sumw':     1.0,

}




TT0l = {

    'tag':      "TT0l",
    'name':     "",
    'xsec':     380.10,
    'color':    ROOT.kAzure+2,
    'sumw':     11597594,
}

TT1l = {

    'tag':      "TT1l",
    'name':     "",
    'xsec':     364.35,
    'color':    ROOT.kAzure+2,
    'sumw':     18048758,
}

TT2l = {

    'tag':      "TT2l",
    'name':     "",
    'xsec':     87.31,
    'color':    ROOT.kAzure+2,
    'sumw':     4959518.0,
}



WW = {

    'tag':      "WW",
    'name':     "",
    'xsec':     12.178,
    'color':    ROOT.kAzure+2,
    'sumw':     1.0,
}

WZ = {

    'tag':      "WZ",
    'name':     "",
    'xsec':     4.42965,
    'color':    ROOT.kAzure+2,
    'sumw':     1.0,
}

ZZ = {

    'tag':      "ZZ",
    'name':     "",
    'xsec':     0.564,
    'color':    ROOT.kAzure+2,
    'sumw':     1.0,
}

WJets0J = {

    'tag':      "WJets0J",
    'name':     "",
    'xsec':     49155.,
    'color':    ROOT.kAzure+2,
    'sumw':     1.0,
}

WJets1J = {

    'tag':      "WJets1J",
    'name':     "",
    'xsec':     8047.00,
    'color':    ROOT.kAzure+2,
    'sumw':     1.0,
}

WJets2J = {

    'tag':      "WJets2J",
    'name':     "",
    'xsec':     3160.00,
    'color':    ROOT.kAzure+2,
    'sumw':     1.0,
}


SingleMuon = {

    'tag':      "singlemuon",
    'name':     "Data",
    'xsec':     1,
    'color':    -1,

}

'''
data_sm:            .*;     SubPath="SingleMuon/",  FillColor=ROOT.kBlue, Label="Data (single Muon)"
data_eg:            .*;     SubPath="HighEGJet/",   FillColor=ROOT.kBlue, Label="Data (EG)"

TT:         .* :    87.31;      FillColor=ROOT.kGreen+2,    Label="t quark", SubPath="TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/"
TT:         .* :    380.10;     FillColor=ROOT.kGreen+2,    Label="t quark", SubPath="TTToHadronic_TuneCP5_13TeV-powheg-pythia8/"
TT:         .* :    364.35;     FillColor=ROOT.kGreen+2,    Label="t quark", SubPath="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/" 

DY:         .* :    6077.22;    FillColor=ROOT.kAzure+2,    Label="DY", SubPath="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/"                   
                     
WJets:      .* :    49155   ;   FillColor=ROOT.kRed+1,      Label="WJets",  SubPath="WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/"
WJets:      .* :    8047.00 ;   FillColor=ROOT.kRed+1,      Label="WJets",  SubPath="WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/"
WJets:      .* :    3160.00 ;   FillColor=ROOT.kRed+1,      Label="WJets",  SubPath="WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/"

Diboson:    .* :    12.178;     FillColor=ROOT.kViolet,     Label="Diboson", SubPath="WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/"
Diboson:    .* :    4.42965;    FillColor=ROOT.kViolet,     Label="Diboson", SubPath="WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/"
Diboson:    .* :    0.564;      FillColor=ROOT.kViolet,     Label="Diboson", SubPath="ZZ_TuneCP5_13TeV-pythia8/"  

'''