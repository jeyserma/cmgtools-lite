#include "TSystem.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <cstdlib> //as stdlib.h                    
#include <cstdio>

#include "../interface/utility.h"

using namespace std;

// root -l -b -q makeLoopTree.C++

void makeLoopTree(const bool isMuon = false, 
		  const string& outfileName = "wmass_varhists.root",
		  const string& usePreFSRvar = "false"  // if "false", use DressedLepton. Pass it as a string to build command below ("1" or "0" work as well)
		  ) 

{

  // load source codes with ++, so that they are always compiled (you never know ...)

  string cmssw_base = getEnvVariable("CMSSW_BASE");
  cout << "CMSSW_BASE = " << cmssw_base << endl;

  string host_name = getEnvVariable("HOSTNAME");
  cout << "HOSTNAME = " << host_name << endl;

 
  cout << "Loading functions.cc" << endl;
  gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/functions.cc+",cmssw_base.c_str())); 
  cout << "Loading functionsWMass.cc" << endl;
  gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/functionsWMass.cc+",cmssw_base.c_str()));
  cout << "Loading loopNtuplesSkeleton.cc" << endl;
  gROOT->ProcessLine(".L loopNtuplesSkeleton.C++");

  //string command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_WSKIM_NEW/\",\"./\",\"" + outfileName + "\")";
  string command = "";
  if (isMuon) {

    if (host_name.find("lxplus") != string::npos) 
      command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_13June2019/\",\"./\",\"" + outfileName + "\",true,"+ usePreFSRvar + ")";

    else if (host_name.find("pccmsrm") != string::npos) 
      command = "loopNtuplesSkeleton(\"/u1/mciprian/trees/muon/TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_13June2019/\",\"./\",\"" + outfileName + "\",true,"+ usePreFSRvar + ")";

  } else {

    if (host_name.find("lxplus") != string::npos) 
      command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/mciprian/TREES_WELE_19Feb2019/\",\"./\",\"" + outfileName + "\",false,"+ usePreFSRvar + ")";
    else if (host_name.find("pccmsrm") != string::npos) 
      command = "loopNtuplesSkeleton(\"/u1/mciprian/trees/electron/TREES_WELE_19Feb2019/\",\"./\",\"" + outfileName + "\",false,"+ usePreFSRvar + ")";
  }

  //string command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_electrons_1l_V6_TINY/\",\"./\",\"" + outfileName + "\")";
  cout << "Executing " << command << endl;
  gROOT->ProcessLine(command.c_str());
  cout << endl;                                          
  cout << "===========================" << endl;                    
  cout << " THE END!" << endl;                          
  cout << "===========================" << endl;         



}
