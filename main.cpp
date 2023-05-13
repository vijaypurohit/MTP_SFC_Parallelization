/*!
 * @author Vijay Purohit
 * @date Created by vijay on 22-02-2023.
 * @name MTP Thesis - Delay Aware Optimal Parallelization of Service Function Chains.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <random> //to generate random number
#include <functional>  // lambda function

using namespace std;
using type_delay = double; ///< data type of delay float /double
using type_wgt = double; ///< data type of distance

int debug = 1;
/**************** Custom Header Files Declaration ****************/
#include "Variables.h"
#include "PhysicalGraph.h" ///< graph structure for physical network, physical node and edge and node capacity
#include "VirtualNetworkFunctions.h" ///< structure for VNFnode, VNFs,
#include "ServiceFunctionChain.h" /// structure for SFC
#include "FileFunctions.h" ///< File reading writing functions
#include "HelperStructures.h" ///< Some of the helper structers
#include "SimulationClass.h" ///< for simulation purpose
#include "DelayCalculationFunctions.h" ///< Various Delay Calculation functions
#include "Algorithms.h" ///< algorithms implemented
#include "Tests.h" ///< Simulation Test implemented

/*!
 * init
 */
bool init(){ //init
    /// Run this if SFC Len is going to be more than 10
    if(maxSFCLength > 10){
        integerCompositionsEnumeration(int(integerCompositions.size()),maxSFCLength);
        find_all_nCk(int(nCk.size()), maxSFCLength);
    }

    std::cout<<"\nInput Directories:  ";
    for(const auto & di : std::filesystem::directory_iterator(input_directory)){  if(di.is_directory()) cout<<di.path().stem()<<",  ";  }

    const vector<std::string> delay_settings = {"PTF", "PFT", "TPF", "TFP", "FTP", "FPT"};
    std::cout<<"\n Delay Settings: ";
    for(int id=0; id<delay_settings.size(); id++){  cout<<"("<<id<<": "<<delay_settings[id]<<") "; }

    const vector<std::string> allSimulationName = {"Basic", "VaryingDeploymentParameter", "ImpactOfChainLength", "ImpactOfFixedChainLength", "ImpactOfNumberOfInstancesPerServer", "ImpactOfHeterogeneousDelays", "ImpactOfVariousDelays"};
    std::cout<<"\n\n id | SimulationName";
    for(int id=0; id<allSimulationName.size(); id++){  cout<<"\n  "<<id<<": "<<allSimulationName[id]; }


    string dirname;  cout<<"\n\nEnter Directory Name (network file should be present as \"network_{DirectoryName}.txt\"): "; cin>>dirname;
    int asid=-1;            cout<<"\nEnter Simulation Id: "; cin>>asid;
    int setting_id = -1;    cout<<"\nEnter Delay Settings Id: "; cin>>setting_id;
    pair<type_delay, type_delay> funExeTimeRange;

    if(asid<0 or asid>allSimulationName.size() or setting_id<0 or setting_id>delay_settings.size() or dirname.empty()){
        cerr<<"\n Invalid inputs.";
        return false;
    }

    if(setDelayParameterSettings(setting_id, funExeTimeRange) == false)
        return false;
//    ifRejectSFC = true;
//    speedOfLight = 300000; velocityFactor = 1.0;
    const string& networkFileName = "network_"+dirname+".txt";
    const string& fileVNF = "VNF"+to_string(maxNumVNFs)+".txt";
    GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, funExeTimeRange, dirname+"/" , fileVNF);

    switch (asid) {
        case 0: //Basic
            ifRejectSFC = true;
            return SimulationTest_Basic(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF, "SFC10_L_R_B.txt", {60,0});
            break;
        case 1: //VaryingDeploymentParameter
            return SimulationTest_VaryingNetworkVNF_DeploymentParameter(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF, "SFC10_L_R_VDP.txt");
            break;
        case 2:  //ImpactOfChainLength
                SimulationTest_ImpactOfChainLength(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF);
        case 3: //ImpactOfFixedChainLength
            return SimulationTest_FixedChainLength(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF);
            break;
        case 4: //ImpactOfNumberOfInstancesPerServer
            return SimulationTest_ImpactOfNumberOfInstancesPerServer(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF);
            break;
        case 5: //ImpactOfHeterogeneousDelays
            return SimulationTest_ImpactOfHeterogeneousDelays(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF);
            break;
        default:
            return SimulationTest_ImpactOfVariousDelays(allSimulationName[asid], dirname, networkFileName);
            return false;
    }
  return false;
}//init

bool AutomateForEachDirectory(){

    for(const string& dirname : {"NewYork"}){

        const string& networkFileName = "network_"+dirname+".txt";

        cout<<"\n"<<"ImpactOfVariousDelays_setid_"<<(0);
        SimulationTest_ImpactOfVariousDelays("ImpactOfVariousDelays_setid_", dirname, networkFileName);

        for(const int setting_id: {0,1,2,3,4,5}){
            const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"ihd.txt";
            pair<type_delay, type_delay> funExeTimeRange;
            setDelayParameterSettings(setting_id, funExeTimeRange);
            GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, funExeTimeRange, dirname+"/" , fileVNF);
            cout<<"\n"<<"ImpactOfHeterogeneousDelays_setid_"<<(setting_id);
            SimulationTest_ImpactOfHeterogeneousDelays("ImpactOfHeterogeneousDelays_setid_"+to_string(setting_id), dirname, networkFileName, fileVNF, {55,1});
            remove((input_directory+dirname+"/"+fileVNF).c_str());
        }
        for(const int setting_id: {0,1,2,3,4,5}){
            const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"nips.txt";
            pair<type_delay, type_delay> funExeTimeRange;
            setDelayParameterSettings(setting_id, funExeTimeRange);
            GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, funExeTimeRange, dirname+"/" , fileVNF);
            cout<<"\n"<<"ImpactOfNumberOfInstancesPerServer_setid_"<<(setting_id);
            SimulationTest_ImpactOfNumberOfInstancesPerServer("ImpactOfNumberOfInstancesPerServer_setid_"+to_string(setting_id), dirname, networkFileName, fileVNF, {55,1});
            remove((input_directory+dirname+"/"+fileVNF).c_str());
        }

        if(dirname == "NewYork"){
            cout<<"\n"<<"SameChainsHeterogeneousDelaysImpactOfNumberOfInstances"<<(0);
            SimulationTest_SameChainsHeterogeneousDelaysImpactOfNumberOfInstancesPerServer(dirname);

            const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"len.txt";
            pair<type_delay, type_delay> funExeTimeRange;
            setDelayParameterSettings(6, funExeTimeRange);
            GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, funExeTimeRange, dirname+"/" , fileVNF);
            SimulationTest_ImpactOfChainLength("ImpactOfChainLength_setid_"+to_string(6), dirname, networkFileName, fileVNF);
            SimulationTest_FixedChainLength("ImpactOfFixedChainLength_setid_"+to_string(6), dirname, networkFileName, fileVNF);
            remove((input_directory+dirname+"/"+fileVNF).c_str());
        }
    } //for directory listing

    return true;
}

int main()
{
    std::cout<<"\n~~~~~~~~~~~~~~~~~~~ [ Delay-Aware Optimal Parallelization of SFC ] ~~~~~~~~~~~~~~~~~~~~~~~~";


int id; cout<<"\n"<<"1. For Individual Testing. \n2. Automate for directories."; cin>>id;
if(id == 1){
    if(init() == false){
        return 1;
    }
}else if(id == 2){
    AutomateForEachDirectory();
}

    cout<<"\n\nsuccess. press any key to exit";
    cout << "\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    return 0;
}


