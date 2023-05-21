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

    const vector<std::string> allSimulationName = {"Basic", "VaryingDeploymentParameter", "ImpactOfChainLength", "ImpactOfFixedChainLength", "ImpactOfNumberOfInstancesPerServer", "ImpactOfHeterogeneousDelays", "ImpactOfVariousDelays", "ComparisonWithExistingPPC"};
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
            SimulationTest_VaryingNetworkVNF_DeploymentParameterComaprisonsBetweenSFCsortting(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF, parallelPairsOpt);
//            SimulationTest_VaryingNetworkVNF_DeploymentParameter(allSimulationName[asid]+"_"+delay_settings[setting_id]+"asc_length", dirname, networkFileName, fileVNF);
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
        case 6: //ImpactOfHeterogeneousDelays
            return SimulationTest_ImpactOfVariousDelays(allSimulationName[asid], dirname, networkFileName);
            break;
        case 7: //ImpactOfHeterogeneousDelays
            return ComparisonWithExistingPPC(allSimulationName[asid]+"_"+delay_settings[setting_id], dirname, networkFileName, fileVNF);
            break;
        default:
            return false;
    }
  return false;
}//init

bool AutomateForEachDirectory(){

    for(const string& dirname : {"NewYork", "India", "Germany", "TA2"}){ //

        const string& networkFileName = "network_"+dirname+".txt";

        cout<<"\n"<<"ImpactOfVariousDelays:";
        SimulationTest_ImpactOfVariousDelays("ImpactOfVariousDelays_setid_", dirname, networkFileName);

        for(const int setting_id: {0,1,2,3,4,5}){
            const string& fileVNF = "VNF"+to_string(maxNumVNFs)+".txt";
            pair<type_delay, type_delay> funExeTimeRange;
            setDelayParameterSettings(setting_id, funExeTimeRange);
            GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, funExeTimeRange, dirname+"/" , fileVNF);

            cout<<"\n"<<"ImpactOfHeterogeneousDelays_setid_"<<(setting_id);
            SimulationTest_ImpactOfHeterogeneousDelays("ImpactOfHeterogeneousDelays_setid_"+to_string(setting_id), dirname, networkFileName, fileVNF, {55,0});

            cout<<"\n"<<"ImpactOfNumberOfInstancesPerServer_setid_"<<(setting_id);
            SimulationTest_ImpactOfNumberOfInstancesPerServer("ImpactOfNumberOfInstancesPerServer_setid_"+to_string(setting_id), dirname, networkFileName, fileVNF,{55,0});
        }

        if(dirname == "NewYork"){
            cout<<"\n"<<"SameChainsHeterogeneousDelaysImpactOfNumberOfInstances"<<(0);
            SimulationTest_SameChainsHeterogeneousDelaysImpactOfNumberOfInstancesPerServer(dirname);
            for(const int set: {3,5}){

                const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"len.txt";
                pair<type_delay, type_delay> funExeTimeRange;
                setDelayParameterSettings(set, funExeTimeRange);
                GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, funExeTimeRange, dirname+"/" , fileVNF);
                SimulationTest_ImpactOfChainLength("ImpactOfChainLength_setid_"+to_string(set), dirname, networkFileName, fileVNF, {55,1});
                SimulationTest_FixedChainLength("ImpactOfFixedChainLength_setid_"+to_string(set), dirname, networkFileName, fileVNF, {55,1});
            }
        }

    } //for directory listing

    return true;
}

int main()
{
    std::cout<<"\n~~~~~~~~~~~~~~~~~~~ [ Delay-Aware Optimal Parallelization of SFC ] ~~~~~~~~~~~~~~~~~~~~~~~~";

//    string dirname = "sample0";
//    Simulations simtest("forReport", dirname);
//    const string& fileNetwork = "network_"+dirname+".txt";
//    const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"_ex.txt";
//    const string& fileSFC = "SFC10_L_R_ex.txt";
//    pair<type_delay, type_delay> funExeTimeRange;
//
//    setDelayParameterSettings(0, funExeTimeRange);
//    readNetwork(simtest.fullDirName,fileNetwork, simtest.PhysicalNetwork,{"fixed-all", 2});
//    GenerateRandomVNFs(4, funServiceRateRange, funExeTimeRange, simtest.fullDirName , fileVNF);
//
//    readVirtualNetworkFunctions(simtest.fullDirName, fileVNF, simtest.VNFNetwork);
//    simtest.findRandomParallelPairs(55,0); /// based on #VNFs
//
//    GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, 20, sfcArrivalRateRange, {true, 4}, simtest.fullDirName, fileSFC);
//    simtest.readGenericServiceFunctionsChains(fileSFC, "asc_length");
//
//    simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
//    simtest.DeploymentVNF_ScoreMethod(0.5,0.75,1,1);
//
//    simtest.showPNsDescription();
//    for(ServiceFunctionChain& sfc: simtest.allSFC){
//        simtest.convert_SeqSFC_to_FullParallel(sfc);
//        simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
//        simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
//        sfc.partialParallelChains = &sfc.subsetPartParSFC;
//        sfc.showSubsetPartialSFC();
//    }if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";



//    return 0;

    int id; cout<<"\n"<<"1. For Individual Testing. \n2. Automate for directories."; cin>>id;
    if(id == 1 and (init() == false)){
            return 1;
    }else if(id == 2 ){
            AutomateForEachDirectory();
    }

    cout<<"\n\nsuccess. press any key to exit";
    cout << "\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    return 0;
}


