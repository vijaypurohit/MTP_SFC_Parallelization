/*!
 * @author Vijay Purohit
 * Created by vijay on 22-02-2023.
 * For MTP Thesis
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
#include "SimulationClass.h" ///< algorithms implemented
#include "DelayCalculationFunctions.h" ///< Various Delay Calculation functions
#include "Algorithms.h" ///< algorithms implemented

bool SimulationTest_VaryingNetworkVNF_DeploymentParameter(){//SimulationTest_VaryingNetworkVNF_DeploymentParameter
    Simulations test("VaryingDeploymentParameter_dsc_rate","sample0");
    test.readDataFromFileInit("network.txt","newVNF.txt","newSFC.txt", {60,0}, "dsc_rate");
    int obs=1;
    for( const float& scale: {0.25f,0.5f,0.75f,1.0f} ){
        for(const float& alpha: {0.25f,0.5f,0.75f}){
            for(const int& pxs: {1,2}){
                for(const int& dist: {1,2}){
                    system("cls");
                    test.DeploymentVNF_ScoreMethod(scale,alpha,pxs,dist);
                    test.TestsResult.clear();

                    auto ft_start1 = std::chrono::steady_clock::now();
                    try{  Heuristic_kShortestPath_InstanceMapping(test);  }
                    catch (std::exception const &e) {
                        std::cerr << "\ncaught: " << e.what();
                        test.showSimulationTestResults(test.TestsResult[name_kshortestpath]);
                    }
                    if(debug) cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start1).count()<<"ms)";
//    test.showSimulationTestResults(test.TestsResult[name_kshortestpath]);

                    auto ft_start2 = std::chrono::steady_clock::now();
                    try {  bruteForce_InstanceMapping(test); }
                    catch (std::exception const &e) {
                        std::cerr << "\ncaught: " << e.what();
                        test.showSimulationTestResults(test.TestsResult[name_bruteForce]);
                    }
                    if (debug) cout << "\nTime:" << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start2).count()<< "ms)";

//    test.showSimulationTestResults(test.TestsResult[name_bruteForce]);
                    string other = to_string(scale)+"_"+to_string(alpha)+"_"+to_string(pxs)+"_"+to_string(dist);
                    if(obs == 1){
                        test.writeInFileSimulationTestResults(obs, ios::out, other);
                    }else{
                        test.writeInFileSimulationTestResults(obs, ios::app, other);
                    }
                    obs++;
                }
            }
        }
    }
    return true;
}//SimulationTest_VaryingNetworkVNF_DeploymentParameter

int main()
{
    /// Run this if SFC Len is going to be more than 10
    if(maxSFCLength > 10){
        integerCompositionsEnumeration(int(integerCompositions.size()),maxSFCLength);
        find_all_nCk(int(nCk.size()), maxSFCLength);
    }
    /// max number of VNFs are 12

    ifRejectSFC = false;
    packetBodySize = 1000; packetHeaderSize = 24; factor_packet = 8;
    bandwidthNW = 1; factor_bandwidth = 1000;

    speedOfLight = 300000; velocityFactor = 1.0;
    read_write_time_per_bit = 0.077e-2;

    return SimulationTest_VaryingNetworkVNF_DeploymentParameter();
 
 
//    unsigned int numPN, numVNF, numOfSFCsNeeded;
//    pair<bool, unsigned int>  lenOfEachSFC;
//    pair<bool, string> WriteToFile;
//    cout<<"\n Number of Physical Nodes: ";
//    cin>>numPN;
//    cout<<"\n Number of VNFs: ";
//    cin>>numVNF;
//    cout<<"\n Number of SFCs to generate: ";
//    cin>>numOfSFCsNeeded;
//    cout<<"\n 0. Random Length of each SFC. ";
//    cout<<"\n 1. Fixed Length of each SFC. ";
//    cin>>lenOfEachSFC.first;
//    if(lenOfEachSFC.first == 1){
//        cout<<"\n\t Length of the Each SFC: ";
//        cin>>lenOfEachSFC.second;
//    }
//
    Simulations test("tt","sample0");

//    GenerateRandomVNFs(numVNF, {1.5,3}, {2,4}, test.fullDirName , "newVNF.txt");
//    GenerateRandomSFCs(numPN, numVNF, numOfSFCsNeeded, {0.3, 0.5}, lenOfEachSFC, test.fullDirName , "newSFC.txt");

    test.readDataFromFileInit("network.txt","newVNF.txt","newSFC.txt", {60,0}, "dsc_rate");
    for(const auto* sfc: test.sortedSFCs){
        (sfc)->showSequentialSFC();
    }
    return 0;
//    test.DeploymentVNF_ScoreMethod(0.25,0.25,1,2);
//    cout<<"\n\n";
//    test.sfccompleted=0;
//    for(const ServiceFunctionChain*const& sfcpointer:test.sortedSFCs) {
//        const ServiceFunctionChain& cSFC = *sfcpointer;
//
//        test.sfccompleted++;
////        cout<<"\r  H-sequential["<<test.sfccompleted<<"/"<<test.sortedSFCs.size()<<"]";
//
//        test.TestsResult[name_kshortestpath].sfcsol[cSFC.index] = SFC_RESULT(); ///< solution for the layer graph algorithm
//        auto sfc_st = std::chrono::steady_clock::now();
//        for(const unsigned int& fn: cSFC.vnfSeq){ /// finding some vnf delays
//            const VNFNode& dstVNFNode = test.VNFNetwork.VNFNodes.at(fn );
//            test.vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF(dstVNFNode);
//            test.vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay(dstVNFNode);
//            for(int fnInst=1; fnInst<=test.finalInstancesCount.at(fn); fnInst++) { ///sequential Queuing Delay
//                test.vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay(cSFC.trafficArrivalRate, dstVNFNode, fnInst, test.TestsResult[name_kshortestpath].seq_utilization);
//            }
//        }
//        //        auto sfc_st = std::chrono::steady_clock::now();
////        if(test.sfccompleted == 1) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = {  {9,1},{8,2},{4,2},{1,2} };
////        }
////        if(test.sfccompleted == 2) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {7,1},{8,1},{6,2},{4,1} };
////        }
////        if(test.sfccompleted == 3) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {6,1},{8,2},{12,1},{1,1},{11,2} };
////        }
////        if(test.sfccompleted == 4) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {2,2},{10,2},{3,1},{12,2},{9,1} };
////        }
////        if(test.sfccompleted == 5) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {3,1},{8,1},{9,2},{4,2},{7,2},{6,1},{1,2},{10,1} };
////        }
////        if(test.sfccompleted == 6) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {12,1},{9,2},{4,1},{1,1},{7,1},{5,1},{10,2},{6,2} };
////        }
////        if(test.sfccompleted == 7) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {7,2},{12,1},{9,2},{5,1},{8,2},{2,1},{11,2},{10,1},{6,1} };
////        }
////        if(test.sfccompleted == 8) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {2,1},{6,3},{1,2},{7,2},{8,1},{11,1},{12,2},{10,1},{4,2} };
////        }
////        if(test.sfccompleted == 9) {
////            test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap = { {8,2},{4,2},{2,1},{1,1},{11,2},{5,1},{7,1},{6,2},{10,2},{9,1} };
////        }
//        if(test.sfccompleted < 10){
////            for(const auto& [fn, fninst]: test.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_fninstmap){
////                test.TestsResult[name_kshortestpath].seq_utilization[fn][fninst] += cSFC.trafficArrivalRate;
////            }
//            kShortestPath_Sequential_Deployement(test, cSFC, showFinal);
////            kShortestPath_FullParallel_Deployment(test, cSFC, showFinal);
//        }else{
////            unordered_map<unsigned int, unsigned int> ppmap ={
////                    {7,2},{4,1},{11,2},{6,3},{1,2},{10,1},{5,1},{9,2},{12,1},{2,2}
////            };
////            cout<<"\nH-d1: "<<calcD_SequentialSFC(cSFC,ppmap , test.TestsResult[name_kshortestpath].seq_utilization, test, showFinal);
////            cout<<"\nH-d2: "<<calcD_ParallelSFC(cSFC, cSFC.subsetPartParSFC.front() ,ppmap , test.TestsResult[name_kshortestpath].seq_utilization, test, showFinal);
//            kShortestPath_Sequential_Deployement(test, cSFC, showFinal);
////            kShortestPath_FullParallel_Deployment(test, cSFC, showFinal);
//        }
//    }

//    test.sfccompleted=0; cout<<"\n\n";
//    for(const ServiceFunctionChain*const& sfcpointer:test.sortedSFCs) {
//        const ServiceFunctionChain& cSFC = *sfcpointer;
//        test.sfccompleted++;
////        cout<<"\r  BF-sequential["<<test.sfccompleted<<"/"<<test.sortedSFCs.size()<<"]";
//
//        test.TestsResult[name_bruteForce].sfcsol[cSFC.index] = SFC_RESULT(); /// using this algo what is optimal result for a single sfc.
//
////        auto sfc_st = std::chrono::steady_clock::now();
////        if(test.sfccompleted == 1) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = {  {9,2},{8,2},{4,2},{1,2} };
////        }
////        if(test.sfccompleted == 2) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {7,1},{8,1},{6,2},{4,1} };
////        }
////        if(test.sfccompleted == 3) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {6,1},{8,2},{12,1},{1,1},{11,2} };
////        }
////        if(test.sfccompleted == 4) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {2,2},{10,2},{3,1},{12,2},{9,1} };
////        }
////        if(test.sfccompleted == 5) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {3,1},{8,1},{9,1},{4,2},{7,2},{6,1},{1,2},{10,1} };
////        }
////        if(test.sfccompleted == 6) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {12,1},{9,2},{4,1},{1,1},{7,2},{5 ,1},{10,1},{6,2} };
////        }
////        if(test.sfccompleted == 7) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {7,1},{12,1},{9,2},{5,1},{8,2},{ 2,1},{11,2},{10,2},{6,1} };
////        }
////        if(test.sfccompleted == 8) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {2,1},{6,3},{1,2},{7,1},{8,1},{11,1},{12,2},{10,2},{4,2} };
////        }
////        if(test.sfccompleted == 9) {
////            test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap = { {8,2},{4,2},{2,1},{1,1},{11,2},{5,1},{7,1},{6,2},{10,1},{9,1} };
////        }
//
//
//        if(test.sfccompleted < 10){
////            for(const auto& [fn, fninst]: test.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_fninstmap){
////                test.TestsResult[name_bruteForce].seq_utilization[fn][fninst] += cSFC.trafficArrivalRate;
////            }
//            bruteForce_Sequential_Deployment(test, cSFC, showFinal);
////            bruteForce_FullParallel_Deployment(test, cSFC, showFinal);
//        }else{
////            unordered_map<unsigned int, unsigned int> ppmap ={
////                    {7,2},{4,2},{11,1},{6,3},{1,2},{10,2},{5,1},{9,1},{12,1},{2,2}
////            };
////            cout<<"\nBF-d1: "<<calcD_SequentialSFC(cSFC,ppmap , test.TestsResult[name_bruteForce].seq_utilization, test, showFinal);
////            cout<<"\nBF-d2: "<<calcD_ParallelSFC(cSFC, cSFC.subsetPartParSFC.front() ,ppmap , test.TestsResult[name_bruteForce].seq_utilization, test, showFinal);
//            bruteForce_Sequential_Deployment(test, cSFC, showFinal);
////            bruteForce_FullParallel_Deployment(test, cSFC, showFinal);
//        }
//    }
///

    auto ft_start1 = std::chrono::steady_clock::now();
    try{  Heuristic_kShortestPath_InstanceMapping(test);  }
    catch (std::exception const &e) {
        std::cerr << "\ncaught: " << e.what();
        test.showSimulationTestResults(test.TestsResult[name_kshortestpath]);
    }
    if(debug) cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start1).count()<<"ms)";

    auto ft_start2 = std::chrono::steady_clock::now();
    try {  bruteForce_InstanceMapping(test); }
    catch (std::exception const &e) {
        std::cerr << "\ncaught: " << e.what();
        test.showSimulationTestResults(test.TestsResult[name_bruteForce]);
    }
    if (debug) cout << "\nTime:" << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start2).count()<< "ms)";

    test.showSimulationTestResults(test.TestsResult[name_kshortestpath]);
    test.showSimulationTestResults(test.TestsResult[name_bruteForce]);
    test.writeInFileSimulationTestResults(1, ios::out);


    cout << "\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    return 0;
}


//    test.PN_2_VNF={
//            {1, {{2,2}, {11,1}, }},
//            {2, {{8,2}, {9,1}, }},
//            {3, {{1,2}, {5,1}, }},
//            {4, {{9,2}, {7,1}, }},
//            {5, {}},
//            {6, {{6,4}, {8,3}, }},
//            {7, {{2,1}, {11,2}, }},
//            {8, {{10,1}, {6,3}, }},
//            {9, {{7,2}, {12,1}, }},
//            {10, {{3,1}, }},
//            {11, {{2,3}, {11,3}, }},
//            {12, {{10,4}, {8,1}, }},
//            {13, {{12,2}, {4,1}, }},
//            {14, {{10,3}, {6,1}, }},
//            {15, {{4,2}, {1,1}, }},
//            {16, {{10,2}, {6,2}, }},
//    };
//    test.finalInstancesCount={{1,2},  {2,3},  {3,1},  {4,2},  {5,1},  {6,4},  {7,2},  {8,3},  {9,2},  {10,4},  {11,3},  {12,2},  };
//
//    test.I_VNFINST_2_PN={
//            {1, {{1,15}, {2,3}, }},
//            {2, {{1,7}, {2,1}, {3,11}, }},
//            {3, {{1,10}, }},
//            {4, {{1,13}, {2,15}, }},
//            {5, {{1,3}, }},
//            {6, {{1,14}, {2,16}, {3,8}, {4,6}, }},
//            {7, {{1,4}, {2,9}, }},
//            {8, {{1,12}, {2,2}, {3,6}, }},
//            {9, {{1,2}, {2,4}, }},
//            {10, {{1,8}, {2,16}, {3,14}, {4,12}, }},
//            {11, {{1,1}, {2,7}, {3,11}, }},
//            {12, {{1,9}, {2,13}, }},
//    };