/*!
 * @author Vijay Purohit
  Created by vijay on 22-02-2023.
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
using type_delay = double;
using type_wgt = double;
int debug = 1;
/**************** Custom Header Files Declaration ****************/
#include "Variables.h"
#include "PhysicalGraph.h" ///< graph structure for physical network, physical node and edge and node capacity
#include "VirtualNetworkFunctions.h" ///< structure for VNFnode, VNFs,
#include "ServiceFunctionChain.h" /// structure for SFC
#include "FileFunctions.h" ///< File reading writing functions
#include "SimulationClass.h" ///< algorithms implemented
#include "DelayCalculationFunctions.h" ///< Various Delay Calculation functions
#include "Algorithms.h" ///< algorithms implemented


//
//template<typename type_wgt, typename type_res>
//void Graph_NewYork_Test(){//NewYork_Test
//
//    string testName = "NewYork";
//    string testDirName = testName + "/"; ///< directory name with slash at the end.
//
//    PhysicalGraph *PhysicalNetwork; ///< Graph Object contains network related functions
//    readNetwork(testDirName, &PhysicalNetwork);
//
//    VirtualMachines *VirtualNetwork; ///< Virtual Machines Object contains Virtual Machines related functions
//    readVirtualMachines(testDirName, &VirtualNetwork);
//
//    VirtualNetworkFunctions *VNFNetwork; ///< VNF Object contains VNF (Virtual Network Function) related code.
//    readVirtualNetworkFunctions(testDirName, &VNFNetwork);
//
//    vector<ServiceFunctionChain *> SFCs; ///< SFCs object contains code related to Service function chains
//    vector<ServiceFunctionChain *> sortedSFCs; ///< SFCs sorted for deployement in order of their priority of traffic rate
//    readGenericServiceFunctionsChains(testDirName, SFCs, sortedSFCs, VNFNetwork);
//
//    { ///< assignment NEWYORK
///// Determination of how many instances of particular VNFs are required.
//        vector<pair<unsigned int, unsigned int>> VNF_TO_InstancesCnt = {
//                {1, 4},{2, 3},{3, 5},{4, 4},{5, 6},{6, 5},{7, 4},{8, 3},{9, 4}
//        };
//        assign_VNF_2_InstancesCnt(VNF_TO_InstancesCnt, VNFNetwork);
//
///// mapping of VNF id, instance id, to VM id
//        vector<vector<unsigned int>> VNF_TO_VM = {
//                {1, 1, 5}, {1, 2, 4}, {1, 3, 12}, {1, 4, 11}, {2, 1, 13}, {2, 2, 9}, {2, 3, 17}, {3, 1, 2}, {3, 2, 15}, {3, 3, 9}, {3, 4, 10}, {3, 5, 11}, {4, 1, 5}, {4, 2, 7}, {4, 3, 8}, {4, 4, 1}, {5, 1, 5}, {5, 2, 6}, {5, 3, 7}, {5, 4, 10}, {5, 5, 14}, {5, 6, 16}, {6, 1, 6}, {6, 2, 15}, {6, 3, 19}, {6, 4, 4}, {6, 5, 18}, {7, 1, 13}, {7, 2, 6}, {7, 3, 18}, {7, 4, 16}, {8, 1, 1}, {8, 2, 19}, {8, 3, 12}, {9, 1, 3}, {9, 2, 4}, {9, 3, 17}, {9, 4, 14}
//        };
//        assign_VNF_2_VM(VNF_TO_VM, VNFNetwork, VirtualNetwork);
//
//        vector<pair<int, int>> VM_TO_PN = {
//                {1, 1},{2, 7},{3, 15},{4, 1},{5, 7},{6, 15},{7, 2},{8, 3},{9, 4},{10, 5},{11, 6},{12, 8},{13, 9},{14, 10},{15, 11},{16, 12},{17, 13},{18, 14},{19, 16}
//        };
//        assign_VM_2_PN<type_wgt, type_res>(VM_TO_PN, VirtualNetwork, PhysicalNetwork);
//    } ///< assignmen
//
//    vector<unsigned int> parallelism_threshold = {40,45,50,55,60,65,70,75};
//    vector<type_delay> pktfactor = { 0.001, 0.005, 0.01, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.1};
//
//    for(const auto& thres: parallelism_threshold){//threshold
//        int observation = 1;
//        parallelPairs.clear();
//        findRandomParallelPairs(testDirName, observation, thres, VNFNetwork);
//        for (const auto &sfc: SFCs) {
//            sfc->vnfBlocksPar.clear();
//            sfc->allPartParSFC.clear();
//            convert_SeqSFC2ParVNFBlocks(sfc, VNFNetwork, (sfc->index == SFCs.size() - 1));
//            assign_Clusters2ParVNFs(sfc);
////                sfc->showParallelSFC(sfc->vnfBlocksPar);
//        }
//
//        observation++;
//        for(const auto& pktfac: pktfactor){ //pktfactor
//            VNFNetwork->seq_utilization.clear();
//            VNFNetwork->fullpar_utilization.clear();
//            VNFNetwork->ppar_utilization.clear();
//            read_write_time_per_bit = pktfac;
//
//            cout<<"\n\n["<<observation<<" | thre:"<<thres<<" | pktfac:"<<pktfac<<"]\n";
//
//            auto ft_start2 = std::chrono::steady_clock::now();
//            try {
//                algo_PartialChains_InstanceMapping<type_wgt, type_res>(sortedSFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork);
//            }
//            catch (std::exception const &e) {
//                std::cerr << "\ncaught: " << e.what();
//                VNFNetwork->showVNFs_Utilization(VNFNetwork->seq_utilization, 0);
//                VNFNetwork->showVNFs_Utilization(VNFNetwork->fullpar_utilization, 1);
//                VNFNetwork->showVNFs_Utilization(VNFNetwork->ppar_utilization, 2);
//            }
//            if (debug)
//                cout << "\nTime:" << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start2).count()<< "ms)";
//
////    showAlgoResults(SFCs, res_layerg);
////    showAlgoResults(SFCs, res_partial);
//
//            ofstream fout;
//            string filepathExt = output_directory+testDirName+"Sol_T_"+to_string(thres)+".csv";;///< path to .gv file without extention
//            fout.open(filepathExt.c_str(), ios::app);
//            if (!fout) {
//                string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
//                fout.clear();
//                throw runtime_error(errorMsg+ __FUNCTION__);
//            }
//            /// observation idx | sfc detial | seq res | full par res | partial par res | seq duration | full par duration | partial par duration | sfc
//            fout << "\n"<<thres<<" | "<< pktfac << ",";
//            for(const ServiceFunctionChain* const sfc: SFCs) {
//                fout<<sfc->index<<","<<sfc->trafficArrivalRate<<","<<sfc->numVNF<<","<<sfc->allPartParSFC.size()<<",";
//                const auto& obj= res_partial.solobj[sfc->index];
//                fout<<obj.seq_delay<<","<<obj.fullpar_delay<<","<<obj.ppar_delay<<",";
//                fout<<obj.seq_duration<<","<<obj.fullpar_duration<<","<<obj.ppar_duration<<", seq";
//
//                if(obj.seq_pid==noResDueToNoStg)fout<<"No Result Obtained for Sequential/Parallel due to no instance combination in one of the stage.";
//                else{
//                    if(obj.seq_pid == noResSeq)fout<<" No Result Obtained for Sequential.)";
//                    else{
//                        for(const auto &blk: sfc->allPartParSFC[obj.seq_pid]) { fout<<" [";  for(const auto& fnid: blk){ fout<<"f"<<fnid<<char(96+obj.seq_fninstmap.at(fnid))<<" "; }   fout<<"]";
//                        }
//                    }
//                    fout<<", fullpar";
//                    if(obj.fullpar_pid == noResSeq)fout<<" No Result Obtained for Full Par.)";
//                    else{
//                        for(const auto &blk: sfc->allPartParSFC[obj.fullpar_pid]) { fout<<" [";  for(const auto& fnid: blk){ fout<<"f"<<fnid<<char(96+obj.fullpar_fninstmap.at(fnid))<<" "; }   fout<<"]";
//                        }
//                    }
//                    fout<<", partpar";
//                    if(obj.ppar_pid == noResPar)fout << " No Result Obtained for Parallel.)";
//                    else{
//                        for(const auto &blk: sfc->allPartParSFC[obj.ppar_pid]) { fout << " [";  for(const auto& fnid: blk){ fout << "f" << fnid << char(96 + obj.ppar_fninstmap.at(fnid)) << " "; }   fout << "]";
//                        }
//                    }
//                }
//                fout<<"\n"<<",";
//            }
//            fout<<"seqdur:"<<res_partial.seq_duration<<" | fullpardur:"<< res_partial.fullpar_duration <<" | ppardur:"<<res_partial.ppar_duration<<"," ;
//            fout.close();
//
//        }//pktfactor
//    }//threshold
//
//    cout << "\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Test ] ~~~~~~~~~~~~~~~~~~~~~~~~";
//
//    for (int ni = 0; ni < SFCs.size(); ni++) {
//        delete SFCs[ni];
//        if (debug) {
//            if (ni == 0) cout << "\n[ SFCs Destructor Completed for sfc[" << ni << "] ";
//            else if (ni == SFCs.size() - 1) cout << "sfc[" << ni << "] ]";
//            else cout << "sfc[" << ni << "] ";
//        }
//    }
//    delete VNFNetwork;
//    delete VirtualNetwork;
//    delete PhysicalNetwork;
//
//}//NewYork_Test

void Sample0Test(){
    Simulations test("test","sample0");
    test.readDataFromFileInit("network.txt","VNFs.txt","SFCs.txt", {50,0});
//    test.DeploymentVNF_ScoreMethod(0.5,0.75,1,1);
    test.PN_2_VNF={
            {1, {{1,3}, {5,1} }},
            {2, {{2,2} }},
            {3, {}},
            {4, {}},
            {5, {{1,4}, {5,2} }},
            {6, {{1,2}, {2,3} }},
            {7, {}},
            {8, {{4,1}, {3,2} }},
            {9, {{6,2}, {7,1} }},
            {10, {{6,1}, {7,2} }},
            {11, {}},
            {12, {{9,1}, {3,1} }},
            {13, {}},
            {14, {{2,4}, {8,1} }},
            {15, {{1,1}, {4,2} }},
            {16, {{2,1}, {9,2} }}
    };
    test.finalInstancesCount={{1,4},  {2,4},  {3,2},  {4,2},  {5,2},  {6,2},  {7,2},  {8,1},  {9,2},  {10,0},};

    test.I_VNFINST_2_PN={
            {1, {{1,15}, {2,6}, {3,1}, {4,5} }},
            {2, {{1,16}, {2,2}, {3,6}, {4,14} }},
            {3, {{1,12}, {2,8} }},
            {4, {{1,8}, {2,15} }},
            {5, {{1,1}, {2,5} }},
            {6, {{1,10}, {2,9} }},
            {7, {{1,9}, {2,10} }},
            {8, {{1,14} }},
            {9, {{1,12}, {2,16} }},
            {10, {}}
    };
//    test.I_VNFINST_2_PN={
//            {1, {{1,15}, {2,6}, {3,1}, {4,5} }},
//            {2, {{1,16}, {2,2}, {3,6}, {4,14} }},
//            {3, {{1,12}, {2,12} }},
//            {4, {{1,8}, {2,15} }},
//            {5, {{1,1}, {2,5} }},
//            {6, {{1,10}, {2,9} }},
//            {7, {{1,9}, {2,10} }},
//            {8, {{1,12} }},
//            {9, {{1,12}, {2,16} }},
//            {10, {}}
//    };

    auto ft_start1 = std::chrono::steady_clock::now();
    try{  algo_LayerGraph_InstanceMapping(test);  }
    catch (std::exception const &e) {
        std::cerr << "\ncaught: " << e.what();
        test.showVNFsUtilization(1, test.TestsResult[name_layerg]);
        test.showVNFsUtilization(2, test.TestsResult[name_layerg]);
        test.showVNFsUtilization(3, test.TestsResult[name_layerg]);
    }
    if(debug) cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start1).count()<<"ms)";
//    test.showSimulationTestResults(test.TestsResult[name_layerg]);

    auto ft_start2 = std::chrono::steady_clock::now();
    try {  algo_PartialChains_InstanceMapping(test); }
    catch (std::exception const &e) {
        std::cerr << "\ncaught: " << e.what();
        test.showVNFsUtilization(1, test.TestsResult[name_partial]);
        test.showVNFsUtilization(2, test.TestsResult[name_partial]);
        test.showVNFsUtilization(3, test.TestsResult[name_partial]);
    }
    if (debug) cout << "\nTime:" << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start2).count()<< "ms)";
//    test.showSimulationTestResults(test.TestsResult[name_partial]);

    ofstream fout;
    string filepathExt = output_directory+test.fullDirName+"S_"+test.sim_name+"_SFCResult.csv";///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    /// observation idx | sfc detial | seq res | full par res | partial par res | seq duration | full par duration | partial par duration | sfc
    fout<<"observation, sfc idx, sfc traffic rate, sfc numVNFs, sfc cntPartialChains, ";
    fout<<"method, Seq Delay, FullPar Delay, PartPar Delay, Seq Duration, FullPar Duration, PartPar Duration, , seq sfc, fullpar sfc, partpar sfc, ";
    fout<<"method, Seq Delay, FullPar Delay, PartPar Delay, Seq Duration, FullPar Duration, PartPar Duration, , seq sfc, fullpar sfc, partpar sfc, ";
    fout << "\n";
    fout<< 1 << ",";
    for(const ServiceFunctionChain& sfc: test.allSFC) {
        fout<<sfc.index<<", "<<sfc.trafficArrivalRate<<", "<<sfc.numVNF<<", "<<sfc.allPartParSFC.size()<<", ";

        for(const auto& obj: test.TestsResult){
            const auto& sfcres= obj.second.sfcsol.at(sfc.index);
            fout<<obj.first<<",";
            fout<<sfcres.seq_delay<<","<<sfcres.fullpar_delay<<","<<sfcres.ppar_delay<<",";
            fout<<sfcres.seq_duration<<","<<sfcres.fullpar_duration<<","<<sfcres.ppar_duration<<", ,seq";
            if(sfcres.seq_pid==noResDueToNoStg)fout<<"No Result Obtained for Sequential/Parallel due to no instance combination in one of the stage.";
            else{
                if(sfcres.seq_pid == noResSeq)fout<<" No Result Obtained for Sequential";
                else{
                    for(const auto &blk: sfc.allPartParSFC[sfcres.seq_pid]) { fout<<" [";  for(const auto& fnid: blk){ fout<<"f"<<fnid<<char(96+sfcres.seq_fninstmap.at(fnid))<<" "; }   fout<<"]";
                    }
                }
                fout<<", fullpar";
                if(sfcres.fullpar_pid == noResSeq)fout<<" No Result Obtained for Full Par";
                else{
                    for(const auto &blk: sfc.allPartParSFC[sfcres.fullpar_pid]) { fout<<" [";  for(const auto& fnid: blk){ fout<<"f"<<fnid<<char(96+sfcres.fullpar_fninstmap.at(fnid))<<" "; }   fout<<"]";
                    }
                }
                fout<<", partpar";
                if(sfcres.ppar_pid == noResPar)fout << " No Result Obtained for Parallel";
                else{
                    for(const auto &blk: sfc.allPartParSFC[sfcres.ppar_pid]) { fout << " [";  for(const auto& fnid: blk){ fout << "f" << fnid << char(96 + sfcres.ppar_fninstmap.at(fnid)) << " "; }   fout << "]";
                    }
                }
            }
            fout<<",";
        }
        fout<<"\n"<<",";
    }
    for(const auto& obj: test.TestsResult){
        fout<<obj.first<<",";
        fout<<""<<obj.second.seq_duration<<" , "<< obj.second.fullpar_duration <<" , "<<obj.second.ppar_duration<<",," ;
    }
    fout.close();
}

// TODO: Finding number of VNF instances, and then VNF Deployement
// TODO: parallelism -> exploring more partial sfc possibilities.
// TODO: existing heuritics compare with bruteforce
// TODO: within same server, packet duplication required or not.
int main()
{
//    clusterSizeEnumeration(int(clusterSz.size()),maxSFClen);
//    all_nCk(int(nCk.size()), maxSFClen);
    Sample0Test();
//    Graph_NewYork_Test<type_wgt_l, type_res_l>();


    cout << "\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    return 0;
}



//vector<vector<unsigned int>> pp1;
//unordered_map<unsigned int, unsigned int> pp1map;
//    pp1 = {{2},{1,7},{6}, {5}};
//    pp1map = {
//            {2,3},
//            {1,2},
//            {7,2},
//            {6,1},
//            {5,2}
//    };
//    cout<<"\n1 ::"<<calcD_ParallelSFC(test.allSFC[5], pp1,pp1map,test.TestsResult[name_partial].fullpar_utilization, test,true);

//    partialChains_Sequential_Deployment(test, 0, test.allSFC[0].allPartParSFC.front());
//    partialChains_Sequential_Deployment(test, 1, test.allSFC[1].allPartParSFC.front());
//    partialChains_Sequential_Deployment(test, 2, test.allSFC[2].allPartParSFC.front());

//    partialChains_FullParallel_Deployment(test, 0, test.allSFC[0].allPartParSFC.back(), true);
//    partialChains_FullParallel_Deployment(test, 1, test.allSFC[1].allPartParSFC.back(), true);
//    partialChains_FullParallel_Deployment(test, 2, test.allSFC[2].allPartParSFC.back(), true);

//    partialChains_PartParallel_Deployment(test, 0);
//    partialChains_PartParallel_Deployment(test, 1,true);
//    partialChains_PartParallel_Deployment(test, 2, true);
