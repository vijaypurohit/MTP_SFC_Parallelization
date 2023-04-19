/*!
 * @author Vijay Purohit
  Created by vijay on 22-02-2023.
 */
#include <iostream>
#include <utility>
#include <vector>
#include <queue>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
//#include <execution>
//#include <limits> // numeric_limits
//#include <chrono>       //to seed random number
#include <random> //to generate random number
#include <functional>  // lambda function


using namespace std;
using type_delay = double;
int debug = 1;
/**************** Custom Header Files Declaration ****************/
#include "Variables.h"
#include "PhysicalGraph.h" ///< graph structure for physical network, physical node and edge and node capacity
#include "VirtualMachines.h" ///<  structure for virtual machines, virtual node
#include "VirtualNetworkFunctions.h" ///< structure for VNFnode, VNFs,
#include "ServiceFunctionChain.h" /// structure for SFC
#include "FileFunctions.h" ///< File reading writing functions
#include "DelayCalculationFunctions.h" ///< Various Delay Calculation functions
#include "AssignmentFunctions.h"
#include "Algorithms.h" ///< algorithms implemented

template<typename type_wgt, typename type_res>
vector<unsigned int> constructPath(unsigned int src, unsigned int dst, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork)
{
    /// If there's no path between node src and dst, simply return an empty array
    if (PhysicalNetwork->nextHop[src][dst] == 0)
        return {};
    // Storing the path in a vector
    vector<unsigned int> path = { src };
    while (src != dst) {
        src = PhysicalNetwork->nextHop[src][dst];
        path.push_back(src);
    }
    return path;
}

template<typename type_wgt, typename type_res>
void solve(const vector<ServiceFunctionChain*>& SFCs, VirtualNetworkFunctions<type_res> *const& VNFNetwork, const VirtualMachines<type_res> *const& VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const& PhysicalNetwork){//solve

//    for(int fn=1; fn<=VNFNetwork->numVNF; fn++){
//        cout<<"\nf: "<<fn<<" = "<<VNFNetwork->cntVNF[fn].size();
//    }cout<<endl;

    const unsigned int numSFCs = SFCs.size();

    pair<unsigned int, unsigned int> access_nodes[numSFCs];
    vector<vector<unsigned int>> bstCandPaths(numSFCs);

    vector<unordered_map<type_delay , unordered_map<type_delay, type_delay>>> deltac(numSFCs);
    vector<unordered_map<unsigned int, unsigned int>> bst_candidate_nodes(numSFCs);

    access_nodes[0] = {2, 6};//shortest path 2,1,5,6
    access_nodes[1] = {1, 8}; //1,2,3,4,8
    access_nodes[2] = {6, 7}; //6,5,3,7

    bstCandPaths[0] = constructPath(access_nodes[0].first, access_nodes[0].second, PhysicalNetwork);
    bstCandPaths[1] = constructPath(access_nodes[1].first, access_nodes[1].second, PhysicalNetwork);
    bstCandPaths[2] = constructPath(access_nodes[2].first, access_nodes[2].second, PhysicalNetwork);

    for(const auto &path: bstCandPaths) {
        for(const auto &bc: path) cout<<bc<<" -> ";cout<<endl;
    }


    for(int sidx=0; sidx<numSFCs; sidx++){/*! for each sfc */
        unordered_map<unsigned int, unsigned int>& fi_bstloc = bst_candidate_nodes[sidx];
        unordered_map<type_delay, unordered_map<type_delay, type_delay>>& sfcdelta_fi_pn = deltac[sidx];

        for(int oi=0; oi<SFCs[sidx]->numVNF; oi++){
            const auto fi= SFCs[sidx]->vnfSeq[oi];
            unsigned int hi_star =   ceil ((oi*bstCandPaths[sidx].size()) / SFCs[sidx]->numVNF);
        cout<<"f"<<fi<<" id:"<<hi_star<<" pn:"<< bstCandPaths[sidx][hi_star]<<endl;
            fi_bstloc[fi] = bstCandPaths[sidx][hi_star];
        }

        for(int oi=0; oi<SFCs[sidx]->numVNF; oi++){
            const unsigned int fi= SFCs[sidx]->vnfSeq[oi];
            unsigned int vhi_star = fi_bstloc[fi];

            unsigned int vhi_1_star ;

            if(oi ==  SFCs[sidx]->numVNF - 1){
                vhi_1_star = access_nodes[sidx].second;
            }else{
                unsigned int fi_1 = SFCs[sidx]->vnfSeq[oi+1];
                vhi_1_star = fi_bstloc[fi_1];
            }

            cout<<"\nf"<<fi<<" vhi*:"<<vhi_star<<" vhi1*:"<< vhi_1_star<<" :: ";
            for(int vh=1; vh<=PhysicalNetwork->numV; vh++)
            {
                const auto dist_vhi_star_2_vh = PhysicalNetwork->dist[vhi_star][vh];
                const auto dist_vh_2_vhi_1_star = PhysicalNetwork->dist[vh][vhi_1_star];
                const auto hopc_vhi_star_2_vh = constructPath(vhi_star, vh, PhysicalNetwork).size()-1;
                const auto hopc_vh_2_vhi_1_star = constructPath(vh, vhi_1_star, PhysicalNetwork).size()-1;

                sfcdelta_fi_pn[fi][vh] = dist_vhi_star_2_vh + dist_vh_2_vhi_1_star;

                cout<<"\n\t";
                cout<<vhi_star<<"->"<<vh<<" = "<<hopc_vhi_star_2_vh<<" ("<<dist_vhi_star_2_vh<<"m)  |  ";
                cout<<vh<<"->"<<vhi_1_star<<" = "<<hopc_vh_2_vhi_1_star<<" ("<<dist_vh_2_vhi_1_star<<"m) ";
            }
        }
    }/*! for each sfc */





}//solve



// TODO: Finding number of VNF instances, and then VNF Deployement
// TODO: parallelism -> exploring more partial sfc possibilities.
// TODO: existing heuritics compare with bruteforce
// TODO: within same server, packet duplication required or not.
int main()
{
    using type_wgt_l = type_delay; ///< Determine Edge Weights Data TYPE (unsigned int or FLOAT or DOUBLE)
    using type_res_l = unsigned int; ///< Determine Resource Data TYPE (unsigned int or FLOAT)

int observation = 1; int totalobservation=50;
//while(observation < 200) {
    // Object Creations

//    clusterSizeEnumeration(int(clusterSz.size()),maxSFClen);
//    all_nCk(int(nCk.size()), maxSFClen);

//    unsigned int testIdx = 0; ///< Test Directory Initialisation.
//    string testName = "sample" + to_string(testIdx); ///< should be without space
    string testName = "sample0";
    string testDirName = testName + "/"; ///< directory name with slash at the end.
//    createDirectory(testDirName);

    PhysicalGraph<type_wgt_l, type_res_l> *PhysicalNetwork; ///< Graph Object contains network related functions
    readNetwork<type_wgt_l, type_res_l>(testDirName, &PhysicalNetwork);
//    PhysicalNetwork->showAdjMatrix();
//    PhysicalNetwork->showAllPairsShortestPath();
//    PhysicalNetwork->showAdjList();
    VirtualMachines<type_res_l> *VirtualNetwork; ///< Virtual Machines Object contains Virtual Machines related functions
    readVirtualMachines<type_res_l>(testDirName, &VirtualNetwork);


    VirtualNetworkFunctions<type_res_l> *VNFNetwork; ///< VNF Object contains VNF (Virtual Network Function) related code.
    readVirtualNetworkFunctions<type_res_l>(testDirName, &VNFNetwork);
    VNFNetwork->findRandomParallelPairs(testDirName, observation);

    vector<ServiceFunctionChain *> SFCs; ///< SFCs object contains code related to Service function chains
    vector<ServiceFunctionChain *> sortedSFCs; ///< SFCs sorted for deployement in order of their priority of traffic rate
    readGenericServiceFunctionsChains<type_res_l>(testDirName, SFCs, sortedSFCs, VNFNetwork);

    for (const auto &sfc: SFCs) {
        convert_SeqSFC2ParVNFBlocks<type_res_l>(sfc, VNFNetwork, (sfc->index == SFCs.size() - 1));
        assign_Clusters2ParVNFs(sfc);
//        sfc->showParallelSFC(sfc->vnfBlocksPar);
    }
//    showSFCsDescriptions(SFCs);

    {
        ///< assignment
        /// Determination of how many instances of particular VNFs are required.
        vector<pair<unsigned int, unsigned int>> VNF_TO_InstancesCnt = {
                {1,  3},{4,  3},{2,  2},{3,  2},{5,  2},{6,  2},{7,  2},{8,  2},{9,  2},{10, 2}
        };
        assign_VNF_2_InstancesCnt<type_res_l>(VNF_TO_InstancesCnt, VNFNetwork);

        /// mapping of VNF id, instance id, to VM id
        vector<vector<unsigned int>> VNF_TO_VM = {
                {1,  1, 3},{1,  2, 1},{1,  3, 12},
                {2,  1, 4},{2,  2, 2},
                {3,  1, 1},{3,  2, 7},
                {4,  1, 2},{4,  2, 13},{4,  3, 12},
                {5,  1, 9},{5,  2, 3},
                {6,  1, 7},{6,  2, 5},
                {7,  1, 8},{7,  2, 6},
                {8,  1, 10},{8,  2, 14},
                {9,  1, 5},{9,  2, 9},
                {10, 1, 11},{10, 2, 11}
        };
        assign_VNF_2_VM<type_res_l>(VNF_TO_VM, VNFNetwork, VirtualNetwork);

        vector<pair<int, int>> VM_TO_PN = {
                {1,  1},{2,  1},{3,  2},{4,  2},{5,  3},{6,  3},{14, 3},{7,  4},{8,  4},{11, 5},{12, 5},{13, 5},{9,  6},{10, 6}
        };
        assign_VM_2_PN<type_wgt_l, type_res_l>(VM_TO_PN, VirtualNetwork, PhysicalNetwork);
    //    vector<vector<pair<int,int>>> VNFType_TO_InstID(SFCs.size());
        //        VNFType_TO_InstID[0] = { {1, 1}, {6, 2}, {4, 1}, {5, 1} };
        //        VNFType_TO_InstID[1] = { {10, 2}, {4, 2}, {8, 1}, {3, 1}, {2, 1}, {1, 1}, {5, 2} };
        //        VNFType_TO_InstID[2] = { {4, 3}, {8, 2}, {7, 1}, {10, 1}, {9, 1} };
        //    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[1], SFCs[1], VNFNetwork);
        //    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[2], SFCs[2], VNFNetwork);
        //    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[3], SFCs[3], VNFNetwork);
    } ///< assignment

//    PhysicalNetwork->showPNs_Description();
//    VirtualNetwork->showVMs_Description();
//    VNFNetwork->showVNFs_Description();
//    for(int ni=0; ni<SFCs.size(); ni++) SFCs[ni]->showSequentialSFC();
//    for(int ni=0; ni<SFCs.size(); ni++) SFCs[ni]->showParallelSFC(sfc->vnfBlocksPar);
    solve(SFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork);
    return 0;

//    vector<vector<unsigned int>> pp1 = {{10},{4,8,3},{2},{1}, {5}};
//    unordered_map<unsigned int, unsigned int> pp1map = {
//            {10,1},
//            {4,3}, {8,1}, {3,1},
//            {2,2}, {1,3},
//            {5,1}
//    };
//    cout<<"\n1 ::"<<calcD_ParallelSFC<type_wgt_l, type_res_l>(pp1, pp1map,VNFNetwork->fullpar_utilization ,SFCs[1], VNFNetwork, VirtualNetwork, PhysicalNetwork,true);
//    pp1map = {
//            {10,1},
//            {4,1}, {8,2}, {3,1},
//            {2,2}, {1,1},
//            {5,2}
//    };
//    cout<<"\n1 ::"<<calcD_ParallelSFC<type_wgt_l, type_res_l>(pp1, pp1map,VNFNetwork->fullpar_utilization ,SFCs[1], VNFNetwork, VirtualNetwork, PhysicalNetwork,true);

//    unordered_map<unsigned int, vnfDelaysPreComputed> vnfDelays;///< pre-calculated VNF delays (processing, execution and queuing delay). Before iterating all the partial par sfc it is better to calculate it for each chain as they will be same for each chain.
//    for(const unsigned int& fn: SFCs[1]->vnfSeq){ /// finding some vnf delays
//        VNFNode<type_res_l> *const dstVNFNode = VNFNetwork->VNFNodes.at(fn );
//        vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF<type_res_l>(dstVNFNode);
//        vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay<type_res_l>(dstVNFNode);
//        for(int fnInst=1; fnInst<=dstVNFNode->numInstances; fnInst++) { ///sequential Queuing Delay
//            vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay<type_res_l>(dstVNFNode, fnInst, VNFNetwork->fullpar_utilization, SFCs[1]);
//        }
//    }
//    layerGraph_FullParallel_Deployment<type_wgt_l, type_res_l>(res_layerg.solobj[SFCs[1]->index], SFCs[1]->allPartParSFC[SFCs[1]->allPartParSFC.size()-1], SFCs[1], vnfDelays, VNFNetwork, VirtualNetwork, PhysicalNetwork, true, true);

//    cout<<"\npartial sfc2::"<<calcD_SequentialSFC<type_wgt_l, type_res_l>(SFCs[1],pp1map, VNFNetwork, VirtualNetwork, PhysicalNetwork,true, true);

    cout<<"\n[observation: "<<observation<<"/"<<totalobservation<<"]";
    auto ft_start1 = std::chrono::steady_clock::now();
    try{
        algo_LayerGraph_InstanceMapping<type_wgt_l, type_res_l>(sortedSFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork);
    }
    catch (std::exception const &e) {
        std::cerr << "\ncaught: " << e.what();
        VNFNetwork->showVNFs_Utilization(VNFNetwork->seq_utilization, 0);
        VNFNetwork->showVNFs_Utilization(VNFNetwork->fullpar_utilization, 1);
        VNFNetwork->showVNFs_Utilization(VNFNetwork->ppar_utilization, 2);
    }
    if(debug)
        cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start1).count()<<"ms)";



    VNFNetwork->seq_utilization.clear();
    VNFNetwork->fullpar_utilization.clear();
    VNFNetwork->ppar_utilization.clear();


    auto ft_start2 = std::chrono::steady_clock::now();
    try {
        algo_PartialChains_InstanceMapping<type_wgt_l, type_res_l>(sortedSFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork);
    }
    catch (std::exception const &e) {
        std::cerr << "\ncaught: " << e.what();
        VNFNetwork->showVNFs_Utilization(VNFNetwork->seq_utilization, 0);
        VNFNetwork->showVNFs_Utilization(VNFNetwork->fullpar_utilization, 1);
        VNFNetwork->showVNFs_Utilization(VNFNetwork->ppar_utilization, 2);
    }
    if (debug)
        cout << "\nTime:" << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start2).count()<< "ms)";

    showAlgoResults(SFCs, res_layerg);
    showAlgoResults(SFCs, res_partial);

//    ofstream fout;
//    string filepathExt = output_directory+testDirName+"all_sfc.csv";///< path to .gv file without extention
//    fout.open(filepathExt.c_str(), ios::app);
//    if (!fout) {
//        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
//        fout.clear();
//        throw runtime_error(errorMsg+ __FUNCTION__);
//    }
    /// observation idx | sfc detial | seq res | full par res | partial par res | seq duration | full par duration | partial par duration | sfc
//    fout << "\n";
//    fout<< observation << ",";
//    for(const ServiceFunctionChain* const sfc: SFCs) {
//        fout<<sfc->index<<"| "<<sfc->trafficArrivalRate<<" | VNFs:"<<sfc->numVNF<<" | partials: "<<sfc->allPartParSFC.size()<<", ,";
//        const auto& obj= res_partial.solobj[sfc->index];
//        fout<<obj.seq_delay<<","<<obj.fullpar_delay<<","<<obj.ppar_delay<<",";
//        fout<<obj.seq_duration<<","<<obj.fullpar_duration<<","<<obj.ppar_duration<<", seq";
//
//        if(obj.seq_pid==noResDueToNoStg)fout<<"No Result Obtained for Sequential/Parallel due to no instance combination in one of the stage.";
//        else{
//            if(obj.seq_pid == noResSeq)fout<<" No Result Obtained for Sequential.)";
//            else{
//                for(const auto &blk: sfc->allPartParSFC[obj.seq_pid]) { fout<<" [";  for(const auto& fnid: blk){ fout<<"f"<<fnid<<char(96+obj.seq_fninstmap.at(fnid))<<" "; }   fout<<"]";
//                }
//            }
//            fout<<", fullpar";
//            if(obj.fullpar_pid == noResSeq)fout<<" No Result Obtained for Full Par.)";
//            else{
//                for(const auto &blk: sfc->allPartParSFC[obj.fullpar_pid]) { fout<<" [";  for(const auto& fnid: blk){ fout<<"f"<<fnid<<char(96+obj.fullpar_fninstmap.at(fnid))<<" "; }   fout<<"]";
//                }
//            }
//            fout<<", partpar";
//            if(obj.ppar_pid == noResPar)fout << " No Result Obtained for Parallel.)";
//            else{
//                for(const auto &blk: sfc->allPartParSFC[obj.ppar_pid]) { fout << " [";  for(const auto& fnid: blk){ fout << "f" << fnid << char(96 + obj.ppar_fninstmap.at(fnid)) << " "; }   fout << "]";
//                }
//            }
//        }
//        fout<<"\n"<<",";
//    }
//    fout<<"seqdur:"<<res_partial.seq_duration<<" | fullpardur:"<< res_partial.fullpar_duration <<" | ppardur:"<<res_partial.ppar_duration<<"," ;
//
//    fout.close();
//    showAlgoResults(SFCs, res_partial);

    //    vector<vector<unsigned int>> pp1 = {{10},{4,8,3},{2,1}, {5}};
//    unordered_map<unsigned int, unsigned int> pp1map = {
//            {10,1},
//            {4,1}, {8,1}, {3,1},
//            {2,2}, {1,2},
//            {5,1}
//    };
//    cout<<"\npartial sfc2::"<<calcD_ParallelSFC<type_wgt_l, type_res_l>(pp1, pp1map,VNFNetwork->par_utilization ,SFCs[2], VNFNetwork, VirtualNetwork, PhysicalNetwork,true);
//      pp1 = {{10},{4,8,3},{2,1}, {5}};
//      pp1map = {
//            {10,1},
//            {4,2}, {8,1}, {3,1},
//            {2,2}, {1,2},
//            {5,2}
//    };
//    cout<<"\nlayergsfc2::"<<calcD_ParallelSFC<type_wgt_l, type_res_l>(pp1, pp1map,VNFNetwork->par_utilization ,SFCs[2], VNFNetwork, VirtualNetwork, PhysicalNetwork,true);


    cout << "\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    for (int ni = 0; ni < SFCs.size(); ni++) {
        delete SFCs[ni];
        if (debug) {
            if (ni == 0) cout << "\n[ SFCs Destructor Completed for sfc[" << ni << "] ";
            else if (ni == SFCs.size() - 1) cout << "sfc[" << ni << "] ]";
            else cout << "sfc[" << ni << "] ";
        }
    }
    delete VNFNetwork;
    delete VirtualNetwork;
    delete PhysicalNetwork;

//observation++;}
    return 0;
}

//#pragma omp parallel default(none)
//{ // Parallel Section
//// Lists number of threads and prints thread number currently executing print command
//printf("There are %d threads. Hello from thread %d\n", omp_get_num_threads(), omp_get_thread_num());
//}

//void debugPrint(const string& msg){
//    if(debug) cout<<msg;
//}

//template<class T>
//typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
//almost_equal(T x, T y, int ulp)
//{
//    // the machine epsilon has to be scaled to the magnitude of the values used and multiplied by the desired precision in ULPs (units in the last place)
//    return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp
//           // unless the result is subnormal
//           || std::fabs(x - y) < std::numeric_limits<T>::min();
////}
//{ ///< assignment
///// Determination of how many instances of particular VNFs are required.
//vector<pair<unsigned int, unsigned int>> VNF_TO_InstancesCnt = {
//        {1,  3},{4,  3},{2,  2},{3,  2},{5,  2},{6,  2},{7,  2},{8,  2},{9,  2},{10, 2}
//};
//assign_VNF_2_InstancesCnt<type_res_l>(VNF_TO_InstancesCnt, VNFNetwork);
//
///// mapping of VNF id, instance id, to VM id
//vector<vector<unsigned int>> VNF_TO_VM = {
//        {1,  1, 3},{1,  2, 1},{1,  3, 12},
//        {2,  1, 4},{2,  2, 2},
//        {3,  1, 1},{3,  2, 7},
//        {4,  1, 2},{4,  2, 13},{4,  3, 12},
//        {5,  1, 9},{5,  2, 3},
//        {6,  1, 7},{6,  2, 5},
//        {7,  1, 8},{7,  2, 6},
//        {8,  1, 10},{8,  2, 14},
//        {9,  1, 5},{9,  2, 9},
//        {10, 1, 11},{10, 2, 11}
//};
//assign_VNF_2_VM<type_res_l>(VNF_TO_VM, VNFNetwork, VirtualNetwork);
//
//vector<pair<int, int>> VM_TO_PN = {
//        {1,  1},{2,  1},{3,  2},{4,  2},{5,  3},{6,  3},{14, 3},{7,  4},{8,  4},{11, 5},{12, 5},{13, 5},{9,  6},{10, 6}
//};
//assign_VM_2_PN<type_wgt_l, type_res_l>(VM_TO_PN, VirtualNetwork, PhysicalNetwork);
///// for SFC -> Mapping of its VNF type to its instance id
//vector<vector<pair<int,int>>> VNFType_TO_InstID(SFCs.size());
////        VNFType_TO_InstID[0] = { {1, 1}, {6, 2}, {4, 1}, {5, 1} };
////        VNFType_TO_InstID[1] = { {10, 2}, {4, 2}, {8, 1}, {3, 1}, {2, 1}, {1, 1}, {5, 2} };
////        VNFType_TO_InstID[2] = { {4, 3}, {8, 2}, {7, 1}, {10, 1}, {9, 1} };
////    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[1], SFCs[1], VNFNetwork);
////    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[2], SFCs[2], VNFNetwork);
////    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[3], SFCs[3], VNFNetwork);
//} ///< assignment


//{ ///< assignment NEWYORK
///// Determination of how many instances of particular VNFs are required.
//vector<pair<unsigned int, unsigned int>> VNF_TO_InstancesCnt = {
//        {1, 4},{2, 3},{3, 5},{4, 4},{5, 6},{6, 5},{7, 4},{8, 3},{9, 4}
//};
//assign_VNF_2_InstancesCnt<type_res_l>(VNF_TO_InstancesCnt, VNFNetwork);
//
///// mapping of VNF id, instance id, to VM id
//vector<vector<unsigned int>> VNF_TO_VM = {
//        {1, 1, 5}, {1, 2, 4}, {1, 3, 12}, {1, 4, 11}, {2, 1, 13}, {2, 2, 9}, {2, 3, 17}, {3, 1, 2}, {3, 2, 15}, {3, 3, 9}, {3, 4, 10}, {3, 5, 11}, {4, 1, 5}, {4, 2, 7}, {4, 3, 8}, {4, 4, 1}, {5, 1, 5}, {5, 2, 6}, {5, 3, 7}, {5, 4, 10}, {5, 5, 14}, {5, 6, 16}, {6, 1, 6}, {6, 2, 15}, {6, 3, 19}, {6, 4, 4}, {6, 5, 18}, {7, 1, 13}, {7, 2, 6}, {7, 3, 18}, {7, 4, 16}, {8, 1, 1}, {8, 2, 19}, {8, 3, 12}, {9, 1, 3}, {9, 2, 4}, {9, 3, 17}, {9, 4, 14}
//};
//assign_VNF_2_VM<type_res_l>(VNF_TO_VM, VNFNetwork, VirtualNetwork);
//
//vector<pair<int, int>> VM_TO_PN = {
//        {1, 1},{2, 7},{3, 15},{4, 1},{5, 7},{6, 15},{7, 2},{8, 3},{9, 4},{10, 5},{11, 6},{12, 8},{13, 9},{14, 10},{15, 11},{16, 12},{17, 13},{18, 14},{19, 16}
//};
//assign_VM_2_PN<type_wgt_l, type_res_l>(VM_TO_PN, VirtualNetwork, PhysicalNetwork);
///// for SFC -> Mapping of its VNF type to its instance id
////        vector<vector<pair<int,int>>> VNFType_TO_InstID(SFCs.size());
////        VNFType_TO_InstID[0] = { {1, 1}, {6, 2}, {4, 1}, {5, 1} };
////        VNFType_TO_InstID[1] = { {10, 2}, {4, 2}, {8, 1}, {3, 1}, {2, 1}, {1, 1}, {5, 2} };
////        VNFType_TO_InstID[2] = { {4, 3}, {8, 2}, {7, 1}, {10, 1}, {9, 1} };
////    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[1], SFCs[1], VNFNetwork);
////    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[2], SFCs[2], VNFNetwork);
////    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[3], SFCs[3], VNFNetwork);
//} ///< assignmen