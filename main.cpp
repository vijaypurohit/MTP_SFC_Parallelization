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

// TODO: Finding number of VNF instances, and then VNF Deployement
// TODO: parallelism -> exploring more partial sfc possibilities.
// TODO: existing heuritics compare with bruteforce
// TODO: can you parallelise some aspects of program??
// TODO: convert seq to par
// TODO: vnfseq in delay calculation, convertion from seq to par.
int main()
{
    using type_wgt_l = type_delay; ///< Determine Edge Weights Data TYPE (unsigned int or FLOAT or DOUBLE)
    using type_res_l = unsigned int; ///< Determine Resource Data TYPE (unsigned int or FLOAT)

    // Object Creations
    PhysicalGraph<type_wgt_l, type_res_l> *PhysicalNetwork; ///< Graph Object contains network related functions
    VirtualMachines<type_res_l> *VirtualNetwork; ///< Virtual Machines Object contains Virtual Machines related functions
    VirtualNetworkFunctions<type_res_l> *VNFNetwork; ///< VNF Object contains VNF (Virtual Network Function) related code.
    vector<ServiceFunctionChain*> SFCs; ///< SFCs object contains code related to Service function chains

//    auto ft_start = std::chrono::steady_clock::now();
//    clusterSizeEnumeration(int(clusterSz.size()),maparVNFs_Cluster_AssignemtnxSFClen);
//    all_nCk(int(nCk.size()), maxSFClen);
//    vector<vector<int>> parSFC_Full = {{1}, {2,3,4,5}};
//    parVNFs_Cluster_Assignment(5, parSFC_Full, false);
//    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";

    unsigned int testIdx = 0; ///< Test Directory Initialisation.
    string testName = "sample" + to_string(testIdx); ///< should be without space
    string testDirName = testName + "/"; ///< directory name with slash at the end.
//    createDirectory(testDirName);

    try{  readNetwork<type_wgt_l, type_res_l>(testDirName, &PhysicalNetwork);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }

    try{   readVirtualMachines<type_res_l>(testDirName, &VirtualNetwork);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }

    try{  readVirtualNetworkFunctions<type_res_l>(testDirName, &VNFNetwork);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }

    try{  readGenericServiceFunctionsChains(testDirName, SFCs);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }
//    PhysicalNetwork->showAdjMatrix();
//    PhysicalNetwork->showAdjList();
//    PhysicalNetwork->showPNs_Description();
//    VirtualNetwork->showVMs_Description();
//    VNFNetwork->showVNFs_Description();
//    for(int ni=1; ni<SFCs.size(); ni++) SFCs[ni]->showSequentialSFC();
//    for(int ni=1; ni<SFCs.size(); ni++) SFCs[ni]->showParallelSFC(SFCs[ni]->vnfBlocksPar)

    { ///< assignment
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
        /// for SFC -> Mapping of its VNF type to its instance id
        vector<vector<pair<int,int>>> VNFType_TO_InstID(SFCs.size());
        VNFType_TO_InstID[1] = { {1, 1}, {6, 2}, {4, 1}, {5, 1} };
        VNFType_TO_InstID[2] = { {10, 2}, {4, 2}, {8, 1}, {3, 1}, {2, 1}, {1, 1}, {5, 2} };
        VNFType_TO_InstID[3] = { {4, 3}, {8, 2}, {7, 1}, {10, 1}, {9, 1} };
//    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[1], SFCs[1], VNFNetwork);
//    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[2], SFCs[2], VNFNetwork);
//    assign_ForSFC_VNFType_2_InstID<type_res_l>(VNFType_TO_InstID[3], SFCs[3], VNFNetwork);
    } ///< assignment

    for(int ni=1; ni<SFCs.size(); ni++) convert_SeqSFC2ParVNFBlocks<type_res_l>(SFCs[ni], VNFNetwork, (ni==SFCs.size()-1));
    for(int ni=1; ni<SFCs.size(); ni++) assign_Clusters2ParVNFs(SFCs[ni]);
    //    find_PartialParalleSFCs( true);

    vector<ServiceFunctionChain*> sortedSFCs;
    priority_queue<ServiceFunctionChain* , vector<ServiceFunctionChain*>, comparator_sfc> pqSortSFC(SFCs.begin()+1, SFCs.end()); ///< sfc[0] is not mapped to any valid sfc, otherwise sort the SFCs
    while(!pqSortSFC.empty()){
        sortedSFCs.push_back(pqSortSFC.top());
        pqSortSFC.pop();
    }

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

    auto ft_start = std::chrono::steady_clock::now();
    algo_LayerGraph_InstanceMapping<type_wgt_l, type_res_l>(sortedSFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork);
    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";


    showAlgoResults(SFCs, res_layerg);

    VNFNetwork->seq_utilization.clear();
    VNFNetwork->par_utilization.clear();

    ft_start = std::chrono::steady_clock::now();
    algo_PartialChains_InstanceMapping<type_wgt_l, type_res_l>(sortedSFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork);
    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";

    showAlgoResults(SFCs, res_partial);



    cout<<"\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    for(int ni=1; ni<SFCs.size(); ni++){
        delete SFCs[ni];
        if (debug) {
            if(ni == 1) cout << "\n[ SFCs Destructor Completed for sfc[" << ni<<"] ";
            else if(ni == SFCs.size()-1) cout << "sfc[" << ni<<"] ]";
            else cout << "sfc[" << ni<<"] ";
        }
    }
    delete VNFNetwork;
    delete VirtualNetwork;
    delete PhysicalNetwork;

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
//}