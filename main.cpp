/*!
 * @author Vijay Purohit
  Created by vijay on 22-02-2023.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <cmath>
//#include <limits> // numeric_limits
//#include <random>       //to generate random number
//#include <chrono>       //to seed random number
//#include <sstream>
#include <random>
#include <functional>

using namespace std;

/**************** Variable Declaration ****************/
int debug = 1;
const string input_directory = "files_input/", output_directory = "files_output/", diagram_directory = "graphs/";
const string filename_constants   = "constants.txt";
const string filename_network   = "network.txt";
const string filename_virtualmachines   = "virtual_machines.txt";
const string filename_vnf   = "VNFs.txt";
const string filename_sfc   = "SFCs.txt", filename_sfc_parallel = "SFCs_Parallel.txt";
const string filename_vnf_parallelpairs = "VNFs_ParallelPairs.txt";
// graphviz colors list
/*const vector<vector<string>> graphviz_colors ={    {"brown", "firebrick4", "crimson", "maroon"},
                                                   {"darkorange", "tomato"},
                                                   {"indigo", "darkslateblue", "deepskyblue4", "darkturquoise"},
                                                   {"gold2"},
                                                   {"forestgreen"},
                                                   {"hotpink", "plum", "thistle"}
                                               };*/

/**************** Custom Header Files Declaration ****************/
#define SFCsrc 0
#define SFCdst (-10)
#define SFCseq (-11)
#define SFCpar (-12)
#define maxSFClen 10

/*!
 * @param factor_packet factor to multiply to convert packet size in bits. 1 Byte is 8 bits
 * @param packetBodySize, packetHeaderSize Size of the Network Packet Body and Header. in Bytes.
 *  Type = unsigned int Range[0,4294967295].
 *  @param factor_bandwidth factor to multiply to convert bandwidth in bits/seconds.
 * @param bandwidthNW Bandwidth of the Network. in Mega bits per second. 1Gb = 1000 Mb. Type = unsigned int Range[0,4294967295].\n
 *  Mb[0,4294967295], Gb_in_Mb[1000 , 4294967.295].
 * @param velocityFactor velocity factor of transmission medium. vaccum = 1.0. copper wise = 0.7
 * @param speedOfLight speed of light in vaccum 3 * 10^8 m/s
 * @param total_SFC total number of SFC in the network
 * @param read_write_time_per_bit 0.077ms (measured by duplicating a large file of 1 MB in a server with Intel i7-8700 core
 */
const unsigned int factor_packet = 8; unsigned int packetBodySize = 1000, packetHeaderSize = 24;
const unsigned int factor_bandwidth = 1000000; unsigned int bandwidthNW = 10;
float velocityFactor = 1.0;
unsigned int  speedOfLight = 300000000;
unsigned int total_SFC = 0;
float read_write_time_per_bit = 0.077e-3;

void debugPrint(const string& msg){
    if(debug) cout<<msg;
}
#include "PhysicalGraph.h" // graph structure for physical network, physical node and edge and node capacity
#include "VirtualMachines.h" //  structure for virtual machines, virtual node
#include "VirtualNetworkFunctions.h" // structure for VNFnode, VNF,
#include "ServiceFunctionChain.h" // structure SFC
#include "FileReadingFunctions.h"

//template<class T>
//typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
//almost_equal(T x, T y, int ulp)
//{
//    // the machine epsilon has to be scaled to the magnitude of the values used and multiplied by the desired precision in ULPs (units in the last place)
//    return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp
//           // unless the result is subnormal
//           || std::fabs(x - y) < std::numeric_limits<T>::min();
//}

#include "TimeCalculationFunctions.h"
#include "AssignmentFunctions.h"
#include "Algorithms.h"



//void dfs_dfs(int levelIdx, unordered_map<int,bool> visited, vector<vector<int>>blks, unordered_map<int,vector<vector<int>>>& levelNodes, vector<vector<vector<int>>>& allPCs){
//    if(levelIdx == levelNodes.size()) {
//        allPCs.push_back(blks);
//        return;
//    };
//    for(const auto& allComb: levelNodes[levelIdx]){
//        bool canBeVisited = true;
//        for(auto node: allComb){
//            if(visited[node]){
//                canBeVisited=false;
//                break;
//            }
//        }
//        if(canBeVisited){
//            for(auto node: allComb)visited[node]=true;
//            blks.push_back(allComb);
//            dfs_dfs(levelIdx+1, visited,blks, levelNodes, allPCs);
//            for(auto node: allComb)visited[node]=false;
//            blks.pop_back();
//        }
//    }
//
//}





int main()
{
// TODO: partial parallel chain?
// TODO: parallel graph data collection
// TODO: Finding number of VNF instances, VNF Deployement
// TODO: existing heuritics compare with bruteforce

//    auto ft_start = std::chrono::steady_clock::now();
//    clusterSizeEnumeration(int(clusterSz.size()),maparVNFs_Cluster_AssignemtnxSFClen);
//    all_nCk(int(nCk.size()), maxSFClen);
//    vector<vector<int>> parSFC_Full = {{1}, {2,3,4,5}};
//    parVNFs_Cluster_Assignment(5, parSFC_Full, false);
//    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";


    using type_wgt_local = float; ///< Determine Edge Weights Data TYPE (unsigned int or FLOAT)
    using type_res_local = unsigned int; ///< Determine Resource Data TYPE (unsigned int or FLOAT)

    unsigned int testIdx = 0; ///< Test Directory Initialisation.
    string testName = "sample" + to_string(testIdx); ///< should be without space
    string testDirName = testName + "/"; ///< directory name with slash at the end.
//    createDirectory(testDirName);

    PhysicalGraph<type_wgt_local, type_res_local> *PhysicalNetwork; ///< Network Creation
    try{  readNetwork<type_wgt_local, type_res_local>(testDirName, &PhysicalNetwork);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }
//    PhysicalNetwork->showAdjMatrix();
//    PhysicalNetwork->showAdjList();
//    PhysicalNetwork->printPhysicalGraph(testDirName);
//    PhysicalNetwork->showAllPairsShortestPath();
    VirtualMachines<type_res_local> *VirtualNetwork;
    try{   readVirtualMachines<type_res_local>(testDirName, &VirtualNetwork);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }

    VirtualNetworkFunctions<type_res_local> *VNFNetwork;
    try{  readVirtualNetworkFunctions<type_res_local>(testDirName, &VNFNetwork);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }

    vector<ServiceFunctionChain*> SFC;
    try{  readGenericServiceFunctionsChains(testDirName, SFC);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }

    for(int ni=1; ni<=total_SFC; ni++) convertSeqSFC_to_FullParallelSFC<type_res_local>(SFC[ni], VNFNetwork);
    if(debug)cout<<"\n\t[SFC converted to Full Parallelised Chain]. Parallel SFC File:"<<output_directory+testDirName+filename_sfc_parallel;
//    for(int ni=1; ni<=total_SFC; ni++) SFC[ni]->showSFC_BlockWise(SFCpar);
    auto ft_start = std::chrono::steady_clock::now();
    parVNFs_Cluster_Assignment(SFC[3], true);
    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";

//    SFC[1]->convertToParallelSFC({{SFCsrc},{1},{6,4},{5},{SFCdst}});
//    SFC[2]->convertToParallelSFC({{SFCsrc},{10},{4,8,3},{2,1},{5},{SFCdst}});
//    SFC[3]->convertToParallelSFC({{SFCsrc},{4,8},{7,10,9},{SFCdst}});
//    for(int ni=1; ni<=total_SFC; ni++) SFC[ni]->showAdjList(SFCpar);
    return 0;
    /// Mapping of VM id to PN id
    vector<pair<int,int>> VM_TO_PN = {
            {1, 1}, {2, 1},
            {3, 2}, {4, 2},
            {5, 3}, {6, 3}, {14, 3},
            {7, 4}, {8, 4},
            {11, 5}, {12, 5}, {13, 5},
            {9, 6}, {10, 6}
    };
    assign_VM_2_PN<type_wgt_local, type_res_local>(VM_TO_PN, VirtualNetwork, PhysicalNetwork);

    /// Determination of how many instances of particular VNFs are required.
    vector<pair<int,int>> VNF_TO_InstancesCnt = {
            {1, 3}, {4, 3},
            {2, 2}, {3, 2}, {5, 2}, {6, 2},  {7, 2}, {8, 2},  {9, 2}, {10, 2}
    };
    assign_VNF_2_InstancesCnt<type_res_local>(VNF_TO_InstancesCnt, VNFNetwork);

    /// mapping of VNF id, instance id, to VM id
    vector<vector<int>> VNF_TO_VM = {
            {1, 1, 3}, {1, 2, 1}, {1, 3, 12},
            {2, 1, 4}, {2, 2, 2},
            {3, 1, 1}, {3, 2, 7},
            {4, 1, 2}, {4, 2, 13}, {4, 3, 12},
            {5, 1, 9}, {5, 2, 3},
            {6, 1, 7}, {6, 2, 5},
            {7, 1, 8}, {7, 2, 6},
            {8, 1, 10}, {8, 2, 14},
            {9, 1, 5}, {9, 2, 9},
            {10, 1, 11}, {10, 2, 11}
    };
    assign_VNF_2_VM<type_res_local>(VNF_TO_VM, VNFNetwork, VirtualNetwork);

    /// for SFC -> Mapping of its VNF type to its instance id
    vector<vector<pair<int,int>>> VNFType_TO_InstID(total_SFC+1);
    VNFType_TO_InstID[1] = { {1, 1}, {6, 2}, {4, 1}, {5, 1} };
    VNFType_TO_InstID[2] = { {10, 2}, {4, 2}, {8, 1}, {3, 1}, {2, 1}, {1, 1}, {5, 2} };
    VNFType_TO_InstID[3] = { {4, 3}, {8, 2}, {7, 1}, {10, 1}, {9, 1} };
    assign_ForSFC_VNFType_2_InstID(VNFType_TO_InstID[1], SFC[1]);
    assign_ForSFC_VNFType_2_InstID(VNFType_TO_InstID[2], SFC[2]);
    assign_ForSFC_VNFType_2_InstID(VNFType_TO_InstID[3], SFC[3]);

//    PhysicalNetwork->showPNs_Description();
//    VirtualNetwork->showVMs_Description();
//    VNFNetwork->showVNFs_Description();
//    for(int ni=1; ni<=total_SFC; ni++) SFC[ni]->showSFC_BlockWise(SFCseq);
//    for(int ni=1; ni<=total_SFC; ni++) SFC[ni]->showSFC_BlockWise(SFCpar);

//    debug=0;
//    for(int ni=1; ni<=total_SFC; ni++){
//        calcObjectiveValueSeq<type_wgt_local, type_res_local>(SFC[ni], SFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);
//        calcObjectiveValuePar<type_wgt_local, type_res_local>(SFC[ni], SFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);
//        calcTime_PacketsDelay<type_res_local>(SFC[ni], VNFNetwork, VirtualNetwork);
//        cout<<"\nSFC["<<ni<<"] E2E: Seq["<<SFC[ni]->distanceSeq[SFCdst]<<"]  Par["<<SFC[ni]->distancePar[SFCdst]<<"]  Pkt["<<SFC[ni]->pktDist[SFC[ni]->vnfBlocksPar.size()-1][SFCdst]<<"]";
//    }
//    debug=1;





//////////////////////////////////////////////calling objective function/////////////////////////////////


//    float objectiveValue = calcObjectiveValuePar<type_wgt_local, type_res_local>(SFC[3], SFC, VNFNetwork, VirtualNetwork, PhysicalNetwork, true);
//    cout << objectiveValue << endl;


    cout<<"\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    for(int ni=1; ni<=total_SFC; ni++) delete SFC[ni];
    delete VNFNetwork;
    delete VirtualNetwork;
    delete PhysicalNetwork;

    return 0;
}