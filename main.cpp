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
//#include <limits> // numeric_limits
//#include <chrono>       //to seed random number
#include <random> //to generate random number
#include <functional>  // lambda function

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
#define maxSFCLength 10
#define maxVNF_Instances 5
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

#include "PhysicalGraph.h" // graph structure for physical network, physical node and edge and node capacity
#include "VirtualMachines.h" //  structure for virtual machines, virtual node
#include "VirtualNetworkFunctions.h" // structure for VNFnode, VNF,
#include "ServiceFunctionChain.h" // structure SFC
#include "FileReadingFunctions.h"
#include "TimeCalculationFunctions.h"
#include "AssignmentFunctions.h"
#include "Algorithms.h"

/*!
* @tparam type_wgt edge weight data type. default=unsigned int.
* @tparam type_res resource data type. default=unsigned int.
 * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
 * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
 * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
 * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] } \n
 *
*/
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
void layerGraphConstruction_and_InstanceSelectionAndRouting(ServiceFunctionChain *SFC,
                                                            VirtualNetworkFunctions<type_res> *VNFNetwork,
                                                            const VirtualMachines<type_res> *VirtualNetwork,
                                                            const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork, bool showInConsole) {

    /*! {{1},{6,4},{5}}; Each Partial SFC is without src and dest block/stage. */
    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}} ;
    unsigned int szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    /*! level to Instances Combinations = set of instance combination in block/stage/level index j. partParSFC = {{1},{6,4},{5}}  \n
     * stg 0 (1 function has 3 instances), B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  }\n
     * stg 1(2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] }\n
     * stg 2(1 function 2 instances), B[2] = {[5a] [5b] [5c]}\n
     * Time to calculte stg2InstCombinations -> if in any block number of parallel functions are 10 and each have\n
        inst = 2 (exe time: 1-2ms) (possibilites: 1024 (2^10)) \n
        inst = 3 (exe time: 38-40ms) (possibilites: 59049 (3^10))\n
        inst = 4 (exe time: 580-600ms) (possibilites: 10 48576 )\n
        inst = 5 (exe time: 5700-5800ms) (possibilites: 197 65625)\n
     */
    unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
    std::function<void(unsigned int, vector<pair<unsigned int,unsigned int>>&, unsigned int&, const vector<unsigned int>&)> findBofGivenBlk = [&findBofGivenBlk, &stg2InstCombinations, &VNFNetwork]
            (unsigned int bi, vector<pair<unsigned int,unsigned int>>& curInstComb, unsigned int& stgid, const vector<unsigned int>& curStg)->void{
        if(bi == curStg.size()){ // all functions in stage iterated.
            stg2InstCombinations[stgid].push_back(curInstComb); // push the one answer into combination stg.
            return;
        }
        unsigned int totInstancs = VNFNetwork->VNFNode[curStg[bi]]->numInstances;
        for(int instid=1; instid<=totInstancs; instid++){
            curInstComb.emplace_back(curStg[bi], instid); // push current instance
            findBofGivenBlk(bi+1, curInstComb, stgid, curStg); // call function for next instance
            curInstComb.pop_back(); // pop curInstComb instance and push next instance of same function.
        }
    };

    for(unsigned int stgid=0; stgid<szStages; stgid++){
        const auto& curStg = partParSFC[stgid];

        if(curStg.size() == 1){ // only one function in the block, then insert all its instance as combination
            for(int instid=1; instid<=VNFNetwork->VNFNode[curStg[0]]->numInstances; instid++)
                stg2InstCombinations[stgid].push_back({{curStg[0], instid}});
        } else{
            vector<pair<unsigned int,unsigned int>> curInstComb;
            findBofGivenBlk(0, curInstComb, stgid, curStg);
        }

    }

//    for(const auto& x: stg2InstCombinations){
//        cout<<"\nSTG["<<x.first<<"]("<<x.second.size()<<") { ";
//        for(const auto& y: x.second){
//            cout<<"[";
//            for(const auto& z: y){
//                cout<<""<<z.first<<char(z.second-1+'a')<<" ";
//            }
//            cout<<"]";
//        }
//        cout<<" }";
//    }

//    std::function<void()> instancesSelection =[&instancesSelection, &partParSFC, &stg2InstCombinations](unsigned int stgid, unordered_map<unsigned int, int>& curVNFType2Inst)->void{
//        if(stgid == partParSFC.size())
//        {
//            return;
//        }
//
//        for(const auto& instComb: stg2InstCombinations[stgid]){
//            for(const auto& [fnType, fnInstId]: instComb){
//                curVNFType2Inst[fnType]=fnType;
//            }
//        }
//
//    };

    struct layerGraphNode{
        pair<unsigned int,unsigned int> lev2InstCombId;
        int time=0;
        layerGraphNode(int _updTime, pair<unsigned int,unsigned int> _givenInstComb){
            this->time = _updTime;
            this->lev2InstCombId = _givenInstComb;
        };
    };
    layerGraphNode srcNode(0,{0,0});
    unordered_map<int,vector<layerGraphNode>> lgAdj;
    lgAdj[0].push_back(srcNode);

    unsigned int stg_nxt_idx = 1;
    /*!
     * stg 0 (1 function has 3 instances), B[0] = 2d{  1d[<1a>] [<1b>] [<1c>]  }\n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] }\n
     * stg 2 (1 function 2 instances), B[2] = {[5a] [5b] [5c}\n
     */

//     for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){
//         cout<<"\nSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
//         for(const auto& instComb: stg2InstCombinations[cur_lvl]){
//             cout<<"[";
//             for(const auto& givenPair: instComb){
//                 cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" ";
//             }
//             cout<<"]";
//         }
//     }


}



int main()
{
// TODO: Finding number of VNF instances, and then VNF Deployement
// TODO: existing heuritics compare with bruteforce
// TODO: try to remove sAdj /pAdj as they are mainly not required.
// TODO: can you parallelise some aspects of program??
// TODO:
// TODO: parallel graph data collection

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
    for(int ni=1; ni<=total_SFC; ni++) convertSeqSFC_to_FullParVNFBlocks<type_res_local>(SFC[ni], VNFNetwork);
//    SFC[1]->convertToParallelSFC({{1},{6,4},{5}});
//    SFC[2]->convertToParallelSFC({{10},{4,8,3},{2,1},{5}});
//    SFC[3]->convertToParallelSFC({{4,8},{7,10,9}});
//    for(int ni=1; ni<=total_SFC; ni++) SFC[ni]->showSFC_BlockWise(SFCpar, SFC[ni]->I_VNFType2Inst);



    /// Determination of how many instances of particular VNFs are required.
    vector<pair<unsigned int,unsigned int>> VNF_TO_InstancesCnt = {
            {1, 3}, {4, 3},
            {2, 2}, {3, 2}, {5, 2}, {6, 2},  {7, 2}, {8, 2},  {9, 2}, {10, 2}
    };
    assign_VNF_2_InstancesCnt<type_res_local>(VNF_TO_InstancesCnt, VNFNetwork);

    /// mapping of VNF id, instance id, to VM id
    vector<vector<unsigned int>> VNF_TO_VM = {
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

    vector<pair<int,int>> VM_TO_PN = {
            {1, 1}, {2, 1},
            {3, 2}, {4, 2},
            {5, 3}, {6, 3}, {14, 3},
            {7, 4}, {8, 4},
            {11, 5}, {12, 5}, {13, 5},
            {9, 6}, {10, 6}
    };
    assign_VM_2_PN<type_wgt_local, type_res_local>(VM_TO_PN, VirtualNetwork, PhysicalNetwork);


    /// for SFC -> Mapping of its VNF type to its instance id
    vector<vector<pair<int,int>>> VNFType_TO_InstID(total_SFC+1);
    VNFType_TO_InstID[1] = { {1, 1}, {6, 2}, {4, 1}, {5, 1} };
    VNFType_TO_InstID[2] = { {10, 2}, {4, 2}, {8, 1}, {3, 1}, {2, 1}, {1, 1}, {5, 2} };
    VNFType_TO_InstID[3] = { {4, 3}, {8, 2}, {7, 1}, {10, 1}, {9, 1} };
    assign_ForSFC_VNFType_2_InstID(VNFType_TO_InstID[1], SFC[1]);
    assign_ForSFC_VNFType_2_InstID(VNFType_TO_InstID[2], SFC[2]);
    assign_ForSFC_VNFType_2_InstID(VNFType_TO_InstID[3], SFC[3]);

    auto ft_start = std::chrono::steady_clock::now();
//    parVNFBlocks_ClusterAssignment_ForSFC(SFC[1]);
    layerGraphConstruction_and_InstanceSelectionAndRouting(SFC[1], VNFNetwork, VirtualNetwork, PhysicalNetwork, true);
    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";

//    PhysicalNetwork->showPNs_Description();
//    VirtualNetwork->showVMs_Description();
//    VNFNetwork->showVNFs_Description();
//    for(int ni=1; ni<=total_SFC; ni++) SFC[ni]->showSFC_BlockWise(SFCseq);
//    for(int ni=1; ni<=total_SFC; ni++) SFC[ni]->showSFC_BlockWise(SFCpar);


//// collecting information for each VNF node that what are its instances and where are they hosted {vmid, pnid}
//    for (const auto& [NFid, NF]: VNFNetwork->VNFNode) {
////        cout<<"\n F["<<NFid<<"] -> { ";
//        for(unsigned int inst_id=1; inst_id<=NF->numInstances; inst_id++){
//            int vm_id = VNFNetwork->I_VNFinst2VM[NFid][inst_id];
//            int pn_id = VirtualNetwork->I_VM2PN[vm_id];
//            NF->inst2nw[inst_id] = {vm_id,pn_id};
////            cout<<"VM["<<vm_id<<"]"<<"PN["<<pn_id<<"] | ";
//        }
////        cout<<" }";
//    }
//    debug=0;
//    for(int ni=1; ni<=total_SFC; ni++){
//        cout<<"\nSFC["<<ni<<"] E2E: ";
//        cout<<"Seq["<<calcObjectiveValueSeq<type_wgt_local, type_res_local>(SFC[ni], SFC[ni]->I_VNFType2Inst, SFC, VNFNetwork, VirtualNetwork, PhysicalNetwork)<<"]  ";
//        cout<<"Par["<<calcObjectiveValuePar<type_wgt_local, type_res_local>(SFC[ni]->vnfBlocksPar,SFC[ni]->I_VNFType2Inst,  SFC[ni], SFC, VNFNetwork, VirtualNetwork, PhysicalNetwork)<<"]  ";
//        cout<<"Pkt["<<calcTime_PacketsDelay<type_res_local>(SFC[ni]->vnfBlocksPar, SFC[ni]->I_VNFType2Inst, SFC[ni], VNFNetwork, VirtualNetwork)<<"]";
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