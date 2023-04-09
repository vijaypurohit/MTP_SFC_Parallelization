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
//#include <execution>
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
void layerGraphConstruction_and_InstanceSelectionAndRouting(ServiceFunctionChain *cSFC, vector<ServiceFunctionChain*> allSFC,
                                                            const VirtualNetworkFunctions<type_res> *VNFNetwork,
                                                            const VirtualMachines<type_res> *VirtualNetwork,
                                                            const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {


//    unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > bst_stg2InstCombinations;
//    bst_stg2InstCombinations = {
//            {0, {      {{1, 1}}, {{1,2}}, {{1,3}}     }        },
//            {1, {      {{6, 1}, {4, 1}}, {{6,1}, {4, 2}}, {{6, 2}, {4, 1}}, {{6,2}, {4, 2}}    }        },
//            {2, {      {{5, 1}}, {{5,2}}     }        }
//    };
//    for(unsigned int stgid=0; stgid<bst_stg2InstCombinations.size(); stgid++) {
//
//        for (const vector<pair<unsigned int, unsigned int>> &instComb: bst_stg2InstCombinations[stgid]) {
//            bool instCanBePicked_UtilizationLess = true;
//            for (const auto &[fnType, fnInstId]: instComb) {
//                // is function is safe to pick ?
//                if (VNFNetwork->utilization.count(fnType) and VNFNetwork->utilization.at(fnType).count(fnInstId) and
//                    ///< old utilization till now of VNF
//                    (VNFNetwork->utilization.at(fnType).at(fnInstId) + cSFC->trafficArrivalRate >
//                     VNFNetwork->VNFNodes.at(fnType)->serviceRate)) {
//                    instCanBePicked_UtilizationLess = false;
//                    break;
//                }
//                curMapping[fnType] = fnInstId;
//            }
//            if (!instCanBePicked_UtilizationLess)continue;
//            instancesEnumeration_And_FindMinTime(stgid + 1, curMapping, curBstMapping, minBstTime, partParSFC,
//                                                 stg2InstCombinations);
//        }
//
//    }//


    /*! Lambda function to find all instances combination of parVNF in that stage. find_stg2IC_ofGivenBlk -> find stage to instnace combination of given block of sfc.
     * @param csfi current stage function index.
     * @param curInstComb current combination in iteration
     * @param stgid stgId/blockId for which we are finding combination of functions in that stg/block and to store in stg2InstCombinations.
     * @param curStg using to iterate all functions in the stage.
     * @param stg2InstCombinations It stores all the stage wise instances combination of all stage in partParSFC. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }}
     * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
     * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
     * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] } \n
     * Time to calculte stg2InstCombinations -> if in any block number of parallel functions are 10 and each have 5 max instances\n
        inst = 2 (exe time: 1-2ms) (possibilites: 1024 (2^10)) \n
        inst = 3 (exe time: 38-40ms) (possibilites: 59 049 (3^10))\n
        inst = 4 (exe time: 580-600ms) (possibilites: 10 48 576 )\n
        inst = 5 (exe time: 5700-5800ms) (possibilites: 97 65 625)\
     */
    std::function<void(unsigned int, vector<pair<unsigned int,unsigned int>>&, unsigned int&, const vector<unsigned int>&, unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>>>&)>
            find_stg2IC_ofGivenBlk = [&find_stg2IC_ofGivenBlk,&cSFC,&VNFNetwork] (unsigned int csfi, vector<pair<unsigned int,unsigned int>>& curInstComb,   unsigned int& stgid,
                                                                             const vector<unsigned int>& curStg, unordered_map<unsigned int,  vector<vector<pair<unsigned int,unsigned int>>>> &stg2InstCombinations)->void{
                if(csfi == curStg.size()){ // all functions in stage iterated. curStg.size()==numOfFunction in that stage.
                        stg2InstCombinations[stgid].push_back(curInstComb); // push the one answer into combination stg.
                    return;
                }
                const unsigned int fnType = curStg[csfi];
                unsigned int totInstancs = VNFNetwork->VNFNodes.at(fnType)->numInstances;
                for(unsigned int fnInstId=1; fnInstId<=totInstancs; fnInstId++){
                    if (VNFNetwork->utilization.count(fnType) and VNFNetwork->utilization.at(fnType).count(fnInstId) and ///< old utilization till now of VNF
                        (VNFNetwork->utilization.at(fnType).at(fnInstId) + cSFC->trafficArrivalRate >  VNFNetwork->VNFNodes.at(fnType)->serviceRate))  {
                        continue; // don't take this instance if its utilisation become more than service rate of function.
                    }
                    curInstComb.emplace_back(fnType, fnInstId); // push current instance
                    find_stg2IC_ofGivenBlk(csfi+1, curInstComb, stgid, curStg, stg2InstCombinations); // call function for next instance
                    curInstComb.pop_back(); // pop curInstComb instance and push next instance of same function.
                }

    };

    /*! For a given stg2InstCombinations, it enumerate all the possible mappings we can give in each stage and calculate time on the go.
     * @param stgid stgId/blockId for which we are enumerating instances.
     * @param curMapping cur function->instance mapping we iterating out of all possibilites.
     * @param bstMapping to save best mapping overall among all partial parallel chain/and its all instances.
     * @param minBstTime min time among all partial parallel chain/and its all instances.
     * @param partParSFC given partial SFC
     * @param stg2InstCombinations It contains all the stage wise instances combination of all stage in partParSFC. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }}
     * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
     * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
     * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] }
     * allMappings are (total 36 = 3*6*2) \n
        0[1a 6a 4a 5a ]         1[1a 6a 4a 5b ]         2[1a 6a 4b 5a ]     3[1a 6a 4b 5b ]         4[1a 6a 4c 5a ]         5[1a 6a 4c 5b ]
        6[1a 6b 4a 5a ]         7[1a 6b 4a 5b ]         8[1a 6b 4b 5a ]     9[1a 6b 4b 5b ]         10[1a 6b 4c 5a ]        11[1a 6b 4c 5b ]
        12[1b 6a 4a 5a ]        13[1b 6a 4a 5b ]        14[1b 6a 4b 5a ]        15[1b 6a 4b 5b ]        16[1b 6a 4c 5a ]        17[1b 6a 4c 5b ]
        18[1b 6b 4a 5a ]        19[1b 6b 4a 5b ]        20[1b 6b 4b 5a ]        21[1b 6b 4b 5b ]        22[1b 6b 4c 5a ]        23[1b 6b 4c 5b ]
        24[1c 6a 4a 5a ]        25[1c 6a 4a 5b ]        26[1c 6a 4b 5a ]        27[1c 6a 4b 5b ]        28[1c 6a 4c 5a ]        29[1c 6a 4c 5b ]
        30[1c 6b 4a 5a ]        31[1c 6b 4a 5b ]        32[1c 6b 4b 5a ]
        3[1c 6b 4b 5b ]        34[1c 6b 4c 5a ]        35[1c 6b 4c 5b ] \n
     * 1 stage -> 10 parallel func each with 5 max instances -> 5^10 possibilities or 97,65,625 instances.
     */
    std::function<void(unsigned int , unordered_map<unsigned int,unsigned int>& , unordered_map<unsigned int,unsigned int>&, float& , const vector<vector<unsigned int>>&,  unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> >&)>
        instancesEnumeration_And_FindMinTime =[&instancesEnumeration_And_FindMinTime, &cSFC, &allSFC, &VNFNetwork, &VirtualNetwork, &PhysicalNetwork, &showInConsole, &showInConsoleDetailed]
                (unsigned int stgid, unordered_map<unsigned int,unsigned int>& curMapping, unordered_map<unsigned int,unsigned int>& curBstMapping,
                float& minBstTime, const vector<vector<unsigned int>>& partParSFC,  unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> >& stg2InstCombinations)->void{
            if(stgid == partParSFC.size()) { // found one mapping then find corresponding time
                        float parallelCost=0, packetCost=0;
                    parallelCost = calcObjectiveValuePar<type_wgt, type_res>(partParSFC,curMapping,  cSFC->index, allSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);
                    // packet delay only in case of parallelism when number of blocks/stages < numOfVNFs (Sequential chain length) AND cost is less than min time among  partial sfc all instances.
                    if( partParSFC.size() < cSFC->numVNF and parallelCost < minBstTime)
                        packetCost = calcTime_PacketsDelay<type_res>(partParSFC, curMapping,  cSFC->index, VNFNetwork, VirtualNetwork);
                   float curMapTime = parallelCost+packetCost; // float_max+flaot+max => inf

                    if( curMapTime < minBstTime){ // current mapping ka time is less than min time among partial sfc all instances.
                        minBstTime =  curMapTime;
                        curBstMapping = curMapping;
                    }
//                    else
//                        return;
                    if(showInConsole and showInConsoleDetailed){
                        cout<<"\n\t"<<"["; for(const auto &blk: partParSFC){ for(const auto& fnid: blk){ cout<<fnid<<char(96+curMapping.at(fnid))<<" ";  }  } cout<<"]";
                        cout<<"["<<curMapTime<<"sec ("<<parallelCost<<"|"<<packetCost<<")]";
                    }
                    return;
            }
            for(const vector<pair<unsigned int,unsigned int>>& instComb: stg2InstCombinations[stgid]){
                bool instCanBePicked_UtilizationLess = true;
                for(const auto& [fnType, fnInstId]: instComb) {
                    // is function is safe to pick ?
                    if (VNFNetwork->utilization.count(fnType) and VNFNetwork->utilization.at(fnType).count(fnInstId) and ///< old utilization till now of VNF
                            (VNFNetwork->utilization.at(fnType).at(fnInstId) + cSFC->trafficArrivalRate >  VNFNetwork->VNFNodes.at(fnType)->serviceRate))  {
                        instCanBePicked_UtilizationLess = false; break;
                    }
                    curMapping[fnType] = fnInstId;
                }
                if(!instCanBePicked_UtilizationLess)continue;
                instancesEnumeration_And_FindMinTime(stgid+1,curMapping, curBstMapping, minBstTime, partParSFC, stg2InstCombinations);
            }
    };

    cSFC->bst_parlen_idx = cSFC->allPartParSFC.size(); // from partParSFC what is the best partial mapping.
    cSFC->bst_seqlen_time = std::numeric_limits<float>::max();
    cSFC->bst_parlen_time = std::numeric_limits<float>::max();
    float bst_partpar_time = std::numeric_limits<float>::max();
    //    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}} ;
    for(int ppsidx=int( cSFC->allPartParSFC.size())-1; ppsidx>=0; --ppsidx){
    /*! {{1},{6,4},{5}}; Each Partial SFC is without src and dest block/stage. */
        const vector<vector<unsigned int>>& partParSFC=cSFC->allPartParSFC.at(ppsidx); ///< for each of the partial parallel SFC of the givenParVNF Blocks
        const unsigned int szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }} */
        unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
        for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
            const auto& curStg = partParSFC[stgid];
                vector<pair<unsigned int,unsigned int>> curInstComb;
                find_stg2IC_ofGivenBlk(0, curInstComb, stgid, curStg, stg2InstCombinations);
        }//stgid<szStages finding instances possibilities of each stage.
        // show stages wise instances combination

        if(showInConsole and showInConsoleDetailed){ // showing partParSFC info
            cout<<"\n partParSFC["<<ppsidx<<"]: "; for(const auto& blks: partParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<") ---------- - --------- - ------";
            for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
                cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
                for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                    cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
                }  cout<<" }";
            }

        }

    /*! finding all the mapping possibilites for the current partParSFC instance combination at each stage.*/
        unordered_map<unsigned int,unsigned int> curMapping; ///< iterating mapping variable
        if(ppsidx == 0){// index 0 partial chain is same as given sfc, so finding best instance combinaton for given sequential sfc.
            instancesEnumeration_And_FindMinTime(0, curMapping,  cSFC->bst_seqlen_mapping, cSFC->bst_seqlen_time, partParSFC ,stg2InstCombinations);
            if(cSFC->bst_seqlen_time < cSFC->bst_parlen_time){ // if sequential is better than parallel than choose sequential.
                cSFC->bst_parlen_time = cSFC->bst_seqlen_time;
                cSFC->bst_parlen_idx = ppsidx;
                cSFC->bst_parlen_mapping = cSFC->bst_seqlen_mapping;
            }
        }
        else{
            instancesEnumeration_And_FindMinTime(0, curMapping,  cSFC->bst_parlen_mapping, bst_partpar_time, partParSFC, stg2InstCombinations);
            if(bst_partpar_time < cSFC->bst_parlen_time){
                cSFC->bst_parlen_idx = ppsidx;
                cSFC->bst_parlen_time = bst_partpar_time;
            }
//            unordered_map<unsigned int,unsigned int> curBstMapping; ///< in a single partPar SFC what is the best mapping which gives minimum time
//            float minCurTime = std::numeric_limits<float>::max(); ///< in a single partPar SFC what is the best minimum time
//            instancesEnumeration_And_FindMinTime(0, curMapping,  curBstMapping, minCurTime, partParSFC, stg2InstCombinations);
//            if(minCurTime < cSFC->bst_parlen_time){
//                cSFC->bst_parlen_idx = ppsidx;
//                cSFC->bst_parlen_time = minCurTime;
//                cSFC->bst_parlen_mapping = curBstMapping;
////                bst_stg2InstCombinations = stg2InstCombinations;
//            }
        }

    }// for each PartParSFC.


    if(cSFC->bst_parlen_idx == cSFC->allPartParSFC.size() ){
        string errorMsg = "Algorithm failed to find best partial parallel mapping for SFC["+to_string(cSFC->index)+ "]. Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    if(cSFC->bst_seqlen_time == std::numeric_limits<float>::max()){
        string errorMsg = "Algorithm failed to find best sequential mapping for SFC["+to_string(cSFC->index)+ "]  Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    if(showInConsole){
        cout<<"\n\n Objective Function Answer For SFC["<<cSFC->index<<"]";
        cout<<"\n\t Sequential: partIdx[0]  time:["<<cSFC->bst_seqlen_time<<"] :(";
        for(const auto &blk: cSFC->allPartParSFC[0]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout<<fnid<<char(96+cSFC->bst_seqlen_mapping.at(fnid))<<" ";
            }   cout<<"]";
        }
        cout<<"\n\t Parallel: partIdx["<<cSFC->bst_parlen_idx<<"]  time:["<<cSFC->bst_parlen_time<<"] :(";
        for(const auto &blk: cSFC->allPartParSFC[cSFC->bst_parlen_idx]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout<<fnid<<char(96+cSFC->bst_parlen_mapping.at(fnid))<<" ";
            }   cout<<"]";
        }

    }


//    for(unsigned int stgid=0; stgid<bst_stg2InstCombinations.size(); stgid++){      // finding instances possibilities of each stage.
//        const auto& curStg = partParSFC[stgid];
//
//    }//stgid<szStages finding instances possibilities of each stage.


}

template<typename type_wgt=unsigned int, typename type_res=unsigned int>
void layerGraphConstruction_and_InstanceSelectionAndRouting1stVersion(ServiceFunctionChain *cSFC, vector<ServiceFunctionChain*> allSFC,
                                                            VirtualNetworkFunctions<type_res> *VNFNetwork,
                                                            const VirtualMachines<type_res> *VirtualNetwork,
                                                            const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {

    /*! Lambda function to find all instances combination of parVNF in that stage.
     * @param csfi current stage function index.
     * @param curInstComb current combination in iteration
     * @param stgid stgId/blockId for which we are finding combination of functions in that stg/block and to store in stg2InstCombinations.
     * @param curStg using to iterate all functions in the stage.
     * @param stg2InstCombinations It stores all the stage wise instances combination of all stage in partParSFC. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }}
     * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
     * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
     * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] } \n
     * Time to calculte stg2InstCombinations -> if in any block number of parallel functions are 10 and each have 5 max instances\n
        inst = 2 (exe time: 1-2ms) (possibilites: 1024 (2^10)) \n
        inst = 3 (exe time: 38-40ms) (possibilites: 59 049 (3^10))\n
        inst = 4 (exe time: 580-600ms) (possibilites: 10 48 576 )\n
        inst = 5 (exe time: 5700-5800ms) (possibilites: 97 65 625)\
     */
    std::function<void(unsigned int, vector<pair<unsigned int,unsigned int>>&, unsigned int&, const vector<unsigned int>&, unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>>>&)>
            find_stg2IC_ofGivenBlk = [&find_stg2IC_ofGivenBlk, &VNFNetwork] (unsigned int csfi, vector<pair<unsigned int,unsigned int>>& curInstComb,   unsigned int& stgid,
                    const vector<unsigned int>& curStg, unordered_map<unsigned int,  vector<vector<pair<unsigned int,unsigned int>>>> &stg2InstCombinations)->void{
                if(csfi == curStg.size()){ // all functions in stage iterated. curStg.size()==numOfFunction in that stage.
                    stg2InstCombinations[stgid].push_back(curInstComb); // push the one answer into combination stg.
                    return;
                }
                unsigned int totInstancs = VNFNetwork->VNFNodes[curStg.at(csfi)]->numInstances;
                for(unsigned int instid=1; instid<=totInstancs; instid++){
                    curInstComb.emplace_back(curStg.at(csfi), instid); // push current instance
                    find_stg2IC_ofGivenBlk(csfi+1, curInstComb, stgid, curStg, stg2InstCombinations); // call function for next instance
                    curInstComb.pop_back(); // pop curInstComb instance and push next instance of same function.
                }
    };
    /*! For a given stg2InstCombinations, it enumerate all the possible mappings we can give in each stage.
     * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
     * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
     * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] }
     * allMappings are (total 36 = 3*6*2) \n
        0[1a 6a 4a 5a ]         1[1a 6a 4a 5b ]         2[1a 6a 4b 5a ]     3[1a 6a 4b 5b ]         4[1a 6a 4c 5a ]         5[1a 6a 4c 5b ]
        6[1a 6b 4a 5a ]         7[1a 6b 4a 5b ]         8[1a 6b 4b 5a ]     9[1a 6b 4b 5b ]         10[1a 6b 4c 5a ]        11[1a 6b 4c 5b ]
        12[1b 6a 4a 5a ]        13[1b 6a 4a 5b ]        14[1b 6a 4b 5a ]        15[1b 6a 4b 5b ]        16[1b 6a 4c 5a ]        17[1b 6a 4c 5b ]
        18[1b 6b 4a 5a ]        19[1b 6b 4a 5b ]        20[1b 6b 4b 5a ]        21[1b 6b 4b 5b ]        22[1b 6b 4c 5a ]        23[1b 6b 4c 5b ]
        24[1c 6a 4a 5a ]        25[1c 6a 4a 5b ]        26[1c 6a 4b 5a ]        27[1c 6a 4b 5b ]        28[1c 6a 4c 5a ]        29[1c 6a 4c 5b ]
        30[1c 6b 4a 5a ]        31[1c 6b 4a 5b ]        32[1c 6b 4b 5a ]
        3[1c 6b 4b 5b ]        34[1c 6b 4c 5a ]        35[1c 6b 4c 5b ] \n
     * 1 stage -> 10 parallel func each with 5 max instances -> 5^10 possibilities or 97,65,625 instances.
     */
    std::function<void(unsigned int, unordered_map<unsigned int,unsigned int>&,  vector< unordered_map<unsigned int,unsigned int>>&,const unsigned int&,unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> >& )>
            instancesEnumerationBackTrack =[&instancesEnumerationBackTrack]  (unsigned int stgid, unordered_map<unsigned int,unsigned int>& curMapping, vector< unordered_map<unsigned int,unsigned int>>& allMappings,
                    const unsigned int& szStages, unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> >& stg2InstCombinations)->void{
                    if(stgid == szStages) {
                        allMappings.push_back(curMapping);
                        return;
                    }
                    for(const vector<pair<unsigned int,unsigned int>>& instComb: stg2InstCombinations[stgid]){
                        for(const auto& [fnType, fnInstId]: instComb){
                            curMapping[fnType]=fnInstId;
                        }
                        instancesEnumerationBackTrack(stgid+1,curMapping, allMappings, szStages, stg2InstCombinations);
                    }
    };


    cSFC->bst_parlen_idx = cSFC->allPartParSFC.size(); // from partParSFC what is the best partial mapping.
    cSFC->bst_seqlen_time = std::numeric_limits<float>::max();
    cSFC->bst_parlen_time = std::numeric_limits<float>::max();
    float bst_partpar_time = std::numeric_limits<float>::max();
    //    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}} ;
    /*! {{1},{6,4},{5}}; Each Partial SFC is without src and dest block/stage. */
    for(unsigned int ppsidx=0; ppsidx<cSFC->allPartParSFC.size(); ppsidx++){
        const vector<vector<unsigned int>>& partParSFC= cSFC->allPartParSFC.at(ppsidx); ///< for each of the partial parallel SFC of the givenParVNF Blocks
        const unsigned int szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

        /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }} */
        unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
        for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
            const auto& curStg = partParSFC[stgid];
            if(curStg.size() == 1){ // only one function idx=0 in the block, then insert all its instance as combination
                for(unsigned int instid_f0=1; instid_f0<=VNFNetwork->VNFNodes[curStg.at(0)]->numInstances; instid_f0++)
                    stg2InstCombinations[stgid].push_back({{curStg.at(0), instid_f0}});
            }else if(curStg.size() == 2){ // two function in the block, then insert all its instance as combination
                for(unsigned int instid_f0=1; instid_f0<=VNFNetwork->VNFNodes[curStg.at(0)]->numInstances; instid_f0++){
                    for(unsigned int instid_f1=1; instid_f1<=VNFNetwork->VNFNodes[curStg.at(1)]->numInstances; instid_f1++){
                        stg2InstCombinations[stgid].push_back({{curStg.at(0), instid_f0},{curStg.at(1), instid_f1}});
                    }
                }
            } else { // if 3 or more func are parallel
                vector<pair<unsigned int,unsigned int>> curInstComb;
                find_stg2IC_ofGivenBlk(0, curInstComb, stgid, curStg, stg2InstCombinations);
            }
        }//stgid<szStages finding instances possibilities of each stage.
        // show stages wise instances combination

        ///finding all the mapping possibilites for the current partParSFC instance combination at each stage.
        vector< unordered_map<unsigned int,unsigned int>> allMappings;
        unordered_map<unsigned int,unsigned int> curMapping;
        instancesEnumerationBackTrack(0, curMapping, allMappings, szStages, stg2InstCombinations);

        if(showInConsole and showInConsoleDetailed){ // showing partParSFC info
            cout<<"\n partParSFC["<<ppsidx<<"]: "; for(const auto& blks: partParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<") ---------- - --------- - ------";
            for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
                cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
                for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                    cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
                }  cout<<" }";
            }

        }

        unsigned int bestMappingId = allMappings.size();
        float minCurTime = std::numeric_limits<float>::max();
        for(unsigned int mapid=0; mapid<allMappings.size(); mapid++){
            const auto& L_VNFType2Inst= allMappings[mapid]; // Local Mapping
            float parallelCost=0, paketCost=0;
                parallelCost = calcObjectiveValuePar<type_wgt, type_res>(partParSFC,L_VNFType2Inst,  cSFC->index, allSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);
            if(szStages < cSFC->numVNF) // packet delay only in case of parallelism when number of blocks/stages < numOfVNFs (Sequential chain length)
                paketCost = calcTime_PacketsDelay<type_res>(partParSFC, L_VNFType2Inst,  cSFC->index, VNFNetwork, VirtualNetwork);

            float curMapTime = parallelCost+paketCost; // float_max+flaot+max => inf

            if( curMapTime < minCurTime - std::numeric_limits<float>::epsilon()){ // current mapping ka time is less than calc time
                bestMappingId = mapid;
                minCurTime =  curMapTime;
            }

            if(showInConsole and showInConsoleDetailed){
                if(mapid%3==0)cout<<"\n";
                cout<<"\t"<<mapid <<"["; for(const auto &blk: partParSFC){ for(const auto& fnid: blk){ cout<<fnid<<char(96+L_VNFType2Inst.at(fnid))<<" ";  }  } cout<<"]";
                cout<<"["<<curMapTime<<"sec ("<<parallelCost<<"|"<<paketCost<<")]";
            }
        } // for each mapping

        if(ppsidx == 0){// index 0 partial chain is same as given sfc
            if(minCurTime < cSFC->bst_seqlen_time){ // current part chain ka time seq time se kam hai.
                cSFC->bst_seqlen_time = minCurTime;
                cSFC->bst_seqlen_mapping = allMappings[bestMappingId];
            }
        }else{ // if it is parallel where chain length is less than sequential
            if(minCurTime < cSFC->bst_parlen_time){ // current part chain ka time seq time se kam hai.
                cSFC->bst_parlen_time = minCurTime;
                cSFC->bst_parlen_idx = ppsidx;
                cSFC->bst_parlen_mapping = allMappings[bestMappingId];
            }
        }

    }// for each PartParSFC.

    if(cSFC->bst_seqlen_time == std::numeric_limits<float>::max()){
        string errorMsg = "Algorithm failed to find best mapping for SFC["+to_string(cSFC->index)+ "]  Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    if(cSFC->bst_parlen_idx == cSFC->allPartParSFC.size() ){
        string errorMsg = "Algorithm failed to find best partial parallel mapping for SFC["+to_string(cSFC->index)+ "]. Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    if(showInConsole){
        cout<<"\n\n Objective Function Answer For SFC["<<cSFC->index<<"]";
        cout<<"\n\t Sequential: partIdx[0]  time:["<<cSFC->bst_seqlen_time<<"] :(";
            for(const auto &blk: cSFC->allPartParSFC[0]) {
                cout<<"[";  for(const auto& fnid: blk){
                            cout<<fnid<<char(96+cSFC->bst_seqlen_mapping.at(fnid))<<" ";
                        }   cout<<"]";
            }
        cout<<"\n\t Parallel: partIdx["<<cSFC->bst_parlen_idx<<"]  time:["<<cSFC->bst_parlen_time<<"] :(";
        for(const auto &blk: cSFC->allPartParSFC[cSFC->bst_parlen_idx]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout<<fnid<<char(96+cSFC->bst_parlen_mapping.at(fnid))<<" ";
            }   cout<<"]";
        }

    }

}//layerGraphConstruction_and_InstanceSelectionAndRouting1stVersion



int main()
{
// TODO: Finding number of VNF instances, and then VNF Deployement
// TODO: existing heuritics compare with bruteforce
// TODO: try to remove sAdj /pAdj as they are mainly not required.
// TODO: can you parallelise some aspects of program??
// TODO:
// TODO: parallel graph data collection
// TODO: parallelism -> fullSFC se sare partial nikalte time, sabka parallely nikal skte. ek partSFC ke stagewise inst combination parallel nikal skte.
// TODO: X ka mapping.

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

    vector<ServiceFunctionChain*> SFCs;
    try{  readGenericServiceFunctionsChains(testDirName, SFCs);
    } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }
    for(int ni=1; ni<SFCs.size(); ni++) convertSeqSFC_to_FullParVNFBlocks<type_res_local>(SFCs[ni], VNFNetwork, (ni==SFCs.size()-1));
//    SFCs[1]->convertToParallelSFC({{1},{6,4},{5}});
//    SFCs[2]->convertToParallelSFC({{10},{4,8,3},{2,1},{5}});
//    SFCs[3]->convertToParallelSFC({{4,8},{7,10,9}});
//    for(int ni=1; ni<=total_SFC; ni++) SFCs[ni]->showSFC_BlockWise(SFCpar, SFCs[ni]->I_VNFType2Inst);
//    for(int ni=1; ni<SFCs.size(); ni++) SFCs[ni]->showSFC_BlockWise(SFCseq, SFCs[ni]->I_VNFType2Inst);

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
    vector<vector<pair<int,int>>> VNFType_TO_InstID(SFCs.size());
    VNFType_TO_InstID[1] = { {1, 1}, {6, 2}, {4, 1}, {5, 1} };
    VNFType_TO_InstID[2] = { {10, 2}, {4, 2}, {8, 1}, {3, 1}, {2, 1}, {1, 1}, {5, 2} };
    VNFType_TO_InstID[3] = { {4, 3}, {8, 2}, {7, 1}, {10, 1}, {9, 1} };
    assign_ForSFC_VNFType_2_InstID<type_res_local>(VNFType_TO_InstID[1], SFCs[1], VNFNetwork);
    assign_ForSFC_VNFType_2_InstID<type_res_local>(VNFType_TO_InstID[2], SFCs[2], VNFNetwork);
    assign_ForSFC_VNFType_2_InstID<type_res_local>(VNFType_TO_InstID[3], SFCs[3], VNFNetwork);


    int indexxxxxx;
    cin>>indexxxxxx;
    for(const auto& fnType: SFCs[indexxxxxx]->vnfSeq){
        if(fnType ==SFCsrc or fnType ==SFCdst)continue;
        VNFNetwork->utilization[fnType][SFCs[indexxxxxx]->I_VNFType2Inst[fnType]] -= SFCs[indexxxxxx]->trafficArrivalRate;
    }
    cin>>VNFNetwork->VNFNodes[4]->serviceRate;
    cout<<endl;
    for(int vnfi=1; vnfi<=10; vnfi++){
        for(int inst=1; inst<=VNFNetwork->VNFNodes[vnfi]->numInstances; inst++){
            cout<<"F["<<vnfi<<"]["<<inst<<"]:"<<VNFNetwork->utilization[vnfi][inst];
            cout<<" ("<<(VNFNetwork->utilization[vnfi][inst]/VNFNetwork->VNFNodes[vnfi]->serviceRate)*100<<")";
            cout<<" | ";
        }
        cout<<endl;
    }
    auto ft_start = std::chrono::steady_clock::now();
    parVNFBlocks_ClusterAssignment_ForSFC(SFCs[indexxxxxx]);
    try{
        layerGraphConstruction_and_InstanceSelectionAndRouting(SFCs[indexxxxxx], SFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork, true, true);
    }catch (std::exception const &e) {  std::cerr << "\ncaught: " << e.what() ; }
    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";

//    PhysicalNetwork->showPNs_Description();
//    VirtualNetwork->showVMs_Description();
//    VNFNetwork->showVNFs_Description();
//    for(int ni=1; ni<=total_SFC; ni++) SFCs[ni]->showSFC_BlockWise(SFCseq);
//    for(int ni=1; ni<=total_SFC; ni++) SFCs[ni]->showSFC_BlockWise(SFCpar);


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
//        cout<<"\nSFCs["<<ni<<"] E2E: ";
//        cout<<"Seq["<<calcObjectiveValueSeq<type_wgt_local, type_res_local>(SFCs[ni], SFCs[ni]->I_VNFType2Inst, SFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork)<<"]  ";
//        cout<<"Par["<<calcObjectiveValuePar<type_wgt_local, type_res_local>(SFCs[ni]->vnfBlocksPar,SFCs[ni]->I_VNFType2Inst,  SFCs[ni]->index, SFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork)<<"]  ";
//        cout<<"Pkt["<<calcTime_PacketsDelay<type_res_local>(SFCs[ni]->vnfBlocksPar, SFCs[ni]->I_VNFType2Inst,  SFCs[ni]->index, VNFNetwork, VirtualNetwork)<<"]";
//    }
//    debug=1;





//////////////////////////////////////////////calling objective function/////////////////////////////////


//    float objectiveValue = calcObjectiveValuePar<type_wgt_local, type_res_local>(SFCs[3], SFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork, true);
//    cout << objectiveValue << endl;


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