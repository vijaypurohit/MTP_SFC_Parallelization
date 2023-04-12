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
using type_delay = float;
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
//#define maxSFCLength 10
//#define maxVNF_Instances 5
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
const unsigned int factor_bandwidth = 1000000; unsigned int bandwidthNW = 10;//10
type_delay velocityFactor = 1.0;
unsigned int  speedOfLight = 300000000;//300000000
//unsigned int total_SFC = 0;
type_delay read_write_time_per_bit = 0.077e-3;

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
                                                            const VirtualNetworkFunctions<type_res> *VNFNetwork, const VirtualMachines<type_res> *VirtualNetwork,
                                                            const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {

    unsigned int mxPathsK = 3;
    std::function<unsigned int(const unsigned int&)>numPathsToTake_Decision = [&mxPathsK](const unsigned int& numOfPathPairs) -> unsigned int{
        if(numOfPathPairs <= 2) return numOfPathPairs; ///< if less than two then take both
        else if(numOfPathPairs <= 4) return 2; // if 2-4 paths then take 2
        return mxPathsK; // if number of paths from src to next stage is more than 2
    };

    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}} ;
    const unsigned int szStages = partParSFC.size();
    struct somePreComputedDelayParameters{
        type_delay exeDelay{};
        type_delay prcDelay{};
        unordered_map<unsigned int, type_delay> queuingDelay;
    };
    unordered_map<unsigned int, somePreComputedDelayParameters> vnfDelays;
    for(const auto&blk : partParSFC){
        for(const auto& fn: blk){
            VNFNode<type_res> *dstVNFNode = VNFNetwork->VNFNodes.at(fn);
            vnfDelays[fn].prcDelay = 10000 *calcTime_MeanProcessingDelayVNF<type_res>(dstVNFNode);
            vnfDelays[fn].exeDelay = 100 *calcTime_FunctionExecutionDelay<type_res>(dstVNFNode);
            for(int fnInst=1; fnInst<=dstVNFNode->numInstances; fnInst++) {
                vnfDelays[fn].queuingDelay[fnInst] = 10000 *calcTime_QueuingDelay<type_res>(dstVNFNode, fnInst,VNFNetwork->utilization, allSFC[cSFC->index]);
            }
        }
    }

    /*! priority queue node to find minimum dist and minimum utilization path from source to destination/current stage.*/
    struct pqNode{
        unsigned int x, y; ///< source and destination pair of previous and current stage lgNode
        type_delay mindist; ///< minimum delay of the path from source to current processing node node
        type_delay utilization; ///< max utilization of the path from source to current processing node node
        vector<unsigned int> path; ///< path constructed till now from source
        pqNode()=default;
        pqNode(unsigned int givenlgSrcId, unsigned int givenlgDstId, type_delay givenDist, type_delay givenUtilization, std::vector<unsigned int> givenPath):x(givenlgSrcId),y(givenlgDstId),mindist(givenDist),utilization(givenUtilization){
                path = std::move(givenPath);
        }
        bool operator<(const struct pqNode& other) const { // overloaded operator for priority queue
            if(mindist == other.mindist){ return utilization > other.utilization; //min heap, return pair of x-y with min utilization
            } else return mindist > other.mindist; //min heap, return pair of x-y with minimum distance.
        }
    };

    /*! Layer Graph Node Vertex */
    struct lgNode{
        unsigned int idx{}; ///< index to detect node uniquely 
        vector<pair<unsigned int,unsigned int>> instCombination; ///< {fnType, instId} pairs showing instance combination at this node
        unordered_map<unsigned int, unsigned int> cntPN; ///< count of physical node,frequency in the instance combination
        type_delay utilization{0}; ///< utilisation percentage of all the instances present in the lgNode
        vector<pair<unsigned int, type_delay>> children; ///<  next stage lgNode index and its distance, that is pair of this->node = {next stg node, min dist}.
        vector<pqNode> kpaths; ///< number of shortest path traverse through this lgNode
        
        lgNode()=default;
        explicit lgNode(unsigned int index):idx(index){};
        lgNode(unsigned int index, const vector<pair<unsigned int,unsigned int>>& givenIC):idx(index), instCombination(givenIC){ } 
    };

    
    unordered_map<unsigned int, lgNode> idx2lgNode;
    idx2lgNode[0] = lgNode(0);
    idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{1, 1}});            idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{1, 2}});        idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{1, 3}});
    idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{6, 1}, {4, 1}});    idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{6,1}, {4, 2}}); idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{6, 2}, {4, 1}});    idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{6,2}, {4, 2}});
    idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{5, 1}});            idx2lgNode[idx2lgNode.size()] = lgNode(idx2lgNode.size(), {{5, 2}});
    unsigned int lgDSTid = idx2lgNode.size();
    idx2lgNode[lgDSTid] = lgNode(lgDSTid);

    unordered_map<unsigned int, vector<unsigned int>> stg2lgNode;
    stg2lgNode = {
            { 0, { 1, 2, 3 } },
            { 1, {  4,  5,  6,  7 } },
            { 2, { 8, 9  } }
    };
    vector<unordered_set<unsigned int>> uniqLgidInStgOfKPaths(szStages); ///< lgNode ids in each stage which is used in k shortest path.
    
    auto processs_min_heap = [&idx2lgNode, &uniqLgidInStgOfKPaths, &numPathsToTake_Decision, &showInConsoleDetailed](unsigned int id_curstg, priority_queue<pqNode>&pq, bool toInsert = true){
        if(showInConsoleDetailed){ cout<<"\nTotalPathPairs:"<<pq.size();}
        unsigned int numOfPathsToConsider= numPathsToTake_Decision(pq.size());
        while(numOfPathsToConsider>0 and !pq.empty()){
            pqNode min_path = pq.top(); /*! pair of src and dst lg node which produce min distance */
            min_path.path.push_back(min_path.y); // create path before inserting into node.
            if(toInsert)uniqLgidInStgOfKPaths[id_curstg].insert(min_path.y);
            idx2lgNode[min_path.y].kpaths.push_back(min_path); ///<consider this node in path
            idx2lgNode[min_path.x].children.push_back({min_path.y, min_path.mindist}); ///< source to child mapping to travers the path
            pq.pop(); numOfPathsToConsider--;
            if(showInConsoleDetailed){ cout<<"\n    p:"<<min_path.mindist<<"sec | "<<min_path.utilization<<"% [";
                for(const auto& kkk: min_path.path)  cout<<kkk<<" -> ";
            }
        }
    };
    
    type_delay  T_tx_init = 100*calcTime_TransmissionDelay();
/********************************************************************************************************/
    /*! Special Case: From Dummy SRC to First Stage(index 0).
     * First Stage consist of Layer Graph Node Indexes (lgnIdy). Each Layer Graph Node consist of instance combinations pairs{vnf type, its instance id}.
     * Calculate maximum delay taken to process that Layer Graph Node, as completion time would be when all instaces in that node finish their execution.
     * inter-duplication + transmission time +  max( "intra-duplication" + "time taken in each server" + "intra-merging")
     * Push the pairs into min priority queue to find pairs which produce minimum delay.
     */ 
    unsigned int id_curstg=0; 
    priority_queue<pqNode> pq; ///< to find minimum delay path in the x-y pairs of current and previous stage.
    
    for(const unsigned int &lgnIdy: stg2lgNode[id_curstg]){
            lgNode& lgy = idx2lgNode[lgnIdy]; ///< layer graph node y, x is dummy source
            if (showInConsoleDetailed) { cout << "\nlg:" << lgy.idx ; cout<<" ["; for(const auto& givenPair: lgy.instCombination){ cout<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  }  cout<<"]";}
            
            type_delay utilization_sum_lgy=0, servicerate_sum_lgy=0; ///< for the lgNode utilization of all instances
            unordered_map<unsigned int, type_delay> T_exe_server; ///< Execution time of server. Maximum delay among all the parallel instances in a server.
        /*! Processing of current lgNode: count of physical servers, max time in each server */
        for (const auto &[d_fnType, d_fnInst]: lgy.instCombination) {
                const auto &d_vm_id = VNFNetwork->I_VNFinst2VM.at(d_fnType).at(d_fnInst); const auto &d_pn_id = VirtualNetwork->I_VM2PN.at(d_vm_id);
                
            lgy.cntPN[d_pn_id] += 1; //< freq of PN in current lgn
            T_exe_server[d_pn_id] = max(T_exe_server[d_pn_id], vnfDelays[d_fnType].prcDelay + vnfDelays[d_fnType].exeDelay + vnfDelays[d_fnType].queuingDelay[d_fnInst]);
            
            if (VNFNetwork->utilization.count(d_fnType) and VNFNetwork->utilization.at(d_fnType).count(d_fnInst)){ ///< old utilization till now of VNF
                utilization_sum_lgy += VNFNetwork->utilization.at(d_fnType).at(d_fnInst);
            } servicerate_sum_lgy += VNFNetwork->VNFNodes.at(d_fnType)->serviceRate;
            
            if (showInConsoleDetailed) { cout << "\n    :F[" << d_fnType << char(96 + d_fnInst) << "]VM[" << d_vm_id << "]PN[" << d_pn_id << "]"
                  <<"  qd[" << vnfDelays[d_fnType].queuingDelay[d_fnInst] << "] | prc[" << vnfDelays[d_fnType].prcDelay << "] | exe[" << vnfDelays[d_fnType].exeDelay << "]";}
        }
        lgy.utilization = (utilization_sum_lgy/servicerate_sum_lgy)*100; //update utilization of the lgy Node

        type_delay mx_delay_x_y = 0;///< maximum delay of the current lgNode pair (x=dummySrc,y=current lgNode).

        /*! Calculation of packet processing time. inter duplication from src to different servers, intra duplication (within same server multiple nodes) and intra merging*/
        for (const auto &[pn_y, pn_y_parallelcnt]: lgy.cntPN) { /*! for each physical server in previous stage*/
            type_delay T_d_hdr = 100*calcTime_IntraDuplicationTime(pn_y_parallelcnt);
            type_delay T_m_hdr = 100*calcTime_IntraMergingTime(pn_y_parallelcnt);
            mx_delay_x_y = max(mx_delay_x_y, T_d_hdr+T_m_hdr + T_exe_server[pn_y]);
        }//curStgPN
        type_delay T_d_pkt = 100*calcTime_InterDuplicationTime( lgy.cntPN.size()); ///< inter duplication time from source to lgy
        mx_delay_x_y +=  T_d_pkt + T_tx_init;
        if (showInConsoleDetailed) {cout << "\n     max:" << mx_delay_x_y; cout<<" | d_pkt:"<<T_d_pkt;}
        pq.emplace(idx2lgNode[0].idx, lgy.idx, mx_delay_x_y, lgy.utilization, vector<unsigned int>{idx2lgNode[0].idx}); /// constructing the dummy src to lgy path.
    }

    /*! Process paths from soruce to next stage. */
    processs_min_heap(id_curstg, pq);

/********************************************************************************************************/   
    /*! From each lgNode (inst combination) in current stage to lgNode(instance combination) in prevous stage. Repeat till last stage.
     *  Process the lgy node first: count physical servers, and maximum time of execution of server.
     *  Then for each lgx node in prev stage. Calculate inter-duplication + transmission + processing time. Take maximum of it for current lgy node.
     *  Out of all servers in lgy take maximum time to be delay for lgy.
     */
    for( id_curstg=1; id_curstg<szStages; id_curstg++) {/*!< Iterating stage ID from 0 to last index */
        unsigned int id_prvstg = id_curstg-1; ///< previous stage id to process
        pq = priority_queue<pqNode>();

    /*****lgnIdy************************************/
        for (const unsigned int &lgnIdy: stg2lgNode[id_curstg]) { /*!< Iterating all the layer graph node index in current stage(destination)  */
                lgNode& lgy = idx2lgNode[lgnIdy]; ///< current layer Layer Graph Node
                        if (showInConsoleDetailed) { cout << "\nlgy:" << lgy.idx ; cout<<" ["; for(const auto& givenPair: lgy.instCombination){ cout<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  }  cout<<"]";}

                /*! Processing of current lgNode: count of physical servers, max time in each server */
                    type_delay utilization_sum_lgy=0, servicerate_sum_lgy=0;
                    unordered_map<unsigned int, type_delay> T_exe_server; ///< Execution time of server. Maximum delay among all the parallel instances in a server.
                for (const auto &[d_fnType, d_fnInst]: lgy.instCombination) { /*!< Iterating all instances combinations in next stage(destination)  */
                        const auto &d_vm_id = VNFNetwork->I_VNFinst2VM.at(d_fnType).at(d_fnInst); const auto &d_pn_id = VirtualNetwork->I_VM2PN.at(d_vm_id);

                        if (VNFNetwork->utilization.count(d_fnType) and VNFNetwork->utilization.at(d_fnType).count(d_fnInst)){ ///< old utilization till now of VNF
                            utilization_sum_lgy += VNFNetwork->utilization.at(d_fnType).at(d_fnInst);
                        } servicerate_sum_lgy += VNFNetwork->VNFNodes.at(d_fnType)->serviceRate;

                        lgy.cntPN[d_pn_id] += 1;
                        type_delay totalDelay_fnType =  vnfDelays[d_fnType].prcDelay + vnfDelays[d_fnType].exeDelay + vnfDelays[d_fnType].queuingDelay[d_fnInst]; ///< total delay of one instance of src and one instance of dst
                        T_exe_server[d_pn_id] = max(T_exe_server[d_pn_id], totalDelay_fnType);

                        if (showInConsoleDetailed) { cout << "\n    :F[" << d_fnType << char(96 + d_fnInst) << "]VM[" << d_vm_id << "]PN[" << d_pn_id << "]"
                                        <<"  qd[" << vnfDelays[d_fnType].queuingDelay[d_fnInst] << "] | prc[" << vnfDelays[d_fnType].prcDelay << "] | exe[" << vnfDelays[d_fnType].exeDelay << "] | sum:" << totalDelay_fnType ;}
                }//d_fnType, d_fnInst
                lgy.utilization = (utilization_sum_lgy/servicerate_sum_lgy)*100;
        /*****lgnIdx************************************/
            for (const unsigned int &lgnIdx: uniqLgidInStgOfKPaths[id_prvstg]) {/*!< Iterating layer graph nodes index (which are in Path) in previous stage(source)  */
                const lgNode &lgx = idx2lgNode[lgnIdx];///< current layer Layer Graph Node

                    if (showInConsoleDetailed) { cout << "\n  lgx:" << lgx.idx; cout << " ["; for (const auto &givenPair: lgx.instCombination) {  cout << givenPair.first << char(givenPair.second - 1 + 'a') << " "; } cout << "]";}

            /*****************************************/
                /*! Once x = lgx, y = lgy are fixed. We will find maximum time for each physical server in y to determine edge (x,y). */
                type_delay mx_delay_x_y = 0;///< maximum delay of the current instance combination of the pair (x,y).
                for(const auto &[pn_y, pn_y_parallelcnt]: lgy.cntPN){ /*! for each physical server in cur node*/
                        /*! Calculation of packet processing time. inter mergring ( different server from cur server), intra duplication (within same server multiple nodes), intra Merging (within same server multiple nodes)*/

                        unsigned int py_px_same=0; if(lgx.cntPN.find(pn_y) != lgx.cntPN.end()) py_px_same =  1;
                        unsigned int cntPrevHopDiffServer = lgx.cntPN.size() - py_px_same;
                        type_delay T_m_pkt = 100*calcTime_InterMergingTime(cntPrevHopDiffServer);

                        type_delay T_d_hdr = 100*calcTime_IntraDuplicationTime(pn_y_parallelcnt);
                        type_delay T_m_hdr = 100*calcTime_IntraMergingTime(pn_y_parallelcnt);

                        type_delay mx_pktPrc =  T_m_pkt  + T_d_hdr + T_m_hdr;///< overall total time spent in packet processing from src to dest.

                        if(showInConsoleDetailed) {
                            cout<<"\n     :"<<pn_y_parallelcnt<<"[py:"<<pn_y<<"]"<<"  prvD:"<<cntPrevHopDiffServer<<" s("<<py_px_same<<")"
                                <<"   [m_pkt:"<<T_m_pkt<<" d_hdr:"<<T_d_hdr<<" m_hdr:"<<T_m_hdr<<"]"<<"   mxServer:"<<T_exe_server[pn_y];
                        }
                    /*! Calculation of inter-duplication time transmission time and  propagation time (we can duplicate the packets right before sending them. */
                    type_delay mx_interdupTxPx_for_y = 0;
                    for (const auto &[pn_x, pn_x_parallelcnt]: lgx.cntPN) { /*! for each physical server in previous node*/
                            unsigned int px_py_same = 0;
                            if (lgy.cntPN.find(pn_x) != lgy.cntPN.end()) px_py_same = 1;
                            unsigned int cntNextHopDiffServer = lgy.cntPN.size() - px_py_same;
                            type_delay T_d_pkt = 100 * calcTime_InterDuplicationTime(cntNextHopDiffServer);
                            type_delay T_tx=0, T_px=0;
                            if(px_py_same == 0){ /// if both server are different then there is transmission and propagation delay
                                T_tx = 100*T_tx_init; T_px = 10000000 * calcTime_PropagationDelay<type_wgt, type_res>(pn_x, pn_y,PhysicalNetwork);
                            }
                            mx_interdupTxPx_for_y = max(mx_interdupTxPx_for_y, T_d_pkt + T_tx + T_px);///< overall total time spent in sending packet from src to dest.
                            if (showInConsoleDetailed) {
                                cout << "\n          " << pn_x_parallelcnt << "(px:" << pn_x << ")  " << "nxtD:" << cntNextHopDiffServer << " s(" << px_py_same << ")"
                                     << "   [d_pkt:" << T_d_pkt << " tx:" << T_tx << " px:" << T_px << "]";
                            }
                    }//prvStgPN
                    mx_delay_x_y = max(mx_delay_x_y, mx_interdupTxPx_for_y + mx_pktPrc + T_exe_server[pn_y]);

                    if(showInConsoleDetailed) {
                        cout<<"   mxTxPx:"<<mx_interdupTxPx_for_y;
                        cout<<"\n       px-py:"<< mx_interdupTxPx_for_y + mx_pktPrc + T_exe_server[pn_y];
                    }

                }//curStgPN
//                if(showInConsoleDetailed) {cout << "\n     max:" <<mx_delay_x_y << "";}

                for(const pqNode& kpq: lgx.kpaths) {/// for each min path in lgx node, push this pair also.
                        pq.emplace(lgnIdx, lgnIdy, kpq.mindist + mx_delay_x_y, max(kpq.utilization, lgy.utilization), kpq.path);
                } //kpq
            }//lgnIdx
        }//lgnIdy

         processs_min_heap(id_curstg, pq);
    }//id_curstg

/********************************************************************************************************/
    /*! From Last Stage to Dummy DST */
    id_curstg=szStages-1;
    pq = priority_queue<pqNode>();
    for(const unsigned int &lgnIdx: uniqLgidInStgOfKPaths[id_curstg]){
        lgNode& lgx = idx2lgNode[lgnIdx];
        if (showInConsoleDetailed) { cout << "\nlg:" << lgx.idx ; cout<<" ["; for(const auto& givenPair: lgx.instCombination){ cout<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  }  cout<<"]";}

        type_delay T_m_pkt = 100*calcTime_InterMergingTime(lgx.cntPN.size());
        type_delay mx_delay_x_y = 100*T_tx_init + T_m_pkt;

        for(const pqNode& kpq: lgx.kpaths) {
            pq.emplace(lgnIdx, lgDSTid, kpq.mindist + mx_delay_x_y, max(kpq.utilization, lgx.utilization), kpq.path);
        }
    }

    processs_min_heap(id_curstg, pq, false);

    for(const auto lgid: idx2lgNode[lgDSTid].kpaths[0].path){
        cout<<" | IC[";
            for(const auto &givenPair: idx2lgNode[lgid].instCombination){
                    cout<<givenPair.first<<char(givenPair.second-1+'a')<<" ";
            }
        cout<<"]";
    }
//    for(const auto &fir: uniqLgidInStgOfKPaths){ cout<<" | v{ "; for(const auto& val: fir){  cout<<val<<" "; }  cout<<" }";}
    if(showInConsole){
        for(const auto& fir: idx2lgNode){ const lgNode& sec= fir.second;
            cout<<"\n "<<fir.first<<"| ";
            cout<<" idx:"<<sec.idx <<" | use:"<<sec.utilization<<" | kpaths:"<<sec.kpaths.size() ;
            cout<<" | child{ "; for(const auto& val: sec.children){ cout<<"("<<val.first<<": "<<val.second<<") "; }  cout<<" }";
            cout<<" | IC["; for(const auto& givenPair: sec.instCombination){ cout<<givenPair.first<<char(givenPair.second-1+'a')<<" ";}  cout<<"]";
        }
    }///show lgNodes
}


int main()
{
// TODO: Finding number of VNF instances, and then VNF Deployement
// TODO: parallelism -> exploring more partial sfc possibilities.
// TODO: existing heuritics compare with bruteforce
// TODO: can you parallelise some aspects of program??

//    auto ft_start = std::chrono::steady_clock::now();
//    clusterSizeEnumeration(int(clusterSz.size()),maparVNFs_Cluster_AssignemtnxSFClen);
//    all_nCk(int(nCk.size()), maxSFClen);
//    vector<vector<int>> parSFC_Full = {{1}, {2,3,4,5}};
//    parVNFs_Cluster_Assignment(5, parSFC_Full, false);
//    if(debug)cout<<"\nTime:"<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";

    using type_wgt_local = type_delay; ///< Determine Edge Weights Data TYPE (unsigned int or FLOAT)
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
//    for(const auto& fnType: SFCs[indexxxxxx]->vnfSeq){
//        if(fnType ==SFCsrc or fnType ==SFCdst)continue;
//        VNFNetwork->utilization[fnType][SFCs[indexxxxxx]->I_VNFType2Inst[fnType]] -= SFCs[indexxxxxx]->trafficArrivalRate;
//    }
    VNFNetwork->utilization[1][2] = 2.5;
    VNFNetwork->utilization[1][3] = 1.5;
//    cin>>VNFNetwork->VNFNodes[4]->serviceRate;
//    cout<<endl;
//    for(int vnfi=1; vnfi<=10; vnfi++){
//        for(int inst=1; inst<=VNFNetwork->VNFNodes[vnfi]->numInstances; inst++){
//            cout<<"F["<<vnfi<<"]["<<inst<<"]:"<<VNFNetwork->utilization[vnfi][inst];
//            cout<<" ("<<(VNFNetwork->utilization[vnfi][inst]/VNFNetwork->VNFNodes[vnfi]->serviceRate)*100<<")";
//            cout<<" | ";
//        }
//        cout<<endl;
//    }
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


//    type_delay objectiveValue = calcObjectiveValuePar<type_wgt_local, type_res_local>(SFCs[3], SFCs, VNFNetwork, VirtualNetwork, PhysicalNetwork, true);
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