//
// Created by vijay on 24-03-2023.
//

#ifndef SFC_PARALLELIZATION_TIMECALCULATIONFUNCTIONS_H
#define SFC_PARALLELIZATION_TIMECALCULATIONFUNCTIONS_H
/*! Calculate Transmission delay between two Physical Nodes n1 and n2
 * (amount of time required to push all the packet's bits into the wire) \n
 * Right now bandwidth in all links is assumed to be same and defined in global variable. \n
 * Formula = (packetSize of the network including header size) / (bandwidth between two nodes n1 and n2). \n
 * @refer bandwidthNW, packetBodySize, packetHeaderSize
 * @param n1 Physical Node n1 @param n2 Physical Node n2. Use in case Bandwidth is different in each link.
 * @param SFC object of class type ServiceFunctionChain. Use in case packet size is different for each SFC.
 * @return transmission time in float. in seconds.\n
 *  random value = 0.0008192 seconds ((1024 * 8 b) / (10 * 10^6 b/s)) \n
 *  min = 0.0000000001 s ( 1 b/ 10000 * 10^6 ). \n
 *  max = 0.1024 s (10240/ 10^6).
 *  @future packet size can be variable according to each SFC.
 */
float calcTime_TransmissionDelay([[maybe_unused]] unsigned int n1=0, [[maybe_unused]] unsigned int n2=0,
                                 [[maybe_unused]] ServiceFunctionChain *SFC = nullptr){
    auto bandwidth_n1_n2 = (float)bandwidthNW * factor_bandwidth; // in Mb/s
    auto packetsize_n1_n2 = ((float)packetBodySize + (float)packetHeaderSize)*factor_packet; // in bits
    return packetsize_n1_n2/(bandwidth_n1_n2);
}

/*! Calculate the propagation delay between two Physical Nodes n1 and n2
 * (time duration taken for a signal to reach its destination)\n
 * Formula = (Distance between two nodes n1 and n2) / Velocity of light in the medium \n
 *        = Dist[n1][n2] / (transmission medium * speed of light in vaccum) \n
 * @refer velocityFactor, speedOfLight and dist of Physical Graph
 * @param n1 Physical Node n1 @param n2 Physical Node n2
 * @param PhysicalNetwork object of class type PhysicalGraph.
 * @return propagation time in float. in seconds. \n
 *  random value =  0.00001 second (3000 m / (1.0 * 3 * 10^8)). \n
 *  min = 3.33e-9 s ( 1 m/ 3*10^8 ). \n
 *  max = 0.04252 s (12756000 m/ 3*10^8).
 */
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
float calcTime_PropagationDelay(unsigned int n1, unsigned int n2, const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork){
    auto len_n1_n2 = (float)PhysicalNetwork->dist[n1][n2];
    auto speed = velocityFactor * (float)speedOfLight;
    return (len_n1_n2 / (speed));
}

/*! Calculate the Queuing delay of a particular VNF (time a job waits in a queue until it can be executed) \n
 * Formula = lambda_c / (mu_f * (mu_f - lambda_c)) \n
 * where lambda_c is traffic arrival rate of chain c or SFC s  and mu_f is the service rate of SFC or say VNF. \n
 * Both in packets per second. \n
 * @param curVNFid VNF Node Id whose arrival rate we have to calculate
 * @param mu_f service rate of current VNF.
 * @param curSFC SFC class whose VNF we are considering
 * @param allSFC all the SFC in the network
 * @return queueing delay in seconds.
 * Sabse phele ham sab SFC me dekhenge ki ye curVNF present hai kya, and curVNF ka jo instance use kr rhe wohi hai to add kr denge us SFC ka arrival rate
 */
float calcTime_QueuingDelay(const int& curVNFid, const float &mu_f, ServiceFunctionChain* curSFC, const vector<ServiceFunctionChain*>& allSFC){
    float lambda = 0;
    for(size_t c=1; c<=total_SFC; c++){
        ServiceFunctionChain *sfc = nullptr;
        if(!(sfc = allSFC[c])){
            string errorMsg = "Error in accessing SFC. Function: ";
            throw runtime_error(errorMsg + __FUNCTION__);
        }
        /// is curVNF of curSFC present in this SFC? and that VNF instance is same as curVNF instance?
        if(sfc->isVNF_Present[curVNFid] and sfc->I_VNFType2Inst[curVNFid] == curSFC->I_VNFType2Inst[curVNFid]){
            lambda += sfc->trafficArrivalRate;
//            cout<<"lambda += SFC["<<sfc->index<<"]";
        }
    }
    if(lambda >= mu_f){
        string errorMsg = "SFC["+ to_string(curSFC->index)+"]->VNF["+ to_string(curVNFid)+"]. Arrival rate ["+ to_string(lambda)+"] >= processing rate of VNF["+ to_string(mu_f)+"]. Function: ";
        throw runtime_error(errorMsg + __FUNCTION__);
    }
    return (lambda / (mu_f * (mu_f - lambda)));
}

/*! Return the time required to complete the execution of the function VNF f. \n
 * @param VNFNode for function execution time.
 * @return execution time in seconds.
 */
template<typename type_res=unsigned int>
float calcTime_FunctionExecutionDelay(VNFNode<type_res> *VNFNode){
    auto exeTime = (float)VNFNode->executionTime;
    return exeTime;
}

/*! Calculate the mean processing delay of particular VNF, time it takes function to process the one packet. \n
 * Fromula = 1/serviceRate of VNF.
 * @param VNFNode VNF Node class for service rate
 * @return mean processing delay in seconds.
 */
template<typename type_res=unsigned int>
float calcTime_MeanProcessingDelayVNF(VNFNode<type_res>* VNFNode){
    auto mu_f = (float)VNFNode->serviceRate;
    return 1/mu_f;
}


/*! The time required by a server v to duplicate the whole packet for instances
 * installed in d different servers of the next step.\n
 * Formula = T_d_pkt * (d-1);
 * @return INTER duplication time in seconds.
 * For example; packetBodySize = 1000 B; packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_pkt=0.630784(d=2),  =1.26157(d=3),  =1.892352(d=4)
 */
float calcTime_InterDuplicationTime(unsigned int cntNextHopDiffServer){
    if( cntNextHopDiffServer <= 1) return 0;
    auto packetsize_n1_n2 = ((float)packetBodySize + (float)packetHeaderSize)*factor_packet;
    float time_to_duplicate_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_duplicate_packet*(float)(cntNextHopDiffServer-1);
}
/*! the time required by a server v to duplicate only the header for its p parallel instances. \n
 * Formula = T_d_hdr * (p-1);
 * @return INTRA duplication time in seconds.
 * For example;packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_hdr=0.014784(p=2),  =0.029568(p=3),  =0.044352(p=4)
 */
float calcTime_IntraDuplicationTime(unsigned int cntParallelServer){
    if(cntParallelServer <= 1) return 0;
    float time_to_duplicate_header = (float)packetHeaderSize *factor_packet * read_write_time_per_bit;
    return time_to_duplicate_header*(float)(cntParallelServer-1);
}

/*! the time required by a server v to merge the whole packets from d different servers of the previous step, \n
 * Formula = T_m_pkt * (d-1);
 * @return INTER merging time in seconds.
 */
float calcTime_InterMergingTime(unsigned int cntPrevHopDiffServer){
    if(cntPrevHopDiffServer <= 1) return 0;
    auto packetsize_n1_n2 = ((float)packetBodySize + (float)packetHeaderSize)*factor_packet;
    float time_to_merging_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_merging_packet*(float)(cntPrevHopDiffServer-1);
}

/*! the time required by a server v to merge the headers from its p parallel instances \n
 * Formula = T_m_hdr * (p-1);
 * @return INTRA merging time in seconds.
 */
float calcTime_IntraMergingTime(unsigned int cntParallelServer){
    if(cntParallelServer <= 1) return 0;
    float time_to_merge_header = (float)packetHeaderSize *factor_packet * read_write_time_per_bit;
    return time_to_merge_header*(float)(cntParallelServer-1);
}

/*! Calculate time it takes for packet duplication, packet merging for each physical SERVER at each level/stage. \n
 * for example: Physical Nodes at each level = { {SFCsrc}, {3, 3, 4}, {3,4,4,5,5,6}, {3,3,4,7,6}, {SFCdst} }; \n
 * SRC -> next level requires only inter duplication AND secLastLevel to LastLevel requires only inter merging time \n
 * at each level requires -> inter merging, intra duplication, intra merging, inter duplication
  * @param SFC given SFC object whose delay we have to colculate
  * @param VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
  * @param VirtualNetwork  VirtualMachines object to find the required Physical Node of the given VM id
  * @param showInConsole to show the calculation debug info in console. Default False.
  * @tparam type_res resource data type. default=unsigned int.
  * @return packet delay of SFC in seconds.
 */
template<typename type_res=unsigned int>
float calcTime_PacketsDelay(ServiceFunctionChain *SFC, VirtualNetworkFunctions<type_res> *VNFNetwork, const VirtualMachines<type_res> *VirtualNetwork, bool showInConsole=false){
//    auto ft_start = std::chrono::steady_clock::now();
    if(debug) cout<<"\n>>>[Function Running: "<<__FUNCTION__<<"]";

//    SFC->vnfBlocksPar = { {SFCsrc}, {3, 3, 4}, {3,4,4,5,5,6}, {3,3,4,7,6}, {SFCdst} };

    int szStages = int(SFC->vnfBlocksPar.size());
    int stg_nxt=1;

    vector<unordered_map<int,int>> stgPN(szStages); ///< all the DIFFERENT physical nodes {nodeid, freq} belonging to each stages.
    unordered_map<int, unordered_map<int, float>> pktdist; ///<stores paket delay for each PHYSICAL NODE
    float totalPktDelay = 0; // total max delay of SFC from src to dst incurred in a process to duplicate and merge packets.

    /// STG[0] --> STG[1] This for-loop is Only For Source Level to Next Level, only duplication required.
    for (const int &vnf_dst_idx: SFC->vnfBlocksPar[stg_nxt]){  //int vnf_dst_inst_idx =0, vm_dst_idx=0,pn_dst_idx =vnf_dst_idx;
              int vnf_dst_inst_idx = SFC->I_VNFType2Inst[vnf_dst_idx]; int vm_dst_idx = VNFNetwork->I_VNFinst2VM[vnf_dst_idx][vnf_dst_inst_idx];int pn_dst_idx = VirtualNetwork->I_VM2PN[vm_dst_idx];
        stgPN[stg_nxt][pn_dst_idx] += 1;
    }
    /// float T_d_pkt = Time to duplicate packets to d diff server from level_i to level_i+1, Here only duplication time required.
    totalPktDelay = pktdist[stg_nxt-1][SFCsrc] = calcTime_InterDuplicationTime(stgPN[stg_nxt].size());
    
    if(showInConsole) {
        cout<<"\n:STG[0]:-->:STG[1]: "<<stgPN[stg_nxt].size()<<"("; for(const auto& pn: stgPN[stg_nxt]) cout<<" PN["<<pn.first<<":"<<pn.second<<"]";  cout<<" )\t[d_pkt:"<<pktdist[stg_nxt-1][SFCsrc]<<"]";
    }

    /*!  From STG = SRC_STG + 1  to szStages-2 Stage, which is prev to DST STG(szStages-1).
     *   each loop first calculates the information of next stage and then process the current stage.
     *   At last, second last stage is info is calculated but not processed.
     */
    for(  stg_nxt=2; stg_nxt <= szStages-1; stg_nxt++) {
        int stg_cur = stg_nxt-1;
        int stg_prv = stg_cur-1;

    // processing the next level stgi+1
        if(stg_nxt != szStages-1){ // not to calulate info for dst vnf stage
        for (const auto &vnf_dst_idx: SFC->vnfBlocksPar[stg_nxt]) { //int vnf_dst_inst_idx =0,vm_dst_idx=0,pn_dst_idx =vnf_dst_idx;
                int vnf_dst_inst_idx = SFC->I_VNFType2Inst[vnf_dst_idx]; int vm_dst_idx = VNFNetwork->I_VNFinst2VM[vnf_dst_idx][vnf_dst_inst_idx];  int pn_dst_idx = VirtualNetwork->I_VM2PN[vm_dst_idx];
            stgPN[stg_nxt][pn_dst_idx] += 1;
        }}

        if(showInConsole) {
            cout<<"\n:STG["<<stg_cur<<"]:-->:STG["<<stg_nxt<<"]: "<<stgPN[stg_nxt].size()<<"("; for(const auto& pn: stgPN[stg_nxt]) cout<<" PN["<<pn.first<<":"<<pn.second<<"]"; cout<<" )";
        }

        unsigned int nxtStgSz = stgPN[stg_nxt].size(), prvStgSz = stgPN[stg_prv].size(); ///< size of stages i.e. number of different servers.
        float stgMaxDelay = 0;
        // processing the current stage for merging and duplication time
        for(const auto &[pn_idx, cntParallelServer]: stgPN[stg_cur]){

            unsigned int curNxtSamePN=0; if(stgPN[stg_nxt].find(pn_idx) != stgPN[stg_nxt].end()) curNxtSamePN =  1;
            unsigned int cntNextHopDiffServer = nxtStgSz - curNxtSamePN;
            float T_d_pkt = calcTime_InterDuplicationTime(cntNextHopDiffServer);
            float T_d_hdr = calcTime_IntraDuplicationTime(cntParallelServer);

            unsigned int curPrvSamePN=0; if(stgPN[stg_prv].find(pn_idx) != stgPN[stg_prv].end()) curPrvSamePN =  1;
            unsigned int cntPrevHopDiffServer = prvStgSz - curPrvSamePN;
            float T_m_pkt = calcTime_InterMergingTime(cntPrevHopDiffServer);
            float T_m_hdr = calcTime_IntraMergingTime(cntParallelServer);

            pktdist[stg_nxt][pn_idx] = T_m_pkt + T_m_hdr + T_d_hdr + T_d_pkt; ///< overall total time spent in packet processing from src to dest.
            stgMaxDelay = max(stgMaxDelay, pktdist[stg_nxt][pn_idx]);
            if(showInConsole) {
                cout<<"\n    "<<cntParallelServer<<"*( PN:"<<pn_idx<<" ) --> "<<"[prvD:"<<cntPrevHopDiffServer<<" s("<<curPrvSamePN<<") ; nxtD:"<<cntNextHopDiffServer<<" s("<<curNxtSamePN<<")]";
                cout<<"\t[m_pkt:"<<T_m_pkt<<" ; m_hdr:"<<T_m_hdr<<"]"; cout<<" [d_pkt:"<<T_d_pkt<<" ; d_hdr:"<<T_d_hdr<<"]";
            }
        }
        totalPktDelay += stgMaxDelay;
    }
//     Processing the last stage. only merging required.
    pktdist[szStages-1][SFCdst] = totalPktDelay + calcTime_InterMergingTime(stgPN[szStages-2].size()); //<float T_m_pkt

    if(showInConsole) {
        cout<<"\n:STG[DST]: "<<"prvD:"<<stgPN[szStages-2].size(); cout<<" \t[E2E:"<<pktdist[szStages-1][SFCdst]<<"]";
    }
    SFC->pktDist = pktdist;
    if(debug)cout<<"\n<<<[Function Completed: "<<__FUNCTION__<<"]";
//    "("<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";
    return pktdist[szStages-1][SFCdst];
}

/*! Another method to find toposort. Although not required but in case to test.
 * @param SFC ServiceFunctionChain object whose SFC graph VNFs we need to topo sort.
 * @param topologicalOrder  store the node index in topological order
 */
void calculateTopologicalOrder(ServiceFunctionChain *SFC, vector<int>& topologicalOrder){
    unordered_map<int, int> indegree; ///<indegree of each node in the SFC
    queue<int> q; ///< node with indegree 0 will be push into queue
    int cntOfNodes=0; ///< variable to store number of nodes processed to detect cycle in the graph.
    indegree[SFCsrc] = 0;
    ///Calculation of indegree of each node of SFC through adjacency list
    for(const auto& node: SFC->pAdj)  {
        for(int w: SFC->pAdj[node.first])
            indegree[w] = indegree[w]+1;
    }

    for(const auto& node: indegree) { //        cout<<node.first<<"->"<<node.second<<endl;
        if(node.second == 0)
            q.push(node.first);
    }

    while(!q.empty()){
        int node = q.front(); q.pop();
        topologicalOrder.push_back(node);
        for(const auto& w: SFC->pAdj[node]){
            if(--indegree[w] == 0)
                q.push(w);
        }
        cntOfNodes++;
    }

    if(cntOfNodes != 2+SFC->numVNF) {
        cout << "There exists a cycle in the graph\n";
        string errorMsg = "There exists a cycle in the SFC " + SFC->name + ". NodesCounted:" + to_string(cntOfNodes) +
                          ". Function: ";
        throw runtime_error(errorMsg + __FUNCTION__);
    }
}

/*! Calculate the end-to-end delay of given PARALLEL SFC with VNF instances assigned. \n
 * End-to-End delay includes - transmission time, propagation time, execution time, queuing delay, processing time.
 * Don't include same vnf again in same SFC.
 * @param SFC given SFC object whose delay we have to colculate
 * @param allSFC all the SFC in the network, required to calculate arrival rate for queuing delay.
 * @param VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
 * @param VirtualNetwork VirtualMachines object to find the required Physical Node of the given VM id
 * @param PhysicalNetwork PhysicalGraph object to find the propagation delay, min dist between two node
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @tparam type_wgt type_wgt edge weight data type. default=unsigned int.
 * @tparam type_res resource data type. default=unsigned int.
 * @return time (in seconds) End to end delay of SFC.
 */
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
float calcObjectiveValuePar( ServiceFunctionChain *SFC, const vector<ServiceFunctionChain*>& allSFC, VirtualNetworkFunctions<type_res> *VNFNetwork, const VirtualMachines<type_res> *VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork, bool showInConsole=false)
{
//    auto ft_start = std::chrono::steady_clock::now();
    if(debug)cout<<"\n>>[Function Running: "<<__FUNCTION__<<"]";

    vector<int> topologicalOrder;
//    calculateTopologicalOrder(SFC, topologicalOrder);
    /// iterate stage wise and push the stage vnfs first before going to next stage.
    for(const auto& blocks: SFC->vnfBlocksPar){
        for(const auto &node: blocks){
            topologicalOrder.push_back(node);
        }
    }
        if(showInConsole) { cout<<"\nTopological Order: ";for (int x: topologicalOrder) { if(x == SFCsrc) cout<<"SRC -> "; else if(x == SFCdst) cout<<" DST"; else cout << x <<char(96+SFC->I_VNFType2Inst[x])<< " --> "; }  }  ///< Print topological order

    float default_min_distance = -1; ///< default min distance to be used in distance(time) calculation.
    unordered_map<int, float> distance; ///<stores longest distance to each node
    for(const auto& vnfid: topologicalOrder)   distance[vnfid] = default_min_distance; ///< initialization of distance
    distance[SFCsrc] = 0;

    size_t sz = topologicalOrder.size(); ///< vnfSeq.size(); 2 + SFC->numVNF; 2+ for source and destination
    if (SFC->vnfSeq.size() != topologicalOrder.size()) {
        string errorMsg = "Error in accessing SFC. num of VNF in SFC != topologicalOrder.size(). Function:";
        throw runtime_error(errorMsg + __FUNCTION__);
    }

    float T_tx, T_tx_init = calcTime_TransmissionDelay();
    float T_px, T_qd, T_prc, T_exe;

/*! Now we will calculate the longest path in the SFC graph. */
    for (const int &vnf_dst_idx: SFC->pAdj[SFCsrc]){
        if (vnf_dst_idx == SFCdst) {
            string errorMsg = "Error in accessing SFC. vnf_dst_idx: " + to_string(vnf_dst_idx) + ". Function:";
            throw runtime_error(errorMsg + __FUNCTION__);
        }
        VNFNode<type_res> *srcVNFNode = VNFNetwork->VNFNode[vnf_dst_idx];
        T_prc = calcTime_MeanProcessingDelayVNF<type_res>(srcVNFNode);
        T_exe = calcTime_FunctionExecutionDelay<type_res>(srcVNFNode);
        try { T_qd = calcTime_QueuingDelay(vnf_dst_idx, srcVNFNode->serviceRate, SFC ,allSFC);
        } catch (std::exception const &e) { std::cerr << "caught: " << e.what() << std::endl; }
        //  T_prc=0.5; T_exe = 1;
        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
        distance[vnf_dst_idx] = distance[SFCsrc] + T_tx_init + T_prc + T_exe + T_qd;
    }

/// from 1 -> SFCsrc(0) already considered, to sz-2 -> sz-1 is destination
    for(size_t idx = 1; idx <= sz-2; idx++) {
        const int& vnf_src_idx = topologicalOrder[idx];
//        if(vnf_src_idx == SFCdst) continue; // because SFCdst is minus10.
//        if(vnf_src_idx == SFCsrc) {
//            //above for loop
//        }
//        else{ //(vnf_src_idx != SFCsrc)
            int vnf_src_inst_idx = SFC->I_VNFType2Inst[vnf_src_idx];
            int vm_src_idx = VNFNetwork->I_VNFinst2VM[vnf_src_idx][vnf_src_inst_idx];
            int pn_src_idx = VirtualNetwork->I_VM2PN[vm_src_idx];
            if(showInConsole) { cout<<"\n:src: F["<<vnf_src_idx<<char(96+vnf_src_inst_idx)<<"],VM["<<vm_src_idx<<"],PN["<<pn_src_idx<<"]"<<"  ";}
            if(distance[vnf_src_idx] != default_min_distance){
                for(const int& vnf_dst_idx: SFC->pAdj[vnf_src_idx]){
                    if(vnf_dst_idx == SFCdst){
//                            T_prc=0.5; //T_tx using same value as at the time of declaration.
                        if(distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx_init )
                            distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx_init ;
                    }else{ //vnfIDdst != SFCdst
                        int vnf_dst_inst_idx = SFC->I_VNFType2Inst[vnf_dst_idx];
                        int vm_dst_idx = VNFNetwork->I_VNFinst2VM[vnf_dst_idx][vnf_dst_inst_idx];
                        int pn_dst_idx = VirtualNetwork->I_VM2PN[vm_dst_idx];
                        VNFNode<type_res> *dstVNFNode = VNFNetwork->VNFNode[vnf_dst_idx];
                        if(pn_src_idx == pn_dst_idx){
                            T_tx=0; T_px=0;
                        }else {
                            T_tx = T_tx_init;
                            T_px = calcTime_PropagationDelay<type_wgt, type_res>(pn_src_idx, pn_dst_idx, PhysicalNetwork);
                        }
                        T_prc = calcTime_MeanProcessingDelayVNF<type_res>(dstVNFNode);
                        T_exe = calcTime_FunctionExecutionDelay<type_res>(dstVNFNode);
                        try{ T_qd=calcTime_QueuingDelay(vnf_dst_idx, dstVNFNode->serviceRate, SFC ,allSFC);
                        } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }
//                    T_tx = 1; T_px = 1; T_qd = 1; T_prc = 1; T_exe = 1;
                        if(showInConsole) {  cout<<" :dst:F["<<vnf_dst_idx<<char(96+vnf_dst_inst_idx)<<"],VM["<<vm_dst_idx<<"],PN["<<pn_dst_idx<<"]"<<"\t" <<"T_tx["<<T_tx<<"] | T_px["<<T_px<<"] | T_qd["<<T_qd<<"]"<<"] | T_prc["<<T_prc<<"]"<<"] | T_exe["<<T_exe<<"]"<<"\n";}
                        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
                        if(distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe)
                            distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe;
                    }  //vnfIDdst != SFCdst
                }   // for source child
            } // if distance != min
//        } // else not SFCsrc
    }//for loop oid each vnf in topological order.

    if(showInConsole) { cout<<endl<<"Distance:\n";  for(const auto& d: topologicalOrder){ cout<<"f["<<d<<"]:"<<distance[d]<<" |  "; }}
    SFC->distancePar = distance;
    if(debug)cout<<"\n<<[Function Completed: "<<__FUNCTION__<<"]";
//    "("<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";
    return distance[SFCdst];
}

/*! Calculate the end-to-end delay of given SEQUENTIAL SFC with VNF instances assigned. \n
 * End-to-End delay includes - transmission time, propagation time, execution time, queuing delay, processing time.
 * Don't include same vnf again in same SFC.
 * @param SFC given SFC object whose delay we have to colculate
 * @param allSFC all the SFC in the network, required to calculate arrival rate for queuing delay.
 * @param VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
 * @param VirtualNetwork VirtualMachines object to find the required Physical Node of the given VM id
 * @param PhysicalNetwork PhysicalGraph object to find the propagation delay, min dist between two node
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @tparam type_wgt type_wgt edge weight data type. default=unsigned int.
 * @tparam type_res resource data type. default=unsigned int.
 * @return time (in seconds) End to end delay of SFC.
 */
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
float calcObjectiveValueSeq( ServiceFunctionChain *SFC, const vector<ServiceFunctionChain*>& allSFC, VirtualNetworkFunctions<type_res> *VNFNetwork, const VirtualMachines<type_res> *VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork, bool showInConsole=false)
{
//    auto ft_start = std::chrono::steady_clock::now();
    if(debug)cout<<"\n>[Function Running: "<<__FUNCTION__<<"]";
//    if(showInConsole) { cout<<"\nTopological Order: ";for (int x: SFC->vnfSeq) { if(x == SFCsrc) cout<<"SRC -> "; else if(x == SFCdst) cout<<" DST"; else cout << x <<char(96+SFC->I_VNFType2Inst[x])<< " --> "; } }  ///< Print topological order

    size_t sz = SFC->vnfSeq.size(); // src + dest = (2)
    float default_min_distance = -1; ///< default min distance to be used in distance(time) calculation.
    unordered_map<int, float> distance; ///<stores longest distance to each node
    for(const auto& vnfid: SFC->vnfSeq)   distance[vnfid] = default_min_distance; ///< initialization of distance
    distance[SFCsrc] = 0;

    float T_tx, T_tx_init = calcTime_TransmissionDelay();
    float T_px, T_qd, T_prc, T_exe;

/*! Now we will calculate the longest path in the SFC graph. */
    for (const int &vnf_dst_idx: SFC->sAdj[SFCsrc]){ // as it is a sequential SFC only one node is there after SFCsrc level but just for precaution
        if (vnf_dst_idx == SFCdst) {
            string errorMsg = "Error in accessing SFC. vnf_dst_idx: " + to_string(vnf_dst_idx) + ". Function:";
            throw runtime_error(errorMsg + __FUNCTION__);
        }
        VNFNode<type_res> *srcVNFNode = VNFNetwork->VNFNode[vnf_dst_idx];
        T_prc = calcTime_MeanProcessingDelayVNF<type_res>(srcVNFNode);
        T_exe = calcTime_FunctionExecutionDelay<type_res>(srcVNFNode);
        try { T_qd = calcTime_QueuingDelay(vnf_dst_idx, srcVNFNode->serviceRate, SFC ,allSFC);
        } catch (std::exception const &e) { std::cerr << "caught: " << e.what() << std::endl; }
//        T_prc=0.5; T_exe = 1;
        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
        distance[vnf_dst_idx] = distance[SFCsrc] + T_tx_init + T_prc + T_exe + T_qd;
    }
    /// from 1 -> SFCsrc(0) already considered, to sz-2 -> sz-1 is destination
    for(size_t idx = 1; idx <= sz-2; idx++) {
        const int& vnf_src_idx = SFC->vnfSeq[idx];
//        if(vnf_src_idx == SFCdst) continue; // because SFCdst is minus10.
//        if(vnf_src_idx == SFCsrc) {
//             //above for loop
//        }
//        else{ //(vnf_src_idx != SFCsrc)
            int vnf_src_inst_idx = SFC->I_VNFType2Inst[vnf_src_idx];
            int vm_src_idx = VNFNetwork->I_VNFinst2VM[vnf_src_idx][vnf_src_inst_idx];
            int pn_src_idx = VirtualNetwork->I_VM2PN[vm_src_idx];
            if(showInConsole) { cout<<"\n:src: F["<<vnf_src_idx<<char(96+vnf_src_inst_idx)<<"],VM["<<vm_src_idx<<"],PN["<<pn_src_idx<<"]"<<" :dst:";}
            if(distance[vnf_src_idx] != default_min_distance){
                for(const int& vnf_dst_idx: SFC->sAdj[vnf_src_idx]){
                    if(vnf_dst_idx == SFCdst){
//                            T_prc=0.5; //T_tx using same value as at the time of declaration.
                        if(distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx_init )
                            distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx_init ;
                    }else{ //vnfIDdst != SFCdst
                        int vnf_dst_inst_idx = SFC->I_VNFType2Inst[vnf_dst_idx];
                        int vm_dst_idx = VNFNetwork->I_VNFinst2VM[vnf_dst_idx][vnf_dst_inst_idx];
                        int pn_dst_idx = VirtualNetwork->I_VM2PN[vm_dst_idx];
                        VNFNode<type_res> *dstVNFNode = VNFNetwork->VNFNode[vnf_dst_idx];
                        if(pn_src_idx == pn_dst_idx){
                            T_tx=0; T_px=0;
                        }else{
                            T_tx = T_tx_init;
                            T_px = calcTime_PropagationDelay<type_wgt, type_res>(pn_src_idx, pn_dst_idx, PhysicalNetwork);
                        }
                        T_prc = calcTime_MeanProcessingDelayVNF<type_res>(dstVNFNode);
                        T_exe = calcTime_FunctionExecutionDelay<type_res>(dstVNFNode);
                        try{ T_qd=calcTime_QueuingDelay(vnf_dst_idx, dstVNFNode->serviceRate, SFC ,allSFC);
                        } catch( std::exception const& e ) { std::cerr << "caught: " << e.what() << std::endl; }
//                    T_tx = 1; T_px = 1; T_qd = 1; T_prc = 1; T_exe = 1;
                        if(showInConsole) {  cout<<" F["<<vnf_dst_idx<<char(96+vnf_dst_inst_idx)<<"],VM["<<vm_dst_idx<<"],PN["<<pn_dst_idx<<"]"<<"\t" <<"T_tx["<<T_tx<<"] | T_px["<<T_px<<"] | T_qd["<<T_qd<<"]"<<"] | T_prc["<<T_prc<<"]"<<"] | T_exe["<<T_exe<<"]";}
                        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
                        if(distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe)
                            distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe;
                    }  //vnfIDdst != SFCdst
                }   // for source child
            } // if distance != min
//        } // else not SFCsrc
    }//for loop oid each vnf in topological order.

    if(showInConsole) { cout<<endl<<"Distance:\n";  for(const auto& d: SFC->vnfSeq){ cout<<"f["<<d<<"]: "<<distance[d]<<" |\t "; }}
    SFC->distanceSeq = distance;
    if(debug)cout<<"\n<[Function Completed: "<<__FUNCTION__<<"]";
//    "("<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ft_start).count()<<"ms)";
    return distance[SFCdst];
}


#endif //SFC_PARALLELIZATION_TIMECALCULATIONFUNCTIONS_H