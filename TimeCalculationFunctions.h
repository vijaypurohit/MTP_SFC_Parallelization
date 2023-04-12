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
 * @return transmission time in type_delay. in seconds.\n
 *  random value = 0.0008192 seconds ((1024 * 8 b) / (10 * 10^6 b/s)) \n
 *  min = 0.0000000001 s ( 1 b/ 10000 * 10^6 ). \n
 *  max = 0.1024 s (10240/ 10^6).
 *  @future packet size can be variable according to each SFC.
 */
type_delay calcTime_TransmissionDelay([[maybe_unused]] unsigned int n1=0, [[maybe_unused]] unsigned int n2=0,
                                 [[maybe_unused]] ServiceFunctionChain *SFC = nullptr){
    auto bandwidth_n1_n2 = (type_delay)bandwidthNW * factor_bandwidth; // in Mb/s
    auto packetsize_n1_n2 = ((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet; // in bits
    return packetsize_n1_n2/(bandwidth_n1_n2);
}

/*! Calculate the propagation delay between two Physical Nodes n1 and n2
 * (time duration taken for a signal to reach its destination)\n
 * Formula = (Distance between two nodes n1 and n2) / Velocity of light in the medium \n
 *        = Dist[n1][n2] / (transmission medium * speed of light in vaccum) \n
 * @refer velocityFactor, speedOfLight and dist of Physical Graph
 * @param n1 Physical Node n1 @param n2 Physical Node n2
 * @param PhysicalNetwork object of class type PhysicalGraph.
 * @return propagation time in type_delay. in seconds. \n
 *  random value =  0.00001 second (3000 m / (1.0 * 3 * 10^8)). \n
 *  min = 3.33e-9 s ( 1 m/ 3*10^8 ). \n
 *  max = 0.04252 s (12756000 m/ 3*10^8).
 */
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
type_delay calcTime_PropagationDelay(unsigned int n1, unsigned int n2, const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork){
    auto len_n1_n2 = (type_delay)PhysicalNetwork->dist[n1][n2];
    auto speed = velocityFactor * (type_delay)speedOfLight;
    return (len_n1_n2 / (speed));
}

/*! Calculate the Queuing delay of a particular VNF (time a job waits in a queue until it can be executed) \n
 * Formula = lambda_c / (mu_f * (mu_f - lambda_c)) \n
 * where lambda_c is traffic arrival rate of chain c or SFC s  and mu_f is the service rate of SFC or say VNF. \n
 * Both in packets per second. \n
 * @param curVNFid VNF Node Id whose queuing delay we have to calculate
 * @param mu_f service rate of current VNF.
 * @param X_VNFType2Instv In curSFC VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param curSFC SFC class whose VNF we are considering
 * @param allSFC all the SFC in the network
 * @return queueing delay in seconds.
 * Sabse phele ham sab SFC me dekhenge ki ye curVNF present hai kya, and curVNF ka jo instance use kr rhe wohi hai to add kr denge us SFC ka arrival rate
 */
type_delay calcTime_QueuingDelay1(const unsigned int& curVNFid, const type_delay &mu_f,    const unordered_map<unsigned int, unsigned int>& X_VNFType2Inst,
                            const unsigned int curSFCindex, const vector<ServiceFunctionChain*>& allSFC){
    type_delay lambda = 0;
    for(unsigned int c=1; c<allSFC.size(); c++){ //0th SFC is un-assigned,
        if(allSFC[c]->index == curSFCindex)
            lambda += allSFC[c]->trafficArrivalRate;
        else if( (allSFC[c]->I_VNFType2Inst.count(curVNFid) != 0) and
                 (allSFC[c]->I_VNFType2Inst.at(curVNFid) == X_VNFType2Inst.at(curVNFid))  ) /// is curVNF of curSFC present in this SFC? and that VNF instance is same as curVNF instance?
        {
                lambda += allSFC[c]->trafficArrivalRate;
        }
    }
//    cout<<"lambda += VNF["<<curVNFid<<char(96+X_VNFType2Inst.at(curVNFid))<<"] = "<<lambda<<"\t";
    if(lambda >= mu_f){
//        return std::numeric_limits<type_delay>::max();
        string errorMsg = "SFC_"+ to_string(curSFCindex)+"->F["+ to_string(curVNFid)+char(96+X_VNFType2Inst.at(curVNFid))+"]. ("+ to_string(lambda)+">="+ to_string(mu_f)+") (SFCs Arrival Rate >= VNF processing rate). Function: ";
        throw runtime_error(errorMsg + __FUNCTION__);
    }
    return (lambda / (mu_f * (mu_f - lambda)));
}

/*! Calculate the Queuing delay of a particular VNF (time a job waits in a queue until it can be executed) \n
 * Formula = lambda_c / (mu_f * (mu_f - lambda_c)) \n
 * where lambda_c is traffic arrival rate of chain c or SFC s  and mu_f is the service rate of SFC or say VNF. \n
 * Both in packets per second. \n
 * @param fnNode VNF Node object whose Queuing delay we have to calculate
 * @param fnInstId In curSFC, VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param oldUtilization old utilisation of the current function type and instance till now.
 * @param cSFC SFC object whose VNF we are considering
 * @return queuing delay in seconds.
 */
template<typename type_res=unsigned int>
type_delay calcTime_QueuingDelay(const VNFNode<type_res> *fnNode, const unsigned int& fnInstId,
                            const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization, const ServiceFunctionChain* cSFC){
    const unsigned int& fnType =  fnNode->index;

    type_delay mu_f = fnNode->serviceRate; ///<service rate of VNF Node.

    type_delay lambda =  cSFC->trafficArrivalRate ; ///< current rate
    if(oldUtilization.count(fnType) and oldUtilization.at(fnType).count(fnInstId)) ///< old utilization till now of VNF
        lambda += oldUtilization.at(fnType).at(fnInstId);

    if(lambda > mu_f){// otherwise queuing delay will be negative
        return std::numeric_limits<type_delay>::max();  // we want to punish this mapping such that with this mapping time is maximum
    }
    else if(lambda == mu_f){// otherwise queuing delay will be infinite,but we can still maximise it to 100% utilization
        return 1;  //
    }
    return (lambda / (mu_f * (mu_f - lambda)));
}

/*! Return the time required to complete the execution of the function VNF f. \n
 * @param VNFNode for function execution time.
 * @return execution time in seconds.
 */
template<typename type_res=unsigned int>
type_delay calcTime_FunctionExecutionDelay(VNFNode<type_res> *VNFNode){
    auto exeTime = (type_delay)VNFNode->executionTime;
    return exeTime;
}

/*! Calculate the mean processing delay of particular VNF, time it takes function to process the one packet. \n
 * Fromula = 1/serviceRate of VNF.
 * @param VNFNode VNF Node class for service rate
 * @return mean processing delay in seconds.
 */
template<typename type_res=unsigned int>
type_delay calcTime_MeanProcessingDelayVNF(VNFNode<type_res>* VNFNode){
    auto mu_f = (type_delay)VNFNode->serviceRate;
    return 1/mu_f;
}


/*! The time required by a server v to duplicate the whole packet for instances
 * installed in d different servers of the next step.\n
 * Formula = T_d_pkt * (d-1);
 * @return INTER duplication time in seconds.
 * For example; packetBodySize = 1000 B; packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_pkt=0.630784(d=2),  =1.26157(d=3),  =1.892352(d=4)
 */
type_delay calcTime_InterDuplicationTime(unsigned int cntNextHopDiffServer){
    if( cntNextHopDiffServer <= 1) return 0;
    auto packetsize_n1_n2 = ((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet;
    type_delay time_to_duplicate_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_duplicate_packet*(type_delay)(cntNextHopDiffServer-1);
}
/*! the time required by a server v to duplicate only the header for its p parallel instances. \n
 * Formula = T_d_hdr * (p-1);
 * @return INTRA duplication time in seconds.
 * For example;packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_hdr=0.014784(p=2),  =0.029568(p=3),  =0.044352(p=4)
 */
type_delay calcTime_IntraDuplicationTime(unsigned int cntParallelServer){
    if(cntParallelServer <= 1) return 0;
    type_delay time_to_duplicate_header = (type_delay)packetHeaderSize *factor_packet * read_write_time_per_bit;
    return time_to_duplicate_header*(type_delay)(cntParallelServer-1);
}
/*! the time required by a server v to merge the whole packets from d different servers of the previous step, \n
 * Formula = T_m_pkt * (d-1);
 * @return INTER merging time in seconds.
 */
type_delay calcTime_InterMergingTime(unsigned int cntPrevHopDiffServer){
    if(cntPrevHopDiffServer <= 1) return 0;
    auto packetsize_n1_n2 = ((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet;
    type_delay time_to_merging_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_merging_packet*(type_delay)(cntPrevHopDiffServer-1);
}
/*! the time required by a server v to merge the headers from its p parallel instances \n
 * Formula = T_m_hdr * (p-1);
 * @return INTRA merging time in seconds.
 */
type_delay calcTime_IntraMergingTime(unsigned int cntParallelServer){
    if(cntParallelServer <= 1) return 0;
    type_delay time_to_merge_header = (type_delay)packetHeaderSize *factor_packet * read_write_time_per_bit;
    return time_to_merge_header*(type_delay)(cntParallelServer-1);
}

/*! Calculate time it takes for packet duplication, packet merging for each physical SERVER at each level/stage. \n
 * for example: Physical Nodes at each level = {SRC--> {3, 3, 4}, {3,4,4,5,5,6}, {3,3,4,7,6},-->DST }; \n
 * SRC -> next level. requires only inter duplication AND secLastLevel --> LastLevel requires only inter merging time \n
 * at each level requires -> inter merging, intra duplication, intra merging, inter duplication
 * @param givenParSFC given enumeration of vnfBlockPar for which we have to calcualte the delay (except src dst block)
 * @param X_VNFType2Instv In curSFC VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param curSFCidx given SFC object whose delay we have to colculate
 * @param VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
 * @param VirtualNetwork  VirtualMachines object to find the required Physical Node of the given VM id
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @tparam type_res resource data type. default=unsigned int.
 * @return packet delay of SFC in seconds. If not reachable (<0) then max time(type_delay value).
 */
template<typename type_res=unsigned int>
type_delay calcTime_PacketsDelay(const vector<vector<unsigned int>>& givenParSFC, const unordered_map<unsigned int, unsigned int>& X_VNFType2Inst,
                            const unsigned int& curSFCidx, const VirtualNetworkFunctions<type_res> *VNFNetwork,
                            const VirtualMachines<type_res> *VirtualNetwork, bool showInConsole=false){
//    givenParSFC = { {3, 3, 4}, {3,4,4,5,5,6}, {3,3,4,7,6} };
    if(debug and showInConsole) cout<<"\n>>>[Function Running: "<<__FUNCTION__<<"]";
    unsigned int szStages = givenParSFC.size();
    /*! cSFCsrc is dummy state, next stage to process is first parallel block  i.e index 0 of parBlockArray\n
     *                      [src] -> [f1] -> [f6 f4] -> [f5] -> [dst]\n
     * arr index(stg_nxt-1)           0        1         2     szStages=3 (excluding src and dst)\n
     *          stg_nxt      0        1        2         3       4
     *                    stg_prv  stg_cur    stg_nxt
     */
    unsigned int stg_nxt=1;

    vector<unordered_map<unsigned int, unsigned int>> stgPN(szStages+2); ///< all the DIFFERENT physical nodes {nodeid, freq} belonging to each stages (+2 for src and dst stage).
    type_delay totalPktDelay; // total max delay of SFC from src to dst incurred in a process to duplicate and merge packets.
//    unordered_map<unsigned int, unordered_map<int, type_delay>> pktdist; ///<stores paket delay for each PHYSICAL NODE {stg -> {pnId, time}}

    /// SFC src --> STG[1] This for-loop is Only For Source Level to Next Level(index 0 of parallel block), only duplication required.
    for (const int &vnf_dst_idx: givenParSFC[stg_nxt-1]){  //int vnf_dst_inst_idx =0, vm_dst_idx=0,pn_dst_idx =vnf_dst_idx;
        const auto & vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx); const auto & vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);const auto & pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
        stgPN[stg_nxt][pn_dst_idx] += 1;
    }
    /// type_delay T_d_pkt = Time to duplicate packets to d diff server from level_i to level_i+1, Here only duplication time required.
    totalPktDelay = calcTime_InterDuplicationTime(stgPN[stg_nxt].size()); //totalPktDelay = pktdist[stg_nxt-1][SFCsrc] = calcTime_InterDuplicationTime(stgPN[stg_nxt].size());

    if(showInConsole) {
        cout<<"\n:STG[0]:-->:STG[1]: "<<stgPN[stg_nxt].size()<<"("; for(const auto& pn: stgPN[stg_nxt]) cout<<" PN["<<pn.first<<":"<<pn.second<<"]";  cout<<" )\t[d_pkt:"<<totalPktDelay<<"]";
    }

    /*!  From STG = 2nd Parallel stage  to szStages-1 2nd Parallel stage,
     *   each loop first calculates the information of next stage and then process the current stage.
     *   At last, second last stage is info is calculated but not processed.
     */
    for(  stg_nxt=2; stg_nxt <= szStages+1; stg_nxt++) {
        unsigned int stg_cur = stg_nxt-1;
        unsigned int stg_prv = stg_cur-1;

    // processing the next level stgi+1, except last level
        if(stg_nxt != szStages+1){
        for (const auto &vnf_dst_idx: givenParSFC[stg_nxt-1]) { //int vnf_dst_inst_idx =0,vm_dst_idx=0,pn_dst_idx =vnf_dst_idx;
            const auto & vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);
            const auto & vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);
            const auto & pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
            stgPN[stg_nxt][pn_dst_idx] += 1;
        }}

        if(showInConsole) {
            cout<<"\n:STG["<<stg_cur<<"]:-->:STG["<<stg_nxt<<"]: "<<stgPN[stg_nxt].size()<<"("; for(const auto& pn: stgPN[stg_nxt]) cout<<" PN["<<pn.first<<":"<<pn.second<<"]"; cout<<" )";
        }

        unsigned int nxtStgSz = stgPN[stg_nxt].size(), prvStgSz = stgPN[stg_prv].size(); ///< size of stages i.e. number of different servers.
        type_delay stgMaxDelay = 0;
        // processing the current stage for merging and duplication time
        for(const auto &[pn_idx, cntParallelServer]: stgPN[stg_cur]){

            unsigned int curNxtSamePN=0; if(stgPN[stg_nxt].find(pn_idx) != stgPN[stg_nxt].end()) curNxtSamePN =  1;
            unsigned int cntNextHopDiffServer = nxtStgSz - curNxtSamePN;
            type_delay T_d_pkt = calcTime_InterDuplicationTime(cntNextHopDiffServer);
            type_delay T_d_hdr = calcTime_IntraDuplicationTime(cntParallelServer);

            unsigned int curPrvSamePN=0; if(stgPN[stg_prv].find(pn_idx) != stgPN[stg_prv].end()) curPrvSamePN =  1;
            unsigned int cntPrevHopDiffServer = prvStgSz - curPrvSamePN;
            type_delay T_m_pkt = calcTime_InterMergingTime(cntPrevHopDiffServer);
            type_delay T_m_hdr = calcTime_IntraMergingTime(cntParallelServer);

            //pktdist[stg_nxt][pn_idx] = T_m_pkt + T_m_hdr + T_d_hdr + T_d_pkt;
            stgMaxDelay = max(stgMaxDelay,  T_m_pkt + T_m_hdr + T_d_hdr + T_d_pkt);///< overall total time spent in packet processing from src to dest.
            if(showInConsole) {
                cout<<"\n    "<<cntParallelServer<<"*( PN:"<<pn_idx<<" ) --> "<<"[prvD:"<<cntPrevHopDiffServer<<" s("<<curPrvSamePN<<") ; nxtD:"<<cntNextHopDiffServer<<" s("<<curNxtSamePN<<")]";
                cout<<"\t[m_pkt:"<<T_m_pkt<<" ; m_hdr:"<<T_m_hdr<<"]"; cout<<" [d_pkt:"<<T_d_pkt<<" ; d_hdr:"<<T_d_hdr<<"]";
            }
        }
        totalPktDelay += stgMaxDelay;
    }

//     Processing the last stage. only merging required.
    totalPktDelay += calcTime_InterMergingTime(stgPN[szStages].size()); //<type_delay T_m_pkt // pktdist[szStages+1][SFCdst] =

    if(showInConsole) {
        cout<<"\n:STG[DST]: "<<"prvD:"<<stgPN[szStages].size(); cout<<" \t[E2E:"<<totalPktDelay<<"]";
    }
    if(debug and showInConsole)cout<<"\n<<<[Function Completed: "<<__FUNCTION__<<"]";

    if ( totalPktDelay < 0) { // if distance obtanied is negative (unreachable or -1 default min distance). delay can be 0 in case of no packet processing.
        return std::numeric_limits<type_delay>::max();
    }
    return totalPktDelay;
}//calcTime_PacketsDelay


/*! Calculate the end-to-end delay of given PARALLEL SFC with VNF instances assigned. \n
 * End-to-End delay includes - transmission time, propagation time, execution time, queuing delay, processing time.
 * Don't include same vnf again in same SFC.
 * @param givenParSFC given enumeration of vnfBlockPar of cSFC for which we have to calcualte the delay (except src dst block)
 * @param X_VNFType2Instv In curSFC VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param curSFCidx given SFC index from which we get the mapping and given SFC detail
 * @param allSFC all the SFC in the network, required to calculate arrival rate for queuing delay.
 * @param VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
 * @param VirtualNetwork VirtualMachines object to find the required Physical Node of the given VM id
 * @param PhysicalNetwork PhysicalGraph object to find the propagation delay, min dist between two node
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @tparam type_wgt type_wgt edge weight data type. default=unsigned int.
 * @tparam type_res resource data type. default=unsigned int.
 * @return time (in seconds) End to end delay of SFC. If not reachable (<=0) then max time(type_delay value).
 */
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
type_delay calcObjectiveValuePar(const vector<vector<unsigned int>>& givenParSFC,   const unordered_map<unsigned int,unsigned int>& X_VNFType2Inst,
                            const unsigned int &curSFCidx, const vector<ServiceFunctionChain*>& allSFC,
                            const VirtualNetworkFunctions<type_res> *VNFNetwork, const VirtualMachines<type_res> *VirtualNetwork,
                            const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork, bool showInConsole=false, bool showInConsoleDetailed=false)
{
    if(debug and showInConsoleDetailed)cout<<"\n>>[Function Running: "<<__FUNCTION__<<"]";

    type_delay default_min_distance = -1; ///< default min distance to be used in distance(time) calculation.
    unordered_map<int, type_delay> distance; ///<stores longest distance to each node
    for(const auto& vnfid: allSFC[curSFCidx]->vnfSeq)   distance[vnfid] = default_min_distance; // jitne bhi function hai unko min dist se intialise kr do
    distance[SFCsrc] = 0;

    type_delay T_tx, T_tx_init = calcTime_TransmissionDelay();
    type_delay T_px, T_qd, T_prc, T_exe;

    unsigned int szStages = givenParSFC.size();
    /*! cSFCsrc is dummy state, next stage to process is first parallel block  i.e index 0 of parBlockArray\n
     *                      [src] -> [f1] -> [f6 f4] -> [f5] -> [dst]\n
     * arr index(stg_nxt-1)           0        1         2     szStages=3 (excluding src and dst)\n
     *          stg_nxt      0        1        2         3       4
     *                    stg_prv  stg_cur    stg_nxt
     */
    unsigned int stg_nxt=1;

/*! Now we will calculate the longest path in the SFC graph. */
    for (const auto &vnf_dst_idx: givenParSFC[stg_nxt-1]){ // process first parallel block, vnfsrc is SFCsrc
        if (vnf_dst_idx == SFCdst) {
            string errorMsg = "Error in accessing SFC. vnf_dst_idx: " + to_string(vnf_dst_idx) + ". Function:";
            throw runtime_error(errorMsg + __FUNCTION__);
        }
        VNFNode<type_res> *dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
        T_prc = calcTime_MeanProcessingDelayVNF<type_res>(dstVNFNode);
        T_exe = calcTime_FunctionExecutionDelay<type_res>(dstVNFNode);
        T_qd = calcTime_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), VNFNetwork->utilization,  allSFC[curSFCidx]);
//        try { T_qd = calcTime_QueuingDelay1(vnf_dst_idx, dstVNFNode->serviceRate, X_VNFType2Inst,curSFCidx,allSFC);
//        } catch (std::exception const &e) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
//          T_prc=0.5; T_exe = 1; T_qd=2;
        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
        distance[vnf_dst_idx] = distance[SFCsrc] + T_tx_init + T_prc + T_exe + T_qd;
    }

    /*!  Stage Wise execution to maintain topological order
     *  for each vnfsrc in cur stage we will process next stage vnfs except last stage (szStages-1).
     */
    for( stg_nxt=1; stg_nxt<=szStages; stg_nxt++){

        for(const auto &vnf_src_idx: givenParSFC[stg_nxt-1]) { //for each vnf in src stage
            const auto & vnf_src_inst_idx = X_VNFType2Inst.at(vnf_src_idx);
            const auto & vm_src_idx = VNFNetwork->I_VNFinst2VM.at(vnf_src_idx).at(vnf_src_inst_idx);
            const auto & pn_src_idx = VirtualNetwork->I_VM2PN.at(vm_src_idx);
            if (showInConsoleDetailed) {
                cout << "\n:src: F[" << vnf_src_idx << char(96 + vnf_src_inst_idx) << "],VM[" << vm_src_idx << "],PN["<< pn_src_idx << "]" << "  ";
            }

           if(stg_nxt == szStages){ // process last stage
               if (distance[SFCdst] < distance[vnf_src_idx] + T_tx_init)
                   distance[SFCdst] = distance[vnf_src_idx] + T_tx_init;
           }  else if (distance[vnf_src_idx] != default_min_distance) { // else process intermediate nodes

                for (const auto &vnf_dst_idx: givenParSFC[stg_nxt]) {  //stg_nxt + 1 - 1, stg_nxt next stage, -1 is the id
                        const auto & vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);
                        const auto & vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);
                        const auto & pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
                        VNFNode<type_res> *dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);

                        if (pn_src_idx == pn_dst_idx) {  T_tx = 0; T_px = 0;
                        } else {  T_tx = T_tx_init; T_px = calcTime_PropagationDelay<type_wgt, type_res>(pn_src_idx, pn_dst_idx, PhysicalNetwork);
                        }
                        T_prc = calcTime_MeanProcessingDelayVNF<type_res>(dstVNFNode);
                        T_exe = calcTime_FunctionExecutionDelay<type_res>(dstVNFNode);
                        T_qd = calcTime_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), VNFNetwork->utilization, allSFC[curSFCidx]);
//                        try {  T_qd = calcTime_QueuingDelay1(vnf_dst_idx, dstVNFNode->serviceRate, X_VNFType2Inst,curSFCidx, allSFC);
//                        } catch (std::exception const &e) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
//                    T_tx = 1; T_px = 1; T_qd = 1; T_prc = 1; T_exe = 1; T_prc=0.5; T_exe = 1; T_qd=2;
                        if (showInConsoleDetailed) {
                            cout << "\n\t :dst:F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "],VM[" << vm_dst_idx
                                 << "],PN[" << pn_dst_idx << "]" << "\t" << "T_tx[" << T_tx << "] | T_px[" << T_px
                                 << "] | T_qd[" << T_qd << "] | T_prc[" << T_prc << "] | T_exe[" << T_exe
                                 << "]" << "\n";
                        }
                        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
                        if (distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe)
                            distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe;

                }   // for each child of source
            } // if distance != min
        } //for each src vnf
    }///<Stage Wise execution to maintain topological order

    if(showInConsole) { cout<<endl<<"Distance:\n";  for(const auto& d: allSFC[curSFCidx]->vnfSeq){ cout<<"f["<<d<<"]:"<<distance[d]<<" |  "; }}

    if(debug and showInConsoleDetailed)cout<<"\n<<[Function Completed: "<<__FUNCTION__<<"]";

    if ( distance[SFCdst] <= 0) { // if distance obtanied is 0 or negative (unreachable or -1 default min distance)
        return std::numeric_limits<type_delay>::max();
    }
    return distance[SFCdst];
}//calcObjectiveValuePar

/*! Calculate the end-to-end delay of given SEQUENTIAL SFC (SFC->vnfSeq) with VNF instances assigned. \n
 * End-to-End delay includes - transmission time, propagation time, execution time, queuing delay, processing time.
 * Don't include same vnf again in same SFC.
 * @param givenSFC given SFC object whose delay we have to colculate. givenSFC->vnfSeq (including src and dest) array used.
 * @param X_VNFType2Instv In curSFC VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param allSFC all the SFC in the network, required to calculate arrival rate for queuing delay.
 * @param VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
 * @param VirtualNetwork VirtualMachines object to find the required Physical Node of the given VM id
 * @param PhysicalNetwork PhysicalGraph object to find the propagation delay, min dist between two node
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @tparam type_wgt type_wgt edge weight data type. default=unsigned int.
 * @tparam type_res resource data type. default=unsigned int.
 * @return time (in seconds) End to end delay of SFC. If not reachable (<=0) then max time(type_delay value).
 */
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
type_delay calcObjectiveValueSeq(const ServiceFunctionChain *givenSFC, const unordered_map<unsigned int, unsigned int>& X_VNFType2Inst,
                            const vector<ServiceFunctionChain*>& allSFC, const VirtualNetworkFunctions<type_res> *VNFNetwork,
                            const VirtualMachines<type_res> *VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *PhysicalNetwork,
                            bool showInConsole=false){

    if(debug and showInConsole)cout<<"\n>[Function Running: "<<__FUNCTION__<<"]";
//    if(showInConsole) { cout<<"\nTopological Order: ";for (int x: givenSFC->vnfSeq) { if(x == SFCsrc) cout<<"SRC -> "; else if(x == SFCdst) cout<<" DST"; else cout << x <<char(96+X_VNFType2Inst[x])<< " --> "; } }  ///< Print topological order

    unsigned int sz = givenSFC->vnfSeq.size(); // src + dest = (2)

    type_delay default_min_distance = -1; ///< default min distance to be used in distance(time) calculation.
    unordered_map< int, type_delay> distance; ///<stores longest distance to each node
    for(const auto& vnfid: givenSFC->vnfSeq)   distance[vnfid] = default_min_distance; ///< initialization of distance
    distance[SFCsrc] = 0;

    type_delay T_tx, T_tx_init = calcTime_TransmissionDelay();
    type_delay T_px, T_qd, T_prc, T_exe;

/*! Now we will calculate the longest path in the SFC graph. */
     // as it is a sequential givenSFC only one node is there after SFCsrc. level(0) is SFCsrc
    {
        const auto &vnf_src_idx = givenSFC->vnfSeq[1]; //first parallel block with index = 1, index 0 is src
        if (vnf_src_idx == SFCdst) {
            string errorMsg = "Error in accessing SFC. vnf_dst_idx: " + to_string(vnf_src_idx) + ". Function:";
            throw runtime_error(errorMsg + __FUNCTION__);
        }
        VNFNode<type_res> *srcVNFNode = VNFNetwork->VNFNodes.at(vnf_src_idx);
        T_prc = calcTime_MeanProcessingDelayVNF<type_res>(srcVNFNode);
        T_exe = calcTime_FunctionExecutionDelay<type_res>(srcVNFNode);

        T_qd = calcTime_QueuingDelay<type_res>(srcVNFNode, X_VNFType2Inst.at(vnf_src_idx), VNFNetwork->utilization, allSFC[givenSFC->index]);
//        try { T_qd = calcTime_QueuingDelay1(vnf_src_idx, srcVNFNode->serviceRate, X_VNFType2Inst,givenSFC->index, allSFC);
//        } catch (std::exception const &e) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
//        T_prc=0.5; T_exe = 1;  T_qd=2;
        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
        distance[vnf_src_idx] = distance[SFCsrc] + T_tx_init + T_prc + T_exe + T_qd;
    }

    /// from 1 as SFCsrc(0) already considered and first parallel block index=1, TO as sz-2 as src -> sz-1 is destination SFCdst
    for(unsigned int idx = 1; idx <= sz-2; idx++) {
        const auto& vnf_src_idx = givenSFC->vnfSeq[idx]; // int qki -10 bhi hai -10 as dst
        const auto& vnf_src_inst_idx = X_VNFType2Inst.at(vnf_src_idx);
        const auto& vm_src_idx = VNFNetwork->I_VNFinst2VM.at(vnf_src_idx).at(vnf_src_inst_idx);
        const auto& pn_src_idx = VirtualNetwork->I_VM2PN.at(vm_src_idx);

            if(showInConsole) { cout<<"\n:src: F["<<vnf_src_idx<<char(96+vnf_src_inst_idx)<<"],VM["<<vm_src_idx<<"],PN["<<pn_src_idx<<"]"<<" :dst:";}
            if(distance[vnf_src_idx] != default_min_distance)
            {
                const int& vnf_dst_idx= givenSFC->vnfSeq[idx+1];
                if(vnf_dst_idx == SFCdst){ //    T_prc=0.5; //T_tx using same value as at the time of declaration.
                        distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx_init ;
                }else{ //vnfIDdst != SFCdst
                    int vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx); int vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx); int pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
                    VNFNode<type_res> *dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
                    if(pn_src_idx == pn_dst_idx){ T_tx=0; T_px=0;
                    }else{ T_tx = T_tx_init; T_px = calcTime_PropagationDelay<type_wgt, type_res>(pn_src_idx, pn_dst_idx, PhysicalNetwork);
                    }
                    T_prc = calcTime_MeanProcessingDelayVNF<type_res>(dstVNFNode);
                    T_exe = calcTime_FunctionExecutionDelay<type_res>(dstVNFNode);
                    T_qd = calcTime_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), VNFNetwork->utilization, allSFC[givenSFC->index]);
//                    try{ T_qd=calcTime_QueuingDelay1(vnf_dst_idx, dstVNFNode->serviceRate, X_VNFType2Inst,givenSFC->index ,allSFC);
//                    } catch( std::exception const& e ) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
//                    T_tx = 1; T_px = 1;  T_prc = 1; T_exe = 1;  T_qd = 1;
                    if(showInConsole) {  cout<<" F["<<vnf_dst_idx<<char(96+vnf_dst_inst_idx)<<"],VM["<<vm_dst_idx<<"],PN["<<pn_dst_idx<<"]"<<"\t" <<"T_tx["<<T_tx<<"] | T_px["<<T_px<<"] | T_qd["<<T_qd<<"] | T_prc["<<T_prc<<"] | T_exe["<<T_exe<<"]";}
                    /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
                    if(distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe)
                        distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe;
                }  //vnfIDdst != SFCdst
            } // if distance != min
    }//for each vnf in topological order.

    if(showInConsole) { cout<<endl<<"Distance:\n";  for(const auto& d: givenSFC->vnfSeq){ cout<<"f["<<d<<"]: "<<distance[d]<<" |  "; }}

    if(debug and showInConsole)cout<<"\n<[Function Completed: "<<__FUNCTION__<<"]";

    if ( distance[SFCdst] <= 0) { // if distance obtanied is 0 or negative (unreachable or -1 default min distance)
        return std::numeric_limits<type_delay>::max();
    }
    return distance[SFCdst];
}//calcObjectiveValueSeq


#endif //SFC_PARALLELIZATION_TIMECALCULATIONFUNCTIONS_H

