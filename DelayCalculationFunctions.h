//
// Created by vijay on 24-03-2023.
//

#ifndef SFC_PARALLELIZATION_DELAYCALCULATIONFUNCTIONS_H
#define SFC_PARALLELIZATION_DELAYCALCULATIONFUNCTIONS_H


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
type_delay calcD_TransmissionDelay([[maybe_unused]] unsigned int n1=0, [[maybe_unused]] unsigned int n2=0,
                                 [[maybe_unused]] ServiceFunctionChain *const SFC = nullptr){
    auto bandwidth_n1_n2 = (type_delay)bandwidthNW * factor_bandwidth; // in Mb/s
    auto packetsize_n1_n2 = ((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet; // in bits
    return packetsize_n1_n2/(bandwidth_n1_n2)*timesfactor_tx;
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
template<typename type_wgt, typename type_res>
type_delay calcD_PropagationDelay(unsigned int n1, unsigned int n2, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork){
    auto len_n1_n2 = (type_delay)PhysicalNetwork->dist[n1][n2];
    auto speed = velocityFactor * (type_delay)speedOfLight;
    return (len_n1_n2 / (speed))*timesfactor_px*timesfactor_px*timesfactor_px;
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
template<typename type_res>
type_delay calcD_QueuingDelay(const VNFNode<type_res> *const fnNode, const unsigned int& fnInstId,
                            const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization, const ServiceFunctionChain*const cSFC){
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
    return (lambda / (mu_f * (mu_f - lambda)))*timesfactor_qd;
}

/*! Return the time required to complete the execution of the function VNF f. \n
 * @param VNFNode for function execution time.
 * @return execution time in seconds.
 */
template<typename type_res>
type_delay calcD_FunctionExecutionDelay(VNFNode<type_res> *const VNFNode){
    auto exeTime = (type_delay)VNFNode->executionTime;
    return exeTime*timesfactor_fnExe;
}

/*! Calculate the mean processing delay of particular VNF, time it takes function to process the one packet. \n
 * Fromula = 1/serviceRate of VNF.
 * @param VNFNode VNF Node class for service rate
 * @return mean processing delay in seconds.
 */
template<typename type_res>
type_delay calcD_MeanProcessingDelayVNF(VNFNode<type_res>*const VNFNode){
    auto mu_f = (type_delay)VNFNode->serviceRate;
    return (1/mu_f)*timesfactor_pkt;
}

/* --------------------------- --------------------------- --------------------------- ------------------------------------------- */
/*! The time required by a server v to duplicate the whole packet for instances
 * installed in d different servers of the next step.\n
 * Formula = T_d_pkt * (d-1);
 * @return INTER duplication time in seconds.
 * For example; packetBodySize = 1000 B; packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_pkt=0.630784(d=2),  =1.26157(d=3),  =1.892352(d=4)
 */
type_delay calcD_InterDuplicationTime(unsigned int cntNextHopDiffServer){
    if( cntNextHopDiffServer <= 1) return 0;
    auto packetsize_n1_n2 = ((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet;
    type_delay time_to_duplicate_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_duplicate_packet*(type_delay)(cntNextHopDiffServer-1)*timesfactor_pkt;
}
/*! the time required by a server v to duplicate only the header for its p parallel instances. \n
 * Formula = T_d_hdr * (p-1);
 * @return INTRA duplication time in seconds.
 * For example;packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_hdr=0.014784(p=2),  =0.029568(p=3),  =0.044352(p=4)
 */
type_delay calcD_IntraDuplicationTime(unsigned int cntParallelServer){
    if(cntParallelServer <= 1) return 0;
    type_delay time_to_duplicate_header = (type_delay)packetHeaderSize *factor_packet * read_write_time_per_bit;
    return time_to_duplicate_header*(type_delay)(cntParallelServer-1)*timesfactor_pkt;
}
/*! the time required by a server v to merge the whole packets from d different servers of the previous step, \n
 * Formula = T_m_pkt * (d-1);
 * @return INTER merging time in seconds.
 */
type_delay calcD_InterMergingTime(unsigned int cntPrevHopDiffServer){
    if(cntPrevHopDiffServer <= 1) return 0;
    auto packetsize_n1_n2 = ((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet;
    type_delay time_to_merging_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_merging_packet*(type_delay)(cntPrevHopDiffServer-1)*timesfactor_pkt;
}
/*! the time required by a server v to merge the headers from its p parallel instances \n
 * Formula = T_m_hdr * (p-1);
 * @return INTRA merging time in seconds.
 */
type_delay calcD_IntraMergingTime(unsigned int cntParallelServer){
    if(cntParallelServer <= 1) return 0;
    type_delay time_to_merge_header = (type_delay)packetHeaderSize *factor_packet * read_write_time_per_bit;
    return time_to_merge_header*(type_delay)(cntParallelServer-1)*timesfactor_pkt;
}

/* --------------------------- --------------------------- --------------------------- ------------------------------------------- */

/*! Calculate the end-to-end delay of given PARALLEL SFC with VNF instances assigned. \n
 * End-to-End delay includes - transmission time, propagation time, execution time, queuing delay, processing time and packet processing time. \n
 * Packet processing time includes inter duplication, intra-duplication, inter-merging and intra merging.
 * @param[in] givenParSFC given enumeration of vnfBlockPar of cSFC for which we have to calcualte the delay (except src dst block)
 * @param[in] X_VNFType2Instv In curSFC VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
 * @param[in] VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
 * @param[in] VirtualNetwork VirtualMachines object to find the required Physical Node of the given VM id
 * @param[in] PhysicalNetwork PhysicalGraph object to find the propagation delay, min dist between two node
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @tparam type_wgt type_wgt edge weight data type.
 * @tparam type_res resource data type.
 * @return time (in seconds) End to end delay of SFC. If not reachable (<=0) then max time(type_delay value).
 * src is dummy state, next stage to process is first parallel block  i.e index 0 of givenParSFC\n
 *                      [src] -> [f1] -> [f6 f4] -> [f5] -> [dst]\n
 * arr index(stg_nxt-1)           0        1         2     szStages=3 (excluding src and dst)
 */
template<typename type_wgt, typename type_res>
type_delay calcD_ParallelSFC(const vector<vector<unsigned int>>& givenParSFC, const unordered_map<unsigned int,unsigned int>& X_VNFType2Inst,
                             const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,  ServiceFunctionChain*const cSFC,
                                 const VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
                                 bool showInConsole=false, bool showInConsoleDetailed=false) {//calcD_ParallelSFC

    unsigned int szStages = givenParSFC.size();
    type_delay  T_tx_init = calcD_TransmissionDelay(); ///< transmission time
    vector<unordered_map<unsigned int, unsigned int>> cntPN(szStages); ///< all the DIFFERENT physical nodes {nodeid, freq} belonging to each stages.
    vector<type_delay> distance(szStages, -1); ///< final dist at each node.

    unsigned int id_curstg=0;
/* ****************************************************************************************************** */
    {//From Dummy SRC to First Stage(index 0).
        unordered_map<unsigned int, type_delay> exePN; ///< maximum execution time of physical node in the instance combination
        if (showInConsoleDetailed) {cout << "\nSTG:" << id_curstg;}

        for (const unsigned int &vnf_dst_idx: givenParSFC[id_curstg]) { ///< processing the current stage.
            const unsigned int& vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);  const unsigned int& vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);
            const unsigned int& pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);

            VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
            type_delay T_prc = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
            type_delay T_exe = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
            type_delay  T_qd = calcD_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), oldUtilization,  cSFC);

            exePN[pn_dst_idx] = max(exePN[pn_dst_idx], T_prc + T_exe + T_qd);  /// maximum execution according to all the instnaces in the server
            cntPN[id_curstg][pn_dst_idx] += 1; ///< count of physical nodes, in current stage
            if (showInConsoleDetailed) {
                cout << "\n   F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "]VM[" << vm_dst_idx<< "]PN[" << pn_dst_idx << "]"
                     << "  [qd: " << T_qd << "]  [prc: " << T_prc << "]  [exe: " << T_exe<< "]" << "  sum: "<<T_prc + T_exe + T_qd;
            }
        }///< processing the current stage.

        const unordered_map<unsigned int, unsigned int>& curStgPN = cntPN[id_curstg];
        type_delay mx_delay_x_y = 0;///< maximum delay of the previous and current stage (x,y).
        for(const auto &[pn_y, pn_y_parallelcnt]: curStgPN){ /*! for each physical server in cur node*/
            /*! Calculation of packet processing time. intra duplication (within same server multiple nodes), intra Merging (within same server multiple nodes)*/
            type_delay T_d_hdr = calcD_IntraDuplicationTime(pn_y_parallelcnt);
            type_delay T_m_hdr = calcD_IntraMergingTime(pn_y_parallelcnt);
            mx_delay_x_y = max(mx_delay_x_y, T_d_hdr+T_m_hdr + exePN[pn_y]);
            if(showInConsoleDetailed) {
                cout<<"\n       :"<<pn_y_parallelcnt<<"[py:"<<pn_y<<"]"<<"  [d_hdr:"<<T_d_hdr<<" m_hdr:"<<T_m_hdr<<"]"<<"   mxServer:"<<exePN.at(pn_y);
            }
        }//curStgPN
        type_delay T_d_pkt =  calcD_InterDuplicationTime(curStgPN.size());
        mx_delay_x_y += T_d_pkt + T_tx_init;
        distance[id_curstg] = mx_delay_x_y;
        if(showInConsoleDetailed) {
            cout<<"\n     px-py:"<< distance[id_curstg]<<"   d_pkt:"<<T_d_pkt;;
        }
    }//From Dummy SRC to First Stage(index 0).
/* ****************************************************************************************************** */
    for( id_curstg=1; id_curstg<szStages; id_curstg++) {/*!< Iterating stage ID from 0 to last index */
        unsigned int id_prvstg = id_curstg-1; ///< previous stage id to process
        if (showInConsoleDetailed) {cout << "\nSTG:" << id_curstg;}
        /*! Processing of current stage: count of physical servers, max time in each server */
            unordered_map<unsigned int, type_delay> exePN; ///< maximum execution time of physical node in the instance combination
            for (const unsigned int &vnf_dst_idx: givenParSFC[id_curstg]) { ///< processing the current stage.
                const unsigned int& vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);  const unsigned int& vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);
                const unsigned int& pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);

                VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
                type_delay T_prc = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
                type_delay T_exe = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
                type_delay T_qd = calcD_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), oldUtilization,  cSFC);
                exePN[pn_dst_idx] = max(exePN[pn_dst_idx], T_prc + T_exe + T_qd);  /// maximum execution according to all the instnaces in the server
                cntPN[id_curstg][pn_dst_idx] += 1; ///< count of physical nodes, in current stage
                if (showInConsoleDetailed) {
                    cout << "\n   F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "]VM[" << vm_dst_idx<< "]PN[" << pn_dst_idx << "]"
                         << "  [qd: " << T_qd << "]  [prc: " << T_prc << "]  [exe: " << T_exe<< "]" << "  sum: "<<T_prc + T_exe + T_qd;
                }
            }///< processing the current stage.

        /*! Once x = prvstage, y = curstage are fixed. We will find maximum time for each physical server in y to determine edge (x,y). */
                const unordered_map<unsigned int, unsigned int>& curStgPN = cntPN[id_curstg];
                const unordered_map<unsigned int, unsigned int>& prvStgPN = cntPN[id_prvstg];

        type_delay mx_delay_x_y = 0;///< maximum delay of the previous and current stage (x,y).
        for(const auto &[pn_y, pn_y_parallelcnt]: curStgPN){ /*! for each physical server in cur node*/
            /*! Calculation of packet processing time. inter mergring ( different server from cur server), intra duplication (within same server multiple nodes), intra Merging (within same server multiple nodes)*/
                unsigned int py_px_same=0; if(prvStgPN.find(pn_y) != prvStgPN.end()) py_px_same =  1;
                unsigned int cntPrevHopDiffServer = prvStgPN.size() - py_px_same;
                type_delay T_m_pkt = calcD_InterMergingTime(cntPrevHopDiffServer);

                type_delay T_d_hdr = calcD_IntraDuplicationTime(pn_y_parallelcnt);
                type_delay T_m_hdr = calcD_IntraMergingTime(pn_y_parallelcnt);

                type_delay mx_pktPrc =  T_m_pkt  + T_d_hdr + T_m_hdr;///< overall total time spent in packet processing from src to dest.

                if(showInConsoleDetailed) {
                    cout<<"\n     :"<<pn_y_parallelcnt<<"[py:"<<pn_y<<"]"<<"  prvD:"<<cntPrevHopDiffServer<<" s("<<py_px_same<<")"
                        <<"   [m_pkt:"<<T_m_pkt<<" d_hdr:"<<T_d_hdr<<" m_hdr:"<<T_m_hdr<<"]"<<"   mxServer:"<<exePN.at(pn_y);
                }
                /*! Calculation of inter-duplication time transmission time and  propagation time (we can duplicate the packets right before sending them. */
                type_delay mx_interdupTxPx_for_y = 0;
                for (const auto &[pn_x, pn_x_parallelcnt]: prvStgPN) { /*! for each physical server in previous node*/
                    unsigned int px_py_same = 0;
                    if (curStgPN.find(pn_x) != curStgPN.end()) px_py_same = 1;
                    unsigned int cntNextHopDiffServer = curStgPN.size() - px_py_same;
                    type_delay T_d_pkt =  calcD_InterDuplicationTime(cntNextHopDiffServer);
                    type_delay T_tx=0, T_px=0;
                    if(pn_y != pn_x){ /// if both server are different then there is transmission and propagation delay
                        T_tx = T_tx_init; T_px =  calcD_PropagationDelay<type_wgt, type_res>(pn_x, pn_y,PhysicalNetwork);
                    }
                    mx_interdupTxPx_for_y = max(mx_interdupTxPx_for_y, T_d_pkt + T_tx + T_px);///< overall total time spent in sending packet from src to dest.
                    if (showInConsoleDetailed) {
                        cout << "\n         " << pn_x_parallelcnt << "(px:" << pn_x << ")  " << "nxtD:" << cntNextHopDiffServer << " s(" << px_py_same << ")"
                             << "   [d_pkt:" << T_d_pkt << " tx:" << T_tx << " px:" << T_px << "]";
                    }
                }/*! for each physical server in previous node*/
                mx_delay_x_y = max(mx_delay_x_y, mx_interdupTxPx_for_y + mx_pktPrc + exePN.at(pn_y));

                if(showInConsoleDetailed) {
                    cout<<"   mxTxPx:"<<mx_interdupTxPx_for_y;
                    cout<<"\n     px-py:"<< mx_interdupTxPx_for_y + mx_pktPrc + exePN.at(pn_y);
                }
            }/// for each physical server in cur node
        distance[id_curstg] = distance[id_prvstg] + mx_delay_x_y;
    } //id_curstg
/* ****************************************************************************************************** */
    {//From last Stage to dummy destination
        id_curstg = szStages-1;
        type_delay T_m_pkt = calcD_InterMergingTime(cntPN[id_curstg].size());
        type_delay mx_delay_x_y = T_tx_init + T_m_pkt;
        distance[id_curstg] = distance[id_curstg] + mx_delay_x_y;
    }//From last Stage to dummy destination
/* ****************************************************************************************************** */
    if(showInConsole) { cout<<endl<<"Delays(Stage wise): ";  for(int i=0; i<szStages; i++){ cout<<"["<<i<<": "<<distance[i]<<"s] | "; }}

    if ( distance[id_curstg] <= 0) { // if distance obtanied is 0 or negative (unreachable or -1 default min distance)
        return std::numeric_limits<type_delay>::max();
    }

    return distance[id_curstg];
}//calcD_ParallelSFC

/*! Calculate the end-to-end delay of given SEQUENTIAL SFC (SFC->vnfSeq) with VNF instances assigned. \n
 * End-to-End delay includes - transmission time, propagation time, execution time, queuing delay, processing time.
 * Don't include same vnf again in same SFC.
 * @param cSFC given SFC object whose delay we have to colculate. cSFC->vnfSeq (including src and dest) array used.
 * @param X_VNFType2Instv In curSFC VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param allSFC all the SFC in the network, required to calculate arrival rate for queuing delay.
 * @param VNFNetwork VirtualNetworkFunctions object to find the required VM of the given vnf id and instance id
 * @param VirtualNetwork VirtualMachines object to find the required Physical Node of the given VM id
 * @param PhysicalNetwork PhysicalGraph object to find the propagation delay, min dist between two node
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @tparam type_wgt type_wgt edge weight data type.
 * @tparam type_res resource data type.
 * @return time (in seconds) End to end delay of SFC. If not reachable (<=0) then max time(type_delay value).
 */
template<typename type_wgt, typename type_res>
type_delay calcD_SequentialSFC(const ServiceFunctionChain *const cSFC, const unordered_map<unsigned int, unsigned int>& X_VNFType2Inst,
                                 const VirtualNetworkFunctions<type_res> *const VNFNetwork,const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
                                 bool showInConsole=false, bool showInConsoleDetailed=false){

    unsigned int Len = cSFC->vnfSeq.size(); // length of the sfc, no src dest stage
    type_delay  T_tx_init = calcD_TransmissionDelay(); ///< transmission time

    vector<unsigned int> whatPN(Len); ///< physical nodes id belonging to each stages, since for sequential it would be only physical server at each stage
    vector<type_delay> distance(Len, -1); ///< final dist at each node.

    for(int id_curstg=0; id_curstg<Len; id_curstg++) {/*!< Iterating stage ID from 0 to last index */
        int id_prvstg = id_curstg-1; ///< previous stage id to process

        /*! Processing of current stage: count of physical servers, max time in each server */
        const auto& vnf_dst_idx = cSFC->vnfSeq[id_curstg];
        const unsigned int& vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);  const unsigned int& vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);
        const unsigned int& pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);

        VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
        type_delay T_prc = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
        type_delay T_exe = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
        type_delay T_qd = calcD_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), VNFNetwork->seq_utilization,  cSFC);
        whatPN[id_curstg] = pn_dst_idx; ///<  physical node, in current stage
        type_delay T_tx=0, T_px=0;
        if(id_prvstg>0 and whatPN[id_prvstg] != pn_dst_idx){ /// if both server are different then there is transmission and propagation delay
            T_tx = T_tx_init; T_px =  calcD_PropagationDelay<type_wgt, type_res>(whatPN[id_prvstg], pn_dst_idx,PhysicalNetwork);
        }
        distance[id_curstg] = (T_qd + T_prc + T_exe);
        if(id_curstg == 0){ ///first stage
            distance[id_curstg] += T_tx_init ; /// maximum execution according to all the instnaces in the server
        } else{
            distance[id_curstg] += (distance[id_prvstg] + T_tx + T_px);
        }

        if (showInConsoleDetailed) {
            cout << "\nSTG:" << id_curstg;
            cout << "\n   F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "]VM[" << vm_dst_idx<< "]PN[" << pn_dst_idx << "]"
                 << "  [qd: " << T_qd << "]  [prc: " << T_prc << "]  [exe: " << T_exe<< "]" << "  sum: "<<T_prc + T_exe + T_qd;
        }
    }//for each vnf in topological order.

    distance[Len-1] += T_tx_init;//for last stage.

    if(showInConsole) { cout<<endl<<"Delays(Stage wise): ";  for(int i=0; i<Len; i++){ cout<<"["<<i<<": "<<distance[i]<<"s] | "; }}

    if ( distance[Len-1] <= 0) { // if distance obtanied is 0 or negative (unreachable or -1 default min distance)
        return std::numeric_limits<type_delay>::max();
    }

    return distance[Len-1];
}//calcObjectiveValueSeq

/* ---------------------------------------------------------------------- */

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
 * @tparam type_wgt type_wgt edge weight data type.
 * @tparam type_res resource data type.
 * @return time (in seconds) End to end delay of SFC. If not reachable (<=0) then max time(type_delay value).
 */
//template<typename type_wgt, typename type_res>
//type_delay calcObjectiveValueSeq(const ServiceFunctionChain *const givenSFC, const unordered_map<unsigned int, unsigned int>& X_VNFType2Inst,
//                            const vector<ServiceFunctionChain*>& allSFC, const VirtualNetworkFunctions<type_res> *const VNFNetwork,
//                            const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
//                            bool showInConsole=false){
//
//    if(debug and showInConsole)cout<<"\n>[Function Running: "<<__FUNCTION__<<"]";
////    if(showInConsole) { cout<<"\nTopological Order: ";for (int x: givenSFC->vnfSeq) { if(x == SFCsrc) cout<<"SRC -> "; else if(x == SFCdst) cout<<" DST"; else cout << x <<char(96+X_VNFType2Inst[x])<< " --> "; } }  ///< Print topological order
//
//    unsigned int sz = givenSFC->vnfSeq.size(); // src + dest = (2)
//
//    type_delay default_min_distance = -1; ///< default min distance to be used in distance(time) calculation.
//    unordered_map< int, type_delay> distance; ///<stores longest distance to each node
//    for(const auto& vnfid: givenSFC->vnfSeq)   distance[vnfid] = default_min_distance; ///< initialization of distance
//    distance[SFCsrc] = 0;
//
//    type_delay T_tx, T_tx_init = calcD_TransmissionDelay();
//    type_delay T_px, T_qd, T_prc, T_exe;
//
///*! Now we will calculate the longest path in the SFC graph. */
//     // as it is a sequential givenSFC only one node is there after SFCsrc. level(0) is SFCsrc
//
//        const auto &vnf_src_idx = givenSFC->vnfSeq[1]; //first parallel block with index = 1, index 0 is src
//        if (vnf_src_idx == SFCdst) {
//            string errorMsg = "Error in accessing SFC. vnf_dst_idx: " + to_string(vnf_src_idx) + ". Function:";
//            throw runtime_error(errorMsg + __FUNCTION__);
//        }
//        VNFNode<type_res> *const srcVNFNode = VNFNetwork->VNFNodes.at(vnf_src_idx);
//        T_prc = calcD_MeanProcessingDelayVNF<type_res>(srcVNFNode);
//        T_exe = calcD_FunctionExecutionDelay<type_res>(srcVNFNode);
//
//        T_qd = calcD_QueuingDelay<type_res>(srcVNFNode, X_VNFType2Inst.at(vnf_src_idx), VNFNetwork->utilization, allSFC[givenSFC->index]);
////        try { T_qd = calcD_QueuingDelay1(vnf_src_idx, srcVNFNode->serviceRate, X_VNFType2Inst,givenSFC->index, allSFC);
////        } catch (std::exception const &e) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
////        T_prc=0.5; T_exe = 1;  T_qd=2;
//        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
//        distance[vnf_src_idx] = distance[SFCsrc] + T_tx_init + T_prc + T_exe + T_qd;
//
//
//    /// from 1 as SFCsrc(0) already considered and first parallel block index=1, TO as sz-2 as src -> sz-1 is destination SFCdst
//    for(unsigned int idx = 1; idx <= sz-2; idx++) {
//        const auto& vnf_src_idx = givenSFC->vnfSeq[idx]; // int qki -10 bhi hai -10 as dst
//        const auto& vnf_src_inst_idx = X_VNFType2Inst.at(vnf_src_idx);
//        const auto& vm_src_idx = VNFNetwork->I_VNFinst2VM.at(vnf_src_idx).at(vnf_src_inst_idx);
//        const auto& pn_src_idx = VirtualNetwork->I_VM2PN.at(vm_src_idx);
//
//            if(showInConsole) { cout<<"\n:src: F["<<vnf_src_idx<<char(96+vnf_src_inst_idx)<<"],VM["<<vm_src_idx<<"],PN["<<pn_src_idx<<"]"<<" :dst:";}
//            if(distance[vnf_src_idx] != default_min_distance)
//
//                const int& vnf_dst_idx= givenSFC->vnfSeq[idx+1];
//                if(vnf_dst_idx == SFCdst){ //    T_prc=0.5; //T_tx using same value as at the time of declaration.
//                        distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx_init ;
//                }else{ //vnfIDdst != SFCdst
//                    int vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx); int vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx); int pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
//                    VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
//                    if(pn_src_idx == pn_dst_idx){ T_tx=0; T_px=0;
//                    }else{ T_tx = T_tx_init; T_px = calcD_PropagationDelay<type_wgt, type_res>(pn_src_idx, pn_dst_idx, PhysicalNetwork);
//                    }
//                    T_prc = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
//                    T_exe = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
//                    T_qd = calcD_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), VNFNetwork->utilization, allSFC[givenSFC->index]);
////                    try{ T_qd=calcD_QueuingDelay1(vnf_dst_idx, dstVNFNode->serviceRate, X_VNFType2Inst,givenSFC->index ,allSFC);
////                    } catch( std::exception const& e ) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
////                    T_tx = 1; T_px = 1;  T_prc = 1; T_exe = 1;  T_qd = 1;
//                    if(showInConsole) {  cout<<" F["<<vnf_dst_idx<<char(96+vnf_dst_inst_idx)<<"],VM["<<vm_dst_idx<<"],PN["<<pn_dst_idx<<"]"<<"\t" <<"T_tx["<<T_tx<<"] | T_px["<<T_px<<"] | T_qd["<<T_qd<<"] | T_prc["<<T_prc<<"] | T_exe["<<T_exe<<"]";}
//                    /// T_tx using same value as at the time of declaration.  from source to all child, add processing delay of packet duplication.
//                    if(distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe)
//                        distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe;
//                }  //vnfIDdst != SFCdst
//              // if distance != min
//    }//for each vnf in topological order.
//
//    if(showInConsole) { cout<<endl<<"Distance:\n";  for(const auto& d: givenSFC->vnfSeq){ cout<<"f["<<d<<"]: "<<distance[d]<<" |  "; }}
//
//    if(debug and showInConsole)cout<<"\n<[Function Completed: "<<__FUNCTION__<<"]";
//
//    if ( distance[SFCdst] <= 0) { // if distance obtanied is 0 or negative (unreachable or -1 default min distance)
//        return std::numeric_limits<type_delay>::max();
//    }
//    return distance[SFCdst];
//}//calcObjectiveValueSeq

//template<typename type_res>
//type_delay calcD_PacketsDelay(const vector<vector<unsigned int>>& givenParSFC, const unordered_map<unsigned int, unsigned int>& X_VNFType2Inst,
//                              const VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork, bool showInConsole= false){
////    givenParSFC = { {3, 3, 4}, {3,4,4,5,5,6}, {3,3,4,7,6} };
//    if(debug and showInConsole) cout<<"\n>>>[Function Running: "<<__FUNCTION__<<"]";
//    unsigned int szStages = givenParSFC.size();
//    /*! cSFCsrc is dummy state, next stage to process is first parallel block   index 0 of parBlockArray\n
//     *                      [src] -> [f1] -> [f6 f4] -> [f5] -> [dst]\n
//     * arr index(stg_nxt-1)           0        1         2     szStages=3 (excluding src and dst)\n
//     *          stg_nxt      0        1        2         3       4
//     *                    stg_prv  stg_cur    stg_nxt
//     */
//    unsigned int stg_nxt=1;
//
//    vector<unordered_map<unsigned int, unsigned int>> stgPN(szStages+2); ///< all the DIFFERENT physical nodes {nodeid, freq} belonging to each   (+2 for src and dst stage).
//    type_delay totalPktDelay; // total max delay of SFC from src to dst incurred in a process to duplicate and merge packets.
////    unordered_map<unsigned int, unordered_map<int, type_delay>> pktdist; ///<stores paket delay for each PHYSICAL NODE {stg -> {pnId, time}}
//
//    /// SFC src --> STG[1] This for-loop is Only For Source Level to Next Level(index 0 of parallel block), only duplication required.
//    for (const int &vnf_dst_idx: givenParSFC[stg_nxt-1]){  //int vnf_dst_inst_idx =0, vm_dst_idx=0,pn_dst_idx =vnf_dst_idx;
//        const auto & vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx); const auto & vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);const auto & pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
//        stgPN[stg_nxt][pn_dst_idx] += 1;
//    }
//    /// type_delay T_d_pkt = Time to duplicate packets to d diff server from level_i to level_i+1, Here only duplication time required.
//    totalPktDelay = calcD_InterDuplicationTime(stgPN[stg_nxt].size()); //totalPktDelay = pktdist[stg_nxt-1][SFCsrc] = calcD_InterDuplicationTime(stgPN[stg_nxt].size());
//
//    if(showInConsole) {
//        cout<<"\n:STG[0]:-->:STG[1]: "<<stgPN[stg_nxt].size()<<"("; for(const auto& pn: stgPN[stg_nxt]) cout<<" PN["<<pn.first<<":"<<pn.second<<"]";  cout<<" )\t[d_pkt:"<<totalPktDelay<<"]";
//    }
//
//    /*!  From STG = 2nd Parallel stage  to szStages-1 2nd Parallel stage,
//     *   each loop first calculates the information of next stage and then process the current stage.
//     *   At last, second last stage is info is calculated but not processed.
//     */
//    for(  stg_nxt=2; stg_nxt <= szStages+1; stg_nxt++) {
//        unsigned int stg_cur = stg_nxt-1;
//        unsigned int stg_prv = stg_cur-1;
//
//        // processing the next level stgi+1, except last level
//        if(stg_nxt != szStages+1){
//            for (const auto &vnf_dst_idx: givenParSFC[stg_nxt-1]) { //int vnf_dst_inst_idx =0,vm_dst_idx=0,pn_dst_idx =vnf_dst_idx;
//                const auto & vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);
//                const auto & vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);
//                const auto & pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
//                stgPN[stg_nxt][pn_dst_idx] += 1;
//            }}
//
//        if(showInConsole) {
//            cout<<"\n:STG["<<stg_cur<<"]:-->:STG["<<stg_nxt<<"]: "<<stgPN[stg_nxt].size()<<"("; for(const auto& pn: stgPN[stg_nxt]) cout<<" PN["<<pn.first<<":"<<pn.second<<"]"; cout<<" )";
//        }
//
//        unsigned int nxtStgSz = stgPN[stg_nxt].size(), prvStgSz = stgPN[stg_prv].size(); ///< size of stages i.e. number of different servers.
//        type_delay stgMaxDelay = 0;
//        // processing the current stage for merging and duplication time
//        for(const auto &[pn_idx, cntParallelServer]: stgPN[stg_cur]){
//
//            unsigned int curNxtSamePN=0; if(stgPN[stg_nxt].find(pn_idx) != stgPN[stg_nxt].end()) curNxtSamePN =  1;
//            unsigned int cntNextHopDiffServer = nxtStgSz - curNxtSamePN;
//            type_delay T_d_pkt = calcD_InterDuplicationTime(cntNextHopDiffServer);
//            type_delay T_d_hdr = calcD_IntraDuplicationTime(cntParallelServer);
//
//            unsigned int curPrvSamePN=0; if(stgPN[stg_prv].find(pn_idx) != stgPN[stg_prv].end()) curPrvSamePN =  1;
//            unsigned int cntPrevHopDiffServer = prvStgSz - curPrvSamePN;
//            type_delay T_m_pkt = calcD_InterMergingTime(cntPrevHopDiffServer);
//            type_delay T_m_hdr = calcD_IntraMergingTime(cntParallelServer);
//
//            //pktdist[stg_nxt][pn_idx] = T_m_pkt + T_m_hdr + T_d_hdr + T_d_pkt;
//            stgMaxDelay = max(stgMaxDelay,  T_m_pkt + T_m_hdr + T_d_hdr + T_d_pkt);///< overall total time spent in packet processing from src to dest.
//            if(showInConsole) {
//                cout<<"\n    "<<cntParallelServer<<"*( PN:"<<pn_idx<<" ) --> "<<"[prvD:"<<cntPrevHopDiffServer<<" s("<<curPrvSamePN<<") ; nxtD:"<<cntNextHopDiffServer<<" s("<<curNxtSamePN<<")]";
//                cout<<"\t[m_pkt:"<<T_m_pkt<<" ; m_hdr:"<<T_m_hdr<<"]"; cout<<" [d_pkt:"<<T_d_pkt<<" ; d_hdr:"<<T_d_hdr<<"]";
//            }
//        }
//        totalPktDelay += stgMaxDelay;
//    }
//
////     Processing the last stage. only merging required.
//    totalPktDelay += calcD_InterMergingTime(stgPN[szStages].size()); //<type_delay T_m_pkt // pktdist[szStages+1][SFCdst] =
//
//    if(showInConsole) {
//        cout<<"\n:STG[DST]: "<<"prvD:"<<stgPN[szStages].size(); cout<<" \t[E2E:"<<totalPktDelay<<"]";
//    }
//    if(debug and showInConsole)cout<<"\n<<<[Function Completed: "<<__FUNCTION__<<"]";
//
//    if ( totalPktDelay < 0) { // if distance obtanied is negative (unreachable or -1 default min distance). delay can be 0 in case of no packet processing.
//        return std::numeric_limits<type_delay>::max();
//    }
//    return totalPktDelay;
//}//calcD_PacketsDelay
//

//template<typename type_wgt, typename type_res>
//type_delay calcObjectiveValuePar(const vector<vector<unsigned int>>& givenParSFC, const unordered_map<unsigned int,unsigned int>& X_VNFType2Inst,
//                                 const unsigned int &curSFCidx, const vector<ServiceFunctionChain*>& allSFC,
//                                 const VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork,
//                                 const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork, bool showInConsole=false, bool showInConsoleDetailed=false){
//    if(debug and showInConsoleDetailed)cout<<"\n>>[Function Running: "<<__FUNCTION__<<"]";
//
//    type_delay default_min_distance = -1; ///< default min distance to be used in distance(time) calculation.
//    unordered_map<int, type_delay> distance; ///<stores longest distance to each node
//    for(const auto& vnfid: allSFC[curSFCidx]->vnfSeq)   distance[vnfid] = default_min_distance; // jitne bhi function hai unko min dist se intialise kr do
//    distance[SFCsrc] = 0;
//
//    type_delay T_tx, T_tx_init = calcD_TransmissionDelay();
//    type_delay T_px, T_qd, T_prc, T_exe;
//
//    unsigned int szStages = givenParSFC.size();
//    /*! cSFCsrc is dummy state, next stage to process is first parallel block    index 0 of parBlockArray\n
//     *                      [src] -> [f1] -> [f6 f4] -> [f5] -> [dst]\n
//     * arr index(stg_nxt-1)           0        1         2     szStages=3 (excluding src and dst)\n
//     *          stg_nxt      0        1        2         3       4
//     *                    stg_prv  stg_cur    stg_nxt
//     */
//    unsigned int stg_nxt=1;
//
///*! Now we will calculate the longest path in the SFC graph. */
//    for (const auto &vnf_dst_idx: givenParSFC[stg_nxt-1]){ // process first parallel block, vnfsrc is SFCsrc
//        if (vnf_dst_idx == SFCdst) {
//            string errorMsg = "Error in accessing SFC. vnf_dst_idx: " + to_string(vnf_dst_idx) + ". Function:";
//            throw runtime_error(errorMsg + __FUNCTION__);
//        }
//        VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
//        T_prc = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
//        T_exe = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
//        T_qd = calcD_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), VNFNetwork->utilization,  allSFC[curSFCidx]);
////        try { T_qd = calcD_QueuingDelay1(vnf_dst_idx, dstVNFNode->serviceRate, X_VNFType2Inst,curSFCidx,allSFC);
////        } catch (std::exception const &e) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
////          T_prc=0.5; T_exe = 1; T_qd=2;
//        /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
//        distance[vnf_dst_idx] = distance[SFCsrc] + T_tx_init + T_prc + T_exe + T_qd;
//    }
//
//    /*!  Stage Wise execution to maintain topological order
//     *  for each vnfsrc in cur stage we will process next stage vnfs except last stage (szStages-1).
//     */
//    for( stg_nxt=1; stg_nxt<=szStages; stg_nxt++){
//
//        for(const auto &vnf_src_idx: givenParSFC[stg_nxt-1]) { //for each vnf in src stage
//            const auto & vnf_src_inst_idx = X_VNFType2Inst.at(vnf_src_idx);
//            const auto & vm_src_idx = VNFNetwork->I_VNFinst2VM.at(vnf_src_idx).at(vnf_src_inst_idx);
//            const auto & pn_src_idx = VirtualNetwork->I_VM2PN.at(vm_src_idx);
//            if (showInConsoleDetailed) {
//                cout << "\n:src: F[" << vnf_src_idx << char(96 + vnf_src_inst_idx) << "],VM[" << vm_src_idx << "],PN["<< pn_src_idx << "]" << "  ";
//            }
//
//            if(stg_nxt == szStages){ // process last stage
//                if (distance[SFCdst] < distance[vnf_src_idx] + T_tx_init)
//                    distance[SFCdst] = distance[vnf_src_idx] + T_tx_init;
//            }  else if (distance[vnf_src_idx] != default_min_distance) { // else process intermediate nodes
//
//                for (const auto &vnf_dst_idx: givenParSFC[stg_nxt]) {  //stg_nxt + 1 - 1, stg_nxt next stage, -1 is the id
//                    const auto & vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);
//                    const auto & vm_dst_idx = VNFNetwork->I_VNFinst2VM.at(vnf_dst_idx).at(vnf_dst_inst_idx);
//                    const auto & pn_dst_idx = VirtualNetwork->I_VM2PN.at(vm_dst_idx);
//                    VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(vnf_dst_idx);
//
//                    if (pn_src_idx == pn_dst_idx) {  T_tx = 0; T_px = 0;
//                    } else {  T_tx = T_tx_init; T_px = calcD_PropagationDelay<type_wgt, type_res>(pn_src_idx, pn_dst_idx, PhysicalNetwork);
//                    }
//                    T_prc = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
//                    T_exe = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
//                    T_qd = calcD_QueuingDelay<type_res>(dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), VNFNetwork->utilization, allSFC[curSFCidx]);
////                        try {  T_qd = calcD_QueuingDelay1(vnf_dst_idx, dstVNFNode->serviceRate, X_VNFType2Inst,curSFCidx, allSFC);
////                        } catch (std::exception const &e) { T_qd = std::numeric_limits<type_delay>::max(); std::cerr << "\ncaught: " << e.what() ; }
////                    T_tx = 1; T_px = 1; T_qd = 1; T_prc = 1; T_exe = 1; T_prc=0.5; T_exe = 1; T_qd=2;
//                    if (showInConsoleDetailed) {
//                        cout << "\n\t :dst:F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "],VM[" << vm_dst_idx
//                             << "],PN[" << pn_dst_idx << "]" << "\t" << "T_tx[" << T_tx << "] | T_px[" << T_px
//                             << "] | T_qd[" << T_qd << "] | T_prc[" << T_prc << "] | T_exe[" << T_exe
//                             << "]" << "\n";
//                    }
//                    /// T_tx using same value as at the time of declaration.  from source to all child, add processing time of packet duplication.
//                    if (distance[vnf_dst_idx] < distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe)
//                        distance[vnf_dst_idx] = distance[vnf_src_idx] + T_tx + T_px + T_qd + T_prc + T_exe;
//
//                }   // for each child of source
//            } // if distance != min
//        } //for each src vnf
//    }///<Stage Wise execution to maintain topological order
//
//    if(showInConsole) { cout<<endl<<"Distance:\n";  for(const auto& d: allSFC[curSFCidx]->vnfSeq){ cout<<"f["<<d<<"]:"<<distance[d]<<" |  "; }}
//
//    if(debug and showInConsoleDetailed)cout<<"\n<<[Function Completed: "<<__FUNCTION__<<"]";
//
//    if ( distance[SFCdst] <= 0) { // if distance obtanied is 0 or negative (unreachable or -1 default min distance)
//        return std::numeric_limits<type_delay>::max();
//    }
//    return distance[SFCdst];
//}//calcObjectiveValuePar


#endif //SFC_PARALLELIZATION_DELAYCALCULATIONFUNCTIONS_H

