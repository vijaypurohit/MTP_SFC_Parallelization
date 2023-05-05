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
type_delay calcD_TransmissionDelay(){
    const auto& bandwidth_n1_n2 = (type_delay)bandwidthNW * factor_bandwidth; // in Mb/s
    const auto& packetsize_n1_n2 = (type_delay)(packetBodySize + packetHeaderSize)*factor_packet; // in bits
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
type_delay calcD_PropagationDelay(const unsigned int& n1, const unsigned int& n2, const PhysicalGraph& PhysicalNetwork){
    const auto& len_n1_n2 = (type_delay)PhysicalNetwork.dist[n1][n2];
    const auto& speed = velocityFactor * (type_delay)speedOfLight;
    return (len_n1_n2 / (speed));
}

/*! Calculate the Queuing delay of a particular VNF (time a job waits in a queue until it can be executed) \n
 * Formula = lambda_c / (mu_f * (mu_f - lambda_c)) \n
 * where lambda_c is traffic arrival rate of chain c or SFC s  and mu_f is the service rate of SFC or say VNF. \n
 * Both in packets per second. \n
 * @param cSFCArrivalRate SFC arrival rate
 * @param fnNode VNF Node object whose Queuing delay we have to calculate
 * @param fnInstId In curSFC, VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param oldUtilization old utilisation of the current function type and instance till now.
 * @return queuing delay in seconds.
 */
type_delay calcD_QueuingDelay( const type_delay& cSFCArrivalRate, const VNFNode& fnNode, const unsigned int& fnInstId,
                            const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization){
    const unsigned int& fnType =  fnNode.index;

    type_delay mu_f = fnNode.serviceRate; ///<service rate of VNF Node.

    type_delay lambda =  cSFCArrivalRate ; ///< current rate
    if(oldUtilization.count(fnType) and oldUtilization.at(fnType).count(fnInstId)) ///< old utilization till now of VNF
        lambda += oldUtilization.at(fnType).at(fnInstId);

    if(lambda >= mu_f){// otherwise queuing delay will be negative
        return lambda/mu_f ; ///< utne factor me delay maximise kr do
//        return std::numeric_limits<type_delay>::max();  // we want to punish this mapping such that with this mapping time is maximum
    }
//    else if(lambda == mu_f){// otherwise queuing delay will be infinite,but we can still maximise it to 100% utilization
//        return 1;  //
//    }
    return (lambda / (mu_f * (mu_f - lambda)));
}

/*! Return the time required to complete the execution of the function VNF f. \n
 * @param VNFNode for function execution time.
 * @return execution time in seconds.
 */
type_delay calcD_FunctionExecutionDelay(const VNFNode& VNFNode){
    const auto& exeTime = (type_delay)VNFNode.executionTime;
    return exeTime;
}

/*! Calculate the mean processing delay of particular VNF, time it takes function to process the one packet. \n
 * Fromula = 1/serviceRate of VNF.
 * @param VNFNode VNF Node class for service rate
 * @return mean processing delay in seconds.
 */
type_delay calcD_MeanProcessingDelayVNF(const VNFNode& VNFNode){
    const auto& mu_f = (type_delay)VNFNode.serviceRate;
    return (1/mu_f);
}

/* --------------------------- --------------------------- --------------------------- ------------------------------------------- */
/*! The time required by a server v to duplicate the whole packet for instances
 * installed in d different servers of the next step.\n
 * Formula = T_d_pkt * (d-1);
 * @return INTER duplication time in seconds.
 * For example; packetBodySize = 1000 B; packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_pkt=0.630784(d=2),  =1.26157(d=3),  =1.892352(d=4)
 */
type_delay calcD_InterDuplicationTime(const unsigned int& cntNextHopDiffServer){
    if( cntNextHopDiffServer <= 1) return 0;
    const auto& packetsize_n1_n2 = (type_delay)(packetBodySize + packetHeaderSize)*factor_packet;
    const type_delay& time_to_duplicate_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_duplicate_packet*(type_delay)(cntNextHopDiffServer-1);
}
/*! the time required by a server v to duplicate only the header for its p parallel instances. \n
 * Formula = T_d_hdr * (p-1);
 * @return INTRA duplication time in seconds.
 * For example;packetHeaderSize = 24 B; read_write_time_per_bit = 0.077e-3;\n
 * T_d_hdr=0.014784(p=2),  =0.029568(p=3),  =0.044352(p=4)
 */
type_delay calcD_IntraDuplicationTime(const unsigned int& cntParallelServerForHDR, const unsigned int& cntParallelServerForPkt){
    type_delay time_to_intra_duplicate = 0;
    if(cntParallelServerForHDR >= 1){
        type_delay time_to_duplicate_header = (type_delay)packetHeaderSize *factor_packet * read_write_time_per_bit;
        time_to_intra_duplicate +=  time_to_duplicate_header*(type_delay)(cntParallelServerForHDR);
    }
    if(cntParallelServerForPkt >= 1){
        type_delay time_to_duplicate_pkt = (type_delay)(packetBodySize + packetHeaderSize)*factor_packet * read_write_time_per_bit;
        time_to_intra_duplicate +=  time_to_duplicate_pkt*(type_delay)(cntParallelServerForPkt);
    }
    return time_to_intra_duplicate;
}
/*! the time required by a server v to merge the whole packets from d different servers of the previous step, \n
 * Formula = T_m_pkt * (d-1);
 * @return INTER merging time in seconds.
 */
type_delay calcD_InterMergingTime(const unsigned int& cntPrevHopDiffServer){
    if(cntPrevHopDiffServer <= 1) return 0;
    const auto& packetsize_n1_n2 = ((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet;
    const type_delay& time_to_merging_packet = packetsize_n1_n2 * read_write_time_per_bit;
    return time_to_merging_packet*(type_delay)(cntPrevHopDiffServer-1);
}
/*! the time required by a server v to merge the headers from its p parallel instances \n
 * Formula = T_m_hdr * (p-1);
 * @return INTRA merging time in seconds.
 */
type_delay calcD_IntraMergingTime(const unsigned int& cntParallelServerForHDR, const unsigned int& cntParallelServerForPkt){
    type_delay time_to_intra_merge = 0;
    if(cntParallelServerForHDR >= 1){
        type_delay time_to_merge_header = (type_delay)packetHeaderSize *factor_packet * read_write_time_per_bit;
        time_to_intra_merge +=  time_to_merge_header*(type_delay)(cntParallelServerForHDR);
    }
    if(cntParallelServerForPkt >= 1){
        type_delay time_to_merge_pkt = (type_delay)((type_delay)packetBodySize + (type_delay)packetHeaderSize)*factor_packet * read_write_time_per_bit;
        time_to_intra_merge +=  time_to_merge_pkt*(type_delay)(cntParallelServerForPkt);
    }
    return time_to_intra_merge;
}

/* --------------------------- --------------------------- --------------------------- ------------------------------------------- */

/*! Calculate the end-to-end delay of given PARALLEL SFC with given VNF instances assigned based on given utilization and sfc arrival rate. \n
 * End-to-End delay includes - transmission time, propagation time, execution time, queuing delay, processing time and packet processing time. \n
 * Packet processing time includes inter duplication, intra-duplication, inter-merging and intra merging.
 * @param[in] cSFC current sfc object
 * @param[in] givenParVersion given enumeration of parallel version of cSFC for which we have to calcualte the delay (except src dst block)
 * @param[in] X_VNFType2Instv In curSFC VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}.
 * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
 * @param[in] Simulate Simulations Class object which contains the necessary details for simulation
 * @param showInConsole to show the calculation debug info in console. Default False.
 * @return End to end delay of SFC. If not reachable (<=0) then max time(type_delay value).
 * src is dummy state, next stage to process is first parallel block  i.e index 0 of givenParVersion\n
 *                      [src] -> [f1] -> [f6 f4] -> [f5] -> [dst]\n
 * arr index(stg_nxt-1)           0        1         2     szStages=3 (excluding src and dst)
 */

type_delay calcD_ParallelSFC(const ServiceFunctionChain& cSFC, const vector<vector<unsigned int>>& givenParVersion,
                             const unordered_map<unsigned int,unsigned int>& X_VNFType2Inst,
                             const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,
                             const Simulations& Simulate,
                             bool showInConsole=false, bool showInConsoleDetailed=false) {//calcD_ParallelSFC

    unsigned int szStages = givenParVersion.size();
    type_delay  T_tx_init = calcD_TransmissionDelay(); ///< transmission time
    vector<unordered_map<unsigned int, vector<unsigned int>>> cntPN(szStages); ///< all the DIFFERENT physical nodes {nodeid, freq} belonging to each stages.
    vector<type_delay> distance(szStages, -1); ///< final dist at each node.

    unsigned int id_curstg=0;
/* ****************************************************************************************************** */
    {//From SFC SRC to First Stage(index 0).
        unordered_map<unsigned int, type_delay> exePN; ///< maximum execution time of physical node in the instance combination
        if (showInConsoleDetailed) {cout << "\nSTG:" << id_curstg;}

        for (const unsigned int &vnf_dst_idx: givenParVersion[id_curstg]) { ///< processing the current stage.
            const unsigned int& vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);
            const unsigned int& pn_dst_idx = Simulate.I_VNFINST_2_PN.at(vnf_dst_idx).at(vnf_dst_inst_idx);

            const VNFNode& dstVNFNode = Simulate.VNFNetwork.VNFNodes.at(vnf_dst_idx);
            type_delay T_prc = calcD_MeanProcessingDelayVNF(dstVNFNode);
            type_delay T_exe = calcD_FunctionExecutionDelay(dstVNFNode);
            type_delay  T_qd = calcD_QueuingDelay(cSFC.trafficArrivalRate, dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), oldUtilization);
            exePN[pn_dst_idx] = max(exePN[pn_dst_idx], T_prc + T_exe + T_qd);  /// maximum execution according to all the instnaces in the server
//            cntPN[id_curstg][pn_dst_idx] += 1; ///< count of physical nodes, in current stage
            cntPN[id_curstg][pn_dst_idx].emplace_back(vnf_dst_idx); ///< vnfs in physical nodes, in current stage
            if (showInConsoleDetailed) {
                cout << "\n   F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "]PN[" << pn_dst_idx << "]"
                     << "  [qd: " << T_qd << "]  [prc: " << T_prc << "]  [exe: " << T_exe<< "]" << "  sum: "<<T_prc + T_exe + T_qd;
            }
        }///< processing the current stage.

        const unordered_map<unsigned int, vector<unsigned int>>& curStgPN = cntPN[id_curstg];
        type_delay mx_delay_x_y = 0;///< maximum delay of the previous and current stage (x,y).
        for(const auto &[pn_y, pn_y_parallel_vnfs]: curStgPN){ /*! for each physical server in cur node*/
            /*! Calculation of packet processing time. intra duplication (within same server multiple nodes), intra Merging (within same server multiple nodes)*/
            type_delay T_d_hdr=0, T_m_hdr=0;
            if(pn_y_parallel_vnfs.size()>1){
                unsigned int needIntraHdrCopy=0, needIntraPktCopy=0;
                for(int j=1; j<pn_y_parallel_vnfs.size(); j++){
                    const auto vnfj = pn_y_parallel_vnfs[j]; bool plzcopy=false;
                    for(int i=0; i<j; i++){
                        if(Simulate.parallelPairs.at(pn_y_parallel_vnfs[i]).at(vnfj) == pktCopy){
                            plzcopy=true;
                            break;
                        }
                    }
                    if(plzcopy)needIntraPktCopy++;
                    else needIntraHdrCopy++;
                }
                T_d_hdr = calcD_IntraDuplicationTime(needIntraHdrCopy, needIntraPktCopy);
                T_m_hdr = calcD_IntraMergingTime(needIntraHdrCopy, needIntraPktCopy);
            }
//            type_delay mx_pktPrc =  T_d_hdr + T_m_hdr;
            const unsigned int& pn_x = cSFC.access_nodes.first; /// sfc source physical node (only one node at previous stage)
            unsigned int px_py_same = 0;
            if (curStgPN.find(pn_x) != curStgPN.end()) px_py_same = 1;
            unsigned int cntNextHopDiffServer = curStgPN.size() - px_py_same;
            type_delay T_d_pkt =  calcD_InterDuplicationTime(cntNextHopDiffServer);
            type_delay T_tx=0, T_px=0;
            if(pn_x != pn_y){ /// if both server are different then there is transmission and propagation delay
                T_tx = T_tx_init; T_px =  calcD_PropagationDelay(pn_x, pn_y, Simulate.PhysicalNetwork);
            }

            mx_delay_x_y = max(mx_delay_x_y, (T_d_pkt + T_tx + T_px)+(T_d_hdr + T_m_hdr) + exePN[pn_y]);
            if(showInConsoleDetailed) {
                cout<<"\n       :"<<pn_y_parallel_vnfs.size()<<"[py:"<<pn_y<<"]"<< "   [d_pkt:" << T_d_pkt << " tx:" << T_tx << " px:" << T_px << "]"<<"  [d_hdr:"<<T_d_hdr<<" m_hdr:"<<T_m_hdr<<"]"<<"   mxServer:"<<exePN.at(pn_y);
            }
        }//curStgPN
        distance[id_curstg] = mx_delay_x_y;
        if(showInConsoleDetailed) {
            cout<<"\n     px-py:"<< distance[id_curstg];
        }
    }//From SFC SRC to First Stage(index 0).
/* ****************************************************************************************************** */
    for( id_curstg=1; id_curstg<szStages; id_curstg++) {/*!< Iterating stage ID from 0 to last index */
        unsigned int id_prvstg = id_curstg-1; ///< previous stage id to process
        if (showInConsoleDetailed) {cout << "\nSTG:" << id_curstg;}
        /*! Processing of current stage: count of physical servers, max time in each server */
            unordered_map<unsigned int, type_delay> exePN; ///< maximum execution time of physical node in the instance combination
            for (const unsigned int &vnf_dst_idx: givenParVersion[id_curstg]) { ///< processing the current stage.
                const unsigned int& vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);
                const unsigned int& pn_dst_idx = Simulate.I_VNFINST_2_PN.at(vnf_dst_idx).at(vnf_dst_inst_idx); 

                const VNFNode& dstVNFNode = Simulate.VNFNetwork.VNFNodes.at(vnf_dst_idx);
                type_delay T_prc = calcD_MeanProcessingDelayVNF(dstVNFNode);
                type_delay T_exe = calcD_FunctionExecutionDelay(dstVNFNode);
                type_delay T_qd = calcD_QueuingDelay(cSFC.trafficArrivalRate, dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), oldUtilization);
                exePN[pn_dst_idx] = max(exePN[pn_dst_idx], T_prc + T_exe + T_qd);  /// maximum execution according to all the instnaces in the server
//                cntPN[id_curstg][pn_dst_idx] += 1; ///< count of physical nodes, in current stage
                cntPN[id_curstg][pn_dst_idx].emplace_back(vnf_dst_idx); ///< count of physical nodes, in current stage
                if (showInConsoleDetailed) {
                    cout << "\n   F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "]PN[" << pn_dst_idx << "]"
                         << "  [qd: " << T_qd << "]  [prc: " << T_prc << "]  [exe: " << T_exe<< "]" << "  sum: "<<T_prc + T_exe + T_qd;
                }
            }///< processing the current stage.

        /*! Once x = prvstage, y = curstage are fixed. We will find maximum time for each physical server in y to determine edge (x,y). */
                const unordered_map<unsigned int, vector<unsigned int>>& curStgPN = cntPN[id_curstg];
                const unordered_map<unsigned int, vector<unsigned int>>& prvStgPN = cntPN[id_prvstg];

        type_delay mx_delay_x_y = 0;///< maximum delay of the previous and current stage (x,y).
        for(const auto &[pn_y, pn_y_parallel_vnfs]: curStgPN){ /*! for each physical server in cur node*/
            /*! Calculation of packet processing time. inter mergring ( different server from cur server), intra duplication (within same server multiple nodes), intra Merging (within same server multiple nodes)*/
                unsigned int py_px_same=0; if(prvStgPN.find(pn_y) != prvStgPN.end()) py_px_same =  1;
                unsigned int cntPrevHopDiffServer = prvStgPN.size() - py_px_same;
                type_delay T_m_pkt = calcD_InterMergingTime(cntPrevHopDiffServer);

                type_delay T_d_hdr=0, T_m_hdr=0;
                if(pn_y_parallel_vnfs.size()>1){
                    unsigned int needIntraHdrCopy=0, needIntraPktCopy=0;
                    for(int j=1; j<pn_y_parallel_vnfs.size(); j++){
                        const auto vnfj = pn_y_parallel_vnfs[j]; bool plzcopy=false;
                        for(int i=0; i<j; i++){
                            if(Simulate.parallelPairs.at(pn_y_parallel_vnfs[i]).at(vnfj) == pktCopy){
                                plzcopy=true; break;
                            }
                        }
                        if(plzcopy)needIntraPktCopy++;
                        else needIntraHdrCopy++;
                    }
                    T_d_hdr = calcD_IntraDuplicationTime(needIntraHdrCopy, needIntraPktCopy);
                    T_m_hdr = calcD_IntraMergingTime(needIntraHdrCopy, needIntraPktCopy);
                }

                type_delay mx_pktPrc =  T_m_pkt  + T_d_hdr + T_m_hdr;///< overall total time spent in packet processing from src to dest.

                if(showInConsoleDetailed) {
                    cout<<"\n     :"<<pn_y_parallel_vnfs.size()<<"[py:"<<pn_y<<"]"<<"  prvD:"<<cntPrevHopDiffServer<<" s("<<py_px_same<<")"
                        <<"   [m_pkt:"<<T_m_pkt<<" d_hdr:"<<T_d_hdr<<" m_hdr:"<<T_m_hdr<<"]"<<"   mxServer:"<<exePN.at(pn_y);
                }
                /*! Calculation of inter-duplication time transmission time and  propagation time (we can duplicate the packets right before sending them. */
                type_delay mx_interdupTxPx_for_y = 0;
                for (const auto &[pn_x, pn_x_parallel_vnfs]: prvStgPN) { /*! for each physical server in previous node*/
                    unsigned int px_py_same = 0;
                    if (curStgPN.find(pn_x) != curStgPN.end()) px_py_same = 1;
                    unsigned int cntNextHopDiffServer = curStgPN.size() - px_py_same;
                    type_delay T_d_pkt =  calcD_InterDuplicationTime(cntNextHopDiffServer);
                    type_delay T_tx=0, T_px=0;
                    if(pn_y != pn_x){ /// if both server are different then there is transmission and propagation delay
                        T_tx = T_tx_init; T_px =  calcD_PropagationDelay(pn_x, pn_y, Simulate.PhysicalNetwork);
                    }
                    mx_interdupTxPx_for_y = max(mx_interdupTxPx_for_y, T_d_pkt + T_tx + T_px);///< overall total time spent in sending packet from src to dest.
                    if (showInConsoleDetailed) {
                        cout << "\n         " << pn_x_parallel_vnfs.size() << "(px:" << pn_x << ")  " << "nxtD:" << cntNextHopDiffServer << " s(" << px_py_same << ")"
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
    {//From last Stage to destination of sfc
        id_curstg = szStages-1;
        const unsigned int& pn_y = cSFC.access_nodes.second; /// sfc destination physical node (only one node at last stage)
        const unordered_map<unsigned int, vector<unsigned int>>& prvStgPN = cntPN[id_curstg];
        type_delay mx_interdupTxPx_for_y=0;
        for (const auto &[pn_x, pn_x_parallel_vnfs]: prvStgPN) { /*! for each physical server in previous node*/
            type_delay T_tx=0, T_px=0;
            if(pn_x != pn_y){ /// if both server are different then there is transmission and propagation delay
                T_tx = T_tx_init; T_px =  calcD_PropagationDelay(pn_x, pn_y, Simulate.PhysicalNetwork);
            }
            mx_interdupTxPx_for_y = max(mx_interdupTxPx_for_y, T_tx + T_px);///< overall total time spent in sending packet from src to dest.
            if (showInConsoleDetailed) {
                cout << "\n         " << pn_x_parallel_vnfs.size() << "(px:" << pn_x << ")  "<< "   [tx:" << T_tx << " px:" << T_px << "]";
            }
        }/*! for each physical server in previous node*/
        unsigned int py_px_same=0; if(prvStgPN.find(pn_y) != prvStgPN.end()) py_px_same =  1;
        unsigned int cntPrevHopDiffServer = prvStgPN.size() - py_px_same;
        type_delay T_m_pkt = calcD_InterMergingTime(cntPrevHopDiffServer);

        distance[id_curstg] += mx_interdupTxPx_for_y+ T_m_pkt;

        if(showInConsoleDetailed) {
            cout<<"   mxTxPx:"<<mx_interdupTxPx_for_y;
            cout<<"\n     px-py:"<< mx_interdupTxPx_for_y +T_m_pkt;
        }
    }//From last Stage to destination of sfc
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
 * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
 * @param[in] Simulate Simulations Class object which contains the necessary details for simulation
 * @param showInConsole to show the calculation debug info in console. Default False.

 * @return End to end delay of SFC. If not reachable (<=0) then max time(type_delay value).
 */

type_delay calcD_SequentialSFC(const ServiceFunctionChain& cSFC, const unordered_map<unsigned int, unsigned int>& X_VNFType2Inst,
                               const unordered_map<unsigned int,  unordered_map<unsigned int, type_delay>>& oldUtilization,
                               const Simulations& Simulate,
                               bool showInConsole=false, bool showInConsoleDetailed=false){

    unsigned int Len = cSFC.vnfSeq.size(); // length of the sfc, no src dest stage
    type_delay  T_tx_init = calcD_TransmissionDelay(); ///< transmission time

    vector<unsigned int> whatPN(Len); ///< physical nodes id belonging to each stages, since for sequential it would be only physical server at each stage
    vector<type_delay> distance(Len, -1); ///< final dist at each node.

    for(int id_curstg=0; id_curstg<Len; id_curstg++) {/*!< Iterating stage ID from 0 to last index */
        int id_prvstg = id_curstg-1; ///< previous stage id to process

        /*! Processing of current stage: count of physical servers, max time in each server */
        const auto& vnf_dst_idx = cSFC.vnfSeq[id_curstg];
        const unsigned int& vnf_dst_inst_idx = X_VNFType2Inst.at(vnf_dst_idx);
        const unsigned int& pn_dst_idx = Simulate.I_VNFINST_2_PN.at(vnf_dst_idx).at(vnf_dst_inst_idx);

        const VNFNode& dstVNFNode = Simulate.VNFNetwork.VNFNodes.at(vnf_dst_idx);
        type_delay T_prc = calcD_MeanProcessingDelayVNF(dstVNFNode);
        type_delay T_exe = calcD_FunctionExecutionDelay(dstVNFNode);
        type_delay T_qd = calcD_QueuingDelay(cSFC.trafficArrivalRate,dstVNFNode, X_VNFType2Inst.at(vnf_dst_idx), oldUtilization);
        whatPN[id_curstg] = pn_dst_idx; ///<  physical node, in current stage

        type_delay T_tx=0, T_px=0;
        if(id_curstg == 0){
            if(cSFC.access_nodes.first != pn_dst_idx){
                T_tx = T_tx_init; T_px =  calcD_PropagationDelay(cSFC.access_nodes.first, pn_dst_idx,Simulate.PhysicalNetwork);
            }
            distance[id_curstg] =  (T_tx + T_px) + (T_qd + T_prc + T_exe);
        }
        else {
            if (whatPN[id_prvstg] != pn_dst_idx) { /// if both server are different then there is transmission and propagation delay
                T_tx = T_tx_init;
                T_px = calcD_PropagationDelay(whatPN[id_prvstg], pn_dst_idx, Simulate.PhysicalNetwork);
            }
            distance[id_curstg] = distance[id_prvstg] + (T_tx + T_px) + (T_qd + T_prc + T_exe);
        }
        if (showInConsoleDetailed) {
            cout << "\nSTG:" << id_curstg;
            cout << "\n   F[" << vnf_dst_idx << char(96 + vnf_dst_inst_idx) << "]PN[" << pn_dst_idx << "]"
                 << "  [qd: " << T_qd << "]  [prc: " << T_prc << "]  [exe: " << T_exe<< "]" << "  sum: "<<T_prc + T_exe + T_qd;
        }
    }//for each vnf in topological order.

    type_delay T_tx=0, T_px=0;
    if(whatPN[Len-1] != cSFC.access_nodes.second){
        T_tx = T_tx_init; T_px =  calcD_PropagationDelay(whatPN[Len-1], cSFC.access_nodes.second,Simulate.PhysicalNetwork);
    }
    distance[Len-1] += (T_tx+T_px);//for last stage.

    if(showInConsole) { cout<<endl<<"Delays(Stage wise): ";  for(int i=0; i<Len; i++){ cout<<"["<<i<<": "<<distance[i]<<"s] | "; }}

    if ( distance[Len-1] <= 0) { // if distance obtanied is 0 or negative (unreachable or -1 default min distance)
        return std::numeric_limits<type_delay>::max();
    }

    return distance[Len-1];
}//calcObjectiveValueSeq

/* ---------------------------------------------------------------------- */


#endif //SFC_PARALLELIZATION_DELAYCALCULATIONFUNCTIONS_H

