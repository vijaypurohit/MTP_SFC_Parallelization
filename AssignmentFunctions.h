//
// Created by vijay on 24-03-2023.
//

#ifndef SFC_PARALLELIZATION_ASSIGNMENTFUNCTIONS_H
#define SFC_PARALLELIZATION_ASSIGNMENTFUNCTIONS_H

/*!  Assign VM id to Physical Node id. \n
 * Also make a mapping of pnid -> list of vm ids.
 * @param VM_TO_PN contains mapping of vnf id,  to pn id.
 * @param VirtualNetwork VirtualMachines class
 * @param PhysicalNetwork PhysicalGraph class
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template<typename type_wgt, typename type_res>
void assign_VM_2_PN(const vector<pair<int,int>>& VM_TO_PN, VirtualMachines<type_res> *VirtualNetwork, PhysicalGraph<type_wgt, type_res> *PhysicalNetwork){
    for(const auto& [vmid, pnid] : VM_TO_PN){
        VirtualNetwork->I_VM2PN[vmid] = pnid;
        PhysicalNetwork->PNode[pnid]->pn2VMs.insert(vmid);
    }
}

/*!  Assign VNF instances to Virtual Machines. \n
 * Also make a mapping of VNF -> list of VNFs id.
 * @param VNF_TO_VM contains mapping of vnf id, instance id,  to vm id.
 * @param VNFNetwork VirtualNetworkFunctions class
 * @param VirtualNetwork VirtualMachines class
 * @tparam type_res resource data type.
 */
template<typename type_res>
void assign_VNF_2_VM(const vector<vector<unsigned int>>& VNF_TO_VM, VirtualNetworkFunctions<type_res> *VNFNetwork, VirtualMachines<type_res> *VirtualNetwork){

    for(const auto& assign : VNF_TO_VM){
        unsigned int fnType = assign[0], instance = assign[1], vmid = assign[2];
        VNFNetwork->I_VNFinst2VM[fnType][instance] = vmid;
        VirtualNetwork->VMNode[vmid]->vm2VNFs.emplace_back(fnType, instance);
    }
}

/*!  Assign total instances required for particular VNF. \n
 * Also make a mapping of VNF -> list of VNFs id.
 * @param VNF_TO_VM contains mapping of vnf id, instance id,  to vm id.
 * @param VirtualNetwork VirtualMachines class
 * @tparam type_res resource data type.
 */
template<typename type_res>
void assign_VNF_2_InstancesCnt(const vector<pair<unsigned int,unsigned int>>& VNF_TO_InstancesCnt, VirtualNetworkFunctions<type_res> *VNFNetwork){
    for(const auto& [fnType, instCnt]: VNF_TO_InstancesCnt){
        VNFNetwork->VNFNodes[fnType]->numInstances = instCnt;
    }
}

template<typename type_res>
void assign_ForSFC_VNFType_2_InstID( const vector<pair<int,int>>& mapping, ServiceFunctionChain* cSFC, VirtualNetworkFunctions<type_res> *VNFNetwork){
    for(const auto& [fnType, fnInstId]: mapping){
        VNFNetwork->ppar_utilization[fnType][fnInstId] += cSFC->trafficArrivalRate;
        VNFNetwork->seq_utilization[fnType][fnInstId] += cSFC->trafficArrivalRate;

    }
}
#endif //SFC_PARALLELIZATION_ASSIGNMENTFUNCTIONS_H
