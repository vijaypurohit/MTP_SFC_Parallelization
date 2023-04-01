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
 * @tparam type_wgt edge weight data type. default=unsigned int.
 * @tparam type_res resource data type. default=unsigned int.
 */
template<typename type_wgt=unsigned int, typename type_res=unsigned int>
void assign_VM_2_PN(const vector<pair<int,int>>& VM_TO_PN, VirtualMachines<type_res> *VirtualNetwork, PhysicalGraph<type_wgt, type_res> *PhysicalNetwork){
    for(const auto& [vmid, pnid] : VM_TO_PN){
        VirtualNetwork->assignVMtoPN(vmid, pnid);
        PhysicalNetwork->PNode[pnid]->virtualMachines.push_back(vmid);
    }
}

/*!  Assign VNF instances to Virtual Machines. \n
 * Also make a mapping of VNF -> list of VNFs id.
 * @param VNF_TO_VM contains mapping of vnf id, instance id,  to vm id.
 * @param VNFNetwork VirtualNetworkFunctions class
 * @param VirtualNetwork VirtualMachines class
 * @tparam type_res resource data type. default=unsigned int.
 */
template<typename type_res=unsigned int>
void assign_VNF_2_VM(const vector<vector<int>>& VNF_TO_VM, VirtualNetworkFunctions<type_res> *VNFNetwork, VirtualMachines<type_res> *VirtualNetwork){

    for(const auto& assign : VNF_TO_VM){
        int vnfid = assign[0], instance = assign[1], vmid = assign[2];
        VNFNetwork->assignVNFinstToVM(vnfid, instance, vmid);
        VirtualNetwork->VMNode[vmid]->VNFs.emplace_back(vnfid,instance);
    }
}

/*!  Assign total instances required for particular VNF. \n
 * Also make a mapping of VNF -> list of VNFs id.
 * @param VNF_TO_VM contains mapping of vnf id, instance id,  to vm id.
 * @param VirtualNetwork VirtualMachines class
 * @tparam type_res resource data type. default=unsigned int.
 */
template<typename type_res=unsigned int>
void assign_VNF_2_InstancesCnt(const vector<pair<int,int>>& VNF_TO_InstancesCnt, VirtualNetworkFunctions<type_res> *VNFNetwork){
    for(const auto& [vnfid, instCnt]: VNF_TO_InstancesCnt){
        VNFNetwork->VNFNode[vnfid]->numInstances = instCnt;
        VNFNetwork->I_VNFinst2VM[vnfid] = vector<int>(instCnt + 1, -1); ///< assign memory according to num of instances.
    }
}

void assign_ForSFC_VNFType_2_InstID( const vector<pair<int,int>>& mapping, ServiceFunctionChain* obj){
    for(const auto& [vnfid, instid]: mapping){
        obj->I_VNFType2Inst[vnfid] = instid;
    }
}
#endif //SFC_PARALLELIZATION_ASSIGNMENTFUNCTIONS_H
