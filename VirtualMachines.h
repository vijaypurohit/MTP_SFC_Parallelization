//
// Created by vijay on 23-02-2023.
//

#ifndef SFC_PARALLELIZATION_VIRTUALMACHINES_H
#define SFC_PARALLELIZATION_VIRTUALMACHINES_H

/*!
 * Virtual Machine Node
 */
class VirtualMachineNode{
public:
    /*! physical node index. Type int.  */
    unsigned int index{};  string name; ///<node name
    NodeCapacity capacity, requirement; /*!< capacity of the VM node to host VNFs. Type NodeCapacity. \n requirements  of the VM node for hosting on physical node. Type NodeCapacity.*/
    vector<pair<unsigned int,unsigned int>> vm2VNFs; ///< vnfs indexes hosted on vm id. {vnf_idx, instance_idx (1-based indexing)}
    /*!
     * @param _index physical node index. Type int.
     * @param _name node name
     * @param _givenCapacity capacity of the VM node to host VNFs. Type NodeCapacity.
     * @param _givenRequirements requirements  of the VM node for hosting on physical node. Type NodeCapacity.
     */
    VirtualMachineNode(unsigned int _index, const string& _name, NodeCapacity _givenCapacity, NodeCapacity _givenRequirements): capacity(_givenCapacity), requirement(_givenRequirements){
        this->index = _index; this->name = _name;
    }
    ~VirtualMachineNode() = default;
};

/*!
 * Virtual Machines Collection Data, link with physical node
 */
class VirtualMachines
{
    unsigned int numVM, srcVM; ///<number of virtual machines . \n srcVM = default 1. Loop starts from 1 to <=numVM
public:
    unordered_map<unsigned int, unsigned int> I_VM2PN; ///< VM index is hosted on which physical machine
    unordered_map<unsigned int, VirtualMachineNode*> VMNode; ///<Index to Virtual machine address
    /*!
     * @param _numVirtualMachines number of virtual machines
     */
    explicit VirtualMachines(unsigned int _numVirtualMachines)
    {
        srcVM = 1; this->numVM = _numVirtualMachines;
    }
    ~VirtualMachines(){
        for(const auto vmInfo: VMNode){
            delete vmInfo.second;
        }
//        for (unsigned int v = srcVM; v <= numVM; ++v) { //deletion of Adj List Edges
//            delete VMNode[v]; // delete all VMs
//        }
        if(debug) cout<<"\n[VirtualMachines Destructor Completed]";
    }

    void showVMs_Description();
    void assignVMtoPN(int vmIndex, int pnIndex);
};

/*! Show all the Virtual Machine nodes and their description
 * @tparam type_res resource data type.
 */
void VirtualMachines::showVMs_Description(){
    cout << "\n\n ----- Virtual Machines Description ::";
    cout<<"\nIndex\t"<<"Name\t"<<" CAP[cores]\t"<<" Req[cores]\t " <<" PNode \t"<<" VNFs";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (unsigned int u = srcVM; u <= numVM; ++u){
        const VirtualMachineNode* vmInfo = VMNode.at(u);
        cout<<"\n"<<vmInfo->index<<" |\t";
        cout<<vmInfo->name<<" |\t[";
        cout<<vmInfo->capacity.cores<<"]\t[";
        cout<<vmInfo->requirement.cores<<"]\t";
        cout << "PN[" << I_VM2PN.at(u) << "]\t";
        for(const auto& vnfid: vmInfo->vm2VNFs)
            cout<<"F["<<vnfid.first<<char(96+vnfid.second)<<"] ";
    }
}

 /*!
  * @param vmIndex Virtual Machine Index
  * @param pnIndex Physical Node Index
  * @tparam type_res resource data type.
  */
void VirtualMachines::assignVMtoPN(int vmIndex, int pnIndex){
     I_VM2PN[vmIndex] = pnIndex;
}


#endif //SFC_PARALLELIZATION_VIRTUALMACHINES_H
