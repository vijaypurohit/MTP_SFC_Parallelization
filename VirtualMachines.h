//
// Created by vijay on 23-02-2023.
//

#ifndef SFC_PARALLELIZATION_VIRTUALMACHINES_H
#define SFC_PARALLELIZATION_VIRTUALMACHINES_H

/*!
 * Virtual Machine Node
 * @tparam type_res resource data type. default=unsigned int.
 */
template <class type_res =unsigned int>
class VirtualMachineNode{
public:
    /*! physical node index. Type int.  */
    int index{};  string name; ///<node name
    NodeCapacity<type_res> capacity, requirement; /*!< capacity of the VM node to host VNFs. Type NodeCapacity. \n requirements  of the VM node for hosting on physical node. Type NodeCapacity.*/
    vector<pair<int,int>> VNFs; ///< vnfs indexes hosted on vm id.
    /*!
     * @param _index physical node index. Type int.
     * @param _name node name
     * @param _givenCapacity capacity of the VM node to host VNFs. Type NodeCapacity.
     * @param _givenRequirements requirements  of the VM node for hosting on physical node. Type NodeCapacity.
     */
    VirtualMachineNode(int _index, const string& _name, NodeCapacity<type_res> _givenCapacity, NodeCapacity<type_res> _givenRequirements): capacity(_givenCapacity), requirement(_givenRequirements){
        this->index = _index; this->name = _name;
    }
    ~VirtualMachineNode() = default;
};

/*!
 * Virtual Machines Collection Data, link with physical node
 * @tparam type_res resource data type. default=unsigned int.
 */
template <class type_res =unsigned int>
class VirtualMachines
{
    unsigned int numVM{}, srcVM; ///<number of virtual machines . \n srcVM = default 1. Loop starts from 1 to <=numVM
public:
    vector<int> I_VM2PN; ///< VM index is hosted on which physical machine
    unordered_map<unsigned int, VirtualMachineNode<type_res>*> VMNode; ///<Index to Virtual machine address
    /*!
     * @param _numVirtualMachines number of virtual machines
     */
    explicit VirtualMachines(unsigned int _numVirtualMachines)
    {
        srcVM = 1;
        this->numVM = _numVirtualMachines;
        I_VM2PN = vector<int>(numVM + 1, -1);
    }
    ~VirtualMachines(){
        for (unsigned int v = srcVM; v <= numVM; ++v) { //deletion of Adj List Edges
            delete VMNode[v]; // delete all VMs
        }
        if(debug) cout<<"\n[VirtualMachines Destructor Completed]";
    }

    void showVMs_Description();
    void assignVMtoPN(int vmIndex, int pnIndex);
};

/*! Show all the Virtual Machine nodes and their description
 * @tparam type_res resource data type. default=unsigned int.
 */
template <class type_res>
void VirtualMachines<type_res>::showVMs_Description(){
    cout << "\n\n ----- Virtual Machines Description ::";
    cout<<"\nIndex\t"<<"Name\t"<<" CAP[cores][memory][disk][speed]\t"<<" Req[cores][memory][disk][speed]\t " <<" PNode \t"<<" VNFs";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (unsigned int u = srcVM; u <= numVM; ++u){
        cout<<"\n"<<VMNode[u]->index<<" |\t";
        cout<<VMNode[u]->name<<" |\t[";
        cout<<VMNode[u]->capacity.cores<<"\t"<<VMNode[u]->capacity.memory<<"\t"<<VMNode[u]->capacity.disk<<"\t"<<VMNode[u]->capacity.cpuSpeed<<"]\t[";
        cout<<VMNode[u]->requirement.cores<<"\t"<<VMNode[u]->requirement.memory<<"\t"<<VMNode[u]->requirement.disk<<"\t"<<VMNode[u]->requirement.cpuSpeed<<"]\t";
        cout << "PN[" << I_VM2PN[u] << "]\t";
        for(const auto& vnfid: VMNode[u]->VNFs)
            cout<<"F["<<vnfid.first<<char(96+vnfid.second)<<"] ";
    }
}

 /*!
  * @param vmIndex Virtual Machine Index
  * @param pnIndex Physical Node Index
  * @tparam type_res resource data type. default=unsigned int.
  */
template <class type_res>
void VirtualMachines<type_res>::assignVMtoPN(int vmIndex, int pnIndex){
     I_VM2PN[vmIndex] = pnIndex;
}


#endif //SFC_PARALLELIZATION_VIRTUALMACHINES_H
