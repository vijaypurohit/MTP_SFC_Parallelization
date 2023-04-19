//
// Created by vijay on 22-03-2023.
//

#ifndef SFC_PARALLELIZATION_FILEFUNCTIONS_H
#define SFC_PARALLELIZATION_FILEFUNCTIONS_H

/*!
 * Create a directory in input and output folder with test name.
 * @param testDirName name of the directory to be created
 * @return true if file does not exist and created now. Else False.
 */
bool createDirectory(const string& testDirName){
    if(!filesystem::exists(input_directory + testDirName) and filesystem::create_directory(input_directory + testDirName)) {
        cout << "Input Directory:[" << input_directory + testDirName << "] created."<<endl;
    }
    if(!filesystem::exists(output_directory + testDirName) and filesystem::create_directory(output_directory + testDirName)) {
        cout<<"Output Directory:["<<output_directory + testDirName<<"] created." <<endl;
    }
    if(!filesystem::exists(output_directory + testDirName+ diagram_directory) and filesystem::create_directory(output_directory + testDirName+ diagram_directory)) {
        cout<<"Output Graph Directory:["<<output_directory + testDirName+ diagram_directory<<"] created."<<endl;
        return true;
    }
    return false;
}

void readConstants(const string& testDirName)
{
    ifstream fin;
    string filepathExt = input_directory + testDirName + filename_constants;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    cout<<"\n Reading File: "<<filepathExt<<endl;
//    fin>>iG_nNodes>>iG_nEdges;

}

/*!
 * Function will read the network data from the input directory and read into graph structure. \n\n
 * Also calculate the all pairs shortest path. Set number of edges read. \n \n
 * First line in file consists of numberofNodes N and NumberOfedges E. Edges can be give approx value. \n
 *      next 2*N Lines consists of for each N nodes --> first row - Name of Node. \n
 * second row - values of index, node capacity(cores, memory, disk, speed). \n
 * Rest of the file consists of edges data (unidrected) in each line - u v wgt. \n
 *      where u is src, v is the dst and wgt is the distance between them \n
 *      distance is in meters. \n
 * For Example \n
 * N E \n
 * ServerName \n
 * index    Capacity[cores  memory  disk    speed] \n
 * u_idx    v_idx   distance_meter
 * @param testDirName path to the current test directory which consists of inputs files
 * @param graph PhysicalGraph class pointer variable to store values.
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template<typename type_wgt, typename type_res>
void readNetwork(const string& testDirName, PhysicalGraph<type_wgt, type_res> **const& graph, bool showinConsole=false)
{
    ifstream fin;
    string filepathExt = input_directory + testDirName + filename_network;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    if(debug)cout<<"\n[Reading Network] Network File: "<<filepathExt<<endl;

//    PhysicalGraph<type_wgt, type_res> *graph = nullptr;
    unsigned int iG_nNodes, iG_nEdges;
    /// if successfully read number of nodes and edges of graph
    if(fin>>iG_nNodes>>iG_nEdges and iG_nNodes){
        *graph = new PhysicalGraph<type_wgt, type_res>(iG_nNodes,iG_nEdges);
        string readNode_name; unsigned int readNode_idx, readNodeCap_cores;
        type_res readNode_memory, readNode_disk, readNode_cpuspeed;
        unsigned int readEdge_u, readEdge_v;
        type_wgt readEdge_wgt;
        for(unsigned int ni=1; ni<=iG_nNodes; ni++) /// For each node there is a row which would be read
        {
            fin.ignore();
            getline(fin, readNode_name);
            if (!(fin>>readNode_idx>>readNodeCap_cores>>readNode_memory>>readNode_disk>>readNode_cpuspeed))
                std::cerr << "\t Node Reading Failed: Node Row["<<ni<<"]\n";
            (*graph)->PNode[readNode_idx] = new PhysicalNode<type_res>(readNode_idx, readNode_name, NodeCapacity<type_res>(readNodeCap_cores, readNode_memory, readNode_disk, readNode_cpuspeed));
            if(showinConsole)cout<<readNode_idx<<"-"<<readNode_name<<"-"<<readNodeCap_cores<<"-"<<readNode_memory<<"-"<<readNode_disk<<"-"<<readNode_cpuspeed<<endl;
        }
        if(debug)cout<<"\tNodes:"<<iG_nNodes<<"\n";
        int ne=0;/// number of edges
        while(!fin.eof()) {
            if(!(fin>>readEdge_u>>readEdge_v>>readEdge_wgt))
                std::cerr << "\t Edge Reading Failed: Edge Row["<<ne<<"]\n";
            ++ne;
            (*graph)->addEdge(new PhysicalEdge<type_wgt>(readEdge_u, readEdge_v, readEdge_wgt));
            if(showinConsole)cout<<readEdge_u<<" - "<<readEdge_v<<" - "<<readEdge_wgt<<endl;
        }
        (*graph)->setNumOfEdges(ne); // set number of edges.
        if(debug)cout<<"\tEdges:"<<ne;

    } else{
        string errorMsg = "Invalid Input File Values V and E. File "+filepathExt+ ". Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fin.close();
    (*graph)->calcAllPairsShortestPath();
//    if(debug)cout<<"\n[Function Completed: "<<__FUNCTION__<<"]";
}

/*! Function will read the Virtual Machine data from the file in the input directory and read into VirtualMachine class.\n\n
 * First Line in file consist of numOfVM M\n
 * from next line upto M times, each group consist\n
 * Name of VM\n
 * VM_index\n
 * Capacity[cores  memory  disk    speed]\n
 * Requirement[cores  memory  disk    speed]
 * @param testDirName  path to the current test directory which consists of inputs files
 * @param ptr VirtualMachines class pointer variable to store value
 * @param showinConsole show the file read values in the console
 * @tparam type_res resource data type.
 */
template<typename type_res>
void readVirtualMachines(const string& testDirName, VirtualMachines<type_res> **const& ptr, bool showinConsole=false) {
    ifstream fin;
    string filepathExt = input_directory + testDirName + filename_virtualmachines;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    if(debug)cout<<"\n["<<__FUNCTION__<<"] VM File:"<<filepathExt<<endl;
    unsigned int i_nVM; // input num of VMs
    if(fin>>i_nVM){
        *ptr = new VirtualMachines<type_res>(i_nVM);
        string readVM_name;   unsigned int readVM_idx, readVM_Cap_cores, readVM_Req_cores;
        type_res readVM_Cap_memory, readVM_Cap_disk, readVM_Cap_cpuspeed, readVM_Req_memory, readVM_Req_disk, readVM_Req_cpuspeed;
        for(int ni=1; ni<=i_nVM; ni++) /// For each node there is a row which would be read
        {
            fin.ignore();
            getline(fin, readVM_name);
            if (!(fin>>readVM_idx>>readVM_Cap_cores>>readVM_Cap_memory>>readVM_Cap_disk>>readVM_Cap_cpuspeed>>readVM_Req_cores>>readVM_Req_memory>>readVM_Req_disk>>readVM_Req_cpuspeed))
                std::cerr << "\t VM Data Reading Failed: VM Row["<<ni<<"]\n";
            (*ptr)->VMNode[readVM_idx] = new VirtualMachineNode<type_res>(readVM_idx, readVM_name,
                                                                          NodeCapacity<type_res>(readVM_Cap_cores, readVM_Cap_memory, readVM_Cap_disk, readVM_Cap_cpuspeed),
                                                                          NodeCapacity<type_res>(readVM_Req_cores, readVM_Req_memory, readVM_Req_disk, readVM_Req_cpuspeed) );
            if(showinConsole)cout<<readVM_idx<<"-"<<readVM_name<<"-["<<readVM_Cap_cores<<"-"<<readVM_Cap_memory<<"-"<<readVM_Cap_disk<<"-"<<readVM_Cap_cpuspeed
                                 <<"]-["<<readVM_Cap_cores<<"-"<<readVM_Cap_memory<<"-"<<readVM_Cap_disk<<"-"<<readVM_Cap_cpuspeed<<"]"<<endl;
        }
        if(debug)cout<<"\t[VM Data Reading Completed] VMs:"<<i_nVM;
    }else{
        string errorMsg = "Invalid Input File Values. File "+filepathExt+ ". Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fin.close();
//    if(debug)cout<<"\n[Function Completed: "<<__FUNCTION__<<"]";
}

/*! Function will read the VNF data from the file in the input directory and read into VNFNode class and VirtualNetworkFunction class.\n\n
 * Also Calculate Parallel Pairs. \n\n
 * First Line in file consist of numOfVNF F\n
 * from next line upto F times, each group consist\n
 * Name of VNF\n
 * VNF_index    NumOfInstances      serviceRate     execTime    Requirement[cores  memory  disk    speed]
 * @param testDirName  path to the current test directory which consists of inputs files
 * @param ptr VirtualNetworkFunctions class pointer variable to store value
 * @param showinConsole show the file read values in the console
 * @param
 * @tparam type_res resource data type.
 */
template<typename type_res>
void readVirtualNetworkFunctions(const string& testDirName, VirtualNetworkFunctions<type_res> **const& ptr, bool showinConsole=false) {
    ifstream fin;
    string filepathExt = input_directory + testDirName + filename_vnf;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    if(debug)cout<<"\n[Reading Virtual Network Functions] VNF File: "<<filepathExt<<endl;

    unsigned int i_nVNF; // input num of VNFs
    if(fin>>i_nVNF){
        *ptr = new VirtualNetworkFunctions<type_res>(i_nVNF);
        string readVNF_name; unsigned int readVNF_idx, readVNFInst; type_delay readVNF_serviceR,readVNF_execTime ; unsigned int readVNF_Req_cores;
        type_res readVNF_Req_memory, readVNF_Req_disk, readVNF_Req_cpuspeed;
        for(int ni=1; ni<=i_nVNF; ni++) /// For each node there is a row which would be read
        {
            fin.ignore();
            getline(fin, readVNF_name);
            if (!(fin>>readVNF_idx>>readVNFInst>>readVNF_serviceR>>readVNF_execTime>>readVNF_Req_cores>>readVNF_Req_memory>>readVNF_Req_disk>>readVNF_Req_cpuspeed))
                std::cerr << "\tVNF Reading Failed: VNF Row["<<ni<<"]\n";
            (*ptr)->VNFNodes[readVNF_idx] = new VNFNode<type_res>(readVNF_idx, readVNF_name, readVNFInst, readVNF_serviceR, readVNF_execTime,
                                                                 NodeCapacity<type_res>(readVNF_Req_cores, readVNF_Req_memory, readVNF_Req_disk, readVNF_Req_cpuspeed) );
            if(showinConsole)cout<<readVNF_idx<<"-"<<readVNF_name<<"-"<<readVNFInst<<"-"<<readVNF_serviceR<<"-"<<readVNF_execTime<<"-["<<readVNF_Req_cores<<"-"<<readVNF_Req_memory<<"-"<<readVNF_Req_disk<<"-"<<readVNF_Req_cpuspeed<<"]"<<endl;
        }
        if(debug)cout<<"\tVNFs:"<<i_nVNF;
    }else{
        string errorMsg = "Invalid Input File Values. File "+filepathExt+ ". Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fin.close();
//    (*ptr)->findRandomParallelPairs(testDirName);
//    if(debug)cout<<"\n[Function Completed: "<<__FUNCTION__<<"]";
}

/*! Function will read the SFC chains data from the file in the input directory and read into ServiceFunctionChain class.\n\n
 * It will simultaneously read vnf data into sequential vector and also construct adj list. \n\n
 * First Line in file consist of numOfSFC S\n
 * from next line upto S times, each s in S group consist\n
 * Name of SFC\n
 * sfc_index    numOfVNFs present except src and dest   arrivalRate of SFC\n
 * This line contains VNFs Type id (for i=1 to numOFVNFs)
 * @param testDirName  path to the current test directory which consists of inputs files
 * @param SFC vector of ServiceFunctionChain objects to store values.
 * @param VNFNetwork VirtualNetworkFunctions class object used to store vnf to sfc mapping.
 * @param showinConsole show the file read values in the console
 * @tparam type_res resource data type.
 */
template<typename type_res>
void readGenericServiceFunctionsChains(const string& testDirName, vector<ServiceFunctionChain*>& allSFC, vector<ServiceFunctionChain *> &sortedSFCs, VirtualNetworkFunctions<type_res> *const& VNFNetwork ,bool showinConsole=false) {
    ifstream fin;
    string filepathExt = input_directory + testDirName + filename_sfc;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    priority_queue<ServiceFunctionChain *, vector<ServiceFunctionChain *>, comparator_sfc> pqSortSFC;
    if(debug)cout<<"\n[Reading Original Sequential SFCs] SFC File: "<<filepathExt<<endl;
    unsigned int i_nSFC; ///< input num of SFCs
    if(fin>>i_nSFC){
        allSFC = vector<ServiceFunctionChain*>(i_nSFC);
        string readSFC_name; unsigned int readSFC_idx, readSFC_totalVNF, vnfid; type_delay readSFC_arrivalRate ;
        for(int ni=0; ni<i_nSFC; ni++) ///< For each node there is a row which would be read
        {
            fin.ignore();
            getline(fin, readSFC_name);
            if (!(fin>>readSFC_idx>>readSFC_totalVNF>>readSFC_arrivalRate))
                std::cerr << "\t SFCs Reading Failed: SFC Row["<<ni<<"]\n";
            allSFC[ni] = new ServiceFunctionChain(readSFC_idx, readSFC_name, readSFC_totalVNF, readSFC_arrivalRate);
            for(int vj = 1; vj<=readSFC_totalVNF; vj++) { ///< read VNFs type ID.
                if (fin >> vnfid) {
                    allSFC[ni]->vnfSeq.push_back(vnfid);
                    VNFNetwork->cntVNF[vnfid].push_back(ni);
                } else std::cerr << "\t SFCs Reading Failed: SFC Row[" << ni << "] fn:" << vnfid << "\n";
            }
            pqSortSFC.push(allSFC[ni]);

            if(showinConsole) {
                 cout << readSFC_idx << "-" << readSFC_name << "-" << readSFC_totalVNF << "-" << readSFC_arrivalRate<<"- [";
                 for(const auto& x: allSFC[ni]->vnfSeq)
                     cout<<x<<" -> "; cout<<"]\n";
             }
        }
        if(debug)cout<<"\tSFCs:"<<allSFC.size();
    }else{
        string errorMsg = "Invalid Input File Values. File "+filepathExt+ ". Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fin.close();

    while (!pqSortSFC.empty()) {
        sortedSFCs.push_back(pqSortSFC.top());
        pqSortSFC.pop();
    }
//    if(debug)cout<<"\n[Function Completed: "<<__FUNCTION__<<"]";
}


#endif //SFC_PARALLELIZATION_FILEFUNCTIONS_H
