//
// Created by vijay on 22-03-2023.
//

#ifndef SFC_PARALLELIZATION_FILEFUNCTIONS_H
#define SFC_PARALLELIZATION_FILEFUNCTIONS_H

/*!
 * Create a directory of output folder with test name.
 * @param fullDirName name of the directory with slash at the end.
 * @return true if file does not exist and created now. Else False.
 */
bool createDirectory(const string& fullDirName){
//    if(!filesystem::exists(input_directory + fullDirName) and filesystem::create_directory(input_directory + fullDirName)) {
//        cout << "Input Directory:[" << input_directory + fullDirName << "] created."<<endl;
//    }
    if(!filesystem::exists(output_directory + fullDirName) and filesystem::create_directory(output_directory + fullDirName)) {
        cout<<"Output Directory:["<<output_directory + fullDirName<<"] created." <<endl;
    }
    if(!filesystem::exists(output_directory + fullDirName+ diagram_directory) and filesystem::create_directory(output_directory + fullDirName+ diagram_directory)) {
        cout<<"Output Graph Directory:["<<output_directory + fullDirName+ diagram_directory<<"] created."<<endl;
        return true;
    }
    return false;
}

/*!Function will read the network data from the input directory and read into graph structure. \n
 * Also calculate the all pairs shortest path. Set number of edges read. \n
 * First line in file consists of numberofNodes N and whether to read data from adjaceny List (L) or Adjacency matrix(M). \n
 * Second line containes number of servers ns to read from it (it should be same as number of nodes as we consider nodes == servers).\n
 * Each next ns line conatins Node_index NumOfCoresInThatServer. (One-Based Indexing) Node Index should be within network range of Num of Nodes.\n
 * If Adj List -> Rest of the file consists of edges data (unidrected) in each line - u v wgt. \n
 *      where u is src, v is the dst and wgt is the distance between them \n
 *      distance is in meters. \n
 * If Adj Matrix -> then N*N matrix to read data with value as weight. Follows One-Based Index.
 * @param dirName name of the directory inside which it contains network file.
 * @param filename_network network data filepath from where we have to read network
 * @param graph PhysicalGraph class pointer variable to save the graph values.
 * @param CoresOpt number of cores in each node. \n
 * {"node-degree", 0} -> according to node degree it will fixed cores \n
 * {"fixed-all", val} -> fixe each node core to value = val \n
 * {"variable-all", maxval} -> scale according to node degree with maxval \n
 * @param showinConsole show reading values in the console to debug.
 */
bool readNetwork(const std::string& dirName, const std::string& filename_network, PhysicalGraph& graph, const pair<string, unsigned int>& CoresOpt = {"fixed-all", 3}, bool showinConsole=false)
{//readNetwork
    if(showMorePrecision) std::cout.precision( numeric_limits<double>::digits10  );

//    auto coresAssignmentLogic = [&](const unsigned int& deg) -> unsigned int{
//        if(CoresOpt.first == "node-degree")
//            return deg;
//        else if(CoresOpt.first == "fixed-all"){
//            return CoresOpt.second;
//        }else if(CoresOpt.first == "variable-all"){
//            if(deg == 1) return deg;
//            else if(deg <= 5) return 2;
//            else if(deg <= 9) return 4;
//            else return CoresOpt.second;
//        }
//        else return 2;
//    };

    ifstream fin;
    string filepathExt = input_directory + dirName + filename_network;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear(); throw runtime_error(errorMsg+ __FUNCTION__);
    }
    if(debug)cout<<"\n[Reading Network] Network File: "<<filepathExt<<endl;

    unsigned int iG_nNodes, iG_nEdges; ///< Input Graph Number of Nodes and Edges in it.
    char edgesReadingAs; ///< reading edges as Adjacency Matrix (M) or List (L)
    unsigned int readEdge_u, readEdge_v; ///< edge source and destination
    type_wgt readEdge_wgt; ///< edge weight
    /// if successfully read number of nodes and edges of graph
    if(fin>>iG_nNodes>>edgesReadingAs and iG_nNodes){
        graph = std::move(PhysicalGraph(iG_nNodes,iG_nEdges));

        unsigned int readNode_idx, readNodeCap_deg;
        unsigned int ns=0; /// number of servers
        fin>>ns;
        for(int pi=1; pi<=ns; pi++) {
            if(!(fin>>readNode_idx>>readNodeCap_deg)){
                std::cerr << "\t Server Reading Failed.\n";
                return false;
            }
//            unsigned int num_cores = coresAssignmentLogic(readNodeCap_deg);
            string readNode_name = "PN"+to_string(readNode_idx);
            graph.PNode[readNode_idx] = PhysicalNode(readNode_idx, readNode_name, NodeCapacity(readNodeCap_deg));
            graph.sum_cores += readNodeCap_deg; ///< summation of all cores in the network
            if(showinConsole)cout<<readNode_idx<<"-"<<readNode_name<<"-"<<readNodeCap_deg<<endl;
        }

        unsigned int ne=0;/// number of edges
//        type_wgt minwgt = std::numeric_limits<type_wgt>::max(), maxwgt = std::numeric_limits<type_wgt>::min();
        if(edgesReadingAs == 'L'){ ///Reading Adj List.
            while(!fin.eof()) {
                if(!(fin>>readEdge_u>>readEdge_v>>readEdge_wgt)){
                    std::cerr << "\t Edge Reading Failed: Edge Row["<<ne<<"]\n";
                    return false;
                } ++ne;
                graph.mat[readEdge_u][readEdge_v] = graph.mat[readEdge_v][readEdge_u] = readEdge_wgt;
                graph.PNode[readEdge_u].neighbours.insert(readEdge_v);
                graph.PNode[readEdge_v].neighbours.insert(readEdge_u);

//                minwgt = min(minwgt, readEdge_wgt);
//                maxwgt = max(maxwgt, readEdge_wgt);
                if(showinConsole)cout<<"("<<readEdge_u<<","<<readEdge_v<<")d="<<readEdge_wgt<<endl;
            }
        } else if(edgesReadingAs == 'M'){ // Reading Adj Matrix
            for(unsigned int src=graph.srcV; src<=iG_nNodes; src++){
                for(unsigned int dst=graph.srcV; dst<=iG_nNodes; dst++){
                    if(!(fin>>readEdge_wgt)){
                        std::cerr << "\t Edge Reading Failed: mat["<<src<<"]["<<dst<<"]\n";
                        return false;
                    }
                    if(readEdge_wgt>0){
                        ++ne;
                        graph.mat[src][dst] = readEdge_wgt;
                        graph.PNode[src].neighbours.insert(dst);
                        if(showinConsole)cout<<"("<<readEdge_u<<","<<readEdge_v<<")d="<<readEdge_wgt<<endl;
                    }
                }//dst
            }//src
            ne = ne/2; //undirected graph each edge counted twice.
        }else return false;
        graph.setNumOfEdges(ne); ///< set number of edges.
        if(debug){ cout<<"\tG("<<iG_nNodes<<","<<ne<<") | Servers: "<<ns; }
//        cout<<"\n minwgt:"<<minwgt<<" maxwgt:"<<maxwgt<<";";
    } else{
        string errorMsg = "Invalid Input Values V and M/L. File "+filepathExt+ ". Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fin.close();
    graph.filename_network = filename_network.substr(0,filename_network.find('.'));
    graph.calcClusteringCoefficient();
    graph.calcAllPairsShortestPath();
    graph.setNodesCores(CoresOpt);
    return true;
}//readNetwork

/*! Function will read the VNF data from the file in the input directory and read into VNFNode class and VirtualNetworkFunction class.\n\n
 * Also Calculate Parallel Pairs. \n\n
 * First Line in file consist of numOfVNF F\n
 * from next line upto F times, each group consist\n
 * Name of VNF\n
 * VNF_index  Requirement[cores]     serviceRate     execTime
 * @param dirName name of the directory inside which it contains vnf file.
 * @param filename_vnf  current VNFs file which consists of VNFs data
 * @param allVNFs VirtualNetworkFunctions class pointer variable to store value of VNFs
 * @param showinConsole show the file reading values in the console
 */
bool readVirtualNetworkFunctions(const std::string& dirName, const std::string& filename_vnf, VirtualNetworkFunctions& allVNFs, bool showinConsole=false) {
    ifstream fin;
    string filepathExt = input_directory + dirName + filename_vnf;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear(); throw runtime_error(errorMsg+ __FUNCTION__);
    }

    if(debug)cout<<"\n[Reading Virtual Network Functions] VNF File: "<<filepathExt<<endl;

    unsigned int i_nVNF; // input num of VNFs
    if(fin>>i_nVNF){
        allVNFs = std::move(VirtualNetworkFunctions(i_nVNF));
        string readVNF_name; unsigned int readVNF_idx, readVNF_Req_cores=0; type_delay readVNF_serviceR=0,readVNF_execTime=0;
        for(int ni=1; ni<=i_nVNF; ni++){ /// For each node there is a row which would be read
            fin.ignore(); getline(fin, readVNF_name);
            if (!(fin>>readVNF_idx>>readVNF_Req_cores>>readVNF_serviceR>>readVNF_execTime)){
                std::cerr << "\tVNF Reading Failed: VNF Row["<<ni<<"]\n";
                return false;
            }
            allVNFs.VNFNodes[readVNF_idx] = VNFNode(readVNF_idx, readVNF_name, readVNF_serviceR, readVNF_execTime,NodeCapacity(readVNF_Req_cores) );
            if(showinConsole)cout<<readVNF_idx<<"-"<<readVNF_name<<"-"<<readVNF_serviceR<<"-"<<readVNF_execTime<<"-["<<readVNF_Req_cores<<"]"<<endl;
        }

    }else{
        string errorMsg = "Invalid Input File Values. File "+filepathExt+ ". Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fin.close();
    allVNFs.filename_vnf = filename_vnf.substr(0,filename_vnf.find('.'));
    if(debug)cout<<"\tVNFs:"<<i_nVNF<<" ("<<allVNFs.filename_vnf<<")";
    return true;
}

#endif //SFC_PARALLELIZATION_FILEFUNCTIONS_H
