//
// Created by vijay on 23-02-2023.
//

#ifndef SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H
#define SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H


/*! @brief Service Function Class */
class ServiceFunctionChain
{
public:
    unsigned int index /*!< SFC index */; string name;  /*!< SFC name */ unsigned int numVNF;  /*!< number of VNFs except source and dest node */
    type_delay trafficArrivalRate;  /*!< traffic arrival rate of SFC s. in packets per second. arrival rate < service rate. */
    pair<unsigned int,unsigned int> access_nodes; ///< source and destination server index from where sfc originates and goes to.
    vector<unsigned int> vnfSeq; ///< vnfids in sequential manner excluding src (SFCsrc=0) and dest (SFCdst=-10).

/// parameters found through algorithms.
    vector<vector<unsigned int>> vnfBlocksPar; ///< vnfids in stage wise, where stg i denotes parallel vnf in ith stage. Except src and dest stg.
    vector<vector<vector<unsigned int>>> allPartParSFC; ///< All the partial parallel VNF blocks of the given sequential SFC, order of parallel vnfs are changed here.
    vector<vector<vector<unsigned int>>> subsetPartParSFC; ///< subset of all the partials chains found which are in order as of vnfSeq.
    vector<vector<vector<unsigned int>>>* partialParallelChains; ///< partial chains in consideration whether allPartParSFC or subset
   /*!
     * @param _index SFC index  @param _name SFC name
     * @param _numVNF number of VNFs except source and dest node
     * @param _trafficRate traffic arrival rate of SFC s
     */
    ServiceFunctionChain(unsigned int _index, string& _name, pair<unsigned int, unsigned int>& _accessnode, unsigned int _numVNF, type_delay _trafficArrivalRate){
        this->index = _index; this->name = std::move(_name);
        this->access_nodes = std::move(_accessnode);
        this->numVNF = _numVNF;
        this->trafficArrivalRate = _trafficArrivalRate;
    }
    ServiceFunctionChain()=default;
    ~ServiceFunctionChain()=default;

    void showSequentialSFC() const;
    void showFullyParallelSFC() const;
    void showAllPartialSFC() const;
    void showSubsetPartialSFC() const;

    [[maybe_unused]] void printOrigSFC(int, const string&);
};

/*!
 * Show Original Sequential SFC in console.
 * @param VNFType2Inst function to instance mapping.
 */
void ServiceFunctionChain::showSequentialSFC() const{
    cout<<"\nSequential SFCid:"<<index<<" | Len: "<<numVNF<<" | ArrivalRate: "<<trafficArrivalRate<<" | [#subsetPartial: "<<subsetPartParSFC.size()<<" #allPartial: "<<allPartParSFC.size()<<"]\n\t";
    cout << "(ps:"<<access_nodes.first<<" -> ";
    for(const auto& fn : vnfSeq){cout <<"f"<< fn << "; ";} cout << " pd:"<<access_nodes.second<<")";
}

/*!
 * Show Fully Parallel SFC in console.
 * @param VNFType2Inst function to instance mapping.
 */
void ServiceFunctionChain::showFullyParallelSFC() const{
    cout<<"\nFully-Parallel SFCid:"<<index<<" | Len: "<<numVNF<<" | ArrivalRate: "<<trafficArrivalRate<<" | [#subsetPartial: "<<subsetPartParSFC.size()<<" #allPartial: "<<allPartParSFC.size()<<"]\n\t";
    cout << "(ps:"<<access_nodes.first<<" -> ";
    for(const auto& blk: vnfBlocksPar){ cout<<" ["; for(int fn: blk){ cout <<"f"<< fn <<"; ";  }  cout<<"] ";  } cout << " pd:"<<access_nodes.second<<")";
}

/*!
 * Show all the partial SFC of the given SFC.
 * @param VNFType2Inst function to instance mapping.
 */
void ServiceFunctionChain::showAllPartialSFC() const{
    cout << "\nAll-Partial-Chains SFCid:"<<index<<" | Len: "<<numVNF<<" | ArrivalRate: "<<trafficArrivalRate<<" | [#subsetPartial: "<<subsetPartParSFC.size()<<" #allPartial: "<<allPartParSFC.size()<<"]\n\t";
    for (int ppid=0; ppid<allPartParSFC.size(); ppid++) {
        cout<<"\n\tpid:"<<ppid<<":: ";
        for (const auto &blk: allPartParSFC[ppid]) { cout << "["; for (int fn: blk) { cout << "f" << fn << "; "; } cout << "] "; }
    }
}

/*!
 * Show all the partial SFC of the given SFC.
 * @param VNFType2Inst function to instance mapping.
 */
void ServiceFunctionChain::showSubsetPartialSFC()const {
    cout << "\nSubset-Partial-Chains SFCid:"<<index<<" | Len: "<<numVNF<<" | ArrivalRate: "<<trafficArrivalRate<<" | [#subsetPartial: "<<subsetPartParSFC.size()<<" #allPartial: "<<allPartParSFC.size()<<"]\n\t";
    for (int spid=0; spid<subsetPartParSFC.size(); spid++) {
        cout<<"\n\tsid:"<<spid<<":: ";
        for (const auto &blk: subsetPartParSFC[spid]) {cout << "[";for (int fn: blk) { cout << "f" << fn << "; "; } cout << "] ";    }
    }
}

/*!
 * Print SFC Graph using Graphviz in a PNG format
 * @param type SFC type (Serial or parallel)
 */
void ServiceFunctionChain::printOrigSFC(const int type, const string& testDirName){
//    unordered_map<int, vector<int>>& adj = typeOfSFCAdj(type);
    unordered_map<int, vector<int>> adj;
    string sfcTypeName = "";//typeOfSFCString(type);
    if (adj.empty()) {
        cout << "SFC Graph is Empty."<<endl;
        return;
    }

    ofstream fout;
    string fileName = "SFC_"+ to_string(index)+"_"+sfcTypeName;///< (SFC + index + type) name without space
    string filepath = output_directory+testDirName+diagram_directory+fileName;///< path to .gv file without extention
    string filepathExt = filepath+".gv";///< filepath with extention of .gv
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
/**********************/
    queue<int> qLoGraph;		// queue level order Graph
    unordered_map<int, bool> visitedBFS; // bfs visited array, size= number of vertexes
    unordered_map<int, int> level; // level array, size= number of vertexes
    unordered_map<int, string> rankSame;// array to keep same rank node in one level, size= number of vertexes

    int max_level=1;
    int u_val, v_val;
    const string& nodeColor = "darkgreen", nodeFontColor = "darkgreen";
    const string& srcDstFillColor = "azure", srcDstFontColor = "darkorange3";
    const string& edgeColor = "darkolivegreen", edgeFontColor = "darkolivegreen";
/**********************/
    fout << "digraph "<<fileName<<"_numVNF_"<<numVNF<<" {" << endl;
    fout << "  label = \""<<name<<" "<<sfcTypeName<<"\" \n";
    fout << "  ranksep=\"equally\"; \n"
            "  rankdir=LR; " << endl << endl;
    fout << "node [fixedsize=shape width=0.45 shape=circle color="<<nodeColor<<" fontcolor="<<nodeFontColor<<"] "<<endl;
    fout << "edge [color="<<edgeColor<<" fontcolor="<<edgeFontColor<<" fontsize=12]"<<endl;
//    fout << SFCsrc << " [label = \"src\", fillcolor="<<srcDstFillColor<<" fontcolor="<<srcDstFontColor<<"]"<< endl;
//    fout << SFCdst << " [label = \"dst\", fillcolor="<<srcDstFillColor<<" fontcolor="<<srcDstFontColor<<"]"<< endl;

    int src = 0;
    qLoGraph.push(src);
    visitedBFS[src]=true;

    level[src]=0;//0
    rankSame[level[src]].append(to_string(src)).append("; ");

    while (!qLoGraph.empty())		// Level order Traversal
    {
        u_val = qLoGraph.front();  qLoGraph.pop();
        if(u_val != 0 and u_val != -1)
            fout << u_val << " [label = <&fnof;<sub>"<<u_val<<"</sub>>]"<< endl;

        for(auto w: adj[u_val]){
            v_val = w;        //adj node value
            if(!visitedBFS[v_val]){
                visitedBFS[v_val]=true;
                qLoGraph.push(v_val);
                level[v_val]=level[u_val]+1;
                if(level[v_val]+1>max_level)max_level++;
                rankSame[level[v_val]].append(to_string(v_val)).append("; ");
            }
            fout <<"\t"<< u_val << " -> " << v_val << endl;
        }//for adj
        fout << endl ;
    }//while qLoGraph

    //printing node with same level together
    for(int i=0; i<max_level; i++){
        if(!rankSame[i].empty())
            fout << "{rank = same; " << rankSame[i] << "}; " << endl;
    }//for i
    fout << "}" << endl;
/**********************/
    fout.close();

    string cmd = "dot -Tpng "+filepathExt+" -o "+filepath+".png";
    system((const char*)cmd.c_str());
    cout <<"\n Graphviz: SFC_"<<index<<" printed. File:"<<filepath<<".png \n";
    if(!debug)remove(filepathExt.c_str());
}//printSFC



#endif //SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H
