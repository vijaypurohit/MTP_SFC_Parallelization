//
// Created by vijay on 23-02-2023.
//

#ifndef SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H
#define SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H

class ServiceFunctionChain
{
public:
// given parameters
    unsigned int index /*!< SFC index */; string name;  /*!< SFC name */ unsigned int numVNF;  /*!< number of VNFs except source and dest node */
    float trafficArrivalRate;  /*!< traffic arrival rate of SFC s. in packets per second. arrival rate < service rate. */
    vector<int> vnfSeq; ///< vnfids in sequential manner including src (SFCsrc=0) and dest (SFCdst=-10).

// parameters found through algorithm.
    vector<vector<unsigned int>> vnfBlocksPar; ///< vnfids in stage wise, where stg i denotes parallel vnf in ith stage. Exvept src and dest stg.
    /*!
     * for final par SFC, VNF type is assigned to its which instance id i.e. {VNFid, instance id (1-based indexing)}. \n
     * value zero means not mapped. \n
     * initially mapped to value 1 at time of reading.
    */
    unordered_map<unsigned int, unsigned int> I_VNFType2Inst; ///< Best mapping for parallel length chain through algorithm
    vector<vector<vector<unsigned int>>> allPartParSFC; ///< All the partial paralle VNF blocks of the given sequential SFC.

    float bst_seqlen_time; ///< Best time for given sequential length chain according to our algorithm.
    unordered_map<unsigned int, unsigned int> bst_seqlen_mapping; ///< Best mapping {fun->its instance taken} for given sequential length chain according to our algorithm.
    unsigned int bst_parlen_idx; ///< idx of allPartParSFC array for which parallelsim give less time.
    float bst_parlen_time; ///< Best time for partial parallel chain with index bst_parlen_idx according to our algorithm.
    unordered_map<unsigned int, unsigned int> bst_parlen_mapping; ///< Best mapping for partial parallel chain with index bst_parlen_idx according to our algorithm.


//    unordered_map<unsigned int, bool> isVNF_Present;  /*!< given VNF index check, if it in SFC or not. {vnfid, true/false}. used in calcualtion of queuing delay */
//    unordered_map<int, vector<int>> pAdj;  /*!<  serial adjaceny list. parallel adjaceny list. */
    unordered_map<int, float> distanceSeq, distancePar; ///<stores longest distance to each node
    unordered_map<unsigned int, unordered_map<int, float>> pktDist;
   /*!
     * @param _index SFC index  @param _name SFC name
     * @param _numVNF number of VNFs except source and dest node
     * @param _trafficRate traffic arrival rate of SFC s
     */
    ServiceFunctionChain(unsigned int _index, const string& _name, unsigned int _numVNF, float _trafficArrivalRate){
        this->index = _index; this->name = _name;  this->numVNF = _numVNF;
        this->trafficArrivalRate = _trafficArrivalRate;
    }
    ~ServiceFunctionChain()  = default;

    unordered_map<int, vector<int>>& typeOfSFCAdj(int);
    static string typeOfSFCString(int) ;

    [[maybe_unused]] void showAdjList(int);
    [[maybe_unused]] void showOrigSFC_BlockWise_UsingAdj(int);
    [[maybe_unused]] void addAllDirectedEdgesToAdj(const vector<pair<int, int>> &);
    [[maybe_unused]] void ConvertSequentialSequenceToAdj();
    void convertToParallelSFC(const vector<vector<unsigned int>>& );

    [[maybe_unused]] void printOrigSFC(int, const string&);

    void showSFC_BlockWise(int, const unordered_map<unsigned int, unsigned int>&);
};

/*!
 * Show SFC (original or parallelised) in block wise order.\n
 * If instances are not mapped then first instance is considered.
 * @param type SFC type (Serial or parallel)
 */
 void ServiceFunctionChain::showSFC_BlockWise(const int type, const unordered_map<unsigned int, unsigned int>& VNFType2Inst){
    if(type == SFCseq){
        cout << "\nSeq SFC:["<<name<<"],nVNF["<< numVNF<<"]::\t";
        for(const int& vnfid : vnfSeq){
            if(vnfid == SFCsrc) cout << "(SRC -> ";
            else if(vnfid == SFCdst) cout << " DST)";
            else cout <<"f"<< vnfid <<char(96+VNFType2Inst.at(vnfid))<< "; -> ";
        }
    }
    else if(type == SFCpar){
        cout << "\nParallel Blocks SFC:["<<name<<"],nVNF["<< numVNF<<"]::\t";
        for(size_t i=0; i<vnfBlocksPar.size(); i++){
            if(i==0) cout << "(SRC -> ";
            cout<<"[ "; for(int vnfid: vnfBlocksPar[i]) cout<<"f"<<vnfid<<char(96+VNFType2Inst.at(vnfid))<<"; ";
            cout<<"] ->";
            if(i==vnfBlocksPar.size()-1)  cout << " DST)";
        }
    }
}






/*! Determine SFC type (Serial or parallel)
 * @param type SFC type (Serial or parallel)
 * @return sAdj, pAdj serial or parallel adjacency reference
 */
//unordered_map<int, vector<int>>& ServiceFunctionChain::typeOfSFCAdj(const int type){
//    return (type==SFCseq) ? sAdj : pAdj;
//}
/*! Determine SFC type (Serial or parallel)
 * @param type SFC type (Serial or parallel)
 * @return serial or parallel name string
 */
string ServiceFunctionChain::typeOfSFCString(const int type) {
    return (type==SFCseq) ? "Sequential" : "Parallel";
}

/*!
 * Show Adjacency list of SFC in console.
 * @param type SFC type (Serial or parallel)
 */
//void ServiceFunctionChain::showAdjList(const int type) {
//    unordered_map<int, vector<int>>& adj = typeOfSFCAdj(type);
//    cout << "\nAdjacency List of "<<name<<", Nodes["<< numVNF<<"] :: \n";
//    for (const auto& src: adj) {
//        if(SFCsrc == src.first)   cout << "SRC:";
//        else cout << "F["<< src.first << "]:";
//        for (auto dest : src.second){
//            if(SFCdst == dest)   cout << " -> " << "DST";
//            else cout << " -> F(" << dest <<")";
//        }
//        printf("\n");
//    }
//}


/*!
 * Show Adjacency list of SFC in console in Stage wise Order using Adjacency List.
 * @param type SFC type (Serial or parallel)
 */
void ServiceFunctionChain::showOrigSFC_BlockWise_UsingAdj(const int type) {
//    unordered_map<int, vector<int>>& adj = typeOfSFCAdj(type);
    unordered_map<int, vector<int>> adj;
    cout << "\nStage wise view of "<<name<<", Nodes["<< numVNF<<"] :: \n";

    unordered_map<int, bool> visitedBFS; ///< bfs visited array, size= number of vertexes
    unordered_map<int, int> level; ///<  level array, size= number of vertexes
    unordered_map<int, string> rankSame;///<  array to keep same rank node in one level, size= number of vertexes
    int max_level=1;

    vector<string> rankPush;
    queue<int> q; q.push(SFCsrc);
    visitedBFS[SFCsrc]=true; level[SFCsrc]=0;
    while(!q.empty()){
        unsigned int sz = q.size();
        while(sz--){
            int u_val = q.front();  q.pop();
            for(auto v_val: adj[u_val]){
                if(v_val == SFCdst)continue;
                if(!visitedBFS[v_val]){
                    visitedBFS[v_val]=true;
                    q.push(v_val);
                    level[v_val]=level[u_val]+1;
                    if(level[v_val]+1>max_level)max_level++;
                    rankSame[level[v_val]].append("f"+to_string(v_val)).append("; ");
                }
            }//for adj
        }
    }

    cout << " SRC -> ";
    //printing node with same level together
    for(int i=0; i<max_level; i++){
        if(!rankSame[i].empty())
            cout << "[ " << rankSame[i] << "] -> ";
    }//for i
    cout << " DST";
}

/*!
 * Print SFC Graph using Graphviz in a PNG format
 * @param type SFC type (Serial or parallel)
 */
void ServiceFunctionChain::printOrigSFC(const int type, const string& testDirName){
//    unordered_map<int, vector<int>>& adj = typeOfSFCAdj(type);
    unordered_map<int, vector<int>> adj;
    string sfcTypeName = typeOfSFCString(type);
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
    fout << SFCsrc << " [label = \"src\", fillcolor="<<srcDstFillColor<<" fontcolor="<<srcDstFontColor<<"]"<< endl;
    fout << SFCdst << " [label = \"dst\", fillcolor="<<srcDstFillColor<<" fontcolor="<<srcDstFontColor<<"]"<< endl;

    int src = SFCsrc;
    qLoGraph.push(src);
    visitedBFS[src]=true;

    level[src]=0;//0
    rankSame[level[src]].append(to_string(src)).append("; ");

    while (!qLoGraph.empty())		// Level order Traversal
    {
        u_val = qLoGraph.front();  qLoGraph.pop();
        if(u_val != SFCsrc and u_val != SFCdst)
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


/*! Explicit Function to Convert the Sequential VNFs id from (vnfSeq) into Adj List variable (sAdj). */
//void ServiceFunctionChain::ConvertSequentialSequenceToAdj(){
//    size_t len = vnfSeq.size(); ///< length of SFC including src and destination
//    for(size_t i=0; i<=len-2; i++ ){ // len-1 is last stage
//        int srcVNFid = vnfSeq[i];
//        int dstVNFid = vnfSeq[i+1];
//        sAdj[srcVNFid].push_back(dstVNFid);
//    }
//}


/*! Convert the Parallel VNFsID in stages from (vnfBlocksPar) into Adj List variable (pAdj).
 * @param blocks stage wise parallel blocks */
//void ServiceFunctionChain::convertToParallelSFC(const vector<vector<unsigned int>>& blocks){
//    this->vnfBlocksPar = std::move(blocks);
//    size_t stages = vnfBlocksPar.size(); ///< stages of SFC.
//    /// for each stage. stages-1 is the last stage
//    for(size_t s=0; s<=stages-2; s++ ){
//        ///for each vnf in source stage
//        for(int srcVNFid : vnfBlocksPar[s]){
//            /// for each vnf in next stage
//            for(int dstVNFid : vnfBlocksPar[s+1]){
//                pAdj[srcVNFid].push_back(dstVNFid);
//            }
//        }
//
//    }
//}

/*!Another Function to make Serial SFC given all Edges. It Add all the directed edges of sfc int Adj. \n
 * it read VNFs data from edges and convert to sequence to array.
 * @param allEdges all the edges {u,v}, including src and dst edges also.
 */
//void ServiceFunctionChain::addAllDirectedEdgesToAdj(const vector<pair<int, int>>& allEdges) {
//    unordered_map<int, vector<int>>& adj = sAdj;
//    for(const auto& edge: allEdges){
//        adj[edge.first].push_back(edge.second);
////        isVNF_Present[edge.second]= true;
//        I_VNFType2Inst[edge.second]= 0;
//    }
//    vnfSeq.push_back(SFCsrc);
//    queue<int> q; q.push(SFCsrc);
//    while(!q.empty()) {
//        int src = q.front(); q.pop();
//        for(const auto& w: adj[src]){
//            vnfSeq.push_back(w);
//            q.push(w);
//        }
//    }
////    cout<<endl<<"test: "; for(auto & x: vnfSeq) cout<<x<<"--";
//}

#endif //SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H
