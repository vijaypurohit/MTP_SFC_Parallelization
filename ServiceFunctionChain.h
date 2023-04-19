//
// Created by vijay on 23-02-2023.
//

#ifndef SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H
#define SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H

/*! @brief For a single SFC, according to the algorithm what is the optimal/best (minimum delay) parameters we have found.*/
struct sfcResult{
    int seq_pid{noResSeq}, ppar_pid{noResPar}, fullpar_pid{noResPar}; ///< idx of allPartParSFC Array for sequential/part parallel/full parallel for which algorithm give optimal answer.
    type_delay seq_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of sequential length chain according to our algorithm.
    type_delay ppar_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of partial parallel chain according to our algorithm.
    type_delay fullpar_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of full parallel chain according to our algorithm.
    unordered_map<unsigned int, unsigned int> seq_fninstmap, ppar_fninstmap, fullpar_fninstmap; ///< Best mapping {fun->its instance taken} for given sequential length chain according to our algorithm.
    type_delay seq_duration, ppar_duration, fullpar_duration;
};
struct finalResults{
    std::string name_sol{}; //< name of the solution
    unordered_map<unsigned int, sfcResult> solobj; //< sfc index to its solution values
    double seq_duration{}, ppar_duration{}, fullpar_duration{}; //< time taken to construct the solution with/without parallelism
    explicit finalResults(string name){
        name_sol = std::move(name);
    }
};
finalResults res_layerg(name_sol_layerg), res_partial(name_sol_partial);
/* ******** Some Structers *************** */
/*! @brief Before finding best mapping of SFC, this will save delays based on VNFs in SFC to pre-calculate in order to avoid multiple computations of same VNF.*/
struct vnfDelaysPreComputed{
    type_delay exeDelay{}; ///< execution delay of the VNF
    type_delay prcDelay{}; ///< processing delay of the VNF
    unordered_map<unsigned int, type_delay> queuingDelay; ///< queuing delay of the VNF and its instance
//    unordered_map<unsigned int, type_delay> seq_queuingDelay; ///< queuing delay of the VNF and its instance, in case of sequential mapping and utilization
//    unordered_map<unsigned int, type_delay> par_queuingDelay; ///< queuing delay of the VNF and its instance, in case of parallel mapping and utilization
};

/*! @brief minheap priority queue node to find minimum dist and minimum utilization path from source to destination/current stage.*/
struct pqNode{
    unsigned int x{}, y{}; ///< source and destination pair of previous and current stage lgNode
    type_delay mindist{}; ///< minimum delay of the path from source to current processing node node
    type_delay utilization{}; ///< max utilization of the path from source to current processing node node
    vector<unsigned int> path; ///< path constructed till now from source
    pqNode()=default;
    pqNode(unsigned int givenlgSrcId, unsigned int givenlgDstId, type_delay givenDist, type_delay givenUtilization, std::vector<unsigned int> givenPath):x(givenlgSrcId),y(givenlgDstId),mindist(givenDist),utilization(givenUtilization){
        path = std::move(givenPath);
    }
//    template<class T>
//    typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
//    almost_equal(T val_x, T val_y, int ulp) const
//    {
//        // the machine epsilon has to be scaled to the magnitude of the values used and multiplied by the desired precision in ULPs (units in the last place)
//        return std::fabs(val_x - val_y) <= std::numeric_limits<T>::epsilon() * std::fabs(val_x + val_y) * ulp
//               // unless the result is subnormal
//               || std::fabs(val_x - val_y) < std::numeric_limits<T>::min();
//    }
    template<class T>
    bool approximatelyEqual(T a, T b, T epsilon) const
    {
        return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
        return (static_cast<int>(a*100.0) == static_cast<int>(b*100.0));
    }

    bool operator<(const struct pqNode& other) const { // overloaded operator for priority queue
//        if(mindist == other.mindist) ///< if distance is same
        if(approximatelyEqual<type_delay>(mindist, other.mindist, 0.0005)) /// upto 2nd decimal digit equal, 38.941, 38.94
        {
            if(utilization>0 or other.utilization>0) /// if some utilizatio to compare
                return utilization > other.utilization; //min heap, return pair of x-y with min utilization
            if(x == other.x)
                return y > other.y; /// sort according to instance id. first come first serve
            return x > other.x; /// sort according to instance id. first come, first served
//            return false;
        }
        else {
            return mindist > other.mindist; //min heap, return pair of x-y with minimum distance.
        }
    }
};

/*! @brief Layer Graph Node Vertex
 * @param idx index to detect node uniquely
 * @param instCombination {fnType, instId} pairs showing instance combination at this node
 * @param cntPN count of physical node,frequency in the instance combination
 * @param exePN maximum execution time of physical node in the instance combination
 * @param utilization utilisation percentage of all the instances present in the lgNode
 * @param children next stage lgNode index and its distance, that is pair of this->node = {next stg node, min dist}.
 * @param kpaths number of shortest path traverse through this lgNode
 */
struct lgNode{
    unsigned int idx{}; ///< index to detect node uniquely
    vector<pair<unsigned int,unsigned int>> instCombination; ///< {fnType, instId} pairs showing instance combination at this node
    unordered_map<unsigned int, unsigned int> cntPN; ///< count of physical node,frequency in the instance combination
    unordered_map<unsigned int, type_delay> exePN; ///< maximum execution time of physical node in the instance combination
    type_delay utilization{0}; ///< utilisation percentage of all the instances present in the lgNode
    vector<pair<unsigned int, type_delay>> children; ///<  next stage lgNode index and its distance, that is pair of this->node = {next stg node, min dist}.
    vector<pqNode> kpaths; ///< number of shortest path traverse through this lgNode

    lgNode()=default;
    explicit lgNode(unsigned int index):idx(index){};
    lgNode(unsigned int index, const vector<pair<unsigned int,unsigned int>>& givenIC):idx(index), instCombination(givenIC){ }
};

/*! @brief Service Function Class */
class ServiceFunctionChain
{
public:
    unsigned int index /*!< SFC index */; string name;  /*!< SFC name */ unsigned int numVNF;  /*!< number of VNFs except source and dest node */
    type_delay trafficArrivalRate;  /*!< traffic arrival rate of SFC s. in packets per second. arrival rate < service rate. */
    vector<unsigned int> vnfSeq; ///< vnfids in sequential manner excluding src (SFCsrc=0) and dest (SFCdst=-10).

// parameters found through algorithms.
    vector<vector<unsigned int>> vnfBlocksPar; ///< vnfids in stage wise, where stg i denotes parallel vnf in ith stage. Exvept src and dest stg.

    vector<vector<vector<unsigned int>>> allPartParSFC; ///< All the partial paralle VNF blocks of the given sequential SFC.

   /*!
     * @param _index SFC index  @param _name SFC name
     * @param _numVNF number of VNFs except source and dest node
     * @param _trafficRate traffic arrival rate of SFC s
     */
    ServiceFunctionChain(unsigned int _index, const string& _name, unsigned int _numVNF, type_delay _trafficArrivalRate){
        this->index = _index; this->name = _name;  this->numVNF = _numVNF;
        this->trafficArrivalRate = _trafficArrivalRate;
    }
    ~ServiceFunctionChain()  = default;

    void showSequentialSFC(const unordered_map<unsigned int, unsigned int>&);
    void showParallelSFC(const vector<vector<unsigned int>>&, const unordered_map<unsigned int, unsigned int>&);

    [[maybe_unused]] void printOrigSFC(int, const string&);
};

/*!
 * @brief It is a comparator function used in max heap to sort the SFCs according to the traffic arrival rate in descending order.
 * If traffic rate are same then sort according to the descending order of the length.
 */
struct comparator_sfc {
    bool operator()(const ServiceFunctionChain  *const sfc1, const ServiceFunctionChain *const sfc2) {
        if(sfc1->trafficArrivalRate == sfc2->trafficArrivalRate){
            return sfc1->numVNF < sfc2->numVNF;
        }
        return sfc1->trafficArrivalRate < sfc2->trafficArrivalRate; /// max heap,
    }
};
/*!
 * Show Sequential SFC in console.
 * @param VNFType2Inst function to instance mapping.
 */
void ServiceFunctionChain::showSequentialSFC(const unordered_map<unsigned int, unsigned int>& VNFType2Inst = unordered_map<unsigned int, unsigned int>()){

    cout<<"\nSeq. SFC:"<<index<<" | TrafficRate: "<<trafficArrivalRate<<" | cntVNFs: "<<numVNF<<" | PartialChains: "<<allPartParSFC.size()<<"\n\t";
    cout << "(SRC -> ";
    for(const auto& fn : vnfSeq){
        cout <<"f"<< fn ;
        if(!VNFType2Inst.empty())cout<<char(96+VNFType2Inst.at(fn));
        cout<< "; -> ";
    } cout << " DST)";
}

/*!
 * Show Parallel VNF SFC in console.
 * @param givenSFC given SFC
 * @param VNFType2Inst function to instance mapping.
 */
void ServiceFunctionChain::showParallelSFC(const vector<vector<unsigned int>>& givenSFC, const unordered_map<unsigned int, unsigned int>& VNFType2Inst = unordered_map<unsigned int, unsigned int>()){
    cout<<"\nPar. SFC:"<<index<<" | TrafficRate: "<<trafficArrivalRate<<" | cntVNFs: "<<numVNF<<" | PartialChains: "<<allPartParSFC.size()<<"\n\t";
    cout << "(SRC ->";
    for(const auto& blk: givenSFC){
        cout<<" ["; for(int fn: blk){
                    cout <<"f"<< fn ;
                    if(!VNFType2Inst.empty())cout<<char(96+VNFType2Inst.at(fn));
                    cout<<"; ";
                }
        cout<<"] ->";
    } cout << " DST)";
}

/*!
 * Show description of all the SFCs in the Network
 */
void showSFCsDescriptions(vector<ServiceFunctionChain *> allSFCs){

    cout << "\n\n ----- Service Functions Chains Description ::";
    cout<<"\nIndex\t"<<"Name\t"<<"cntVNFs\t"<<"TrafficRate\t"<<"PartialChains\t"<<"Sequential SFC\t"<<"Parallel SFC";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (const auto &sfc: allSFCs){
        cout<<"\n"<<sfc->index<<" | ";
        cout<<sfc->name<<" | ";
        cout<<sfc->numVNF<<" | ";
        cout<<sfc->trafficArrivalRate<<" | ";
        cout<<sfc->allPartParSFC.size()<<" | ";
        cout << "(SRC->";
        for(const auto& fn : sfc->vnfSeq){ cout <<"f"<< fn ;   cout<< ";->";
        } cout << " DST)\t";
        cout << "(SRC->";
        for(const auto& blk: sfc->allPartParSFC[sfc->allPartParSFC.size()-1]){
            cout<<" ["; for(int fn: blk){
                cout <<"f"<< fn ;  cout<<"; ";
            }
            cout<<"]->";
        } cout << "DST)";
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



#endif //SFC_PARALLELIZATION_SERVICEFUNCTIONCHAIN_H
