// Created by vijay on 22-02-2023.
//
#ifndef SFC_PARALLELIZATION_PHYSICALGRAPH_H
#define SFC_PARALLELIZATION_PHYSICALGRAPH_H

/*!
 * Store Resources Values, number of physical+virtual cores,  memory in GB, disk space, CPU speed
 * @tparam type_res resource data type. default=unsigned int.
 */
template <class type_res>
class NodeCapacity
{
public:
    unsigned int cores{}; ///<number of physical+virtual cores
    type_res memory,disk,cpuSpeed; ///< memory in GB, disk space, CPU speed
    /*!
     * @param _cores number of physical+virtual cores
     * @param _memory default 0. in GB.
     * @param _disk default 0.
     * @param _cpuSpeed default 0.
     */
    explicit NodeCapacity(unsigned int _cores, type_res _memory=0, type_res _disk=0, type_res _cpuSpeed=0){
        this->cores = _cores;
        this->memory = _memory; this->disk = _disk; this->cpuSpeed = _cpuSpeed;
    }
    ~NodeCapacity() = default;
};

/*!
 * @tparam type_res resource data type.
 */
template <class type_res >
class PhysicalNode{
public:
    /// \brief physical node index. Type int.
    unsigned int index{}; string name{}; ///< physical node name
    NodeCapacity<type_res> capacity; /*!< capacity of the node. Type NodeCapacity. */
    unordered_set<unsigned int> pn2VMs; /*!< list of VM index present in the physical machine. */
    /*!
     * @param _index physical node index (int)
     * @param _name node name
     * @param _givenCapacity capacity of the node
     */
    PhysicalNode(unsigned int& _index, const string& _name, NodeCapacity<type_res> _givenCapacity): capacity(_givenCapacity){
        this->index = std::move(_index);
        this->name = std::move(_name);
    }
    ~PhysicalNode() = default;
};

/*!
 * @tparam type_wgt edge weight data type.
 */
template <class type_wgt = unsigned int>
struct PhysicalEdge
{
    /*! @param u edge source physical node index,  @param v edge destination physical node index  */
    unsigned int u, v;   type_wgt wt;///<edge(u,v) weight
    /*! Physical Edge Constructor
     * @param  _source soruce node
     * @param _destination destination node
     * @param _weight distance between _src to _destination
     */
    PhysicalEdge(unsigned int _source, unsigned int _destination, type_wgt _weight){
        this->u = _source; this->v = _destination;
        this->wt = _weight;
    }
    ~PhysicalEdge() = default;
};

/*!
 * Physical Graph Class for Network Creation
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
class PhysicalGraph {
    unsigned int numV, numE; /*!< numV=number of vertices, numE=number of edges  */
    unsigned int srcV; ///<source vertex. default 1, loop from 1 to <= numVertexes
public:
    unordered_map<unsigned int, PhysicalNode<type_res>*> PNode; ///< index to Physical Node address
    vector<PhysicalEdge<type_wgt> *> *adj; ///< adjacency list[u] = {v1,wt1}->{v2,wt2}
    vector<vector<type_wgt>> mat /*!< adjacency Matrix original distances. in Meters. max distance 2*Radius of Earth = 12756000 meters,*/;
    vector<vector<type_wgt>> dist; ///<distance Matrix calculated using all pairs shortest path
    vector<vector<unsigned int>> nextHop; ///< nextHop of distance Matrix
    type_wgt TypeMaxValue = std::numeric_limits<type_wgt>::max(); /*!< max Value of the data type. \n float: 3.40282e+38 or 0x1.fffffep+127 \n size_t: 18446744073709551615 or 0xffffffffffffffff.*/
    type_wgt EPS = std::numeric_limits<type_wgt>::epsilon();/*!< EPS Returns the machine epsilon, that is, the difference between 1.0 and the next value representable by the floating-point type T. \n It is only meaningful if std::numeric_limits<T>::is_integer == false. \n double	= DBL_EPSILON, float = FLT_EPSILON, unsigned int = 0 */

    /*! Physical Graph Constructor
     * @param _numVertexes Number of vertexes in network (1-based Indexing)
     * @param _numEdges Total Edges in network
     * @param _sourceVertex default 1, loop from 1 to <= numVertexes
     */
    PhysicalGraph(unsigned int _numVertexes, unsigned int _numEdges, unsigned int _sourceVertex = 1) {
        srcV = _sourceVertex; numV = _numVertexes; numE = _numEdges;
        adj = new vector<PhysicalEdge<type_wgt> *>[numV + 1];
        mat = vector<vector<type_wgt>>(numV + 1, vector<type_wgt>(numV + 1, 0));
        dist = vector<vector<type_wgt>>(numV + 1, vector<type_wgt>(numV + 1, 0));
        nextHop = vector<vector<unsigned int>>(numV + 1, vector<unsigned int>(numV + 1, 0));
    }

    ~PhysicalGraph() {
        for (unsigned int v = srcV; v <= numV; ++v) { ///deletion of Adj List Edges
            for (const auto edge: adj[v]) delete edge; /// delete all edges
            delete PNode.at(v); /// delete nodes
        }
//        for(const auto pnInfo: PNode){
//            for (const auto edge: adj[pnInfo.second->index]) delete edge; /// delete all edges
//            delete pnInfo.second;
//        }
        delete[] adj; // delete entire adj list
        if (debug) cout << "\n[PhysicalGraph Destructor Completed for G(" << numV<<", "<<numE<<")]";
    }

    [[maybe_unused]] void showAdjList() const;

    [[maybe_unused]] void showAdjMatrix() const;

    [[maybe_unused]] void showAllPairsShortestPath() const;

    [[maybe_unused]] void printPhysicalGraph(const string&) const;

    void setNumOfEdges(const unsigned int val){
        numE = val;
    }

//    void setNumOfNodes(const unsigned int val){
//        numV = val;
//    }
    void addEdge(PhysicalEdge<type_wgt> *);
    void calcAllPairsShortestPath();

    [[maybe_unused]] void showPNs_Description();
};

/*!
 * Show Network Adjacency List in console.
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
void PhysicalGraph<type_wgt,type_res>::showAdjList()const {
    cout << "\n\n Adjacency List G(" << numV<<", "<<numE<<") ::";
    for (unsigned int v = srcV; v <= numV; ++v) {
        cout << "\n V["<< v << "]:";
        for (const auto edge : adj[v])
            cout << " -> (" << edge->v<<", "<<edge->wt<<")";
    }
}

/*!
 * Show Network Adjacency Matrix in console.
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
void PhysicalGraph<type_wgt,type_res>::showAdjMatrix() const{
    cout << "\n\n ----- Adjacency Matrix of G(" << numV<<", "<<numE<<") ::";
    cout<<"\nV:\t";  for (unsigned int v = srcV; v <= numV; ++v)  cout <<v<< "\t"; cout<<endl;
    for (unsigned int v = 0; v <= numV; ++v) cout <<"---"<< "\t";
    for (unsigned int u = srcV; u <= numV; ++u){
        cout<<"\n"<<u<<"|\t";
        for (unsigned int v = srcV; v <= numV; ++v)
            cout<<mat[u][v]<<"\t";
    }
}

/*! Show all the Nodes and their description
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
void PhysicalGraph<type_wgt,type_res>::showPNs_Description(){
    cout << "\n\n ----- Nodes Description of G(" << numV<<", "<<numE<<") ::";
    cout<<"\nIndex\t"<<"Name\t\t"<<" CAP[cores][memory][disk][speed]\t"<<"VirtualMachines";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (unsigned int u = srcV; u <= numV; ++u){
        const PhysicalNode<type_res> * pn = PNode.at(u);
        cout<<"\n"<<pn->index<<" |\t";
        cout<<pn->name<<" |\t";
        cout<<pn->capacity.cores<<"\t"<<pn->capacity.memory<<"\t"<<pn->capacity.disk<<"\t"<<pn->capacity.cpuSpeed<<"\t";
        for(const unsigned int& vmid: pn->pn2VMs)
            cout<<"VM["<<vmid<<"] ";
    }
}
/*!
 * Show Network Adjacency Matrix in console After calculating All Pairs Shortest Path.
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
void PhysicalGraph<type_wgt,type_res>::showAllPairsShortestPath() const{
    cout << "\n\n All Pairs Shortest Path G(" << numV<<", "<<numE<<") ::";
    cout<<"\nV:\t";  for (unsigned int v = srcV; v <= numV; ++v)  cout <<v<< "\t"; cout<<endl;
    for (unsigned int v = 0; v <= numV; ++v) cout <<"---"<< "\t";
    for (unsigned int u = srcV; u <= numV; ++u){
        cout<<"\n"<<u<<"|\t";
        for (unsigned int v = srcV; v <= numV; ++v)
        {
            if(dist[u][v] == TypeMaxValue) cout<<"inf\t";
            else cout<<dist[u][v]<<"\t";
        }
    }
    cout<<endl;
    cout << "\n\t Next Hop G(" << numV<<", "<<numE<<") ::";
    cout<<"\nV:\t";  for (unsigned int v = srcV; v <= numV; ++v)  cout <<v<< "\t"; cout<<endl;
    for (unsigned int v = 0; v <= numV; ++v) cout <<"---"<< "\t";
    for (unsigned int u = srcV; u <= numV; ++u){
        cout<<"\n"<<u<<"|\t";
        for (unsigned int v = srcV; v <= numV; ++v) {
            if(nextHop[u][v] == TypeMaxValue) cout<<"inf\t";
            else cout<<nextHop[u][v]<<"\t";
        }
    }
}

/*!
 * Print Network Graph using Graphviz in a PNG format
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
[[maybe_unused]] void PhysicalGraph<type_wgt,type_res>::printPhysicalGraph(const string& testDirName) const{
    if (mat.empty()) {
        cout << "Network Graph is Empty."<<endl;
        return;
    }

    ofstream fout;
    string fileName = "NetGraph_G_"+to_string(numV)+"_"+to_string(numE); ///< (Network Graph + V + E) name without space
    string filepath = output_directory+testDirName+diagram_directory+fileName; ///<path to .gv file without extention
    string filepathExt = filepath+".gv"; ///<filepath with extention of .gv
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    if(debug)cout<<"\n[Function Running: "<<__FUNCTION__<<"] File:"<<filepathExt<<endl;
/**********************/
    queue<unsigned int> qLoGraph;		// queue level order Graph
    vector<bool> visitedNodes(numV+1, false); // bfs visited array, size= number of vertexes + 1 (1 based indexing)
    vector<int> level(numV+1); // level array, size= number of vertexes
    vector<string> rankSame(numV+1);// array to keep same rank node in one level, size= number of vertexes
    vector<vector<bool>> visitedEdges(numV+1, vector<bool>(numV+1, false));

    int max_level=1;
    unsigned int u_val, v_val;
    const string& nodeColor = "darkslateblue", nodeFontColor = "firebrick4";
    const string& edgeColor = "darkslateblue", edgeFontColor = "firebrick4";
/**********************/
    fout << "graph "<<fileName<<" {" << endl;
    fout << "  label = \"Network Graph G ("<<numV<<","<<numE<<")\" \n";
    fout << "  ranksep=\"equally\"; \n"
            "  rankdir=LR; " << endl << endl;
    fout << "node [fixedsize=shape margin=0.3 width=0.6 shape=circle color="<<nodeColor<<" fontcolor="<<nodeFontColor<<"] "<<endl;
    fout << "edge [color="<<edgeColor<<" fontcolor="<<edgeFontColor<<" fontsize=12]"<<endl;

    unsigned int src = srcV;
    qLoGraph.push(src);
    visitedNodes[src]=true;
    level[src]=0;
    rankSame[level[src]].append(to_string(src)).append("; ");

    while (!qLoGraph.empty())		// Level order Traversal
    {
        u_val = qLoGraph.front();  qLoGraph.pop();
        fout << u_val << " [label = <&eta;<sub>"<<u_val<<"</sub>>];"<< endl;
        for(const auto edge: adj[u_val]){
            v_val = edge->v;        //adj node value
            if(!visitedNodes[v_val] ){
                visitedNodes[v_val]=true;
                qLoGraph.push(v_val);

                level[v_val]=level[u_val]+1;
                if(level[v_val]+1>max_level)max_level++;
                rankSame[level[v_val]].append(to_string(v_val)).append("; ");
            }
            if(!visitedEdges[u_val][v_val]) // add edges for graphviz
            {
                visitedEdges[u_val][v_val] = visitedEdges[v_val][u_val] = true;
                fout <<"\t"<< u_val << " -- " << v_val <<" [label ="<<mat[u_val][v_val]<<"]"<< endl;
            }
        }//for adj
    }//while qLoGraph

/*  // Printing Edges and Nodes using Adjacency Matrix
    for(u_val = srcV; u_val <= numV; u_val++)
    {
        fout << u_val << " [label = <&eta;<sub>"<<u_val<<"</sub>>];"<< endl;
        for (v_val = u_val+1; v_val<=numV ; ++v_val) {
            if(mat[u_val][v_val] != (float)0){
                fout <<"\t"<< u_val << " -- " << v_val <<" [label ="<<mat[u_val][v_val]<<"]"<< endl;
            }
        }
        fout << endl ;
    }
*/

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
    if(debug)cout<<"\n[Function Completed: "<<__FUNCTION__<<"]"<<" Graphviz: Original Network GraphG ("<<numV<<","<<numE<<") printed. File:"<<filepath<<".png \n";
    if(!debug)remove(filepathExt.c_str());
}

/*!
 * Add undirected edge to graph adjacency list and matrix.
 * @param edge pointer of structure Physical Edge.
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
void PhysicalGraph<type_wgt,type_res>::addEdge(PhysicalEdge<type_wgt>* edge) {
    adj[edge->u].emplace_back(edge);
    adj[edge->v].emplace_back(new PhysicalEdge<type_wgt>(edge->v, edge->u, edge->wt));
    mat[edge->u][edge->v] =  mat[edge->v][edge->u] = edge->wt;
}

/*!
 * Calculates All Pairs Shortest Path using adjacency matrix
 * @param mat adjacency matrix
 * @param dist updates distance matrix and @param nextHop updates matrix to take next vertex for shortest path.
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
template <class type_wgt, class type_res>
void PhysicalGraph<type_wgt,type_res>::calcAllPairsShortestPath()
{
    //preprocessing
//    cout<<TypeMaxValue<<"------------"<<EPS<<endl;
    for(unsigned int i=srcV; i<=numV; i++) // source
    {
        for(unsigned int j=srcV; j<=numV; j++) // destination
        {
            if(i==j){
                dist[i][j] = 0;
                nextHop[i][j] = j;
            }
            else if(mat[i][j] != 0){ //mat[0][0] is zero 0.
                dist[i][j] = mat[i][j];
                nextHop[i][j] = j;
            }
            else
                dist[i][j] = TypeMaxValue;
        }
    }
    // calculation
    for(unsigned int k=srcV; k<=numV; k++) // intermediate vertices
    {
        for(unsigned int i=srcV; i<=numV; i++) // source
        {
            for(unsigned int j=srcV; j<=numV; j++) // destination
            {
                if(i==k || j==k || dist[i][k]==TypeMaxValue  || dist[k][j]==TypeMaxValue || i==j)
                    continue;
//                if(dist[i][j] > dist[i][k] + dist[k][j])  {
//                    dist[i][j] = dist[i][k] + dist[k][j];
//                    nextHop[i][j] = nextHop[i][k];
//                }
                if (dist[i][j] - EPS > dist[i][k] + dist[k][j]) ///EPS - machine epsilon
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    nextHop[i][j] = nextHop[i][k];
                }
            }
        }
    }
    if(debug)cout<<"\n\t[All Pairs Shortest Path Completed]";
}


#endif //SFC_PARALLELIZATION_PHYSICALGRAPH_H
