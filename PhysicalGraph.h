// Created by vijay on 22-02-2023.
//
#ifndef SFC_PARALLELIZATION_PHYSICALGRAPH_H
#define SFC_PARALLELIZATION_PHYSICALGRAPH_H

/*! Store Resources Values, number of physical+virtual cores,  memory in GB, disk space, CPU speed */
class NodeCapacity
{ public:
    unsigned int cores{0}; ///<number of physical+virtual cores
//    type_res memory{},disk{},cpuSpeed{}; ///< memory in GB, disk space, CPU speed
    /*!  @param _cores number of physical+virtual cores */
    explicit NodeCapacity(unsigned int _cores=0):cores(_cores){ }
    ~NodeCapacity() = default;
};

/*!  Hold Node Data */
class PhysicalNode{
public:
    /// \brief physical node index. Type int.
    unsigned int index{}; string name{}; ///< physical node name
    NodeCapacity capacity; /*!< capacity of the node. Type NodeCapacity. */
    unordered_set<unsigned int> neighbours; ///< neighbours of the node
    type_delay clusteringcoeff_local; ///< local clustering coefficient.
    /*!
     * @param _index physical node index (int)
     * @param _name node name
     * @param _givenCapacity capacity of the node
     */
    PhysicalNode(unsigned int& _index, const string& _name, NodeCapacity _givenCapacity): capacity(_givenCapacity){
        this->index = _index;
        this->name = std::move(_name);
    }
    PhysicalNode() = default;
    ~PhysicalNode() = default;
};


/*!
 * Physical Graph Class for Network Creation
 */
class PhysicalGraph {
public:
    string filename_network;
    unsigned int numV, numE; /*!< numV=number of vertices, numE=number of edges  */
    unsigned int srcV; ///<source vertex. default 1, loop from 1 to <= numVertexes
    unordered_map<unsigned int, PhysicalNode> PNode; ///< index to Physical Node/Server address

    vector<vector<type_wgt>> mat /*!< adjacency Matrix original distances. in Meters. max distance 2*Radius of Earth = 12756000 meters,*/;
    vector<vector<type_wgt>> dist; ///<distance Matrix calculated using all pairs shortest path
    vector<vector<unsigned int>> nextHop; ///< nextHop of distance Matrix
    type_wgt TypeMaxValue = std::numeric_limits<type_wgt>::max(); /*!< max Value of the data type. \n float: 3.40282e+38 or 0x1.fffffep+127 \n size_t: 18446744073709551615 or 0xffffffffffffffff.*/
    type_wgt EPS = std::numeric_limits<type_wgt>::epsilon();/*!< EPS Returns the machine epsilon, that is, the difference between 1.0 and the next value representable by the floating-point type T. \n It is only meaningful if std::numeric_limits<T>::is_integer == false. \n double	= DBL_EPSILON, float = FLT_EPSILON, unsigned int = 0 */

    unsigned long sum_cores = 0; ///< total cores in the entire network
    double clusteringcoeff_mean = 0; ///< mean clustering coefficient of the graph.
    double clusteringcoeff_global = 0; ///< global clustering coefficient of the graph.

    /*! Physical Graph Constructor
     * @param _numVertexes Number of vertexes in network (1-based Indexing)
     * @param _numEdges Total Edges in network
     * @param _sourceVertex default 1, loop from 1 to <= numVertexes
     */
    PhysicalGraph(unsigned int _numVertexes, unsigned int _numEdges) {
        srcV = 1; numV = _numVertexes; numE = _numEdges;
        mat = vector<vector<type_wgt>>(numV + 1, vector<type_wgt>(numV + 1, 0));
        dist = vector<vector<type_wgt>>(numV + 1, vector<type_wgt>(numV + 1, 0));
        nextHop = vector<vector<unsigned int>>(numV + 1, vector<unsigned int>(numV + 1, 0));
    }
    PhysicalGraph() = default;
    ~PhysicalGraph() = default;

    [[maybe_unused]] void showPNs_Description() const;
    [[maybe_unused]] void showAdjMatrix() const;
    [[maybe_unused]] void showAllPairsShortestPath() const;

    void setNumOfEdges(const unsigned int& val){ this->numE = val; }
    void calcAllPairsShortestPath();
    vector<unsigned int> constructShortestPath(const unsigned int &FromSrc, const unsigned int &ToDst);
    bool calcClusteringCoefficient(bool);
};


/*!  Show Network Adjacency Matrix in console.*/
void PhysicalGraph::showAdjMatrix() const{
    cout << "\n\n ----- Adjacency Matrix of G(" << numV<<", "<<numE<<") :: (in meters)";
    cout<<std::setw(2)<<"\nV:\t";  for (unsigned int v = srcV; v <= numV; ++v)  cout <<std::setw(7)<<v<< "\t"; cout<<endl;
    for (unsigned int v = 0; v <= numV; ++v) cout <<std::setw(7)<<"------"<< "\t";
    for (unsigned int u = srcV; u <= numV; ++u){
        cout<<"\n"<<std::setw(2)<<u<<"|\t";
        for (unsigned int v = srcV; v <= numV; ++v)
            cout<<std::setw(7)<<mat[u][v]<<"\t";
    }
}

/*! Show all the Nodes and their description*/
void PhysicalGraph::showPNs_Description() const{
    cout << "\n\n ----- Nodes Description of G(" << numV<<", "<<numE<<") ::";
    cout<<"\nIdx |\t"<<"Name | "<<"Cores | "<<"Deg"<<" | "<<"CoeffL";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----;";
    for (unsigned int u = srcV; u <= numV; ++u){
//        if(PNode.find(u) == PNode.end()) continue;
        const PhysicalNode& pn = PNode.at(u);
        cout<<"\n"<<std::setw(2)<<pn.index<<" | ";
        cout<<std::setw(5)<<pn.name<<" | ";
        cout<<std::setw(2)<<pn.capacity.cores<<" | ";
        cout<<std::setw(2)<<pn.neighbours.size()<<" | ";
        cout<<std::setw(5)<<pn.clusteringcoeff_local<<" | ";
    }
}

/*!  Show Network Adjacency Matrix in console After calculating All Pairs Shortest Path. */
void PhysicalGraph::showAllPairsShortestPath() const{
    cout << "\n\n All Pairs Shortest Path G(" << numV<<", "<<numE<<") :: (in meters)";
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
 * Calculates All Pairs Shortest Path using adjacency matrix
 * @param mat adjacency matrix
 * @param dist updates distance matrix and @param nextHop updates matrix to take next vertex for shortest path.
 * @tparam type_wgt edge weight data type.
 * @tparam type_res resource data type.
 */
void PhysicalGraph::calcAllPairsShortestPath()
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
    if(debug)cout<<"\n\t[Calculated All Pairs Shortest Path]";
}

vector<unsigned int> PhysicalGraph::constructShortestPath(const unsigned int& FromSrc, const unsigned int& ToDst)
{
    if (nextHop[FromSrc][ToDst] == 0) /// If there's no path between node src and dst, simply return an empty array
        return {};
    unsigned int src = FromSrc;
    vector<unsigned int> path = { src };/// Storing the path in a vector
    while (src != ToDst) {
        src = nextHop[src][ToDst];
        path.push_back(src);
    }
    return path;
}

/*! Calculating Clustering Coefficient -> local, mean, gloabal coefficient (refer wikipedia)
 * The local clustering coefficient of a vertex (node) in a graph quantifies how close its neighbours are to being a clique (complete graph). \n
 *  Duncan J. Watts and Steven Strogatz introduced the measure in 1998 to determine whether a graph is a small-world network. \n
 *  A graph G=(V,E) formally consists of a set of vertices V and a set of edges E between them. An edge e_{ij} connects vertex v_{i} with vertex v_{j}. \n
 *  The neighbourhood N_{i} for a vertex v_{i} is defined as its immediately connected neighbours as follows:  N_{i} = \{v_{j} : e_{ij} in E lor e_{ji} in E\} \n
 *  k_{i} as the number of vertices, in the neighbourhood,  N_{i}, of a vertex. \n
 *  The local clustering coefficient C_{i} for a vertex v_{i} is then given by a proportion of the number of links between the vertices within its neighbourhood \n
 *  divided by the number of links that could possibly exist between them. \n
 *  These measures are 1 if every neighbour connected to v_{i} is also connected to every other vertex within the neighbourhood, and 0 if no vertex that is connected to \n
 *  v_{i} connects to any other vertex that is connected to v_{i}.
 */
bool PhysicalGraph::calcClusteringCoefficient(bool showInConsoleDetailed = false) {//calcClusteringCoefficient
    unsigned long sum_num = 0;
    unsigned long sum_den = 0;
    for(unsigned int pi=srcV; pi<=numV; pi++){
        PhysicalNode& piServer = PNode[pi];
        unsigned long denominator;
        if(piServer.neighbours.size() <= 1) ///< clustering coefficient would be zero
            continue;
        else
            denominator = (piServer.neighbours.size()*(piServer.neighbours.size()-1));

        unsigned long numerator=0;
        for(unsigned int pj: piServer.neighbours){ /// all neighbours of the node pi.
            for(unsigned int pk: piServer.neighbours){
                if(pj!=pk)
                    numerator += ((mat[pi][pj]>0)and(mat[pj][pk]>0)and(mat[pk][pi]>0));
            }
        }
        sum_num +=  numerator;
        sum_den +=  denominator;
        clusteringcoeff_mean += piServer.clusteringcoeff_local = numerator*(1.0/denominator);
        if(showInConsoleDetailed){
            cout<<"\npi:"<<pi<<"  CoeffL["<<numerator<<"/"<<denominator<<" = "<<piServer.clusteringcoeff_local<<"]  ";
            cout<<"Neigh("<<piServer.neighbours.size()<<") -> [";
            for(const auto& n: piServer.neighbours)  cout<<n<<","; cout<<"]";
        }
    }//for each server
    clusteringcoeff_mean = clusteringcoeff_mean/numV;
    clusteringcoeff_global = (double)1.0*sum_num/sum_den;
    if(debug)cout<<"\n\tClustering Coeff mean:"<<clusteringcoeff_mean<<" | global:"<<clusteringcoeff_global;
    return true;
}//calcClusteringCoefficient

#endif //SFC_PARALLELIZATION_PHYSICALGRAPH_H
