//
// Created by vijay on 20-04-2023.
//

#ifndef SFC_PARALLELIZATION_SIMULATIONCLASS_H
#define SFC_PARALLELIZATION_SIMULATIONCLASS_H

/* ******** Some Structers *************** */
/*! @brief Before finding best mapping of SFC, this will save delays based on VNFs in SFC to pre-calculate in order to avoid multiple computations of same VNF.*/
struct vnfDelaysPreComputed{
    type_delay exeDelay{}; ///< execution delay of the VNF
    type_delay prcDelay{}; ///< processing delay of the VNF
    unordered_map<unsigned int, type_delay> queuingDelay; ///< queuing delay of the VNF and its instance
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
//        return (static_cast<int>(a*100.0) == static_cast<int>(b*100.0));
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
    unordered_map<unsigned int, vector<unsigned int>> cntPN; ///< count of physical node,frequency in the instance combination
    unordered_map<unsigned int, type_delay> exePN; ///< maximum execution time of physical node in the instance combination
    type_delay utilization{0}; ///< utilisation percentage of all the instances present in the lgNode
    vector<pair<unsigned int, type_delay>> children; ///<  next stage lgNode index and its distance, that is pair of this->node = {next stg node, min dist}.
    vector<pqNode> kpaths; ///< number of shortest path traverse through this lgNode

    lgNode()=default;
    explicit lgNode(unsigned int index):idx(index){};
    lgNode(unsigned int index, const vector<pair<unsigned int,unsigned int>>& givenIC):idx(index), instCombination(givenIC){ }
};

/*! @brief For a single SFC, according to the algorithm what is the optimal/best (minimum delay) parameters we have found.*/
struct SFC_RESULT{
    int seq_pid{noResSeq}, ppar_pid{noResPar}, fullpar_pid{noResPar}; ///< idx of allPartParSFC Array for sequential/part parallel/full parallel for which algorithm give optimal answer.
    type_delay seq_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of sequential length chain according to our algorithm.
    type_delay ppar_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of partial parallel chain according to our algorithm.
    type_delay fullpar_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of full parallel chain according to our algorithm.
    double seq_duration{0}, ppar_duration{0}, fullpar_duration{0};
    unordered_map<unsigned int, unsigned int> seq_fninstmap, ppar_fninstmap, fullpar_fninstmap; ///< Best mapping {fun->its instance taken} for given sequential length chain according to our algorithm.
};

/*! @brief Final Output produced by the algorithm.*/
class SimTEST{
public:
    std::string name{}; ///< name of the solution
    unordered_map<unsigned int, SFC_RESULT> sfcsol; ///< sfc index to its solution values
    double seq_duration{}, ppar_duration{}, fullpar_duration{}; ///< time taken to construct the solution with/without parallelism
    unordered_map<unsigned int, unordered_map<unsigned int, type_delay>> seq_utilization, ppar_utilization, fullpar_utilization; ///< utilization of the VNF_Inst. {VNFid -> {instid -> utilization}} for sequential and parallel chain


//    unsigned int packetBodySize = 1000, packetHeaderSize = 24; ///<Size of the Network Packet Body and Header. in Bytes. Type = unsigned int Range[0,4294967295].
//    unsigned int factor_packet = 8; ///<factor to multiply in order to convert packet size in bits. 1 Byte is 8 bits
//
//    unsigned int bandwidthNW = 10; ///<Bandwidth of the Network. in Mega bits per second. 1Gb = 1000 Mb. Type = unsigned int Range[0,4294967295]. *Mb[0,4294967295], Gb_in_Mb[1000 , 4294967.295].
//    unsigned int factor_bandwidth = 1000000;  ///< factor to multiply to convert bandwidth in bits/seconds.
//
//    unsigned int speedOfLight = 300000000;///<speed of light in vaccum 3 * 10^8 m/s
//    type_delay velocityFactor = 1.0; //<velocity factor of transmission medium. vaccum = 1.0. copper wise = 0.7
//    type_delay read_write_time_per_bit = 0.077e-3; ///<0.077ms (measured by duplicating a large file of 1 MB in a server with Intel i7-8700 core
//    type_delay timesfactor = 1;
//    type_delay timesfactor_pkt=1;
//    type_delay timesfactor_tx=10;
//    type_delay timesfactor_px=10;
//    type_delay timesfactor_fnExe=1;
//    type_delay timesfactor_qd=10;

    explicit SimTEST(string name){
        this->name = std::move(name);
    }
    SimTEST() = default;


};

class Simulations{
public:
    string sim_name{};
    string dirName{}/*!< directory name */, fullDirName{} /*!< directory name with slash at the end*/;
    /// Common Object
    PhysicalGraph PhysicalNetwork; ///< Graph Object contains network related functions
    VirtualNetworkFunctions VNFNetwork; ///< VNF Object contains VNF (Virtual Network Function) related code.

    /// SFC Related
    string filename_sfc{};///< name of the file from where SFCs are read;
    vector<ServiceFunctionChain> allSFC; ///< SFCs object contains code related to Service function chains
    vector<ServiceFunctionChain*> sortedSFCs; ///< SFCs sorted for deployement in order of their priority of traffic rate and length
    unordered_map<unsigned int, unordered_set<unsigned int>> VNF_2_SFC; ///< VNFs Type are present in what SFCs Requests (sfc vector idx).

    /// VNF related
    unsigned int numParallelPairs{0}, totalPairs{0};
    unordered_map<unsigned int, unordered_map<unsigned int, unsigned int>> parallelPairs; ///< {i_vnfid ->{j_vnfid it is parallel with copy/without copy}} pairs identifying which are parallel
    type_delay likelihood_mean; ///< Mean parallelism likelihood of all funtion pairs (fi,fj);
    unordered_map<unsigned int, unordered_map<unsigned int, double>> likelihood;/*!< Likelihood of parallelizing two functions fi and fj defined as L_{fi,fj} = |num of sfc in which fi and fj are present and are parallelizable| divided by |numOfSFCS|*/

    ///VNF Deployement Mappings
    unordered_map<unsigned int, unsigned int> finalInstancesCount; ///< count of instances of each function
    unordered_map<unsigned int, vector<pair<unsigned int,unsigned int>>> PN_2_VNF;
    unordered_map<unsigned int, unordered_map<unsigned int, unsigned int>> I_VNFINST_2_PN; ///< VNF {type,inst} is hosted on which PN. {VNFid -> {instid -> PN id}} ie. arr[vnf][inst]=pnid;, inst->1based indexing

    /// Test
    unordered_map<string, SimTEST> TestsResult;
    unordered_map<unsigned int, vnfDelaysPreComputed> vnfDelays;///< pre-calculated VNF delays (processing, execution and queuing delay). Before iterating all the partial par sfc it is better to calculate it for each chain as they will be same for each chain.


    Simulations(string simulation_name, string directory){
        this->sim_name = std::move(simulation_name);
        this->dirName = std::move(directory);
        this->fullDirName = dirName + "/";
        createDirectory();
    };

    bool createDirectory();
    bool readGenericServiceFunctionsChains(const string &, bool = false);
    bool findRandomParallelPairs(const float &, int = 0);

    bool readDataFromFileInit(const string &, const string &, const string &, const pair<float, int>& = {50,2}, const int& = 1);

    void convert_SeqSFC_to_FullParallel_vs1(ServiceFunctionChain &);

    void convert_parVNFBlk_to_PartialChains_vs1(ServiceFunctionChain &, bool=false, bool=false);

    void showSFCsDescriptions();

    bool calcLikelihoodOfTwoFunctions(int=1, bool=false);

    bool DeploymentVNF_ScoreMethod(float , float , int , int , bool , bool );

    void showSimulationTestResults(const SimTEST &result);

    void showVNFsUtilization(const int &, const SimTEST&) const;

    void getUtilization(const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>&) const;

    void showPNsDescription() const;

    void showVNFsDescription() const;
};


/*!
 * Create a directory in input and output folder with test name.
 * @param testDirName name of the directory to be created
 * @return true if file does not exist and created now. Else False.
 */
bool Simulations::createDirectory(){
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

/*! Function will read the SFC chains data from the file in the input directory and read into ServiceFunctionChain class.\n\n
 * It will simultaneously read vnf data into sequential vector and also construct adj list. \n\n
 * First Line in file consist of numOfSFC S\n
 * from next line upto S times, each s group consist\n
 * sfc_index  (arrivalRate of SFC) (access nodes src dst) (numOfVNFs present except src and dest)   \n
 * This line contains VNFs Type id (for i=1 to numOFVNFs)
 * @param filename_sfc Name of the SFCs file which consist of the data
 * @param sample vnfs to sfc mapping.
 * @param showinConsole show the file read values in the console
 */
bool Simulations::readGenericServiceFunctionsChains(const std::string& filename, bool showinConsole) {
    allSFC.clear(), sortedSFCs.clear(), VNF_2_SFC.clear();
    filename_sfc = filename;
    ifstream fin;
    string filepathExt = input_directory + fullDirName + filename;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    priority_queue<ServiceFunctionChain*, vector<ServiceFunctionChain*>, comparator_sfc> pqSortSFC;
    if(debug)cout<<"\n[Reading Original Sequential SFCs] SFC File: "<<filepathExt<<endl;
    unsigned int i_nSFC; ///< input num of SFCs
    if(fin>>i_nSFC){
        allSFC = vector<ServiceFunctionChain>(i_nSFC);
        unsigned int readSFC_idx, readSFC_totalVNF, vnfid; type_delay readSFC_arrivalRate ;
        for(int ni=0; ni<i_nSFC; ni++) ///< For each node there is a row which would be read
        {
            pair<unsigned int, unsigned int> readSFC_accessnode;
            if (!(fin>>readSFC_idx>>readSFC_arrivalRate>>readSFC_accessnode.first>>readSFC_accessnode.second>>readSFC_totalVNF)){
                std::cerr << "\t SFCs Reading Failed: SFC Row["<<ni<<"]\n";
                return false;
            }
            string readSFC_name = "SFC"+to_string(readSFC_idx);
            allSFC[ni] =  ServiceFunctionChain(readSFC_idx, readSFC_name,readSFC_accessnode, readSFC_totalVNF, readSFC_arrivalRate);
            for(int vj = 1; vj<=readSFC_totalVNF; vj++) { ///< read VNFs type ID.
                if (fin >> vnfid) {
                    allSFC[ni].vnfSeq.push_back(vnfid);
                    VNF_2_SFC[vnfid].insert(ni);
                }else {
                    std::cerr << "\t SFCs Reading Failed: SFC Row[" << ni << "] fn:" << vnfid << "\n";
                    return false;
                }
            }
            pqSortSFC.push(&allSFC[ni]);
            if(showinConsole) {
                cout << readSFC_idx << "-" << readSFC_name << "-" << readSFC_totalVNF << "-" << readSFC_arrivalRate<<"- [";
                for(const auto& x: allSFC[ni].vnfSeq) cout<<x<<" -> "; cout<<"]\n"; }
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
    return true;
}

/*!
 * Reading Data First Time
 * @param fileNetwork FileName of the Network present in the directory.
 * @param fileVNF FileName of the VNFs to read VNFs data.
 * @param fileSFC fileName of the SFCs to read SFCs data.
 * @param parallelPairsOpt how you want parallel pairs to be generated
 * @param likelihoodOpt how you want likelihood mean to be generated
 * @return status
 */
bool Simulations::readDataFromFileInit(const string& fileNetwork, const string& fileVNF, const string& fileSFC, const pair<float, int>& parallelPairsOpt, const int& likelihoodOpt){
    if( readNetwork(fullDirName,fileNetwork, PhysicalNetwork) and
        readVirtualNetworkFunctions(fullDirName,fileVNF, VNFNetwork) and
        readGenericServiceFunctionsChains(fileSFC) and
        findRandomParallelPairs(parallelPairsOpt.first,parallelPairsOpt.second) and /// based on VNFs
        calcLikelihoodOfTwoFunctions(likelihoodOpt) ///based on parallel pairs and VNFs in SFC
        ){
        for(ServiceFunctionChain& sfc: allSFC){
            convert_SeqSFC_to_FullParallel_vs1(sfc);
            convert_parVNFBlk_to_PartialChains_vs1(sfc);
        }
        if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";
    }
    else return false;
    return true;
}

/*!
 * Show description of all the SFCs in the Network
 */
void Simulations::showSFCsDescriptions(){
    cout << "\n\n ----- Service Functions Chains Description ::\n";
    cout<<std::setw(2)<<"Id"<<std::setw(7)<<"Name\t"<<"len\t"<<"TrafficRate\t"<<"PartialChains\t"<<"(ps->pd)\t"<<"Sequential SFC\t"<<"Fully Parallel SFC";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (const auto &sfc: allSFC){
        cout<<"\n"<<std::setw(2)<<sfc.index<<" | ";
        cout<<std::setw(7)<<sfc.name<<" | ";
        cout<<std::setw(2)<<sfc.numVNF<<" | ";
        cout<<std::setw(4)<<sfc.trafficArrivalRate<<"ms | ";
        cout<<std::setw(2)<<sfc.allPartParSFC.size()<<" | ";
        cout<<std::setw(7)<<"("<<sfc.access_nodes.first<<"->"<<sfc.access_nodes.second<<") | ";
        cout << "(";
        for(const auto& fn : sfc.vnfSeq){ cout <<"f"<< fn ;   cout<< ";->";
        } cout << ")\t";
        cout << "(";
        for(const auto& blk: sfc.allPartParSFC[sfc.allPartParSFC.size()-1]){
            cout<<" ["; for(int fn: blk){
                cout <<"f"<< fn ;  cout<<"; ";
            }
            cout<<"]->";
        } cout << ")";
    }
}

/*!
 * Make VNF pairs (uniform (approx 50%)) out of total pairs to be parallelizable/independent so that they can execute parallel. Using Uniform Int Distribution.  \n
 * Total Pairs = numVNF*(numVNF-1) ordered (1,2)(2,1) except (1,1)(2,2)... \n
 * It then writes the pairs to the filename_vnf_parallelpairs {i_vnf -> {set of j_vnf it is parallel to}}
 * @param threshold percentage of pairs we want to parallelize. (1-100).
 * @param option way to obtain parallel pairs, 0 -> fixed, 1->approx around threshold, 2->exact to threshold.
 */
bool Simulations::findRandomParallelPairs(const float& threshold, int option){
    if(threshold<=0 or threshold>=101) return false;
    const unsigned int numVNF = VNFNetwork.numVNF;
    totalPairs =  numVNF*(numVNF-1) - numVNF;
    if(option == 0) {
//        numParallelPairs = 43;
//        parallelPairs={
//                {1,{{9,pktNoCopy},{7,pktNoCopy},{6,pktNoCopy},{5,pktNoCopy},{4,pktNoCopy},{2,pktCopy}}},
//                {2,{{10,pktCopy},{7,pktNoCopy},{5,pktCopy},{4,pktNoCopy}}},
//                {3,{{8,pktNoCopy},{7,pktCopy},{6,pktNoCopy},{2,pktCopy}}},
//                {4,{{9,pktNoCopy},{8,pktCopy},{2,pktNoCopy},{1,pktNoCopy}}},
//                {5,{{10,pktNoCopy},{9,pktCopy},{8,pktCopy},{4,pktNoCopy},{3,pktNoCopy}}},
//                {6,{{9,pktCopy},{5,pktNoCopy},{3,pktNoCopy}}},
//                {7,{{9,pktCopy},{5,pktCopy},{4,pktCopy},{3,pktCopy}}},
//                {8,{{10,pktCopy},{7,pktNoCopy},{2,pktNoCopy},{1,pktNoCopy}}},
//                {9,{{10,pktCopy},{8,pktCopy},{7,pktNoCopy},{4,pktCopy},{1,pktNoCopy}}},
//                {10,{{8,pktNoCopy},{6,pktCopy},{3,pktCopy},{2,pktNoCopy}}}
//        };
        numParallelPairs = 40;
        parallelPairs={
            {1,{{7,1},{10,1},{9,1},{5,1},{6,2}}},
            {2,{{10,2},}},
            {3,{{5,1},{2,2},{8,2},{7,2}}},
            {4,{{8,2},{1,2}}},
            {5,{{3,2}}},
            {6,{{7,2},{1,2},{9,2},{8,2},{10,1},{2,1}}},
            {7,{{8,1},{4,2},{5,2},{10,2},{3,1},{9,2}}},
            {8,{{10,1},{2,1},{6,1},{3,1},{4,2}}},
            {9,{{10,2},{5,2},{8,1},{6,2},{4,1},{3,2},{1,1}}},
            {10,{{5,1},{1,2},{9,2}}}
        };
        if(debug)cout<<"\n\tRandom Parallel Pairs: "<<numParallelPairs<<" | "<<(100.0*numParallelPairs/totalPairs)<<" % out of total:"<<totalPairs<<" unique pairs.";
        return true;
    }

    parallelPairs.clear();
    numParallelPairs=0;
    const unsigned int numOfIndependentPairsNeeded = round(totalPairs*(threshold/100.0));

    unsigned seed = chrono::system_clock::now().time_since_epoch().count(); ///< (present time) and clock's epoch
    std::mt19937 rd_generator(seed); ///< mt19937 is a standard mersenne_twister_engine
    std::uniform_int_distribution<unsigned int> withCopyOrWithout_distribution(pktNoCopy, pktCopy);      ///< for selecting operation Parallel/Not Parallel

    if(option == 1){
        std::uniform_int_distribution<unsigned int> rd_threshold_distribution(1, totalPairs);      ///< for selecting operation Parallel/Not Parallel
        for(unsigned int i_vnf=VNFNetwork.srcVNF; i_vnf<=numVNF; i_vnf++) {
            for (unsigned int j_vnf = VNFNetwork.srcVNF; j_vnf <= numVNF; j_vnf++) {
                if(i_vnf == j_vnf)continue;
                if( rd_threshold_distribution(rd_generator) <= numOfIndependentPairsNeeded){
//                    if(parallelPairs.count(j_vnf)>0 and parallelPairs[j_vnf].count(i_vnf)>0){ ///< if fj->fi value exist then take that value
//                        parallelPairs[i_vnf][j_vnf] = parallelPairs[j_vnf][i_vnf];
//                    }else{ /// else generate new value of fi fj pair
                        parallelPairs[i_vnf][j_vnf] = withCopyOrWithout_distribution(rd_generator);
//                    }
                    numParallelPairs++;
                }// threshold %pairs are generated
            }//j
        }// i
    }else if (option == 2){
        std::uniform_int_distribution<unsigned int> fn_distribution(VNFNetwork.srcVNF, numVNF);
        for(unsigned int i=1; i<=numOfIndependentPairsNeeded; i++) {
            unsigned int fi = fn_distribution(rd_generator);
            unsigned int fj = fn_distribution(rd_generator);
            while( (fi == fj) or parallelPairs[fi].count(fj)>0 ){ // fj same as fi or fi already paired with fj
                fj=fn_distribution(rd_generator);
            }
//            if(parallelPairs.count(fj)>0 and parallelPairs[fj].count(fi)>0){ ///< if fj->fi value exist then take that value
//                parallelPairs[fi][fj] = parallelPairs[fj][fi];
//            }else{ /// else generate new value of fi fj pair
                parallelPairs[fi][fj] = withCopyOrWithout_distribution(rd_generator);
//            }
            numParallelPairs++;
        }
    }//opt 2

    if(debug)cout<<"\n\tRandom Parallel Pairs: "<<(100.0*numParallelPairs/totalPairs)<<"% | "<<numParallelPairs<<" out of "<<totalPairs<<" pairs.";
    /// writing to file.
    ofstream fout;
    string filepathExt = output_directory+fullDirName+"Sim_"+sim_name+"_Th_"+to_string((100*numParallelPairs/totalPairs))+"_ParallelPairs_"+VNFNetwork.filename_vnf;///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fout<<"\n{";
    for(int i_vnf=VNFNetwork.srcVNF; i_vnf<=numVNF; i_vnf++) {
        fout<<"\n  {"<<i_vnf<<",{";
        for(const auto& j_vnf: parallelPairs[i_vnf]){;
            fout<<"{"<<j_vnf.first<<","<<j_vnf.second<<"},";
        }
        fout<<"}},";
    }// i
    fout<<"\n};";
    fout.close();
    return true;
}

/*! Likelihood of parallelizing two functions fi and fj defined as
 * L_{fi,fj} = |num of sfc in which fi and fj are present and are parallelizable| divided by |numOfSFCS|
 * @param opt opt=1 mean calculated using SumLikelihood/(totalpairs - numvnfs) unique pairs, IF opt=2 then  SumLikelihood/(numParallelPairs)
 */
bool Simulations::calcLikelihoodOfTwoFunctions(int opt, bool showInConsoleDetailed){
    likelihood.clear();

    auto intersect = [&](const unordered_set<unsigned int>& SFC_OF_FI, const unordered_set<unsigned int>& SFC_OF_FJ)->unsigned int{
        unsigned int count=0;
        for(const auto& sfcidx: SFC_OF_FI){
            if(SFC_OF_FJ.count(sfcidx) > 0) count++;
        }
        return count;
    };

    likelihood_mean = 0;
    for(const auto& [fi, IsParallelWith]: parallelPairs){
        for(const auto& [fj,_]: IsParallelWith){
            likelihood_mean+= likelihood[fi][fj] = (double)intersect(VNF_2_SFC[fi], VNF_2_SFC[fj])/allSFC.size();
        }
    }
    if(showInConsoleDetailed){
        cout<<"\nLikelihood:";
        for(unsigned int fi=VNFNetwork.srcVNF; fi<=VNFNetwork.numVNF; fi++){
            cout<<"\n\tf"<<fi<<" with ";
            for(const auto& [fj, val]: likelihood[fi]){
                cout<<"[f"<<fj<<"= "<<likelihood[fi][fj]<<"] ";
            }
        }
        cout<<"\n   likelihood_mean (/uniquepairs): "<<likelihood_mean/(totalPairs);
        cout<<"\n   likelihood_mean (/parallelpairs): "<<likelihood_mean/(numParallelPairs);
    }

    if(opt == 1){
        likelihood_mean = likelihood_mean/(totalPairs); // total unique pairs.
    }else if(opt == 2){
        likelihood_mean = likelihood_mean/(numParallelPairs); // total unique pairs.
    }
    if(debug)cout<<"\n\tLikelhood_fi_fj mean: "<<likelihood_mean;
    return true;
}

/*! It reads data from sequential vector and then convert into full Parallel SFC (parallel VNF blocks) chain by detecting
 * whether two function can be parallelised or not.\n
 * Each block denote fully parallel VNFs in that step.
 * @param[in, out] SFC ServiceFunctionChain object to convert sequential into parallel.
 */
void Simulations::convert_SeqSFC_to_FullParallel_vs1(ServiceFunctionChain& csfc){//convert_SeqSFC_to_FullParallel_vs1
    unsigned int sz = csfc.vnfSeq.size(); // total vnfs including src and dest

    unsigned int cid=0;
    csfc.vnfBlocksPar.push_back({csfc.vnfSeq[cid]}); // stage 0, pushing src node
    for(cid=1; cid<sz; cid++){ // from SFCsrc+1 stg to SFCdst-1 stage
        unsigned int prv_vnf = csfc.vnfSeq[cid-1], cur_vnf = csfc.vnfSeq[cid];  // Checking prv_vnf --> cur_vnf pairs
        // NOT PARALLEL, if prv_vnf does not exist as first vnf in pair or if prv_vnf exist it is not parallel to cur_vnf, then it is not parallel pair
        if(parallelPairs.count(prv_vnf)==0 or parallelPairs.at(prv_vnf).find(cur_vnf) == parallelPairs.at(prv_vnf).end()){
            csfc.vnfBlocksPar.push_back({cur_vnf});// push cur_vnf as separte new stage
        }else if(csfc.vnfBlocksPar.back().size() == 1) { //PARALLEL Pairs STG Size=1 and size of previous stg is just one, then we can directly push into that stg.
            csfc.vnfBlocksPar.back().push_back(cur_vnf);
        }else{ // PARALLEL Pairs and STG size > 1, from previous stg we have to check if cur_vnf is parallel to all prv_vnf in previous stage
            bool pushInLstStg = true;
            for(const unsigned int& lst_stg_vnf: csfc.vnfBlocksPar.back()){
                if(parallelPairs.count(lst_stg_vnf)==0 or parallelPairs.at(lst_stg_vnf).find(cur_vnf) == parallelPairs.at(lst_stg_vnf).end()){
                    pushInLstStg = false;   break;
                }
            }
            if(pushInLstStg) csfc.vnfBlocksPar.back().push_back(cur_vnf);
            else csfc.vnfBlocksPar.push_back({cur_vnf});
        }//lst stg size>1
    }// for cid
}//convert_SeqSFC_to_FullParallel_vs1

/*! generate all the feasible partial parallel SFC for the given full parallel SFC.
 * @param nVNFs is number of VNFs except src and dest in fully parallel SFC.
 * @param clusterSz for each of the cluster enumeration for size[nVNFs].
 * @param nCk n=block size and k={cluster_i[l] value i.e. l(level) index value dentoes number of function(k) to be chosen as parallel out of all functions(n) in blocks.
 * @param fullParVNFBlocks is fully parallel VNF Blocks in sequence where each block/step denotes all the parallelizable functions in that block/step.
 * @param[out] allPartParSFC All the Partial parallel Clusters of the fully parallel VNF Blocks. Each SFC is without src and dest.
 * @example: fullParVNFBlocks={{1},{2,3,4,5}}, nVNFs:5 (f1,f2,f3,f4,f5).\n
 * Blk(0) = {1} only 1 parallel function, Blk(1) = {2,3,4,5} all 4 are parallel.\n
 * clusterSz[5(size=nVNFs)] = {{1, 2,1,1},{1,4},{3,1,1} ... } so on. cluster_i[l] denotes number of function parallel in that level. \n
 * where cluster_i {1, 2, 1, 1} means in level[0] only one function runs, level[2] = 2 functions run together, and level[3] and level[4] one-one function are there \n
 * {[1], [2,1,1]} is mapped to {[f1], [f2,f3,f4,f4]} such that [ 1-c combination of f1 ]-> [2-c combination of f2,f3,f4,f5] -> [1-c combination of f2,f3,f4,f5] -> [1-c combination of f2,f3,f4,f5] \n
 */
void Simulations::convert_parVNFBlk_to_PartialChains_vs1(ServiceFunctionChain& csfc, bool showInConsole, bool showInConsoleDetailed){
    const vector<vector<unsigned int>>& fullParVNFBlocks = csfc.vnfBlocksPar;
    const unsigned int& nBlk = fullParVNFBlocks.size(); ///< number of blocks of the fully parallel SFC, (including src and dst block)
    vector<vector<vector<unsigned int>>>& allPartParSFC = csfc.allPartParSFC; ///< All the Partial parallel Clusters of the fully parallel VNF Blocks. Each SFC is without src and dest.
//    const vector<vector<unsigned int>> SK = { {1,  2,1,1},{1, 2,2}, {1,4}, {2,2,1}, {3,1,1} };

    /// lambda backtrack function to find all partial sfc corresponding to cluster_i and parSFC_Full.
    std::function<void(unsigned int,vector<vector<unsigned int>>&,unordered_map<unsigned int,vector<unsigned int>>&, unsigned int&)> findAllPartSFC_Backtrack
            =[&findAllPartSFC_Backtrack, &allPartParSFC, &fullParVNFBlocks]
                    (unsigned int cur_level_idx, vector<vector<unsigned int>>&partSFC, unordered_map<unsigned int,vector<unsigned int>>& levelInfo, unsigned int& mask) ->void
            {
                // if current level is equal to total level in cluster/SFC.
                if(cur_level_idx == levelInfo.size()) {
                    allPartParSFC.push_back(partSFC);
                    return;
                }

                unsigned int blkid = levelInfo[cur_level_idx][0], nl = levelInfo[cur_level_idx][1], kl= levelInfo[cur_level_idx][2];
                for(const vector<unsigned int>& curCombination: nCk[nl][kl]){ //n=4,k=1 {{1}, {2}, {3}, {4}} } | k=2 {{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}} }
                    bool thisCombinationCanBeVisited = true;

                    vector<unsigned int> curBlkFunc; //< if we can visit this combination then this vector will be current blk
                    unsigned int localMask=0;
                    for(const unsigned int& idx: curCombination){ // check if the curCombination idx{2,3} mapped to block func {fw (id 1-1), fx (id 2-1), fy (id 3-1)} --> fx, fy can be visited or its node already visited.
                        int fn_id = fullParVNFBlocks[blkid][idx-1];
                        if((mask & (1<<fn_id)) != 0){ thisCombinationCanBeVisited=false; break; }
                        curBlkFunc.push_back(fn_id);
                        localMask |= (1<<fn_id);
                    }
                    if(thisCombinationCanBeVisited){
                        mask |= localMask;
                        partSFC.push_back(std::move(curBlkFunc));
                        findAllPartSFC_Backtrack(cur_level_idx+1, partSFC, levelInfo, mask);
                        mask ^= localMask;
                        partSFC.pop_back();
                    }
                }
            };

    for(const vector<unsigned int>& cluster: clusterSz[csfc.numVNF]){ //
        if(cluster[0] > fullParVNFBlocks[0].size()) // if in first block(after src blk) 2 function is there, but cluster saying 3 needs to be parallel then continue next cluster.
            continue;

        unsigned int cur_level=0; ///< current level at which we have to insert all the nodes
        unordered_map<unsigned int,vector<unsigned int>> levelInfo;///< this stores level wise info {level i -> {0 -> blkId, 1->n, 2-> k}} to find nCk for block blkId of parSFC_Full


        if(showInConsoleDetailed){
            cout<<"\ncluster["; for(const auto& x: cluster) cout<<x<<" "; cout<<"]";
        }
        bool allBlksOfSFCDone = true;
        // iterate through the blocks except src (1) and dst(nBlk-1) and map it to cur cluster.
        for(unsigned int blk_id=0; blk_id<nBlk; blk_id++){
            const unsigned int& curBlk_size = fullParVNFBlocks[blk_id].size();  ///< current block size, that is num of VNFs present in it.

            /*!
             * It checks whether cluster_i is a feasible vector size and we can find same number of parallel VNFs specified by the cluster. cluster_i[l] number of func can run in parallel. \n
             * If we cannot find the same number of parallel VNFs specified by cluster_i, then it is not feasible. e.g. [2,2] is not feasible for {f1,f2,f3} block. \n
             * feasible [1, 2,1,1] and {{f1},{f2,f3,f4,f5}} where [1] is mapped to {f1} and  entire [2+1+1] is mapped to {f2,f3,f4,f5}
             */
            bool foundMapping=false; unsigned int delta=0; //< delta is range of level [endLevel-curLevel] upto which we can parallelize current block.
            for(unsigned int li=cur_level, sum_s=0; li<cluster.size() and sum_s < curBlk_size; li++){
                sum_s += cluster[li]; delta++;
                if(sum_s == curBlk_size){
                    foundMapping = true; break;
                }
            }
            if(!foundMapping) {
                if(showInConsoleDetailed){ cout<<"\tNot feasible for block{"; for(const auto& x: fullParVNFBlocks[blk_id]) cout<<"f"<<x<<","; cout<<"}";}
                allBlksOfSFCDone = false; // break lag gya isliye dfs call mt krna
                break;
            }

            for(unsigned int li=cur_level; li<cur_level+delta; li++){
                if(showInConsoleDetailed){
                    cout<<"\n\tL["<<li+1<<"]: \t";
                    for(const auto& allComb: nCk[curBlk_size][cluster[li]]) {
                        cout<<"{"; for(auto node: allComb) cout<<"f"<<fullParVNFBlocks[blk_id][node-1]<<",";  cout<<"}";
                    }
                }
                levelInfo[li] = {blk_id, curBlk_size, cluster[li]};
            }
            cur_level = cur_level+delta;
        }
        if(allBlksOfSFCDone) { //            cout<<"dfs called.";
            vector<vector<unsigned int>>partSFC; unsigned int mask=0;
            findAllPartSFC_Backtrack(0, partSFC, levelInfo, mask);
        }
    }

    if(showInConsole){
        cout<<"\nTotal PartSFC:"<<allPartParSFC.size();
        for(int idx=0; idx<allPartParSFC.size(); idx++){
            const auto& PCs = allPartParSFC[idx];
            cout<<"\nPC["<<idx+1<<"] ( "; unordered_map<unsigned int,unsigned int> freq;
            for(const auto& blks: PCs){
                cout<<"["; for(auto fn_id: blks){
                    cout<<"f"<<fn_id<<" ";
                    if(++freq[fn_id]>1)  throw runtime_error("Error in calculation of allPartSFC. Some fn_id repeated");
                } cout<<"]";
            } cout<<")";
        }
    }

}

/*!
 * @param fnInstScalingFactor (0.5-1] how much scaling (<=1 ==1 full scaling) in function instance we wanted as compared to actual instances needed (val=0).
 * @param alphaDecayRate (0-1) weightage of giving more importance to far away nodes
 * @param choice_pxs opt=1 choose bst candidate node for fcur, opt=2 choose bst candidate node for fcur-1(prev)
 * @param choice_dist_pref opt=1 distance prefernce based on hop count, opt=2 based on actual nearest distance
 * @param showInConsole final output in console
 * @param showInConsoleDetailed detailed ouput in console
 */
bool Simulations::DeploymentVNF_ScoreMethod(float fnInstScalingFactor=0.5,float alphaDecayRate=0.75, int choice_pxs = 1, int choice_dist_pref = 1, bool showInConsole = false, bool showInConsoleDetailed=false){//DeploymentVNF_ScoreMethod

    if((fnInstScalingFactor<0 or fnInstScalingFactor>1) or (alphaDecayRate<=0 or alphaDecayRate>=1) or (choice_pxs<=0 or choice_pxs >=3) or (choice_dist_pref<=0 or choice_dist_pref >=3))
        return false;

    finalInstancesCount.clear();
    PN_2_VNF.clear();
    I_VNFINST_2_PN.clear();

    const unsigned int& numVNFs = VNFNetwork.numVNF;
    const unsigned int& numPNs =  PhysicalNetwork.numV;
    const unsigned int& numSFCs = allSFC.size();

    vector<unordered_map<unsigned int , unordered_map<unsigned int, double>>> dist_sfc_f_p(numSFCs); ///< score/priority of sfc s to put function fi on physical node pn.
    ///< Calculation of priority of sfc to place f on p based on its shortest path
    for(unsigned int sid=0; sid<numSFCs; sid++){ /*! for each sfc */
        const ServiceFunctionChain& sfc = allSFC[sid];
        const unsigned int& Len = sfc.numVNF;

        unordered_map<unsigned int, unsigned int> f_bstloc;
        unordered_map<unsigned int, unordered_map<unsigned int, double>>& priority_fi_pn = dist_sfc_f_p[sid];
        vector<unsigned int> access_nodes_shortest_paths = PhysicalNetwork.constructShortestPath(sfc.access_nodes.first, sfc.access_nodes.second);

        for(unsigned int vi=0; vi<Len; vi++){ ///finding best physical nodes for each function in sfc
            unsigned int path_id_star =   ceil((vi*access_nodes_shortest_paths.size()) / sfc.numVNF);
            f_bstloc[sfc.vnfSeq[vi]] = access_nodes_shortest_paths.at(path_id_star);
//            if(showInConsoleDetailed){ cout<<"\nf"<<sfc.vnfSeq[vi]<<" id:"<<path_id_star<<" pn:"<< access_nodes_shortest_paths[sid][path_id_star]; }
        }

        for(unsigned int i=0; i<Len; i++){ ///< priority for other nodes other than physical
            const unsigned int& fcur = sfc.vnfSeq[i];
            unsigned int pxs, pys; ///! px* -> pi -> py*,  (px* prev fun bst pNode) (pi current fun pNode) (py* next fun bst pNode)

            if(choice_pxs == 1){ ///< one way of choice pxs is the best cancidate of fcur
                pxs = f_bstloc[fcur];
            } else if(choice_pxs == 2){ /// othewise chose the previous f bst candidate to be pxs
                if(i==0) {  pxs = sfc.access_nodes.first;
                }else{      pxs = f_bstloc[sfc.vnfSeq.at(i-1)];  }
            }
            if(i ==  sfc.numVNF - 1){  pys = sfc.access_nodes.second;
            }else{                     pys = f_bstloc[sfc.vnfSeq.at(i+1)];  }

//            if(showInConsoleDetailed){ cout<<"\n[SFC:"<<sid<<"]"<<" f"<<fcur<<" (bst_p:"<<f_bstloc[fcur]<<") (px*:"<<pxs<<" py*:"<< pys<<") :: "; }
            for(unsigned int pi=PhysicalNetwork.srcV; pi<=numPNs; pi++)
            {
                double val_pxs_pi, val_pi_pys;
                if(choice_dist_pref == 1){
                    val_pxs_pi = PhysicalNetwork.constructShortestPath(pxs, pi).size()-1;
                    val_pi_pys = PhysicalNetwork.constructShortestPath(pi, pys).size()-1;
                    priority_fi_pn[fcur][pi] = (val_pxs_pi + val_pi_pys);
                }else if(choice_dist_pref == 2){
                     val_pxs_pi = PhysicalNetwork.dist[pxs][pi] ;// /(PhysicalNetwork.constructShortestPath(pxs, pi).size());
                     val_pi_pys = PhysicalNetwork.dist[pi][pys] ;// /(PhysicalNetwork.constructShortestPath(pi, pys).size());
                    priority_fi_pn[fcur][pi] = (val_pxs_pi + val_pi_pys);
                }
//                if(showInConsoleDetailed){ cout<<"\n\t"<<pxs<<"->"<<pi<<"->"<<pys<<" ("<<val_pxs_pi<<" + "<<val_pi_pys<<") | sum: "<<priority_fi_pn[fcur][pi] ;}
            }
        }///< priority for other nodes other than physical
    }/*! for each sfc */

    unordered_map<unsigned int, unordered_map<unsigned int, double>> score_f_p; ///< final score of placing vnf f on physical node p
    ///< Calculation of Final Score based on VNF clustering and Score Order.
    for (unsigned int f = VNFNetwork.srcVNF; f<=numVNFs; ++f) {///for each f
        for(unsigned int p = PhysicalNetwork.srcV; p<=numPNs; p++){///for each p
            double score_order_f_p=0; ///< score of placing vnf f on physical node p based on order of fn in sfc
            double score_cluster_f_p=0; ///< score of placing vnf f on physical node p based on independent functions clusters.

            for(int sidx=0; sidx<numSFCs; sidx++){
                if(dist_sfc_f_p[sidx].count(f) > 0){ ///< if the chain contains the function, otherwise zero
                    if(dist_sfc_f_p[sidx][f][p] != 0) /// chain contains function and dist is not zero
                        score_order_f_p += (1/dist_sfc_f_p[sidx][f][p]);
                    else
                        score_order_f_p += 1;
                }
            }

            for (unsigned int fj = VNFNetwork.srcVNF; fj <= numVNFs; ++fj){
                if(f == fj)continue;
                score_cluster_f_p += (likelihood[f][fj]-likelihood_mean) ;
            }
            score_cluster_f_p *= (PhysicalNetwork.PNode.at(p).clusteringcoeff_local - PhysicalNetwork.clusteringcoeff_mean);

            score_f_p[f][p] = score_order_f_p * (1 + score_cluster_f_p);  /// final score
            if(showInConsoleDetailed){cout<<"\nf"<<f<<" on p"<<p<<"\t is order: "<< score_order_f_p<<" + cluster: "<<score_cluster_f_p<<" = "<< score_order_f_p+score_cluster_f_p << "\t order(1+cluster)=" <<score_f_p[f][p];}
        }///for each p
    }///for each f


    unordered_map<unsigned int, unsigned int> f_inst_needed; ///< function instances needed according to usage in sfc.
    unsigned long sum_instances = 0;
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        type_delay arrival_rate_sfc = 0;
        for(const auto& sfcid: VNF_2_SFC[f]){
            arrival_rate_sfc += allSFC[sfcid].trafficArrivalRate;
        }
        sum_instances += f_inst_needed[f] = std::ceil(arrival_rate_sfc/ std::floor(VNFNetwork.VNFNodes[f].serviceRate));
        if(showInConsoleDetailed){
            cout<<"\n f:"<<f<<" is on #sfc("<<VNF_2_SFC[f].size()<<")";
//            for(const auto& sfcid: VNF_2_SFC[f]) cout<<"{"<<sfcid<<","<<allSFC[sfcid].trafficArrivalRate<<"}";
            cout<<" needed instances = "<<arrival_rate_sfc<<"/"<<VNFNetwork.VNFNodes[f].serviceRate<<" -> "<<f_inst_needed[f];
        }
    }

    vector<pair<unsigned int, unsigned int>> f_inst_final; ///< finall {function ->instance count mapping}
    unsigned int used_cores = 0; ///< number of cores actually used
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        unsigned int f_inst_on_scaled=0, instcnt=0;
        f_inst_on_scaled = floor(1.0*f_inst_needed[f] * PhysicalNetwork.sum_cores/sum_instances);
        if(PhysicalNetwork.sum_cores >= sum_instances){
            instcnt = (1-fnInstScalingFactor)*(f_inst_needed[f]) + fnInstScalingFactor*((float)f_inst_on_scaled);
        }else{
            instcnt = f_inst_on_scaled;
            if(f_inst_needed[f] > 0 and instcnt == 0 and used_cores < PhysicalNetwork.sum_cores){
                instcnt = 1;
            }
        }
        used_cores += instcnt;
        f_inst_final.emplace_back(f, instcnt);
        if(showInConsoleDetailed){  cout<<"\n f:"<<f<<" "<<PhysicalNetwork.sum_cores<<"/"<<sum_instances<<" scaled instances = "<<f_inst_on_scaled<<" final instances = "<<instcnt;}
    }

    ///< sorting function in descending order based on their original instance demand. scaling would be in similar manner.
    sort(f_inst_final.begin(), f_inst_final.end(), [&](const auto& left, const auto& right){
        const auto& first_instcnt = f_inst_needed[left.first];
        const auto& second_instcnt = f_inst_needed[right.first];
        if(first_instcnt == second_instcnt){
            return VNFNetwork.VNFNodes[left.first].serviceRate > VNFNetwork.VNFNodes[right.first].serviceRate;
        }
        else return first_instcnt > second_instcnt;
    });

    ///< in case of low cores, atleast use all cores OR if full utilization of all nodes is the request
    if((PhysicalNetwork.sum_cores < sum_instances and used_cores < PhysicalNetwork.sum_cores) or fnInstScalingFactor==1){
        int idx=0;
        while(used_cores < PhysicalNetwork.sum_cores){ ///< if cores are still unused and demand is there then fill it up
            f_inst_final[idx].second += 1;
            used_cores++;
            idx++;
        }
    }

    unordered_map<unsigned int,unsigned int> nodesToDeploy;
    for(unsigned int p = PhysicalNetwork.srcV; p<=numPNs; p++){
        nodesToDeploy[p] = PhysicalNetwork.PNode[p].capacity.cores;
    }

    for(const auto&[f, instcnt]: f_inst_final){
        finalInstancesCount[f] = instcnt;
        for(unsigned int i=1; i<=instcnt; i++){
            unsigned int pstar = 0; type_delay mxscore = 0;

            for(const auto& [p,_]: nodesToDeploy){
                if(score_f_p[f][p] > mxscore){
                    pstar = p;
                    mxscore = score_f_p[f][p];
                }
            }

            if(pstar == 0) break;
            I_VNFINST_2_PN[f][i] = pstar;  ///< saving ans
            PN_2_VNF[pstar].push_back({f, i});
            nodesToDeploy[pstar] -= 1;
            if(nodesToDeploy[pstar] <= 0) nodesToDeploy.erase(pstar);


            for(auto& [p, _]: nodesToDeploy){
                unsigned int path = PhysicalNetwork.constructShortestPath(p, pstar).size() - 1;
                double ratio = 1.0/path;
                score_f_p[f][p] = score_f_p[f][p] * pow(alphaDecayRate, ratio); //exp(ratio);

                for(const auto [fj,__]:parallelPairs[f]){
                    score_f_p[fj][p] = score_f_p[fj][p]*pow( 1 + likelihood[f][fj] , ratio);
                }
            }


        }
    }

    ofstream fout;
    string filepathExt = output_directory+fullDirName+"Sim_"+sim_name+"_VNFDeploymentScore_"+VNFNetwork.filename_vnf;///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fout<<"VNFs Deployement on the Network.";
    fout<<"\n   * opt pxs: "<<choice_pxs<<" Choice of pxs(opt=1 choose bst candidate physical node for vnfcur, opt=2 choose bst candidate node for fcur-1(prev))";
    fout<<"\n   * opt dist: "<<choice_dist_pref<<" Choice of distance in sfc preferece(opt=1 distance prefernce based on hop count, opt=2 based on actual nearest distance";
    fout<<"\n   * opt scale: "<<fnInstScalingFactor <<" Function Instance Scaling Factor (0.5-1] (how much scaling (<=1 ==1 full scaling) in function instance we wanted as compared to actual instances needed (val=0))";
    fout<<"\n   * opt rate: "<<alphaDecayRate <<" Alpha Decay Rate((0-1) weightage of giving more importance to far away nodes)";

    /// for readability
    fout<<"\n\nVNF to Node mapping:";
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        fout<<"\n   VNF["<<f<<"] [cnt: "<<f_inst_needed[f]<<"->"<<finalInstancesCount[f]<<"]  ";
        for(unsigned int i=1; i<=finalInstancesCount[f]; i++){
            fout<<"f"<<f<<char(i+96)<<" -> p:"<<I_VNFINST_2_PN[f][i]<<" | \t";
        }
    }
    fout<<"\n\nPhysical Node to VNFs mapping:";
    for(unsigned int p=PhysicalNetwork.srcV; p<=numPNs; p++){
        fout<<"\n   p:"<<p<<" [deg:"<<PhysicalNetwork.PNode[p].neighbours.size()<<"][Cl:"<<PhysicalNetwork.PNode[p].clusteringcoeff_local<<"] \t-> { ";
        for(const auto&[f, i]: PN_2_VNF[p])
            fout<<"f"<<f<<char(i+96)<<", ";
        fout<<"}";
    }
    fout<<"\n\tClustering Coeff mean:"<<PhysicalNetwork.clusteringcoeff_mean<<" | global:"<<PhysicalNetwork.clusteringcoeff_global;

    /// for resuse
    fout<<"\n\nfinalInstancesCount={";
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        fout<<"{"<<f<<","<<finalInstancesCount[f]<<"},  ";
    }
    fout<<"};";

    fout<<"\nI_VNFINST_2_PN={";
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        fout<<"\n{"<<f<<", {";
        for(unsigned int i=1; i<=finalInstancesCount[f]; i++){
            fout<<"{"<<i<<","<<I_VNFINST_2_PN[f][i]<<"}, ";
        } fout<<"}},";
    }
    fout<<"\n};";

    fout<<"\nPN_2_VNF={";
    for(unsigned int p=PhysicalNetwork.srcV; p<=numPNs; p++){
        fout<<"\n  {"<<p<<", {";
        for(const auto&[f, i]: PN_2_VNF[p])
            fout<<"{"<<f<<","<<i<<"}, ";
        fout<<"}},";
    }
    fout<<"\n};";

    fout.close();

//    if(showInConsole){
//        cout<<"\n VNF to Node mapping:";
//        for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
//            cout<<"\n";
//            for(unsigned int i=1; i<=finalInstancesCount[f]; i++){
//                cout<<"f"<<f<<char(i+96)<<" -> p:"<<I_VNFINST_2_PN[f][i]<<"\t";
//            }
//        }
//        cout<<"\n Node to VNF mapping:";
//        for(unsigned int p=PhysicalNetwork.srcV; p<=numPNs; p++){
//            cout<<"\np:"<<p<<" -> { ";
//            for(const auto&[f, i]: PN_2_VNF[p])
//                cout<<"f"<<f<<char(i+96)<<", ";
//            cout<<"}";
//        }
//    }
    return true;
}//DeploymentVNF_ScoreMethod

void Simulations::showSimulationTestResults(const SimTEST& result){

    cout<<"\n  "<<result.name<<" : ";
    for(const ServiceFunctionChain& sfc: allSFC){
        cout<<"\nSFC:"<<sfc.index<<" | TrafficRate: "<<sfc.trafficArrivalRate<<" | cntVNFs: "<<sfc.numVNF<<" | PartialChains: "<<sfc.allPartParSFC.size();
        const auto& obj= result.sfcsol.at(sfc.index);
        if(obj.seq_pid==noResDueToNoStg)cout<<"\nNo Result Obtained for Sequential/Parallel due to no instance combination in one of the stage.";
        else{
            cout<<"\n\tSeq. id(";
            if(obj.seq_pid == noResSeq)cout<<" No Result Obtained for Sequential.)";
            else{ cout<<obj.seq_pid<<"):\tD: "<<obj.seq_delay<<"ms  \t";
                for(const auto &blk: sfc.allPartParSFC[obj.seq_pid]) { cout<<" [";  for(const auto& fnid: blk){ cout<<"f"<<fnid<<char(96+obj.seq_fninstmap.at(fnid))<<" "; }   cout<<"]";
                }
            }
            cout<<"\n\tFull Par. id(";
            if(obj.fullpar_pid == noResPar)cout << " No Result Obtained for Partial Parallel.)";
            else{ cout << obj.fullpar_pid << "):\tD: " << obj.fullpar_delay << "ms \t";
                for(const auto &blk: sfc.allPartParSFC[obj.fullpar_pid]) { cout << " [";  for(const auto& fnid: blk){ cout << "f" << fnid << char(96 + obj.fullpar_fninstmap.at(fnid)) << " "; }   cout << "]";
                }
            }
            cout<<"\n\tPartial Par. id(";
            if(obj.ppar_pid == noResPar)cout << " No Result Obtained for Partial Parallel.)";
            else{ cout << obj.ppar_pid << "):\tD: " << obj.ppar_delay << "ms \t";
                for(const auto &blk: sfc.allPartParSFC[obj.ppar_pid]) { cout << " [";  for(const auto& fnid: blk){ cout << "f" << fnid << char(96 + obj.ppar_fninstmap.at(fnid)) << " "; }   cout << "]";
                }
            }
        }

    }
    cout<<"\nDeployement Duration: Seq: "<<result.seq_duration<<"ms." << " | Full Par: " << result.fullpar_duration << "ms." << " | Partial Par: " << result.ppar_duration << "ms.";
    cout<<"\n---------------------------------------------------------";
}

/*! Show Utiliztion of the VNFs
 * @param type 1 for seq, 2 for fully parallel, 3 for partial parallel
 * @param result SimTest Object from which we have to see the utilization.
 */
void Simulations::showVNFsUtilization(const int& type, const SimTEST& result) const{
    if(type == 1){
        cout<<"\nSequential SFC Utilization::\n";
        getUtilization(result.seq_utilization) ;
    } else if(type == 2){
        cout<<"\nFully Parallel SFC Utilization::\n";
        getUtilization(result.fullpar_utilization) ;
    }
    else if(type == 3){
        cout<<"\nPartial Parallel SFC Utilization::\n";
        getUtilization(result.ppar_utilization) ;
    }
}
void Simulations::getUtilization(const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& util)const {
    for (unsigned int f = VNFNetwork.srcVNF; f <=VNFNetwork.numVNF; ++f){
        const VNFNode& vnfInfo = VNFNetwork.VNFNodes.at(f);
        cout<<"f"<<vnfInfo.index<<" ("<<vnfInfo.serviceRate<<"ms):  ";
        for(unsigned int inst = 1; inst <=finalInstancesCount.at(f); inst++){
            cout<<"\t[f"<<vnfInfo.index<<char(96+inst)<<": ";
            if(util.count(vnfInfo.index) and util.at(vnfInfo.index).count(inst))
                cout<< util.at(vnfInfo.index).at(inst) <<" ("<<(util.at(vnfInfo.index).at(inst)/vnfInfo.serviceRate)*100<<" %)";
            else cout<<"0 (0 %)";

            cout<<"]\t";
        }cout<<"\n";
    }
}

/*! Show all the Physical Nodes and their description*/
void Simulations::showPNsDescription() const{
    cout << "\n\n ----- Nodes Description of G(" << PhysicalNetwork.numV<<", "<<PhysicalNetwork.numE<<") ::";
    cout<<"\nIdx |\t"<<"Name | "<<"Cores | "<<"Deg"<<" | "<<"CoeffL"<<" | "<<"VNFs";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----;";
    for (unsigned int u = PhysicalNetwork.srcV; u <= PhysicalNetwork.numV; ++u){
//        if(PNode.find(u) == PNode.end()) continue;
        const PhysicalNode& pn = PhysicalNetwork.PNode.at(u);
        cout<<"\n"<<std::setw(2)<<pn.index<<" | ";
        cout<<std::setw(5)<<pn.name<<" | ";
        cout<<std::setw(2)<<pn.capacity.cores<<" | ";
        cout<<std::setw(2)<<pn.neighbours.size()<<" | ";
        cout<<std::setw(5)<<pn.clusteringcoeff_local<<" | [";
        for(const auto& [f, finst]: PN_2_VNF.at(pn.index))
            cout<<"f"<<f<<char(finst+96)<<", ";
        cout<<"]";
    }
}

/*! Show all the Virtual Network Function nodes and their description */
void Simulations::showVNFsDescription() const{
    cout << "\n\n ----- Virtual Network Functions Description ::";
    cout<<"\nIdx |"<<"Name|\t"<<"ServiceRate | "<<"ExeTime | "<<"Cores | "<<"PNode";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (unsigned int f = VNFNetwork.srcVNF; f <= VNFNetwork.numVNF; ++f){
        const VNFNode& vnfInfo = VNFNetwork.VNFNodes.at(f);
        cout<<std::setw(3)<<"\nf"<<vnfInfo.index<<" | ";
        cout<<std::setw(12)<<vnfInfo.name<<" | ";
//        cout<<std::setw(2)<<"("<<vnfInfo.numInstances<<") | ";
        cout<<std::setw(4)<<vnfInfo.serviceRate<<"ps | ";
        cout<<std::setw(4)<<vnfInfo.executionTime<<"ms | ";
        cout<<std::setw(2)<<"["<<vnfInfo.requirement.cores<<"]";
        for(unsigned int inst = 1; inst <=finalInstancesCount.at(f); inst++){
            if(I_VNFINST_2_PN.at(f).empty()) cout << "PN[-1] ";
            else cout << "PN[" << I_VNFINST_2_PN.at(f).at(inst) << "] ";
        }
    }
}
#endif //SFC_PARALLELIZATION_SIMULATIONCLASS_H
