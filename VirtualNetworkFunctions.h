 //
// Created by vijay on 25-02-2023.
//

#ifndef SFC_PARALLELIZATION_VIRTUALNETWORKFUNCTIONS_H
#define SFC_PARALLELIZATION_VIRTUALNETWORKFUNCTIONS_H
/*!
 * Virtual Network Function Node
 * @tparam type_res resource data type.
 */
template <class type_res >
class VNFNode
{
public:
    unsigned int index  /*! physical node index. Type int.*/;  string name /*! node name  */;
    unsigned int numInstances; ///< number of instances of vnf.
    type_delay serviceRate; ///< rate of service of vnf. in packets per second. arrival rate < service rate.
    type_delay executionTime; ///< time taken to execute the particular function
    NodeCapacity<type_res> requirement; ///<  requirements  of the VM node. Type NodeCapacity.

    VNFNode(unsigned int _index, const string& _name,unsigned int _instance, type_delay _serviceRate, type_delay _execTime, NodeCapacity<type_res> _givenRequirements): requirement(_givenRequirements){
        this->index = _index; this->name = _name; this->numInstances = _instance;
        this->serviceRate = _serviceRate; this->executionTime = _execTime;
    }
    ~VNFNode() = default;
};

/*!
 * Virtual Network Function Collection Data
 * @tparam type_res resource data type.
 */
template <class type_res >
class VirtualNetworkFunctions
{
    unsigned int numParallelPairs{}; //< count of number of pairs which are parallel
    unsigned int numVNF{}  /*!<  number of unique type of VNFs  */; unsigned int srcVNF /*!< starting VNF to iterate the loop*/;
public:
    unordered_map<unsigned int, VNFNode<type_res>*> VNFNodes; ///<Index to VNF structure
    unordered_map<unsigned int,unordered_set<unsigned int>> parallelPairs; ///< {i_vnfid ->{j_vnfid it is parallel to}} pairs identifying which are parallel

    unordered_map<unsigned int, unordered_map<unsigned int, unsigned int>> I_VNFinst2VM; ///< VNF {type,inst} is hosted on which VM. {VNFid -> {instid -> VM id}} ie. arr[vnf][inst]=vmid;, inst->1based indexing
//    unordered_map<unsigned int, unordered_map<unsigned int, type_delay>> utilization; ///< utilization of the VNF_Inst. {VNFid -> {instid -> utilization}}
    unordered_map<unsigned int, unordered_map<unsigned int, type_delay>> seq_utilization, par_utilization; ///< utilization of the VNF_Inst. {VNFid -> {instid -> utilization}} for sequential and parallel chain

    /*! @param _numVirtualNetworkFunctions number of VNFs */
    explicit VirtualNetworkFunctions(unsigned int _numVirtualNetworkFunctions) {
        srcVNF = 1;
        this->numVNF = _numVirtualNetworkFunctions;
    }

    ~VirtualNetworkFunctions(){
//        for (unsigned int v = srcVNF; v <= numVNF; ++v) { //deletion of Adj List Edges
//            delete VNFNodes[v]; // delete mapping
////            cout<<VNFNodes[v]->index<<"....";
//        }
        for(const auto vnfInfo: VNFNodes){
            delete vnfInfo.second;
        }
        if(debug) cout<<"\n[VirtualNetworkFunctions Destructor Completed]";
    }

    void showVNFs_Description();
    void assignVNFinstToVM(int, int, int);

    void Algorithm_NF_Parallelism_Identification();

    void findRandomParallelPairs(const string&);

    void showVNFs_Utilization(const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>> &);
};

/*! Show all the Virtual Machine nodes and their description
 * @tparam type_res resource data type.
 */
template <class type_res>
void VirtualNetworkFunctions<type_res>::showVNFs_Description(){
    cout << "\n\n ----- Virtual Network Function Description ::";
    cout<<"\nIndex\t"<<"Name\t"<<"Instances\t"<<"ServiceRate\t"<<"ExeTime\t"<<" Req[cores][memory][disk][speed]\t"<<" VMNode";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (unsigned int u = srcVNF; u <= numVNF; ++u){
        const VNFNode<type_res> *vnfInfo = VNFNodes.at(u);
        cout<<"\n"<<vnfInfo->index<<" |\t";
        cout<<vnfInfo->name<<" |\t";
        cout<<vnfInfo->numInstances<<" |\t";
        cout<<vnfInfo->serviceRate<<" |\t";
        cout<<vnfInfo->executionTime<<" |\t[";
        cout<<vnfInfo->requirement.cores<<"\t"<<vnfInfo->requirement.memory<<"\t"<<vnfInfo->requirement.disk<<"\t"<<vnfInfo->requirement.cpuSpeed<<"]\t";
        for(unsigned int inst = 1; inst <=vnfInfo->numInstances; inst++){
            if(I_VNFinst2VM.at(u).empty()) cout << "VM[-1] ";
            else cout << "VM[" << I_VNFinst2VM.at(u).at(inst) << "] ";
        }
    }
}

/*! Show Utiliztion of the VNFs
 * @param util sequential or parallel
 */
 template <class type_res>
 void VirtualNetworkFunctions<type_res>::showVNFs_Utilization(const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& util){

     for (unsigned int u = srcVNF; u <= numVNF; ++u){
         const VNFNode<type_res> *vnfInfo = VNFNodes.at(u);
         for(unsigned int inst = 1; inst <=vnfInfo->numInstances; inst++){
             cout<<"[f"<<vnfInfo->index<<char(96+inst)<<": ";
             if(util.count(vnfInfo->index) and util.at(vnfInfo->index).count(inst))
                cout<< util.at(vnfInfo->index).at(inst) <<" ("<<(util.at(vnfInfo->index).at(inst)/vnfInfo->serviceRate)*100<<")";
             else cout<<"0 (0)";
             cout<<"]\t";
         }cout<<"\n";
     }
 }
 /*
 * @param vnfIndex VNF index
 * @param vmIndex Virtual Machine Index
 * @param instance Index of the VNF type instance 1,2,3 (1-based indexing)
 * @tparam type_res resource data type. .
 */
template <class type_res>
void VirtualNetworkFunctions<type_res>::assignVNFinstToVM(int vnfIndex, int instance, int vmIndex){
    I_VNFinst2VM[vnfIndex][instance] = vmIndex;
}

/*!
 * Make VNF pairs (uniform (approx 50%)) out of total pairs to be parallelizable/independent so that they can execute parallel. Using Uniform Int Distribution.  \n
 * Total Pairs = numVNF*(numVNF-1) ordered (1,2)(2,1) except (1,1)(2,2)... \n
 * It then writes the pairs to the filename_vnf_parallelpairs {i_vnf -> {set of j_vnf it is parallel to}}
 * @param testDirName  path to the current test directory which consists of inputs files
 */
template <class type_res>
void VirtualNetworkFunctions<type_res>::findRandomParallelPairs(const string& testDirName){

     {
         numParallelPairs = 44;
         parallelPairs={{1,{2,3,7,9,10}},
         {2,{1,4,6,7,10}},
         {3,{4,5,6,7,8}},
         {4,{1,6,3,8,9}},
         {5,{1,2,6,7,8}},
         {6,{2,4}},
         {7,{2,4,6,9,10}},
         {8,{3,7,9}},
         {9,{1,3,7,10}},
         {10,{3,6,9}}};
         if(debug)cout<<"\n\t[Random Parallel Pairs: "<<numParallelPairs<<" | "<<(100.0*numParallelPairs/(numVNF*(numVNF-1)))<<" % out of Total:"<<numVNF*(numVNF-1)<<" pairs]";
     }
     return;

     ofstream fout;
     string filepathExt = output_directory+testDirName+filename_vnf_parallelpairs;///< path to .gv file without extention
     fout.open(filepathExt.c_str(), ios::out);
     if (!fout) {
         string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
         fout.clear();
         throw runtime_error(errorMsg+ __FUNCTION__);
     }

    unsigned seed = chrono::system_clock::now().time_since_epoch().count(); ///< (present time) and clock's epoch
    std::mt19937 rd_generator(seed); ///< mt19937 is a standard mersenne_twister_engine
    std::uniform_int_distribution<unsigned int> rd_distribution(0, 1);      ///< for selecting operation Parallel/Not Parallel
    // iterating all possible pairs

    numParallelPairs=0;
    for(int i_vnf=1; i_vnf<=numVNF; i_vnf++) {
        fout<<"\n{"<<i_vnf<<",{";
        for (int j_vnf = 1; j_vnf <= numVNF; j_vnf++) {
            if(i_vnf == j_vnf)continue;
            if( rd_distribution(rd_generator)>0){
                parallelPairs[i_vnf].insert(j_vnf);
                numParallelPairs++;
                fout<<j_vnf<<",";
            }// 50%
        }//j
        fout<<"}},";
    }// i
    fout.close();
     if(debug)cout<<"\n\t[Random Parallel Pairs: "<<numParallelPairs<<" | "<<(100.0*numParallelPairs/(numVNF*(numVNF-1)))<<" % out of Total:"<<numVNF*(numVNF-1)<<" pairs]";
     /* // to output in console in given map structure
      * cout<<"\n{";
     for(const auto& i_vnf: parallelPairs){
         cout<<"\n  {"<<i_vnf.first<<",{"; bool first = false;
         for(const auto& j_vnf: i_vnf.second){
             if(!first)  { first = true;}
             else cout << ",";
             cout<<j_vnf;
         }
         cout<<"}},";
     }
      cout<<"\n}";*/
}

 /*! : For Order(NF1, before, NF2), whether the two NFs are parallelizable, and whether we need to copy packets if the two NFs can be executed in parallel.
     * Green blocks denote parallelizable, no need to copy. \n
     * Orange blocks denote parallelizable, need copy pkts. \n
     * Gray blocks denote not parallelizable situations. \n
     * For read-write or write-write case, we need not copy packets if two NFs modify different fields.
     * @refer Table 2-3, pg 46. Paper NFP: Enabling Network Function Parallelism in NFV
     */
template <class type_res>
void VirtualNetworkFunctions<type_res>::Algorithm_NF_Parallelism_Identification(){

     struct pktFields{
         unordered_map<int, int> info{};
         explicit pktFields(unordered_map<int, int> x){  info = std::move(x);  }
     };

    enum {R=1, W=2, RW=2, T=3};
    enum {SIP=1, DIP, SPORT, DPORT, Payload, AddRm, Drop};
    enum {NOT_PARALLELIZABLE=-1, PARALLELIZABLE_NO_COPY=1, PARALLELIZABLE_WITH_COPY=2};

    const unordered_map<int, unordered_map<int,int>> DepTable ={
            {R,     {{R, PARALLELIZABLE_NO_COPY},{W, PARALLELIZABLE_NO_COPY}, {AddRm, PARALLELIZABLE_NO_COPY},  {Drop, PARALLELIZABLE_NO_COPY}} },
            {W,     {{R, NOT_PARALLELIZABLE},   {W, PARALLELIZABLE_NO_COPY},  {AddRm, PARALLELIZABLE_NO_COPY},  {Drop, PARALLELIZABLE_NO_COPY}} },
            {AddRm, {{R, NOT_PARALLELIZABLE},   {W, NOT_PARALLELIZABLE},        {AddRm, PARALLELIZABLE_NO_COPY},  {Drop, PARALLELIZABLE_NO_COPY}} },
            {Drop,  {{R, NOT_PARALLELIZABLE},   {W, NOT_PARALLELIZABLE},        {AddRm, NOT_PARALLELIZABLE},        {Drop, PARALLELIZABLE_NO_COPY}} }
    };

    int n = 12; // number of VNF
    vector<pktFields*> vnfPktInfo(n+1);
    vnfPktInfo[1] = new pktFields({{SIP, R}, {DIP, R}, {SPORT, R}, {DPORT, R}, {Drop, T}});
    vnfPktInfo[2] = new pktFields({{SIP, R}, {DIP, R}, {SPORT, R}, {DPORT, R}, {Payload, R}});
    vnfPktInfo[3] = new pktFields({{SIP, R}, {DIP, R}});
    vnfPktInfo[4] = new pktFields({{SIP, RW}, {DIP, RW}, {SPORT, R}, {DPORT, R}});
    vnfPktInfo[5] = new pktFields({{SIP,     R}, {SPORT,   R}, {Payload, R}});
    vnfPktInfo[6] = new pktFields({{SIP,     R}, {DIP,     R}, {Payload, RW}, {AddRm,   T}});
    vnfPktInfo[7] = new pktFields({{SIP,   RW},  {DIP,   RW}, {SPORT, RW}, {DPORT, RW}});
    vnfPktInfo[8] = new pktFields({{SIP, RW}, {DIP, RW}});
    vnfPktInfo[9] = new pktFields({{Payload, RW}});
    vnfPktInfo[10] = new pktFields({});
    vnfPktInfo[11] = new pktFields({{SIP,   R}, {DIP,   R}, {SPORT, R}, {DPORT, R}});
    vnfPktInfo[12] = new pktFields({{SIP,   R}, {DIP,   W}, {SPORT, R}, {DPORT, W}});
 
    // iterating all possible pairs
    int cnt=0;
    for(int i_vnf=1; i_vnf<=n; i_vnf++){
        for(int j_vnf=1; j_vnf<=n; j_vnf++){
            if(i_vnf==j_vnf)continue; // same vnf id
            bool canBeParallelized = true;
            // for each field in packet 1
            for(auto &[al1_field_name, al1_field_val]:  vnfPktInfo[i_vnf]->info){
                // if field exist in 2nd packet
                if(vnfPktInfo[j_vnf]->info.count(al1_field_name)){
                    int al2_field_val = vnfPktInfo[j_vnf]->info[al1_field_name];
                        switch (DepTable[al1_field_val][al2_field_val]) {
                            case NOT_PARALLELIZABLE: canBeParallelized = false; break;
                            case PARALLELIZABLE_NO_COPY: continue;
                            case PARALLELIZABLE_WITH_COPY:
                                canBeParallelized = false; break;
                                break;
                        }
                }
                if(!canBeParallelized) break;
            } // foreach
            if(canBeParallelized)cnt++;
            cout<<"\n"<<cnt<<"["<<i_vnf<<":"<<j_vnf<<"] : canBeParallelized:"<<canBeParallelized;
        }// for j
    }// for i
}

#endif //SFC_PARALLELIZATION_VIRTUALNETWORKFUNCTIONS_H
