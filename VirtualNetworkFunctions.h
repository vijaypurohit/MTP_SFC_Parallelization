 //
// Created by vijay on 25-02-2023.
//

#ifndef SFC_PARALLELIZATION_VIRTUALNETWORKFUNCTIONS_H
#define SFC_PARALLELIZATION_VIRTUALNETWORKFUNCTIONS_H

/*!Virtual Network Function Node*/
class VNFNode
{   public:
    unsigned int index  /*! physical node index. Type int.*/;  string name /*! node name  */;
    type_delay serviceRate; ///< rate of service of vnf. in packets per second. arrival rate < service rate.
    type_delay executionTime; ///< time taken to execute the particular function
    NodeCapacity requirement; ///<  requirements  of the VM node. Type NodeCapacity.

    VNFNode(unsigned int& _index, const string& _name, type_delay& _serviceRate, type_delay& _execTime, NodeCapacity _givenRequirements): requirement(_givenRequirements){
        this->index = _index; this->name = std::move(_name);
        this->serviceRate = _serviceRate; this->executionTime = _execTime;
    }
    VNFNode() = default;
    ~VNFNode() = default;
};

/*!Virtual Network Function Collection Data*/
class VirtualNetworkFunctions
{
public:
    string filename_vnf;
    unsigned int numVNF{}  /*!<  number of unique type of VNFs  */; unsigned int srcVNF /*!< starting VNF to iterate the loop*/;
    unordered_map<unsigned int, VNFNode> VNFNodes; ///<Index to VNF structure

     /*! @param _numVirtualNetworkFunctions number of VNFs */
    explicit VirtualNetworkFunctions(unsigned int _numVirtualNetworkFunctions) {
        srcVNF = 1;
        this->numVNF = _numVirtualNetworkFunctions;
    }
    VirtualNetworkFunctions()=default;
    ~VirtualNetworkFunctions()=default;

    void showVNFs_Description() const;

    void Algorithm_NF_Parallelism_Identification();

    void setVNFFnx(const type_delay &T_fnx);
    void setVNFFnxFromRange(const pair<type_delay, type_delay>& funExeTimeRange);
};

void VirtualNetworkFunctions::setVNFFnx(const type_delay& T_fnx){
    for (auto& vnfInfo: VNFNodes){
        vnfInfo.second.executionTime = T_fnx;
    }
}

 void VirtualNetworkFunctions::setVNFFnxFromRange(const pair<type_delay, type_delay>& funExeTimeRange){
     std::mt19937_64 rd_generator(std::chrono::system_clock::now().time_since_epoch().count());
     std::uniform_real_distribution<type_delay> rd_fnexetime_distribution(funExeTimeRange.first, funExeTimeRange.second);
     for (auto& vnfInfo: VNFNodes){
         vnfInfo.second.executionTime = rd_fnexetime_distribution(rd_generator);
     }
 }

/*! Show all the Virtual Machine nodes and their description */
void VirtualNetworkFunctions::showVNFs_Description() const{
    cout << "\n\n ----- Virtual Network Functions Description ::";
    cout<<"\nIdx |"<<"Name|\t"<<"ServiceRate | "<<"ExeTime | "<<"Cores | "<<"PNode";
    cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
    for (unsigned int u = srcVNF; u <= numVNF; ++u){
        const VNFNode& vnfInfo = VNFNodes.at(u);
        cout<<std::setw(3)<<"\nf"<<vnfInfo.index<<" | ";
        cout<<std::setw(12)<<vnfInfo.name<<" | ";
        cout<<std::setw(4)<<vnfInfo.serviceRate<<"ps | ";
        cout<<std::setw(4)<<vnfInfo.executionTime<<"ms | ";
        cout<<std::setw(2)<<"["<<vnfInfo.requirement.cores<<"]";
    }
}


 /*! : For Order(NF1, before, NF2), whether the two NFs are parallelizable, and whether we need to copy packets if the two NFs can be executed in parallel.
     * Green blocks denote parallelizable, no need to copy. \n
     * Orange blocks denote parallelizable, need copy pkts. \n
     * Gray blocks denote not parallelizable situations. \n
     * For read-write or write-write case, we need not copy packets if two NFs modify different fields.
     * @refer Table 2-3, pg 46. Paper NFP: Enabling Network Function Parallelism in NFV
     */
void VirtualNetworkFunctions::Algorithm_NF_Parallelism_Identification(){

     struct pktFields{
         vector<pair<int, int>> info{};
         explicit pktFields(vector<pair<int, int>> x){  info = std::move(x);  }
     };

    enum {R=10, W, RW, T};
    enum {SIP=20, DIP, SPORT, DPORT, Payload, AddRm, Drop};
    enum {NOT_PARALLELIZABLE=30, PARALLELIZABLE_NO_COPY, PARALLELIZABLE_WITH_COPY};

     unordered_map<int, unordered_map<int,int>> DepTable ={
            {R,     {{R, PARALLELIZABLE_NO_COPY},{W, PARALLELIZABLE_WITH_COPY}, {AddRm, PARALLELIZABLE_WITH_COPY},  {Drop, PARALLELIZABLE_NO_COPY}} },
            {W,     {{R, NOT_PARALLELIZABLE},   {W, PARALLELIZABLE_WITH_COPY},  {AddRm, PARALLELIZABLE_WITH_COPY},  {Drop, PARALLELIZABLE_NO_COPY}} },
            {AddRm, {{R, NOT_PARALLELIZABLE},   {W, NOT_PARALLELIZABLE},        {AddRm, PARALLELIZABLE_WITH_COPY},  {Drop, PARALLELIZABLE_NO_COPY}} },
            {Drop,  {{R, NOT_PARALLELIZABLE},   {W, NOT_PARALLELIZABLE},        {AddRm, NOT_PARALLELIZABLE},        {Drop, PARALLELIZABLE_NO_COPY}} }
    };

    int n = 12; // number of VNF
    vector<pktFields*> vnfPktInfo(n+1);
    vnfPktInfo[1] = new pktFields({{SIP, R}, {DIP, R}, {SPORT, R}, {DPORT, R}, {Drop, Drop}});
    vnfPktInfo[2] = new pktFields({{SIP, R}, {DIP, R}, {SPORT, R}, {DPORT, R}, {Payload, R}});
    vnfPktInfo[3] = new pktFields({{SIP, R}, {DIP, R}});
    vnfPktInfo[4] = new pktFields({{SIP, R}, {SIP, W},{DIP, R},{DIP, W}, {SPORT, R}, {DPORT, R}});
    vnfPktInfo[5] = new pktFields({{SIP,     R}, {SPORT,   R}, {Payload, R}});
    vnfPktInfo[6] = new pktFields({{SIP,     R}, {DIP,     R}, {Payload, R},{Payload, W}, {AddRm,   AddRm}});
    vnfPktInfo[7] = new pktFields({{SIP,   R},{SIP,   W},  {DIP,   R},{DIP,   W}, {SPORT, R},{SPORT, W}, {DPORT, R},{DPORT, W}});
    vnfPktInfo[8] = new pktFields({{SIP, R}, {DIP, R},{SIP, W}, {DIP, W}});
    vnfPktInfo[9] = new pktFields({{Payload, R},{Payload, W}});
    vnfPktInfo[10] = new pktFields({});
    vnfPktInfo[11] = new pktFields({{SIP,   R}, {DIP,   R}, {SPORT, R}, {DPORT, R}});
    vnfPktInfo[12] = new pktFields({{SIP,   R}, {DIP,   W}, {SPORT, R}, {DPORT, W}});
 
    // iterating all possible pairs
     unordered_map<unsigned int, unordered_map<unsigned int, unsigned int>> pairs;
    int cnt=0;
    for(int i_vnf=1; i_vnf<=n; i_vnf++){
        for(int j_vnf=1; j_vnf<=n; j_vnf++){
            if(i_vnf==j_vnf)continue; // same vnf id
            bool canBeParallelized = true;
            // for each field in packet 1
            int conflicting_actions = 0;
            for(auto &[al1_field_name, al1_field_val]:  vnfPktInfo[i_vnf]->info){
                for(auto &[al2_field_name, al2_field_val]:  vnfPktInfo[j_vnf]->info) {

                    if((al1_field_val == R and  al2_field_val == W) or (al1_field_val == W and  al2_field_val == W)){
                        if(al1_field_name == al2_field_name)
                            conflicting_actions++;
                        continue; //canBeParallelized unchanged
                    }
                    // if field exist in 2nd packet
                    switch (DepTable[al1_field_val][al2_field_val]) {
                        case NOT_PARALLELIZABLE:
                            canBeParallelized = false;
                            break;
                        case PARALLELIZABLE_NO_COPY:  //do nothing
                            break;
                        case PARALLELIZABLE_WITH_COPY:
                            conflicting_actions++;
                            break;
                    }
                    if(!canBeParallelized) break;
                }
                if(!canBeParallelized) break;
            } // foreach
            if(canBeParallelized){
                cnt++;
                pairs[i_vnf][j_vnf] = 1;
                cout<<"\n"<<cnt<<"["<<i_vnf<<":"<<j_vnf<<"] : canBeParallelized:"<<canBeParallelized<< " : "<<conflicting_actions;
                if(pairs.count(j_vnf) and pairs[j_vnf].count(i_vnf) and pairs[j_vnf][i_vnf]==0)
                    cout<<"  -- found.";
            }else
                pairs[i_vnf][j_vnf] = 0;
        }// for j
    }// for i
}

#endif //SFC_PARALLELIZATION_VIRTUALNETWORKFUNCTIONS_H
