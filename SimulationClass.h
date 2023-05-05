//
// Created by vijay on 20-04-2023.
//

#ifndef SFC_PARALLELIZATION_SIMULATIONCLASS_H
#define SFC_PARALLELIZATION_SIMULATIONCLASS_H

class Simulations{
public:
    string sim_name{};
    string dirName{}/*!< directory name */, fullDirName{} /*!< directory name with slash at the end*/;
    /// Common Object
    PhysicalGraph PhysicalNetwork; ///< Graph Object contains network related functions
    VirtualNetworkFunctions VNFNetwork; ///< VNF Object contains VNF (Virtual Network Function) related code.

    /// SFC Related
    string filename_sfc{};///< name of the file from where SFCs are read;
//    unordered_map<unsigned int, unsigned int> sfcIndex_2_arrIndex;
    vector<ServiceFunctionChain> allSFC; ///< SFCs object contains code related to Service function chains
    vector<ServiceFunctionChain*> sortedSFCs; ///< SFCs sorted for deployement in order of their priority of traffic rate and length
    unordered_map<unsigned int, unordered_set<unsigned int>> VNF_2_SFC; ///< VNFs Type are present in what SFCs Requests (sfc vector idx). unordered set for each check

    /// VNF related
    unsigned int numParallelPairs{0}, totalPairs{0};
    unordered_map<unsigned int, unordered_map<unsigned int, unsigned int>> parallelPairs; ///< {i_vnfid ->{j_vnfid it is parallel with copy/without copy}} pairs identifying which are parallel
    type_delay likelihood_mean; ///< Mean parallelism likelihood of all funtion pairs (fi,fj);
    unordered_map<unsigned int, unordered_map<unsigned int, double>> likelihood;/*!< Likelihood of parallelizing two functions fi and fj defined as L_{fi,fj} = |num of sfc in which fi and fj are present and are parallelizable| divided by |numOfSFCS|*/

    ///VNF Deployement Mappings
    type_delay total_service_capacity = 0;
    unordered_map<unsigned int, unsigned int> finalInstancesCount; ///< count of instances of each function
    unordered_map<unsigned int, vector<pair<unsigned int,unsigned int>>> PN_2_VNF;
    unordered_map<unsigned int, unordered_map<unsigned int, unsigned int>> I_VNFINST_2_PN; ///< VNF {type,inst} is hosted on which PN. {VNFid -> {instid -> PN id}} ie. arr[vnf][inst]=pnid;, inst->1based indexing

    /// Test
    unordered_map<string, SimTEST> TestsResult;
    unordered_map<unsigned int, vnfDelaysPreComputed> vnfDelays;///< pre-calculated VNF delays (processing, execution and queuing delay). Before iterating all the partial par sfc it is better to calculate it for each chain as they will be same for each chain.
    int sfccompleted{1};
    /*! @param simulation_name Name to identify the simulation instance
     * @param directory directory on which it is running.
     */
    Simulations(string simulation_name, string directory){
        this->sim_name = std::move(simulation_name);
        this->dirName = std::move(directory);
        this->fullDirName = dirName + "/";
        createDirectory(this->fullDirName);
    };

    void showPNsDescription() const;
    void showVNFsDescription() const;
    void showSFCsDescriptions() const;
    void showVNFsUtilization(const int &, const SimTEST&) const;
    void getUtilization(const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>&) const;
    void showSimulationTestResults(const SimTEST &result);

    bool readGenericServiceFunctionsChains(const string &, const std::string& sortOpt, const int&);
    bool readDataFromFileInit(const string &, const string &, const string &, const pair<float, int>& = {50,2}, const string& = "asc_length");
    bool writeInFileSimulationTestResults(const int& observation , const _Ios_Openmode& filemode, const string& other);

    bool findRandomParallelPairs(const float &, int = 0);
    bool calcLikelihoodOfTwoFunctions(const int& opt=1, const int &showInConsole=dontShow);
    bool DeploymentVNF_ScoreMethod(const float& fnInstScalingFactor=0.5,const float& alphaDecayRate=0.75, const int& choice_pxs = 1, const int& choice_dist_pref = 1, const int& showInConsole = dontShow);


    void findAllPartialSFC_Backtrack(unsigned int cur_level_idx,unordered_map<unsigned int, vector<unsigned int>> &levelInfo,vector<vector<unsigned int>> &partSFC, unsigned int &mask, ServiceFunctionChain& csfc);
    vector<unsigned int> findParallelSubset(const unsigned int &seqid, const unsigned int &SZ, const vector<unsigned int> &vnfSeq);
    void convert_SeqSFC_to_FullParallel(ServiceFunctionChain &csfc);
    void convert_fullParVNFBlk_to_AllPartialChains(ServiceFunctionChain &csfc, const int &showInConsole = dontShow);
    void convert_SeqSFC_to_SubsetPartialChains(ServiceFunctionChain &csfc, const int &showInConsole = dontShow);

    void printGenericNetworkGraph() const;
};

/*  **************************************************************************************************************** */
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
        if(finalInstancesCount.count(f)){
            for(unsigned int inst = 1; inst <=finalInstancesCount.at(f); inst++){
                if(I_VNFINST_2_PN.at(f).empty()) cout << "PN[-1] ";
                else cout << "PN[" << I_VNFINST_2_PN.at(f).at(inst) << "] ";
            }
        }
    }
}
/*!
 * Show description of all the SFCs in the Network
 */
void Simulations::showSFCsDescriptions() const{
    cout << "\n\n ----- Service Functions Chains Description ::\n";
    cout<<"Id  | "<<"Name |\t"<<"Len |\t"<<"ArrivalRate |/t"
        <<"[#subsetPartial "<<"#allPartial] |\t"
        <<"(ps->pd)\t"<<"Sequential SFC |\t"<<"Fully Parallel SFC";
    cout<<"\n-----\t -----\t -----\t -----\t -----\t -----\t -----\t -----\t -----\t -----";
    for (const auto &sfc: allSFC){
        cout<<"\n"<<std::setw(2)<<sfc.index<<" | " <<std::setw(7)<<sfc.name<<" | " <<std::setw(2)<<sfc.numVNF<<" | " <<std::setw(4)<<sfc.trafficArrivalRate<<" | [";
        cout<<std::setw(2)<<sfc.subsetPartParSFC.size()<<" | "<<std::setw(2)<<sfc.allPartParSFC.size()<<"] | ";
        cout<<std::setw(7)<<"("<<sfc.access_nodes.first<<"->"<<sfc.access_nodes.second<<") | ";
        cout << "("; for(const auto& fn : sfc.vnfSeq){ cout <<"f"<< fn ;   cout<< "; "; } cout << ")\t";
        cout << "("; for(const auto& blk: sfc.vnfBlocksPar){ cout<<" ["; for(int fn: blk){ cout <<"f"<< fn ;  cout<<"; "; }  cout<<"] "; } cout << ")";
    }
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
    type_delay tcap=0;
    type_delay trate=0;
    for (unsigned int f = VNFNetwork.srcVNF; f <=VNFNetwork.numVNF; ++f){
        const VNFNode& vnfInfo = VNFNetwork.VNFNodes.at(f);
        cout<<"f"<<vnfInfo.index<<" ("<<vnfInfo.serviceRate<<"ms):  ";
        for(unsigned int inst = 1; inst <=finalInstancesCount.at(f); inst++){
            cout<<"\t[f"<<vnfInfo.index<<char(96+inst)<<": ";
            if(util.count(vnfInfo.index) and util.at(vnfInfo.index).count(inst))
                cout<< util.at(vnfInfo.index).at(inst) <<" ("<<(util.at(vnfInfo.index).at(inst)/vnfInfo.serviceRate)*100<<" %)";
            else cout<<"0 (0 %)";

            cout<<"]\t";
            trate += util.at(vnfInfo.index).at(inst);
            tcap += vnfInfo.serviceRate;
        }cout<<"\n";
    }
    cout<<"System Load: "<<(trate/tcap)*100.0<<"%";
}

void Simulations::showSimulationTestResults(const SimTEST& result){
    cout<<"\n["<<result.name<<"]";
    for(const ServiceFunctionChain& sfc: allSFC){
        cout<<"\nSFCid:"<<sfc.index<<" | Len: "<<sfc.numVNF<<" | ArrivalRate: "<<sfc.trafficArrivalRate
                << " | #subsetC: " << sfc.subsetPartParSFC.size()
                <<"  | #allC: "<<sfc.allPartParSFC.size()
                << " | (ps:"<<sfc.access_nodes.first<<" -> "<< " pd:"<<sfc.access_nodes.second<<")";

        if(result.sfcsol.count(sfc.index) == 0) continue;
        const auto& sfcres= result.sfcsol.at(sfc.index);

        cout<<"\n\tSequential:";
        if(sfcres.seq_status==noResDueToNoStg)cout<<" No Result obtained due to no instance combination in one of the stage.";
        else if(!sfcres.seq_fninstmap.empty()){
            cout<<"\tD: "<<sfcres.seq_delay<<" ("<<sfcres.seq_load<<"%) \t(";
                for(const auto& fnid : sfc.vnfSeq){ cout<<"f"<<fnid<<char(96+sfcres.seq_fninstmap.at(fnid))<<"; "; } cout<<")";
        }//if(obj.seq_status==noResDueToNoStg)
        else cout<<" -- ";

        cout<<"\n\tFully Parallel:";
        if(sfcres.fullpar_status==noResDueToNoStg)cout<<" No Result obtained due to no instance combination in one of the stage.";
        else if(!sfcres.fullpar_fninstmap.empty()){
            cout<<"\tD: " << sfcres.fullpar_delay<<" ("<<sfcres.fullpar_load<<"%) \t";
            for(const auto& blk: sfc.vnfBlocksPar){  cout<<" ["; for(int fnid: blk){ cout << "f" << fnid << char(96 + sfcres.fullpar_fninstmap.at(fnid)) << " "; } cout<<"]"; }
        }else cout<<"  -- ";

        cout<<"\n\tPartial-Parallel";
        if(sfcres.ppar_pid==noResDueToNoStg)cout<<" No Result obtained due to no instance combination in one of the stage.";
        else if(!sfcres.ppar_fninstmap.empty()){
            cout<<"["<< sfcres.ppar_pid << "]:\tD: " << sfcres.ppar_delay<<" ("<<sfcres.ppar_load<<"%) \t";
            for(const auto &blk: (*sfc.partialParallelChains)[sfcres.ppar_pid]) { cout << " [";  for(const auto& fnid: blk){ cout << "f" << fnid << char(96 + sfcres.ppar_fninstmap.at(fnid)) << " "; }   cout << "]";}
        }else cout<<":  -- ";
    }//all SFC
    cout<<"\nTotal Delay: [Seq: "<<result.total_seq_delay<<"] [Fully-Par: " << result.total_fullpar_delay << "] [Partial-Par: " << result.total_ppar_delay<<"]";
    cout<<"\nDeployement Duration: [Seq: "<<result.total_seq_duration<<"] [Fully-Par: " << result.total_fullpar_duration << "] [Partial-Par: " << result.total_ppar_duration<<"]";
    cout<<"\nSysyem Load: [Seq: "<<(result.total_seq_load/total_service_capacity)*100.0<<"%] [Fully-Par: " <<(result.total_fullpar_load/total_service_capacity)*100.0 << "%] [Partial-Par: " << (result.total_ppar_load/total_service_capacity)*100.0<<"%]";
    cout<<"\n---------------------------------------------------------";
}

/*   **************************************************************************************************************** */

/*! Function will read the SFC chains data from the file in the input directory and read into ServiceFunctionChain class.\n\n
 * It will simultaneously read vnf data into sequential vector and also construct adj list. \n\n
 * First Line in file consist of numOfSFC S\n
 * from next line upto S times, each s group consist\n
 * sfc_index  (arrivalRate of SFC) (access nodes src dst) (numOfVNFs present except src and dest)   \n
 * This line contains VNFs Type id (for i=1 to numOFVNFs)
 * @param filename_sfc Name of the SFCs file which consist of the data
 * @param sortOpt how to sort the SFCs, asc_length, dsc_length, dsc_rate.
 * @param showinConsole output in console (0 don't show, 1 final output, 2 detailed output)
 */
bool Simulations::readGenericServiceFunctionsChains(const std::string& filename, const std::string& sortOpt="asc_length", const int& showInConsole = dontShow) {
    allSFC.clear(), sortedSFCs.clear(), VNF_2_SFC.clear();
    filename_sfc = filename.substr(0,filename.find('.'));
    ifstream fin;
    string filepathExt = input_directory + fullDirName + filename;
    fin.open(filepathExt.c_str(), ios::in);
    if(!fin) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fin.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }


    if(debug)cout<<"\n[Reading Original Sequential SFCs] SFC File: "<<filepathExt<<endl;
    unsigned int i_nSFC; ///< input num of SFCs
    if(fin>>i_nSFC){
        allSFC = std::move(vector<ServiceFunctionChain>(i_nSFC));
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
            sortedSFCs.push_back(&allSFC[ni]);
            if(showInConsole == showDetailed) {
                cout << readSFC_idx << "-" << readSFC_name << "-" << readSFC_totalVNF << "-" << readSFC_arrivalRate<<"- [";
                for(const auto& x: allSFC[ni].vnfSeq) cout<<x<<" -> "; cout<<"]\n"; }
        }
        if(debug)cout<<"\tSFCs:"<<allSFC.size();
    }else{
        string errorMsg = "Invalid Input File Values. File "+filepathExt+ ". Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fin.close();

    if(sortOpt == "asc_length")
        sort(sortedSFCs.begin(), sortedSFCs.end(), comparator_sfc_asc_length);
    else if(sortOpt == "dsc_length")
        sort(sortedSFCs.begin(), sortedSFCs.end(), comparator_sfc_dsc_length);
    else if(sortOpt == "dsc_rate")
        sort(sortedSFCs.begin(), sortedSFCs.end(), comparator_sfc_dsc_rate);
    else {
        cerr<<"\n Invalid sort option.";
        return false;
    }

    return true;
}

/*!
 * Reading Data First Time
 * @param fileNetwork FileName of the Network present in the directory.
 * @param fileVNF FileName of the VNFs to read VNFs data.
 * @param fileSFC fileName of the SFCs to read SFCs data.
 * @param parallelPairsOpt how you want parallel pairs to be generated. {threshold, Option}.
 * threshold (1-99):percentage of pairs we want to parallelize. (1-99).
 * Option:way to obtain parallel pairs, 0 -> fixed, 1->approx around threshold, 2->exact to threshold (1-85%)
 * @param sfc_sort_opt how to sort the SFCs, asc_length, dsc_length, dsc_rate.
 * @return status
 */
bool Simulations::readDataFromFileInit(const string& fileNetwork, const string& fileVNF, const string& fileSFC, const pair<float, int>& parallelPairsOpt, const string& sfc_sort_opt){
    if( readNetwork(fullDirName,fileNetwork, PhysicalNetwork) and
        readVirtualNetworkFunctions(fullDirName,fileVNF, VNFNetwork) and
        readGenericServiceFunctionsChains(fileSFC, sfc_sort_opt) and
        findRandomParallelPairs(parallelPairsOpt.first,parallelPairsOpt.second) and /// based on VNFs
        calcLikelihoodOfTwoFunctions() ///based on parallel pairs and VNFs in SFC
        ){
        for(ServiceFunctionChain& sfc: allSFC){
            convert_SeqSFC_to_FullParallel(sfc);
            convert_fullParVNFBlk_to_AllPartialChains(sfc);
            convert_SeqSFC_to_SubsetPartialChains(sfc);
            sfc.partialParallelChains = &sfc.subsetPartParSFC;
        }
        if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";
    }
    else return false;
    return true;
}

/*!
 * Write the Simulation Result into csv file.
 */
bool Simulations::writeInFileSimulationTestResults(const int& observation = 1, const _Ios_Openmode& filemode = ios::out, const string& otherval = ""){
    bool allowHeader =false;
    ofstream fout;
    string filepathExt = output_directory+fullDirName+"ST_"+sim_name+"-Result.csv";///< path to .gv file without extention
    fout.open(filepathExt.c_str(), filemode);
    if (!fout) {
        fout.clear();
        filepathExt += "-Result-Renamed-"+to_string(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count())+".csv";///< path to .gv file without extention
        fout.open(filepathExt.c_str(), ios::out);
        if(!fout){
            fout.clear();
            std::cerr <<"\n"<<__FUNCTION__ <<" caught: "+filepathExt+ " failed to open.";
            return false;
        }
        allowHeader = true;
    }
    if(filemode != ios::app or allowHeader){ /// if not appened mode
        /// observation idx | sfc detial | seq res | full par res | partial par res | seq duration | full par duration | partial par duration | sfc
        fout<<"Observation, sfcID, Len, ArrivalRate, numSubsetChains, numPartialChains"; ///< HEADER LINE
        for(const auto& _: TestsResult){
            fout<<","<<_.first+", SeqDelay, FullyParDelay, PartParDelay,  SeqLoad, FullyParLoad, PartParLoad, SeqDuration, FullyParDuration, PartParDuration,  SeqRes, FullyParRes, PartParRes, FullyParDeg, PartParDeg";
        }
    }
    if(!otherval.empty())fout<<"\n"<<otherval<<",";
    fout << "\n" << observation<<",";
    for(const ServiceFunctionChain*const& sfcpointer: sortedSFCs) {
        const ServiceFunctionChain& sfc = *sfcpointer;
        fout<<sfc.index<<","<<sfc.numVNF<<", "<<sfc.trafficArrivalRate<<","<<sfc.subsetPartParSFC.size()<<","<<sfc.allPartParSFC.size();

        for(auto&[_,objST]: TestsResult){
            if(objST.sfcsol.count(sfc.index) == 0) continue;
            const auto& sfcres= objST.sfcsol.at(sfc.index);
            fout<<",->,";
            fout<<sfcres.seq_delay<<","<<sfcres.fullpar_delay<<","<<sfcres.ppar_delay<<",";
            fout<<sfcres.seq_load<<","<<sfcres.fullpar_load<<","<<sfcres.ppar_load<<",";
            fout<<sfcres.seq_duration<<","<<sfcres.fullpar_duration<<","<<sfcres.ppar_duration;

            if(sfcres.seq_fninstmap.empty()) fout<<",";
            else{ fout<<",Seq:("; for(const auto& fnid : sfc.vnfSeq){ fout<<"f"<<fnid<<char(96+sfcres.seq_fninstmap.at(fnid))<<"; "; } fout<<")"; }

            unsigned int degOfParallelism_fullpar = 0; ///< the maximal number of parallelized NFs in any step
            if(sfcres.fullpar_fninstmap.empty()) fout<<",";
            else{ fout<<",FullPar:";
                 for(const auto& blk: sfc.vnfBlocksPar){  fout<<" ["; for(int fnid: blk){ fout << "f" << fnid << char(96 + sfcres.fullpar_fninstmap.at(fnid)) << " "; } fout<<"]";
                     if(blk.size() > degOfParallelism_fullpar) degOfParallelism_fullpar = blk.size(); }
            }

            unsigned int degOfParallelism_ppar = 0; ///< the maximal number of parallelized NFs in any step
            if(sfcres.ppar_fninstmap.empty()) fout<<",";
            else{ fout<<",PartPar("<< sfcres.ppar_pid << "):";
                for(const auto &blk: (*sfc.partialParallelChains)[sfcres.ppar_pid]) { fout << " [";  for(const auto& fnid: blk){ fout << "f" << fnid << char(96 + sfcres.ppar_fninstmap.at(fnid)) << " "; }   fout << "]";
                    if(blk.size() > degOfParallelism_ppar) degOfParallelism_ppar = blk.size(); }
            }
            fout<<","<<degOfParallelism_fullpar<<","<<degOfParallelism_ppar;
        }
        fout<<"\n"<<",";
    }
    fout<<",,,,,";
    for(const auto& [nameST,objST]: TestsResult){
        fout<<nameST<<","<<objST.total_seq_delay<<","<<objST.total_fullpar_delay<<","<<objST.total_ppar_delay<<",";
        fout<<(objST.total_seq_load/total_service_capacity)*100.0<<","<<(objST.total_fullpar_load/total_service_capacity)*100.0<<","<<(objST.total_ppar_load/total_service_capacity)*100.0<<",";
        fout<<objST.total_seq_duration<<", "<< objST.total_fullpar_duration <<", "<<objST.total_ppar_duration<<", , , , , ,";
    }

    fout.close();
    return true;
}

/* SFC Related **************************************************************************************************************** */

/*! Helper backtrack function to find all partial sfc corresponding to cluster_i and fullParVNFBlocks.
 * @param cur_level_idx current level on which we are checking function combination
 * @param levelInfo  all the level wise info, what functions present and its combination present in the level
 * @param partSFC partial sfc in current iteration
 * @param mask visited mask of the function used
 * @param allPartParSFC all the partial parallel sfc to be collected
 * @param fullParVNFBlocks full parallel version of the sfc.
 */
void Simulations::findAllPartialSFC_Backtrack(unsigned int cur_level_idx, unordered_map<unsigned int,vector<unsigned int>>& levelInfo, vector<vector<unsigned int>>&partSFC, unsigned int& mask, ServiceFunctionChain& csfc)
{
    // if current level is equal to total level in cluster/SFC.
    if(cur_level_idx == levelInfo.size()) {
        csfc.allPartParSFC.push_back(partSFC);
        return;
    }// base case

    unsigned int blkid = levelInfo[cur_level_idx][0], nl = levelInfo[cur_level_idx][1], kl= levelInfo[cur_level_idx][2];
    for(const vector<unsigned int>& curCombination: nCk[nl][kl]){ //n=4,k=1 {{1}, {2}, {3}, {4}} } | k=2 {{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}} }

        bool thisCombinationCanBeVisited = true;
        vector<unsigned int> curBlkFunc; //< if we can visit this combination then this vector will be current blk
        unsigned int localMask=0; /// how many visited in this stage

        for(const unsigned int& idx: curCombination){ // check if the curCombination idx{2,3} mapped to block func {fw (id 1-1), fx (id 2-1), fy (id 3-1)} --> fx, fy can be visited or its node already visited.
            int fn_id = csfc.vnfBlocksPar[blkid][idx-1];
            if((mask & (1<<fn_id)) != 0){ thisCombinationCanBeVisited=false; break; }
            curBlkFunc.push_back(fn_id);
            localMask |= (1<<fn_id);
        }
        if(thisCombinationCanBeVisited){
            mask |= localMask;
            partSFC.push_back(std::move(curBlkFunc));
            findAllPartialSFC_Backtrack(cur_level_idx+1, levelInfo, partSFC, mask, csfc);
            mask ^= localMask;
            partSFC.pop_back();
        }
    }//recursive case
};

/*! It reads data from sequential vector (csfc.vnfSeq) and then convert into full Parallel SFC (parallel VNF blocks) chain by detecting
 * whether two function can be parallelised or not  (parallelPairs).\n
 * Each block denote fully parallel VNFs in that step.
 * Please not that functions are not mutually parallel f1->f2, f2->f3 it is not necessary that f1->f3 would be parallel.
 * @param[in] csfc SFC object whose fully parallel version we have to find.
 */
void Simulations::convert_SeqSFC_to_FullParallel(ServiceFunctionChain& csfc){//convert_SeqSFC_to_FullParallel
    const unsigned int& SZ = csfc.vnfSeq.size(); // total vnfs including src and dest

    unsigned int cid=0;
    csfc.vnfBlocksPar.push_back({csfc.vnfSeq[cid]}); // stage 0, pushing src node
    for(cid=1; cid<SZ; cid++){ // from SFCsrc+1 stg to SFCdst-1 stage
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
}//convert_SeqSFC_to_FullParallel

/*! generate all the feasible partial parallel SFC for the given full parallel SFC.
 * @param[in] csfc SFC whose partial version we have to find.
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 * @example fullParVNFBlocks={{1},{2,3,4,5}}, nVNFs:5 (f1,f2,f3,f4,f5).\n
 * Blk(0) = {1} only 1 parallel function, Blk(1) = {2,3,4,5} all 4 are parallel.\n
 * clusterSz[5(size=nVNFs)] = {{1, 2,1,1},{1,4},{3,1,1} ... } so on. cluster_i[l] denotes number of function parallel in that level. \n
 * where cluster_i {1, 2, 1, 1} means in level[0] only one function runs, level[2] = 2 functions run together, and level[3] and level[4] one-one function are there \n
 * {[1], [2,1,1]} is mapped to {[f1], [f2,f3,f4,f4]} such that [ 1-c combination of f1 ]-> [2-c combination of f2,f3,f4,f5] -> [1-c combination of f2,f3,f4,f5] -> [1-c combination of f2,f3,f4,f5] \n
 */
void Simulations::convert_fullParVNFBlk_to_AllPartialChains(ServiceFunctionChain& csfc, const int& showInConsole){
    const vector<vector<unsigned int>>& fullParVNFBlocks = csfc.vnfBlocksPar;
    const unsigned int& nBlk = fullParVNFBlocks.size(); ///< number of blocks of the fully parallel SFC, (including src and dst block)
//    const vector<vector<unsigned int>> SK = { {1,  2,1,1},{1, 2,2}, {1,4}, {2,2,1}, {3,1,1} };

    for(const vector<unsigned int>& cluster: integerCompositions[csfc.numVNF]){ //
        if(cluster[0] > fullParVNFBlocks[0].size()) // if in first block(after src blk) 2 function is there, but cluster saying 3 needs to be parallel then continue next cluster.
            continue;

        unsigned int cur_level=0; ///< current level at which we have to insert all the nodes
        unordered_map<unsigned int,vector<unsigned int>> levelInfo;///< this stores level wise info {level i -> {0 -> blkId, 1->n, 2-> k}} to find nCk for block blkId of parSFC_Full

        if(showInConsole == showDetailed){
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
                if(showInConsole == showDetailed){ cout<<"\tNot feasible for block{"; for(const auto& x: fullParVNFBlocks[blk_id]) cout<<"f"<<x<<","; cout<<"}";}
                allBlksOfSFCDone = false; // break lag gya isliye dfs call mt krna
                break;
            }

            for(unsigned int li=cur_level; li<cur_level+delta; li++){
                if(showInConsole == showDetailed){
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
            findAllPartialSFC_Backtrack(0, levelInfo, partSFC, mask, csfc);
        }
    }

    if(showInConsole == showFinal){
        cout<<"\nTotal PartSFC:"<<csfc.allPartParSFC.size();
        for(int idx=0; idx<csfc.allPartParSFC.size(); idx++){
            const auto& PCs = csfc.allPartParSFC[idx];
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

/*! Helper Function. It check if we can parallelise the function from index (seqid) upto SZ and return that array.
 * @param seqid fun index in vnfSeq array
 * @param SZ cluster size, this many we want parallel in current stage
 * @param vnfSeq orginal sequential sfc
 * @return vector of fn of size SZ, which are all parallel
 */
vector<unsigned int> Simulations::findParallelSubset(const unsigned int& seqid, const unsigned int& SZ, const vector<unsigned int>& vnfSeq){
    unsigned int qi = seqid;
    vector<unsigned int> part = {vnfSeq[qi]};
    for(qi = seqid+1; qi<seqid+SZ; qi++){/// vnfseq me seqid se leke seqid+sz tk check kro ki parallel hai ya nhi
        for(const unsigned int& prv_vnf: part){
            if(parallelPairs.count(prv_vnf)==0 or parallelPairs.at(prv_vnf).find(vnfSeq[qi]) == parallelPairs.at(prv_vnf).end()){
                return {};
            }
        }
        part.push_back(vnfSeq[qi]);
    }//qi
    return part;
}//findSubsetPartialSFC

/*!
 * It converts Sequential SFC into partial chains based on the Cluster Enumeration. If they are parallel according to cluster then it constructs one of the partial sfc
 * @param[in] csfc SFC object whose subset partial parallel version we have to find.
 * @param[in] parallelPairs to check whether two functions are parallel or not.
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 * @example vnfSeq={{1,2,3}, Len:3 \n
 * Blk(0) = {1} only 1 parallel function, Blk(1) = {2,3,4,5} all 4 are parallel.\n
 * clusterSz[3(size=Len)] = {{1,1,1},{2,1},{1,2},{3} } \n
 * where cluster_i {1,2} means in level[0] only one function runs, level[1] = check 2 functions can run together. if it so than it is one of the parallel chain.
 */
void Simulations::convert_SeqSFC_to_SubsetPartialChains(ServiceFunctionChain& csfc, const int& showInConsole){//convert_SeqSFC_to_SubsetPartialChains
    vector<vector<vector<unsigned int>>>& subsetPartParSFC = csfc.subsetPartParSFC;

    for(const vector<unsigned int>& cluster: integerCompositions[csfc.vnfSeq.size()]){ /// how we can arrange the seq sfc
        vector<vector<unsigned int>>partSFC; ///< partial sfc in current iteration
        unsigned int vnfdone = 0; ///< in a vnfSeq how many functions are completed from start
        for(int cid = 0; cid<cluster.size(); cid++){
            vector<unsigned int> part =  std::move(findParallelSubset(vnfdone, cluster[cid], csfc.vnfSeq));
            if(part.empty()){ break; }
            vnfdone += part.size();
            partSFC.push_back(std::move(part));
        }
        if(vnfdone == csfc.vnfSeq.size()) subsetPartParSFC.push_back(partSFC);
        else continue;
        if(showInConsole == showDetailed){
            cout<<"\ncluster["; for(const auto& x: cluster) cout<<x<<" "; cout<<"]\t(";
            for(const auto& blks: partSFC){
                cout<<"["; for(auto fn_id: blks){ cout<<"f"<<fn_id<<" ";  } cout<<"]";
            } cout<<")";
        }
    }// each cluster

    if(showInConsole == showFinal){
        cout<<"\nTotal PartSFC:"<<subsetPartParSFC.size();
        for(int idx=0; idx<subsetPartParSFC.size(); idx++){
            const auto& PCs = subsetPartParSFC[idx];
            cout<<"\nPC["<<idx+1<<"] ( "; unordered_map<unsigned int,unsigned int> freq;
            for(const auto& blks: PCs){
                cout<<"["; for(auto fn_id: blks){
                    cout<<"f"<<fn_id<<" ";
                    if(++freq[fn_id]>1)  throw runtime_error("Error in calculation of allPartSFC. Some fn_id repeated");
                } cout<<"]";
            } cout<<")";
        }
    }
}//convert_SeqSFC_to_SubsetPartialChains

/* **************************************************************************************************************** */
/*!
 * Make VNF pairs (uniform (approx 50%)) out of total pairs to be parallelizable/independent so that they can execute parallel. Using Uniform Int Distribution.  \n
 * Total Pairs = numVNF*(numVNF-1) ordered (1,2)(2,1) except (1,1)(2,2)... \n
 * It then writes the pairs to the filename_vnf_parallelpairs {i_vnf -> {set of j_vnf it is parallel to}}
 * @param threshold percentage of pairs we want to parallelize. (1-99).
 * @param option way to obtain parallel pairs, 0 -> fixed, 1->approx around threshold, 2->exact to threshold (1-85%).
 */
bool Simulations::findRandomParallelPairs(const float& threshold, int option){
    if(threshold<=0 or threshold>=100) {
        cerr<<"\n Threshold (1-99) value for random parallel pairs is invalid: "<<threshold;
        return false;
    }
    const unsigned int numVNF = VNFNetwork.numVNF;
    totalPairs =  numVNF*(numVNF-1) - numVNF; ///< total unique pairs
    if(option == 0) {
        if(threshold <= 50){
            numParallelPairs = 60;
            totalPairs = 120;
            parallelPairs={
                    {1,{{12,pktNoCopy},{9,pktCopy},{7,pktCopy},{3,pktCopy},{2,pktNoCopy},{4,pktNoCopy}}},
                    {2,{{10,pktCopy},{4,pktNoCopy},{5,pktCopy},{8,pktNoCopy},{12,pktNoCopy},{6,pktNoCopy},{7,pktNoCopy}}},
                    {3,{{12,pktNoCopy}}},
                    {4,{{3,pktNoCopy},{2,pktCopy},{11,pktNoCopy}}},
                    {5,{{3,pktNoCopy},{12,pktCopy},{8,pktNoCopy},{1,pktNoCopy},{7,pktCopy},{6,pktNoCopy},{10,pktCopy}}},
                    {6,{{8,pktCopy},{9,pktCopy}}},
                    {7,{{12,pktCopy},{4,pktCopy},{10,pktNoCopy},{11,pktNoCopy},{6,pktNoCopy},{9,pktCopy},{3,pktNoCopy},{8,pktNoCopy},{5,pktNoCopy}}},
                    {8,{{5,pktNoCopy},{6,pktNoCopy}}},
                    {9,{{7,pktCopy},{2,pktCopy},{8,pktNoCopy},{11,pktNoCopy},{10,pktCopy},{4,pktNoCopy}}},
                    {10,{{1,pktNoCopy},{7,pktCopy},{12,pktNoCopy},{4,pktNoCopy}}},
                    {11,{{12,pktCopy},{9,pktCopy},{4,pktCopy},{2,pktNoCopy},{7,pktNoCopy},{1,pktCopy},{8,pktCopy},{6,pktCopy}}},
                    {12,{{6,pktCopy},{3,pktCopy},{9,pktNoCopy},{8,pktNoCopy},{7,pktCopy}}}
            };
        }//50
        else if (threshold <= 60){
                numParallelPairs = 72;
                totalPairs = 120;
                parallelPairs={
                        {1,{    {6,pktCopy},    {10,pktCopy},{12,pktNoCopy},{8,pktCopy},{4,pktCopy }}},
                        {2,{    {10,pktCopy},   {9,pktNoCopy},{11,pktCopy},{3,pktNoCopy},{12,pktCopy},{4,pktCopy},{1,pktCopy},{5,pktCopy},{8,pktNoCopy},{7,pktNoCopy},{6,pktCopy }}},
                        {3,{    {5,pktNoCopy},  {11,pktCopy},{12,pktNoCopy }}},
                        {4,{    {9,pktCopy},    {7,pktCopy},{6,pktNoCopy},{11,pktNoCopy},{12,pktNoCopy }}},
                        {5,{    {10,pktNoCopy}, {3,pktCopy},{11,pktNoCopy},{12,pktNoCopy},{7,pktNoCopy},{6,pktNoCopy},{1,pktNoCopy},{9,pktCopy},{2,pktCopy }}},
                        {6,{    {7,pktNoCopy},  {10,pktCopy},{1,pktCopy }}},
                        {7,{    {5,pktCopy},    {2,pktCopy},{11,pktNoCopy},{4,pktNoCopy},{8,pktCopy},{9,pktCopy},{6,pktNoCopy},{10,pktNoCopy},{3,pktNoCopy }}},
                        {8,{    {4,pktNoCopy},  {6,pktCopy},{12,pktCopy},{9,pktNoCopy},{11,pktCopy }}},
                        {9,{    {1,pktNoCopy},  {7,pktCopy},{10,pktNoCopy},{2,pktCopy }}},
                        {10,{   {3,pktNoCopy},  {4,pktNoCopy},{8,pktCopy},{9,pktNoCopy }}},
                        {11,{   {9,pktCopy},    {10,pktNoCopy},{1,pktCopy},{2,pktNoCopy},{6,pktCopy},{12,pktNoCopy},{7,pktCopy }}},
                        {12,{   {3,pktCopy},    {11,pktNoCopy},{9,pktCopy},{5,pktNoCopy},{4,pktCopy},{2,pktCopy},{6,pktNoCopy }}}
                };
        }



        if(debug)cout<<"\n\tRandom Parallel Pairs: "<<numParallelPairs<<" | "<<(100.0*numParallelPairs/totalPairs)<<" % out of total:"<<totalPairs<<" unique pairs.";
        return true;
    }//if(option == 0)

    parallelPairs.clear();
    numParallelPairs=0;
    const unsigned int numOfIndependentPairsNeeded = round(totalPairs*(threshold/100.0));

//    std::uniform_int_distribution<unsigned int> withCopyOrWithout_distribution(pktNoCopy, pktCopy);      ///< for selecting operation Parallel/Not Parallel

    if(option == 1 or threshold >=85){
        std::mt19937_64 rd_generator(chrono::system_clock::now().time_since_epoch().count()); ///< mt19937 is a standard mersenne_twister_engine
        std::uniform_int_distribution<unsigned int> rd_threshold_distribution(1, totalPairs);      ///< for selecting operation Parallel/Not Parallel
        for(unsigned int i_vnf=VNFNetwork.srcVNF; i_vnf<=numVNF; i_vnf++) {
            for (unsigned int j_vnf = VNFNetwork.srcVNF; j_vnf <= numVNF; j_vnf++) {
                if(i_vnf == j_vnf)continue;
                if( rd_threshold_distribution(rd_generator) <= numOfIndependentPairsNeeded){
                    parallelPairs[i_vnf][j_vnf] = boolean_distribution(rd_generator);
                    numParallelPairs++;
                }// threshold %pairs are generated
            }//j
        }// i
    }else if (option == 2){
        std::mt19937 rd_generator(chrono::system_clock::now().time_since_epoch().count()); ///< mt19937 is a standard mersenne_twister_engine
        std::uniform_int_distribution<unsigned int> fn_distribution(VNFNetwork.srcVNF, numVNF);
        for(unsigned int i=1; i<=numOfIndependentPairsNeeded; i++) {
            unsigned int fi = fn_distribution(rd_generator);
            unsigned int fj = fn_distribution(rd_generator);
            while( (fi == fj) or parallelPairs[fi].count(fj)>0 ){ // fj same as fi or fi already paired with fj
                fj=fn_distribution(rd_generator);
            }

            parallelPairs[fi][fj] = boolean_distribution(rd_generator);
            numParallelPairs++;
        }
    }//opt 2

    if(debug)cout<<"\n\tRandom Parallel Pairs: "<<(100.0*numParallelPairs/totalPairs)<<"% | "<<numParallelPairs<<" out of "<<totalPairs<<" unique pairs.";
    /// writing to file.
    ofstream fout;
    string filepathExt = output_directory+fullDirName+"ST_"+sim_name+"-Th"+to_string((100*numParallelPairs/totalPairs))+"_ParallelPairs.txt";
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fout<<"\n\tRandom Parallel Pairs: "<<(100.0*numParallelPairs/totalPairs)<<"% | "<<numParallelPairs<<" out of "<<totalPairs<<" unique pairs.";
    ///for reuse
    fout<<"\nparallelPairs={";
    for(int i_vnf=VNFNetwork.srcVNF; i_vnf<=numVNF; i_vnf++) {
        fout<<"\n  {"<<i_vnf<<",{"; int cnt = 1;
        for(const auto& j_vnf: parallelPairs[i_vnf]){

            fout<<"{"<<j_vnf.first<<",";
            if(j_vnf.second == pktNoCopy) fout<<"pktNoCopy}";
            else fout<<"pktCopy}";
            if(cnt != parallelPairs[i_vnf].size()) fout<<", ";
            cnt++;
        }
        fout<<"}}";
        if(i_vnf != numVNF) fout<<", ";
    }// i
    fout<<"\n};";
    fout.close();
    return true;
}

/*! Likelihood of parallelizing two functions fi and fj defined as
 * L_{fi,fj} = |num of sfc in which fi and fj are present and are parallelizable| divided by |numOfSFCS|
 * @param opt opt=1 mean calculated using SumLikelihood/(totalpairs - numvnfs) unique pairs, IF opt=2 then  SumLikelihood/(numParallelPairs)
 */
bool Simulations::calcLikelihoodOfTwoFunctions(const int& opt, const int& showInConsole){
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
    if(showInConsole == showDetailed){
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

/*!
 * @param fnInstScalingFactor (0.5-1] how much scaling (<=1 ==1 full scaling) in function instance we wanted as compared to actual instances needed (val=0).
 * @param alphaDecayRate (0-1) weightage of giving more importance to far away nodes
 * @param choice_pxs opt=1 choose bst candidate node for fcur, opt=2 choose bst candidate node for fcur-1(prev)
 * @param choice_dist_pref opt=1 distance prefernce based on hop count, opt=2 based on actual nearest distance
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 */
bool Simulations::DeploymentVNF_ScoreMethod(const float& fnInstScalingFactor, const float& alphaDecayRate, const int& choice_pxs, const int& choice_dist_pref, const int& showInConsole){//DeploymentVNF_ScoreMethod

    if((fnInstScalingFactor<0 or fnInstScalingFactor>1) or (alphaDecayRate<=0 or alphaDecayRate>=1) or (choice_pxs<=0 or choice_pxs >=3) or (choice_dist_pref<=0 or choice_dist_pref >=3))
        return false;

    total_service_capacity = 0;
    finalInstancesCount.clear();
    PN_2_VNF.clear();
    I_VNFINST_2_PN.clear();

    const unsigned int& numVNFs = VNFNetwork.numVNF;
    const unsigned int& numPNs =  PhysicalNetwork.numV;
    const unsigned int& numSFCs = allSFC.size();

    vector<unordered_map<unsigned int , unordered_map<unsigned int, double>>> dist_sfc_f_p(numSFCs); ///< score/priority of sfc s to put function fi on physical node pn.
    ///< Calculation of priority of sfc to place f on p based on its shortest path
    for(unsigned int aid=0; aid<numSFCs; aid++){ /*! for each sfc */
        const ServiceFunctionChain& sfc = allSFC[aid];
        const unsigned int& Len = sfc.numVNF;

        unordered_map<unsigned int, unsigned int> f_bstloc;
        unordered_map<unsigned int, unordered_map<unsigned int, double>>& priority_fi_pn = dist_sfc_f_p[aid];
        vector<unsigned int> access_nodes_shortest_paths = PhysicalNetwork.constructShortestPath(sfc.access_nodes.first, sfc.access_nodes.second);

        for(unsigned int vi=0; vi<Len; vi++){ ///finding best physical nodes for each function in sfc
            unsigned int path_id_star =   ceil((vi*access_nodes_shortest_paths.size()) / sfc.numVNF);
            f_bstloc[sfc.vnfSeq[vi]] = access_nodes_shortest_paths.at(path_id_star);
//            if(showInConsole == showDetailed){ cout<<"\nf"<<sfc.vnfSeq[vi]<<" id:"<<path_id_star<<" pn:"<< access_nodes_shortest_paths[aid][path_id_star]; }
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

//            if(showInConsole == showDetailed){ cout<<"\n[SFC:"<<aid<<"]"<<" f"<<fcur<<" (bst_p:"<<f_bstloc[fcur]<<") (px*:"<<pxs<<" py*:"<< pys<<") :: "; }
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
//                if(showInConsole == showDetailed){ cout<<"\n\t"<<pxs<<"->"<<pi<<"->"<<pys<<" ("<<val_pxs_pi<<" + "<<val_pi_pys<<") | sum: "<<priority_fi_pn[fcur][pi] ;}
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
            if(showInConsole == showDetailed){cout<<"\nf"<<f<<" on p"<<p<<"\t is order: "<< score_order_f_p<<" + cluster: "<<score_cluster_f_p<<"\t order(1+cluster)=" <<score_f_p[f][p];}
        }///for each p
    }///for each f

/////////////////////////////////////////////////////
    unordered_map<unsigned int, unsigned int> f_inst_needed; ///< approx minimum function instances needed according to arrival rate of sfc.
    unsigned long demand = 0;
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        vector<type_delay> items;
        for(const auto& sfcArrid: VNF_2_SFC[f]){  ///< push the arrival rate into items vector
            items.push_back(allSFC[sfcArrid].trafficArrivalRate);
        }
        f_inst_needed[f] = best_fit_decreasing<type_delay>(VNFNetwork.VNFNodes[f].serviceRate, items).size();
        demand += f_inst_needed[f]*VNFNetwork.VNFNodes[f].requirement.cores;
        if(showInConsole == showDetailed){
            cout<<"\n f:"<<f<<" is on #sfc("<<VNF_2_SFC[f].size()<<")";
            type_delay arrival_rate_sfc = 0;
            for(const auto& sfcArrid: VNF_2_SFC[f]){
//                cout<<"{"<<sfcid<<","<<allSFC[sfcArrid].trafficArrivalRate<<"}";
                arrival_rate_sfc +=  allSFC[sfcArrid].trafficArrivalRate;
            }
            cout<<" needed instances = "<<arrival_rate_sfc<<"/"<<VNFNetwork.VNFNodes[f].serviceRate<<" -> div:"<<arrival_rate_sfc/VNFNetwork.VNFNodes[f].serviceRate<<" -> BF:"<<f_inst_needed[f]<< " -> FF:";
//            vector<Bin<type_delay>> bestfit_bins = best_fit_decreasing<type_delay>(VNFNetwork.VNFNodes[f].serviceRate, items);
            vector<Bin<type_delay>> firstfit_bins = first_fit_decreasing<type_delay>(VNFNetwork.VNFNodes[f].serviceRate, items);
            cout<<firstfit_bins.size();
//            for (int i = 0; i < bestfit_bins.size(); ++i) {
//                cout << "\n\tbestfit_bin " << i + 1 << ": Capacity = " << bestfit_bins[i].get_capacity() << ", Space Left = " << bestfit_bins[i].get_space_left();
////            cout<<" \t["; for(const auto& x: bins[i].items) cout<<x<<"; "; cout<<"]";
//            }
        }//showInConsole
    }/// min instances needed

    vector<pair<unsigned int, unsigned int>> f_inst_final; ///< final {function -> instance count mapping}
    unsigned int used_cores = 0; ///< number of cores actually used
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        unsigned int instcnt=0;
        float f_inst_on_scaled=0;
        f_inst_on_scaled = (1.0*f_inst_needed[f] * PhysicalNetwork.sum_cores/demand);
        if(PhysicalNetwork.sum_cores >= demand){
            instcnt = (1.0-fnInstScalingFactor)*(f_inst_needed[f]) + fnInstScalingFactor*(f_inst_on_scaled);
        }else{
            instcnt = floor(f_inst_on_scaled);
        }
        /// phele kuch instance need tha but abhi zero ho gya and cores bhi hai use krne ko,
        if(f_inst_needed[f] > 0 and instcnt == 0 and (used_cores + VNFNetwork.VNFNodes.at(f).requirement.cores <= PhysicalNetwork.sum_cores)){
            instcnt = 1;
        }
        used_cores += instcnt*VNFNetwork.VNFNodes.at(f).requirement.cores;
        f_inst_final.emplace_back(f, instcnt);
        if(showInConsole == showDetailed){  cout<<"\n f:"<<f<<" "<<PhysicalNetwork.sum_cores<<"/"<<demand<<" scaled instances = ("<<f_inst_on_scaled<<" -> "<<instcnt<<")";}
    }

    ///< sorting function in descending order based on their original instance demand. scaling would be in similar manner.
    sort(f_inst_final.begin(), f_inst_final.end(), [&](const auto& left, const auto& right){
        const auto& first_instcnt = f_inst_needed[left.first];/// original instance cnt needed by left fn
        const auto& second_instcnt = f_inst_needed[right.first]; /// original instance cnt needed by right fn
        if(first_instcnt == second_instcnt){
            return VNFNetwork.VNFNodes[left.first].serviceRate > VNFNetwork.VNFNodes[right.first].serviceRate;
        }
        else return first_instcnt > second_instcnt;
    });

///////////////////////////////////////////////////////////////////////
    /*unordered_map<unsigned int,unsigned int> nodesToDeploy;
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
    }*/

//// DEPLOYMENT
    unsigned int totalInstDeployed = 0, totalCoresUsed = 0;
    unordered_map<unsigned int, Bin<unsigned int>> nodesToDeploy;
    for(unsigned int p = PhysicalNetwork.srcV; p<=numPNs; p++){
        nodesToDeploy[p] = Bin<unsigned int>(PhysicalNetwork.PNode[p].capacity.cores);
    }

    auto ForEachFnInst = [&](const unsigned int& f, const unsigned int& instfirstval, const unsigned int& instcnt){
        const unsigned int& reqcore = VNFNetwork.VNFNodes.at(f).requirement.cores; ///< this VNfunction required how many cores
        for(unsigned int i=instfirstval; i<=instcnt; i++){
            unsigned int pstar = 0; type_delay mxscore = 0;

            for(const auto& [p, bin]: nodesToDeploy){ /// among all the nodes
                if(bin.get_space_left() >= reqcore and score_f_p[f][p] > mxscore){ /// if that  node has the capacity then find the max score in that
                    pstar = p;
                    mxscore = score_f_p[f][p];
                }
            }

            if(pstar == 0) break;
            I_VNFINST_2_PN[f][i] = pstar;  ///< saving ans
            PN_2_VNF[pstar].push_back({f, i});
            nodesToDeploy[pstar].pack(reqcore);
            if(nodesToDeploy[pstar].get_space_left() == 0) nodesToDeploy.erase(pstar);

            for(auto& [p, _]: nodesToDeploy){
                double hopCnt;
                if(choice_dist_pref == 1){
                    hopCnt = PhysicalNetwork.constructShortestPath(p, pstar).size() - 1;
                }else if(choice_dist_pref == 2){
                    hopCnt = PhysicalNetwork.dist[p][pstar] ;// /(PhysicalNetwork.constructShortestPath(pxs, pi).size());
                }
                const double ratio = 1.0/hopCnt;
                /// update score for the same function type, based on alphaDecayRate, node far from pstar has a higher priority to install an instance of the same function f
                score_f_p[f][p] = score_f_p[f][p] * pow(alphaDecayRate, ratio); //exp(ratio);
                /// update score for the different independent function type fj, such that deploy them near to f
                for(const auto [fj,__]:parallelPairs[f]){
                    score_f_p[fj][p] = score_f_p[fj][p]*pow( 1 + likelihood[f][fj] , ratio);
                }
            }

            finalInstancesCount[f]++; /// update the final instance count of each VNF
            total_service_capacity += VNFNetwork.VNFNodes[f].serviceRate;
            totalInstDeployed++;
            totalCoresUsed+=reqcore;
        } /// for each inst
    };

    for(const auto&[f, instcnt]: f_inst_final){
        ForEachFnInst(f, 1, instcnt);
    }

    if((PhysicalNetwork.sum_cores < demand and used_cores < PhysicalNetwork.sum_cores) or fnInstScalingFactor==1){
        int idx=0;
        while((idx<f_inst_final.size()) and used_cores + VNFNetwork.VNFNodes.at(f_inst_final[idx].first).requirement.cores <= PhysicalNetwork.sum_cores){
            f_inst_final[idx].second += 1; /// increase instance cnt
            used_cores += VNFNetwork.VNFNodes.at(f_inst_final[idx].first).requirement.cores;
            ForEachFnInst(f_inst_final[idx].first, finalInstancesCount[f_inst_final[idx].first]+1, finalInstancesCount[f_inst_final[idx].first]+1);
            idx++;
        }
    }



    ofstream fout;
    string filepathExt = output_directory+fullDirName+"ST_"+sim_name+"-VNFDeployment-"+filename_sfc+".txt";

    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
    fout<<"With Options:";
    fout<<"\n   * opt scale: "<<fnInstScalingFactor <<" Function Instance Scaling Factor (0.5-1] (how much scaling (<=1 ==1 full scaling) in function instance we wanted as compared to actual instances needed (val=0))";
    fout<<"\n   * opt rate: "<<alphaDecayRate <<" Alpha Decay Rate((0-1) weightage of giving more importance to far away nodes)";
    fout<<"\n   * opt pxs: "<<choice_pxs<<" Choice of pxs(opt=1 choose bst candidate physical node for vnfcur, opt=2 choose bst candidate node for fcur-1(prev))";
    fout<<"\n   * opt dist: "<<choice_dist_pref<<" Choice of distance in sfc preferece(opt=1 distance prefernce based on hop count, opt=2 based on actual nearest distance";

    fout<<"\n\nVNFs ("<<VNFNetwork.filename_vnf<<") Deployement on the Network ("<<PhysicalNetwork.filename_network<<") Based on SFCs ("<<filename_sfc<<")";
    fout<<"\n\tRandom Parallel Pairs: "<<(100.0*numParallelPairs/totalPairs)<<"% | "<<numParallelPairs<<" out of "<<totalPairs<<" unique pairs.";
    fout<<"\n\tLikelihood mean: "<<likelihood_mean;

/// for readability --------------------------------
    fout<<"\n\nVNF to Node mapping: Total [NW Cores: "<<PhysicalNetwork.sum_cores<<"] [Instances Deployed: "<<totalInstDeployed<<" | Cores Used: "<<totalCoresUsed<<" | Demand: "<<demand<<"]";
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        fout<<"\n   VNF:"<<setw(2)<<f<<" [cnt: "<<f_inst_needed[f]<<"->"<<finalInstancesCount[f]<<"]  ";
        for(unsigned int i=1; i<=finalInstancesCount[f]; i++){
            fout<<"f"<<f<<char(i+96)<<" -> p:"<<I_VNFINST_2_PN[f][i]<<" | \t";
        }
    }

    fout<<"\n\nPhysical Nodes to VNFs mapping: G("<<PhysicalNetwork.numV<<","<<PhysicalNetwork.numE<<"): ";
    for(unsigned int p=PhysicalNetwork.srcV; p<=numPNs; p++){
        fout<<"\n   p:"<<setw(2)<<p;
        fout<<" [cores:"<<PhysicalNetwork.PNode[p].capacity.cores<<"]\t-> { ";
        for(const auto&[f, i]: PN_2_VNF[p]) fout<<setw(2)<<"f"<<f<<char(i+96)<<","; fout<<"}";
        fout<<"\t[deg:"<<PhysicalNetwork.PNode[p].neighbours.size()<<"][Cl:"<<PhysicalNetwork.PNode[p].clusteringcoeff_local<<"]";
    }  fout<<"\n\tClustering Coeff mean:"<<PhysicalNetwork.clusteringcoeff_mean<<" | global:"<<PhysicalNetwork.clusteringcoeff_global;

/// for resuse --------------------------------------------
    fout<<"\n--------------------------------------------";
    fout<<"\n\nfinalInstancesCount={";
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        fout<<"{"<<setw(2)<<f<<","<<setw(2)<<finalInstancesCount[f]<<"}";
        if(f != numVNFs) fout<<", ";
    }
    fout<<"};";

    fout<<"\nI_VNFINST_2_PN={";
    for(unsigned int f=VNFNetwork.srcVNF; f<=numVNFs; f++){
        fout<<"\n  {"<<setw(2)<<f<<", {";
        for(unsigned int i=1; i<=finalInstancesCount[f]; i++){
            fout<<"{"<<setw(2)<<i<<","<<setw(2)<<I_VNFINST_2_PN[f][i]<<"}";
            if(i != finalInstancesCount[f]) fout<<", ";
        } fout<<"}}";
        if(f != numVNFs) fout<<",";
    }
    fout<<"\n};";

    fout<<"\nPN_2_VNF={";
    for(unsigned int p=PhysicalNetwork.srcV; p<=numPNs; p++){
        fout<<"\n  {"<<setw(2)<<p<<", {"; int cnt=1;
        for(const auto&[f, i]: PN_2_VNF[p]){
            fout<<"{"<<setw(2)<<f<<","<<setw(2)<<i<<"}";
            if(cnt != PN_2_VNF[p].size()) fout<<", "; cnt++;
        }
        fout<<"}}";
        if(p != numPNs) fout<<",";
    }
    fout<<"\n};";

    fout.close();

//    if(showInConsole == showFinal){
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


void Simulations::printGenericNetworkGraph() const{
    if (PhysicalNetwork.mat.empty()) {
        cout << "Network Graph is Empty."<<endl;
        return;
    }
    ofstream fout;
    const string& fileName = "Network_G_"+to_string(PhysicalNetwork.numV)+"_"+to_string(PhysicalNetwork.numE); ///< (Network Graph + V + E) name without space
    const string& filepath = output_directory+fullDirName+diagram_directory+fileName; ///<path to .gv file without extention
    string filepathExt = filepath+".gv"; ///<filepath with extention of .gv
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

/**********************/
    const string& nodeColor = "darkslateblue", nodeFontColor = "firebrick4";
    const string& edgeColor = "darkslateblue", edgeFontColor = "firebrick4";

    fout << "graph "<<fileName<<" {" << endl;
    fout << "  label = \"Network Graph G("<<PhysicalNetwork.numV<<","<<PhysicalNetwork.numE<<")\" \n";
    fout << "  ranksep=\"equally\"; \n"
            "  rankdir=LR; " << endl << endl;
    fout << "node [fixedsize=shape margin=0.3 width=0.6 shape=circle color="<<nodeColor<<" fontcolor="<<nodeFontColor<<"] "<<endl;
    fout << "edge [color="<<edgeColor<<" fontcolor="<<edgeFontColor<<" fontsize=12]"<<endl;

    // Printing Edges and Nodes using Adjacency Matrix
    for(unsigned int u_val = PhysicalNetwork.srcV; u_val <= PhysicalNetwork.numV; u_val++)
    {
        fout << u_val << " [label = <&eta;<sub>"<<u_val<<"</sub>>];"<< endl;
        for (unsigned int v_val = u_val+1; v_val<=PhysicalNetwork.numV ; ++v_val) {
            if(PhysicalNetwork.mat[u_val][v_val] != (type_delay)0){
                fout <<"\t"<< u_val << " -- " << v_val <<" [label ="<<PhysicalNetwork.mat[u_val][v_val]<<"]"<< endl;
            }
        }
        fout << endl ;
    }
    fout << "}" << endl;
/**********************/
    fout.close();
    const string cmd = "dot -Tpng "+filepathExt+" -o "+filepath+".png";
    system((const char*)cmd.c_str());
    if(debug)cout<<"\nGraphviz: Original Network Graph G("<<PhysicalNetwork.numV<<","<<PhysicalNetwork.numE<<") printed. File:"<<filepath<<".png \n";
    if(!debug)remove(filepathExt.c_str());
}

#endif //SFC_PARALLELIZATION_SIMULATIONCLASS_H
