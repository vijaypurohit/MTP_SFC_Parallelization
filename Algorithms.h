//
// Created by vijay on 01-04-2023.
//


#ifndef SFC_PARALLELIZATION_ALGORITHMS_H
#define SFC_PARALLELIZATION_ALGORITHMS_H


/*! Dynamic Program to search all the available COMPOSITION for a SFC of length K.\n
 * A COMPOSITION of an integer n is a tuple (ordered list) of positive integers whose elements sum to n.
 * (sometimes also called integer composition, ordered partition or ordered integer partition). This is an additive representation of n.  \n
 * A part in a composition is sometimes also called a summand.
 * A composition of n into k (positive) parts is an ordered k-tuple (x1, . . . , xk) with each xi ∈ N and x1 + · · · + xk = n.
 * If the order is not taken into account then the sum is a partition. \n
 * c(n) = 2^(n-1) number of composition of size n. \n
 *  <a href="https://oeis.org/wiki/Integer_compositions">wiki/Integer_compositions</a>
 * @param k small k is starting point, from which you have to calculate, as some will be precalculated and saved.
 * @param K chain length which is maxSFClen, that is upto what you want to calculate all possible enumeration
 * @param showInConsole to show the output in console. Time will increase to 2000-2500 ms
 * @return updated clusterSz[i] for k <= i <=K.\n
 * @example
 *  K[1] total(1), {{1}} \n
 *  K[2] total(2), {{1,1},{2}} \n
 *  K[3] total(4), {{1,1,1},{2,1},{1,2},{3}} \n
 *  K[4] total(8), {{1,1,1,1},{2,1,1},{1,2,1},{3,1},{1,1,2},{2,2},{1,3},{4}}
 */
void integerCompositionsEnumeration(unsigned int k, unsigned int K, bool showInConsole=false) {
    if(k>K) return;
    if(integerCompositions.find(k - 1) == integerCompositions.end()){
        string errorMsg = "Previous Cluster size k-1="+to_string(k-1)+ " does not exist. Caclculate that first. Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
//    vector<vector<vector<int>>> clusterSz(K+1);
    for(unsigned int ki = k; ki<=K; ki++){
        vector<vector<unsigned int>> ans;
        for(unsigned int i=1; i<=ki; i++){
            for(auto s_dash: integerCompositions[ki - i]){
                s_dash.push_back(i);
                ans.push_back(s_dash);
            }
        }
        integerCompositions[ki] = std::move(ans);
    }

    if(showInConsole == showFinal){
        cout<<"\nintegerCompositions={";
        for(unsigned int ki = k; ki<=K; ki++){
            cout << "\n\tK[" << ki << "] total(" << integerCompositions[ki].size() << "), {";
            for(const auto& x: integerCompositions[ki]) {
                cout<<"{"; for(unsigned int i=0; i<x.size()-1; i++) cout<<x[i]<<","; cout<<x[x.size()-1]<<"},";
            }  cout<<"}";
        } cout<<"\n}end;";
    }
}


/*! Calculate all possible combination vectors (n Choose k) for given n and k.
 * @param n starting point from which you have to calculate new value, as some will be precalculated and saved.
 * @param N maximum value upto which you want to calculate all possible combination
 * @param showInConsole  to show the output in console. Time will increase to 2000-2500 ms
 * @return updated nCk[i] for n <= i <=N.\n
 *      n=1, {k=1, [1] }\n
 *      n=2, {k=1, [{1},{2}]  |         k=2, [{1,2}]   }\n
 *      n=3, {k=1, [{1},{2},{3} |       k=2, [{1,2},{1,3},{2,3}] | {k=3, [{1,2,3}] }\n
 *      n=4, {k=1, [{1},{2},{3},{4}] |  k=2, [{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}] | {k=3, [{1,2,3},{1,2,4},{1,3,4}, {2,3,4}] | k=4, [{1,2,3,4}] }
 */
void find_all_nCk(unsigned int n, unsigned int N, bool showInConsole=false){
    if(n>N) return;
    // Lambda Function to calculate nCk using backtracking application.
    std::function<void(unsigned int,unsigned int,vector<unsigned int>&,unsigned int&,unsigned int&)> combineHelper = [&combineHelper]
    (unsigned int st, unsigned int k, vector<unsigned int>& cur, unsigned int &x, unsigned int &y)->void
    {
        if(k==0) {
            nCk[x][y].push_back(cur);
            return;
        }
        for(unsigned int value=st; value <= x-k+1; value++) {
            cur.push_back(value);
            combineHelper(value+1, k-1, cur, x, y);
            cur.pop_back();
        }
    };

    for(unsigned int ni=n; ni<=N; ni++){
        for(unsigned int kr=1; kr<=ni; kr++){
            vector<unsigned int> cur;
            combineHelper(1, kr, cur, ni, kr);
        }
    }
    if(showInConsole == showFinal){
        cout<<"\nnCk={";
        for(unsigned int ni=n; ni<=N; ni++){
            cout<<"\n\tN["<<ni<<"]"<<", {";
            for(unsigned int kr=1; kr<=ni; kr++){
                cout<<"\n\t\tk["<<kr<<"], { ";
                for(auto vec: nCk[ni][kr] ) {
                    cout<<"{"; for(unsigned int i=0; i<vec.size()-1; i++) cout<<vec[i]<<","; cout<<vec[vec.size()-1]<<"},";
                }
                cout<<"}";
            }
            cout<<"\n\t}";
        }cout<<"\n}end;";
    }
}

/*!
 * @param numPN number of Physical Nodes in the network
 * @param numVNF  number of Virtual Network Function in the network
 * @param numOfSFCsNeeded number of sfc to generate
 * @param arrivalRateRange arrival rate range of each sfc
 * @param lenOfEachSFC len of each sfc (random or fixed length)
 * @param fullDirName directory name with slash on which to generate file
 * @param filename_sfc sfc file to be saved
 * @return status
 */
bool GenerateRandomSFCs(const unsigned int& numPN, const unsigned int& numVNF, const unsigned int& numOfSFCsNeeded, const pair<type_delay, type_delay>& arrivalRateRange, const pair<bool, unsigned int>& lenOfEachSFC, const string& fullDirName, const string& filename_sfc){//GenerateRandomSFCs

    const unsigned int srcVNF = 1, srcPN = 1;
    unsigned int minLenSFC = std::floor(0.3*numVNF); ///< if vnf = 10 or 12 then min length 3
    unsigned int maxLenSFC = std::floor(0.9*numVNF); ///< if max len of sfc not given then take according to number of VNFs
    if(maxLenSFC > maxSFCLength) maxLenSFC = maxSFCLength;

    std::mt19937_64 rd_generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<type_delay> rd_arrivalrate_distribution(arrivalRateRange.first, arrivalRateRange.second);

    /*const pair<bool, unsigned int> givenMinLenSFC; const pair<bool, unsigned int> givenMaxLenSFC;
    if(givenMinLenSFC.first and givenMinLenSFC.second > numVNF){ ///< if max len of sfc given then it should be less than number of VNF{
        cerr<<"\nGiven Minimum Length of SFC ("<<givenMinLenSFC.second<<") is more than number of VNFs ("<<numVNF<<")";
        return false;
    }else{
        minLenSFC = givenMinLenSFC.second;
    }
    if(givenMaxLenSFC.first and givenMaxLenSFC.second > numVNF){ ///< if max len of sfc given then it should be less than number of VNF{
        cerr<<"\nGiven Maximum Length of SFC ("<<givenMaxLenSFC.second<<") is more than number of VNFs ("<<numVNF<<")";
        return false;
    }else{
        maxLenSFC = givenMaxLenSFC.second;
    }*/

    std::vector<unsigned int> vnfseq(numVNF); ///All the VNFS index present in array index from 0 to numVNF-1 (size = numVNF)
    std::iota(begin(vnfseq), end(vnfseq), srcVNF); ///< fill with sequence from 1 to numVNF

    std::mt19937 rd_length_generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<unsigned int> rd_length_distribution(minLenSFC, maxLenSFC); ///< random distribution for sfc length

//    if(showInConsole == showFinal){
//        cout<<"\nPN Range: ["<<srcPN<<" - "<<numPN<<"]";
//        cout<<"\nVNF Range: ["<<srcVNF<<" - "<<numVNF<<"]";
//        cout<<"\nSFC Range: ["<<minLenSFC<<" - "<<maxLenSFC<<"] : cnt:"<<numOfSFCsNeeded;
//    }

    std::knuth_b rd_accessnode_generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<unsigned int> rd_access_nodes_distribution(srcPN, numPN); ///< sfc access nodes distribution according to physical network

    ofstream fout;
    string filepathExt = input_directory+fullDirName+filename_sfc;///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    fout<<numOfSFCsNeeded<<"\n\n";

    for(int idx=0; idx<numOfSFCsNeeded; idx++){//cnt
        std::shuffle(begin(vnfseq), end(vnfseq), rd_generator); ///< random shuffle the vnfs in array

        unsigned int psrc = rd_access_nodes_distribution(rd_accessnode_generator); ///< sfc src node
        unsigned int pdst = rd_access_nodes_distribution(rd_accessnode_generator); ///< sfc dst node
        while(psrc == pdst){  pdst = rd_access_nodes_distribution(rd_accessnode_generator); } ///< if dst node same as src node

        fout<<idx<<'\n';
        fout<<rd_arrivalrate_distribution(rd_generator)<<'\t'<<psrc<<' '<<pdst<<'\n';
        if(lenOfEachSFC.first){ ///< if length of each sfc is constant, take first len vnfs
            fout<<lenOfEachSFC.second<<'\n';
            for (int x = 0; x < lenOfEachSFC.second; ++x) {
                fout << vnfseq[x] << " ";
            }
        }else{
            unsigned int lenOfSFC = rd_length_distribution(rd_length_generator); ///< generate random length
            fout<<lenOfSFC<<'\n';
            for (int i = 0; i < lenOfSFC; ++i){
                fout << vnfseq[i] << ' ';
            }
        }
        fout<<"\n\n";
    }//cnt<=numOfSFCsNeeded

    return true;
}//GenerateRandomSFCs

/*!
 * @param numOfVNFsNeeded number of vnfs to be generated
 * @param serviceRateRange range of service rate for each vnf
 * @param funExeTimeRange range of function time for each vnf
 * @param fullDirName full directory name
 * @param filename_vnfs filename of vnf to be saved
 * @return status
 */
bool GenerateRandomVNFs(unsigned int& numOfVNFsNeeded, const pair<type_delay, type_delay>& serviceRateRange, const pair<type_delay, type_delay>& funExeTimeRange, const string& fullDirName, const string& filename_vnfs){//GenerateRandomVNFs

    if(numOfVNFsNeeded > maxNumVNFs) numOfVNFsNeeded = maxNumVNFs;

    const unsigned int srcVNF = 1;

    std::mt19937_64 rd_generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<type_delay> rd_servicerate_distribution(serviceRateRange.first, serviceRateRange.second);
    std::uniform_real_distribution<type_delay> rd_fnexetime_distribution(funExeTimeRange.first, funExeTimeRange.second);


    ofstream fout;
    string filepathExt = input_directory+fullDirName+filename_vnfs;///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::out);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    fout<<numOfVNFsNeeded;


    for(int idx=1; idx<=numOfVNFsNeeded; idx++){//cnt
        fout<<"\nVNF"<<idx<<'\n';
        fout<<idx<<'\t'<<1<<'\t';
        fout<<rd_servicerate_distribution(rd_generator)<<'\t'<<rd_fnexetime_distribution(rd_generator);
    }//cnt<=numOfSFCsNeeded

    return true;
}//GenerateRandomVNFs


/* **************************************************************************************************************** */ 
/*! Function to find all instances combination of parVNF in that stage. . find stage to instnace combination of given block of sfc.
     * @param csfi function index of current stage.
     * @param curInstComb current combination in iteration
     * @param stgid stgId/blockId for which we are finding combination of functions in that stg/block and to store in stg2InstCombinations.
     * @param curStg using to iterate all functions in the stage.
     * @param[out] stg2InstCombinations It stores all the stage wise instances combination of all stage in partParSFC. {stgid . 2d{ 1d instances combinations{pair<fun, inst>}  }}
     * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
     * @param[in] SimTest Simulation Object which contains SFC on which to perform the test.
     * @param[in]  cSFC  given SFC object for which we have to find minimum delay mapping
     * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
     * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
     * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] } \n
     * Time to calculte stg2InstCombinations . if in any block number of parallel functions are 10 and each have 5 max instances\n
        inst = 2 (exe time: 1-2ms) (possibilites: 1024 (2^10)) \n
        inst = 3 (exe time: 38-40ms) (possibilites: 59 049 (3^10))\n
        inst = 4 (exe time: 580-600ms) (possibilites: 10 48 576 )\n
        inst = 5 (exe time: 5700-5800ms) (possibilites: 97 65 625)\
     */
void construct_stageInstanceCombination(unsigned int csfi, vector<pair<unsigned int,unsigned int>>& curInstComb,  const unsigned int& stgid, const vector<unsigned int>& curStg,
                                        unordered_map<unsigned int,  vector<vector<pair<unsigned int,unsigned int>>>> &stg2InstCombinations,
                                        const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,
                                        const Simulations& SimTest, const ServiceFunctionChain& cSFC ){ //construct_stageInstanceCombination
    
    if(csfi == curStg.size()){ // all functions in stage iterated. curStg.size()==numOfFunction in that stage.
        stg2InstCombinations[stgid].push_back(curInstComb); // push the one answer into combination stg.
        return;
    }
    const unsigned int fnType = curStg[csfi];
    unsigned int totInstancs = SimTest.finalInstancesCount.at(fnType);
    for(unsigned int fnInstId=1; fnInstId<=totInstancs; fnInstId++){
        if (ifRejectSFC and oldUtilization.count(fnType) and oldUtilization.at(fnType).count(fnInstId) and ///< old utilization till now of VNF
            (oldUtilization.at(fnType).at(fnInstId) + cSFC.trafficArrivalRate >  SimTest.VNFNetwork.VNFNodes.at(fnType).serviceRate))  {
            continue; // don't take this instance if its utilization become more than service rate of function.
        }
        curInstComb.emplace_back(fnType, fnInstId); // push current instance
        construct_stageInstanceCombination(csfi+1, curInstComb, stgid, curStg, stg2InstCombinations, oldUtilization, SimTest, cSFC ); // call function for next instance
        curInstComb.pop_back(); // pop curInstComb instance and push next instance of same function.
    }
}//construct_stageInstanceCombination

/*! For a given stg2InstCombinations, it enumerate all the possible mappings we can give in each stage and calculate delay on the go.
     * @param stgid stgId/blockId for which we are enumerating instances.
     * @param prv_combination_utilization max utilization of the current path
     * @param curMapping cur function.instance mapping we iterating out of all possibilites.
     * @param bstMapping to save best mapping overall among all partial parallel chain/and its all instances.
     * @param minBstDelay min delay among all partial parallel chain/and its all instances.
     * @param minBstUtil min utilization of the path chosen with minBstDelay and among all paths with minBstDelay choose minUtil
     * @param partParSFC given partial SFC
     * @param stg2InstCombinations It contains all the stage wise instances combination of all stage in partParSFC. {stgid . 2d{ 1d instances combinations{pair<fun, inst>}  }}
     * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
     * @param[in] SimTest Simulation Object which contains SFC on which to perform the test.
     * @param[in]  cSFC  given SFC object for which we have to find minimum delay mapping
     * @param showInConsoleDetailed
     * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
     * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
     * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] }
     * allMappings are (total 36 = 3*6*2) \n
        0[1a 6a 4a 5a ]         1[1a 6a 4a 5b ]         2[1a 6a 4b 5a ]     3[1a 6a 4b 5b ]         4[1a 6a 4c 5a ]         5[1a 6a 4c 5b ]
        6[1a 6b 4a 5a ]         7[1a 6b 4a 5b ]         8[1a 6b 4b 5a ]     9[1a 6b 4b 5b ]         10[1a 6b 4c 5a ]        11[1a 6b 4c 5b ]
        12[1b 6a 4a 5a ]        13[1b 6a 4a 5b ]        14[1b 6a 4b 5a ]        15[1b 6a 4b 5b ]        16[1b 6a 4c 5a ]        17[1b 6a 4c 5b ]
        18[1b 6b 4a 5a ]        19[1b 6b 4a 5b ]        20[1b 6b 4b 5a ]        21[1b 6b 4b 5b ]        22[1b 6b 4c 5a ]        23[1b 6b 4c 5b ]
        24[1c 6a 4a 5a ]        25[1c 6a 4a 5b ]        26[1c 6a 4b 5a ]        27[1c 6a 4b 5b ]        28[1c 6a 4c 5a ]        29[1c 6a 4c 5b ]
        30[1c 6b 4a 5a ]        31[1c 6b 4a 5b ]        32[1c 6b 4b 5a ]
        3[1c 6b 4b 5b ]        34[1c 6b 4c 5a ]        35[1c 6b 4c 5b ] \n
     * 1 stage . 10 parallel func each with 5 max instances . 5^10 possibilities or 97,65,625 instances.
     */
void enumerate_InstanceCombination_CalcDelay(unsigned int stgid,  type_delay prv_combination_utilization, unordered_map<unsigned int,unsigned int>& curMapping, unordered_map<unsigned int,unsigned int>& curBstMapping, type_delay& minBstDelay, type_delay& minBstUtil,
                                             const vector<vector<unsigned int>>& partParSFC,  const unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> >& stg2InstCombinations,
                                             const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,
                                             const Simulations& SimTest, const ServiceFunctionChain& cSFC, const int& showInConsole = dontShow){//enumerate_InstanceCombination_CalcDelay
    if(stgid == partParSFC.size()) { // found one mapping then find corresponding delay 
        type_delay parSfcCost = calcD_ParallelSFC(cSFC, partParSFC, curMapping, oldUtilization, SimTest);

        if((parSfcCost < minBstDelay) or (approximatelyEqual<type_delay>(parSfcCost,minBstDelay,0.0005) and prv_combination_utilization < minBstUtil)){ // current mapping ka delay is less than min delay among partial sfc all instances.
            minBstDelay =  parSfcCost;
            curBstMapping = curMapping;
            minBstUtil = prv_combination_utilization;
        }
        else return;
        if(showInConsole >= showDetailed){ /// showing instance and its delay
            cout<<"\n\t"<<"["; for(const auto &blk: partParSFC){ for(const auto& fnid: blk){ cout<<fnid<<char(96+curMapping.at(fnid))<<" ";  }  } cout<<"]";
            cout<<"["<<parSfcCost<<" ] ("<<prv_combination_utilization<<" %)";
        }
        return;
    }
    for(const vector<pair<unsigned int,unsigned int>>& instComb: stg2InstCombinations.at(stgid)){
        type_delay utilization_sum_lgy=0, servicerate_sum_lgy=0;
        for(const auto& [fnType, fnInstId]: instComb) {
            curMapping[fnType] = fnInstId;
            if (oldUtilization.count(fnType) and oldUtilization.at(fnType).count(fnInstId)){ ///< old utilization till now of VNF
                utilization_sum_lgy += oldUtilization.at(fnType).at(fnInstId);
            } servicerate_sum_lgy += SimTest.VNFNetwork.VNFNodes.at(fnType).serviceRate;
        }
        enumerate_InstanceCombination_CalcDelay(stgid+1, max((utilization_sum_lgy/servicerate_sum_lgy)*100, prv_combination_utilization), curMapping, curBstMapping, minBstDelay, minBstUtil, partParSFC, stg2InstCombinations, oldUtilization, SimTest, cSFC, showInConsole);
    }
}//enumerate_InstanceCombination_CalcDelay

/*! FINDING THE SFC's VNF Instance Assignment IF PARALLELISM is DISABLED IN SFC
 * @param SimTest Simulation Object which contains SFC on which to perform the test.
 * @param cSFC  given SFC object for which we have to find minimum delay mapping
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 * @return status of the algorithm.
 */
int bruteForce_Sequential_Deployment(Simulations& SimTest, const ServiceFunctionChain& cSFC, const int& showInConsole = dontShow) {//bruteForce_Sequential_Deployment
    SFC_RESULT& singleSfcRes = SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index];
    const vector<vector<unsigned int>>& seqSFC = (*cSFC.partialParallelChains).front();
    const unsigned int& szStages = seqSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid . 2d{ 1d instances combinations{pair<fun, inst>}  }} */
    unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
    for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
        const auto& curStg = seqSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_stageInstanceCombination(0, curInstComb, stgid, curStg, stg2InstCombinations, SimTest.TestsResult[name_bruteForce].seq_utilization, SimTest, cSFC );
        if(stg2InstCombinations.find(stgid) == stg2InstCombinations.end()){
            return singleSfcRes.seq_status = noResDueToNoStg; ///< if there is no instance combinations in current stage then we can't proceed
        }
    }//stgid<szStages finding instances possibilities of each stage.
         
    if(showInConsole >= showDetailed){ // showing partParSFC info
        cout<<"\n Sequential-SFC:"; for(const auto& blks: seqSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; }  
        for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
            for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
            }  cout<<" }";
        }
    }// show stages wise instances combination

    /*! finding all the mapping possibilites for the current partParSFC instance combination at each stage.*/
    unordered_map<unsigned int,unsigned int> curMapping; ///< iterating mapping variable
    enumerate_InstanceCombination_CalcDelay(0, 0, curMapping, singleSfcRes.seq_fninstmap, singleSfcRes.seq_delay, singleSfcRes.seq_load, seqSFC , stg2InstCombinations, SimTest.TestsResult[name_bruteForce].seq_utilization, SimTest, cSFC, showInConsole);//, showInConsoleDetailed

    for(const auto& [fn, fninst]: singleSfcRes.seq_fninstmap){
        SimTest.TestsResult[name_bruteForce].seq_utilization[fn][fninst] += cSFC.trafficArrivalRate;
//        cout<<"{"<<fn<<","<<fninst<<"},";
    }
    if(showInConsole >= showFinal){
        cout<<"\n BruteForce-Instance-Assignment:: Sequential: SFCid:"<<cSFC.index<<" | delay:["<<singleSfcRes.seq_delay<<"] ("<<singleSfcRes.seq_load<<"%) :(";
        for(const auto& fnid : cSFC.vnfSeq){
            cout<<"f"<<fnid<<char(96+singleSfcRes.seq_fninstmap.at(fnid))<<"; ";
        } cout<<")";
    }
    return singleSfcRes.seq_status = algosuccess;
}//bruteForce_Sequential_Deployment

/*! FINDING THE SFC's VNF Instance Assignment IF Full PARALLELISM Enabled IN SFC
 * @param SimTest Simulation Object which contains SFC on which to perform the test. 
 * @param cSFC  given SFC object for which we have to find minimum delay mapping
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output) 
 * @return status of the algorithm.
 */
int bruteForce_FullParallel_Deployment(Simulations& SimTest, const ServiceFunctionChain& cSFC, const int& showInConsole = dontShow) {//bruteForce_FullParallel_Deployment
    SFC_RESULT& singleSfcRes = SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index];
    const vector<vector<unsigned int>>& fullParSFC = cSFC.vnfBlocksPar;
    const unsigned int& szStages = fullParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    unsigned long long totalInstComb = 1;
    /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid . 2d{ 1d instances combinations{pair<fun, inst>}  }} */
    unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
    for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
        const auto& curStg = fullParSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_stageInstanceCombination(0, curInstComb, stgid, curStg, stg2InstCombinations, SimTest.TestsResult[name_bruteForce].fullpar_utilization, SimTest, cSFC );
        if(stg2InstCombinations.find(stgid) == stg2InstCombinations.end()){
            return singleSfcRes.fullpar_status = noResDueToNoStg; ///< if there is no instance combinations in current stage then we can't proceed
        }
        totalInstComb *= stg2InstCombinations[stgid].size();
    }//stgid<szStages finding instances possibilities of each stage.

    cout<<"\n  BF-fully-parallel["<<SimTest.sfccompleted<<"/"<<SimTest.sortedSFCs.size()<<"]:"<<totalInstComb;

    if(showInConsole >= showDetailed){ // showing partParSFC info
        cout<<"\n Fully-Parallel-SFC:"; for(const auto& blks: fullParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; }  
        for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
            for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
            }  cout<<" }";
        }
    }// show stages wise instances combination

    /*! finding all the mapping possibilites for the current partParSFC instance combination at each stage.*/
    unordered_map<unsigned int,unsigned int> curMapping; ///< iterating mapping variable
    enumerate_InstanceCombination_CalcDelay(0,0, curMapping, singleSfcRes.fullpar_fninstmap, singleSfcRes.fullpar_delay, singleSfcRes.fullpar_load, fullParSFC , stg2InstCombinations, SimTest.TestsResult[name_bruteForce].fullpar_utilization, SimTest, cSFC, showInConsole);//, showInConsoleDetailed

    for(const auto& [fn, fninst]: singleSfcRes.fullpar_fninstmap){
        SimTest.TestsResult[name_bruteForce].fullpar_utilization[fn][fninst] += cSFC.trafficArrivalRate;
    }
    if(showInConsole == showFinal){
        cout<<"\n BruteForce-Instance-Assignment:: Fully-Parallel: SFCid:"<<cSFC.index<<" |  delay:["<<singleSfcRes.fullpar_delay<<"] ("<<singleSfcRes.fullpar_load<<"%) :";
        for(const auto& blk: cSFC.vnfBlocksPar){
            cout<<" ["; for(int fnid: blk){
                cout << "f" << fnid << char(96 + singleSfcRes.fullpar_fninstmap.at(fnid)) << " ";
            } cout<<"]";
        }
    }
    return singleSfcRes.fullpar_status=algosuccess;
}//bruteForce_FullParallel_Deployment

/*! FINDING THE Optimal SFC's VNF Instance Assignment in case IF PARALLELISM Enabled IN SFC.
 * Intance deployment is based on min delay of the path.
 * For a partial parallel chain, find all its instance combinations in each stage, then out of all instance mappings find min delay.
 * @param[in] SimTest Simulation Object which contains SFC on which to perform the test. 
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 */
bool bruteForce_PartParallel_Deployment(Simulations& SimTest, const ServiceFunctionChain& cSFC, const int& showInConsole = dontShow) {//bruteForce_PartParallel_Deployment
    const vector<vector<vector<unsigned int>>>& allPartParSFC = *cSFC.partialParallelChains;
    SFC_RESULT& singleSfcRes = SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index];

    type_delay bst_partpar_delay = std::numeric_limits<type_delay>::max();
    //    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}} ;
    for(int ppsidx=int( allPartParSFC.size())-1; ppsidx>=0; --ppsidx){ /*! For each partial parallel SFC without src and dest block/stage.*/
        /*! {{1},{6,4},{5}}; Each Partial SFC is without src and dest block/stage. */
//        cout<<"\r  ppar["<<SimTest.sfccompleted<<"/"<<SimTest.sortedSFCs.size()<<"]("<<allPartParSFC.size()-ppsidx<<"/"<<allPartParSFC.size()<<")";
//        cout<<"\r\t\t\t("<<allPartParSFC.size()-ppsidx<<"/"<<allPartParSFC.size()<<")"; //showing progress

        const vector<vector<unsigned int>>& partParSFC=allPartParSFC.at(ppsidx); ///< for each of the partial parallel SFC of the givenParVNF Blocks
        const unsigned int szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

        unsigned long long totalInstComb = 1;
        /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid . 2d{ 1d instances combinations{pair<fun, inst>}  }} */
        unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
        for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
            const auto& curStg = partParSFC[stgid];
            vector<pair<unsigned int,unsigned int>> curInstComb;
            construct_stageInstanceCombination(0, curInstComb, stgid, curStg, stg2InstCombinations, SimTest.TestsResult[name_bruteForce].ppar_utilization, SimTest, cSFC );
            if(stg2InstCombinations.find(stgid) == stg2InstCombinations.end()){
                return singleSfcRes.ppar_pid = noResDueToNoStg; ///< if there is no instance combinations in current stage then we can't proceed
            }
            totalInstComb *= stg2InstCombinations[stgid].size();
        }//stgid<szStages finding instances possibilities of each stage.
        if(stg2InstCombinations.size() != szStages)continue; ///< if number of stages here are not same as sfc then no need to proceed.

        cout<<"\r  BF-partial-parallel["<<SimTest.sfccompleted<<"/"<<SimTest.sortedSFCs.size()<<"]("<<allPartParSFC.size()-ppsidx<<"/"<<allPartParSFC.size()<<"):"<<totalInstComb<<"       ";

        if(showInConsole >= showDetailed){ // showing partParSFC info
            cout<<"\n Partial-Parallel-SFC["<<ppsidx<<"]: "; for(const auto& blks: partParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<") ---------- - --------- - ------";
            for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
                cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
                for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                    cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
                }  cout<<" }";
            }
        }// show stages wise instances combination

        /*! finding all the mapping possibilites for the current partParSFC instance combination at each stage.*/
        unordered_map<unsigned int,unsigned int> curMapping; ///< iterating mapping variable

        enumerate_InstanceCombination_CalcDelay(0, 0, curMapping, singleSfcRes.ppar_fninstmap, bst_partpar_delay, singleSfcRes.ppar_load, partParSFC, stg2InstCombinations, SimTest.TestsResult[name_bruteForce].ppar_utilization, SimTest, cSFC, showInConsole);
        if(bst_partpar_delay < singleSfcRes.ppar_delay ){
            singleSfcRes.ppar_pid = ppsidx;
            singleSfcRes.ppar_delay = bst_partpar_delay;
        } 

    }// for each allPartParSFC ppsidx.
 
//    if(singleSfcRes.ppar_pid == noResPar){
//        return algostopped;
//    }
    
    // update
    for(const auto& [fn, fninst]: singleSfcRes.ppar_fninstmap){
        SimTest.TestsResult[name_bruteForce].ppar_utilization[fn][fninst] += cSFC.trafficArrivalRate;
    }
    
    if(showInConsole == showFinal){
        cout << "\n BruteForce-Instance-Assignment:: Partial-Parallel: SFCid:" << cSFC.index <<"  ppId:" << singleSfcRes.ppar_pid << " | delay:[" << singleSfcRes.ppar_delay << "] ("<<singleSfcRes.ppar_load<<"%) :";
        for(const auto &blk: allPartParSFC[singleSfcRes.ppar_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout << fnid << char(96+singleSfcRes.ppar_fninstmap.at(fnid)) << " ";
            }   cout<<"]";
        }
    }
    return algosuccess;
}//bruteForce_PartParallel_Deployment



/* **************************************************************************************************************** */
/*! Recursive Function to construct all Layer Graph Nodes with their instances combination of a given stage (which consists of parallel VNFs).
 * Once a single instance combination is found, then construct lgNode with that instance combination and process that stage to precompute some of the delays, frequency of PNs .
 * It may be the case that some combination is not feasible then don't construct.
 * @param csfi function index of current stage.
 * @param curStg using to iterate all functions in the stage.
 * @param[in,out] curInstComb current combination construction in progress
 * @param[in,out] lgnids It stores all the stage wise lgNode indexes {stgid . all lgNodes id}
 * @param[in,out] idx2lgNode It stores mapping lgNode index to lgNode structure
 * @param[in] oldUtilization old utilization of the System on that basis we have to calculate the current combination feasibility.
 * @param[in] SimTest Simulation Object which contains SFC on which to perform the test.
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
 * stg 0 (1 function has 3 instances),    3 lgNodes  =  lgNode.instCombination (index,1d[pair<1a>])  lgNode.instCombination (index,[<1b>]) lgNode.instCombination (index,[<1c>])  \n
 * stg 1 (2 par function 2 & 3 instances),6 lgNodes  =  [<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c]  \n
 * stg 2 (1 function 2 instances),        2 lgNodes  =  [5a] [5b]
*/
void construct_layerGraphNodesIC(unsigned int csfi, const vector<unsigned int>& curStg, vector<pair<unsigned int,unsigned int>>& curInstComb,
                                 vector<unsigned int>& lgnids, unordered_map<unsigned int, lgNode>& idx2lgNode,
                                 const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,
                                 const Simulations& SimTest, const ServiceFunctionChain& cSFC
                                 ){//construct_layerGraphNodesIC

    if(csfi == curStg.size()){ // all functions in stage iterated. curStg.size()==numOfFunction in that stage.
        const unsigned int& lgNid = idx2lgNode.size(); ///< index to assign to new lgNode
        idx2lgNode[lgNid] = lgNode(lgNid, curInstComb); ///< create a node and mapping of index to lgNode
        lgnids.push_back(lgNid);
        lgNode& lgn = idx2lgNode[lgNid];
        /*! Processing of current lgNode: count of physical servers, max time in each server */
        type_delay utilization_sum_lgy=0, servicerate_sum_lgy=0;
        for (const auto &[d_fnType, d_fnInst]: lgn.instCombination) {
            const auto &d_pn_id = SimTest.I_VNFINST_2_PN.at(d_fnType).at(d_fnInst);
            const vnfDelaysPreComputed& fndelay = SimTest.vnfDelays.at(d_fnType);

            lgn.cntPN[d_pn_id].push_back(d_fnType); //< freq of PN in current lgn
            lgn.exePN[d_pn_id] = max(lgn.exePN[d_pn_id], fndelay.prcDelay + fndelay.exeDelay + fndelay.queuingDelay.at(d_fnInst));

            if (oldUtilization.count(d_fnType) and oldUtilization.at(d_fnType).count(d_fnInst)){ ///< old utilization till now of VNF
                utilization_sum_lgy += oldUtilization.at(d_fnType).at(d_fnInst);
            } servicerate_sum_lgy += SimTest.VNFNetwork.VNFNodes.at(d_fnType).serviceRate;
        }//for d_fnType, d_fnInst
        lgn.utilization = (utilization_sum_lgy/servicerate_sum_lgy)*100;
        return;
    }// base case

    const unsigned int& fn = curStg[csfi];
    const unsigned int& totInstancs = SimTest.finalInstancesCount.at(fn);
    for(unsigned int fninst=1; fninst<=totInstancs; fninst++){
        if (ifRejectSFC and oldUtilization.count(fn) and oldUtilization.at(fn).count(fninst)){ ///< old utilization till now of VNF
            if(oldUtilization.at(fn).at(fninst) + cSFC.trafficArrivalRate >  SimTest.VNFNetwork.VNFNodes.at(fn).serviceRate) {
                continue; // don't take this instance if its utilization become more than service rate of function.
            }
        }
        curInstComb.emplace_back(fn, fninst); // push current instance
        construct_layerGraphNodesIC(csfi+1, curStg, curInstComb, lgnids, idx2lgNode, oldUtilization, SimTest, cSFC); // call function for next instance
        curInstComb.pop_back(); // pop last instance and push next instance of same function for another combinations.
    }//recursive case
}//construct_layerGraphNodesIC

/*! Given partial SFC and Layer Graph (stg2lgnids, idx2lgNode), it traverse the layer graph to find the best mapping for given sfc with minimum time and min utilization.
 * In each stage we traverse atmost 3 (mxPathsK) paths from previous stage to current stage. Then mxPathsK*(stg size=num of lg nodes in stage) pairs of new paths are generated, out of which we again select atmost mxPathsK paths for next iteration.
 * Selection of paths are based on minimum time and min utilization of the path. (utilization of the path is max of any lgNode in that path)
 * @param partParSFC given partial parallel SFC
 * @param stg2lgnids mapping of stage wise lgNode indexes {stgid . all lgNodes id}
 * @param idx2lgNode mapping lgNode index to lgNode structure
 * @param[in,out] minBstDelay min delay found till now among all partial parallel chains of a given sfc.
 * @param minBstUtil min utilization of the path chosen with minBstDelay and among all paths with minBstDelay choose minUtil
 * @param[in,out] bstMapping to save best mapping overall among all partial parallel chain/and its all instances.
 * @param[in] SimTest Simulation Object which contains SFC on which to perform the test.
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 *
 */
void traverse_layerGraph(const vector<vector<unsigned int>>& partParSFC,  const vector<vector<unsigned int>>& stg2lgnids, unordered_map<unsigned int, lgNode>& idx2lgNode,
                          type_delay& minBstDelay, type_delay& minBstUtil, unordered_map<unsigned int,unsigned int>& curBstMapping,
                         const Simulations& SimTest, const ServiceFunctionChain& cSFC, const int& showInConsole = dontShow) {//traverse_layerGraph

    const unsigned int& szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.
    vector<unordered_set<unsigned int>> uniqLgidInStgOfKPaths(szStages); ///< unique lgNode ids in each stage which is used in k shortest path. So that we don't have to iterate all lgNode in stages again and itearte only which is necessary.
    type_delay  T_tx_init = calcD_TransmissionDelay(); ///< transmission time
    unsigned int lgDSTid = idx2lgNode.size()-1; ///< destination lgNode index
    const unsigned int mxPathsK = 5; ///< maximum number of paths to consider in any stage
/* ************* function to find min pair path ************************************************************ */
    std::function<unsigned int(unsigned int)> findPathsToConsider = [&](unsigned int pqsize)->unsigned int{
        if(pqsize <= 2) return pqsize; ///< if less than two then take both
        else if(pqsize <= 6) return 4; // if 3-6 paths then take 4
        else return mxPathsK; // from 7 to ---
//        return pqsize;
    };
    /*! to process the min heap. select min time and min utilization paths of all the paths present in min heap.
     * Early stopping criteria . At any point/stage if path's mindist is more than what we found the min delay for sfc till now among all partial chains,
        then don't consider this min path. Since next min path will also be more than minBst delay then no point in checking rest of the paths.
     * @param id_curstd cur stage id which we are processing
     * @param pq priority queue structure consist of all the paths of processed previous-current stage.
     * @param lastStage if we are processing last stage, then we just need one path and assign it to SFC.
     */
    std::function<int(const unsigned int, priority_queue<pqNode>&, unsigned int)>   processs_min_heap = [&idx2lgNode, &uniqLgidInStgOfKPaths, &minBstDelay, &showInConsole](const unsigned int& id_curstg, priority_queue<pqNode>&pq, unsigned int numOfPathsToConsider)->int{
        if(showInConsole == showFinal){ cout<<"\n TotalPathPairs:"<<pq.size();}
        int numOfPathsInserted = 0;
        type_delay prv_mindist_taken = 0;
        while(numOfPathsToConsider>0 and !pq.empty()){
            pqNode min_path = pq.top(); /*! pair of src and dst lg node which produce min distance */

            /*! Early stopping criteria. At any point/stage if path's mindist is more than what we found the min delay for sfc till now among all partial chains,
             * then don't consider this min path. Since next min path will also be more than minBst delay then no point in checking rest of the paths. */
            if(min_path.mindist > minBstDelay) {
                return numOfPathsInserted;
            }

            // min 3 path insert kar diye and hamare pass abhi kuch paths consider krne ko hai and pq me aur bhi paths hai, par abhi wale ka dist same hai prv wale se
             if(numOfPathsInserted>2 and (pq.size() > numOfPathsToConsider) and  fabs(prv_mindist_taken - min_path.mindist)<1){
                 pq.pop(); // is path ko consider hi nhi kar rha, aage ke 2 path ko consider kr rha
                 continue;
             }
            prv_mindist_taken = min_path.mindist;


            numOfPathsInserted++;
            min_path.path.push_back(min_path.y); ///< create path by inserting new destination for which it is minimum.

            idx2lgNode[min_path.y].kpaths.push_back(min_path); ///<consider this path traversing through this y index lgNode
            idx2lgNode[min_path.x].children.emplace_back(min_path.y, min_path.mindist); ///< from lgNode x source, where the paths traverse to lgNode y and its distance.

            uniqLgidInStgOfKPaths[id_curstg].insert(min_path.y); ///< if(!lastStage) last stage don't require insert but doesn't create a difference

            pq.pop(); numOfPathsToConsider--;
            if(showInConsole == showFinal){ cout<<"\n     p:"<<min_path.mindist<<"sec | "<<min_path.utilization<<"% ["; for(const auto& kkk: min_path.path)  cout<<kkk<<" -> ";}
        }
        return numOfPathsInserted;
    };
/* ************* stage wise calculation *************************************************************** */
    /*! Special Case: From Dummy SRC to First Stage(index 0).
     * First Stage consist of Layer Graph Node Indexes (lgnIdy). Each Layer Graph Node consist of instance combinations pairs{vnf type, its instance id}.
     * Calculate maximum delay taken to process that Layer Graph Node, as completion time would be when all instaces in that node finish their execution.
     * inter-duplication + transmission time +  max( "intra-duplication" + "time taken in each server" + "intra-merging")
     * Push the pairs into min priority queue to find pairs which produce minimum delay.
     */
    unsigned int id_curstg=0;
    { ///From SFC SRC to First Stage(index 0).
        priority_queue<pqNode> pqs; ///< to find minimum delay path in the x-y pairs of first and dummy source stage.

        for (const unsigned int &lgnIdy: stg2lgnids[id_curstg]) {
            const lgNode &lgy = idx2lgNode[lgnIdy]; ///< layer graph node y, x is dummy source
            if (showInConsole >= showDetailed) { cout << "\n lg:" << lgy.idx << " ["; for (const auto &givenPair: lgy.instCombination) { cout << givenPair.first << char(givenPair.second - 1 + 'a') << " "; } cout << "]";}

            type_delay mx_delay_x_y = 0;///< maximum delay of the current lgNode pair (x=dummySrc,y=current lgNode).
            /*! Calculation of packet processing time. inter duplication from src to different servers, intra duplication (within same server multiple nodes) and intra merging*/
            for (const auto &[pn_y, pn_y_parallel_vnfs]: lgy.cntPN) {
                type_delay T_d_hdr=0, T_m_hdr=0;
                if(pn_y_parallel_vnfs.size()>1){
                    unsigned int needIntraHdrCopy=0, needIntraPktCopy=0;
                    for(int j=1; j<pn_y_parallel_vnfs.size(); j++){
                        const auto vnfj = pn_y_parallel_vnfs[j]; bool plzcopy=false;
                        for(int i=0; i<j; i++){
                            if(SimTest.parallelPairs.at(pn_y_parallel_vnfs[i]).at(vnfj) == pktCopy){
                                plzcopy=true;
                                break;
                            }
                        }
                        if(plzcopy)needIntraPktCopy++;
                        else needIntraHdrCopy++;
                    }
                    T_d_hdr = calcD_IntraDuplicationTime(needIntraHdrCopy, needIntraPktCopy);
                    T_m_hdr = calcD_IntraMergingTime(needIntraHdrCopy, needIntraPktCopy);
                }
                const unsigned int pn_x = cSFC.access_nodes.first; /// sfc source physical node (only one node at previous stage)
                unsigned int px_py_same = 0;
                if (lgy.cntPN.find(pn_x) != lgy.cntPN.end()) px_py_same = 1;
                unsigned int cntNextHopDiffServer = lgy.cntPN.size() - px_py_same;
                type_delay T_d_pkt =  calcD_InterDuplicationTime(cntNextHopDiffServer); ///< inter duplication time from source to lgy
                type_delay T_tx=0, T_px=0;
                if(pn_x != pn_y){ /// if both server are different then there is transmission and propagation delay
                    T_tx = T_tx_init; T_px =  calcD_PropagationDelay(pn_x, pn_y, SimTest.PhysicalNetwork);
                }
                mx_delay_x_y = max(mx_delay_x_y, (T_d_pkt + T_tx + T_px)+(T_d_hdr + T_m_hdr) + lgy.exePN.at(pn_y));
            }//curStgPN

            if (showInConsole >= showDetailed) { cout << "\n      max:" << mx_delay_x_y; }
            pqs.emplace(idx2lgNode[0].idx, lgy.idx, mx_delay_x_y, lgy.utilization, vector<unsigned int>{idx2lgNode[0].idx}); /// constructing the dummy src to lgy path.
        }

        /*! Process paths from soruce to next stage. */
        if(processs_min_heap(id_curstg, pqs, findPathsToConsider((unsigned int)pqs.size())) <= 0){
            return; /// if min heap process stopped mid-way because early stopping crteira and not a single path taken for consideration then return. No sense of processing it further
        }
    }//From SFC SRC to First Stage(index 0).
/* ****************************************************************************************************** */
    /*! From each lgNode (inst combination) in current stage to lgNode(instance combination) in prevous stage. Repeat till last stage.
     *  Process the lgy node first: count physical servers, and maximum time of execution of server.
     *  Then for each lgx node in prev stage. Calculate inter-duplication + transmission + processing time. Take maximum of it for current lgy node.
     *  Out of all servers in lgy take maximum time to be delay for lgy.
     */
    for( id_curstg=1; id_curstg<szStages; id_curstg++) {/*!< Iterating stage ID from 0 to last index */
        unsigned int id_prvstg = id_curstg-1; ///< previous stage id to process
        priority_queue<pqNode> pq; ///< to find minimum delay path in the x-y pairs of previous and current source stage.

        /*  lgnIdy ********************************* */
        for (const unsigned int &lgnIdy: stg2lgnids[id_curstg]) { /*!< Iterating all the layer graph node index in current stage(destination)  */
            const lgNode& lgy = idx2lgNode[lgnIdy]; ///< current layer Layer Graph Node
            if (showInConsole >= showDetailed) { cout << "\n lgy:" << lgy.idx ; cout<<" ["; for(const auto& givenPair: lgy.instCombination){ cout<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  }  cout<<"]";}

            /*  lgnIdx*********************************** */
            for (const unsigned int &lgnIdx: uniqLgidInStgOfKPaths[id_prvstg]) {/*!< Iterating layer graph nodes index (which are in Path) in previous stage(source)  */
                const lgNode &lgx = idx2lgNode[lgnIdx];///< current layer Layer Graph Node
                if (showInConsole >= showDetailed) { cout << "\n   lgx:" << lgx.idx; cout << " ["; for (const auto &givenPair: lgx.instCombination) {  cout << givenPair.first << char(givenPair.second - 1 + 'a') << " "; } cout << "]";}

                /* *************************************** */
                /*! Once x = lgx, y = lgy are fixed. We will find maximum time for each physical server in y to determine edge (x,y). */
                type_delay mx_delay_x_y = 0;///< maximum delay of the current instance combination of the pair (x,y).

                for(const auto &[pn_y, pn_y_parallel_vnfs]: lgy.cntPN){ /*! for each physical server in cur node*/

                    /*! Calculation of packet processing time. inter mergring ( different server from cur server), intra duplication (within same server multiple nodes), intra Merging (within same server multiple nodes)*/
                    unsigned int py_px_same=0; if(lgx.cntPN.find(pn_y) != lgx.cntPN.end()) py_px_same =  1;
                    unsigned int cntPrevHopDiffServer = lgx.cntPN.size() - py_px_same;
                    type_delay T_m_pkt = calcD_InterMergingTime(cntPrevHopDiffServer);

                    type_delay T_d_hdr=0, T_m_hdr=0;
                    if(pn_y_parallel_vnfs.size()>1){
                        unsigned int needIntraHdrCopy=0, needIntraPktCopy=0;
                        for(int j=1; j<pn_y_parallel_vnfs.size(); j++){
                            const auto vnfj = pn_y_parallel_vnfs[j]; bool plzcopy=false;
                            for(int i=0; i<j; i++){
                                if(SimTest.parallelPairs.at(pn_y_parallel_vnfs[i]).at(vnfj) == pktCopy){
                                    plzcopy=true; break;
                                }
                            }
                            if(plzcopy)needIntraPktCopy++;
                            else needIntraHdrCopy++;
                        }
                        T_d_hdr = calcD_IntraDuplicationTime(needIntraHdrCopy, needIntraPktCopy);
                        T_m_hdr = calcD_IntraMergingTime(needIntraHdrCopy, needIntraPktCopy);
                    }

                    type_delay mx_pktPrc =  T_m_pkt  + T_d_hdr + T_m_hdr;///< overall total time spent in packet processing from src to dest.

                    if(showInConsole >= showDetailed) {
                        cout<<"\n      :"<<pn_y_parallel_vnfs.size()<<"[py:"<<pn_y<<"]"<<"  prvD:"<<cntPrevHopDiffServer<<" s("<<py_px_same<<")"
                            <<"   [m_pkt:"<<T_m_pkt<<" d_hdr:"<<T_d_hdr<<" m_hdr:"<<T_m_hdr<<"]"<<"   mxServer:"<<lgy.exePN.at(pn_y);
                    }
                    /*! Calculation of inter-duplication time transmission time and  propagation time (we can duplicate the packets right before sending them. */
                    type_delay mx_interdupTxPx_for_y = 0;
                    for (const auto &[pn_x, pn_x_parallel_vnfs]: lgx.cntPN) { /*! for each physical server in previous node*/
                        unsigned int px_py_same = 0;
                        if (lgy.cntPN.find(pn_x) != lgy.cntPN.end()) px_py_same = 1;
                        unsigned int cntNextHopDiffServer = lgy.cntPN.size() - px_py_same;
                        type_delay T_d_pkt =  calcD_InterDuplicationTime(cntNextHopDiffServer);
                        type_delay T_tx=0, T_px=0;
                        if(pn_y != pn_x){ /// if both server are different then there is transmission and propagation delay
                            T_tx = T_tx_init; T_px =  calcD_PropagationDelay(pn_x, pn_y,SimTest.PhysicalNetwork);
                        }
                        mx_interdupTxPx_for_y = max(mx_interdupTxPx_for_y, T_d_pkt + T_tx + T_px);///< overall total time spent in sending packet from src to dest.
                        if (showInConsole >= showDetailed) {
                            cout << "\n           " << pn_x_parallel_vnfs.size() << "(px:" << pn_x << ")  " << "nxtD:" << cntNextHopDiffServer << " s(" << px_py_same << ")"
                                 << "   [d_pkt:" << T_d_pkt << " tx:" << T_tx << " px:" << T_px << "]";
                        }
                    }//lgx prvStgPN
                    mx_delay_x_y = max(mx_delay_x_y, mx_interdupTxPx_for_y + mx_pktPrc + lgy.exePN.at(pn_y));

                    if(showInConsole >= showDetailed) {
                        cout<<"   mxTxPx:"<<mx_interdupTxPx_for_y;
                        cout<<"\n        px-py:"<< mx_interdupTxPx_for_y + mx_pktPrc + lgy.exePN.at(pn_y);
                    }

                }//lgy curStgPN

                for(const pqNode& kpq: lgx.kpaths) {/// for each min path in lgx node, push this pair also.
                    pq.emplace(lgnIdx, lgnIdy, kpq.mindist + mx_delay_x_y, max(kpq.utilization, lgy.utilization), kpq.path);
                } //kpq
            }//lgnIdx
        }//lgnIdy
        if(processs_min_heap(id_curstg, pq, findPathsToConsider((unsigned int)pq.size())) <= 0){
            return; /// if min heap process stopped mid-way because early stopping crteira and not a single path taken for consideration then return. No sense of processing it further
        }
    }//id_curstg

/*    From last Stage to destination of sfc **** ******************************************************* */
    { //From last Stage to destination of sfc
        id_curstg = szStages - 1;
        priority_queue<pqNode> pqd; ///< to find minimum delay path in the x-y pairs of last stage and dummy dst stage.
        const unsigned int& pn_y = cSFC.access_nodes.second; /// sfc destination physical node (only one node at last stage)
        for (const unsigned int &lgnIdx: uniqLgidInStgOfKPaths[id_curstg]) { //lgx
            const lgNode &lgx = idx2lgNode[lgnIdx];
            if (showInConsole >= showDetailed) {
                cout << "\n lg:" << lgx.idx << " ["; for (const auto &givenPair: lgx.instCombination) {  cout << givenPair.first << char(givenPair.second - 1 + 'a') << " ";  } cout << "]";
            }
            type_delay mx_interdupTxPx_for_y = 0;
            for (const auto &[pn_x, pn_x_parallel_vnfs]: lgx.cntPN) { /*! for each physical server in previous node*/

                type_delay T_tx=0, T_px=0;
                if(pn_x != pn_y){ /// if both server are different then there is transmission and propagation delay
                    T_tx = T_tx_init; T_px =  calcD_PropagationDelay(pn_x, pn_y,SimTest.PhysicalNetwork);
                }
                mx_interdupTxPx_for_y = max(mx_interdupTxPx_for_y, T_tx + T_px);///< overall total time spent in sending packet from src to dest.
                if (showInConsole >= showDetailed) {
                    cout << "\n           " << pn_x_parallel_vnfs.size() << "(px:" << pn_x << ")  "<< "   [tx:" << T_tx << " px:" << T_px << "]";
                }
            }//lgx prvStgPN
            unsigned int py_px_same=0; if(lgx.cntPN.find(pn_y) != lgx.cntPN.end()) py_px_same =  1;
            unsigned int cntPrevHopDiffServer = lgx.cntPN.size() - py_px_same;
            type_delay T_m_pkt = calcD_InterMergingTime(cntPrevHopDiffServer);

            type_delay mx_delay_x_y = mx_interdupTxPx_for_y + T_m_pkt;

            for (const pqNode &kpq: lgx.kpaths) {
                pqd.emplace(lgnIdx, lgDSTid, (kpq.mindist + mx_delay_x_y), max(kpq.utilization, lgx.utilization), kpq.path);
            }
        }//lgx
        if(processs_min_heap(id_curstg, pqd, 1) <= 0){
            return; /// if min heap process stopped mid-way because early stopping crteira and not a single path taken for consideration then return. No sense of processing it further
        }
    }//From last Stage to destination of sfc

    if(idx2lgNode[lgDSTid].kpaths[0].mindist < minBstDelay){ ///< 0th index is shortest one
        minBstDelay =  idx2lgNode[lgDSTid].kpaths[0].mindist;
        minBstUtil = idx2lgNode[lgDSTid].kpaths[0].utilization;
        for(const auto lgid: idx2lgNode[lgDSTid].kpaths[0].path){
//            if(idx2lgNode[lgid].instCombination.empty())continue;
            for(const auto &[fn, fninst]: idx2lgNode[lgid].instCombination){
                curBstMapping[fn]=fninst;
            }
        }
    }
    else return;

    if(showInConsole >= showDetailed){
        cout << "\n-- Layer Graph Nodes ::";
        cout<<"\nId\t"<<"kp\t"<<"Utlz\t\t"<<"maxExePN\t"<<"child";
        cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
        for(int i=0; i<idx2lgNode.size(); i++){ const lgNode& sec= idx2lgNode[i];
            cout<<"\n"<<sec.idx<<"  | "<<sec.kpaths.size()<<" | " <<sec.utilization<<" |\t [";
            for(const auto& val: sec.exePN){ cout<<"p"<<val.first<<":"<<val.second<<"s ";}cout<<"]\t[";
            for(const auto& val: sec.children){ cout<<"id"<<val.first<<": "<<val.second<<"s; "; } cout<<"] ";;
        }
    }///show lgNodes

    if(showInConsole == showFinal){ /// showing instance and its delay
        cout<<"    "<<"["; for(const auto &blk: partParSFC){ for(const auto& fnid: blk){ cout<<fnid<<char(96+curBstMapping.at(fnid))<<" ";  }  } cout<<"]";
        cout<<"["<<minBstDelay<<"] ("<<minBstUtil<<"%)";
    }

}//traverse_layerGraph

/*! FINDING THE SFC's VNF Instance Assignment IF PARALLELISM is DISABLED IN SFC\n
 * For a given sequential sfc, construct the layer graph and find the final instance mapping/deployment based on previous utilization.\n
 * Intance mapping is based on min delay of the path with min utilization (lightly loaded).
 * @param[in] SimTest Simulation Object which contains SFC on which to perform the test.
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * @param showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 * @return status of the algorithm.
 */
bool kShortestPath_Sequential_Deployement(Simulations& SimTest, const ServiceFunctionChain& cSFC, const int& showInConsole = dontShow) {//kShortestPath_Sequential_Deployement

    SFC_RESULT& singleSfcRes = SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index];
    const vector<vector<unsigned int>>& seqSFC = (*cSFC.partialParallelChains).front();
    const unsigned int& Len = seqSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    vector<vector<unsigned int>> stg2lgnids(Len); ///< stage wise lgNode indexes in order to process/traverse it.
    unordered_map<unsigned int, lgNode> idx2lgNode;///< given index it is mapped to actual layer graph node structer so that we can process it uniquely.
    idx2lgNode[0] = lgNode(0); ///< source lgNode.

    for(unsigned int stgid=0; stgid<Len; stgid++){      // finding instances possibilities of each stage and constructing lgNodes.
        const auto& curStg = seqSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_layerGraphNodesIC(0, curStg, curInstComb, stg2lgnids[stgid], idx2lgNode, SimTest.TestsResult[name_kshortestpath].seq_utilization, SimTest, cSFC );
        if(stg2lgnids[stgid].empty()){
            return singleSfcRes.seq_status = noResDueToNoStg; ///< if there is no instance combinations in current stage then we can't proceed
        }
    }//stgid<szStages finding instances possibilities of each stage.

    unsigned int lgDSTid = idx2lgNode.size(); ///< destination lgNode
    idx2lgNode[lgDSTid] = lgNode(lgDSTid);

    if(showInConsole >= showDetailed){ // showing partParSFC info
        cout<<"\nSequential-SFC: ";
        for(const auto& blks: seqSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<")";
        for(int cur_lvl=0; cur_lvl<Len; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2lgnids[cur_lvl].size()<<"):";
            for(const auto& lgid: stg2lgnids[cur_lvl]){
                cout<<lgid<<"["; for(const auto &givenPair: idx2lgNode[lgid].instCombination){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  } cout<<"] ";
            }
        }
    }

    traverse_layerGraph(seqSFC, stg2lgnids, idx2lgNode, singleSfcRes.seq_delay, singleSfcRes.seq_load,singleSfcRes.seq_fninstmap, SimTest, cSFC);
  
    /// if we found the mapping for sequential chain
    for(const auto& [fn, fninst]: singleSfcRes.seq_fninstmap){
        SimTest.TestsResult[name_kshortestpath].seq_utilization[fn][fninst] += cSFC.trafficArrivalRate;
//        cout<<"{"<<fn<<","<<fninst<<"},";
    }

    if(showInConsole >= showFinal) {
        cout << "\n kShortestPath-Instance-Assignment:: Sequential: SFCid:" << cSFC.index << " | delay:["<<singleSfcRes.seq_delay<<"] ("<<singleSfcRes.seq_load<<"%) :(";
        for(const auto& fnid : cSFC.vnfSeq){
            cout<<"f"<<fnid<<char(96+singleSfcRes.seq_fninstmap.at(fnid))<<"; ";
        } cout<<")";
    }
    return singleSfcRes.seq_status = algosuccess;
}//kShortestPath_Sequential_Deployement


/*! FINDING THE SFC's VNF Instance Assignment IF Full PARALLELISM Enabled IN SFC \n
 * For a given partial sfc and its full parallel vnf, construct the layer graph and find the final instance mapping/deployment based on previous utilization.\n
 * Intance mapping is based on min delay of the path with min utilization (lightly loaded).\n 
 * @param[in] SimTest Simulation Object which contains SFC on which to perform the test.
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * @param[in] showInConsole output in console (0 don't show, 1 final output, 2 detailed output) 
 * @return status of the algorithm.
 */
bool kShortestPath_FullParallel_Deployment(Simulations& SimTest, const ServiceFunctionChain& cSFC,
                                            const int& showInConsole = dontShow) {//kShortestPath_FullParallel_Deployment
    const vector<vector<unsigned int>>& fullSFC = cSFC.vnfBlocksPar;
    SFC_RESULT& singleSfcRes = SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index];

    const unsigned int Len = fullSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    vector<vector<unsigned int>> stg2lgnids(Len); ///< stage wise lgNode indexes in order to process/traverse it.
    unordered_map<unsigned int, lgNode> idx2lgNode;///< given index it is mapped to actual layer graph node structer so that we can process it uniquely.
    idx2lgNode[0] = lgNode(0); ///< source lgNode.

    unsigned long long totalInstComb = 0;
    for(unsigned int stgid=0; stgid<Len; stgid++){      // finding instances possibilities of each stage and constructing lgNodes.
        const auto& curStg = fullSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_layerGraphNodesIC(0, curStg, curInstComb, stg2lgnids[stgid], idx2lgNode, SimTest.TestsResult[name_kshortestpath].fullpar_utilization, SimTest, cSFC );
        if(stg2lgnids[stgid].empty()){
            return singleSfcRes.fullpar_status = noResDueToNoStg; ///< if there is no instance combinations in current stage then we can't proceed
        }
        if(stgid != 0){
            totalInstComb += min(stg2lgnids[stgid-1].size(), 5ULL)*stg2lgnids[stgid].size();
        }else totalInstComb = stg2lgnids[stgid].size();
    }//stgid<szStages finding instances possibilities of each stage.

    cout<<"\n  H-fully-parallel["<<SimTest.sfccompleted<<"/"<<SimTest.sortedSFCs.size()<<"]:"<<totalInstComb;

    unsigned int lgDSTid = idx2lgNode.size(); ///< destination lgNode
    idx2lgNode[lgDSTid] = lgNode(lgDSTid);

    if(showInConsole == showDetailed){ // showing partParSFC info
        cout<<"\nFully-Parallel-SFC: ";
        for(const auto& blks: fullSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<")";
        for(int cur_lvl=0; cur_lvl<Len; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2lgnids[cur_lvl].size()<<"):";
            for(const auto& lgid: stg2lgnids[cur_lvl]){
                cout<<lgid<<"["; for(const auto &givenPair: idx2lgNode[lgid].instCombination){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  } cout<<"] ";
            }
        }
    }

    traverse_layerGraph(fullSFC, stg2lgnids, idx2lgNode, singleSfcRes.fullpar_delay, singleSfcRes.fullpar_load, singleSfcRes.fullpar_fninstmap, SimTest, cSFC);
 

    /// if we found the mapping for sequential chain
    for(const auto& [fn, fninst]: singleSfcRes.fullpar_fninstmap){
        SimTest.TestsResult[name_kshortestpath].fullpar_utilization[fn][fninst] += cSFC.trafficArrivalRate;
    }

    if(showInConsole == showFinal) {
        cout << "\n kShortestPath-Instance-Assignment: Fully-Parallel: SFCid:" << cSFC.index << " |  delay:["<<singleSfcRes.fullpar_delay<<"] ("<<singleSfcRes.fullpar_load<<"%) :";
        for(const auto& blk: cSFC.vnfBlocksPar){
            cout<<" ["; for(int fnid: blk){
                cout << "f" << fnid << char(96 + singleSfcRes.fullpar_fninstmap.at(fnid)) << " ";
            } cout<<"]";
        }
    }
    return singleSfcRes.fullpar_status = algosuccess;
}//kShortestPath_FullParallel_Deployment

/*! FINDING THE Optimal SFC's VNF Instance Assignment in case IF PARALLELISM Enabled IN SFC.\n
 * For a given SFC and its all partial parallel chains, construct the layer graph and find the final instance mapping/deployment based on previous utilization.\n
 * Intance mapping is based on min delay of the path with min utilization (lightly loaded).\n
 * For a partial parallel chain, and its instances combination in stages, it enumerates only k shortest paths in each stage from previous to current stage. 
 * @param[in] SimTest Simulation Object which contains SFC on which to perform the test. 
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * @param[in] showInConsole output in console (0 don't show, 1 final output, 2 detailed output)
 * @return status of the algorithm.
 */
bool kShortestPath_PartParallel_Deployement(Simulations& SimTest, const ServiceFunctionChain& cSFC,
                                            const int& showInConsole = dontShow) {//kShortestPath_PartParallel_Deployement

    SFC_RESULT& singleSfcRes = SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index];
    const vector<vector<vector<unsigned int>>>& allPartParSFC = *cSFC.partialParallelChains;
    
    type_delay bst_partpar_delay = std::numeric_limits<type_delay>::max();
    //    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}}; ///Testing purpose
    for(int ppsidx=int( allPartParSFC.size())-1; ppsidx>=0; --ppsidx){/*! For each partial parallel SFC without src and dest block/stage.*/

//      cout<<"\r\t\t\t("<<allPartParSFC.size()-ppsidx<<"/"<<allPartParSFC.size()<<")"; //showing progress
        const vector<vector<unsigned int>>& partParSFC=allPartParSFC[ppsidx]; ///< for each of the partial parallel SFC of the givenParVNF Blocks
        const unsigned int szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

        vector<vector<unsigned int>> stg2lgnids(szStages); ///< stage wise lgNode indexes in order to process/traverse it.
        unordered_map<unsigned int, lgNode> idx2lgNode;///< given index it is mapped to actual layer graph node structer so that we can process it uniquely.
        idx2lgNode[0] = lgNode(0); ///< source lgNode.

        unsigned long long totalInstComb = 0;
//        unsigned int mxNumOflgNodesInStg = 0; // maximum number of lgNodes present in any stage;
        for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage and constructing lgNodes.
            const auto& curStg = partParSFC[stgid];
            vector<pair<unsigned int,unsigned int>> curInstComb;
            construct_layerGraphNodesIC(0, curStg, curInstComb, stg2lgnids[stgid], idx2lgNode, SimTest.TestsResult[name_kshortestpath].ppar_utilization, SimTest, cSFC );
//            mxNumOflgNodesInStg=max(mxNumOflgNodesInStg, (unsigned int)(stg2lgnids[stgid].size()));

            if(stg2lgnids[stgid].empty()){
                return singleSfcRes.ppar_pid = noResDueToNoStg; ///< if there is no instance combinations in current stage then we can't proceed
            }
            if(stgid != 0){
                totalInstComb += min(stg2lgnids[stgid-1].size(), 5ULL)*stg2lgnids[stgid].size();
            } else totalInstComb = stg2lgnids[stgid].size();
        }//stgid<szStages finding instances possibilities of each stage.
        if(stg2lgnids.size() != szStages)continue; ///< if number of stages here are not same as sfc then no need to proceed.

        cout<<"\r  H-partial-parallel["<<SimTest.sfccompleted<<"/"<<SimTest.sortedSFCs.size()<<"]("<<allPartParSFC.size()-ppsidx<<"/"<<allPartParSFC.size()<<"):"<<totalInstComb<<"       ";

        unsigned int lgDSTid = idx2lgNode.size(); ///< destination lgNode
        idx2lgNode[lgDSTid] = lgNode(lgDSTid);

        if(showInConsole >= showDetailed){ // showing partParSFC info
            cout<<"\nPartial-Parallel-SFC["<<ppsidx<<"]: ";
                for(const auto& blks: partParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<")";
                for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
                    cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2lgnids[cur_lvl].size()<<"):";
                    for(const auto& lgid: stg2lgnids[cur_lvl]){
                        cout<<lgid<<"["; for(const auto &givenPair: idx2lgNode[lgid].instCombination){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  } cout<<"] ";
                    }
                }
        }

        traverse_layerGraph(partParSFC, stg2lgnids, idx2lgNode, bst_partpar_delay, singleSfcRes.ppar_load, singleSfcRes.ppar_fninstmap, SimTest, cSFC, showInConsole);
        if(bst_partpar_delay < singleSfcRes.ppar_delay ){
            singleSfcRes.ppar_pid = ppsidx;
            singleSfcRes.ppar_delay = bst_partpar_delay;
        }
    }//for each allPartParSFC ppsidx.
//    if(singleSfcRes.ppar_pid == noResPar){
//        return algostopped;
//    }
    /// if we found the mapping for parallel chain
    for(const auto& [fn, fninst]: singleSfcRes.ppar_fninstmap){
        SimTest.TestsResult[name_kshortestpath].ppar_utilization[fn][fninst] += cSFC.trafficArrivalRate;
    }

    if(showInConsole == showFinal) {
        cout << "\n kShortestPath-Instance-Assignment:: Partial-Parallel: SFCid:" << cSFC.index << "  ppId:" << singleSfcRes.ppar_pid << " | delay:[" << singleSfcRes.ppar_delay << "] ("<<singleSfcRes.ppar_load<<"%) :";
        for(const auto &blk: allPartParSFC[singleSfcRes.ppar_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout << fnid << char(96+singleSfcRes.ppar_fninstmap.at(fnid)) << " ";
            }   cout<<"]";
        }
    }

    return algosuccess;
}//kShortestPath_PartParallel_Deployement

/* **************************************************************************************************************** */
/*! Find SFC's VNF deployment for sequential sfc only, for fully parallel sfc only, for partial-parallel sfc only.
 * @param SimTest Simulation Object which contains SFC on which to perform the test.
 */
void Heuristic_kShortestPath_InstanceMapping(Simulations& SimTest) {//Heuristic_kShortestPath_InstanceMapping

    SimTest.TestsResult[name_kshortestpath] =  SimTEST(name_kshortestpath); ///< simulation test object
    cout<<"\n ["<<name_kshortestpath<<"]\n";
    
// SEQ ----------------------------------------------------------
    SimTest.sfccompleted=0;
    for(const ServiceFunctionChain*const& sfcpointer:SimTest.sortedSFCs) {
        const ServiceFunctionChain& cSFC = *sfcpointer;
        SimTest.sfccompleted++;
        cout<<"\r  H-sequential["<<SimTest.sfccompleted<<"/"<<SimTest.sortedSFCs.size()<<"]";

        SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index] = SFC_RESULT(); ///< solution for the layer graph algorithm
        auto sfc_st = std::chrono::steady_clock::now();
        for(const unsigned int& fn: cSFC.vnfSeq){ /// finding some vnf delays
            const VNFNode& dstVNFNode = SimTest.VNFNetwork.VNFNodes.at(fn );
            SimTest.vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF(dstVNFNode);
            SimTest.vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay(dstVNFNode);
            for(int fnInst=1; fnInst<=SimTest.finalInstancesCount.at(fn); fnInst++) { ///sequential Queuing Delay
                SimTest.vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay(cSFC.trafficArrivalRate, dstVNFNode, fnInst, SimTest.TestsResult[name_kshortestpath].seq_utilization);
            }
        }

        kShortestPath_Sequential_Deployement(SimTest, cSFC);
        SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        SimTest.TestsResult[name_kshortestpath].total_seq_delay += SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_delay;
        SimTest.TestsResult[name_kshortestpath].total_seq_duration += SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].seq_duration;
    }

// FULL ----------------------------------------------------------
    SimTest.sfccompleted=0;
    for(const ServiceFunctionChain*const& sfcpointer:SimTest.sortedSFCs) {
        const ServiceFunctionChain& cSFC = *sfcpointer;
        SimTest.sfccompleted++; 
        
        SimTest.vnfDelays.clear();
        auto sfc_st = std::chrono::steady_clock::now();
        for(const unsigned int& fn: cSFC.vnfSeq){ /// finding some vnf delays
            const VNFNode& dstVNFNode = SimTest.VNFNetwork.VNFNodes.at(fn);
            SimTest.vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF(dstVNFNode);
            SimTest.vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay(dstVNFNode);
            for(int fnInst=1; fnInst<=SimTest.finalInstancesCount.at(fn); fnInst++) { ///sequential Queuing Delay
                SimTest.vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay(cSFC.trafficArrivalRate, dstVNFNode, fnInst, SimTest.TestsResult[name_kshortestpath].fullpar_utilization);
            }
        }
        kShortestPath_FullParallel_Deployment(SimTest, cSFC);
        SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].fullpar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        SimTest.TestsResult[name_kshortestpath].total_fullpar_delay += SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].fullpar_delay;
        SimTest.TestsResult[name_kshortestpath].total_fullpar_duration += SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].fullpar_duration;
    }

// AllPART ----------------------------------------------------------
    SimTest.sfccompleted=0; cout<<"\n";
    for(const ServiceFunctionChain*const& sfcpointer:SimTest.sortedSFCs) {
        const ServiceFunctionChain& cSFC = *sfcpointer;
        SimTest.sfccompleted++;  cout<<"\n";
        
        SimTest.vnfDelays.clear();
        auto sfc_st = std::chrono::steady_clock::now();
        for(const unsigned int& fn: cSFC.vnfSeq){ /// finding some vnf delays
            const VNFNode & dstVNFNode = SimTest.VNFNetwork.VNFNodes.at(fn);
            SimTest.vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF(dstVNFNode);
            SimTest.vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay(dstVNFNode);
            for(int fnInst=1; fnInst<=SimTest.finalInstancesCount.at(fn); fnInst++) { ///Parallel Queuing Delay
                SimTest.vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay(cSFC.trafficArrivalRate, dstVNFNode, fnInst, SimTest.TestsResult[name_kshortestpath].ppar_utilization);
            }
        }

        kShortestPath_PartParallel_Deployement(SimTest, cSFC);
        SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].ppar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        SimTest.TestsResult[name_kshortestpath].total_ppar_delay += SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].ppar_delay;
        SimTest.TestsResult[name_kshortestpath].total_ppar_duration += SimTest.TestsResult[name_kshortestpath].sfcsol[cSFC.index].ppar_duration;
    }
}//Heuristic_kShortestPath_InstanceMapping

/*! Find SFC's VNF deployment for sequential sfc only, for fully parallel sfc only, for partial-parallel sfc only.
 * @param SimTest Simulation Object which contains SFC on which to perform the test. 
 */
void bruteForce_InstanceMapping(Simulations& SimTest) {//bruteForce_InstanceMapping

    SimTest.TestsResult[name_bruteForce] = SimTEST(name_bruteForce); ///< simulation test object
    cout<<"\n ["<<name_bruteForce<<"]\n";

// SEQ ----------------------------------------------------------
    SimTest.sfccompleted=0;
    for(const ServiceFunctionChain*const& sfcpointer:SimTest.sortedSFCs) {
        const ServiceFunctionChain& cSFC = *sfcpointer;
        SimTest.sfccompleted++;
        cout<<"\r  BF-sequential["<<SimTest.sfccompleted<<"/"<<SimTest.sortedSFCs.size()<<"]";

        SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index] = SFC_RESULT(); /// using this algo what is optimal result for a single sfc.

        auto sfc_st = std::chrono::steady_clock::now();
        bruteForce_Sequential_Deployment(SimTest, cSFC);
        SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        SimTest.TestsResult[name_bruteForce].total_seq_delay +=  SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_delay;
        SimTest.TestsResult[name_bruteForce].total_seq_duration += SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].seq_duration;
    }

// FULL ----------------------------------------------------------
    SimTest.sfccompleted=0;
    for(const ServiceFunctionChain*const& sfcpointer:SimTest.sortedSFCs) {
        const ServiceFunctionChain& cSFC = *sfcpointer;
        SimTest.sfccompleted++;

        auto sfc_st = std::chrono::steady_clock::now();
        bruteForce_FullParallel_Deployment(SimTest, cSFC);
        SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].fullpar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        SimTest.TestsResult[name_bruteForce].total_fullpar_delay +=  SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].fullpar_delay;
        SimTest.TestsResult[name_bruteForce].total_fullpar_duration += SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].fullpar_duration;
    }

// ALL PART----------------------------------------------------------
    SimTest.sfccompleted=0; cout<<"\n";
    for(const ServiceFunctionChain*const& sfcpointer:SimTest.sortedSFCs) {
        const ServiceFunctionChain& cSFC = *sfcpointer;
        SimTest.sfccompleted++;  cout<<"\n";

        auto sfc_st = std::chrono::steady_clock::now();
        bruteForce_PartParallel_Deployment(SimTest, cSFC);
        SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].ppar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        SimTest.TestsResult[name_bruteForce].total_ppar_delay +=  SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].ppar_delay;
        SimTest.TestsResult[name_bruteForce].total_ppar_duration +=  SimTest.TestsResult[name_bruteForce].sfcsol[cSFC.index].ppar_duration;
    }

}//bruteForce_InstanceMapping

/* **************************************************************************************************************** */

#endif //SFC_PARALLELIZATION_ALGORITHMS_H

