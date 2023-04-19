//
// Created by vijay on 01-04-2023.
//

#ifndef SFC_PARALLELIZATION_ALGORITHMS_H
#define SFC_PARALLELIZATION_ALGORITHMS_H

/*!
 * Dynamic Program to search all the available sequential cluster sizes for a Chain C of length K
 * @param k small k is starting point, from which you have to calculate, as some will be precalculated and saved.
 * @param K chain length which is maxSFClen, that is upto what you want to calculate all possible enumeration
 * @param showInConsole to show the output in console. Time will increase to 2000-2500 ms
 * @return updated clusterSz[i] for k <= i <=K.\n
 *  K[1] total(1), {{1}} \n
 *  K[2] total(2), {{1,1},{2}} \n
 *  K[3] total(4), {{1,1,1},{2,1},{1,2},{3}} \n
 *  K[4] total(8), {{1,1,1,1},{2,1,1},{1,2,1},{3,1},{1,1,2},{2,2},{1,3},{4}}
 */
void clusterSizeEnumeration(unsigned int k, unsigned int K, bool showInConsole=false) {
    if(k>K) return;
    if(clusterSz.find(k-1) == clusterSz.end()){
        string errorMsg = "Previous Cluster size k-1="+to_string(k-1)+ " does not exist. Caclculate that first. Function: ";
        throw runtime_error(errorMsg+ __FUNCTION__);
    }
//    vector<vector<vector<int>>> clusterSz(K+1);
    for(unsigned int ki = k; ki<=K; ki++){
        vector<vector<unsigned int>> ans;
        for(unsigned int i=1; i<=ki; i++){
            for(auto s_dash: clusterSz[ki-i]){
                s_dash.push_back(i);
                ans.push_back(s_dash);
            }
        }
        clusterSz[ki] = std::move(ans);
    }

    if(showInConsole){
        cout<<"\nclusterSize={";
        for(unsigned int ki = k; ki<=K; ki++){
            cout<<"\n\tK["<<ki<<"] total("<<clusterSz[ki].size()<<"), {";
            for(const auto& x: clusterSz[ki]) {
                cout<<"{"; for(unsigned int i=0; i<x.size()-1; i++) cout<<x[i]<<","; cout<<x[x.size()-1]<<"},";
            }  cout<<"}";
        } cout<<"\n}end;";
    }
}

/*! Calculate all possible combination vectors n Choose k for given n and k.
 * @param n starting point from which you have to calculate new value, as some will be precalculated and saved.
 * @param N maximum value upto which you want to calculate all possible combination
 * @param showInConsole  to show the output in console. Time will increase to 2000-2500 ms
 * @return updated nCk[i] for n <= i <=N.\n
 *      n=1, {k=1, [1] }\n
 *      n=2, {k=1, [{1},{2}]  |         k=2, [{1,2}]   }\n
 *      n=3, {k=1, [{1},{2},{3} |       k=2, [{1,2},{1,3},{2,3}] | {k=3, [{1,2,3}] }\n
 *      n=4, {k=1, [{1},{2},{3},{4}] |  k=2, [{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}] | {k=3, [{1,2,3},{1,2,4},{1,3,4}, {2,3,4}] | k=4, [{1,2,3,4}] }
 */
void all_nCk(unsigned int n, unsigned int N, bool showInConsole=false){
    if(n>N) return;
    // Lambda Function to calculate nCk using backtracking application.
    std::function<void(unsigned int,unsigned int,vector<unsigned int>&,unsigned int&,unsigned int&)> combineHelper = [&combineHelper]
    (unsigned int st, unsigned int k, vector<unsigned int>& cur, unsigned int &x, unsigned int &y) ->void
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
    if(showInConsole){
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

/* **************************************************************************************************************** */

/*! It reads data from sequential vector and then convert into full Parallel SFC (parallel VNF blocks) chain by detecting
 * whether two function can be parallelised or not.\n
 * Each block denote fully parallel VNFs in that step.
 * @param[in, out] SFC ServiceFunctionChain object to convert sequential into parallel.
 * @param VNFNetwork object of VirtualNetworkFunctions class required to check two vnf are parallelizable or not.
 */
template<typename type_res>
void convert_SeqSFC2ParVNFBlocks(ServiceFunctionChain* const SFC, const VirtualNetworkFunctions<type_res> * const VNFNetwork, bool isLast=false){
    unsigned int sz = SFC->vnfSeq.size(); // total vnfs including src and dest
    const unordered_map<unsigned int, unordered_set<unsigned int>>& parallelPairs = VNFNetwork->parallelPairs;
    unsigned int cid=0;
    SFC->vnfBlocksPar.push_back({SFC->vnfSeq[cid]}); // stage 0, pushing src node
    for(cid=1; cid<sz; cid++){ // from SFCsrc+1 stg to SFCdst-1 stage
        unsigned int prv_vnf = SFC->vnfSeq[cid-1], cur_vnf = SFC->vnfSeq[cid];  // Checking prv_vnf --> cur_vnf pairs
        // NOT PARALLEL, if prv_vnf does not exist as first vnf in pair or if prv_vnf exist it is not parallel to cur_vnf, then it is not parallel pair
        if(parallelPairs.count(prv_vnf)==0 or parallelPairs.at(prv_vnf).find(cur_vnf) == parallelPairs.at(prv_vnf).end()){
            SFC->vnfBlocksPar.push_back({cur_vnf});// push cur_vnf as separte new stage
        }else if(SFC->vnfBlocksPar.back().size() == 1) { //PARALLEL Pairs STG Size=1 and size of previous stg is just one, then we can directly push into that stg.
            SFC->vnfBlocksPar.back().push_back(cur_vnf);
        }else{ // PARALLEL Pairs and STG size > 1, from previous stg we have to check if cur_vnf is parallel to all prv_vnf in previous stage
            bool pushInLstStg = true;
            for(const unsigned int& lst_stg_vnf: SFC->vnfBlocksPar.back()){
                if(parallelPairs.count(lst_stg_vnf)==0 or parallelPairs.at(lst_stg_vnf).find(cur_vnf) == parallelPairs.at(lst_stg_vnf).end()){
                    pushInLstStg = false;   break;
                }
            }
            if(pushInLstStg) SFC->vnfBlocksPar.back().push_back(cur_vnf);
            else SFC->vnfBlocksPar.push_back({cur_vnf});
        }//lst stg size>1
    }// for cid
    if(debug and isLast)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";
}

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
void assign_Clusters2ParVNFs(ServiceFunctionChain* const SFC, bool showInConsole = false, bool showInConsoleDetailed = false){
    const vector<vector<unsigned int>>& fullParVNFBlocks = SFC->vnfBlocksPar;
    const unsigned int& nBlk = fullParVNFBlocks.size(); ///< number of blocks of the fully parallel SFC, (including src and dst block)
    vector<vector<vector<unsigned int>>>& allPartParSFC = SFC->allPartParSFC; ///< All the Partial parallel Clusters of the fully parallel VNF Blocks. Each SFC is without src and dest.
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


    if(showInConsoleDetailed){
        SFC->showParallelSFC(fullParVNFBlocks);
    }
    for(const vector<unsigned int>& cluster: clusterSz[SFC->numVNF]){ //
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

/*! generate all the feasible partial parallel SFC for the given full parallel SFC.
 * @param nVNFs is number of VNFs except src and dest in fully parallel SFC.
 * @param clusterSz for each of the cluster enumeration for size[nVNFs].
 * @param nCk n=block size and k={cluster_i[l] value i.e. l(level) index value dentoes number of function(k) to be chosen as parallel out of all functions(n) in blocks.
 * @param fullParVNFBlocks is fully parallel VNF Blocks in sequence where each block/step denotes all the parallelizable functions in that block/step.
 * @return allPartParSFC All the Partial parallel Clusters of the fully parallel VNF Blocks. Each SFC is without src and dest.
 * @ForExample: fullParVNFBlocks={{1},{2,3,4,5}}, nVNFs:5 (f1,f2,f3,f4,f5).\n
 * Blk(0) = {1} only 1 parallel function, Blk(1) = {2,3,4,5} all 4 are parallel.\n
 * clusterSz[5(size=nVNFs)] = {{1, 2,1,1},{1,4},{3,1,1} ... } so on. cluster_i[l] denotes number of function parallel in that level. \n
 * where cluster_i {1, 2, 1, 1} means in level[0] only one function runs, level[2] = 2 functions run together, and level[3] and level[4] one-one function are there \n
 * {[1], [2,1,1]} is mapped to {[f1], [f2,f3,f4,f4]} such that [ 1-c combination of f1 ]-> [2-c combination of f2,f3,f4,f5] -> [1-c combination of f2,f3,f4,f5] -> [1-c combination of f2,f3,f4,f5] \n
 */
void find_PartialParalleSFCs(bool showInConsole = false, bool showInConsoleDetailed = false){
    auto numVNF = 4;
    const vector<vector<unsigned int>>& fullParVNFBlocks = {{1,2,3},{4}};// {{1,2,3},{4}}
    const unsigned int& nBlk = fullParVNFBlocks.size(); ///< number of blocks of the fully parallel SFC, (including src and dst block)
    vector<vector<vector<unsigned int>>> allPartParSFC; ///< All the Partial parallel Clusters of the fully parallel VNF Blocks. Each SFC is without src and dest.
//    const vector<vector<unsigned int>> SK = { {1,  2,1,1},{1, 2,2}, {1,4}, {2,2,1}, {3,1,1} };
    if(showInConsoleDetailed){
        cout << "(SRC ->";
        for(const auto& blk: fullParVNFBlocks){
            cout<<" ["; for(int fn: blk){ cout <<"f"<< fn <<"; "; } cout<<"] ->";
        } cout << " DST)";
    }

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



    for(const vector<unsigned int>& cluster: clusterSz[numVNF]){ //

//        if(cluster[0] > fullParVNFBlocks[0].size()) // if in first block(after src blk) 2 function is there, but cluster saying 3 needs to be parallel then continue next cluster.
//            continue;

        unsigned int cur_level=0; ///< current level at which we have to insert all the nodes
        unordered_map<unsigned int,vector<unsigned int>> levelInfo;///< this stores level wise info {level i -> {0 -> blkId, 1->n, 2-> k}} to find nCk for block blkId of parSFC_Full


        if(showInConsoleDetailed){
            cout<<"\ncluster["; for(const auto& x: cluster) cout<<x<<" "; cout<<"]";
        }
        bool allBlksOfSFCDone = true;
        // iterate through the blocks and map it to cur cluster.
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

/* **************************************************************************************************************** */ 
/*! Function to find all instances combination of parVNF in that stage. -> find stage to instnace combination of given block of sfc.
     * @param csfi function index of current stage.
     * @param curInstComb current combination in iteration
     * @param stgid stgId/blockId for which we are finding combination of functions in that stg/block and to store in stg2InstCombinations.
     * @param curStg using to iterate all functions in the stage.
     * @param[out] stg2InstCombinations It stores all the stage wise instances combination of all stage in partParSFC. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }}
     * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
     * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
     * stg 0 (1 function has 3 instances),     B[0] = 2d{  1d[ pair<1a> ] [<1b>] [<1c>]  } \n
     * stg 1 (2 par function 2 & 3 instances), B[1] = 2d{ 1d[<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c] } \n
     * stg 2 (1 function 2 instances),         B[2] = 2d{ 1d[5a] [5b] [5c] } \n
     * Time to calculte stg2InstCombinations -> if in any block number of parallel functions are 10 and each have 5 max instances\n
        inst = 2 (exe time: 1-2ms) (possibilites: 1024 (2^10)) \n
        inst = 3 (exe time: 38-40ms) (possibilites: 59 049 (3^10))\n
        inst = 4 (exe time: 580-600ms) (possibilites: 10 48 576 )\n
        inst = 5 (exe time: 5700-5800ms) (possibilites: 97 65 625)\
     */
template<typename type_res>
void construct_stageInstanceCombination(unsigned int csfi, vector<pair<unsigned int,unsigned int>>& curInstComb,  unsigned int& stgid, const vector<unsigned int>& curStg,
                                           unordered_map<unsigned int,  vector<vector<pair<unsigned int,unsigned int>>>> &stg2InstCombinations, const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,
                                           ServiceFunctionChain *const cSFC, const VirtualNetworkFunctions<type_res> *const VNFNetwork){ //construct_stageInstanceCombination
    
    if(csfi == curStg.size()){ // all functions in stage iterated. curStg.size()==numOfFunction in that stage.
        stg2InstCombinations[stgid].push_back(curInstComb); // push the one answer into combination stg.
        return;
    }
    const unsigned int fnType = curStg[csfi];
    unsigned int totInstancs = VNFNetwork->VNFNodes.at(fnType)->numInstances;
    for(unsigned int fnInstId=1; fnInstId<=totInstancs; fnInstId++){
        if (oldUtilization.count(fnType) and oldUtilization.at(fnType).count(fnInstId) and ///< old utilization till now of VNF
            (oldUtilization.at(fnType).at(fnInstId) + cSFC->trafficArrivalRate >  VNFNetwork->VNFNodes.at(fnType)->serviceRate))  {
            continue; // don't take this instance if its utilisation become more than service rate of function.
        }
        curInstComb.emplace_back(fnType, fnInstId); // push current instance
        construct_stageInstanceCombination(csfi+1, curInstComb, stgid, curStg, stg2InstCombinations, oldUtilization, cSFC, VNFNetwork); // call function for next instance
        curInstComb.pop_back(); // pop curInstComb instance and push next instance of same function.
    }
}//construct_stageInstanceCombination

/*! For a given stg2InstCombinations, it enumerate all the possible mappings we can give in each stage and calculate delay on the go.
     * @param stgid stgId/blockId for which we are enumerating instances.
     * @param curMapping cur function->instance mapping we iterating out of all possibilites.
     * @param bstMapping to save best mapping overall among all partial parallel chain/and its all instances.
     * @param minBstDelay min delay among all partial parallel chain/and its all instances.
     * @param partParSFC given partial SFC
     * @param stg2InstCombinations It contains all the stage wise instances combination of all stage in partParSFC. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }}
     * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
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
     * 1 stage -> 10 parallel func each with 5 max instances -> 5^10 possibilities or 97,65,625 instances.
     */
template<typename type_wgt, typename type_res>
void enumerate_InstanceCombination_CalcDelay(unsigned int stgid, unordered_map<unsigned int,unsigned int>& curMapping, unordered_map<unsigned int,unsigned int>& curBstMapping, type_delay& minBstDelay, 
                                 const vector<vector<unsigned int>>& partParSFC,  unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> >& stg2InstCombinations,
                                 const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,
                                 ServiceFunctionChain * const cSFC, VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
                                 bool showInConsoleDetailed = false){//enumerate_InstanceCombination_CalcDelay
    if(stgid == partParSFC.size()) { // found one mapping then find corresponding delay 
        type_delay parSfcCost =  parSfcCost = calcD_ParallelSFC<type_wgt, type_res>(partParSFC, curMapping, oldUtilization,  cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);
        if( parSfcCost < minBstDelay){ // current mapping ka delay is less than min delay among partial sfc all instances.
            minBstDelay =  parSfcCost;
            curBstMapping = curMapping;
        }
        else
            return;
        if(showInConsoleDetailed){ /// showing instance and its delay
            cout<<"\n\t"<<"["; for(const auto &blk: partParSFC){ for(const auto& fnid: blk){ cout<<fnid<<char(96+curMapping.at(fnid))<<" ";  }  } cout<<"]";
            cout<<"["<<parSfcCost<<"sec]";
        }
        return;
    }
    for(const vector<pair<unsigned int,unsigned int>>& instComb: stg2InstCombinations[stgid]){
        for(const auto& [fnType, fnInstId]: instComb) {
            curMapping[fnType] = fnInstId;
        }
        enumerate_InstanceCombination_CalcDelay<type_wgt, type_res>(stgid+1,curMapping, curBstMapping, minBstDelay, partParSFC, stg2InstCombinations, oldUtilization, cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork, showInConsoleDetailed);
    }
}//enumerate_InstanceCombination_CalcDelay

/*! For a given SFC and its all partial parallel function, find all the possibile instance mappings (brute force).
 * Intance mapping is based on min delay of the path.
 * For a partial parallel chain, it first find all its instance combinations in each stage, then out of all instance mappings it find min delay.
 * @param[in,out] cSFC given SFC for which we have to find minimum delay mapping
 * @param[in,out] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * @param PhysicalNetwork Physical Network class consists of network graph
 * @param showInConsole print in the console.
 * @param showInConsoleDetailed print of each step in the console
 */
template<typename type_wgt, typename type_res>
int partialChains_Sequential_Deployment(sfcResult& sol_partial, const vector<vector<unsigned int>>& seqSFC, ServiceFunctionChain * const cSFC, VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
                                        bool showInConsole = false, bool showInConsoleDetailed = false) {//partialChains_Sequential_Deployment
        
    const unsigned int szStages = seqSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }} */
    unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
    for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
        const auto& curStg = seqSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_stageInstanceCombination<type_res>(0, curInstComb, stgid, curStg, stg2InstCombinations, VNFNetwork->seq_utilization, cSFC, VNFNetwork);
        if(stg2InstCombinations.find(stgid) == stg2InstCombinations.end()){
            sol_partial.seq_pid = noResDueToNoStg;
            return algostopped; ///< if there is no instance combinations in current stage then we can't proceed
        }
    }//stgid<szStages finding instances possibilities of each stage.
         
    if(showInConsoleDetailed){ // showing partParSFC info
        cout<<"\n seqSFC:"; for(const auto& blks: seqSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<") ---------- - --------- - ------";
        for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
            for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
            }  cout<<" }";
        }
    }// show stages wise instances combination

    /*! finding all the mapping possibilites for the current partParSFC instance combination at each stage.*/
    unordered_map<unsigned int,unsigned int> curMapping; ///< iterating mapping variable
    enumerate_InstanceCombination_CalcDelay<type_wgt, type_res>(0, curMapping,  sol_partial.seq_fninstmap, sol_partial.seq_delay, seqSFC ,stg2InstCombinations, VNFNetwork->seq_utilization, cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);//, showInConsoleDetailed
    sol_partial.seq_pid = 0;
    for(const auto& [fn, fninst]: sol_partial.seq_fninstmap){
        VNFNetwork->seq_utilization[fn][fninst] += cSFC->trafficArrivalRate;
    }
    if(showInConsole){
        cout<<"\n PartialEnumeration: SFC["<<cSFC->index<<"]   Seq: idx["<<sol_partial.seq_pid<<"]  delay:["<<sol_partial.seq_delay<<"] :";
        for(const auto &blk: cSFC->allPartParSFC[sol_partial.seq_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout<<fnid<<char(96+sol_partial.seq_fninstmap.at(fnid))<<" ";
            }   cout<<"]";
        }
    }
    return algosuccess;
}//partialChains_Sequential_Deployment

template<typename type_wgt, typename type_res>
int partialChains_FullParallel_Deployment(sfcResult& sol_partial, const vector<vector<unsigned int>>& fullParSFC, ServiceFunctionChain * const cSFC, VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
                                        bool showInConsole = false, bool showInConsoleDetailed = false) {//partialChains_FullParallel_Deployment

    const unsigned int szStages = fullParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }} */
    unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
    for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
        const auto& curStg = fullParSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_stageInstanceCombination<type_res>(0, curInstComb, stgid, curStg, stg2InstCombinations, VNFNetwork->fullpar_utilization, cSFC, VNFNetwork);
        if(stg2InstCombinations.find(stgid) == stg2InstCombinations.end()){
            sol_partial.fullpar_pid = noResDueToNoStg;
            return algostopped; ///< if there is no instance combinations in current stage then we can't proceed
        }
    }//stgid<szStages finding instances possibilities of each stage.

    if(showInConsoleDetailed){ // showing partParSFC info
        cout<<"\n fullParSFC:"; for(const auto& blks: fullParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<") ---------- - --------- - ------";
        for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
            for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
            }  cout<<" }";
        }
    }// show stages wise instances combination

    /*! finding all the mapping possibilites for the current partParSFC instance combination at each stage.*/
    unordered_map<unsigned int,unsigned int> curMapping; ///< iterating mapping variable
    enumerate_InstanceCombination_CalcDelay<type_wgt, type_res>(0, curMapping,  sol_partial.fullpar_fninstmap, sol_partial.fullpar_delay, fullParSFC ,stg2InstCombinations, VNFNetwork->fullpar_utilization, cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);//, showInConsoleDetailed
    sol_partial.fullpar_pid = cSFC->allPartParSFC.size()-1;
    for(const auto& [fn, fninst]: sol_partial.fullpar_fninstmap){
        VNFNetwork->fullpar_utilization[fn][fninst] += cSFC->trafficArrivalRate;
    }
    if(showInConsole){
        cout<<"\n PartialEnumeration: SFC["<<cSFC->index<<"]   fullpar: idx["<<sol_partial.fullpar_pid<<"]  delay:["<<sol_partial.fullpar_delay<<"] :";
        for(const auto &blk: cSFC->allPartParSFC[sol_partial.fullpar_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout<<fnid<<char(96+sol_partial.fullpar_fninstmap.at(fnid))<<" ";
            }   cout<<"]";
        }
    }
    return algosuccess;
}//partialChains_FullParallel_Deployment

/*! For a given SFC and its all partial parallel function, find all the possibile instances (brute force).
 * Intance deployment is based on min delay of the path.
 * For a partial parallel chain, find all its instance combinations in each stage, then out of all instance mappings find min delay.
 * @param[in,out] cSFC given SFC for which we have to find minimum delay mapping
 * @param[in,out] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * @param PhysicalNetwork Physical Network class consists of network graph
 * @param showInConsole print in the console.
 * @param showInConsoleDetailed print of each step in the console
 */
template<typename type_wgt, typename type_res>
bool partialChains_PartParallel_Deployment(sfcResult& sol_partial, ServiceFunctionChain * const cSFC, VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
                                       bool showInConsole = false, bool showInConsoleDetailed = false) {//partialChains_PartParallel_Deployment
 
    type_delay bst_partpar_delay = std::numeric_limits<type_delay>::max();
    //    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}} ;
    for(int ppsidx=int( cSFC->allPartParSFC.size())-1; ppsidx>=0; --ppsidx){ /*! For each partial parallel SFC without src and dest block/stage.*/
        /*! {{1},{6,4},{5}}; Each Partial SFC is without src and dest block/stage. */
        cout<<"\r\t\t\t("<<cSFC->allPartParSFC.size()-ppsidx<<"/"<<cSFC->allPartParSFC.size()<<")"; //showing progress

        const vector<vector<unsigned int>>& partParSFC=cSFC->allPartParSFC.at(ppsidx); ///< for each of the partial parallel SFC of the givenParVNF Blocks
        const unsigned int szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

        /*! level to Instances Combinations = set of instance combination in block/stage/level index j. {stgid -> 2d{ 1d instances combinations{pair<fun, inst>}  }} */
        unordered_map<unsigned int, vector<vector<pair<unsigned int,unsigned int>>> > stg2InstCombinations;
        for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage.
            const auto& curStg = partParSFC[stgid];
            vector<pair<unsigned int,unsigned int>> curInstComb;
            construct_stageInstanceCombination<type_res>(0, curInstComb, stgid, curStg, stg2InstCombinations, VNFNetwork->ppar_utilization, cSFC, VNFNetwork);
            if(stg2InstCombinations.find(stgid) == stg2InstCombinations.end()){
                sol_partial.ppar_pid = noResDueToNoStg;
                return algostopped; ///< if there is no instance combinations in current stage then we can't proceed
            }
        }//stgid<szStages finding instances possibilities of each stage.
        if(stg2InstCombinations.size() != szStages)continue; ///< if number of stages here are not same as sfc then no need to proceed.

        if(showInConsoleDetailed){ // showing partParSFC info
            cout<<"\n partParSFC["<<ppsidx<<"]: "; for(const auto& blks: partParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<") ---------- - --------- - ------";
            for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
                cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2InstCombinations[cur_lvl].size()<<") { ";
                for(const auto& instComb: stg2InstCombinations[cur_lvl]){
                    cout<<"[";  for(const auto& givenPair: instComb){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" "; } cout<<"]";
                }  cout<<" }";
            }
        }// show stages wise instances combination

        /*! finding all the mapping possibilites for the current partParSFC instance combination at each stage.*/
        unordered_map<unsigned int,unsigned int> curMapping; ///< iterating mapping variable 
        enumerate_InstanceCombination_CalcDelay<type_wgt, type_res>(0, curMapping, sol_partial.ppar_fninstmap, bst_partpar_delay, partParSFC, stg2InstCombinations, VNFNetwork->ppar_utilization, cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork);
        if(bst_partpar_delay < sol_partial.ppar_delay ){
            sol_partial.ppar_pid = ppsidx;
            sol_partial.ppar_delay = bst_partpar_delay;
        } 

    }// for each PartParSFC ppsidx.
 
    if(sol_partial.ppar_pid == noResPar){
        return algostopped;
    }

    for(const auto& [fn, fninst]: sol_partial.ppar_fninstmap){
        VNFNetwork->ppar_utilization[fn][fninst] += cSFC->trafficArrivalRate;
    }
    if(showInConsole){
        cout << "\n PartialEnumeration: SFC[" << cSFC->index << "]   Par: idx[" << sol_partial.ppar_pid << "]  delay:[" << sol_partial.ppar_delay << "] :";
        for(const auto &blk: cSFC->allPartParSFC[sol_partial.ppar_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout << fnid << char(96+sol_partial.ppar_fninstmap.at(fnid)) << " ";
            }   cout<<"]";
        }
    }
    return algosuccess;
}//partialChains_PartParallel_Deployment



/* **************************************************************************************************************** */
/*! Recursive Function to construct all Layer Graph Nodes with their instances combination of a given stage (which consists of parallel VNFs).
 * Once a single instance combination is found, then construct lgNode with that instance combination and process that stage to precompute some of the delays, frequency of PNs .
 * It may be the case that some combination is not feasible then don't construct.
 * @param csfi function index of current stage.
 * @param curStg using to iterate all functions in the stage.
 * @param[in,out] curInstComb current combination construction in progress
 * @param[in,out] lgnids It stores all the stage wise lgNode indexes {stgid -> all lgNodes id}
 * @param[in,out] idx2lgNode It stores mapping lgNode index to lgNode structure
 * @param[in] vnfDelays It consists of the delays (processing, execution, queuing) of the vnfs in that chain.
 * @param[in] oldUtilization utilization of the vnfs till now based on previous deployment.
 * @param[in] cSFC current SFC object for which we have to find minimum delay mapping
 * @param[in] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param[in] VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * For example:  partParSFC = { {1}, {6,4}, {5} }  \n
 * stg 0 (1 function has 3 instances),    3 lgNodes  =  lgNode->instCombination (index,1d[pair<1a>])  lgNode->instCombination (index,[<1b>]) lgNode->instCombination (index,[<1c>])  \n
 * stg 1 (2 par function 2 & 3 instances),6 lgNodes  =  [<6a> <4a>], [<6a> <4b>], [6a 4c], [6b 4a], [6b 4b], [6b 4c]  \n
 * stg 2 (1 function 2 instances),        2 lgNodes  =  [5a] [5b]
*/
template<typename type_res>
void construct_layerGraphNodesIC(unsigned int csfi, const vector<unsigned int>& curStg, vector<pair<unsigned int,unsigned int>>& curInstComb,
                                 vector<unsigned int>& lgnids, unordered_map<unsigned int, lgNode>& idx2lgNode,
                                 const unordered_map<unsigned int, vnfDelaysPreComputed>& vnfDelays,  const unordered_map<unsigned int, unordered_map<unsigned int, type_delay>>& oldUtilization,
                                 const ServiceFunctionChain *const cSFC, const VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork
                                 ){//construct_layerGraphNodesIC

    if(csfi == curStg.size()){ // all functions in stage iterated. curStg.size()==numOfFunction in that stage.
        const unsigned int& lgNid = idx2lgNode.size(); ///< index to assign to new lgNode
        idx2lgNode[lgNid] = lgNode(lgNid, curInstComb); ///< create a node and mapping of index to lgNode
        lgnids.push_back(lgNid);
        lgNode& lgn = idx2lgNode[lgNid];
        /*! Processing of current lgNode: count of physical servers, max time in each server */
        type_delay utilization_sum_lgy=0, servicerate_sum_lgy=0;
        for (const auto &[d_fnType, d_fnInst]: lgn.instCombination) {
            const auto &d_vm_id = VNFNetwork->I_VNFinst2VM.at(d_fnType).at(d_fnInst); const auto &d_pn_id = VirtualNetwork->I_VM2PN.at(d_vm_id);
            const vnfDelaysPreComputed& fndelay = vnfDelays.at(d_fnType);

            lgn.cntPN[d_pn_id] += 1; //< freq of PN in current lgn
            lgn.exePN[d_pn_id] = max(lgn.exePN[d_pn_id], fndelay.prcDelay + fndelay.exeDelay + fndelay.queuingDelay.at(d_fnInst));

            if (oldUtilization.count(d_fnType) and oldUtilization.at(d_fnType).count(d_fnInst)){ ///< old utilization till now of VNF
                utilization_sum_lgy += oldUtilization.at(d_fnType).at(d_fnInst);
            } servicerate_sum_lgy += VNFNetwork->VNFNodes.at(d_fnType)->serviceRate;
        }//for d_fnType, d_fnInst
        lgn.utilization = (utilization_sum_lgy/servicerate_sum_lgy)*100;
        return;
    }// base case

    const unsigned int& fn = curStg[csfi];
    const unsigned int& totInstancs = VNFNetwork->VNFNodes.at(fn)->numInstances;
    for(unsigned int fninst=1; fninst<=totInstancs; fninst++){
        if (oldUtilization.count(fn) and oldUtilization.at(fn).count(fninst)){ ///< old utilization till now of VNF
            if(oldUtilization.at(fn).at(fninst) + cSFC->trafficArrivalRate >  VNFNetwork->VNFNodes.at(fn)->serviceRate) {
                continue; // don't take this instance if its utilisation become more than service rate of function.
            }
        }
        curInstComb.emplace_back(fn, fninst); // push current instance
        construct_layerGraphNodesIC<type_res>(csfi+1, curStg, curInstComb, lgnids, idx2lgNode,  vnfDelays, oldUtilization, cSFC, VNFNetwork, VirtualNetwork); // call function for next instance
        curInstComb.pop_back(); // pop last instance and push next instance of same function for another combinations.
    }//recursive case
}//construct_layerGraphNodesIC


/*! Given partial SFC and Layer Graph (stg2lgnids, idx2lgNode), it traverse the layer graph to find the best mapping for given sfc with minimum time and min utilization.
 * In each stage we traverse atmost 3 (mxPathsK) paths from previous stage to current stage. Then mxPathsK*(stg size=num of lg nodes in stage) pairs of new paths are generated, out of which we again select atmost mxPathsK paths for next iteration.
 * Selection of paths are based on minimum time and min utilization of the path. (utilization of the path is max of any lgNode in that path)
 * @param partParSFC given partial parallel SFC
 * @param stg2lgnids mapping of stage wise lgNode indexes {stgid -> all lgNodes id}
 * @param idx2lgNode mapping lgNode index to lgNode structure
 * @param[in,out] minBstDelay min delay found till now among all partial parallel chains of a given sfc.
 * @param[in,out] bstMapping to save best mapping overall among all partial parallel chain/and its all instances.
 * @param cSFC current SFC class.
 * @param PhysicalNetwork Physical Network class consists of network graph
 * @param showInConsole print in the console.
 * @param showInConsoleDetailed print of each step in the console
 */
template<typename type_wgt, typename type_res>
void traverse_layerGraph(const vector<vector<unsigned int>>& partParSFC,  vector<vector<unsigned int>>& stg2lgnids, unordered_map<unsigned int, lgNode>& idx2lgNode,
                          type_delay& minBstDelay, unordered_map<unsigned int,unsigned int>& curBstMapping,
                          ServiceFunctionChain *const cSFC, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {//traverse_layerGraph

    const unsigned int& szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.
    vector<unordered_set<unsigned int>> uniqLgidInStgOfKPaths(szStages); ///< unique lgNode ids in each stage which is used in k shortest path. So that we don't have to iterate all lgNode in stages again and itearte only which is necessary.
    type_delay  T_tx_init = calcD_TransmissionDelay(); ///< transmission time
    unsigned int lgDSTid = idx2lgNode.size()-1; ///< destination lgNode index
    const unsigned int mxPathsK = 5; ///< maximum number of paths to consider in any stage
/* ************* function to find min pair path ************************************************************ */
    std::function<unsigned int(unsigned int)> findPathsToConsider = [&](unsigned int pqsize)->unsigned int{
        if(pqsize <= 2) return pqsize; ///< if less than two then take both
        else if(pqsize <= 4) return 2; // if 2-4 paths then take 2
        else return mxPathsK;
//        return 16;
    };
    /*! to process the min heap. select min time and min utilization paths of all the paths present in min heap.
     * Early stopping criteria -> At any point/stage if path's mindist is more than what we found the min delay for sfc till now among all partial chains,
        then don't consider this min path. Since next min path will also be more than minBst delay then no point in checking rest of the paths.
     * @param id_curstd cur stage id which we are processing
     * @param pq priority queue structure consist of all the paths of processed previous-current stage.
     * @param lastStage if we are processing last stage, then we just need one path and assign it to SFC.
     */
    std::function<int(const unsigned int, priority_queue<pqNode>&, unsigned int)> processs_min_heap = [&idx2lgNode, &uniqLgidInStgOfKPaths, &minBstDelay, &showInConsole](const unsigned int& id_curstg, priority_queue<pqNode>&pq, unsigned int numOfPathsToConsider)->int{
        if(showInConsole){ cout<<"\n TotalPathPairs:"<<pq.size();}
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
            if(showInConsole){ cout<<"\n     p:"<<min_path.mindist<<"sec | "<<min_path.utilization<<"% ["; for(const auto& kkk: min_path.path)  cout<<kkk<<" -> ";}
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
    { ///From Dummy SRC to First Stage(index 0).
        priority_queue<pqNode> pqs; ///< to find minimum delay path in the x-y pairs of first and dummy source stage.

        for (const unsigned int &lgnIdy: stg2lgnids[id_curstg]) {
            const lgNode &lgy = idx2lgNode[lgnIdy]; ///< layer graph node y, x is dummy source
            if (showInConsoleDetailed) { cout << "\n lg:" << lgy.idx << " ["; for (const auto &givenPair: lgy.instCombination) { cout << givenPair.first << char(givenPair.second - 1 + 'a') << " "; } cout << "]";}

            type_delay mx_delay_x_y = 0;///< maximum delay of the current lgNode pair (x=dummySrc,y=current lgNode).
            /*! Calculation of packet processing time. inter duplication from src to different servers, intra duplication (within same server multiple nodes) and intra merging*/
            for (const auto &[pn_y, pn_y_parallelcnt]: lgy.cntPN) { /*! for each physical server in previous stage*/
                type_delay T_d_hdr = calcD_IntraDuplicationTime(pn_y_parallelcnt);
                type_delay T_m_hdr = calcD_IntraMergingTime(pn_y_parallelcnt);
                mx_delay_x_y = max(mx_delay_x_y, T_d_hdr + T_m_hdr + lgy.exePN.at(pn_y));
            }//curStgPN
            type_delay T_d_pkt = calcD_InterDuplicationTime( lgy.cntPN.size()); ///< inter duplication time from source to lgy
            mx_delay_x_y += T_d_pkt + T_tx_init;
            if (showInConsoleDetailed) { cout << "\n      max:" << mx_delay_x_y << " | d_pkt:" << T_d_pkt; }
            pqs.emplace(idx2lgNode[0].idx, lgy.idx, mx_delay_x_y, lgy.utilization,
                        vector<unsigned int>{idx2lgNode[0].idx}); /// constructing the dummy src to lgy path.
        }

        /*! Process paths from soruce to next stage. */
        if(processs_min_heap(id_curstg, pqs, findPathsToConsider((unsigned int)pqs.size())) <= 0){
            return; /// if min heap process stopped mid-way because early stopping crteira and not a single path taken for consideration then return. No sense of processing it further
        }
    }//From Dummy SRC to First Stage(index 0).
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
            if (showInConsoleDetailed) { cout << "\n lgy:" << lgy.idx ; cout<<" ["; for(const auto& givenPair: lgy.instCombination){ cout<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  }  cout<<"]";}

            /*  lgnIdx*********************************** */
            for (const unsigned int &lgnIdx: uniqLgidInStgOfKPaths[id_prvstg]) {/*!< Iterating layer graph nodes index (which are in Path) in previous stage(source)  */
                const lgNode &lgx = idx2lgNode[lgnIdx];///< current layer Layer Graph Node
                if (showInConsoleDetailed) { cout << "\n   lgx:" << lgx.idx; cout << " ["; for (const auto &givenPair: lgx.instCombination) {  cout << givenPair.first << char(givenPair.second - 1 + 'a') << " "; } cout << "]";}

                /* *************************************** */
                /*! Once x = lgx, y = lgy are fixed. We will find maximum time for each physical server in y to determine edge (x,y). */
                type_delay mx_delay_x_y = 0;///< maximum delay of the current instance combination of the pair (x,y).

                for(const auto &[pn_y, pn_y_parallelcnt]: lgy.cntPN){ /*! for each physical server in cur node*/

                    /*! Calculation of packet processing time. inter mergring ( different server from cur server), intra duplication (within same server multiple nodes), intra Merging (within same server multiple nodes)*/
                    unsigned int py_px_same=0; if(lgx.cntPN.find(pn_y) != lgx.cntPN.end()) py_px_same =  1;
                    unsigned int cntPrevHopDiffServer = lgx.cntPN.size() - py_px_same;
                    type_delay T_m_pkt = calcD_InterMergingTime(cntPrevHopDiffServer);

                    type_delay T_d_hdr = calcD_IntraDuplicationTime(pn_y_parallelcnt);
                    type_delay T_m_hdr = calcD_IntraMergingTime(pn_y_parallelcnt);

                    type_delay mx_pktPrc =  T_m_pkt  + T_d_hdr + T_m_hdr;///< overall total time spent in packet processing from src to dest.

                    if(showInConsoleDetailed) {
                        cout<<"\n      :"<<pn_y_parallelcnt<<"[py:"<<pn_y<<"]"<<"  prvD:"<<cntPrevHopDiffServer<<" s("<<py_px_same<<")"
                            <<"   [m_pkt:"<<T_m_pkt<<" d_hdr:"<<T_d_hdr<<" m_hdr:"<<T_m_hdr<<"]"<<"   mxServer:"<<lgy.exePN.at(pn_y);
                    }
                    /*! Calculation of inter-duplication time transmission time and  propagation time (we can duplicate the packets right before sending them. */
                    type_delay mx_interdupTxPx_for_y = 0;
                    for (const auto &[pn_x, pn_x_parallelcnt]: lgx.cntPN) { /*! for each physical server in previous node*/
                        unsigned int px_py_same = 0;
                        if (lgy.cntPN.find(pn_x) != lgy.cntPN.end()) px_py_same = 1;
                        unsigned int cntNextHopDiffServer = lgy.cntPN.size() - px_py_same;
                        type_delay T_d_pkt =  calcD_InterDuplicationTime(cntNextHopDiffServer);
                        type_delay T_tx=0, T_px=0;
                        if(pn_y != pn_x){ /// if both server are different then there is transmission and propagation delay
                            T_tx = T_tx_init; T_px =  calcD_PropagationDelay<type_wgt, type_res>(pn_x, pn_y,PhysicalNetwork);
                        }
                        mx_interdupTxPx_for_y = max(mx_interdupTxPx_for_y, T_d_pkt + T_tx + T_px);///< overall total time spent in sending packet from src to dest.
                        if (showInConsoleDetailed) {
                            cout << "\n           " << pn_x_parallelcnt << "(px:" << pn_x << ")  " << "nxtD:" << cntNextHopDiffServer << " s(" << px_py_same << ")"
                                 << "   [d_pkt:" << T_d_pkt << " tx:" << T_tx << " px:" << T_px << "]";
                        }
                    }//lgx prvStgPN
                    mx_delay_x_y = max(mx_delay_x_y, mx_interdupTxPx_for_y + mx_pktPrc + lgy.exePN.at(pn_y));

                    if(showInConsoleDetailed) {
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

/*    From Last Stage to Dummy DST **** ******************************************************* */
    { //From Last Stage to Dummy DST
        id_curstg = szStages - 1;
        priority_queue<pqNode> pqd; ///< to find minimum delay path in the x-y pairs of last stage and dummy dst stage.
        for (const unsigned int &lgnIdx: uniqLgidInStgOfKPaths[id_curstg]) {
            lgNode &lgx = idx2lgNode[lgnIdx];
            if (showInConsoleDetailed) {
                cout << "\n lg:" << lgx.idx << " ["; for (const auto &givenPair: lgx.instCombination) {  cout << givenPair.first << char(givenPair.second - 1 + 'a') << " ";  } cout << "]";
            }

            type_delay T_m_pkt = calcD_InterMergingTime(lgx.cntPN.size());
            type_delay mx_delay_x_y = T_tx_init + T_m_pkt;

            for (const pqNode &kpq: lgx.kpaths) {
                pqd.emplace(lgnIdx, lgDSTid, (kpq.mindist + mx_delay_x_y), max(kpq.utilization, lgx.utilization), kpq.path);
            }
        }
        if(processs_min_heap(id_curstg, pqd, 1) <= 0){
            return; /// if min heap process stopped mid-way because early stopping crteira and not a single path taken for consideration then return. No sense of processing it further
        }
    }//From Last Stage to Dummy DST

    if(idx2lgNode[lgDSTid].kpaths[0].mindist < minBstDelay){ ///< 0th index is shortest one
        minBstDelay =  idx2lgNode[lgDSTid].kpaths[0].mindist;
        for(const auto lgid: idx2lgNode[lgDSTid].kpaths[0].path){
            if(idx2lgNode[lgid].instCombination.empty())continue;
            for(const auto &[fn, fninst]: idx2lgNode[lgid].instCombination){
                curBstMapping[fn]=fninst;
            }
        }
    }
    else return;

    if(showInConsoleDetailed){
        cout << "\n-- Layer Graph Nodes ::";
        cout<<"\nId\t"<<"kp\t"<<"Utlz\t\t"<<"maxExePN\t"<<"child";
        cout<<"\n-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----\t"<<"-----";
        for(int i=0; i<idx2lgNode.size(); i++){ const lgNode& sec= idx2lgNode[i];
            cout<<"\n"<<sec.idx<<"  | "<<sec.kpaths.size()<<" | " <<sec.utilization<<" |\t [";
            for(const auto& val: sec.exePN){ cout<<"p"<<val.first<<":"<<val.second<<"s ";}cout<<"]\t[";
            for(const auto& val: sec.children){ cout<<"id"<<val.first<<": "<<val.second<<"s; "; } cout<<"] ";;
        }
    }///show lgNodes

    if(showInConsole){ /// showing instance and its delay
        cout<<"    "<<"["; for(const auto &blk: partParSFC){ for(const auto& fnid: blk){ cout<<fnid<<char(96+curBstMapping.at(fnid))<<" ";  }  } cout<<"]";
        cout<<"["<<minBstDelay<<"sec]";
    }

}//traverse_layerGraph

/*! FINDING THE VNF DEPLOYMENT by Layer Graph Construction IF PARALLELISM is DISABLED IN SFC\n
 * For a given partial sfc which is same as sequential sfc, construct the layer graph and find the final instance mapping/deployment based on previous utilization.\n
 * Intance mapping is based on min delay of the path with min utilization (lightly loaded).\n
 * For a partial parallel chain, and its instances combination in stages, it enumerates only k shortest paths in each stage from previous to current stage.
 * @param[out] sol_layerg solution/deployement obtained by algorithm
 * @param[in] seqSFC sequential chain, this partial chain is same as given sfc.
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * @param[in] vnfDelays pre-computed vnfs delay for the given sfc
 * @param[in] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param[in] VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * @param[in] PhysicalNetwork Physical Network class consists of network graph
 * @param[in] showInConsole print in the console.
 * @param[in] showInConsoleDetailed print of each step in the console
 * @return status of the algorithm.
 */
template<typename type_wgt, typename type_res>
bool layerGraph_Sequential_Deployement(sfcResult& sol_layerg, const vector<vector<unsigned int>>& seqSFC, ServiceFunctionChain *const cSFC, const unordered_map<unsigned int, vnfDelaysPreComputed>& vnfDelays,
                                       VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork,
                                       const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {

    const unsigned int Len = seqSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    vector<vector<unsigned int>> stg2lgnids(Len); ///< stage wise lgNode indexes in order to process/traverse it.
    unordered_map<unsigned int, lgNode> idx2lgNode;///< given index it is mapped to actual layer graph node structer so that we can process it uniquely.
    idx2lgNode[0] = lgNode(0); ///< source lgNode.

    for(unsigned int stgid=0; stgid<Len; stgid++){      // finding instances possibilities of each stage and constructing lgNodes.
        const auto& curStg = seqSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_layerGraphNodesIC<type_res>(0, curStg, curInstComb, stg2lgnids[stgid], idx2lgNode, vnfDelays, VNFNetwork->seq_utilization, cSFC, VNFNetwork, VirtualNetwork);
        if(stg2lgnids[stgid].empty()){
            sol_layerg.seq_pid = noResDueToNoStg;
            return algostopped; ///< if there is no instance combinations in current stage then we can't proceed
        }
    }//stgid<szStages finding instances possibilities of each stage.

    unsigned int lgDSTid = idx2lgNode.size(); ///< destination lgNode
    idx2lgNode[lgDSTid] = lgNode(lgDSTid);

    if(showInConsoleDetailed){ // showing partParSFC info
        cout<<"\nseqSFC: ";
        for(const auto& blks: seqSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<")";
        for(int cur_lvl=0; cur_lvl<Len; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2lgnids[cur_lvl].size()<<"):";
            for(const auto& lgid: stg2lgnids[cur_lvl]){
                cout<<lgid<<"["; for(const auto &givenPair: idx2lgNode[lgid].instCombination){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  } cout<<"] ";
            }
        }
    }

    traverse_layerGraph<type_wgt, type_res>(seqSFC, stg2lgnids, idx2lgNode, sol_layerg.seq_delay, sol_layerg.seq_fninstmap,cSFC, PhysicalNetwork);
    sol_layerg.seq_pid = 0;

    /// if we found the mapping for sequential chain
    for(const auto& [fn, fninst]: sol_layerg.seq_fninstmap){
        VNFNetwork->seq_utilization[fn][fninst] += cSFC->trafficArrivalRate;
    }

    if(showInConsole) {
        cout << "\n LayerGraph: SFC[" << cSFC->index << "]"<<"\t Sequential partIdx["<<sol_layerg.seq_pid<<"]  delay:["<<sol_layerg.seq_delay<<"]:";
        for(const auto &blk: cSFC->allPartParSFC[sol_layerg.seq_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout<<fnid<<char(96+sol_layerg.seq_fninstmap.at(fnid))<<" ";
            }   cout<<"]";
        }
    }
    return algosuccess;
}


/*! FINDING THE VNF DEPLOYMENT by Layer Graph Construction IF Full PARALLELISM Enabled IN SFC\n
 * For a given partial sfc and its full parallel vnf, construct the layer graph and find the final instance mapping/deployment based on previous utilization.\n
 * Intance mapping is based on min delay of the path with min utilization (lightly loaded).\n
 * For a partial parallel chain, and its instances combination in stages, it enumerates only k shortest paths in each stage from previous to current stage.
 * @param[out] sol_layerg solution/deployement obtained by algorithm
 * @param[in] fullSFC sequential chain, this partial chain is same as given sfc.
 * @param[in] cSFC  given SFC object for which we have to find minimum delay mapping
 * @param[in] vnfDelays pre-computed vnfs delay for the given sfc
 * @param[in] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param[in] VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * @param[in] PhysicalNetwork Physical Network class consists of network graph
 * @param[in] showInConsole print in the console.
 * @param[in] showInConsoleDetailed print of each step in the console
 * @return status of the algorithm.
 */
template<typename type_wgt, typename type_res>
bool layerGraph_FullParallel_Deployment(sfcResult& sol_layerg, const vector<vector<unsigned int>>& fullSFC, ServiceFunctionChain *const cSFC, const unordered_map<unsigned int, vnfDelaysPreComputed>& vnfDelays,
                                       VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork,
                                       const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {

    const unsigned int Len = fullSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

    vector<vector<unsigned int>> stg2lgnids(Len); ///< stage wise lgNode indexes in order to process/traverse it.
    unordered_map<unsigned int, lgNode> idx2lgNode;///< given index it is mapped to actual layer graph node structer so that we can process it uniquely.
    idx2lgNode[0] = lgNode(0); ///< source lgNode.

    for(unsigned int stgid=0; stgid<Len; stgid++){      // finding instances possibilities of each stage and constructing lgNodes.
        const auto& curStg = fullSFC[stgid];
        vector<pair<unsigned int,unsigned int>> curInstComb;
        construct_layerGraphNodesIC<type_res>(0, curStg, curInstComb, stg2lgnids[stgid], idx2lgNode, vnfDelays, VNFNetwork->fullpar_utilization, cSFC, VNFNetwork, VirtualNetwork);
        if(stg2lgnids[stgid].empty()){
            sol_layerg.fullpar_pid = noResDueToNoStg;
            return algostopped; ///< if there is no instance combinations in current stage then we can't proceed
        }
    }//stgid<szStages finding instances possibilities of each stage.

    unsigned int lgDSTid = idx2lgNode.size(); ///< destination lgNode
    idx2lgNode[lgDSTid] = lgNode(lgDSTid);

    if(showInConsoleDetailed){ // showing partParSFC info
        cout<<"\nfullSFC: ";
        for(const auto& blks: fullSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<")";
        for(int cur_lvl=0; cur_lvl<Len; cur_lvl++){  // showing stage wise combination
            cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2lgnids[cur_lvl].size()<<"):";
            for(const auto& lgid: stg2lgnids[cur_lvl]){
                cout<<lgid<<"["; for(const auto &givenPair: idx2lgNode[lgid].instCombination){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  } cout<<"] ";
            }
        }
    }

    traverse_layerGraph<type_wgt, type_res>(fullSFC, stg2lgnids, idx2lgNode, sol_layerg.fullpar_delay, sol_layerg.fullpar_fninstmap,cSFC, PhysicalNetwork, showInConsole);
    sol_layerg.fullpar_pid = cSFC->allPartParSFC.size()-1;

    /// if we found the mapping for sequential chain
    for(const auto& [fn, fninst]: sol_layerg.fullpar_fninstmap){
        VNFNetwork->fullpar_utilization[fn][fninst] += cSFC->trafficArrivalRate;
    }

    if(showInConsole) {
        cout << "\n LayerGraph: SFC[" << cSFC->index << "]"<<"\t Sequential partIdx["<<sol_layerg.fullpar_pid<<"]  delay:["<<sol_layerg.fullpar_delay<<"]:";
        for(const auto &blk: cSFC->allPartParSFC[sol_layerg.fullpar_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout<<fnid<<char(96+sol_layerg.fullpar_fninstmap.at(fnid))<<" ";
            }   cout<<"]";
        }
    }
    return algosuccess;
}

/*! FINDING THE VNF DEPLOYMENT by Layer Graph Construction IF PARALLELISM ENABLED IN SFC.\n
 * For a given SFC and its all partial parallel chains, construct the layer graph and find the final instance mapping/deployment based on previous utilization.\n
 * Intance mapping is based on min delay of the path with min utilization (lightly loaded).\n
 * For a partial parallel chain, and its instances combination in stages, it enumerates only k shortest paths in each stage from previous to current stage.
 * @param[out] sol_layerg solution/deployement obtained by algorithm
 * @param[in] cSFC  given SFC for which we have to find minimum delay mapping
 * @param[in] vnfDelays pre-computed vnfs delay for the given sfc
 * @param[in] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param[in] VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * @param[in] PhysicalNetwork Physical Network class consists of network graph
 * @param[in] showInConsole print in the console.
 * @param[in] showInConsoleDetailed print of each step in the console
 * @return status of the algorithm.
 */
template<typename type_wgt, typename type_res>
bool layerGraph_PartParallel_Deployement(sfcResult& sol_layerg, ServiceFunctionChain *const cSFC, const unordered_map<unsigned int, vnfDelaysPreComputed>& vnfDelays,
                                     VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork,
                                     const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {//layerGraph_PartParallel_Deployement

    type_delay bst_partpar_delay = std::numeric_limits<type_delay>::max();
    //    vector<vector<unsigned int>> partParSFC = {{1},{6,4},{5}}; ///Testing purpose
    for(int ppsidx=int( cSFC->allPartParSFC.size())-1; ppsidx>=0; --ppsidx){/*! For each partial parallel SFC without src and dest block/stage.*/
        cout<<"\r\t\t\t("<<cSFC->allPartParSFC.size()-ppsidx<<"/"<<cSFC->allPartParSFC.size()<<")"; //showing progress
        const vector<vector<unsigned int>>& partParSFC=cSFC->allPartParSFC[ppsidx]; ///< for each of the partial parallel SFC of the givenParVNF Blocks
        const unsigned int szStages = partParSFC.size(); ///< number of block/stage/level of the partParSFC without src and dst block/stage.

        vector<vector<unsigned int>> stg2lgnids(szStages); ///< stage wise lgNode indexes in order to process/traverse it.
        unordered_map<unsigned int, lgNode> idx2lgNode;///< given index it is mapped to actual layer graph node structer so that we can process it uniquely.
        idx2lgNode[0] = lgNode(0); ///< source lgNode.

//        unsigned int mxNumOflgNodesInStg = 0; // maximum number of lgNodes present in any stage;
        for(unsigned int stgid=0; stgid<szStages; stgid++){      // finding instances possibilities of each stage and constructing lgNodes.
            const auto& curStg = partParSFC[stgid];
            vector<pair<unsigned int,unsigned int>> curInstComb;
            construct_layerGraphNodesIC<type_res>(0, curStg, curInstComb, stg2lgnids[stgid], idx2lgNode, vnfDelays, VNFNetwork->ppar_utilization, cSFC, VNFNetwork, VirtualNetwork);
//            mxNumOflgNodesInStg=max(mxNumOflgNodesInStg, (unsigned int)(stg2lgnids[stgid].size()));
            if(stg2lgnids[stgid].empty()){
                sol_layerg.ppar_pid = noResDueToNoStg;
                return algostopped; ///< if there is no instance combinations in current stage then we can't proceed
            }
        }//stgid<szStages finding instances possibilities of each stage.
        if(stg2lgnids.size() != szStages)continue; ///< if number of stages here are not same as sfc then no need to proceed.

        unsigned int lgDSTid = idx2lgNode.size(); ///< destination lgNode
        idx2lgNode[lgDSTid] = lgNode(lgDSTid);

        if(showInConsoleDetailed){ // showing partParSFC info
            cout<<"\npartParSFC["<<ppsidx<<"]: ";
                for(const auto& blks: partParSFC){ cout<<"["; for(auto fn_id: blks){  cout<<"f"<<fn_id<<" ";  } cout<<"]"; } cout<<")";
                for(int cur_lvl=0; cur_lvl<szStages; cur_lvl++){  // showing stage wise combination
                    cout<<"\n\tSTG["<<cur_lvl<<"]("<<stg2lgnids[cur_lvl].size()<<"):";
                    for(const auto& lgid: stg2lgnids[cur_lvl]){
                        cout<<lgid<<"["; for(const auto &givenPair: idx2lgNode[lgid].instCombination){  cout<<""<<givenPair.first<<char(givenPair.second-1+'a')<<" ";  } cout<<"] ";
                    }
                }
        }

        traverse_layerGraph<type_wgt, type_res>(partParSFC, stg2lgnids, idx2lgNode, bst_partpar_delay, sol_layerg.ppar_fninstmap, cSFC, PhysicalNetwork);
        if(bst_partpar_delay < sol_layerg.ppar_delay ){
            sol_layerg.ppar_pid = ppsidx;
            sol_layerg.ppar_delay = bst_partpar_delay;
        }
    }//for each PartParSFC ppsidx.
    if(sol_layerg.ppar_pid == noResPar){
        return algostopped;
    }

    /// if we found the mapping for parallel chain
    for(const auto& [fn, fninst]: sol_layerg.ppar_fninstmap){
        VNFNetwork->ppar_utilization[fn][fninst] += cSFC->trafficArrivalRate;
    }

    if(showInConsole) {
        cout << "\n [LayerGraph]  SFC[" << cSFC->index << "]" << "\t Parallel: partIdx[" << sol_layerg.ppar_pid << "] delay:[" << sol_layerg.ppar_delay << "] :";
        for(const auto &blk: cSFC->allPartParSFC[sol_layerg.ppar_pid]) {
            cout<<"[";  for(const auto& fnid: blk){
                cout << fnid << char(96+sol_layerg.ppar_fninstmap.at(fnid)) << " ";
            }   cout<<"]";
        }
    }

    return algosuccess;
}//layerGraph_PartParallel_Deployement

/* **************************************************************************************************************** */
/*! For a given SFC and its all partial parallel chains, construct the layer graph and find the final instance mapping/deployment based on previous utilization.\n
 * Intance mapping is based on min delay of the path with min utilization (lightly loaded).\n
 * For a partial parallel chain, and its instances combination in stages, it enumerates only k shortest paths in each stage from previous to current stage.
 * partial chain with index 0 is same as sequential sfc. Therefore find deployment in case of no parallelism also.
 * @param[in,out] cSFC  given SFC for which we have to find minimum delay mapping
 * @param[in,out] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param[in] VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * @param[in] PhysicalNetwork Physical Network class consists of network graph
 * @param[in] showInConsole print in the console.
 * @param[in] showInConsoleDetailed print of each step in the console
 * @return status of the algorithm.
 */
template<typename type_wgt, typename type_res>
void algo_LayerGraph_InstanceMapping(const vector<ServiceFunctionChain*>& sortedSFCs, VirtualNetworkFunctions<type_res>*const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork,
                                     const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork, bool showInConsole = false, bool showInConsoleDetailed = false) {//algo_LayerGraph_InstanceMapping

    int sfccompleted=1, sfctotal = sortedSFCs.size(); cout<<"\n [Layer Graph]\n";
    auto deploy_st = std::chrono::steady_clock::now();
    for(ServiceFunctionChain *cSFC:sortedSFCs) {
        cout<<"\r  seq["<<sfccompleted<<"/"<<sfctotal<<"]"; sfccompleted++;
//        cout<<"\n------------------------------------------------------";
//        cout<<"\n[LayerG] Seq | sfc:"<<cSFC->index<<" :  lambda: "<<cSFC->trafficArrivalRate<<" :  vnfs: "<<cSFC->numVNF;

        sfcResult sol_layerg; ///< solution for the layer graph algorithm

        auto sfc_st = std::chrono::steady_clock::now();
        unordered_map<unsigned int, vnfDelaysPreComputed> vnfDelays;///< pre-calculated VNF delays (processing, execution and queuing delay). Before iterating all the partial par sfc it is better to calculate it for each chain as they will be same for each chain.
        for(const unsigned int& fn: cSFC->vnfSeq){ /// finding some vnf delays
            VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(fn );
            vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
            vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
            for(int fnInst=1; fnInst<=dstVNFNode->numInstances; fnInst++) { ///sequential Queuing Delay
                vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay<type_res>(dstVNFNode, fnInst, VNFNetwork->seq_utilization, cSFC);
            }
        }

        layerGraph_Sequential_Deployement<type_wgt, type_res>(sol_layerg, cSFC->allPartParSFC[0], cSFC, vnfDelays, VNFNetwork, VirtualNetwork, PhysicalNetwork);
        sol_layerg.seq_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        res_layerg.solobj[cSFC->index] = sol_layerg; // save result.
    }
    res_layerg.seq_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - deploy_st).count();

    sfccompleted=1; cout<<"\n";
    deploy_st = std::chrono::steady_clock::now();
    for(ServiceFunctionChain *cSFC:sortedSFCs) {
        cout<<"\r  fullpar["<<sfccompleted<<"/"<<sfctotal<<"]"; sfccompleted++;
//        cout<<"\n------------------------------------------------------";
//        cout<<"\n[LayerG] FullPar | sfc:"<<cSFC->index<<" :  lambda: "<<cSFC->trafficArrivalRate<<" :  vnfs: "<<cSFC->numVNF<<" : PartialChains: "<<cSFC->allPartParSFC.size();

        auto sfc_st = std::chrono::steady_clock::now();
        unordered_map<unsigned int, vnfDelaysPreComputed> vnfDelays;///< pre-calculated VNF delays (processing, execution and queuing delay). Before iterating all the partial par sfc it is better to calculate it for each chain as they will be same for each chain.
        for(const unsigned int& fn: cSFC->vnfSeq){ /// finding some vnf delays
            VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(fn );
            vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
            vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
            for(int fnInst=1; fnInst<=dstVNFNode->numInstances; fnInst++) { ///sequential Queuing Delay
                vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay<type_res>(dstVNFNode, fnInst, VNFNetwork->fullpar_utilization, cSFC);
            }
        }
        layerGraph_FullParallel_Deployment<type_wgt, type_res>(res_layerg.solobj[cSFC->index], cSFC->allPartParSFC[cSFC->allPartParSFC.size()-1], cSFC, vnfDelays, VNFNetwork, VirtualNetwork, PhysicalNetwork);
        res_layerg.solobj[cSFC->index].fullpar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
    }
    res_layerg.fullpar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - deploy_st).count();


    sfccompleted=1;
    deploy_st = std::chrono::steady_clock::now();
    for(ServiceFunctionChain *cSFC:sortedSFCs) {
        cout<<"\n  ppar["<<sfccompleted<<"/"<<sfctotal<<"]"; sfccompleted++;
//        cout<<"\n---------------------------------------------------------";
//        cout<<"\n[LayerG] PartPar | sfc:"<<cSFC->index<<" :  lambda: "<<cSFC->trafficArrivalRate<<" :  vnfs: "<<cSFC->numVNF<<" : PartialChains: "<<cSFC->allPartParSFC.size();
        auto sfc_st = std::chrono::steady_clock::now();
        unordered_map<unsigned int, vnfDelaysPreComputed> vnfDelays;///< pre-calculated VNF delays (processing, execution and queuing delay). Before iterating all the partial par sfc it is better to calculate it for each chain as they will be same for each chain.
        for(const unsigned int& fn: cSFC->vnfSeq){ /// finding some vnf delays
            VNFNode<type_res> *const dstVNFNode = VNFNetwork->VNFNodes.at(fn );
            vnfDelays[fn].prcDelay = calcD_MeanProcessingDelayVNF<type_res>(dstVNFNode);
            vnfDelays[fn].exeDelay = calcD_FunctionExecutionDelay<type_res>(dstVNFNode);
            for(int fnInst=1; fnInst<=dstVNFNode->numInstances; fnInst++) { ///Parallel Queuing Delay
                vnfDelays[fn].queuingDelay[fnInst] = calcD_QueuingDelay<type_res>(dstVNFNode, fnInst, VNFNetwork->ppar_utilization, cSFC);
            }
        }

        layerGraph_PartParallel_Deployement<type_wgt, type_res>(res_layerg.solobj[cSFC->index], cSFC, vnfDelays, VNFNetwork, VirtualNetwork, PhysicalNetwork);
        res_layerg.solobj[cSFC->index].ppar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
    }
    res_layerg.ppar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - deploy_st).count();

}//algo_LayerGraph_InstanceMapping

/*! For a given SFC and its all partial parallel function, find all the possibile instance mappings (brute force).
 * Intance mapping is based on min delay of the path.
 * For a partial parallel chain, it first find all its instance combinations in each stage, then out of all instance mappings it find min delay.
 * @param[in,out] cSFC given SFC for which we have to find minimum delay mapping
 * @param[in,out] VNFNetwork Virtual Network function class consist of VNF nodes
 * @param VirtualNetwork  Virtual Network class consists of Virtual machine nodes
 * @param PhysicalNetwork Physical Network class consists of network graph
 * @param showInConsole print in the console.
 * @param showInConsoleDetailed print of each step in the console
 */
template<typename type_wgt, typename type_res>
void algo_PartialChains_InstanceMapping(const vector<ServiceFunctionChain*>& sortedSFCs, VirtualNetworkFunctions<type_res> *const VNFNetwork, const VirtualMachines<type_res> *const VirtualNetwork, const PhysicalGraph<type_wgt, type_res> *const PhysicalNetwork,
                                        bool showInConsole = false, bool showInConsoleDetailed = false) {//algo_PartialChains_InstanceMapping

    int sfccompleted=1, sfctotal = sortedSFCs.size(); cout<<"\n[Partial Enumeration]\n";
    auto deploy_st = std::chrono::steady_clock::now();
    for(ServiceFunctionChain *cSFC:sortedSFCs) {
        cout<<"\r  seq["<<sfccompleted<<"/"<<sfctotal<<"]"; sfccompleted++;
//        cout<<"\n------------------------------------------------------";
//        cout<<"\n[PartialCh] Seq | sfc:"<<cSFC->index<<" :  lambda: "<<cSFC->trafficArrivalRate<<" :  vnfs: "<<cSFC->numVNF;

        sfcResult sol_partial; // using this algo what is optimal result for a single sfc.

        auto sfc_st = std::chrono::steady_clock::now();
        partialChains_Sequential_Deployment<type_wgt, type_res>(sol_partial, cSFC->allPartParSFC[0], cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork, showInConsole, showInConsoleDetailed);
        sol_partial.seq_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        res_partial.solobj[cSFC->index] = sol_partial; // save result
    }
    res_partial.seq_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - deploy_st).count();


    sfccompleted=1; cout<<"\n";
    deploy_st = std::chrono::steady_clock::now();
    for(ServiceFunctionChain *cSFC:sortedSFCs) {
        cout<<"\r  fullpar["<<sfccompleted<<"/"<<sfctotal<<"]"; sfccompleted++;
//        cout<<"\n------------------------------------------------------";
//        cout<<"\n[PartialCh] FullPar | sfc:"<<cSFC->index<<" :  lambda: "<<cSFC->trafficArrivalRate<<" :  vnfs: "<<cSFC->numVNF<<" : PartialChains: "<<cSFC->allPartParSFC.size();

        auto sfc_st = std::chrono::steady_clock::now();
        partialChains_FullParallel_Deployment<type_wgt, type_res>(res_partial.solobj[cSFC->index], cSFC->allPartParSFC[cSFC->allPartParSFC.size()-1], cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork, showInConsole, showInConsoleDetailed);
        res_partial.solobj[cSFC->index].fullpar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
    }
    res_partial.fullpar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - deploy_st).count();


    sfccompleted=1;
    deploy_st = std::chrono::steady_clock::now();
    for(ServiceFunctionChain *cSFC:sortedSFCs) {
        cout<<"\n  ppar["<<sfccompleted<<"/"<<sfctotal<<"]"; sfccompleted++;
//        cout<<"\n------------------------------------------------------";
//        cout<<"\n[PartialCh] PartPar | sfc:"<<cSFC->index<<" :  lambda: "<<cSFC->trafficArrivalRate<<" :  vnfs: "<<cSFC->numVNF<<" : PartialChains: "<<cSFC->allPartParSFC.size();
        auto sfc_st = std::chrono::steady_clock::now();
        partialChains_PartParallel_Deployment<type_wgt, type_res>(res_partial.solobj[cSFC->index], cSFC, VNFNetwork, VirtualNetwork, PhysicalNetwork, showInConsole, showInConsoleDetailed);
        res_partial.solobj[cSFC->index].ppar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
    }
    res_partial.ppar_duration = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - deploy_st).count();


}//algo_PartialChains_InstanceMapping

/* **************************************************************************************************************** */

void showAlgoResults(vector<ServiceFunctionChain*>& SFCs, finalResults fr){

    cout<<"\n  "<<fr.name_sol<<" : ";
    for(const ServiceFunctionChain* const sfc: SFCs){
        cout<<"\nSFC:"<<sfc->index<<" | TrafficRate: "<<sfc->trafficArrivalRate<<" | cntVNFs: "<<sfc->numVNF<<" | PartialChains: "<<sfc->allPartParSFC.size();
        const auto& obj= fr.solobj[sfc->index];
            if(obj.seq_pid==noResDueToNoStg)cout<<"No Result Obtained for Sequential/Parallel due to no instance combination in one of the stage.";
            else{
                cout<<"\n\tSeq. id(";
                if(obj.seq_pid == noResSeq)cout<<" No Result Obtained for Sequential.)";
                else{ cout<<obj.seq_pid<<"):\tD: "<<obj.seq_delay<<"sec  \t";
                    for(const auto &blk: sfc->allPartParSFC[obj.seq_pid]) { cout<<" [";  for(const auto& fnid: blk){ cout<<"f"<<fnid<<char(96+obj.seq_fninstmap.at(fnid))<<" "; }   cout<<"]";
                    }
                }
                cout<<"\n\tFull Par. id(";
                if(obj.fullpar_pid == noResPar)cout << " No Result Obtained for Partial Parallel.)";
                else{ cout << obj.fullpar_pid << "):\tD: " << obj.fullpar_delay << "sec \t";
                    for(const auto &blk: sfc->allPartParSFC[obj.fullpar_pid]) { cout << " [";  for(const auto& fnid: blk){ cout << "f" << fnid << char(96 + obj.fullpar_fninstmap.at(fnid)) << " "; }   cout << "]";
                    }
                }
                cout<<"\n\tPartial Par. id(";
                if(obj.ppar_pid == noResPar)cout << " No Result Obtained for Partial Parallel.)";
                else{ cout << obj.ppar_pid << "):\tD: " << obj.ppar_delay << "sec \t";
                    for(const auto &blk: sfc->allPartParSFC[obj.ppar_pid]) { cout << " [";  for(const auto& fnid: blk){ cout << "f" << fnid << char(96 + obj.ppar_fninstmap.at(fnid)) << " "; }   cout << "]";
                    }
                }
            }

    }
    cout<<"\nDeployement Duration: Seq: "<<fr.seq_duration<<"ms." << " | Full Par: " << fr.fullpar_duration << "ms." << " | Partial Par: " << fr.ppar_duration << "ms.";
    cout<<"\n---------------------------------------------------------";
}

#endif //SFC_PARALLELIZATION_ALGORITHMS_H
