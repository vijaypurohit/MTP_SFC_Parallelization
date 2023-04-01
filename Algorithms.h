//
// Created by vijay on 01-04-2023.
//

#ifndef SFC_PARALLELIZATION_ALGORITHMS_H
#define SFC_PARALLELIZATION_ALGORITHMS_H

/*! It reads data from sequential vector and then convert into full parallel service chain by detecting whether two function can be parallelised or not. \n
 * Also convert into parallel Adj list. \n
 * In the end write the parallel sfc into file.
 * @param testDirName to save the parallel conversion to file.
 * @param SFC ServiceFunctionChain object to convert sequential into parallel.
 * @param VNFNetwork object of VirtualNetworkFunctions class required to check two vnf are parallelizable or not.
 */
template<typename type_res=unsigned int>
void convertSeqSFC_to_FullParallelSFC(const string& testDirName, ServiceFunctionChain* SFC, VirtualNetworkFunctions<type_res> *VNFNetwork){
    size_t sz = SFC->vnfSeq.size(); // total vnfs including src and dest
    SFC->vnfBlocksPar.push_back({SFCsrc}); // stage 0, pushing src node
    for(int idx=1; idx<sz; idx++){
        int prv_vnf = SFC->vnfSeq[idx-1], cur_vnf = SFC->vnfSeq[idx];  // Checking prv_vnf --> cur_vnf pairs
        // NOT PARALLEL, if prv_vnf does not exist as first vnf in pair or if prv_vnf exist it is not parallel to cur_vnf, then it is not parallel pair
        if(VNFNetwork->parallelPairs.count(prv_vnf)==0 or VNFNetwork->parallelPairs[prv_vnf].find(cur_vnf) == VNFNetwork->parallelPairs[prv_vnf].end()){
            SFC->vnfBlocksPar.push_back({cur_vnf});// push cur_vnf as separte new stage
        }else if(SFC->vnfBlocksPar.back().size() == 1) { //PARALLEL Pairs STG Size=1 and size of previous stg is just one, then we can directly push into that stg.
            SFC->vnfBlocksPar.back().push_back(cur_vnf);
        }else{ // PARALLEL Pairs and STG size > 1, from previous stg we have to check if cur_vnf is parallel to all prv_vnf in previous stage
            bool pushInLstStg = true;
            for(int lst_stg_vnf: SFC->vnfBlocksPar.back()){
                if(VNFNetwork->parallelPairs.count(lst_stg_vnf)==0 or VNFNetwork->parallelPairs[lst_stg_vnf].find(cur_vnf) == VNFNetwork->parallelPairs[lst_stg_vnf].end()){
                    pushInLstStg = false;   break;
                }
            }
            if(pushInLstStg) SFC->vnfBlocksPar.back().push_back(cur_vnf);
            else SFC->vnfBlocksPar.push_back({cur_vnf});
        }//lst stg size>1
    }// for idx

    // Inserting into Adj List.
    size_t stages = SFC->vnfBlocksPar.size(); ///< stages of SFC.
    /// for each stage. stages-1 is the last stage
    for(size_t s=0; s<=stages-2; s++ ){
        ///for each vnf in source stage
        for(int srcVNFid : SFC->vnfBlocksPar[s]){
            /// for each vnf in next stage
            for(int dstVNFid : SFC->vnfBlocksPar[s+1]){
                SFC->pAdj[srcVNFid].push_back(dstVNFid);
            }
        }

    }

    // Writing to File.
    ofstream fout;
    string filepathExt = output_directory+testDirName+filename_sfc_parallel;///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::app);
    if (!fout) {
        string errorMsg = "File "+filepathExt+ " failed to open. Function: ";
        fout.clear();
        throw runtime_error(errorMsg+ __FUNCTION__);
    }

    fout << "\nParallelised SFC:["<<SFC->name<<"],nVNF["<< SFC->numVNF<<"]::\t";
    for(size_t i=0; i<SFC->vnfBlocksPar.size(); i++){
        if(i==0 and SFC->vnfBlocksPar[i][0] == SFCsrc) fout << "(SRC -> ";
        else if(i==SFC->vnfBlocksPar.size()-1 and SFC->vnfBlocksPar[i][0] == SFCdst)  fout << " DST)";
        else{
            fout<<"[ ";
            for(int vnfid: SFC->vnfBlocksPar[i])
                fout<<"f"<<vnfid<<char(96+SFC->I_VNFType2Inst[vnfid])<<"; ";
            fout<<"] -> ";
        }
    }
    fout.close();
}

#endif //SFC_PARALLELIZATION_ALGORITHMS_H
