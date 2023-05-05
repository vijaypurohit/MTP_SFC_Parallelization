//
// Created by vijay on 05-05-2023.
//

#ifndef SFC_PARALLELIZATION_TESTS_H
#define SFC_PARALLELIZATION_TESTS_H

/*! Based on same set of SFC, vary the VNF deployment parameters */
bool SimulationTest_Basic(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const string& fileSFC, const pair<float, int>& parallelPairsOpt, const string& sfc_sort_opt){//SimulationTest_VaryingNetworkVNF_DeploymentParameter

    Simulations simtest(testname, dirname);

    readNetwork(simtest.fullDirName,fileNetwork, simtest.PhysicalNetwork);
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOpt.first,parallelPairsOpt.second); /// based on VNFs

    GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, 10, {0.3, 0.5}, {false, 0}, simtest.fullDirName, fileSFC);
    simtest.readGenericServiceFunctionsChains(fileSFC, sfc_sort_opt);

    simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
    simtest.DeploymentVNF_ScoreMethod(0.5,0.75,1,1);

    for(ServiceFunctionChain& sfc: simtest.allSFC){
        simtest.convert_SeqSFC_to_FullParallel(sfc);
        simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
        simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
        sfc.partialParallelChains = &sfc.subsetPartParSFC;
    }if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";

    try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
    catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
        simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
    }

    try {  bruteForce_InstanceMapping(simtest); }
    catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
        simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
    }
    simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
    simtest.showVNFsUtilization(2,simtest.TestsResult[name_kshortestpath]);

    simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
    simtest.showVNFsUtilization(1,simtest.TestsResult[name_bruteForce]);
    simtest.writeInFileSimulationTestResults(1, ios::out);
    return true;
}//SimulationTest_VaryingNetworkVNF_DeploymentParameter

/*! Based on same set of SFC, vary the VNF deployment parameters */
bool SimulationTest_VaryingNetworkVNF_DeploymentParameter(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const string& fileSFC, const pair<float, int>& parallelPairsOpt, const string& sfc_sort_opt){//SimulationTest_VaryingNetworkVNF_DeploymentParameter

    Simulations simtest(testname,dirname);

    readNetwork(simtest.fullDirName,fileNetwork, simtest.PhysicalNetwork);
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOpt.first,parallelPairsOpt.second); /// based on VNFs

    GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, 10, {0.3, 0.5}, {false, 0}, simtest.fullDirName, fileSFC);
    simtest.readGenericServiceFunctionsChains(fileSFC, sfc_sort_opt);

    simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC

    for(ServiceFunctionChain& sfc: simtest.allSFC){
        simtest.convert_SeqSFC_to_FullParallel(sfc);
        simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
        simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
        sfc.partialParallelChains = &sfc.subsetPartParSFC;
    }if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";

    int obs=1;
    for( const float& scale: {0.25f,0.5f,0.75f,1.0f} ){
        for(const float& alpha: {0.25f,0.5f,0.75f}){
            for(const int& pxs: {1,2}){
                for(const int& dist: {1,2}){
                    string other = to_string(scale)+"_"+to_string(alpha)+"_"+to_string(pxs)+"_"+to_string(dist);
                    cout<<"\n[parameters: "<<other<<"]------------";

                    simtest.DeploymentVNF_ScoreMethod(scale,alpha,pxs,dist);

                    simtest.TestsResult.clear();
                    try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
                    catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
                        simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
                    }

                    try {  bruteForce_InstanceMapping(simtest); }
                    catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
                        simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
                    }

                    if(obs == 1){  simtest.writeInFileSimulationTestResults(obs, ios::out, other);
                    }else{ simtest.writeInFileSimulationTestResults(obs, ios::app, other);
                    } obs++;
                    system(CLR);
                }
            }
        }
    }
    return true;
}//SimulationTest_VaryingNetworkVNF_DeploymentParameter

/*! Comparisorn of Execution Time for brute force vs Heuristic, for each sfclength from 2 to 10, Generate 20 random SFCs, Deploy VNFs and then Deploy SFC. */
bool SimulationTest_ExecutionTimeComparison(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const pair<float, int>& parallelPairsOpt, const string& sfc_sort_opt){//SimulationTest_ExecutionTimeComparison
    Simulations simtest(testname, dirname);
    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);/// reading network
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOpt.first, parallelPairsOpt.second); /// based on VNFs

    int obs=1;
    for(const auto& sfclen: {2,3,4,5,6,7,8,9,10}){
        const string& filesfc = "SFC20_L_"+to_string(sfclen)+"_ETC.txt";
        cout<<"\n[sfclen:"<<sfclen<<": "<<filesfc<<"]------------";
        GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, 20, {0.3, 0.5}, {true, sfclen}, simtest.fullDirName, filesfc);
        simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);

        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_FullParallel(sfc);
            simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
            simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
            sfc.partialParallelChains = &sfc.subsetPartParSFC;
        }
        if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";

        simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
        simtest.DeploymentVNF_ScoreMethod(1,0.5,1,1);

        simtest.TestsResult.clear();
        try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
        }

        try {  bruteForce_InstanceMapping(simtest); }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
        }

        if(obs == 1){ simtest.writeInFileSimulationTestResults(obs, ios::out, filesfc);
        }else{  simtest.writeInFileSimulationTestResults(obs, ios::app, filesfc);
        } obs++;
        system(CLR);
    }//sfclen
    return true;
}//SimulationTest_ExecutionTimeComparison

/*! Comparisorn of Execution Time for brute force vs Heuristic, On same VNF deployement (based on previous last SFC), for each sfclength from 2 to 10, take previous 20 generated SFCs and Deploy SFC.  */
bool SimulationTest_ExecutionTimeComparisonSameDeployement(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const string& fileSFC, const pair<float, int>& parallelPairsOpt, const string& sfc_sort_opt){//SimulationTest_ExecutionTimeComparison

    Simulations simtest(testname+"SameDeployement_"+sfc_sort_opt, dirname);

    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);/// reading network
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.readGenericServiceFunctionsChains(fileSFC, sfc_sort_opt);

    simtest.findRandomParallelPairs(parallelPairsOpt.first, parallelPairsOpt.second); /// based on VNFs
    simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
    simtest.DeploymentVNF_ScoreMethod(1,0.5,1,1); /// based on above sfc deploy vnfs

    int obs=1;
    for(const auto& sfclen: {2,3,4,5,6,7,8,9,10}){
        const string& filesfc = "SFC20_L_"+to_string(sfclen)+"-exetime.txt";
        cout<<"\n[sfclen:"<<sfclen<<": "<<filesfc<<"]------------";
        simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);

        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_FullParallel(sfc);
            simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
            simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
            sfc.partialParallelChains = &sfc.subsetPartParSFC;
        }
        if(debug)cout<<"\n\tSFCs converted to Full Parallel VNFs Blocks";

        simtest.TestsResult.clear();

        try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
        }

        try {  bruteForce_InstanceMapping(simtest); }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
        }

        if(obs == 1){ simtest.writeInFileSimulationTestResults(obs, ios::out, filesfc);
        }else{  simtest.writeInFileSimulationTestResults(obs, ios::app, filesfc);
        } obs++;
        system(CLR);
    }//sfclen
    return true;
}//SimulationTest_ExecutionTimeComparison

/*! Comparisorn of Execution Time for brute force vs Heuristic, for each sfclength from 2 to 10, Generate 20 random SFCs, Deploy VNFs and then Deploy SFC. */
bool SimulationTest_ExecutionTimeComparisonWithFixedLength(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const pair<float, int>& parallelPairsOpt, const string& sfc_sort_opt){//SimulationTest_ExecutionTimeComparison
    Simulations simtest(testname, dirname);
    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);/// reading network
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOpt.first, parallelPairsOpt.second); /// based on VNFs

    int sfclen = 10;
    for(int obs=1; obs<=20; obs++){
        const string& filesfc = "SFC20_L_"+ to_string(sfclen)+"_ETCFL.txt";
        cout<<"\n[obs:"<<obs<<": "<<filesfc<<"]------------";
        GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, 20, {0.3, 0.5}, {true, sfclen}, simtest.fullDirName, filesfc);
        simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);

        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_FullParallel(sfc);
            simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
            simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
            sfc.partialParallelChains = &sfc.subsetPartParSFC;
        }
        if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";

        simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
        simtest.DeploymentVNF_ScoreMethod(1,0.5,1,1);

        simtest.TestsResult.clear();
        try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
        }

        try {  bruteForce_InstanceMapping(simtest); }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
        }

        if(obs == 1){ simtest.writeInFileSimulationTestResults(obs, ios::out);
        }else{  simtest.writeInFileSimulationTestResults(obs, ios::app);
        }
        system(CLR);
    }//sfclen
    return true;
}//SimulationTest_ExecutionTimeComparison


#endif //SFC_PARALLELIZATION_TESTS_H
