//
// Created by vijay on 05-05-2023.
//

#ifndef SFC_PARALLELIZATION_TESTS_H
#define SFC_PARALLELIZATION_TESTS_H

bool setDelayParameterSettings(const int& setting_id, pair<type_delay, type_delay>& funExeTimeRange){
    switch(setting_id){
        case 0:         /*! PTF               rwtpb = 0.00099 ;    pktsz = 1040 B;  bw = 2 Mbps
                         *  t_d_pkt(1 packet) = 8.2368ms   t_tx = 4.16 ms         t_fnx = 1.1-1.5 ms
                         */
            packetBodySize = 1000; packetHeaderSize = 40; //factor_packet = 8;
            read_write_time_per_bit = 0.00099;
            bandwidthNW = 2; //factor_bandwidth = 1000;
            funExeTimeRange ={1.1,1.51};
            break;
        case 1:
            /*! PFT                  rwtpb = 0.00099 ;   pktsz = 1040 B;   bw = 7.5 Mbps
             *  t_d_pkt(1 packet) = 8.2368 ms    t_fnx = 4.0-5.1 ms      t_tx = 1.10933 ms
             */
            packetBodySize = 1000; packetHeaderSize = 40; //factor_packet = 8;
            read_write_time_per_bit = 0.00099;
            bandwidthNW = 7.5; //factor_bandwidth = 1000;
            funExeTimeRange ={4.0,5.1};
            break;
        case 2:
            /*! TPF                 rwtpb = 0.00055 ;  pktsz = 1040 B;      bw = 1 Mbps
             *  t_tx = 8.32 ms      t_d_pkt(1 packet) = 3.652 ms    t_fnx = 1.1-1.5 ms
             */
            packetBodySize = 1000; packetHeaderSize = 40; //factor_packet = 8;
            read_write_time_per_bit = 0.00055;
            bandwidthNW = 1; //factor_bandwidth = 1000;
            funExeTimeRange = {1.1,1.51};
            break;
        case 3:
            /*! TFP                   rwtpb = 7.7e-05 ;   pktsz = 1040 B;   bw = 1 Mbps
             *  t_tx = 8.32 ms    t_fnx = 4.0-5.1 ms      t_d_pkt(1 packet) = 0.64064 ms
             */
            packetBodySize = 1000; packetHeaderSize = 40; //factor_packet = 8;
            read_write_time_per_bit = 7.7e-05;
            bandwidthNW = 1; //factor_bandwidth = 1000;
            funExeTimeRange = {4.0,5.1};
            break;
        case 4:
            /*!  FTP                       pktsz = 528 B;  rwtpb = 7.7e-05 ;  bw = 1 Mbps
            *  t_fnx = 8.0-9.0 ms      t_tx = 4.224 ms    t_d_pkt(1 packet) = 0.325248 ms
            */
            packetBodySize = 500; packetHeaderSize = 28; //factor_packet = 8;
            read_write_time_per_bit = 7.7e-05;
            bandwidthNW = 1; //factor_bandwidth = 1000;
            funExeTimeRange = {8.0,9.1};
            break;
        case 5:
            /*!  FPT                       pktsz = 830 B;  rwtpb = 0.00055 ;  bw = 7.5 Mbps
             *  t_fnx = 8.0-9.0 ms      t_d_pkt(1 packet) = 3.652 ms          t_tx = 0.885333 ms
             */
            packetBodySize = 800; packetHeaderSize = 30; //factor_packet = 8;
            read_write_time_per_bit = 0.00055;
            bandwidthNW = 7.5; //factor_bandwidth = 1000;
            funExeTimeRange = {8.0,9.1};
            break;
        default:
            cerr<<"\nInvalid Delay Settings.";
            return false;
    }
    return true;
}

unsigned int setnumOfSFCsNeeded(const unsigned int& numV){
    if(numV >= 10 and numV <= 25){
        return 20;
    }else if(numV <= 35){
        return 30;
    }else if(numV <= 45){
        return 40;
    }else if(numV <= 55){
        return 50;
    }else{
        return 100;
    }

}

/*! Based on same set of SFC, vary the VNF deployment parameters */
bool SimulationTest_Basic(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const string& fileSFC, const pair<float, int>& parallelPairsOption = parallelPairsOpt, const string& sfc_sort_opt= "asc_length"){//SimulationTest_VaryingNetworkVNF_DeploymentParameter

    Simulations simtest(testname, dirname);

    readNetwork(simtest.fullDirName,fileNetwork, simtest.PhysicalNetwork);
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOption.first,parallelPairsOption.second); /// based on #VNFs

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
//    simtest.showVNFsUtilization(2,simtest.TestsResult[name_kshortestpath]);
    simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
//    simtest.showVNFsUtilization(1,simtest.TestsResult[name_bruteForce]);
    simtest.writeInFileSimulationTestResults(1, ios::out);
    return true;
}//SimulationTest_VaryingNetworkVNF_DeploymentParameter

/*! Based on samet of SFC, vary the VNF deployment parameters */
bool SimulationTest_VaryingNetworkVNF_DeploymentParameter(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const string& fileSFC, const pair<float, int>& parallelPairsOption = parallelPairsOpt, const string& sfc_sort_opt= "asc_length"){//SimulationTest_VaryingNetworkVNF_DeploymentParameter

    Simulations simtest(testname,dirname);

    readNetwork(simtest.fullDirName,fileNetwork, simtest.PhysicalNetwork);
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOption.first,parallelPairsOption.second); /// based on VNFs

    unsigned int numOfSFCsNeeded = setnumOfSFCsNeeded(simtest.PhysicalNetwork.numV);

    GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, numOfSFCsNeeded, sfcArrivalRateRange, {false, 0}, simtest.fullDirName, fileSFC);
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
bool SimulationTest_ImpactOfChainLength(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const pair<float, int>& parallelPairsOption = parallelPairsOpt, const string& sfc_sort_opt= "asc_length"){//SimulationTest_ImpactOfChainLength
    Simulations simtest(testname, dirname);
    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);/// reading network
    readVirtualNetworkFunctions(simtest.fullDirName, fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOption.first, parallelPairsOption.second); /// based on VNFs

    unsigned int numOfSFCsNeeded = setnumOfSFCsNeeded(simtest.PhysicalNetwork.numV);

    int obs=1;
    for(const auto& sfclen: {2,3,4,5,6,7,8,9,10}){
        const string& filesfc = "SFC"+to_string(numOfSFCsNeeded)+"_L_"+to_string(sfclen)+"_ETC.txt";
        cout<<"\n[sfclen:"<<sfclen<<": "<<filesfc<<"]------------";
        GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, numOfSFCsNeeded, sfcArrivalRateRange, {true, sfclen}, simtest.fullDirName, filesfc);
        simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);

        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_FullParallel(sfc);
            simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
            simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
            sfc.partialParallelChains = &sfc.subsetPartParSFC;
        } if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";

        simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
        simtest.DeploymentVNF_ScoreMethod(1,0.75,1,1);

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
}//SimulationTest_ImpactOfChainLength

/*! Comparisorn of Execution Time for brute force vs Heuristic, for each fixed sfclength, Generate 20 random SFCs, Deploy VNFs and then Deploy SFC. */
bool SimulationTest_FixedChainLength(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const pair<float, int>& parallelPairsOption = parallelPairsOpt, const string& sfc_sort_opt= "asc_length"){//SimulationTest_FixedChainLength
    Simulations simtest(testname, dirname);
    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);/// reading network
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOption.first, parallelPairsOption.second); /// based on VNFs

    unsigned int numOfSFCsNeeded =  setnumOfSFCsNeeded(simtest.PhysicalNetwork.numV);

    int sfclen = 7;
    for(int obs=1; obs<=15; obs++){
        const string& filesfc = "SFC"+ to_string(numOfSFCsNeeded)+"_L_"+ to_string(sfclen)+"_FCL.txt";
        cout<<"\n[obs:"<<obs<<": "<<filesfc<<"]------------";
        GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, numOfSFCsNeeded, sfcArrivalRateRange, {true, sfclen}, simtest.fullDirName, filesfc);
        simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);

        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_FullParallel(sfc);
            simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
            simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
            sfc.partialParallelChains = &sfc.subsetPartParSFC;
        } if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";
        simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
        simtest.DeploymentVNF_ScoreMethod(1,0.75,1,1);

        simtest.TestsResult.clear();
        try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
        }

        try {  bruteForce_InstanceMapping(simtest); }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_bruteForce]);
        }

        if(obs == 1){ simtest.writeInFileSimulationTestResults(obs, ios::out, filesfc+ " obs:"+to_string(obs));
        }else{  simtest.writeInFileSimulationTestResults(obs, ios::app, filesfc + " obs:"+to_string(obs));
        }
        system(CLR);
    }//sfclen
    return true;
}//SimulationTest_FixedChainLength

/*! Comparisorn of impact of instances per server from 1 to 7 for Heuristic */
bool SimulationTest_ImpactOfNumberOfInstancesPerServer(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const pair<float, int>& parallelPairsOption = parallelPairsOpt, const string& sfc_sort_opt= "asc_length"){//SimulationTest_ImpactOfNumberOfInstancesPerServer
    Simulations simtest(testname, dirname);
    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);/// reading network
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOption.first, parallelPairsOption.second); /// based on VNFs

    unsigned int numOfSFCsNeeded =  setnumOfSFCsNeeded(simtest.PhysicalNetwork.numV);

    int sfclen = 7;
    const string& filesfc = "SFC"+ to_string(numOfSFCsNeeded)+"_L_"+ to_string(sfclen)+"_INIPS.txt";
    GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, numOfSFCsNeeded, sfcArrivalRateRange, {true, sfclen}, simtest.fullDirName, filesfc);
    simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);

    for(ServiceFunctionChain& sfc: simtest.allSFC){
        simtest.convert_SeqSFC_to_FullParallel(sfc);
        simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
        simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
        sfc.partialParallelChains = &sfc.subsetPartParSFC;
    }if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";
    simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC

    for(int cores=1; cores<=7; cores++){
        cout<<"\n["<<filesfc<<" cores:"<<cores<<"]------------";

        simtest.PhysicalNetwork.setNodesCores({"fixed-all", cores});
        simtest.DeploymentVNF_ScoreMethod(1,0.75,1,1);

        simtest.TestsResult.clear();
        try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
        }

        if(cores == 1){ simtest.writeInFileSimulationTestResults(cores, ios::out, filesfc+" cores:"+ to_string(cores));
        }else{  simtest.writeInFileSimulationTestResults(cores, ios::app, filesfc+" cores:"+ to_string(cores));
        }
        system(CLR);
    }//
    return true;
}//SimulationTest_ImpactOfNumberOfInstancesPerServer

/*! Comparisorn of Heterogenous delay (when one time dominates the other) */
bool SimulationTest_ImpactOfHeterogeneousDelays(const string& testname, const string& dirname, const string& fileNetwork, const string& fileVNF, const pair<float, int>& parallelPairsOption = parallelPairsOpt, const string& sfc_sort_opt= "asc_length"){//SimulationTest_ImpactOfHeterogeneousDelays
    Simulations simtest(testname, dirname);
    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);/// reading network
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOption.first, parallelPairsOption.second); /// based on VNFs

    int obs = 1;
    for(const unsigned int& numOfSFCsNeeded: {10,20,30,40,50,60,70}){
        const string& filesfc = "SFC"+ to_string(numOfSFCsNeeded)+"_R_IHD.txt";
        cout<<"\n["<<filesfc<<" obs:"<<obs<<"]------------";
        GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, numOfSFCsNeeded, sfcArrivalRateRange, {false, 0}, simtest.fullDirName, filesfc);
        simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);
        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_FullParallel(sfc);
            simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
            simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
            sfc.partialParallelChains = &sfc.subsetPartParSFC;
        }if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";
        simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC

        simtest.DeploymentVNF_ScoreMethod(1,0.75,1,1);

        simtest.TestsResult.clear();
        try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
        catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
            simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
        }

        if(obs == 1){ simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::out, filesfc+" obs:"+ to_string(obs));
        }else{  simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::app, filesfc+" obs:"+ to_string(obs));
        }
        obs++;
        system(CLR);
    }//numOfSFCsNeeded

    return true;
}//SimulationTest_ImpactOfHeterogeneousDelays

/*! Comparisorn of impact of instances per server from 1 to 7 for Heuristic */
bool SimulationTest_SameChainsHeterogeneousDelaysImpactOfNumberOfInstancesPerServer(const string& dirname, const string& sfc_sort_opt = "asc_length"){//SimulationTest_ImpactOfNumberOfInstancesPerServer
    const string& networkFileName = "network_"+dirname+".txt";
    const string& fileVNF = "VNF"+to_string(maxNumVNFs)+".txt";
    Simulations simtest("SameChainsHeterogeneousDelaysImpactOfNumberOfInstancesPerServer", dirname);
    readNetwork(simtest.fullDirName,networkFileName,simtest.PhysicalNetwork);/// reading network first time, we will read again.
    GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, {0,0}, dirname+"/" , fileVNF); ///first time, we will generate exe time belo
    readVirtualNetworkFunctions(simtest.fullDirName,fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOpt.first, parallelPairsOpt.second); /// based on #VNFs

    unsigned int numOfSFCsNeeded = setnumOfSFCsNeeded(simtest.PhysicalNetwork.numV);

    int sfclen = 7;
    const string& filesfc = "SFC"+ to_string(numOfSFCsNeeded)+"_L_"+ to_string(sfclen)+"_INIPS.txt";
    GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, numOfSFCsNeeded, sfcArrivalRateRange, {true, sfclen}, simtest.fullDirName, filesfc);
    simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);

    for(ServiceFunctionChain& sfc: simtest.allSFC){
        simtest.convert_SeqSFC_to_FullParallel(sfc);
        simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
        simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
        sfc.partialParallelChains = &sfc.subsetPartParSFC;
    }if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";
    simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC

    for(const int& setting_id: {0,1,2,3,4,5}){
        simtest.sim_name = "SameChainsHeterogeneousDelaysImpactOfNumberOfInstancesPerServer_settingid_" + to_string(setting_id);

        pair<type_delay, type_delay> funExeTimeRange;
        setDelayParameterSettings(setting_id, funExeTimeRange);
        simtest.VNFNetwork.setVNFFnxFromRange(funExeTimeRange);

        for(int cores=1; cores<=7; cores++){
            cout<<"\n["<<filesfc<<" cores:"<<cores<<"]------------";

            simtest.PhysicalNetwork.setNodesCores({"fixed-all", cores});
            simtest.DeploymentVNF_ScoreMethod(1,0.75,1,1);

            simtest.TestsResult.clear();
            try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
            catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
                simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
            }

            if(cores == 1){ simtest.writeInFileSimulationTestResults(cores, ios::out, " setting_id:"+ to_string(setting_id)+" cores:"+ to_string(cores));
            }else{  simtest.writeInFileSimulationTestResults(cores, ios::app, " setting_id:"+ to_string(setting_id)+" cores:"+ to_string(cores));
            }
            system(CLR);
        }//
    }
    return true;
}//SimulationTest_ImpactOfNumberOfInstancesPerServer

bool SimulationTest_ImpactOfVariousDelays(const string& testname, const string& dirname, const string& fileNetwork, const pair<float, int>& parallelPairsOption = parallelPairsOpt, const string& sfc_sort_opt= "asc_length"){//SimulationTest_ImpactOfFunctionExeDelays
    Simulations simtest(testname, dirname);
    readNetwork(simtest.fullDirName,fileNetwork,simtest.PhysicalNetwork);        /// reading network

    const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"ivd.txt";
    GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, {0,0}, dirname+"/" , fileVNF); //baad me set kr rhe hai
    readVirtualNetworkFunctions(simtest.fullDirName, fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOption.first, parallelPairsOption.second); /// based on VNFs

    unsigned int numOfSFCsNeeded = setnumOfSFCsNeeded(simtest.PhysicalNetwork.numV);

    int sfclen = 7;
    const string& filesfc = "SFC"+ to_string(numOfSFCsNeeded)+"_L_"+ to_string(sfclen)+"_IVD.txt";
    GenerateRandomSFCs(simtest.PhysicalNetwork.numV, simtest.VNFNetwork.numVNF, numOfSFCsNeeded, sfcArrivalRateRange, {true, sfclen}, simtest.fullDirName, filesfc);
    simtest.readGenericServiceFunctionsChains(filesfc, sfc_sort_opt);
    for(ServiceFunctionChain& sfc: simtest.allSFC){
        simtest.convert_SeqSFC_to_FullParallel(sfc);
        simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
        simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
        sfc.partialParallelChains = &sfc.subsetPartParSFC;
    }if(debug)cout<<"\n\t[SFCs converted to Full Parallel VNFs Blocks]";

    simtest.calcLikelihoodOfTwoFunctions(); ///based on parallel pairs and VNFs in SFC
    simtest.DeploymentVNF_ScoreMethod(1,0.75,1,1);

    {   int obs = 1;
        pair<type_delay, type_delay> dummy;
        setDelayParameterSettings(0, dummy);
        simtest.sim_name = testname+"_Fx_PTF";
        vector<pair<type_delay, type_delay>> funExeTimeRange = {{0.1, 0.51}, {0.7, 1.01}, {3.0, 5.0}, {8.0,10.0}, {13.0,15.0}, {18.0,20.0}, {20.0,23.0}, {30.0,33.0}, {40.0,43.0}, {50.0,53.0}};
        for(const pair<type_delay, type_delay>& funExeTime: funExeTimeRange){
            simtest.VNFNetwork.setVNFFnxFromRange(funExeTime);
            cout<<"\n["<<filesfc<<" obs:"<<obs<<"]------------";

            simtest.TestsResult.clear();
            try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
            catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
                simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
            }

            if(obs == 1){ simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::out, filesfc+" range:"+ to_string(funExeTime.first)+"_"+ to_string(funExeTime.second));
            }else{  simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::app, filesfc+" range:"+ to_string(funExeTime.first)+"_"+ to_string(funExeTime.second));
            }
            obs++;
            system(CLR);
        }//funExeTime
    }//funExeTime

    {   int obs = 1;
        pair<type_delay, type_delay> funExeTimeRange ;
        setDelayParameterSettings(1, funExeTimeRange);
        simtest.VNFNetwork.setVNFFnxFromRange(funExeTimeRange);
        read_write_time_per_bit = 0.00055;
        simtest.sim_name = testname+"_Tx_PFT";
//        vector<tuple<unsigned int,unsigned int, type_delay>> txRange = {{400,25,7.5},{400,25,3.5},{500,28,3},{500,28,1.5},{800,30,2},{1000,40,2},{1000,40,1.5},{800,30,1},{1000,40,1}};
//        for(const auto& [pktbody, pkthdr, bw]: txRange){
//            packetBodySize = pktbody;
//            packetHeaderSize = pkthdr;
//            bandwidthNW = bw;
//            type_delay tx = calcD_TransmissionDelay();
        fixTx = true;
        for(const type_delay& tx: {0.5,1.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0}){
            fixTxVal = tx;
            cout<<"\n["<<filesfc<<" obs:"<<obs<<"]------------";
            simtest.TestsResult.clear();
            try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
            catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
                simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
            }

            if(obs == 1){ simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::out, filesfc+" tx:"+ to_string(tx));
            }else{  simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::app, filesfc+" tx:"+ to_string(tx));
            }
            obs++;
            system(CLR);
        }//tx
        fixTx = false;
    }//tx

    {   int obs = 1;
        pair<type_delay, type_delay> funExeTimeRange ;
        setDelayParameterSettings(4, funExeTimeRange);
        simtest.VNFNetwork.setVNFFnxFromRange(funExeTimeRange);
        simtest.sim_name = testname+"_Pkt_FTP";
        vector<tuple<unsigned int,unsigned int, type_delay>> txRange = {{400,25,0.00099},{500,28,0.00099},{800,30,0.00099},{1000,40,0.00099},{1500,40,0.00099},{2000,60,0.00099},{2500,50,0.00099},{3000,80,0.00099},{5000,140,0.00099}};
        for(const auto& [pktbody, pkthdr, rw]: txRange){
            packetBodySize = pktbody;
            packetHeaderSize = pkthdr;
            read_write_time_per_bit = rw;
            type_delay pktx = calcD_InterDuplicationTime(2);
            cout<<"\n["<<filesfc<<" obs:"<<obs<<"]------------";
            simtest.TestsResult.clear();
            try{  Heuristic_kShortestPath_InstanceMapping(simtest);  }
            catch (std::exception const &e) { std::cerr << "\ncaught: " << e.what();
                simtest.showSimulationTestResults(simtest.TestsResult[name_kshortestpath]);
            }

            if(obs == 1){ simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::out, filesfc+" pktx:"+ to_string(pktx));
            }else{  simtest.writeInFileSimulationTestResults(numOfSFCsNeeded, ios::app, filesfc+" pktx:"+ to_string(pktx));
            }
            obs++;
            system(CLR);
        }//Pkt
    }//Pkt

    return true;
}//SimulationTest_ImpactOfFunctionExeDelays


bool TestTime_SFCConvertions(){
    const string& fileVNF = "VNF"+to_string(maxNumVNFs)+".txt";
    Simulations simtest("TimeSFCConvertions", "sample0");
    GenerateRandomVNFs(maxNumVNFs, funServiceRateRange, {1,1}, simtest.fullDirName , fileVNF);
    readVirtualNetworkFunctions(simtest.fullDirName, fileVNF, simtest.VNFNetwork);
    simtest.findRandomParallelPairs(parallelPairsOpt.first, parallelPairsOpt.second); /// based on VNFs
    remove(fileVNF.c_str());

    ofstream fout;
    string filepathExt = output_directory+"TT_SFCConvertionsTime.csv";///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::out);
    fout<<"SFCLen, numPartialChains, TimeToConvertToAllPart, numSubsetPartialChains, TimeToConvertToSubset,"; ///< HEADER LINE


    unsigned int numOfSFCsNeeded = 50;
    for(int sfclen = 3; sfclen<=10; sfclen++){
        const string& filesfc = "SFC"+ to_string(numOfSFCsNeeded)+"_L_"+ to_string(sfclen)+"_IntComp.txt";
        GenerateRandomSFCs(50, maxNumVNFs, numOfSFCsNeeded, sfcArrivalRateRange, {true, sfclen}, simtest.fullDirName, filesfc);
        simtest.readGenericServiceFunctionsChains(filesfc);
        fout<<"\n"<<sfclen<<",";
        int sumChains=0;
        auto sfc_st = std::chrono::steady_clock::now();
        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_FullParallel(sfc);
            simtest.convert_fullParVNFBlk_to_AllPartialChains(sfc);
            sumChains+=sfc.allPartParSFC.size();
        }
        fout<<sumChains<<","<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count()<<",";

        sumChains=0;
        sfc_st = std::chrono::steady_clock::now();
        for(ServiceFunctionChain& sfc: simtest.allSFC){
            simtest.convert_SeqSFC_to_SubsetPartialChains(sfc);
            sumChains+=sfc.subsetPartParSFC.size();
        }
        fout<<sumChains<<","<<std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        remove(filesfc.c_str());
    }
    fout.close();
    return true;
}

bool TestTime_IntegerComposition(){
    std::cout.precision( numeric_limits<double>::digits10  );
    ofstream fout;
    string filepathExt = output_directory+"TT_IntegerCompositionTime.txt";///< path to .gv file without extention
    fout.open(filepathExt.c_str(), ios::out);

    fout<<"Length, Time(ms), size"; ///< HEADER LINE
    for(int len = 1; len<=20; len++){
        std::unordered_map<unsigned int, std::vector<std::vector<unsigned int>>> ic = {{0 , {{}} }};
        auto sfc_st = std::chrono::steady_clock::now();
        for(unsigned int ki = 1; ki<=len; ki++){
            vector<vector<unsigned int>> ans;
            for(unsigned int i=1; i<=ki; i++){
                for(auto s_dash: ic[ki - i]){
                    s_dash.push_back(i);
                    ans.push_back(s_dash);
                }
            }
            ic[ki] = std::move(ans);
        }
        double end=std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - sfc_st).count();
        fout<<"\n"<<len<<","<<end<<","<<ic[len].size();
    }
    fout.close();
    return true;
}

#endif //SFC_PARALLELIZATION_TESTS_H
