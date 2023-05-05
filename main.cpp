/*!
 * @author Vijay Purohit
 * @date Created by vijay on 22-02-2023.
 * @name MTP Thesis - Delay Aware Optimal Parallelization of Service Function Chains.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <random> //to generate random number
#include <functional>  // lambda function

using namespace std;
using type_delay = double; ///< data type of delay float /double
using type_wgt = double; ///< data type of distance

int debug = 1;
/**************** Custom Header Files Declaration ****************/
#include "Variables.h"
#include "PhysicalGraph.h" ///< graph structure for physical network, physical node and edge and node capacity
#include "VirtualNetworkFunctions.h" ///< structure for VNFnode, VNFs,
#include "ServiceFunctionChain.h" /// structure for SFC
#include "FileFunctions.h" ///< File reading writing functions
#include "HelperStructures.h" ///< Some of the helper structers
#include "SimulationClass.h" ///< for simulation purpose
#include "DelayCalculationFunctions.h" ///< Various Delay Calculation functions
#include "Algorithms.h" ///< algorithms implemented
#include "Tests.h" ///< Simulation Test implemented


int main()
{
    /// Run this if SFC Len is going to be more than 10
    if(maxSFCLength > 10){
        integerCompositionsEnumeration(int(integerCompositions.size()),maxSFCLength);
        find_all_nCk(int(nCk.size()), maxSFCLength);
    }

    std::cout<<"\n~~~~~~~~~~~~~~~~~~~ [ Delay-Aware Optimal Parallelization of SFC ] ~~~~~~~~~~~~~~~~~~~~~~~~";

    std::cout<<"\nInput Directories: ";
    for(const auto & di : std::filesystem::directory_iterator(input_directory)){
        if(di.is_directory()) cout<<"\n\t"<<di.path().stem();
    }

    const vector<std::string> allSimulationName = {"Basic", "VaryingDeploymentParameter", "ExecutionTimeComaprison", "ExecutionTimeComparisonWithFixedLength"};
    std::cout<<"\n\n id | SimulationName";
    for(int id=0; id<allSimulationName.size(); id++){
        cout<<"\n  "<<id<<": "<<allSimulationName[id];
    }


    string dirname;  cout<<"\n\nEnter Directory Name (network file should be present as \"network_{DirectoryName}.txt\"): "; cin>>dirname;
    int asid; cout<<"\nEnter Simulation Id: "; cin>>asid;


    if(asid == 0)
    {
        ifRejectSFC = false;
        packetBodySize = 1000; packetHeaderSize = 24; factor_packet = 8;
        bandwidthNW = 1; factor_bandwidth = 1000;

        speedOfLight = 300000; velocityFactor = 1.0;
        read_write_time_per_bit = 0.077e-2;

        const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"_B.txt";
        GenerateRandomVNFs(maxNumVNFs, {1.5,3.2}, {2.0,4.1}, dirname+"/" , fileVNF);
        SimulationTest_Basic(allSimulationName[asid], dirname, "network_"+dirname+".txt", fileVNF, "SFC10_L_R_B.txt", {60,0},"asc_length");
    }

    if(asid == 1)
    {
        ifRejectSFC = false;
        packetBodySize = 1000; packetHeaderSize = 24; factor_packet = 8;
        bandwidthNW = 1; factor_bandwidth = 1000;

        speedOfLight = 300000; velocityFactor = 1.0;
        read_write_time_per_bit = 0.077e-2;

        const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"_VDP.txt";
        GenerateRandomVNFs(maxNumVNFs, {1.5,3.2}, {2.0,4.1}, dirname+"/" , fileVNF);
//        SimulationTest_VaryingNetworkVNF_DeploymentParameter(allSimulationName[asid], dirname, "network_"+dirname+".txt", fileVNF, "SFC10_L_R_VDP.txt", {60,2},"asc_length");
    }

    if(asid == 2)
    {
        ifRejectSFC = false;
        packetBodySize = 1000; packetHeaderSize = 24; factor_packet = 8; /// 1000+24B * 8 = bits
        read_write_time_per_bit = 0.077e-2; /// for packet processing.

        bandwidthNW = 1; factor_bandwidth = 1000; /// 1Mb/ms
        speedOfLight = 300000; velocityFactor = 1.0; ///

        const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"_ETC.txt";
        GenerateRandomVNFs(maxNumVNFs, {1.5,3.2}, {2.0,4.1}, dirname+"/" , fileVNF);
        SimulationTest_ExecutionTimeComparison(allSimulationName[asid], dirname, "network_"+dirname+".txt", fileVNF, {55, 2}, "asc_length");
//        SimulationTest_ExecutionTimeComparisonSameDeployement(allSimulationName[asid], dirname, "network-"+dirname+".txt", fileVNF, "SFC20_L_10_ETC.txt", {55, 2}, "asc_length");
    }

    if(asid == 3)
    {
        ifRejectSFC = false;
        packetBodySize = 1000; packetHeaderSize = 24; factor_packet = 8; /// 1000+24B * 8 = bits
        read_write_time_per_bit = 0.077e-2; /// for packet processing.

        bandwidthNW = 1; factor_bandwidth = 1000; /// 1Mb/ms
        speedOfLight = 300000; velocityFactor = 1.0; ///

        const string& fileVNF = "VNF"+to_string(maxNumVNFs)+"_ETCFL.txt";
        GenerateRandomVNFs(maxNumVNFs, {1.5,3.2}, {2.0,4.1}, dirname+"/" , fileVNF);
        SimulationTest_ExecutionTimeComparisonWithFixedLength(allSimulationName[asid], dirname, "network_"+dirname+".txt", fileVNF, {55, 2}, "asc_length");
    }

    cout << "\n\n~~~~~~~~~~~~~~~~~~~ [ Ending Program ] ~~~~~~~~~~~~~~~~~~~~~~~~";
    return 0;
}


