//
// Created by vijay on 02-05-2023.
//
#ifndef SFC_PARALLELIZATION_HELPERSTRUCTURES_H
#define SFC_PARALLELIZATION_HELPERSTRUCTURES_H
/* ********  ********   ***************  ********   ***************  ********   ***************   *************** */
/*! Bin Class for bin-packing algorithm
 * @tparam type_bin data type for the content of the bin
 */
template<class type_bin>
class Bin {
private:
    type_bin capacity_; ///< capacity of the bin
    type_bin space_left_; ///< space left in the bin
public:
//    vector<type_delay> items;
    Bin() = default;
    explicit Bin(const type_bin& givenCapacity) : capacity_(givenCapacity), space_left_(givenCapacity) {}
    type_bin get_capacity() const { return capacity_; }
    type_bin get_space_left() const { return space_left_; }
    void pack(const type_bin& item_weight) { space_left_ -= item_weight;
//        items.push_back(item_weight);
    }//pack
};

/*! The First-Fit Algorithm for bin packing. It keeps all bins open, in the order in which they were opened.
 * It attempts to place each new item into the first bin in which it fits.
 * Its approximation ratio is FF(L) <= floor(1.7OPT) and there is family of input lists L for which FF(L) matches this bound.
 * @param bin_capacity capacity of the bin which will accomadate the items
 * @param items weight of each items needed to be accomadate
 * @return all the bins with items
 */
template <class type_bin>
vector<Bin<type_bin>> first_fit(const type_bin& bin_capacity, const vector<type_bin>& items) {
    vector<Bin<type_bin>> bins; ///< bins which are required to fit the items
    for (int i = 0; i < items.size(); ++i) {
        type_bin item_weight = items[i];
        bool packed = false;
        for (int j = 0; j < bins.size(); ++j) {
            if (bins[j].get_space_left() >= item_weight) {
                bins[j].pack(item_weight);
                packed = true;
                break;
            }
        }
        if (!packed) {
            Bin<type_bin> newBin(bin_capacity);
            newBin.pack(item_weight);
            bins.push_back(newBin);
        }
    }
    return bins;
}

/*! It orders the items by descending size, then calls First-Fit. Its approximation ratio is FFD(I) = 11/9 OPT(I) + 6/9\n
 * Offline algorithm has an asymptotic approximation ratio of 1.22 ~ 11/9 <= R_A^{inf} <= 5/4 = 1.25. \n
 * <a href="https://en.wikipedia.org/wiki/Bin_packing_problem">wiki/Bin_packing_problem</a>
 * @param bin_capacity capacity of the bin which will accomadate the items
 * @param items weight of each items needed to be accomadate
 * @return all the bins with items
 */
template <class type_bin>
vector<Bin<type_bin>> first_fit_decreasing(const type_delay& bin_capacity, vector<type_delay>& items){
    sort(begin(items), end(items), greater<type_bin>());
    return first_fit<type_bin>(bin_capacity, items);
}

/*! The Best-Fit Algorithm for bin packing.
 * It keeps all bins open, but attempts to place each new item into the bin with the maximum load in which it fits.
 * Its approximation ratio is identical to that of FF,
 * that is BF(L) <= floor(1.7OPT) and there is family of input lists L for which BF(L) matches this bound.
 * @param bin_capacity capacity of the bin which will accomadate the items
 * @param items weight of each items needed to be accomadate
 * @return all the bins with items
 */
template <class type_bin>
vector<Bin<type_bin>> best_fit(const type_bin& bin_capacity, const vector<type_bin>& items) {
    vector<Bin<type_bin>> bins;
    for (int i = 0; i < items.size(); ++i) {
        type_bin item_size = items[i];
        type_bin min_space_left = bin_capacity + 1;
        int best_bin_index = -1;
        for (int j = 0; j < bins.size(); ++j) {
            if (bins[j].get_space_left() >= item_size && bins[j].get_space_left() < min_space_left) {
                best_bin_index = j;
                min_space_left = bins[j].get_space_left();
            }
        }
        if (best_bin_index != -1) {
            bins[best_bin_index].pack(item_size);
        } else {
            Bin<type_bin> newBin(bin_capacity);
            newBin.pack(item_size);
            bins.push_back(newBin);
        }
    }
    return bins;
}
/*! It orders the items by descending size, then calls Best-Fit. Its approximation ratio is BFD(I) = 11/9 OPT(I) + 6/9.\n
 * Offline algorithm has an asymptotic approximation ratio of 1.22 ~ 11/9 <= R_A^{inf} <= 5/4 = 1.25. \n
 * <a href="https://en.wikipedia.org/wiki/Bin_packing_problem">wiki/Bin_packing_problem</a>
 * @param bin_capacity capacity of the bin which will accomadate the items
 * @param items weight of each items needed to be accomadate
 * @return all the bins with items
 */
template <class type_bin>
vector<Bin<type_bin>> best_fit_decreasing(const type_bin& bin_capacity, vector<type_bin>& items){
    sort(begin(items), end(items), greater<type_bin>());
    return best_fit<type_bin>(bin_capacity, items);
}

/* ********  ********   ***************  ********   ***************  ********   ***************   *************** */
/*! @brief Before finding best mapping of SFC, this will save delays based on VNFs in SFC to pre-calculate in order to avoid multiple computations of same VNF.*/
struct vnfDelaysPreComputed{
    type_delay exeDelay{}; ///< execution delay of the VNF
    type_delay prcDelay{}; ///< processing delay of the VNF
    unordered_map<unsigned int, type_delay> queuingDelay; ///< queuing delay of the VNF and its instance
};

/*! @brief It is a comparator function used to sort the SFCs according to the length of sfc in ascending order.
 * If length of sfc are same then sort according to the acending order of the arrival rate.
 */
bool comparator_sfc_asc_length (const ServiceFunctionChain*const& sfc1, const ServiceFunctionChain*const& sfc2) {
    if(sfc1->numVNF == sfc2->numVNF){ /// sort by asc order of len of vnf
        return sfc1->trafficArrivalRate < sfc2->trafficArrivalRate;
    }
    return sfc1->numVNF < sfc2->numVNF; /// min heap,
} //
/*! @brief It is a comparator function used to sort the SFCs according  to the length of sfc in descending order.
 * If length of sfc are same then sort according to the descending order of the arrival rate.
 */
bool comparator_sfc_dsc_length (const ServiceFunctionChain*const& sfc1, const ServiceFunctionChain*const& sfc2) {
    if(sfc1->numVNF == sfc2->numVNF){ /// sort by asc order of len of vnf
        return sfc1->trafficArrivalRate > sfc2->trafficArrivalRate;
    }
    return sfc1->numVNF > sfc2->numVNF; /// min heap,
} //
/*! @brief It is a comparator function used to sort the SFCs according to the traffic arrival rate in descending order.
 * If traffic rate are same then sort according to the descending order of the length.
 */
bool comparator_sfc_dsc_rate (const ServiceFunctionChain*const& sfc1, const ServiceFunctionChain*const& sfc2){
    if(sfc1->trafficArrivalRate == sfc2->trafficArrivalRate){
        return sfc1->numVNF > sfc2->numVNF;
    }
    return sfc1->trafficArrivalRate > sfc2->trafficArrivalRate; /// < max heap, min heap >
}//
/* ********  ********   ***************  ********   ***************  ********   ***************   *************** */
/*! @brief For a single SFC, according to the algorithm what is the optimal/best (minimum delay) parameters we have found.*/
struct SFC_RESULT{
    int seq_status{-1}, fullpar_status{-1}, ppar_pid{-1} ; ///< idx of allPartParSFC Array for sequential/part parallel/full parallel for which algorithm give optimal answer.
    type_delay seq_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of sequential length chain according to our algorithm.
    type_delay ppar_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of partial parallel chain according to our algorithm.
    type_delay fullpar_delay{std::numeric_limits<type_delay>::max()}; ///< Best time of full parallel chain according to our algorithm.
    double seq_duration{0}, ppar_duration{0}, fullpar_duration{0};
    unordered_map<unsigned int, unsigned int> seq_fninstmap, ppar_fninstmap, fullpar_fninstmap; ///< Best mapping {fun->its instance taken} for given sequential length chain according to our algorithm.
    type_delay seq_load{1000}, fullpar_load{1000}, ppar_load{1000};
};

/*! @brief Final Output produced by the algorithm.*/
class SimTEST{
public:
    std::string name{}; ///< name of the solution
    unordered_map<unsigned int, SFC_RESULT> sfcsol; ///< sfc index to its solution values
    double total_seq_duration{0}, total_ppar_duration{0}, total_fullpar_duration{0}; ///< time taken to construct the solution with/without parallelism
    unordered_map<unsigned int, unordered_map<unsigned int, type_delay>> seq_utilization, ppar_utilization, fullpar_utilization; ///< utilization of the VNF_Inst. {VNFid -> {instid -> utilization}} for sequential and parallel chain
    type_delay total_seq_delay{0}, total_ppar_delay{0}, total_fullpar_delay{0};///< sum of delays
    type_delay total_seq_load{0}, total_fullpar_load{0}, total_ppar_load{0};
    explicit SimTEST(string name){
        this->name = std::move(name);
    }
    SimTEST() = default;
};

/* ********  ********   ***************  ********   ***************  ********   ***************   *************** */
/*! @brief minheap priority queue node to find minimum dist and minimum utilization path from source to destination/current stage.*/
template<class T>
bool approximatelyEqual(const T& a, const T& b, const T& epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
//        return (static_cast<int>(a*100.0) == static_cast<int>(b*100.0));
}
struct pqNode{
    unsigned int x{}, y{}; ///< source and destination pair of previous and current stage lgNode
    type_delay mindist{}; ///< minimum delay of the path from source to current processing node node
    type_delay utilization{}; ///< max utilization of the path from source to current processing node node
    vector<unsigned int> path; ///< path constructed till now from source
    pqNode()=default;
    pqNode(unsigned int givenlgSrcId, unsigned int givenlgDstId, type_delay givenDist, type_delay givenUtilization, std::vector<unsigned int> givenPath):x(givenlgSrcId),y(givenlgDstId),mindist(givenDist),utilization(givenUtilization){
        path = std::move(givenPath);
    }

    bool operator<(const struct pqNode& other) const { // overloaded operator for priority queue
        if(approximatelyEqual<type_delay>(mindist, other.mindist, 0.0005)) /// upto 2nd decimal digit equal, 38.941, 38.94
        {
            if(static_cast<int>(utilization*100.0) != static_cast<int>(other.utilization*100.0)){ /// if util to compare
                return utilization > other.utilization; //min heap, return pair of x-y with min utilization
            }
            return boolean_distribution(generator); /// random sampling
        } else {
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
#endif //SFC_PARALLELIZATION_HELPERSTRUCTURES_H
