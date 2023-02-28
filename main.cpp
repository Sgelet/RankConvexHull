#include "ConvexHull.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include "random"
#include "list"
//#include "toplevel_tree.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

bool verify(CHTree<int>& ch, std::vector<int> data, int lb, int ub) {
    std::vector<Point_2> out = {};
    std::vector<int> sorted_data;
    for(int i=lb; i<ub; ++i){
        sorted_data.push_back(data[i]);
    }
    std::vector<Point_2> rank_data;
    std::sort(sorted_data.begin(),sorted_data.end());
    for(int i=0; i<sorted_data.size(); ++i){
        rank_data.push_back({i+1,sorted_data[i]});
    }
    CGAL::upper_hull_points_2(rank_data.begin(), rank_data.end(),std::back_inserter(out));
    std::reverse(out.begin(), out.end());
    auto dyn = ch.UpperPoints();
    for (int i = 0; i < out.size(); ++i) {
        if (dyn.at(i+1).first != out.at(i).x() || dyn.at(i+1).second != out.at(i).y()){
            return false;
        }
    }

    out.clear();

    CGAL::lower_hull_points_2(rank_data.begin(), rank_data.end(),std::back_inserter(out));
    //std::reverse(out.begin(), out.end());
    dyn = ch.LowerPoints();
    for (int i = 0; i < out.size(); ++i) {
        if (dyn.at(i).first != out.at(i).x() || dyn.at(i).second != out.at(i).y()){
            return false;
        }
    }
    return true;
}

struct BucketDynHull{
    std::vector<std::list<Point_2>> buckets={std::list<Point_2>()};
    std::vector<std::list<Point_2>> gbuckets={std::list<Point_2>()};
    std::vector<std::list<Point_2>> hulls={std::list<Point_2>()};
    std::vector<std::list<Point_2>> ghulls={std::list<Point_2>()};

    void Insert(Point_2 p){
        buckets[0].emplace_back(p); // TODO: Recompute hull
        if(buckets[0].size() + gbuckets[0].size() == 2<<4){
            for(int i=1;;++i ){
                if(i == buckets.size()){
                    buckets.emplace_back();
                    gbuckets.emplace_back();
                    hulls.emplace_back();
                    ghulls.emplace_back();
                }
                hulls[i].clear();
                ghulls[i].clear();
                if(buckets[i].size() == 0){
                    for(int j=0; j<i; ++j){
                        buckets[i].splice(buckets[i].cbegin(),buckets[j]);
                        gbuckets[i].splice(gbuckets[i].cbegin(),gbuckets[j]);
                    }
                    // Clear matching tombstones
                    bool removed;
                    for(auto iter = gbuckets[i].begin(); iter != gbuckets[i].end();iter++){
                        for(auto inner = buckets[i].begin(); inner != buckets[i].end(); inner++){
                            if(*iter == *inner){
                                buckets[i].erase(inner);
                                iter = --gbuckets[i].erase(iter);
                                break;
                            }
                        }
                    }


                }
                // Move buckets to closest power of two
                // Yes, doing bit-twiddling gives constant time computation of rounding to nearest power of two but also unreadable
                for(int j=i;j>0;j--){
                    if(buckets[i].size() <= 2<<(1+j)) continue;// If you fit into previous bucket
                    std::swap(buckets[i],buckets[j]);
                    std::swap(gbuckets[i],gbuckets[j]);
                    i = j;
                    break;
                }
                // Construct hulls
                CGAL::convex_hull_2(buckets[i].begin(),buckets[i].end(), std::back_inserter(hulls[i]));
                CGAL::convex_hull_2(gbuckets[i].begin(),gbuckets[i].end(), std::back_inserter(ghulls[i]));
                break;
            }
        }
    }
};

bool verificationTest(int verify_step, bool shuffle){
    std::random_device rd;
    std::mt19937 g(rd());
    std::vector<int> data;
    auto CH = CHTree<int>();
    auto BH = BucketDynHull();
    int acc = 0;
    int VERIFY_AFTER = 1;
    while(std::cin >> acc){
        data.emplace_back(acc);
    }
    if(shuffle) std::shuffle(data.begin(),data.end(),g);
    // Insert things
    acc = 0;
    for(auto iter = data.begin(); iter != data.end(); ++iter) {
        CH.Insert(*iter);
        if (++acc % VERIFY_AFTER == 0) {
            if (!verify(CH, data, 0, acc)) return false;
        }
    }
    // Now remove things
    if(shuffle) std::shuffle(data.begin(),data.end(),g);

    acc = 0;
    for(auto iter = data.begin(); iter != data.end(); ++iter) {
        CH.Remove(*iter);
        if (++acc % VERIFY_AFTER == 0) {
            if (!verify(CH, data, acc, data.size())) return false;
        }
    }
    return true;
}

template <typename T>
struct as {};
template <typename T>
std::vector<T> generate_data(size_t n) {
    std::vector<T> data(n);
    std::mt19937 engine(42);

    using RandomFunction = std::function<T()>;
    if constexpr (std::is_floating_point<T>()) {
        RandomFunction lognormal = std::bind(std::lognormal_distribution<T>(0, 0.5), engine);
        RandomFunction exponential = std::bind(std::exponential_distribution<T>(1.2), engine);
        auto rand = GENERATE_COPY(as<RandomFunction>{}, lognormal, exponential);
        std::generate(data.begin(), data.end(), rand);
    } else {
        T min = 0;
        if constexpr (std::is_signed_v<T>)
            //min = -10000;
        RandomFunction uniform_dense = std::bind(std::uniform_int_distribution<T>(min, 10000), engine);
        RandomFunction uniform_sparse = std::bind(std::uniform_int_distribution<T>(min, 10000000), engine);
        RandomFunction binomial = std::bind(std::binomial_distribution<T>(50000), engine);
        RandomFunction geometric = std::bind(std::geometric_distribution<T>(0.8), engine);
        std::generate(data.begin(), data.end(), uniform_sparse);
    }

    std::sort(data.begin(), data.end());
    return data;
}
int main(int argc, char* argv[]){
    if(argc < 3){
        std::cout << "You must supply test-identifier and window size as arguments"<<std::endl;
        return -1;
    }
    if(!std::strcmp(argv[1],"VER")){
        std::cout << "Running verification against CGAL convex hull: ";
        if(verificationTest(atoi(argv[2]),true)) std::cout << "SUCCESS";
        else std::cout << "FAILURE";
        std::cout << std::endl;
    }else if(!std::strcmp(argv[1],"GEN")){
        auto data = generate_data<int>(atoi(argv[2]));
        for(auto e: data){
            std::cout << e << std::endl;
        }
    }
    // TODO: Compare outputs to bucket
    // TODO: Measure time
}