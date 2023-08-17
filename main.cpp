#include <chrono>
#include <cstring>
#include "ConvexHull.h"
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/convex_hull_2.h>
#include "random"
#include "list"
//#include "toplevel_tree.h"


//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef K::Point_2 Point_2;

using hrc = std::chrono::high_resolution_clock;
/*
struct BucketDynHull{
    std::vector<std::list<Point_2>> buckets={std::list<Point_2>()};
    std::vector<std::list<Point_2>> gbuckets={std::list<Point_2>()};
    std::vector<std::vector<Point_2>> uhulls={std::vector<Point_2>()};
    std::vector<std::vector<Point_2>> lhulls={std::vector<Point_2>()};
    std::vector<std::vector<Point_2>> ghulls={std::vector<Point_2>()};

    void Update(Point_2 p, bool deletion){
        if(deletion) gbuckets[0].emplace_back(p);
        else buckets[0].emplace_back(p);
        if(buckets[0].size() + gbuckets[0].size() == 2){
            int i=1;
            for(;;++i) {
                if (i == buckets.size()) {
                    buckets.emplace_back();
                    gbuckets.emplace_back();
                    uhulls.emplace_back();
                    lhulls.emplace_back();
                    ghulls.emplace_back();
                }
                uhulls[i].clear();
                lhulls[i].clear();
                ghulls[i].clear();
                if (buckets[i].size() == 0) {
                    for (int j = 0; j < i; ++j) {
                        buckets[i].splice(buckets[i].cbegin(), buckets[j]);
                        gbuckets[i].splice(gbuckets[i].cbegin(), gbuckets[j]);
                        lhulls[j].clear();
                        uhulls[j].clear();
                    }
                    // Clear matching tombstones
                    bool removed;
                    for (auto iter = gbuckets[i].begin(); iter != gbuckets[i].end(); iter++) {
                        for (auto inner = buckets[i].begin(); inner != buckets[i].end(); inner++) {
                            if (iter->y() == inner->y()) {
                                buckets[i].erase(inner);
                                iter = --gbuckets[i].erase(iter);
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            // Move buckets to closest power of two
            // Yes, doing bit-twiddling gives constant time computation of rounding to nearest power of two but also unreadable
            for(int j=i-1;j>0;j--){
                if(buckets[i].size() <= 2<<(j)) continue;// If you fit into previous bucket
                std::swap(buckets[i],buckets[j]);
                std::swap(gbuckets[i],gbuckets[j]);
                i = j;
            }
            // Construct hulls
            std::vector<int> temp;
            for(auto e: buckets[i]){
                temp.emplace_back(e.y());
            }
            std::sort(temp.begin(),temp.end());
            buckets[i].clear();
            for(int j=1; j<=temp.size();++j){
                buckets[i].push_back({j,temp[j-1]});
            }
            CGAL::lower_hull_points_2(buckets[i].begin(),buckets[i].end(), std::back_inserter(lhulls[i]));
            //std::reverse(lhulls[i].begin(), lhulls[i].end());
            CGAL::upper_hull_points_2(buckets[i].begin(),buckets[i].end(), std::back_inserter(uhulls[i]));
            std::reverse(uhulls[i].begin(), uhulls[i].end());
            uhulls[i].insert(uhulls[i].begin(),lhulls[i].front());
            lhulls[i].emplace_back(uhulls[i].back());
        } else {
            uhulls[0].clear();
            lhulls[0].clear();
            std::vector<int> temp;
            for(auto e: buckets[0]){
                temp.emplace_back(e.y());
            }
            std::sort(temp.begin(),temp.end());
            buckets[0].clear();
            for(int i=1; i<=temp.size();++i){
                buckets[0].push_back({i,temp[i-1]});
            }
            CGAL::lower_hull_points_2(buckets[0].begin(),buckets[0].end(), std::back_inserter(lhulls[0]));
            CGAL::upper_hull_points_2(buckets[0].begin(),buckets[0].end(), std::back_inserter(uhulls[0]));
            std::reverse(uhulls[0].begin(), uhulls[0].end());
            uhulls[0].insert(uhulls[0].begin(),lhulls[0].front());
            lhulls[0].emplace_back(uhulls[0].back());
        }
    }

    void Insert(Point_2 p){
        Update(p,false);
    }

    void Remove(Point_2 p){
        Update(p, true);
    }

    std::vector<std::pair<int,int>> UpperPoints(){
        for(auto e: uhulls){
            if(!e.empty()){
                auto res = std::vector<std::pair<int,int>>();
                for(auto p : e){
                    res.emplace_back(p.x(),p.y());
                }
                return res;
            }
        }
    }
    std::vector<std::pair<int,int>> LowerPoints(){
        for(auto e: lhulls){
            if(!e.empty()){
                auto res = std::vector<std::pair<int,int>>();
                for(auto p : e){
                    res.emplace_back(p.x(),p.y());
                }
                return res;
            }
        }
    }
    std::pair<double, double> ExtractLine(std::vector<Point_2>& hull, double x){
        int i = 0;
        int j = hull.size();
        while(i<j){
            auto m = i + (j-i)/2;
            if(hull.at(m).x() >= x) j = m;
            else i = m+1;
        }
        auto ub = hull[i];
        auto lb = hull[i-1];
        auto slope = (double)(ub.y() - lb.y())/(ub.x() - lb.x());
        auto offset = ub.y() - slope * ub.x();
        return {slope, offset};
    }

    bool Covers(double x, double y){
        for(int i=0; i<buckets.size(); i++){
            if(buckets[i].empty() || x <= lhulls[i].front().x() || x > lhulls[i].back().x()) continue;
            auto[uslope, uoffset] = ExtractLine(uhulls[i],x);
            auto[lslop, loffset] = ExtractLine(lhulls[i],x);
            if(y <= uslope*x+uoffset && y >= lslop * x + loffset) return true;
        }
    }
};

bool verify(CHTree<int>& ch, BucketDynHull& bh, std::vector<int> data, int lb, int ub) {
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
    auto bucket = bh.UpperPoints();
    for (int i = 0; i < out.size(); ++i) {
        if (dyn.at(i+1).first != out.at(i).x() || dyn.at(i+1).second != out.at(i).y()){
            return false;
        }
        if (bucket.at(i+1).first != out.at(i).x() || bucket.at(i+1).second != out.at(i).y()){
            return false;
        }
    }

    out.clear();

    CGAL::lower_hull_points_2(rank_data.begin(), rank_data.end(),std::back_inserter(out));
    //std::reverse(out.begin(), out.end());
    dyn = ch.LowerPoints();
    bucket = bh.LowerPoints();
    for (int i = 0; i < out.size(); ++i) {
        if (dyn.at(i).first != out.at(i).x() || dyn.at(i).second != out.at(i).y()){
            return false;
        }
        if (bucket.at(i).first != out.at(i).x() || bucket.at(i).second != out.at(i).y()){
            return false;
        }
    }
    return true;
}
bool verificationTest(int verify_step, bool shuffle){
    std::random_device rd;
    std::mt19937 g(rd());
    std::vector<int> data;
    auto CH = CHTree<int>();
    auto BH = BucketDynHull();
    int acc = 0;
    while(std::cin >> acc){
        data.emplace_back(acc);
    }
    if(shuffle) std::shuffle(data.begin(),data.end(),g);
    // Insert things
    acc = 0;
    for(auto iter = data.begin(); iter != data.end(); ++iter) {
        acc++;
        CH.Insert(*iter);
        BH.Insert({acc,*iter});
        if (acc && !(acc & (acc - 1))) {
            if (!verify(CH, BH, data, 0, acc)) return false;
        }
    }

    if(!CH.Covers(2,3) || !BH.Covers(2,3)){
        std::cout << "Failed cover 1" << std::endl;
    }
    if(!CH.Covers(3,4) || !BH.Covers(3,4)){
        std::cout << "Failed cover 2" << std::endl;
    }
    if(!CH.Covers(1.5,1.5) || !BH.Covers(1.5,1.5)){
        std::cout << "Failed cover on edge" << std::endl;
    }
    if(CH.Covers(4,1.5) || CH.Covers(4,4)  ){
        std::cout << "Failed to fail cover" << std::endl;
    }
    if(BH.Covers(4,1.5) || BH.Covers(4,4)  ){
        std::cout << "BH Failed to fail cover" << std::endl;
    }
    // Now remove things
    if(shuffle) std::shuffle(data.begin(),data.end(),g);

    acc = 0;
    int remaining = data.size();
    for(auto iter = data.begin(); iter != data.end(); ++iter) {
        acc++;
        remaining--;
        CH.Remove(*iter);
        BH.Remove({acc, *iter});
        if (remaining && !(remaining & (remaining - 1))) {
            if (!verify(CH, BH, data, acc, data.size())) return false;
        }
    }
    return true;
}
*/
void runtimeTest(int window_size){
    // Read data
    std::random_device rd;
    std::mt19937 g(rd());
    std::vector<int> data;
    int x = 0;
    while(std::cin >> x){
        data.emplace_back(x);
    }
    std::shuffle(data.begin(),data.end(),g);

    auto CH = CHTree<int>();

    auto t0 = hrc::now();
    auto t1 = hrc::now();
    auto acc = t1-t0;

    std::cout << "Inserting"<<std::endl;
    for(int i=0; i<data.size()/window_size;i++){
        t0 = hrc::now();
        for(int j=0; j<window_size; j++){
            CH.Insert(data[i*window_size+j]);
        }
        t1 = hrc::now();
        acc += (t1-t0);
        std::cout << (i+1)*window_size << " " << (t1 - t0).count() * 1e-9 << " "<< acc.count() * 1e-9 << " 0 0" <<std::endl;
    }
/*
    std::cout << "Deleting ROCH"<<std::endl;
    for(int i=0; i<data.size()/window_size;i++){
        t0 = hrc::now();
        for(int j=0; j<window_size; j++){
            CH.Remove(data[i*window_size+j]);
        }
        t1 = hrc::now();
        std::cout << (i+1)*window_size << ", " << (t1 - t0).count() * 1e-9 << std::endl;
    }
*/
}

template <typename T>
struct as {};
template <typename T>
std::vector<T> generate_data(size_t n) {
    std::vector<T> data(1<<n);
    std::mt19937 engine(42);

    std::uniform_int_distribution<int> d;
    std::generate(data.begin(),data.end(),[&](){return d(engine);});

    return data;
}
int main(int argc, char* argv[]){
    if(argc < 3){
        std::cout << "You must supply test-identifier and window size as arguments"<<std::endl;
        return -1;
    }
    if(!std::strcmp(argv[1],"VER")){
        std::cout << "Running verification against CGAL convex hull: ";
        //if(verificationTest(atoi(argv[2]),true)) std::cout << "SUCCESS";
        //else std::cout << "FAILURE";
        std::cout << std::endl;
    }else if(!std::strcmp(argv[1],"GEN")){
        auto data = generate_data<int>(atoi(argv[2]));
        for(auto e: data){
            std::cout << e << std::endl;
        }
    }else if(!std::strcmp(argv[1],"RUN")){
        runtimeTest(atoi(argv[2]));
    }
    // TODO: Compare outputs to bucket
    // TODO: Measure time
}