#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

struct Option { double p; int cost; int weight; int index; };
using Pmf = std::vector<double>;

static const int N = 14;
static const double P[N][4] = {
    {0,.90,.93,.94},{0,.95,.94,.96},{0,.85,.90,.92},{0,.83,.87,.89},
    {0,.94,.93,.95},{0,.99,.98,.99},{0,.91,.92,.94},{0,.81,.90,.92},
    {0,.97,.99,.99},{0,.83,.85,.88},{0,.94,.95,.96},{0,.79,.82,.85},
    {0,.98,.99,.99},{0,.90,.92,.95}
};
static const int C[N][4] = {
    {0,1,1,2},{0,2,1,3},{0,2,3,4},{0,3,4,5},{0,2,2,3},{0,3,3,4},
    {0,4,4,5},{0,3,5,6},{0,2,3,4},{0,4,4,5},{0,3,4,5},{0,2,3,4},
    {0,2,3,4},{0,4,4,5}
};
static const int W[N][4] = {
    {0,3,4,5},{0,8,10,11},{0,7,5,6},{0,5,6,7},{0,4,3,5},{0,5,4,6},
    {0,7,8,9},{0,4,7,8},{0,8,9,10},{0,6,5,6},{0,5,6,7},{0,4,5,6},
    {0,5,5,6},{0,6,7,8}
};

struct Search {
    std::vector<std::vector<Option>> options;
    std::vector<std::vector<int>> suffix_cost;
    std::vector<int> suffix_weight;
    int threshold;
    double target;
    long long nodes=0, leaves=0;
    int best_cost=std::numeric_limits<int>::max();
    double best_rel=-1;
    std::vector<int> best_chrom;

    Pmf convolve(const Pmf& pmf, double p, int weight) const {
        Pmf next(threshold+1,0.0);
        for(int s=0;s<=threshold;s++) {
            double mass=pmf[s];
            next[s] += mass*(1-p);
            next[std::min(threshold,s+weight)] += mass*p;
        }
        return next;
    }

    double upper(const Pmf& pmf, int start) const {
        Pmf bound=pmf;
        for(int i=start;i<N;i++) {
            double p=0; int w=0;
            for(const auto& o: options[i]) { p=std::max(p,o.p); w=std::max(w,o.weight); }
            bound=convolve(bound,p,w);
        }
        return bound[threshold];
    }

    void phase1(int i,int cost,int capacity,const Pmf& pmf,std::vector<int>& chrom) {
        nodes++;
        if(cost>best_cost) return;
        if(capacity+suffix_weight[i]<threshold) return;
        int need=std::max(0,threshold-capacity);
        if(cost+suffix_cost[i][need]>best_cost) return;
        if(upper(pmf,i)+1e-12<target) return;
        if(i==N) {
            leaves++;
            if(pmf[threshold]+1e-12>=target) best_cost=std::min(best_cost,cost);
            return;
        }
        for(const auto& o: options[i]) {
            chrom.push_back(o.index);
            phase1(i+1,cost+o.cost,capacity+o.weight,convolve(pmf,o.p,o.weight),chrom);
            chrom.pop_back();
        }
    }

    void phase2(int i,int cost,int capacity,const Pmf& pmf,std::vector<int>& chrom) {
        nodes++;
        if(cost>best_cost) return;
        if(capacity+suffix_weight[i]<threshold) return;
        int need=std::max(0,threshold-capacity);
        if(cost+suffix_cost[i][need]>best_cost) return;
        if(upper(pmf,i)+1e-12<target) return;
        if(i==N) {
            leaves++;
            if(cost==best_cost && pmf[threshold]+1e-12>=target && pmf[threshold]>best_rel+1e-12) {
                best_rel=pmf[threshold]; best_chrom=chrom;
            }
            return;
        }
        for(const auto& o: options[i]) {
            chrom.push_back(o.index);
            phase2(i+1,cost+o.cost,capacity+o.weight,convolve(pmf,o.p,o.weight),chrom);
            chrom.pop_back();
        }
    }

    void prepare() {
        suffix_cost.assign(N+1,std::vector<int>(threshold+1,1000000000));
        suffix_cost[N][0]=0;
        suffix_weight.assign(N+1,0);
        for(int i=N-1;i>=0;i--) {
            suffix_weight[i]=suffix_weight[i+1]+options[i].back().weight;
            for(int need=0;need<=threshold;need++) {
                for(const auto& o: options[i]) {
                    int rest=std::max(0,need-o.weight);
                    suffix_cost[i][need]=std::min(suffix_cost[i][need],o.cost+suffix_cost[i+1][rest]);
                }
            }
        }
    }

    double evaluate(const std::vector<int>& chromosome) const {
        Pmf pmf(threshold+1,0.0); pmf[0]=1.0;
        for(int i=0;i<N;i++) {
            const Option* selected=nullptr;
            for(const auto& option: options[i]) if(option.index==chromosome[i]) { selected=&option; break; }
            pmf=convolve(pmf,selected->p,selected->weight);
        }
        return pmf[threshold];
    }

    void seed_greedy() {
        std::vector<int> chromosome(N,0);
        for(int i=0;i<N;i++) chromosome[i]=std::max_element(options[i].begin(),options[i].end(),[](const Option&a,const Option&b){return a.weight<b.weight;})->index;
        double reliability=evaluate(chromosome);
        if(reliability+1e-12<target) return;
        bool improved=true;
        while(improved) {
            improved=false;
            int current_cost=0; for(int i=0;i<N;i++) for(const auto&o:options[i]) if(o.index==chromosome[i]) current_cost+=o.cost;
            for(int i=0;i<N && !improved;i++) {
                Option current; for(const auto&o:options[i]) if(o.index==chromosome[i]) current=o;
                for(const auto& candidate:options[i]) if(candidate.cost<current.cost) {
                    auto trial=chromosome; trial[i]=candidate.index;
                    double trial_r=evaluate(trial); int trial_cost=current_cost-current.cost+candidate.cost;
                    if(trial_cost<current_cost && trial_r+1e-12>=target) { chromosome=trial; reliability=trial_r; improved=true; break; }
                }
            }
        }
        best_cost=0; for(int i=0;i<N;i++) for(const auto&o:options[i]) if(o.index==chromosome[i]) best_cost+=o.cost;
    }
};

int main(int argc,char** argv) {
    std::ofstream out(argc>1?argv[1]:"fyffe_campaign_cpp.csv");
    out << "alpha,gamma,best_cost,best_reliability,nodes,leaves\n";
    for(double alpha: {0.0,.3,.5,.7,.9,1.0}) for(double gamma: {0.0,.3,.5,.7,1.0}) {
        double k1=40,k2=50,k3=60;
        double kl=k1+alpha*(k2-k1), ku=k3-alpha*(k3-k2);
        int threshold=(int)std::llround((ku-gamma*(ku-kl))*100.0);
        Search s; s.threshold=threshold; s.target=.90; s.options.resize(N);
        for(int i=0;i<N;i++) for(int j=0;j<4;j++) {
            double lower=W[i][j]*(.9+.1*alpha), upper=W[i][j]*(1.1-.1*alpha);
            int weight=(int)std::llround((lower+gamma*(upper-lower))*100.0);
            s.options[i].push_back({P[i][j],C[i][j],weight,j});
        }
        int common=threshold;
        for(const auto& row:s.options) for(const auto& option:row) common=std::gcd(common,option.weight);
        if(common>1) {
            threshold/=common;
            for(auto& row:s.options) for(auto& option:row) option.weight/=common;
        }
        for(auto& row:s.options) std::sort(row.begin(),row.end(),[](const Option&a,const Option&b){return a.cost<b.cost;});
        s.prepare(); s.seed_greedy(); Pmf initial(threshold+1,0.0); initial[0]=1.0; std::vector<int> chrom;
        s.phase1(0,0,0,initial,chrom);
        if(s.best_cost==std::numeric_limits<int>::max()) { out<<alpha<<","<<gamma<<",-1,0,"<<s.nodes<<","<<s.leaves<<"\n"; continue; }
        s.phase2(0,0,0,initial,chrom);
        out<<std::fixed<<std::setprecision(1)<<alpha<<","<<gamma<<","<<s.best_cost<<","<<std::setprecision(12)<<s.best_rel<<","<<s.nodes<<","<<s.leaves<<"\n";
        std::cerr<<"alpha="<<alpha<<" gamma="<<gamma<<" cost="<<s.best_cost<<" R="<<s.best_rel<<" nodes="<<s.nodes<<"\n";
    }
}
