// Generic two-phase exact solver for the finite-catalog weighted k-out-of-n:G
// design problem.  It is a faithful C++ port of `fuzzy_wkon.solver.solve_exact`
// and is used as a *faster backend*, never as a different algorithm.
//
// Division of labour: the audited Python layer performs the exact rational
// scalarisation and hands this program an already-integral problem, so no
// fuzzy arithmetic is duplicated (and therefore cannot diverge) here.
//
// stdin format
//   <N> <target> <threshold> <time_limit_seconds>
//   for each subsystem i:
//     <m_i>
//     <choice_index> <p> <cost> <weight>      (m_i lines)
//
// stdout: one JSON object.
//
// Build:  g++ -O3 -std=c++17 exact_solver.cpp -o exact_solver

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

namespace {

struct Option {
    int choice_index;
    double reliability;
    int cost;
    int weight;
};

using Pmf = std::vector<double>;
using Clock = std::chrono::steady_clock;

struct TimeLimit {};

class Search {
public:
    std::vector<std::vector<Option>> options;
    int threshold = 0;
    double target = 0.0;
    double time_limit = 0.0;

    long long nodes = 0;
    long long leaves = 0;
    long long pruned_cost = 0;
    long long pruned_capacity = 0;
    long long pruned_reliability = 0;

    int best_cost = std::numeric_limits<int>::max();
    double best_reliability = -1.0;
    std::vector<int> best_chromosome;

    void prepare() {
        const int n = static_cast<int>(options.size());
        for (auto& row : options) {
            std::sort(row.begin(), row.end(), [](const Option& a, const Option& b) {
                if (a.cost != b.cost) return a.cost < b.cost;
                if (a.reliability != b.reliability) return a.reliability > b.reliability;
                if (a.weight != b.weight) return a.weight > b.weight;
                return a.choice_index < b.choice_index;
            });
        }
        const int infinity = 1000000000;
        suffix_cost_.assign(n + 1, std::vector<int>(threshold + 1, infinity));
        suffix_cost_[n][0] = 0;
        suffix_weight_.assign(n + 1, 0);
        suffix_max_p_.assign(n + 1, 0.0);
        suffix_max_w_.assign(n + 1, 0);
        for (int i = n - 1; i >= 0; --i) {
            int max_weight = 0;
            double max_p = 0.0;
            for (const Option& option : options[i]) {
                max_weight = std::max(max_weight, option.weight);
                max_p = std::max(max_p, option.reliability);
            }
            suffix_weight_[i] = suffix_weight_[i + 1] + max_weight;
            suffix_max_p_[i] = max_p;
            suffix_max_w_[i] = max_weight;
            for (int need = 0; need <= threshold; ++need) {
                int best = infinity;
                for (const Option& option : options[i]) {
                    const int rest = std::max(0, need - option.weight);
                    best = std::min(best, option.cost + suffix_cost_[i + 1][rest]);
                }
                suffix_cost_[i][need] = best;
            }
        }
        started_ = Clock::now();
    }

    // Greedy feasible incumbent: a pruning upper bound only; it never changes
    // the certified optimum.
    void seed_incumbent() {
        const int n = static_cast<int>(options.size());
        std::vector<int> chromosome(n);
        for (int i = 0; i < n; ++i) {
            const Option* heaviest = &options[i].front();
            for (const Option& option : options[i])
                if (option.weight > heaviest->weight) heaviest = &option;
            chromosome[i] = heaviest->choice_index;
        }
        if (evaluate(chromosome) + 1e-12 < target) return;

        bool improved = true;
        while (improved) {
            improved = false;
            for (int i = 0; i < n && !improved; ++i) {
                const Option* current = find(i, chromosome[i]);
                for (const Option& candidate : options[i]) {
                    if (candidate.cost >= current->cost) continue;
                    std::vector<int> trial = chromosome;
                    trial[i] = candidate.choice_index;
                    if (evaluate(trial) + 1e-12 >= target) {
                        chromosome = trial;
                        improved = true;
                        break;
                    }
                }
            }
        }
        int cost = 0;
        for (int i = 0; i < n; ++i) cost += find(i, chromosome[i])->cost;
        best_cost = cost;
    }

    void phase(int phase_index) {
        Pmf initial(threshold + 1, 0.0);
        initial[0] = 1.0;
        std::vector<int> chromosome;
        chromosome.reserve(options.size());
        recurse(phase_index, 0, 0, 0, initial, chromosome);
    }

    double evaluate(const std::vector<int>& chromosome) const {
        Pmf pmf(threshold + 1, 0.0);
        pmf[0] = 1.0;
        for (std::size_t i = 0; i < options.size(); ++i) {
            const Option* option = find(static_cast<int>(i), chromosome[i]);
            pmf = convolve(pmf, option->reliability, option->weight);
        }
        return pmf[threshold];
    }

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - started_).count();
    }

private:
    std::vector<std::vector<int>> suffix_cost_;
    std::vector<int> suffix_weight_;
    std::vector<double> suffix_max_p_;
    std::vector<int> suffix_max_w_;
    Clock::time_point started_;

    const Option* find(int subsystem, int choice_index) const {
        for (const Option& option : options[subsystem])
            if (option.choice_index == choice_index) return &option;
        return nullptr;
    }

    Pmf convolve(const Pmf& pmf, double probability, int weight) const {
        Pmf next(threshold + 1, 0.0);
        const double complement = 1.0 - probability;
        for (int state = 0; state <= threshold; ++state) {
            const double mass = pmf[state];
            if (mass == 0.0) continue;
            next[state] += mass * complement;
            next[std::min(threshold, state + weight)] += mass * probability;
        }
        return next;
    }

    // Optimistic bound: replace every remaining subsystem by a virtual item
    // carrying the separate catalog maxima of reliability and capacity.  The
    // coupling argument in the manuscript shows this dominates every
    // admissible completion.
    double optimistic(const Pmf& pmf, int start) const {
        Pmf bound = pmf;
        for (std::size_t i = start; i < options.size(); ++i)
            bound = convolve(bound, suffix_max_p_[i], suffix_max_w_[i]);
        return bound[threshold];
    }

    void check_clock() {
        if ((nodes & 0xFFF) == 0 && time_limit > 0.0 && elapsed() > time_limit)
            throw TimeLimit{};
    }

    void recurse(int phase_index, int i, int cost, int capacity, const Pmf& pmf,
                 std::vector<int>& chromosome) {
        ++nodes;
        check_clock();

        // Equal-cost branches are retained: strict '>' is what the reviewer
        // counterexample required.
        if (cost > best_cost) { ++pruned_cost; return; }
        if (capacity + suffix_weight_[i] < threshold) { ++pruned_capacity; return; }
        const int need = std::max(0, threshold - capacity);
        if (cost + suffix_cost_[i][need] > best_cost) { ++pruned_cost; return; }
        if (optimistic(pmf, i) + 1e-12 < target) { ++pruned_reliability; return; }

        if (i == static_cast<int>(options.size())) {
            ++leaves;
            const double reliability = pmf[threshold];
            if (reliability + 1e-12 < target) return;
            if (phase_index == 1) {
                best_cost = std::min(best_cost, cost);
            } else if (cost == best_cost && reliability > best_reliability + 1e-12) {
                best_reliability = reliability;
                best_chromosome = chromosome;
            }
            return;
        }

        for (const Option& option : options[i]) {
            chromosome.push_back(option.choice_index);
            recurse(phase_index, i + 1, cost + option.cost, capacity + option.weight,
                    convolve(pmf, option.reliability, option.weight), chromosome);
            chromosome.pop_back();
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    (void)argc; (void)argv;
    std::ios::sync_with_stdio(false);

    int n = 0;
    double target = 0.0;
    long long threshold = 0;
    double time_limit = 0.0;
    if (!(std::cin >> n >> target >> threshold >> time_limit)) {
        std::cout << "{\"status\":\"error\",\"message\":\"bad header\"}\n";
        return 1;
    }

    Search search;
    search.threshold = static_cast<int>(threshold);
    search.target = target;
    search.time_limit = time_limit;
    search.options.resize(n);
    for (int i = 0; i < n; ++i) {
        int m = 0;
        std::cin >> m;
        for (int j = 0; j < m; ++j) {
            Option option{};
            std::cin >> option.choice_index >> option.reliability >> option.cost >> option.weight;
            search.options[i].push_back(option);
        }
    }

    search.prepare();

    std::ostringstream out;
    out << std::setprecision(17);
    try {
        search.seed_incumbent();
        search.phase(1);
        if (search.best_cost == std::numeric_limits<int>::max()) {
            out << "{\"status\":\"completed\",\"feasible\":false,\"nodes\":" << search.nodes
                << ",\"runtime_seconds\":" << search.elapsed() << "}";
        } else {
            search.phase(2);
            out << "{\"status\":\"completed\",\"feasible\":true,\"best_cost\":" << search.best_cost
                << ",\"best_reliability\":" << search.best_reliability << ",\"best_chromosome\":[";
            for (std::size_t i = 0; i < search.best_chromosome.size(); ++i)
                out << (i ? "," : "") << search.best_chromosome[i];
            out << "],\"nodes\":" << search.nodes << ",\"leaves\":" << search.leaves
                << ",\"pruned_cost\":" << search.pruned_cost
                << ",\"pruned_capacity\":" << search.pruned_capacity
                << ",\"pruned_reliability\":" << search.pruned_reliability
                << ",\"runtime_seconds\":" << search.elapsed() << "}";
        }
    } catch (const TimeLimit&) {
        out << "{\"status\":\"timeout\",\"time_limit_seconds\":" << time_limit
            << ",\"nodes\":" << search.nodes
            << ",\"runtime_seconds\":" << search.elapsed() << "}";
    }
    std::cout << out.str() << "\n";
    return 0;
}
