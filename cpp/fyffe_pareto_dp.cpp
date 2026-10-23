#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

// Exact cost-bounded dynamic program for the adapted Fyffe catalog.
// At a fixed layer and cost, first-order stochastically dominated capped
// capacity distributions are discarded. Common non-negative future
// convolutions preserve that order, so the pruning is lossless.

struct Option {
    double p;
    int cost;
    int nominal_weight;
    int choice;
};

struct State {
    std::vector<double> pmf;
    std::vector<double> tail;
    std::vector<int> chromosome;
};

static const std::array<std::array<Option, 4>, 14> CATALOG = {{
    {{{0.00,0,0,0},{0.90,1,3,1},{0.93,1,4,2},{0.94,2,5,3}}},
    {{{0.00,0,0,0},{0.95,2,8,1},{0.94,1,10,2},{0.96,3,11,3}}},
    {{{0.00,0,0,0},{0.85,2,7,1},{0.90,3,5,2},{0.92,4,6,3}}},
    {{{0.00,0,0,0},{0.83,3,5,1},{0.87,4,6,2},{0.89,5,7,3}}},
    {{{0.00,0,0,0},{0.94,2,4,1},{0.93,2,3,2},{0.95,3,5,3}}},
    {{{0.00,0,0,0},{0.99,3,5,1},{0.98,3,4,2},{0.99,4,6,3}}},
    {{{0.00,0,0,0},{0.91,4,7,1},{0.92,4,8,2},{0.94,5,9,3}}},
    {{{0.00,0,0,0},{0.81,3,4,1},{0.90,5,7,2},{0.92,6,8,3}}},
    {{{0.00,0,0,0},{0.97,2,8,1},{0.99,3,9,2},{0.99,4,10,3}}},
    {{{0.00,0,0,0},{0.83,4,6,1},{0.85,4,5,2},{0.88,5,6,3}}},
    {{{0.00,0,0,0},{0.94,3,5,1},{0.95,4,6,2},{0.96,5,7,3}}},
    {{{0.00,0,0,0},{0.79,2,4,1},{0.82,3,5,2},{0.85,4,6,3}}},
    {{{0.00,0,0,0},{0.98,2,5,1},{0.99,3,5,2},{0.99,4,6,3}}},
    {{{0.00,0,0,0},{0.90,4,6,1},{0.92,4,7,2},{0.95,5,8,3}}}
}};

static const std::array<double, 6> ALPHAS = {0.0, 0.3, 0.5, 0.7, 0.9, 1.0};
static const std::array<double, 5> GAMMAS = {0.0, 0.3, 0.5, 0.7, 1.0};
// Independently re-evaluated feasible V1 designs. These are only safe search
// caps; the DP proves whether a lower cost exists and optimizes reliability.
static const int COST_CAP[6][5] = {
    {36,28,23,20,14}, {32,26,23,20,17}, {28,26,23,21,18},
    {26,24,23,22,20}, {24,24,23,23,22}, {23,23,23,23,23}
};

static std::vector<double> make_tail(const std::vector<double>& pmf) {
    std::vector<double> tail(pmf.size());
    double total = 0.0;
    for (int i = static_cast<int>(pmf.size()) - 1; i >= 0; --i) {
        total += pmf[i];
        tail[i] = total;
    }
    return tail;
}

static bool dominates(const State& a, const State& b) {
    constexpr double eps = 1e-12;
    for (std::size_t i = 0; i < a.tail.size(); ++i) {
        if (a.tail[i] + eps < b.tail[i]) return false;
    }
    return true;
}

static State convolve(const State& state, double p, int weight, int choice) {
    const int threshold = static_cast<int>(state.pmf.size()) - 1;
    State out;
    out.pmf.assign(state.pmf.size(), 0.0);
    for (int current = 0; current <= threshold; ++current) {
        const double mass = state.pmf[current];
        if (mass == 0.0) continue;
        out.pmf[current] += mass * (1.0 - p);
        out.pmf[std::min(threshold, current + weight)] += mass * p;
    }
    out.tail = make_tail(out.pmf);
    out.chromosome = state.chromosome;
    out.chromosome.push_back(choice);
    return out;
}

static void insert_nondominated(std::vector<State>& frontier, State&& candidate) {
    for (const State& existing : frontier) {
        if (dominates(existing, candidate)) return;
    }
    frontier.erase(std::remove_if(frontier.begin(), frontier.end(),
        [&](const State& existing) { return dominates(candidate, existing); }), frontier.end());
    frontier.push_back(std::move(candidate));
}

static std::string chromosome_text(const std::vector<int>& chromosome) {
    std::ostringstream stream;
    for (std::size_t i = 0; i < chromosome.size(); ++i) {
        if (i) stream << '-';
        stream << chromosome[i];
    }
    return stream.str();
}

int main(int argc, char** argv) {
    const std::string output = argc > 1 ? argv[1] : "../results/generated/fyffe_pareto_dp.csv";
    std::ofstream csv(output);
    if (!csv) {
        std::cerr << "Cannot open output: " << output << '\n';
        return 2;
    }
    csv << "alpha,gamma,best_cost,best_reliability,chromosome,peak_states,seconds\n";
    csv << std::setprecision(17);

    for (std::size_t ai = 0; ai < ALPHAS.size(); ++ai) {
        for (std::size_t gi = 0; gi < GAMMAS.size(); ++gi) {
            const double alpha = ALPHAS[ai];
            const double gamma = GAMMAS[gi];
            // All configured decimals have at most two effective digits after
            // scalarization. Reduce by gcd to match Python's rational scaling.
            const int factor100 = static_cast<int>(std::llround(
                100.0 * (0.9 + 0.1 * alpha + gamma * (0.2 - 0.2 * alpha))));
            const int threshold100 = static_cast<int>(std::llround(
                100.0 * ((60.0 - 10.0 * alpha) - gamma * (20.0 - 20.0 * alpha))));
            int divisor = threshold100;
            for (const auto& row : CATALOG) {
                for (const auto& option : row) divisor = std::gcd(divisor, option.nominal_weight * factor100);
            }
            const int threshold = threshold100 / divisor;
            const int cap = COST_CAP[ai][gi];
            std::vector<std::vector<State>> states(cap + 1);
            State initial;
            initial.pmf.assign(threshold + 1, 0.0);
            initial.pmf[0] = 1.0;
            initial.tail = make_tail(initial.pmf);
            states[0].push_back(std::move(initial));
            std::size_t peak_states = 1;
            const auto started = std::chrono::steady_clock::now();

            for (const auto& row : CATALOG) {
                std::vector<std::vector<State>> expanded(cap + 1);
                for (int cost = 0; cost <= cap; ++cost) {
                    for (const State& state : states[cost]) {
                        for (const Option& option : row) {
                            const int next_cost = cost + option.cost;
                            if (next_cost > cap) continue;
                            const int weight = option.nominal_weight * factor100 / divisor;
                            insert_nondominated(expanded[next_cost], convolve(state, option.p, weight, option.choice));
                        }
                    }
                }
                states = std::move(expanded);
                std::size_t count = 0;
                for (const auto& bucket : states) count += bucket.size();
                peak_states = std::max(peak_states, count);
            }

            int best_cost = -1;
            double best_reliability = -1.0;
            std::vector<int> best_chromosome;
            for (int cost = 0; cost <= cap && best_cost < 0; ++cost) {
                for (const State& state : states[cost]) {
                    const double reliability = state.pmf[threshold];
                    if (reliability + 1e-12 >= 0.90 && reliability > best_reliability + 1e-12) {
                        best_cost = cost;
                        best_reliability = reliability;
                        best_chromosome = state.chromosome;
                    }
                }
            }
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            if (best_cost < 0) {
                std::cerr << "No feasible state within verified cap at alpha=" << alpha
                          << " gamma=" << gamma << '\n';
                return 3;
            }
            csv << alpha << ',' << gamma << ',' << best_cost << ',' << best_reliability
                << ',' << chromosome_text(best_chromosome) << ',' << peak_states << ',' << seconds << '\n';
            csv.flush();
            std::cout << "alpha=" << alpha << " gamma=" << gamma << " cost=" << best_cost
                      << " R=" << std::setprecision(12) << best_reliability
                      << " states=" << peak_states << " seconds=" << seconds << '\n';
        }
    }
    return 0;
}
