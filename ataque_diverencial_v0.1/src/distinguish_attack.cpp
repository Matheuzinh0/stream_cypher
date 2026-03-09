#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using namespace std;
using namespace NTL;
using json = nlohmann::json;

GF2X parsePolynomial(const string& s) {
    GF2X result;
    if (s.empty()) return result;
    istringstream ss(s);
    long exponent;
    while (ss >> exponent) SetCoeff(result, exponent);
    return result;
}

int bitAt(const GF2E& value, long idx) {
    return IsOne(coeff(rep(value), idx)) ? 1 : 0;
}

int filter_function(const GF2E& state) {
    GF2X poly = rep(state);
    long s0 = IsOne(coeff(poly, 0));
    long s10 = IsOne(coeff(poly, 10));
    long s20 = IsOne(coeff(poly, 20));
    return static_cast<int>(s20 ^ (s0 & s10));
}

GF2E fromMask(uint64_t mask) {
    GF2X p;
    for (int i = 0; i < 64; ++i) {
        if ((mask >> i) & 1ULL) SetCoeff(p, i);
    }
    return conv<GF2E>(p);
}

struct SeedResult {
    long eligible = 0;
    long matches = 0;
    long mismatches = 0;
};

SeedResult run_target_pair(const GF2E& seed, const GF2E& a, const GF2E& c, long rounds) {
    GF2E x1 = seed;
    GF2E x2 = seed + fromMask(1ULL);  // diferença inicial Delta0 = 1

    GF2E delta = fromMask(1ULL);

    SeedResult result;

    for (long t = 0; t < rounds; ++t) {
        int z1 = filter_function(x1);
        int z2 = filter_function(x2);
        int dz = z1 ^ z2;

        int d0 = bitAt(delta, 0);
        int d10 = bitAt(delta, 10);
        int d20 = bitAt(delta, 20);

        if (d0 == 0 && d10 == 0) {
            result.eligible++;
            if (dz == d20) {
                result.matches++;
            } else {
                result.mismatches++;
            }
        }

        x1 = a * x1 + c;
        x2 = a * x2 + c;
        delta = a * delta;
    }

    return result;
}

SeedResult run_random_pair(const GF2E& a, long rounds, mt19937_64& rng) {
    uniform_int_distribution<int> bit(0, 1);

    GF2E delta = fromMask(1ULL);
    SeedResult result;

    for (long t = 0; t < rounds; ++t) {
        int dz = bit(rng) ^ bit(rng);

        int d0 = bitAt(delta, 0);
        int d10 = bitAt(delta, 10);
        int d20 = bitAt(delta, 20);

        if (d0 == 0 && d10 == 0) {
            result.eligible++;
            if (dz == d20) {
                result.matches++;
            } else {
                result.mismatches++;
            }
        }

        delta = a * delta;
    }

    return result;
}

int main() {
    ifstream file("unified_config.json");
    if (!file.is_open()) {
        cerr << "[ERRO] unified_config.json nao encontrado.\n";
        return 1;
    }

    json configs;
    file >> configs;

    mt19937_64 rng(0xC0FFEEULL);

    cout << "=== Distinguishing Attack Check ===\n";

    bool all_configs_consistent = true;
    long cfg_id = 1;

    for (const auto& config : configs) {
        string h_str = config["base_params"]["h"];
        string a_str = config["base_params"]["a"];
        string c_str = config["base_params"]["c"];

        GF2X h = parsePolynomial(h_str);
        GF2E::init(h);

        GF2E a = conv<GF2E>(parsePolynomial(a_str));
        GF2E c = conv<GF2E>(parsePolynomial(c_str));

        long rounds = max(128L, 4L * deg(h));

        long target_eligible = 0;
        long target_mismatches = 0;

        for (const auto& seed_val : config["seeds"]) {
            GF2E seed = conv<GF2E>(parsePolynomial(seed_val.get<string>()));
            SeedResult r = run_target_pair(seed, a, c, rounds);
            target_eligible += r.eligible;
            target_mismatches += r.mismatches;
        }

        long random_trials = 200;
        long random_eligible = 0;
        long random_mismatches = 0;

        for (long i = 0; i < random_trials; ++i) {
            SeedResult rr = run_random_pair(a, rounds, rng);
            random_eligible += rr.eligible;
            random_mismatches += rr.mismatches;
        }

        double target_mismatch_rate = (target_eligible > 0)
            ? static_cast<double>(target_mismatches) / target_eligible
            : 0.0;

        double random_mismatch_rate = (random_eligible > 0)
            ? static_cast<double>(random_mismatches) / random_eligible
            : 0.0;

        bool target_is_target = (target_mismatches == 0);
        bool random_not_target = (random_mismatch_rate > 0.30);
        bool cfg_ok = target_is_target && random_not_target;

        cout << "\n[Config " << cfg_id << "] GF(2^" << deg(h) << ")\n";
        cout << "  target: eligible=" << target_eligible
             << " mismatches=" << target_mismatches
             << " rate=" << target_mismatch_rate << "\n";
        cout << "  random: eligible=" << random_eligible
             << " mismatches=" << random_mismatches
             << " rate=" << random_mismatch_rate << "\n";
        cout << "  decision: " << (cfg_ok ? "CONSISTENTE" : "INCONSISTENTE") << "\n";

        all_configs_consistent = all_configs_consistent && cfg_ok;
        cfg_id++;
    }

    cout << "\n=== RESULTADO GLOBAL: "
         << (all_configs_consistent ? "CONSISTENTE" : "INCONSISTENTE")
         << " ===\n";

    return all_configs_consistent ? 0 : 2;
}
