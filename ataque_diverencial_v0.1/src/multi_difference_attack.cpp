#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <fstream>
#include <iostream>
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

GF2E deltaFromExponent(long exponent) {
    GF2X p;
    SetCoeff(p, exponent);
    return conv<GF2E>(p);
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

struct DeltaStats {
    long equation_checks = 0;
    long equation_matches = 0;
    long recovered_s0 = 0;
    long recovered_s10 = 0;
    long recovered_both = 0;
};

int main() {
    ifstream file("unified_config.json");
    if (!file.is_open()) {
        cerr << "[ERRO] unified_config.json nao encontrado.\n";
        return 1;
    }

    json configs;
    file >> configs;

    ofstream report("relatorio_multi_difference.txt");
    report << "=== Relatorio: Multi-Difference Differential Attack ===\n";

    // Conjunto de diferencas iniciais (multiplas caracteristicas)
    const vector<long> delta_exponents = {0, 1, 2, 5, 10};

    cout << "=== Multi-Difference Differential Attack ===\n";
    report << "=== Multi-Difference Differential Attack ===\n";
    cout << "Deltas iniciais (expoentes): ";
    report << "Deltas iniciais (expoentes): ";
    for (size_t i = 0; i < delta_exponents.size(); ++i) {
        cout << delta_exponents[i] << (i + 1 == delta_exponents.size() ? "\n" : ", ");
        report << delta_exponents[i] << (i + 1 == delta_exponents.size() ? "\n" : ", ");
    }

    bool all_configs_consistent = true;
    long config_id = 1;

    for (const auto& config : configs) {
        string h_str = config["base_params"]["h"];
        string a_str = config["base_params"]["a"];
        string c_str = config["base_params"]["c"];

        GF2X h = parsePolynomial(h_str);
        GF2E::init(h);

        GF2E a = conv<GF2E>(parsePolynomial(a_str));
        GF2E c = conv<GF2E>(parsePolynomial(c_str));

        const long rounds = max(128L, 4L * deg(h));

        DeltaStats stats;

        for (const auto& seed_val : config["seeds"]) {
            GF2E x = conv<GF2E>(parsePolynomial(seed_val.get<string>()));

            vector<GF2E> deltas;
            deltas.reserve(delta_exponents.size());
            for (long exp : delta_exponents) deltas.push_back(deltaFromExponent(exp));

            for (long t = 0; t < rounds; ++t) {
                int z = filter_function(x);
                int s0_real = bitAt(x, 0);
                int s10_real = bitAt(x, 10);

                bool found_s0 = false;
                bool found_s10 = false;
                int s0_est = 0;
                int s10_est = 0;

                for (size_t i = 0; i < deltas.size(); ++i) {
                    const GF2E& d = deltas[i];
                    int d0 = bitAt(d, 0);
                    int d10 = bitAt(d, 10);
                    int d20 = bitAt(d, 20);

                    int zp = filter_function(x + d);
                    int dz = z ^ zp;

                    int predicted = d20 ^ (s0_real & d10) ^ (s10_real & d0) ^ (d0 & d10);

                    stats.equation_checks++;
                    if (dz == predicted) stats.equation_matches++;

                    // Casos lineares uteis para recuperar bits por multiplas diferencas
                    // d0=0,d10=1 => dz ^ d20 = s0
                    if (d0 == 0 && d10 == 1 && !found_s0) {
                        s0_est = dz ^ d20;
                        found_s0 = true;
                        if (s0_est == s0_real) stats.recovered_s0++;
                    }

                    // d0=1,d10=0 => dz ^ d20 = s10
                    if (d0 == 1 && d10 == 0 && !found_s10) {
                        s10_est = dz ^ d20;
                        found_s10 = true;
                        if (s10_est == s10_real) stats.recovered_s10++;
                    }
                }

                if (found_s0 && found_s10) {
                    stats.recovered_both++;
                }

                x = a * x + c;
                for (auto& d : deltas) d = a * d;
            }
        }

        double eq_rate = (stats.equation_checks > 0)
            ? static_cast<double>(stats.equation_matches) / stats.equation_checks
            : 0.0;

        bool cfg_ok = (eq_rate == 1.0);
        all_configs_consistent = all_configs_consistent && cfg_ok;

        cout << "\n[Config " << config_id << "] GF(2^" << deg(h) << ")\n";
           report << "\n[Config " << config_id << "] GF(2^" << deg(h) << ")\n";
        cout << "  equation matches : " << stats.equation_matches << " / " << stats.equation_checks
             << " (rate=" << eq_rate << ")\n";
           report << "  equation matches : " << stats.equation_matches << " / " << stats.equation_checks
                << " (rate=" << eq_rate << ")\n";
        cout << "  recovered s0     : " << stats.recovered_s0 << "\n";
           report << "  recovered s0     : " << stats.recovered_s0 << "\n";
        cout << "  recovered s10    : " << stats.recovered_s10 << "\n";
           report << "  recovered s10    : " << stats.recovered_s10 << "\n";
        cout << "  recovered both   : " << stats.recovered_both << "\n";
           report << "  recovered both   : " << stats.recovered_both << "\n";
        cout << "  decision         : " << (cfg_ok ? "CONSISTENTE" : "INCONSISTENTE") << "\n";
           report << "  decision         : " << (cfg_ok ? "CONSISTENTE" : "INCONSISTENTE") << "\n";

        config_id++;
    }

    cout << "\n=== RESULTADO GLOBAL: "
         << (all_configs_consistent ? "CONSISTENTE" : "INCONSISTENTE")
         << " ===\n";
    report << "\n=== RESULTADO GLOBAL: "
           << (all_configs_consistent ? "CONSISTENTE" : "INCONSISTENTE")
           << " ===\n";
    report.close();
    cout << "Relatorio salvo em: relatorio_multi_difference.txt\n";

    return all_configs_consistent ? 0 : 2;
}
