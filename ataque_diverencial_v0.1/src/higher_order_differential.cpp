#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
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

struct Stats {
    long checks = 0;
    long matches = 0;
    long random_checks = 0;
    long random_matches = 0;
};

int main() {
    ifstream file("unified_config.json");
    if (!file.is_open()) {
        cerr << "[ERRO] unified_config.json nao encontrado.\n";
        return 1;
    }

    json configs;
    file >> configs;

    ofstream report("relatorio_higher_order.txt");
    report << "=== Relatorio: Higher-Order Differential (2nd order) ===\n";

    // Diferenciais de 2a ordem: D_{a,b} f(x) = f(x) ^ f(x^a) ^ f(x^b) ^ f(x^a^b)
    // Para f = s20 ^ (s0&s10), o termo de 2a ordem depende de a,b nos bits (0,10):
    // D2 = (a0 & b10) ^ (a10 & b0), independente de x.
    const long alpha_exp = 0;   // alpha0 = 1
    const long beta_exp = 10;   // beta10 = 1

    mt19937_64 rng(0xBADC0FFEEULL);
    uniform_int_distribution<int> bit(0, 1);

    cout << "=== Higher-Order Differential Check (2nd order) ===\n";
    cout << "alpha = x^" << alpha_exp << ", beta = x^" << beta_exp << "\n";
    report << "=== Higher-Order Differential Check (2nd order) ===\n";
    report << "alpha = x^" << alpha_exp << ", beta = x^" << beta_exp << "\n";

    bool all_consistent = true;
    long cfg_id = 1;

    for (const auto& config : configs) {
        string h_str = config["base_params"]["h"];
        string a_str = config["base_params"]["a"];
        string c_str = config["base_params"]["c"];

        GF2X h = parsePolynomial(h_str);
        GF2E::init(h);

        GF2E a = conv<GF2E>(parsePolynomial(a_str));
        GF2E c = conv<GF2E>(parsePolynomial(c_str));

        const long rounds = max(128L, 4L * deg(h));

        Stats stats;

        for (const auto& seed_val : config["seeds"]) {
            GF2E x = conv<GF2E>(parsePolynomial(seed_val.get<string>()));
            GF2E alpha = deltaFromExponent(alpha_exp);
            GF2E beta = deltaFromExponent(beta_exp);

            for (long t = 0; t < rounds; ++t) {
                int f0 = filter_function(x);
                int fa = filter_function(x + alpha);
                int fb = filter_function(x + beta);
                int fab = filter_function(x + alpha + beta);
                int d2_obs = f0 ^ fa ^ fb ^ fab;

                int a0 = bitAt(alpha, 0);
                int a10 = bitAt(alpha, 10);
                int b0 = bitAt(beta, 0);
                int b10 = bitAt(beta, 10);
                int d2_pred = (a0 & b10) ^ (a10 & b0);

                stats.checks++;
                if (d2_obs == d2_pred) stats.matches++;

                // Baseline aleatorio para contraste
                int r_obs = bit(rng) ^ bit(rng) ^ bit(rng) ^ bit(rng);
                stats.random_checks++;
                if (r_obs == d2_pred) stats.random_matches++;

                x = a * x + c;
                alpha = a * alpha;
                beta = a * beta;
            }
        }

        double hit_rate = (stats.checks > 0)
            ? static_cast<double>(stats.matches) / stats.checks
            : 0.0;
        double random_rate = (stats.random_checks > 0)
            ? static_cast<double>(stats.random_matches) / stats.random_checks
            : 0.0;

        bool cfg_ok = (hit_rate == 1.0) && (random_rate > 0.40 && random_rate < 0.60);
        all_consistent = all_consistent && cfg_ok;

        cout << "\n[Config " << cfg_id << "] GF(2^" << deg(h) << ")\n";
           report << "\n[Config " << cfg_id << "] GF(2^" << deg(h) << ")\n";
        cout << "  observed vs predicted : " << stats.matches << " / " << stats.checks
             << " (rate=" << hit_rate << ")\n";
           report << "  observed vs predicted : " << stats.matches << " / " << stats.checks
                << " (rate=" << hit_rate << ")\n";
        cout << "  random baseline rate  : " << random_rate << "\n";
           report << "  random baseline rate  : " << random_rate << "\n";
        cout << "  decision              : " << (cfg_ok ? "CONSISTENTE" : "INCONSISTENTE") << "\n";
           report << "  decision              : " << (cfg_ok ? "CONSISTENTE" : "INCONSISTENTE") << "\n";

        cfg_id++;
    }

    cout << "\n=== RESULTADO GLOBAL: "
         << (all_consistent ? "CONSISTENTE" : "INCONSISTENTE")
         << " ===\n";
    report << "\n=== RESULTADO GLOBAL: "
           << (all_consistent ? "CONSISTENTE" : "INCONSISTENTE")
           << " ===\n";
    report.close();
    cout << "Relatorio salvo em: relatorio_higher_order.txt\n";

    return all_consistent ? 0 : 2;
}
