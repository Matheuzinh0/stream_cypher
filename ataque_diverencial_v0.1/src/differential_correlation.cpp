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

struct CorrAccumulator {
    long n = 0;
    long n11 = 0;
    long n10 = 0;
    long n01 = 0;
    long n00 = 0;

    void add(int x, int y) {
        n++;
        if (x == 1 && y == 1) n11++;
        else if (x == 1 && y == 0) n10++;
        else if (x == 0 && y == 1) n01++;
        else n00++;
    }

    double phi() const {
        double a = static_cast<double>(n11);
        double b = static_cast<double>(n10);
        double c = static_cast<double>(n01);
        double d = static_cast<double>(n00);
        double den = (a + b) * (c + d) * (a + c) * (b + d);
        if (den <= 0.0) return 0.0;
        return (a * d - b * c) / sqrt(den);
    }
};

int main() {
    ifstream file("unified_config.json");
    if (!file.is_open()) {
        cerr << "[ERRO] unified_config.json nao encontrado.\n";
        return 1;
    }

    json configs;
    file >> configs;

    ofstream report("relatorio_differential_correlation.txt");
    report << "=== Relatorio: Differential Correlation Attack ===\n";

    mt19937_64 rng(0x123456789ULL);
    uniform_int_distribution<int> bit(0, 1);

    cout << "=== Differential Correlation Attack Check ===\n";
    report << "=== Differential Correlation Attack Check ===\n";

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

        CorrAccumulator passive_corr_target;      // corr(z, s20)
        CorrAccumulator passive_corr_random;      // corr(rand, s20)
        CorrAccumulator differential_corr_target; // corr(dz(delta=e0), s10)
        CorrAccumulator differential_corr_random; // corr(rand, s10)

        for (const auto& seed_val : config["seeds"]) {
            GF2E x = conv<GF2E>(parsePolynomial(seed_val.get<string>()));
            GF2E e0 = deltaFromExponent(0); // probe diferencial fixo por clock: x_t xor e0

            for (long t = 0; t < rounds; ++t) {
                int z = filter_function(x);
                int s20 = bitAt(x, 20);
                int s10 = bitAt(x, 10);

                int zp = filter_function(x + e0);
                int dz = z ^ zp;

                passive_corr_target.add(z, s20);
                differential_corr_target.add(dz, s10);

                passive_corr_random.add(bit(rng), s20);
                differential_corr_random.add(bit(rng), s10);

                x = a * x + c;
            }
        }

        double phi_passive_target = passive_corr_target.phi();
        double phi_passive_random = passive_corr_random.phi();
        double phi_diff_target = differential_corr_target.phi();
        double phi_diff_random = differential_corr_random.phi();

        bool passive_ok = (phi_passive_target > 0.35) && (abs(phi_passive_random) < 0.08);
        bool diff_ok = (phi_diff_target > 0.95) && (abs(phi_diff_random) < 0.08);
        bool cfg_ok = passive_ok && diff_ok;
        all_consistent = all_consistent && cfg_ok;

        cout << "\n[Config " << cfg_id << "] GF(2^" << deg(h) << ")\n";
        report << "\n[Config " << cfg_id << "] GF(2^" << deg(h) << ")\n";
        cout << "  phi(z, s20) target      : " << phi_passive_target << "\n";
        report << "  phi(z, s20) target      : " << phi_passive_target << "\n";
        cout << "  phi(z, s20) random      : " << phi_passive_random << "\n";
        report << "  phi(z, s20) random      : " << phi_passive_random << "\n";
        cout << "  phi(dz, s10) target     : " << phi_diff_target << "\n";
        report << "  phi(dz, s10) target     : " << phi_diff_target << "\n";
        cout << "  phi(dz, s10) random     : " << phi_diff_random << "\n";
        report << "  phi(dz, s10) random     : " << phi_diff_random << "\n";
        cout << "  decision                : " << (cfg_ok ? "CONSISTENTE" : "INCONSISTENTE") << "\n";
        report << "  decision                : " << (cfg_ok ? "CONSISTENTE" : "INCONSISTENTE") << "\n";

        cfg_id++;
    }

    cout << "\n=== RESULTADO GLOBAL: "
         << (all_consistent ? "CONSISTENTE" : "INCONSISTENTE")
         << " ===\n";
    report << "\n=== RESULTADO GLOBAL: "
           << (all_consistent ? "CONSISTENTE" : "INCONSISTENTE")
           << " ===\n";
    report.close();
    cout << "Relatorio salvo em: relatorio_differential_correlation.txt\n";

    return all_consistent ? 0 : 2;
}
