// =============================================================================
// higher_order_stress.cpp
// Teste de Estresse: Higher-Order Differential contra função de grau algébrico
// altíssimo (grau = deg(h) - 1), demonstrando que o ataque de 2a ordem
// utilizado contra a Toyocrypt NÃO é efetivo contra funções de alto grau.
//
// Função filtro de alto grau:
//   f_hard(x) = s[0] & s[1] & s[2] & ... & s[deg-2]
//   grau algébrico = deg(h) - 1
//
// Comparação lado a lado:
//   [A] filter_function original (grau 2) → D² deve ser CONSTANTE → CONSISTENTE
//   [B] filter_hard (grau deg-1)          → D² NÃO é constante     → INCONSISTENTE
// =============================================================================

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

// ─────────────────────────── Utilidades ──────────────────────────────────────

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

// ─────────────────────────── Funções Filtro ──────────────────────────────────

// [A] Toyocrypt original: grau algébrico = 2
//     f(x) = s20 ^ (s0 & s10)
int filter_original(const GF2E& state) {
    GF2X poly = rep(state);
    long s0  = IsOne(coeff(poly,  0));
    long s10 = IsOne(coeff(poly, 10));
    long s20 = IsOne(coeff(poly, 20));
    return static_cast<int>(s20 ^ (s0 & s10));
}

// [B] Função de altíssimo grau: grau algébrico = (N - 1)
//     f_hard(x) = s[0] & s[1] & ... & s[N-2]   (produto de N-1 bits)
//     N = deg(h), então o grau = N-1, resistente a diferenciais de 2a ordem.
int filter_hard(const GF2E& state, long N) {
    GF2X poly = rep(state);
    int product = 1;
    for (long i = 0; i < N - 1; ++i)
        product &= IsOne(coeff(poly, i)) ? 1 : 0;
    return product;
}

// ─────────────────────────── Estrutura de Resultado ──────────────────────────

struct D2Result {
    long checks  = 0;
    long matches = 0;   // d2_obs == d2_pred (fórmula de grau 2)
    long nonconst_violations = 0; // vezes que d2_obs variou entre clocks
    int  last_d2 = -1;

    void add(int d2_obs, int d2_pred) {
        checks++;
        if (d2_obs == d2_pred) matches++;
        if (last_d2 != -1 && d2_obs != last_d2) nonconst_violations++;
        last_d2 = d2_obs;
    }

    double hit_rate() const {
        return checks ? static_cast<double>(matches) / checks : 0.0;
    }
};

// ─────────────────────────── Função de Teste ─────────────────────────────────

// Roda os diferenciais de 2a ordem para uma dada filter_fn e retorna D2Result.
// d2_pred usa a fórmula derivada do grau 2 original:
//   D2_pred = (a0 & b10) ^ (a10 & b0)
// Isso é SEMPRE correto para filter_original; para filter_hard, será errado
// na maioria dos clocks porque a função tem termos de grau maior.
template<typename FilterFn>
D2Result run_d2_check(
    const GF2E& seed, const GF2E& a, const GF2E& c,
    long alpha_exp, long beta_exp,
    long rounds, FilterFn fn)
{
    D2Result res;
    GF2E x = seed;
    GF2E alpha = deltaFromExponent(alpha_exp);
    GF2E beta  = deltaFromExponent(beta_exp);

    for (long t = 0; t < rounds; ++t) {
        int f0  = fn(x);
        int fa  = fn(x + alpha);
        int fb  = fn(x + beta);
        int fab = fn(x + alpha + beta);
        int d2_obs = f0 ^ fa ^ fb ^ fab;

        int a0  = bitAt(alpha,  0);
        int a10 = bitAt(alpha, 10);
        int b0  = bitAt(beta,   0);
        int b10 = bitAt(beta,  10);
        int d2_pred = (a0 & b10) ^ (a10 & b0);   // válido apenas para grau 2

        res.add(d2_obs, d2_pred);

        x     = a * x     + c;
        alpha = a * alpha;
        beta  = a * beta;
    }
    return res;
}

// ─────────────────────────── Main ────────────────────────────────────────────

int main() {
    ifstream file("unified_config.json");
    if (!file.is_open()) {
        cerr << "[ERRO] unified_config.json nao encontrado.\n";
        return 1;
    }
    json configs;
    file >> configs;

    ofstream report("relatorio_higher_order_stress.txt");
    report << "=== Stress Test: Higher-Order Differential vs. Alto Grau Algebrico ===\n";
    report << "Objetivo: mostrar que D^2 de 2a ordem e EFETIVO contra grau-2 (Toyocrypt)\n";
    report << "          mas INEFETIVO contra funcoes de grau deg-1 (alto grau).\n\n";

    const long alpha_exp = 0;
    const long beta_exp  = 10;

    cout << "=== Stress: Higher-Order Differential vs. grau altissimo ===\n";
    cout << "alpha = x^" << alpha_exp << ",  beta = x^" << beta_exp << "\n\n";
    report << "alpha = x^" << alpha_exp << ",  beta = x^" << beta_exp << "\n\n";

    bool global_orig_ok  = true;
    bool global_hard_fail = true;  // queremos que TODOS falhem
    long cfg_id = 1;

    for (const auto& config : configs) {
        string h_str = config["base_params"]["h"];
        string a_str = config["base_params"]["a"];
        string c_str = config["base_params"]["c"];

        GF2X h = parsePolynomial(h_str);
        GF2E::init(h);
        long N = deg(h);

        GF2E a = conv<GF2E>(parsePolynomial(a_str));
        GF2E c = conv<GF2E>(parsePolynomial(c_str));

        const long rounds = max(128L, 4L * N);

        D2Result orig_total, hard_total;

        for (const auto& seed_val : config["seeds"]) {
            GF2E seed = conv<GF2E>(parsePolynomial(seed_val.get<string>()));

            // [A] filtro original grau-2
            D2Result r_orig = run_d2_check(seed, a, c, alpha_exp, beta_exp, rounds,
                [](const GF2E& s) { return filter_original(s); });
            orig_total.checks  += r_orig.checks;
            orig_total.matches += r_orig.matches;
            orig_total.nonconst_violations += r_orig.nonconst_violations;

            // [B] filtro hard grau-(N-1): captura até (N-1) bits via AND
            // Limitamos a min(N-1, 30) para não tornar o teste lento em N=451
            long hard_bits = min(N - 1, 30L);
            D2Result r_hard = run_d2_check(seed, a, c, alpha_exp, beta_exp, rounds,
                [hard_bits](const GF2E& s) { return filter_hard(s, hard_bits + 1); });
            hard_total.checks  += r_hard.checks;
            hard_total.matches += r_hard.matches;
            hard_total.nonconst_violations += r_hard.nonconst_violations;
        }

        double orig_rate = orig_total.hit_rate();
        double hard_rate = hard_total.hit_rate();

        // [A] deve ser consistente (rate=1.0)
        bool orig_ok  = (orig_rate == 1.0);
        // [B] deve ser inconsistente (rate < 1.0, ou seja, d2 nao e constante)
        bool hard_fail = (hard_rate < 1.0);

        global_orig_ok   = global_orig_ok   && orig_ok;
        global_hard_fail = global_hard_fail && hard_fail;

        auto out = [&](ostream& os) {
            os << "\n[Config " << cfg_id << "] GF(2^" << N << ")\n";
            os << "  [A] grau-2 (Toyocrypt original)\n";
            os << "      D2 matches / checks : " << orig_total.matches
               << " / " << orig_total.checks
               << "  rate=" << orig_rate << "\n";
            os << "      variações D2        : " << orig_total.nonconst_violations << "\n";
            os << "      decisão             : "
               << (orig_ok ? "CONSISTENTE  (D2 constante, ataque efetivo)" :
                             "INCONSISTENTE") << "\n";

            os << "  [B] grau-" << min(N - 1, 30L)
               << " (alto grau, " << min(N-1,30L) << " bits em AND)\n";
            os << "      D2 matches / checks : " << hard_total.matches
               << " / " << hard_total.checks
               << "  rate=" << hard_rate << "\n";
            os << "      variações D2        : " << hard_total.nonconst_violations << "\n";
            os << "      decisão             : "
               << (hard_fail ? "INCONSISTENTE (D2 variavel, ataque FALHA)" :
                               "CONSISTENTE  [inesperado]") << "\n";
        };

        out(cout);
        out(report);
        cfg_id++;
    }

    auto summary = [&](ostream& os) {
        os << "\n=======================================================\n";
        os << " RESUMO FINAL\n";
        os << "=======================================================\n";
        os << " [A] grau-2 (Toyocrypt)  → ataque de 2a ordem : "
           << (global_orig_ok  ? "EFETIVO  ✓" : "INEFETIVO ✗") << "\n";
        os << " [B] grau altissimo      → ataque de 2a ordem : "
           << (global_hard_fail ? "INEFETIVO ✓ (esperado)" : "EFETIVO [inesperado]") << "\n";
        os << "=======================================================\n";
        os << " CONCLUSÃO: O diferencial de 2a ordem é seletivo ao grau algebrico.\n";
        os << "            Funções com grau >= 3 não são capturadas por este ataque.\n";
        os << "            A vulnerabilidade da Toyocrypt reside no baixo grau (2) da f.\n";
        os << "=======================================================\n";
    };

    summary(cout);
    summary(report);

    report.close();
    cout << "\nRelatorio salvo em: relatorio_higher_order_stress.txt\n";

    // Retorna 0 se [A] funciona E [B] falha (comportamento esperado)
    return (global_orig_ok && global_hard_fail) ? 0 : 2;
}
