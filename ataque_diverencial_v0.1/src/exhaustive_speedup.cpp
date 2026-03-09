#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace NTL;

GF2X parsePolynomial(const string& s) {
    GF2X result;
    if (s.empty()) return result;
    long val = 0;
    bool in_num = false;
    bool neg = false;
    for (char ch : s) {
        if (ch == '-') {
            neg = true;
        } else if (ch >= '0' && ch <= '9') {
            in_num = true;
            val = val * 10 + (ch - '0');
        } else {
            if (in_num) {
                long e = neg ? -val : val;
                if (e >= 0) SetCoeff(result, e);
                val = 0;
                in_num = false;
                neg = false;
            }
        }
    }
    if (in_num) {
        long e = neg ? -val : val;
        if (e >= 0) SetCoeff(result, e);
    }
    return result;
}

GF2E fromMask(uint64_t mask) {
    GF2X p;
    for (int i = 0; i < 64; ++i) {
        if ((mask >> i) & 1ULL) SetCoeff(p, i);
    }
    return conv<GF2E>(p);
}

int filter_function(const GF2E& state) {
    GF2X poly = rep(state);
    long s0 = IsOne(coeff(poly, 0));
    long s10 = IsOne(coeff(poly, 10));
    long s20 = IsOne(coeff(poly, 20));
    return static_cast<int>(s20 ^ (s0 & s10));
}

vector<int> generate_stream(GF2E x, const GF2E& a, const GF2E& c, long rounds) {
    vector<int> out;
    out.reserve(rounds);
    for (long t = 0; t < rounds; ++t) {
        out.push_back(filter_function(x));
        x = a * x + c;
    }
    return out;
}

bool eq_vec(const vector<int>& a, const vector<int>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return false;
    return true;
}

int main() {
    // Instância pequena para demonstrar speedup de busca exaustiva com chave relacionada K xor Delta
    // Reusa o campo GF(2^31) já usado no projeto, mas restringe o espaço de busca a 2^k candidatos.
    string h_str = "31 13 8 3 0";
    string a_str = "4 0";
    string c_str = "0";

    GF2X h = parsePolynomial(h_str);
    GF2E::init(h);

    GF2E a = conv<GF2E>(parsePolynomial(a_str));
    GF2E c = conv<GF2E>(parsePolynomial(c_str));

    const long k = 12;          // espaço de busca: 2^12
    const long rounds = 64;     // bits observados por stream
    const uint64_t delta = 1ULL; // K' = K xor 1

    // segredo de teste dentro do espaço 2^k
    const uint64_t secret = 0b101101001011ULL;

    GF2E secret_key = fromMask(secret);
    GF2E secret_rel = fromMask(secret ^ delta);

    vector<int> obs1 = generate_stream(secret_key, a, c, rounds);
    vector<int> obs2 = generate_stream(secret_rel, a, c, rounds);

    // Baseline: testa todos os candidatos (K) em [0, 2^k)
    long baseline_trials = 0;
    uint64_t baseline_found = UINT64_MAX;

    for (uint64_t cand = 0; cand < (1ULL << k); ++cand) {
        baseline_trials++;
        vector<int> s1 = generate_stream(fromMask(cand), a, c, rounds);
        vector<int> s2 = generate_stream(fromMask(cand ^ delta), a, c, rounds);
        if (eq_vec(s1, obs1) && eq_vec(s2, obs2)) {
            baseline_found = cand;
            break;
        }
    }

    // Otimizado: testa apenas representantes de pares {K, K xor Delta}
    // Com Delta=1, os pares são (par, ímpar). Testamos apenas os pares e aceitamos qualquer ordem das streams.
    long opt_trials = 0;
    uint64_t opt_found_repr = UINT64_MAX;

    for (uint64_t cand = 0; cand < (1ULL << k); cand += 2) {
        opt_trials++;
        vector<int> s1 = generate_stream(fromMask(cand), a, c, rounds);
        vector<int> s2 = generate_stream(fromMask(cand ^ delta), a, c, rounds);

        bool order1 = eq_vec(s1, obs1) && eq_vec(s2, obs2);
        bool order2 = eq_vec(s1, obs2) && eq_vec(s2, obs1);

        if (order1 || order2) {
            opt_found_repr = cand;
            break;
        }
    }

    bool same_pair = false;
    if (baseline_found != UINT64_MAX && opt_found_repr != UINT64_MAX) {
        uint64_t bf0 = baseline_found & ~1ULL;
        uint64_t of0 = opt_found_repr & ~1ULL;
        same_pair = (bf0 == of0);
    }

    double speedup = (opt_trials > 0)
        ? static_cast<double>(baseline_trials) / static_cast<double>(opt_trials)
        : 0.0;

    cout << "=== Exhaustive Search Speedup Check ===\n";
    cout << "field: GF(2^31), search-space: 2^" << k << "\n";
    cout << "secret key      : " << secret << "\n";
    cout << "baseline trials : " << baseline_trials << "\n";
    cout << "opt trials      : " << opt_trials << "\n";
    cout << "speedup         : " << speedup << "x\n";
    cout << "baseline found  : " << (baseline_found == UINT64_MAX ? -1 : static_cast<long long>(baseline_found)) << "\n";
    cout << "opt repr found  : " << (opt_found_repr == UINT64_MAX ? -1 : static_cast<long long>(opt_found_repr)) << "\n";
    cout << "consistency     : " << (same_pair ? "CONSISTENTE" : "INCONSISTENTE") << "\n";

    return same_pair ? 0 : 2;
}
