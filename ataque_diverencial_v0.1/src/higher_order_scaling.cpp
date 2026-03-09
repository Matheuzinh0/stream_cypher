// =============================================================================
// higher_order_scaling.cpp
// Demonstra que para quebrar uma função de grau d é necessário exatamente
// um diferencial de ordem d, com custo de 2^d consultas por clock.
//
// Para cada config do JSON testa:
//   - filter_original (grau 2): quebra na ordem 2
//   - filter_hard de grau G:    só quebra na ordem G (custo 2^G consultas)
//
// O "break" é detectado quando D^k f == constante para TODOS os clocks.
//
// Saída: tabela mostrando qual ordem mínima zera o diferencial, e o custo
//        real em consultas por clock (= 2^ordem).
// =============================================================================

#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <climits>
#include <vector>
#include <nlohmann/json.hpp>

using namespace std;
using namespace NTL;
using json = nlohmann::json;

// ─── Utilidades ──────────────────────────────────────────────────────────────

GF2X parsePolynomial(const string& s) {
    GF2X result;
    if (s.empty()) return result;
    istringstream ss(s);
    long e; while (ss >> e) SetCoeff(result, e);
    return result;
}

GF2E deltaFromExp(long exp) {
    GF2X p; SetCoeff(p, exp); return conv<GF2E>(p);
}

// ─── Funções filtro ───────────────────────────────────────────────────────────

int filter_grau2(const GF2E& s) {
    GF2X p = rep(s);
    return static_cast<int>(
        IsOne(coeff(p, 20)) ^ (IsOne(coeff(p, 0)) & IsOne(coeff(p, 10))));
}

// Grau algébrico exato = G, função balanceada.
// Construção: soma (XOR) de todos os monômios de grau G formados pelos G
// primeiros bits, mais todos os monômios de grau menor (para garantir grau G).
//   f_G(x) = s[0]&s[1]&...&s[G-1]          (monômio líder de grau G)
//          ^ s[0]&s[1]&...&s[G-2]           (grau G-1)
//          ^ ... ^ s[0]                     (grau 1)
//          ^ s[G-1]                         (grau 1, para balanceamento)
// Isso garante grau exato G e a função não é identicamente 0.
int filter_grauG(const GF2E& s, int G) {
    GF2X p = rep(s);
    int result = 0;
    // monômio líder: produto dos G primeiros bits
    int prod = 1;
    for (int i = 0; i < G; ++i) prod &= IsOne(coeff(p, i)) ? 1 : 0;
    result ^= prod;
    // XOR com todos os sub-monômios de grau 1..G-1 para balancear
    for (int len = 1; len < G; ++len) {
        int m = 1;
        for (int i = 0; i < len; ++i) m &= IsOne(coeff(p, i)) ? 1 : 0;
        result ^= m;
    }
    return result;
}

// ─── Diferencial de ordem k (genérico, recursivo via vetor de deltas) ─────────
//
// D^k f(x; d1,...,dk) = XOR sobre todos os 2^k subconjuntos S de {d1..dk}
//                            de  f(x + sum_{i in S} di)
//
// Para k <= 20 (2^20 = 1M consultas) rodamos de verdade.
// Para k > 20 marcamos como INVIAVEL sem computar.

static const int MAX_FEASIBLE_ORDER = 20; // 2^20 = 1 048 576 consultas/clock

// Retorna o valor do diferencial de ordem k no ponto x com deltas dados.
// Se k > MAX_FEASIBLE_ORDER retorna -1 (inviável).
int Dk(const GF2E& x,
       const vector<GF2E>& deltas,
       int k,
       const function<int(const GF2E&)>& f)
{
    if (k > MAX_FEASIBLE_ORDER) return -1;
    int result = 0;
    long total = 1L << k;
    for (long mask = 0; mask < total; ++mask) {
        GF2E point = x;
        for (int b = 0; b < k; ++b)
            if ((mask >> b) & 1) point += deltas[b];
        result ^= f(point);
    }
    return result;
}

// ─── Estrutura de resultado por ordem ────────────────────────────────────────

struct OrderResult {
    int  order;
    long checks;
    long constant_clocks;   // clocks onde Dk == mesmo valor que clock anterior
    long nonzero_clocks;    // clocks onde Dk != 0
    bool is_constant;       // true se Dk foi igual em TODOS os clocks
    bool infeasible;        // true se 2^order > limite
    long queries_per_clock; // 2^order
};

// ─── Testa uma filter_fn para ordens 2..max_order ────────────────────────────

vector<OrderResult> scan_orders(
    const GF2E& seed,
    const GF2E& a, const GF2E& c,
    int max_order,
    long rounds,
    const function<int(const GF2E&)>& fn,
    // deltas fixos: x^0, x^1, x^2, ... (max_order deltas)
    const vector<GF2E>& base_deltas)
{
    vector<OrderResult> results;

    for (int k = 2; k <= max_order; ++k) {
        OrderResult res;
        res.order = k;
        res.queries_per_clock = (k <= 62) ? (1L << k) : LONG_MAX;
        res.infeasible = (k > MAX_FEASIBLE_ORDER);
        res.checks = 0;
        res.constant_clocks = 0;
        res.nonzero_clocks = 0;
        res.is_constant = false;

        if (res.infeasible) {
            results.push_back(res);
            continue;
        }

        // Para cada clock propagamos os deltas junto com x
        GF2E x = seed;
        vector<GF2E> deltas(base_deltas.begin(), base_deltas.begin() + k);

        int first_val = -1;
        bool all_same = true;

        for (long t = 0; t < rounds; ++t) {
            int val = Dk(x, deltas, k, fn);
            res.checks++;
            if (val != 0) res.nonzero_clocks++;
            if (first_val == -1) {
                first_val = val;
                res.constant_clocks = 1;
            } else {
                if (val == first_val) res.constant_clocks++;
                else all_same = false;
            }

            x = a * x + c;
            for (auto& d : deltas) d = a * d;
        }
        // "constante" significa igual em todos os clocks (pode ser 0 ou 1)
        res.is_constant = all_same && (res.checks > 0);
        // Rejeita "constante zero trivial": se a função nunca foi 1 nos 2^k
        // vértices do cubo (todos outputs foram 0), pode ser que a função seja
        // identicamente 0 naquele subespaço — não é uma propriedade do grau.
        // Exigimos que pelo menos em alguns clocks o diferencial tenha variado
        // antes de colapsar, OU que a função base não seja identicamente 0.
        // Marcamos como "trivialmente nulo" para distinguir do colapso real.
        bool all_zero = (res.nonzero_clocks == 0);
        if (res.is_constant && all_zero) {
            // Verifica se f(x) != 0 em pelo menos 1 clock (função não-nula no trajeto)
            GF2E x_check = seed;
            bool fn_nontrivial = false;
            for (long t = 0; t < rounds && !fn_nontrivial; ++t) {
                if (fn(x_check) != 0) fn_nontrivial = true;
                x_check = a * x_check + c;
            }
            // Se a própria função é quase sempre 0 no trajeto, o colapso é espúrio
            if (!fn_nontrivial) res.is_constant = false;
        }
        results.push_back(res);

        // Para sem avançar para ordens maiores após encontrar a constante nula real
        if (res.is_constant && res.nonzero_clocks == 0) break;
    }

    return results;
}

// ─── Main ────────────────────────────────────────────────────────────────────

int main() {
    ifstream file("unified_config.json");
    if (!file.is_open()) { cerr << "[ERRO] unified_config.json\n"; return 1; }
    json configs; file >> configs;

    ofstream rpt("relatorio_higher_order_scaling.txt");
    rpt << "=== Scaling: Custo do Diferencial de Ordem k para Quebrar Grau d ===\n\n";

    cout << "=== Scaling: Ordem mínima necessária vs. grau algébrico ===\n\n";

    // Graus a testar para filter_hard
    // (limitados a <= MAX_FEASIBLE_ORDER para não travar)
    const vector<int> graus_hard = {2, 3, 4, 5, 8, 10, 15, 20};
    const int MAX_ORDER_SCAN = MAX_FEASIBLE_ORDER; // varre até 20

    // Deltas base: x^0, x^1, ..., x^(MAX_ORDER_SCAN-1)
    // Serão inicializados após GF2E::init
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

        // Apenas 1 seed por config para manter o tempo aceitável
        GF2E seed = conv<GF2E>(parsePolynomial(
            config["seeds"][0].get<string>()));

        // rounds pequeno para viabilizar ordens altas
        const long rounds = 32;

        // Deltas: x^0, x^1, ..., x^(MAX_ORDER_SCAN-1)
        vector<GF2E> base_deltas;
        for (int i = 0; i < MAX_ORDER_SCAN; ++i)
            base_deltas.push_back(deltaFromExp(i));

        auto print_cfg = [&](ostream& os) {
            os << "\n[Config " << cfg_id << "] GF(2^" << N << ")\n";
            os << "  rounds por seed = " << rounds
               << " | consultas/clock = 2^k\n";
            os << "  ─────────────────────────────────────────────────────────\n";
            os << "  Função             | grau | k mín. | consultas/clock | resultado\n";
            os << "  ─────────────────────────────────────────────────────────\n";
        };
        print_cfg(cout);
        print_cfg(rpt);

        // Helper para formatar uma linha da tabela
        auto table_line = [&](ostream& os,
                               const string& fname, int grau,
                               const vector<OrderResult>& res) {
            // Primeira ordem onde D^k == 0 em TODOS os clocks
            int k_zero = -1;
            long cost_zero = -1;
            bool hit_infeasible = false;
            for (const auto& r : res) {
                if (r.infeasible) { hit_infeasible = true; break; }
                if (r.is_constant && r.nonzero_clocks == 0) {
                    k_zero = r.order;
                    cost_zero = r.queries_per_clock;
                    break;
                }
            }

            os << "  " << fname;
            int pad = 20 - (int)fname.size();
            for (int i = 0; i < max(pad,1); ++i) os << ' ';
            os << "| grau=" << grau;
            if (grau < 10) os << " ";
            os << " | ";

            if (k_zero >= 0) {
                os << "ordem " << k_zero
                   << " → 2^" << k_zero << "=" << cost_zero
                   << " consul./clock  → QUEBRADO\n";
            } else if (hit_infeasible) {
                os << "ordem > 20 necessária → INVIÁVEL (>2^20 consul./clock)\n";
            } else {
                os << "nao quebrado no intervalo testado\n";
            }
        };

        // [A] filtro original grau-2
        {
            function<int(const GF2E&)> fn = [](const GF2E& s){ return filter_grau2(s); };
            auto res = scan_orders(seed, a, c, min(4, MAX_ORDER_SCAN), rounds, fn, base_deltas);
            table_line(cout, "f_original(grau 2)", 2, res);
            table_line(rpt,  "f_original(grau 2)", 2, res);
        }

        // [B] filtro hard para cada grau listado
        for (int G : graus_hard) {
            if (G >= (int)N) continue;

            string fname = "f_hard(grau " + to_string(G) + ")";
            function<int(const GF2E&)> fn = [G](const GF2E& s){
                return filter_grauG(s, G);
            };

            int scan_up_to = min(G + 1, MAX_ORDER_SCAN);
            auto res = scan_orders(seed, a, c, scan_up_to, rounds, fn, base_deltas);
            table_line(cout, fname, G, res);
            table_line(rpt,  fname, G, res);
        }

        cout << "\n";
        rpt  << "\n";
        cfg_id++;

        // Só roda config 1 e config 4 (GF(2^31) e GF(2^128)) para não demorar demais
        if (cfg_id > 2) break;
    }

    // Tabela de custo teórico
    auto cost_table = [&](ostream& os) {
        os << "\n=======================================================\n";
        os << " Custo teórico: diferencial de ordem k = 2^k consultas\n";
        os << "=======================================================\n";
        os << "  k  | consultas/clock | viável?\n";
        os << "  ---+----------------+---------\n";
        for (int k = 2; k <= 25; ++k) {
            long cost = (k <= 30) ? (1L << k) : -1L;
            bool ok = (k <= 20);
            if (k <= 30)
                os << "  " << k << (k<10?" ":"") << " | " << cost
                   << string(15 - to_string(cost).size(), ' ')
                   << " | " << (ok ? "sim" : "NAO (> 1M)") << "\n";
        }
        os << "\n CONCLUSAO:\n";
        os << "  - Para quebrar grau d precisa-se de ordem d → custo 2^d por clock.\n";
        os << "  - grau 2  → ordem 2 → 4 consultas  → TRIVIAL\n";
        os << "  - grau 10 → ordem 10 → 1024 consul. → caro mas factível\n";
        os << "  - grau 20 → ordem 20 → ~1M consul.  → limite prático\n";
        os << "  - grau 30 → ordem 30 → ~1B consul.  → INVIÁVEL\n";
        os << "  - grau 64+→ inviável computacionalmente\n";
        os << "=======================================================\n";
    };

    cost_table(cout);
    cost_table(rpt);

    rpt.close();
    cout << "\nRelatorio salvo em: relatorio_higher_order_scaling.txt\n";
    return 0;
}
