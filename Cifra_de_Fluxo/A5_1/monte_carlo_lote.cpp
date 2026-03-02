#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <nlohmann/json.hpp>

using namespace std;
using namespace NTL;
using json = nlohmann::json;

// --- Funções Auxiliares ---
GF2X parsePolynomial(const string& s) {
    GF2X result;
    if (s.empty()) return result;
    istringstream ss(s);
    long coeff;
    while (ss >> coeff) SetCoeff(result, coeff);
    return result;
}

// Calcula o Peso de Hamming do estado polinomial
long getHammingWeight(const GF2E& element) {
    GF2X poly = rep(element);
    long weight = 0;
    for (long i = 0; i <= deg(poly); ++i) {
        if (IsOne(coeff(poly, i))) weight++;
    }
    return weight;
}

struct PRNGState {
    GF2E x, a, c;
};

// --- MOTOR DE CLOCK IRREGULAR COM TAPS DINÂMICOS ---
// Usa a regra da maioria entre 3 bits definidos por c1, c2 e c3
bool decide_clock(const GF2E& state, long c1, long c2, long c3) {
    GF2X poly = rep(state);
    long b1 = IsOne(coeff(poly, c1));
    long b2 = IsOne(coeff(poly, c2));
    long b3 = IsOne(coeff(poly, c3));
    
    long maioria = (b1 & b2) ^ (b1 & b3) ^ (b2 & b3);
    return (maioria == 1);
}

int main() {
    // 1. Leitura do JSON Unificado
    ifstream file("unified_config.json");
    if(!file.is_open()) {
        cerr << "Erro: Arquivo 'unified_config.json' não encontrado!" << endl;
        return 1;
    }
    
    json configs;
    file >> configs;
    cout << "JSON carregado! " << configs.size() << " configurações encontradas." << endl;

    // 2. Arquivo de Log
    ofstream log_file("relatorio_monte_carlo_irregular.txt");
    log_file << "========================================================\n";
    log_file << " SIMULAÇÃO MONTE CARLO: EFEITO ESPACIAL DOS BITS DE CLOCK\n";
    log_file << "========================================================\n\n";

    // Parâmetros da Simulação
    int num_samples = 1000; // 1000 amostras por teste (suficiente para significância estatística rápida)
    int blanking_rounds = 100; // Rodadas de "blindagem" iniciais
    int config_id = 1;

    // 3. Loop sobre as Arquiteturas de LCG
    for (const auto& config : configs) {
        string h_str = config["base_params"]["h"];
        string a_str = config["base_params"]["a"];
        string c_str = config["base_params"]["c"];

        GF2X h = parsePolynomial(h_str);
        GF2E::init(h); 
        long degree = deg(h);

        GF2E a = conv<GF2E>(parsePolynomial(a_str));
        GF2E c = conv<GF2E>(parsePolynomial(c_str));

        log_file << "--------------------------------------------------------\n";
        log_file << "CONFIGURAÇÃO " << config_id++ << " | Grau N = " << degree << "\n";
        log_file << " h(x) = " << h_str << "\n";
        log_file << " a(x) = " << a_str << "\n";
        log_file << "--------------------------------------------------------\n";

        // 4. Definição Esparsa dos Bits de Controle (Taps)
        // Criamos 3 baterias de teste variando onde ficam os bits de clock
        vector<vector<long>> tap_sets = {
            {1, 2, 3},                                     // Conjunto 1: Colados no S_0
            {degree / 4, degree / 2, (3 * degree) / 4},    // Conjunto 2: Espalhados (A5/1 style)
            {degree - 3, degree - 2, degree - 1}           // Conjunto 3: No fim do registrador
        };

        vector<string> desc_taps = {
            "CONCENTRADOS NO INICIO (A diferença atinge o relógio quase instantaneamente)",
            "ESPALHADOS ESTRATEGICAMENTE (Difusão regular)",
            "CONCENTRADOS NO FINAL (A diferença viaja o máximo possível antes do impacto)"
        };

        // 5. Roda a Simulação para cada conjunto de Taps
        for (size_t t_idx = 0; t_idx < tap_sets.size(); ++t_idx) {
            long c1 = tap_sets[t_idx][0];
            long c2 = tap_sets[t_idx][1];
            long c3 = tap_sets[t_idx][2];

            // Trava de segurança caso o polinômio seja menor que 3 bits
            if(c1 >= degree) c1 = degree - 1;
            if(c2 >= degree) c2 = degree - 1;
            if(c3 >= degree) c3 = degree - 1;

            log_file << "\n  [>] TAPS DE CONTROLE: {" << c1 << ", " << c2 << ", " << c3 << "}\n";
            log_file << "      Topologia: " << desc_taps[t_idx] << "\n";

            int sobrevivencias_perfeitas = 0;
            vector<double> avg_hw_over_time(blanking_rounds, 0.0);

            // SIMULAÇÃO MONTE CARLO
            for (int sample = 0; sample < num_samples; ++sample) {
                // Estado original aleatório
                PRNGState secret_state;
                secret_state.x = random_GF2E();
                secret_state.a = a;
                secret_state.c = c;

                // Estado com diferença (Injetando Delta = 1 na posição S_0)
                PRNGState dist_state = secret_state;
                GF2E delta_inicial = conv<GF2E>(parsePolynomial("0"));
                dist_state.x = secret_state.x + delta_inicial;

                bool perdeu_sincronismo = false;

                // Evolução temporal (Blanking rounds)
                for (int t = 0; t < blanking_rounds; ++t) {
                    GF2E diferenca_atual = secret_state.x + dist_state.x;
                    avg_hw_over_time[t] += getHammingWeight(diferenca_atual);

                    // Verifica se a diferença corrompeu a decisão de clock
                    bool clock_original = decide_clock(secret_state.x, c1, c2, c3);
                    bool clock_modificado = decide_clock(dist_state.x, c1, c2, c3);

                    if (clock_original != clock_modificado) {
                        perdeu_sincronismo = true;
                    }

                    // Se não perdeu sincronismo até aqui, avança os LCGs
                    if (clock_original) secret_state.x = secret_state.a * secret_state.x + secret_state.c;
                    if (clock_modificado) dist_state.x = dist_state.a * dist_state.x + dist_state.c;
                }

                if (!perdeu_sincronismo) {
                    sobrevivencias_perfeitas++;
                }
            }

            // Exibição dos Dados do Lote
            double taxa_sobrevivencia = (double)sobrevivencias_perfeitas / num_samples * 100.0;
            log_file << "      -> TAXA DE SOBREVIVÊNCIA DA DIFERENÇA: " << fixed << setprecision(2) << taxa_sobrevivencia << "%\n";
            log_file << "      -> EVOLUÇÃO DO PESO DE HAMMING (Difusão Média):\n";
            
            for (int t = 0; t < blanking_rounds; t += 10) {
                double avg_hw = avg_hw_over_time[t] / num_samples;
                log_file << "         T=" << setw(3) << t << " | Média de " << fixed << setprecision(2) << avg_hw << " bits diferentes\n";
            }
        }
        log_file << "\n";
    }

    log_file.close();
    cout << "Processamento finalizado! Leia o 'relatorio_monte_carlo_irregular.txt'." << endl;
    return 0;
}
