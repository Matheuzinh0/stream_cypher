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

GF2X parsePolynomial(const string& s) {
    GF2X result;
    if (s.empty()) return result;
    istringstream ss(s);
    long coeff;
    while (ss >> coeff) SetCoeff(result, coeff);
    return result;
}

struct PRNGState {
    GF2E x, a, c;
};

// Motor Irregular (Regra da Maioria)
bool decide_clock(const GF2E& state, long c1, long c2, long c3) {
    GF2X poly = rep(state);
    long b1 = IsOne(coeff(poly, c1));
    long b2 = IsOne(coeff(poly, c2));
    long b3 = IsOne(coeff(poly, c3));
    
    long maioria = (b1 & b2) ^ (b1 & b3) ^ (b2 & b3);
    return (maioria == 1);
}

int main() {
    ifstream file("unified_config.json");
    if(!file.is_open()) {
        cerr << "Erro: Arquivo 'unified_config.json' não encontrado!" << endl;
        return 1;
    }
    
    json configs;
    file >> configs;
    
    // Gerando um CSV para facilitar a análise de dados e plotagem de gráficos
    ofstream csv_file("mapa_vulnerabilidade_taps.csv");
    csv_file << "ConfigID,Grau_N,Polinomio_A,C1,C2,C3,Taxa_Sobrevivencia_%\n";

    // --- PARÂMETROS ---
    int num_samples = 100; // Reduzido para 100 para viabilizar a combinatória gigante
    int blanking_rounds = 100; 

    int config_id = 1;

    for (const auto& config : configs) {
        string h_str = config["base_params"]["h"];
        string a_str = config["base_params"]["a"];
        string c_str = config["base_params"]["c"];

        GF2X h = parsePolynomial(h_str);
        GF2E::init(h); 
        long N = deg(h);

        GF2E a = conv<GF2E>(parsePolynomial(a_str));
        GF2E c = conv<GF2E>(parsePolynomial(c_str));

        // Calcula total de combinacoes para log
        long total_combinations = (N - 1) * (N - 2) * (N - 3) / 6;
        
        cout << "\n========================================================" << endl;
        cout << "CONFIGURAÇÃO " << config_id << " | Grau N = " << N << endl;
        cout << "Total de Permutações a testar: " << total_combinations << endl;
        cout << "========================================================" << endl;

        long combinacoes_feitas = 0;

        // Loop combinatório: Testar todos os C1 < C2 < C3 possíveis.
        // Começamos do 1 porque o bit 0 é onde injetamos a falha diferencial.
        for (long c1 = 1; c1 < N - 2; ++c1) {
            for (long c2 = c1 + 1; c2 < N - 1; ++c2) {
                for (long c3 = c2 + 1; c3 < N; ++c3) {
                    
                    int sobrevivencias_perfeitas = 0;

                    // SIMULAÇÃO MONTE CARLO (100 amostras por permutação)
                    for (int sample = 0; sample < num_samples; ++sample) {
                        PRNGState secret_state;
                        secret_state.x = random_GF2E();
                        secret_state.a = a;
                        secret_state.c = c;

                        PRNGState dist_state = secret_state;
                        GF2E delta_inicial = conv<GF2E>(parsePolynomial("0"));
                        dist_state.x = secret_state.x + delta_inicial;

                        bool perdeu_sincronismo = false;

                        for (int t = 0; t < blanking_rounds; ++t) {
                            bool clock_original = decide_clock(secret_state.x, c1, c2, c3);
                            bool clock_modificado = decide_clock(dist_state.x, c1, c2, c3);

                            if (clock_original != clock_modificado) {
                                perdeu_sincronismo = true;
                                break; // Otimização: se perdeu, não precisa continuar os 100 clocks
                            }

                            if (clock_original) secret_state.x = secret_state.a * secret_state.x + secret_state.c;
                            if (clock_modificado) dist_state.x = dist_state.a * dist_state.x + dist_state.c;
                        }

                        if (!perdeu_sincronismo) {
                            sobrevivencias_perfeitas++;
                        }
                    }

                    // Calcula a taxa e salva diretamente no CSV
                    double taxa = (double)sobrevivencias_perfeitas / num_samples * 100.0;
                    csv_file << config_id << "," << N << ",\"" << a_str << "\"," 
                             << c1 << "," << c2 << "," << c3 << "," 
                             << fixed << setprecision(2) << taxa << "\n";

                    // Log de progresso a cada 10.000 combinações
                    combinacoes_feitas++;
                    if (combinacoes_feitas % 10000 == 0) {
                        cout << "Progresso Config " << config_id << ": " 
                             << combinacoes_feitas << " / " << total_combinations 
                             << " permutacoes avaliadas...\r" << flush;
                    }
                }
            }
        }
        cout << "\nConfiguracao " << config_id << " concluida!" << endl;
        config_id++;
    }

    csv_file.close();
    cout << "\nVarredura total finalizada! Dados salvos em 'mapa_vulnerabilidade_taps.csv'." << endl;
    return 0;
}
