#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <nlohmann/json.hpp>
using namespace std;
using namespace NTL;

// Funções Auxiliares
GF2X parsePolynomial(const string& s) {
    GF2X result;
    if (s.empty()) return result;
    istringstream ss(s);
    long coeff;
    while (ss >> coeff) SetCoeff(result, coeff);
    return result;
}

// Calcula o Peso de Hamming (Quantidade de bits 1 no polinômio)
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

// --- O MOTOR IRREGULAR (STOP-AND-GO) ---
// Retorna 'true' se o gerador deve avançar, 'false' se deve congelar
bool decide_clock(const GF2E& state) {
    GF2X poly = rep(state);
    
    // Escolhemos 3 bits arbitrários espaçados para controlar o clock (como no A5/1)
    long b1 = IsOne(coeff(poly, 12));
    long b2 = IsOne(coeff(poly, 34));
    long b3 = IsOne(coeff(poly, 56));
    
    // Regra da Maioria: Retorna 1 se pelo menos dois bits forem 1
    long maioria = (b1 & b2) ^ (b1 & b3) ^ (b2 & b3);
    return (maioria == 1);
}

int main() {
    // 1. Configuração do Corpo de Galois (Usando um polinômio de 64 bits)
    string h_str = "64 42 33 22 0";
    GF2X h = parsePolynomial(h_str);
    GF2E::init(h);

    // Usando um a(x) multiplicador simples para vermos o HW crescer gradualmente
    GF2E a = conv<GF2E>(parsePolynomial("5 0")); 
    GF2E c = conv<GF2E>(parsePolynomial("0"));

    int num_samples = 10000;
    int blanking_rounds = 100; // Os 100 clocks iniciais de inicialização
    
    int sobrevivencias_perfeitas = 0;
    vector<double> avg_hw_over_time(blanking_rounds, 0.0);

    cout << "========================================================" << endl;
    cout << " SIMULAÇÃO MONTE CARLO: LCG COM CLOCK IRREGULAR" << endl;
    cout << " Amostras: " << num_samples << " | Clocks Iniciais: " << blanking_rounds << endl;
    cout << "========================================================\n" << endl;

    for (int sample = 0; sample < num_samples; ++sample) {
        // Gera um estado interno aleatório (Semente Real)
        PRNGState secret_state;
        secret_state.x = random_GF2E();
        secret_state.a = a;
        secret_state.c = c;

        // Injeta a diferença Delta = 1 no bit 0
        PRNGState dist_state = secret_state;
        GF2E delta_inicial = conv<GF2E>(parsePolynomial("0"));
        dist_state.x = secret_state.x + delta_inicial;

        bool perdeu_sincronismo = false;

        // Roda a fase de "Blanking" (Os primeiros clocks onde não vemos a saída)
        for (int t = 0; t < blanking_rounds; ++t) {
            // Registra o HW da diferença neste instante
            GF2E diferenca_atual = secret_state.x + dist_state.x;
            avg_hw_over_time[t] += getHammingWeight(diferenca_atual);

            // Cada universo calcula o seu próprio clock independentemente!
            bool clock_original = decide_clock(secret_state.x);
            bool clock_modificado = decide_clock(dist_state.x);

            if (clock_original != clock_modificado) {
                perdeu_sincronismo = true;
            }

            // Avança os geradores apenas se o clock deles permitir
            if (clock_original) {
                secret_state.x = secret_state.a * secret_state.x + secret_state.c;
            }
            if (clock_modificado) {
                dist_state.x = dist_state.a * dist_state.x + dist_state.c;
            }
        }

        if (!perdeu_sincronismo) {
            sobrevivencias_perfeitas++;
        }
    }

    // Processar e exibir resultados
    cout << "--- RESULTADOS DA SIMULAÇÃO ---" << endl;
    double taxa_sobrevivencia = (double)sobrevivencias_perfeitas / num_samples * 100.0;
    cout << "> Sobrevivências Perfeitas: " << sobrevivencias_perfeitas << " / " << num_samples << endl;
    cout << "> Probabilidade da Diferença Sobreviver: " << fixed << setprecision(2) << taxa_sobrevivencia << "%" << endl;
    
    cout << "\n--- EVOLUÇÃO DO PESO DE HAMMING (Difusão) ---" << endl;
    cout << "Tempo T | Peso de Hamming Médio da Diferença" << endl;
    cout << "--------------------------------------------" << endl;
    for (int t = 0; t < blanking_rounds; t += 10) { // Imprimindo a cada 10 clocks
        double avg_hw = avg_hw_over_time[t] / num_samples;
        cout << " T = " << setw(3) << t << " | HW = " << fixed << setprecision(2) << avg_hw << " bits contaminados" << endl;
    }
    cout << "========================================================" << endl;

    return 0;
}
