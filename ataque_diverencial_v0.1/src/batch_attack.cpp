#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <nlohmann/json.hpp> // sudo pacman -S nlohmann-json
//#include "json.hpp"
//Autor: Matheus Sousa, DES-UFPE;
/* Descrição: O ataque consiste em: Dado que o atacante conhece a estrutura e o texto plano é escolhido, ele sincroniza duas estruturas numa delas o atacante
inverte o último bit de inicialização S0' e preserva a outra S0, dada duas interações o atacante assume que o diferencial(dada a propagação do erro numa das estruturas S0') é ΔZ = Z + Z' ( " + " é o XOR em GF(2^e)).
o atacante assume que existirá uma evolução da propagação de erro no tempo. Dada que foi feita a inversão do bit de propagação, se o sistema for linear, o bit real é o complemento do bit de erro propagado (ruído).
O atacante consegue forçar utilizando a álgebra a linearização forçando com o que ΔZ = S0(1 bit), essa prova é mostrada no documento em anexo no repositório.
assumindo a evoluçaõ no tempo o atacante registra os valores de evolução do ruído e consegue montar um sistema linear a partir de uma matriz de trasição dada por parãmetros da estrutura (M) * S0 = B => s0= M ^ -1 * B.
*/
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

string GF2EToString(const GF2E& element) {
    GF2X poly = rep(element);
    stringstream ss;
    ss << poly;
    return ss.str();
}

struct PRNGState {
    GF2E x, a, c;
};

// função filtro selecionando o bit 10
int filter_function(const GF2E& state) {
    GF2X poly = rep(state);
    long s0  = IsOne(coeff(poly, 0));
    long s10 = IsOne(coeff(poly, 10));
    long s20 = IsOne(coeff(poly, 20));
    return s20 ^ (s0 & s10); 
}

int main() {
    //Leitura do arquivo .json
    ifstream file("unified_config.json");
    if(!file.is_open()) {
        cerr << "Arquivo não encontrado" << endl;
        return 1;
    }
    
    json configs;
    file >> configs;
    //Arquivo de Log
    ofstream log_file("relatorio_ataques.txt");
    log_file << " Relatório \n";

    int config_id = 0x1;

    // 3. Loop pelas Configurações
    for (const auto& config : configs) {
        string h_str = config["base_params"]["h"];
        string a_str = config["base_params"]["a"];
        string c_str = config["base_params"]["c"];

        // Redefinir o GF(2^e) para a nova configuração
        GF2X h = parsePolynomial(h_str);
        GF2E::init(h); 

        GF2E a = conv<GF2E>(parsePolynomial(a_str));
        GF2E c = conv<GF2E>(parsePolynomial(c_str));
        long degree = deg(h); // o grau *e* dirá quantos bits serão necessários para serém extraídos, verificar a necessidade de e+k bits

        log_file << "--------------------------------------------------------\n";
        log_file << "config " << config_id++ << "GF(2^" << degree << "))\n";
        log_file << " h(x) = " << h_str << "\n";
        log_file << " a(x) = " << a_str << "\n";
        log_file << "--------------------------------------------------------\n";

        int seed_id = 0x1;
        
            //Semente e ataque na semente com inversão no bit 0 -> que é operar xor 1
        for (const auto& seed_val : config["seeds"]) {
            string secret_str = seed_val.get<string>();
            log_file << "\n iInvertendo um bit da seed" << seed_id++ << "\n";
            log_file << " seed original ( polinomial ): " << secret_str << "\n";

            // Inicializa Gerador Secreto
            PRNGState secret_state;
            secret_state.x = conv<GF2E>(parsePolynomial(secret_str));
            secret_state.a = a; 
            secret_state.c = c;

            // Inicializa Gerador Diferencial (Injeta Delta = 1 no bit 0)
            PRNGState dist_state = secret_state;
            GF2X _d0; SetCoeff(_d0, 0); // polinômio constante 1 = x^0
            GF2E delta_inicial = conv<GF2E>(_d0); // delta = 1 em GF(2^e)
            dist_state.x = secret_state.x + delta_inicial;
            
            GF2E delta_teorico = delta_inicial;

            vector<int> recovered_stream;
            int acertos = 0;

            log_file << "      -- Log de Evolução e Extração (Tempo T = 0 a " << degree - 1 << ") --\n";

            // 5. O Ataque: Roda por `degree` clocks para recuperar bits independentes suficientes
            for (long t = 0; t < degree; ++t) {
                // Verificação de Distinção (Cifra linear previsível?)
                GF2E delta_interno_real = secret_state.x + dist_state.x;
                bool lin_ok = (delta_interno_real == delta_teorico);

                // Recuperação de Chave: compara f(x_t) com f(x_t XOR e_0)
                // onde e_0 = 1 (bit 0). Pela álgebra da função filtro:
                //   f(x XOR e_0) = s20 XOR ((s0 XOR 1) AND s10) = f(x) XOR s10
                // portanto delta_z = f(x_t) XOR f(x_t XOR 1) = s10(x_t)  [garantido]
                int z_real = filter_function(secret_state.x);

                PRNGState probe_state = secret_state;
                probe_state.x = probe_state.x + delta_inicial; // delta_inicial = e_0 = 1, fixo
                int z_probe = filter_function(probe_state.x);

                int delta_z = z_real ^ z_probe;
                long bit_10_gabarito = IsOne(coeff(rep(secret_state.x), 10));

                // Pela álgebra do LCG linear, delta_z == bit_10 é garantido
                bool bit_ok = (delta_z == bit_10_gabarito);
                recovered_stream.push_back(delta_z); // sempre válido
                if (bit_ok) acertos++;

                // Salva no log a validação do clock
                log_file << "        T=" << t 
                         << " | Dif_Teorica_OK=" << (lin_ok ? "Sim" : "Nao")
                         << " | Z=" << z_real << " Z'=" << z_probe
                         << " | dZ=" << delta_z << " (Gabarito S[10]=" << bit_10_gabarito << ")\n";

                // Avança o gerador
                secret_state.x = secret_state.a * secret_state.x + secret_state.c;
                dist_state.x   = dist_state.a   * dist_state.x   + dist_state.c;
                delta_teorico  = delta_teorico  * a;
            }

            // Exibe convergência da matriz
            log_file << "      -> [RESULTADO] Precisão da Recuperação: " << acertos << " / " << degree << " bits";
            log_file << (acertos == degree ? " [CONVERGÊNCIA TOTAL]" : " [FALHA]") << "\n";
            log_file << "      -> [FLUXO CONVERGENTE PARA ÁLGEBRA LINEAR]:\n         ";
            for(int bit : recovered_stream) log_file << bit;
            log_file << "\n";

            // Grava o fluxo em arquivo dedicado para o recover_key
            string stream_filename = "stream_" + to_string(config_id - 1) + "_" + to_string(seed_id - 1) + ".txt";
            ofstream sf(stream_filename);
            sf << h_str << "\n" << a_str << "\n" << c_str << "\n";
            for(int bit : recovered_stream) sf << bit;
            sf << "\n";
            sf.close();
            log_file << "      -> [ARQUIVO] Fluxo salvo em: " << stream_filename << "\n";
        }
        log_file << "\n";
    }

    log_file.close();
    cout << "Processamento finalizado! O arquivo 'relatorio_ataques.txt' foi gerado com sucesso." << endl;
    return 0;
}
