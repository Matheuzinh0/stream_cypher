#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <deque>
using namespace std;
using namespace NTL;

GF2X parsePolynomial(const string& s) {
    GF2X result;
    istringstream ss(s);
    long coeff;
    while (ss >> coeff) {
        SetCoeff(result, coeff);
    }
    return result;
}

string GF2EToString(const GF2E& element) {
    GF2X poly = rep(element);
    stringstream ss;
    ss << poly;
    return ss.str();
}

struct PRNGState {
    GF2E x, x1, a, c;
    long fieldDegree;
    deque<bool> bit_buffer;
};

static PRNGState prng_state;

unsigned int prng_generator(PRNGState& s) {
    while (s.bit_buffer.size() < 32) {
        s.x = s.a * s.x + s.c;
        GF2X xPoly = rep(s.x);
        for (long j = 0; j < s.fieldDegree; ++j) {
            bool bit = IsOne(coeff(xPoly, j));
            s.bit_buffer.push_back(bit);
        }
    }

    unsigned int output = 0;
    for (int i = 0; i < 32; ++i) {
        output <<= 1;
        output |= s.bit_buffer.front() ? 1 : 0;
        s.bit_buffer.pop_front();
    }
    return output;
}

double prng_generator_out(PRNGState& s) {
    unsigned int v = prng_generator(s);
    return v / (double)UINT32_MAX;
}
//g++ -O3 -march=znver2 -pg -flto -fopenmp run_tests.cpp -o run -ltestu01 -lprobdist -lntl -lgmp -pthread && ./run
//g++ -O3 -march=znver2 prng.cpp -o prng -lntl -lgmp -lprobdist -ltestu01 -lm

int main(int argc, char* argv[]) {
    
    string h_str = "31 13 8 3 0";
    string a_str = "4 0";
    string c_str = "0";
    string x_str = "31 29 28 25 23 22 21 19 17 14 13 11 10 9 7 5 0";
    GF2X h = parsePolynomial(h_str);
    GF2E::init(h);

    PRNGState original_state;
    original_state.fieldDegree = deg(h);
    original_state.a = conv<GF2E>(parsePolynomial(a_str));
    original_state.c = conv<GF2E>(parsePolynomial(c_str));
    original_state.x = conv<GF2E>(parsePolynomial(x_str));

    PRNGState delta_prng;
    delta_prng.fieldDegree = deg(h);
    delta_prng.a = conv<GF2E>(parsePolynomial(a_str));
    delta_prng.c = conv<GF2E>(parsePolynomial(c_str));
    string x1_str = "31 29 28 25 23 22 21 19 17 14 13 11 10 9 7 5"; // x^0 xor 1, então 1 xor 1 = 0, então só mudamos o bit menos significativo
    delta_prng.x = conv<GF2E>(parsePolynomial(x1_str));
/*
    // Loop sincronizado: ambos avançam um passo por vez
    for (int i = 0; i < 100; ++i) {
        // Avançar estados manualmente (um passo)
        GF2E x_new_orig = original_state.a * original_state.x + original_state.c;
        original_state.x = x_new_orig;

        GF2E x_new_delta = delta_prng.a * delta_prng.x + delta_prng.c;
        delta_prng.x = x_new_delta;me explique então o que eu fiz. Como sei que o que eu fiz no meu código está correto? como posso validar?

        // Calcular deltas
        GF2E delta = original_state.x + delta_prng.x;  // XOR em GF(2)

        cout << "Iteração " << i + 1 << ": \n original x = " << GF2EToString(original_state.x) 
            // << "\n delta x = " << GF2EToString(delta_prng.x) 
             << "\n    delta   = " << GF2EToString(delta_prng.x) << endl;

        if (delta_prng.x == original_state.x) {
            cout << "Estados iguais na iteração " << i + 1 << endl;
            break;
        }
    }
*/
// --- PREPARAÇÃO PARA VALIDAÇÃO ---
    // Delta Teórico começa igual à diferença das sementes (1)
    // O polinômio "0" em NTL representa x^0, que é o valor 1.
    GF2E delta_teorico = conv<GF2E>(parsePolynomial("0")); 
    
    // Precisamos do multiplicador 'a' para a previsão
    GF2E a_poly = original_state.a;

    cout << "--- INICIANDO VALIDAÇÃO DO ATAQUE ---" << endl;

    for (int i = 0; i < 100; ++i) {
        // 1. Avançar os Estados Reais (Simulação da Cifra)
        original_state.x = original_state.a * original_state.x + original_state.c;
        delta_prng.x     = delta_prng.a * delta_prng.x + delta_prng.c;

        // 2. Calcular o Delta Real (O que observamos na saída)
        GF2E delta_real = original_state.x + delta_prng.x;

        // 3. Calcular o Delta Teórico (Previsão Matemática do Atacante)
        // O atacante não sabe o estado 'x', ele só multiplica o delta anterior por 'a'
        delta_teorico = delta_teorico * a_poly;

        // 4. Comparação e Impressão
        cout << "Iteração " << i + 1 << ":" << endl;
        // Corrigido: Agora imprimimos o delta_real, não o estado delta_prng.x
        cout << "   Delta Real (Obs): " << GF2EToString(delta_real) << endl;
        cout << "   Delta Teorico:    " << GF2EToString(delta_teorico) << endl;

        // Validação
        if (delta_real == delta_teorico) {
            cout << "   [OK] Previsão confirmada." << endl;
        } else {
            cout << "   [ERRO] A previsão falhou!" << endl;
            return 1; // Para o programa se houver erro
        }
        cout << "----------------------------------------" << endl;
    }

    cout << "SUCESSO: O sistema é perfeitamente linear." << endl;
    return 0;
}
