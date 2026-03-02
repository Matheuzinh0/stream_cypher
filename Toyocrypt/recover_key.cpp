#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;
using namespace NTL;

// Funções Auxiliares
GF2X parsePolynomial(const string& s) {
    GF2X result;
    istringstream ss(s);
    long c;
    while (ss >> c) SetCoeff(result, c);
    return result;
}

// Imprime a chave secreta exatamente no formato de texto
void printPolynomialExponents(const GF2X& poly) {
    bool first = true;
    for (long i = deg(poly); i >= 0; --i) {
        if (IsOne(coeff(poly, i))) {
            if (!first) cout << " ";
            cout << i;
            first = false;
        }
    }
    cout << endl;
}

int main() {
    //Parametros da estrutura
    string h_str = "31 13 8 3 0";
    string a_str = "4 0";
    string c_str = "0";
    
    // fluxo convergente dos diferenciais
    string recovered_stream = "1001000110100101111011010000100"; 
    
    long target_bit = 10; // O bit s[k] 

    // Inicialização do GF(q)
    GF2X h = parsePolynomial(h_str);
    GF2E::init(h);
    long N = deg(h); // N é o grau, ou seja, 31.

    GF2E a = conv<GF2E>(parsePolynomial(a_str));
    GF2E c = conv<GF2E>(parsePolynomial(c_str));

    cout << "Sistema Linear " << N << "x" << N << "..." << endl;

    //Motagem da matriz de transiç ão (M) e (Y)
    mat_GF2 M;
    M.SetDims(N, N);
    vec_GF2 Y;
    Y.SetLength(N);

    GF2E C_i = conv<GF2E>(GF2X::zero()); // Contador da constante c, por definição c=0 para classe analisada
    GF2E a_pow = conv<GF2E>(GF2X(0, 1)); // Multiplicador a^i, começa em a^0 = 1

    for (long i = 0; i < N; ++i) {
        // Bit 'B_i' observado no ataque diferencial
        long B_i = recovered_stream[i] - '0'; 
        
        // C_i: Deslocamento caso o LCG tenha 'c' diferente de zero
        long c_bit = IsOne(coeff(rep(C_i), target_bit));
        
        // Lado direito da equação (Vetor Coluna Y)
        Y[i] = B_i ^ c_bit;

        // Monta a Linha 'i' da matriz (Como S_0 influencia o bit S_10 no clock 'i')
        for (long j = 0; j < N; ++j) {
            GF2X xj;
            SetCoeff(xj, j); // Polinômio x^j (Representa 1 bit específico de S_0)
            
            // Avança aquele bit específico no tempo: a^i * x^j
            GF2E term = a_pow * conv<GF2E>(xj);
            
            // Se ele parar no 'target_bit', a contribuição na matriz é 1
            M[i][j] = IsOne(coeff(rep(term), target_bit));
        }

        // Avança o modelo matemático de 'a' e 'c' para o clock seguinte
        C_i = a * C_i + c;
        a_pow = a * a_pow;
    }

    // ======================================================================
    // 3. RESOLVENDO A MATRIZ: S_0 = M^(-1) * Y
    // ======================================================================
    mat_GF2 M_inv;
    try {
        inv(M_inv, M); // Inversão de matriz em GF(2) pela NTL
    } catch (...) {
        cerr << "[ERRO] A matriz de transicao nao e inversivel. Ataque falhou." << endl;
        return 1;
    }

    vec_GF2 S_0 = M_inv * Y;

    // ======================================================================
    // 4. RESULTADO (A Chave Mestra)
    // ======================================================================
    GF2X secret_poly;
    for (long j = 0; j < N; ++j) {
        if (S_0[j] == 1) {
            SetCoeff(secret_poly, j);
        }
    }

    cout << "\n========================================================" << endl;
    cout << " [VITÓRIA ABSOLUTA] CHAVE SECRETA RECUPERADA COM SUCESSO! " << endl;
    cout << "========================================================" << endl;
    cout << "> Forma Polinomial: " << secret_poly << endl;
    cout << "> Chave (Semente):  ";
    printPolynomialExponents(secret_poly);
    cout << "========================================================" << endl;

    return 0;
}