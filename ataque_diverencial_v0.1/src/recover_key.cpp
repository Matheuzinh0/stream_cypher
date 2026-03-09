#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <iostream>
#include <fstream>
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

// Imprime os expoentes do polinômio no formato "e1 e2 e3 ..."
string polynomialToExponents(const GF2X& poly) {
    ostringstream oss;
    bool first = true;
    for (long i = deg(poly); i >= 0; --i) {
        if (IsOne(coeff(poly, i))) {
            if (!first) oss << " ";
            oss << i;
            first = false;
        }
    }
    return oss.str();
}

int main(int argc, char* argv[]) {
    // Arquivo de fluxo gerado pelo batch_attack (padrão: stream_1_1.txt)
    string stream_file = (argc > 1) ? argv[1] : "stream_1_1.txt";

    // Lê o arquivo: linha 1=h, linha 2=a, linha 3=c, linha 4=fluxo de bits
    ifstream sf(stream_file);
    if (!sf.is_open()) {
        cerr << "[ERRO] Arquivo '" << stream_file << "' nao encontrado.\n"
             << "       Execute batch_attack primeiro para gerar os arquivos de fluxo.\n";
        return 1;
    }
    string h_str, a_str, c_str, stream_bits;
    getline(sf, h_str);
    getline(sf, a_str);
    getline(sf, c_str);
    getline(sf, stream_bits);
    sf.close();

    long target_bit = 10; // bit s[k] extraído pelo ataque diferencial

    // Inicialização do GF(2^e)
    GF2X h = parsePolynomial(h_str);
    GF2E::init(h);
    long N = deg(h);

    if ((long)stream_bits.size() < N) {
        cerr << "[ERRO] Fluxo com " << stream_bits.size() << " bits, esperado >= " << N << ".\n";
        return 1;
    }

    GF2E a = conv<GF2E>(parsePolynomial(a_str));
    GF2E c = conv<GF2E>(parsePolynomial(c_str));

    cout << "GF(2^" << N << ") | Sistema Linear " << N << "x" << N << "...\n";

    // ======================================================================
    // 1. MONTAGEM DA MATRIZ M e VETOR Y
    //    M[i][j] = bit target_bit de (a^i * x^j)
    //    Y[i]    = B_i XOR bit target_bit de C_i
    // ======================================================================
    mat_GF2 M;
    M.SetDims(N, N);
    vec_GF2 Y;
    Y.SetLength(N);

    GF2E C_i;                              // acumulador de c: C_0=0, C_i = a*C_{i-1}+c
    GF2X one_poly; SetCoeff(one_poly, 0);  // polinômio 1 = x^0
    GF2E a_pow = conv<GF2E>(one_poly);     // a^0 = 1 (correto, sem ambiguidade)

    for (long i = 0; i < N; ++i) {
        long B_i   = stream_bits[i] - '0';
        long c_bit = IsOne(coeff(rep(C_i), target_bit));
        Y[i] = B_i ^ c_bit;

        for (long j = 0; j < N; ++j) {
            GF2X xj; SetCoeff(xj, j);          // monômio x^j
            GF2E term = a_pow * conv<GF2E>(xj); // a^i * x^j em GF(2^e)
            M[i][j] = IsOne(coeff(rep(term), target_bit));
        }

        C_i   = a * C_i + c;
        a_pow = a * a_pow;
    }

    // ======================================================================
    // 2. RESOLUÇÃO: S_0 = M^(-1) * Y
    // ======================================================================
    mat_GF2 M_inv;
    try {
        inv(M_inv, M);
    } catch (...) {
        cerr << "[ERRO] A matriz de transicao nao e inversivel. Verifique o fluxo.\n";
        return 1;
    }

    vec_GF2 S_0_vec = M_inv * Y;

    // ======================================================================
    // 3. RESULTADO
    // ======================================================================
    GF2X recovered_poly;
    for (long j = 0; j < N; ++j)
        if (S_0_vec[j] == 1) SetCoeff(recovered_poly, j);

    string recovered_str = polynomialToExponents(recovered_poly);

    cout << "\n========================================================\n";
    cout << " [CHAVE RECUPERADA]\n";
    cout << "========================================================\n";
    cout << "> Arquivo de fluxo : " << stream_file << "\n";
    cout << "> GF(2^" << N << ") | h(x) = " << h_str << "\n";
    cout << "> Chave recuperada : " << recovered_str << "\n";
    cout << "  (Nota: valor reduzido mod h(x), que e o estado real do gerador.\n";
    cout << "   Sementes do JSON com grau = e sao reduzidas pela NTL ao inicializar.)\n";
    cout << "========================================================\n";

    return 0;
}