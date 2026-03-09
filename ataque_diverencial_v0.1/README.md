# Ataques Diferenciais em Cifra de Fluxo

Implementação e validação de ataques criptanalíticos contra um gerador pseudoaleatório (PRNG) baseado num **LCG sobre GF(2ᵉ)** com função filtro $f(x) = s_{20} \oplus (s_0 \wedge s_{10})$.

> Referência teórica: `docs/Differential Cryptanalysis in Stream Ciphers.pdf`

---

## Estrutura do Projeto

```
ataque_diferencial_v0.1/
├── src/                         # Código-fonte dos ataques
│   ├── batch_attack.cpp         # Ataque diferencial em lote + geração de fluxos
│   ├── recover_key.cpp          # Recuperação de chave via sistema linear (M·S₀ = B)
│   ├── differential_correlation.cpp  # Mede φ(z,s₂₀) e φ(ΔZ,s₁₀)
│   ├── distinguish_attack.cpp   # Distinguidor: dZ=Δ₂₀ quando Δ₀=Δ₁₀=0
│   ├── higher_order_differential.cpp # D²f(x; α, β) = const. quando α·β relevantes
│   ├── multi_difference_attack.cpp   # Múltiplas características → recupera s₀ e s₁₀
│   └── exhaustive_speedup.cpp   # Speedup 2× via simetria de pares {K, K⊕Δ}
├── config/
│   └── unified_config.json      # Parâmetros GF(2ᵉ): h(x), a(x), c(x), seeds
├── docs/
│   └── Differential Cryptanalysis in Stream Ciphers.pdf
├── output/
│   ├── streams/                 # Fluxos de bits gerados por batch_attack
│   └── reports/                 # Relatórios .txt de cada ataque
├── bin/                         # Binários compilados (gerado pelo make)
├── Makefile
└── README.md
```

---

## Dependências

| Biblioteca | Arch/Manjaro | Debian/Ubuntu |
|---|---|---|
| NTL | `sudo pacman -S ntl` | `sudo apt install libntl-dev` |
| GMP | `sudo pacman -S gmp` | `sudo apt install libgmp-dev` |
| nlohmann/json | `sudo pacman -S nlohmann-json` | `sudo apt install nlohmann-json3-dev` |

---

## Compilação e Execução

```bash
# Compila todos os executáveis
make

# Executa os 4 ataques de validação (lêem config/unified_config.json)
make run-validation

# Gera fluxos de bits e recupera chave
make run-batch

# Teste de speedup de busca exaustiva (independente do JSON)
make exhaustive

# Limpa binários
make clean
```

---

## Descrição dos Ataques

### 1. `batch_attack` + `recover_key` — Recuperação de Chave via Álgebra Linear

**Ideia:** injeta diferença $\Delta = 1$ (bit 0) no estado inicial. A linearidade do LCG garante que $\Delta_t = a^t \cdot \Delta_0$. Da função filtro:

$$\Delta Z_t = f(x_t) \oplus f(x_t \oplus \Delta_0) = s_{10}(x_t)$$

Com $N = \deg(h)$ observações de $s_{10}$ monta-se o sistema linear $M \cdot S_0 = B$ e resolve-se por inversão de matriz em GF(2).

---

### 2. `differential_correlation` — Correlação Diferencial (Coeficiente de Phi)

Mede a correlação de Matthews entre:
- $\phi(z, s_{20})$: correlação passiva (viés da saída)
- $\phi(\Delta Z, s_{10})$: correlação diferencial esperada = **1.0**

---

### 3. `distinguish_attack` — Distinguidor Determinístico

Quando $\Delta_0 = \Delta_{10} = 0$ (i.e. $t$ tal que $a^t \cdot \Delta$ anula os bits 0 e 10):

$$\Delta Z_t = \Delta_{20}(t)$$

A taxa de discordância no alvo é **0%**; no aleatório é ~50%.

---

### 4. `higher_order_differential` — Diferencial de 2ª Ordem

Para $\alpha = x^0$, $\beta = x^{10}$ o diferencial de segunda ordem:

$$D^2 f(x; \alpha, \beta) = f(x) \oplus f(x\oplus\alpha) \oplus f(x\oplus\beta) \oplus f(x\oplus\alpha\oplus\beta)$$

é **constante** (independente de $x$), provando que $f$ tem grau algébrico 2.

---

### 5. `multi_difference_attack` — Múltiplas Características

Usa 5 deltas simultâneos $\{x^0, x^1, x^2, x^5, x^{10}\}$ para recuperar $s_0$ e $s_{10}$ diretamente por combinações lineares da diferença de saída.

---

### 6. `exhaustive_speedup` — Speedup na Busca Exaustiva

Demonstra que, dado que os pares $\{K, K \oplus \Delta\}$ produzem os mesmos dois fluxos (em ordens trocadas), a busca exaustiva pode ser reduzida à metade testando apenas representantes pares.

---

## Resultados Esperados

Todos os ataques retornam **CONSISTENTE** em todas as 8 configurações do JSON (campos GF(2³¹) a GF(2⁴⁵¹)).

| Executável | Saída esperada |
|---|---|
| `differential_correlation` | `φ(ΔZ, s₁₀) = 1.0` em todas as configs |
| `distinguish_attack` | `mismatches = 0` no alvo, ~50% no aleatório |
| `higher_order_differential` | `rate = 1.0` em todos os checks |
| `multi_difference_attack` | `equation rate = 1.0` + s₀ e s₁₀ recuperados |
| `exhaustive_speedup` | `speedup = 2.0x`, consistente |
| `batch_attack` + `recover_key` | chave recuperada coincide com semente original |
