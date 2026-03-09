# =============================================================================
# Ataques Diferenciais em Cifra de Fluxo LCG sobre GF(2^e)
# Dependências: NTL, GMP, nlohmann/json
#   Arch/Manjaro : sudo pacman -S ntl gmp nlohmann-json
#   Debian/Ubuntu: sudo apt install libntl-dev libgmp-dev nlohmann-json3-dev
# =============================================================================

CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall
LIBS     := -lntl -lgmp

SRCDIR   := src
BINDIR   := bin
CFGDIR   := config
OUTDIR   := output

# Targets que lêem unified_config.json via caminho relativo
JSON_READERS := differential_correlation distinguish_attack \
                higher_order_differential higher_order_stress \
                higher_order_scaling \
                multi_difference_attack batch_attack

# Targets independentes do JSON
STANDALONE := exhaustive_speedup recover_key

ALL_TARGETS := $(JSON_READERS) $(STANDALONE)
BINS := $(patsubst %, $(BINDIR)/%, $(ALL_TARGETS))

.PHONY: all clean run-all run-batch run-validation

all: $(BINS)

$(BINDIR)/%: $(SRCDIR)/%.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@

$(BINDIR):
	mkdir -p $(BINDIR)

# Executa todos os ataques de validação (lêem config do diretório atual via symlink/cópia)
run-validation: all
	@echo ">>> [1/4] Differential Correlation"
	cd $(CFGDIR) && ../$(BINDIR)/differential_correlation
	@echo ">>> [2/4] Distinguishing Attack"
	cd $(CFGDIR) && ../$(BINDIR)/distinguish_attack
	@echo ">>> [3/4] Higher-Order Differential"
	cd $(CFGDIR) && ../$(BINDIR)/higher_order_differential
	@echo ">>> [4/4] Multi-Difference Attack"
	cd $(CFGDIR) && ../$(BINDIR)/multi_difference_attack

stress-higher-order: all
	@echo ">>> Stress Test: Higher-Order Differential vs. alto grau algebrico"
	cd $(CFGDIR) && ../$(BINDIR)/higher_order_stress
	@mv $(CFGDIR)/relatorio_higher_order_stress.txt $(OUTDIR)/reports/ 2>/dev/null || true

scaling-higher-order: all
	@echo ">>> Scaling: ordem mínima para quebrar grau d (custo 2^k)"
	cd $(CFGDIR) && ../$(BINDIR)/higher_order_scaling
	@mv $(CFGDIR)/relatorio_higher_order_scaling.txt $(OUTDIR)/reports/ 2>/dev/null || true

# Gera os fluxos de bits via batch_attack e recupera a chave do primeiro stream
run-batch: all
	cd $(CFGDIR) && ../$(BINDIR)/batch_attack
	@mkdir -p $(OUTDIR)/streams
	@mv $(CFGDIR)/stream_*.txt $(OUTDIR)/streams/ 2>/dev/null || true
	@mv $(CFGDIR)/relatorio_ataques.txt $(OUTDIR)/reports/ 2>/dev/null || true
	@echo ">>> Recuperando chave de output/streams/stream_1_1.txt"
	$(BINDIR)/recover_key $(OUTDIR)/streams/stream_1_1.txt

run-all: run-validation run-batch

exhaustive:
	$(BINDIR)/exhaustive_speedup

clean:
	rm -rf $(BINDIR)
	rm -f $(CFGDIR)/stream_*.txt $(CFGDIR)/relatorio_*.txt
