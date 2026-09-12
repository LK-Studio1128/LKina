# LKina — one-command build, test, and paper reproduction
#
#   make build         compile LKina locally (Linux/macOS)
#   make test          17-item automated regression suite
#   make docker-build  build the reproducible container
#   make docker-test   run the regression suite inside the container
#   make paper         regenerate paper tables/figures from archived benchmark JSONs
#   make clean         remove local build artifacts
#
# The full raw-PDB cohort pipelines (110-mode coverage, 104-system redocking,
# reactive/covalent suites, parameter sweeps) are single scripts:
#   benchmarks/run_metal_coverage_v3.py
#   benchmarks/metallocomplex_redocking_benchmark.py
#   benchmarks/run_reactive_tests.py  /  run_covalent_full.py
#   benchmarks/run_parameter_sweeps.py
# Run them (hours of CPU) before `make paper` to regenerate from raw inputs.

.DEFAULT_GOAL := help
SHELL := /bin/bash
REPO  := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
BIN   := $(REPO)/build/linux/release/LKina
ifeq ($(shell uname),Darwin)
BIN   := $(REPO)/build/mac/release/LKina
endif

.PHONY: help build test docker-build docker-test paper clean

help:
	@grep -E '^#   make' $(firstword $(MAKEFILE_LIST)) | sed 's/^#   //'

build:
	@cd $(REPO) && ./build_LKina_linux.sh

test:
	@cd $(REPO) && BIN=$(BIN) ./tests/reactive_regression.sh

docker-build:
	docker build -t lkina $(REPO)

docker-test:
	docker run --rm --entrypoint /src/tests/reactive_regression.sh lkina

# Regenerate CSV exports and all paper figures from the archived JSON results.
paper:
	@cd $(REPO)/benchmarks && python3 export_benchmark_data.py
	@cd $(REPO)/benchmarks && python3 make_supplementary_figures.py
	@cd $(REPO)/paper && for s in make_figures_v3.py make_redock_figures.py make_figS3S4.py make_graphical_abstract.py; do \
		if [ -f $$s ]; then echo "== $$s =="; python3 $$s; fi; done
	@echo "paper assets regenerated (figures/ and benchmarks/*.csv)"

clean:
	rm -rf $(REPO)/build/linux/release $(REPO)/build/mac/release
