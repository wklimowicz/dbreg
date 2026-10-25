export OPENBLAS_NUM_THREADS ?= 1
export OMP_NUM_THREADS ?= 2

.PHONY: test test-databricks test-spark check document install website help

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

test: ## Run tinytest::test_all()
	Rscript -e "pkgload::load_all('.'); tinytest::test_all('.')"

test-databricks: ## Run native synthetic Databricks tests (requires configured credentials/driver)
	Rscript .github/scripts/run-databricks-tests.R

test-spark: ## Run small real Spark 4.2.0 syntax/algorithm checks (requires sparklyr)
	Rscript .github/scripts/run-spark-tests.R

testall: ## Run all tests including NYC taxi and binsreg
	DBREG_TEST_NYC=TRUE Rscript -e "pkgload::load_all('.'); tinytest::test_all('.')"

check: ## Full R CMD check (no manual)
	R CMD build . && R CMD check --no-manual *.tar.gz

document: ## Regenerate man pages
	Rscript -e "devtools::document('.')"

install: ## Install package locally
	R CMD INSTALL .

website: ## Build docs website (altdoc)
	Rscript -e "altdoc::render_docs(freeze = TRUE)"
