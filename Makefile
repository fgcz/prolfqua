.DEFAULT_GOAL := check

PKG_NAME := $(shell awk '/^Package:/ {print $$2; exit}' DESCRIPTION)
PKG_VERSION := $(shell awk '/^Version:/ {print $$2; exit}' DESCRIPTION)
TARBALL := ../$(PKG_NAME)_$(PKG_VERSION).tar.gz

# Load package-specific command overrides and targets before defining the
# shared rules. Every package commits this file, even when it is empty.
include make/package.mk

DOCUMENT_CMD ?= Rscript -e "devtools::document()"
BUILD_CMD ?= Rscript -e "devtools::build(vignettes = TRUE)"
CHECK_CMD ?= Rscript -e "devtools::check()"
CHECK_FAST_CMD ?= Rscript -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-vignettes', vignettes = FALSE)"
CHECK_BIOC_CMD ?= Rscript -e "BiocCheck::BiocCheck()"
BUILD_VIGNETTES_CMD ?= Rscript -e "devtools::build_vignettes()"
TEST_CMD ?= Rscript -e "devtools::test()"
COVERAGE_CMD ?= Rscript -e "covr::package_coverage() |> print()"
INSTALL_CMD ?= R CMD INSTALL $(TARBALL)
INSTALL_BUILD_TARGET ?= build
LINT_CMD ?= Rscript -e "lints <- lintr::lint_package(); print(lints); if (length(lints) > 0L) quit(status = 1L)"
SITE_CMD ?= Rscript -e "pkgdown::build_site(install = FALSE)"
NEW_VERSION_CMD ?= Rscript -e "d <- read.dcf('DESCRIPTION'); old <- d[1, 'Version']; parts <- as.integer(strsplit(old, '.', fixed = TRUE)[[1]]); if (length(parts) < 3) parts <- c(parts, rep(0L, 3L - length(parts))); parts[3] <- parts[3] + 1L; new <- paste(parts, collapse = '.'); x <- readLines('DESCRIPTION'); x <- sub('^Version: .*', paste0('Version: ', new), x); writeLines(x, 'DESCRIPTION'); cat(new)"

.PHONY: all help help-common help-package hooks precommit-prepare document build build-vignettes install test check-fast check-bioc check coverage lint format site clean clean-package new-version vignette

all: check

help: help-common help-package

help-common:
	@echo "$(PKG_NAME) development targets:"
	@echo ""
	@echo "Common:"
	@echo "  make install          - build and install the package"
	@echo "  make test             - run testthat tests"
	@echo "  make check            - run the full R CMD check"
	@echo "  make lint             - run lintr as a hard gate"
	@echo "  make format           - format with Air"
	@echo "  make site             - build the documentation site"
	@echo "  make clean            - remove generated artifacts"
	@echo ""
	@echo "Advanced:"
	@echo "  make document         - generate roxygen2 documentation"
	@echo "  make build            - build the source package with vignettes"
	@echo "  make build-vignettes  - build vignettes into inst/doc"
	@echo "  make check-fast       - run R CMD check without vignettes"
	@echo "  make check-bioc       - run BiocCheck"
	@echo "  make coverage         - generate a code coverage report"
	@echo "  make vignette V=Name  - render one R Markdown vignette"
	@echo "  make hooks            - activate the repository pre-commit hook"
	@echo "  make new-version      - bump, commit, tag, and push a patch release"

help-package:

hooks:
	git config core.hooksPath .githooks
	@echo "Installed pre-commit hook (core.hooksPath = .githooks)."

precommit-prepare:

document:
	$(DOCUMENT_CMD)

build: document
	$(BUILD_CMD)

build-vignettes: document
	rm -rf doc inst/doc
	$(BUILD_VIGNETTES_CMD)
	mkdir -p inst/doc
	cp doc/*.html doc/*.Rmd doc/*.R inst/doc/ 2>/dev/null || true

install: $(INSTALL_BUILD_TARGET)
	$(INSTALL_CMD)

test: document
	$(TEST_CMD)

check-fast: document
	$(CHECK_FAST_CMD)

check-bioc:
	$(CHECK_BIOC_CMD)

check: build
	$(CHECK_CMD)

coverage: document
	$(COVERAGE_CMD)

lint:
	$(LINT_CMD)

format:
	air format .

site: install
	$(SITE_CMD)

vignette:
ifndef V
	$(error Usage: make vignette V=<vignette_name>, e.g. make vignette V=Benchmark_prolfqua)
endif
	Rscript -e "rmarkdown::render('vignettes/$(V).Rmd')"

new-version:
	@NEW_VERSION="$$( $(NEW_VERSION_CMD) )"; \
	echo "Bumped version to $$NEW_VERSION"; \
	git add DESCRIPTION; \
	git commit -m "new version $$NEW_VERSION"; \
	git tag "$$NEW_VERSION"; \
	git push && git push --tags; \
	echo "Released $$NEW_VERSION"

clean: clean-package
	rm -rf *.Rcheck
	rm -rf *.BiocCheck
	rm -f Rplots.pdf
	rm -rf inst/doc doc Meta
	rm -f vignettes/*.html vignettes/*.R

clean-package:
