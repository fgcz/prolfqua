.DEFAULT_GOAL := check

.PHONY: check check-fast build build-vignettes document coverage install lint format clean help site deploy

DOCUMENT_CMD = Rscript -e "devtools::document()"
BUILD_CMD = Rscript -e "devtools::build()"
CHECK_CMD = Rscript -e "devtools::check()"
CHECK_FAST_CMD = Rscript -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-vignettes')"
BUILD_VIGNETTES_CMD = Rscript -e "devtools::build_vignettes()"
TEST_CMD = Rscript -e "devtools::test()"
COVERAGE_CMD = Rscript -e "covr::package_coverage() |> print()"
INSTALL_CMD = Rscript -e "devtools::install()"
LINT_CMD = Rscript -e "lintr::lint_package()"
SITE_CMD = Rscript -e "pkgdown::build_site()"
DEPLOY_CMD = Rscript -e "pkgdown::deploy_to_branch()"

help:
	@echo "prolfqua development targets:"
	@echo "  make check     - R CMD check (runs document, build first)"
	@echo "  make check-fast - R CMD check without rebuilding vignettes during check"
	@echo "  make install   - install package locally"
	@echo "  make lint      - run lintr"
	@echo "  make format    - format package with air"
	@echo "  make clean     - remove build artifacts"
	@echo ""
	@echo "Advanced:"
	@echo "  make document  - generate roxygen2 docs"
	@echo "  make build     - build tarball"
	@echo "  make build-vignettes - build vignettes into inst/doc"
	@echo "  make coverage  - code coverage report"
	@echo "  make site      - build pkgdown site locally"
	@echo "  make deploy    - build pkgdown site and push to gh-pages"
document:
	$(DOCUMENT_CMD)

build: document
	$(BUILD_CMD)

check: build
	$(CHECK_CMD)

build-vignettes: document
	$(BUILD_VIGNETTES_CMD)
	mkdir -p inst/doc
	cp doc/*.html doc/*.Rmd doc/*.R inst/doc/ 2>/dev/null || true

check-fast: document
	$(CHECK_FAST_CMD)

coverage: document
	$(COVERAGE_CMD)

install: document
	$(INSTALL_CMD)

lint:
	$(LINT_CMD)

format:
	air format .

site: document
	$(SITE_CMD)

deploy: document
	$(DEPLOY_CMD)

clean:
	rm -rf *.Rcheck
	rm -f Rplots.pdf
	rm -rf inst/doc doc Meta
