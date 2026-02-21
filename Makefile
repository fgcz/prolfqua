.PHONY: all check check-fast test build build-vignettes document coverage install lint format clean help site deploy renv-init renv-restore renv-snapshot

all: check

help:
	@echo "prolfqua development targets:"
	@echo "  make all       - full pipeline: document -> build -> check (default)"
	@echo "  make check     - R CMD check (runs document, build first)"
	@echo "  make check-fast - R CMD check without rebuilding vignettes during check"
	@echo "  make build-vignettes - build vignettes into inst/doc"
	@echo "  make test      - run testthat tests (runs document first)"
	@echo "  make build     - build tarball (runs document first)"
	@echo "  make document  - generate roxygen2 docs"
	@echo "  make coverage  - code coverage report"
	@echo "  make install   - install package locally"
	@echo "  make lint      - run lintr"
	@echo "  make format    - format with air"
	@echo "  make clean     - remove build artifacts"
	@echo "  make site      - build pkgdown site locally"
	@echo "  make deploy    - build pkgdown site and push to gh-pages"
	@echo ""
	@echo "  Environment (renv):"
	@echo "  make renv-init       - initialize renv and install all deps (first time)"
	@echo "  make renv-restore    - restore environment from renv.lock"
	@echo "  make renv-snapshot   - update renv.lock after installing new packages"

document:
	Rscript -e "devtools::document()"

build: document
	Rscript -e "devtools::build()"

check: build
	Rscript -e "devtools::check()"

build-vignettes: document
	Rscript -e "devtools::build_vignettes()"
	mkdir -p inst/doc
	cp doc/*.html doc/*.Rmd doc/*.R inst/doc/ 2>/dev/null || true

check-fast: document
	Rscript -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-vignettes')"

test: document
	Rscript -e "devtools::test()"

coverage: document
	Rscript -e "covr::package_coverage() |> print()"

install: document
	Rscript -e "devtools::install()"

lint:
	Rscript -e "lintr::lint_package()"

format:
	air format .

site: document
	Rscript -e "pkgdown::build_site()"

deploy: document
	Rscript -e "pkgdown::deploy_to_branch()"

renv-init:
	Rscript -e "renv::init(bioconductor = TRUE)"

renv-restore:
	Rscript -e "renv::restore()"

renv-snapshot:
	Rscript -e "renv::snapshot()"

clean:
	rm -rf *.Rcheck
	rm -f Rplots.pdf
	rm -rf inst/doc doc Meta
