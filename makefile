default:
	@echo targets: roxy install

roxy:
	R -e "devtools::document()"

install:
	R CMD install .
