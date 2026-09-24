.PHONY: setup data demo test

setup:
	conda env create -f environment.yml

data:
	bash scripts/download_pbmc_data.sh

demo:
	python scripts/demo.py "$(QUESTION)"

test:
	pytest tests/ -v
