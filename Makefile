notebooks := $(wildcard ./*.ipynb)

format: $(notebooks)
	black -t py27 $(notebooks)
