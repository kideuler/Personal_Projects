.SUFFIXES:
FC = python3
main: 
	
.PHONY: clean
clean:
	@rm -rf __pycache__/*.pyc
.PHONY: interior mesh
interior mesh:
	@$(FC) main.py 1
.PHONY: exterior mesh
exterior mesh:
	@$(FC) main.py 0
