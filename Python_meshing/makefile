.SUFFIXES:
FC = python3
main: 
	@$(FC) main.py $(deg)
.PHONY: clean
clean:
	@rm -rf __pycache__/*.pyc
.PHONY: meshpoints
meshpoints:
	@$(FC) main_mesh.py $(deg) $(n) 1
.PHONY: meshcurve
meshcurve:
	@$(FC) main_mesh.py $(deg) $(n) 2
.PHONY: test
test:
	@$(FC) main_mesh.py 3 50 1