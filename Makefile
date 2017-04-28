exec:
	gcc shock.c -o shock.out -lm
	gcc sedov.c -o sedov.out -lm
shock:
	./shock.out
sedov:
	./sedov.out
plotshock:
	python shock_plot.py
plotsedov:
	python sedov_plot.py
clean: 
	rm -f density.txt pressure.txt velocity.txt
	rm -f r10.txt r60.txt r120.txt t10.txt t60.txt t120.txt
	rm -f sedov.out
	rm -f shock.out