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

