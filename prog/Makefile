ResultsDir = "/home/sinyamada/results/kappa_013700.013640/ts64/spin00_raw"
#BSUB = bsub
SRC_SPIN = 0
SINK_SPIN = 0
LOG= ../log/$@
all:
	g++  matrix.cpp  noise_redection.cpp -o NR.out
#	g++ -DNDEBUG  matrix.cpp  noise_redection.cpp -o NR.out
clean:
	rm  ./*.out ./*~ ./core*
runnr:
	$(BSUB) ./NR.out $(ResultsDir) $(SINK_SPIN) $(SRC_SPIN)
