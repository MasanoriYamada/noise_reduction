ResultsDir = "/work/data2/yamada/OmgOmg/Set1/ts32/spin00"
BSUB = bsub
SRC_SPIN = 0
SINK_SPIN = 0
#LOG= -o ../log/$@
all:
#	g++  matrix.cpp  noise_redection.cpp -o NR.out
	g++ -DNDEBUG  matrix.cpp  noise_redection.cpp -o NR.out
clean:
	rm  ./*.out ./*~ ./core*
runnr:
	$(BSUB) $(LOG) ./NR.out $(ResultsDir) $(SINK_SPIN) $(SRC_SPIN)
