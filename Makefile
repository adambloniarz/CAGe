export CC = g++
export BAMTOOLS_PATH = <PATH TO BAMTOOLS>

all:
	$(MAKE) -C cage_src
	$(MAKE) -C bamdump_src
	mkdir -p bin
	mv cage_src/cage bin
	mv bamdump_src/bamdump bin

clean:
	$(MAKE) clean -C cage_src
	$(MAKE) clean -C bamdump_src
	rm bin/cage
	rm bin/bamdump
	rmdir bin 
