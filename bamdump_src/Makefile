CFLAGS = -Wall -O3 -std=c++11
INCLUDES = -I$(BAMTOOLS_PATH)include -I../lib/tclap/include
LFLAGS = -L$(BAMTOOLS_PATH)lib -lbamtools

bamdump: bamdump.cpp
	$(CC) -o bamdump bamdump.cpp $(CFLAGS) $(LFLAGS) $(INCLUDES)

clean:
	rm -f bamdump
