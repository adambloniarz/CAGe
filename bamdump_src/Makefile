CFLAGS = -Wall -O3 -I$(BAMTOOLS_PATH)include -std=c++11
LFLAGS = -L$(BAMTOOLS_PATH)lib -lbamtools

bamdump: bamdump.cpp
	$(CC) -o bamdump bamdump.cpp $(CFLAGS) $(LFLAGS)

clean:
	rm -f bamdump
