CFLAGS = -Wall -O3 -std=c++11
LFLAGS = -lsqlite3
INCLUDES = -I../lib/tclap/include
OBJECTS = cage.o cost.o dataloader.o variantcaller.o alignment.o suffstats.o pileslist.o snpdb.o pelt.o

cage: $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(CFLAGS) $(LFLAGS) 

-include $(OBJECTS:.o=.d)

%.o: %.cpp
	$(CC) -c $< $(CFLAGS) $(INCLUDES)
	$(CC) -MM $< $(CFLAGS) $(INCLUDES) > $*.d

clean:
	rm -f *.o *.d cage
