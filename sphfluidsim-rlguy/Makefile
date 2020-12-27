CC=g++
CFLAGS=
LIBS=
OBJ =  \
     gridcell.o \
     cellhash.o \
     main.o \
     quaternion.o \
     spatialgrid.o \
     sphfluidsimulation.o \
     utils.o
     
%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

sphfluidsim.exe: $(OBJ)
	$(CC) -o $@ $^  $(LIBS) 

clean:
	rm -f *.o 
	rm -f *.exe
