CC=g++
CFLAGS=
LIBS=
OBJ = main.o 
     
%.o: %.cpp
	$(CC) -c -o $@ $< 

mueller3d.exe: main.o
	$(CC) -o $@ $^  $(LIBS) 

clean:
	rm -f *.o 
	rm -f *.exe
