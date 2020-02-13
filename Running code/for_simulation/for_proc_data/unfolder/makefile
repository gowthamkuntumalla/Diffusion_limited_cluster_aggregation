unfolder: unfolder.o
	g++ -o unfolder unfolder.o 
unfolder.o: unfolder.cpp 
	g++ -c -Wno-deprecated unfolder.cpp   

clean:
	rm *.o *~ unfolder