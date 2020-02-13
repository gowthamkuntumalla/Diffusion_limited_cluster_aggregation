cluster_in_a_box: cluster_in_a_box.o
	g++ -o cluster_in_a_box cluster_in_a_box.o 
cluster_in_a_box.o: cluster_in_a_box.cpp 
	g++ -c -Wno-deprecated cluster_in_a_box.cpp   

clean:
	rm *.o *~ cluster_in_a_box