CFLAGS = -I/usr/local/include -L/usr/local/lib
LIBGSL = -lgsl -lgslcblas -lm -O2 -larmadillo
OBJ = data_processing.o decomposition.o mapping_BFS.o mapping_SVD.o fct_utils.o

all : main clean

test : test_sp.cpp
	g++ -Wall -o test_sp test_sp.cpp $(LIBGSL)

main : main.o $(OBJ)
	g++ -Wall -o main main.o lib/*.o $(LIBGSL)

main.o : main.cpp
	g++ -Wall -g -o main.o -c main.cpp

data_processing.o : lib/data_processing.cpp
	g++ -Wall -g -o lib/data_processing.o -c lib/data_processing.cpp

decomposition.o : lib/decomposition.cpp
	g++ -Wall -g -o lib/decomposition.o -c lib/decomposition.cpp

mapping_BFS.o : lib/mapping_BFS.cpp
	g++ -Wall -g -o lib/mapping_BFS.o -c lib/mapping_BFS.cpp

mapping_SVD.o : lib/mapping_SVD.cpp
	g++ -Wall -g -o lib/mapping_SVD.o -c lib/mapping_SVD.cpp

fct_utils.o : lib/fct_utils.cpp
	g++ -Wall -o lib/fct_utils.o -c lib/fct_utils.cpp

clean :
	rm lib/*.o
	rm *.o