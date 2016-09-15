all: graph.o maxflow.o EasyBMP.o
	echo 'building maxflow lib and easybmp lib'
graph.o:
	g++ -O2 -c maxflow/graph.cpp
maxflow.o:
	g++ -O2 -c maxflow/maxflow.cpp
EasyBMP.o:
	g++ -O2 -c EasyBMP/EasyBMP.cpp
clean:
	rm *.o
