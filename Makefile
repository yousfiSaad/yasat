
main: build/yasat

clean: 
	rm build/* build/objects/*

build/objects/CDCL_solver.o: ./src/implementations/CDCL_solver.cpp
	g++ -O3 -c -o $@ $<

build/objects/main.o: ./src/main.cpp
	g++ -O3 -c -o $@ $<

build/yasat: ./build/objects/main.o ./build/objects/CDCL_solver.o
	g++ -o $@ $^

