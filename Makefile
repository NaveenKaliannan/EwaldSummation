
main:main.cpp ewald.cpp
	g++ main.cpp ewald.cpp -o main 

clean:
	rm -f main
	rm -f output/data/*
	rm -f output/graph/*

