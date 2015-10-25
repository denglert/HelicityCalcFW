all : build 
	@ ./bin/main

test : build
	@ ./bin/main_test

####################

build : build_main build_main_test

build_main : 
	@ cd src; make ../bin/main;

build_main_test : 
	@ cd src; make ../bin/main_test;

clean :
	rm ./lib/*.o
