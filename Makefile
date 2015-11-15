all : build 
	@ ./bin/main

test : build
	@ ./bin/main_test

main02 : build_main02
	@ ./bin/main02

####################

build : build_main build_main_test

build_main : 
	@ cd src; make ../bin/main;

build_main_test : 
	@ cd src; make ../bin/main_test;

build_main02 : 
	@ cd src; make ../bin/main02;

clean :
	rm ./lib/*.o
