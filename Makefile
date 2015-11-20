all : build 

test_rh_tautau : build 
	@ ./bin/test_rh_tautau

test : build
	@ ./bin/main_test

test_rh_6f : build
	@ ./bin/test_rh_6f

cuba : build
	@ ./bin/main_cuba_test

####################

build : build_allbinaries

build_allbinaries : 
	@ cd src; make binaries;

clean :
	rm -f ./lib/*.o
	find ./bin/ -type f -not -name 'dummy' | xargs rm -f
	rm -f ./src/.depend_cpp
	touch ./src/.depend_cpp
