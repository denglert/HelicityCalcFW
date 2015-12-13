all : build 

test_rh_tautau : build 
	@ ./bin/test_rh_tautau | tee ./log/test_rh_tautau

test : build
	@ ./bin/main_test

test_rh_6f_k0_1001 : build
	@ ./bin/test_rh_6f_k0_1001 | tee ./log/test_rh_6f_k0_1001

test_rh_6f_k0_physical : build
	@ ./bin/test_rh_6f_k0_physical | tee ./log/test_rh_6f_k0_physical

test_HelicityTools : build
	@ ./bin/test_HelicityTools

test_FortranUtils : build
	@ ./bin/test_FortranUtils

test_PhaseSpaceTools : build
	@ ./bin/test_PhaseSpaceTools

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
