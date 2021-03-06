all : build 
test_rh_tautau : build @ ./bin/test_rh_tautau | tee ./log/test_rh_tautau

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

test_PhaseSpaceTools2 : build
	@ ./bin/test_PhaseSpaceTools2

test_PhaseSpaceTools3 : build
	@ ./bin/test_PhaseSpaceTools3

test_cuba : build
	@ ./bin/test_cuba

LifeTime_Muon_Polarized : build
	@ ./bin/LifeTime_Muon_Polarized | tee ./log/LifeTime_Muon_Polarized

SpinPolarization : build
	@ ./bin/SpinPolarization | tee ./log/SpinPolarization

test_PhaseSpaceIntegration : build
	@ ./bin/test_PhaseSpaceIntegration | tee ./log/test_PhaseSpaceIntegration

test_PhaseSpaceIntegration2 : build
	@ ./bin/test_PhaseSpaceIntegration2 | tee ./log/test_PhaseSpaceIntegration2

test_PhaseSpaceIntegration3 : build
	@ ./bin/test_PhaseSpaceIntegration3 | tee ./log/test_PhaseSpaceIntegration3

test_PhaseSpaceIntegration4 : build
	@ ./bin/test_PhaseSpaceIntegration4 | tee ./log/test_PhaseSpaceIntegration4

Scalar_2Fermion_6Fermion_SingleBranchIntegration : build
	@ ./bin/Scalar_2Fermion_6Fermion_SingleBranchIntegration | tee ./log/Scalar_2Fermion_6Fermion_SingleBranchIntegration

Scalar_2Fermion_6Fermion_FullIntegration : build
	@ ./bin/Scalar_2Fermion_6Fermion_FullIntegration | tee ./log/Scalar_2Fermion_6Fermion_FullIntegration

Analyzer-h-tautau-6f : build
	@ ./bin/Analyzer-h-tautau-6f | tee ./log/Analyzer-h-tautau-6f


####################

build : build_allbinaries

build_allbinaries : 
	@ cd src; make binaries;

clean :
	rm -f ./lib/*.o
	find ./bin/ -type f -not -name 'dummy' | xargs rm -f
	rm -f ./src/.depend_cpp
	touch ./src/.depend_cpp
