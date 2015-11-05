#include "UtilFunctions.h"

void displayTLorentzVector(TLorentzVector *v)
{
	std::cout << Form("(E,px,py,pz): (%6.2f, %6.2f, %6.2f, %6.2f)\n", v->E(), v->Px(), v->Py(), v->Pz() );
};