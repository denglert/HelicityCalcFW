#ifndef PHYSCONS_H
#define PHYSCONS_H

///////////////////////////////
// --- Physical constans --- //
///////////////////////////////

const double hbar_c  = 0.1973269718; // GeV*fm
const double hbar    = 6.58211927*1E-25; // GeV*s
const double c = 299792458.0; // m/s

const double c_tau_tau = 87.03*1E9; // fm
const double c_tau_muon = 658.6384*1E15; // fm

const double G_Fermi = 1.16*1E-5; // GeV^{-2}

//////////////////////////////////////
// --- Particle masses in [GeV] --- //
//////////////////////////////////////

// -- Scalars -- //
const double m_higgs = 125.0; 

// -- Leptons -- //
//const double m_ele = 0.000511; 
//const double m_muo = 0.10569; 
const double m_ele = 0.0; 
const double m_muo = 0.0; 
const double m_tau = 1.77582; 

// -- Lepton Neutrinos -- //
const double m_nu_ele = 0.0; 
const double m_nu_muo = 0.0; 
const double m_nu_tau = 0.0; 


#endif
