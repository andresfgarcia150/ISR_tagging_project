/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

Other functions
*/

#include "ROOTFunctions.h"
const Double_t PI = TMath::Pi();

/*
Function that calculates the angular difference between two angles
*/
Double_t deltaAng(Double_t phi1, Double_t phi2)
{
	Double_t DeltaAng =TMath::Abs(phi1-phi2);
	if (DeltaAng > PI) DeltaAng = 2*PI - DeltaAng;
	return DeltaAng;
}
