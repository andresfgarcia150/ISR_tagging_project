// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

/*
-------     Universidad de los Andes      -------
-------      Departamento de Fisica       -------
-------   Proyecto Joven Investigador     -------
-------  Andres Felipe Garcia Albarracin  -------
-------    Juan Carlos Sanabria Arenas    -------

This code develops pythia hadronization. Takes as
parameter a .cmnd file, where a .lhe file from MadGraph
and other parameters are specified. Then the code
produces .hep files after making the hadronization

Obs: The class MyUserHooks is written in order to
veto all the ISR emissions produced after the 
first ISR parton. It is an extension of the code
hadronization01

run as ./hadronization02 input.cmnd [output.hep]

The MakeFile has been also modified to compile 
this file
*/

#include "Pythia8/Pythia.h"
#include "stdhep.h"
#include "stdcnt.h"
#include "stdhep_mcfio.h"
#include <string.h>

using namespace Pythia8;
void fill_stdhep(int i, Event &e);

// Write own derived UserHooks class.

class MyUserHooks : public UserHooks {

public:

	// Constructor.
	MyUserHooks() { }

	// Destructor.
	~MyUserHooks() { }

	// Allow a veto of ISR emissions
	virtual bool canVetoISREmission(){
		return true; 	// Interrupts the initial shower emission after each emission 
				// and allow the emission to be vetoed by the next method.
	}
	
	// Analize each emissionand asks for the number of the ISR emissions so far, in order 
	// to allow just 1 ISR parton per event
	virtual bool doVetoISREmission(int sizeOld, const Event& event, int iSys){
		// counts the number of ISR partons (i.e. the numer of particles with status 43)
		int ISR_part = 0;
		for( int i = 0; i < event.size(); i++){
			if (event[i].status() == 43 || event[i].status() == -43)
				ISR_part ++;
		}
		if (ISR_part > 1)
			return true;
		else
			return false;
	}
};

//====================================================


int main(int argc, char** argv) {
	
	// Interface for conversion from Pythia8::Event to HepMC event.
	char fileout[500], title[100];
	strcpy(title,"output_pythia8\0");	

        // Set up generation.
	// Declare Pythia object
        Pythia pythia;

        // Set simulation configurations. Read the file as parameter. If none, it reads hadro_input.cmnd
        if (argc > 1 ) pythia.readFile(argv[1]);
        else {
		cout << "ERROR: \n No  parameters file has passed as parameter. Abort " << endl;
		return 1;
	}

	// Specify the name of the output file
	if (argc > 2 ) strcpy(fileout,argv[2]);
	else strcpy(fileout,"output_pythia8.hep\0");

	// Especify the number of events
        int nEvent = pythia.mode("Main:numberOfEvents"); // For reading only
	int nAbort = 10; // Maximum number of failures accepted
	int iAbort = 0; // Abortions counter

	// Necessary stdhep functions
	int istr(0);
	int ierr = StdHepXdrWriteOpen(fileout, title, nEvent, istr);
	
	// Set up to do a user veto and send it in.
	MyUserHooks* myUserHooks = new MyUserHooks();
	pythia.setUserHooksPtr( myUserHooks);
	
	// Initialize simulation
	pythia.init();

	// Begin event loop; generate until none left in input file.
	for (int iEvent = 0; iEvent < nEvent ; ++iEvent) {
		// Generate events, and check whether generation failed.
		if (!pythia.next()) {
			// If failure because reached end of file then exit event loop.
			if (pythia.info.atEndOfFile()) break;
			// First few failures write off as "acceptable" errors, then quit.
			if (++iAbort < nAbort) continue;
			break;
		}

		// Fill stdhep file
		fill_stdhep(iEvent+1,pythia.event);
		ierr = StdHepXdrWrite(1,istr);		
	}
	
	StdHepXdrEnd(istr);
	pythia.stat();
	cout << ierr;
	delete myUserHooks;
	return 0;

}

// This functions writes in stdhep format. It was written by Steve Mrenna
void fill_stdhep(int i, Event &e)
{
   int num = e.size();
   hepevt_.nevhep = i;
   hepevt_.nhep = num;
   for (int j = 0; j < num; j++) {
      hepevt_.idhep[j] = e[j].id();
      hepevt_.isthep[j] = e[j].statusHepMC();
      hepevt_.jmohep[j][0] = (e[j].mother1()>0) ? e[j].mother1()+1 : 0;
      hepevt_.jmohep[j][1] = (e[j].mother2()>0) ? e[j].mother2()+1 : 0;
      hepevt_.jdahep[j][0] = (e[j].daughter1()>0) ? e[j].daughter1()+1 : 0;
      hepevt_.jdahep[j][1] = (e[j].daughter2()>0) ? e[j].daughter2()+1 : 0;
      hepevt_.phep[j][0] = e[j].px();
      hepevt_.phep[j][1] = e[j].py();
      hepevt_.phep[j][2] = e[j].pz();
      hepevt_.phep[j][3] = e[j].e();
      hepevt_.phep[j][4] = e[j].m();
      hepevt_.vhep[j][0] = e[j].xProd();
      hepevt_.vhep[j][1] = e[j].yProd();
      hepevt_.vhep[j][2] = e[j].zProd();
      hepevt_.vhep[j][3] = e[j].tProd();
   }
}

