/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

*/

#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"

using namespace std;
// Global Variables
const Double_t PI = TMath::Pi();


int main(int argc, char **argv){
	std::cout.precision(4);
	// Counting time
	Double_t initialTime = clock();

	// Create chains of root trees
	TChain chain_Pythia("STDHEP");
	TChain chain_Delphes("Delphes");

	// Loading simulations from Delphes
	Char_t *local_path;
	local_path = (Char_t*) malloc(512*sizeof(Char_t));
	strcpy(local_path,"/home/af.garcia1214/PhenoMCsamples/Simulations/MG_pythia8_delphes_parallel/"); // At the server

	Char_t head_folder[] = "_Tops_Events_WI/";

	Char_t current_folder[] = "_Tops_MG_1K_AG_WI_003/";

	Char_t unidad = '3'; Char_t decena = '0'; Char_t centena = '0';

	Char_t *file_pythia;
	file_pythia = (Char_t*) malloc(512*sizeof(Char_t));
	strcpy(file_pythia,local_path);
	strcat(file_pythia,head_folder);
	strcat(file_pythia,current_folder);
	strcat(file_pythia,"Events/run_01/output_pythia8.root");

	Char_t *file_delphes;
	file_delphes = (Char_t*) malloc(512*sizeof(Char_t));
	strcpy(file_delphes,local_path);
	strcat(file_delphes,head_folder);
	strcat(file_delphes,current_folder);
	strcat(file_delphes,"Events/run_01/output_delphes.root");

	cout << "\nReading the default files: \nPythia8: " << file_pythia << "\nDelphes: " << file_delphes << endl;

	chain_Pythia.Add(file_pythia);
	chain_Delphes.Add(file_delphes);

	// Objects of class ExRootTreeReader for reading the information
	ExRootTreeReader *treeReader_Pythia = new ExRootTreeReader(&chain_Pythia);
	ExRootTreeReader *treeReader_Delphes = new ExRootTreeReader(&chain_Delphes);

	Long64_t numberOfEntries = treeReader_Pythia->GetEntries();
	Long64_t numberOfEntries_Delphes = treeReader_Delphes->GetEntries();

	// Get pointers to branches used in this analysis
	TClonesArray *branchParticlePythia = treeReader_Pythia->UseBranch("GenParticle");
	TClonesArray *branchJet = treeReader_Delphes->UseBranch("Jet");
	TClonesArray *branchMissingET = treeReader_Delphes->UseBranch("MissingET");

	cout << endl;
	cout << " Number of Entries Pythia = " << numberOfEntries << endl;
	cout << " Number of Entries Delphes = " << numberOfEntries_Delphes << endl;
	cout << endl;

	// particles, jets and vectors
	TRootGenParticle *particle_pythia;
	TRootGenParticle *ISR_particle;
    MissingET *METpointer;
    TLorentzVector *vect_ISR_particle = new TLorentzVector;

	/*
	 * Main cycle of the program
	 */
	numberOfEntries = 10;
	for (Int_t entry = 0; entry < numberOfEntries; ++entry){
		// Load selected branches with data from specified event
		treeReader_Pythia->ReadEntry(entry);
		treeReader_Delphes->ReadEntry(entry);

		// MET
		METpointer = (MissingET*) branchMissingET->At(0);
		Double_t MET = METpointer->MET;
		cout << "Event: " << entry << " MET: "<< MET << endl;
	}

	return 0;
}
