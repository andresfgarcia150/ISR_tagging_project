/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Fisica       -------
-------        Joven Investigador         -------
-------  Andres Felipe Garcia Albarracin  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

This algorithm looks for the ISR parton into the
pythia8 simulation file and then finds the
corresponding ISR jet

It also stores in a binary file the matching
results
*/

#include <iostream>
#include "ROOTFunctions.h"
#include "graphs_Funcs.h"
#include "functions.h"
#include "DelphesFunctions.h"

using namespace std;
// Global Variables
const Double_t PI = TMath::Pi();

// Other simulations parameters
const Char_t channel = 's'; // 's' for sTops and '_' for Tops
const Char_t ISR_or_NOT[] = "WI"; // "WI" with ISR, "SI" without (Here it does not make any sense), "bb" bjets production
const Bool_t atServer = true; // True if it is run at the server, false at the university's pc
const Bool_t Matching = true; // True if a matching has been done between MG and Pythia, false otherwise

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
	if (atServer)
		strcpy(local_path,"/home/af.garcia1214/PhenoMCsamples/Simulations/MG_pythia8_delphes_parallel/"); // At the server
	else
		strcpy(local_path,"/home/afgarcia1214/Documentos/Simulations/"); // At the University's pc

	Char_t *head_folder;
	head_folder = (Char_t*) malloc(512*sizeof(Char_t));
	if (Matching)
		strcpy(head_folder,"_Tops_Events_WI_Matching/");
	else
		strcpy(head_folder,"_Tops_Events_WI/");
	head_folder[0] = channel;
	head_folder[13] = ISR_or_NOT[0];
	head_folder[14] = ISR_or_NOT[1];

	Char_t current_folder[] = "_Tops_MG_1K_AG_WI_003/";
	current_folder[0] = channel;
	current_folder[15] = ISR_or_NOT[0];
	current_folder[16] = ISR_or_NOT[1];

	Char_t unidad = '3'; Char_t decena = '0'; Char_t centena = '0';

	if (argc > 1){
		cout << "The number of the simulation should consist of 3 digits" << endl;
		centena = argv[1][0];
		decena = argv[1][1];
		unidad = argv[1][2];
		current_folder[18] = centena;
		current_folder[19] = decena;
		current_folder[20] = unidad;
	}

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

	if (argc > 1){
		cout << "\nReading the files: \nPythia8: " << file_pythia << "\nDelphes: " << file_delphes << endl;
	}
	else
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

	// Temporary variables
	Bool_t ISR_parton_found = false; // true if the initial ISR_parton (with status 43) was found
	Int_t pos_ISR = -1; // position of the ISR_parton into the branchParticlePythia array
	Double_t MET = 0.0; // Missing transverse energy

	/*
	 * Some variables used through the code
	 */
	Int_t NumEvents1ISRJet = 0;     // Number of events where the number of ISR jets is 1
	Int_t NumMatches = 0;           // Number of matches
	Int_t NumJets = 0;
	Int_t ISR_match_index = -1;
	Double_t Cut_matching_DPT = 50.0;
	Double_t Cut_matching_DEta = 0.4;
	Double_t Cut_matching_DPhi = 0.4;
	Double_t Cut_matching_Dy = 0.4;
	Int_t ISR_jets[numberOfEntries];

	/*
	 * Main cycle of the program. Cycle over the events
	 */
	numberOfEntries = 100000;
	for (Int_t entry = 0; entry < numberOfEntries; ++entry){
		// Progress
		if(numberOfEntries>10 && (entry%((int)numberOfEntries/10))==0.0){
			cout<<"progress = "<<(entry*100/numberOfEntries)<<"%\t";
			cout<< "Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;
		}

		// Load selected branches with data from specified event
		treeReader_Pythia->ReadEntry(entry);
		treeReader_Delphes->ReadEntry(entry);

		// By default, the ISR jet was not matched
		ISR_jets[entry] = -1;

		// MET
		METpointer = (MissingET*) branchMissingET->At(0);
		MET = METpointer->MET;

		// Finding the ISR parton
		ISR_parton_found = false;
		pos_ISR = -1;
		for(Int_t iPart = 0; iPart < branchParticlePythia->GetEntries(); iPart++){
			particle_pythia = (TRootGenParticle*) branchParticlePythia->At(iPart);
			if( abs(particle_pythia->Status) == 43){
				pos_ISR = iPart;
				ISR_particle = (TRootGenParticle*) branchParticlePythia->At(pos_ISR);
				ISR_parton_found = true;
//				The following lines were used to check that everything was going well
//				cout << pos_ISR << "\t\t" << ISR_particle->Status << "\t\t" << ISR_particle->PID
//					<< "\t\t" << ISR_particle->M1 << "\t\t" << ISR_particle->M2
//					<< "\t\t" << ISR_particle->D1 << "\t\t" << ISR_particle->D2 << endl;
			}
		}

		// If there is not ISR parton, pass to the next event
		if (ISR_parton_found == false){
			continue;
		}

		// Finding the last copy of the ISR_parton
		ISR_parton_found = false;
		while (!ISR_parton_found){
			if (ISR_particle->D1 != ISR_particle->D2)
				ISR_parton_found = true;
			else{
				pos_ISR = ISR_particle->D1;
				if(pos_ISR != -1) // To avoid an incoherent event
					ISR_particle = (TRootGenParticle*) branchParticlePythia->At(pos_ISR);
				else
					ISR_parton_found = true; // To end up the while loop
			}
		}

		if (pos_ISR == -1) // End the incoherent events
			continue;

		// Matching algorithm
		// Matching between the ISR parton and a jet
		// Auxiliary variables
		Double_t R_min = 2.0;
		Double_t r; // Current deltaR
		ISR_match_index = -1;
		Int_t mixJets = 0;
		TLorentzVector *vect_Jet1 = new TLorentzVector();       // Four-momentum of the jet of the 1st for
		TLorentzVector *vect_Jetc = new TLorentzVector();       // Four-momentum of the jet of the 2nd, 3rd ... for
		TLorentzVector *vect_Jets = new TLorentzVector();       // Four-momentum of the sum of jets
		TLorentzVector *vect_Jeto = new TLorentzVector();       // Four-momentum of the optimal combination
		Jet *jet = new Jet();
		Jet *jet2 = new Jet();

		NumJets = branchJet->GetEntries();
		vect_ISR_particle->SetPtEtaPhiE(ISR_particle->PT,ISR_particle->Eta,ISR_particle->Phi,ISR_particle->E);

		if (NumJets < 3) // Minimun 3 jets per event
			continue;

		// Finding the jet with the minimum R to the ISR parton
		for ( Int_t j = 0; j < NumJets; j++ ) {     // Loop over jets finding the one with the minimum R
				jet = (Jet*) branchJet->At(j);
				vect_Jet1->SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
				r = vect_ISR_particle->DeltaR(*vect_Jet1);
				if ( r < R_min ) {
						R_min = r;
						ISR_match_index = j;
						mixJets = 1;
						*vect_Jeto = *vect_Jet1;
				}
				// Checking if there are two jets mixed
				for ( Int_t k = j+1; k<NumJets; k++){
					jet2 = (Jet*) branchJet->At(k);
					vect_Jetc->SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
					*vect_Jets = *vect_Jet1 + *vect_Jetc;
					r = vect_ISR_particle->DeltaR(*vect_Jets);
					if ( r < R_min ) {
						R_min = r;
						ISR_match_index = j;
						mixJets = 2;
						*vect_Jeto = *vect_Jets;
					}
				// Checking if there are three jets mixed
				for (Int_t m = k+1; m<NumJets; m++){
					jet2 = (Jet*) branchJet->At(m);
					vect_Jetc->SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
					*vect_Jets = *vect_Jets + *vect_Jetc;
					r = vect_ISR_particle->DeltaR(*vect_Jets);
					if ( r < R_min ) {
						R_min = r;
						ISR_match_index = j;
						mixJets = 3;
						*vect_Jeto = *vect_Jets;
        	                }
				// Checking if there are four jets mixed
				for (Int_t n = m+1; n<NumJets; n++){
					jet2 = (Jet*) branchJet->At(n);
					vect_Jetc->SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
					*vect_Jets = *vect_Jets + *vect_Jetc;
					r = vect_ISR_particle->DeltaR(*vect_Jets);
					if ( r < R_min ) {
						R_min = r;
						ISR_match_index = j;
						mixJets = 4;
						*vect_Jeto = *vect_Jets;
					}
				}
                        	}
                		}
		}     // Loop over jets finding the one with the minimum R

		if( (mixJets == 1) && (ISR_match_index >= 0) && (ISR_match_index < NumJets) ) {
        	        NumEvents1ISRJet++;
        	        Double_t Delta_PT = TMath::Abs(vect_Jeto->Pt() - vect_ISR_particle->Pt());
        	        Double_t Delta_Eta = TMath::Abs(vect_Jeto->Eta() - vect_ISR_particle->Eta());
        	        Double_t Delta_Phi = vect_Jeto->DeltaPhi(*vect_ISR_particle);
        	        Double_t Delta_y = TMath::Abs(vect_Jeto->Rapidity() - vect_ISR_particle->Rapidity());

        	        if ( (Delta_PT > Cut_matching_DPT) || (Delta_Eta > Cut_matching_DEta) || (Delta_Phi > Cut_matching_DPhi ) || (Delta_y > Cut_matching_Dy) ) {
        	                ISR_jets[entry] = -1;
        	        }
        	        else {
        	                NumMatches++;
        	                ISR_jets[entry] = ISR_match_index;
        	        }
        	}

        	if (ISR_jets[entry] >= NumJets){
        		cout << "Error en el matching" << endl;
        		return 1;
        	}
	} // End of the cycle over the events

	cout<<"progress = 100%\t";
	cout<< "Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;

	/*
	 * Writing results
	 */
	Char_t *local_path_results;
	local_path_results = (Char_t*) malloc(512*sizeof(Char_t));
	if (atServer)
		strcpy(local_path_results,"/home/af.garcia1214/PhenoMCsamples/Results/matching_Results/"); // At the server
	else
	strcpy(local_path_results,"/home/afgarcia1214/Documentos/Results_and_data/matching_Results/"); // At the University's pc

	Char_t *head_folder_results;
	head_folder_results = (Char_t*) malloc(512*sizeof(Char_t));
	if (Matching)
		strcpy(head_folder_results,"_Tops_matchs_WI_Matching/");
	else
		strcpy(head_folder_results,"_Tops_matchs_WI/");
	head_folder_results[0] = channel;
	head_folder_results[13] = ISR_or_NOT[0];
	head_folder_results[14] = ISR_or_NOT[1];

	Char_t matching_name[] = "ISR_jets_Tops_WI_003.bn";
	matching_name[8] = channel;
	matching_name[14] = ISR_or_NOT[0];
	matching_name[15] = ISR_or_NOT[1];

	if (argc > 1){
		matching_name[17] = centena;
    		matching_name[18] = decena;
    		matching_name[19] = unidad;
        }

	Char_t * fileName;
	fileName = (Char_t*) malloc(512*sizeof(Char_t));
	strcpy(fileName,local_path_results);
	strcat(fileName,head_folder_results);
	strcat(fileName,matching_name);

	if (argc > 1)
		cout << "*** Writing the binary file...:" << fileName << endl;
    	else
        	cout<<"*** Writing the default binary file...:" << fileName << endl;

	ofstream ofs(fileName,ios::out|ios::binary);
	if (!ofs){
		cout << "Problemas al escribir el archivo" << endl;
	}
	else{
		for(Int_t j = 0; j<numberOfEntries; j++){
			ofs.write((Char_t *) (ISR_jets+j),sizeof(Int_t));
		}
	}
	ofs.close();

	cout << endl;
	cout << "Number of events with a single ISR jet = " << NumEvents1ISRJet <<endl;
	cout << "Number of matches = " << NumMatches << endl;
	cout << endl;

	return 0;
}
