/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

This algorithm looks for the ISR parton into the
pythia8 simulation file and then finds the
corresponding ISR jet

It also stores in a binary file the matching
results

To run, type

./ISR_matching_improved [config.txt] [000]

where [config.txt] is the configuration file and
[000] is the seed of the simulation under analysis
*/

#include <iostream>
#include "ROOTFunctions.h"
#include "graphs_Funcs.h"
#include "functions.h"
#include "DelphesFunctions.h"

using namespace std;
// Global Variables
const Double_t PI = TMath::Pi();

int main(int argc, char **argv){

	std::cout.precision(4);
	// Counting time
	Double_t initialTime = clock();

	// Folder variables
	string head_folder = "/home/af.garcia1214/PhenoMCsamples/Simulations/MG_pythia8_delphes_parallel/_Tops_Events_WI_Matching/";
	string current_folder = "_Tops_MG_1K_AG_WI_003/";

	string head_folder_results = "/home/af.garcia1214/PhenoMCsamples/Results_Improved_Codes/matching_Results/_Tops_matchs_WI_Matching/";
	string matching_name = "ISR_jets_Tops_WI_003.bn";

	// Checking input parameters
	string config_file_name = "Debug/config_file.txt";
	// Reading the file as first parameter
	if (argc>1){
		config_file_name = argv[1];
	}
	else{
		cout << "It is necessary to type a configuration file as parameter. Execute as ./ISR_matching config_file.txt [000]" << endl;
		return 1;
	}
	cout << "Reading input parameters" << endl;
	cout << "\tUsing as parameters' file: " << config_file_name << endl;

	ifstream config_file (config_file_name);
	if (config_file.is_open()){
		cout << "\tReading file" << endl;
		string line;
		int number_line = 1;
		while (getline(config_file,line)){
			// Skipping commented lines
			if (line[0] == '!')
				continue;

			// Finding the position of the equal sign
			int pos_equal = -1;
			pos_equal = line.find('=');

			if (pos_equal == -1){
				cout << "\tLine " << number_line << " is incorrect" << endl;
				continue;
			}

			// Splitting the line according to the position of equal sign
			string var_name = line.substr(0,pos_equal);
			string var_value = line.substr(pos_equal+1);

			// Reading head folder
			if(var_name.compare("head_folder") == 0){
				head_folder = var_value;
				cout << "\tVariable head folder set as: " << head_folder << endl;
			}
			// Reading current folder
			else if (var_name.compare("current_folder") == 0){
				current_folder = var_value;
				cout << "\tVariable current folder set as: " << current_folder <<endl;
			}
			// Reading head folder results
			else if (var_name.compare("head_folder_results") == 0){
				head_folder_results = var_value;
				cout << "\tVariable head folder results set as: " << head_folder_results << endl;
			}
			// Reading matching name
			else if (var_name.compare("matching_name") == 0){
				matching_name = var_value;
				cout << "\tVariable matching_name set as: " << matching_name << endl;
			}

			number_line ++;
		}
	}
	else
	{
		cout << "ERROR: File " << config_file_name << " does not exist. Terminating program" << endl;
		return 0;
	}

	// Reading the seed of the simulation. This parameter is optional and is the second of argv

	Char_t unidad = '3'; Char_t decena = '0'; Char_t centena = '0';
	if (argc > 2){
		cout << "\tRemember: The number of the simulation should consist of 3 digits" << endl;
		centena = argv[2][0];
		decena = argv[2][1];
		unidad = argv[2][2];
		current_folder[current_folder.size()-4] = centena;
		current_folder[current_folder.size()-3] = decena;
		current_folder[current_folder.size()-2] = unidad;
		matching_name[matching_name.size()-6] = centena;
		matching_name[matching_name.size()-5] = decena;
		matching_name[matching_name.size()-4] = unidad;
	}

	cout << "\tThe seed of the simulation is: " << centena << decena << unidad << endl;

	// Full path name of pythia and Delphes simulations
	string file_pythia_str = head_folder + current_folder + "Events/run_01/output_pythia8.root";
	Char_t *file_pythia = (Char_t *) file_pythia_str.c_str(); //Pass string to char_t *

	string file_delphes_str = head_folder + current_folder + "Events/run_01/output_delphes.root";
	Char_t *file_delphes = (Char_t *) file_delphes_str.c_str();

	if (argc > 2){
		cout << "\n\tReading the files: \n\tPythia8: " << file_pythia << "\n\tDelphes: " << file_delphes << endl;
	}
	else
		cout << "\n\tReading the default files: \n\tPythia8: " << file_pythia << "\n\tDelphes: " << file_delphes << endl;


	// Loading simulations of Pythia and Delphes
	cout << "\nLoading simulations of Pythia and Delphes" << endl;
	// Create chains of root trees
	TChain chain_Pythia("STDHEP");
	TChain chain_Delphes("Delphes");

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
	cout << "\tNumber of Entries Pythia = " << numberOfEntries << endl;
	cout << "\tNumber of Entries Delphes = " << numberOfEntries_Delphes << endl;
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
	 * Main cycle of the program
	 */
    cout << "Running the matching algorithm" << endl;
	numberOfEntries = 100000;
	for (Int_t entry = 0; entry < numberOfEntries; ++entry){
		// Progress
		if(numberOfEntries>10 && (entry%((int)numberOfEntries/10))==0.0){
			cout<<"\tprogress = "<<(entry*100/numberOfEntries)<<"%\t";
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
        	cout << "Error en el matching. Terminating program" << endl;
        	return 1;
        }
	}

    cout<<"\tprogress = 100%\t";
    cout<< "Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;

    /*
     * Writing results
     */
    cout << "\nWriting files" << endl;
    string fileName_str = head_folder_results + matching_name;

    Char_t * fileName = (Char_t *) fileName_str.c_str();

    if (argc > 2)
		cout << "\t Writing the binary file...:" << fileName << endl;
    else
        cout<<"\t Writing the default binary file...:" << fileName << endl;

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

    cout << "\nSome overal results: " << endl;
    cout << "\tNumber of events with a single ISR jet = " << NumEvents1ISRJet <<endl;
    cout << "\tNumber of matches = " << NumMatches << endl;
    cout << endl;

	return 0;
}
