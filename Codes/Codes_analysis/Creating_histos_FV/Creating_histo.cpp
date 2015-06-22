/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

This algorithm fills 2 N-dimensional histograms.
The histograms contain kinematic variables of ISR
jets and non ISR jets.

The user can choose 3 of 8 variables for filling
 the histograms:
1. PT
2. Abs(Eta) // Eta is a pair function
3. Delta Phi_MET
4. PT_ratio
5. Delta Eta_aver
6. Delta Phi_MET_others
7. Delta PT_others
8. Delta Eta_others

In order to choose them, the code should be run as
./Creating_histo N1 N2 N3, where N1 N2 and N3 are
the index of the 3 variables. If no parameter is
passed as parameter, N1 N2 and N3 will be 0,1 and 3
by default.
*/


#include "ROOTFunctions.h"
#include "graphs_Funcs.h"
#include "functions.h"
#include "histoN.h"
#include "DelphesFunctions.h"

// Global Variables
const Double_t PI = TMath::Pi();

// Other simulations parameters
const Char_t channel_histo = '_'; // 's' for sTops and '_' for Tops
const Char_t ISR_or_NOT_histo[] = "WI"; // "WI" with ISR, "SI" without (Here it does not make any sense), "bb" bjets production
const Bool_t atServer = true; // True if it is run at the server, false at the university's pc
const Bool_t Matching = true; // True if a matching has been done between MG and Pythia, false otherwise

int main(int argc, char **argv){
	std::cout.precision(4);
	// Counting time
	Double_t initialTime = clock();

	cout << "\n *** Creating histograms *** \n" << endl;

	// Variables for initializing histograms
    Int_t dims = 3;
    Double_t min_Values[3] = {0,-5.2,0};
    Double_t max_Values[3] = {800,5.2,PI};

	/*
	 * Read inputs and set variables for analysis
	 */
	Int_t var_index[3] = {0,1,2}; // Index of the 3 variables for analysis. By default 0, 1 and 2
	string variables[8] = {"PT","Abs(Eta)","Delta Phi_MET","PT_ratio","Delta Eta_aver","Delta Phi_MET_others","Delta PT_leading","Delta Eta_leading"};
	Double_t var_values[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Vector with the values of the 8 variables
	// Min and maximun values of the eight variables
	Double_t var_min_values[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	Double_t var_max_values[8] = {800,5.2,PI,8.0,7.0,PI,500,6.5};

	if (argc == 4){
		cout << "Filling histograms with the variables:" << endl;
		for (Int_t ind = 0; ind < 3; ind ++){
			var_index[ind] = atoi(argv[ind+1]);
		}
		cout << endl;
	}
	else if (argc == 1) {
		cout << "Filling histograms with the default variables:" << endl;
	}
	else {
		cout << "Error at calling this algorithm. Use as:" << endl;
		cout << "\t ./Creating_histo N1 N2 N3 or ./Creating_histo" << endl;
		cout << "Read the documentation at the beginning of the code for further information\n" << endl;
		return 1;
	}

	cout << "Var \t\t min_Value \t max_Value" << endl;
	for (Int_t ind = 0; ind < 3; ind ++){
		min_Values[ind] = var_min_values[var_index[ind]];
		max_Values[ind] = var_max_values[var_index[ind]];
		cout << var_index[ind] << ". " << variables[var_index[ind]] <<
				"\t" << min_Values[ind] << "\t" << max_Values[ind] << endl;
	}
	cout << endl;


	/*
	 * Initializing the 3-dimensional histogram
	 */
    Int_t bins[3] = {20,20,20};
    histoN* histoISR = new histoN(dims,min_Values,max_Values,bins);
    histoN* histoNonISR = new histoN(dims,min_Values,max_Values,bins);
    // Input variables of each histogram
	Double_t values[3] = {0.0,0.0,0.0};

	for(int iRun = 100; iRun < 261; iRun ++){
		// Create chains of root trees
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
		head_folder[0] = channel_histo;
		head_folder[13] = ISR_or_NOT_histo[0];
		head_folder[14] = ISR_or_NOT_histo[1];

		Char_t current_folder[] = "_Tops_MG_1K_AG_WI_003/";
		current_folder[0] = channel_histo;
		current_folder[15] = ISR_or_NOT_histo[0];
		current_folder[16] = ISR_or_NOT_histo[1];

        Char_t unidad = 0x30 + iRun%10;
        Char_t decena = 0x30 + int(iRun/10)%10;
        Char_t centena = 0x30 + int(iRun/100)%10;

		current_folder[18] = centena;
		current_folder[19] = decena;
		current_folder[20] = unidad;

		Char_t *file_delphes;
		file_delphes = (Char_t*) malloc(512*sizeof(Char_t));
		strcpy(file_delphes,local_path);
		strcat(file_delphes,head_folder);
		strcat(file_delphes,current_folder);
		strcat(file_delphes,"Events/run_01/output_delphes.root");

        cout << "Writing run: "<<centena<<decena<<unidad<<endl;
		cout << "\nReading the file: \nDelphes: " << file_delphes << endl;

		chain_Delphes.Add(file_delphes);
		// Objects of class ExRootTreeReader for reading the information
		ExRootTreeReader *treeReader_Delphes = new ExRootTreeReader(&chain_Delphes);

		Long64_t numberOfEntries = treeReader_Delphes->GetEntries();

		// Get pointers to branches used in this analysis
		TClonesArray *branchJet = treeReader_Delphes->UseBranch("Jet");
		TClonesArray *branchMissingET = treeReader_Delphes->UseBranch("MissingET");

		cout << endl;
		cout << " Number of Entries Delphes = " << numberOfEntries << endl;
		cout << endl;

		// particles, jets and vectors
		MissingET *METpointer;
		TLorentzVector *vect_currentJet = new TLorentzVector;
		TLorentzVector *vect_auxJet = new TLorentzVector;
		TLorentzVector *vect_leading = new TLorentzVector;
		Jet *currentJet = new Jet;
		Jet *auxJet = new Jet;

		// Temporary variables
		Double_t MET = 0.0; // Missing transverse energy
		Double_t delta_phi = 0.0; // difference between the phi angle of MET and the jet
		Double_t transverse_mass = 0.0; // Transverse mass
		Int_t numMatches = 0; // Number of matched jets
		Double_t delta_PT_jet = 0.0; // |PT-<PT>|
		Double_t PT_sum = 0.0; // sum(PT)
		Double_t PT_aver = 0.0; // <PT>
		Double_t Delta_eta_aver = 0.0; // sum_i|eta-eta_i|/(Nj-1)
		Double_t Delta_phi_sum = 0.0; // sum delta_phi
		Double_t Delta_phi_other_jets = 0.0; // Average of delta phi of other jets
		Double_t PT_ratio = 0.0; // PT/PT_others
		Double_t Delta_PT_leading = 0.0; // PT - PT_leading
		Double_t Delta_Eta_leading = 0.0; // |Eta - Eta_leading|

		// Jet with greatest PT
		Double_t PT_max = 0;
		Int_t posLeadingPT = -1;
		Int_t ISR_greatest_PT = 0;
		Double_t MT_leading_jet = 0.0; // Transverse mass

		/*
		 * Some variables used through the code
		 */
		Int_t ISR_jets[numberOfEntries];
		Int_t NumJets = 0;

	    Char_t *local_path_binary;
	    local_path_binary = (Char_t*) malloc(512*sizeof(Char_t));
	    if (atServer)
	    	strcpy(local_path_binary,"/home/af.garcia1214/PhenoMCsamples/Results/matching_Results/"); // At the server
	    else
	    	strcpy(local_path_binary,"/home/afgarcia1214/Documentos/Results_and_data/matching_Results/"); // At the University's pc

		Char_t *head_folder_binary;
		head_folder_binary = (Char_t*) malloc(512*sizeof(Char_t));
		if (Matching)
			strcpy(head_folder_binary,"_Tops_matchs_WI_Matching/");
		else
			strcpy(head_folder_binary,"_Tops_matchs_WI/");
	    head_folder_binary[0] = channel_histo;
	    head_folder_binary[13] = ISR_or_NOT_histo[0];
	    head_folder_binary[14] = ISR_or_NOT_histo[1];

	    Char_t matching_name[] = "ISR_jets_Tops_WI_003.bn";
	    matching_name[8] = channel_histo;
	    matching_name[14] = ISR_or_NOT_histo[0];
	    matching_name[15] = ISR_or_NOT_histo[1];

		matching_name[17] = centena;
		matching_name[18] = decena;
		matching_name[19] = unidad;

	    Char_t * fileName;
	    fileName = (Char_t*) malloc(512*sizeof(Char_t));
	    strcpy(fileName,local_path_binary);
	    strcat(fileName,head_folder_binary);
	    strcat(fileName,matching_name);

		ifstream ifs(fileName,ios::in | ios::binary);

		for (Int_t j = 0; j<numberOfEntries; j++){
				ifs.read((Char_t *) (ISR_jets+j),sizeof(Int_t));
		}
		ifs.close();

		/*
		 * Main cycle of the program
		 */
		numberOfEntries = 100000;
		for (Int_t entry = 0; entry < numberOfEntries; ++entry){
			// Progress
			if(numberOfEntries>10 && (entry%((int)numberOfEntries/10))==0.0){
				cout<<"progress = "<<(entry*100/numberOfEntries)<<"%\t";
				cout<< "Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;
			}

			// Load selected branches with data from specified event
			treeReader_Delphes->ReadEntry(entry);

			// MET
			METpointer = (MissingET*) branchMissingET->At(0);
			MET = METpointer->MET;

			NumJets=branchJet->GetEntries();

			// checking the ISR
			if (ISR_jets[entry] == -1 || NumJets < 3)
				continue;

			if (ISR_jets[entry] >= NumJets){
				cout << "Error en el matching" << endl;
				return 1;
			}

			// 3 PT ratio
			PT_aver = 0.0;
			PT_sum = 0.0;
			PT_ratio = 0.0;

			// 4 Delta Eta aver
			Delta_eta_aver = 0.0;

			// 5 Delta Phi others
			Delta_phi_sum = 0.0;
			Delta_phi_other_jets = 0.0;

			// 6 Delta PT leading
			PT_max = 0.0;
			Delta_PT_leading = 0.0;
			delta_PT_jet = 0.0; // If needed

			// 7 Delta Eta leading
			Delta_Eta_leading = 0.0;

			// Reset Var_values (Not necessary)
			for(Int_t ind = 0; ind < 8; ind++){
				var_values[ind] = 0.0;
				if (ind < dims) values[ind] = 0.0;
			}

			// Preliminary for. It is used to calculate PT_aver and Delta_phi_sum
			for (Int_t iJet = 0; iJet<NumJets; iJet++){
				currentJet = (Jet*) branchJet->At(iJet);
				vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
				PT_sum += vect_currentJet->Pt();
				delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);
				Delta_phi_sum += delta_phi;
				// PT Leading jet
				if(PT_max < vect_currentJet->Pt()){
					PT_max = vect_currentJet->Pt();
					posLeadingPT = iJet;
				}
			}

			numMatches++;

			//PT_aver
			PT_aver = PT_sum/NumJets;

			// Leading PT
			currentJet = (Jet*) branchJet->At(posLeadingPT);
			vect_leading->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);

			for (Int_t iJet = 0; iJet<NumJets; iJet++){
				currentJet = (Jet*) branchJet->At(iJet);
				vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);

				// 2 Delta Phi MET
				delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);

				// PT ratio
				PT_ratio = vect_currentJet->Pt()*(NumJets-1)/(PT_sum-vect_currentJet->Pt());

				// 4 Delta Eta Aver
				Delta_eta_aver = 0.0;
				// For cycle used to calculate Delta_eta_aver
				for(Int_t iJet2 = 0; iJet2<NumJets; iJet2++){
					auxJet = (Jet*) branchJet->At(iJet2);
					vect_auxJet->SetPtEtaPhiM(auxJet->PT,auxJet->Eta,auxJet->Phi,auxJet->Mass);
					if (iJet2 != iJet) Delta_eta_aver += TMath::Abs(vect_auxJet->Eta()-vect_currentJet->Eta());
				}
				Delta_eta_aver = Delta_eta_aver/(NumJets-1);

				// 5 Delta Phi MET Others
				Delta_phi_other_jets = (Delta_phi_sum-delta_phi)/(NumJets-1);

				// 6 Delta PT leading
				Delta_PT_leading = vect_leading->Pt()-vect_currentJet->Pt();

				// 7 Delta Eta leading
				Delta_Eta_leading = TMath::Abs(vect_currentJet->Eta()-vect_leading->Eta());

				// Other variables
				delta_PT_jet = TMath::Abs(vect_currentJet->Pt()-PT_aver);
				transverse_mass = sqrt(2*vect_currentJet->Pt()*MET*(1-cos(delta_phi)));

				// Filling the array with the variables' values
				var_values[0] = vect_currentJet->Pt();
				var_values[1] = TMath::Abs(vect_currentJet->Eta());
				var_values[2] = delta_phi;
				var_values[3] = PT_ratio;
				var_values[4] = Delta_eta_aver;
				var_values[5] = Delta_phi_other_jets;
				var_values[6] = Delta_PT_leading;
				var_values[7] = Delta_Eta_leading;

				for (Int_t ind = 0; ind < dims; ind++){
					int pos = *(var_index+ind);
					values[ind] = *(var_values+pos);
				}

				if (iJet != ISR_jets[entry]){
					// Non ISR jet
					histoNonISR->fill(values);
				}
				else{
					// ISR jet
					histoISR->fill(values);
				}

			}
		}

		cout<<"progress = 100%\t";
		cout<<"Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;
		cout<<"\nNumber of Written Events: "<<numMatches<<endl;
	} // End run's for cicle

	/*
	 * Writing the histogram
	 */
    // Counting time
	Double_t partialTime = clock();
	cout<< "Time building the histogram:"<< (partialTime-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;

	// Writing the histogram
	cout<<"Min value 1: "<<min_Values[0]<<endl;
	Int_t* freq;
	for(Int_t j = 0; j<dims; j++){
			cout<<"ISR Jets - Events of the dimension:\t"<<j<<endl;
			freq = histoISR->getHistDim(j);
			for(Int_t i = 0; i<bins[j];i++){
					cout<<"\t"<<freq[i];
					if(i>0 && ((i+1)%10 == 0.0))
							cout<<endl;
			}
			cout<<endl;
	}

	cout<<endl<<"\t\t ***"<<endl<<endl;
	for(Int_t j = 0; j<dims; j++){
			cout<<"Non ISR Jets - Events of the dimension:\t"<<j<<endl;
			freq = histoNonISR->getHistDim(j);
			for(Int_t i = 0; i<bins[j];i++){
					cout<<"\t"<<freq[i];
					if(i>0 && ((i+1)%10 == 0.0))
							cout<<endl;
			}
			cout<<endl;
	}

	cout<<"Entries: "<<histoISR->getEntries()<<endl;

	/*
	 * Creating histograms
	 */
	cout<<"\nWriting..."<<endl;

	// Defining the names of the files
	Char_t combination[] = "______"; // Combination of variables
	for (Int_t ind = 0; ind < dims; ind ++){
		*(combination+(ind*2)+1) = (Char_t) (0x30 + var_index[ind]); // Int to char
	}

    Char_t *local_path_results;
    local_path_results = (Char_t*) malloc(512*sizeof(Char_t));
    if (atServer)
    	strcpy(local_path_results,"/home/af.garcia1214/PhenoMCsamples/Results/histo_folder/"); // At the server
    else
    	strcpy(local_path_results,"/home/afgarcia1214/Documentos/Results_and_data/histo_folder/"); // At the University's pc

    Char_t *head_folder_results;
	head_folder_results = (Char_t*) malloc(512*sizeof(Char_t));
	if (Matching)
		strcpy(head_folder_results,"_Tops_histos_WI_Matching/");
	else
		strcpy(head_folder_results,"_Tops_histos_WI/");
    head_folder_results[0] = channel_histo;
    head_folder_results[13] = ISR_or_NOT_histo[0];
    head_folder_results[14] = ISR_or_NOT_histo[1];

	Char_t *info_ISR_name;
	info_ISR_name = (Char_t*) malloc(sizeof(char)*512);
	strcpy(info_ISR_name,local_path_results);
	strcat(info_ISR_name,head_folder_results);
	strcat(info_ISR_name,"info_histo_ISR");
	strcat(info_ISR_name,combination);
	strcat(info_ISR_name,".txt");

	Char_t *array_ISR_name;
	array_ISR_name = (Char_t*) malloc(sizeof(char)*512);
	strcpy(array_ISR_name,local_path_results);
	strcat(array_ISR_name,head_folder_results);
	strcat(array_ISR_name,"array_histo_ISR");
	strcat(array_ISR_name,combination);
	strcat(array_ISR_name,".bn");

	Char_t *info_Non_ISR_name;
	info_Non_ISR_name = (Char_t*) malloc(sizeof(char)*512);
	strcpy(info_Non_ISR_name,local_path_results);
	strcat(info_Non_ISR_name,head_folder_results);
	strcat(info_Non_ISR_name,"info_histo_Non_ISR");
	strcat(info_Non_ISR_name,combination);
	strcat(info_Non_ISR_name,".txt");

	Char_t *array_Non_ISR_name;
	array_Non_ISR_name = (Char_t*) malloc(sizeof(char)*512);
	strcpy(array_Non_ISR_name,local_path_results);
	strcat(array_Non_ISR_name,head_folder_results);
	strcat(array_Non_ISR_name,"array_histo_Non_ISR");
	strcat(array_Non_ISR_name,combination);
	strcat(array_Non_ISR_name,".bn");

	cout << "Output files:\n\t" << info_ISR_name << "\n\t" << array_ISR_name << "\n\t" << info_Non_ISR_name << "\n\t" << array_Non_ISR_name << endl;

	histoISR->writeClass((Char_t*) info_ISR_name,(Char_t*) array_ISR_name);
	histoNonISR->writeClass((Char_t*) info_Non_ISR_name,(Char_t*)array_Non_ISR_name);
	cout<< "Time writing the file:"<< (clock()-partialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;

	cout<<"Fin :)"<<endl;

	return 0;
}
