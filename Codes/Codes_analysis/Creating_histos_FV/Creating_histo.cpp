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
./Creating_histo config_file.txt [N1 N2 N3]

where [config_file.txt] is a configuration file with
all parameters needed for the simulation.

N1 N2 and N3 are the index of the 3 variables.
If no parameter is passed as parameter, N1 N2 and N3
will be 0,1 and 2 by default.
*/


#include "ROOTFunctions.h"
#include "graphs_Funcs.h"
#include "functions.h"
#include "histoN.h"
#include "DelphesFunctions.h"

// Global Variables
const Double_t PI = TMath::Pi();

int main(int argc, char **argv){
	std::cout.precision(4);
	// Counting time
	Double_t initialTime = clock();

	// Folder variables
	string head_folder = "/home/af.garcia1214/PhenoMCsamples/Simulations/MG_pythia8_delphes_parallel/_Tops_Events_WI_Matching/";
	string current_folder = "_Tops_MG_1K_AG_WI_003/";

	string head_folder_binary = "/home/af.garcia1214/PhenoMCsamples/Results_Improved_Codes/matching_Results/_Tops_matchs_WI_Matching/";
	string matching_name = "ISR_jets_Tops_WI_003.bn";

	string head_folder_results = "/home/af.garcia1214/PhenoMCsamples/Results_Improved_Codes/histo_folder/_Tops_histos_WI_Matching/";

	// Checking input parameters
	string config_file_name = "Debug/config_file.txt";
	// Reading the file as first parameter
	if (argc>1){
		config_file_name = argv[1];
	}
	else{
		cout << "It is necessary to type a configuration file as parameter. Execute as ./Creating_histo config_file.txt [N1 N2 N3]" << endl;
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
			// Reading head folder binary
			else if (var_name.compare("head_folder_binary") == 0){
				head_folder_binary = var_value;
				cout << "\tVariable head folder binary set as: " << head_folder_binary << endl;
			}
			// Reading matching name
			else if (var_name.compare("matching_name") == 0){
				matching_name = var_value;
				cout << "\tVariable matching_name set as: " << matching_name << endl;
			}
			// Reading head folder results
			else if (var_name.compare("head_folder_results") == 0){
				head_folder_results = var_value;
				cout << "\tVariable head folder results set as: " << head_folder_results << endl;
			}

			number_line ++;
		}
	}
	else
	{
		cout << "ERROR: File " << config_file_name << " does not exist. Terminating program" << endl;
		return 0;
	}


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

	if (argc == 5){
		cout << "Filling histograms with the variables:" << endl;
		for (Int_t ind = 0; ind < 3; ind ++){
			var_index[ind] = atoi(argv[ind+2]);
		}
		cout << endl;
	}
	else if (argc == 1 || argc == 2) {
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

	// For loop over several simulations. iRun is the seed of the current simulation
	for(int iRun = 11; iRun < 261; iRun ++){
		// Create chains of root trees
		TChain chain_Delphes("Delphes");

		Char_t unidad = 0x30 + iRun%10;
        Char_t decena = 0x30 + int(iRun/10)%10;
        Char_t centena = 0x30 + int(iRun/100)%10;

		current_folder[current_folder.size()-4] = centena;
		current_folder[current_folder.size()-3] = decena;
		current_folder[current_folder.size()-2] = unidad;
		matching_name[matching_name.size()-6] = centena;
		matching_name[matching_name.size()-5] = decena;
		matching_name[matching_name.size()-4] = unidad;

		string file_delphes_str = head_folder + current_folder + "Events/run_01/output_delphes.root";

		Char_t *file_delphes = (Char_t *) file_delphes_str.c_str();

        cout << "\nWriting run: "<<centena<<decena<<unidad<<endl;
		cout << "\tReading the file: \n\tDelphes: " << file_delphes << endl;

		chain_Delphes.Add(file_delphes);
		// Objects of class ExRootTreeReader for reading the information
		ExRootTreeReader *treeReader_Delphes = new ExRootTreeReader(&chain_Delphes);

		Long64_t numberOfEntries = treeReader_Delphes->GetEntries();

		// Get pointers to branches used in this analysis
		TClonesArray *branchJet = treeReader_Delphes->UseBranch("Jet");
		TClonesArray *branchMissingET = treeReader_Delphes->UseBranch("MissingET");

		cout << "\tNumber of Entries Delphes = " << numberOfEntries << endl;
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

		string fileName_str = head_folder_binary + matching_name;

	    Char_t * fileName = (Char_t *) fileName_str.c_str();

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
				cout<<"\tprogress = "<<(entry*100/numberOfEntries)<<"%\t";
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

		cout<<"\tprogress = 100%\t";
		cout<<"Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;
		cout<<"\n\tNumber of Written Events: "<<numMatches<<endl;
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
	string combination = "______"; // Combination of variables
	for (Int_t ind = 0; ind < dims; ind ++){
		combination[(ind*2)+1] = (Char_t) (0x30 + var_index[ind]); // Int to char
	}

	string info_ISR_name_str = head_folder_results + "info_histo_ISR" + combination + ".txt";
	Char_t *info_ISR_name = (Char_t *) info_ISR_name_str.c_str();

	string array_ISR_name_str = head_folder_results + "array_histo_ISR" + combination + ".bn";
	Char_t *array_ISR_name = (Char_t *) array_ISR_name_str.c_str();

	string info_Non_ISR_name_str = head_folder_results + "info_histo_Non_ISR" + combination + ".txt";
	Char_t *info_Non_ISR_name = (Char_t *) info_Non_ISR_name_str.c_str();

	string array_Non_ISR_name_str = head_folder_results + "array_histo_Non_ISR" + combination + ".bn";
	Char_t *array_Non_ISR_name = (Char_t *) array_Non_ISR_name_str.c_str();

	cout << "Output files:\n\t" << info_ISR_name << "\n\t" << array_ISR_name << "\n\t" << info_Non_ISR_name << "\n\t" << array_Non_ISR_name << endl;

	histoISR->writeClass((Char_t*) info_ISR_name,(Char_t*) array_ISR_name);
	histoNonISR->writeClass((Char_t*) info_Non_ISR_name,(Char_t*)array_Non_ISR_name);
	cout<< "Time writing the file:"<< (clock()-partialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;

	cout<<"Fin :)"<<endl;

	return 0;
}
