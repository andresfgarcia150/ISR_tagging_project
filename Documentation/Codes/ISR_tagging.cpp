/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Fisica       -------
-------        Joven Investigador         -------
-------  Andres Felipe Garcia Albarracin  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

This algorithm tags ISR jet in a certain sample.
It takes 2 N-dimensional histograms which contain
information about ISR and Non ISR Jets as input
and developes the ISR tagging in another sample.

The user can choose 3 of 8 variables for
developing the algorithm
1. PT
2. Abs(Eta) // Eta is a pair function
3. Delta Phi_MET
4. PT_ratio
5. Delta Eta_aver
6. Delta Phi_MET_others
7. Delta PT_others
8. Delta Eta_others

In order to choose them, the code should be run as:

./ISR tagging config_file.txt [N1 N2 N3]

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

	string head_folder_histos = "/home/af.garcia1214/PhenoMCsamples/Results_Improved_Codes/histo_folder/_Tops_histos_WI_Matching/";
	string head_folder_results = "/home/af.garcia1214/PhenoMCsamples/Results_Improved_Codes/resultsTagging/_Tops_histos_WI_Matching/";

	Bool_t ISR_OR_NOT = true;

	// Variables for analysis
	Double_t pt_cut = 0.0; // ISR jet pt cut
	Double_t Jet_cut = 2;  // Ptobability cut 'K'
	Bool_t do_pt_cut = false;
	Bool_t do_jet_cut = false;

	// Checking input parameters
	string config_file_name = "Debug/config_file.txt";
	// Reading the file as first parameter
	if (argc>1){
		config_file_name = argv[1];
	}
	else{
		cout << "It is necessary to type a configuration file as parameter. Execute as ./ISR tagging config_file.txt [N1 N2 N3]" << endl;
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
			// Reading head folder histos
			else if (var_name.compare("head_folder_histos") == 0){
				head_folder_histos = var_value;
				cout << "\tVariable head folder histos set as: " << head_folder_histos << endl;
			}
			// Reading head folder results
			else if (var_name.compare("head_folder_results") == 0){
				head_folder_results = var_value;
				cout << "\tVariable head folder results set as: " << head_folder_results << endl;
			}
			// Reading pt_cut
			else if (var_name.compare("pt_cut") == 0){
				pt_cut = atof((Char_t *) var_value.c_str());
				do_pt_cut = true;
			}
			// Reading jet_cut
			else if (var_name.compare("Jet_cut") == 0){
				Jet_cut = atof((Char_t *) var_value.c_str());
				do_jet_cut = true;
			}
			// Reading ISR_OR_NOT
			else if (var_name.compare("ISR_OR_NOT") == 0){
				if (var_value.compare("1") == 0)
					ISR_OR_NOT = true;
				else
					ISR_OR_NOT = false;
			}

			number_line ++;
		}
	}
	else
	{
		cout << "ERROR: File " << config_file_name << " does not exist. Terminating program" << endl;
		return 0;
	}

	cout << "\n *** Running the tagging Algorithm *** \n" << endl;

	// Variables for initializing histograms
	Int_t dims = 3;

	/*
	 * Read inputs and set variables for analysis
	 */
	Int_t var_index[3] = {0,1,2}; // Index of the 3 variables for analysis. By default 0, 1 and 2
	string variables[8] = {"PT","Abs(Eta)","Delta Phi_MET","PT_ratio","Delta Eta_aver","Delta Phi_MET_others","Delta PT_leading","Delta Eta_leading"};
	Double_t var_values[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Vector with the values of the 8 variables


	if (argc == 2) {
		cout << "\tRunning the algorithm with the default variables:" << endl;
	}
	else if (argc >= 5){
		cout << "\tRunning the algorithm with the variables:" << endl;
		for (Int_t ind = 0; ind < 3; ind ++){
			var_index[ind] = atoi(argv[ind+2]);
		}
	}
	else {
		cout << "\tError at calling this algorithm. Use as:" << endl;
		cout << "\t ./ISR_tagging config_file.txt,  ./ISR_tagging config_file.txt N1 N2 N3 or just ./ISR_tagging" << endl;
		cout << "\tRead the documentation at the beginning of the code for further information\n" << endl;
		return 1;
	}

	cout << "\tVar \t\t min_Value \t max_Value" << endl;
	for (Int_t ind = 0; ind < 3; ind ++){

		cout << "\t" << var_index[ind] << ". " << variables[var_index[ind]] << endl;
	}
	cout << endl;

	cout << "\tTransverse momentum of the ISR: " << pt_cut << endl;

	/*
	 * Initializing the 3-dimensional histogram
	 */
	// Defining the names of the files
	string combination = "______"; // Combination of variables
	for (Int_t ind = 0; ind < dims; ind ++){
		combination[(ind*2)+1] = (Char_t) (0x30 + var_index[ind]); // Int to char
	}

	string info_ISR_name_str = head_folder_histos + "info_histo_ISR" + combination + ".txt";
	Char_t *info_ISR_name = (Char_t *) info_ISR_name_str.c_str();

	string array_ISR_name_str = head_folder_histos + "array_histo_ISR" + combination + ".bn";
	Char_t *array_ISR_name = (Char_t *) array_ISR_name_str.c_str();

	string info_Non_ISR_name_str = head_folder_histos + "info_histo_Non_ISR" + combination + ".txt";
	Char_t *info_Non_ISR_name = (Char_t *) info_Non_ISR_name_str.c_str();

	string array_Non_ISR_name_str = head_folder_histos + "array_histo_Non_ISR" + combination + ".bn";
	Char_t *array_Non_ISR_name = (Char_t *) array_Non_ISR_name_str.c_str();

	histoN* histoISR = new histoN(info_ISR_name,array_ISR_name);
	histoN* histoNonISR = new histoN(info_Non_ISR_name,array_Non_ISR_name);

	cout << "\tEntradas ISR: " << histoISR->getEntries() << endl;
	cout << "\tEntradas FSR: " << histoNonISR->getEntries() << endl;

	// Input variables of each histogram
	Double_t values[3] = {0.0,0.0,0.0};

	/*
	 * MET histograms
	 */
	TH1 *h_MET = new TH1F("Missing ET","All events",300,0,1000);
	Char_t *name_histo_MET;
	name_histo_MET = (Char_t*) malloc(sizeof(char)*512);
	strcpy(name_histo_MET,"ISR jet PT > ");
	Char_t pt_str[] = "   ";
	pt_str[0] = 0x30 + int(pt_cut/100)%10;
	pt_str[1] = 0x30 + int(pt_cut/10)%10;
	pt_str[2] = 0x30 + int(pt_cut)%10;
	strcat(name_histo_MET,pt_str);
	strcat(name_histo_MET,"-k = ");
	Char_t k_str[] = "   ";
	k_str[0] = 0x30 + int(Jet_cut)%10;
	k_str[1] = '.';
	k_str[2] = 0x30 + int(Jet_cut*10)%10;
	strcat(name_histo_MET,k_str);
	TH1 *h_MET_hpt1 = new TH1F(name_histo_MET,"Missing ET high_ISR_pt-1",300,0.0,1000);

	if (do_jet_cut && do_pt_cut)
		cout << "\tThe algorithm will evaluate the MET for a sample with PT > " << pt_str << " at k = " << k_str << endl;
	/*
	 * Tagging variables
	 */

	cout << "\tJet cut, k = " << Jet_cut << endl;

	// Arrays with the number of tags, Misstags and events rejected
	// Probability cut
	Double_t Prob_cut = 0;
	Double_t k_min = 1.2; // Minimum probability cut = k_min/num_jets
	Double_t k_max = 3.0; // Maximum probability cut = k_max/num_jets
	Int_t k_bins = 100; // Number of values of k between k_min and k_max
	Double_t k_step = (Double_t) (k_max-k_min)/k_bins;
	Double_t k_values[k_bins];
	for(Int_t ind = 0; ind < k_bins; ind ++){
		k_values[ind] = k_min + k_step*ind;
	}

	// Tagging results
	Int_t Num_Tags = 0;
	Int_t Num_MissTags = 0;
	Int_t Num_Rejected = 0;

	Double_t Num_Tags_array[k_bins];
	Double_t Num_MissTags_array[k_bins];
	Double_t Num_Rejected_array[k_bins];
	Double_t Num_Total_Jets[k_bins];

	Double_t Num_Tags_array_hpt[k_bins];
	Double_t Num_MissTags_array_hpt[k_bins];
//	Double_t Num_Rejected_array_hpt[k_bins];
	Double_t Num_Total_Jets_hpt[k_bins];


	for (Int_t ind = 0; ind < k_bins; ind ++){
		Num_Tags_array[ind] = 0;
		Num_MissTags_array[ind] = 0;
		Num_Rejected_array[ind] = 0;
		Num_Total_Jets[ind] = 0;
		Num_Tags_array_hpt[ind] = 0;
		Num_MissTags_array_hpt[ind] = 0;
//		Num_Rejected_array_hpt[ind] = 0;
		Num_Total_Jets_hpt[ind] = 0;
	}

	// Variables of the ISR tagging algorithm
	Double_t H_ISR, H_Non_ISR, alpha;
	Double_t prob_max = 0;
	Double_t probISR = 0;
	Double_t k_ISR = 0;
	Double_t k_ISR_pos = 0; // Position of the ISR in the vector
	Int_t ISR_tag_index = -1;

	// Cycle over several runs . iRun corresponds to the seed of the current run
	for(int iRun = 261; iRun < 270; iRun ++){
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

		cout << "\n\tStudying run: "<<centena<<decena<<unidad<<endl;
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

		if (ISR_OR_NOT == true){
			ifstream ifs(fileName,ios::in | ios::binary);

			for (Int_t j = 0; j<numberOfEntries; j++){
				ifs.read((Char_t *) (ISR_jets+j),sizeof(Int_t));
			}
			ifs.close();
		}
		else{
			for (Int_t j = 0; j<numberOfEntries; j++){
				ISR_jets[j] = -2; // There is not ISR jet but also there is not matching
			}
		}

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
			if (NumJets < 3 || ISR_jets[entry] == -1)
				continue;

			h_MET->Fill(MET);

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

			//PT_aver
			PT_aver = PT_sum/NumJets;

			// Leading PT
			currentJet = (Jet*) branchJet->At(posLeadingPT);
			vect_leading->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);

			// The best ISR candidate
			TLorentzVector *vect_optimum = new TLorentzVector;

			// Reset variables
			probISR = 0.0;
			k_ISR = 0.0;
			prob_max = 0;
			ISR_tag_index = -1;

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

				// Comparing with histos
				H_ISR = histoISR->getProbVal(values);
				H_Non_ISR = histoNonISR->getProbVal(values);

				if (H_ISR >3e-7 || H_Non_ISR>3e-7){
					alpha = NumJets/(H_Non_ISR*(NumJets-1)+H_ISR);
					probISR = alpha*H_ISR/NumJets;

					if(probISR > (1.0 + 1.0e-10)){
						cout << setprecision(20) << "\n\t *** ERROR: La probabilidad no puede ser mayor a 1 ***" << endl;
						return 1;
					}

					if (probISR >= prob_max){
						prob_max = probISR;
						vect_optimum->SetPtEtaPhiM(vect_currentJet->Pt(),vect_currentJet->Eta(),vect_currentJet->Phi(),vect_currentJet->M());
						ISR_tag_index = iJet;
					}
				}
			}

			k_ISR = prob_max*NumJets;

			// Check the tagging results
			k_ISR_pos = findPosition(k_min,k_max,k_bins,k_ISR);

			if(k_ISR == 0.0) k_ISR_pos = -1;

			if (ISR_jets[entry] != -1 && ISR_OR_NOT == true){
				// A comparison can be handled
				for (Int_t ind = 0; ind < k_ISR_pos + 1; ind++){
					if (ISR_tag_index == ISR_jets[entry])
						Num_Tags_array[ind]++;
					else
						Num_MissTags_array[ind]++;
				}
				for (Int_t ind = k_ISR_pos+1; ind < k_bins; ind++){
					Num_Rejected_array[ind]++;
				}
			}
			else if (ISR_jets[entry] == -2 && ISR_OR_NOT == false){
				for (Int_t ind = 0; ind < k_ISR_pos + 1; ind++){
					Num_MissTags_array[ind]++;
				}
				for (Int_t ind = k_ISR_pos+1; ind < k_bins; ind++){
					Num_Rejected_array[ind]++;
				}
			}

			if (ISR_tag_index != -1 && vect_optimum->Pt()>pt_cut && ISR_OR_NOT == true){  // != S means bb or WI
				for (Int_t ind = 0; ind < k_ISR_pos + 1; ind++){
					if (ISR_tag_index == ISR_jets[entry])
						Num_Tags_array_hpt[ind]++;
					else
						Num_MissTags_array_hpt[ind]++;
				}
				for (Int_t ind = k_ISR_pos+1; ind < k_bins; ind++){
//					Num_Rejected_array_hpt[ind]++;		// Under a certain k_cut, there are not rejected events
				}
			}
		
			Prob_cut = Jet_cut/NumJets;
			if(prob_max >= Prob_cut){
				if (ISR_tag_index == ISR_jets[entry] && ISR_OR_NOT == true)  // != S means bb or WI
					Num_Tags++;
				else
					Num_MissTags++;

				// Cheching MET boosting
				if(vect_optimum->Pt()>pt_cut){
					h_MET_hpt1->Fill(MET);
				}
			}
			else
				Num_Rejected++;

		}

		cout<<"\tprogress = 100%\t";
		cout<<"Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;

	} // End run's for cicle

	/*
	 * Tagging results
	 */

	Int_t Num_Studied = Num_Tags + Num_MissTags + Num_Rejected;
	cout << "\nOverall tagging results" << endl;
	cout << "\tNumber of compared events (between the matching and tagging algorithms) : " << Num_Studied << endl;
	cout << "\tPer. Tags: \t" << ((Double_t)Num_Tags/Num_Studied)*100 << "%" << endl;
	cout << "\tPer. MissTags: \t" << ((Double_t) Num_MissTags/Num_Studied)*100 << "%" << endl;
	cout << "\tPer. Rejected: \t" << ((Double_t) Num_Rejected/Num_Studied)*100 << "%" << endl;

	// Calculating percentages
	for (Int_t ind=0; ind < k_bins; ind++){
		Num_Total_Jets[ind] = Num_Tags_array[ind] + Num_MissTags_array[ind] + Num_Rejected_array[ind];
		Num_Tags_array[ind] = Num_Tags_array[ind]/Num_Total_Jets[ind];
		Num_MissTags_array[ind] = Num_MissTags_array[ind]/Num_Total_Jets[ind];
		Num_Rejected_array[ind] = Num_Rejected_array[ind]/Num_Total_Jets[ind];
		Num_Total_Jets_hpt[ind] = Num_Tags_array_hpt[ind] + Num_MissTags_array_hpt[ind]; // + Num_Rejected_array_hpt[ind];
		Num_Tags_array_hpt[ind] = Num_Tags_array_hpt[ind]/Num_Total_Jets_hpt[ind];
		Num_MissTags_array_hpt[ind] = Num_MissTags_array_hpt[ind]/Num_Total_Jets_hpt[ind];
//		Num_Rejected_array_hpt[ind] = Num_Rejected_array_hpt[ind]/Num_Total_Jets_hpt[ind];
	}

	/*
	 * Writing results
	 */
	Bool_t archivoExiste = false;

	Char_t outNamept[] = "Percenta_hpt-100";
	outNamept[13] = 0x30 + int(pt_cut/100)%10;
	outNamept[14] = 0x30 + int(pt_cut/10)%10;
	outNamept[15] = 0x30 + int(pt_cut)%10;

	string outFileTotal_str = head_folder_results + "Overall_percs" + combination + ".txt";
	Char_t *outFileTotal = (Char_t *) outFileTotal_str.c_str();

	string outFileTotalpt_str = head_folder_results + outNamept + combination + ".txt";
	Char_t *outFileTotalpt = (Char_t *) outFileTotalpt_str.c_str();

	ifstream my_file(outFileTotal);
	if(my_file.good()){
		archivoExiste = true;
	}
	my_file.close();

	ofstream ofs_over(outFileTotal,ios::out);
	if(!archivoExiste){
		// If file already exists
	}

	ofs_over << "# Number of Tags, Misstags and Rejected as a function of k" << endl;
	ofs_over << "# Number of Events " << Num_Total_Jets[0] << endl;
	ofs_over << "# k_cut \t Tags \t MissTags \t Rejected \t Total_Events " << endl;


	for (Int_t ind = 0; ind < k_bins; ind ++){
    	ofs_over << setiosflags(ios::fixed) << setprecision(6) << setw(6) << k_values[ind]
				<< "\t" << Num_Tags_array[ind] << "\t" << Num_MissTags_array[ind] << "\t" << Num_Rejected_array[ind]
				<< "\t" << setprecision(0) << Num_Total_Jets[ind] << endl;
	}

	if (do_pt_cut){
	    ofstream ofs_pt(outFileTotalpt,ios::out);
		ofs_pt << "# Number of Tags, Misstags and Rejected as a function of k. The ISR has pt > " << pt_cut << endl;
		ofs_pt << "# Number of Events " << Num_Total_Jets_hpt[0] << endl;
		ofs_pt << "# k_cut \t Tags \t MissTags \t Total_Events " << endl;
		for (Int_t ind = 0; ind < k_bins; ind ++){
	    	ofs_pt << setiosflags(ios::fixed) << setprecision(6) << setw(6) << k_values[ind]
					<< "\t" << Num_Tags_array_hpt[ind] << "\t" << Num_MissTags_array_hpt[ind] // << "\t" << Num_Rejected_array_hpt[ind]
					<< "\t" << setprecision(0) << Num_Total_Jets_hpt[ind] << endl;
	    }
	    ofs_pt.close();
	}

	if (do_jet_cut){
		Char_t outNameMET[] = "AbsValue_MET_pt_000_k_2.0";
		outNameMET[16] = pt_str[0];
		outNameMET[17] = pt_str[1];
		outNameMET[18] = pt_str[2];
		outNameMET[22] = k_str[0];
		outNameMET[23] = k_str[1];
		outNameMET[24] = k_str[2];

		string outFileMET_str = head_folder_results + outNameMET + combination;
		Char_t *outFileMET = (Char_t *) outFileMET_str.c_str();

		Char_t *outFilehist;
		outFilehist = (Char_t*) malloc(sizeof(char)*512);
		strcpy(outFilehist,outFileMET);
		strcat(outFilehist,".root");

		TFile* hfile = new TFile("histos.root", "RECREATE");
		TCanvas *C = new TCanvas(outFileMET,"MET in a sample with high PT ISR jets",1280,720);
		Present(h_MET,h_MET_hpt1,C,2,"MET [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();
		hfile->Close();

	}

	ofs_over.close();

	cout<<"\nEND :)"<<endl;

	return 0;
}
