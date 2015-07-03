/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Fisica       -------
-------        Joven Investigador         -------
-------  Andres Felipe Garcia Albarracin  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

This algorithm studies the kinematic properties
of the ISR jets. It reads the results of the
matching algorithm
*/


#include "ROOTFunctions.h"
#include "graphs_Funcs.h"
#include "functions.h"
#include "Rtypes.h"
#include "DelphesFunctions.h"

// Global Variables
const Double_t PI = TMath::Pi();

// Other simulations parameters
const Char_t channel = '_'; // 's' for sTops and '_' for Tops
const Char_t ISR_or_NOT[] = "WI"; // "WI" with ISR, "SI" without (Here it does not make any sense), "bb" bjets production
const Bool_t atServer = true; // True if it is run at the server, false at the university's pc
const Bool_t Matching = true; // True if a matching has been done between MG and Pythia, false otherwise


int main(int argc, char **argv){
	std::cout.precision(4);
	// Counting time
	Double_t initialTime = clock();

	/*
	 * Histograms
	 */
	// All jets
	TH1 *h_numberJet = new TH1F("Number Jets","Number Jets",11,-0.5,10.5);

	// Non Isr jets
	TH1 *h_jet_PT = new TH1F("Jet PT","Jet PT", 201,0.0,600.0);
	TH1 *h_jet_Eta = new TH1F("Jet Eta","Jet Eta", 171,-5.0,5.0);
	TH1 *h_jet_Phi = new TH1F("Jet Phi","Jet Phi", 375,-3.5,3.5);
	TH1 *h_jet_DPhi_MET = new TH1F("Jet - MET Delta_Phi","Jet - MET Delta_Phi",300,0.0,4.0);
	TH1 *h_jet_DPhi_MET_hpt = new TH1F("Jet - MET Delta_Phi_hpt","Jet - MET Delta_Phi_hpt",300,0.0,4.0);
	TH1 *h_jet_MT = new TH1F("Jet Transverse mass","Jet Transverse Mass",201,0.0,600.0);
	TH1 *h_jet_Delta_PT = new TH1F("Jet Delta-PT","Non ISR Delta-PT", 201,0.0,300.0);
	TH1 *h_jet_PT_HT = new TH1F("Jet PT-HT ratio","Jet PT-HT ratio",201,-0.0025,1.0025);
	TH1 *h_jet_PT_over_PT_others = new TH1F("Jet PT/PT_others","Jet PT/PT_others",401,-0.0025,2.0025);
	TH1 *h_jet_Eta_over_Eta_others = new TH1F("Jet Eta/Eta_others","Jet Eta/Eta_others",401,-0.0025,2.0025);
	TH1 *h_jet_DPhi_over_Phi_others = new TH1F("Jet Phi/Phi_others","Jet Phi/Phi_others",401,-0.0025,2.0025);
	TH1 *h_jet_Delta_Eta = new TH1F("Jet Delta-Eta","Jet Delta-Eta", 171,0.0,5.0);
	TH1 *h_jet_DPhi_MET_other = new TH1F("Jet - MET Delta_Phi other","Jet - MET Delta_Phi other",300,0.0,4.0);
	TH1 *h_jet_multiplicity = new TH1F("Jet - Multiplicity","Jet - Multiplicity",101,-0.5,100.5);
	TH1 *h_jet_DeltaR = new TH1F ("Jet - Delta_R","Jet - Delta_R",201,-0.0025,0.8025);
	TH1 *h_jet_Delta_PT_leading = new TH1F("Delta PT: leading - Jet","Delta PT: leading - Jet", 201,0.0,600.0);
	TH1 *h_jet_Delta_Eta_leading = new TH1F("Delta Eta: Jet - leading","Delta Eta: Jet - leading", 171,0.0,8.0);

	TH2 *h2_jet_PTEta=new TH2F("Non_ISR_Jet_PT_Eta","Non ISR Jet PT Vs. Eta",201,-1.25,501.25,201,-4.02,4.02);

	// ISR jets
	TH1 *h_ISR_PT = new TH1F("ISR PT","ISR PT", 201,0.0,600.0);
	TH1 *h_ISR_Eta = new TH1F("ISR Eta","ISR Eta", 171,-5.0,5.0);
	TH1 *h_ISR_Phi = new TH1F("ISR Phi","ISR Phi", 375,-3.5,3.5);
	TH1 *h_ISR_DPhi_MET = new TH1F("ISR - MET Delta_Phi","ISR - MET Delta_Phi",300,0.0,4.0);
	TH1 *h_ISR_DPhi_MET_hpt = new TH1F("ISR - MET Delta_Phi_hpt","ISR - MET Delta_Phi_hpt",300,0.0,4.0);
	TH1 *h_ISR_MT = new TH1F("ISR Transverse mass","ISR Transverse Mass",201,0.0,600.0);
	TH1 *h_ISR_Delta_PT = new TH1F("ISR Delta-PT","ISR Delta-PT", 201,0.0,300.0);
	TH1 *h_ISR_PT_HT = new TH1F("ISR PT-HT ratio","ISR PT-HT ratio",201,-0.0025,1.0025);
	TH1 *h_ISR_PT_over_PT_others = new TH1F("ISR PT/PT_others","ISR PT/PT_others",401,-0.0025,2.0025);
	TH1 *h_ISR_Eta_over_Eta_others = new TH1F("ISR Eta/Eta_others","ISR Eta/Eta_others",401,-0.0025,2.0025);
	TH1 *h_ISR_DPhi_over_Phi_others = new TH1F("ISR Phi/Phi_others","ISR Phi/Phi_others",401,-0.0025,2.0025);
	TH1 *h_ISR_Delta_Eta = new TH1F("ISR Delta-Eta","ISR Delta-Eta", 171,0.0,5.0);
	TH1 *h_ISR_DPhi_MET_other = new TH1F("ISR - MET Delta_Phi other","ISR - MET Delta_Phi other",300,0.0,4.0);
	TH1 *h_ISR_multiplicity = new TH1F("ISR - Multiplicity","ISR - Multiplicity",101,-0.5,100.5);
	TH1 *h_ISR_DeltaR = new TH1F ("ISR - Delta_R","ISR - Delta_R",201,-0.0025,0.8025);
	TH1 *h_ISR_Delta_PT_leading = new TH1F("Delta PT: leading - ISR","Delta PT: leading - ISR", 201,0.0,600.0);
	TH1 *h_ISR_Delta_Eta_leading = new TH1F("Delta Eta: ISR - leading","Delta Eta: ISR - leading", 171,0.0,8.0);

	TH2 *h2_ISR_PTEta=new TH2F("ISR_Jet_PT_Eta","ISR Jet PT Vs. Eta",201,-1.25,501.25,201,-4.02,4.02);

	// MET
	TH1 *h_MET = new TH1F("Missing ET","Missing ET",200,0,600);
	TH1 *h_MET_hpt1 = new TH1F("Missing ET high_ISR_pt-1","Missing ET high_ISR_pt-1",200,0.0,600.0);
	TH1 *h_MET_hpt2 = new TH1F("Missing ET high_ISR_pt-2","Missing ET high_ISR_pt-2",200,0.0,600.0);
	TH1 *h_MET_hpt3 = new TH1F("Missing ET high_ISR_pt-3","Missing ET high_ISR_pt-3",200,0.0,600.0);
	TH1 *h_MET_hpt4 = new TH1F("Missing ET high_ISR_pt-4","Missing ET high_ISR_pt-4",200,0.0,600.0);

	TH2 *h2_dif_PTEta=new TH2F("FSR_ISR_Jet_PT_Eta_Difference","Difference between FSR and ISR Jet PT Vs. Eta distributions",201,-1.25,501.25,201,-4.02,4.02);
	TH2 *h2_dif_lead_PTEta=new TH2F("Lead_ISR_Jet_PT_Eta_Difference","Difference between Lead and ISR Jet PT Vs. Eta distributions",201,-1.25,501.25,201,-4.02,4.02);

	// Leading PT
	TH1 *h_leading_PT = new TH1F("Leading PT","Leading PT", 201,0.0,600.0);
	TH1 *h_leading_MT = new TH1F("Leading Transverse mass","Leading Transverse Mass",201,0.0,600.0);
	TH1 *h_leading_Eta = new TH1F("Leading Eta","Leading Eta", 171,-5.0,5.0);
	TH1 *h_leading_DPhi_MET = new TH1F("Leading - MET Delta_Phi","Leading - MET Delta_Phi",300,0.0,4.0);

	TH2 *h2_leading_PTEta=new TH2F("Leading_Jet_PT_Eta","Leading Jet PT Vs. Eta",201,-1.25,501.25,201,-4.02,4.02);

	// Other variables
	TH1 *h_HT = new TH1F("HT","HT",201,0.0,600.0);
	TH1 *h_HT_R1 = new TH1F("HT_R1","HT_R1",51,-0.01,1.01);
	TH1 *h_HT_R2 = new TH1F("HT_R2","HT_R2",51,-0.01,1.01);

	// B tagging
	TH1 *h_BTag = new TH1F("BTag","BTag",5,-0.5,4.5);
	TH1 *h_BTag_PT = new TH1F("BTag PT","BTag PT", 201,0.0,600.0);
	TH1 *h_BTag_Eta = new TH1F("BTag Eta","BTag Eta", 171,-5.0,5.0);
	TH1 *h_BTag_DPhi_MET = new TH1F("BTag - MET Delta_Phi","BTag - MET Delta_Phi",300,0.0,4.0);
	TH1 *h_BTags_per_Event = new TH1F("BTags per event","BTags per event",5,-0.5,4.5);

	// Further analysis
	TH1 *h_ISR_PT_comp = new TH1F("ISR PT for comparison","ISR PT for comparison with histo", 20,0.0,800.0);
	TH1 *h_ISR_Eta_comp = new TH1F("ISR Eta for comparison","ISR Eta for comparison with histo", 20,-4.2,4.2);
	TH1 *h_ISR_DPhi_MET_comp = new TH1F("ISR Phi for comparison","ISR Phi for comparison with histo", 20,0,PI);

	// To check the histograms' creation
	TH1 *hist_ISR_PT = new TH1F("ISR PT comp","ISR PT comp", 20,0.0,800.0);
	TH1 *hist_ISR_Abs_Eta = new TH1F("ISR Abs Eta comp","ISR Abs Eta comp", 20,0.0,5.2);
	TH1 *hist_ISR_DPhi_MET = new TH1F("ISR Delta Phi comp","ISR Delta Phi comp", 20,0.0,PI);
	TH1 *hist_ISR_PT_ratio = new TH1F("ISR PT/PT_others comp","ISR PT/PT_others comp",20,0.0,8.0);
	TH1 *hist_ISR_Delta_Eta = new TH1F("ISR Delta-Eta comp","ISR Delta-Eta comp", 20,0.0,7.0);
	TH1 *hist_ISR_DPhi_MET_other = new TH1F("ISR - MET Delta_Phi other comp","ISR - MET Delta_Phi other comp",20,0.0,PI);
	TH1 *hist_ISR_Delta_PT_leading = new TH1F("Delta PT: leading - ISR comp","Delta PT: leading - ISR comp", 20,0.0,500.0);
	TH1 *hist_ISR_Delta_Eta_leading = new TH1F("Delta Eta: ISR - leading comp","Delta Eta: ISR - leading comp", 20,0.0,6.5);
	TH1 *hist_jet_PT = new TH1F("Jet PT comp","Jet PT comp", 20,0.0,800.0);
	TH1 *hist_jet_Abs_Eta = new TH1F("Jet Abs Eta comp","Jet Abs Eta comp", 20,0.0,5.2);
	TH1 *hist_jet_DPhi_MET = new TH1F("Jet Delta Phi comp","Jet Delta Phi comp", 20,0.0,PI);
	TH1 *hist_jet_PT_ratio = new TH1F("Jet PT/PT_others comp","Jet PT/PT_others comp",20,0.0,7.0);
	TH1 *hist_jet_Delta_Eta = new TH1F("Jet Delta-Eta comp","Jet Delta-Eta comp", 20,0.0,8.0);
	TH1 *hist_jet_DPhi_MET_other = new TH1F("Jet - MET Delta_Phi other comp","Jet - MET Delta_Phi other comp",20,0.0,PI);
	TH1 *hist_jet_Delta_PT_leading = new TH1F("Delta PT: leading - Jet comp","Delta PT: leading - Jet comp", 20,0.0,500.0);
	TH1 *hist_jet_Delta_Eta_leading = new TH1F("Delta Eta: Jet - leading comp","Delta Eta: Jet - leading comp", 20,0.0,6.5);

	for(int iRun = 1; iRun < 11; iRun ++){
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
		head_folder[0] = channel;
		head_folder[13] = ISR_or_NOT[0];
		head_folder[14] = ISR_or_NOT[1];

		Char_t current_folder[] = "_Tops_MG_1K_AG_WI_003/";
		current_folder[0] = channel;
		current_folder[15] = ISR_or_NOT[0];
		current_folder[16] = ISR_or_NOT[1];

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
		TRefArray array_temp;

		// Temporary variables
		Double_t MET = 0.0; // Missing transverse energy
		Double_t delta_phi = 0.0; // difference between the phi angle of MET and the jet
		Double_t transverse_mass = 0.0; // Transverse mass
		Double_t HT = 0.0; // Sum of jets' PT
		Double_t HT_R1 = 0.0; // Sum of jets' PT which are in the same hemisphere of the ISR jet hemisphere
		Double_t HT_R2 = 0.0; // Sum of jets' PT which are in the opposite hemisphere of the ISR jet hemisphere
		Double_t ISR_Eta = 0.0; // Pseudorapidity of the ISR jet
		Int_t number_Btags = 0; // Number of B jets per event
		Int_t ISR_Btags = 0; // Number of BTags which are also ISR jets
		Double_t delta_PT_jet = 0.0; // |PT-<PT>|
		Double_t PT_sum = 0.0; // sum(PT)
		Double_t PT_aver = 0.0; // <PT>
		Double_t Delta_eta_aver = 0.0; // sum_i|eta-eta_i|/(Nj-1)
		Double_t Delta_phi_sum = 0.0; // sum delta_phi
		Double_t Delta_phi_other_jets = 0.0; // Average of delta phi of other jets
		Double_t PT_ratio = 0.0; // PT/PT_others
		Double_t Eta_ratio = 0.0; // Eta/Eta_others
		Double_t Eta_sum = 0.0; // sum(Eta)
		Double_t Delta_R = 0.0; // Size of the jet
		Double_t Delta_phi_ratio = 0.0; // Delta_phi/Delta_phi_others
		Double_t Delta_PT_leading = 0.0; // PT - PT_leading
		Double_t Delta_Eta_leading = 0.0; // |Eta - Eta_leading|

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
		head_folder_binary[0] = channel;
		head_folder_binary[13] = ISR_or_NOT[0];
		head_folder_binary[14] = ISR_or_NOT[1];

		Char_t matching_name[] = "ISR_jets_Tops_WI_003.bn";
		matching_name[8] = channel;
		matching_name[14] = ISR_or_NOT[0];
		matching_name[15] = ISR_or_NOT[1];

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

		// Jet with greatest PT
		Double_t PT_max = 0;
		Int_t posLeadingPT = -1;
		Int_t ISR_greatest_PT = 0;
		Double_t MT_leading_jet = 0.0; // Transverse mass

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
			h_MET->Fill(MET);

			NumJets=branchJet->GetEntries();
			h_numberJet->Fill(NumJets);

			// checking the ISR
			if (ISR_jets[entry] == -1 || NumJets < 3)
				continue;

			PT_max = 0;
			posLeadingPT = -1;
			HT = 0;
			HT_R1 = 0;
			HT_R2 = 0;
			number_Btags = 0;

			delta_PT_jet = 0.0;
			PT_aver = 0.0;
			PT_sum = 0.0;
			Delta_eta_aver = 0.0;
			Delta_phi_sum = 0.0;
			Delta_phi_other_jets = 0.0;
			Delta_phi_ratio = 0.0;
			Delta_PT_leading = 0.0;
			Delta_Eta_leading = 0.0;

			PT_ratio = 0.0;
			Eta_ratio = 0.0;
			Eta_sum = 0.0;

			Delta_R = 0.0;

			if (ISR_jets[entry] >= NumJets){
				cout << "Error en el matching" << endl;
				return 1;
			}

			// Preliminary for. It is used to calculate PT_aver and Delta_phi_sum
			for (Int_t iJet = 0; iJet<NumJets; iJet++){
				currentJet = (Jet*) branchJet->At(iJet);
				vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
				delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);
				PT_sum += vect_currentJet->Pt();
				Eta_sum += vect_currentJet->Eta();
				Delta_phi_sum += delta_phi;
				// HT
				HT += vect_currentJet->Pt();
				// HT ratios
				if((vect_currentJet->Eta()*ISR_Eta) > 0)
					HT_R1 += vect_currentJet->Pt();
				else
					HT_R2 += vect_currentJet->Pt();
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

			// ISR jet
			currentJet = (Jet*) branchJet->At(ISR_jets[entry]);
			vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
			ISR_Eta = vect_currentJet->Eta();

			for (Int_t iJet = 0; iJet<NumJets; iJet++){
				currentJet = (Jet*) branchJet->At(iJet);
				vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
				delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);
				transverse_mass = sqrt(2*vect_currentJet->Pt()*MET*(1-cos(delta_phi)));

				// Correlated variables
				delta_PT_jet = TMath::Abs(vect_currentJet->Pt()-PT_aver);
				Delta_phi_other_jets = (Delta_phi_sum-delta_phi)/(NumJets-1);
				PT_ratio = vect_currentJet->Pt()*(NumJets-1)/(PT_sum-vect_currentJet->Pt());
				Eta_ratio = vect_currentJet->Eta()*(NumJets-1)/(Eta_sum-vect_currentJet->Eta());
				Delta_phi_ratio = delta_phi*(NumJets-1)/(Delta_phi_sum-delta_phi);

				Delta_Eta_leading = TMath::Abs(vect_currentJet->Eta()-vect_leading->Eta());
				Delta_PT_leading = vect_leading->Pt()-vect_currentJet->Pt();

				Delta_eta_aver = 0.0;
				// For cycle used to calculate Delta_eta_aver
				for(Int_t iJet2 = 0; iJet2<NumJets; iJet2++){
					auxJet = (Jet*) branchJet->At(iJet2);
					vect_auxJet->SetPtEtaPhiM(auxJet->PT,auxJet->Eta,auxJet->Phi,auxJet->Mass);
					if (iJet2 != iJet) Delta_eta_aver += TMath::Abs(vect_auxJet->Eta()-vect_currentJet->Eta());
				}
				Delta_eta_aver = Delta_eta_aver/(NumJets-1);
				Delta_R = sqrt(pow(currentJet->DeltaEta,2)+pow(currentJet->DeltaPhi,2));

				// Multiplicity
				array_temp = (TRefArray) currentJet->Constituents;

				if (iJet != ISR_jets[entry]){ // Non ISR
					h_jet_PT->Fill(vect_currentJet->Pt());
					h_jet_Eta->Fill(vect_currentJet->Eta());
					h_jet_Phi->Fill(vect_currentJet->Phi());
					h_jet_DPhi_MET->Fill(delta_phi);
					h_jet_MT->Fill(transverse_mass);
					h_jet_Delta_PT->Fill(delta_PT_jet);
					h_jet_Delta_Eta->Fill(Delta_eta_aver);
					h_jet_DPhi_MET_other->Fill(Delta_phi_other_jets);
					h_jet_PT_HT->Fill(vect_currentJet->Pt()/HT);
					h_jet_multiplicity->Fill(array_temp.GetEntries());
					h_jet_PT_over_PT_others->Fill(PT_ratio);
					h_jet_Eta_over_Eta_others->Fill(Eta_ratio);
					h_jet_DeltaR->Fill(Delta_R);
					h_jet_DPhi_over_Phi_others->Fill(Delta_phi_ratio);
					h_jet_Delta_PT_leading->Fill(Delta_PT_leading);
					h_jet_Delta_Eta_leading->Fill(Delta_Eta_leading);
					if (vect_currentJet->Pt()>240)
						h_jet_DPhi_MET_hpt->Fill(delta_phi);
					h2_jet_PTEta->Fill(vect_currentJet->Pt(),vect_currentJet->Eta());

					// For testing creating histo
					hist_jet_PT->Fill(vect_currentJet->Pt());
					hist_jet_Abs_Eta->Fill(TMath::Abs(vect_currentJet->Eta()));
					hist_jet_DPhi_MET->Fill(delta_phi);
					hist_jet_PT_ratio->Fill(PT_ratio);
					hist_jet_Delta_Eta->Fill(Delta_eta_aver);
					hist_jet_DPhi_MET_other->Fill(Delta_phi_other_jets);
					hist_jet_Delta_PT_leading->Fill(Delta_PT_leading);
					hist_jet_Delta_Eta_leading->Fill(Delta_Eta_leading);
				}
				else{ //ISR
					h_ISR_PT->Fill(vect_currentJet->Pt());
					h_ISR_Eta->Fill(vect_currentJet->Eta());
					h_ISR_Phi->Fill(vect_currentJet->Phi());
					h_ISR_DPhi_MET->Fill(delta_phi);
					h_ISR_Eta_comp->Fill(vect_currentJet->Eta());
					h_ISR_PT_comp->Fill(vect_currentJet->Pt());
					h_ISR_DPhi_MET_comp->Fill(delta_phi);
					h_ISR_Delta_PT->Fill(delta_PT_jet);
					h_ISR_Delta_Eta->Fill(Delta_eta_aver);
					h_ISR_DPhi_MET_other->Fill(Delta_phi_other_jets);
					h_ISR_PT_HT->Fill(vect_currentJet->Pt()/HT);
					h_ISR_multiplicity->Fill(array_temp.GetEntries());
					h_ISR_PT_over_PT_others->Fill(PT_ratio);
					h_ISR_Eta_over_Eta_others->Fill(Eta_ratio);
					h_ISR_DeltaR->Fill(Delta_R);
					h_ISR_DPhi_over_Phi_others->Fill(Delta_phi_ratio);
					h_ISR_Delta_PT_leading->Fill(Delta_PT_leading);
					h_ISR_Delta_Eta_leading->Fill(Delta_Eta_leading);
					if (vect_currentJet->Pt()>120)
						h_MET_hpt1->Fill(MET);
					if (vect_currentJet->Pt()>200)
						h_MET_hpt2->Fill(MET);
					if (vect_currentJet->Pt()>240){
						h_MET_hpt3->Fill(MET);
						h_ISR_DPhi_MET_hpt->Fill(delta_phi);
					}
					if (vect_currentJet->Pt()>300)
						h_MET_hpt4->Fill(MET);
					h2_ISR_PTEta->Fill(vect_currentJet->Pt(),vect_currentJet->Eta());
					// Transverse mass
					h_ISR_MT->Fill(transverse_mass);

					// For testing creating histo
					hist_ISR_PT->Fill(vect_currentJet->Pt());
					hist_ISR_Abs_Eta->Fill(TMath::Abs(vect_currentJet->Eta()));
					hist_ISR_DPhi_MET->Fill(delta_phi);
					hist_ISR_PT_ratio->Fill(PT_ratio);
					hist_ISR_Delta_Eta->Fill(Delta_eta_aver);
					hist_ISR_DPhi_MET_other->Fill(Delta_phi_other_jets);
					hist_ISR_Delta_PT_leading->Fill(Delta_PT_leading);
					hist_ISR_Delta_Eta_leading->Fill(Delta_Eta_leading);
				}

				// BTag
				h_BTag->Fill(currentJet->BTag);
				if (currentJet->BTag == 1){ // The current jet is B Tagged
					h_BTag_PT->Fill(vect_currentJet->Pt());
					h_BTag_Eta->Fill(vect_currentJet->Eta());
					h_BTag_DPhi_MET->Fill(delta_phi);
					number_Btags++;

					if (iJet == ISR_jets[entry]){ // If the ISR jet is also a B jet
						ISR_Btags++;
					}
				}
			}

			// Jet with greatest PT
			if (posLeadingPT != -1){
				h_leading_PT->Fill(PT_max);
				if(posLeadingPT == ISR_jets[entry]) ISR_greatest_PT++;

				currentJet = (Jet*) branchJet->At(posLeadingPT);
				vect_currentJet->SetPtEtaPhiM(currentJet->PT,currentJet->Eta,currentJet->Phi,currentJet->Mass);
				delta_phi = deltaAng(vect_currentJet->Phi(), METpointer->Phi);
				MT_leading_jet = sqrt(2*vect_currentJet->Pt()*MET*(1-cos(delta_phi)));
				h_leading_MT->Fill(MT_leading_jet);

				h_leading_Eta->Fill(vect_currentJet->Eta());
				h_leading_DPhi_MET->Fill(delta_phi);

				h2_leading_PTEta->Fill(vect_currentJet->Pt(),vect_currentJet->Eta());
			}

			// HT
			if (1 < HT_R1/HT || 1 < HT_R2/HT){
				cout << "Error en el evento: " << entry << endl;
				cout << "HT: " << HT << "\tHT_R1: " << HT_R1 << "\tHT_R2: " << HT_R2 << endl;
				return 1;
			}

			h_HT->Fill(HT);
			h_HT_R1->Fill(HT_R1/HT);
			h_HT_R2->Fill(HT_R2/HT);
			h_BTags_per_Event->Fill(number_Btags);

		}

		cout<<"progress = 100%\t";
		cout<< "Time :"<< (clock()-initialTime)/double_t(CLOCKS_PER_SEC)<<"s"<<endl;
		cout<< "Percentage of events where the ISR jet is the jet with greatest PT: " << (Double_t) (ISR_greatest_PT*100)/numberOfEntries << "%\n";
		cout<< "Percentage of events where the ISR jet is tagged as Bjet: " << (Double_t) (ISR_Btags*100)/numberOfEntries << "%\n";

	} // End run's for cicle

	TFile* hfile = new TFile("./histos/histos.root", "RECREATE");
	h_jet_DPhi_MET->Write();
	h_jet_Eta->Write();
	h_jet_PT->Write();
	h_jet_Phi->Write();
	h_jet_MT->Write();
	h_jet_Delta_PT->Write();
	h_jet_Delta_Eta->Write();
	h_jet_DPhi_MET_other->Write();
	h_jet_PT_HT->Write();
	h_jet_multiplicity->Write();
	h_jet_PT_over_PT_others->Write();
	h_jet_Eta_over_Eta_others->Write();
	h_jet_DeltaR->Write();
	h_jet_DPhi_over_Phi_others->Write();
	h_jet_Delta_Eta_leading->Write();
	h_jet_Delta_PT_leading->Write();

	h_ISR_DPhi_MET->Write();
	h_ISR_Eta->Write();
	h_ISR_PT->Write();
	h_ISR_Phi->Write();
	h_ISR_MT->Write();
	h_ISR_Delta_PT->Write();
	h_ISR_Delta_Eta->Write();
	h_ISR_DPhi_MET_other->Write();
	h_ISR_PT_HT->Write();
	h_ISR_multiplicity->Write();
	h_ISR_PT_over_PT_others->Write();
	h_ISR_Eta_over_Eta_others->Write();
	h_ISR_DeltaR->Write();
	h_ISR_DPhi_over_Phi_others->Write();
	h_ISR_Delta_Eta_leading->Write();
	h_ISR_Delta_PT_leading->Write();

	h_MET->Write();
	h_MET_hpt1->Write();
	h_MET_hpt2->Write();
	h_MET_hpt3->Write();

	h_leading_MT->Write();
	h_leading_PT->Write();
	h_leading_Eta->Write();
	h_leading_DPhi_MET->Write();

	h_HT->Write();
	h_HT_R1->Write();
	h_HT_R2->Write();

	h_numberJet->Write();

	h_BTag->Write();
	h_BTag_PT->Write();
	h_BTag_Eta->Write();
	h_BTag_DPhi_MET->Write();
	h_BTags_per_Event->Write();

	h2_ISR_PTEta->Write();
	h2_jet_PTEta->Write();
	h2_dif_PTEta->Add(h2_ISR_PTEta,h2_jet_PTEta,1,-1);
	h2_dif_PTEta->Write();

	h2_dif_lead_PTEta->Add(h2_ISR_PTEta,h2_leading_PTEta,1,-1);
	h2_dif_lead_PTEta->Write();

	{
		TCanvas *C = new TCanvas("Eta","Pseudorapidity",1280,720);
		Present(h_ISR_Eta,h_jet_Eta,C,1,"h","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Eta ISR vs BTag","Pseudorapidity ISR vs BTag",1280,720);
		Present(h_ISR_Eta,h_BTag_Eta,C,1,"h","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Eta ISR vs Leading","Pseudorapidity ISR vs Leading",1280,720);
		Present(h_ISR_Eta,h_leading_Eta,C,1,"h","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Transverse momentum","Transverse momentum",1280,720);
		Present(h_ISR_PT,h_jet_PT,C,2,"PT [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Transverse momentum ISR vs Leading","Transverse momentum ISR vs Leading",1280,720);
		Present(h_ISR_PT,h_leading_PT,C,2,"PT [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Transverse momentum ISR vs B_Tag","Transverse momentum ISR vs B_Tag",1280,720);
		Present(h_ISR_PT,h_BTag_PT,C,2,"PT [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Transverse momentum ISR, B_Tag, Leading","Transverse momentum ISR, B_Tag, Leading",1280,720);
		Present_3(h_ISR_PT,h_BTag_PT,h_leading_PT,C,2,"PT [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Transverse momentum ISR, B_Tag, Leading LOG","Transverse momentum ISR, B_Tag, Leading LOG",1280,720);
		Present_3(h_ISR_PT,h_BTag_PT,h_leading_PT,C,2,"PT [GeV]","Num. Jets / Total",12,12,true);
		C->Write();
		C->Close();

		C = new TCanvas("Transverse mass Leading vs ISR Jet","Transverse mass Leading vs ISR Jet",1280,720);
		Present(h_ISR_MT,h_leading_MT,C,2,"MT [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Transverse mass ISR vs Jet","Transverse mass ISR vs Jet",1280,720);
		Present(h_ISR_MT,h_jet_MT,C,2,"MT [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Phi","Phi",1280,720);
		Present(h_ISR_Phi,h_jet_Phi,C,3,"f","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Delta Phi - Jet - MET","Delta Phi - Jet - MET",1280,720);
		Present(h_ISR_DPhi_MET,h_jet_DPhi_MET,C,3,"Df","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Delta Phi - Jet - MET - Btag","Delta Phi - Jet - MET - Btag",1280,720);
		Present(h_ISR_DPhi_MET,h_BTag_DPhi_MET,C,3,"Df","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Delta Phi - Jet - MET - leading","Delta Phi - Jet - MET - leading",1280,720);
		Present(h_ISR_DPhi_MET,h_leading_DPhi_MET,C,1,"Df","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("MET > 120","MET > 120",1280,720);
		Present(h_MET,h_MET_hpt1,C,2,"MET","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("MET > 200","MET > 200",1280,720);
		Present(h_MET,h_MET_hpt2,C,2,"MET","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("MET > 240","MET > 240",1280,720);
		Present(h_MET,h_MET_hpt3,C,2,"MET","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("HT ratio comparison","HT ratio comparison",1280,720);
		Present(h_HT_R1,h_HT_R2,C,2,"HT","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("PT vs ETA - ISR","PT vs ETA - ISR",1280,720);
		Plot_Single_2D(h2_ISR_PTEta,C,2, "PT [GeV]", "h", 12, 122);
		C->Write();
		C->Close();

		C = new TCanvas("PT vs ETA - Jet","PT vs ETA - Jet",1280,720);
		Plot_Single_2D(h2_jet_PTEta,C,2, "PT [GeV]", "h", 12, 122);
		C->Write();
		C->Close();

		C = new TCanvas("PT vs ETA - Diff with any jet","PT vs ETA - Diff with any jet",1280,720);
		Plot_Single_2D(h2_dif_PTEta,C,2, "PT [GeV]", "h", 12, 122);
		C->Write();
		C->Close();

		C = new TCanvas("PT vs ETA - leading","PT vs ETA - leading",1280,720);
		Plot_Single_2D(h2_leading_PTEta,C,2, "PT [GeV]", "h", 12, 122);
		C->Write();
		C->Close();

		C = new TCanvas("PT vs ETA - Diff with leading","PT vs ETA - Diff with leading",1280,720);
		Plot_Single_2D(h2_dif_lead_PTEta,C,2, "PT [GeV]", "h", 12, 122);
		C->Write();
		C->Close();

		C = new TCanvas("HT","HT",1280,720);
		Plot_Single(h_HT,C,2, "HT [GeV]", "Num. Jets / Total", 12, 12);
		C->Write();
		C->Close();

		C = new TCanvas("Number_of_B_Tags","Number of B Tags",1280,720);
		Plot_Single(h_BTags_per_Event,C,2, "B Tags / event", "Num. Jets / Total", 12, 12);
		C->Write();
		C->Close();

		C = new TCanvas("Jet_multiplitcity","Jet multiplicity",1280,720);
		Present(h_ISR_multiplicity,h_jet_multiplicity,C,2,"Tracks","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Delta_R_-_Jet_size","Delta R - Jet Size",1280,720);
		Present(h_ISR_DeltaR,h_jet_DeltaR,C,1,"Delta_R","Num. Jets / Total");
		C->Write();
		C->Close();

		// Correlated variables
		C = new TCanvas("Cor_Delta_PT_Jet", "Delta PT jet",1280,720);
		Present(h_ISR_Delta_PT,h_jet_Delta_PT,C,2,"PT [GeV]","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Cor_PT_proportion","PT proportion",1280,720);
		Present(h_ISR_PT_HT,h_jet_PT_HT,C,2,"PT/HT","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Cor_Delta_Eta_Average","Delta Eta Average",1280,720);
		Present(h_ISR_Delta_Eta,h_jet_Delta_Eta,C,2,"Dh","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Cor_Delta_Phi_Jet_MET_other_jets","Delta Phi - Jet MET - other jets",1280,720);
		Present(h_ISR_DPhi_MET_other,h_jet_DPhi_MET_other,C,2,"Df","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Cor_PT_over_<PT_other>","PT/<PT_other>",1280,720);
		Present(h_ISR_PT_over_PT_others,h_jet_PT_over_PT_others,C,2,"PT/<PT>","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Cor_Eta_over_<Eta_other>","Eta/<Eta_other>",1280,720);
		Present(h_ISR_Eta_over_Eta_others,h_jet_Eta_over_Eta_others,C,3,"h/<h>","Num. Jets / Total",122);
		C->Write();
		C->Close();

		C = new TCanvas("Cor_Delta_Phi_over_<Delta_Phi_other>","Delta_Phi/<Delta_Phi_other>",1280,720);
		Present(h_ISR_DPhi_over_Phi_others,h_jet_DPhi_over_Phi_others,C,3,"Df/<Df>","Num. Jets / Total",122);
		C->Write();
		C->Close();

		// Comparison with the leading Jet
		C = new TCanvas("Leading_Delta_PT","Delta PT: PT_leading-PT",1280,720);
		Present(h_ISR_Delta_PT_leading,h_jet_Delta_PT_leading,C,2,"(PT_leading - PT)","Num. Jets / Total");
		C->Write();
		C->Close();

		C = new TCanvas("Leading_Delta_Eta","Delta Eta: |Eta-Eta_leading|",1280,720);
		Present(h_ISR_Delta_Eta_leading,h_jet_Delta_Eta_leading,C,2,"|Eta - Eta_leading|","Num. Jets / Total");
		C->Write();
		C->Close();

	}

	hfile->Close();

	TFile* hfile2 = new TFile("./histos/histos2.root", "RECREATE");
	h_ISR_PT_comp->Write();
	h_ISR_Eta_comp->Write();
	h_ISR_DPhi_MET_comp->Write();

	hist_ISR_PT->Write();
	hist_ISR_Abs_Eta->Write();
	hist_ISR_DPhi_MET->Write();
	hist_ISR_PT_ratio->Write();
	hist_ISR_Delta_Eta->Write();
	hist_ISR_DPhi_MET_other->Write();
	hist_ISR_Delta_PT_leading->Write();
	hist_ISR_Delta_Eta_leading->Write();

	hist_jet_PT->Write();
	hist_jet_Abs_Eta->Write();
	hist_jet_DPhi_MET->Write();
	hist_jet_PT_ratio->Write();
	hist_jet_Delta_Eta->Write();
	hist_jet_DPhi_MET_other->Write();
	hist_jet_Delta_PT_leading->Write();
	hist_jet_Delta_Eta_leading->Write();

	hfile2->Close();

	return 0;
}
