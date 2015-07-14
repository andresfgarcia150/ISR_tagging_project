/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

Functions for plotting charts
*/

#include "ROOTFunctions.h"

/*
 * Folder where the pictures will be contained
 * The name of the picture will be FOLDER_GRAPHS + Canvas' name
 * Observe that you can leave FOLDER_GRAPS empty ("") and specify the full path name
 * of the picture in the Canvas' name
 */
#define FOLDER_GRAPHS "";

/*
 * Function to draw two histograms in the same canvas
 * Taken from: Sample code to show Two Statistics Boxes for Two Histograms Side by Side
 * http://hep.uchicago.edu/~okumura/pukiwiki/latest/index.php?Sample%20code%20to%20show%20Two%20Statistics%20Boxes%20for%20Two%20Histograms%20Side%20by%20Side
 */
void Show(TH1* h1, TH1* h2, TCanvas* C){
	Double_t scale1 = 1.0/h1->Integral();
	h1->Scale(scale1);
	h2->SetLineColor(kRed);
	Double_t scale2 = 1.0/h2->Integral();
	h2->Scale(scale2);
	Double_t max1 = h1->GetMaximum();
	Double_t max2 = h2->GetMaximum();
	if ( max1 > max2 ) {
		h1->Draw();
		C->Update();
		TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");

		h2->Draw();
		C->Update();
		TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
		ps2->SetTextColor(kRed);

		//ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.6);
		ps1->SetX1(-4.0); ps1->SetX2(0.5);
		ps2->SetX1NDC(0.65); ps2->SetX2NDC(0.85);

		h1->Draw("");
		h2->Draw("same");
	}
	else {
		h2->Draw();
		C->Update();
		TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
		ps2->SetTextColor(kRed);

		h1->Draw();
		C->Update();
		TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");

		ps1->SetX1NDC(0.15); ps1->SetX2NDC(0.3);
		//ps1->SetX1(-4.0); ps1->SetX2(-2.5);
		ps2->SetX1NDC(0.65); ps2->SetX2NDC(0.85);

		h2->Draw("");
		h1->Draw("same");
	}
}

/*
 * Another function to plot two histograms. It's better than show() and uses THStack class
 * It writes the picture in .png format
 * @param: h1 The first histogram
 * @param: h2 The second histogram
 * @param: C The canvas which stores the information
 * @param: pos1 Position of the statistics of the histograms. pos = {1,2,3,4} => Corners of the screen
 * @param: xaxis Name of the x axis
 * @param: yaxis Name of the y axis
 * @param: Font of the x axis. Recall that fontX = 122 are greek letters
 * @param: Font of the y axis. Recall that fontY = 122 are greek letters
 * @param: Scale of the y axis. (Logarithmic scale or not)
 *
 */
void Present(TH1* h1, TH1* h2, TCanvas* C, Int_t pos, const Char_t * xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false){
	// THStack is a collection of TH1 objects
	THStack* hs = new THStack("hs"," ");

	// Normalizing both histograms
	h1->SetLineColor(kBlue);
	Double_t scale1 = 1.0/h1->Integral();
	h1->Scale(scale1);
	h2->SetLineColor(kRed);
	Double_t scale2 = 1.0/h2->Integral();
	h2->Scale(scale2);

	// Setting the stats boxes
	h1->Draw();
	C->Update();
	TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
	h2->Draw();
	C->Update();
	TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");

	// Log function
	if (log){
		C->SetLogy();
	}

	// Position of the stats
	if (pos == 1){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
		ps2->SetX1NDC(0.28); ps2->SetX2NDC(0.43);
		ps1->SetY1NDC(0.72); ps1->SetY2NDC(0.87);
		ps2->SetY1NDC(0.72); ps2->SetY2NDC(0.87);
	}
	else if (pos == 2){
		ps1->SetX1NDC(0.58); ps1->SetX2NDC(0.73);
		ps2->SetX1NDC(0.74); ps2->SetX2NDC(0.89);
		ps1->SetY1NDC(0.72); ps1->SetY2NDC(0.87);
		ps2->SetY1NDC(0.72); ps2->SetY2NDC(0.87);
	}
	else if (pos == 3){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
		ps2->SetX1NDC(0.28); ps2->SetX2NDC(0.43);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
		ps2->SetY1NDC(0.12); ps2->SetY2NDC(0.27);
	}
	else if (pos == 4){
		ps1->SetX1NDC(0.58); ps1->SetX2NDC(0.73);
		ps2->SetX1NDC(0.74); ps2->SetX2NDC(0.89);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
		ps2->SetY1NDC(0.12); ps2->SetY2NDC(0.27);
	}
	else{
		ps1->SetX1NDC(0.85); ps1->SetX2NDC(1);
		ps2->SetX1NDC(0.85); ps2->SetX2NDC(1);
		ps1->SetY1NDC(0.75); ps1->SetY2NDC(0.90);
		ps2->SetY1NDC(0.59); ps2->SetY2NDC(0.74);
	}
	ps1->SetTextColor(kBlue);
	ps2->SetTextColor(kRed);


	// Adding the two histograms to the stack
	hs->Add(h1);
	hs->Add(h2);

	// Draw both histograms on screen. It's important "nostack" as parameter for avoiding changes in the scale
	hs->Draw("nostack");

	// Setting the label of the axis
	hs->GetXaxis()->SetTitle(xaxis);
	hs->GetYaxis()->SetTitle(yaxis);
	hs->GetXaxis()->SetTitleFont(fontX);
	hs->GetYaxis()->SetTitleFont(fontY);
	hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitleSize(0.05);

	// Setting the title
   TPaveText *pt = new TPaveText(0.4212321,0.94,0.5787679,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(22);
   TText *text = pt->AddText(C->GetTitle());
   text->SetTextSize(0.07);
   pt->Draw();
   C->Modified();

   Char_t nameFile[256] = FOLDER_GRAPHS;
   strcat(nameFile,C->GetName());
   strcat(nameFile,".png");
   C->SaveAs(nameFile);
}

/*
 * Another function to plot three histograms. It's better than show() and uses THStack class
 * It writes the picture in .png format
 * @param: h1 The first histogram
 * @param: h2 The second histogram
 * @param: h3 The third histogram
 * @param: C The canvas which stores the information
 * @param: pos1 Position of the statistics of the histograms. pos = {1,2,3,4} => Corners of the screen
 * @param: xaxis Name of the x axis
 * @param: yaxis Name of the y axis
 * @param: Font of the x axis. Recall that fontX = 122 are greek letters
 * @param: Font of the y axis. Recall that fontY = 122 are greek letters
 * @param: Scale of the y axis. (Logarithmic scale or not)
 */
void Present_3(TH1* h1, TH1* h2, TH1* h3, TCanvas* C, Int_t pos, const Char_t * xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false){
	// THStack is a collection of TH1 objects
	THStack* hs = new THStack("hs"," ");

	// Normalizing both histograms
	h1->SetLineColor(kBlue);
	Double_t scale1 = 1.0/h1->Integral();
	h1->Scale(scale1);
	h2->SetLineColor(kRed);
	Double_t scale2 = 1.0/h2->Integral();
	h2->Scale(scale2);
	h3->SetLineColor(kGreen);
	Double_t scale3 = 1.0/h3->Integral();
	h3->Scale(scale3);

	// Setting the stats boxes
	h1->Draw();
	C->Update();
	TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
	h2->Draw();
	C->Update();
	TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
	h3->Draw();
	C->Update();
	TPaveStats *ps3 = (TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");

	// Log function
	if (log == true){
		C->SetLogy();
	}

	// Position of the stats
	if (pos == 1){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
		ps2->SetX1NDC(0.28); ps2->SetX2NDC(0.43);
		ps3->SetX1NDC(0.45); ps3->SetX2NDC(0.60);
	}
	else if (pos == 2){
		ps1->SetX1NDC(0.42); ps1->SetX2NDC(0.57);
		ps2->SetX1NDC(0.58); ps2->SetX2NDC(0.73);
		ps3->SetX1NDC(0.74); ps3->SetX2NDC(0.89);
	}
	else if (pos == 3){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
		ps2->SetX1NDC(0.28); ps2->SetX2NDC(0.43);
		ps3->SetX1NDC(0.45); ps3->SetX2NDC(0.60);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
		ps2->SetY1NDC(0.12); ps2->SetY2NDC(0.27);
		ps3->SetY1NDC(0.12); ps3->SetY2NDC(0.27);
	}
	else if (pos == 4){
		ps1->SetX1NDC(0.42); ps1->SetX2NDC(0.57);
		ps2->SetX1NDC(0.58); ps2->SetX2NDC(0.73);
		ps3->SetX1NDC(0.74); ps3->SetX2NDC(0.89);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
		ps2->SetY1NDC(0.12); ps2->SetY2NDC(0.27);
		ps3->SetY1NDC(0.12); ps3->SetY2NDC(0.27);
	}
	else{
		ps1->SetX1NDC(0.85); ps1->SetX2NDC(1);
		ps2->SetX1NDC(0.85); ps2->SetX2NDC(1);
		ps3->SetX1NDC(0.85); ps3->SetX2NDC(1);
		ps1->SetY1NDC(0.75); ps1->SetY2NDC(0.90);
		ps2->SetY1NDC(0.59); ps2->SetY2NDC(0.74);
		ps3->SetY1NDC(0.43); ps3->SetY2NDC(0.58);
	}
	ps1->SetTextColor(kBlue);
	ps2->SetTextColor(kRed);
	ps3->SetTextColor(kGreen);

	// Adding the two histograms to the stack
	hs->Add(h1);
	hs->Add(h2);
	hs->Add(h3);

	// Draw both histograms on screen. It's important "nostack" as parameter for avoiding changes in the scale
	hs->Draw("nostack");

	// Setting the label of the axis
	hs->GetXaxis()->SetTitle(xaxis);
	hs->GetYaxis()->SetTitle(yaxis);
	hs->GetXaxis()->SetTitleFont(fontX);
	hs->GetYaxis()->SetTitleFont(fontY);
	hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitleSize(0.05);

	// Setting the title
   TPaveText *pt = new TPaveText(0.4212321,0.94,0.5787679,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(22);
   TText *text = pt->AddText(C->GetTitle());
   text->SetTextSize(0.07);
   pt->Draw();
   C->Modified();

   Char_t nameFile[256] = FOLDER_GRAPHS;
   strcat(nameFile,C->GetTitle());
   strcat(nameFile,".png");
   C->SaveAs(nameFile);
}

/*
 * A Function to plot single histograms
 * It writes the picture in .png format
 * @param: h1 The histogram
 * @param: C The canvas which stores the information
 * @param: pos1 Position of the statistics of the histograms. pos = {1,2,3,4} => Corners of the screen
 * @param: xaxis Name of the x axis
 * @param: yaxis Name of the y axis
 * @param: Font of the x axis. Recall that fontX = 122 are greek letters
 * @param: Font of the y axis. Recall that fontY = 122 are greek letters
 * @param: Scale of the y axis. (Logarithmic scale or not)
 */
void Plot_Single(TH1* h1, TCanvas* C, Int_t pos, const Char_t *xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false){
	// THStack is a collection of TH1 objects
	THStack* hs = new THStack("hs"," ");

	// Normalizing both histograms
	h1->SetLineColor(kBlue);
	Double_t scale1 = 1.0/h1->Integral();
	h1->Scale(scale1);

	// Setting the stats boxes
	h1->Draw();
	C->Update();
	TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");

	// Log function
	if (log){
		C->SetLogy();
	}

	// Position of the stats
	if (pos == 1){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
	}
	else if (pos == 2){
		ps1->SetX1NDC(0.74); ps1->SetX2NDC(0.89);
	}
	else if (pos == 3){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
	}
	else if (pos == 4){
		ps1->SetX1NDC(0.74); ps1->SetX2NDC(0.89);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
	}
	else{
		ps1->SetX1NDC(0.85); ps1->SetX2NDC(1);
		ps1->SetY1NDC(0.75); ps1->SetY2NDC(0.90);
	}
	ps1->SetTextColor(kBlue);


	// Adding the two histograms to the stack
	hs->Add(h1);

	// Draw both histograms on screen. It's important "nostack" as parameter for avoiding changes in the scale
	hs->Draw("nostack");

	// Setting the label of the axis
	hs->GetXaxis()->SetTitle(xaxis);
	hs->GetYaxis()->SetTitle(yaxis);
	hs->GetXaxis()->SetTitleFont(fontX);
	hs->GetYaxis()->SetTitleFont(fontY);
	hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetYaxis()->SetTitleSize(0.05);

	// Setting the title
   TPaveText *pt = new TPaveText(0.4212321,0.94,0.5787679,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(22);
   TText *text = pt->AddText(C->GetTitle());
   text->SetTextSize(0.07);
   pt->Draw();
   C->Modified();

   Char_t nameFile[256] = FOLDER_GRAPHS;
   strcat(nameFile,C->GetTitle());
   strcat(nameFile,".png");
   C->SaveAs(nameFile);
}


/*
 * A Function to plot single 2D histograms
 * It writes the picture in .png format
 * @param: h1 The two dimensional histogram
 * @param: C The canvas which stores the information
 * @param: pos1 Position of the statistics of the histograms. pos = {1,2,3,4} => Corners of the screen
 * @param: xaxis Name of the x axis
 * @param: yaxis Name of the y axis
 * @param: Font of the x axis. Recall that fontX = 122 are greek letters
 * @param: Font of the y axis. Recall that fontY = 122 are greek letters
 * @param: Scale of the y axis. (Logarithmic scale or not)
 */
void Plot_Single_2D(TH2* h1, TCanvas* C, Int_t pos, const Char_t *xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false){

	// Normalizing both histograms
	Double_t scale1 = 1.0/h1->Integral();
	h1->Scale(scale1);
    h1->SetOption("COLZ");

	// Setting the stats boxes
	h1->Draw();
	C->Update();
	TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");

	// Log function
	if (log){
		C->SetLogy();
	}

	// Position of the stats
	if (pos == 1){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
	}
	else if (pos == 2){
		ps1->SetX1NDC(0.74); ps1->SetX2NDC(0.89);
	}
	else if (pos == 3){
		ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.27);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
	}
	else if (pos == 4){
		ps1->SetX1NDC(0.74); ps1->SetX2NDC(0.89);
		ps1->SetY1NDC(0.12); ps1->SetY2NDC(0.27);
	}
	else{
		ps1->SetX1NDC(0.85); ps1->SetX2NDC(1);
		ps1->SetY1NDC(0.75); ps1->SetY2NDC(0.90);
	}
	ps1->SetTextColor(kBlue);


	// Adding the two histograms to the stack
	//h1->Add(h1);

	// Draw both histograms on screen. It's important "nostack" as parameter for avoiding changes in the scale
	//h1->Draw("nostack");

	// Setting the label of the axis
    h1->GetXaxis()->SetTitle(xaxis);
    h1->GetYaxis()->SetTitle(yaxis);
    h1->GetXaxis()->SetTitleFont(fontX);
    h1->GetYaxis()->SetTitleFont(fontY);
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->SetTitleFont(12);
    // Setting the title
    TPaveText *pt = new TPaveText(0.4212321,0.94,0.5787679,0.995,"blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(22);
    //TText *text = pt->AddText(C->GetTitle());
    //text->SetTextSize(0.07);
    pt->Draw();
    C->Modified();

    Char_t nameFile[256] = FOLDER_GRAPHS;
    strcat(nameFile,C->GetTitle());
    strcat(nameFile,".png");
    C->SaveAs(nameFile);
}

