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

/*
 * Function to draw two histograms in the same canvas
 * Taken from: Sample code to show Two Statistics Boxes for Two Histograms Side by Side
 * http://hep.uchicago.edu/~okumura/pukiwiki/latest/index.php?Sample%20code%20to%20show%20Two%20Statistics%20Boxes%20for%20Two%20Histograms%20Side%20by%20Side
 */
void Show(TH1* h1, TH1* h2, TCanvas* C);

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
 */
void Present(TH1* h1, TH1* h2, TCanvas* C, Int_t pos, const Char_t * xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false);

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
void Plot_Single(TH1* h1, TCanvas* C, Int_t pos, const Char_t *xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false);


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
void Plot_Single_2D(TH2* h1, TCanvas* C, Int_t pos, const Char_t *xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false);

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
void Present_3(TH1* h1, TH1* h2, TH1* h3, TCanvas* C, Int_t pos, const Char_t * xaxis, const Char_t * yaxis, Int_t fontX = 12, Int_t fontY = 12, Bool_t log = false);
