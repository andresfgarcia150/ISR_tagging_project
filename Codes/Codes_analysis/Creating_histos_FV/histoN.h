/*
 * histoN.h
 *
 *  Created on: Sep 7, 2014
 *      Author: Andrés Felipe García Albarracín
 *
 *      -------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

Class of N-dimensional histograms
 */

#ifndef HISTON_H_
#define HISTON_H_
#include "ROOTFunctions.h"

class histoN{

private:
	Int_t dimensions;		// dimensions of the histogram
	Int_t entries;			// Number of entries of histo
	Double_t* minValues;	// Array with the minimum values in each dimension
	Double_t* maxValues;	// Array with the maximum values in each dimension
	Int_t* bins;			// Array with the number of bins in each dimension
	Double_t* steps;		// Array with the step size or bin size in each dimension
	Int_t* histo;			// Array whose entries are the frequencies of the bins
	Int_t sizeHisto;		// Size of histo

public:
	histoN(Int_t dims, Double_t* minValues_array, Double_t* maxValues_array, Int_t* bins_array);
	histoN(Char_t* textFile,Char_t* binFile);
	virtual ~histoN();
	Int_t findPosition(Double_t* values);
	bool fill(Double_t* values);
	Int_t getFreqPos(Int_t dim, Int_t entryDim);
	Int_t* getHistDim(Int_t dim);
	Double_t* getProbDim(Int_t dim);
	Double_t getProbVal(Double_t* values);
	Int_t* getHisto();
	Int_t getDims();
	Int_t getEntries();
	Double_t* getMinValues();
	Double_t* getMaxValues();
	Int_t* getBins();
	Double_t* getSteps();
	Int_t getSizeHisto();
	Bool_t writeClass(Char_t* textFile,Char_t* binFile);
};


#endif /* HISTON_H_ */
