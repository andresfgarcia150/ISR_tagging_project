/*
 * histoN.cpp
 *
 *  Created on: Sep 7, 2014
 *      Author: felipe-dell

-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

Class of N-dimensional histograms
*/

#include "histoN.h"
#include "ROOTFunctions.h"

/*
 * Constructor of the class histoN
 * @param dims: Histogram dimensions (number of variables)
 * @param minValues_array: Array with the minimum values in each dimension
 * @param maxValues_array: Array with the maximum values in each dimension
 * @param bins_array: Array with the number of bins in each dimension
 * The attributes of the class are initialized and the memory associated to them is allocated
 */
histoN::histoN(Int_t dims, Double_t* minValues_array, Double_t* maxValues_array, Int_t* bins_array) {
	// Auto-generated constructor stub
	dimensions = dims;
	entries = 0;
	// The idea is that minValues doesn't point to the same memory direction that minValues_array, but to other direction with the same information
	minValues = (Double_t*) malloc(dimensions*sizeof(Double_t));
	maxValues = (Double_t*) malloc(dimensions*sizeof(Double_t));
	bins = (Int_t*) malloc(dimensions*sizeof(Int_t));
	steps = (Double_t*) malloc(dimensions*sizeof(Double_t));
	sizeHisto = 1;
	for(Int_t j = 0; j<dimensions; j++){
		*(minValues+j) = *(minValues_array+j); //The information pointed by minValues is the same that the info. pointed by minValues_array (The directions, however, aren't the same)
		*(maxValues+j) = (*(maxValues_array+j));
		*(bins+j) = *(bins_array+j);
		*(steps+j) = (*(maxValues+j)-*(minValues+j))/(*(bins+j));
		sizeHisto = (Int_t) sizeHisto*(*(bins+j));
	}
	histo = (Int_t*) calloc(sizeHisto,sizeof(Int_t));
}

histoN::histoN(Char_t* textFile,Char_t* binFile){
	FILE* inFile;
	inFile = fopen(textFile,"r");
	Int_t cont = 0;
	// Two lines of comments
	while(cont !=2 ){
		if(fgetc(inFile) == '\n') cont++;
	}

	fscanf(inFile,"%d %d %d\n",&dimensions,&entries,&sizeHisto);
	//cout<<dimensions<<"\t"<<entries<<"\t"<<sizeHisto<<endl;

	// Another line of comments
	cont = 0;
	while(cont !=1 ){
		if(fgetc(inFile) == '\n') cont++;
	}

	minValues = (Double_t*) malloc(dimensions*sizeof(Double_t));
	maxValues = (Double_t*) malloc(dimensions*sizeof(Double_t));
	bins = (Int_t*) malloc(dimensions*sizeof(Int_t));
	steps = (Double_t*) malloc(dimensions*sizeof(Double_t));
	for (Int_t j = 0; j<dimensions; j++){
		fscanf(inFile,"%lf %lf %d %lf\n",(minValues+j),(maxValues+j),(bins+j),(steps+j));
		//cout<<minValues[j]<<"\t"<<maxValues[j]<<"\t"<<bins[j]<<"\t"<<steps[j]<<endl;
	}

	fclose(inFile);

	ifstream ifs(binFile,std::ios::in | std::ios::binary);
	histo = (Int_t*) calloc(sizeHisto,sizeof(Int_t));
	for (Int_t j = 0; j<sizeHisto; j++){
		ifs.read((Char_t *) (histo+j),sizeof(Int_t));
	}
	ifs.close();

}

/*
 * Destructor of the class histoN
 * The memory allocated is released
 */
histoN::~histoN() {
	// Auto-generated destructor stub
	free(minValues);
	free(maxValues);
	free(bins);
	free(histo);
	free(steps);
}

/*
 * findPosition: Finds the position of a set of values within the array histo
 * @param values: Array with the variable values
 * @return Position within the array histo of the input set of variables.
 * 			-1 if any value doesn't fit within the histo
 */
Int_t histoN::findPosition(Double_t* values){
	Int_t position = 0;
	Int_t relPosition = 0;
	Int_t product = 1;
	for (Int_t j = dimensions-1; j>=0; j--){
		if ((*(values+j))<(*(minValues+j))||(*(values+j))>(*(maxValues+j))){
			position = -1;
			break;
		}
		relPosition = (Int_t)((*(values+j)-*(minValues+j))/(*(steps+j)));
		if (relPosition == *(bins+j)) relPosition = relPosition -1; 		// In case that value = maxValue
		position += relPosition*product;
		product *= *(bins+j);
	}
	return position;
}

/*
 * fill: Fills the histogram according to a set of values
 * @param values: Array with the variable values
 * @return: True if the histo is filled, False otherwise
 */
bool histoN::fill(Double_t* values){
	Int_t pos = findPosition(values);
	if (pos == -1) return false;
	else{
		histo[pos] += 1;
		entries ++;
		return true;
	}
}

/*
 * getFreqPos: Return the number of data of all the bins whose entry in dimension dim is entryDim
 * @param dim: Int_t dimension considered (dim<=0, dim < dimensions)
 * @param entryDim: Int_t Entry of the dimension considered (entryDim<=0, entryDim < *(bins+dim))
 */
Int_t histoN::getFreqPos(Int_t dim, Int_t entryDim){
	Int_t period = *(bins+dim);
	Int_t offset = entryDim;
	Int_t subEnt = 1;
	Int_t repet;
	Int_t data = 0;
	for(Int_t j=dim+1;j<dimensions;j++){
		period *= *(bins+j);
		offset *= *(bins+j);
		subEnt *= *(bins+j);
	}
	// Number of times a block of entryDim is found in the histogram
	repet = sizeHisto/period;

	// Loop over the histogram array
	Int_t pos = 0;
	Int_t posT = 0;
	for(Int_t j = 0; j<repet; j++){
		pos = j*period + offset;
		// Loop over the blocks of entryDim
		for (Int_t i = 0; i<subEnt; i++){
			posT = pos + i;
			data += *(histo+posT);
		}
	}
	return data;
}

/*
 * getProbDim: Gets an array with the probabilities in dimension dim
 * @param dim: Int_t dimension considered (dim<=0, dim < dimensions)
 */
Double_t* histoN::getProbDim(Int_t dim){
	Double_t* prob = (Double_t*) calloc(*(bins+dim),sizeof(Double_t));
	Int_t* freq = getHistDim(dim);
	for(Int_t j=0; j<*(bins+dim);j++){
		*(prob+j) = ((Double_t) *(freq+j))/entries;
	}
	return prob;
}

/*
 * getHistDim: Gets an array with the frequencies of the entries in dimension dim
 * @param dim: Int_t dimension considered (dim>=0, dim < dimensions)
 */
Int_t* histoN::getHistDim(Int_t dim){
	Int_t* freq = (Int_t*) calloc(*(bins+dim),sizeof(Int_t));
	for(Int_t j=0; j<*(bins+dim);j++){
		*(freq+j) = getFreqPos(dim,j);
	}
	return freq;
}

Double_t histoN::getProbVal(Double_t* values){
	Int_t pos = findPosition(values);
	return ((Double_t)histo[pos])/entries;
}

/*
 *	Gets the pointer to the histogram
 *	@return: Pointer to the histogram array
 */
Int_t* histoN::getHisto(){
	return histo;
}

/*
 * Gets the dimensions (number of variables) of the histogram
 * @return dimensions
 */
Int_t histoN::getDims(){
	return dimensions;
}

/*
 * Gets the total number of entries of the histogram
 * @return entries
 */
Int_t histoN::getEntries(){
	return entries;
}

/*
 * Gets the pointer to the array of minValues
 * @return Pointer to minValues
 */
Double_t* histoN::getMinValues(){
	return minValues;
}

/*
 * Gets the pointer to the array of maxValues
 * @return Pointer to maxValues
 */
Double_t* histoN::getMaxValues(){
	return maxValues;
}

/*
 * Gets the pointer to the array of number of bins
 * @return Pointer to bins
 */
Int_t* histoN::getBins(){
	return bins;
}

/*
 * Gets the pointer to the array of step sizes
 * @return Pointer to steps
 */
Double_t* histoN::getSteps(){
	return steps;
}

/*
 * Gets the size of the histogram (number of elements of the array)
 * @return sizeHisto
 */
Int_t histoN::getSizeHisto(){
	return sizeHisto;
}

/*
 * Write the information of the class in two files
 * textFile: All the information except the array histo
 * binFile: Array histo
 * @param textFile: Name of the textFile
 * @param binFile: Name of the binFile
 * @return True: If the information was written
 */
Bool_t histoN::writeClass(Char_t* textFile,Char_t* binFile){
	ofstream ofs(textFile, std::ios::out);
	if (!ofs){
		return false;
	}
	ofs<<"# Archivo con los datos de la clase histoN\n";
	ofs<<"# dimensions\tentries\tsizeHisto\n";
	ofs<<"\t\t"<<dimensions<<"\t\t"<<entries<<"\t\t"<<sizeHisto<<std::endl;
	ofs<<"#minValues\tmaxValues\tbins\tsteps"<<std::endl;
	for(Int_t j = 0;j<dimensions;j++){
		ofs<<"\t\t"<<minValues[j]<<"\t\t"<<maxValues[j]<<"\t\t"<<bins[j]<<"\t\t"<<steps[j]<<std::endl;
	}
	ofs<<"#Archivo binario"<<std::endl;
	ofs<<binFile<<std::endl;
	ofs.close();

	ofs.open(binFile,std::ios::out|std::ios::binary);
	if (!ofs){
		return false;
	}
	for(Int_t j = 0; j<sizeHisto; j++){
		ofs.write((Char_t *) (histo+j),sizeof(Int_t));
	}
	ofs.close();

	return true;
}
