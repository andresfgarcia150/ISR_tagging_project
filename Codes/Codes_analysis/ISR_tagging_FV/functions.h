/*
-------------------------------------------------
-------     Universidad de los Andes      -------
-------      Departamento de Física       -------
-------        Joven Investigador         -------
-------  Andrés Felipe García Albarracín  -------
-------    Juan Carlos Sanabria Arenas    -------
-------------------------------------------------

Other functions
*/

/*
 * Function that calculates the angular difference between two angles
 * @param: phi1 First angle
 * @param: phi2 Second angle
 * @return: The angular difference between phi1 and phi2
 */
Double_t deltaAng(Double_t phi1, Double_t phi2);

/*
 * Function that given a min_Value, a max_Value, the number of bins and a desired
 * value, returns the position of such value in the array
 * @param: min_Value The minimum value of the array
 * @param: max_Value The maximum value of the array
 * @param: bins		 The number of bins in the array
 * @param: value	 The value whose position in the array will be found
 * @return: The position in the array of the desired value
 * In case of underflow, it returns -1.
 * In case of overflow, it returns (bins -1)
 */
Int_t findPosition(Double_t min_Value, Double_t max_Value, Int_t bins, Double_t value);
