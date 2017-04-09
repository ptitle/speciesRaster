#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

// return list of species per cell
// [[Rcpp::export(name = spListPerCell, rng = false)]]
List spListPerCell(NumericMatrix input) {

	int nrows = input.nrow();
	int ncols = input.ncol();
	List dimnames = input.attr( "dimnames" );
	CharacterVector colnames = dimnames[1];
	List out(nrows);

	for (int i = 0; i < nrows; i++) { 
		CharacterVector vec(ncols);
		vec = colnames[input(i,_) == 1];
		if (vec.size() > 0) {
			out[i] = vec.sort();
		} else {
			out[i] = NA_REAL;
		}
	}

	return out;
}

// Given a numeric vector and a cutoff, return the index
// of the [cutoff] largest values
// [[Rcpp::export(name = sort_, rng = false)]]
NumericVector sort_(NumericVector x, int cutoff) {

	Range idx(0, cutoff - 1);
	NumericVector sub(x.size());
	NumericVector sorted = clone(x).sort(true);
	sub = match(sorted, x);
	NumericVector res = sub[idx];
	return res;

}

// for a matrix of values and a vector of cutoff integers,
//return the index of the [cutoff] largest values per row
// [[Rcpp::export(name = returnTopIndices, rng = false)]]
List returnTopIndices(NumericMatrix input, IntegerVector cutoff) {

	int nrows = input.nrow();
	int ncols = input.ncol();

	List out(nrows);

	for (int i = 0; i < nrows; i++) {
		NumericVector x(ncols);
		x = input(i,_);

		if (cutoff[i] != 0) {

			Range idx(0, cutoff[i] - 1);
			NumericVector sub(x.size());
			NumericVector sorted = clone(x).sort(true);
			sub = match(sorted, x);
			NumericVector res = sub[idx];

			out[i] = res;
		} else {
			out[i] = NA_REAL;
		}
	}

	return out;

}

// return either mean or median value of trait per cell
// [[Rcpp::export(name = cellAvg, rng = false)]]
NumericVector cellAvg(List input, NumericVector trait, String stat) {
		
	int n = input.size();
	NumericVector out(n);
	
	for (int i = 0; i < n; i++) {
		
		CharacterVector ind = as< CharacterVector >(input[i]);
		if (ind[0] != "NA") {
			NumericVector vals = trait[ind];
		
			if (stat == "mean") {
				out[i] = double(mean(vals));
			} else if (stat == "median") {
				out[i] = double(median(vals));				
			}
		} else {
			out[i] = NA_REAL;
		}
	}
	
	return out;	
}



// given a list of character vectors, and a character vector, 
// return the list of vectors, intersected with the second vector
// [[Rcpp::export(name = intersectList, rng = false)]]
List intersectList(List input, StringVector vec) {

	int n = input.size();
	List out(n);

	for (int i = 0; i < n; i++) {
		
		StringVector sp = as< StringVector > (input[i]);

		if (all(!is_na(sp))) {

			StringVector res = intersect(sp, vec);
	
			if (res.size() > 0) {
				out[i] = intersect(sp, vec);
			} else {
				out[i] = NA_REAL;
			}
		} else {
			out[i] = NA_REAL;
		}
	}

	return out;

}

// [[Rcpp::export(name = flattenMatrix, rng = false)]]
NumericVector flattenMatrix(NumericMatrix mat) {

	std::vector<double> out = as< std::vector<double> >(mat);
	NumericVector ret = wrap(out);
	return ret;

}



// Legendre difference index
// Carvalho et al. 2013, Legendre 2014 
// [[Rcpp::export(name = LegendreDiffIndex, rng = false)]]
double LegendreDiffIndex(StringVector commA, StringVector commB) {
	if (all(is_na(commA)) | all(is_na(commB))) {
		return NA_REAL;
	} else {
		double A = union_(commA, commB).size();
		double B = setdiff(commA, commB).size();
		double C = setdiff(commB, commA).size();
		return std::abs(B - C) / (2 * A + B + C);
	}
}


// turnover component of the Jaccard’s dissimilarity index
// [[Rcpp::export(name = jaccardDissimilarity, rng = false)]]
double jaccardDissimilarity(StringVector commA, StringVector commB) {
	if (all(is_na(commA)) | all(is_na(commB))) {
		return NA_REAL;
	} else {
		double A = union_(commA, commB).size();
		double B = setdiff(commA, commB).size();
		double C = setdiff(commB, commA).size();
		return (B + C) / (A + B + C);
	}
}

// turnover component of the Jaccard’s dissimilarity index
// [[Rcpp::export(name = jaccardTurnover, rng = false)]]
double jaccardTurnover(StringVector commA, StringVector commB) {
	if (all(is_na(commA)) | all(is_na(commB))) {
		return NA_REAL;
	} else {
		double A = union_(commA, commB).size();
		double B = setdiff(commA, commB).size();
		double C = setdiff(commB, commA).size();
		return (2 * std::min(B, C)) / (A + 2 * std::min(B, C));
	}
}

// turnover component of the Jaccard’s dissimilarity index
// [[Rcpp::export(name = jaccardNestedness, rng = false)]]
double jaccardNestedness(StringVector commA, StringVector commB) {
	if (all(is_na(commA)) | all(is_na(commB))) {
		return NA_REAL;
	} else {
		double A = union_(commA, commB).size();
		double B = setdiff(commA, commB).size();
		double C = setdiff(commB, commA).size();
		return (std::abs(B - C) / (A + B + C)) * (A / (A + 2 * std::min(B, C)));
	}
}


// Function will return index of list entries that are only NA
// a list will be returned, list[[0]] will be index for empty
// list[[1]] will be index for non-empty
// [[Rcpp::export(name = ListIsEmpty, rng = false)]]
std::vector<std::vector<int> > ListIsEmpty(List input) {
	
	std::vector<std::vector<int> > out;

	int nx = input.size();
	
	//std::vector<int> y;
	//y.reserve(nx);

	std::vector<int> empty;
	std::vector<int> nonempty;

	for(int i = 0; i < nx; i++) {
		CharacterVector vec = as< CharacterVector >(input[i]);
		if (vec.size() == 1) {
			if (vec[0] == "NA") {
				//y.push_back(i);
				empty.push_back(i);
			} else {
				nonempty.push_back(i);
			}
		} else {
			nonempty.push_back(i);
		}
	}

	//IntegerVector empty = wrap(y);

	// empty is indices for elements that are only NA
	//IntegerVector v = Rcpp::seq(0, input.size() - 1);
	//IntegerVector nonEmpty = setdiff(v, empty);

	//List out(2);
	//out[0] = empty;
	//out[1] = nonEmpty;
	out.push_back(empty);
	out.push_back(nonempty);

	return out;
}

// [[Rcpp::export(name = ListIsEmptyR, rng = false)]]
List ListIsEmptyR(List input) {
	
	List out(2);

	int nx = input.size();
	
	//std::vector<int> y;
	//y.reserve(nx);

	std::vector<int> empty;
	std::vector<int> nonempty;

	for(int i = 0; i < nx; i++) {
		CharacterVector vec = as< CharacterVector >(input[i]);
		if (vec.size() == 1) {
			if (vec[0] == "NA") {
				//y.push_back(i);
				empty.push_back(i);
			} else {
				nonempty.push_back(i);
			}
		} else {
			nonempty.push_back(i);
		}
	}

	IntegerVector empty2 = wrap(empty);
	IntegerVector nonempty2 = wrap(nonempty);

	out[0] = empty2;
	out[1] = nonempty2;

	return out;
}

// http://stackoverflow.com/questions/23849354/equivalent-of-which-function-in-rcpp
// [[Rcpp::export(name = whichCpp, rng = false)]]
IntegerVector whichCpp(NumericVector x) {

	int nx = x.size();
	std::vector<int> y;
	y.reserve(nx);

	for (int i = 0; i < nx; i++) {
		if (R_IsNA(x[i])) y.push_back(i);
	}

	return wrap(y);
}




// http://stackoverflow.com/questions/30175104/how-to-effectively-combine-a-list-of-numericvectors-into-one-large-numericvector
// [[Rcpp::export(name = characterUnlist, rng = false)]]
StringVector characterUnlist(const List& list)
{
   std::size_t n = list.size();

   // Figure out the length of the output vector
   std::size_t total_length = 0;
   for (std::size_t i = 0; i < n; ++i)
      total_length += Rf_length(list[i]);

   // Allocate the vector
   StringVector output = no_init(total_length);

   // Loop and fill
   std::size_t index = 0;
   for (std::size_t i = 0; i < n; ++i)
   {
      StringVector el = list[i];
      std::copy(el.begin(), el.end(), output.begin() + index);

      // Update the index
      index += el.size();
   }

   return output;

}

// Drop NA from std::vector
// [[Rcpp::export(name = naomit, rng = false)]]
std::vector<int> naomit(std::vector<int> x){
	std::vector<int> r(x.size());
	int k=0;
	for (int i = 0; i < x.size(); ++i) {
		if (x[i] == x[i]) {
			r[k] = x[i];
			k++;
		}
	}
	r.resize(k);
	return r;
}    

// Drop negatives from std::vector
// [[Rcpp::export(name = dropNeg, rng = false)]]
std::vector<int> dropNeg(std::vector<int> x){
	std::vector<int> r(x.size());
	int k=0;
	for (int i = 0; i < x.size(); ++i) {
		if (x[i] >= 0) {
			r[k] = x[i];
			k++;
		}
	}
	r.resize(k);
	return r;
}    


// calculate multi site values for beta diversity metrics
// [[Rcpp::export(name = multiPrepCpp, rng = false)]]
std::vector<double> multiPrepCpp(List allComm) {
	
	// Si = total number of species in site i
	// S_T = total number of species in all sites
	// b_ij = number of species exclusive to site i relative to j
	// b_ji = number of species exclusive to site j relative to i

	std::vector<double> out(3);
	
	int nComm = allComm.size();
	
	// get all pairwise minimums and maximums
	std::vector<int> minMat(nComm * nComm);
	std::vector<int> maxMat(nComm * nComm);
	for (int i = 0; i < nComm; i++) {
		for (int j = 0; j < nComm; j++) {
			if (i < j) {
				std::vector<std::string> commI = as< std::vector<std::string> >(allComm[i]);
				std::vector<std::string> commJ = as< std::vector<std::string> >(allComm[j]);
				sort(commI.begin(), commI.end());
				sort(commJ.begin(), commJ.end());
				std::vector<std::string> diffIJ;
				set_difference(
					commI.begin(),
				    commI.end(),
				    commJ.begin(),
				    commJ.end(),
				    back_inserter(diffIJ));
				int ij = diffIJ.size();

				std::vector<std::string> diffJI;
				set_difference(
					commJ.begin(),
				    commJ.end(),
				    commI.begin(),
				    commI.end(),
				    back_inserter(diffJI));
				int ji = diffJI.size();

				minMat[nComm * i + j] = std::min(ij, ji);
				maxMat[nComm * i + j] = std::max(ij, ji);

			} else {
				minMat[nComm * i + j] = -1;
				maxMat[nComm * i + j] = -1;
			}
		}
	}
		
	// unique species
	CharacterVector allSp = characterUnlist(allComm);
	CharacterVector uniqueSp = unique(allSp);
	int sumSi = allSp.size();
	int S_T = uniqueSp.size();
	int Sdiff = sumSi - S_T;

	std::vector<int> minmat2 = dropNeg(minMat);
	std::vector<int> maxmat2 = dropNeg(maxMat);
	int sumMin = std::accumulate(minmat2.begin(), minmat2.end(), 0);
	int sumMax = std::accumulate(maxmat2.begin(), maxmat2.end(), 0);
	
	out[0] = double(sumMin);
	out[1] = double(sumMax);
	out[2] = double(Sdiff);
	
	return out;
}


// calculate multi site values for beta diversity metrics
// [[Rcpp::export(name = multiPrep, rng = false)]]
std::vector<double> multiPrep(List allComm) {
	
	// Si = total number of species in site i
	// S_T = total number of species in all sites
	// b_ij = number of species exclusive to site i relative to j
	// b_ji = number of species exclusive to site j relative to i

	std::vector<double> out(3);
	
	int nComm = allComm.size();
	
	// get all pairwise minimums
	NumericMatrix minMat(nComm, nComm);
	for (int i = 0; i < nComm; i++) {
		for (int j = 0; j < nComm; j++) {
			if (i < j) {
				StringVector commI = as< StringVector >(allComm[i]);
				StringVector commJ = as< StringVector >(allComm[j]);
				int ij = setdiff(commI, commJ).size();
				int ji = setdiff(commJ, commI).size();
				minMat(i,j) = std::min(ij, ji);
			} else {
				minMat(i,j) = NA_REAL;
			}
		}
	}
	
	// get all pairwise maximums
	NumericMatrix maxMat(nComm, nComm);
	for (int i = 0; i < nComm; i++) {
		for (int j = 0; j < nComm; j++) {
			if (i < j) {
				StringVector commI = as< StringVector >(allComm[i]);
				StringVector commJ = as< StringVector >(allComm[j]);
				int ij = setdiff(commI, commJ).size();
				int ji = setdiff(commJ, commI).size();
				maxMat(i,j) = std::max(ij, ji);
			} else {
				maxMat(i,j) = NA_REAL;
			}
		}
	}
	
	// unique species
	CharacterVector allSp = characterUnlist(allComm);
	CharacterVector uniqueSp = unique(allSp);
	int sumSi = allSp.size();
	int S_T = uniqueSp.size();
	int Sdiff = sumSi - S_T;

	NumericVector minMatFlat = flattenMatrix(minMat);
	NumericVector maxMatFlat = flattenMatrix(maxMat);
	int sumMin = sum(na_omit(minMatFlat));
	int sumMax = sum(na_omit(maxMatFlat));
	
	out[0] = sumMin;
	out[1] = sumMax;
	out[2] = Sdiff;
	
	return out;
}


// Jaccard dissimilarity for multiple sites
// [[Rcpp::export(name = betaJAC, rng = false)]]
double betaJAC(std::vector<double> vec) {
	return (vec[0] + vec[1]) / (vec[2] + vec[0] + vec[1]);
}

// Jaccard dissimilarity turnover component for multiple sites
// [[Rcpp::export(name = betaJTU, rng = false)]]
double betaJTU(std::vector<double> vec) {
	return (2 * vec[0]) / (vec[2] + 2 * vec[0]);
}

// Jaccard dissimilarity nestedness component for multiple sites
// [[Rcpp::export(name = betaJNE, rng = false)]]
double betaJNE(std::vector<double> vec) {
	return ((vec[1] - vec[0]) / (vec[2] + vec[0] + vec[1])) * (vec[2] / (vec[2] + 2 * vec[0]));
}

// Sorensen dissimilarity BETA SOR for multiple sites
// [[Rcpp::export(name = betaSOR, rng = false)]]
double betaSOR(std::vector<double> vec) {
	return (vec[0] + vec[1]) / (2 * vec[2] + vec[0] + vec[1]);
}

// Simpson dissimilarity; Sorensen dissimilarity for multiple sites
// [[Rcpp::export(name = betaSIM, rng = false)]]
double betaSIM(std::vector<double> vec) {
	return vec[0] / (vec[2] + vec[0]);
}

// nestedness resultant component of Sorensen dissimilarity for multiple sites
// [[Rcpp::export(name = betaSNE, rng = false)]]
double betaSNE(std::vector<double> vec) {
	return ((vec[1] - vec[0]) / (2 * vec[2] + vec[0] + vec[1])) * (vec[2] / (vec[2] + vec[0]));
}

// Calculate beta diversity distance for all cells, using mean pairwise metrics
// [[Rcpp::export(name = calcBetaPairwise, rng = false)]]
NumericVector calcBetaPairwise(List spByCell, List nbList, String metric) {
	
	int n = spByCell.size();
	NumericVector out(n);

	if (n != nbList.size()) {
		throw std::range_error("Input lists of unequal length");
	}
	
	for (int i = 0; i < n; i++) {

		// focal species comm
		StringVector commA = as< StringVector >(spByCell[i]);

	 	if (all(!is_na(commA))) {

	 		// get neighbor cells
	 		NumericVector nbCells = as< NumericVector >(nbList[i]);
		
	 		NumericVector cellDists(nbCells.size());
	 		for (int j = 0; j < nbCells.size(); j++) {
				
	 			int nb = nbCells[j];
	 			if (nb > n) {
	 				throw std::range_error("Requested cells are out of range.");
	 			}
	 			StringVector commB = as< StringVector >(spByCell[nb - 1]);

	 			if (metric == "jaccardTurnover") {
	 				cellDists[j] = jaccardTurnover(commA, commB);
	 			} else if (metric == "jaccardNestedness") {
	 				cellDists[j] = jaccardNestedness(commA, commB);
	 			} else if (metric == "jaccardDissimilarity") {
	 				cellDists[j] = jaccardDissimilarity(commA, commB);
	 			} else if (metric == "LegendreDiffIndex") {
		 			cellDists[j] = LegendreDiffIndex(commA, commB);
				} else {
					throw std::range_error("Metric not recognized.");
				}
	 		}
			
	 		out[i] = mean(na_omit(cellDists));
	 	} else {
	 		out[i] = NA_REAL;
	 	}
	}
	
	return out;
}




// Calculate beta diversity distance for all cells, using multi-site metrics
// [[Rcpp::export(name = calcBetaMultiSite, rng = false)]]
NumericVector calcBetaMultiSite(List spByCell, List nbList, String metric) {
	
	NumericVector out(spByCell.size());

	if (spByCell.size() != nbList.size()) {
		throw std::range_error("Input lists of unequal length");
	}

	Rcout << "Starting indexing...\n";
	// drop all cells that are empty
	std::vector<std::vector<int> > emptyInd = ListIsEmpty(spByCell);
	Rcout << "...Finished indexing...\n";

	IntegerVector emptyCells = wrap( emptyInd[0] );
	IntegerVector nonEmptyCells = wrap( emptyInd[1] );
	Rcout << "...Converting done...\n";

	int n = nonEmptyCells.size();

	Rcout << "Starting loop...\n";

	for (int i = 0; i < n; i++) {

		// 	subset spByCell list to relevant cells
		IntegerVector nbCells = as< IntegerVector >(nbList[nonEmptyCells[i]]);
		//std::vector<int> nbCells = as< std::vector<int> >(nbList[nonEmptyCells[i]]);
	
		// fix indexing as passed from R
		nbCells = nbCells - 1;
		//for (int j = 0; j < n; j++) {
		//	nbCells[j] = nbCells[j] - 1;
		//}
	
		// drop cells that have been identified as empty
		// IntegerVector goodCells = setdiff(nbCells, emptyCells);
		// determine which cells are empty in this set
		std::vector<std::vector<int> > allInd = ListIsEmpty(spByCell[nbCells]);
		IntegerVector ind = wrap(allInd[1]);
		IntegerVector goodCells = nbCells[ind];

	 	List subList = spByCell[goodCells];

	 	// calculate values that will be needed
	 	std::vector<double> prepVec = multiPrepCpp(subList);

		if (metric == "betaJAC") {
				out[nonEmptyCells[i]] = betaJAC(prepVec);
			} else if (metric == "betaJTU") {
				out[nonEmptyCells[i]] = betaJTU(prepVec);
			} else if (metric == "betaJNE") {
				out[nonEmptyCells[i]] = betaJNE(prepVec);
			} else if (metric == "betaSOR") {
				out[nonEmptyCells[i]] = betaSOR(prepVec);
			} else if (metric == "betaSIM") {
				out[nonEmptyCells[i]] = betaSIM(prepVec);
			} else if (metric == "betaSNE") {
				out[nonEmptyCells[i]] = betaSNE(prepVec);
			} else {
			throw std::range_error("Metric not recognized.");
		}
	 	//out[nonEmptyCells[i]] = wrap(prepVec);
	}


	out[emptyCells] = NA_REAL;
	
	return out;
}




// Calculate beta diversity distance for all cells, using multi-site metrics
// [[Rcpp::export(name = calcBetaMultiSiteBlock, rng = false)]]
List calcBetaMultiSiteBlock(List spByCell, List nbList, String metric) {
	
	List out(2);
	NumericVector cellVals(nbList.size());
	List spVals(nbList.size());

	//Rcout << "\tStarting indexing...\n";
	// drop all cells that are empty
	//std::vector<std::vector<int> > emptyInd = ListIsEmpty(spByCell);
	//Rcout << "\t...Finished indexing...\n";

	//IntegerVector emptyCells = wrap( emptyInd[0] );
	//IntegerVector nonEmptyCells = wrap( emptyInd[1] );
	//Rcout << "\t...Converting done...\n";

	int n = nbList.size();

	Rcout << "\tStarting loop...\n";

	for (int i = 0; i < n; i++) {

		// 	subset spByCell list to relevant cells
		IntegerVector nbCells = as< IntegerVector >(nbList[i]);
	
		// fix indexing as passed from R
		nbCells = nbCells - 1;
	
		// drop cells that have been identified as empty
		// determine which cells are empty in this set
		std::vector<std::vector<int> > allInd = ListIsEmpty(spByCell[nbCells]);
		IntegerVector ind = wrap(allInd[1]);
		IntegerVector goodCells = nbCells[ind];

		if (goodCells.size() > 0) {

		 	List subList = spByCell[goodCells];
		 	CharacterVector allSp = characterUnlist(subList);
			CharacterVector uniqueSp = unique(allSp);
			spVals[i] = "hello";

		 	// calculate values that will be needed
		 	std::vector<double> prepVec = multiPrepCpp(subList);

			if (metric == "betaJAC") {
					cellVals[i] = betaJAC(prepVec);
				} else if (metric == "betaJTU") {
					cellVals[i] = betaJTU(prepVec);
				} else if (metric == "betaJNE") {
					cellVals[i] = betaJNE(prepVec);
				} else if (metric == "betaSOR") {
					cellVals[i] = betaSOR(prepVec);
				} else if (metric == "betaSIM") {
					cellVals[i] = betaSIM(prepVec);
				} else if (metric == "betaSNE") {
					cellVals[i] = betaSNE(prepVec);
				} else {
				throw std::range_error("Metric not recognized.");
			}
	 	} else {
	 		cellVals[i] = NA_REAL;
	 		spVals[i] = "NA";
	 	}
	}

	out[0] = cellVals;
	out[1] = spVals;

	return out;
}
