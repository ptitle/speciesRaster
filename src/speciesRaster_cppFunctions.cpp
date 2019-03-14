// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <progress.hpp>
#include <progress_bar.hpp>
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

// from vector of values, create pairwise matrix and take minimum
// [[Rcpp::export(name = meanNNdist, rng = false)]]
double meanNNdist(NumericVector input) {

	int n = input.size();

	if (n == 1) {
		return 0;
	} else {
	
		// get all pairwise distances
		NumericVector minVals(n);
		for (int i = 0; i < n; i++) {
			NumericVector vec(n, NumericVector::get_na());
			for (int j = 0; j < n; j++) {
				if (i != j) {
					vec[j] = std::abs(input[i] - input[j]);
				}
			}
			minVals[i] = min(na_omit(vec));
		}

		return mean(na_omit(minVals));
	}
}

// from vector of values, create pairwise matrix and take minimum
// [[Rcpp::export(name = minNNdist, rng = false)]]
double minNNdist(NumericVector input) {

	int n = input.size();

	if (n == 1) {
		return 0;
	} else {
		
		// get all pairwise distances
		NumericVector minVals(n);
		for (int i = 0; i < n; i++) {
			NumericVector vec(n, NumericVector::get_na());
			for (int j = 0; j < n; j++) {
				if (i != j) {
					vec[j] = std::abs(input[i] - input[j]);
				}
			}
			minVals[i] = min(na_omit(vec));
		}

		return min(na_omit(minVals));
	}
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
			} else if (stat == "variance") {
				out[i] = double(var(vals));
			} else if (stat == "mean_NN_dist") {
				out[i] = double(meanNNdist(vals));
			} else if (stat == "min_NN_dist") {
				out[i] = double(minNNdist(vals));
			} else if (stat == "range") {
				out[i] = double(max(vals) - min(vals));
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

		Rcpp::checkUserInterrupt();		
		StringVector sp = as< StringVector > (input[i]);

		if (all(!is_na(sp))) {

			StringVector res = intersect(sp, vec);
	
			if (res.size() > 0) {
				out[i] = res;
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
	
	std::vector<int> empty;
	std::vector<int> nonempty;

	for(int i = 0; i < nx; i++) {
		std::vector<std::string> vec = as< std::vector<std::string> >(input[i]);
		if (vec.size() == 1) {
			if (vec[0] == "empty") {
				//y.push_back(i);
				empty.push_back(i);
			} else {
				nonempty.push_back(i);
			}
		} else {
			nonempty.push_back(i);
		}
	}

	out.push_back(empty);
	out.push_back(nonempty);

	return out;
}

// [[Rcpp::export(name = ListIsEmptyR, rng = false)]]
List ListIsEmptyR(List input) {
	
	List out(2);

	int nx = input.size();
	
	std::vector<int> empty;
	std::vector<int> nonempty;

	for(int i = 0; i < nx; i++) {
		std::vector<std::string> vec = as< std::vector<std::string> >(input[i]);
		if (vec.size() == 1) {
			if (vec[0] == "empty") {
				//y.push_back(i);
				empty.push_back(i + 1);
			} else {
				nonempty.push_back(i + 1);
			}
		} else {
			nonempty.push_back(i + 1);
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

// convert values assigned to unique cell communities, to vector of all grid cells
// [[Rcpp::export(name = uniqueCommResToFullList, rng = false)]]
NumericVector uniqueCommResToFullList(NumericVector resVal, IntegerVector cellCommInd) {

	NumericVector out(cellCommInd.size());

	for (int i = 0; i < resVal.size(); i++) {
		for (int j = 0; j < cellCommInd.size(); j++) {
			if (cellCommInd[j] == (i + 1)) {
				out[j] = resVal[i];
			}
		}
	}

	return out;
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

	Rcout << "\tIndexing cells...";
	// drop all cells that are empty
	std::vector<std::vector<int> > emptyInd = ListIsEmpty(spByCell);
	Rcout << "done\n";

	IntegerVector emptyCells = wrap( emptyInd[0] );
	IntegerVector nonEmptyCells = wrap( emptyInd[1] );
	//Rcout << "...Converting done...\n";

	int n = nonEmptyCells.size();

	Rcout << "\tStarting metric calculation...";

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

	Rcout << "done\n";


	out[emptyCells] = NA_REAL;
	
	return out;
}




// Calculate beta diversity distance for all cells, using multi-site metrics
// [[Rcpp::export(name = calcBetaMultiSiteBlock, rng = false)]]
NumericVector calcBetaMultiSiteBlock(List spByCell, List nbList, String metric) {
	
	//NumericVector out(2);
	NumericVector cellVals(nbList.size());
	//List spVals(nbList.size());

	//Rcout << "\tStarting indexing...\n";
	// drop all cells that are empty
	//std::vector<std::vector<int> > emptyInd = ListIsEmpty(spByCell);
	//Rcout << "\t...Finished indexing...\n";

	//IntegerVector emptyCells = wrap( emptyInd[0] );
	//IntegerVector nonEmptyCells = wrap( emptyInd[1] );
	//Rcout << "\t...Converting done...\n";

	int n = nbList.size();

	//Rcout << "\t...starting metric calculation...";

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
		 	//CharacterVector allSp = characterUnlist(subList);
			//CharacterVector uniqueSp = unique(allSp);
			//spVals[i] = "hello";

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
	 		//spVals[i] = "NA";
	 	}
	}
	//Rcout << "done\n";

	//out[0] = cellVals;
	//out[1] = spVals;

	return cellVals;
}


// Return index of cells that map to each unique community set
// [[Rcpp::export(name = mapComm, rng = false)]]
List mapComm(CharacterVector uniqueCommLabels, CharacterVector allComm) {

	int n = uniqueCommLabels.size();
	List out(n);
	
	std::vector<std::string> uniqueCommLabels2 = as< std::vector<std::string> >(uniqueCommLabels);
	std::vector<std::string> allComm2 = as< std::vector<std::string> >(allComm);
	for (int i = 0; i < n; i++) {

		std::vector<int> tmp;
		for (int j = 0; j < allComm.size(); j++) {
			if (allComm2[j] == uniqueCommLabels2[i]) {
				tmp.push_back(j + 1);
			}
		}

		out[i] = wrap(tmp);
	}

	return out;

}



// return index of x in integer vec
// [[Rcpp::export(name = c_which_int, rng = false)]]
int c_which_int(std::vector<int> vec, int x) {
	int nx = vec.size();
	for (int i = 0; i < nx; i++) {
		if (vec[i] == x) {
			return i;
		}
	}
	return -1;
}


// return edge indices
// [[Rcpp::export(name = getRootToTipEdges, rng = false)]]
List getRootToTipEdges(List phylo) {

	// extract components from phylo list
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge1a = edge(_, 0);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge1 = as< std::vector<int> >(edge1a);
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);
	int rootnode = tipLabels.size() + 1;

	List out = tipLabels.size();
	
	for (int i = 0; i < tipLabels.size(); i++) {
		std::vector<int> nodes;
		int childnode = i + 1;
		while (childnode != rootnode) {
			int parentnode = edge1[c_which_int(edge2, childnode)];
			nodes.push_back(childnode);
			childnode = parentnode;
		}

		std::vector<int> edgesInd(nodes.size());
		for (int j = 0; j < nodes.size(); j++) {
			edgesInd[j] = c_which_int(edge2, nodes[j]);
		}

		out[i] = edgesInd;
	}


	return out;
}


// function to determine the number of geographic cells for every branch in phylogeny
// [[Rcpp::export(name = phyloBranchRanges, rng = false)]]
List phyloBranchRanges(List phylo, List speciesList, List tipEdges) {

	// extract components from phylo list
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);
	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge1a = edge(_, 0);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge1 = as< std::vector<int> >(edge1a);
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);

	std::vector<int> allnodes = edge2;

	// tipEdges are edge indices. Convert to child node number
	// tipNodes and tipEdges are in order of tiplabels
	List tipNodes(tipEdges.size());
	for (int i = 0; i < tipEdges.size(); i++) {
		std::vector<int> spTipEdges = as< std::vector<int> >(tipEdges[i]);
		for (int j = 0; j < spTipEdges.size(); j++) {
			spTipEdges[j] = edge2[spTipEdges[j]];
		}
		tipNodes[i] = spTipEdges;
	}

	Rcpp::checkUserInterrupt();

	// For each node number, find which species have it listed
	// speciesFromNode is in order of edge2 nodes
	List speciesFromNode(allnodes.size());
	for (int i = 0; i < allnodes.size(); i++) {
		Rcpp::checkUserInterrupt();
		std::vector<int> isPresent;
		if (allnodes[i] <= tipLabels.size()) {
			isPresent.push_back(allnodes[i]);
		} else {
			for (int j = 0; j < tipNodes.size(); j++) { 
				std::vector<int> nodesForSp = as< std::vector<int> >(tipNodes[j]);
				for (int k = 0; k < nodesForSp.size(); k++) { 
					if (nodesForSp[k] == allnodes[i]) {
						isPresent.push_back(j+1);
						break;
					}
				}
			}
		}
		speciesFromNode[i] = isPresent;
	}

//	convert tip indices to names
//	allLeaves is in order of edge2 nodes
	List allLeaves(speciesFromNode.size());
	for (int i = 0; i < speciesFromNode.size(); i++) {
		Rcpp::checkUserInterrupt();
		std::vector<int> tipIndices = as< std::vector<int> >(speciesFromNode[i]);
		std::vector<std::string> tipNames(tipIndices.size());
		for (int j = 0; j < tipIndices.size(); j++) { 
			tipNames[j] = tipLabels[tipIndices[j] - 1];
		}
		allLeaves[i] = tipNames;
	}

	// for each branch (= child node), count the number of cells that contain any of the
	// tip taxa
	std::vector<int> branchCellCount(allLeaves.size());
	for (int i = 0; i < allLeaves.size(); i++) {
		Rcpp::checkUserInterrupt();
		int counter = 0;
		std::vector<std::string> nodeSp = as< std::vector<std::string> >(allLeaves[i]);
		for (int j = 0; j < speciesList.size(); j++) {
			int counter2 = 0;
			std::vector<std::string> cellSp = as< std::vector<std::string> >(speciesList[j]);
			for (int k = 0; k < nodeSp.size(); k++) {
				if (std::find(cellSp.begin(), cellSp.end(), nodeSp[k]) != cellSp.end()) {
					counter2 = counter2 + 1;
				}
			}
			if (counter2 > 0) {
				counter = counter + 1;
			}
		}
		branchCellCount[i] = counter;
	}

	List out(2);
	out[0] = edgeLengths;
	out[1] = branchCellCount;

	return out;
}


// for each species in vec, count how many cells it is found in
// [[Rcpp::export(name = countCells, rng = false)]]
NumericVector countCells(List cellList, StringVector vec) {
	
	std::vector<std::string> uniqueSp = as< std::vector<std::string> >(vec);
	std::vector<int> out(uniqueSp.size());
	
	for (int i = 0; i < cellList.size(); i++) {
		
		std::vector<std::string> cell = as< std::vector<std::string> >(cellList[i]);
		
		for (int j = 0; j < uniqueSp.size(); j++) {
		
			if (std::find(cell.begin(), cell.end(), uniqueSp[j]) != cell.end()) {
				out[j] = out[j] + 1;
			}
		}
	}
	
	return wrap(out);
}

// after creating speciesRaster in blocks, merge separate lists into 1
// [[Rcpp::export(name = mergeLists, rng = false)]]
List mergeLists(List input) {

	List firstList = as< List >(input[0]);
	
	int n = firstList.size();
	List out(n);
	
	int iterN = input.size();
		
	// first loop: blocks
	for (int i = 0; i < n; i++) {

		std::vector<std::string> cellSp;
	
		// across cells
		for (int j = 0; j < iterN; j++) {
			List subList = as< List >(input[j]);
			std::vector<std::string> cell = as< std::vector<std::string> >(subList[i]);
			for (int k = 0; k < cell.size(); k++) {	
				if (cell[k] != "empty") {
					cellSp.push_back(cell[k]);
				}
			}
		}
		out[i] = cellSp;
	}

	return out;
}


// return intersect of two species vectors
// [[Rcpp::export(name = getComponentA, rng = false)]]
std::vector<std::string> getComponentA(std::vector<std::string> commI, std::vector<std::string> commJ) {

	// intersect(cellI, cellJ)
	std::vector<std::string> a;
	// for each species in cellI, is it present in cellJ?
	for (int k = 0; k < commI.size(); k++) { 
		if (std::find(commJ.begin(), commJ.end(), commI[k]) != commJ.end()) {
			a.push_back(commI[k]);
		}
	}

	return a;
}


// return species in commJ but not commI
// [[Rcpp::export(name = getComponentB, rng = false)]]
std::vector<std::string> getComponentB(std::vector<std::string> commI, std::vector<std::string> commJ) {

	// setdiff(cellI, cellJ)
	std::vector<std::string> diffIJ;

	sort(commI.begin(), commI.end());
	sort(commJ.begin(), commJ.end());
	set_difference(
		commI.begin(),
	    commI.end(),
	    commJ.begin(),
	    commJ.end(),
	    back_inserter(diffIJ));
	
	return diffIJ;
}


// return species in commI but not commJ
// [[Rcpp::export(name = getComponentC, rng = false)]]
std::vector<std::string> getComponentC(std::vector<std::string> commI, std::vector<std::string> commJ) {

	// setdiff(cellJ, cellI)
	std::vector<std::string> diffJI;

	sort(commI.begin(), commI.end());
	sort(commJ.begin(), commJ.end());
	set_difference(
		commJ.begin(),
	    commJ.end(),
	    commI.begin(),
	    commI.end(),
	    back_inserter(diffJI));
	
	return diffJI;
}


// return index of x in character vec
// [[Rcpp::export(name = c_which_char, rng = false)]]
int c_which_char(std::vector<std::string> vec, std::string x) {
	int nx = vec.size();
	for (int i = 0; i < nx; i++) {
		if (vec[i] == x) {
			return i;
		}
	}
	return -1;
}



// return geog area of phylo branches, given species
// [[Rcpp::export(name = weightedPhylo, rng = false)]]
double weightedPhylo(std::vector<std::string> a, std::vector<std::string> tipLabels, List spEdges, std::vector<double> edgeArea1, std::vector<double> edgeArea2) {
		
	std::vector<int> branchIndices;
	
	// get index of species labels and get branch indices from spEdges
	for (int i = 0; i < a.size(); i++) {
		
		int tmpInd = c_which_char(tipLabels, a[i]);
		std::vector<int> branchInd = as< std::vector<int> >(spEdges[tmpInd]);
		
		for (int j = 0; j < branchInd.size(); j++) {
			if (std::find(branchIndices.begin(), branchIndices.end(), branchInd[j]) == branchIndices.end()) {
				branchIndices.push_back(branchInd[j]);
			}
		}
	}

	// use branchIndices to look up branch and branch area in edgeArea table
	double summedRes = 0;

	for (int i = 0; i < branchIndices.size(); i++) {
		summedRes = summedRes + edgeArea1[branchIndices[i]] / edgeArea2[branchIndices[i]];
	}

	return summedRes;
	
}




// Calculate taxonomic beta diversity
// [[Rcpp::export(name = calcRWTurnover_taxonomic_old, rng = false)]]
NumericVector calcRWTurnover_taxonomic_old(List spByCell, List nbList) {
	
	NumericVector out(spByCell.size());

	if (spByCell.size() != nbList.size()) {
		throw std::range_error("Input lists of unequal length");
	}

	// for each cell, identify neighborhood cells
	for (int i = 0; i < spByCell.size(); i++) {

		std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[i]);

 		// 	pull out the neighborhood cells for the i'th cell
 		std::vector<int> cellNeighbors = as< std::vector<int> >(nbList[i]);
 		NumericVector cellVec(cellNeighbors.size());
		
 		for (int j = 0; j < cellNeighbors.size(); j++) {

			// for the focal cell and each neighbor, calculate:
			// a = intersect(cellI, cellJ)
			// b = setdiff(cellI, cellJ)
			// c = setdiff(cellJ, cellI)

			std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[cellNeighbors[j] - 1]);
			std::vector<std::string> a = getComponentA(commI, commJ);

			if (commI.size() == commJ.size() && commI.size() == a.size()) {
				cellVec[j] = 0.0;
			} else if (a.size() == 0) {
				cellVec[j] = 1.0;
			} else {

				std::vector<std::string> b = getComponentB(commI, commJ);
				std::vector<std::string> c = getComponentC(commI, commJ);

				// Simpson's beta diversity index
				// cellVec[j] = 1.0 - double(a.size()) / (double(a.size()) + std::min(double(b.size()), double(c.size())));
				
				// Sorenson metric
				cellVec[j] = 1.0 - 2 * double(a.size()) / (2 * double(a.size()) + double(b.size()) + double(c.size()));
			
				// Jaccard metric
				//cellVec[j] = 1.0 - double(a.size()) / (double(a.size()) + double(b.size() + double(c.size())));
			}
		}
		out[i] = mean(cellVec);
	}

	return out;
}


// Calculate range-weighted beta diversity
// [[Rcpp::export(name = calcRWTurnover_rangeWeighted_old, rng = false)]]
NumericVector calcRWTurnover_rangeWeighted_old(List spByCell, List nbList, NumericVector cellCountsR) {
	
	NumericVector out(spByCell.size());

	if (spByCell.size() != nbList.size()) {
		throw std::range_error("Input lists of unequal length");
	}

	// pull out cell count information
	std::vector<std::string> cellCountNames = as< std::vector<std::string> >(cellCountsR.attr("names"));
	std::vector<double> cellCounts = as< std::vector<double> >(cellCountsR);

	// for each cell, identify neighborhood cells
	for (int i = 0; i < spByCell.size(); i++) {

		std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[i]);

 		// 	pull out the neighborhood cells for the i'th cell
 		std::vector<int> cellNeighbors = as< std::vector<int> >(nbList[i]);
 		NumericVector cellVec(cellNeighbors.size());
		
 		for (int j = 0; j < cellNeighbors.size(); j++) {

			// for the focal cell and each neighbor, calculate:
			// a = intersect(cellI, cellJ)
			// b = setdiff(cellI, cellJ)
			// c = setdiff(cellJ, cellI)

			std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[cellNeighbors[j] - 1]);
			std::vector<std::string> a = getComponentA(commI, commJ);

			if (commI.size() == commJ.size() && commI.size() == a.size()) {
				cellVec[j] = 0.0;
			} else if (a.size() == 0) {
				cellVec[j] = 1.0;
			} else {

				std::vector<std::string> b = getComponentB(commI, commJ);
				std::vector<std::string> c = getComponentC(commI, commJ);

				double rwA = 0.0;
				for (int k = 0; k < a.size(); k++) {
					int ind = c_which_char(cellCountNames, a[k]);
					rwA = rwA + cellCounts[ind];
				}

				double rwB = 0.0;
				for (int k = 0; k < b.size(); k++) {
					int ind = c_which_char(cellCountNames, b[k]);
					rwB = rwB + cellCounts[ind];
				}

				double rwC = 0.0;
				for (int k = 0; k < c.size(); k++) {
					int ind = c_which_char(cellCountNames, c[k]);
					rwC = rwC + cellCounts[ind];
				}

				cellVec[j] = 1.0 - rwA / (rwA + rwB + rwC);
			}
		}
		out[i] = mean(cellVec);
	}

	return out;

}





// Calculate range-weighted phylogenetic turnover
// [[Rcpp::export(name = calcRWTurnover_phyloRangeWeighted_old, rng = false)]]
NumericVector calcRWTurnover_phyloRangeWeighted_old(List spByCell, List nbList, List phylo, List spEdges, NumericMatrix edgeArea) {
	
	NumericVector out(spByCell.size());

	if (spByCell.size() != nbList.size()) {
		throw std::range_error("Input lists of unequal length");
	}

	// extract relevant info from input data
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	
	NumericVector edgeArea1a = edgeArea(_, 0);
	NumericVector edgeArea2a = edgeArea(_, 1);
  
	std::vector<double> edgeArea1 = as< std::vector<double> >(edgeArea1a);
	std::vector<double> edgeArea2 = as< std::vector<double> >(edgeArea2a);

	// for each cell, identify neighborhood cells
	for (int i = 0; i < spByCell.size(); i++) {

		std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[i]);

 		// 	pull out the neighborhood cells for the i'th cell
 		std::vector<int> cellNeighbors = as< std::vector<int> >(nbList[i]);
 		NumericVector cellVec(cellNeighbors.size());
		
 		for (int j = 0; j < cellNeighbors.size(); j++) {

			// for the focal cell and each neighbor, calculate:
			// a = intersect(cellI, cellJ)
			// b = setdiff(cellI, cellJ)
			// c = setdiff(cellJ, cellI)

			std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[cellNeighbors[j] - 1]);
			std::vector<std::string> a = getComponentA(commI, commJ);

			if (commI.size() == commJ.size() && commI.size() == a.size()) {
				cellVec[j] = 0.0;
			} else if (a.size() == 0) {
				cellVec[j] = 1.0;
			} else {

				std::vector<std::string> b = getComponentB(commI, commJ);
				std::vector<std::string> c = getComponentC(commI, commJ);

				double wpA = weightedPhylo(a, tipLabels, spEdges, edgeArea1, edgeArea2);
				double wpB = weightedPhylo(b, tipLabels, spEdges, edgeArea1, edgeArea2);
				double wpC = weightedPhylo(c, tipLabels, spEdges, edgeArea1, edgeArea2);


				cellVec[j] = 1.0 - wpA / (wpA + wpB + wpC);
			}
		}
		out[i] = mean(cellVec);
	}

	return out;

}



// Calculate beta diversity distance for one cell, averaged across its neighbors
// [[Rcpp::export(name = calcRWTurnover_taxonomic_singleCell, rng = false)]]
NumericVector calcRWTurnover_taxonomic_singleCell(StringVector focalCell, List nbList) {

	std::vector<std::string> commI = as< std::vector<std::string> >(focalCell);

	std::vector<double> resVec(nbList.size());
			
 	for (int i = 0; i < nbList.size(); i++) {

		// for the focal cell and each neighbor, calculate:
		// a = intersect(cellI, cellJ)
		// b = setdiff(cellI, cellJ)
		// c = setdiff(cellJ, cellI)

		std::vector<std::string> commJ = as< std::vector<std::string> >(nbList[i]);
		std::vector<std::string> a = getComponentA(commI, commJ);

		if (commI.size() == commJ.size() && commI.size() == a.size()) {
			resVec[i] = 0.0;
		} else if (a.size() == 0) {
			resVec[i] = 1.0;
		} else {

			std::vector<std::string> b = getComponentB(commI, commJ);
			std::vector<std::string> c = getComponentC(commI, commJ);

			// Simpson's beta diversity index
			// resVec[i] = 1.0 - double(a.size()) / (double(a.size()) + std::min(double(b.size()), double(c.size())));
			
			// Sorenson metric
			resVec[i] = 1.0 - 2.0 * double(a.size()) / (2 * double(a.size()) + double(b.size()) + double(c.size()));
			
			// Jaccard metric
			//resVec[i] = 1.0 - double(a.size()) / (double(a.size()) + double(b.size()) + double(c.size()));
		}
	}
		
	double res = std::accumulate(resVec.begin(), resVec.end(), 0.0);
	double out = res / resVec.size();

	return (wrap(out));
}




// Calculate beta diversity distance for one cell, averaged across its neighbors
// [[Rcpp::export(name = calcRWTurnover_rangeWeighted_singleCell, rng = false)]]
NumericVector calcRWTurnover_rangeWeighted_singleCell(StringVector focalCell, List nbList, NumericVector cellCountsR) {

	std::vector<std::string> commI = as< std::vector<std::string> >(focalCell);
	
	// pull out cell count information
	std::vector<std::string> cellCountNames = as< std::vector<std::string> >(cellCountsR.attr("names"));
	std::vector<double> cellCounts = as< std::vector<double> >(cellCountsR);

	std::vector<double> resVec(nbList.size());
			
 	for (int i = 0; i < nbList.size(); i++) {

		// for the focal cell and each neighbor, calculate:
		// a = intersect(cellI, cellJ)
		// b = setdiff(cellI, cellJ)
		// c = setdiff(cellJ, cellI)

		std::vector<std::string> commJ = as< std::vector<std::string> >(nbList[i]);
		std::vector<std::string> a = getComponentA(commI, commJ);

		if (commI.size() == commJ.size() && commI.size() == a.size()) {
			resVec[i] = 0.0;
		} else if (a.size() == 0) {
			resVec[i] = 1.0;
		} else {

			std::vector<std::string> b = getComponentB(commI, commJ);
			std::vector<std::string> c = getComponentC(commI, commJ);

			double rwA = 0.0;
			for (int k = 0; k < a.size(); k++) {
				int ind = c_which_char(cellCountNames, a[k]);
				rwA = rwA + cellCounts[ind];
			}

			double rwB = 0.0;
			for (int k = 0; k < b.size(); k++) {
				int ind = c_which_char(cellCountNames, b[k]);
				rwB = rwB + cellCounts[ind];
			}

			double rwC = 0.0;
			for (int k = 0; k < c.size(); k++) {
				int ind = c_which_char(cellCountNames, c[k]);
				rwC = rwC + cellCounts[ind];
			}

			resVec[i] = 1.0 - rwA / (rwA + rwB + rwC);
		}
	}

		
	double res = std::accumulate(resVec.begin(), resVec.end(), 0.0);
	double out = res / resVec.size();

	return (wrap(out));
}


// Calculate beta diversity distance for one cell, averaged across its neighbors
// [[Rcpp::export(name = calcRWTurnover_phyloRangeWeighted_singleCell, rng = false)]]
NumericVector calcRWTurnover_phyloRangeWeighted_singleCell(StringVector focalCell, List nbList, StringVector phyloTipLabels, List spEdges, NumericMatrix edgeArea) {

	// extract relevant info from input data
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phyloTipLabels);
	
	NumericVector edgeArea1a = edgeArea(_, 0);
	NumericVector edgeArea2a = edgeArea(_, 1);
  
	std::vector<double> edgeArea1 = as< std::vector<double> >(edgeArea1a);
	std::vector<double> edgeArea2 = as< std::vector<double> >(edgeArea2a);

	std::vector<std::string> commI = as< std::vector<std::string> >(focalCell);

	std::vector<double> resVec(nbList.size());
			
 	for (int i = 0; i < nbList.size(); i++) {

		// for the focal cell and each neighbor, calculate:
		// a = intersect(cellI, cellJ)
		// b = setdiff(cellI, cellJ)
		// c = setdiff(cellJ, cellI)

		std::vector<std::string> commJ = as< std::vector<std::string> >(nbList[i]);
		std::vector<std::string> a = getComponentA(commI, commJ);

		if (commI.size() == commJ.size() && commI.size() == a.size()) {
			resVec[i] = 0.0;
		} else if (a.size() == 0) {
			resVec[i] = 1.0;
		} else {

			std::vector<std::string> b = getComponentB(commI, commJ);
			std::vector<std::string> c = getComponentC(commI, commJ);

			double wpA = weightedPhylo(a, tipLabels, spEdges, edgeArea1, edgeArea2);
			double wpB = weightedPhylo(b, tipLabels, spEdges, edgeArea1, edgeArea2);
			double wpC = weightedPhylo(c, tipLabels, spEdges, edgeArea1, edgeArea2);

			resVec[i] = 1.0 - wpA / (wpA + wpB + wpC);
		}	

	}
		
	double res = std::accumulate(resVec.begin(), resVec.end(), 0.0);
	double out = res / resVec.size();

	return (wrap(out));
}




// function to return cell number, given row and column
// [[Rcpp::export(name = getCellFromRowCol, rng = false)]]
int getCellFromRowCol(int rowInd, int colInd, int nCol) {
	return (rowInd) * nCol + colInd;
}



// return cells within window, where raster values are returned if not negative
// this way, bad cells can be coded as -1
// [[Rcpp::export(name = getMovingWindowCells, rng = false)]]
std::vector<int> getMovingWindowCells(int nRow, int nCol, int focalCell, int radius, std::vector<int>  rasterValues) {

	// convert indexing from R
	int cell = focalCell - 1;

	// get row and column from cell number
	int cellRow = std::floor(cell / nCol - 1) + 1;
	int cellCol = cell - nCol * cellRow;

//	Rcout << "cellrow: " << cellRow << std::endl;
//	Rcout << "cellcol: " << cellCol << std::endl;

	// get neighborhood row and column bounds
	int leftBoundCol = std::max(0, cellCol - radius);
	int rightBoundCol = std::min(nCol - 1, cellCol + radius);

	int topBoundRow = std::max(0, cellRow - radius);
	int bottomBoundRow = std::min(nRow - 1, cellRow + radius);

	std::vector<int> allRows((bottomBoundRow - topBoundRow) + 1);
	std::iota(allRows.begin(), allRows.end(), topBoundRow);

	std::vector<int> nbCells;
	for (int i = 0; i < allRows.size(); i++) {
		// for each row, pull out cells in window
		int seqStart = getCellFromRowCol(allRows[i], leftBoundCol, nCol);
		int seqEnd = getCellFromRowCol(allRows[i], rightBoundCol, nCol);
		std::vector<int> tmpCells((seqEnd - seqStart) + 1);
		std::iota(tmpCells.begin(), tmpCells.end(), seqStart);

		// identify which cells are valid, and store
		for (int j = 0; j < tmpCells.size(); j++) {
			if (rasterValues[tmpCells[j]] > 0 && tmpCells[j] != cell) {
				nbCells.push_back(rasterValues[tmpCells[j]]);
			}
		}
	}

	return nbCells;
}


// Calculate taxonomic beta diversity
// [[Rcpp::export(name = calcRWTurnover_taxonomic, rng = false)]]
NumericVector calcRWTurnover_taxonomic(List spByCell, int radius, int rasterNRow, int rasterNCol, NumericVector rasterValuesR, NumericVector nonNAcellsR, bool showProgress) {
	
	NumericVector out(spByCell.size());

	std::vector<int> rasterValues = as< std::vector<int> >(rasterValuesR);
	std::vector<int> nonNAcells = as< std::vector<int> >(nonNAcellsR);

	int n = spByCell.size();
	Progress p(n, showProgress);

	// for each cell, identify neighborhood cells
	for (int i = 0; i < n; i++) {

		Rcpp::checkUserInterrupt();
		p.increment(); 
		// if (fmod(i / double(n), 0.05) < (1 / double(n))) {
		// 	Rcpp::checkUserInterrupt();
		// 	//p.increment(); 
		// 	//Rcout << i / double(n) << "% completed..." << std::endl;
		// }

		std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[i]);

 		// 	pull out the neighborhood cells for the i'th cell
 		std::vector<int> cellNeighbors = getMovingWindowCells(rasterNRow, rasterNCol, nonNAcells[i], radius, rasterValues);

 		NumericVector cellVec(cellNeighbors.size());
		
 		for (int j = 0; j < cellNeighbors.size(); j++) {

			// for the focal cell and each neighbor, calculate:
			// a = intersect(cellI, cellJ)
			// b = setdiff(cellI, cellJ)
			// c = setdiff(cellJ, cellI)

			std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[cellNeighbors[j]]);
			std::vector<std::string> a = getComponentA(commI, commJ);

			if (commI.size() == commJ.size() && commI.size() == a.size()) {
				cellVec[j] = 0.0;
			} else if (a.size() == 0) {
				cellVec[j] = 1.0;
			} else {

				std::vector<std::string> b = getComponentB(commI, commJ);
				std::vector<std::string> c = getComponentC(commI, commJ);

				// Simpson's beta diversity index
				// cellVec[j] = 1.0 - double(a.size()) / (double(a.size()) + std::min(double(b.size()), double(c.size())));
				
				// Sorenson metric
				cellVec[j] = 1.0 - 2 * double(a.size()) / (2 * double(a.size()) + double(b.size()) + double(c.size()));
			
				// Jaccard metric
				//cellVec[j] = 1.0 - double(a.size()) / (double(a.size()) + double(b.size() + double(c.size())));
			}
		}
		out[i] = mean(cellVec);
	}

	return out;
}



// Calculate range-weighted beta diversity
// [[Rcpp::export(name = calcRWTurnover_rangeWeighted, rng = false)]]
NumericVector calcRWTurnover_rangeWeighted(List spByCell, int radius, int rasterNRow, int rasterNCol, NumericVector rasterValuesR, NumericVector nonNAcellsR, NumericVector cellCountsR, bool showProgress) {
	
	NumericVector out(spByCell.size());

	std::vector<int> rasterValues = as< std::vector<int> >(rasterValuesR);
	std::vector<int> nonNAcells = as< std::vector<int> >(nonNAcellsR);

	// pull out cell count information
	std::vector<std::string> cellCountNames = as< std::vector<std::string> >(cellCountsR.attr("names"));
	std::vector<double> cellCounts = as< std::vector<double> >(cellCountsR);

	int n = spByCell.size();
	Progress p(n, showProgress);

	// for each cell, identify neighborhood cells
	for (int i = 0; i < n; i++) {

		Rcpp::checkUserInterrupt();
		p.increment(); 
		// if (fmod(i / double(n), 0.05) < (1 / double(n))) {
		// 	Rcpp::checkUserInterrupt();
		// 	//p.increment(); 
		// 	//Rcout << i / double(n) << "% completed..." << std::endl;
		// }

		std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[i]);

 		// 	pull out the neighborhood cells for the i'th cell
 		std::vector<int> cellNeighbors = getMovingWindowCells(rasterNRow, rasterNCol, nonNAcells[i], radius, rasterValues);

 		NumericVector cellVec(cellNeighbors.size());
		
 		for (int j = 0; j < cellNeighbors.size(); j++) {

			// for the focal cell and each neighbor, calculate:
			// a = intersect(cellI, cellJ)
			// b = setdiff(cellI, cellJ)
			// c = setdiff(cellJ, cellI)

			std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[cellNeighbors[j]]);
			std::vector<std::string> a = getComponentA(commI, commJ);

			if (commI.size() == commJ.size() && commI.size() == a.size()) {
				cellVec[j] = 0.0;
			} else if (a.size() == 0) {
				cellVec[j] = 1.0;
			} else {

				std::vector<std::string> b = getComponentB(commI, commJ);
				std::vector<std::string> c = getComponentC(commI, commJ);

				double rwA = 0.0;
				for (int k = 0; k < a.size(); k++) {
					int ind = c_which_char(cellCountNames, a[k]);
					rwA = rwA + cellCounts[ind];
				}

				double rwB = 0.0;
				for (int k = 0; k < b.size(); k++) {
					int ind = c_which_char(cellCountNames, b[k]);
					rwB = rwB + cellCounts[ind];
				}

				double rwC = 0.0;
				for (int k = 0; k < c.size(); k++) {
					int ind = c_which_char(cellCountNames, c[k]);
					rwC = rwC + cellCounts[ind];
				}

				cellVec[j] = 1.0 - rwA / (rwA + rwB + rwC);
			}
		}
		out[i] = mean(cellVec);
	}

	return out;
}



// Calculate range-weighted phylogenetic turnover
// [[Rcpp::export(name = calcRWTurnover_phyloRangeWeighted, rng = false)]]
NumericVector calcRWTurnover_phyloRangeWeighted(List spByCell, int radius, int rasterNRow, int rasterNCol, NumericVector rasterValuesR, NumericVector nonNAcellsR, List phylo, List spEdges, NumericMatrix edgeArea, bool showProgress) {
	
	NumericVector out(spByCell.size());

	std::vector<int> rasterValues = as< std::vector<int> >(rasterValuesR);
	std::vector<int> nonNAcells = as< std::vector<int> >(nonNAcellsR);

	// extract relevant info from input data
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	
	NumericVector edgeArea1a = edgeArea(_, 0);
	NumericVector edgeArea2a = edgeArea(_, 1);
  
	std::vector<double> edgeArea1 = as< std::vector<double> >(edgeArea1a);
	std::vector<double> edgeArea2 = as< std::vector<double> >(edgeArea2a);

	int n = spByCell.size();
	Progress p(n, showProgress);

	// for each cell, identify neighborhood cells
	for (int i = 0; i < n; i++) {

		Rcpp::checkUserInterrupt();
		p.increment(); 
		// if (fmod(i / double(n), 0.05) < (1 / double(n))) {
		// 	Rcpp::checkUserInterrupt();
		// 	//p.increment(); 
		// 	//Rcout << i / double(n) << "% completed..." << std::endl;
		// }

		std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[i]);

 		// 	pull out the neighborhood cells for the i'th cell
 		std::vector<int> cellNeighbors = getMovingWindowCells(rasterNRow, rasterNCol, nonNAcells[i], radius, rasterValues);

 		NumericVector cellVec(cellNeighbors.size());
		
 		for (int j = 0; j < cellNeighbors.size(); j++) {

			// for the focal cell and each neighbor, calculate:
			// a = intersect(cellI, cellJ)
			// b = setdiff(cellI, cellJ)
			// c = setdiff(cellJ, cellI)

			std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[cellNeighbors[j]]);
			std::vector<std::string> a = getComponentA(commI, commJ);

			if (commI.size() == commJ.size() && commI.size() == a.size()) {
				cellVec[j] = 0.0;
			} else if (a.size() == 0) {
				cellVec[j] = 1.0;
			} else {

				std::vector<std::string> b = getComponentB(commI, commJ);
				std::vector<std::string> c = getComponentC(commI, commJ);

				double wpA = weightedPhylo(a, tipLabels, spEdges, edgeArea1, edgeArea2);
				double wpB = weightedPhylo(b, tipLabels, spEdges, edgeArea1, edgeArea2);
				double wpC = weightedPhylo(c, tipLabels, spEdges, edgeArea1, edgeArea2);


				cellVec[j] = 1.0 - wpA / (wpA + wpB + wpC);
			}
		}
		out[i] = mean(cellVec);
	}

	return out;

}










// union two std string vectors
// [[Rcpp::export(name = getUnion, rng = false)]]
std::vector<std::string> getUnion(std::vector<std::string> vec1, std::vector<std::string> vec2) {

	std::vector<std::string> out;

	for (int i = 0; i < vec1.size(); i++) {
		if (std::find(out.begin(), out.end(), vec1[i]) == out.end()) {
			out.push_back(vec1[i]);
		}
	}

	for (int i = 0; i < vec2.size(); i++) {
		if (std::find(out.begin(), out.end(), vec2[i]) == out.end()) {
			out.push_back(vec2[i]);
		}
	}

	return out;
}


// return a list of all nodes, containing tiplabels that are descendant from each
// [[Rcpp::export(name = getLeavesForNodes, rng = false)]]
List getLeavesForNodes(List phylo) {

	// extract components from phylo list
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge1a = edge(_, 0);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge1 = as< std::vector<int> >(edge1a);
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);
	int rootnode = tipLabels.size() + 1;
	int nodemax = *std::max_element(std::begin(edge2), std::end(edge2));
	std::vector<int> allnodes(nodemax);
	std::iota(allnodes.begin(), allnodes.end(), 1);

	// create list of tips, containing root-to-tip nodes
	List rootToTipNodes(tipLabels.size());
	for (int i = 0; i < tipLabels.size(); i++) {
		std::vector<int> nodes;
		int childnode = i + 1;
		while (childnode != rootnode) {
			int parentnode = edge1[c_which_int(edge2, childnode)];
			nodes.push_back(childnode);
			childnode = parentnode;
		}
		nodes.push_back(rootnode);
		rootToTipNodes[i] = nodes;
	}

	// create list of all nodes, and for each node, list the descendant taxa
	List nodeLeaves(allnodes.size());
	for (int i = 0; i < allnodes.size(); i++) {
		std::vector<std::string> foundLeaves;
		for (int j = 0; j < rootToTipNodes.size(); j++) {
			std::vector<int> tmp = as< std::vector<int> >(rootToTipNodes[j]);
			if (std::find(tmp.begin(), tmp.end(), allnodes[i]) != tmp.end()) {
				foundLeaves.push_back(tipLabels[j]);
			}
		}
		nodeLeaves[i] = foundLeaves;
	}


	return nodeLeaves;
}


// function to receive the nodeLeaves list and a set of tiplabels, to return the MRCA
// [[Rcpp::export(name = getMRCA_from_nodeLeaves, rng = false)]]
int getMRCA_from_nodeLeaves(List nodeLeaves, std::vector<std::string> taxa) {

//	std::vector<std::string> taxa = as< std::vector<std::string> >(tip);

	// for each node, does it contain all target taxa downstream?
	std::vector<int> commonNodes;
	for (int i = 0; i < nodeLeaves.size(); i++) {
		std::vector<std::string> node = as< std::vector<std::string> >(nodeLeaves[i]);
		int counter = 0;

		for (int j = 0; j < taxa.size(); j++) {
			if (std::find(node.begin(), node.end(), taxa[j]) != node.end()) {
				counter = counter + 1;
			}
		}

		// if node contained all taxa downstream, then counter should be equal to number of taxa
		if (counter == taxa.size()) {
			commonNodes.push_back(i + 1);
		}
	}

	// most recent common ancestor will be max node number
	// if just one tip, then it will be the minimum

	int mrca;

	if (taxa.size() > 1) {
	 	mrca = *std::max_element(std::begin(commonNodes), std::end(commonNodes));
	} else {
		mrca = *std::min_element(std::begin(commonNodes), std::end(commonNodes));
	}

	return mrca;
}

// // function to return the terminal branch length, given a species name
// // [[Rcpp::export(name = returnTerminalBranch, rng = false)]]
// std::vector<int> returnTerminalBranch(std::vector<std::string> tip, List phylo) {

// 	//std::vector<std::string> tip = as< std::vector<std::string> >(a);
	
// 	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
// 	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);

// 	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
// //	NumericVector edge1a = edge(_, 0);
// 	NumericVector edge2a = edge(_, 1);
  
// //	std::vector<int> edge1 = as< std::vector<int> >(edge1a);
// 	std::vector<int> edge2 = as< std::vector<int> >(edge2a);

// 	std::vector<int> branchIndices;
// 	int spInd = c_which_char(tipLabels, tip[0]);
// 	branchIndices.push_back(c_which_int(edge2, spInd + 1));

// 	return branchIndices;
// }


// return geog area of phylo branches, given species
// [[Rcpp::export(name = FaithPD_branchIndices, rng = false)]]
std::vector<int> FaithPD_branchIndices(std::vector<std::string> a, List phylo, List nodeLeaves, List spEdges) {
	
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);

	NumericMatrix edge = as<NumericMatrix>(phylo["edge"]);
	NumericVector edge2a = edge(_, 1);
  
	std::vector<int> edge2 = as< std::vector<int> >(edge2a);

	std::vector<int> branchIndices;

	// get mrca node of taxon set a
	if (a.size() > 1) {
		int mrca = getMRCA_from_nodeLeaves(nodeLeaves, a);

		// get mrca index in edge2 (= edge index)
		int mrcaInd = c_which_int(edge2, mrca);
		
		// get index of species labels and get branch indices from spEdges
		// only keep branch if branch index is greater than MRCA index
		for (int i = 0; i < a.size(); i++) {
			
			int tmpInd = c_which_char(tipLabels, a[i]);
			std::vector<int> branchInd = as< std::vector<int> >(spEdges[tmpInd]);
			
			for (int j = 0; j < branchInd.size(); j++) {
				if (std::find(branchIndices.begin(), branchIndices.end(), branchInd[j]) == branchIndices.end()) {
					if (branchInd[j] > mrcaInd) {
						branchIndices.push_back(branchInd[j]);
					}
				}
			}
		}

	} else {
		// if 1 species, return terminal branch length
		int spInd = c_which_char(tipLabels, a[0]);
		branchIndices.push_back(c_which_int(edge2, spInd + 1));
	}

	return branchIndices;
}






// Calculate phylosor phylogenetic turnover
// [[Rcpp::export(name = calcPhylosor, rng = false)]]
NumericVector calcPhylosor(List spByCell, int radius, int rasterNRow, int rasterNCol, NumericVector rasterValuesR, NumericVector nonNAcellsR, List phylo, bool showProgress) {
	
	NumericVector out(spByCell.size());

	std::vector<int> rasterValues = as< std::vector<int> >(rasterValuesR);
	std::vector<int> nonNAcells = as< std::vector<int> >(nonNAcellsR);

	// extract relevant info from input data
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);
	
	// get leaves for all nodes as well as tip edges
	List nodeLeaves = getLeavesForNodes(phylo);
	List spEdges = getRootToTipEdges(phylo);

	int n = spByCell.size();
	Progress p(n, showProgress);

	// for each cell, identify neighborhood cells
	for (int i = 0; i < n; i++) {

		Rcpp::checkUserInterrupt();
		p.increment(); 

		std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[i]);

		// get Faith's PD for commI 
		std::vector<int> edgesCommI = FaithPD_branchIndices(commI, phylo, nodeLeaves, spEdges);
		double pdCommI = 0;
		for (int k = 0; k < edgesCommI.size(); k++) {
			pdCommI = pdCommI + edgeLengths[edgesCommI[k]];
		}

 		// 	pull out the neighborhood cells for the i'th cell
 		std::vector<int> cellNeighbors = getMovingWindowCells(rasterNRow, rasterNCol, nonNAcells[i], radius, rasterValues);

 		NumericVector cellVec(cellNeighbors.size());
		
 		for (int j = 0; j < cellNeighbors.size(); j++) {

			std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[cellNeighbors[j]]);

			// get intersection of commI and J
			std::vector<std::string> a = getComponentA(commI, commJ);

			// if there is complete overlap in I and J, then turnover is 1 and no need to
			// do calculations
			if (commI.size() == commJ.size() && commI.size() == a.size()) {
				cellVec[j] = 1.0;
			} else {

				double pdCombined = 0;
				double pdCommJ = 0;

				// get Faith's PD for commJ
				std::vector<int> edgesCommJ = FaithPD_branchIndices(commJ, phylo, nodeLeaves, spEdges);
				for (int k = 0; k < edgesCommJ.size(); k++) {
					pdCommJ = pdCommJ + edgeLengths[edgesCommJ[k]];
				}

				// numerator is sum of branches shared by commI and commJ
				std::vector<int> sharedIJ;
				// for each branch in cellI, is it present in cellJ?
				for (int k = 0; k < edgesCommI.size(); k++) { 
					if (std::find(edgesCommJ.begin(), edgesCommJ.end(), edgesCommI[k]) != edgesCommJ.end()) {
						sharedIJ.push_back(edgesCommI[k]);
					}
				}

				if (sharedIJ.size() > 0) {
					for (int i = 0; i < sharedIJ.size(); i++) {
						pdCombined = pdCombined + edgeLengths[sharedIJ[i]];
					}
				}

				cellVec[j] = pdCombined / (0.5 * (pdCommI + pdCommJ));
			}
		}
		out[i] = mean(cellVec);
	}

	return out;

}



// Calculate phylosor phylogenetic turnover
// [[Rcpp::export(name = calcPhylosor2, rng = false)]]
NumericVector calcPhylosor2(List spByCell, List phylo, bool showProgress) {
	
	// NumericVector out(spByCell.size());

	// std::vector<int> rasterValues = as< std::vector<int> >(rasterValuesR);
	// std::vector<int> nonNAcells = as< std::vector<int> >(nonNAcellsR);

	// extract relevant info from input data
	std::vector<std::string> tipLabels = as< std::vector<std::string> >(phylo["tip.label"]);
	std::vector<double> edgeLengths = as< std::vector<double> >(phylo["edge.length"]);
	
	// get leaves for all nodes as well as tip edges
	List nodeLeaves = getLeavesForNodes(phylo);
	List spEdges = getRootToTipEdges(phylo);

	// int n = spByCell.size();
	// Progress p(n, showProgress);
	
	std::vector<std::string> commI = as< std::vector<std::string> >(spByCell[0]);
	std::vector<std::string> commJ = as< std::vector<std::string> >(spByCell[1]);

	NumericVector cellVec(1);

	// get intersection of commI and J
	std::vector<std::string> a = getComponentA(commI, commJ);

	// if there is complete overlap in I and J, then turnover is 1
	if (commI.size() == commJ.size() && commI.size() == a.size()) {
		cellVec[0] = 1.0;
	} else {

		// // get mrca of commI
		// int mrcaI = getMRCA_from_nodeLeaves(nodeLeaves, commI);

		// // get mrca of commJ
		// int mrcaJ = getMRCA_from_nodeLeaves(nodeLeaves, commJ);

		// get mrca of combined commI and commJ			
		// std::vector<std::string> commIJ = getUnion(commI, commJ);
		// int mrcaIJ = getMRCA_from_nodeLeaves(nodeLeaves, commIJ);

		double pdCombined = 0;
		double pdCommI = 0;
		double pdCommJ = 0;

		// get Faith's PD for union(commI, commJ)
		// std::vector<int> edgeIndicesCombined = FaithPD_branchIndices(a, phylo, nodeLeaves, spEdges);

		// for (int i = 0; i < edgeIndicesCombined.size(); i++) {
		// 	pdCombined = pdCombined + edgeLengths[edgeIndicesCombined[i]];
		// }

		// get Faith's PD for commI and commJ separately
		std::vector<int> edgesCommI = FaithPD_branchIndices(commI, phylo, nodeLeaves, spEdges);
		std::vector<int> edgesCommJ = FaithPD_branchIndices(commJ, phylo, nodeLeaves, spEdges);

		// numerator is sum of branches shared by commI and commJ
		std::vector<int> sharedIJ;
		// for each species in cellI, is it present in cellJ?
		for (int k = 0; k < edgesCommI.size(); k++) { 
			if (std::find(edgesCommJ.begin(), edgesCommJ.end(), edgesCommI[k]) != edgesCommJ.end()) {
				sharedIJ.push_back(edgesCommI[k]);
			}
		}

		if (sharedIJ.size() > 0) {
			for (int i = 0; i < sharedIJ.size(); i++) {
				pdCombined = pdCombined + edgeLengths[sharedIJ[i]];
			}
		}

		for (int i = 0; i < edgesCommI.size(); i++) {
			pdCommI = pdCommI + edgeLengths[edgesCommI[i]];
		}

		for (int i = 0; i < edgesCommJ.size(); i++) {
			pdCommJ = pdCommJ + edgeLengths[edgesCommJ[i]];
		}


		// Rcout << "PD for commI: " << pdCommI << std::endl;
		// Rcout << "PD for commJ: " << pdCommJ << std::endl;
		// Rcout << "PD for commIJ: " << pdCombined << std::endl;

		cellVec[0] = pdCombined / (0.5 * (pdCommI + pdCommJ));
	}

	return cellVec;

}