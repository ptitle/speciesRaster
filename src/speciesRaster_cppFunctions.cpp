#include <Rcpp.h>
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





