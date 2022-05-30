#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cstdarg>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cum_min(NumericVector x){
    int n = x.size();
    if (n <2 ) return x;
    NumericVector out(n);
    out[0] = x[0];
    bool na_found = false;
    for(int i =1; i<n;  ++i) {
        if(na_found || R_IsNA(x[i])){
         na_found = true;
         out[i] = NA_REAL;
        }else{
        out[i] = std::min(out[i-1], x[i]);
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cum_max(NumericVector x){
    int n = x.size();
    if (n <2 ) return x;
    NumericVector out(n);
    out[0] = x[0];
    for(int i =1; i<n;  ++i) {
        out[i] = std::max(out[i-1], x[i]);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cum_prod(NumericVector x){
    int n = x.size();
    NumericVector out(n);
    out[0] = x[0];
    for(int i =1; i<n;  ++i) {
        out[i] = out[i-1] *  x[i];
    }
    return out;
}

// [[Rcpp::export]]
NumericVector diffC(NumericVector x, int lag = 1, int differences = 1)
{
    int len = x.size();
    if (lag >= len) return NumericVector(0);
    NumericVector out(len - lag);
    for(int i = 0; i< len - lag;  ++i ){
        out[i] = x[i + lag] - x[i];
    }
    if (differences == 1){
        return out;
    }else{
        differences = differences - 1;
        return diffC(out, lag, differences);
    }
}

// [[Rcpp::export]]
NumericVector rangeC(NumericVector x, bool na_rm)
{
   if (x.size() == 1) return NumericVector::create(x[0], x[0]);
   NumericVector::iterator it;
// Initialize with biggest, smallest possible values as min and max
   NumericVector out = NumericVector::create(R_PosInf, R_NegInf);
   double cur;
   for(it = x.begin(); it != x.end(); ++it){
        cur = *it;
        if(cur == NA_REAL){
            if(na_rm){
                continue;
            }else{
                return NumericVector::create(NA_REAL, NA_REAL);
            }
        }
        out[0] = std::min(out[0], cur);
        out[1] = std::max(out[1], cur);
   }
       return out;
}

// [[Rcpp::export]]
double varC(NumericVector x, bool na_rm = false)
{

   int n = x.size();
   if (n < 2) return NA_REAL;
   double xbar_last = x[0];
   double xbar_cur = NA_REAL;
   double M2_last = 0;
   double M2_cur;
   int nas = 0;

   for(int i =0; i <n;  ++i){
       if(R_IsNA(x[i])){
           if(na_rm){
               nas = nas + 1;
               continue;
           }else{
               return NA_REAL;
           }
       }
       //First non-NA
       if(R_IsNA(xbar_last)){
        xbar_last = x[i];
        continue;
       }

       xbar_cur = xbar_last + ( x[i] - xbar_last )/(i - nas + 1); //subtract number of NAs skipped
       M2_cur = M2_last +(x[i] - xbar_last) * (x[i] - xbar_cur);
       xbar_last = xbar_cur;
       M2_last = M2_cur;
   }
   int valid = n - 1 - nas;
   double out;
   // If 1 valid value or less, return NA, since variance undefined
   if (valid > 0){
    out = M2_last/valid;
   }else{
    out = NA_REAL;
   }
   return out;
}

// [[Rcpp::export]]
double medianC( std::vector<double> x, bool na_rm = false){
    int n = x.size();
    std::vector<double> filtered;
    // See https://stackoverflow.com/questions/21204676/modern-way-to-filter-stl-container
    std::copy_if (x.begin(), x.end(), std::back_inserter(filtered), [](double i){return !R_IsNA(i);} );
    if(!na_rm && x.size() > filtered.size()) return NA_REAL;
    int which = filtered.size() / 2;
    std::partial_sort(filtered.begin(), filtered.begin() + which, filtered.end());
    if(n % 2 == 0){
        return (filtered[which -1] + filtered[which] ) / 2;
    }else{
        return filtered[which];

    }
}


// [[Rcpp::export]]
LogicalVector inC(std::vector<double> x, std::vector<double> y){
    int n = x.size();
    std::unordered_set<double> table(y.begin(), y.end());
    LogicalVector out(n);


    for(int i = 0; i < n;  ++i){
        out[i] = bool(table.count(x[i]));
        //R_isnancpp
        //if (na_present && R_IsNA(x[i] )) out[i] = true;
    }
    return out;
}

//template <int RTYPE>
// [[Rcpp::export]]
std::vector<double> UniqueCC(std::vector<double> x, std::vector<double> incomparables ){
    std::unordered_set<double> inc(incomparables.begin(), incomparables.end());
    std::unordered_set<double> seen;
    std::vector<double> out;
    double cur;
    for(int i = 0; i < x.size(); ++ i  ){
        cur = x[i];
        if(seen.insert(cur).second || bool(inc.count(cur))){
out.push_back(cur); //Add if new or marked incomparable
}
    }
    return out;

}
// [[Rcpp::export]]
int whichmaxC(std::vector<double> x){
    return std::distance(x.begin(), std::max_element(x.begin(), x.end())) + 1;
}


// [[Rcpp::export]]
std::vector<int> intersectC(std::vector<int> x, std::vector<int> y){
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    std::vector<int> out;
    std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(out));
    return out;
}

// [[Rcpp::export]]
std::vector<double> unionC(std::vector<double> x, std::vector<double> y){
    //Actually this is how R implements it
    x.insert(x.end(), y.begin(), y.end());
    return UniqueCC(x, std::vector<double>());
}

// [[Rcpp::export]]
std::vector<int> setdiffC(std::vector<int> x, std::vector<int> y){
    std::unordered_set<int> table(y.begin(), y.end());
    std::vector<int> out;
    for(int i=0; i < x.size();  ++ i){
        if(!bool(table.count(x[i]))) out.push_back(x[i]);
    }

    return out;
}

// [[Rcpp::export]]
// Variadic args in C++ are a nightmare, sue me
double maxC_impl(std::list< std::vector<double> > x, bool na_rm = false){
    double cur_max = R_NegInf;
    std::list< std::vector<double> >::iterator it;
    for(it = x.begin(); it != x.end();  ++it){
    if(!na_rm && *std::min_element(it->begin(), it->end()) == -2147483648){
        return NA_REAL;
    }
    //double res = *std::min_element(it->begin(), it->end());
    //std::cout << res << "\n";
        double this_max = *std::max_element(it->begin(), it->end());

        cur_max = std::max(cur_max, this_max);
    }
    return cur_max;
}
