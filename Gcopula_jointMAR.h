#ifndef Gcopula_jointMAR_H
#define Gcopula_jointMAR_H

// Second Paper //
//Rcpp::List integrate_condnormal(double r12, double upper, double lower, double upper2, double lower2);
//double integrate_condnormal_double(double r12, double upper, double lower, double upper2, double lower2);

double doubleint(double & r12, double & upper1, double & upper2);

double survival_dist(arma::vec & gamma, arma::vec & lambda, arma::vec & beta,arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X, int & k,int & p,int & n);
double betalog(double & beta_i, arma::vec & betaold, int & ind, arma::vec & gamma, arma::vec & lambda, arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X, int & k,int & p,int & n, arma::vec & beta0, arma::mat & betaSig0);
double lambdalog(double & lambda_k, arma::vec & lambdaold, int & ind, arma::vec & gamma, arma::vec & beta,arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X, int & k,int & p,int & n, double & shape0, double & rate0);
double thetalog(double & theta_k, arma::vec & thetaold, int & ind, arma::vec & lambda, arma::vec & beta,arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X, int & k,int & p,int & n, double & shape0, double & rate0);
double corr_rlog(arma::mat & Y,arma::mat & delta, arma::mat & X, int & k,int & p,int & n, double & corr_r_i, arma::vec & corr_r, int & ind, arma::vec & gamma, arma::vec & lambda, arma::vec & beta,double & sig20);
arma::vec int_v_out(double & r12, double & var1_r2, arma::mat & upper, int & k, int & n);

// First Paper //
arma::mat Gammacal(arma::vec & r12, int & m);
arma::mat AR1cal(double & ar1, int & m, double & Gammastr_ar1);
//double survcdflog_cpp(       arma::mat & Xcenm, arma::vec & Xcent, int & m, int & n, int & M, int & J, arma::mat & Gamma, arma::vec & nobs, arma::uvec & uvec0,arma::uvec & IMcolnum,arma::uvec & IMcenidx,arma::uvec & IMcenfullyobserved);
//double obsliklog_cpp(        arma::mat & Xobs, 					   int & m, 		 int & M, int & J, arma::mat & Gamma, arma::vec & nobs, 				   arma::uvec & IMcolnum0,arma::uvec & IMobsidx,arma::uvec & IMobsfullyobserved);
double survcdflog_cpp(       arma::mat & Xcenm, arma::vec & Xcent, int & m, int & n, int & M, int & J, arma::mat & Gamma, Rcpp::List & idx); //arma::vec & nobs, arma::uvec & uvec0,arma::uvec & IMcolnum,arma::uvec & IMcenidx,arma::uvec & IMcenfullyobserved);
double obsliklog_cpp(        arma::mat & Xobs, 					   int & m, 		 int & M, int & J, arma::mat & Gamma, Rcpp::List & idx); //arma::vec & nobs, 				   arma::uvec & IMcolnum0,arma::uvec & IMobsidx,arma::uvec & IMobsfullyobserved);
double censliklog_cpp(       arma::mat & Xcenm, 				   int & m, 		 int & M, int & J, arma::mat & Gamma, Rcpp::List & idx);
Rcpp::List r12postdist_cpp(  arma::mat & Y, int & n, int & m, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, double & r, arma::vec & r_old, int & ind, double & sigr2, double & mh, 																		  													 Rcpp::List & idx, arma::uvec & cenIM, arma::uvec & obsIM);
Rcpp::List Gammapostdist_cpp(arma::mat & Y, arma::vec & v, int & n, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, double & r, arma::vec & r_old, int & ind, double & sigr2, double & ar1, arma::vec & alpha, double & mh, 								  													 Rcpp::List & idx, arma::uvec & cenIM, arma::uvec & obsIM, double & Gammastr_ar1 );
double missingliklog_cpp(    arma::mat & Yorigin, double & ymiss, int & ind, arma::vec & missind, arma::vec & Ximp, arma::vec & v, arma::vec & Xt, arma::vec & Xcent, arma::cube & W, int & n, int & m, int & M, int & J, arma::mat & Gamma, arma::vec & mu, arma::vec & sig2,  						  													 Rcpp::List & idx);
double ar1postdist_cpp(      arma::mat & Y, arma::vec & v, int & n, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, arma::vec & r12, double & ar11, double & ar12, double & ar1, arma::vec & alpha, 														  													 Rcpp::List & idx, double & ar1priorBETA, arma::uvec & cenIM, arma::uvec & obsIM, double & Gammastr_ar1);
double alphapostdist_cpp(    arma::mat & Y, arma::vec & v, int & n, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, arma::vec & r12, double & ar1, double & alphaa, arma::vec & alpha_old, int & ind, arma::vec & alpha0, arma::mat & Sigalpha0, double & mh, 													 Rcpp::List & idx, arma::uvec & cenIM, arma::uvec & obsIM , double & Gammastr_ar1);
double mupostdist_cpp(       arma::mat & Y, arma::vec & v, int & n, int & m, int & M, int & J, arma::vec & Xt, double & muu, arma::vec & mu_old, int & ind, arma::vec & mu0, arma::mat & Sig0, arma::vec & sig2, arma::mat & Gamma, arma::cube & W, double & mh, 										  													 Rcpp::List & idx);
double sig2postdist_cpp(     arma::mat & Y, arma::vec & v, int & n, int & m, int & M, int & J, arma::vec & Xt, double & sig2, double & c0, double & d0, arma::vec & mu, arma::mat & Gamma, arma::cube & W, 														 										  													 Rcpp::List & idx);
double sig2multpostdist_cpp( arma::mat & Y, arma::vec & v, int & n, int & m, int & M, int & J, arma::vec & Xt, double & sig22, arma::vec & sig2_old, int & ind, double & c0, double & d0, arma::vec & mu, arma::mat & Gamma, arma::cube & W, 																												 Rcpp::List & idx);
double lambdapostdist_cpp(   arma::vec & Time, arma::mat & Z, arma::vec & v, arma::vec & delta, int & n, int & m, int & M, int & J, double & Xtmax, double & Xtmin, arma::mat & Xm, arma::mat & Xcenm, double & lambdaa, arma::vec & lambda_old, int & ind, arma::vec & timeint, double & a0, double & b0, arma::vec & beta, arma::mat & Gamma, 			 Rcpp::List & idx);
double betapostdist_cpp(     arma::vec & Time, arma::mat & Z, arma::vec & v, arma::vec & delta, int & n, int & m, int & M, int & J, double & Xtmax, double & Xtmin, arma::mat & Xm, arma::mat & Xcenm, double & betaa, arma::vec & beta_old, int & ind, arma::vec & timeint, arma::vec & beta0, arma::mat & Sigbeta0, arma::vec & lambda, arma::mat & Gamma, Rcpp::List & idx);



arma::vec slice_sample_cpp(Rcpp::String & param_type, arma::vec & x0, double & w, double & slice_m, double & lower, double & upper, Rcpp::List & params);

#endif 