# define ARMA_DONT_USE_WRAPPER
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>   
#include <vector>
#include <ctime>
#include "Gcopula_jointMAR.h"

using namespace Rcpp;
using namespace Numer;

//[[Rcpp::export]]
List MCMC_curerate(List & priors, arma::mat & Y, arma::mat & delta, arma::mat & X, arma::mat & dat, List & inits, List & slice_param, int & M, int & burnin, int & MHmcmc, int & thin, int & saving, double & accepttol, String & filename) {

  Rcout << "Importing priors, ";
  // ---------------------- PRIORS ------------------------- //
  double lambda_s0, lambda_r0, theta_s0, theta_r0, corr_var0;
  lambda_s0 = priors["lambda_s0"] ; lambda_r0 = priors["lambda_r0"] ; 
  theta_s0  = priors["theta_s0"]  ; theta_r0  = priors["theta_r0"]  ; 
  corr_var0 = priors["corr_var0"] ;
  arma::vec beta0 = priors["beta0"] ; arma::mat betaSig0 = priors["betaSig0"]  ;
  // ------------------------------------------------------- //

  Rcout << "data, ";
  // ---------------------------- DATA ------------------------------ //
  int k,n,p,pk,k_1;
  k = Y.n_cols; n = Y.n_rows; pk = X.n_cols; p = pk / k; k_1=k-1;
  // ---------------------------------------------------------------- //

  Rcout << "and initial values...\n";
  // ------------------------------------------- INITS ---------------------------------------------------- //
  arma::vec corr_r = inits["corr_r"] ; 
  arma::mat GammaCorrmat = inits["GammaCorrmat"];
  arma::vec beta = inits["beta"] ; arma::vec lambda = inits["lambda"] ; arma::vec theta = inits["theta"] ; 
  // ------------------------------------------------------------------------------------------------------ //
  
  Rcout << "Making bags for the posterior samples...\n";
  // ------------------ BAGS ------------------------ //
  arma::mat BETA(pk,M, arma::fill::zeros);
  arma::mat LAMBDA(k,M, arma::fill::zeros);
  arma::mat THETA(k,M, arma::fill::zeros);
  arma::mat CORR_R(1,M, arma::fill::zeros);

  arma::mat BETAMH(pk,MHmcmc, arma::fill::zeros);
  int mhind, mhindmu; mhind = mhindmu = 0;
  // ------------------------------------------------ //

  double wlong = slice_param["w"]; double slice_w_big = 20.0; double slice_w_small = 0.3; 
  double slice_m = slice_param["m"];

  int mcmcind=0; 
  int totalMCMC = M*thin + burnin + 2*MHmcmc+1; int burnmh1 = burnin + MHmcmc; int burnmh = burnin+2*MHmcmc; 

  double accpt,accptmu,accptind ; accpt=accptmu=accptind=0 ; 

  // ------------------------------- COMPUTING TIMES -------------------------------- //
  arma::vec timebeta(totalMCMC, arma::fill::zeros);
  arma::vec timelambda(totalMCMC, arma::fill::zeros);
  arma::vec timetheta(totalMCMC, arma::fill::zeros);
  arma::vec timecorr(totalMCMC, arma::fill::zeros);
  arma::vec timeNsupp(totalMCMC, arma::fill::zeros);
  // -------------------------------------------------------------------------------- //
  
  //-------- Lists for log-likelihoods ------//
  List params_beta   = List::create(Named("Y")=Y, _["delta"]=delta, _["X"]=X, _["gamma"]=theta, _["lambda"]=lambda, _["GammaCorrmat"]=GammaCorrmat, _["k"]=k, _["p"]=p, _["n"]=n, _["beta0"]=beta0, _["betaSig0"]=betaSig0);
  List params_lambda = List::create(Named("Y")=Y, _["delta"]=delta, _["X"]=X, _["gamma"]=theta, _["beta"]=beta,_["GammaCorrmat"]=GammaCorrmat, _["k"]=k, _["p"]=p, _["n"]=n, _["shape0"]=lambda_s0, _["rate0"]=lambda_r0);
  List params_theta  = List::create(Named("Y")=Y, _["delta"]=delta, _["X"]=X, _["lambda"]=lambda, _["beta"]=beta,_["GammaCorrmat"]=GammaCorrmat, _["k"]=k, _["p"]=p, _["n"]=n, _["shape0"]=theta_s0, _["rate0"]=theta_r0);
  List params_corr   = List::create(Named("Y")=Y, _["delta"]=delta, _["X"]=X, _["k"]=k, _["p"]=p, _["n"]=n,_["gamma"]=theta, _["lambda"]=lambda, _["beta"]=beta,_["sig20"]=corr_var0);
  
  String param_type_beta="beta_curerate"; String param_type_lambda="lambda_curerate"; String param_type_theta="theta_curerate"; String param_type_corr="corr_curerate";
  double lower_inf = -arma::datum::inf; double lower0 = 0.001 ; double lower_neg1 = -0.99 ; double upper_inf = arma::datum::inf ; double upper1 = 0.99;
  // --------------------------------------- //

  //-- Metropolis-Hastings Algorithm --//
  arma::mat VV(pk,pk); arma::mat MHvarbeta(pk,pk); arma::mat BetaMH(pk,MHmcmc,arma::fill::zeros);
  double accptbeta=0; double constbeta=1.0;
  int betaidx = 1; int mhindbeta = 0;
  // --------------------------------- //  

  Rcout << "Start MCMC\n\n";
  // MCMC SAMPLING
  for (int mcmc=0; mcmc<totalMCMC; mcmc++) {

    if (mcmc % saving == 0) {
      Rcout << mcmc << "\n"; 
      Rcout << "beta: "<<beta<<"\n"<< "lambda: "<<lambda<<"\n"<< "gamma: "<<theta<<"\n"<< "corr: "<<corr_r<<"\n\n";
      if (mcmc > 1) Rcout << "timebeta: " << timebeta(mcmc-1)<<"  timelambda: "<<timelambda(mcmc-1)<< "  timetheta: "<<timetheta(mcmc-1)<< "  timeNsupp: "<<timeNsupp(mcmc-1)<< "  timecorr: "<<timecorr(mcmc-1)<< "\n";
    }

     

    // --- BETA --- //
    time_t startb = time(0);

    if (MHmcmc==0 || mcmc<=burnmh1) {
      params_beta["lambda"]=lambda; params_beta["gamma"]=theta; params_beta["GammaCorrmat"]=GammaCorrmat;
      beta = slice_sample_cpp(param_type_beta,beta,wlong,slice_m,lower_inf,upper_inf,params_beta);
    } else {
      
        arma::vec mm = beta;
        // arma::mat vv(mm.n_elem,mm.n_elem); vv = VV(arma::span(start, end),arma::span(start, end));
        arma::vec betaparnew = arma::mvnrnd(mm,VV,1);
        arma::vec betanew    = betaparnew;
        double betai = beta(betaidx); double betanewi = betanew(betaidx);
        double lpold = betalog(betai, beta, betaidx, theta,lambda,GammaCorrmat,Y,delta,X,k,p,n,beta0,betaSig0);
        double lpnew = betalog(betanewi, betanew, betaidx, theta,lambda,GammaCorrmat,Y,delta,X,k,p,n,beta0,betaSig0);        
        arma::vec tmp1 = arma::ones(2); tmp1(1) = exp(lpnew-lpold); 
        double accprob = arma::min(tmp1);
        double runif = arma::randu<double>();
        if (runif < accprob) {
          accptbeta = accptbeta + 1;
          beta = betanew;
        }
      
     
      if (mcmc == burnmh) {
        double tmp;
        if (accptbeta/burnmh < 0.01)
          tmp = 1e-30;
        else 
          tmp = accptbeta/burnmh;
        constbeta = (constbeta*R::qnorm(accepttol/2, 0.0,1.0,TRUE,FALSE)) / R::qnorm(0.5*tmp, 0.0,1.0,TRUE,FALSE);

        Rcout << "accprob:" << accptbeta/burnmh << "\n" ;
        Rcout << "const:" << constbeta << "\n" ;

        VV = (constbeta*constbeta) * MHvarbeta;
      }      
    }
    if (MHmcmc > 0) {
      if ((mcmc > burnin) && (mcmc <= burnmh1)) {
        BetaMH(arma::span::all,mhindbeta) = beta;
        mhindbeta +=1;
      }
      if (mcmc==burnmh1) {
        MHvarbeta  = arma::cov(BetaMH.t());
        Rcout << "MHvarbeta" << MHvarbeta << "\n";
        VV = MHvarbeta;
      }
    }
    
  
    time_t endb = time(0);
    double timeb = endb-startb;


    // --- LAMBDA --- //
    time_t startl = time(0);

    params_lambda["GammaCorrmat"]=GammaCorrmat; params_lambda["gamma"]=theta; params_lambda["beta"]=beta;   
    lambda = slice_sample_cpp(param_type_lambda,lambda,wlong,slice_m,lower0,upper_inf,params_lambda);
    time_t endl = time(0);
    double timel = endl-startl;

    // --- THETA --- //
    params_theta["lambda"]=lambda; params_theta["GammaCorrmat"]=GammaCorrmat; params_theta["beta"]=beta;
    time_t startt = time(0);
    theta = slice_sample_cpp(param_type_theta,theta,wlong,slice_m,lower0,upper_inf,params_theta);
    
    time_t endt = time(0);
    double timet = endt-startt;

    // --- CORR --- //
    time_t startc = time(0);
    params_corr["lambda"]=lambda; params_corr["gamma"]=theta; params_corr["beta"]=beta;
    corr_r = slice_sample_cpp(param_type_corr,corr_r,wlong,slice_m,lower_inf,upper_inf,params_corr);
    GammaCorrmat = Gammacal(corr_r,k_1);

    time_t endc = time(0);
    double timec = endc-startc;
    
//    timebeta(mcmc)   = timeb ;
//    timeNsupp(mcmc)  = timeN ;
//    timelambda(mcmc) = timel ;
//    timetheta(mcmc) = timet ;
//    timecorr(mcmc)   = timec ;

    // --- Saving Data --- //
    int mcmcburnmh = mcmc-burnin-2*MHmcmc;
    if ((mcmc > burnmh) && (mcmcburnmh % thin == 0)) {
      BETA(arma::span::all,mcmcind) = beta ; LAMBDA(arma::span::all,mcmcind) = lambda ; THETA(arma::span::all,mcmcind) = theta ; CORR_R(arma::span::all,mcmcind) = corr_r;
      mcmcind +=1;
    }
  }

  List Result = List::create(Named("M")=M, 
    _["BETA"]=BETA.t(), _["LAMBDA"]=LAMBDA.t(), _["THETA"]=THETA.t(), _["CORR_R"]=CORR_R.t(), 
    _["dat"]=dat);   
  return Result;
}
