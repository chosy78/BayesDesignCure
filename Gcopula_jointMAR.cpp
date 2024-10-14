# define ARMA_DONT_USE_WRAPPER
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>   
#include <vector>

using namespace Rcpp;
using namespace Numer;

/******************/
/* BASE FUNCTIONS */
/******************/
//[[Rcpp::export]]
arma::mat Gammacal(arma::vec & r12, int & m) {
  
  arma::mat R = arma::mat(m+1,m+1, arma::fill::eye);
  
  int ind = 0;
  for (int j = 0; j < m+1; j++)
  {
    for (int i = 0; i < m+1; i++) 
    {
      if (i < j)
      {
        R(i,j) = r12[ind];
        ind ++;
      }
    } 
  }
  
  arma::mat Sigma = inv(R.t() * R);
  arma::vec diag  = 1/sqrt(Sigma.diag());
  arma::mat Diag  = diagmat(diag);
  arma::mat Gamma = Diag * Sigma * Diag;
  
  // arma::mat Gammainv    = inv(Gamma);
  // arma::mat Gammamm     = Gamma(arma::span(1,m),arma::span(1,m));
  // arma::mat Gammamminv = inv(Gammamm);
  
  //List Result = List::create(Named("Gamma")=Gamma, _["Gammainv"]=Gammainv, _["Gammamm"]=Gammamm, _["Gammamminv"]=Gammamminv);
  
  return Gamma;
}

//[[Rcpp::export]]
arma::mat AR1cal(double & ar1, int & m, double & Gammastr_ar1) {
  arma::mat R = arma::mat(m,m,arma::fill::zeros);
  arma::mat RR = arma::mat(m,m,arma::fill::zeros);
  arma::mat Rho = arma::mat(m,m) ; Rho.fill(ar1);
  arma::mat exponent;
  
  if (Gammastr_ar1==1.0) {
    for (int i=1; i<m; i++) {
      arma::vec a(m) ; arma::rowvec b(m) ; a.fill(i); b.fill(i);
      R.col(i) = a;
      RR.row(i) = b;
    }
    exponent = exp(abs(R-RR)%log(Rho));
  } else {
    exponent = Rho; exponent.diag().ones();
  }
  return(exponent);
}

/*********************************/
/* SECOND PAPER RELATED FUNCTION */
/*********************************/
class conditionalnormalIntegrand: public Func
{
private:
  double r12;
  double upper;
  double lower;
public:
  conditionalnormalIntegrand(double r12_, double upper_) : r12(r12_), upper(upper_) {}

  double operator()(const double & x) const
  {
    return 1/(sqrt(2*arma::datum::pi)) * exp(-0.5*x*x) * (R::pnorm((upper-r12*x)/(sqrt(1-r12*r12)), 0.0,1.0,TRUE,FALSE));  
    //return (R::pnorm((upper-r12*x)/(sqrt(1-r12*r12)), 0.0,1.0,TRUE,FALSE) - R::pnorm((lower-r12*x)/(sqrt(1-r12*r12)), 0.0,1.0,TRUE,FALSE));  
  }
};

class intcont
{
public:
  double res;
  double err_est;
  int err_code;

  intcont(double res_, double err_est_, int err_code_) : res(res_), err_est(err_est_), err_code(err_code_) {}
};
 
intcont integrate_condnormal_class(double r12, double upper, double upper2)
{
  const double ul = upper2, ll=R_NegInf; 
 
  conditionalnormalIntegrand f(r12,upper);
  double err_est;
  int err_code;
  double res = integrate(f, ll, ul, err_est, err_code);
 
  intcont xxx(res,err_est,err_code);         
  return xxx;
}   

//[[Rcpp::export]]
arma::vec int_v_out(double & r12, double & var1_r2, arma::mat & upper, int & k, int & n) {

  arma::vec result(n,arma::fill::zeros); 
  arma::vec tmp = arma::regspace(-5,6); arma::mat intsum(tmp.n_elem,1); intsum.col(0) = tmp; 

  for (int i=0; i<n; i++){  
    double c,it,fpc,gpc,u_rc,l_rc,k0,k1,k2,k3,k4,m0,m1,m2,m3,m4;
    
    if (upper(i,k-1)<-5.0) {
      result(i)=0.0;
    } else {
      arma::mat iintsum = intsum.rows(arma::find(intsum<upper(i,k-1))); //int tmpn = tmp.n_elem+1;
      
      if (upper(i,k-1)<intsum.max()) {
        arma::mat ttmp(1,1); ttmp.fill(upper(i,k-1));
        iintsum = join_cols(iintsum,ttmp); 
      } 
      
      for(int j=0; j<iintsum.size(); j++) {
        it = iintsum(j); c = it-0.5; 
        if (it==iintsum(iintsum.n_rows-1,0)) c = 0.5*(iintsum(iintsum.n_rows-1,0)+iintsum(iintsum.n_rows-2,0));
        if (j==0) {
          l_rc = -6.0;
        } else {
          l_rc = iintsum(j-1);
        }
      
        m0 = R::pnorm(it,0.0,1.0,TRUE,FALSE) - R::pnorm(l_rc,0.0,1.0,TRUE,FALSE); 
        m1 = -(R::dnorm(it,0.0,1.0,FALSE) - R::dnorm(l_rc,0.0,1.0,FALSE)); 
        m2 = -(R::dnorm(it,0.0,1.0,FALSE)*it-R::dnorm(l_rc,0.0,1.0,FALSE)*l_rc) + m0;
        m3 = -(R::dnorm(it,0.0,1.0,FALSE)*pow(it,2)-R::dnorm(l_rc,0.0,1.0,FALSE)*pow(l_rc,2)) + 2.0*m1;   
        m4 = -(R::dnorm(it,0.0,1.0,FALSE)*pow(it,3)-R::dnorm(l_rc,0.0,1.0,FALSE)*pow(l_rc,3)) + 3.0*m2;      
        
        fpc = -r12/sqrt(2*arma::datum::pi*var1_r2)*exp(-(pow(upper(i,0)-r12*c,2))/(2*var1_r2));
        u_rc = upper(i,0)-r12*c; 

        k0 = R::pnorm(u_rc/sqrt(var1_r2),0.0,1.0,TRUE,FALSE);
        k1 = fpc;
        k2 = r12/(2*var1_r2)*(u_rc*fpc);
        k3 = r12*r12/(6*var1_r2)*((u_rc*u_rc-var1_r2)*fpc/var1_r2);
        k4 = -pow(r12/var1_r2,3)/24*(3*var1_r2-pow(u_rc,2))*u_rc*fpc;
      
        result(i) += k0*m0 + k1*(m1-c*m0) + k2*(m2-2*c*m1+pow(c,2)*m0) + k3*(m3-3*c*(m2-c*m1)-pow(c,3)*m0) + k4*(m4-4*c*m3+6*pow(c,2)*m2-4*pow(c,3)*m1+pow(c,4)*m0);
      }

    }
  }
  return result;
}

//[[Rcpp::export]]
double survival_dist(arma::vec & gamma, arma::vec & lambda, arma::vec & beta, arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X,int & k,int & p,int & n) {
  
/*  if (X.n_rows!=Y.n_rows) {
    Rcout << "X:"<< X.n_rows << "Y:" << Y.n_rows<<"\n";
    stop("CSY: The data size is different!");
  } */

  double r12 = GammaCorrmat(0,k-1); double var1_r2 = 1-r12*r12; 
  arma::mat theta(n,k); 
  for (int kk=0; kk<k; kk++)
    theta.col(kk) = exp(X(arma::span::all,arma::span(kk*p,((kk+1)*p)-1))*beta(arma::span(kk*p,((kk+1)*p)-1)));

  double loglik;

  //Rcout  << "before calculation \n";
  
  arma::mat Fcdfuncured(n,k);      Fcdfuncured      = 1 - exp(-arma::repelem(gamma.t(),n,1) % pow(Y.each_row(),lambda.t()));
  arma::mat fdensityuncured(n,k); fdensityuncured = arma::repelem((gamma % lambda).t(),n,1) % (pow(Y.each_row(),lambda.t()-1) % (1-Fcdfuncured));
  arma::mat Fcdf(n,k); Fcdf = 1-(exp(-theta % Fcdfuncured)); 
  arma::mat ul(n,k); ul = 1-exp(-theta); 
  arma::mat fdensity(n,k); fdensity = theta % fdensityuncured % (1-Fcdf);
  arma::vec CopulaCDF(n); arma::vec logcopulaftn(n);
  arma::mat uppernum(n,k);

  arma::vec condpnorm0(n), condpnorm1(n); 
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<k; j++) {
      //pnormupper(i,j) = R::pnorm(ul(i,j), 0.0,1.0,TRUE,FALSE); plnormupper(i,j) = R::pnorm(ul(i,j),0.0,1.0,TRUE,TRUE);
      uppernum(i,j) = R::qnorm(Fcdf(i,j), 0.0,1.0,TRUE,FALSE); //- R::qnorm(0.5 * (1-exp(-theta(i,j))), 0.0,1.0,TRUE,FALSE);
      //upper(i,j) = R::qnorm(ul(i,j)*pnormupper(i,j), 0.0,1.0,TRUE,FALSE);
    }
    condpnorm0(i) = R::pnorm((uppernum(i,1)-r12*uppernum(i,0))/sqrt(var1_r2),0.0,1.0,TRUE,FALSE); 
    condpnorm1(i) = R::pnorm((uppernum(i,0)-r12*uppernum(i,1))/sqrt(var1_r2),0.0,1.0,TRUE,FALSE);
  }
  
  arma::vec numdoubleint(n); numdoubleint = int_v_out(r12,var1_r2,uppernum,k,n);
//  arma::vec dendoubleint(n); dendoubleint = int_v_out(r12,var1_r2,upper,k,n);
  
  if (!numdoubleint.is_finite()) {
      arma::uvec idx = arma::find_nonfinite(numdoubleint);
    for (int i=0;i<idx.n_elem;i++) {
      numdoubleint(idx(i)) = integrate_condnormal_class(r12,uppernum(idx(i),0),uppernum(idx(i),k-1)).res;
    }
  }
//  if (!dendoubleint.is_finite()) {
//      arma::uvec idx = arma::find_nonfinite(dendoubleint);
//    for (int i=0;i<idx.n_elem;i++) {
//      dendoubleint(idx(i)) = integrate_condnormal_class(r12,upper(idx(i),0),upper(idx(i),k-1)).res;
//    }
//  }

  CopulaCDF = numdoubleint;///dendoubleint;

  arma::vec copuladens(n); copuladens = -0.5*(1/var1_r2)*(r12*r12*sum(pow(uppernum,2),1)-2*r12*prod(uppernum,1));
  logcopulaftn = -0.5*log(var1_r2) + copuladens + sum(log(fdensity),1);// + sum(plnormupper,1) - log(dendoubleint);
  //logcopulaftn = log(2*arma::datum::pi) + copuladens + sum(log(fdensity),1) - log(dendoubleint);

  arma::uvec obsobs,cenobs,obscen,cencen; 
  obsobs = arma::find(sum(delta,1)==k); cencen = arma::find(sum(delta,1)==0); 
  cenobs = arma::find(sum(delta,1)==1 && delta(arma::span::all,0)==0) ; obscen = arma::find(sum(delta,1)==1 && delta(arma::span::all,0)==1);
  arma::vec obsobsvec, cenobsvec, obscenvec, cencenvec; arma::vec tmp; arma::vec tmpcencenvec;
  
  // obsobs
  obsobsvec = logcopulaftn(obsobs);
  
  //cenobs
  //tmp = log(fdensity(arma::span::all,k-1)) - log(dendoubleint) + log(dendoubleint - pnormupper.col(k-1)%condpnorm1); cenobsvec = tmp(cenobs);
  tmp = log(fdensity(arma::span::all,k-1)) + log(1 - condpnorm1); cenobsvec = tmp(cenobs);
  
  //obscen
  //tmp = log(fdensity(arma::span::all,0)) - log(dendoubleint) + log(dendoubleint - pnormupper.col(0)%condpnorm0); obscenvec = tmp(obscen);
  tmp = log(fdensity(arma::span::all,0)) + log(1 - condpnorm0); obscenvec = tmp(obscen);

  //cencen
  tmp = -sum(Fcdf,1) + 1 + CopulaCDF; cencenvec = tmp(cencen); //tmp(tmpidx).fill(1.0); 
  
  //if (!cencenvec.is_finite()) Rcout<<"D:";//<<(tmp(cencennegidx));;
  
  loglik = accu(obsobsvec(arma::find_finite(obsobsvec))) + accu(cenobsvec(arma::find_finite(cenobsvec))) + accu(obscenvec(arma::find_finite(obscenvec))) + accu(log(cencenvec(arma::find(cencenvec>0))));

  return loglik;
}

//[[Rcpp::export]]
double betalog(double & beta_i, arma::vec & betaold, int & ind, arma::vec & gamma, arma::vec & lambda, arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X,int & k,int & p,int & n, arma::vec & beta0, arma::mat & betaSig0) {

  arma::vec beta = betaold ; beta(ind-1) = beta_i ;
  // Prior
  double prior = - 0.5 * accu((beta-beta0).t()*betaSig0*(beta-beta0)) ;
  // Likelihood
  double logLik = survival_dist(gamma, lambda, beta, GammaCorrmat, Y, delta, X,k,p,n);
  double result = logLik + prior;
  return result;
}

//[[Rcpp::export]]
double lambdalog(double & lambda_k, arma::vec & lambdaold, int & ind, arma::vec & gamma, arma::vec & beta, arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X,int & k,int & p,int & n, double & shape0, double & rate0) {
  arma::vec lambda = lambdaold ; lambda(ind-1) = lambda_k ; 
  // Prior
  double prior = (shape0-1)*log(lambda_k) - rate0*lambda_k;
  // Likelihood
  double loglik = survival_dist(gamma, lambda, beta, GammaCorrmat, Y, delta, X,k,p,n);
  
  double result = prior + loglik;
  return result ;
}

//[[Rcpp::export]]
double thetalog(double & theta_k, arma::vec & thetaold, int & ind, arma::vec & lambda, arma::vec & beta, arma::mat & GammaCorrmat, arma::mat & Y, arma::mat & delta, arma::mat & X,int & k,int & p,int & n, double & shape0, double & rate0) {
  arma::vec theta = thetaold ; theta(ind-1) = theta_k ; 
  // Prior
  double prior = (shape0-1)*log(theta_k) - rate0*theta_k;
  // Likelihood
  double loglik = survival_dist(theta, lambda, beta, GammaCorrmat, Y, delta, X,k,p,n);
  
  double result = prior + loglik;
  return result ;
}

//[[Rcpp::export]]
double corr_rlog(arma::mat & Y,arma::mat & delta,arma::mat & X, int & k,int & p,int & n, double & corr_r_i, arma::vec & corr_r, int & ind, arma::vec & gamma, arma::vec & lambda, arma::vec & beta, double & sig20) {
  // Prior
  arma::vec r = corr_r; r(ind-1) = corr_r_i; int kk = k-1;
  double prior = -accu(corr_r % corr_r) / (2*sig20);
  arma::mat GammaCorrmat = Gammacal(r, kk);
  // Likelihood
  double loglik = survival_dist(gamma, lambda, beta, GammaCorrmat, Y, delta, X,k,p,n);

  double result = prior + loglik;
  return result;
}

//[[Rcpp::export]]
arma::vec slice_sample_cpp(String & param_type, arma::vec & x0, double & w, double & slice_m, double & lower, double & upper, List & params) { 
  
  double J, K;
  arma::vec L(x0.size()), R(x0.size()), param_new(x0.size());
  arma::vec paramvec = param_new = x0;

  //List idx = params["idx"]; 
  //arma::vec nobsobs = idx["nobsobs"] ; arma::vec nobscen = idx["nobscen"]; arma::uvec uvec0=idx["uvec0"];arma::uvec IMcolnum=idx["IMcolnum"];arma::uvec IMuvec=idx["IMuvec"];arma::uvec IMuvecfullyobserved=idx["IMuvecfullyobserved"];   

  for (int j=0; j < x0.size(); j++) {

    double logy=0.0; double logz=0.0;

    double param_old = paramvec(j); 
    int jj = j + 1;

    // Calculate initial horizontal interval;
    L(j) = param_old - R::runif(0.0, w);
    R(j) = L(j) + w;
    // Truncate bounds to support of the parameter space;
    L(j) = std::max(L(j),lower);
    R(j) = std::min(R(j),upper);
    // Step out;
    J = floor(slice_m * R::runif(0.0,1.0));
    K = (slice_m-1)-J;

    double funL, funR, funnew;
    
    if (param_type=="beta_curerate") { // Second Paper from here //

      arma::vec gamma=params["gamma"]; arma::vec lambda=params["lambda"]; arma::mat GammaCorrmat=params["GammaCorrmat"];
      arma::mat Y=params["Y"]; arma::mat delta=params["delta"];  arma::mat X=params["X"];
      int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
      arma::vec beta0=params["beta0"]; arma::mat betaSig0=params["betaSig0"]; 
     
      logy = betalog(param_old, paramvec, jj, gamma,lambda,GammaCorrmat,Y,delta,X,k,p,n,beta0,betaSig0);
      funL = betalog(L(j), paramvec, jj, gamma,lambda,GammaCorrmat,Y,delta,X,k,p,n,beta0,betaSig0);
      funR = betalog(R(j), paramvec, jj, gamma,lambda,GammaCorrmat,Y,delta,X,k,p,n,beta0,betaSig0);
     
     } else if (param_type=="lambda_curerate") {
      arma::vec gamma=params["gamma"]; arma::vec beta=params["beta"]; arma::mat GammaCorrmat=params["GammaCorrmat"];
      arma::mat Y=params["Y"]; arma::mat delta=params["delta"];  arma::mat X=params["X"];
      int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
      double shape0=params["shape0"]; double rate0=params["rate0"]; 

      logy = lambdalog(param_old, paramvec, jj, gamma,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);
      funL = lambdalog(L(j), paramvec, jj, gamma,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);
      funR = lambdalog(R(j), paramvec, jj, gamma,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);  

     } else if (param_type=="theta_curerate") {
      arma::vec lambda=params["lambda"]; arma::vec beta=params["beta"]; arma::mat GammaCorrmat=params["GammaCorrmat"];
      arma::mat Y=params["Y"]; arma::mat delta=params["delta"];  arma::mat X=params["X"];
      int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
      double shape0=params["shape0"]; double rate0=params["rate0"]; 

      logy = thetalog(param_old, paramvec, jj, lambda,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);
      funL = thetalog(L(j), paramvec, jj, lambda,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);
      funR = thetalog(R(j), paramvec, jj, lambda,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);  

     } else if (param_type=="corr_curerate") {
      arma::mat Y=params["Y"]; arma::mat delta=params["delta"]; arma::mat X=params["X"];
      int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
      arma::vec lambda=params["lambda"]; arma::vec gamma=params["gamma"]; arma::vec beta=params["beta"]; 
      double sig20=params["sig20"];

      logy = corr_rlog(Y,delta,X,k,p,n,param_old,paramvec, jj, gamma,lambda,beta,sig20);
      funL = corr_rlog(Y,delta,X,k,p,n,L(j),paramvec, jj, gamma,lambda,beta,sig20);
      funR = corr_rlog(Y,delta,X,k,p,n,R(j),paramvec, jj, gamma,lambda,beta,sig20);
     }

    // draw uniformly from [0, y]
    logz = logy - R::rexp(1.0);
      
    while ( J > 0 && L(j) > lower && funL > logz )
    {
      L(j) = L(j) - w; 
      if (L(j) <= lower) 
        L(j) = lower;
      J = J-1;
    }
    while ( K > 0 && R(j) < upper && funR > logz )
    {
      R(j) = R(j) + w;
      if (R(j) >= upper) 
        R(j) = upper;
      K = K-1;
    }

    // shrinkage procedure
    int cnt = 0;
    do {
      cnt++;
      paramvec(j) = R::runif(L(j),R(j));

      
      if (param_type=="beta_curerate") { // Second Paper from here //

        arma::vec gamma=params["gamma"]; arma::vec lambda=params["lambda"]; arma::mat GammaCorrmat=params["GammaCorrmat"];
        arma::mat Y=params["Y"]; arma::mat delta=params["delta"];  arma::mat X=params["X"];
        int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
        arma::vec beta0=params["beta0"]; arma::mat betaSig0=params["betaSig0"]; 
       
        funnew = betalog(paramvec(j), paramvec, jj, gamma,lambda,GammaCorrmat,Y,delta,X,k,p,n,beta0,betaSig0);
      
     } else if (param_type=="lambda_curerate") {

      arma::vec gamma=params["gamma"]; arma::vec beta=params["beta"]; arma::mat GammaCorrmat=params["GammaCorrmat"];
      arma::mat Y=params["Y"]; arma::mat delta=params["delta"];  arma::mat X=params["X"];
      int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
      double shape0=params["shape0"]; double rate0=params["rate0"]; 

      funnew = lambdalog(paramvec(j), paramvec, jj, gamma,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);

     } else if (param_type=="theta_curerate") {

      arma::vec lambda=params["lambda"]; arma::vec beta=params["beta"]; arma::mat GammaCorrmat=params["GammaCorrmat"];
      arma::mat Y=params["Y"]; arma::mat delta=params["delta"];  arma::mat X=params["X"];
      int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
      double shape0=params["shape0"]; double rate0=params["rate0"]; 

      funnew = thetalog(paramvec(j), paramvec, jj, lambda,beta,GammaCorrmat,Y,delta,X,k,p,n,shape0,rate0);

     } else if (param_type=="corr_curerate") {
      arma::mat Y=params["Y"]; arma::mat delta=params["delta"]; arma::mat X=params["X"];
      int k=params["k"]; int p=params["p"] ; int n=params["n"]; 
      arma::vec lambda=params["lambda"]; arma::vec gamma=params["gamma"]; arma::vec beta=params["beta"]; 
      double sig20=params["sig20"];

      funnew = corr_rlog(Y,delta,X,k,p,n,paramvec(j),paramvec, jj, gamma,lambda,beta,sig20);
     }
      
      if ( funnew > logz )
        break;
      if ( paramvec(j) < param_old )
        L(j) = paramvec(j);
      else
        R(j) = paramvec(j);
    } while (cnt<1e4);
    if (cnt==1e4) {
        std::string print = param_type; std::cout << print << "  funnew: "<<funnew<<"  logz:"<<logz<<"\n";
        ::Rf_error("slice_sample_cpp loop did not finish");
      }

      if (-0.0000000001 <= L(j) - R(j) && L(j) - R(j) <= 0.0000000001)
        paramvec(j) = 0.5 * (L(j) + R(j));


    }

  return paramvec;
}