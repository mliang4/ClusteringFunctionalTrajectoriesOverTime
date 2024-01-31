#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Mathematical constants
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

namespace help{

////////////////////////////////////
// Helper functions from PGBVS /////
////////////////////////////////////

// Generate exponential distribution random variates
double exprnd(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables 
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = help::exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss(double z, double t)
{
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / help::truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = help::randinvg(mu);
    }
  }    
  return X;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1) 
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + help::exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = help::tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = help::aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1) 
    {
      Sn = Sn + asgn * help::aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// Simulate MVT normal data
arma::mat mvrnormArma( int n, arma::vec mu, arma::mat sigma ) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn( n, ncols );
  return arma::repmat( mu, 1, n ).t()  + Y * arma::chol( sigma );
}

// LSE trick to help with underflow, eg: when calculating p=exp(-1000,-1000,-1000), without LSE p would become (0,0,0) but it should be (1/3,1/3,1/3)
double logsumexp(arma::vec x){
  double c = max(x);
  return c + log(sum(exp(x-c)));
}


// // Algorithm to sample from truncated normal
// /// Check if simpler subalgorithm is appropriate.
// inline bool CheckSimple(const double low, ///< lower bound of distribution
//                         const double high ///< upper bound of distribution
// ) {
//   // Init Values Used in Inequality of Interest
//   double val1 = (2 * sqrt(exp(1.0))) / (low + sqrt(pow(low, 2.0) + 4));
//   double val2 = exp((pow(low, 2.0) - low * sqrt(pow(low, 2.0) + 4)) / (4)) ;
//   //
//   
//   // Test if Simple is Preferred
//   if (high > low + val1 * val2) {
//     return true ;
//   } else {
//     return false ;
//   }
// }
// 
// /// Draw using algorithm 1.
// 
// /// 
// /// Naive Accept-Reject algorithm.
// /// 
// inline double UseAlg1(const double low, ///< lower bound of distribution
//                       const double high ///< upper bound of distribution
// ) {
//   // Init Valid Flag
//   int valid = 0 ;
//   //
//   
//   // Init Draw Storage
//   double z = 0.0 ;
//   //
//   
//   // Loop Until Valid Draw
//   while (valid == 0) {
//     z = Rf_rnorm(0.0, 1.0) ;
//     
//     if (z <= high && z >= low) {
//       valid = 1 ;
//     }
//   }
//   //
//   
//   // Returns
//   return z ;
//   //
// }
// 
// /// Draw using algorithm 2.
// 
// /// 
// ///  Accept-Reject Algorithm
// ///
// 
// inline double UseAlg2(const double low ///< lower bound of distribution
// ) {
//   // Init Values
//   const double alphastar = (low +
//                             sqrt(pow(low, 2.0) + 4.0)
//   ) / (2.0) ;
//   const double alpha = alphastar ;
//   double e = 0 ;
//   double z = 0 ;
//   double rho = 0 ;
//   double u = 0 ;
//   //
//   
//   // Init Valid Flag
//   int valid = 0 ;
//   //
//   
//   // Loop Until Valid Draw
//   while (valid == 0) {
//     e = Rf_rexp(1.0) ;
//     z = low + e / alpha ;
//     
//     rho = exp(-pow(alpha - z, 2.0) / 2) ;
//     u = Rf_runif(0, 1) ;
//     if (u <= rho) {
//       // Keep Successes
//       valid = 1 ;
//     }
//   }
//   //
//   
//   // Returns
//   return z ;
//   //
// }
// 
// /// Draw using algorithm 3.
// 
// /// 
// /// Accept-Reject Algorithm
// /// 
// 
// inline double UseAlg3(const double low, ///< lower bound of distribution
//                       const double high ///< upper bound of distribution
// ) {
//   // Init Valid Flag
//   int valid = 0 ;
//   //
//   
//   // Declare Qtys
//   double rho = 0 ;
//   double z = 0 ;
//   double u = 0 ;
//   //
//   
//   // Loop Until Valid Draw
//   while (valid == 0) {
//     z = Rf_runif(low, high) ;
//     if (0 < low) {
//       rho = exp((pow(low, 2.0) - pow(z, 2.0)) / 2) ;
//     } else if (high < 0) {
//       rho = exp((pow(high, 2.0) - pow(z, 2.0)) / 2) ;
//     } else if (0 < high && low < 0) {
//       rho = exp(- pow(z, 2.0) / 2) ;
//     }
//     
//     u = Rf_runif(0, 1) ;
//     if (u <= rho) {
//       valid = 1 ;
//     }
//   }
//   //
//   
//   // Returns
//   return z ;
//   //
// }


/// Draw from an arbitrary truncated normal distribution.

///
/// See Robert (1995): <br />
/// Reference Type: Journal Article <br />
/// Author: Robert, Christian P. <br />
/// Primary Title: Simulation of truncated normal variables <br />
/// Journal Name: Statistics and Computing <br />
/// Cover Date: 1995-06-01 <br />
/// Publisher: Springer Netherlands <br />
/// Issn: 0960-3174 <br />
/// Subject: Mathematics and Statistics <br />
// Start Page: 121 <br />
// End Page: 125 <br />
/// Volume: 5 <br />
/// Issue: 2 <br />
/// Url: http://dx.doi.org/10.1007/BF00143942 <br />
/// Doi: 10.1007/BF00143942 <br />
///

// double rtn1(const double mean,
//             const double sd,
//             const double low,
//             const double high
// ) {
//   // Namespace
//   using namespace Rcpp ;
//   //
//   
//   // Init Useful Values
//   double draw = 0;
//   int type = 0 ;
//   int valid = 0 ; // used only when switching to a simplified version
//   // of Alg 2 within Type 4 instead of the less
//   // efficient Alg 3
//   //
//   
//   // Set Current Distributional Parameters
//   const double c_mean = mean ;
//   double c_sd = sd ;
//   const double c_low = low ;
//   const double c_high = high ;
//   double c_stdlow = (c_low - c_mean) / c_sd ;
//   double c_stdhigh = (c_high - c_mean) / c_sd ; // bounds are standardized
//   //
//   
//   // Map Conceptual Cases to Algorithm Cases
//   // Case 1 (Simple Deterministic AR)
//   // mu \in [low, high]
//   if (0 <= c_stdhigh &&
//       0 >= c_stdlow
//   ) {
//     type = 1 ;
//   }
//   
//   // Case 2 (Robert 2009 AR)
//   // mu < low, high = Inf
//   if (0 < c_stdlow &&
//       c_stdhigh == INFINITY
//   ) {
//     type = 2 ;
//   }
//   
//   // Case 3 (Robert 2009 AR)
//   // high < mu, low = -Inf
//   if (0 > c_stdhigh &&
//       c_stdlow == -INFINITY
//   ) {
//     type = 3 ;
//   }
//   
//   // Case 4 (Robert 2009 AR)
//   // mu -\in [low, high] & (abs(low) =\= Inf =\= high)
//   if ((0 > c_stdhigh || 0 < c_stdlow) &&
//       !(c_stdhigh == INFINITY || c_stdlow == -INFINITY)
//   ) {
//     type = 4 ;
//   }
//   
//   ////////////
//   // Type 1 //
//   ////////////
//   if (type == 1) {
//     draw = UseAlg1(c_stdlow, c_stdhigh) ;
//   }
//   
//   ////////////
//   // Type 3 //
//   ////////////
//   if (type == 3) {
//     c_stdlow = -1 * c_stdhigh ;
//     c_stdhigh = INFINITY ;
//     c_sd = -1 * c_sd ; // hack to get two negative signs to cancel out
//     
//     // Use Algorithm #2 Post-Adjustments
//     type = 2 ;
//   }
//   
//   ////////////
//   // Type 2 //
//   ////////////
//   if (type == 2) {
//     draw = UseAlg2(c_stdlow) ;
//   }
//   
//   ////////////
//   // Type 4 //
//   ////////////
//   if (type == 4) {
//     if (CheckSimple(c_stdlow, c_stdhigh)) {
//       while (valid == 0) {
//         draw = UseAlg2(c_stdlow) ;
//         // use the simple
//         // algorithm if it is more
//         // efficient
//         if (draw <= c_stdhigh) {
//           valid = 1 ;
//         }
//       }
//     } else {
//       draw = UseAlg3(c_stdlow, c_stdhigh) ; // use the complex
//       // algorithm if the simple
//       // is less efficient
//     }
//   }
//   
//   
//   
//   // Returns
//   return  c_mean + c_sd * draw ;
//   //
// }
// 
// // Sample from truncated (from below) normal distribution (at a)
// double rtnorm( double mu, double sigma, double a ){
//   double test_out = help::rtn1(mu, sigma, a, R_PosInf);
//   return test_out;
// }

////////////////////////////////////////////
// Helper functions created for MCMC   /////
////////////////////////////////////////////

// Sample from an integer vector 
int sample_cpp( arma::vec x ){
  // Calling sample()
  Function f( "sample" );
  IntegerVector sampled = f( x, Named( "size" ) = 1 );
  return sampled[ 0 ];
}

// Sample according to a vector of probability
int sample_prob_cpp( arma::vec x, arma::vec prob ){
  // Calling sample()
  Function f( "sample" );
  IntegerVector sampled = f( x, Named( "size" ) = 1, Named( "prob" ) = prob );
  return sampled[ 0 ];
}

// Function :: sampler for inverse-gamma distribution with shape (a) and scale (b) parameters
double rinvgamma_cpp( double a, double b ){
  double invg_return = 1 / Rcpp::rgamma( 1, a, 1/b )[0];
  return invg_return;
}

// Function :: Calculate normal log-density ( univariate )
double log_normal_cpp( double value, double mu, double sigma2 ){
  double log_normal = -0.50*log( 2*M_PI*sigma2 ) - 1/( 2*sigma2 )*pow( value - mu ,2 );
  
  // Return output
  return log_normal; 
}


// Function :: Calculate Bernoulli log-density ( univariate ) ( logit(P(value=1)) = exp(psi)/(1+exp(psi)) )
double log_bern_cpp( double value, double psi ){
  double log_bern = value*psi - log(1+exp(psi)); 
  
  // Return output
  return log_bern; 
}


// Function :: Calculate Beta log-density ( univariate )
double log_beta_cpp( double value, double a, double b ){
  double log_beta = (a-1) * log(value) + (b-1) * log(1-value) - lgamma(a) - lgamma(b) + lgamma(a+b);
  
  // Return output
  return log_beta; 
}

// Function :: Calculate multivariate normal log-density only with diagonal covariance matrix
double log_mvnorm_cpp( arma::vec value, arma::vec mu, arma::vec diag_sig2 ){
  int n = value.size();
  double d = 0;
  
  for (int i = 0; i < n; ++i){
    d += help::log_normal_cpp( value[i], mu[i], diag_sig2[i] );
  }
  
  return d;
}

// Function :: calculate (standardized) CDF for normal distribution for matrix
arma::mat pnorm_cpp(
    arma::mat Aprime
){
  int N = Aprime.n_rows;
  int P = Aprime.n_cols;
  arma::mat out( N, P, fill::zeros);
  double tmp = 0;
  
  for(int i = 0; i < N; i++ ){
    for( int j = 0; j < P; j++){
      tmp = Aprime( i, j );
      out( i, j ) = R::pnorm( tmp, 0.0, 1.0, true, false );
    }
  }
  
  return out;
}

// Function :: convert partition C1 to canonical form (position that this number first appears)
arma::vec canonical_form(arma::vec C1){
  int l = C1.size();
  arma::vec C2(l, fill::zeros);
  for(int i = 0; i < l; i++){
    uvec pos = find(C1==C1[i]);
    C2[i] = pos[0];
  }
  return C2;
}

// Function :: decide whether two partitions are the same
bool sameC(
    arma::vec C1,
    arma::vec C2
){
  // the following partitions are the same: (1,1,1,2,2,3,3,3) and (2,2,2,3,3,1,1,1)
  
  C1 = help::canonical_form(C1);
  C2 = help::canonical_form(C2);
  int l = C1.size();
  for(int i = 0; i < l; i++){
    if(C1[i] != C2[i]){
      return FALSE;
    }
  }
  
  return TRUE;
}

// Function :: return the clustering indicator for the reduced partition C_j^{R_j}
arma::vec getCjRj(
    int i,
    int Ctype, // 1 = add i, -1 = minus i, other = without adding/erasing i
    arma::vec gammaj,
    arma::vec Cj
){
  
  if( Ctype == 1 ){
    gammaj[ i ] = 1;
  }else if( Ctype == -1 ){
    gammaj[ i ] = 0;
  }
  
  uvec Rj = find( gammaj == 1 ); // The vector of all subjects in cluster k
  
  arma::vec CjRj = Cj( Rj );
  return CjRj;
}

// Function :: decide whether C1 and C2 are compatible given gammaj
bool compatible(
    int i,
    int Ctype, // can force i to be added/erased
    arma::vec C1,
    arma::vec C2,
    arma::vec gammaj
){
  arma::vec C1Rj = help::getCjRj(i, Ctype, gammaj, C1);
  arma::vec C2Rj = help::getCjRj(i, Ctype, gammaj, C2);
  bool comp = help::sameC(C1Rj, C2Rj);
  return comp;
}

// Function :: calculate prob of c_ij = k given some existing partition C0
double probcij( 
    int k,
    arma::vec C0,
    double a
){
  int N = C0.size(); // total number of units
  if( N == 0 ){
    return 1;
  }
  arma::vec uniqC0 = unique( C0 ); // unique clusters
  int kj = uniqC0.size();
  double denom = N + a;
  
  double p = a / denom; // if k doesn't equal to any existing labels in C0
  uvec tmp;
  double tsize = 0;
  
  for( int ktmp = 0; ktmp < kj; ktmp++ ){
    if( uniqC0[ ktmp ] == k ){
      tmp = find( C0 == uniqC0[ ktmp ] );
      tsize = tmp.size();
      p = tsize / denom;
    }
  }
  
  return p;
}

// return a vector v with element i removed
arma::vec remove_i(arma::vec v, int i){
  int n = v.size();
  arma::vec out( n-1, fill::zeros );
  
  if( i != 0 ){
    out.subvec( 0, i-1 ) = v.subvec( 0, i-1 );
  }
  
  if( i != (n-1) ){
    out.subvec( i, n-2 ) = v.subvec( i+1, n-1 );
  }
  
  return out;
}

// return a matrix m with row i removed
arma::mat remove_rowi(arma::mat m, int i){
  int n = m.n_rows;
  int p = m.n_cols;
  arma::mat out( n-1, p, fill::zeros );
  
  if( i != 0 ){
    out.rows( 0, i-1 ) = m.rows( 0, i-1 );
  }
  
  if( i != (n-1) ){
    out.rows( i, n-2 ) = m.rows( i+1, n-1 );
  }
  
  return out;
}

////////////////////////////////////////////////
//              Main functions             /////
////////////////////////////////////////////////

List update_omega(
  List T,
  List Z,
  arma::cube Beta,
  arma::mat Theta,
  arma::mat IDstart,
  arma::mat IDend
){
  
  int J = Beta.n_slices;
  List updated_omega(J);
  
  for(int j = 0; j < J; j++){
    arma::vec idS = IDstart.col(j);
    arma::vec idE = IDend.col(j);
    
    arma::mat Tj = T[j];
    arma::mat Zj = Z[j];
    
    arma::vec Thetaj = Theta.col( j );
    
    int nobs = Tj.n_rows;
    arma::vec omegaj( nobs );
    
    arma::mat betaj = Beta.slice(j);
    int N = betaj.n_rows;
    
    int pos = 0;
    
    for(int i = 0; i < N; i++){
      for(int k = idS[i]; k <= idE[i]; k++){
        arma::mat psi(1,1);
        psi = Tj.row( pos ) * betaj.row( i ).t() + Zj.row( pos ) * Thetaj;
        omegaj[ pos ] = help::samplepg( psi[ 0 ] );
        
        pos += 1;
      }
    } // end for i
    
    updated_omega[j] = omegaj;
    
  } // end for j
  
  return updated_omega;
}

// Function to update each gamma element as described in section 1.1
arma::mat update_gamma_ij(
    int i,
    int j,
    arma::cube X,
    arma::mat Gamma,
    arma::mat C,
    arma::mat Eta,
    double a_loc
){
  arma::mat Alpha, Xj, XEta;
  arma::vec gammaj, Cj, CjRjmi, Cjm1;
  double phi_ij, post_cij, post_gammaij;
  int gnew;
  bool comp;
  
  if( j > 0 ){ // only update gamma for j>0
    
    // Check compatibility
    gammaj = Gamma.col( j );
    Cjm1 = C.col( j-1 );
    Cj = C.col( j );
    comp = help::compatible( i, 1, Cjm1, Cj, gammaj); // compatibility of fixed subjects between j-1 and j with i added
    
    if(!comp){
      
      gnew = 0;
      
    }else{
      
      Xj = X.slice(j);
      XEta = Xj.row(i) * Eta.col(j);
      phi_ij = exp(XEta[0]) / ( 1 + exp(XEta[0]) );
      CjRjmi = help::getCjRj(i, -1, gammaj, Cj); // get the clustering indicator for the reduced partition C_j^{R_j (-i)}
      post_cij = help::probcij(Cj[i], CjRjmi, a_loc); // p( cij = k | rho_j^Rj (-i) )
      post_gammaij = phi_ij / ( phi_ij + (1-phi_ij) * post_cij ); // p( gammaij = 1 | - )
      gnew = rbinom( 1, 1, post_gammaij )[0];
      
    }
    
    Gamma( i, j ) = gnew;
    
  } // end if j>0
  
  return Gamma;
}

arma::cube update_beta_jk_PG(
    int j, // time period
    int k, // the cluster indicator
    List Y,
    List T,
    List omega,
    List Z,
    arma::mat C,
    arma::cube Beta,
    arma::mat Theta,
    arma::cube Lambda2,
    arma::mat Tau2,
    arma::mat IDstart,
    arma::mat IDend
){
  
  // Create some temporary variables for the polya-gamma sampler

  arma::mat Zj = Z[j];
  arma::vec muij = Zj * Theta.col(j);
    
  arma::vec Yj = Y[j];
  arma::vec Kj = Yj - 0.5;
  
  arma::vec omegaj = omega[j];
  arma::vec zbarj = Kj/omegaj - muij;
  
  arma::mat Tj = T[j];
  
  // Create some temporary variables for the prior
  
  arma::mat Lambda2j = Lambda2.slice(j);
  arma::vec Tau2j = Tau2.col(j);
  
  // Find subjects and the number of subjects belong to kth cluster
  
  arma::vec Cj = C.col(j);
  arma::vec idS = IDstart.col(j);
  arma::vec idE = IDend.col(j);
  
  uvec pos = find( Cj == k ); // The vector of all subjects in cluster k
  int Cjk_size = pos.size();
  int i = pos[ 0 ]; // one of the subjects in cluster k
  
  // Make the prior variance matrix
  
  int D = Lambda2j.n_cols;
  arma::vec s2_diag( D );
  for(int d = 0; d < D; d++){
    s2_diag[ d ] = Lambda2j( i, d ) * Tau2j[ i ];
  }
  arma::mat Sigma_star_inv = diagmat( 1 / s2_diag );
  arma::mat Sigma_bjk_star = Sigma_star_inv; // Starting point of calculating Sigma_{beta^jk}*^-1
  arma::vec mu_bjk_star( D, fill::zeros ); // Starting point of calculating mu_{beta^jk}*
  
  arma::mat Tij, omega_mat_ij;
  arma::vec Zij, omegaij;
  
  for(int idx = 0; idx < Cjk_size; ++idx){
    
    i = pos[ idx ]; // subject id
    int is = idS[i]; // Y_ij = Yj[is:ie]
    int ie = idE[i];

    Tij = Tj.rows( is, ie );
    Zij = zbarj.subvec( is, ie );
    omegaij = omegaj.subvec( is, ie );
    omega_mat_ij = diagmat( omegaij );
    
    Sigma_bjk_star += Tij.t() * omega_mat_ij * Tij;
    mu_bjk_star += Tij.t() * omega_mat_ij * Zij;
  }
  
  Sigma_bjk_star = inv( Sigma_bjk_star );
  mu_bjk_star = Sigma_bjk_star * mu_bjk_star;
  
  arma::mat beta_new = help::mvrnormArma( 1, mu_bjk_star, Sigma_bjk_star );
  arma::mat beta_temp = Beta.slice( j );
  
  // put updated beta back into every member within the same cluster
  for(int idx = 0; idx < Cjk_size; ++idx){
    i = pos[ idx ]; // subject id
    beta_temp.row( i ) = beta_new;
  }
  
  Beta.slice( j ) = beta_temp;
  
  return Beta;
}

arma::mat update_theta_PG(
    List Y, 
    List T, 
    List Z,
    List omega1, 
    arma::cube Beta,
    arma::mat IDstart, 
    arma::mat IDend, 
    arma::mat Theta,
    arma::mat SigTheta
){
  
  
  // int dz = SigTheta.n_cols;
  int J = Beta.n_slices;
  int N = Beta.n_rows;
  
  // arma::mat Obarj, Tj, invSigTheta, VTheta_inv, VTheta, Betaj;
  // arma::vec MTheta, Yj, Kj, omega1j, Vj, ztildej, phiij;
  // arma::vec kappa( N, fill::zeros );
  // arma::mat Xj( N, dx, fill::zeros);
  
  arma::mat Betaj, Tij, Obarj, invSigTheta, VTheta_inv, VTheta;
  arma::vec Kj, Ztildej, idS, idE, phiij, MTheta;
  int is, ie;
  
  for(int j = 0; j < J; j++){
    
    // Create some temporary variables for the polya-gamma sampler
    arma::vec Yj = Y[j];
    Kj = Yj - 0.5;
    
    // int nobs = Yj.size();
    
    arma::mat Tj = T[j];
    arma::mat Zprimej = Z[j]; // nobs*dz matrix
    
    arma::vec omega1j = omega1[j];
    Ztildej = Kj/omega1j;
    
    Betaj = Beta.slice( j );
    
    idS = IDstart.col(j);
    idE = IDend.col(j);
    
    for(int i = 0; i < N; i++){
      is = idS[i]; // Y_ij = Yj[is:ie]
      ie = idE[i];
      
      Tij = Tj.rows( is, ie );
      phiij = Tij * Betaj.row(i).t();
      Ztildej.subvec( is, ie ) -= phiij;
    }
    
    Obarj = diagmat( omega1j );

    invSigTheta = inv(SigTheta);
    VTheta_inv = invSigTheta;
    VTheta_inv += Zprimej.t() * Obarj * Zprimej;

    VTheta = inv(VTheta_inv);
    MTheta = VTheta * ( Zprimej.t() * Obarj * Ztildej );
    // MTheta = VTheta * ( Zprimej.t() * Obarj * Ztildej + invSigTheta * MuTheta );

    Theta.col(j) = help::mvrnormArma( 1, MTheta, VTheta).t();

  }
  
  return Theta;
}



arma::cube update_lambda2_jk( 
    int j, // time period
    int k, // the cluster indicator
    arma::mat C,
    arma::cube Beta,
    arma::cube Lambda2,
    arma::mat Tau2,
    arma::cube Nulam
){
  
  arma::mat Betaj = Beta.slice(j);
  arma::mat lambda2_temp = Lambda2.slice(j);
  arma::vec Tau2j = Tau2.col(j);
  arma::mat Nulamj = Nulam.slice(j);
  
  arma::vec Cj = C.col(j);
  
  // Find subjects and the number of subjects belong to kth cluster
  uvec pos = find( Cj == k ); // The vector of all subjects in cluster k
  int Cjk_size = pos.size();
  int i = pos[ 0 ]; // one of the subjects in cluster k
  int D = Betaj.n_cols;
  
  double beta_jkd = 0;
  double tau2_jk = 0;
  double nulam_jkd = 0;
  double bterm = 0;
  double lambda2_new = 0;
  
  for(int d = 0; d < D; d++){
    beta_jkd = Betaj( i, d );
    tau2_jk = Tau2j[i];
    nulam_jkd = Nulamj( i, d );
    bterm = beta_jkd * beta_jkd / tau2_jk / 2 + 1 / nulam_jkd;
    lambda2_new = help::rinvgamma_cpp(1, bterm);
    
    for(int idx = 0; idx < Cjk_size; ++idx){
      i = pos[ idx ]; // subject id
      lambda2_temp( i, d ) = lambda2_new;
    }
  }
  
  Lambda2.slice( j ) = lambda2_temp;
  
  return Lambda2;
}

arma::mat update_tau2_jk( 
    int j, // time period
    int k, // the cluster indicator
    arma::mat C,
    arma::cube Beta,
    arma::cube Lambda2,
    arma::mat Tau2,
    arma::mat Nutau
){
  
  arma::mat Betaj = Beta.slice(j);
  arma::mat Lambda2j = Lambda2.slice(j);
  arma::vec tau2_temp = Tau2.col(j);
  arma::vec Nutauj = Nutau.col(j);
  
  arma::vec Cj = C.col(j);
  
  // Find subjects and the number of subjects belong to kth cluster
  uvec pos = find( Cj == k ); // The vector of all subjects in cluster k
  int Cjk_size = pos.size();
  int i = pos[ 0 ]; // one of the subjects in cluster k
  int D = Betaj.n_cols;
  
  double beta_jkd = 0;
  double lambda2_jkd = 0;
  
  double aterm = 0.5 * (D+1);
  double bterm = 1/Nutauj[i]; // start with 1/nu_tau_jk*
  
  for(int d = 0; d < D; d++){
    beta_jkd = Betaj( i, d );
    lambda2_jkd = Lambda2j( i, d );
    bterm += beta_jkd * beta_jkd / lambda2_jkd / 2;
  }
  
  double tau2_new = help::rinvgamma_cpp(aterm, bterm);
  
  for(int idx = 0; idx < Cjk_size; ++idx){
    i = pos[ idx ]; // subject id
    tau2_temp[ i ] = tau2_new;
  }
  
  Tau2.col( j ) = tau2_temp;
  
  return Tau2;
}


arma::cube update_nu_lambda_jk( 
    int j, // time period
    int k, // the cluster indicator
    arma::mat C,
    arma::cube Lambda2,
    arma::cube Nulam
){
  
  arma::mat Lambda2j = Lambda2.slice(j);
  arma::mat nulam_temp = Nulam.slice(j);
  
  arma::vec Cj = C.col(j);
  
  // Find subjects and the number of subjects belong to kth cluster
  uvec pos = find( Cj == k ); // The vector of all subjects in cluster k
  int Cjk_size = pos.size();
  int i = pos[ 0 ]; // one of the subjects in cluster k
  int D = Lambda2.n_cols;
  
  double lambda2_jkd = 0;
  double bterm = 0;
  double nulam_new = 0;
  
  for(int d = 0; d < D; d++){
    lambda2_jkd = Lambda2j( i, d );
    bterm = 1 + 1 / lambda2_jkd;
    nulam_new = help::rinvgamma_cpp(1, bterm);
    
    for(int idx = 0; idx < Cjk_size; ++idx){
      i = pos[ idx ]; // subject id
      nulam_temp( i, d ) = nulam_new;
    }
  }
  
  Nulam.slice( j ) = nulam_temp;
  
  return Nulam;
}

arma::mat update_nu_tau_jk( 
    int j, // time period
    int k, // the cluster indicator
    arma::mat C,
    arma::mat Tau2,
    arma::mat Nutau 
){
  
  arma::vec Tau2j = Tau2.col(j);
  arma::vec nutau_temp = Nutau.col(j);
  
  arma::vec Cj = C.col(j);
  
  // Find subjects and the number of subjects belong to kth cluster
  uvec pos = find( Cj == k ); // The vector of all subjects in cluster k
  int Cjk_size = pos.size();
  int i = pos[ 0 ]; // one of the subjects in cluster k
  
  double nutau_new = help::rinvgamma_cpp(1, Tau2j[ i ]);
  
  for(int idx = 0; idx < Cjk_size; ++idx){
    i = pos[ idx ]; // subject id
    nutau_temp[ i ] = nutau_new;
  }
  
  Nutau.col( j ) = nutau_temp;
  
  return Nutau;
}



// List dp_update_C_ij(
//     int i,
//     int j,
//     List Y,
//     List T,
//     arma::mat C,
//     arma::cube Beta,
//     arma::cube Lambda2,
//     arma::mat Tau2,
//     arma::mat IDstart,
//     arma::mat IDend,
//     double a
// ){
// 
//   arma::vec Yj = Y[j];
//   arma::vec Yij = Yj.subvec( IDstart(i,j), IDend(i,j) );
//   arma::mat Tj = T[j];
//   arma::mat Tij = Tj.rows( IDstart(i,j), IDend(i,j) );
//   int N = Yj.size();
//   int ni = Yij.size();
//   
//   arma::mat Betaj = Beta.slice(j);
//   arma::mat Betajmi = help::remove_rowi(Betaj, i);
//   arma::mat Lambda2j = Lambda2.slice(j);
//   int D = Lambda2j.n_cols;
//   arma::mat Lambda2jmi = help::remove_rowi(Lambda2j, i);
//   arma::vec Tau2j = Tau2.col(j);
//   arma::vec Tau2jmi = help::remove_i(Tau2j, i);
//   
//   arma::vec Cj = C.col(j);
//   arma::vec Ctemp = Cj; // a temporary vector to hold cluster indicators to check compatibility after change
//   arma::vec Cjmi = help::remove_i(Cj, i);
//   
//   arma::vec uniqCjmi = unique( Cjmi ); // unique cluster indicators in S after remove unit i
//   int kjmi = uniqCjmi.size(); // number of unique cluster indicators in S after remove unit i
//   
//   arma::mat uniqBeta( kjmi + 1, D, fill::zeros); // matrix of unique beta for each cluster Sjk^{-i}
//   arma::mat uniqLambda2(kjmi + 1, D, fill::zeros);
//   arma::vec uniqTau2(kjmi+1, fill::zeros);
//   arma::vec postProb_cij_vec( kjmi+1, fill::zeros ); // probabilities of moving to each cluster
//   arma::vec candidate_cij( kjmi+1, fill::zeros ); // candidates for new cluster indicators
//   arma::vec idvec = linspace(0, kjmi, kjmi+1); // a sequence from 0 to kjmi
//   
//   
//   candidate_cij.subvec( 0, kjmi - 1 ) = uniqCjmi;
//   double postProb_cij = 0;
//   
//   arma::vec diagtmp(D, fill::zeros); // a vector to hold the diagnoal elements of Sigma for new cluster
//   arma::mat newSigma; // a matrix to hold Sigma for new cluster
//   arma::vec newmu(D, fill::zeros); // a vector to hold mu for new cluster
//   
//   uvec Sjk;
//   int Sjk0, Sjk_size, new_pos;
//   arma::mat meantmp; // a place to hold T*beta
// 
//     
//   // prob of assigning to one of the existing clusters
//   for(int k = 0; k < kjmi; k++){
//     
//     // Find subjects and the number of subjects belong to kth cluster
//     Sjk = find( Cjmi == uniqCjmi[ k ] ); // The vector of all subjects in cluster S_jk
//     Sjk_size = Sjk.size();
//     Sjk0 = Sjk[ 0 ]; // pick one of the subjects in cluster k
//     
//     uniqBeta.row( k ) = Betajmi.row( Sjk0 ); //
//     uniqTau2[ k ] = Tau2jmi[ Sjk0 ];
//     uniqLambda2.row( k ) = Lambda2jmi.row( Sjk0 );
//     
//     postProb_cij = 0;
//     for(int t = 0; t < ni; t++ ){
//       meantmp = Tij.row(t) * uniqBeta.row( k ).t();
//       postProb_cij += help::log_bern_cpp( Yij[t], meantmp[0] );
//     }
//     
//     postProb_cij += log( Sjk_size / (N + a - 1) );
//     postProb_cij_vec[k] = postProb_cij;
//     
//   } // end for k
//   
//   // prob of assigning to a new cluster
//   candidate_cij[ kjmi ] = max( uniqCjmi ) + 1; // new cluster label
//   
//   uniqTau2[ kjmi ] = Rcpp::rgamma( 1, 2, 1/2.0 )[0]; // shape, scale parameters
//   
//   for(int d = 0; d < D; d++){
//     uniqLambda2( kjmi, d ) = Rcpp::rgamma( 1, 2, 1/2.0 )[0];
//     diagtmp[ d ] = uniqLambda2( kjmi, d ) * uniqTau2[ kjmi ];
//   }
//   
//   newSigma = diagmat( diagtmp );
//   uniqBeta.row( kjmi ) = help::mvrnormArma(1, newmu, newSigma);
//   
//   postProb_cij = 0;
//   for(int t = 0; t < ni; t++ ){
//     meantmp = Tij.row(t) * uniqBeta.row( kjmi ).t();
//     postProb_cij += help::log_bern_cpp( Yij[t], meantmp[0] );
//   }
//   postProb_cij += log( a / (N + a - 1) );
//   postProb_cij_vec[ kjmi ] = postProb_cij;
//   
//   postProb_cij_vec = exp( postProb_cij_vec - help::logsumexp( postProb_cij_vec ) ); // convert log-probabilities to probabilities
//   
//   new_pos = help::sample_prob_cpp( idvec, postProb_cij_vec );
//   C( i, j ) = candidate_cij[ new_pos ];
//   Betaj.row( i ) = uniqBeta.row( new_pos );
//   Lambda2j.row( i ) = uniqLambda2.row( new_pos );
//   Tau2( i, j ) = uniqTau2[ new_pos ];
//   
//   Beta.slice( j ) = Betaj;
//   Lambda2.slice( j ) = Lambda2j;
//   
//   List output( 4 );
//   output[ 0 ] = C;
//   output[ 1 ] = Beta;
//   output[ 2 ] = Lambda2;
//   output[ 3 ] = Tau2;
//   
//   return output;
// }


List update_C_ij(
    bool DFP,
    int i,
    int j,
    List Y,
    List T,
    List Z,
    arma::mat Gamma,
    arma::mat C,
    arma::cube Beta,
    arma::mat Theta,
    arma::cube Lambda2,
    arma::mat Tau2,
    arma::mat IDstart,
    arma::mat IDend,
    double a
){
  int N = IDstart.n_rows;
  int J = C.n_cols;
  arma::vec Yj = Y[j];
  arma::vec Yij = Yj.subvec( IDstart(i,j), IDend(i,j) );
  arma::mat Tj = T[j];
  arma::mat Tij = Tj.rows( IDstart(i,j), IDend(i,j) );
  arma::mat Zj = Z[j];
  arma::mat Zij = Zj.rows( IDstart(i,j), IDend(i,j) );
  
  int ni = Yij.size();
  
  arma::mat Betaj = Beta.slice(j);
  arma::mat Betajmi = help::remove_rowi(Betaj, i);
  arma::mat Lambda2j = Lambda2.slice(j);
  int D = Lambda2j.n_cols;
  arma::mat Lambda2jmi = help::remove_rowi(Lambda2j, i);
  arma::vec Tau2j = Tau2.col(j);
  arma::vec Tau2jmi = help::remove_i(Tau2j, i);

  arma::vec Cj = C.col(j);
  arma::vec Ctemp = Cj; // a temporary vector to hold cluster indicators to check compatibility after change
  arma::vec Cjmi = help::remove_i(Cj, i);

  arma::vec uniqCjmi = unique( Cjmi ); // unique cluster indicators in S after remove unit i
  int kjmi = uniqCjmi.size(); // number of unique cluster indicators in S after remove unit i

  arma::mat uniqBeta( kjmi + 1, D, fill::zeros); // matrix of unique beta for each cluster Sjk^{-i}
  arma::mat uniqLambda2(kjmi + 1, D, fill::zeros);
  arma::vec uniqTau2(kjmi+1, fill::zeros);
  arma::vec postProb_cij_vec( kjmi+1, fill::zeros ); // probabilities of moving to each cluster
  arma::vec candidate_cij( kjmi+1, fill::zeros ); // candidates for new cluster indicators
  arma::vec idvec = linspace(0, kjmi, kjmi+1); // a sequence from 0 to kjmi


  candidate_cij.subvec( 0, kjmi - 1 ) = uniqCjmi;
  double postProb_cij = 0;

  arma::vec diagtmp(D, fill::zeros); // a vector to hold the diagnoal elements of Sigma for new cluster
  arma::mat newSigma; // a matrix to hold Sigma for new cluster
  arma::vec newmu(D, fill::zeros); // a vector to hold mu for new cluster

  uvec Sjk;
  int Sjk0, Sjk_size, new_pos;
  arma::mat meantmp; // a place to hold T*beta
  bool comp;
  
  
  if( (DFP && (Gamma(i,j) == 0)) || !DFP ){

    // prob of assigning to one of the existing clusters
    for(int k = 0; k < kjmi; k++){

      // Find subjects and the number of subjects belong to kth cluster
      Sjk = find( Cjmi == uniqCjmi[ k ] ); // The vector of all subjects in cluster S_jk
      Sjk_size = Sjk.size();
      Sjk0 = Sjk[ 0 ]; // pick one of the subjects in cluster k

      uniqBeta.row( k ) = Betajmi.row( Sjk0 ); //
      uniqTau2[ k ] = Tau2jmi[ Sjk0 ];
      uniqLambda2.row( k ) = Lambda2jmi.row( Sjk0 );
      
      postProb_cij = 0;
      for(int t = 0; t < ni; t++ ){
        meantmp = Tij.row(t) * uniqBeta.row( k ).t() + Zij.row(t) * Theta.col(j);
        postProb_cij += help::log_bern_cpp( Yij[t], meantmp[0] );
      }

      postProb_cij += log( Sjk_size / (N + a - 1) );
      postProb_cij_vec[k] = postProb_cij;

    } // end for k

    // prob of assigning to a new cluster
    candidate_cij[ kjmi ] = max( uniqCjmi ) + 1; // new cluster label

    uniqTau2[ kjmi ] = Rcpp::rgamma( 1, 2, 1/2.0 )[0]; // shape, scale parameters
    // Rcout << "uniqTau2: " << uniqTau2 << endl;
    
    for(int d = 0; d < D; d++){
      uniqLambda2( kjmi, d ) = Rcpp::rgamma( 1, 2, 1/2.0 )[0];
      diagtmp[ d ] = uniqLambda2( kjmi, d ) * uniqTau2[ kjmi ];
    }
    // Rcout << "diagtmp: " << diagtmp << endl;
    
    
    newSigma = diagmat( diagtmp );
    // Rcout << "newSigma: " << newSigma << endl;
    
    uniqBeta.row( kjmi ) = help::mvrnormArma(1, newmu, newSigma);

    postProb_cij = 0;
    for(int t = 0; t < ni; t++ ){
      meantmp = Tij.row(t) * uniqBeta.row( kjmi ).t() + Zij.row(t) * Theta.col(j);
      postProb_cij += help::log_bern_cpp( Yij[t], meantmp[0] );
    }
    postProb_cij += log( a / (N + a - 1) );
    postProb_cij_vec[ kjmi ] = postProb_cij;

    postProb_cij_vec = exp( postProb_cij_vec - help::logsumexp( postProb_cij_vec ) ); // convert log-probabilities to probabilities

    if( DFP && (j < (J-1)) ){
      if( Gamma( i, j + 1 ) == 1 ){ // check compatibility
        for(int k = 0; k < (kjmi+1); k++){
          Ctemp[ i ] = candidate_cij[ k ]; // propose a temp partition
          comp = help::compatible( 0, 0, Ctemp, C.col(j + 1), Gamma.col( j + 1 ) );
          // set first 2 arguments to 0 to check compatibility of Ctemp and C only, without touching unit i
          if( !comp ){
            postProb_cij_vec[ k ] = 0;
          }
        } // end for
      } // end if gamma
    } // end if j < J-1

    new_pos = help::sample_prob_cpp( idvec, postProb_cij_vec );
    C( i, j ) = candidate_cij[ new_pos ];
    Betaj.row( i ) = uniqBeta.row( new_pos );
    Lambda2j.row( i ) = uniqLambda2.row( new_pos );
    Tau2( i, j ) = uniqTau2[ new_pos ];

    Beta.slice( j ) = Betaj;
    Lambda2.slice( j ) = Lambda2j;
  } // end if
  
  
  List output( 4 );
  output[ 0 ] = C;
  output[ 1 ] = Beta;
  output[ 2 ] = Lambda2;
  output[ 3 ] = Tau2;
  
  return output;
}



arma::mat update_eta_PG(
    arma::cube X,
    arma::mat Gamma,
    arma::mat Eta,
    arma::vec MuEta,
    arma::mat SigEta,
    arma::mat w
){
  
  int dx = X.n_cols;
  int J = X.n_slices;
  int N = X.n_rows;
  
  arma::mat wj, invSigEta, VEta_inv, VEta;
  arma::vec MEta;
  arma::vec kappa( N, fill::zeros );
  arma::mat Xj( N, dx, fill::zeros);
  
  for(int j = 1; j < J; j++){
    kappa = Gamma.col( j ) - 0.5;
    
    Xj = X.slice( j );
    wj = diagmat( w.col(j) );
    
    invSigEta = inv(SigEta);
    VEta_inv = invSigEta;
    VEta_inv += Xj.t() * wj * Xj; /// 
    
    VEta = inv(VEta_inv);
    MEta = VEta * ( Xj.t() * kappa + invSigEta * MuEta ); 
    
    Eta.col(j) = help::mvrnormArma( 1, MEta, VEta).t();
  }
  
  return Eta;
}


arma::mat update_w( 
    arma::cube X, 
    arma::mat Eta
){
  
  int dx = X.n_cols;
  int J = X.n_slices;
  int N = X.n_rows;
  
  arma::mat Xmat( N, dx, fill::zeros );
  arma::mat updated_w( N, J, fill::zeros );
  arma::vec psi( N, fill::zeros );
  
  for(int j = 1; j < J; j++){
    Xmat = X.slice( j );
    psi = Xmat * Eta.col(j);
    for(int i = 0; i < N; i++ ){
      updated_w( i, j ) = help::samplepg( psi[ i ] );
    }
  }
  return updated_w;
  
}

}  // For namespace 'help'

// Function :: MCMC algorithm
// [[Rcpp::export]]
List DFPmcmc(
    int iterations,
    int thin,
    List Y,     // We use list to store the observation because the number of Ys at each period could be difference. Each slice is a vector of length sum(i,t): [y1,...,yN]
    List T,     // For the same reason, we use a list to store the basis functions T for each period. Each slice is a matrix of size sum(i,t)*D
    List Z,     // For the ease of coding, we require the baseline covariate Z to have the same length as Y, even it's technically only N*dz. This could be done by a wrapper function on the outside. 
    arma::cube X, // N*dx*J cube
    arma::cube Beta, // N*D*J cube
    arma::cube Lambda2, // N*D*J cube
    arma::cube Nulam, // N*D*J cube
    arma::mat Gamma, // N*J matrix
    arma::mat C, // N*J matrix
    arma::mat Tau2, // N*J matrix
    arma::mat Nutau, // N*J matrix
    arma::mat Eta, // dx*J matrix
    arma::mat Theta, // dz*J matrix
    arma::mat IDstart, // N*J
    arma::mat IDend, // N*J
    arma::mat SigEta,
    arma::mat SigTheta,
    arma::vec MuEta,
    double a
){
  int N = Tau2.n_rows;
  int J = Tau2.n_cols;
  int dx = X.n_cols;
  int dz = SigTheta.n_cols;
  
  List omega = Y; // The auxiliary parameters for updating the beta. Each slice j of the list, same shape (sum_i sum_t) and layout as Y_j.
  List omega1 = Y; // The auxiliary parameters for updating the theta. Each slice also accommodates all observations Y_j
  arma::mat w( N, J, fill::zeros ); // The auxiliary parameters for updating Eta.
  
  // Create these to store the output:
  arma::cube GammaCube( N, J, iterations/thin, fill::zeros );
  arma::cube CCube( N, J, iterations/thin, fill::zeros );
  List BetaList( iterations/thin );
  List Lambda2List( iterations/thin );
  arma::cube Tau2Cube( N, J, iterations/thin, fill::zeros );
  List NulamList( iterations/thin );
  arma::cube NutauCube( N, J, iterations/thin, fill::zeros );
  arma::cube EtaCube( dx, J, iterations/thin, fill::zeros );
  arma::cube ThetaCube( dz, J, iterations/thin, fill::zeros );
  arma::vec uniqCj;
  
  List clist( 4 ); // to store update_C_ij outputs
    
  int kj = 0; // number of unique cluster at period j
  
  for( int iter = 0; iter < iterations; iter++ ){
    // Rcout << "iter" << iter << endl;
    
    omega = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    omega1 = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    
    for( int j = 0; j < J; j++ ){
      
      for( int i = 0; i < N; i++ ){
        Gamma = help::update_gamma_ij( i, j, X, Gamma, C, Eta, a );
      }

      for( int i = 0; i < N; i++ ){
        clist = help::update_C_ij( TRUE, i, j, Y, T, Z, Gamma, C, Beta, Theta, Lambda2, Tau2, IDstart, IDend, a);
        C = as<arma::mat>( clist[ 0 ]);
        Beta = as<arma::cube>( clist[ 1 ] );
        Lambda2 = as<arma::cube>( clist[ 2 ] );
        Tau2 = as<arma::mat>( clist[ 3 ] );
      }

      uniqCj = unique( C.col( j ) );
      kj = uniqCj.size(); // size of all clusters: 1+q_p

      for( int k = 0; k < kj; k++ ){

        Beta = help::update_beta_jk_PG( j, uniqCj[ k ], Y, T, omega, Z, C, Beta, Theta, Lambda2, Tau2, IDstart, IDend );
        
        Lambda2 = help::update_lambda2_jk( j, uniqCj[ k ], C, Beta, Lambda2, Tau2, Nulam );
        Tau2 = help::update_tau2_jk( j, uniqCj[ k ], C, Beta, Lambda2, Tau2, Nutau );
        Nulam = help::update_nu_lambda_jk( j, uniqCj[ k ], C, Lambda2, Nulam );
        Nutau = help::update_nu_tau_jk( j, uniqCj[ k ], C, Tau2, Nutau );

      } // end for k
      
    } // end for j

    w = help::update_w(X, Eta);
    Eta = help::update_eta_PG(X, Gamma, Eta, MuEta, SigEta, w);

    Theta = help::update_theta_PG(Y, T, Z, omega1, Beta, IDstart, IDend, Theta, SigTheta );
    
    if( ( iter + 1 ) % thin == 0 ){
      BetaList[ ( iter + 1 )/thin - 1  ] = Beta;
      Lambda2List[ ( iter + 1 )/thin - 1  ] = Lambda2;
      NulamList[ ( iter + 1 )/thin - 1  ] = Nulam;
      GammaCube.slice( ( iter + 1 )/thin - 1  ) = Gamma;
      CCube.slice( ( iter + 1 )/thin - 1  ) = C;
      Tau2Cube.slice( ( iter + 1 )/thin - 1  ) = Tau2;
      NutauCube.slice( ( iter + 1 )/thin - 1  ) = Nutau;
      EtaCube.slice( ( iter + 1 )/thin - 1  ) = Eta;
      ThetaCube.slice( ( iter + 1 )/thin - 1  ) = Theta;
    }
    
    // Print out progress
    double printer = iter % 50;
    if( printer == 0 ){
      Rcout << "Iteration = " << iter << endl;
    }
    
  } // end for iter
  
  List output( 9 );
  output[ 0 ] = GammaCube;
  output[ 1 ] = CCube;
  output[ 2 ] = BetaList;
  output[ 3 ] = Lambda2List;
  output[ 4 ] = Tau2Cube;
  output[ 5 ] = NulamList;
  output[ 6 ] = NutauCube;
  output[ 7 ] = EtaCube;
  output[ 8 ] = ThetaCube;

  return output;
}


// [[Rcpp::export]]
List DPmcmc(
    int iterations,
    int thin,
    List Y,     // We use list to store the observation because the number of Ys at each period could be difference. Each slice is a vector of length sum(i,t): [y1,...,yN]
    List T,     // For the same reason, we use a list to store the basis functions T for each period. Each slice is a matrix of size sum(i,t)*D
    List Z,     // For the ease of coding, we require the baseline covariate Z to have the same length as Y, even it's technically only N*dz. This could be done by a wrapper function on the outside. 
    arma::cube Beta,
    arma::cube Lambda2,
    arma::cube Nulam,
    arma::mat C,
    arma::mat Tau2,
    arma::mat Nutau,
    arma::mat Theta, // dz*J matrix
    arma::mat IDstart,
    arma::mat IDend,
    arma::mat SigTheta,
    double a
){
  int N = Tau2.n_rows;
  int J = Tau2.n_cols;
  int dz = SigTheta.n_cols;

  List omega = Y; // The auxiliary parameters for updating the beta. Each slice j of the list, same shape (sum_i sum_t) and layout as Y_j.
  List omega1 = Y; // The auxiliary parameters for updating the theta. Each slice also accommodates all observations Y_j
  
  arma::mat Gamma( N, J, fill::zeros ); // not required, just for functions to run
  arma::cube CCube( N, J, iterations/thin, fill::zeros );
  List OmegaList( iterations/thin );
  List BetaList( iterations/thin );
  List Lambda2List( iterations/thin );
  arma::cube Tau2Cube( N, J, iterations/thin, fill::zeros );
  List NulamList( iterations/thin );
  arma::cube NutauCube( N, J, iterations/thin, fill::zeros );
  arma::cube ThetaCube( dz, J, iterations/thin, fill::zeros );
  arma::vec uniqCj;
  
  List clist( 4 ); // to store update_C_ij outputs
  
  int kj = 0; // number of unique cluster at period j
  
  for( int iter = 0; iter < iterations; iter++ ){
    
    omega = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    omega1 = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    
    for( int j = 0; j < J; j++ ){

      for( int i = 0; i < N; i++ ){
        clist = help::update_C_ij( FALSE, i, j, Y, T, Z, Gamma, C, Beta, Theta, Lambda2, Tau2, IDstart, IDend, a);
        C = as<arma::mat>( clist[ 0 ]);
        Beta = as<arma::cube>( clist[ 1 ] );
        Lambda2 = as<arma::cube>( clist[ 2 ] );
        Tau2 = as<arma::mat>( clist[ 3 ] );
      }

      uniqCj = unique( C.col( j ) );
      kj = uniqCj.size(); // size of all clusters: 1+q_p
      
      for( int k = 0; k < kj; k++ ){
        Beta = help::update_beta_jk_PG( j, uniqCj[ k ], Y, T, omega, Z, C, Beta, Theta, Lambda2, Tau2, IDstart, IDend );
        
        Lambda2 = help::update_lambda2_jk( j, uniqCj[ k ], C, Beta, Lambda2, Tau2, Nulam );
        Tau2 = help::update_tau2_jk( j, uniqCj[ k ], C, Beta, Lambda2, Tau2, Nutau );
        Nulam = help::update_nu_lambda_jk( j, uniqCj[ k ], C, Lambda2, Nulam );
        Nutau = help::update_nu_tau_jk( j, uniqCj[ k ], C, Tau2, Nutau );
        
      } // end for k
      
    } // end for j
    
    Theta = help::update_theta_PG(Y, T, Z, omega1, Beta, IDstart, IDend, Theta, SigTheta );
    
    if( ( iter + 1 ) % thin == 0 ){
      BetaList[ ( iter + 1 )/thin - 1  ] = Beta;
      Lambda2List[ ( iter + 1 )/thin - 1  ] = Lambda2;
      NulamList[ ( iter + 1 )/thin - 1  ] = Nulam;
      CCube.slice( ( iter + 1 )/thin - 1  ) = C;
      Tau2Cube.slice( ( iter + 1 )/thin - 1  ) = Tau2;
      NutauCube.slice( ( iter + 1 )/thin - 1  ) = Nutau;
      ThetaCube.slice( ( iter + 1 )/thin - 1  ) = Theta;
    }
    
    // Print out progress
    double printer = iter % 50;
    if( printer == 0 ){
      Rcout << "Iteration = " << iter << endl;
    }
    
  } // end for iter
  
  List output( 7 );
  output[ 0 ] = CCube;
  output[ 1 ] = BetaList;
  output[ 2 ] = Lambda2List;
  output[ 3 ] = Tau2Cube;
  output[ 4 ] = NulamList;
  output[ 5 ] = NutauCube;
  output[ 6 ] = ThetaCube;
  
  return output;
}
