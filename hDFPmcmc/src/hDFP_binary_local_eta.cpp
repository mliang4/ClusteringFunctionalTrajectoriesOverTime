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

// Function :: Calculate multivariate normal log-density (diagonal covariance matrix only)
double log_mvnorm_cpp( arma::vec value, arma::vec mu, arma::vec diag_sig2 ){
  int n = value.size();
  double d = 0;
  
  for (int i = 0; i < n; ++i){
    d += help::log_normal_cpp( value[i], mu[i], diag_sig2[i] );
  }
  
  return d;
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


// // return a matrix with column j removed
// arma::mat remove_colj(arma::mat m, int j){
//   int N = m.n_rows;
//   int J = m.n_cols;
//   arma::mat out( N, J-1, fill::zeros );
//   
//   if( j != 0 ){
//     out.cols( 0, j-1 ) = m.cols( 0, j-1 );
//   }
//   
//   if( j != (J-1) ){
//     out.cols( j, J-2 ) = m.cols( j+1, J-1 );
//   }
//   
//   return out;
// }


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

// return position in a vector given the position ij in a matrix and the number of rows in the matrix N
int id_m2v(int N, int i, int j){
  return j*N + i;
}

// return position ij given the position pos in a vector and the number of rows in the matrix N
arma::vec id_v2m(int N, int pos){
  arma::vec ij( 2, fill::zeros ); 
  ij[1] = pos/N; // j
  ij[0] = pos - N*ij[1]; // i
  return ij;
}

// return the number of clusters eating dish d, with options to remove subject ij/cluster jk
int get_ldotd( 
    arma::mat C, 
    arma::vec D, 
    int mi, // option to remove subject i in period j
    int mj, // need to specify j if removing subject/cluster
    int mk, // option to remove cluster k in period j
    int l // the dish
){
  int N = C.n_rows;
  int J = C.n_cols;
  int ldotd = 0;
  for (int j = 0; j < J; j++){
    arma::vec Cj = C.col(j);
    arma::vec Dj = D.subvec( j*N, (j+1)*N-1 );
    
    if(j == mj && mi >= 0){
      Cj = help::remove_i(Cj, mi); // remove the i^th subject
      Dj = help::remove_i(Dj, mi); // to stay consistent!
    }
    arma::uvec uniqCj_id = find_unique(Cj); // subject rep position of unique clusters
    
    // Rcout << "uniqCj_id \n" << uniqCj_id << endl;
    
    arma::vec clus_dish = Dj( uniqCj_id );
    
    // Rcout << "clus_dish \n" << clus_dish << endl;
    
    if(j == mj && mk >= 0){
      clus_dish = help::remove_i(clus_dish, mk); // remove the k^th clus_dish
    }
    arma::uvec clus_dish_l = find(clus_dish == l);
    ldotd += clus_dish_l.size();
    
    // Rcout << "j \n" << j << endl;
    // Rcout << "ldotd \n" << ldotd << endl;
  }
  
  return ldotd;
  
}

// Function :: remove the gaps within the integer vector (e.g. [0, 1, 3, 5] --> [0, 1, 2, 3])
// Also let the vector start from 0 
arma::vec tighten(arma::vec somevec) {
  arma::vec uni = unique(somevec); // unique elements in the vector
  int unisize = uni.size();
  int start = 0;
  
  for(int i = 0; i < unisize; ++i){
    uvec pos = find( somevec == uni(i) ); // position in the vector that equals to value uni[i]
    arma::vec rpl( pos.size() );
    rpl.fill( start + i );
    somevec( pos ) = rpl;
  }
  
  return somevec;
}


////////////////////////////////////////////////
//              Main functions             /////
////////////////////////////////////////////////

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

// Function to jointly update k and l as described in section 1.2.1
List joint_update_cd_ij(
    int i,
    int j,
    List Y,
    List T,
    List Z,
    arma::cube Beta,
    arma::cube Lambda2,
    arma::mat Gamma,
    arma::mat C,
    arma::mat Theta,
    arma::mat Tau2,
    arma::mat IDstart,
    arma::mat IDend,
    arma::vec D,
    double a_glo,
    double a_loc
){
  
  // THINGS NEEDED FOR THIS UPDATE:
  // item1 - k^{-i}_{jc_old}: number of subjects in cluster c_old after removing i. This is needed for all c_old at period j.
  // item2 - l^{-i}_{d_old}: number of clusters eating dish d_old after removing i. This is needed for all d_old.
  // item3 - l^{-i}: number of clusters in total after removing i. Equal to the sum of all item2.
  
  // Rcout << "--- i " << i << endl;
  // Rcout << "--- j " << j << endl;
  
  int N = IDstart.n_rows;
  int J = IDstart.n_cols;
  arma::vec Yj = Y[j];
  arma::vec Yij = Yj.subvec( IDstart(i,j), IDend(i,j) );
  arma::mat Tj = T[j];
  arma::mat Tij = Tj.rows( IDstart(i,j), IDend(i,j) );
  arma::mat Zj = Z[j];
  arma::mat Zij = Zj.rows( IDstart(i,j), IDend(i,j) );
  int dt = Tj.n_cols;
  int ni = Yij.size();
  
  arma::vec Cj = C.col(j); // table assignment at period j
  arma::vec Cjmi = help::remove_i(Cj, i); // table assignment at period j (i removed)
  arma::vec uniqCjmi = unique( Cjmi ); // unique cluster indicators in C after remove unit i
  int kjmi = uniqCjmi.size(); // number of unique cluster indicators in C after remove unit i
  
  int Didx_ij = help::id_m2v( N, i, j ); // find position ij in vector D
  arma::vec Dmi = help::remove_i(D, Didx_ij); // dish assignment after remove i
  arma::vec uniqDmi = unique( Dmi ); // unique dishes after remove i
  int lmi = uniqDmi.size(); // number of unique dishes after remove i
  
  int l_ = 0; // item3: total number of existing clusters (without subject i)
  for(int jj = 0; jj < J; jj++){
    arma::vec Ctemp;
    if(jj == j){
      Ctemp = Cjmi;
    }else{
      Ctemp = C.col(jj); // make sure the temporary cluster vector leave out subject i
    }
    
    arma::vec uniqCtemp = unique( Ctemp ); // unique clusters
    int uniqCtemp_size = uniqCtemp.size(); // number of unique clusters at period j
    l_ += uniqCtemp_size; // add to d_dotdot for the total number of clusters
  } // end for jj: get item 3
  
  
  arma::vec l_d(lmi, fill::zeros); // number of clusters eating dish l
  for(int d = 0; d < lmi; d++){
    l_d[d] = help::get_ldotd(C, D, i, j, -1, uniqDmi[d]); 
  }
  
  arma::vec C_old = Cj; // a temporary vector to hold cluster indicators to check compatibility after change
  
  // Initialize the candidates for c,d,beta,lambda,tau
  arma::vec post_prob_cdij_vec( kjmi+lmi+1, fill::zeros ); // probabilities of moving to each cluster
  arma::vec candidate_cij( kjmi+lmi+1, fill::zeros ); // candidates for new cluster indicators
  arma::vec candidate_dij( kjmi+lmi+1, fill::zeros ); // candidates for new dish indicators (fill in loop for easier implimentation)
  arma::mat candidate_beta( kjmi+lmi+1, dt, fill::zeros); // candidates for beta 
  arma::mat candidate_lambda2( kjmi+lmi+1, dt, fill::zeros); // candidates for lambda2 
  arma::vec candidate_tau2( kjmi+lmi+1, fill::zeros); // candidates for tau2 
  
  candidate_cij.subvec( 0, kjmi - 1 ) = uniqCjmi; // first kjmi candidate are the previous cluster indicators
  candidate_cij.subvec( kjmi, kjmi+lmi ).fill( 1 + max(uniqCjmi) ); // after kjmi, the candidate cluster indicators are all equal to max+1, indicating a new cluster assignment
  
  arma::vec diagtmp(dt, fill::zeros); // a vector to hold the diagnoal elements of Sigma for new cluster
  arma::mat newSigma; // a matrix to hold Sigma for new cluster
  arma::vec newmu(dt, fill::zeros); // a vector to hold mu for new cluster
  
  double post_prob_cdij;
  arma::mat meantmp;
  
  if(Gamma(i,j) == 0){
    
    // Rcout << "moving" << endl;
    
    // assign to one of the existing clusters
    arma::mat Betaj = Beta.slice(j);
    arma::mat Betajmi = help::remove_rowi(Betaj, i);
    arma::mat Lambda2j = Lambda2.slice(j);
    arma::mat Lambda2jmi = help::remove_rowi(Lambda2j, i);
    arma::vec Tau2j = Tau2.col(j);
    arma::vec Tau2jmi = help::remove_i(Tau2j, i);
    arma::vec Dj = D.subvec( j*N, (j+1)*N-1 ); // dish assignment for cluster j
    arma::vec Djmi = help::remove_i(Dj, i);
    
    for(int cd = 0; cd < kjmi; cd++){ 
      
      // Rcout << "moving to --- uniqCjmi[ cd ] \n" << candidate_cij[ cd ] << endl;
      
      // Find subjects and the number of subjects belong to the proposed cluster
      arma::uvec Sjk = find( Cjmi == uniqCjmi[ cd ] ); // The vector of all subjects in cluster S_jk
      int Sjk_size = Sjk.size(); // number of subjects in this cluster (item1)
      
      // Rcout << "Sjk_size \n" << Sjk_size << endl;
      
      int Sjk0 = Sjk[ 0 ]; // pick one of the subjects in cluster k
      
      candidate_beta.row( cd ) = Betajmi.row( Sjk0 );
      candidate_tau2[ cd ] = Tau2jmi[ Sjk0 ];
      candidate_lambda2.row( cd ) = Lambda2jmi.row( Sjk0 );
      candidate_dij[ cd ] = Djmi[ Sjk0 ];
      
      post_prob_cdij = 0;
      
      for(int t = 0; t < ni; t++ ){
        meantmp = Tij.row(t) * candidate_beta.row( cd ).t() + Zij.row(t) * Theta.col(j);
        post_prob_cdij += help::log_bern_cpp( Yij[t], meantmp[0] );
      }
      
      // Rcout << "candidate_beta_row \n" << candidate_beta.row( cd ) << endl;
      // Rcout << "Likelihood \n" << post_prob_cdij << endl;
      
      post_prob_cdij += log( Sjk_size / (N + a_loc - 1) );
      
      // Rcout << "Likelihood + attract \n" << post_prob_cdij << endl;
      
      post_prob_cdij_vec[ cd ] = post_prob_cdij;
      
    } // end of assign to existing clusters
    
    // new cluster, existing dish
    
    arma::vec ij( 2, fill::zeros ); // to store i and j from vector indices
    
    for(int cd = kjmi; cd < kjmi + lmi; cd++){
      
      // Rcout << "moving to --- new cluster " << candidate_cij[ cd ] << " old dish --- " << uniqDmi[ cd - kjmi ] << endl;
      
      // find unique dishes
      arma::uvec dl = find( D == uniqDmi[ cd - kjmi ] ); // unique dish l (search among all d so that we don't have that -i problem in the indexing)
      int dl_size = l_d[ cd - kjmi ]; // number of table eating this dish (d_dot_l)
      int dl0 = dl[ 0 ]; // one vector index
      ij = id_v2m( N, dl0 );
      
      // set the new cluster to have one of the existing dishes
      arma::mat Beta_tmp = Beta.slice( ij[1] );
      arma::mat Lambda2_tmp = Lambda2.slice( ij[1] );
      arma::vec Tau2_tmp = Tau2.col( ij[1] );
      
      candidate_beta.row( cd ) = Beta_tmp.row( ij[0] );
      candidate_tau2[ cd ] = Tau2_tmp[ ij[0] ];
      candidate_lambda2.row( cd ) = Lambda2_tmp.row( ij[0] );
      candidate_dij[ cd ] = uniqDmi[ cd - kjmi ];
      
      post_prob_cdij = 0;
      for(int t = 0; t < ni; t++ ){
        meantmp = Tij.row(t) * candidate_beta.row( cd ).t() + Zij.row(t) * Theta.col(j);
        post_prob_cdij += help::log_bern_cpp( Yij[t], meantmp[0] );
      }
      // Rcout << "candidate_beta_row \n" << candidate_beta.row( cd ) << endl;
      // Rcout << "Likelihood \n" << post_prob_cdij << endl;
      
      post_prob_cdij += ( log(a_loc / (N + a_loc - 1)) + log(dl_size/(l_ + a_glo)) );
      
      // Rcout << "Likelihood + attract \n" << post_prob_cdij << endl;
      
      post_prob_cdij_vec[ cd ] = post_prob_cdij;
      
    }
    
    // new cluster, new dish
    
    int cd = kjmi + lmi;
    
    // Rcout << "moving to --- new cluster new dish \n" << cd << endl;
    
    candidate_dij[ cd ] = max( uniqDmi ) + 1; // new dish
    candidate_tau2[ cd ] = Rcpp::rgamma( 1, 2, 1/2.0 )[0]; // shape, scale parameters
    for(int d = 0; d < dt; d++){
      candidate_lambda2( cd, d ) = Rcpp::rgamma( 1, 2, 1/2.0 )[0];
      diagtmp[ d ] = candidate_lambda2( cd, d ) * candidate_tau2[ cd ];
    }
    newSigma = diagmat( diagtmp );
    candidate_beta.row( cd ) = help::mvrnormArma(1, newmu, newSigma);
    
    post_prob_cdij = 0;
    for(int t = 0; t < ni; t++ ){
      meantmp = Tij.row(t) * candidate_beta.row( cd ).t() + Zij.row(t) * Theta.col(j);
      post_prob_cdij += help::log_bern_cpp( Yij[t], meantmp[0] );
    }
    // Rcout << "candidate_beta_row \n" << candidate_beta.row( cd ) << endl;
    // Rcout << "Likelihood \n" << post_prob_cdij << endl;
    
    post_prob_cdij += ( log(a_loc / (N + a_loc - 1)) + log(a_glo/(l_ + a_glo)) );
    
    // Rcout << "Likelihood + attract \n" << post_prob_cdij << endl;
    
    post_prob_cdij_vec[ cd ] = post_prob_cdij;
    
    post_prob_cdij_vec = exp( post_prob_cdij_vec - help::logsumexp( post_prob_cdij_vec ) ); // convert log-probabilities to probabilities
    
    // Rcout << "post_prob_cdij_vec \n" << post_prob_cdij_vec << endl;
    
    if( j < (J-1) ){
      // if( DFP && (j < (J-1)) ){
      if( Gamma( i, j + 1 ) == 1 ){ // check compatibility
        for(int k = 0; k < (kjmi+lmi+1); k++){
          C_old[ i ] = candidate_cij[ k ]; // propose a temp partition
          
          // Rcout << "proposed C" << C_old << endl;
          
          bool comp = help::compatible( 0, 0, C_old, C.col(j + 1), Gamma.col( j + 1 ) );
          // set first 2 arguments to 0 to check compatibility of Ctemp and C only, without touching unit i
          if( !comp ){
            post_prob_cdij_vec[ k ] = 0;
          }
        } // end for
      } // end if gamma
    } // end if j < J-1
    
    
    // Rcout << "post_prob_cdij_vec_checked \n" << post_prob_cdij_vec << endl;
    
    // if(sum(post_prob_cdij_vec)>0){
    arma::vec idvec = linspace(0, kjmi+lmi, kjmi+lmi+1); // a sequence from 0 to kjmi
    int new_pos = help::sample_prob_cpp( idvec, post_prob_cdij_vec );
    
    // Rcout << "idvec \n" << idvec << endl;
    // Rcout << "new_pos \n" << new_pos << endl;
    
    C( i, j ) = candidate_cij[ new_pos ];
    Betaj.row( i ) = candidate_beta.row( new_pos );
    Lambda2j.row( i ) = candidate_lambda2.row( new_pos );
    Tau2( i, j ) = candidate_tau2[ new_pos ];
    int d_id_ij = help::id_m2v(N,i,j);
    
    // Rcout << "d_id_ij: " << d_id_ij << endl;
    
    D[d_id_ij] = candidate_dij[ new_pos ];
    
    // Rcout << "new C \n" << candidate_cij[ new_pos ] << endl;
    // Rcout << "new D \n" << candidate_dij[ new_pos ] << endl;
    
    // Rcout << "Betaj \n" << Betaj << endl;
    
    Beta.slice( j ) = Betaj;
    Lambda2.slice( j ) = Lambda2j;
    
    // }else if(sum(post_prob_cdij_vec)>0 && j>0){
    //   
    //   C(i,j) = C(i,j-1);
    //   
    //   int d_id_ij = help::id_m2v(N,i,j);
    //   int d_id_ij1 = help::id_m2v(N,i,j-1);
    //   D[d_id_ij] = D[d_id_ij1];
    //   
    //   Beta( span(i,i), span(0,dt-1), span(j,j) ) = Beta( span(i,i), span(0,dt-1), span(j-1,j-1) );
    //   Lambda2( span(i,i), span(0,dt-1), span(j,j) ) = Lambda2( span(i,i), span(0,dt-1), span(j-1,j-1) );
    //   Tau2(i,j) = Tau2(i,j-1);
    // }
    
    
  }else{
    
    // Rcout << "cannot move" << endl;
    
    // C(i,j) = C(i,j-1);
    // 
    // int d_id_ij = help::id_m2v(N,i,j);
    // int d_id_ij1 = help::id_m2v(N,i,j-1);
    // D[d_id_ij] = D[d_id_ij1];
    // 
    // Beta( span(i,i), span(0,dt-1), span(j,j) ) = Beta( span(i,i), span(0,dt-1), span(j-1,j-1) );
    // Lambda2( span(i,i), span(0,dt-1), span(j,j) ) = Lambda2( span(i,i), span(0,dt-1), span(j-1,j-1) );
    // Tau2(i,j) = Tau2(i,j-1);
    
  } // end if
  
  List output( 5 );
  output[ 0 ] = C;
  output[ 1 ] = D;
  output[ 2 ] = Beta;
  output[ 3 ] = Lambda2;
  output[ 4 ] = Tau2;
  
  return output;
}

// Function to update the dish for the c^th cluster at period j with C being fixed
List cond_update_d_jk(
    int j,
    int k, // the cluster index, not the actual dish itself
    List Y,
    List T,
    List Z,
    arma::cube Beta,
    arma::cube Lambda2,
    arma::mat C,
    arma::mat Theta,
    arma::mat Tau2,
    arma::mat IDstart,
    arma::mat IDend,
    arma::vec D,
    double a_glo
){
  int N = IDstart.n_rows;
  arma::vec Yj = Y[j];
  arma::mat Tj = T[j];
  arma::mat Zj = Z[j];
  int dt = Tj.n_cols;
  
  // get total number of existing dishes (without cluster k)
  arma::vec Cj = C.col(j); // table assignment at period j
  arma::vec uniqCj = unique( Cj ); // should be a length bigger than k
  
  // Rcout << "uniqCj \n" << uniqCj << endl;
  
  arma::uvec pos_cjk = find(Cj == uniqCj[k]); // index i of the subjects in cluster Cjk
  int cjk_size = pos_cjk.size(); // number of subjects in cluster Cjk
  D( j*N + pos_cjk ).fill(-999); // set dish indices for cluster Cjk to -999
  
  // Rcout << "D \n" << D << endl;
  
  arma::vec Dmjk = D(find(D>=0)); // vector D (without cluster Cjk)
  
  // Rcout << "Dmjk \n" << Dmjk << endl;
  
  arma::vec uniq_dish = unique( Dmjk ); // unique dish (without cluster Cjk)
  int n_dish = uniq_dish.size(); // number of unique dishes (without cluster Cjk)
  
  // Rcout << "n_dish \n" << n_dish << endl;
  
  // initialize candidates (beta, lambda2, tau2, )
  arma::vec post_prob_ljk_vec( n_dish+1, fill::zeros ); // probabilities of moving to each cluster
  arma::vec candidate_ljk( n_dish+1, fill::zeros ); // candidates for new dish indicators
  arma::mat candidate_beta( n_dish+1, dt, fill::zeros); // candidates for beta 
  arma::mat candidate_lambda2( n_dish+1, dt, fill::zeros); // candidates for lambda2 
  arma::vec candidate_tau2( n_dish+1, fill::zeros); // candidates for tau2 
  
  // for each existing dishes, get total number of tables eating that dish (with jk^th cluster removed)
  arma::vec ldotd( n_dish, fill::zeros );
  for(int l = 0; l < n_dish; l++){
    candidate_ljk[ l ] = uniq_dish[ l ];
    ldotd[l] = help::get_ldotd(C,D,-1,j,k,uniq_dish[l]);
  }
  
  // Rcout << "ldotd \n" << ldotd << endl;
  
  int ddotdot = sum(ldotd); // total number of clusters after removing the jk^th cluster
  
  // Rcout << "ddotdot \n" << ddotdot << endl;
  
  arma::vec diagtmp(dt, fill::zeros); // a vector to hold the diagnoal elements of Sigma for new cluster
  arma::mat newSigma; // a matrix to hold Sigma for new cluster
  arma::vec newmu(dt, fill::zeros); // a vector to hold mu for new cluster
  
  double post_prob_ljk;
  arma::mat meantmp;
  
  arma::mat Betaj = Beta.slice( j );
  arma::mat Lambda2j = Lambda2.slice( j );
  arma::vec Tau2j = Tau2.col( j );
  
  // calculate the likelihood by summing over all subjects in cluster jk
  for(int l = 0; l < n_dish; l++){
    
    // Find 1 subject eating the proposed dish
    arma::uvec dl = find( D == uniq_dish[ l ] ); // The vector of all subjects that eating dish l
    int dl0 = dl[ 0 ]; // pick one of the subjects
    arma::vec dl0_ij = help::id_v2m(N, dl0);
    
    arma::mat Beta_tmpj = Beta.slice( dl0_ij[1] ); // tmpj is the period subject dl0 is at
    arma::mat Lambda2_tmpj = Lambda2.slice( dl0_ij[1] );
    arma::vec Tau2_tmpj = Tau2.col( dl0_ij[1] );
    
    candidate_beta.row( l ) = Beta_tmpj.row( dl0_ij[0] );
    candidate_tau2[ l ] = Tau2_tmpj[ dl0_ij[0] ];
    candidate_lambda2.row( l ) = Lambda2_tmpj.row( dl0_ij[0] );
    
    post_prob_ljk = 0;
    // find corresponding parameters, let everyone in CLUSTER JK to have this coefficients, calculate likelihood
    // arma::uvec djk = find( D == -999 ); // D index for everyone in cluster jk
    for(int tmpi = 0; tmpi < cjk_size; tmpi++){
      arma::vec Yij = Yj.subvec( IDstart(pos_cjk[tmpi],j), IDend(pos_cjk[tmpi],j) );
      arma::mat Tij = Tj.rows( IDstart(pos_cjk[tmpi],j), IDend(pos_cjk[tmpi],j) );
      arma::mat Zij = Zj.rows( IDstart(pos_cjk[tmpi],j), IDend(pos_cjk[tmpi],j) );
      int ni = Yij.size();
      for(int tmpii = 0; tmpii < ni; tmpii++){ // calculating the likelihood for one subject
        meantmp = Tij.row(tmpii) * candidate_beta.row( l ).t() + Zij.row(tmpii) * Theta.col(j);
        post_prob_ljk += help::log_bern_cpp( Yij[tmpii], meantmp[0] );
      } // end of calculating the likelihood for one subject
    } // end of summing the likelihood for all subjects in cluster jk
    
    post_prob_ljk += log( ldotd[l] / (ddotdot + a_glo) );
    post_prob_ljk_vec[ l ] = post_prob_ljk;
  } // end of each existing dish l
  
  // propose new dish + calculate likelihood
  candidate_ljk[ n_dish ] = max(uniq_dish) + 1; // new dish label
  
  candidate_tau2[ n_dish ] = Rcpp::rgamma( 1, 2, 1/2.0 )[0];
  for(int d = 0; d < dt; d++){
    candidate_lambda2( n_dish, d ) = Rcpp::rgamma( 1, 2, 1/2.0 )[0];
    diagtmp[ d ] = candidate_lambda2( n_dish, d ) * candidate_tau2[ n_dish ];
  }
  newSigma = diagmat( diagtmp );
  candidate_beta.row( n_dish ) = help::mvrnormArma(1, newmu, newSigma);
  
  post_prob_ljk = 0;
  for(int tmpi = 0; tmpi < cjk_size; tmpi++){
    arma::vec Yij = Yj.subvec( IDstart(pos_cjk[tmpi],j), IDend(pos_cjk[tmpi],j) );
    arma::mat Tij = Tj.rows( IDstart(pos_cjk[tmpi],j), IDend(pos_cjk[tmpi],j) );
    arma::mat Zij = Zj.rows( IDstart(pos_cjk[tmpi],j), IDend(pos_cjk[tmpi],j) );
    int ni = Yij.size();
    for(int tmpii = 0; tmpii < ni; tmpii++){ // calculating the likelihood for one subject
      meantmp = Tij.row(tmpii) * candidate_beta.row( n_dish ).t() + Zij.row(tmpii) * Theta.col(j);
      post_prob_ljk += help::log_bern_cpp( Yij[tmpii], meantmp[0] );
    } // end of calculating the likelihood for one subject
  } // end of summing the likelihood for all subjects in cluster jk
  
  post_prob_ljk += log( a_glo / (ddotdot + a_glo) );
  post_prob_ljk_vec[ n_dish ] = post_prob_ljk;
  
  // sample and save decision
  post_prob_ljk_vec = exp( post_prob_ljk_vec - help::logsumexp( post_prob_ljk_vec ) ); // convert log-probabilities to probabilities
  
  arma::vec idvec = linspace(0, n_dish, n_dish+1); // a sequence from 0 to n_dish
  int new_pos = help::sample_prob_cpp( idvec, post_prob_ljk_vec ); // this is the new D index
  D( j*N + pos_cjk ).fill( candidate_ljk[ new_pos ] ); // set dish indices for cluster Cjk to this new dish
  for(int tmpi = 0; tmpi < cjk_size; tmpi++){
    Betaj.row( pos_cjk[tmpi] ) = candidate_beta.row( new_pos );
    Lambda2j.row( pos_cjk[tmpi] ) = candidate_lambda2.row( new_pos );
    Tau2( pos_cjk[tmpi], j ) = candidate_tau2[ new_pos ];
  } // end for all subjects in cluster jk
  
  Beta.slice( j ) = Betaj;
  Lambda2.slice( j ) = Lambda2j;
  
  List output( 4 );
  output[ 0 ] = D;
  output[ 1 ] = Beta;
  output[ 2 ] = Lambda2;
  output[ 3 ] = Tau2;
  
  return output;
  
  
}


// function to update the beta for the l^th dish
arma::cube update_beta_l(
    int dl, // the actual dish (not the index)
    List Y,
    List T,
    List omega_beta,
    List Z,
    arma::mat C,
    arma::vec D,
    arma::cube Beta,
    arma::mat Theta,
    arma::cube Lambda2,
    arma::mat Tau2,
    arma::mat IDstart,
    arma::mat IDend
){
  
  // In vector D, find the position that are eating dish dl, those are the ij's we need to look at
  
  int N = IDstart.n_rows;
  uvec pos_l = find( D == dl ); // subjects (D index) that eat dish l
  
  // Initialize posterior Sig and mean
  
  int pos_l0 = pos_l[0];
  arma::vec pos_ij0 = help::id_v2m(N, pos_l0);
  // -- create some temporary variables for the prior
  arma::mat Lambda2j = Lambda2.slice( pos_ij0[1] );
  arma::vec Tau2j = Tau2.col( pos_ij0[1] );
  int DT = Lambda2j.n_cols;
  arma::vec s2_diag( DT );
  for(int dT = 0; dT < DT; dT++){
    s2_diag[ dT ] = Lambda2j( pos_ij0[0], dT ) * Tau2j[ pos_ij0[0] ];
  }
  arma::mat Sigma_star_inv = diagmat( 1 / s2_diag );
  arma::mat Sigma_bl_cs = Sigma_star_inv; // prior variance as the starting point
  
  // Rcout << Sigma_bl_cs << endl;
  
  arma::vec mu_bl_cs( DT, fill::zeros ); // empty vector as the starting point
  arma::mat Tij, omega_mat_ij;
  arma::vec Zij, omegaij;
  
  // For each ij, we add together the t(T)OmegaT and t(T)OmegaZbar for the post Sig and mean
  int dl_size = pos_l.size(); // number of subjects eating dish l
  for(int l_id=0; l_id<dl_size; l_id++){
    int l = pos_l[l_id]; // D index of one subject
    arma::vec ij = help::id_v2m(N, l); // C index of one subject
    
    // todo: revise this
    int is = IDstart( ij[0], ij[1] );
    int ie = IDend( ij[0], ij[1] );
    
    arma::mat Tj = T[ ij[1] ];
    arma::mat Zj = Z[ ij[1] ];
    arma::vec muij = Zj * Theta.col(ij[1]);
    arma::vec Yj = Y[ ij[1] ];
    arma::vec Kj = Yj - 0.5;
    
    arma::vec omegaj = omega_beta[ ij[1] ];
    arma::vec zbarj = Kj/omegaj - muij;
    
    Tij = Tj.rows( is, ie );
    Zij = zbarj.subvec( is, ie );
    omegaij = omegaj.subvec( is, ie );
    omega_mat_ij = diagmat( omegaij );
    
    Sigma_bl_cs += Tij.t() * omega_mat_ij * Tij;
    mu_bl_cs += Tij.t() * omega_mat_ij * Zij;
  }
  
  Sigma_bl_cs = inv( Sigma_bl_cs );
  mu_bl_cs = Sigma_bl_cs * mu_bl_cs;
  
  // Rcout << Sigma_bl_cs << endl;
  
  arma::mat beta_new = help::mvrnormArma( 1, mu_bl_cs, Sigma_bl_cs );  // Sample posterior beta
  
  // Replace the original beta for the corresponding ij's
  for(int l_id=0; l_id<dl_size; l_id++){
    int l = pos_l[l_id]; // D index of one subject
    arma::vec ij = help::id_v2m(N, l); // C index of one subject
    Beta( span(ij[0],ij[0]), span(0,DT-1), span(ij[1],ij[1]) ) = beta_new;
  }
  
  return Beta;
}


// function to update the auxiliary parameters for Beta
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

arma::mat update_theta_PG(
    List Y, 
    List T, 
    List Z,
    List omega_theta, 
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
    
    arma::vec omega_thetaj = omega_theta[j];
    Ztildej = Kj/omega_thetaj;
    
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
    
    Obarj = diagmat( omega_thetaj );
    
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

arma::cube update_lambda2_l( 
    int dl,
    arma::vec D,
    arma::cube Beta,
    arma::cube Lambda2,
    arma::mat Tau2,
    arma::cube Nulam
){
  // Find subjects and the number of subjects belong to kth cluster
  uvec pos_l = find( D == dl ); // subjects (D index) that eat dish l
  int Dl_size = pos_l.size();
  int dl0 = pos_l[ 0 ]; // one of the subjects in cluster k
  int N = Beta.n_rows;
  int DT = Beta.n_cols;
  arma::vec ij0 = help::id_v2m(N, dl0);
  
  double beta_ld = 0;
  double tau2_l = 0;
  double nulam_ld = 0;
  double bterm = 0;
  double lambda2_new = 0;
  
  for(int d = 0; d < DT; d++){
    
    beta_ld = Beta( ij0[0], d, ij0[1] );
    tau2_l = Tau2( ij0[0], ij0[1] );
    nulam_ld = Nulam( ij0[0], d, ij0[1] );
    bterm = beta_ld * beta_ld / tau2_l / 2 + 1 / nulam_ld;
    lambda2_new = help::rinvgamma_cpp(1, bterm);
    
    for(int idx = 0; idx < Dl_size; ++idx){ // update lambda2_new for every subject
      int dl = pos_l[ idx ]; // one of the subjects in cluster k
      arma::vec ij = help::id_v2m(N, dl);
      Lambda2( ij[0], d, ij[1] ) = lambda2_new;
    }
  }
  
  return Lambda2;
}

arma::mat update_tau2_l( 
    int dl,
    arma::mat D,
    arma::cube Beta,
    arma::cube Lambda2,
    arma::mat Tau2,
    arma::mat Nutau
){
  uvec pos_l = find( D == dl ); 
  int Dl_size = pos_l.size();
  int dl0 = pos_l[ 0 ]; // one of the subjects in cluster k
  int N = Beta.n_rows;
  int DT = Beta.n_cols;
  arma::vec ij0 = help::id_v2m(N, dl0);
  
  double beta_ld = 0;
  double lambda2_ld = 0;
  
  double aterm = 0.5 * (DT+1);
  double bterm = 1/Nutau(ij0[0], ij0[1]); // start with 1/nu_tau_l*
  
  for(int d = 0; d < DT; d++){
    beta_ld = Beta( ij0[0], d, ij0[1] );
    lambda2_ld = Lambda2( ij0[0], d, ij0[1] );
    bterm += beta_ld * beta_ld / lambda2_ld / 2;
  }
  
  double tau2_new = help::rinvgamma_cpp(aterm, bterm);
  
  for(int idx = 0; idx < Dl_size; ++idx){
    int dl = pos_l[ idx ]; // one of the subjects in cluster k
    arma::vec ij = help::id_v2m(N, dl);
    Tau2( ij[0], ij[1] ) = tau2_new;
  }
  
  return Tau2;
}


arma::cube update_nu_lambda_l( 
    int dl,
    arma::mat D,
    arma::cube Lambda2,
    arma::cube Nulam
){  
  uvec pos_l = find( D == dl ); 
  int Dl_size = pos_l.size();
  int dl0 = pos_l[ 0 ]; // one of the subjects in cluster k
  int N = Lambda2.n_rows;
  int DT = Lambda2.n_cols;
  arma::vec ij0 = help::id_v2m(N, dl0);
  
  double lambda2_ld = 0;
  double bterm = 0;
  double nulam_new = 0;
  
  for(int d = 0; d < DT; d++){
    lambda2_ld = Lambda2( ij0[0], d, ij0[1] );
    bterm = 1 + 1 / lambda2_ld;
    nulam_new = help::rinvgamma_cpp(1, bterm);
    
    for(int idx = 0; idx < Dl_size; ++idx){
      int dl = pos_l[ idx ]; // one of the subjects in cluster k
      arma::vec ij = help::id_v2m(N, dl);
      Nulam( ij[0], d, ij[1] ) = nulam_new;
    }
  }
  
  return Nulam;
}

arma::mat update_nu_tau_l( 
    int dl,
    arma::mat D,
    arma::mat Tau2,
    arma::mat Nutau 
){
  uvec pos_l = find( D == dl ); 
  int Dl_size = pos_l.size();
  int dl0 = pos_l[ 0 ]; // one of the subjects in cluster k
  int N = Tau2.n_rows;
  arma::vec ij0 = help::id_v2m(N, dl0);
  double nutau_new = help::rinvgamma_cpp( 1, Tau2( ij0[0], ij0[1]) );
  
  for(int idx = 0; idx < Dl_size; ++idx){
    int dl = pos_l[ idx ]; // one of the subjects in cluster k
    arma::vec ij = help::id_v2m(N, dl);
    Nutau( ij[0], ij[1] ) = nutau_new;
  }
  return Nutau;
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

// function to update the auxiliary parameters for Eta
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
List hDFPmcmc(
    int iterations,
    int thin,
    List Y,     // We use list to store the observation because the number of Ys at each period could be different. Each slice is a vector of length sum(i,t): [y1,...,yN]
    List T,     // For the same reason, we use a list to store the basis functions T for each period. Each slice is a matrix of size sum(i,t)*dt
    List Z,     // For the ease of coding, we require the baseline covariate Z to have the same length as Y, even it's technically only N*dz. This could be done by a wrapper function on the outside. 
    arma::cube X, // N*dx*J cube
    arma::cube Beta, // N*dt*J cube
    arma::cube Lambda2, // N*dt*J cube
    arma::cube Nulam, // N*dt*J cube
    arma::mat Gamma, // N*J matrix
    arma::mat C, // N*J matrix (table assignment)
    arma::mat Tau2, // N*J matrix
    arma::mat Nutau, // N*J matrix
    arma::mat Eta, // dx*J matrix
    arma::mat Theta, // dz*J matrix
    arma::mat IDstart, // N*J
    arma::mat IDend, // N*J
    arma::mat SigEta,
    arma::mat SigTheta,
    arma::vec D, // N*J vector (dish assignment)
    arma::vec MuEta,
    double a_glo, // concentration parameter at the global (dish) level
    double a_loc // concentration parameter at the local (table) level
){
  int N = Beta.n_rows;
  int J = Beta.n_slices;
  int dx = X.n_cols;
  int dz = Theta.n_rows;
  
  List omega_beta = Y; // The auxiliary parameters for updating the beta. Each slice j of the list has the same shape (sum_i sum_t) and layout as Y_j.
  List omega_theta = Y; // The auxiliary parameters for updating the theta. Same shape/layout as omega_beta
  arma::mat w( N, J, fill::zeros ); // The auxiliary parameters for updating Eta
  
  // Create these to store the output:
  List BetaList( iterations/thin );
  List Lambda2List( iterations/thin );
  List NulamList( iterations/thin );
  arma::cube GammaCube( N, J, iterations/thin, fill::zeros );
  arma::cube CCube( N, J, iterations/thin, fill::zeros );
  arma::cube Tau2Cube( N, J, iterations/thin, fill::zeros );
  arma::cube NutauCube( N, J, iterations/thin, fill::zeros );
  arma::cube EtaCube( dx, J, iterations/thin, fill::zeros );
  arma::cube ThetaCube( dz, J, iterations/thin, fill::zeros );
  arma::mat DMat( iterations/thin, N*J, fill::zeros );
  
  List cdlist( 5 ); // to store joint_update_cd_ij outputs
  List dlist( 4 ); 
  
  for( int iter = 0; iter < iterations; iter++ ){
    
    omega_beta = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    omega_theta = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    
    for( int j = 0; j < J; j++ ){
      
      for( int i = 0; i < N; i++ ){
        Gamma = help::update_gamma_ij( i, j, X, Gamma, C, Eta, a_loc );
      }
      
      // Rcout << "gamma" << Gamma << endl;
      
      // Rcout << "update cd" << endl;
      
      for( int i = 0; i < N; i++ ){
        cdlist = help::joint_update_cd_ij( i, j, Y, T, Z, Beta, Lambda2, Gamma, C, Theta, Tau2, IDstart, IDend, D, a_glo, a_loc );
        C = as<arma::mat>( cdlist[ 0 ] );
        D = as<arma::vec>( cdlist[ 1 ] );
        Beta = as<arma::cube>( cdlist[ 2 ] );
        Lambda2 = as<arma::cube>( cdlist[ 3 ] );
        Tau2 = as<arma::mat>( cdlist[ 4 ] );
        
        // if(j < (J-1)){
        //   Rcout << "C_ij \n" << C(i,j) << endl;
        //   Rcout << "C_ij+1 \n" << C(i,j+1) << endl;
        //   Rcout << "D_ij \n" << D[j*N + i] << endl;
        //   Rcout << "D_ij+1 \n" << D[(j+1)*N + i] << endl;
        // }
        
        
      }
      
      C.col(j) = help::tighten(C.col(j));
      D = help::tighten(D);
      
      
      
      // Rcout << "update d" << endl;
      
      arma::vec uniqCj = unique( C.col( j ) );
      int kj = uniqCj.size(); // number of clusters at period j
      
      for( int k = 0; k < kj; k++ ){
        dlist = help::cond_update_d_jk( j, k, Y, T, Z, Beta, Lambda2, C, Theta, Tau2, IDstart, IDend, D, a_glo );
        D = as<arma::vec>( dlist[ 0 ] );
        Beta = as<arma::cube>( dlist[ 1 ] );
        Lambda2 = as<arma::cube>( dlist[ 2 ] );
        Tau2 = as<arma::mat>( dlist[ 3 ] );
        
        
      } // end for k
      
      D = help::tighten(D);
      // 
      // Rcout << "C \n" << C << endl;
      // Rcout << "D \n" << D << endl;
      
    } // end for j
    
    arma::vec uniqD = unique( D );
    int L = uniqD.size();
    for( int l = 0; l < L; l++ ){
      Beta = help::update_beta_l( uniqD[l], Y,T,omega_beta,Z,C,D,Beta,Theta,Lambda2,Tau2,IDstart,IDend);
      Lambda2 = help::update_lambda2_l( uniqD[l], D, Beta, Lambda2, Tau2, Nulam );
      Tau2 = help::update_tau2_l( uniqD[l], D, Beta, Lambda2, Tau2, Nutau );
      Nulam = help::update_nu_lambda_l( uniqD[l], D, Lambda2, Nulam );
      Nutau = help::update_nu_tau_l( uniqD[l], D, Tau2, Nutau );
    }  
    
    w = help::update_w(X, Eta);
    Eta = help::update_eta_PG(X, Gamma, Eta, MuEta, SigEta, w);
    Theta = help::update_theta_PG(Y, T, Z, omega_theta, Beta, IDstart, IDend, Theta, SigTheta );
    
    
    
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
      DMat.row( ( iter + 1 )/thin - 1  ) = D.t();
    }
    
    
    // Print out progress
    double printer = iter % 50;
    if( printer == 0 ){
      Rcout << "Iteration = " << iter << endl;
    }
    
  } // end for iter
  
  List output( 10 );
  output[ 0 ] = BetaList;
  output[ 1 ] = Lambda2List;
  output[ 2 ] = NulamList;
  output[ 3 ] = GammaCube;
  output[ 4 ] = CCube;
  output[ 5 ] = Tau2Cube;
  output[ 6 ] = NutauCube;
  output[ 7 ] = EtaCube;
  output[ 8 ] = ThetaCube;
  output[ 9 ] = DMat;
  
  return output;
}


// [[Rcpp::export]]
List hDPmcmc(
    int iterations,
    int thin,
    List Y,     // We use list to store the observation because the number of Ys at each period could be different. Each slice is a vector of length sum(i,t): [y1,...,yN]
    List T,     // For the same reason, we use a list to store the basis functions T for each period. Each slice is a matrix of size sum(i,t)*dt
    List Z,     // For the ease of coding, we require the baseline covariate Z to have the same length as Y, even it's technically only N*dz. This could be done by a wrapper function on the outside. 
    arma::cube Beta, // N*dt*J cube
    arma::cube Lambda2, // N*dt*J cube
    arma::cube Nulam, // N*dt*J cube
    arma::mat C, // N*J matrix (table assignment)
    arma::mat Tau2, // N*J matrix
    arma::mat Nutau, // N*J matrix
    arma::mat Theta, // dz*J matrix
    arma::mat IDstart, // N*J
    arma::mat IDend, // N*J
    arma::mat SigTheta,
    arma::vec D, // N*J vector (dish assignment)
    double a_glo, // concentration parameter at the global (dish) level
    double a_loc // concentration parameter at the local (table) level
){
  int N = Beta.n_rows;
  int J = Beta.n_slices;
  
  int dz = Theta.n_rows;
  
  List omega_beta = Y; // The auxiliary parameters for updating the beta. Each slice j of the list has the same shape (sum_i sum_t) and layout as Y_j.
  List omega_theta = Y; // The auxiliary parameters for updating the theta. Same shape/layout as omega_beta
  
  arma::mat Gamma( N, J, fill::zeros ); // not required, just for functions to run
  
  // Create these to store the output:
  List BetaList( iterations/thin );
  List Lambda2List( iterations/thin );
  List NulamList( iterations/thin );
  arma::cube CCube( N, J, iterations/thin, fill::zeros );
  arma::cube Tau2Cube( N, J, iterations/thin, fill::zeros );
  arma::cube NutauCube( N, J, iterations/thin, fill::zeros );
  arma::cube ThetaCube( dz, J, iterations/thin, fill::zeros );
  arma::mat DMat( iterations/thin, N*J, fill::zeros );
  
  List cdlist( 5 ); // to store joint_update_cd_ij outputs
  List dlist( 4 ); 
  
  for( int iter = 0; iter < iterations; iter++ ){
    
    omega_beta = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    omega_theta = help::update_omega(T, Z, Beta, Theta, IDstart, IDend);
    
    for( int j = 0; j < J; j++ ){
      
      for( int i = 0; i < N; i++ ){
        cdlist = help::joint_update_cd_ij( i, j, Y, T, Z, Beta, Lambda2, Gamma, C, Theta, Tau2, IDstart, IDend, D, a_glo, a_loc );
        // cdlist = help::joint_update_cd_ij( FALSE, i, j, Y, T, Z, Beta, Lambda2, Gamma, C, Theta, Tau2, IDstart, IDend, D, a_glo, a_loc );
        C = as<arma::mat>( cdlist[ 0 ] );
        D = as<arma::vec>( cdlist[ 1 ] );
        Beta = as<arma::cube>( cdlist[ 2 ] );
        Lambda2 = as<arma::cube>( cdlist[ 3 ] );
        Tau2 = as<arma::mat>( cdlist[ 4 ] );
      }
      
      arma::vec uniqCj = unique( C.col( j ) );
      int kj = uniqCj.size(); // number of clusters at period j
      
      for( int k = 0; k < kj; k++ ){
        dlist = help::cond_update_d_jk( j, k, Y, T, Z, Beta, Lambda2, C, Theta, Tau2, IDstart, IDend, D, a_glo );
        D = as<arma::vec>( dlist[ 0 ] );
        Beta = as<arma::cube>( dlist[ 1 ] );
        Lambda2 = as<arma::cube>( dlist[ 2 ] );
        Tau2 = as<arma::mat>( dlist[ 3 ] );
      } // end for k
      
      D = help::tighten(D);
      
    } // end for j
    
    arma::vec uniqD = unique( D );
    int L = uniqD.size();
    for( int l = 0; l < L; l++ ){
      Beta = help::update_beta_l( uniqD[l], Y,T,omega_beta,Z,C,D,Beta,Theta,Lambda2,Tau2,IDstart,IDend);
      Lambda2 = help::update_lambda2_l( uniqD[l], D, Beta, Lambda2, Tau2, Nulam );
      Tau2 = help::update_tau2_l( uniqD[l], D, Beta, Lambda2, Tau2, Nutau );
      Nulam = help::update_nu_lambda_l( uniqD[l], D, Lambda2, Nulam );
      Nutau = help::update_nu_tau_l( uniqD[l], D, Tau2, Nutau );
    }  
    
    Theta = help::update_theta_PG(Y, T, Z, omega_theta, Beta, IDstart, IDend, Theta, SigTheta );
    
    if( ( iter + 1 ) % thin == 0 ){
      BetaList[ ( iter + 1 )/thin - 1  ] = Beta;
      Lambda2List[ ( iter + 1 )/thin - 1  ] = Lambda2;
      NulamList[ ( iter + 1 )/thin - 1  ] = Nulam;
      CCube.slice( ( iter + 1 )/thin - 1  ) = C;
      Tau2Cube.slice( ( iter + 1 )/thin - 1  ) = Tau2;
      NutauCube.slice( ( iter + 1 )/thin - 1  ) = Nutau;
      ThetaCube.slice( ( iter + 1 )/thin - 1  ) = Theta;
      DMat.row( ( iter + 1 )/thin - 1  ) = D.t();
    }
    
    // Print out progress
    double printer = iter % 50;
    if( printer == 0 ){
      Rcout << "Iteration = " << iter << endl;
    }
    
  } // end for iter
  
  List output( 8 );
  output[ 0 ] = BetaList;
  output[ 1 ] = Lambda2List;
  output[ 2 ] = NulamList;
  output[ 3 ] = CCube;
  output[ 4 ] = Tau2Cube;
  output[ 5 ] = NutauCube;
  output[ 6 ] = ThetaCube;
  output[ 7 ] = DMat;
  
  return output;
}
