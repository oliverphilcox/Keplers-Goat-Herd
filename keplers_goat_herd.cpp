// keplers_goat_herd.cpp -- Oliver Philcox, 2021.
// Solve Kepler's equation with a variety of numerical techniques for an array of mean anomalies
//
// Methods implemented are:
//   - Newton-Raphson: The quadratic Newton-Raphson root finder.
//   - Danby: The quartic root finder described in Danby (1988).
//   - Series: An elliptic series method, as described in Murray & Dermott.
//   - Contour: A method based on contour integration, described in Philcox et al. (2021).
//
// Given an array of mean anomalies, an eccentricity and a desired precision, the code will estimate the eccentric anomaly using each method.
// We generate a grid of mean anomalies from a uniformly spaced grid in eccentric anomaly.
// The hyperparameter of each approach is increased until the desired precision is reached, relative to the true result.
// Timing is performed using the `chrono' package.

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <complex>
#include <chrono>
#include <cmath>
using namespace std::chrono;

// Could swap between single and double precision here.
typedef double Float;

// ========================== Class to hold basic attributes ================

class Approximations {
  private:
    int N_ell; // number of ell points
    Float e; // eccentricity
    Float* ell_arr; // array of mean anomalies
  public:
    Approximations(int _N_ell, Float _e, Float* _ell_arr){
      N_ell = _N_ell;
      e = _e;

      // Create and copy ell array
      ell_arr = new Float[N_ell];
      for(int i=0;i<N_ell;i++) ell_arr[i] = _ell_arr[i];

    }
    ~Approximations(){};

// ========================== Approximation Methods ================

  void compute_newton_raphson(int N_it, Float *output){
    /// Newton-Raphson method 
    //
    // Hanno Rein: - Improved initial guess. 
    //             - Combined operations in iteration.
    //             - sincos used to speed up sin/cos calculations

    Float this_ell, old_E, sinE, cosE;

    old_E = 0;

    for(int i=0;i<N_ell;i++){
      this_ell = ell_arr[i];
      
      // Original initial estimate
      //if((sin(this_ell))<0) old_E = this_ell - 0.85*e;
      //else old_E = this_ell + 0.85*e;

      // Define initial estimate (fourth order in e)
      sincos(this_ell,&sinE,&cosE);
      old_E = this_ell + e*sinE/sqrt(1.-2.*e*cosE+e*e);

      // Perform Newton-Raphson estimate
      for(int j=0;j<N_it;j++) {

        // Compute sin and cos in one step.
        // Note: the compiler is usually smart enough 
        // to notice that optimization automatically and
        // this does therefore not improve the speed 
        // for most compiler options.
        // Might need to use __sincos() instead of sincos().
        sincos(old_E,&sinE,&cosE);

        // Combine update in one step
        Float new_E = (this_ell - e*(old_E*cosE-sinE))/(1.-e*cosE);

        old_E = new_E;
      }

      // Add to array
      output[i] = old_E;
    }
  }

  void compute_danby(int N_it, Float *output){
    /// Compute Danby (1988) fourth-order root-finding 
    ///  method with given number of steps
    //
    // Hanno Rein: - Improved initial guess. 
    //             - Combined operations in iteration.
    //             - sincos used to speed up sin/cos calculations
    // This has quartic convergence.
    // The initial step is defined as E_0 = ell + sgn(sin(ell))*e*k following Danby (1988)

    Float f_E, fP_E, fPP_E, fPPP_E, this_ell, old_E, delta_i2, delta_i3, esinE, ecosE;

    for(int i=0;i<N_ell;i++){
      this_ell = ell_arr[i];

      // Define initial estimate (second order in e)
      old_E = this_ell + e*sin(this_ell);

      // Perform Newton-Raphson estimate
      for(int j=0;j<N_it;j++) {

        // Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
        sincos(old_E,&esinE,&ecosE);
        esinE = e*esinE;
        ecosE = e*ecosE;
        f_E = old_E - esinE-this_ell;
        fP_E = 1. - ecosE;
        fPP_E = esinE/2.;
        fPPP_E = ecosE;

        delta_i2 = -f_E*fP_E/(fP_E*fP_E-f_E*fPP_E);
        delta_i3 = -f_E/(fP_E + delta_i2*(fPP_E + fPPP_E*delta_i2/6.));

        // Update E
        old_E += delta_i3;
      }

      // Add to array
      output[i] = old_E;
    }
  }

  void compute_series(int N_it, Float *output){
    // Solve Kepler's equation via the series method described in Murray & Dermot.
    // We use a specifies maximum number of iterations and precompute the coupling coefficients.

    Float coeff;

    // Take initial guess
    for(int i=0; i<N_ell;i++){
      output[i] = ell_arr[i];
    }

    // Iterate over number of iterations
    for(int s=1; s<=N_it; s++){
        coeff = 2*jn(s,s*e)/s; // define coefficient

        for(int i=0;i<N_ell;i++){
          output[i] += coeff*sin(s*ell_arr[i]);
        }
    }
  }

void compute_contour(int N_it, Float *output){
    // Solve Kepler's equation via the contour integration method of Philcox et al. (2021)
    // This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.
    // N_it specifies the number of grid-points.

    Float ft_gx2, ft_gx1, this_ell, freq, zR, zI, cosC, sinC, esinRadius, ecosRadius, center;
    Float fxR, fxI, ftmp, tmpcosh, tmpsinh, tmpcos, tmpsin;

    // Define sampling points (actually use one more than this)
    int N_points = N_it-2;
    int N_fft = (N_it-1)*2;

    // Define contour radius
    Float radius = e/2;

    // Generate e^{ikx} sampling points and precompute real and imaginary parts
    Float exp2R[N_points], exp2I[N_points], exp4R[N_points], exp4I[N_points], coshI[N_points], sinhI[N_points], ecosR[N_points], esinR[N_points];
    for(int jj=0;jj<N_points;jj++){
      // NB: j = jj+1
      freq = 2.0*M_PI*(jj+1)/N_fft;
      exp2R[jj] = cos(freq);
      exp2I[jj] = sin(freq);
      exp4R[jj] = cos(2.0*freq);
      exp4I[jj] = sin(2.0*freq);
      coshI[jj] = cosh(radius*exp2I[jj]);
      sinhI[jj] = sinh(radius*exp2I[jj]);
      ecosR[jj] = e*cos(radius*exp2R[jj]);
      esinR[jj] = e*sin(radius*exp2R[jj]);
    }

    // Precompute e sin(e/2) and e cos(e/2)
    esinRadius = e*sin(radius);
    ecosRadius = e*cos(radius);

    // Iterate over array of mean anomalies
    for(int i=0;i<N_ell;i++){
      this_ell = ell_arr[i];

      // Define contour center for each ell and precompute sin(center), cos(center)
      if(this_ell<M_PI) center = this_ell+e/2;
      else center = this_ell-e/2;
      sinC = sin(center);
      cosC = cos(center);
      output[i] = center;

      // Accumulate Fourier coefficients
      // NB: we halve the range by symmetry, absorbing factor of 2 into ratio

      ///////////////
      // Separate out j = 0 piece, which is simpler

      // Compute z in real and imaginary parts (zI = 0 here)
      zR = center + radius;

      // Compute e*sin(zR) from precomputed quantities
      tmpsin = sinC*ecosRadius+cosC*esinRadius; // sin(zR)

      // Compute f(z(x)) in real and imaginary parts (fxI = 0)
      fxR = zR - tmpsin - this_ell;

      // Add to array, with factor of 1/2 since an edge
      ft_gx2 = 0.5/fxR;
      ft_gx1 = 0.5/fxR;

      ///////////////
      // Compute for j = 1 to N_points
      // NB: j = jj+1
      for(int jj=0;jj<N_points;jj++){

        // Compute z in real and imaginary parts
        zR = center + radius*exp2R[jj];
        zI = radius*exp2I[jj];

        // Compute f(z(x)) in real and imaginary parts
        // can use precomputed cosh / sinh / cos / sin for this!
        tmpcosh = coshI[jj]; // cosh(zI)
        tmpsinh = sinhI[jj]; // sinh(zI)
        tmpsin = sinC*ecosR[jj]+cosC*esinR[jj]; // e sin(zR)
        tmpcos = cosC*ecosR[jj]-sinC*esinR[jj]; // e cos(zR)

        fxR = zR - tmpsin*tmpcosh-this_ell;
        fxI = zI - tmpcos*tmpsinh;

        // Compute 1/f(z) and append to array
        ftmp = fxR*fxR+fxI*fxI;
        fxR /= ftmp;
        fxI /= ftmp;

        ft_gx2 += (exp4R[jj]*fxR+exp4I[jj]*fxI);
        ft_gx1 += (exp2R[jj]*fxR+exp2I[jj]*fxI);
      }

      ///////////////
      // Separate out j = N_it piece, which is simpler

      // Compute z in real and imaginary parts (zI = 0 here)
      zR = center - radius;

      // Compute sin(zR) from precomputed quantities
      tmpsin = sinC*ecosRadius-cosC*esinRadius; // sin(zR)

      // Compute f(z(x)) in real and imaginary parts (fxI = 0 here)
      fxR = zR - tmpsin-this_ell;

      // Add to sum, with 1/2 factor for edges
      ft_gx2 += 0.5/fxR;
      ft_gx1 += -0.5/fxR;

      ///////////////
      // Compute E(ell)
      output[i] += radius*ft_gx2/ft_gx1;
      }
    }
};

// ========================== Timing Comparison ================

int main(int argc, char *argv[]) {

  // PARAMETERS
  int N_ell = 1000000; // ell array size
  Float e = 0.1; // Eccentricity
  Float tol = 1e-12; // tolerance for error acceptance

  // Print parameters
  printf("N_ell = %d\n",N_ell);
  printf("e = %.2f\n",e);
  printf("tolerance = %.2e\n",tol);

  // Define ell array from a linearly spaced grid of E
  Float* E_exact = new Float[N_ell];
  Float* ell_arr = new Float[N_ell];
  for(int i=0;i<N_ell;i++){
    E_exact[i] = 2.0*M_PI*(i+0.5)/N_ell;
    ell_arr[i] = E_exact[i]-e*sin(E_exact[i]);
  }

  // Create output class to hold methods
  Approximations approx(N_ell, e, ell_arr);

  // Initialize timers
  auto start = high_resolution_clock::now();
  auto stop = high_resolution_clock::now();

  // Output estimates
  Float* E_newton_raphson = new Float[N_ell];
  Float* E_Danby = new Float[N_ell];
  Float* E_series = new Float[N_ell];
  Float* E_contour = new Float[N_ell];

  // Compute Newton-Raphson quadratic estimate
  int N_NR = 0; // Newton-Raphson iterations
  Float err_NR;
  long long int duration_NR;

  // Increase N_NR until we reach tolerance!
  while (N_NR<100){ // max limit!
    start = high_resolution_clock::now(); // starting time
    approx.compute_newton_raphson(N_NR,E_newton_raphson);
    stop = high_resolution_clock::now(); // ending time
    duration_NR = duration_cast<microseconds>(stop - start).count(); // duration

    err_NR = 0;
    for(int i=0;i<N_ell;i++) err_NR += abs(E_exact[i]-E_newton_raphson[i])/N_ell;
    if(err_NR<tol) break;
    N_NR ++;

  }
  printf("Computed Newton-Raphson estimate in %d steps after %.1f ms with mean-error %.2e\n",N_NR,float(duration_NR/1000.),err_NR);

  // Compute Danby quartic estimate
  int N_Danby = 0; // Danby iterations
  Float err_Danby;
  long long int duration_Danby;

  while (N_Danby<100){ // max limit!
    start = high_resolution_clock::now(); // starting time
    approx.compute_danby(N_Danby,E_Danby);
    stop = high_resolution_clock::now(); // ending time
    duration_Danby = duration_cast<microseconds>(stop - start).count(); // duration

    err_Danby = 0;
    for(int i=0;i<N_ell;i++) err_Danby += abs(E_exact[i]-E_Danby[i])/N_ell;
    if(err_Danby<tol) break;
    N_Danby++;

  }
  printf("Computed Danby estimate in %d steps after %.1f ms with mean-error %.2e\n",N_Danby,float(duration_Danby/1000.),err_Danby);

  // Compute series estimate
  int N_series = 0; // Series iterations
  Float err_series;
  long long int duration_series;

  if(e>0.6627434){
    printf("### Series method is non-convergent; skipping!\n");
  }
  else{
    while (N_series<100){ // max limit!
      start = high_resolution_clock::now(); // starting time
      approx.compute_series(N_series,E_series);
      stop = high_resolution_clock::now(); // ending time
      duration_series = duration_cast<microseconds>(stop - start).count(); // duration

      err_series = 0;
      for(int i=0;i<N_ell;i++) err_series += abs(E_exact[i]-E_series[i])/N_ell;
      if(err_series<tol) break;
      N_series++;
    }
    printf("Computed series estimate in %d steps after %.1f ms with mean-error %.2e\n",N_series,float(duration_series/1000.),err_series);
  }

  // Compute contour estimate
  int N_contour = 2; // number of integration steps
  Float err_contour;
  long long int duration_contour;

  while (N_contour<256){ // max limit!
    start = high_resolution_clock::now(); // starting time
    approx.compute_contour(N_contour,E_contour);
    stop = high_resolution_clock::now(); // ending time
    duration_contour = duration_cast<microseconds>(stop - start).count(); // duration

    err_contour = 0;
    for(int i=0;i<N_ell;i++) err_contour += abs(E_exact[i]-E_contour[i])/N_ell;
    if(err_contour<tol) break;
    N_contour+=1;
  }
  printf("Computed contour estimate in %d steps after %.1f ms with mean-error %.2e\n",N_contour,float(duration_contour/1000.),err_contour);

}
