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

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <utility>
#include <vector>

// Could swap between single and double precision here.
using Float = double;

#ifdef __cpp_lib_math_constants
#include <numbers>
constexpr inline Float pi = std::numbers::pi_v<Float>;
#else
constexpr inline Float pi{3.14159265358979323846};
#endif

// ========================== Class to hold basic attributes ================

constexpr inline std::size_t max_contour_iters = 256;

class Approximations {
  private:
    constexpr static std::size_t max_points = max_contour_iters - 2;
    std::size_t N_ell;
    Float e; // eccentricity
    std::vector<Float> ell_arr;
    // storage space for contour integration
    std::vector<Float> exp2R = std::vector<Float>(max_points);
    std::vector<Float> exp2I = std::vector<Float>(max_points);
    std::vector<Float> exp4R = std::vector<Float>(max_points);
    std::vector<Float> exp4I = std::vector<Float>(max_points);
    std::vector<Float> coshI = std::vector<Float>(max_points);
    std::vector<Float> sinhI = std::vector<Float>(max_points);
    std::vector<Float> ecosR = std::vector<Float>(max_points);
    std::vector<Float> esinR = std::vector<Float>(max_points);
  public:
    Approximations(Float _e, std::vector<Float> _ell_arr)
      : e{_e}, ell_arr{std::move(_ell_arr)}
    {
      N_ell = ell_arr.size();
    }

// ========================== Approximation Methods ================

  void compute_newton_raphson(std::size_t N_it, Float *output){
    /// Compute Newton-Raphson method with given number of steps
    // This has quadratic convergence.
    // The initial step is defined as E_0 = ell + sgn(sin(ell))*e*k following Danby (1988)

    using std::cos, std::sin;

    Float k = 0.85;
    Float f_E, fP_E, this_ell, old_E;

    for(std::size_t i=0;i<N_ell;i++){
      this_ell = ell_arr[i];

      // Define initial estimate
      if((sin(this_ell))<0) old_E = this_ell - k*e;
      else old_E = this_ell + k*e;

      // Perform Newton-Raphson estimate
      for(std::size_t j=0;j<N_it;j++) {

        // Compute f(E) and f'(E)
        f_E = old_E - e*sin(old_E)-this_ell;
        fP_E = 1. - e*cos(old_E);

        // Update E
        old_E -= f_E/fP_E;
      }

      // Add to array
      output[i] = old_E;
    }
  }

  void compute_danby(std::size_t N_it, Float *output){
    /// Compute Danby (1988) fourth-order root-finding method with given number of steps
    // This has quartic convergence.
    // The initial step is defined as E_0 = ell + sgn(sin(ell))*e*k following Danby (1988)

    using std::cos, std::sin;

    Float k = 0.85;
    Float f_E, fP_E, fPP_E, fPPP_E, this_ell, old_E, delta_i1, delta_i2, delta_i3, esinE, ecosE;

    for(std::size_t i=0;i<N_ell;i++){
      this_ell = ell_arr[i];

      // Define initial estimate
      if((sin(this_ell))<0) old_E = this_ell - k*e;
      else old_E = this_ell + k*e;

      // Perform Newton-Raphson estimate
      for(std::size_t j=0;j<N_it;j++) {

        // Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
        esinE = e*sin(old_E);
        ecosE = e*cos(old_E);
        f_E = old_E - esinE-this_ell;
        fP_E = 1. - ecosE;
        fPP_E = esinE;
        fPPP_E = ecosE;

        delta_i1 = -f_E/fP_E;
        delta_i2 = -f_E/(fP_E+1./2.*delta_i1*fPP_E);
        delta_i3 = -f_E/(fP_E+1./2.*delta_i2*fPP_E+1./6.*fPPP_E*delta_i2*delta_i2);

        // Update E
        old_E += delta_i3;
      }

      // Add to array
      output[i] = old_E;
    }
  }

  void compute_series(std::size_t N_it, Float *output){
    // Solve Kepler's equation via the series method described in Murray & Dermot.
    // We use a specifies maximum number of iterations and precompute the coupling coefficients.

    using std::cyl_bessel_j, std::sin;

    Float coeff;

    // Take initial guess
    for(std::size_t i=0; i<N_ell;i++){
      output[i] = ell_arr[i];
    }

    // Iterate over number of iterations
    for(std::size_t s=1; s<=N_it; s++){
        coeff = 2*cyl_bessel_j(s,s*e)/s; // define coefficient

        for(std::size_t i=0;i<N_ell;i++){
          output[i] += coeff*sin(s*ell_arr[i]);
        }
    }
  }

void compute_contour(std::size_t N_it, Float *output){
    // Solve Kepler's equation via the contour integration method of Philcox et al. (2021)
    // This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.
    // N_it specifies the number of grid-points.

    using std::cos, std::cosh, std::sin, std::sinh;

    Float ft_gx2, ft_gx1, this_ell, freq, zR, zI, cosC, sinC, esinRadius, ecosRadius, center;
    Float fxR, fxI, ftmp, tmpcosh, tmpsinh, tmpcos, tmpsin;

    // Define sampling points (actually use one more than this)
    std::size_t N_points = N_it-2;
    std::size_t N_fft = (N_it-1)*2;

    // Define contour radius
    Float radius = e/2;

    // Generate e^{ikx} sampling points and precompute real and imaginary parts
    Float cf, sf;
    for(std::size_t jj=0;jj<N_points;jj++){
      // NB: j = jj+1
      freq = 2.0*pi*(jj+1)/N_fft;
      cf = cos(freq);
      sf = sin(freq);
      exp2R[jj] = cf;
      exp2I[jj] = sf;
      exp4R[jj] = cf*cf-sf*sf;
      exp4I[jj] = 2.0*cf*sf;
      coshI[jj] = cosh(radius*exp2I[jj]);
      sinhI[jj] = sinh(radius*exp2I[jj]);
      ecosR[jj] = e*cos(radius*exp2R[jj]);
      esinR[jj] = e*sin(radius*exp2R[jj]);
    }

    // Precompute e sin(e/2) and e cos(e/2)
    esinRadius = e*sin(radius);
    ecosRadius = e*cos(radius);

    // Iterate over array of mean anomalies
    for(std::size_t i=0;i<N_ell;i++){
      this_ell = ell_arr[i];

      // Define contour center for each ell and precompute sin(center), cos(center)
      if(this_ell<pi) center = this_ell+e/2;
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
      for(std::size_t jj=0;jj<N_points;jj++){

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

  using namespace std::chrono;
  using std::abs, std::sin;

  // PARAMETERS
  constexpr std::size_t N_ell = 1000000; // ell array size
  Float e = 0.5; // Eccentricity
  Float tol = 1e-12; // tolerance for error acceptance

  // Print parameters
  std::printf("N_ell = %zu\n",N_ell);
  std::printf("e = %.2f\n",e);
  std::printf("tolerance = %.2e\n",tol);

  // Define ell array from a linearly spaced grid of E
  std::vector<Float> E_exact(N_ell);
  std::vector<Float> ell_arr(N_ell);
  for(std::size_t i=0;i<N_ell;i++){
    E_exact[i] = 2.0*pi*(i+0.5)/N_ell;
    ell_arr[i] = E_exact[i]-e*sin(E_exact[i]);
  }

  // Create output class to hold methods
  Approximations approx(e, ell_arr);

  // Initialize timers
  auto start = high_resolution_clock::now();
  auto stop = high_resolution_clock::now();

  // Output estimates
  std::vector<Float> E_newton_raphson(N_ell);
  std::vector<Float> E_Danby(N_ell);
  std::vector<Float> E_series(N_ell);
  std::vector<Float> E_contour(N_ell);

  // Compute Newton-Raphson quadratic estimate
  std::size_t N_NR = 0; // Newton-Raphson iterations
  Float err_NR;
  long long int duration_NR;

  // Increase N_NR until we reach tolerance!
  while (N_NR<100){ // max limit!
    start = high_resolution_clock::now(); // starting time
    approx.compute_newton_raphson(N_NR, E_newton_raphson.data());
    stop = high_resolution_clock::now(); // ending time
    duration_NR = duration_cast<microseconds>(stop - start).count(); // duration

    err_NR = 0;
    for(std::size_t i=0;i<N_ell;i++) err_NR += abs(E_exact[i]-E_newton_raphson[i])/N_ell;
    if(err_NR<tol) break;
    N_NR ++;

  }
  std::printf("Computed Newton-Raphson estimate in %zu steps after %.1f ms with mean-error %.2e\n",N_NR,float(duration_NR/1000.),err_NR);

  // Compute Danby quartic estimate
  std::size_t N_Danby = 0; // Danby iterations
  Float err_Danby;
  long long int duration_Danby;

  while (N_Danby<100){ // max limit!
    start = high_resolution_clock::now(); // starting time
    approx.compute_danby(N_Danby, E_Danby.data());
    stop = high_resolution_clock::now(); // ending time
    duration_Danby = duration_cast<microseconds>(stop - start).count(); // duration

    err_Danby = 0;
    for(std::size_t i=0;i<N_ell;i++) err_Danby += abs(E_exact[i]-E_Danby[i])/N_ell;
    if(err_Danby<tol) break;
    N_Danby++;

  }
  std::printf("Computed Danby estimate in %zu steps after %.1f ms with mean-error %.2e\n",N_Danby,float(duration_Danby/1000.),err_Danby);

  // Compute series estimate
  std::size_t N_series = 0; // Series iterations
  Float err_series;
  long long int duration_series;

  if(e>0.6627434){
    std::printf("### Series method is non-convergent; skipping!\n");
  }
  else{
    while (N_series<100){ // max limit!
      start = high_resolution_clock::now(); // starting time
      approx.compute_series(N_series, E_series.data());
      stop = high_resolution_clock::now(); // ending time
      duration_series = duration_cast<microseconds>(stop - start).count(); // duration

      err_series = 0;
      for(std::size_t i=0;i<N_ell;i++) err_series += abs(E_exact[i]-E_series[i])/N_ell;
      if(err_series<tol) break;
      N_series++;
    }
    std::printf("Computed series estimate in %zu steps after %.1f ms with mean-error %.2e\n",N_series,float(duration_series/1000.),err_series);
  }

  // Compute contour estimate
  std::size_t N_contour = 2; // number of integration steps
  Float err_contour;
  long long int duration_contour;

  while (N_contour<max_contour_iters){ // max limit!
    start = high_resolution_clock::now(); // starting time
    approx.compute_contour(N_contour,E_contour.data());
    stop = high_resolution_clock::now(); // ending time
    duration_contour = duration_cast<microseconds>(stop - start).count(); // duration

    err_contour = 0;
    for(std::size_t i=0;i<N_ell;i++) err_contour += abs(E_exact[i]-E_contour[i])/N_ell;
    if(err_contour<tol) break;
    N_contour+=1;
  }
  std::printf("Computed contour estimate in %zu steps after %.1f ms with mean-error %.2e\n",N_contour,float(duration_contour/1000.),err_contour);

}
