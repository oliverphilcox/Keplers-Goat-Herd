import numpy as np, time

def compute_contour(ell_array, eccentricity, N_it):
    """Solve Kepler's equation, E - e sin E = ell, via the contour integration method of Philcox et al. (2021)
    This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.

    Args:
        ell_array (np.ndarray): Array of mean anomalies, ell, in the range (0,2 pi).
        eccentricity (float): Eccentricity. Must be in the range 0<e<1.
        N_it (float): Number of grid-points.

    Returns:
        np.ndarray: Array of eccentric anomalies, E.
    """

    # Check inputs
    if eccentricity<=0.:
        raise Exception("Eccentricity must be greater than zero!")
    elif eccentricity>=1:
        raise Exception("Eccentricity must be less than unity!")
    if np.max(ell_array)>2.*np.pi:
        raise Exception("Mean anomaly should be in the range (0, 2 pi)")
    if np.min(ell_array)<0:
        raise Exception("Mean anomaly should be in the range (0, 2 pi)")
    if N_it<2:
        raise Exception("Need at least two sampling points!")

    # Define sampling points
    N_points = N_it - 2
    N_fft = (N_it-1)*2

    # Define contour radius
    radius = eccentricity/2

    # Generate e^{ikx} sampling points and precompute real and imaginary parts
    j_arr = np.arange(N_points)
    freq = (2*np.pi*(j_arr+1.)/N_fft)[:,np.newaxis]
    exp2R = np.cos(freq)
    exp2I = np.sin(freq)
    ecosR= eccentricity*np.cos(radius*exp2R)
    esinR = eccentricity*np.sin(radius*exp2R)
    exp4R = exp2R*exp2R-exp2I*exp2I
    exp4I = 2.*exp2R*exp2I
    coshI = np.cosh(radius*exp2I)
    sinhI = np.sinh(radius*exp2I)

    # Precompute e sin(e/2) and e cos(e/2)
    esinRadius = eccentricity*np.sin(radius);
    ecosRadius = eccentricity*np.cos(radius);

    # Define contour center for each ell and precompute sin(center), cos(center)
    filt = ell_array<np.pi
    center = ell_array-eccentricity/2.
    center[filt] += eccentricity
    sinC = np.sin(center)
    cosC = np.cos(center)
    output = center

    ## Accumulate Fourier coefficients
    # NB: we halve the integration range by symmetry, absorbing factor of 2 into ratio

    ## Separate out j = 0 piece, which is simpler

    # Compute z in real and imaginary parts (zI = 0 here)
    zR = center + radius

    # Compute e*sin(zR) from precomputed quantities
    tmpsin = sinC*ecosRadius+cosC*esinRadius

    # Compute f(z(x)) in real and imaginary parts (fxI = 0)
    fxR = zR - tmpsin - ell_array

     # Add to arrays, with factor of 1/2 since an edge
    ft_gx2 = 0.5/fxR
    ft_gx1 = 0.5/fxR

    ## Compute j = 1 to N_points pieces

    # Compute z in real and imaginary parts
    zR = center + radius*exp2R
    zI = radius*exp2I

    # Compute f(z(x)) in real and imaginary parts
    # can use precomputed cosh / sinh / cos / sin for this!
    tmpsin = sinC*ecosR+cosC*esinR # e sin(zR)
    tmpcos = cosC*ecosR-sinC*esinR # e cos(zR)

    fxR = zR - tmpsin*coshI-ell_array
    fxI = zI - tmpcos*sinhI

    # Compute 1/f(z) and append to array
    ftmp = fxR*fxR+fxI*fxI;
    fxR /= ftmp;
    fxI /= ftmp;

    ft_gx2 += np.sum(exp4R*fxR+exp4I*fxI,axis=0)
    ft_gx1 += np.sum(exp2R*fxR+exp2I*fxI,axis=0)

    ## Separate out j = N_it piece, which is simpler

    # Compute z in real and imaginary parts (zI = 0 here)
    zR = center - radius

    # Compute sin(zR) from precomputed quantities
    tmpsin = sinC*ecosRadius-cosC*esinRadius

    # Compute f(z(x)) in real and imaginary parts (fxI = 0 here)
    fxR = zR - tmpsin-ell_array

    # Add to sum, with 1/2 factor for edges
    ft_gx2 += 0.5/fxR;
    ft_gx1 += -0.5/fxR;

    ### Compute and return the solution E(ell,e)
    output += radius*ft_gx2/ft_gx1;

    return output

if __name__=="__main__":
    """Test the Python function above with a simple example"""

    # Parameters
    N_ell = 1000000
    eccentricity = 0.5
    N_it = 10

    print("\n##### PARAMETERS #####")
    print("# N_ell = %d"%N_ell)
    print("# Eccentricity = %.2f"%eccentricity)
    print("# Iterations: %d"%N_it)
    print("######################")

    # Create ell array from E
    E_true = (2.0*np.pi*(np.arange(N_ell)+0.5))/N_ell
    ell_input = E_true - eccentricity*np.sin(E_true)

    # Time the function
    init = time.time()
    E_out = compute_contour(ell_input,eccentricity,N_it)
    runtime = time.time()-init

    print("\nEstimation complete after %.1f millseconds, achieving mean error %.2e.\n"%(runtime*1000.,np.mean(np.abs(E_out-E_true))))
