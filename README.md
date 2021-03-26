# Keplers-Goat-Herd

Solve Kepler's equation using contour integration.

The C++ packagwith a variety of numerical techniques for an array of mean anomalies
//
// Methods implemented are:
//   - Newton-Raphson: The quadratic Newton-Raphson root finder.
//   - Danby: The quartic root finder described in Danby (1988).
//   - Series: An elliptic series method, as described in Murray & Dermott.
//   - Contour: A method based on contour integration, described in Philcox et al. (2021).
//
// Given an array of mean anomalies, an eccentricity and a desired precision, the code will estimate the eccentricity using each method.
// The accuracy of each approach is increased until the desired precision is reached, relative to an (overconverged) Danby estimate with 100 steps.
// Timing is performed using the `chrono` package.

**Authors**:
- Oliver Philcox (Princeton, [ohep2@cantab.ac.uk](mailto:ohep2@cantab.ac.uk))
