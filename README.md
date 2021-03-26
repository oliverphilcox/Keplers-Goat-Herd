# Kepler's Goat Herd

C++ code for solving Kepler's equation using contour integration, following Philcox et al. (2021, in prep.). This uses a method originally proposed by Ullisch (2020) to solve the "geometric goat problem".

The code contains implementations of a variety of solution methods:
- Newton-Raphson: The quadratic Newton-Raphson root finder.
- Danby: The quartic root finder described in Danby (1988).
- Series: An elliptical series method, as described in Murray & Dermott.
- Contour: A new method based on contour integration.

Given an array of mean anomalies, an eccentricity and a desired precision, the code will estimate the eccentricity using each method. The accuracy of each approach is increased until the desired precision is reached, relative to an (overconverged) Danby estimate with 100 steps, and timing is performed using the C++ `chrono` package.

To compile the code using g++, simply run ```rm kepler; g++ -o kepler keplers_goat_herd.cpp -std=c++17 -ffast-math -Wall -O3; ./kepler```.

**Authors**:
- Oliver Philcox (Princeton, [ohep2@cantab.ac.uk](mailto:ohep2@cantab.ac.uk))
