# Kepler's Goat Herd

Code for solving Kepler's equation using contour integration, following Philcox et al. (2021, [arXiv](https://arxiv.org/abs/2103.15829)). This uses a method originally proposed by Ullisch (2020) to solve the "geometric goat problem".

The C++ code contains implementations of a variety of solution methods:
- Newton-Raphson: The quadratic Newton-Raphson root finder.
- Danby: The quartic root finder described in Danby (1988).
- Series: An elliptical series method, as described in Murray & Dermott.
- Contour: A new method based on contour integration.

Given an array of mean anomalies, an eccentricity and a desired precision, the code will estimate the eccentric anomaly using each method. The accuracy of each approach is increased until the desired precision is reached, and timing is performed using the C++ `chrono` package.

To compile the code using g++, simply run ```g++ -o kepler keplers_goat_herd.cpp -std=c++17 -ffast-math -Wall -O3```.  The code can be run using ```./kepler```. The individual functions, e.g. ```compute_contour``` can also be used outside of this script, given an input array of mean anomalies and an eccentricity.

For non-gcc compilers that do not support the GNU variable length arrays extension, the file ```keplers_goat_herd_std17.cpp``` can be used instead.

We also provide a pure numpy version of the contour integration function in ```keplers_goat_herd.py```. This is around 9 times slower than the C++ code. Python bindings for the C++ code can be added if these would be of use.

**Authors**:
- Oliver Philcox (Princeton, [ohep2@cantab.ac.uk](mailto:ohep2@cantab.ac.uk))
