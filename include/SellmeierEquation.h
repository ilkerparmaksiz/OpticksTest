// ----------------------------------------------------------------------------
// nexus | SellmeierEquation.h
//
// The Sellmeier equation is an empirical relationship between refractive
// index and wavelength for a dielectric transparent medium.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef SELLMEIER_EQUATION_H
#define SELLMEIER_EQUATION_H

#include "math.h"
namespace nexus {
  class SellmeierEquation
  {
  public:
    SellmeierEquation(double* B, double* C);
    ~SellmeierEquation();

    double RefractiveIndex(double wavelength);

  private:
    double B_[3];
    double C_[3];
  };

  // Inline definitions ///////////////////////////////////

  inline SellmeierEquation::SellmeierEquation(double* B, double* C)
  {
    for (unsigned int i=0; i<3; i++) {
      B_[i] = B[i]; C_[i] = C[i];
    }
  }

  inline SellmeierEquation::~SellmeierEquation() {}

  inline double SellmeierEquation::RefractiveIndex(double wavelength)
  {
    double n2 = 1.;
    double wl2 = wavelength * wavelength;

    for (unsigned int i=0; i<3; i++)
      n2 += (B_[i] * wl2) / (wl2 - C_[i]);

    return sqrt(n2);
  }
}
#endif
