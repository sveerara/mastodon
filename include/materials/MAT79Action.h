/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MAT79ACTION_H
#define MAT79ACTION_H

#include "Action.h"

class MAT79Action: public Action
{
public:
  MAT79Action(const InputParameters & params);

  virtual void act();

private:
  const std::string _data_file;
  const Real _Poissons;
  const Real _b_exp;
  const Real _a0;
  const Real _a1;
  const Real _a2;
  const Real _p0;
  const Real _p1;
  const Real _p2;
  const Real _p3;
  const Real _K0;
  const Real _G0;
  const Real _p_ref;
  const Real _OCR;
  const Real _PI;
  const Real _fit_2_theta1;
  const Real _fit_2_theta2;
  const Real _fit_2_theta3;
  const Real _fit_2_theta4;
  const Real _fit_2_theta5;
  const Real _taumax;
  const int  _numberofpoints;
  const int  _soiltype;
  const std::vector<SubdomainName> _block;
  unsigned int _number;
  
  bool parseNextLineReals( std::ifstream & ifs, std::vector<Real> & myvec);
  void parseColumns( std::vector<Real> & x, std::vector<Real> & y);
};

template<>
InputParameters validParams<MAT79Action>();

#endif //MAT79ACTION_H
