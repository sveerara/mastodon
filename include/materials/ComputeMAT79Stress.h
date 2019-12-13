/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEMAT79STRESS_H
#define COMPUTEMAT79STRESS_H

#include "ComputeFiniteStrainElasticStress.h"


class ComputeMAT79Stress : public ComputeFiniteStrainElasticStress
{
public:
  ComputeMAT79Stress(const InputParameters & parameters);

protected:

  virtual void computeQpStress();


  virtual void computeStress();

  virtual void initQpStatefulProperties();

  virtual void computeNonmasing();

  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _total_strain;
  const MaterialProperty<RankTwoTensor> & _total_strain_old;

  std::vector<std::string> _base_models;

  std::vector<MaterialProperty<RankTwoTensor> *>  _stress_model;
  std::vector<MaterialProperty<RankTwoTensor> *>  _stress_model_old;

  std::vector<Real> _backbone_stress;
  std::vector<Real> _backbone_strain;
  std::vector<Real> _sudo_backbone_stress;
  std::vector<Real> _sudo_backbone_strain;
  std::vector<Real> _mapped_backbone_stress;
  std::vector<Real> _mapped_backbone_strain;
  Real _reduction_factor;
  Real _Gm;
  Real _b_exp;
  Real _Poissons;
  Real _a0;
  Real _a1;
  Real _a2;
  Real _p0;
  Real _p1;
  Real _p2;
  Real _p3;
  Real _K0;
  Real _G0;
  Real _p_ref;


  std::vector<Real> _youngs;
  std::vector<Real> _yield_stress;
  // Used in stress calculation
  RankTwoTensor _stress_new;
  RankTwoTensor _individual_stress_increment;
  RankTwoTensor _deviatoric_trial_stress;
  Real _strength_pressure_correction;
  Real _stiffness_pressure_correction;
  Real _tangent_modulus;
  private:
  	Real _initial_reversal_check;


  MaterialProperty<Real> &_reversal_check;
  const MaterialProperty<Real> &_reversal_check_old;

  private:
  	Real _initial_strain_memory;


  MaterialProperty<Real> &_strain_memory;
  const MaterialProperty<Real> &_strain_memory_old;
  private:
  	unsigned int _number_n;

};

#endif //COMPUTEMAT79STRESS_H
