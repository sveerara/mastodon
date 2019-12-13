/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FREQUENCYINDEPENDENTDAMPER_H
#define FREQUENCYINDEPENDENTDAMPER_H

#include "ElementUserObject.h"
#include "RankFourTensor.h"
#include "ColumnMajorMatrix.h"

class FrequencyIndependentDamper;

template <>
InputParameters validParams<FrequencyIndependentDamper>();

/**
 * Works on top of NodalNormalsPreprocessor
 */
class FrequencyIndependentDamper : public ElementUserObject
{
public:
  FrequencyIndependentDamper(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & /*uo*/) override{};
  virtual void finalize() override{};

  Real getResidual(const unsigned int i, const std::vector<Real> & vel) const;
  Real getValue(const unsigned int i, const unsigned int j) const;
  mutable std::map<dof_id_type, ColumnMajorMatrix> _C_map;
protected:
  const std::string _base_name;
  const Real _damping_coeff;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;
  const MaterialProperty<Real> & _density;
  const Real _beta;
  const bool _second;
  FEType _fe_type;
  const VariablePhiValue & _phi;
  const VariableTestValue & _test;
  const VariablePhiGradient & _grad_phi;
  const VariableTestGradient & _grad_test;
};

#endif /* FREQUENCYINDEPENDENTDAMPER_H */
