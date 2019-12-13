/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "Kernel.h"
#include "Material.h"
#include "FrequencyIndependentDamper.h"

// Forward Declarations
class DampingForce;

template <>
InputParameters validParams<DampingForce>();

class DampingForce : public Kernel
{
public:
  DampingForce(const InputParameters & parameters);

  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

protected:
  virtual void computeResidual() override;
  virtual Real computeQpResidual() override {return 0.0; };

private:
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  std::vector<unsigned int> _disp_num;
  std::vector<unsigned int> _vel_num;
  std::vector<unsigned int> _accel_num;
  const Real _beta;
  const Real _gamma;
  const FrequencyIndependentDamper * const _damping_matrix_calculator;
  const unsigned int _component;
};

#endif // DAMPINGFORCE_H
