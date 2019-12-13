/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FrequencyIndependentDamper.h"
#include "libmesh/quadrature.h"
#include "MooseMesh.h"
#include "Assembly.h"
#include "ElasticityTensorTools.h"

template <>
InputParameters
validParams<FrequencyIndependentDamper>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredParam<Real>("damping_coefficient", "Damping coefficient for frequency independent damping model.");
  params.addRequiredParam<Real>("beta", "Newmark time integration parameter.");
  params.set<MultiMooseEnum>("execute_on") = "INITIAL";
  params.suppressParameter<MultiMooseEnum>("execute_on");
  return params;
}

FrequencyIndependentDamper::FrequencyIndependentDamper(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _damping_coeff(getParam<Real>("damping_coefficient")),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _density(getMaterialPropertyByName<Real>("density")),
    _beta(getParam<Real>("beta")),
    _second(_mesh.hasSecondOrderElements()),
    _fe_type(_second ? SECOND : FIRST, LAGRANGE),
    _phi(_assembly.fePhi(_fe_type)),
    _test(_assembly.fePhi(_fe_type)),
    _grad_phi(_assembly.feGradPhi(_fe_type)),
    _grad_test(_assembly.feGradPhi(_fe_type))
{
}

void
FrequencyIndependentDamper::initialize()
{
  ColumnMajorMatrix C;
  C.reshape(0,0);
  ConstElemRange & elem_range = *_mesh.getActiveLocalElementRange();
  for (const auto & elem : elem_range)
    _C_map[elem->id()] = C;
}
void
FrequencyIndependentDamper::execute()
{
  if (_C_map[_current_elem->id()].n() == 0 && _C_map[_current_elem->id()].m() == 0)
  {
    _C_map[_current_elem->id()].reshape(_test.size() * _mesh.dimension(), _test.size() * _mesh.dimension());
    _C_map[_current_elem->id()].zero();

    ColumnMajorMatrix K(_test.size() * _mesh.dimension(), _phi.size() * _mesh.dimension());
    K.zero();
    for (unsigned int component = 0; component < _mesh.dimension(); ++component)
      for (unsigned int coupled_component = 0; coupled_component < _mesh.dimension(); ++coupled_component)
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
              K(i + component * _test.size(), j + coupled_component * _test.size()) += _JxW[qp] * ElasticityTensorTools::elasticJacobian(_Jacobian_mult[qp],
                                                       component,
                                                       coupled_component,
                                                       _grad_test[i][qp],
                                                       _grad_phi[j][qp]);

    ColumnMajorMatrix M(_test.size() * _mesh.dimension(), _phi.size() * _mesh.dimension());
    M.zero();
    for (unsigned int component = 0; component < _mesh.dimension(); ++component)
      for (unsigned int coupled_component = 0; coupled_component < _mesh.dimension(); ++coupled_component)
        for (unsigned int i = 0; i < _test.size(); ++i)
          for (unsigned int j = 0; j < _phi.size(); ++j)
            for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
            {
              if (component == coupled_component)
                M(i + component * _test.size(), j + coupled_component * _test.size()) += _JxW[qp] * _test[i][qp] * _density[qp] * _phi[j][qp];
              else
                M(i + component * _test.size(), j + coupled_component * _test.size()) += 0.0;
            }


    if (M.n() > 0 && M.m() > 0)
    {
      ColumnMajorMatrix B(_test.size() * _mesh.dimension(), _test.size() * _mesh.dimension());
      M.inverse(B);
      B = (B + B.transpose()) * 0.5;

      ColumnMajorMatrix A = B * K;
      A = (A + A.transpose()) * 0.5;

      // calculate eigenvalues and eigenvectors of Inv(M) * K
      ColumnMajorMatrix eval, evec;
      A.eigen(eval, evec);
      eval.print();
      ColumnMajorMatrix sqeval(_test.size() * _mesh.dimension(), _test.size() * _mesh.dimension());

      for (unsigned int i = 0; i < _test.size() * _mesh.dimension(); ++i)
        sqeval(i, i) = std::sqrt(std::abs(eval(i)));

      ColumnMajorMatrix inv_evec;
      evec.inverse(inv_evec);

      _C_map[_current_elem->id()] = (M * evec * sqeval * inv_evec) * 2.0 * _damping_coeff;
    }
  }
}

Real
FrequencyIndependentDamper::getResidual(const unsigned int i, const std::vector<Real> & vel) const
{
  Real res = 0.0;

  if (_C_map[_current_elem->id()].n() > 0 && _C_map[_current_elem->id()].m() > 0)
  {
    for (unsigned int j = 0; j < _phi.size() * _mesh.dimension(); ++j)
      res += _C_map[_current_elem->id()](i, j) * vel[j];
  }
  else
  {
    printf("Shouldn't be here!\n");
  }
  return res;
}

Real
FrequencyIndependentDamper::getValue(const unsigned int i, const unsigned int j) const
{
  Real value = 0.0;

  if (_C_map[_current_elem->id()].n() > 0 && _C_map[_current_elem->id()].m() > 0)
    value = _C_map[_current_elem->id()](i,j);

  return value;
}
