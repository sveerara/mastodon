/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DampingForce.h"
#include "SubProblem.h"
#include "NonlinearSystem.h"
#include "AuxiliarySystem.h"
#include "MooseVariable.h"
#include "Assembly.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<DampingForce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for the frequency independent damping.");
  params.set<bool>("use_displaced_mesh") = true;
  params.addRequiredCoupledVar("velocities", "velocity variable");
  params.addRequiredCoupledVar("displacements", "displacement variable");
  params.addRequiredCoupledVar("accelerations", "aceeleration variable");
  params.addRequiredParam<Real>("beta", "beta parameter for Newmark Time integration");
  params.addRequiredParam<Real>("gamma", "gamma parameter for Newmark Time integration");
  params.addRequiredParam<UserObjectName>("damping_matrix_userobject", "User object that contains the damping matrix");
  params.addRequiredParam<unsigned int>("component", "The component for which damping force is being calculated.");
  return params;
}

DampingForce::DampingForce(const InputParameters & parameters)
  : Kernel(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp_var(3),
    _disp_num(3),
    _vel_num(3),
    _accel_num(3),
    _beta(getParam<Real>("beta")),
    _gamma(getParam<Real>("gamma")),
    _damping_matrix_calculator(&getUserObject<FrequencyIndependentDamper>("damping_matrix_userobject")),
    _component(getParam<unsigned int>("component"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp_var[i] = coupled("displacements", i);

    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number();

    MooseVariable * vel_variable = getVar("velocities", i);
    _vel_num[i] = vel_variable->number();

    MooseVariable * accel_variable = getVar("accelerations", i);
    _accel_num[i] = accel_variable->number();
  }
}

void
DampingForce::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  mooseAssert(re.size() == _test.size(), "Residual vector should have same size as the number of nodes in the element.");
  _local_re.resize(re.size());
  _local_re.zero();

  if (_dt != 0)
  {
    // fetch the two end nodes for _current_elem
    std::vector<Node *> node;
    for (unsigned int i = 0; i < _test.size(); ++i)
      node.push_back(_current_elem->get_node(i));

    // Fetch the solution for the nodes of the element
    NonlinearSystemBase & nonlinear_sys = _fe_problem.getNonlinearSystemBase();
    const NumericVector<Number> & sol = *nonlinear_sys.currentSolution();
    const NumericVector<Number> & sol_old = nonlinear_sys.solutionOld();

    AuxiliarySystem & aux = _fe_problem.getAuxiliarySystem();
    const NumericVector<Number> & aux_sol_old = aux.solutionOld();
    std::vector<Real> vel(_test.size() * _mesh.dimension());

    for (unsigned int j = 0; j < _mesh.dimension(); ++j)
    {
      for (_i = 0; _i < _test.size(); ++_i)
      {
        Real disp = sol(node[_i]->dof_number(nonlinear_sys.number(), _disp_num[j], 0));
        Real disp_old= sol_old(node[_i]->dof_number(nonlinear_sys.number(), _disp_num[j], 0));

        Real vel_old = aux_sol_old(node[_i]->dof_number(aux.number(), _vel_num[j], 0));

        Real accel_old = aux_sol_old(node[_i]->dof_number(aux.number(), _accel_num[j], 0));

        Real accel = 1. / _beta * (((disp - disp_old) / (_dt * _dt)) - vel_old / _dt -
                                   accel_old * (0.5 - _beta));

        vel[_i + j * _test.size()] = vel_old + (_dt * (1 - _gamma)) * accel_old + _gamma * _dt * accel;
      }
    }


    for (_i = 0; _i < _test.size(); ++_i)
      _local_re(_i) = _damping_matrix_calculator->getResidual(_i + _component * _test.size(), vel);

  }

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _save_in.size(); ++i)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}

void
DampingForce::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  if (_dt > 0)
  {
    for (unsigned int i = 0; i < _test.size(); ++i)
      for (unsigned int j = 0; j < _phi.size(); ++j)
        _local_ke(i, j) = _damping_matrix_calculator->getValue(i +  _component * _test.size(), j + _component * _phi.size());

    _local_ke *= _gamma / _beta / _dt;
  }

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i = 0; i < rows; ++i)
      diag(i) = _local_ke(i, i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}

void
DampingForce::computeOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _var.number())
    computeJacobian();
  else
  {
    unsigned int coupled_component = 0;
    bool disp_coupled = false;

    for (unsigned int i = 0; i < _ndisp; ++i)
      if (jvar == _disp_var[i])
      {
        coupled_component = i;
        disp_coupled = true;
        break;
      }

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);

    if (disp_coupled && _dt > 0)
    {
      for (unsigned int i = 0; i < _test.size(); ++i)
        for (unsigned int j = 0; j < _phi.size(); ++j)
          ke(i, j) += _damping_matrix_calculator->getValue(i + _component * _test.size(), j + coupled_component * _phi.size()) * _gamma / _beta / _dt;
    }
  }
}
