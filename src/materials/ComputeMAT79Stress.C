/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeMAT79Stress.h"
#include"ColumnMajorMatrix.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<ComputeMAT79Stress>()
{
  InputParameters params = validParams<ComputeFiniteStrainElasticStress>();
  params.addClassDescription("Compute stress using a set of elastic perfectly plastic stress-strain curves");
  params.addRequiredParam<Real>("b_exp", "The exponent for pressure dependent stiffness");
  params.addRequiredParam<Real>("poissons", "The poissons ratio for the material");
  params.addRequiredParam<Real>("a0", "pressure dependent strength constant");
  params.addRequiredParam<Real>("a1", "pressure dependent strength constant");
  params.addRequiredParam<Real>("a2", "pressure dependent strength constant");
  params.addRequiredParam<Real>("p0", "tension cutoff");
  params.addRequiredParam<Real>("p1", "MRDF type nonmasing parameter");
  params.addRequiredParam<Real>("p2", "MRDF type nonmasing parameter");
  params.addRequiredParam<Real>("p3", "MRDF type nonmasing parameter");
  params.addRequiredParam<Real>("K0", "Initial bulk modulus");
  params.addRequiredParam<Real>("G0", "Initial Shear Modulus");
  params.addRequiredParam<Real>("p_ref", "Reference Pressure");
  params.addRequiredParam<std::vector<Real> >("backbone_stress", "Backbone stress values");
  params.addRequiredParam<std::vector<Real> >("backbone_strain", "Backbone strain values");
  params.addRequiredParam<std::vector<std::string> >("base_models", "Base name for each elastic perfectly plastic model.");
  params.addParam<Real>("initial_reversal_check",0.0,"The initial value for reversal check");
  params.addParam<Real>("initial_strain_memory",0.0,"The initial value for strain memory");
  return params;
}

ComputeMAT79Stress::ComputeMAT79Stress(const InputParameters & parameters) :
    ComputeFiniteStrainElasticStress(parameters),
    _strain_increment(getMaterialProperty<RankTwoTensor>(_base_name + "strain_increment")),
    _total_strain(getMaterialProperty<RankTwoTensor>(_base_name + "total_strain")),
    _total_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "total_strain")),
    _base_models(getParam<std::vector<std::string> >("base_models")),
    _stress_model(_base_models.size()),
    _stress_model_old(_base_models.size()),
    _backbone_stress(getParam<std::vector<Real> >("backbone_stress")),
    _backbone_strain(getParam<std::vector<Real> >("backbone_strain")),
    _sudo_backbone_stress(_backbone_stress.size()),
    _sudo_backbone_strain(_backbone_stress.size()),
    _mapped_backbone_stress(_backbone_stress.size()),
    _mapped_backbone_strain(_backbone_stress.size()),
    _reduction_factor(1.0),
    _Gm(0.0),
    _b_exp(getParam<Real>("b_exp")),
	_Poissons(getParam<Real>("poissons")),
    _a0(getParam<Real>("a0")),
    _a1(getParam<Real>("a1")),
    _a2(getParam<Real>("a2")),
    _p0(getParam<Real>("p0")),
    _p1(getParam<Real>("p1")),
    _p2(getParam<Real>("p2")),
    _p3(getParam<Real>("p3")),
    _K0(getParam<Real>("K0")),
    _G0(getParam<Real>("G0")),
    _p_ref(getParam<Real>("p_ref")),
    _youngs(_backbone_stress.size()),
    _yield_stress(_backbone_stress.size()),
    _strength_pressure_correction(0.0),
    _stiffness_pressure_correction(0.0),
    _tangent_modulus(0.0),
	_initial_reversal_check(getParam<Real>("initial_reversal_check")),
	_reversal_check(declareProperty<Real>("_reversal_check")),
	_reversal_check_old(declarePropertyOld<Real>("_reversal_check")),
	_initial_strain_memory(getParam<Real>("initial_strain_memory")),
	_strain_memory(declareProperty<Real>("_strain_memory")),
	_strain_memory_old(declarePropertyOld<Real>("_strain_memory"))
{
   _number_n = 1;
   for (unsigned int i = 0; i < _base_models.size(); i++)
   {
     _stress_model[i] = &declareProperty<RankTwoTensor>(_base_models[i] + "_stress_model");
     _stress_model_old[i] = &declarePropertyOld<RankTwoTensor>(_base_models[i] + "_stress_model");
   }
}
void
ComputeMAT79Stress::computeNonmasing()
{

//   if ((_youngs.size() != _yield_stress.size()) || (_youngs.size() != _base_models.size()))
//     mooseError("Youngs modulus, yield stress and base_models should be of the same size");


   _stress_new.zero();
   _individual_stress_increment.zero();
   _deviatoric_trial_stress.zero();
  _number_n = _base_models.size();
  // Calculating the Youngs modulus for each of the n curves
  ColumnMajorMatrix A( _number_n, _number_n);
  ColumnMajorMatrix InvA(_number_n, _number_n);

  for (unsigned int i = 0; i < _base_models.size(); i++)
  {
    for (unsigned int j = 0; j < _base_models.size(); j++)
    {
      if (i <= j)
        A(i, j) = _backbone_strain[i];
      else
        A(i, j) = _backbone_strain[j];
    }
     _sudo_backbone_stress[i] = _backbone_stress[i];
     _sudo_backbone_strain[i] = _backbone_strain[i];
  }

  // create upper triangular matrix and modified stress
//  std::vector<Real> _youngs(_number_n);
  std::vector<Real> G0_component(_number_n);
//  std::vector<Real> _yield_stress(_number_n);
  std::vector<Real> modified_stress(_number_n);  InvA = A;

  for (unsigned int i = 1; i < _base_models.size(); i++)
  {
    for (unsigned int j = 0; j < _base_models.size(); j++)
    {
      InvA(i, j) = A(i, j) - A(i-1, j);
      modified_stress[i] = _backbone_stress[i] - _backbone_stress[i-1];
    }
  }


  modified_stress[0] = _backbone_stress[0];

  // backward substitution
  G0_component[_base_models.size()-1] = modified_stress[_base_models.size() - 1]/InvA(_base_models.size() -1, _base_models.size() - 1);
  _yield_stress[_base_models.size()-1] = _backbone_strain[_base_models.size()-1];
  int i = _base_models.size()- 2;
  while (i >= 0)
  {
    G0_component[i] = 0.0;
    for (unsigned int j = i+1; j < _base_models.size(); j++)
       G0_component[i] += InvA(i, j)* G0_component[j];

    G0_component[i] = (modified_stress[i] - G0_component[i]) / InvA(i,i);
    _yield_stress[i] = _backbone_strain[i];
    i = i-1;
  }
 // printf("p1, p2, p3: %e, %e, %e \n", _p1, _p2, _p3);
 // printf("G0: %e \n", _G0);

  for (unsigned int i = 0; i <_base_models.size(); i++)
  {
 //   printf("G0 component, yield stress0: %e, %e \n", G0_component[i], _yield_stress[i]);
 //   printf("strain, stress %e, %e \n", _backbone_strain[i], _backbone_stress[i]);
  }
  //scaling
  for (unsigned int i = 0; i < _base_models.size(); i++)
  {
    _youngs[i] = G0_component[i] * 2.0 * (1.0 + _Poissons);
    _yield_stress[i] = _yield_stress[i]*std::sqrt(3.0)/ (2.0 * (1.0 + _Poissons));
  }
}

void
ComputeMAT79Stress::initQpStatefulProperties()
{
  computeNonmasing();

  for (unsigned int i = 0; i < _base_models.size(); i++)
  {
    (*_stress_model[i])[_qp].zero();
    (*_stress_model_old[i])[_qp].zero();
  }

   _reversal_check[_qp] = _initial_reversal_check;
   _strain_memory[_qp] = _initial_strain_memory;
  // determine the lateral and vertical stresses
   Real total_vertical_stress = _stress[_qp](2,2);
   Real total_lateral_xx = _stress[_qp](0,0);
   Real total_lateral_yy = _stress[_qp](1,1);
   Real residual_vertical =  total_vertical_stress;
   Real residual_xx = total_lateral_xx;
   Real residual_yy = total_lateral_yy;
   Real mean_stress = _stress[_qp].trace()/(-3.0);
   _stiffness_pressure_correction = pow((mean_stress - _p0) / _p_ref, _b_exp);
   _strength_pressure_correction = std::sqrt(_a0+_a1*(mean_stress - _p0)+ _a2*(mean_stress - _p0)*(mean_stress - _p0))/std::sqrt(_a0+_a1*(_p_ref)+ _a2*(_p_ref)*(_p_ref));
   Real total_stress = 0.0;
   // Calculate the K0 consistent stress distribution
   RankTwoTensor dev_model;
   for (unsigned int i = 0; i < _base_models.size(); i++)
   {
     Real mean_pressure = 0.0;
     if (residual_vertical != 0.0)
     {
       Real sum_youngs = 0.0;

       for (unsigned int j = i; j < _base_models.size(); j++)
         sum_youngs += _youngs[j];

      (*_stress_model[i])[_qp](2,2) = residual_vertical * _youngs[i]/sum_youngs;
      (*_stress_model[i])[_qp](0,0) = residual_xx * _youngs[i]/sum_youngs;
      (*_stress_model[i])[_qp](1,1) = residual_yy * _youngs[i]/sum_youngs;
      dev_model =  ((*_stress_model[i])[_qp]).deviatoric() / _youngs[i];
      mean_pressure = (*_stress_model[i])[_qp].trace() / 3.0;
      Real J2_model = dev_model.doubleContraction(dev_model);
      Real dev_stress_model = std::sqrt(3.0 / 2.0 * J2_model);
      if (dev_stress_model > _yield_stress[i] * _strength_pressure_correction)
        dev_model *= (_yield_stress[i] * _strength_pressure_correction) / dev_stress_model;

      (*_stress_model[i])[_qp] = dev_model * _youngs[i]; // stress_model contains only the deviatoric part of the stress
     }
     residual_vertical = residual_vertical - (*_stress_model[i])[_qp](2,2) - mean_pressure;
     residual_xx = residual_xx - (*_stress_model[i])[_qp](0,0) - mean_pressure;
     residual_yy = residual_yy - (*_stress_model[i])[_qp](1,1) - mean_pressure;
     total_stress += (*_stress_model[i])[_qp](2,2)+ mean_pressure;
  }
  ComputeStressBase::initQpStatefulProperties();
  _stress_old[_qp] = _stress[_qp];
}

void
ComputeMAT79Stress::computeQpStress()
{
  // Nothing to update during the first time step, return immediately
  if (_t_step == 0)
    return;

  _stress_new.zero();
  computeStress();
  _stress[_qp] = _rotation_increment[_qp] * _stress_new * _rotation_increment[_qp].transpose();

  //Compute dstress_dstrain
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp] * _tangent_modulus; //This is NOT the exact jacobian
}

void
ComputeMAT79Stress::computeStress()
{
  if (_t_step == 0)
    return;

  _individual_stress_increment.zero();
  _deviatoric_trial_stress.zero();
  _tangent_modulus = 0.0;

  // current pressure calculation
  //Real mean_stress = pow(-_K0 / pow(_p_ref, _b_exp) * log (1.0 + _strain_increment[_qp].trace()) * (1-_b_exp) + pow(_stress_old[_qp].trace()/(-3.0) - _p0, 1.0 - _b_exp), 1.0/(1.0-_b_exp)) + _p0;
  Real mean_stress = _stress_old[_qp].trace() / (-3.0);
  if (mean_stress < _p0)
    mean_stress = _p0;

  _stiffness_pressure_correction = pow((mean_stress - _p0) / _p_ref, _b_exp);
  _strength_pressure_correction = std::sqrt(_a0 + _a1 * (mean_stress - _p0) + _a2 * (mean_stress - _p0) * (mean_stress - _p0)) / std::sqrt(_a0 + _a1 * (_p_ref) + _a2 * (_p_ref) * (_p_ref));
  Real mean_pressure = 0.0;
//    printf("strain_increment:  %e \n",_strain_increment[_qp](0,2));
  for (unsigned int i = 0; i < _base_models.size(); i++)
  {
//    printf("stress:  %e \n", (*_stress_model[i])[_qp](0,2) );
    // compute trial stress increment - note that _elasticity_tensor here assumes youngs modulus = 1
    _individual_stress_increment = _elasticity_tensor[_qp] * (_strain_increment[_qp]);
    // calculate pressure for each element
    mean_pressure += _individual_stress_increment.trace() / 3.0 * _youngs[i] * _stiffness_pressure_correction;

    // compute the deviatoric trial stress normalized by non pressure dependent youngs modulus.
    _deviatoric_trial_stress = _individual_stress_increment.deviatoric() * _stiffness_pressure_correction + (*_stress_model_old[i])[_qp] / (_youngs[i]);

    // compute the effective trial stress
    Real dev_trial_stress_squared = _deviatoric_trial_stress.doubleContraction(_deviatoric_trial_stress);
    Real effective_trial_stress = std::sqrt(3.0 / 2.0 * dev_trial_stress_squared);

    // check yield condition and calculate plastic strain
    Real yield_condition = effective_trial_stress - _yield_stress[i] * _strength_pressure_correction;

    if (yield_condition >= 0.0)
        _deviatoric_trial_stress *= _yield_stress[i] * _strength_pressure_correction / effective_trial_stress;
    else
      _tangent_modulus += _youngs[i];

    (*_stress_model[i])[_qp] = _youngs[i] * (_deviatoric_trial_stress);
//    printf("zx:  %e \n", (*_stress_model[i])[_qp](0,2));

    _stress_new += (*_stress_model[i])[_qp];

//   printf("stress:  %e \n", (*_stress_model[i])[_qp](0,2) );
  }

  _stress_new(0,0) += mean_pressure - mean_stress;
  _stress_new(1,1) += mean_pressure - mean_stress;
  _stress_new(2,2) += mean_pressure - mean_stress;

   if (_p1 != 1.0 && _p2 != 0.0)
    {
   // check if unloading - reloading occured
    RankTwoTensor dev_stress;
    RankTwoTensor dev_stress_old;
    RankTwoTensor dev_strain;
    RankTwoTensor dev_strain_old;
    dev_stress = _stress_new.deviatoric();
    Real dev_stress_squared = 0.5*(dev_stress(0,0)*dev_stress(0,0)+dev_stress(1,1)*dev_stress(1,1)+dev_stress(2,2)*dev_stress(2,2)) + dev_stress(0,1)*dev_stress(0,1) + dev_stress(0,2)*dev_stress(0,2)+dev_stress(1,2)*dev_stress(1,2);
    Real eff_stress = std::sqrt(3.0 / 2.0 * dev_stress_squared);
    dev_stress_old = _stress_old[_qp].deviatoric();
    Real dev_stress_old_squared = 0.5*(dev_stress_old(0,0)*dev_stress_old(0,0)+dev_stress_old(1,1)*dev_stress_old(1,1)+dev_stress_old(2,2)*dev_stress_old(2,2)) + dev_stress_old(0,1)*dev_stress_old(0,1) + dev_stress_old(0,2)*dev_stress_old(0,2)+dev_stress_old(1,2)*dev_stress_old(1,2);
    Real eff_stress_old = std::sqrt(3.0 / 2.0 * dev_stress_old_squared);
    dev_strain = _total_strain[_qp].deviatoric();
    dev_strain_old = _total_strain_old[_qp].deviatoric();
    Real dev_strain_squared = pow((2.0*dev_strain(0,1)), 2.0) + pow((2.0*dev_strain(0,2)),2.0) +pow((2.0*dev_strain(1,2)), 2.0) ;
    Real dev_strain_squared_old = pow((2.0*dev_strain_old(0,1)), 2.0) + pow((2.0*dev_strain_old(0,2)),2.0) +pow((2.0*dev_strain_old(1,2)), 2.0) ;
	Real eff_shear_strain = std::sqrt(dev_strain_squared);
	Real eff_shear_strain_old = std::sqrt(dev_strain_squared_old);

	if ((eff_shear_strain < (eff_shear_strain_old-_backbone_strain[0])) && (eff_shear_strain > _strain_memory_old[_qp]) && (eff_shear_strain > _backbone_strain[0]))
	  {
		_reversal_check[_qp] = _reversal_check_old[_qp] + 1.0;
		_strain_memory[_qp] = eff_shear_strain;
//		printf("reversal_check:  %e \n", _reversal_check[_qp]);
//		printf("eff_shear_strain:  %e \n",eff_shear_strain);
//		printf("eff_shear_strain_old:  %e \n",eff_shear_strain_old);
//		printf("total_strain:  %e \n",_total_strain[_qp](0,2));
//		printf("total_strain_old:  %e \n",_total_strain_old[_qp](0,2));
	  }
	else
          {
	    _reversal_check[_qp] = _reversal_check_old[_qp]+0.0;
	    _strain_memory[_qp] = _strain_memory_old[_qp]+0.0;
	  }
//	printf("reversal_check:  %e \n", _reversal_check[_qp]);
//	printf("total_strain_zx:  %e \n", _total_strain[_qp](0,2));

//	if  (_reversal_check[_qp] > _reversal_check_old[_qp])
//	  {

//		_reduction_factor = 1.0 ;
//		printf("reduction factor:  %e \n", _reduction_factor);
//	  }


   if (_reversal_check[_qp] > _reversal_check_old[_qp])
    { 
	int record = 0;
	Real correction = _stiffness_pressure_correction/_strength_pressure_correction;

	for (unsigned int i = 0; i < _base_models.size(); i++)

// interpolate for intermediate strains
	  {	
//		printf("_backbone_strain:  %e \n", _sudo_backbone_strain[i]);
//		printf("_backbone_stress:  %e \n", _sudo_backbone_stress[i]);
		if ((i > 0))
		 {
		//printf("_backbone_strain_i-1:  %e \n", _backbone_strain[i-1]);
		_sudo_backbone_strain[i] = _backbone_strain[i];
		_sudo_backbone_stress[i] = _backbone_stress[i];    
		if ((_sudo_backbone_strain[i-1] < correction*eff_shear_strain) && (abs(_sudo_backbone_strain[i-1] - correction*eff_shear_strain) > 1e-6 ))
		 {
			if  (_sudo_backbone_strain[i] > correction*eff_shear_strain)
			 {
				record = i;
            		//	printf("record:  %i \n", record);			
		 	 }
		 }
		 }
	  }
	if (record != 0)
	 {

		_sudo_backbone_stress[record] = _sudo_backbone_stress[record-1] + (_sudo_backbone_stress[record] - _sudo_backbone_stress[record - 1])* (correction*eff_shear_strain - _sudo_backbone_strain[record-1])/(_sudo_backbone_strain[record] - _sudo_backbone_strain[record-1]);
		_sudo_backbone_strain[record] = eff_shear_strain*correction;
		_Gm = _sudo_backbone_stress[record] / _sudo_backbone_strain[record];
//		printf("dev_stress:  %e \n", std::sqrt(dev_stress_old_squared));
		printf("Gm:  %e \n", _Gm);
	_reduction_factor =  _p1 - _p2*pow((1.0 - _Gm / (_sudo_backbone_stress[0]/_sudo_backbone_strain[0])), _p3);
		printf("reduction_factor:  %e \n", _reduction_factor);

	for (unsigned int i = 0; i < _base_models.size(); i++)
	  {	
		printf("sudo_backbone_strain:  %e \n", _sudo_backbone_strain[i]);
		printf("sudo_backbone_stress:  %e \n", _sudo_backbone_stress[i]);
		if (_sudo_backbone_strain[i]<= correction*eff_shear_strain)
		 {
			_mapped_backbone_stress[i] = _sudo_backbone_stress[i]*_reduction_factor + (1.0-_reduction_factor)*_Gm*_sudo_backbone_strain[i];
			_mapped_backbone_strain[i] = _sudo_backbone_strain[i];
		 }
		else
		 {
			_mapped_backbone_stress[i] = _sudo_backbone_stress[i];
			_mapped_backbone_strain[i] = _sudo_backbone_strain[i];
		 }
		printf("mapped_backbone_strain:  %e \n", _mapped_backbone_strain[i]);
		printf("mapped_backbone_stress:  %e \n", _mapped_backbone_stress[i]);		
	  }

	   _stress_new.zero();
	   _individual_stress_increment.zero();
	   _deviatoric_trial_stress.zero();
	  _number_n = _base_models.size();
	  // Calculating the Youngs modulus for each of the n curves
	  ColumnMajorMatrix A( _number_n, _number_n);
	  ColumnMajorMatrix InvA(_number_n, _number_n);

	  for (unsigned int i = 0; i < _base_models.size(); i++)
	  {
		for (unsigned int j = 0; j < _base_models.size(); j++)
		{
		  if (i <= j)
		    A(i, j) = _mapped_backbone_strain[i];
		  else
		    A(i, j) = _mapped_backbone_strain[j];
		}
	  }

	  // create upper triangular matrix and modified stress
	//  std::vector<Real> _youngs(_number_n);
	  std::vector<Real> G0_component(_number_n);
	//  std::vector<Real> _yield_stress(_number_n);
	  std::vector<Real> modified_stress(_number_n);  InvA = A;

	  for (unsigned int i = 1; i < _base_models.size(); i++)
	  {
		for (unsigned int j = 0; j < _base_models.size(); j++)
		{
		  InvA(i, j) = A(i, j) - A(i-1, j);
		  modified_stress[i] = _mapped_backbone_stress[i] - _mapped_backbone_stress[i-1];
		}
	  }


	  modified_stress[0] = _mapped_backbone_stress[0];

	  // backward substitution
	  G0_component[_base_models.size()-1] = modified_stress[_base_models.size() - 1]/InvA(_base_models.size() -1, _base_models.size() - 1);
	  _yield_stress[_base_models.size()-1] = _mapped_backbone_strain[_base_models.size()-1];
	  int i = _base_models.size()- 2;
	  while (i >= 0)
	  {
		G0_component[i] = 0.0;
		for (unsigned int j = i+1; j < _base_models.size(); j++)
		   G0_component[i] += InvA(i, j)* G0_component[j];

		G0_component[i] = (modified_stress[i] - G0_component[i]) / InvA(i,i);
		_yield_stress[i] = _mapped_backbone_strain[i];
		i = i-1;
	  }
	//  printf("p1, p2, p3: %e, %e, %e \n", _p1, _p2, _p3);
	//  printf("G0: %e \n", _G0);

//	  for (unsigned int i = 0; i <_base_models.size(); i++)
//	  {
//		printf("G0 component, yield stress0: %e, %e \n", G0_component[i], _yield_stress[i]);
//
//	  }
//	printf("zx:  %e \n", _stress_new(0,2)); 
	  //scaling
	  for (unsigned int i = 0; i < _base_models.size(); i++)
	  {
		_youngs[i] = G0_component[i] * 2.0 * (1.0 + _Poissons);
		_yield_stress[i] = _yield_stress[i]*std::sqrt(3.0)/ (2.0 * (1.0 + _Poissons));
	  }
   	  Real dev_model;
   	  Real residual_vertical_stress = 0.0;
     	  Real residual_lateral_xx = 0.0;
      	  Real residual_lateral_yy = 0.0;
      	  Real residual_zx = 0.0;
	  Real residual_xy = 0.0;
	  Real residual_yz = 0.0;
	  for (unsigned int i = 0; i < _base_models.size(); i++)
	  {

  		Real mean_stress = _stress_old[_qp].trace() / (-3.0);
   	//	printf(" zz %e \n",mean_stress);

	RankTwoTensor component_stress = (*_stress_model_old[i])[_qp];
//	printf("zx:  %e \n", (*_stress_model_old[i])[_qp](0,2));
        (*_stress_model_old[i])[_qp](0,0) +=  residual_lateral_xx;
        (*_stress_model_old[i])[_qp](1,1) +=  residual_lateral_yy;
        (*_stress_model_old[i])[_qp](2,2) +=  residual_vertical_stress;
        (*_stress_model_old[i])[_qp](0,1) +=  residual_xy;
        (*_stress_model_old[i])[_qp](0,2) +=  residual_zx;
        (*_stress_model_old[i])[_qp](1,2) +=  residual_yz;
	(*_stress_model_old[i])[_qp](1,0) = (*_stress_model_old[i])[_qp](0,1);
	(*_stress_model_old[i])[_qp](2,0) = (*_stress_model_old[i])[_qp](0,2);
	(*_stress_model_old[i])[_qp](2,1) = (*_stress_model_old[i])[_qp](1,2);

      	RankTwoTensor component_deviatoric =  ((*_stress_model_old[i])[_qp]).deviatoric() / _youngs[i];
//    	printf(" youngs %e \n", _youngs[i]);
      	Real component_mean_stress = (*_stress_model_old[i])[_qp].trace() / 3.0;
//    	printf(" component mean sterss %e \n", component_mean_stress);
      	Real J2_model = component_deviatoric.doubleContraction(component_deviatoric);
      	Real dev_stress_model = std::sqrt(3.0 / 2.0 * J2_model);
//  		_stiffness_pressure_correction = pow((mean_stress - _p0) / _p_ref, _b_exp);
//  		_strength_pressure_correction = std::sqrt(_a0 + _a1 * (mean_stress - _p0) + _a2 * (mean_stress - _p0) * (mean_stress - _p0)) / std::sqrt(_a0 + _a1 * (_p_ref) + _a2 * (_p_ref) * (_p_ref));
//	printf("i:  %i \n", i);
//	printf("G0 component, yield stress0: %e, %e \n", G0_component[i]/2.0 * (1.0 + _Poissons), _yield_stress[i]);
//	printf("zx:  %e \n", (*_stress_model_old[i])[_qp](0,2));
//	printf("dev_stress_model:  %e \n", dev_stress_model);	
//	printf("yield_stress:  %e \n", _yield_stress[i] * _strength_pressure_correction);	      	
		if (dev_stress_model >= _yield_stress[i] * _strength_pressure_correction)
		{ 
        	component_deviatoric *= (_yield_stress[i] * _strength_pressure_correction)/dev_stress_model ;
     		(*_stress_model_old[i])[_qp] = component_deviatoric * _youngs[i]; // stress_model contains only the deviatoric part of the stress
		(*_stress_model_old[i])[_qp](0,0) += component_mean_stress; 
		(*_stress_model_old[i])[_qp](1,1) += component_mean_stress; 
		(*_stress_model_old[i])[_qp](2,2) += component_mean_stress;
   	  	residual_vertical_stress += component_stress(2,2)-(*_stress_model_old[i])[_qp](2,2);
      		residual_lateral_xx += component_stress(0,0)-(*_stress_model_old[i])[_qp](0,0);
      		residual_lateral_yy += component_stress(1,1)-(*_stress_model_old[i])[_qp](1,1);
      		residual_zx += component_stress(0,2)-(*_stress_model_old[i])[_qp](0,2);
	  	residual_xy += component_stress(0,1)-(*_stress_model_old[i])[_qp](0,1);
	  	residual_yz += component_stress(1,2)-(*_stress_model_old[i])[_qp](1,2); 
		}
 	    	else
		{
     	    	residual_vertical_stress = 0.0;
      		residual_lateral_xx = 0.0;
      		residual_lateral_yy = 0.0;
      		residual_zx = 0.0;
	        residual_xy = 0.0;
	        residual_yz = 0.0;   
 		}

//	printf("residual_zx:  %e \n", residual_zx);
	  } 
//	printf("zx:  %e \n", _stress_old[_qp](0,2));
//	printf("zx:  %e \n", _stress_new(0,2)); 

  _individual_stress_increment.zero();
  _deviatoric_trial_stress.zero();
  _tangent_modulus = 0.0;

  // current pressure calculation
  //Real mean_stress = pow(-_K0 / pow(_p_ref, _b_exp) * log (1.0 + _strain_increment[_qp].trace()) * (1-_b_exp) + pow(_stress_old[_qp].trace()/(-3.0) - _p0, 1.0 - _b_exp), 1.0/(1.0-_b_exp)) + _p0;
  Real mean_stress = _stress_old[_qp].trace() / (-3.0);
  if (mean_stress < _p0)
    mean_stress = _p0;

  _stiffness_pressure_correction = pow((mean_stress - _p0) / _p_ref, _b_exp);
  _strength_pressure_correction = std::sqrt(_a0 + _a1 * (mean_stress - _p0) + _a2 * (mean_stress - _p0) * (mean_stress - _p0)) / std::sqrt(_a0 + _a1 * (_p_ref) + _a2 * (_p_ref) * (_p_ref));
  Real mean_pressure = 0.0;
  for (unsigned int i = 0; i < _base_models.size(); i++)
  {
//    printf("stress:  %e \n", (*_stress_model[i])[_qp](0,2) );
    // compute trial stress increment - note that _elasticity_tensor here assumes youngs modulus = 1
    _individual_stress_increment = _elasticity_tensor[_qp] * (_strain_increment[_qp]);

    // calculate pressure for each element
    mean_pressure += _individual_stress_increment.trace() / 3.0 * _youngs[i] * _stiffness_pressure_correction;

    // compute the deviatoric trial stress normalized by non pressure dependent youngs modulus.
    _deviatoric_trial_stress = _individual_stress_increment.deviatoric() * _stiffness_pressure_correction + (*_stress_model_old[i])[_qp] / (_youngs[i]);

    // compute the effective trial stress
    Real dev_trial_stress_squared = _deviatoric_trial_stress.doubleContraction(_deviatoric_trial_stress);
    Real effective_trial_stress = std::sqrt(3.0 / 2.0 * dev_trial_stress_squared);

    // check yield condition and calculate plastic strain
    Real yield_condition = effective_trial_stress - _yield_stress[i] * _strength_pressure_correction;

    if (yield_condition >= 0.0)
       {
        _deviatoric_trial_stress *= _yield_stress[i] * _strength_pressure_correction / effective_trial_stress;
//	printf("i:  %i \n", i);
	}
    else
       {
      _tangent_modulus += _youngs[i];
	}
    (*_stress_model[i])[_qp] = _youngs[i] * (_deviatoric_trial_stress);

    _stress_new += (*_stress_model[i])[_qp];


  }
//   printf("stress:  %e \n", _stress_new(0,2) );
  _stress_new(0,0) += mean_pressure - mean_stress;
  _stress_new(1,1) += mean_pressure - mean_stress;
  _stress_new(2,2) += mean_pressure - mean_stress;
	 }

    }
    }
  _tangent_modulus *= _stiffness_pressure_correction;
  //_stress_new(0,0) += -mean_stress;
  //_stress_new(1,1) += -mean_stress;
  //_stress_new(2,2) += -mean_stress;
}
