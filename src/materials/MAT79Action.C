/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "MAT79Action.h"

#include "Factory.h"
#include "FEProblem.h"
#include "Parser.h"
#include "Conversion.h"

template<>
InputParameters validParams<MAT79Action>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Set up MAT79 material model");

  params.addRequiredParam<std::vector<NonlinearVariableName> >("displacements", "displacements in the problem");
  params.addParam<FileName>("data_file"," ", "File name of the file containing stress-strain curve");
  params.addRequiredParam<Real>("poissons", "The poissons ratio for the material");
  params.addRequiredParam<Real>("b_exp", "exponential factor for pressure dependent stiffness");
  params.addRequiredParam<Real>("a0", "pressure dependent strength parameter 0");
  params.addRequiredParam<Real>("a1", "pressure dependent strength parameter 1");
  params.addRequiredParam<Real>("a2", "pressure dependent strength parameter 2");
  params.addRequiredParam<Real>("p0", "tension cut-off");
  params.addRequiredParam<Real>("p1", "MRDF type nonmasing parameter");
  params.addRequiredParam<Real>("p2", "MRDF type nonmasing parameter");
  params.addRequiredParam<Real>("p3", "MRDF type nonmasing parameter");
  params.addParam<Real>("K0", -1.0, "Initial Bulk Modulus");
  params.addParam<Real>("G0",-1.0, "Initial Shear Modulus");
  params.addRequiredParam<Real>("p_ref", "Reference pressure at which the parameters are defined");
  params.addParam<Real>("OCR",-1.0, "over consolidation ratio");
  params.addParam<Real>("PI",-1.0, "plasticity index");
  params.addParam<Real>("theta_1",1.0, "curve fit coefficient");
  params.addParam<Real>("theta_2",1.0, "curve fit coefficient");
  params.addParam<Real>("theta_3",1.0, "curve fit coefficient");
  params.addParam<Real>("theta_4",1.0, "curve fit coefficient");
  params.addParam<Real>("theta_5",1.0, "curve fit coefficient");
  params.addParam<Real>("taumax",-1.0, "ultimate shear strength");
  params.addParam<int>("numberofpoints",-1, "total number of data points to define backbone curve");
  params.addRequiredParam<int>("soiltype", "determined whether it is user defined soil or Darendeli formulation");
  params.addRequiredParam<std::vector<SubdomainName> >("block", "The blocks where this material model is applied.");
  params.addParam<std::vector<FunctionName> >("initial_stress", "The function values for initial_stress");
  return params;
}

MAT79Action::MAT79Action(const InputParameters & params) :
  Action(params),
  _data_file(getParam<FileName>("data_file")),
  _Poissons(getParam<Real>("poissons")),
  _b_exp(getParam<Real>("b_exp")),
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
  _OCR(getParam<Real>("OCR")),
  _PI(getParam<Real>("PI")),
  _fit_2_theta1(getParam<Real>("theta_1")),
  _fit_2_theta2(getParam<Real>("theta_2")),
  _fit_2_theta3(getParam<Real>("theta_3")),
  _fit_2_theta4(getParam<Real>("theta_4")),
  _fit_2_theta5(getParam<Real>("theta_5")),
  _taumax(getParam<Real>("taumax")),
  _numberofpoints(getParam<int>("numberofpoints")),
  _soiltype(getParam<int>("soiltype")),
  _block(getParam<std::vector<SubdomainName> >("block"))
{
  _number = 1;
}

void
MAT79Action::act()
{
  std::vector<Real> stress;
  std::vector<Real> strain;
  std::vector<Real> modulus;


   if (_soiltype == 0)
    {
     if (_data_file == " ")
      {
        mooseError("Invalid stress strain curve. Data file should be provided");
      }
     parseColumns(strain, stress);

    }
   else if (_soiltype == 1)
    {
     if (_OCR <= 0.0 || _PI <= 0.0  || _numberofpoints <= 0 || _K0 <= 0.0 || _G0 <= 0.0 )
      {
        mooseError("OCR, PI,numberofpoints,K0, G0 have to be provided. Please provide positive values");
      }
     strain.resize(_numberofpoints);
     stress.resize(_numberofpoints);
     modulus.resize(_numberofpoints);
     Real phi1 = 0.0352;
     Real phi2 = 0.001;
     Real phi3 = 0.3246;
     Real phi4 = 0.3483;
     Real phi5 = 0.9190;
     Real ref_strain = (phi1+phi2*_PI*pow(_OCR,phi3))*pow(_p_ref*0.00987,phi4);
     printf("ref_strain:  %e \n", ref_strain);
     printf("numberofpoints:  %d \n", _numberofpoints);
     for (unsigned int i = 0; i < _numberofpoints; i++)
        {
        strain[i] = pow(10.0,-6.0+5.0/(_numberofpoints-1)*i);
        printf("strain:  %e \n", strain[i]);
        stress[i] = _G0*(1.0/(1.0+pow(100.0*strain[i]/ref_strain,phi5)))*strain[i];
        }
     }
   else if (_soiltype == 2)
    {
//     if (  _fit_2_theta1+_fit_2_theta2+_fit_2_theta3+_fit_2_theta4+_fit_2_theta5 >1.0 || //_numberofpoints <= 0  || _taumax <= 0.0 || _K0<= 0.0 || _G0<= 0.0)
//      {
//        mooseError("theta_1 through 5, K0, G0, taumax numberofpoints have to be provided. Please provide the appropriate values e.g. taumax, K0,G0, numberpoints should be positive. Sum of theta 1 through 5 should be smaller or equal to 1  ");
//      }

     strain.resize(_numberofpoints);
     stress.resize(_numberofpoints);
     modulus.resize(_numberofpoints);
     for (unsigned int i = 0; i < _numberofpoints; ++i)
     {
            strain[i] = pow(10.0,-6.0+5.0/(_numberofpoints-1)*i);
           Real theta = _fit_2_theta1 + _fit_2_theta2 * pow(strain[i]/(_taumax/_G0),_fit_2_theta5) /(pow(_fit_2_theta3, _fit_2_theta5) + pow(strain[i]/(_taumax/_G0), _fit_2_theta5));
            modulus[i] = (theta == 0) ? 1.0 / (1.0 + strain[i]/(_taumax/_G0)) : (1.0 / (strain[i]/(_taumax/_G0)))*(0.5 / theta*(1.0 + strain[i]/(_taumax/_G0) - sqrt(pow(1.0 + strain[i]/(_taumax/_G0), 2.0) - 4.0*theta*strain[i]/(_taumax/_G0))));
            stress[i] = modulus[i]*_G0*strain[i];
     }
    }


  _number = stress.size();
  // create small incremental strain block
  std::vector<NonlinearVariableName> displacements = getParam<std::vector<NonlinearVariableName> >("displacements");
  unsigned int ndisp = displacements.size();

  std::vector<VariableName> coupled_displacements;
  for (unsigned int i = 0; i < ndisp; ++i)
    coupled_displacements.push_back(displacements[i]);

  InputParameters params = _factory.getValidParams("ComputeIncrementalSmallStrain");
  std::string unique_strain_name = "strain_"+ _block[0];
  params.set<std::vector<SubdomainName> >("block") = _block;
  params.set<std::vector<VariableName> >("displacements") = coupled_displacements;
  params.set<bool>("stateful_displacements") = true;
  _problem->addMaterial("ComputeIncrementalSmallStrain", unique_strain_name, params);

  // create Elasticty tensor with E = 1 and poissons ratio taken as input
  params = _factory.getValidParams("ComputeIsotropicElasticityTensor");

  std::string unique_elasticity_name = "elasticity_" + _block[0];
  params.set<std::vector<SubdomainName> >("block") = _block;
  params.set<Real>("youngs_modulus") = 1.0;
  params.set<Real>("poissons_ratio") = _Poissons;
  _problem->addMaterial("ComputeIsotropicElasticityTensor", unique_elasticity_name, params);

  // Stress calculation
  params = _factory.getValidParams("ComputeMAT79Stress");

  std::vector<std::string> base_names(_number);
  for (unsigned int i = 0; i < _number; i++)
    base_names[i] = "MAT79_" + Moose::stringify(i);

  params.set<std::vector<std::string> >("base_models") = base_names;
  params.set<std::vector<Real> >("backbone_stress") = stress;
  params.set<std::vector<Real> >("backbone_strain") = strain;
  params.set<std::vector<SubdomainName> >("block") = _block;
  params.set<Real>("poissons") = _Poissons;
  params.set<Real>("a0") = _a0;
  params.set<Real>("a1") = _a1;
  params.set<Real>("a2") = _a2;
  params.set<Real>("b_exp") = _b_exp;
  params.set<Real>("p0") = _p0;
  params.set<Real>("p1") = _p1;
  params.set<Real>("p2") = _p2;
  params.set<Real>("p3") = _p3;
  params.set<Real>("K0") = _K0;
  params.set<Real>("G0") = _G0;
  params.set<Real>("p_ref") = _p_ref;
  params.set<bool>("store_stress_old") = true;
  params.set<std::vector<FunctionName> >("initial_stress") = getParam<std::vector<FunctionName> >("initial_stress");
  _problem->addMaterial("ComputeMAT79Stress", "stress_" + _block[0], params);


/*  std::vector<MaterialName> names(_number);
  std::vector<std::string> base_names(_number);
  for (unsigned int i = 0; i < _number; ++i)
  {
    // create unique material name for each material
    std::string unique_material_name = "MAT79_mat" + Moose::stringify(i);
    std::string base = "MAT79_" + Moose::stringify(i);

    InputParameters params = _factory.getValidParams("RecomputeRadialReturnIsotropicPlasticity");

    params.set<Real>("yield_stress") = scaled_yield0_component[i];
    params.set<Real>("hardening_constant") = 0.0;
    params.set<std::vector<SubdomainName> >("block") = _block;
    params.set<std::string>("model_name") = base;
    _problem->addMaterial("RecomputeRadialReturnIsotropicPlasticity", unique_material_name, params);

    names[i] = unique_material_name;
    base_names[i] = base;
  }

  // create Elasticty tensor with E = 1 and possions ratio taken as input
  InputParameters params = _factory.getValidParams("ComputeIsotropicElasticityTensor");

  std::string unique_elasticity_name = "elasticity" + Moose::stringify(i);
  params.set<std::vector<SubdomainName> >("block") = _block;
  params.set<Real>("youngs_modulus") = 1.0;
  params.set<Real>("poissons_ratio") = _Poissons;
  _problem->addMaterial("ComputeIsotropicElasticityTensor", unique_elasticity_name, params);

  // create small incremental strain blocks
  params = _factory.getValidParams("ComputeIncrementalSmallStrain");
  std::string unique_strain_name = "strain" + Moose::stringify(i);
  params.set<std::vector<SubdomainName> >("block") = _block;
  params.set<std::vector<VariableName> >("displacements") = coupled_displacements;
  params.set<bool>("stateful_displacements") = true;
  _problem->addMaterial("ComputeIncrementalSmallStrain", unique_strain_name, params);


  // Stress calculation
  params = _factory.getValidParams("ComputeMAT79Stress");

  params.set<std::vector<Real> >("youngs_modulus") = E0_component;
  params.set<std::vector<SubdomainName> >("block") = _block;
  params.set<Real>("poissons_ratio") = _Poissons;
  params.set<Real>("a0") = _a0;
  params.set<Real>("a1") = _a1;
  params.set<Real>("a2") = _a2;
  params.set<Real>("b_exp") = _b_exp;
  params.set<Real>("p0") = _p0;
  params.set<Real>("p_ref") = _p_ref;
  params.set<std::vector<MaterialName> >("return_mapping_models") = names;
  params.set<std::vector<std::string> >("base_models") = base_names;

  _problem->addMaterial("ComputeMAT79Stress", "stress", params);
*/}

bool
MAT79Action::parseNextLineReals(std::ifstream & ifs, std::vector<Real> &myvec)
{
  std::string line;
  myvec.clear();
  bool gotline(false);
  if (getline(ifs,line))
  {
    gotline=true;

    //Replace all commas with spaces
    while (size_t pos=line.find(','))
    {
      if (pos == line.npos)
        break;
      line.replace(pos,1,1,' ');
    }

    //Harvest floats separated by whitespace
    std::istringstream iss(line);
    Real f;
    while (iss>>f)
    {
      myvec.push_back(f);
    }
  }
  return gotline;
}

void
MAT79Action::parseColumns( std::vector<Real> & x, std::vector<Real> & y )
{
  std::ifstream file(_data_file.c_str());
  if (!file.good())
    mooseError("In MAT79Action " + _name + ": Error opening file '" + _data_file + "'.");

  std::vector<Real> scratch;
  unsigned int x_index = 0;
  unsigned int y_index = 1;

  unsigned int line_index = 0;
  while (parseNextLineReals(file, scratch))
  {
    if (scratch.size() > 0)
    {
      if (x_index < scratch.size())
        x.push_back(scratch[x_index]);
      else
        mooseError("In MAT79Action " + _name + ": column " + std::to_string(x_index) + " for x does not exist on line " + std::to_string(line_index));

      if (y_index < scratch.size())
        y.push_back(scratch[y_index]);
      else
        mooseError("In MAT79Action " + _name + ": column " + std::to_string(y_index) + " for y does not exist on line " + std::to_string(line_index));

      if (scratch.size() != 2)
        mooseError("In MAT79Action " + _name + ": Read more than 2 columns of data from file '" + _data_file + "'.");
    }

    line_index++;
  }
}
