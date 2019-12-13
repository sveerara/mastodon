[Mesh]
  type = FileMesh #Read in mesh from file
  file = 'Model.e'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
#  use_displaced_mesh = false
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[AuxVariables]
  [./vel_x]
  [../]
  [./accel_x]
  [../]
  [./vel_y]
  [../]
  [./accel_y]
  [../]
  [./vel_z]
  [../]
  [./accel_z]
  [../]
[]

[Kernels]
  [./DynamicTensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = false
  [../]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.25
    gamma = 0.5
    use_displaced_mesh = false
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.25
    gamma = 0.5
    use_displaced_mesh = false
  [../]
  [./inertia_z]
    type = InertialForce
    variable = disp_z
    velocity = vel_z
    acceleration = accel_z
    beta = 0.25
    gamma = 0.5
    use_displaced_mesh = false
  [../]
  [./damping_force_x]
    type = DampingForce
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    beta = 0.25
    gamma = 0.5
    component = 0
    damping_matrix_userobject = fid_uo
    use_displaced_mesh = false
  [../]
  [./damping_force_y]
    type = DampingForce
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    beta = 0.25
    gamma = 0.5
    component = 1
    damping_matrix_userobject = fid_uo
    use_displaced_mesh = false
  [../]
  [./damping_force_z]
    type = DampingForce
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    beta = 0.25
    gamma = 0.5
    component = 2
    damping_matrix_userobject = fid_uo
    use_displaced_mesh = false
  [../]
  [./gravity_z]
    type = Gravity
    variable = disp_z
    value = -9.81
    use_displaced_mesh = false
  [../]
[]

[UserObjects]
  [./fid_uo]
    type = FrequencyIndependentDamper
    damping_coefficient = 0.03
    beta = 0.25
  [../]
[]

[AuxKernels]
  [./accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.25
    execute_on = timestep_end
  [../]
  [./vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = 0.25
    execute_on = timestep_end
  [../]
  [./vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./accel_z]
    type = NewmarkAccelAux
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = 0.25
    execute_on = timestep_end
  [../]
  [./vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = 0.5
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./bottom_z]
    type = PresetBC
    variable = disp_z
    boundary = 101
    value = 0.0
  [../]
  [./bottom_y]
    type = PresetBC
    variable = disp_y
    boundary = 101
    value = 0.0
  [../]
  [./lateral_1]
    type = PresetBC
    variable = disp_y
    boundary = 1024
    value = 0.0
  [../]
  [./lateral_2]
    type = PresetBC
    variable = disp_y
    boundary = 1022
    value = 0.0
  [../]
  [./bottom_accel]
	type = PresetAcceleration
    boundary = 101
    function = acceleration_bottom
    variable = disp_x
    beta = 0.25
	acceleration = accel_x
	velocity = vel_x
  [../]
  [./Periodic]
    [./y_dir]
      variable = 'disp_x disp_y disp_z'
      primary = 1024
      secondary = 1022
      translation = '0.0 1.0 0.0'
    [../]
    [./x_dir]
      variable = 'disp_x disp_y disp_z'
      primary = 1025
      secondary = 1023
      translation = '1.0 0.0 0.0'
    [../]
  [../]
[]

[Functions]
  [./acceleration_bottom]
    type = PiecewiseLinear
    data_file = chichi_bc_mpss.csv
    format = columns
  [../]
[]

[Materials]
  [./strain1]
    #Computes the strain, assuming small strains
    type = ComputeIncrementalSmallStrain
    block = '10001 10002 10003 10004 10005 10006 10007 10008 10009 10010 10011 10012 10013 10014 10015 10016 10017 10018 10019 10020'
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress1]
    #Computes the stress, using linear elasticity
    type = ComputeFiniteStrainElasticStress
    store_stress_old = true
    block = '10001 10002 10003 10004 10005 10006 10007 10008 10009 10010 10011 10012 10013 10014 10015 10016 10017 10018 10019 10020'
  [../]
  [./mat1]
    youngs_modulus = 325e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10001
  [../]
  [./mat2]
    youngs_modulus = 307e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10002
  [../]
  [./mat3]
    youngs_modulus = 290e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10003
  [../]
  [./mat4]
    youngs_modulus = 270e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10004
  [../]
  [./mat5]
    youngs_modulus = 252e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10005
  [../]
  [./mat6]
    youngs_modulus = 234e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10006
  [../]
  [./mat7]
    youngs_modulus = 216e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10007
  [../]
  [./mat8]
    youngs_modulus = 200e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10008
  [../]
  [./mat9]
    youngs_modulus = 184e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10009
  [../]
  [./mat10]
    youngs_modulus = 168e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10010
  [../]
  [./mat11]
    youngs_modulus = 154e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10011
  [../]
  [./mat12]
    youngs_modulus = 140e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10012
  [../]
  [./mat13]
    youngs_modulus = 127e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10013
  [../]
  [./mat14]
    youngs_modulus = 114e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10014
  [../]
  [./mat15]
    youngs_modulus = 102e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10015
  [../]
  [./mat16]
    youngs_modulus = 91e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10016
  [../]
  [./mat17]
    youngs_modulus = 80e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10017
  [../]
  [./mat18]
    youngs_modulus = 70e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10018
  [../]
  [./mat19]
    youngs_modulus = 61e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10019
  [../]
  [./mat20]
    youngs_modulus = 52e6 #Pa
    poissons_ratio = 0.3
    type = ComputeIsotropicElasticityTensor
    block = 10020
  [../]
  [./den1]
    type = GenericConstantMaterial
    block = '10001 10002 10003 10004 10005 10006 10007 10008 10009 10010 10011 10012 10013 10014 10015 10016 10017 10018 10019 10020'
    prop_names = density
    prop_values = 2000 #kg/m3
  [../]
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./accel_bot]
    type = PointValue
    variable = accel_x
    point = '0.0 0.0 -0.5'
  [../]
  [./accel_mid]
    type = PointValue
    variable = accel_x
    point = '0.0 0.0 9.5'
  [../]
  [./accel_top_1]
    type = PointValue
    variable = accel_x
    point = '0.0 0.0 18.5'
  [../]
  [./accel_top]
    type = PointValue
    variable = accel_x
    point = '0.0 0.0 19.5'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-8
  l_tol = 1e-8
  start_time = 0
  end_time = 85.4
  dt = 0.005
  dtmin = 0.005
  timestep_tolerance = 1e-8
#  petsc_options = '-snes_check_jacobian'
#  petsc_options = '-snes_ksp_ew'
#  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
#  petsc_options_value = '201                hypre    boomeramg      4'
#  line_search = 'none'
[]

[Outputs]
   exodus = true
   csv = true
   print_perf_log = true
   print_linear_residuals = true
   [./screen]
       type = Console
       max_rows = 1
       interval = 100
  [../]
[]
