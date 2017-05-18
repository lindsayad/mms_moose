[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  elem_type = EDGE3
[]

[Variables]
  [./u]
    order = SECOND
    family = LAGRANGE
  [../]
  [./v]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./u_diff]
    type = Diffusion
    variable = u
  [../]
  [./v_diff]
    type = Diffusion
    variable = v
  [../]
  [./u_coupled_gradient_source]
    type = CoupledGradientSource
    variable = u
    v = v
  [../]
  [./u_source]
    type = UserForcingFunction
    function = u_source_func
    variable = u
  [../]
  [./v_source]
    type = UserForcingFunction
    function = v_source_func
    variable = v
  [../]
[]

[BCs]
  [./u]
    type = FunctionDirichletBC
    boundary = 'left right'
    function = u_func
    variable = u
  [../]
  [./v]
    type = FunctionDirichletBC
    function = v_func
    boundary = 'left right'
    variable = v
  [../]
[]

[Functions]
  [./u_source_func]
    type = ParsedFunction
    value = ''
  [../]
  [./v_source_func]
    type = ParsedFunction
    value = ''
  [../]
  [./u_func]
    type = ParsedFunction
    value = ''
  [../]
  [./v_func]
    type = ParsedFunction
    value = ''
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = NEWTON
  [../]
[]

[Executioner]
  type = Steady
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels'
  petsc_options_value = '300                bjacobi  ilu          4'
  # line_search = none
  nl_rel_tol = 1e-12
  nl_max_its = 10
  l_tol = 1e-6
  l_max_its = 50
[]

[Outputs]
  [./exodus]
    type = Exodus
    file_base = 'forward_test'
  [../]
  [./csv]
    type = CSV
    file_base = 'forward_test'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Postprocessors]
  [./L2u]
    type = ElementL2Error
    variable = u
    function = u_func
    outputs = 'console csv'
  [../]
  [./L2v]
    variable = v
    function = v_func
    type = ElementL2Error
    outputs = 'console csv'
  [../]
  [./L2nu]
    type = NodalL2Error
    variable = u
    function = u_func
    outputs = 'console csv'
  [../]
  [./L2nv]
    variable = v
    function = v_func
    type = NodalL2Error
    outputs = 'console csv'
  [../]
[]
