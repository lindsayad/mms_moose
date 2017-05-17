[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
[]

[Variables]
  [./u]
  [../]
  [./v]
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
    type = DirichletBC
    boundary = 'left'
    value = 0
    variable = u
  [../]
  [./v]
    type = FunctionDirichletBC
    function = v_func
    boundary = 'left right'
    variable = v
  [../]
  [./u_test]
    type = CoupledGradientBC
    variable = u
    boundary = 'right'
    v = v
  [../]
[]

[Functions]
  [./u_source_func]
    type = ParsedFunction
    value = '-4'
  [../]
  [./v_source_func]
    type = ParsedFunction
    value = '-6'
  [../]
  [./u_func]
    type = ParsedFunction
    value = '2*x^2 + 4*x'
  [../]
  [./v_func]
    type = ParsedFunction
    value = '3*x^2 + 2*x + 1'
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = PJFNK
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
[]
