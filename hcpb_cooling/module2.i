[Mesh]
  file = hcpb_homogenised.e
[]

[Variables]
  [T]
    initial_condition = 300.0
  []
[]

[AuxVariables]
  [heat_source]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  [initial_heat_source]
    type = ConstantIC
    variable = heat_source
    value = 0
  []
[]

[Functions]

  [helium_k_8_MPa]
    type = PiecewiseLinear
    xy_data = '
      280 0.154
      300 0.1611
      350 0.1784
      400 0.1950
      500 0.2267
      600 0.2566
      700 0.2851
      800 0.3125
      900 0.3329
    '
  []

  [helium_C_8_MPa]
    type = PiecewiseLinear
    xy_data = '
      280 5.189e3
      300 5.188e3
      350 5.187e3
      400 5.186e3
      500 5.187e3
      600 5.187e3
      700 5.187e3
      800 5.188e3
      900 5.189e3
    '
  []

  [helium_rho_8_MPa]
    type = PiecewiseLinear
    xy_data = '
      280 13.23
      300 12.38
      350 10.67
      400 9.380
      500 7.548
      600 6.314
      700 5.427
      800 4.758
      900 4.235
    '
  []

  [helium_k_200_KPa]
    type = PiecewiseLinear
    xy_data = '
      280 0.1488
      300 0.1561
      350 0.1736
      400 0.1904
      500 0.2224
      600 0.2525
      700 0.2811
      800 0.3086
      900 0.3351
    '
  []

  [helium_C_200_KPa]
    type = PiecewiseLinear
    xy_data = '
      280 5.193e3
      300 5.193e3
      350 5.193e3
      400 5.193e3
      500 5.193e3
      600 5.193e3
      700 5.193e3
      800 5.193e3
      900 5.193e3
    '
  []

  [helium_rho_200_KPa]
    type = PiecewiseLinear
    xy_data = '
      280 0.3435
      300 0.3206
      350 0.2749
      400 0.2405
      500 0.1925
      600 0.1604
      700 0.1375
      800 0.1203
      900 0.1070
    '
  []

  #[breeder_rho]
  #  type = ParsedFunction
  #  expression = '2.278e3 * exp(-3*(18.8*T + 0.5 * 1.66e-2 * (T)^2) * 1e-6)'
  #[]

  #[breeder_k]
  #  type = ParsedFunction
  #  expression = '(1.98 + 850 / T)* (0.498634 / (1 + (1 - 0.498634) * (2.14 - 7e-4*T)))'
  #[]

  #[breeder_c]
  #  type = ParsedFunction
  #  expression = '939.9 + 1.4577*(T) - 4.011e7/(T)^2'
  #[]

  #[multiplier_rho]
  #  type = ParsedFunction
  #  expression = '1.884e3 * exp(-3*(8.4305*(T) + 0.5 * 1.1464e-2 * (T)^2 - 1/3 * 2.9752e-6*(T)^3) * 1e-6)'
  #[]

  #[multiplier_k]
  #  type = ParsedFunction
  #  expression = '430.35 - 1.1674*(T) + 1.6044e-3*(T)^2 - 1.0097e-6*(T)^3 + 2.3642e-10*(T)^4'
  #[]

  #[multiplier_c]
  #  type = ParsedFunction
  #  expression = '606.91 - 5.3382*(T) - 4.1726e-3*(T)^2 + 1.2723e-6*(T)^3'
  #[]

  [fw_h]
    type = ConstantFunction
    value = 1.0e4
  []

  [fw_Tinf]
    type = ConstantFunction
    value = 300
  []

  [cp_h]
    type = ConstantFunction
    value = 1.0e4
  []

  [cp_Tinf]
    type = ConstantFunction
    value = 300
  []

  [bp_h]
    type = ConstantFunction
    value = 100
  []

  [bp_Tinf]
    type = ConstantFunction
    value = 300.0
  []

  [helium_8_MPa_rho]
    type = ParsedFunction
    expression = '2.989887e1 - 8.521529e-2*(x) + 1.029465e-04*(x)^2 - 4.444487e-08*(x)^3'
  []

  [tungsten_k]
    type = ParsedFunction
    expression = '240.51 - 0.2899*(x) + 2.5403e-4*(x)^2 + 1.0263e-7*(x)^3 + 1.5238e-11*(x)^4'
  []

  [tungsten_C]
    type = ParsedFunction
    expression = '116.37 - 7.1119e-2*(x) - 6.5828e-5*(x)^2 + 3.2396e-8*(x)^3 - 5.4523e-12*(x)^4'
  []

  [tungsten_rho]
    type = ParsedFunction
    expression = '19.3e3 * exp(-3*(5.0777*(x) + 0.5 * 5.6862e-4* (x)^2) * 1e-6)'
  []

  [eurofer_k]
    type = ParsedFunction
    expression = '22.867 - 1.4546e-2*(x-273.15) - 2.3056e-5*(x-273.15)^2 + 1.4815e-8*(x-273.15)^3'
  []

  [eurofer_C]
    type = ParsedFunction
    expression = '441.12 + 0.44049*(x - 273.15) -5.5848e-4*(x-273.15)^2 + 1.427e-6*(x-273.15)^3'
  []

  [eurofer_rho]
    type = ParsedFunction
    expression = '7.78e3 * exp(-3*(9.861 * (x - 273.15) + 0.5 * 6.6532e-3 * (x - 273.15)^2 + 1/3 * -3.7695e-6 * (x - 273.15)^3) * 1e-6)'
  []

  [breeder_k]
    type = ParsedFunction
    expression = '-1.165068e-09*(x)^3 + 2.684793e-06*(x)^2 - 2.035117e-03*(x) + 1.084675e+00'
  []

  [breeder_C]
    type = ParsedFunction
    expression = '1.912038e-11*(x)^5 - 6.511734e-08*(x)^4 + 8.888057e-05*(x)^3 - 6.142962e-02*(x)^2 + 2.343275e+01*(x) - 2.488073e+03'
  []

  [breeder_rho]
    type = ParsedFunction
    expression = '1.13638189e+03 - 6.59934269e-02*(x) - 2.28793141e-05*(x)^2 - 1.63595520e-09*(x)^3'
  []

  [multiplier_k]
    type = ParsedFunction
    expression = '1.219087e-10*(x)^4 - 5.206525e-07*(x)^3 + 8.272984e-04*(x)^2 - 6.018222e-01*(x) + 2.219560e+02'
  []

  [multiplier_C]
    type = ParsedFunction
    expression = '4.344442e-12*(x)^4 + 1.259689e-06*(x)^3 - 4.158392e-03*(x)^2 - 5.345847e+00*(x) + 6.092845e+02'
  []

  [multiplier_rho]
    type = ParsedFunction
    expression = '8.781650e-13*(x)^4 + 7.092810e-10*(x)^3 - 1.337280e-05*(x)^2 - 2.625565e-02*(x) + 9.720660e+02'
  []

  [dummy_k]
    type = ParsedFunction
    expression = '1e-10'
  []

  [dummy_C]
    type = ParsedFunction
    expression = '1e-2'
  []

  [dummy_rho]
    type = ParsedFunction
    expression = '1e-3'
  []
[]

[Materials]
  # all Eurofer regions

  [eurofer_conductivity_heat] 
    type = HeatConductionMaterial
    block = 'fw backplate'
    temp = T
    thermal_conductivity_temperature_function = eurofer_k
    specific_heat_temperature_function = eurofer_C
    min_T = 300
  []

  [eurofer_density]
    type = DerivativeParsedMaterial
    block = 'fw backplate'
    property_name = 'density'
    coupled_variables = 'T'
    expression = '7.78e3 * exp(-3*(9.861 * (T - 273.15) + 0.5 * 6.6532e-3 * (T - 273.15)^2 + 1/3 * -3.7695e-6 * (T - 273.15)^3) * 1e-6)'
    derivative_order = 1
  []

  # Tungtsen armour

  [tungsten_conductivity_heat] 
    type = HeatConductionMaterial
    block = 'armour'
    temp = T
    thermal_conductivity_temperature_function = tungsten_k
    specific_heat_temperature_function = tungsten_C
    min_T = 300
  []

  [tungsten_density] 
    type = DerivativeParsedMaterial
    block = 'armour'       
    property_name = 'density'
    coupled_variables = 'T'
    expression = '19.3e3 * exp(-3*(5.0777*(T) + 0.5 * 5.6862e-4* (T)^2) * 1e-6)'
    derivative_order = 1
  []

  # Helium coolant
  [helium_coolant_heat] 
    type = HeatConductionMaterial
    block = 'coolant_pipe_1 coolant_pipe_2 coolant_pipe_3 coolant_pipe_4 cp_bottom cp_top'
    temp = T
    thermal_conductivity_temperature_function = helium_k_8_MPa
    specific_heat_temperature_function = helium_C_8_MPa
    min_T = 300
  []

  [helium_coolant_density]
    type = DerivativeParsedMaterial
    block = 'coolant_pipe_1 coolant_pipe_2 coolant_pipe_3 coolant_pipe_4 cp_bottom cp_top'
    property_name = 'density'
    coupled_variables = 'T'
    expression = '2.989887e1 - 8.521529e-2*(T) + 1.029465e-04*(T)^2 - 4.444487e-08*(T)^3'
    derivative_order = 1
  []

  [breeder_conductivity_heat] 
    type = HeatConductionMaterial
    block = 'breeder'
    temp = T
    thermal_conductivity_temperature_function = breeder_k
    specific_heat_temperature_function = breeder_C
    min_T = 300
  []

  [breeder_density]
    type = DerivativeParsedMaterial
    block = 'breeder'
    property_name = 'density'
    coupled_variables = 'T'
    expression = '1.13638189e+03 - 6.59934269e-02*(T) - 2.28793141e-05*(T)^2 - 1.63595520e-09*(T)^3'
    derivative_order = 1
  []

  [multiplier_density]
    type = DerivativeParsedMaterial
    block = 'be_top be_bottom pre_slab_be_TETRA 34471'
    property_name = 'density'
    coupled_variables = 'T'
    expression = '8.781650e-13*(T)^4 + 7.092810e-10*(T)^3 - 1.337280e-05*(T)^2 - 2.625565e-02*(T) + 9.720660e+02'
    derivative_order = 1
  []

  [multiplier_conductivity_heat] 
    type = HeatConductionMaterial
    block = 'be_top be_bottom pre_slab_be_TETRA 34471'
    temp = T
    thermal_conductivity_temperature_function = multiplier_k
    specific_heat_temperature_function = multiplier_C
    min_T = 300
  []
[]


[Kernels]
  [time_derivative]
    type = HeatConductionTimeDerivative
    variable = T
  []

  [heat_conduction]
    type = HeatConduction
    variable = T
  []

  [heat_source_kernel]
    type = CoupledForce
    variable = T
    v = heat_source
  []
[]

[BCs]

  [source_heat]
    type = NeumannBC
    variable = T
    boundary = 'fw_front'
    value = 1e3
  []

  [fw_cooling]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'coolant_pipe_int_1 coolant_pipe_int_2 coolant_pipe_int_3 coolant_pipe_int_4'
    coefficient = fw_h
    T_infinity = fw_Tinf
  []

  [coolant_plate_cooling]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'bottom_be_cp_int top_be_cp_int bottom_cp_bp_int bottom_cp_bb_int bottom_cp_slab_int top_cp_bb_int top_cp_bp_int top_cp_slab_int'
    coefficient = cp_h
    T_infinity = cp_Tinf
  []

  [backplate_cooling]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'backplate_rear'
    coefficient = bp_h
    T_infinity = bp_Tinf
  []
[]

[Postprocessors]
  [source_integral]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    execute_on = transfer
  []
[]

[MultiApps]
  [openmc]
    type = TransientMultiApp
    input_files = 'openmc_cooling.i'
    execute_on = timestep_begin
  []
[]


[Transfers]
  [heat_source_from_openmc]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = openmc
    variable = heat_source
    source_variable = heating_local
  []

  [temp_to_openmc]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = openmc
    variable = temp
    source_variable = T
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-20
  l_tol = 1e-3
  nl_max_its = 10
  end_time = 10
  dt = 1
  accept_on_max_fixed_point_iteration = true
  steady_state_detection = true
  steady_state_tolerance = 1e-6
[]

[Outputs]
  exodus = true
  csv = true
  problem_summary = true
[]
