# cardinal_hcpb_fixed_source.i
# Cardinal OpenMC cell-average problem configured for an HCPB-like blanket.
# Mesh file should exist: hcpb_homogenised.e (created from your VTK -> Exodus step)
# MultiApp "conduction" should be a MOOSE input file (module.i) that solves heat conduction.
#
# Notes:
# - Tallies: H3_production, heating_local, flux
# - Uses fixed-source normalization (source_strength) — suitable for fusion D-T neutrons
# - temperature_blocks must match Exodus element-block names (see `ncdump`/meshio output)

[Mesh]
  [file]
    type = FileMeshGenerator
    file = hcpb_homogenised.e
  []
[]

#######################
# Temperature aux var  #
#######################
[AuxVariables]
  [cell_temperature]
    family = MONOMIAL
    order = CONSTANT
  []
  [cell_density]
    family = MONOMIAL
    order = CONSTANT
  []
  [material_id]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [cell_temperature]
    type = CellTemperatureAux
    variable = cell_temperature
  []

  [cell_density]
    type = CellDensityAux
    variable = cell_density
  []

  [material_id]
    type = CellMaterialIDAux
    variable = material_id
  []  
[]


########################
# OpenMC / Cardinal Job #
########################
[Problem]
  type = OpenMCCellAverageProblem
  verbose = true
  source_strength = 2e18
  xml_directory = model.xml
  # scaling factor applied to OpenMC-derived tallies before transfer (optional)
  initial_properties = xml
  temperature_blocks = 'fw armour coolant_pipe_1 coolant_pipe_2 coolant_pipe_3 coolant_pipe_4 pre_slab_be_TETRA be_top cp_top breeder cp_bottom be_bottom backplate 34471'
  density_blocks = 'fw armour coolant_pipe_1 coolant_pipe_2 coolant_pipe_3 coolant_pipe_4 pre_slab_be_TETRA be_top cp_top breeder cp_bottom be_bottom backplate 34471'
  # hierarchy level for cell averaging (0 => element block level; 1 => cell/universe level)
  cell_level = 0

  particles = 1000

  # turn off relaxation while you debug (relaxation cant be used if the mapping changes)
  relaxation = none

  # Skinner disabled here: enable only if MoabSkinner is registered in your build.
  #skinner = moab

  [Tallies]
    [hcpb_mesh]
      type = MeshTally
      mesh_template = hcpb_homogenised.e
      score = 'H3_production heating_local flux'
      normalize_by_global_tally = true
      output = unrelaxed_tally_std_dev
    []
  []
[]

#####################
# Skinner / Skinsmap #
#####################
# NOTE: MoabSkinner is commented out because it was not registered in your Cardinal build.
# Uncomment and configure only if MoabSkinner is available.
#[UserObjects]
 # [moab]
 #   type = MoabSkinner
 #   temperature = temp
 #   temperature_min = 300
 #   temperature_max = 2000
 #   n_temperature_bins = 40
#    build_graveyard = true
#    output_skins = true
##  []
#[]

#####################
# Postprocessors    #
#####################
[Postprocessors]
  [heating_integral]
    type = ElementIntegralVariablePostprocessor
    variable = heating_local
  []
  [tritium_integral]
    type = ElementIntegralVariablePostprocessor
    variable = H3_production
  []
  [flux_integral]
    type = ElementIntegralVariablePostprocessor
    variable = flux
  []
  [tritium_error]
    type = TallyRelativeError
    tally_score = H3_production
    value_type = average
  []
  [heating_error]
    type = TallyRelativeError
    tally_score = heating_local
    value_type = average
  []
  [flux_error]
    type = TallyRelativeError
    tally_score = flux
    value_type = average
  []
[]


####################
# Initial / ICs    #
####################
[ICs]
  [H3_init]
    type = ConstantIC
    variable = H3_production
    value = 0.0
  []

  [heat_init]
    type = ConstantIC
    variable = heating_local
    value = 0
  []

  [flux_init]
    type = ConstantIC
    variable = flux
    value = 0
  []

  [temp_init]
    type = ConstantIC
    variable = temp
    value = 300
  []
  
  [eurofer_init]
    type = ConstantIC
    variable = density
    value = 7.78e3
    block = 'fw backplate'
  []

  [tungsten_init]
    type = ConstantIC
    variable = density
    value = 19.3e3
    block = 'armour'
  []

  [coolant_init]
    type = ConstantIC
    variable = density
    value = 6.3114002161460885
    block = 'coolant_pipe_1 coolant_pipe_2 coolant_pipe_3 coolant_pipe_4 cp_bottom cp_top'
  []

  [breeder_init]
    type = ConstantIC
    variable = density
    value = 1.136e3
    block = 'breeder'
  []

  [multiplier_init]
    type = ConstantIC
    variable = density
    value = 0.9721144468361395e3
    block = 'be_top be_bottom pre_slab_be_TETRA 34471'
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