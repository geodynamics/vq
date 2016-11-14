sim.time.end_year                 = 1000
sim.system.progress_period        = 5
sim.greens.use_normal             = true


# -===== FILE I/O =========
# == If calculating Green's functions for this model for the first time, use
sim.greens.method                 = standard
sim.greens.output                  = /path/to/greens/functions.h5

# == If using previously-calculated Green's functions, use ==
# sim.greens.method               = file
# sim.greens.input                = /path/to/greens/functions.h5

sim.file.input                    = /path/to/fault/model/file.h5
sim.file.input_type               = hdf5
sim.file.output_event             = /path/to/event/save/file.h5
sim.file.output_event_type        = hdf5


# -===== FRICTION PARAMETERS =========
sim.friction.dynamic              = 0.2
sim.friction.dynamic_stress_drops = true


# -===== UCERF3 GREENS LIMITS =========
sim.greens.shear_offdiag_min      = -745348
sim.greens.shear_offdiag_max      =  3590676
sim.greens.normal_offdiag_min     = -197265.50322689
sim.greens.normal_offdiag_max     =  202035.55337244


# -===== UCERF2 GREENS LIMITS =========
# sim.greens.shear_offdiag_min      = -745348
# sim.greens.shear_offdiag_max      =  3590676
# sim.greens.normal_offdiag_min     = -837457
# sim.greens.normal_offdiag_max     =  953668


# -===== Other possible parameters (with default values where given) =========
# sim.bass.max_generations = 0
# sim.bass.b = 1
# sim.bass.c = 0.1
# sim.bass.d = 300
# sim.bass.dm = 1.25
# sim.bass.max_generations = 0
# sim.bass.mm = 4
# sim.bass.p = 1.25
# sim.bass.q = 1.35
# sim.file.output_stress =
# sim.file.output_stress_index =
# sim.file.output_stress_type =
# sim.greens.bh_theta = 0
# sim.greens.input =
# sim.greens.kill_distance = 0
# sim.greens.output =
# sim.greens.sample_distance = 1000
# sim.system.sanity_check = false
# sim.system.checkpoint_period = 0
# sim.system.checkpoint_prefix = sim_state_
# sim.system.transpose_matrix = 1
