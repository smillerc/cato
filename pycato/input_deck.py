# -*- coding: utf-8 -*-
class CatoInput(object):
    def __init__(self, ini_file: None):
        # general
        self.title = None  # Name of the simulation
        self.unit_system = None  # Type of units to output (CGS, ICF, MKS, etc...)

        # reference state
        self.reference_pressure = 1.0
        self.reference_density = 1.0
        self.reference_mach = 1.0

        # grid
        self.grid_type = "2d_regular"  # Structure/layout of the grid, e.g. '2d_regular'
        self.xmin = 0.0
        # Minimum extent of the grid in x (ignored for .h5 initial grids)
        self.xmax = 0.0
        # Maximum extent of the grid in x (ignored for .h5 initial grids)
        self.ymin = 0.0
        # Minimum extent of the grid in y (ignored for .h5 initial grids)
        self.ymax = 0.0
        # Maximum extent of the grid in y (ignored for .h5 initial grids)
        self.ni_nodes = 0
        # # of i nodes (not including ghost) (ignored for .h5 initial grids)
        self.nj_nodes = 0
        # # of j nodes (not including ghost) (ignored for .h5 initial grids)
        self.n_ghost_layers = 2
        # # of ghost layers to use (based on spatial reconstruction order)

        # initial conditions
        self.initial_condition_file = None
        self.read_init_cond_from_file = False

        self.init_x_velocity = 0.0
        self.init_y_velocity = 0.0
        self.init_density = 0.0
        self.init_pressure = 0.0

        # restart files
        self.restart_from_file = False
        self.restart_file = None

        # boundary conditions
        self.bc_pressure_input_file = None
        self.apply_constant_bc_pressure = False
        self.constant_bc_pressure_value = 0.0
        self.bc_density = 5e-3  # density input for a pressure boundary condition
        self.bc_pressure_scale_factor = 1.0
        self.plus_x_bc = "periodic"  # Boundary condition at +x
        self.minus_x_bc = "periodic"  # Boundary condition at -x
        self.plus_y_bc = "periodic"  # Boundary condition at +y
        self.minus_y_bc = "periodic"  # Boundary condition at -y

        # source inputs (i.e. energy source injection at a grid location)
        self.enable_source_terms = False
        self.source_term_type = "energy"
        self.apply_constant_source = False
        self.source_file = None
        self.constant_source_value = 0.0
        self.source_scale_factor = 1.0
        self.source_ilo = 0
        self.source_ihi = 0
        self.source_jlo = 0
        self.source_jhi = 0

        # io
        self.contour_io_format = "xdmf"
        self.append_date_to_result_folder = False
        self.plot_reconstruction_states = False
        self.plot_reference_states = False
        self.plot_evolved_states = False
        self.plot_64bit = True
        self.plot_ghost_cells = True

        # timing
        self.max_time = 1.0
        self.cfl = 0.1  # Courant–Friedrichs–Lewy condition
        self.initial_delta_t = 0.0
        self.use_constant_delta_t = False
        self.contour_interval_dt = 0.5
        self.max_iterations = huge(1)
        self.time_integration_strategy  # How is time integration handled? e.g. 'rk2', 'rk4', etc.
        self.smooth_residuals = True

        # physics
        self.polytropic_index = 5.0 / 3.0  # e.g. gamma for the simulated gas

        # finite volume scheme specifics
        self.flux_solver = "AUSMPW+"  # flux solver, ('FVLEG', 'AUSM+-up', etc.)
        self.spatial_reconstruction = (
            "MUSCL"  # (MUSCL, green_gauss) How are the edge values interpolated?
        )
        self.apply_low_mach_fix = True  # some flux solvers have this option
        self.limiter = (
            "minmod"  # Flux limiter, e.g. minmod, superbee, TVD3, MLP3, MLP5, etc.
        )

        # AUSM solver specifics, see the AUSM solver packages for more details.
        # These are only read in if the flux solver is from the AUSM family
        self.ausm_beta = 1.0 / 8.0  # beta parameter
        self.ausm_pressure_diffusion_coeff = 0.25  # K_p; pressure diffusion coefficient
        self.ausm_pressure_flux_coeff = 0.75  # K_u; pressure flux coefficient
        self.ausm_sonic_point_sigma = (
            1.0  # sigma; another pressure diffusion coefficient
        )

        self.tau = 1.0e-5  # time increment for FVEG and FVLEG schemes

        # debug
        self.plot_limiters = False
        self.plot_gradients = False
