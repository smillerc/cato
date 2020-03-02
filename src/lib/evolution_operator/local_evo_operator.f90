module mod_local_evo_operator
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, &
                                                                              std_err => error_unit, std_out => output_unit
  use mod_globals, only: debug_print, print_recon_data

  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: near_zero, equal
  use math_constants, only: pi, rad2deg
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
  use mod_cone, only: cone_t, new_cone
  use mod_grid, only: grid_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_reconstruction_factory, only: reconstruction_factory

  implicit none

  private
  public :: local_evo_operator_t

  type, extends(abstract_evo_operator_t) :: local_evo_operator_t
    !< Local Evolution Operator (E0)
  contains
    procedure, public :: initialize
    procedure, public :: evolve_leftright_midpoints
    procedure, public :: evolve_downup_midpoints
    procedure, public :: evolve_corners
    procedure, public :: e0_operator
    procedure, public :: copy
    final :: finalize
  end type local_evo_operator_t

contains

  subroutine initialize(self, input, grid_target, recon_operator_target)
    !< Constructor for the FVLEG operator
    class(local_evo_operator_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    class(grid_t), intent(in), target :: grid_target
    class(abstract_reconstruction_t), intent(in), target :: recon_operator_target

    self%name = "FVLEG"
    call self%set_tau(input%tau)
    self%grid => grid_target
    self%reconstruction_operator => recon_operator_target

  end subroutine

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self

    call debug_print('Running local_evo_operator_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%name)) deallocate(self%name)
    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%reconstructed_state)) nullify(self%reconstructed_state)
    if(associated(self%reconstruction_operator)) nullify(self%reconstruction_operator)
  end subroutine finalize

  subroutine copy(out_evo, in_evo)
    class(abstract_evo_operator_t), intent(in) :: in_evo
    class(local_evo_operator_t), intent(inout) :: out_evo

    call debug_print('Running local_evo_operator_t%copy()', __FILE__, __LINE__)

    if(associated(out_evo%grid)) nullify(out_evo%grid)
    out_evo%grid => in_evo%grid

    if(associated(out_evo%reconstructed_state)) nullify(out_evo%reconstructed_state)
    out_evo%reconstructed_state => in_evo%reconstructed_state

    if(associated(out_evo%reconstruction_operator)) nullify(out_evo%reconstruction_operator)
    out_evo%reconstruction_operator => in_evo%reconstruction_operator

    if(allocated(out_evo%name)) deallocate(out_evo%name)
    allocate(out_evo%name, source=in_evo%name)

  end subroutine

  subroutine evolve_leftright_midpoints(self, evolved_state, lbounds, error_code)
    !< Create Mach cones and evolve the state at all of the left/right midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the left and right.

    class(local_evo_operator_t), intent(in) :: self

    integer(ik), dimension(3), intent(in) :: lbounds
    integer(ik), intent(out) :: error_code
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(out) :: evolved_state
    !< ((rho, u, v, p), i, j); Reconstructed U at each midpoint on the left/right edges

    ! Locals
    type(cone_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik), parameter :: midpoint_idx = 2 !< used to select the midpoint in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the edge in the reconstructed_state
    integer(ik), dimension(2, 2) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 2) :: midpoint_edge_vectors
    !< ((x,y), head/tail, vector_id); set of vectors that define the midpoint

    real(rk), dimension(4, 2) :: reconstructed_midpoint_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the midpoint with respect to each cell

    integer(ik) :: ilo, ihi, jlo, jhi
    real(rk), dimension(4) :: evolved_primitive_vars

    error_code = 0

    !       cell 1
    !      (i-1,j)
    !  N4----M3----N3
    !  |            |
    !  M4    C1     M2
    !  |            |
    !  N1----M1----N2
    !
    !  P1----O----P2  edge vector
    !
    !  N4----M3----N3
    !  |           |
    !  M4    C2    M2
    !  |           |
    !  N1----M1----N2
    !       cell 2
    !      (i,j-1)
    !
    ! For left/right midpoints, the edge vectors go left (O-P1) then right (O-P2).
    ! The neighboring cells are up (i, j+1) and down (i, j-1)

    ! We only reconstruct the real domain
    ilo = lbound(evolved_state, dim=2) ! in the i-direction, # left/right midpoints = # cells
    ihi = ubound(evolved_state, dim=2) ! in the i-direction, # left/right midpoints = # cells
    jlo = lbound(evolved_state, dim=3) ! in the j-direction, # left/right midpoints = # nodes
    jhi = ubound(evolved_state, dim=3) ! in the j-direction, # left/right midpoints = # nodes

    do j = jlo, jhi
      do i = ilo, ihi
        ! print*, i, j
        neighbor_cell_indices(:, 1) = [i, j]     ! above
        neighbor_cell_indices(:, 2) = [i, j - 1] ! below

        ! <----M---->  (left and right vectors) these have (x_tail, y_tail) and (x_head, y_head)
        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='bottom')

        ! reconstructed_state indexing
        !< ((rho, u ,v, p), point, node/midpoint, i, j);

        ! Cell 1 -> use M1 from the cell above
        point_idx = 1
        reconstructed_midpoint_state(:, 1) = self%reconstructed_state(:, point_idx, midpoint_idx, i, j)

        ! Cell 2 -> use M3 from the cell below
        point_idx = 3
        reconstructed_midpoint_state(:, 2) = self%reconstructed_state(:, point_idx, midpoint_idx, i, j - 1)

        ! if (minval(reconstructed_midpoint_state(4,:)) < 0.0_rk) then
        !   call print_recon_data('p', i, j, self%reconstructed_state)
        ! end if

        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             cell_indices=neighbor_cell_indices, &
                             cone_location='left/right midpoint')

        ! Set the evolved state at the midpoint
        call self%e0_operator(mach_cone, evolved_primitive_vars, 'evolve_leftright_midpoints', error_code)
        ! if(i == 601 .and. j == 1) then
        !   print *, mach_cone
        !   print *, 'evolved_primitive_vars (601,1): ', evolved_primitive_vars
        ! end if
        if(error_code /= 0) then
          write(std_out, '(a)') 'The E0 operator returned an error in '// &
            'local_evo_operator_t%evolve_leftright_midpoints(), check cato.error for details'
          write(std_err, '(2(a,i0))') 'Error in local_evo_operator_t%evolve_leftright_midpoints() at i=', i, ', j=', j
          write(std_err, *) mach_cone
          return
        end if
        evolved_state(:, i, j) = evolved_primitive_vars
        ! evolved_state(:, i, j) = self%e0_operator(mach_cone)

      end do
    end do
  end subroutine evolve_leftright_midpoints

  subroutine evolve_downup_midpoints(self, evolved_state, lbounds, error_code)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    integer(ik), dimension(3), intent(in) :: lbounds
    integer(ik), intent(out) :: error_code
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(out) :: evolved_state
    !< ((rho, u, v, p), i, j); Evolved U at each midpoint on the left/right edges

    ! Locals
    type(cone_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik), parameter :: midpoint_idx = 2 !< used to select the midpoint in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the edge in the reconstructed_state
    integer(ik), dimension(2, 2) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 2) :: midpoint_edge_vectors !< ((x,y), (tail,head), (vector1:vector2))

    real(rk), dimension(4, 2) :: reconstructed_midpoint_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the midpoint with respect to each cell

    integer(ik) :: ilo, ihi, jlo, jhi
    real(rk), dimension(4) :: evolved_primitive_vars

    error_code = 0

    !       cell 1          edge          cell 2
    !      (i-1,j)         vector          (i,j)
    !  N4----M3----N3       P2       N4----M3----N3
    !  |           |        |       |            |
    !  M4    C1    M2       O       M4    C2     M2
    !  |           |        |       |            |
    !  N1----M1----N2      P1       N1----M1----N2
    !
    ! For down/up midpoints, the edge vectors go down (O-P1) then up (O-P2).
    ! The neighboring cells are left (i, j) and right (i-1, j)

    ilo = lbound(evolved_state, dim=2)  ! in the i-direction, # down/up midpoints = # nodes
    ihi = ubound(evolved_state, dim=2)  ! in the i-direction, # down/up midpoints = # nodes
    jlo = lbound(evolved_state, dim=3)  ! in the j-direction, # down/up midpoints = # cells
    jhi = ubound(evolved_state, dim=3)  ! in the j-direction, # down/up midpoints = # cells

    do j = jlo, jhi
      do i = ilo, ihi
        ! print*, i,j

        ! cell ordering is 1) left, 2) right
        neighbor_cell_indices(:, 1) = [i - 1, j] ! left
        neighbor_cell_indices(:, 2) = [i, j]     ! right

        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='left')

        ! Cell 1: use M2 from the left cell
        point_idx = 2
        reconstructed_midpoint_state(:, 1) = self%reconstructed_state(:, point_idx, midpoint_idx, i - 1, j)

        ! Cell 2: cell to the right -> use M4 from the right cell
        point_idx = 4
        reconstructed_midpoint_state(:, 2) = self%reconstructed_state(:, point_idx, midpoint_idx, i, j)
        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             cell_indices=neighbor_cell_indices, &
                             cone_location='down/up midpoint')

        ! Set the evolved state at the midpoint
        call self%e0_operator(mach_cone, evolved_primitive_vars, 'evolve_downup_midpoints', error_code)
        ! if(i == 200 .and. j == 1) then
        !   print *, mach_cone
        !   print *, 'evolved_primitive_vars down/up (601,1): ', evolved_primitive_vars
        ! !   ! error stop
        ! end if
        if(error_code /= 0) then
          write(std_out, '(a)') 'The E0 operator returned an error in '// &
            'local_evo_operator_t%evolve_downup_midpoints(), check cato.error for details'
          write(std_err, '(2(a,i0))') 'Error in local_evo_operator_t%evolve_downup_midpoints() at i=', i, ', j=', j
          write(std_err, *) mach_cone
          return
        end if
        evolved_state(:, i, j) = evolved_primitive_vars
        ! evolved_state(:, i, j) = self%e0_operator(mach_cone)
      end do
    end do
  end subroutine evolve_downup_midpoints

  subroutine evolve_corners(self, evolved_state, lbounds, error_code)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    integer(ik), dimension(3), intent(in) :: lbounds
    integer(ik), intent(out) :: error_code
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(out) :: evolved_state
    !< ((rho, u, v, p), i, j); Reconstructed U at each corner

    ! Locals
    type(cone_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik), parameter :: corner_idx = 1 !< used to select the corner in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the point in the reconstructed_state
    integer(ik), dimension(2, 4) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 4) :: corner_edge_vectors !< ((x,y), (tail,head), (vector1:vector4))
    !< ((x,y), head/tail, vector_id); set of vectors that define the corner

    real(rk), dimension(4, 4) :: reconstructed_corner_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the corner with respect to each cell

    integer(ik) :: ilo, ihi, jlo, jhi
    real(rk), dimension(4) :: evolved_primitive_vars

    error_code = 0

    !       cell 4                   cell 3
    !      (i-1,j)                   (i,j)
    !  N4----M3----N3    P3   N4----M3----N3
    !  |            |    |    |            |
    !  M4    C4    M2    |    M4    C3    M2
    !  |            |    |    |            |
    !  N1----M1----N2    |    N1----M1----N2
    !                    |
    !  P4----------------O-----------------P2
    !                    |
    !  N4----M3----N3    |    N4----M3----N3
    !  |            |    |    |            |
    !  M4    C1    M2    |    M4    C2    M2
    !  |            |    |    |            |
    !  N1----M1----N2    P1   N1----M1----N2
    !      cell 1                  cell 2
    !     (i-1,j-1)               (i,j-1)

    ! For corners, the edge vectors go down (0-P1), right(O-P2), up(O-P3), right (O-P4)
    ! O is the corner point shared by the 4 cells

    ! Unlike the midpoints, the corner ref state dimensions are the same in i,j as the # nodes
    ilo = lbound(evolved_state, dim=2)
    ihi = ubound(evolved_state, dim=2)
    jlo = lbound(evolved_state, dim=3)
    jhi = ubound(evolved_state, dim=3)

    do j = jlo, jhi
      do i = ilo, ihi
        ! cell ordering is 1) lower left, 2) lower right, 3) upper right, 4) upper left
        neighbor_cell_indices(:, 1) = [i - 1, j - 1] ! lower left
        neighbor_cell_indices(:, 2) = [i, j - 1]     ! lower right
        neighbor_cell_indices(:, 3) = [i, j]         ! upper right
        neighbor_cell_indices(:, 4) = [i - 1, j]     ! upper left

        corner_edge_vectors = self%grid%get_corner_vectors(cell_ij=[i, j], corner='lower-left')

        ! reconstructed_state indexing
        !< ((rho, u ,v, p), point, node/midpoint, i, j);

        ! Cell 1: lower left cell -> corner is in the upper right (N3) of its parent cell
        reconstructed_corner_state(:, 1) = self%reconstructed_state(:, 3, corner_idx, i - 1, j - 1)

        ! Cell 2: lower right cell -> corner is in the upper left (N4) of its parent cell
        reconstructed_corner_state(:, 2) = self%reconstructed_state(:, 4, corner_idx, i, j - 1)

        ! Cell 3: upper right cell-> corner is in the lower left (N1) of its parent cell
        reconstructed_corner_state(:, 3) = self%reconstructed_state(:, 1, corner_idx, i, j)

        ! Cell 4: upper left cell -> corner is in the lower right (N2) of its parent cell
        reconstructed_corner_state(:, 4) = self%reconstructed_state(:, 2, corner_idx, i - 1, j)

        mach_cone = new_cone(tau=self%tau, edge_vectors=corner_edge_vectors, &
                             reconstructed_state=reconstructed_corner_state, &
                             cell_indices=neighbor_cell_indices, &
                             cone_location='corner')

        ! Set the evolved state at the midpoint
        call self%e0_operator(mach_cone, evolved_primitive_vars, 'evolve_corners', error_code)
        if(error_code /= 0) then
          write(std_out, '(a)') 'The E0 operator returned an error in '// &
            'local_evo_operator_t%evolve_corners(), check cato.error for details'
          write(std_err, '(2(a,i0))') 'Error in local_evo_operator_t%evolve_corners() at i=', i, ', j=', j
          write(std_err, *) mach_cone
          return
        end if
        evolved_state(:, i, j) = evolved_primitive_vars
      end do
    end do
  end subroutine evolve_corners

  subroutine e0_operator(self, mach_cone, primitive_variables, calling_routine, error_code)
    use mod_eos, only: eos
    class(local_evo_operator_t), intent(in) :: self
    class(cone_t), intent(in) :: mach_cone
    character(len=*), intent(in) :: calling_routine
    integer(ik), intent(out) :: error_code
    real(rk), dimension(4), intent(out) :: primitive_variables

    real(rk), dimension(4) :: p_prime_prim_vars !< [rho, u, v, p] at P'
    real(rk) :: pressure, density, u, v
    real(rk) :: eps
    integer(ik) :: i

    error_code = 0 ! All ok at first

    do i = 1, 4
      if(mach_cone%p_prime_in_cell(i)) then
        p_prime_prim_vars = self%reconstruction_operator%reconstruct_point( &
                            xy=mach_cone%p_prime_xy, &
                            cell_ij=mach_cone%p_prime_ij(:, i))
      end if
    end do

    associate(theta_ie=>mach_cone%theta_ie, &
              theta_ib=>mach_cone%theta_ib, &
              rho_tilde=>mach_cone%reference_density, &
              a_tilde=>mach_cone%reference_sound_speed, &
              rho_p_prime=>p_prime_prim_vars(1), &
              pressure_p_prime=>p_prime_prim_vars(4), &
              u_i=>mach_cone%arc_primitive_vars(2, :, :), &
              v_i=>mach_cone%arc_primitive_vars(3, :, :), &
              p_i=>mach_cone%arc_primitive_vars(4, :, :))

      pressure = sum(p_i * (theta_ie - theta_ib) &
                     - rho_tilde * a_tilde * u_i * (sin(theta_ie) - sin(theta_ib)) &
                     + rho_tilde * a_tilde * v_i * (cos(theta_ie) - cos(theta_ib))) / (2.0_rk * pi)

      density = rho_p_prime + ((pressure - pressure_p_prime) / a_tilde**2)
      ! if(abs(pressure - pressure_p_prime) > 0.0_rk) then
      !   ! print*, 'supersonic!'
      !   print *, "rho(P')", rho_p_prime
      !   print *, "rho(P)", density
      !   print *, "p(P)", pressure
      !   print *, "p(P')", pressure_p_prime
      !   print *, "p(P) - p(P')", pressure - pressure_p_prime
      !   print *, "p_i / (rho_tilde * a_tilde)", p_i / (rho_tilde * a_tilde)
      !   print *
      !   ! error stop
      ! end if

      u = sum((-p_i / (rho_tilde * a_tilde)) * (sin(theta_ie) - sin(theta_ib)) &
              + u_i * ( &
              ((theta_ie - theta_ib) / 2.0_rk) &
              + ((sin(2.0_rk * theta_ie) - sin(2.0_rk * theta_ib)) / 4.0_rk) &
              ) &
              - v_i * ((cos(2.0_rk * theta_ie) - cos(2.0_rk * theta_ib)) / 4.0_rk) &
              ) / pi

      v = sum((p_i / (rho_tilde * a_tilde)) * (cos(theta_ie) - cos(theta_ib)) &
              - u_i * ((cos(2.0_rk * theta_ie) - cos(2.0_rk * theta_ib)) / 4.0_rk) &
              + v_i * ( &
              ((theta_ie - theta_ib) / 2.0_rk) &
              - ((sin(2.0_rk * theta_ie) - sin(2.0_rk * theta_ib)) / 4.0_rk) &
              ) &
              ) / pi

      if(abs(u) < 1.0e-15_rk) u = 0.0_rk
      if(abs(v) < 1.0e-15_rk) v = 0.0_rk
      primitive_variables = [density, u, v, pressure]

      if(density < 0.0_rk .or. pressure < 0.0_rk) then

        print *, "Density or pressure is negative, called by: ", trim(calling_routine)
        write(*, *) mach_cone
        write(*, *) 'radius', mach_cone%radius
        print *, 'Primitive variables ', primitive_variables

        write(std_err, '(a)') 'Error in local_evo_operator_t%e0_operator():'
        write(std_err, '(a)') 'This may just be a problem during a stage in '// &
          'time-integration, but just to be explicit, here it is...'
        write(std_err, '(a, es15.8)') 'Density(P):      ', density
        write(std_err, '(a, es15.8)') "Density(P'):     ", rho_p_prime
        write(std_err, '(a, es15.8)') "p(P)/a_tilde^2:  ", pressure / a_tilde**2
        write(std_err, '(a, es15.8)') "Pressure(P'):    ", pressure_p_prime
        write(std_err, '(a, es15.8)') 'Pressure(P):     ', pressure
        write(std_err, '(a, es15.8)') "a_tilde:         ", a_tilde
        write(std_err, '(a, es15.8)') "p(P) - p(P')/a^2:",((pressure - pressure_p_prime) / a_tilde**2)
        write(std_err, '(a, es15.8)') 'X Velocity:      ', u
        write(std_err, '(a, es15.8)') 'Y Velocity:      ', v
        write(std_err, '(a, es15.8)') '=========================='
        error stop "Density or pressure is negative"
      end if

    end associate
  end subroutine e0_operator

end module mod_local_evo_operator
