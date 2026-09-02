!----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Provides an implementation of the Psy layer for physics

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_phys_mod

  use constants_mod,         only : i_def, r_def
  use field_mod,             only : field_type, field_proxy_type
  use integer_field_mod,     only : integer_field_type, integer_field_proxy_type
  use mesh_mod,              only : mesh_type
  
  implicit none
  public

contains
  !---------------------------------------------------------------------
  !> LFRic and PSyclone currently do not have a mechanism to loop over a subset
  !> of cells in a horizontal domain. Furthermore this is a psykal_lite file which is not PSycloned.
  !> PSyclone development to this aim can be tracked at https://github.com/stfc/PSyclone/issues/3193
  !> Future refactoring may seek to push these loops into a kernel layer and then to apply PSyclone.
  !>
  !> The orographic drag kernel only needs to be applied to a subset of land
  !> points where the standard deviation of subgrid orography is more than zero.
  !>
  !> invoke_orographic_drag_kernel: Invokes the kernel which calls the UM
  !> orographic drag scheme only on points where the standard deviation
  !> of subgrid orography is more than zero.
  subroutine invoke_orographic_drag_kernel(                              &
                      du_orog_blk, dv_orog_blk, du_orog_gwd, dv_orog_gwd,&
                      dtemp_orog_blk, dtemp_orog_gwd, u_in_w3, v_in_w3,&
                      wetrho_in_w3, theta_in_wth, exner_in_w3, sd_orog,  &
                      grad_xx_orog, grad_xy_orog, grad_yy_orog,          &
                      mr_v, mr_cl, mr_ci,                                &
                      height_w3, height_wth,                             &
                      taux_orog_blk, tauy_orog_blk,                      &
                      taux_orog_gwd, tauy_orog_gwd )

    use orographic_drag_kernel_mod, only: orographic_drag_kernel_code
    use mesh_mod, only: mesh_type
    use physics_config_mod,  only : gw_segment
    implicit none

    ! Increments from orographic drag
    type(field_type), intent(inout) :: du_orog_blk, dv_orog_blk, & ! Winds
                                       du_orog_gwd, dv_orog_gwd, &
                                       dtemp_orog_blk, dtemp_orog_gwd ! Temperature

    ! Inputs to orographic drag scheme
    type(field_type), intent(in) :: u_in_w3, v_in_w3,           & ! Winds
                                    wetrho_in_w3, theta_in_wth, & ! Density, Temperature
                                    exner_in_w3,                & ! Exner pressure
                                    sd_orog, grad_xx_orog,      & ! Orography ancils
                                    grad_xy_orog, grad_yy_orog, & !
                                    mr_v, mr_cl, mr_ci,         & ! mixing ratios
                                    height_w3, height_wth         ! Heights
    ! Diagnostics from orographic drag
    ! Stress from ...
    type(field_type), intent(inout) :: taux_orog_blk, tauy_orog_blk, & ! ... orographic flow blocking drag
                                       taux_orog_gwd, tauy_orog_gwd    ! ... orographic gravity wave drag

    integer :: cell

    ! Number of degrees of freedom
    integer :: ndf_w3, undf_w3, ndf_wtheta, undf_wtheta

    ! Integers for segmentation
    integer :: applicable_points, nlayers, loop_upper_bound, segment, seg_len, &
               l_bound, u_bound, n_segments , seg_target
    
    ! These are in ANY_DISCONTINUOUS_SPACE_1
    integer :: ndf_adspc1_sd_orog, undf_adspc1_sd_orog

    integer, allocatable ::    cell_index(:) !dhc record the points with orography

    
    type(field_proxy_type) :: du_blk_proxy, dv_blk_proxy,             &
                              du_orog_gwd_proxy, dv_orog_gwd_proxy,   &
                              dtemp_blk_proxy, dtemp_orog_gwd_proxy,  &
                              u_in_w3_proxy, v_in_w3_proxy,           &
                              wetrho_in_w3_proxy, theta_in_wth_proxy, &
                              exner_in_w3_proxy,                      &
                              sd_orog_proxy, grad_xx_orog_proxy,      &
                              grad_xy_orog_proxy, grad_yy_orog_proxy, &
                              mr_v_proxy, mr_cl_proxy, mr_ci_proxy,   &
                              height_w3_proxy, height_wth_proxy,      &
                              taux_blk_proxy, tauy_blk_proxy,         &
                              taux_gwd_proxy, tauy_gwd_proxy

    integer, pointer :: map_adspc1_sd_orog(:,:) => null(), &
                        map_w3(:,:) => null(),                  &
                        map_wtheta(:,:) => null()

    TYPE(mesh_type), pointer :: mesh => null()

    ! Initialise field and/or operator proxies
    du_blk_proxy = du_orog_blk%get_proxy()
    dv_blk_proxy = dv_orog_blk%get_proxy()
    du_orog_gwd_proxy = du_orog_gwd%get_proxy()
    dv_orog_gwd_proxy = dv_orog_gwd%get_proxy()
    dtemp_blk_proxy = dtemp_orog_blk%get_proxy()
    dtemp_orog_gwd_proxy = dtemp_orog_gwd%get_proxy()
    u_in_w3_proxy = u_in_w3%get_proxy()
    v_in_w3_proxy = v_in_w3%get_proxy()
    wetrho_in_w3_proxy = wetrho_in_w3%get_proxy()
    theta_in_wth_proxy = theta_in_wth%get_proxy()
    exner_in_w3_proxy = exner_in_w3%get_proxy()
    sd_orog_proxy = sd_orog%get_proxy()
    grad_xx_orog_proxy = grad_xx_orog%get_proxy()
    grad_xy_orog_proxy = grad_xy_orog%get_proxy()
    grad_yy_orog_proxy = grad_yy_orog%get_proxy()
    mr_v_proxy = mr_v%get_proxy()
    mr_cl_proxy = mr_cl%get_proxy()
    mr_ci_proxy = mr_ci%get_proxy()
    height_w3_proxy = height_w3%get_proxy()
    height_wth_proxy = height_wth%get_proxy()

    taux_blk_proxy = taux_orog_blk%get_proxy()
    tauy_blk_proxy = tauy_orog_blk%get_proxy()
    taux_gwd_proxy = taux_orog_gwd%get_proxy()
    tauy_gwd_proxy = tauy_orog_gwd%get_proxy()

    ! Initialise number of layers
    nlayers = du_blk_proxy%vspace%get_nlayers()

    ! Create a mesh object
    mesh => du_orog_blk%get_mesh()

    ! Look-up dofmaps for each function space
    map_w3 => du_blk_proxy%vspace%get_whole_dofmap()
    map_wtheta => dtemp_blk_proxy%vspace%get_whole_dofmap()
    map_adspc1_sd_orog => sd_orog_proxy%vspace%get_whole_dofmap()

    ! Initialise number of DoFs for w3
    ndf_w3 = du_blk_proxy%vspace%get_ndf()
    undf_w3 = du_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for wtheta
    ndf_wtheta = dtemp_blk_proxy%vspace%get_ndf()
    undf_wtheta = dtemp_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for adspc1_sd_orog
    ndf_adspc1_sd_orog = sd_orog_proxy%vspace%get_ndf()
    undf_adspc1_sd_orog = sd_orog_proxy%vspace%get_undf()

    ! Temporary variable for loop bound helps OpenMP
    loop_upper_bound = mesh%get_last_edge_cell()    
  
    applicable_points = 0   
    !cell index padded with 0s beyond applicable_points
    allocate(cell_index(loop_upper_bound))
    cell_index = 0
    do cell=1, loop_upper_bound
      ! Only call orographic_drag_kernel_code at points where the
      ! standard deviation of the subgrid orography is more than zero.
      if ( sd_orog_proxy%data(map_adspc1_sd_orog(1, cell)) > 0.0_r_def ) then

         ! We want to select the cells which have orographic drag
         applicable_points = applicable_points+1
         cell_index(applicable_points) = cell

      end if ! sd_orog_proxy%data(map_adspc1_sd_orog(1, cell)) > 0.0_r_def

    end do

    ! Check that the value of gw_segment is acceptable as a target length, otherwise use 1 for safety
    if (gw_segment .le. 0 .or. gw_segment .gt. applicable_points) then
      seg_target = 1
    else
      seg_target = gw_segment
    end if

    n_segments = (applicable_points + seg_target - 1) / seg_target

    ! Call orographic_drag_kernel_code if you have applicable points
    ! seg_len will be seg_target unless you have too few points at the end - l and u bound are of the segment iterated over
    if (applicable_points .gt. 0) then
      !$omp parallel do default(shared), private(segment, seg_len, l_bound, u_bound), schedule(dynamic)
      do segment=1, n_segments
        seg_len = min(seg_target, applicable_points - (segment - 1) * seg_target)
        l_bound = (segment - 1) * seg_target + 1
        u_bound = l_bound + seg_len - 1
        call orographic_drag_kernel_code(                                                &
                        loop_upper_bound, nlayers, du_blk_proxy%data, dv_blk_proxy%data, &
                        du_orog_gwd_proxy%data, dv_orog_gwd_proxy%data,                  &
                        dtemp_blk_proxy%data, dtemp_orog_gwd_proxy%data,                 &
                        u_in_w3_proxy%data, v_in_w3_proxy%data,                          &
                        wetrho_in_w3_proxy%data, theta_in_wth_proxy%data,                &
                        exner_in_w3_proxy%data,                                          &
                        sd_orog_proxy%data, grad_xx_orog_proxy%data,                     &
                        grad_xy_orog_proxy%data, grad_yy_orog_proxy%data,                &
                        mr_v_proxy%data, mr_cl_proxy%data, mr_ci_proxy%data,             &
                        height_w3_proxy%data, height_wth_proxy%data,                     &
                        taux_blk_proxy%data, tauy_blk_proxy%data,                        &
                        taux_gwd_proxy%data, tauy_gwd_proxy%data,                        &
                        ndf_w3, undf_w3, map_w3,                                         &
                        ndf_wtheta, undf_wtheta, map_wtheta,                             &
                        ndf_adspc1_sd_orog, undf_adspc1_sd_orog,                         &
                        map_adspc1_sd_orog, seg_len, cell_index(l_bound:u_bound))
      end do
      !$omp end parallel do
    end if
    deallocate(cell_index)

    ! Set halos dirty/clean for fields modified in the above loop
    call du_blk_proxy%set_dirty()
    call dv_blk_proxy%set_dirty()
    call du_orog_gwd_proxy%set_dirty()
    call dv_orog_gwd_proxy%set_dirty()
    call dtemp_blk_proxy%set_dirty()
    call dtemp_orog_gwd_proxy%set_dirty()
    call taux_blk_proxy%set_dirty()
    call tauy_blk_proxy%set_dirty()
    call taux_gwd_proxy%set_dirty()
    call tauy_gwd_proxy%set_dirty()

  end subroutine invoke_orographic_drag_kernel
  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the stochastic physics
  !> forcing pattern. At the moment it requires to pass the arrays
  !> my_coeff_rad and my_phi_stph via the kernel argument.
  !> PSyclone does not recognize arrays in the argument yet
  !> this functionality is being developed in PSyclone ticket 1312
  !> at https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
  subroutine invoke_spectral_2_cs_kernel_type(fp, longitude, pnm_star,         &
                                              coeffc_phase, coeffs_phase,      &
                                              stph_level_bottom,               &
                                              stph_level_top,                  &
                                              stph_n_min, stph_n_max)

  use spectral_2_cs_kernel_mod, ONLY: spectral_2_cs_code
  use mesh_mod, ONLY: mesh_type

  implicit none

  integer(KIND=i_def), intent(in) :: stph_level_bottom, stph_level_top,stph_n_min, stph_n_max
  type(field_type), intent(in) :: fp, longitude, pnm_star
  integer(KIND=i_def) cell
  integer(KIND=i_def) nlayers
  type(field_proxy_type) fp_proxy, longitude_proxy, pnm_star_proxy
  integer(KIND=i_def), pointer :: map_adspc1_longitude(:,:) => null(), map_adspc2_pnm_star(:,:) => null(), &
  &map_wspace(:,:) => null()
  integer(KIND=i_def) ndf_wspace, undf_wspace, ndf_adspc1_longitude, undf_adspc1_longitude, ndf_adspc2_pnm_star, &
  &undf_adspc2_pnm_star
  type(mesh_type), pointer :: mesh => null()

  ! Add arrays my_coeff_rad & my_phi_stph by hand
  real(kind=r_def), intent(in), dimension(:,:) :: coeffc_phase
  real(kind=r_def), intent(in), dimension(:,:) :: coeffs_phase
  integer(kind=i_def), parameter :: nranks_array = 2
  integer(kind=i_def), dimension(nranks_array) :: dims_array

  ! Get the upper bound for each rank of each scalar array
  ! Do dims_array for coeffc and coeffs?
  dims_array = shape(coeffc_phase)

  !
  ! Initialise field and/or operator proxies
  !
  fp_proxy     = fp%get_proxy()
  longitude_proxy  = longitude%get_proxy()
  pnm_star_proxy   = pnm_star%get_proxy()
  !
  ! Initialise number of layers
  !
  nlayers = fp_proxy%vspace%get_nlayers()
  !
  ! Create a mesh object
  !
  mesh => fp_proxy%vspace%get_mesh()
  !
  ! Look-up dofmaps for each function space
  !
  map_wspace           => fp_proxy%vspace%get_whole_dofmap()
  map_adspc1_longitude => longitude_proxy%vspace%get_whole_dofmap()
  map_adspc2_pnm_star  => pnm_star_proxy%vspace%get_whole_dofmap()
  !
  ! Initialise number of DoFs for wtheta
  !
  ndf_wspace  = fp_proxy%vspace%get_ndf()
  undf_wspace = fp_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc1_longitude
  !
  ndf_adspc1_longitude  = longitude_proxy%vspace%get_ndf()
  undf_adspc1_longitude = longitude_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc2_pnm_star
  !
  ndf_adspc2_pnm_star  = pnm_star_proxy%vspace%get_ndf()
  undf_adspc2_pnm_star = pnm_star_proxy%vspace%get_undf()
  !
  ! Call kernels and communication routines
  !
  do cell=1,mesh%get_last_edge_cell()
      call spectral_2_cs_code(nlayers, &
                              ! Add fields
                              fp_proxy%data, longitude_proxy%data, pnm_star_proxy%data, &
                              ! Add arrays
                              nranks_array, dims_array, coeffc_phase, coeffs_phase, &
                              ! Add SPT scalars
                              stph_level_bottom, stph_level_top, stph_n_min, stph_n_max, &
                              ! Add fields' assoc. space variables
                              ndf_wspace, undf_wspace, map_wspace(:,cell), &
                              ndf_adspc1_longitude, undf_adspc1_longitude, map_adspc1_longitude(:,cell), &
                              ndf_adspc2_pnm_star, undf_adspc2_pnm_star, map_adspc2_pnm_star(:,cell))
  end do
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  call fp_proxy%set_dirty()
  !
  !
  end subroutine invoke_spectral_2_cs_kernel_type
  !---------------------------------------------------------------------
  ! PSyclone currently doesn't work for DOMAIN kernels with stencil fields
  ! See PSyclone #1948

  ! Dumb first try, but does it work without this
  ! or do i need to modify the kernel / algorithm as well ?

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_heaviside_kernel_type(exner_wth, nplev, plevs, plev_heaviside, p_zero, kappa)
      USE heaviside_kernel_mod, ONLY: heaviside_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: exner_wth, plev_heaviside
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_heaviside_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_wth_data => null()
      TYPE(field_proxy_type) exner_wth_proxy, plev_heaviside_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_exner_wth(:,:) => null(), map_adspc2_plev_heaviside(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_exner_wth, undf_adspc1_exner_wth, ndf_adspc2_plev_heaviside, undf_adspc2_plev_heaviside
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      exner_wth_proxy = exner_wth%get_proxy()
      exner_wth_data => exner_wth_proxy%data
      plev_heaviside_proxy = plev_heaviside%get_proxy()
      plev_heaviside_data => plev_heaviside_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = exner_wth_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => exner_wth_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_exner_wth => exner_wth_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_heaviside => plev_heaviside_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_exner_wth
      !
      ndf_adspc1_exner_wth = exner_wth_proxy%vspace%get_ndf()
      undf_adspc1_exner_wth = exner_wth_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_heaviside
      !
      ndf_adspc2_plev_heaviside = plev_heaviside_proxy%vspace%get_ndf()
      undf_adspc2_plev_heaviside = plev_heaviside_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop0_start,loop0_stop
        !
        CALL heaviside_code(nlayers, exner_wth_data, nplev, plevs, plev_heaviside_data, p_zero, kappa, ndf_adspc1_exner_wth, &
&undf_adspc1_exner_wth, map_adspc1_exner_wth(:,cell), ndf_adspc2_plev_heaviside, undf_adspc2_plev_heaviside, &
&map_adspc2_plev_heaviside(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_heaviside_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_heaviside_kernel_type

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_temp_on_pres_kernel_type(temp, exner_wth, &
         height_wth, nplev, plevs, plev_temp, p_zero, kappa, ex_power)
      USE temp_on_pres_kernel_mod, ONLY: temp_on_pres_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa, ex_power
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: temp, exner_wth, height_wth, plev_temp
      INTEGER(KIND=i_def) cell
      INTEGER df
      INTEGER(KIND=i_def) loop1_start, loop1_stop
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_temp_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: height_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: temp_data => null()
      TYPE(field_proxy_type) temp_proxy, exner_wth_proxy, height_wth_proxy, plev_temp_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_temp(:,:) => null(), map_adspc2_plev_temp(:,:) => null(), map_wtheta(:,:) => null()
      INTEGER(KIND=i_def) ndf_aspc1_temp, undf_aspc1_temp, ndf_adspc1_temp, undf_adspc1_temp, ndf_wtheta, undf_wtheta, &
&ndf_adspc2_plev_temp, undf_adspc2_plev_temp
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      temp_proxy = temp%get_proxy()
      temp_data => temp_proxy%data
      exner_wth_proxy = exner_wth%get_proxy()
      exner_wth_data => exner_wth_proxy%data
      height_wth_proxy = height_wth%get_proxy()
      height_wth_data => height_wth_proxy%data
      plev_temp_proxy = plev_temp%get_proxy()
      plev_temp_data => plev_temp_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = temp_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => temp_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_temp => temp_proxy%vspace%get_whole_dofmap()
      map_wtheta => height_wth_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_temp => plev_temp_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for aspc1_temp
      !
      ndf_aspc1_temp = temp_proxy%vspace%get_ndf()
      undf_aspc1_temp = temp_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc1_temp
      !
      ndf_adspc1_temp = temp_proxy%vspace%get_ndf()
      undf_adspc1_temp = temp_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for wtheta
      !
      ndf_wtheta = height_wth_proxy%vspace%get_ndf()
      undf_wtheta = height_wth_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_temp
      !
      ndf_adspc2_plev_temp = plev_temp_proxy%vspace%get_ndf()
      undf_adspc2_plev_temp = plev_temp_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop1_start = 1
      loop1_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop1_start,loop1_stop
        !
        CALL temp_on_pres_code(nlayers, temp_data, exner_wth_data, height_wth_data, nplev, plevs, plev_temp_data, p_zero, kappa, &
&ex_power, ndf_adspc1_temp, undf_adspc1_temp, map_adspc1_temp(:,cell), ndf_wtheta, undf_wtheta, map_wtheta(:,cell), &
&ndf_adspc2_plev_temp, undf_adspc2_plev_temp, map_adspc2_plev_temp(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_temp_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_temp_on_pres_kernel_type

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_pres_interp_kernel_type(u_in_w3, exner_w3, nplev, plevs, plev_u, p_zero, kappa)
      USE pres_interp_kernel_mod, ONLY: pres_interp_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: u_in_w3, exner_w3, plev_u
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_u_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_w3_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: u_in_w3_data => null()
      TYPE(field_proxy_type) u_in_w3_proxy, exner_w3_proxy, plev_u_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_u_in_w3(:,:) => null(), map_adspc2_plev_u(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_u_in_w3, undf_adspc1_u_in_w3, ndf_adspc2_plev_u, undf_adspc2_plev_u
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      u_in_w3_proxy = u_in_w3%get_proxy()
      u_in_w3_data => u_in_w3_proxy%data
      exner_w3_proxy = exner_w3%get_proxy()
      exner_w3_data => exner_w3_proxy%data
      plev_u_proxy = plev_u%get_proxy()
      plev_u_data => plev_u_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = u_in_w3_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => u_in_w3_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_u_in_w3 => u_in_w3_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_u => plev_u_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_u_in_w3
      !
      ndf_adspc1_u_in_w3 = u_in_w3_proxy%vspace%get_ndf()
      undf_adspc1_u_in_w3 = u_in_w3_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_u
      !
      ndf_adspc2_plev_u = plev_u_proxy%vspace%get_ndf()
      undf_adspc2_plev_u = plev_u_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop0_start,loop0_stop
        !
        CALL pres_interp_code(nlayers, u_in_w3_data, exner_w3_data, nplev, plevs, plev_u_data, p_zero, kappa, ndf_adspc1_u_in_w3, &
&undf_adspc1_u_in_w3, map_adspc1_u_in_w3(:,cell), ndf_adspc2_plev_u, undf_adspc2_plev_u, map_adspc2_plev_u(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_u_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_pres_interp_kernel_type

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_geo_on_pres_kernel_type(height_w3, exner_w3, theta_wth, height_wth, exner_wth, nplev, plevs, plev_geopot, &
&p_zero, kappa, cp, gravity, ex_power)
      USE geo_on_pres_kernel_mod, ONLY: geo_on_pres_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa, cp, gravity, ex_power
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: height_w3, exner_w3, theta_wth, height_wth, exner_wth, plev_geopot
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_geopot_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: height_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: theta_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_w3_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: height_w3_data => null()
      TYPE(field_proxy_type) height_w3_proxy, exner_w3_proxy, theta_wth_proxy, height_wth_proxy, exner_wth_proxy, plev_geopot_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_height_w3(:,:) => null(), map_adspc2_plev_geopot(:,:) => null(), &
&map_wtheta(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_height_w3, undf_adspc1_height_w3, ndf_wtheta, undf_wtheta, ndf_adspc2_plev_geopot, &
&undf_adspc2_plev_geopot
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      height_w3_proxy = height_w3%get_proxy()
      height_w3_data => height_w3_proxy%data
      exner_w3_proxy = exner_w3%get_proxy()
      exner_w3_data => exner_w3_proxy%data
      theta_wth_proxy = theta_wth%get_proxy()
      theta_wth_data => theta_wth_proxy%data
      height_wth_proxy = height_wth%get_proxy()
      height_wth_data => height_wth_proxy%data
      exner_wth_proxy = exner_wth%get_proxy()
      exner_wth_data => exner_wth_proxy%data
      plev_geopot_proxy = plev_geopot%get_proxy()
      plev_geopot_data => plev_geopot_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = height_w3_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => height_w3_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_height_w3 => height_w3_proxy%vspace%get_whole_dofmap()
      map_wtheta => theta_wth_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_geopot => plev_geopot_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_height_w3
      !
      ndf_adspc1_height_w3 = height_w3_proxy%vspace%get_ndf()
      undf_adspc1_height_w3 = height_w3_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for wtheta
      !
      ndf_wtheta = theta_wth_proxy%vspace%get_ndf()
      undf_wtheta = theta_wth_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_geopot
      !
      ndf_adspc2_plev_geopot = plev_geopot_proxy%vspace%get_ndf()
      undf_adspc2_plev_geopot = plev_geopot_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop0_start,loop0_stop
        !
        CALL geo_on_pres_code(nlayers, height_w3_data, exner_w3_data, theta_wth_data, height_wth_data, exner_wth_data, nplev, &
&plevs, plev_geopot_data, p_zero, kappa, cp, gravity, ex_power, ndf_adspc1_height_w3, undf_adspc1_height_w3, &
&map_adspc1_height_w3(:,cell), ndf_wtheta, undf_wtheta, map_wtheta(:,cell), ndf_adspc2_plev_geopot, undf_adspc2_plev_geopot, &
&map_adspc2_plev_geopot(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_geopot_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_geo_on_pres_kernel_type
  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_thetaw_kernel_type(plev_thetaw, plev_temp, plev_qv, nplev, plevs)
      USE thetaw_kernel_mod, ONLY: thetaw_code
      USE mesh_mod, ONLY: mesh_type
      INTEGER(KIND=i_def), intent(in) :: nplev
      TYPE(field_type), intent(in) :: plev_temp, plev_qv, plev_thetaw
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers_plev_thetaw
      REAL(KIND=r_def), pointer, dimension(:) :: plev_thetaw_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: plev_qv_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: plev_temp_data => null()
      TYPE(field_proxy_type) plev_qv_proxy, plev_thetaw_proxy, plev_temp_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_plev_thetaw(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_plev_thetaw, undf_adspc1_plev_thetaw
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      plev_qv_proxy = plev_qv%get_proxy()
      plev_qv_data => plev_qv_proxy%data
      plev_temp_proxy = plev_temp%get_proxy()
      plev_temp_data => plev_temp_proxy%data
      plev_thetaw_proxy = plev_thetaw%get_proxy()
      plev_thetaw_data => plev_thetaw_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_plev_thetaw = plev_thetaw_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => plev_thetaw_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_plev_thetaw => plev_thetaw_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_plev_thetaw
      !
      ndf_adspc1_plev_thetaw = plev_thetaw_proxy%vspace%get_ndf()
      undf_adspc1_plev_thetaw = plev_thetaw_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell = loop0_start, loop0_stop, 1
        CALL thetaw_code(nlayers_plev_thetaw, plev_thetaw_data, plev_temp_data, plev_qv_data, nplev, plevs, ndf_adspc1_plev_thetaw, &
&undf_adspc1_plev_thetaw, map_adspc1_plev_thetaw(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_thetaw_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_thetaw_kernel_type
end module psykal_lite_phys_mod
