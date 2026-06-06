! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: mercurial -*-
!=============================================================================!
!                          XRD_SF MODULE                                      !
!=============================================================================!
! This module calculates the X-ray structure factors for the specified model. !
!-----------------------------------------------------------------------------!
! Written by Benjamin Shi, version 1.0, 15/09/2020                            !
!=============================================================================!

module xrd_sf
  use constants, only : dp, file_maxpath, periodic_table_name
  use trace, only: trace_entry, trace_exit
  use cell,  only : current_cell
  use model, only: model_state
  use density, only: electron_density
  use io, only: io_abort,io_allocate_abort, io_clock, io_finalise, io_initialise, io_initialise_farm, &
  & io_initialise_master, io_success, master_stdout, stdout, stderr, io_open_file, &
  & io_read_memory_usage, io_uppercase, io_code_real, io_code_present, io_code_integer, io_code_block, &
  & io_code_string, io_code_logical, io_getenv, io_atomic_to_unit, io_unit_label, io_flush
!Margo changed this (for calculating the atomic electron scattering factors)
  use otfg_scf
  use otfg_basis
 
  
  implicit none
  private

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: xrd_sf_calculate
contains
  subroutine xrd_sf_calculate(mdl)
    !=========================================================================!
    ! Main routine for calculate X-ray structure factors                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl (in) : The model to be used.                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Benjamin Shi, v1, 15/09/2020
    !=========================================================================!
    use io, only: io_atomic_to_unit, io_flush, io_unit_label, io_open_file,&
    io_code_real, io_code_present, io_code_integer, io_code_block, io_clock,&
    & io_code_string, io_code_logical, io_uppercase, io_allocate_abort, io_abort
    use constants, only : dp, cmplx_0, cmplx_1, cmplx_i, pi, two_pi
    use density, only: density_calculate_soft, density_augment, density_calculate_radial, &
    density_symmetrise, density_write_formatted, electron_density, density_zero,&
    density_allocate,density_deallocate
    use ion, only: ps_gmax, kkbeta, ps_num_points, rlog, rab, max_mesh, mesh, qfunc_bl, &
    max_bl_projectors, num_bl_projectors,bl_projector_l,bl_projector_m, rs_core_radial_charge, &
    rs_iam_radial_charge, ps_max_points,  atomic_radial_charge
    use basis, only : current_basis, basis_real_to_recip_grid, &
                      basis_radial_to_recip_grid,basis_radial_to_recip_grid
    use cell,  only : current_cell
    use model, only: model_state
    use parameters, only: current_params
    use comms!, only: on_root_node,num_in_gv_group,my_node_id
    implicit none

    ! Arguments
    type(model_state), intent(inout) :: mdl  ! The model to check for the grounds state
    type(ae_basis)     :: ab !Margo added this part for electron scattering factors

    real(kind=dp) :: time0, time1, time2 ! timing clocks
    integer :: i,j,k,l,i1,j1,k1,nsp,nsp1,ni,n,m,n1, m1, nx,ny,nz,ind,ierr,nsf, dummy_int, l_min, l_max, ng, nr, &
    sum_min, sum_max, m3                     ! Various dummy variables for loops
    real(kind=dp) :: dummy_core, dummy_dw, dummy_dw_exp, dummy_iam, dummy_phase
    complex(kind=dp) :: dummy_aug, dummy_sph, dummy_sph1, dummy_sph2, dummy1, dummy2, dummy_sf

    complex(kind=dp), dimension(:), allocatable :: sf_grid_dummy
    complex(kind=dp), dimension(:), allocatable :: total_sf, total_sf_soft_val, total_sf_aug_val, total_sf_core, &
    total_sf_iam, total_sf_hard_val ! Total contributions to SF
    complex(kind=dp), dimension(:), allocatable :: loc_total_sf, loc_total_sf_soft_val, loc_total_sf_aug_val, loc_total_sf_core, &
    loc_total_sf_iam, loc_total_sf_hard_val ! Total contributions to SF
    complex(kind=dp), dimension(:), allocatable :: total_sf_req, total_sf_soft_val_req, total_sf_aug_val_req, total_sf_core_req, &
    total_sf_iam_req, total_sf_hard_val_req ! Total contributions to SF
    real(kind=dp), dimension(3) :: H_vector, pos ! The (hkl) vector and ion position in a.u.


    ! Variables for the core and IAM charge density
    complex(kind=dp), dimension(:,:), allocatable ::  sf_grid_core, sf_grid_iam
    ! Reciprocal core and IAM charge density for each atom species on fine grid
    real(kind=dp), dimension(:,:), allocatable :: recip_core_radial_charge, recip_iam_radial_charge, recip_val_radial_charge 
    !Margo changed this
    real(kind=dp), dimension(:,:), allocatable :: recip_core_potential, &
     recip_val_potential, recip_ae_potential, core_potential, &
     val_potential, ae_potential, nuc_potential,&
     rad_grid
    real(kind=dp) :: ehar ! hartree energy
    integer                     :: ZZ, indx
    integer                     :: ZZZ
    character(len=8)            :: cname
    integer, dimension(:), allocatable         :: npts_arr(:)
    ! Fourier-Bessel transform of real space radial core and IAM charge density


    ! Variables for the soft and hard charge. By hard charge, I mean soft charge + pseudised augmentation charge
    type(electron_density)::dummy_charge      ! Soft charge on a FINE grid
    complex(kind=dp), dimension(:,:,:), allocatable ::  sf_grid    ! Reciprocal soft charge density for each atom and ion species
    real(kind=dp), dimension(:,:,:), allocatable :: hirshfeld_atomic_weighting ! Partitions the density into different atomic contributions
    complex(kind=dp), dimension(:,:), allocatable ::  sf_grid_soft, sf_grid_hard ! Reciprocal charge density for each atom species

    ! Variables for all-electron augmentation charge
    real(kind=dp), dimension(:,:,:,:), allocatable :: qfuncL_bl
    real(kind=dp), dimension(:,:,:,:,:), allocatable :: qfuncL_bl_recip
    complex(kind=dp), dimension(:,:,:,:,:), allocatable :: qfuncL_bl_fine_recip
    integer, dimension(:,:), allocatable :: blochl_l, blochl_l_index
    integer, dimension(current_cell%num_species) :: num_blochl_l
    real(kind=dp), dimension(0:3,-3:3,0:3,-3:3,0:6,-6:6) :: aug_coefficient
    complex(kind=dp), dimension(:,:,:,:), allocatable :: rho, rho_factor !rho is D_{u1,u2} and rho_factor is the factor that multiplies it
    logical                                      :: n_neq_m
    real(kind=dp) :: S, x, y, z, theta, phi
    ! S is the norm of the H_vector, theta and phi are the azimuthal and polar angles arising from cartesian coordinates x, y and z
    complex(kind=dp), dimension(:,:), allocatable ::  sf_grid_aug_val

    integer :: node, num_data, iostat
    logical :: in_first_gv_group

    ! Variables from reading the input file
    character(len=:),allocatable :: sf_string, dw_string ! The string of structure factor values taken from devel_code in seed.param
    integer, dimension(:), allocatable :: sf_list ! sf_string split into an integer array
    integer, dimension(:,:), allocatable :: sf_list_formatted, loc_sf_list_formatted, sf_list_formatted_req ! sf_list formatted correctly
    integer :: num_sf, num_sf_numbers, len_sf_string, len_dw_string, num_dw ! The number of structure factor values
    real(kind=dp), dimension(3,3,current_cell%num_species) :: dw_matrix ! 3x3 matrix containing the DW factors
    real(kind=dp), dimension(:,:), allocatable :: dw_list   ! Unformatted dw_list taken straight from .param file
    character(len=8), dimension(:), allocatable:: dw_list_species_symbol ! The species symbol taken from the .param file
    logical :: dw_aniso, dw_factor_set, sf_calc_all ! Logical determining whether a 3x3 matrix of debye waller factors is given or just the isotropic value
    ! Variables for writing the output files
    character(len=file_maxpath) :: filename, fmt_name ! On the root node, the output filename
    character(len=20), dimension(3) :: flags ! On the root node, the flags to pass to io_open_file
    integer :: output_file, fmt_file         ! On the root node, the handle for the output file
    integer :: status                        ! The status flag returned by io_open_file
    character(len=256)                                :: iomsg
    ! Reads in the list of structure factor values listed between the CALCULATE_XRD_SF block in devel_code in the .param file
    len_sf_string = len(trim(io_code_block(current_params%devel_code,'CALCULATE_XRD_SF')))
    allocate(character(len=len_sf_string) :: sf_string,stat=ierr)
    if (ierr/=0) call io_allocate_abort('sf_string','xrd_sf_calculate')
    sf_string = trim(io_code_block(current_params%devel_code,'CALCULATE_XRD_SF'))
    !
    ! Oh-so-clever way of counting blanks provokes bug in gfortran compiler debug build. Replace with
    ! easier to understand version which is also robust against input variations eg sequences of blanks. KR
    ! 
    !           num_sf_numbers = count([(sf_string(i:i),i=1,len_sf_string)].eq.' ')+1
    num_sf_numbers=0
    do i=1,len_sf_string
       if( i == 1 ) then
          if( sf_string(i:i) /= ' ' ) num_sf_numbers=num_sf_numbers + 1
       else
          if( sf_string(i:i) /= ' ' .and. sf_string(i-1:i-1) == ' ' ) num_sf_numbers=num_sf_numbers + 1
       end if
    end do
    
    time0=io_clock()

    ! We need the following until the end of the routine. So allocate them upfront
    allocate(sf_grid_core(current_basis%max_fine_grid_points, current_cell%num_species), stat=ierr)
    if (ierr/=0) call io_allocate_abort('sf_grid_core','xrd_sf_calculate')
    allocate(sf_grid_iam(current_basis%max_fine_grid_points, current_cell%num_species), stat=ierr)
    if (ierr/=0) call io_allocate_abort('sf_grid_iam','xrd_sf_calculate')



    ! Determines what structure factors to write out. If the block is empty, then the entire fine grid will be written out
    if (num_sf_numbers .eq. 1) then
      sf_calc_all = .true.
      num_sf = current_basis%num_fine_recip_points
      allocate(sf_list_formatted(current_basis%max_fine_grid_points,3), stat=ierr)
      if (ierr/=0) call io_allocate_abort('sf_list_formatted','xrd_sf_calculate')
      sf_list_formatted = 0
      ind = 0
      do nsf=1,current_basis%num_fine_recip_columns
        ny = current_basis%fine_recip_column_y(nsf)
        nz = current_basis%fine_recip_column_z(nsf)
        do nx = 1, current_basis%ngx_fine
          ind = ind + 1
          sf_list_formatted(ind,1) = current_basis%fine_recip_x_grid_coord(nx)
          sf_list_formatted(ind,2) = current_basis%fine_recip_y_grid_coord(ny)
          sf_list_formatted(ind,3) = current_basis%fine_recip_z_grid_coord(nz)
        end do
      end do
    else
      sf_calc_all = .false.
      ! Ensure that the values stated in parameter file makes sense
      if (modulo(num_sf_numbers,3) .ne. 0) call io_abort('Error: (hkl) vectors not stated correctly')
      num_sf = num_sf_numbers/3
      allocate(sf_list(num_sf_numbers),stat=ierr)
      if (ierr/=0) call io_allocate_abort('sf_list','xrd_sf_calculate')
      allocate(sf_list_formatted(num_sf,3), stat=ierr)
      if (ierr/=0) call io_allocate_abort('sf_list_formatted','xrd_sf_calculate')
      sf_list = 0
      sf_list_formatted = 0
      read(sf_string,*) sf_list
      ind = 0
      do nsf = 1, num_sf
        sf_list_formatted(nsf,1) = sf_list(ind+1)
        sf_list_formatted(nsf,2) = sf_list(ind+2)
        sf_list_formatted(nsf,3) = sf_list(ind+3)
        ind = ind + 3
      end do
    end if

    ! Reads in the list of Debye-Waller factors from devel_code block, DW factors can either be isotropic or anistropic
    dw_matrix = 0.0_dp
    if (io_code_present(current_params%devel_code,'DW_FACTOR')) then
      len_dw_string = len(trim(io_code_block(current_params%devel_code,'DW_FACTOR')))
      allocate(character(len=len_dw_string) :: dw_string)
      dw_string = trim(io_code_block(current_params%devel_code,'DW_FACTOR'))
      !
      ! Oh-so-clever way of counting blanks provokes bug in gfortran compiler debug build. Replace with
      ! easier to understand version which is also robust against input variations eg sequences of blanks. KR
      ! 
      !                 num_dw = count([(dw_string(i:i),i=1,len_dw_string)].eq.' ')+1
      num_dw=0
      do i=1,len_dw_string
         if( i == 1 ) then
            if( dw_string(i:i) /= ' ' ) num_dw=num_dw + 1
         else
            if( dw_string(i:i) /= ' ' .and. dw_string(i-1:i-1) == ' ' ) num_dw=num_dw + 1
         end if
      end do

      if (num_dw .eq. current_cell%num_species*2) then
        dw_aniso = .false.
      else if (num_dw .eq. current_cell%num_species*10) then
        dw_aniso = .true.
      else
        call io_abort('Error: Wrong number of Debye-Waller factors stated in .param file')
      end if

      if (dw_aniso) then
        allocate(dw_list(9,current_cell%num_species),stat=ierr)
      else
        allocate(dw_list(1,current_cell%num_species),stat=ierr)
      end if
      if (ierr/=0) call io_allocate_abort('dw_list','xrd_sf_calculate')
      allocate(dw_list_species_symbol(current_cell%num_species),stat=ierr)
      if (ierr/=0) call io_allocate_abort('dw_list_species_symbol','xrd_sf_calculate')
      dw_list = 0.0_dp
      if (dw_aniso) then
        read(dw_string,*) (dw_list_species_symbol(nsp),(dw_list(i,nsp),i=1,9),nsp=1,current_cell%num_species)
      else
        read(dw_string,*) (dw_list_species_symbol(nsp),dw_list(1,nsp),nsp=1,current_cell%num_species)
      end if

      do nsp = 1, current_cell%num_species
        dw_factor_set = .false.
        do nsp1 = 1, current_cell%num_species
          if (trim(io_uppercase(current_cell%species_symbol(nsp))) .eq. &
          trim(io_uppercase(dw_list_species_symbol(nsp1)))) then
          dw_factor_set = .true.
          if (dw_aniso) then
            ind = 0
            do i = 1, 3
              do j = 1, 3
                ind = ind + 1
                dw_matrix(i,j,nsp) = dw_list(ind, nsp1)
              end do
            end do
          else
            dw_matrix(1,1,nsp) = real(dw_list(1,nsp1),dp)
            dw_matrix(2,2,nsp) = real(dw_list(1,nsp1),dp)
            dw_matrix(3,3,nsp) = real(dw_list(1,nsp1),dp)
          end if
        else
          cycle
        end if
      end do
      if (.not. dw_factor_set) call io_abort('Error: Debye-Waller factor for '// trim(current_cell%species_symbol(nsp)) //&
      &' not stated in .param file')
    end do
  end if

  allocate(recip_core_radial_charge(ps_max_points,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('recip_core_radial_charge','xrd_sf_calculate')
  allocate(recip_iam_radial_charge(ps_max_points,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('recip_iam_radial_charge','xrd_sf_calculate')
  allocate(recip_val_radial_charge(ps_max_points,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('recip_val_radial_charge','xrd_sf_calculate')
  !again, Margo changed this
  allocate(core_potential(max_mesh,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('core_potential','xrd_sf_calculate')
  allocate(val_potential(max_mesh,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('val_potential','xrd_sf_calculate')
  allocate(ae_potential(max_mesh,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('ae_potential','xrd_sf_calculate') 
  allocate(nuc_potential(max_mesh,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('nuc_potential','xrd_sf_calculate')   
  allocate(rad_grid(max_mesh,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('rad_grid','xrd_sf_calculate')
  
  allocate(recip_core_potential(ps_max_points,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('recip_core_potential','xrd_sf_calculate')
  allocate(recip_val_potential(ps_max_points,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('recip_val_potential','xrd_sf_calculate')
  allocate(recip_ae_potential(ps_max_points,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('recip_ae_potential','xrd_sf_calculate')
    
  ! Get the core electron density and IAM charge on the fine grid
  ! Zero the arrays
  sf_grid_core = cmplx_0
  sf_grid_iam = cmplx_0
  recip_core_radial_charge = 0.0_dp
  core_potential = 0.0_dp
  val_potential = 0.0_dp
  ae_potential = 0.0_dp
  nuc_potential = 0.0_dp
  
  recip_core_potential = 0.0_dp
  recip_val_potential = 0.0_dp
  recip_ae_potential = 0.0_dp
  
  recip_iam_radial_charge = 0.0_dp
  recip_val_radial_charge = 0.0_dp
  !Again, Morgo's changes, reading in the potentials maybe easier
  open(unit=400,file='one.proj_pot',status='old',action='read')
  read(400,*) nsp
  allocate(npts_arr(nsp))
  do i=1,nsp
     read(400,*) npts_arr(i)
  end do
  close(400)
  open(unit=401,file='two.proj_pot',status='old',action='read')
 
 
  do i=1,nsp
     read(401,*) ZZZ
     do j=1,npts_arr(i)
        read(401,*) rad_grid(j,i)
     end do
     do j=1,npts_arr(i)
        read(401,*) core_potential(j,i)
     end do
     do j=1,npts_arr(i)
        read(401,*) ae_potential(j,i)
     end do
     do j=1,npts_arr(i)
        read(401,*) nuc_potential(j,i)
     end do
  end do
  close(401)
  ! Generate the radial reciprocal space grids
  do nsp = 1, current_cell%num_species
  !Margo changed this for electron scattering factors
    ZZ = -1 ! Set atomic number to -1 to check later
    do i=1,size(periodic_table_name)
       ! Get the atomic_number for the symbol
       indx  = index(current_cell%species_symbol(nsp),':')
       if(indx>0) then
          cname = current_cell%species_symbol(nsp)(:indx-1)
       else
          cname = current_cell%species_symbol(nsp)
       end if
       if(cname.eq.periodic_table_name(i)) ZZ=i
    end do
    if(ZZ<1) call io_abort('Error: otfg_inquire_psp: current_cell%species_symbol not known')

    !call otfg_init_ae_basis(ab, ZZ) !Margo - comment to myself: this part works fine (ie calculating the potential) but maybe it is easier just to read in the potentials from the two. file written in the otfg_generate step 

    !call otfg_ae_hartree(ab,rs_core_radial_charge(:,nsp),core_potential(:,nsp),ehar)
    !call otfg_ae_hartree(ab,(rs_iam_radial_charge(:,nsp)-rs_core_radial_charge(:,nsp)),val_potential(:,nsp),ehar)
    !call otfg_ae_hartree(ab,rs_iam_radial_charge(:,nsp),ae_potential(:,nsp),ehar)
    
    !End of the electron scattering factor changes
    call xrd_sf_radial_transform(0,0,0,mesh(nsp),rlog(:,nsp),rab(:,nsp),ps_num_points(nsp), &
    ps_gmax(nsp),rs_core_radial_charge(:,nsp),recip_core_radial_charge(:,nsp))
    call xrd_sf_radial_transform(0,0,0,mesh(nsp),rlog(:,nsp),rab(:,nsp),ps_num_points(nsp), &
    ps_gmax(nsp),rs_iam_radial_charge(:,nsp),recip_iam_radial_charge(:,nsp))
    !Margo changed this, supposed to give the valence scattering factor
    call xrd_sf_radial_transform(0,0,0,mesh(nsp),rlog(:,nsp),rab(:,nsp),ps_num_points(nsp), &
    ps_gmax(nsp),(rs_iam_radial_charge(:,nsp)-rs_core_radial_charge(:,nsp)),recip_val_radial_charge(:,nsp))
    
    !call xrd_sf_radial_transform(0,0,0,mesh(nsp),rlog(:,nsp),rab(:,nsp),ps_num_points(nsp), &
    !ps_gmax(nsp),core_potential(:,nsp),recip_core_potential(:,nsp))
    !call xrd_sf_radial_transform(0,0,0,mesh(nsp),rlog(:,nsp),rab(:,nsp),ps_num_points(nsp), &
    !ps_gmax(nsp),val_potential(:,nsp),recip_val_potential(:,nsp))
    call xrd_sf_radial_transform(0,0,0,mesh(nsp),rad_grid(:,nsp),rab(:,nsp),ps_num_points(nsp), &
    ps_gmax(nsp),(ae_potential(:,nsp)+nuc_potential(:,nsp)),recip_ae_potential(:,nsp))

    open(unit=2005, file="additional_data.out", status="unknown")
    write(2005,*) mesh(nsp)
    write(2005, *) rab(:,nsp)
    write(2005, *) ps_num_points(nsp)
    write(2005, *) ps_gmax(nsp)
  end do
  time1=io_clock()


  !Margo changed this 
  open(unit=2001,file="recip_iam_radial_charge.out",status="unknown")
  write(2001,*) recip_iam_radial_charge
  close(2001)

  open(unit=2002,file="radial_grid.out",status="unknown")
  write(2002,*) ps_num_points, ps_gmax
  close(2002)
  
  open(unit=2003,file="recip_core_radial_charge.out",status="unknown")
  write(2003,*) recip_core_radial_charge
  close(2003)

  open(unit=2004,file="recip_val_radial_charge.out",status="unknown")
  write(2004,*) recip_val_radial_charge
  close(2004)
     
  !open(unit=2005,file="recip_core_potential.out",status="unknown")
  !write(2005,*) recip_core_potential
  !close(2005)
  
  !open(unit=2006,file="recip_val_potential.out",status="unknown")
  !write(2006,*) recip_val_potential
  !close(2006)
  
  open(unit=2007,file="recip_ae_potential.out",status="unknown")
  write(2007,*) recip_ae_potential
  close(2007)
  
  
  open(unit=2008,file="debug_ae_potential.out",status="unknown")
  write(2008,*) ae_potential
  close(2008)
  
  open(unit=2009,file="debug_nuc_potential.out",status="unknown")
  write(2009,*) nuc_potential
  close(2009)
  !END OF MARGO'S CHANGES
  
!!  if(on_root_node)write(stdout,*) 'set up:', time1-time0

  allocate(sf_grid_dummy(current_basis%max_fine_grid_points), stat=ierr)
  if (ierr/=0) call io_allocate_abort('sf_grid_dummy','xrd_sf_calculate')

  ! Core structure factors
  dummy_sf = cmplx_0
  ! Loop over species
  do nsp=1,mdl%cell%num_species
    sf_grid_dummy = cmplx_0
    ! Recip space radial charge for this species -> sf_grid_dummy
    call basis_radial_to_recip_grid(ps_num_points(nsp),ps_gmax(nsp), &
    recip_core_radial_charge(:,nsp),'EVEN',sf_grid_dummy(:),'FINE')
    ! Loop over atoms of species number nsp
    do ni=1,mdl%cell%num_ions_in_species(nsp)
      ! Get Cartesian position of atom nsp/ni
      pos = mdl%cell%ionic_positions(1,ni,nsp)*mdl%cell%real_lattice(1,:) + &
      mdl%cell%ionic_positions(2,ni,nsp)*mdl%cell%real_lattice(2,:) + &
      mdl%cell%ionic_positions(3,ni,nsp)*mdl%cell%real_lattice(3,:)
      do i=1,current_basis%num_fine_recip_points
        ! Calculate phase factor
        dummy_sf = exp(1.0_dp*cmplx(0.0_dp,current_basis%g_fine_recip_grid(1,i)*pos(1) + &
        current_basis%g_fine_recip_grid(2,i)*pos(2) + &
        current_basis%g_fine_recip_grid(3,i)*pos(3),dp))
        sf_grid_core(i,nsp)= sf_grid_core(i,nsp) + dummy_sf*sf_grid_dummy(i)
      end do
    end do
  end do

  ! IAM structure factors
  dummy_sf = cmplx_0

  ! Loop over species
  do nsp=1,mdl%cell%num_species
    sf_grid_dummy = cmplx_0
    ! Recip space radial charge for this species -> spec_charge FINE grid
    call basis_radial_to_recip_grid(ps_num_points(nsp),ps_gmax(nsp), &
    recip_iam_radial_charge(:,nsp),'EVEN',sf_grid_dummy,'FINE')
    ! Loop over atoms of species number nsp
    do ni=1,mdl%cell%num_ions_in_species(nsp)
      ! Get Cartesian position of atom nsp/ni
      pos = mdl%cell%ionic_positions(1,ni,nsp)*mdl%cell%real_lattice(1,:) + &
      mdl%cell%ionic_positions(2,ni,nsp)*mdl%cell%real_lattice(2,:) + &
      mdl%cell%ionic_positions(3,ni,nsp)*mdl%cell%real_lattice(3,:)
      do i=1,current_basis%num_fine_recip_points
        ! Calculate structure factor exp(-ig.pos)
        dummy_sf = exp(1.0_dp*cmplx(0.0_dp,current_basis%g_fine_recip_grid(1,i)*pos(1) + &
        current_basis%g_fine_recip_grid(2,i)*pos(2) + &
        current_basis%g_fine_recip_grid(3,i)*pos(3),dp))
        sf_grid_iam(i,nsp)= sf_grid_iam(i,nsp) + dummy_sf*sf_grid_dummy(i)
      end do
    end do
  end do

  deallocate(sf_grid_dummy,stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating sf_grid_dummy in xrd_sf_calculate')
  deallocate(recip_core_radial_charge, stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating recip_core_radial_charge in xrd_sf_calculate')
  deallocate(recip_iam_radial_charge, stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating recip_iam_radial_charge in xrd_sf_calculate')

  time2=io_clock()
!!  if(on_root_node)write(stdout,*) 'core IAM:', time2-time1

  ! Allocate arrays
  allocate(hirshfeld_atomic_weighting(current_basis%max_fine_grid_points,current_cell%max_ions_in_species,&
  current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('hirshfeld_atomic_weighting','xrd_sf_calculate')

  allocate(sf_grid(current_basis%max_fine_grid_points,current_cell%max_ions_in_species,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('sf_grid','xrd_sf_calculate')
  allocate(sf_grid_soft(current_basis%max_fine_grid_points,current_cell%num_species), stat=ierr)
  if (ierr/=0) call io_allocate_abort('sf_grid_hard','xrd_sf_calculate')
  allocate(sf_grid_hard(current_basis%max_fine_grid_points,current_cell%num_species), stat=ierr)

  ! Zero the arrays
  sf_grid = cmplx_0
  sf_grid_soft = cmplx_0
  sf_grid_hard = cmplx_0

  call density_allocate(dummy_charge,den_type='D',do_ked=.false.)
  call density_zero(dummy_charge)
  call density_calculate_soft(mdl%wvfn,mdl%occ,dummy_charge, do_ked=.false.)
  call density_symmetrise(dummy_charge)
  call xrd_sf_hirshfeld_weight(mdl,hirshfeld_atomic_weighting)
  call xrd_sf_density_to_sf(dummy_charge,sf_grid,hirshfeld_atomic_weighting,'C')  ! soft valence charge in reciprocal space

  do nsp = 1, current_cell%num_species
    do ni = 1, current_cell%num_ions_in_species(nsp)
      sf_grid_soft(:,nsp) = sf_grid_soft(:,nsp) + sf_grid(:,ni,nsp)
    end do
  end do

  ! Calculate the structure factor from hard valence change = soft + pseudised augmentation charge
  sf_grid = cmplx_0
  call density_zero(dummy_charge)
  call density_calculate_soft(mdl%wvfn,mdl%occ,dummy_charge, do_ked=.false.)
  call density_symmetrise(dummy_charge)
  call density_augment(mdl%wvfn,mdl%occ,dummy_charge,do_ked=.false.)
  call xrd_sf_density_to_sf(dummy_charge,sf_grid,hirshfeld_atomic_weighting,'C')  ! hard valence charge in reciprocal space

  do nsp = 1, current_cell%num_species
    do ni = 1, current_cell%num_ions_in_species(nsp)
      sf_grid_hard(:,nsp) = sf_grid_hard(:,nsp) + sf_grid(:,ni,nsp)
    end do
  end do

  deallocate(sf_grid, stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating sf_grid in xrd_sf_calculate')



  time1=io_clock()
!!  if(on_root_node)write(stdout,*) 'soft:   ', time1-time2

  ! Calculate the all-electron augmentation charge contribution to structure factor
  ! Generate indexing for blochl projectors for l angular momentum only,
  ! no m angular momentum because it doesn't affect radial components
  num_blochl_l = 0

  do nsp = 1, current_cell%num_species
    dummy_int = 0
    do while (dummy_int .lt. num_bl_projectors(nsp))
      dummy_int = dummy_int + 2*bl_projector_l(dummy_int+1,nsp) + 1
      num_blochl_l(nsp) = num_blochl_l(nsp) + 1
    end do
  end do

  allocate(blochl_l(maxval(num_blochl_l),current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('hirshfeld_atomic_weighting','xrd_sf_calculate')
  allocate(blochl_l_index(max_bl_projectors,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('hirshfeld_atomic_weighting','xrd_sf_calculate')
  allocate(qfuncL_bl_recip(maxval(ps_num_points),max_bl_projectors,max_bl_projectors,0:6,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('hirshfeld_atomic_weighting','xrd_sf_calculate')
  allocate(qfuncL_bl(maxval(kkbeta),max_bl_projectors,max_bl_projectors,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('hirshfeld_atomic_weighting','xrd_sf_calculate')
  allocate(qfuncL_bl_fine_recip(current_basis%max_fine_grid_points,max_bl_projectors,max_bl_projectors,&
  0:6,current_cell%num_species),stat=ierr)

  blochl_l = 0
  blochl_l_index = 0
  qfuncL_bl = 0.0_dp
  qfuncL_bl_recip = 0.0_dp
  qfuncL_bl_fine_recip = cmplx_0


  do nsp = 1, current_cell%num_species
    m = 0
    do n=1, num_blochl_l(nsp)
      m = m + 2*bl_projector_l(m+1,nsp) + 1 ! 2l+1
      blochl_l(n,nsp) = bl_projector_l(m,nsp)
      blochl_l_index(m-2*bl_projector_l(m,nsp):m,nsp) = n
    end do
  end do

  ! Set up a small qfuncL_bl
  do nsp = 1, current_cell%num_species
    n1 = 0
    do n = 1, num_blochl_l(nsp)
      n1 = n1 + 2*blochl_l(n,nsp) + 1
      m1 = 0
      do m = 1, num_blochl_l(nsp)
        m1 = m1 + 2*blochl_l(m,nsp) + 1
        qfuncL_bl(1:kkbeta(nsp),n,m,nsp) = qfunc_bl(1:kkbeta(nsp),n1,m1,nsp)
      end do
    end do
  end do

  ! Generate the deltaR (which has variable name qfunc_bl) on the fine reciprocal grid
  do nsp = 1, current_cell%num_species
    do n = 1, num_blochl_l(nsp)
      do m = n, num_blochl_l(nsp)
        l_min = abs(blochl_l(n,nsp) - blochl_l(m,nsp))
        l_max = blochl_l(n,nsp) + blochl_l(m,nsp)
        do l = l_min, l_max
          call xrd_sf_radial_transform(l,-1,0,kkbeta(nsp),&
          rlog(:,nsp),rab(:,nsp),ps_num_points(nsp),ps_gmax(nsp), qfuncL_bl(:,n,m,nsp), &
          qfuncL_bl_recip(:,n,m,l,nsp))

          call basis_radial_to_recip_grid(ps_num_points(nsp),ps_gmax(nsp), &
          qfuncL_bl_recip(:,n,m,l,nsp),'S',qfuncL_bl_fine_recip(:,n,m,l,nsp),'FINE')

          if (n .ne. m) then
            qfuncL_bl_recip(:,m,n,l,nsp) = qfuncL_bl_recip(:,n,m,l,nsp)
            qfuncL_bl_fine_recip(:,m,n,l,nsp) = qfuncL_bl_fine_recip(:,n,m,l,nsp)
          end if
        end do
      end do
    end do
  end do
  qfuncL_bl_fine_recip = sqrt(4.0_dp*pi)*qfuncL_bl_fine_recip ! Normalise

  deallocate(qfuncL_bl,stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating qfuncL_bl in xrd_sf_calculate')
  deallocate(qfuncL_bl_recip,stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating qfuncL_bl_recip in xrd_sf_calculate')
  allocate(rho(current_cell%max_ions_in_species,max_bl_projectors,max_bl_projectors,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('rho','xrd_sf_calculate')
  allocate(rho_factor(current_basis%max_fine_grid_points,&
  max_bl_projectors,max_bl_projectors,current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('rho_factor','xrd_sf_calculate')

  time1=io_clock()
!!  if(on_root_node)write(stdout,*) 'blochl   ', time1-time2

  rho = cmplx_0
  rho_factor = cmplx_0
  ! Calculate the rho matrix and gaunt_coefficients
  call xrd_sf_rho(mdl%wvfn,mdl%wvfn%nkpts,mdl%wvfn%nspins,mdl%occ,current_cell%kpoint_weights_on_node,rho)
  aug_coefficient = 0.0_dp
  call xrd_sf_gaunt_coefficients(aug_coefficient)

  do nsf = 1, current_basis%num_fine_recip_points
    S = 0.0_dp
    H_vector = 0.0_dp
    H_vector = current_basis%g_fine_recip_grid(1:3,nsf)
    S=sqrt(sum(H_vector(:)**2))
    x = H_vector(1)
    y = H_vector(2)
    z = H_vector(3)
    if( S > 0.0_dp ) then
       theta = atan2(sqrt(real(x,dp)**2+real(y,dp)**2),real(z,dp))
    else
       theta = 0.0_dp
    end if
    if( abs(x) > 0.0_dp .or. abs(y) .gt. 0.0_dp ) then
       phi = atan2(real(y,dp),real(x,dp))
    else
       phi=0.0_dp
    end if
    do nsp = 1, current_cell%num_species
      do n = 1,num_bl_projectors(nsp)
        do m = n,num_bl_projectors(nsp)
          n_neq_m = (n .ne. m)
          sum_min = 0
          sum_max = 0
          sum_min = abs(bl_projector_l(n,nsp) - bl_projector_l(m,nsp))
          sum_max = bl_projector_l(n,nsp) + bl_projector_l(m,nsp)
          m3 = -bl_projector_m(n,nsp) + bl_projector_m(m,nsp)
          dummy1 = cmplx_0
          dummy2 = cmplx_0
          do l = sum_min,sum_max
            if (abs(m3) .gt. l) then
              cycle
            else if (aug_coefficient(bl_projector_l(n,nsp), bl_projector_m(n,nsp), &
              bl_projector_l(m,nsp),bl_projector_m(m,nsp),l,m3) .eq. 0.0_dp) then
              cycle
            else
              dummy1 = dummy1 + 4.0_dp*pi*(cmplx_i**l)*aug_coefficient(bl_projector_l(n,nsp), bl_projector_m(n,nsp), &
              bl_projector_l(m,nsp),bl_projector_m(m,nsp),l,m3)*xrd_sf_apply_ylm_cmplx(theta,phi,l,m3)*&
              qfuncL_bl_fine_recip(nsf,blochl_l_index(n,nsp),blochl_l_index(m,nsp),l,nsp)
              if (n_neq_m) then
                dummy2 = dummy2 + 4.0_dp*pi*(cmplx_i**l)*aug_coefficient(bl_projector_l(m,nsp), bl_projector_m(m,nsp), &
                bl_projector_l(n,nsp),bl_projector_m(n,nsp),l,-m3)*xrd_sf_apply_ylm_cmplx(theta,phi,l,-m3)*&
                qfuncL_bl_fine_recip(nsf,blochl_l_index(n,nsp),blochl_l_index(m,nsp),l,nsp)
              end if
            end if
          end do
          rho_factor(nsf,n,m,nsp) = dummy1

          if (n_neq_m) then
            rho_factor(nsf,m,n,nsp) = dummy2
          end if
        end do
      end do
    end do
  end do

  allocate(sf_grid_aug_val(current_basis%max_fine_grid_points, current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('sf_grid_aug_val','xrd_sf_calculate')

  sf_grid_aug_val = cmplx_0
  do nsf = 1, current_basis%num_fine_recip_points
    H_vector = 0.0_dp
    H_vector = current_basis%g_fine_recip_grid(1:3,nsf)
    dummy_sph = cmplx_0
    do nsp = 1, current_cell%num_species
      do ni = 1, current_cell%num_ions_in_species(nsp)
        pos = 0.0_dp
        pos(:) = current_cell%ionic_positions(1,ni,nsp)*current_cell%real_lattice(1,:)   &
        & + current_cell%ionic_positions(2,ni,nsp)*current_cell%real_lattice(2,:) &
        & + current_cell%ionic_positions(3,ni,nsp)*current_cell%real_lattice(3,:)
        dummy_aug = cmplx_0
        do n = 1,num_bl_projectors(nsp)
          do m = n,num_bl_projectors(nsp)
            dummy_aug = dummy_aug + rho_factor(nsf,n,m,nsp)*rho(ni,n,m,nsp)
            if (n .ne. m) then
              dummy_aug = dummy_aug + rho_factor(nsf,m,n,nsp)*conjg(rho(ni,m,n,nsp))
            end if
          end do
        end do
        sf_grid_aug_val(nsf,nsp) = sf_grid_aug_val(nsf,nsp) + dummy_aug*exp(1.0_dp*cmplx(0.0_dp,H_vector(1)*pos(1) &
        + H_vector(2)*pos(2) + H_vector(3)*pos(3),dp))
      end do
    end do
  end do

  if(num_in_gv_group == 1)then
    ! this case everything we need is already on the root node
    ! and we just need to do the final calculation and write out what we need

    allocate(total_sf(num_sf),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rho','xrd_sf_calculate')
    allocate(total_sf_iam(num_sf),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rho','xrd_sf_calculate')
    allocate(total_sf_core(num_sf),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rho','xrd_sf_calculate')
    allocate(total_sf_soft_val(num_sf),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rho','xrd_sf_calculate')
    allocate(total_sf_hard_val(num_sf),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rho','xrd_sf_calculate')
    allocate(total_sf_aug_val(num_sf),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rho','xrd_sf_calculate')

    total_sf_iam = cmplx_0
    total_sf_core = cmplx_0
    total_sf_soft_val = cmplx_0
    total_sf_hard_val = cmplx_0
    total_sf_aug_val = cmplx_0
    total_sf = cmplx_0

    do nsf = 1, num_sf
      i = sf_list_formatted(nsf, 1)
      j = sf_list_formatted(nsf, 2)
      k = sf_list_formatted(nsf, 3)
      i1 = i
      j1 = j
      k1 = k
      if (i .lt. 0) i1 = current_basis%ngx_fine + i
      if (j .lt. 0) j1 = current_basis%ngy_fine + j
      if (k .lt. 0) k1 = current_basis%ngz_fine + k
      ind = i1 + j1*(current_basis%ngx_fine) + k1*(current_basis%ngx_fine*current_basis%ngy_fine) + 1
      H_vector = 0.0_dp
      H_vector = i*current_cell%recip_lattice(1,:) + j*current_cell%recip_lattice(2,:) + k*current_cell%recip_lattice(3,:)
      H_vector = H_vector*1.8897261258369282_dp/two_pi ! Convert to 1/Angstrom
      do nsp = 1, current_cell%num_species
        dummy_dw = 0.0_dp
        dummy_dw_exp = 0.0_dp
        do n = 1, 3
          do m = 1, 3
            dummy_dw = dummy_dw + H_vector(n)*H_vector(m)*dw_matrix(n,m,nsp)
          end do
        end do
        dummy_dw_exp = exp(-0.25_dp*dummy_dw)
        total_sf_iam(nsf) = total_sf_iam(nsf) + sf_grid_iam(ind,nsp)*dummy_dw_exp
        total_sf_core(nsf) = total_sf_core(nsf) + sf_grid_core(ind,nsp)*dummy_dw_exp
        total_sf_aug_val(nsf) = total_sf_aug_val(nsf) + sf_grid_aug_val(ind,nsp)*dummy_dw_exp
        total_sf_soft_val(nsf) = total_sf_soft_val(nsf) + sf_grid_soft(ind,nsp)*dummy_dw_exp
        total_sf_hard_val(nsf) = total_sf_hard_val(nsf) + sf_grid_hard(ind,nsp)*dummy_dw_exp
      end do
      ! Symmetrise the total_sf phase using the phase from total_sf_hard_augment, which has been symmetrised
      dummy_sf = cmplx_0
      dummy_sf = total_sf_hard_val(nsf) + total_sf_core(nsf)
      if (abs(dummy_sf) .lt. 0.00001) then
        dummy_phase = 0.0_dp
      else if (abs(real(dummy_sf)) .lt. 0.00001 .and. aimag(dummy_sf) .gt. 0.00001) then
        dummy_phase = pi/2.0_dp
      else if (abs(real(dummy_sf)) .lt. 0.00001 .and. aimag(dummy_sf) .lt. -0.00001) then
        dummy_phase = -pi/2.0_dp
      else
        dummy_phase = atan2(aimag(dummy_sf),real(dummy_sf,dp))
      end if
      total_sf(nsf) = abs(total_sf_soft_val(nsf) + total_sf_aug_val(nsf) + total_sf_core(nsf))*exp(cmplx(0.0_dp,dummy_phase,dp))
      total_sf_aug_val(nsf) = total_sf(nsf) - total_sf_soft_val(nsf) - total_sf_core(nsf)
    end do

    filename = trim(current_params%seedname)//'.xrd_sf'
    ! Open the file
    flags(1) = 'write'
    flags(2) = 'formatted'
    flags(3) = 'replace'
    ! The root node must do the writing
    if (on_root_node) then
      call io_open_file(filename,3,flags,output_file,status,iomsg=iomsg)
      if (status.ne.io_success) call io_abort('Error write_local_potential: Failed to open .xrd_sf file: '//trim(iomsg))
      write(output_file,1) 'test1', 'k', 'l', 'Re(F_PAW)', 'Im(F_PAW)', 'Re(F_IAM)', 'Im(F_IAM)', 'Re(F_PP)', 'Im(F_PP)' , &
      &'Re(F_core)', 'Im(F_core)', 'Re(F_soft)', 'Im(F_soft)', 'Re(F_aug)', 'Im(F_aug)'
      do nsf = 1, num_sf
        i = sf_list_formatted(nsf, 1)
        j = sf_list_formatted(nsf, 2)
        k = sf_list_formatted(nsf, 3)
        write(output_file,2) i, j, k, real(total_sf(nsf)), aimag(total_sf(nsf)), real(total_sf_iam(nsf)),&
        aimag(total_sf_iam(nsf)), real(total_sf_hard_val(nsf) + total_sf_core(nsf)), &
        aimag(total_sf_hard_val(nsf) + total_sf_core(nsf)), real(total_sf_core(nsf)), &
        aimag(total_sf_core(nsf)), real(total_sf_soft_val(nsf)), aimag(total_sf_soft_val(nsf)), &
        real(total_sf_aug_val(nsf)), aimag(total_sf_aug_val(nsf))
      end do
    end if

    deallocate(sf_grid_core, stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating sf_grid_core in xrd_sf_calculate')
    deallocate(sf_grid_iam, stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating sf_grid_iam in xrd_sf_calculate')
    deallocate(sf_grid_soft, stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating sf_grid_soft in xrd_sf_calculate')
    deallocate(sf_grid_hard, stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating sf_grid_hard in xrd_sf_calculate')
    deallocate(sf_grid_aug_val, stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating sf_grid_aug_val in xrd_sf_calculate')


  else
    ! G-vector parallel case

    ! now the data is spread across the nodes.
    ! to save bandwidth and RAM we'll compute all SF on the local nodes
    ! We'll then send all the SF to the root node in turn
    ! If we need all then we'll just write to disk
    ! If we want a subset we'll look for these and write out
    ! This part of the code is written by JRY Aug 2021

    if(.not.sf_calc_all) then
      allocate(total_sf_req(num_sf),stat=ierr)
      if (ierr/=0) call io_allocate_abort('total_sf_req','xrd_sf_calculate')
      allocate(total_sf_iam_req(num_sf),stat=ierr)
      if (ierr/=0) call io_allocate_abort('total_sf_iam_req','xrd_sf_calculate')
      allocate(total_sf_core_req(num_sf),stat=ierr)
      if (ierr/=0) call io_allocate_abort('total_sf_core_req','xrd_sf_calculate')
      allocate(total_sf_soft_val_req(num_sf),stat=ierr)
      if (ierr/=0) call io_allocate_abort('total_sf_soft_val_req','xrd_sf_calculate')
      allocate(total_sf_hard_val_req(num_sf),stat=ierr)
      if (ierr/=0) call io_allocate_abort('total_sf_hard_val_req','xrd_sf_calculate')
      allocate(total_sf_aug_val_req(num_sf),stat=ierr)
      if (ierr/=0) call io_allocate_abort('total_sf_aug_val_req','xrd_sf_calculate')
      total_sf_iam_req = cmplx_0
      total_sf_core_req = cmplx_0
      total_sf_soft_val_req = cmplx_0
      total_sf_hard_val_req = cmplx_0
      total_sf_aug_val_req = cmplx_0
      total_sf_req = cmplx_0
      allocate(sf_list_formatted_req(num_sf,3), stat=ierr)
      if (ierr/=0) call io_allocate_abort('sf_list_formatted_req','xrd_sf_calculate')
      sf_list_formatted_req=sf_list_formatted
      deallocate(sf_list_formatted)
      allocate(sf_list_formatted(current_basis%max_fine_grid_points,3), stat=ierr)
      if (ierr/=0) call io_allocate_abort('sf_list_formatted','xrd_sf_calculate')
      sf_list_formatted=0
      ind = 0
      do nsf=1,current_basis%num_fine_recip_columns
        ny = current_basis%fine_recip_column_y(nsf)
        nz = current_basis%fine_recip_column_z(nsf)
        do nx = 1, current_basis%ngx_fine
          ind = ind + 1
          sf_list_formatted(ind,1) = current_basis%fine_recip_x_grid_coord(nx)
          sf_list_formatted(ind,2) = current_basis%fine_recip_y_grid_coord(ny)
          sf_list_formatted(ind,3) = current_basis%fine_recip_z_grid_coord(nz)
        end do
      end do
    end if


    allocate(total_sf(current_basis%max_fine_grid_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('total_sf','xrd_sf_calculate')
    allocate(total_sf_iam(current_basis%max_fine_grid_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('total_sf_iam','xrd_sf_calculate')
    allocate(total_sf_core(current_basis%max_fine_grid_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('total_sf_core','xrd_sf_calculate')
    allocate(total_sf_soft_val(current_basis%max_fine_grid_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('total_sf_soft_val','xrd_sf_calculate')
    allocate(total_sf_hard_val(current_basis%max_fine_grid_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('total_sf_hard_val','xrd_sf_calculate')
    allocate(total_sf_aug_val(current_basis%max_fine_grid_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('total_sf_aug_val','xrd_sf_calculate')

    total_sf_iam = cmplx_0
    total_sf_core = cmplx_0
    total_sf_soft_val = cmplx_0
    total_sf_hard_val = cmplx_0
    total_sf_aug_val = cmplx_0
    total_sf = cmplx_0

    do nsf = 1,current_basis%num_fine_recip_points
      H_vector = current_basis%g_fine_recip_grid(1:3,nsf)*1.8897261258369282_dp/two_pi ! Convert to 1/Angstrom
      do nsp = 1, current_cell%num_species
        dummy_dw = 0.0_dp
        dummy_dw_exp = 0.0_dp
        do n = 1, 3
          do m = 1, 3
            dummy_dw = dummy_dw + H_vector(n)*H_vector(m)*dw_matrix(n,m,nsp)
          end do
        end do
        dummy_dw_exp = exp(-0.25_dp*dummy_dw)
        total_sf_iam(nsf) = total_sf_iam(nsf) + sf_grid_iam(nsf,nsp)*dummy_dw_exp
        total_sf_core(nsf) = total_sf_core(nsf) + sf_grid_core(nsf,nsp)*dummy_dw_exp
        total_sf_aug_val(nsf) = total_sf_aug_val(nsf) + sf_grid_aug_val(nsf,nsp)*dummy_dw_exp
        total_sf_soft_val(nsf) = total_sf_soft_val(nsf) + sf_grid_soft(nsf,nsp)*dummy_dw_exp
        total_sf_hard_val(nsf) = total_sf_hard_val(nsf) + sf_grid_hard(nsf,nsp)*dummy_dw_exp
      end do

      ! Symmetrise the total_sf phase using the phase from total_sf_hard_augment, which has been symmetrised
      dummy_sf = cmplx_0
      dummy_sf = total_sf_hard_val(nsf) + total_sf_core(nsf)
      if (abs(dummy_sf) .lt. 0.00001_dp) then
        dummy_phase = 0.0_dp
      else if (abs(real(dummy_sf)) .lt. 0.00001_dp .and. aimag(dummy_sf) .gt. 0.00001_dp) then
        dummy_phase = pi/2.0_dp
      else if (abs(real(dummy_sf)) .lt. 0.00001_dp .and. aimag(dummy_sf) .lt. -0.00001_dp) then
        dummy_phase = -pi/2.0_dp
      else
        dummy_phase = atan2(aimag(dummy_sf),real(dummy_sf,dp))
      end if
      total_sf(nsf) = abs(total_sf_soft_val(nsf) + total_sf_aug_val(nsf) + total_sf_core(nsf))*exp(cmplx(0.0_dp,dummy_phase,dp))
      total_sf_aug_val(nsf) = total_sf(nsf) - total_sf_soft_val(nsf) - total_sf_core(nsf)
    end do



    ! Now all the nodes have their structure factors

    ! Now work out if this node is in the first g-vector group
    in_first_gv_group = .false.
    if(id_in_kp_group.eq.master_in_kp_group)in_first_gv_group=.true.





    if(.not.on_root_node)then

      ! If not on root node, but in first g-vector group send data to root node
      if(in_first_gv_group)then
        call comms_send(root_node_id,current_basis%num_fine_recip_points,1)
        call comms_send(root_node_id,sf_list_formatted(1,1),current_basis%max_fine_grid_points*3)
        call comms_send(root_node_id,total_sf(1),current_basis%max_fine_grid_points)
        call comms_send(root_node_id,total_sf_hard_val(1),current_basis%max_fine_grid_points)
        call comms_send(root_node_id,total_sf_iam(1),current_basis%max_fine_grid_points)
        call comms_send(root_node_id,total_sf_core(1),current_basis%max_fine_grid_points)
        call comms_send(root_node_id,total_sf_soft_val(1),current_basis%max_fine_grid_points)
        call comms_send(root_node_id,total_sf_aug_val(1),current_basis%max_fine_grid_points)
      endif


    else

      ! Allocate work space

      allocate(loc_total_sf(current_basis%max_fine_grid_points),stat=ierr)
      if (ierr/=0) call io_allocate_abort('loc_total_sf','xrd_sf_calculate')
      allocate(loc_total_sf_iam(current_basis%max_fine_grid_points),stat=ierr)
      if (ierr/=0) call io_allocate_abort('loc_total_sf_iam','xrd_sf_calculate')
      allocate(loc_total_sf_core(current_basis%max_fine_grid_points),stat=ierr)
      if (ierr/=0) call io_allocate_abort('loc_total_sf_core','xrd_sf_calculate')
      allocate(loc_total_sf_soft_val(current_basis%max_fine_grid_points),stat=ierr)
      if (ierr/=0) call io_allocate_abort('loc_total_sf_soft_val','xrd_sf_calculate')
      allocate(loc_total_sf_hard_val(current_basis%max_fine_grid_points),stat=ierr)
      if (ierr/=0) call io_allocate_abort('loc_total_sf_hard_val','xrd_sf_calculate')
      allocate(loc_total_sf_aug_val(current_basis%max_fine_grid_points),stat=ierr)
      if (ierr/=0) call io_allocate_abort('loc_total_sf_aug_val','xrd_sf_calculate')
      allocate(loc_sf_list_formatted(current_basis%max_fine_grid_points,3), stat=ierr)
      if (ierr/=0) call io_allocate_abort('loc_sf_list_formatted','xrd_sf_calculate')

      loc_sf_list_formatted=0
      loc_total_sf_iam = cmplx_0
      loc_total_sf_core = cmplx_0
      loc_total_sf_soft_val = cmplx_0
      loc_total_sf_hard_val = cmplx_0
      loc_total_sf_aug_val = cmplx_0
      loc_total_sf = cmplx_0

      filename = trim(current_params%seedname)//'.xrd_sf'
      ! Open the file
      flags(1) = 'write'
      flags(2) = 'formatted'
      flags(3) = 'replace'
      ! The root node must do the writing
      call io_open_file(filename,3,flags,output_file,status,iomsg=iomsg)
      if (status.ne.io_success) call io_abort('Error write_local_potential: Failed to open .xrd_sf file: '//trim(iomsg))
      write(output_file,1) 'test2', 'k', 'l', 'Re(F_PAW)', 'Im(F_PAW)', 'Re(F_IAM)', 'Im(F_IAM)', 'Re(F_PP)', 'Im(F_PP)' , &
      &'Re(F_core)', 'Im(F_core)', 'Re(F_soft)', 'Im(F_soft)', 'Re(F_aug)', 'Im(F_aug)'


      if(sf_calc_all) then
        ! We are on the root node, so first of all write out local current
        do nsf = 1,current_basis%num_fine_recip_points
          !  write(myunit,'(4i4,2f13.6)') nsf,sf_list_formatted(nsf,1),sf_list_formatted(nsf,2),sf_list_formatted(nsf,3),total_sf(nsf)
          write(output_file,2,iostat=iostat) sf_list_formatted(nsf,1),sf_list_formatted(nsf,2),sf_list_formatted(nsf,3), &
          real(total_sf(nsf),dp), aimag(total_sf(nsf)), real(total_sf_iam(nsf),dp),&
          aimag(total_sf_iam(nsf)), real(total_sf_hard_val(nsf) + total_sf_core(nsf),dp), &
          aimag(total_sf_hard_val(nsf) + total_sf_core(nsf)), real(total_sf_core(nsf),dp), &
          aimag(total_sf_core(nsf)), real(total_sf_soft_val(nsf)), aimag(total_sf_soft_val(nsf)), &
          real(total_sf_aug_val(nsf),dp), aimag(total_sf_aug_val(nsf))
        end do
        if(iostat.ne.0)call io_abort('ERROR: xrd_sf_calculate: error writing .xrd_sf file')
      else
        ! look for any matches
        do nsf = 1,current_basis%num_fine_recip_points
          do i=1,num_sf
            if (sf_list_formatted(nsf,1)==sf_list_formatted_req(i,1).and.sf_list_formatted(nsf,2)==sf_list_formatted_req(i,2)&
            .and.sf_list_formatted(nsf,3)==sf_list_formatted_req(i,3)) then
            total_sf_iam_req(i) = total_sf_iam(nsf)
            total_sf_core_req(i) = total_sf_core(nsf)
            total_sf_soft_val_req(i) = total_sf_soft_val(nsf)
            total_sf_hard_val_req(i) = total_sf_hard_val(nsf)
            total_sf_aug_val_req(i) = total_sf_aug_val(nsf)
            total_sf_req(i) = total_sf(nsf)
            cycle
          end if
        end do
      end do

    end if



    ! Now loop over the rest of the nodes indices and receive their current in node order
    do node=1,num_in_gv_group

      if(nodes_in_gv_group(node).ne.root_node_id)then
        call comms_recv(nodes_in_gv_group(node),num_data,1)
        call comms_recv(nodes_in_gv_group(node),loc_sf_list_formatted(1,1),current_basis%max_fine_grid_points*3)
        call comms_recv(nodes_in_gv_group(node),loc_total_sf(1),current_basis%max_fine_grid_points)
        call comms_recv(nodes_in_gv_group(node),loc_total_sf_hard_val(1),current_basis%max_fine_grid_points)
        call comms_recv(nodes_in_gv_group(node),loc_total_sf_iam(1),current_basis%max_fine_grid_points)
        call comms_recv(nodes_in_gv_group(node),loc_total_sf_core(1),current_basis%max_fine_grid_points)
        call comms_recv(nodes_in_gv_group(node),loc_total_sf_soft_val(1),current_basis%max_fine_grid_points)
        call comms_recv(nodes_in_gv_group(node),loc_total_sf_aug_val(1),current_basis%max_fine_grid_points)

        if(sf_calc_all) then
          do nsf=1,num_data
            write(output_file,2,iostat=iostat) loc_sf_list_formatted(nsf,1),loc_sf_list_formatted(nsf,2),&
            loc_sf_list_formatted(nsf,3), &
            real(loc_total_sf(nsf),dp), aimag(loc_total_sf(nsf)), real(loc_total_sf_iam(nsf),dp),&
            aimag(loc_total_sf_iam(nsf)), real(loc_total_sf_hard_val(nsf) + loc_total_sf_core(nsf),dp), &
            aimag(loc_total_sf_hard_val(nsf) + loc_total_sf_core(nsf)), real(loc_total_sf_core(nsf),dp), &
            aimag(loc_total_sf_core(nsf)), real(loc_total_sf_soft_val(nsf)), aimag(loc_total_sf_soft_val(nsf)), &
            real(loc_total_sf_aug_val(nsf),dp), aimag(loc_total_sf_aug_val(nsf))
          end do
          if(iostat.ne.0)call io_abort('ERROR: xrd_sf_calculate: error writing .xrd_sf file')

        else

          do nsf = 1,num_data
            do i=1,num_sf
              if (loc_sf_list_formatted(nsf,1)==sf_list_formatted_req(i,1).and.&
                loc_sf_list_formatted(nsf,2)==sf_list_formatted_req(i,2)&
                .and.loc_sf_list_formatted(nsf,3)==sf_list_formatted_req(i,3)) then
                  total_sf_iam_req(i) = loc_total_sf_iam(nsf)
                  total_sf_core_req(i) = loc_total_sf_core(nsf)
                  total_sf_soft_val_req(i) = loc_total_sf_soft_val(nsf)
                  total_sf_hard_val_req(i) = loc_total_sf_hard_val(nsf)
                  total_sf_aug_val_req(i) = loc_total_sf_aug_val(nsf)
                  total_sf_req(i) = loc_total_sf(nsf)
                  cycle
                end if
              end do
            end do
          end if


      endif ! if(node id is in 1st gv-group, but not root node)

    end do ! node loop

    if(.not.sf_calc_all) then
      do nsf=1,num_sf
        write(output_file,2,iostat=iostat) sf_list_formatted_req(nsf,1),sf_list_formatted_req(nsf,2),&
        sf_list_formatted_req(nsf,3), &
        real(total_sf_req(nsf),dp), aimag(total_sf_req(nsf)), real(total_sf_iam_req(nsf),dp),&
        aimag(total_sf_iam_req(nsf)), real(total_sf_hard_val_req(nsf) + total_sf_core_req(nsf),dp), &
        aimag(total_sf_hard_val_req(nsf) + total_sf_core_req(nsf)), real(total_sf_core_req(nsf),dp), &
        aimag(total_sf_core_req(nsf)), real(total_sf_soft_val_req(nsf)), aimag(total_sf_soft_val_req(nsf)), &
        real(total_sf_aug_val_req(nsf),dp), aimag(total_sf_aug_val_req(nsf))
      end do
      if(iostat.ne.0)call io_abort('ERROR: xrd_sf_calculate: error writing .xrd_sf file')
    end if

  endif ! on_root?

end if


  time2=io_clock()
  ! if(on_root_node)write(stdout,*) 'end:   ', time2-time1

  ! close the output file
  if (on_root_node) close(output_file)

  1 format (3A4, *(A12))
  2 format (3I4, *(F12.6))
  return
end subroutine xrd_sf_calculate

subroutine xrd_sf_density_to_sf(den,atomic_recip_grid,Hirshfeld_w,cmp)

  !--------------------------------------------------------------------------!
  ! Purpose:                                                                 !
  ! This subroutine transforms on component of a density into structure      !
  ! factor. The optional argument cmp determines which component of the      !
  ! density is to be transformed: 'C'-charge, 'S'-spin. If cmp is not present!
  ! then the default is cmp='C'.                                             !
  !                                                                          !
  ! In the case of NCM the magnitude of the spin is returned
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  ! den(in): The density to be transformed                                   !
  ! recip_grid(out): The transformed density on the FINE grid                !
  ! cmp(in, optional): The flag indicating which component is transformed    !
  !--------------------------------------------------------------------------!
  ! Parent module variables used: type(electron_density)                     !
  !--------------------------------------------------------------------------!
  ! Modules used: constants, io, basis, cell, trace                          !
  !--------------------------------------------------------------------------!
  ! Key internal variables: none                                             !
  !--------------------------------------------------------------------------!
  ! Necessary conditions: cmp is either 'C' or 'S'                           !
  !--------------------------------------------------------------------------!
  ! Modified from density_to_recip by Benjamin Shi, 15/09/2020               !
  !--------------------------------------------------------------------------!
  use constants, only : cmplx_0
  use io,        only : io_abort,io_allocate_abort
  use basis!,     only : basis_recip_grid_to_real,current_basis
  use cell,      only : current_cell
  use trace,     only : trace_entry, trace_exit
  use comms

  implicit none

  type(electron_density), intent(in):: den  ! Input density
  complex(kind=dp), dimension(current_basis%max_fine_grid_points,current_cell%max_ions_in_species,&
  current_cell%num_species), intent(out):: atomic_recip_grid ! Output transformed C or S density
  character, optional, intent(in)::cmp
  real(kind=dp), dimension(:,:,:), intent(in) :: Hirshfeld_w


  complex(kind=dp), dimension(current_basis%max_fine_grid_points) :: recip_grid ! Output transformed C or S density
  complex(kind=dp), dimension(current_basis%max_fine_grid_points) :: total_recip_grid ! Output transformed C or S density

  logical::transf_ch ! .true. if charge is transformed, .false. if spin is transformed

  integer                         :: ng, nsp, ni
  integer                         :: status

  call trace_entry('xrd_sf_density_to_sf',status)

  if(present(cmp))then
    if(.not.(cmp.eq.'C'.or.cmp.eq.'S'))call io_abort('ERROR: density_to_recip - cmp must be C or S')
  endif

  !  ! Check for optional cmp
  !  if(.not.present(cmp))then
  !     transf_ch=.true.
  !  else
  !     if(cmp.eq.'C')transf_ch=.true.
  !     if(cmp.eq.'S')transf_ch=.false.
  !  endif

  !  ! If spin density, make sure it exits
  !  if(.not.transf_ch.and.den%nspins.eq.1)then
  !     if(den%has_ncm) then
  !        allocate(magnitude_of_spin(1:current_basis%num_fine_real_points), stat=status)
  !        if (status/=0) call io_allocate_abort('magnitude_of_spin','density_to_recip')
  !        call density_spin_magnitude(den,magnitude_of_spin)
  !     else
  !        call io_abort('ERROR: density_to_recip - transform spin but den%nspins=1')
  !     end if
  !  end if

  ! transf_ch=.true.

  ! Copy either charge or spin density onto recip_grid
  ! if(transf_ch)then
  if(den%have_cmplx_den)then
    recip_grid(1:current_basis%num_fine_real_points)=den%charge(1:current_basis%num_fine_real_points)
  else
    recip_grid(1:current_basis%num_fine_real_points)=cmplx(den%real_charge(1:current_basis%num_fine_real_points),0.0_dp,dp)
  endif
  ! endif
  ! else
  !   if (den%has_ncm) then
  !     recip_grid(1:current_basis%num_fine_real_points)=magnitude_of_spin(1:current_basis%num_fine_real_points)
  !   else if(den%have_cmplx_den)then
  !     recip_grid(1:current_basis%num_fine_real_points)=den%spin(1:current_basis%num_fine_real_points)
  !   else
  !     recip_grid(1:current_basis%num_fine_real_points)=cmplx(den%real_spin(1:current_basis%num_fine_real_points),0.0_dp,dp)
  !   endif
  ! endif

  do ng=1+current_basis%num_fine_real_points,current_basis%max_fine_grid_points
    recip_grid(ng) = cmplx_0
  end do

  total_recip_grid = cmplx_0
  total_recip_grid = recip_grid

  do nsp = 1, current_cell%num_species
    do ni = 1, current_cell%num_ions_in_species(nsp)
      recip_grid = cmplx_0
      recip_grid(:) = total_recip_grid(:)*Hirshfeld_w(:,ni,nsp)
      call basis_real_to_recip_grid(recip_grid,'FINE')
      !          call basis_recip_grid_to_real(recip_grid,'FINE')
      atomic_recip_grid(:,ni,nsp) = conjg(recip_grid(:)) !bxs21 thinks this is needed due to diff defn of FT see email 17aug21
    end do
  end do

  ! if (den%has_ncm .and. (.not. transf_ch)) then
  !   deallocate(magnitude_of_spin, stat=status)
  !   if (status/=0) call io_abort('Error, deallocating magnitude_of_spin in xrd_sf_density_to_sf')
  ! end if

  call trace_exit('xrd_sf_density_to_sf',status)

  return

end subroutine xrd_sf_density_to_sf

subroutine xrd_sf_hirshfeld_weight(mdl, Hirshfeld_w)
  !=========================================================================!
  ! This subroutine calculates the Hirshfeld weight on the FINE grid        !
  !-------------------------------------------------------------------------!
  ! Arguments:                                                              !
  !   mdl (in) : The model to be used.                                      !
  !   Hirshfeld_w (out) : The Hirshfeld weights                             !
  !-------------------------------------------------------------------------!
  ! Modules used:                                                           !
  !-------------------------------------------------------------------------!
  ! Key Internal Variables:                                                 !
  !   spec_charge : charge for a species                                    !
  !-------------------------------------------------------------------------!
  ! Necessary conditions:                                                   !
  !   None                                                                  !
  !-------------------------------------------------------------------------!
  ! Written by Benjamin Shi, Version 1, 15/09/2020                          !
  !=========================================================================!

  use constants, only : pi
  use io, only : io_open_file, io_atomic_to_unit, io_unit_label
  use basis, only : current_basis,   &
  basis_radial_to_recip_grid,basis_recip_grid_to_real
  use density, only : electron_density,density_zero,density_allocate,density_deallocate, &
  density_symmetrise
  use ion, only : atomic_radial_charge, ps_num_points, ps_gmax
  use parameters, only: current_params
  use cell, only: cell_cart_lattice_to_abc
  use trace, only: trace_entry, trace_exit
  use comms

  implicit none

  ! Arguments

  type(model_state), intent(in) :: mdl ! The current model
  real(kind=dp), dimension(:,:,:), intent(out) :: Hirshfeld_w
  ! Local variables

  complex(kind=dp), dimension(:), allocatable :: spec_charge, atom_charge
  type(electron_density)::lcao_charge
  real(kind=dp), dimension(:,:,:), allocatable :: atomic_charges

  integer :: status ! The status flag returned by io_open_file
  integer :: i,j,nsp,ni
  real(kind=dp), dimension(3)::pos  ! Position of atom
  complex(kind=dp)::sf       ! Structure factor for one grid point

  ! allocate and calculate the density
  call trace_entry('xrd_sf_hirshfeld_weight',status)

  allocate(spec_charge(current_basis%max_fine_grid_points),stat=status)
  if (status/=0) call io_allocate_abort('spec_charge','xrd_sf_hirshfeld_weight')
  allocate(atom_charge(current_basis%max_fine_grid_points),stat=status)
  if (status/=0) call io_allocate_abort('atom_charge','xrd_sf_hirshfeld_weight')
  allocate(atomic_charges(current_basis%max_fine_grid_points,current_cell%max_ions_in_species,&
  current_cell%num_species),stat=status)
  if (status/=0) call io_allocate_abort('atomic_charges','xrd_sf_hirshfeld_weight')

  call density_allocate(lcao_charge,1,'R',do_ked=.false.)
  call density_zero(lcao_charge)
  atomic_charges = 0.0_dp
  Hirshfeld_w = 0.0_dp

  ! Loop over species
  do nsp=1,mdl%cell%num_species

    spec_charge=(0.0_dp,0.0_dp)
    ! Recip space radial charge for this species -> spec_charge FINE grid
    call basis_radial_to_recip_grid(ps_num_points(nsp),ps_gmax(nsp), &
    atomic_radial_charge(:,nsp),'EVEN',spec_charge,'FINE')

    ! Loop over atoms of species number nsp
    do ni=1,mdl%cell%num_ions_in_species(nsp)

      ! Get Cartesian position of atom nsp/ni
      pos = mdl%cell%ionic_positions(1,ni,nsp)*mdl%cell%real_lattice(1,:) + &
      mdl%cell%ionic_positions(2,ni,nsp)*mdl%cell%real_lattice(2,:) + &
      mdl%cell%ionic_positions(3,ni,nsp)*mdl%cell%real_lattice(3,:)
      atom_charge = cmplx(0.0_dp,0.0_dp,dp)

      do i=1,current_basis%num_fine_recip_points
        ! Calculate structure factor exp(-ig.pos)
        sf = exp(-1.0_dp*cmplx(0.0_dp,current_basis%g_fine_recip_grid(1,i)*pos(1) + &
        current_basis%g_fine_recip_grid(2,i)*pos(2) + &
        current_basis%g_fine_recip_grid(3,i)*pos(3),dp))
        atom_charge(i) = sf * spec_charge(i) * mdl%cell%mixture_weight(ni,nsp)
      end do
      call basis_recip_grid_to_real(atom_charge,'FINE')
      do j=1,current_basis%num_fine_real_points
        if ( real(atom_charge(j)) .lt. 0.0_dp ) then
          atom_charge(j) = cmplx(0.0_dp,0.0_dp,dp)
        end if
      end do
      atomic_charges(:,ni,nsp) = real(atom_charge(:),dp)
      lcao_charge%real_charge(:) = lcao_charge%real_charge(:) + real(atom_charge(:),dp)
    end do
  end do

  do nsp = 1, mdl%cell%num_species
    do ni = 1, mdl%cell%num_ions_in_species(nsp)
      Hirshfeld_w(:,ni,nsp) = atomic_charges(:,ni,nsp)/lcao_charge%real_charge(:)
    end do
  end do

  ! deallocate the array

  deallocate(atomic_charges,stat=status)
  if (status/=0) call io_abort('Error in deallocating atom_positions in xrd_sf_hirshfeld_weight')
  deallocate(atom_charge,stat=status)
  if (status/=0) call io_abort('Error in deallocating atom_charge in xrd_sf_hirshfeld_weight')
  deallocate(spec_charge,stat=status)
  if (status/=0) call io_abort('Error in deallocating spec_charge in xrd_sf_hirshfeld_weight')
  call density_deallocate(lcao_charge)

  call trace_exit('xrd_sf_hirshfeld_weight',status)

  return
end subroutine xrd_sf_hirshfeld_weight

subroutine xrd_sf_radial_transform(l,lmax,power,nrpts,rlog,rab,nqpts,qmax,rfunc,qfunc)

  use io, only : io_abort,io_allocate_abort

  !=========================================================================!
  ! Transforms a radial real space function on a logarithmic grid to a      !
  ! regular reciprocal space grid                                           !
  !-------------------------------------------------------------------------!
  ! Arguments:                                                              !
  !                                                                         !
  ! 1) l      : input  : which bessel transform required                    !
  ! 2) lmax   : input  : maximum l for which the lookup table is constructed!
  ! 3) power  : input  : power of r needed in integral                      !
  ! 4) nrpts  : input  : number of real space points                        !
  ! 5) rlog   : input  : the logarithmic grid                               !
  ! 6) rab    : input  : the "spacings" of the logarithmic grid             !
  ! 7) nqpts  : input  : number of points in the regular reciprocal grid    !
  ! 8) qmax   : input  : the maximum q-vector for the transform             !
  ! 9) rfunc  : input  : the real space radial function                     !
  ! 10) qfunc : output : the reciprocal space function                      !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Parent module variables used:                                           !
  !                                                                         !
  ! 0) None                                                                 !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Modules used:                                                           !
  !                                                                         !
  ! 1) io                                                                   !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Key Internal Variables:                                                 !
  !                                                                         !
  ! 1) check  : a saved variable that is used to check if the lookup table  !
  !              can still be used.                                         !
  ! 2) lookup : the lookup table, contains the weighted Bessel functions.   !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Necessary conditions:                                                   !
  !                                                                         !
  ! 1) nrpts should be odd                                                  !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Chris J. Pickard, V0.0, 8th June 2000                                   !
  ! This is essentially a clone of the routine from otfg                    !
  ! B.X. Shi reported that using the OTFG version sometimes gave NaNs       !
  ! This version is modified slightly. Eventually, this routine             !
  ! should be replaced with the more efficient multi version                !
  ! at such time we should test the OTFG version again                      !
  !=========================================================================!

  use trace, only: trace_entry, trace_exit

  implicit none

  integer,                          intent(in) :: l
  integer,                          intent(in) :: lmax
  integer,                          intent(in) :: power
  integer,                          intent(in) :: nrpts
  real(kind=dp), dimension(:),      intent(in) :: rlog
  real(kind=dp), dimension(:),      intent(in) :: rab
  integer,                          intent(in) :: nqpts
  real(kind=dp),                    intent(in) :: qmax
  real(kind=dp), dimension(:),      intent(in) :: rfunc
  real(kind=dp), dimension(:),      intent(out) :: qfunc

  !-------------------------------------------------------------------------!

  integer       :: nq,nr,ll
  integer       :: ierr
  real(kind=dp) :: q,qr
  real(kind=dp) :: rwgt(nrpts)
  real(kind=dp) :: current_check(7)

  real(kind=dp), save, dimension(:,:,:), allocatable :: lookup
  real(kind=dp) :: check(7) = 0.0_dp ! This will be saved

  integer                         :: status

  call trace_entry('xrd_sf_radial_transform',status)

  ! Check conditions
  if(nrpts==2*(nrpts/2)) call io_abort('Error: xrd_sf_radial_transform: nrpts should be odd')
  if(size(rlog)<nrpts) call io_abort('Error: xrd_sf_radial_transform: rlog should be of size nrpts')
  if(size(rab)<nrpts) call io_abort('Error: xrd_sf_radial_transform: rab should be of size nrpts')
  if(size(rfunc)<nrpts) call io_abort('Error: xrd_sf_radial_transform: rfunc should be of size nrpts')
  if(size(qfunc)<nqpts) call io_abort('Error: xrd_sf_radial_transform: qfunc should be of size nqpts')


  ! ** If lmax = -999, reset xrd_sf_radial_transform

  if(lmax==-999) then
    if (allocated(lookup)) then
      deallocate(lookup,stat=ierr)
      if (ierr/=0) call io_abort('Error in deallocating lookup in xrd_sf_radial_transform')
    end if
    check = 0.0_dp
    call trace_exit('xrd_sf_radial_transform',status)

    return
  end if

  ! ** If lmax == -1, only store the requested angular momentum channel

  if(lmax==-1) then

    current_check(1) = real(nqpts,dp)
    current_check(2) = qmax
    current_check(3) = real(nrpts,dp)
    current_check(4) = rlog(1)
    current_check(5) = rab(1)
    current_check(6) = real(-l,dp)
    current_check(7) = real(power,dp)

    if((.not.allocated(lookup)).or.(.not.all(check==current_check))) then
      if (allocated(lookup)) then
        deallocate(lookup,stat=ierr)
        if (ierr/=0) call io_abort('Error in deallocating lookup in xrd_sf_radial_transform')
      end if
      allocate(lookup(nrpts,nqpts,1),stat=ierr)
      if (ierr/=0) call io_allocate_abort('lookup','xrd_sf_radial_transform')

      ! ** Set the integration weights

      rwgt(1) = 4.0_dp/12.0_dp

      ! Optimisation: trap special cases of power being 0 or 1, which are extremely common and allow
      ! significant simplifications
      if(power==0) then
        do nr=2,nrpts-3,2
          rwgt(nr)=rab(nr)*16.0_dp/12.0_dp
          rwgt(nr+1)=rab(nr+1)*8.0_dp/12.0_dp
        end do

        rwgt(nrpts-1) = rab(nrpts-1)*16.0_dp/12.0_dp
        rwgt(nrpts) = 4.0_dp/12.0_dp
      else if(power==1) then
        do nr=2,nrpts-3,2
          rwgt(nr)=rab(nr)*16.0_dp/12.0_dp*rlog(nr)
          rwgt(nr+1)=rab(nr+1)*8.0_dp/12.0_dp*rlog(nr+1)
        end do

        rwgt(nrpts-1) = rab(nrpts-1)*16.0_dp/12.0_dp*rlog(nrpts-1)
        rwgt(nrpts) = 4.0_dp/12.0_dp*rlog(nrpts)
      else
        do nr=2,nrpts-3,2
          rwgt(nr)=rab(nr)*16.0_dp/12.0_dp*rlog(nr)**power
          rwgt(nr+1)=rab(nr+1)*8.0_dp/12.0_dp*rlog(nr+1)**power
        end do

        rwgt(nrpts-1) = rab(nrpts-1)*16.0_dp/12.0_dp*rlog(nrpts-1)**power
        rwgt(nrpts) = 4.0_dp/12.0_dp*rlog(nrpts)**power
      end if

      ! ** Loop over the radial reciprocal space grid

      do nq=1,nqpts
        q = real(nq-1,dp)*qmax/real(nqpts-1,dp)
        lookup(1,nq,1) = xrd_sf_bessfn_r(0.0_dp,l)*rwgt(1)
        do nr=2,nrpts
          qr = q*rlog(nr)
          lookup(nr,nq,1) = xrd_sf_bessfn_r(qr,l)*rwgt(nr)
        end do
      end do

      check(1) = real(nqpts,dp)
      check(2) = qmax
      check(3) = real(nrpts,dp)
      check(4) = rlog(1)
      check(5) = rab(1)
      check(6) = real(-l,dp)
      check(7) = real(power,dp)

    end if

    ! ** Do the matrix-vector multiplication

    call dgemv('T',nrpts,nqpts,1.0_dp,lookup(1,1,1),nrpts,rfunc,1,0.0_dp,qfunc,1)

  else

    current_check(1) = real(nqpts,dp)
    current_check(2) = qmax
    current_check(3) = real(nrpts,dp)
    current_check(4) = rlog(1)
    current_check(5) = rab(1)
    current_check(6) = real(lmax,dp)
    current_check(7) = real(power,dp)

    if((.not.allocated(lookup)).or.(.not.all(check==current_check))) then
      if (allocated(lookup)) then
        deallocate(lookup,stat=ierr)
        if (ierr/=0) call io_abort('Error in deallocating lookup in xrd_sf_radial_transform')
      end if
      allocate(lookup(nrpts,nqpts,0:lmax),stat=ierr)
      if (ierr/=0) call io_allocate_abort('lookup','xrd_sf_radial_transform')

      ! ** Set the integration weights

      rwgt(1) = 4.0_dp/12.0_dp

      ! Optimisation: trap special cases of power being 0 or 1, which are extremely common and allow
      ! significant simplifications
      if(power==0) then
        do nr=2,nrpts-3,2
          rwgt(nr)=rab(nr)*16.0_dp/12.0_dp
          rwgt(nr+1)=rab(nr+1)*8.0_dp/12.0_dp
        end do

        rwgt(nrpts-1) = rab(nrpts-1)*16.0_dp/12.0_dp
        rwgt(nrpts) = 4.0_dp/12.0_dp
      else if(power==1) then
        do nr=2,nrpts-3,2
          rwgt(nr)=rab(nr)*16.0_dp/12.0_dp*rlog(nr)
          rwgt(nr+1)=rab(nr+1)*8.0_dp/12.0_dp*rlog(nr+1)
        end do

        rwgt(nrpts-1) = rab(nrpts-1)*16.0_dp/12.0_dp*rlog(nrpts-1)
        rwgt(nrpts) = 4.0_dp/12.0_dp*rlog(nrpts)
      else
        do nr=2,nrpts-3,2
          rwgt(nr)=rab(nr)*16.0_dp/12.0_dp*rlog(nr)**power
          rwgt(nr+1)=rab(nr+1)*8.0_dp/12.0_dp*rlog(nr+1)**power
        end do

        rwgt(nrpts-1) = rab(nrpts-1)*16.0_dp/12.0_dp*rlog(nrpts-1)**power
        rwgt(nrpts) = 4.0_dp/12.0_dp*rlog(nrpts)**power
      end if

      ! ** Loop over the radial reciprocal space grid

      do nq=1,nqpts
        q = real(nq-1,dp)*qmax/real(nqpts-1,dp)
        do ll=0,lmax
          lookup(1,nq,ll) = xrd_sf_bessfn_r(0.0_dp,ll)*rwgt(1)
        end do
        do nr=2,nrpts
          qr = q*rlog(nr)
          do ll=0,lmax
            lookup(nr,nq,ll) = xrd_sf_bessfn_R(qr,ll)*rwgt(nr)
          end do
        end do
      end do

      check(1) = real(nqpts,dp)
      check(2) = qmax
      check(3) = real(nrpts,dp)
      check(4) = rlog(1)
      check(5) = rab(1)
      check(6) = real(lmax,dp)
      check(7) = real(power,dp)

    end if

    ! ** Do the matrix-vector multiplication

    call dgemv('T',nrpts,nqpts,1.0_dp,lookup(1,1,l),nrpts,rfunc,1,0.0_dp,qfunc,1)

  end if

  call trace_exit('xrd_sf_radial_transform',status)

end subroutine xrd_sf_radial_transform

function xrd_sf_bessfn_r(x,l)

  !=========================================================================!
  ! Returns the value of a spherical bessel function                        !
  !-------------------------------------------------------------------------!
  ! Arguments:                                                              !
  !                                                                         !
  ! 1) x : input : the value at which the bessel function is to be obtained !
  ! 2) l : input : which spherical harmonic is required                     !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Parent module variables used:                                           !
  !                                                                         !
  ! 0) None                                                                 !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Modules used:                                                           !
  !                                                                         !
  ! 0) None                                                                 !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Key Internal Variables:                                                 !
  !                                                                         !
  ! 1) IACC : determines accuracy of the result                             !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Necessary conditions:                                                   !
  !                                                                         !
  ! 0) None                                                                 !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Chris J. Pickard, V0.0, 8th June 2000                                   !
  ! This is essentially a clone of the routine from otfg                    !
  ! B.X. Shi reported that using the OTFG version sometimes gave NaNs       !
  ! This version is modified slightly. Eventually, this routine             !
  ! should be replaced with the more efficient multi version                !
  ! at such time we should test the OTFG version again                      !
  !=========================================================================!

  use trace, only: trace_entry, trace_exit

  implicit none

  integer,       intent(in) :: l
  real(kind=dp), intent(in) :: x

  real(kind=dp)             :: xrd_sf_bessfn_r
  real(kind=dp)             :: xrd_sf_bessfn_series
  integer, parameter :: bess_smax=1
  integer, parameter :: fact_lmax=10

  !-------------------------------------------------------------------------!

  integer                  :: i, s, m, isum
  real(kind=dp)            :: j,jm,jp,oox, summ
  real(kind=dp), parameter :: xsplit = 0.1_dp
  integer, parameter       :: IACC=40

  integer                  :: status

  call trace_entry('xrd_sf_bessfn_r',status)

  j = 0.0_dp
  jm = 0.0_dp
  jp = 0.0_dp
  oox = 0.0_dp

  if(x<10.0_dp*tiny(x)) then
    if(l.eq.0) then
      xrd_sf_bessfn_r=1.0_dp
    else
      xrd_sf_bessfn_r=0.0_dp
    end if
  else if(x>xsplit*l) then
    oox = 1.0_dp/x
    jm=sin(x)*oox
    j=(jm-cos(x))*oox
    if(l.eq.0) then
      xrd_sf_bessfn_r=jm
    else if (l.eq.1) then
      xrd_sf_bessfn_r=j
    else
      do i=1,l-1
        jp=real(2*i+1,dp)*oox*j-jm
        jm=j
        j=jp
      end do
      xrd_sf_bessfn_r=j
    end if
  else
    oox = 1.0_dp/x
    m=2*((l+int(sqrt(real(IACC*l,dp))))/2)
    xrd_sf_bessfn_series = 0.0_dp
    isum=0
    summ=0
    jp=0
    j=1
    do i=m,1,-1
      jm=real(2*i+1,dp)*oox*j-jp
      jp=j
      j=jm
      if(isum/=0) summ=summ+j
      isum=1-isum
      if(i==l) xrd_sf_bessfn_series=jp
    end do
    summ=2.0_dp*summ-j
    xrd_sf_bessfn_series=xrd_sf_bessfn_series/summ

    xrd_sf_bessfn_r = xrd_sf_bessfn_series
  end if

  call trace_exit('xrd_sf_bessfn_r',status)

end function xrd_sf_bessfn_r
subroutine xrd_sf_gaunt_coefficients(aug_coefficient)

  use trace, only: trace_entry, trace_exit

  implicit none

  real(kind=dp), dimension(0:3,-3:3,0:3,-3:3,0:6,-6:6), intent(inout) :: aug_coefficient

  !=========================================================================!
  ! This subroutine contains is a look-up table for Gaunt Coefficients.     !
  ! familiar with these beasts, this would be the integral over spherical   !
  ! harmonics YlmYkqYlpmp. These particular coefficients were generated     !
  ! using the Gaunt and threej routines from SymPy in Python.               !
  !-------------------------------------------------------------------------!
  ! Arguments:                                                              !
  !   Table is of the form Coefficient = GC_table(l1, m1, l2, m2, l3, m3)   !
  !-------------------------------------------------------------------------!
  ! Modules used:                                                           !
  !-------------------------------------------------------------------------!
  ! Key Internal Variables:                                                 !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Necessary conditions:                                                   !
  !   None                                                                  !
  !-------------------------------------------------------------------------!
  ! Written by Benjamin Shi, Version 1, 15/09/2020                          !
  ! NOTE: They have thusfar only been generated for l1 = 0,3, l2 = 0,3      !
  ! and l3 = 0,6.                                                           !
  ! This is for an integral over Y*l1m1Yl2m2Y*l3m3, which is different      !
  ! from paw_recon_gaunt as here, the complex conjugate is on the first and !
  ! third spherical harmonics                                               !
  !=========================================================================!

  integer                         :: status

  call trace_entry('xrd_sf_gaunt_coefficients',status)

  aug_coefficient=0.0_dp
  aug_coefficient( 0, 0, 0, 0, 0, 0) =  0.282094791774d0
  aug_coefficient( 0, 0, 1,-1, 1,-1) =  0.282094791774d0
  aug_coefficient( 0, 0, 1, 0, 1, 0) =  0.282094791774d0
  aug_coefficient( 0, 0, 1, 1, 1, 1) =  0.282094791774d0
  aug_coefficient( 0, 0, 2,-2, 2,-2) =  0.282094791774d0
  aug_coefficient( 0, 0, 2,-1, 2,-1) =  0.282094791774d0
  aug_coefficient( 0, 0, 2, 0, 2, 0) =  0.282094791774d0
  aug_coefficient( 0, 0, 2, 1, 2, 1) =  0.282094791774d0
  aug_coefficient( 0, 0, 2, 2, 2, 2) =  0.282094791774d0
  aug_coefficient( 0, 0, 3,-3, 3,-3) =  0.282094791774d0
  aug_coefficient( 0, 0, 3,-2, 3,-2) =  0.282094791774d0
  aug_coefficient( 0, 0, 3,-1, 3,-1) =  0.282094791774d0
  aug_coefficient( 0, 0, 3, 0, 3, 0) =  0.282094791774d0
  aug_coefficient( 0, 0, 3, 1, 3, 1) =  0.282094791774d0
  aug_coefficient( 0, 0, 3, 2, 3, 2) =  0.282094791774d0
  aug_coefficient( 0, 0, 3, 3, 3, 3) =  0.282094791774d0
  aug_coefficient( 1,-1, 0, 0, 1, 1) = -0.282094791774d0
  aug_coefficient( 1,-1, 1,-1, 0, 0) =  0.282094791774d0
  aug_coefficient( 1,-1, 1,-1, 2, 0) = -0.126156626101d0
  aug_coefficient( 1,-1, 1, 0, 2, 1) = -0.218509686118d0
  aug_coefficient( 1,-1, 1, 1, 2, 2) = -0.309019361619d0
  aug_coefficient( 1,-1, 2,-2, 1,-1) =  0.309019361619d0
  aug_coefficient( 1,-1, 2,-2, 3,-1) = -0.082588898361d0
  aug_coefficient( 1,-1, 2,-1, 1, 0) =  0.218509686118d0
  aug_coefficient( 1,-1, 2,-1, 3, 0) = -0.143048168103d0
  aug_coefficient( 1,-1, 2, 0, 1, 1) =  0.126156626101d0
  aug_coefficient( 1,-1, 2, 0, 3, 1) = -0.202300659403d0
  aug_coefficient( 1,-1, 2, 1, 3, 2) = -0.261169028265d0
  aug_coefficient( 1,-1, 2, 2, 3, 3) = -0.319865427934d0
  aug_coefficient( 1,-1, 3,-3, 2,-2) =  0.319865427934d0
  aug_coefficient( 1,-1, 3,-3, 4,-2) = -0.061558130307d0
  aug_coefficient( 1,-1, 3,-2, 2,-1) =  0.261169028265d0
  aug_coefficient( 1,-1, 3,-2, 4,-1) = -0.106621809311d0
  aug_coefficient( 1,-1, 3,-1, 2, 0) =  0.202300659403d0
  aug_coefficient( 1,-1, 3,-1, 4, 0) = -0.150786008773d0
  aug_coefficient( 1,-1, 3, 0, 2, 1) =  0.143048168103d0
  aug_coefficient( 1,-1, 3, 0, 4, 1) = -0.194663900273d0
  aug_coefficient( 1,-1, 3, 1, 2, 2) =  0.082588898361d0
  aug_coefficient( 1,-1, 3, 1, 4, 2) = -0.238413613504d0
  aug_coefficient( 1,-1, 3, 2, 4, 3) = -0.282094791774d0
  aug_coefficient( 1,-1, 3, 3, 4, 4) = -0.325735007935d0
  aug_coefficient( 1, 0, 0, 0, 1, 0) =  0.282094791774d0
  aug_coefficient( 1, 0, 1,-1, 2,-1) =  0.218509686118d0
  aug_coefficient( 1, 0, 1, 0, 0, 0) =  0.282094791774d0
  aug_coefficient( 1, 0, 1, 0, 2, 0) =  0.252313252202d0
  aug_coefficient( 1, 0, 1, 1, 2, 1) =  0.218509686118d0
  aug_coefficient( 1, 0, 2,-2, 3,-2) =  0.184674390922d0
  aug_coefficient( 1, 0, 2,-1, 1,-1) =  0.218509686118d0
  aug_coefficient( 1, 0, 2,-1, 3,-1) =  0.233596680328d0
  aug_coefficient( 1, 0, 2, 0, 1, 0) =  0.252313252202d0
  aug_coefficient( 1, 0, 2, 0, 3, 0) =  0.247766695083d0
  aug_coefficient( 1, 0, 2, 1, 1, 1) =  0.218509686118d0
  aug_coefficient( 1, 0, 2, 1, 3, 1) =  0.233596680328d0
  aug_coefficient( 1, 0, 2, 2, 3, 2) =  0.184674390922d0
  aug_coefficient( 1, 0, 3,-3, 4,-3) =  0.162867503968d0
  aug_coefficient( 1, 0, 3,-2, 2,-2) =  0.184674390922d0
  aug_coefficient( 1, 0, 3,-2, 4,-2) =  0.213243618623d0
  aug_coefficient( 1, 0, 3,-1, 2,-1) =  0.233596680328d0
  aug_coefficient( 1, 0, 3,-1, 4,-1) =  0.238413613504d0
  aug_coefficient( 1, 0, 3, 0, 2, 0) =  0.247766695083d0
  aug_coefficient( 1, 0, 3, 0, 4, 0) =  0.246232521230d0
  aug_coefficient( 1, 0, 3, 1, 2, 1) =  0.233596680328d0
  aug_coefficient( 1, 0, 3, 1, 4, 1) =  0.238413613504d0
  aug_coefficient( 1, 0, 3, 2, 2, 2) =  0.184674390922d0
  aug_coefficient( 1, 0, 3, 2, 4, 2) =  0.213243618623d0
  aug_coefficient( 1, 0, 3, 3, 4, 3) =  0.162867503968d0
  aug_coefficient( 1, 1, 0, 0, 1,-1) = -0.282094791774d0
  aug_coefficient( 1, 1, 1,-1, 2,-2) = -0.309019361619d0
  aug_coefficient( 1, 1, 1, 0, 2,-1) = -0.218509686118d0
  aug_coefficient( 1, 1, 1, 1, 0, 0) =  0.282094791774d0
  aug_coefficient( 1, 1, 1, 1, 2, 0) = -0.126156626101d0
  aug_coefficient( 1, 1, 2,-2, 3,-3) = -0.319865427934d0
  aug_coefficient( 1, 1, 2,-1, 3,-2) = -0.261169028265d0
  aug_coefficient( 1, 1, 2, 0, 1,-1) =  0.126156626101d0
  aug_coefficient( 1, 1, 2, 0, 3,-1) = -0.202300659403d0
  aug_coefficient( 1, 1, 2, 1, 1, 0) =  0.218509686118d0
  aug_coefficient( 1, 1, 2, 1, 3, 0) = -0.143048168103d0
  aug_coefficient( 1, 1, 2, 2, 1, 1) =  0.309019361619d0
  aug_coefficient( 1, 1, 2, 2, 3, 1) = -0.082588898361d0
  aug_coefficient( 1, 1, 3,-3, 4,-4) = -0.325735007935d0
  aug_coefficient( 1, 1, 3,-2, 4,-3) = -0.282094791774d0
  aug_coefficient( 1, 1, 3,-1, 2,-2) =  0.082588898361d0
  aug_coefficient( 1, 1, 3,-1, 4,-2) = -0.238413613504d0
  aug_coefficient( 1, 1, 3, 0, 2,-1) =  0.143048168103d0
  aug_coefficient( 1, 1, 3, 0, 4,-1) = -0.194663900273d0
  aug_coefficient( 1, 1, 3, 1, 2, 0) =  0.202300659403d0
  aug_coefficient( 1, 1, 3, 1, 4, 0) = -0.150786008773d0
  aug_coefficient( 1, 1, 3, 2, 2, 1) =  0.261169028265d0
  aug_coefficient( 1, 1, 3, 2, 4, 1) = -0.106621809311d0
  aug_coefficient( 1, 1, 3, 3, 2, 2) =  0.319865427934d0
  aug_coefficient( 1, 1, 3, 3, 4, 2) = -0.061558130307d0
  aug_coefficient( 2,-2, 0, 0, 2, 2) =  0.282094791774d0
  aug_coefficient( 2,-2, 1,-1, 1, 1) = -0.309019361619d0
  aug_coefficient( 2,-2, 1,-1, 3, 1) =  0.082588898361d0
  aug_coefficient( 2,-2, 1, 0, 3, 2) =  0.184674390922d0
  aug_coefficient( 2,-2, 1, 1, 3, 3) =  0.319865427934d0
  aug_coefficient( 2,-2, 2,-2, 0, 0) =  0.282094791774d0
  aug_coefficient( 2,-2, 2,-2, 2, 0) = -0.180223751573d0
  aug_coefficient( 2,-2, 2,-2, 4, 0) =  0.040299255968d0
  aug_coefficient( 2,-2, 2,-1, 2, 1) = -0.220728115442d0
  aug_coefficient( 2,-2, 2,-1, 4, 1) =  0.090111875786d0
  aug_coefficient( 2,-2, 2, 0, 2, 2) = -0.180223751573d0
  aug_coefficient( 2,-2, 2, 0, 4, 2) =  0.156078347227d0
  aug_coefficient( 2,-2, 2, 1, 4, 3) =  0.238413613504d0
  aug_coefficient( 2,-2, 2, 2, 4, 4) =  0.337167765672d0
  aug_coefficient( 2,-2, 3,-3, 1,-1) =  0.319865427934d0
  aug_coefficient( 2,-2, 3,-3, 3,-1) = -0.132980760134d0
  aug_coefficient( 2,-2, 3,-3, 5,-1) =  0.023961469724d0
  aug_coefficient( 2,-2, 3,-2, 1, 0) =  0.184674390922d0
  aug_coefficient( 2,-2, 3,-2, 3, 0) = -0.188063194516d0
  aug_coefficient( 2,-2, 3,-2, 5, 0) =  0.053579475145d0
  aug_coefficient( 2,-2, 3,-1, 1, 1) =  0.082588898361d0
  aug_coefficient( 2,-2, 3,-1, 3, 1) = -0.206012907746d0
  aug_coefficient( 2,-2, 3,-1, 5, 1) =  0.092802373193d0
  aug_coefficient( 2,-2, 3, 0, 3, 2) = -0.188063194516d0
  aug_coefficient( 2,-2, 3, 0, 5, 2) =  0.141757966610d0
  aug_coefficient( 2,-2, 3, 1, 3, 3) = -0.132980760134d0
  aug_coefficient( 2,-2, 3, 1, 5, 3) =  0.200476038955d0
  aug_coefficient( 2,-2, 3, 2, 5, 4) =  0.268966830577d0
  aug_coefficient( 2,-2, 3, 3, 5, 5) =  0.347234685170d0
  aug_coefficient( 2,-1, 0, 0, 2, 1) = -0.282094791774d0
  aug_coefficient( 2,-1, 1,-1, 1, 0) =  0.218509686118d0
  aug_coefficient( 2,-1, 1,-1, 3, 0) = -0.143048168103d0
  aug_coefficient( 2,-1, 1, 0, 1, 1) = -0.218509686118d0
  aug_coefficient( 2,-1, 1, 0, 3, 1) = -0.233596680328d0
  aug_coefficient( 2,-1, 1, 1, 3, 2) = -0.261169028265d0
  aug_coefficient( 2,-1, 2,-2, 2,-1) =  0.220728115442d0
  aug_coefficient( 2,-1, 2,-2, 4,-1) = -0.090111875786d0
  aug_coefficient( 2,-1, 2,-1, 0, 0) =  0.282094791774d0
  aug_coefficient( 2,-1, 2,-1, 2, 0) =  0.090111875786d0
  aug_coefficient( 2,-1, 2,-1, 4, 0) = -0.161197023871d0
  aug_coefficient( 2,-1, 2, 0, 2, 1) = -0.090111875786d0
  aug_coefficient( 2,-1, 2, 0, 4, 1) = -0.220728115442d0
  aug_coefficient( 2,-1, 2, 1, 2, 2) = -0.220728115442d0
  aug_coefficient( 2,-1, 2, 1, 4, 2) = -0.254874873736d0
  aug_coefficient( 2,-1, 2, 2, 4, 3) = -0.238413613504d0
  aug_coefficient( 2,-1, 3,-3, 3,-2) =  0.210261043502d0
  aug_coefficient( 2,-1, 3,-3, 5,-2) = -0.063396089939d0
  aug_coefficient( 2,-1, 3,-2, 1,-1) =  0.261169028265d0
  aug_coefficient( 2,-1, 3,-2, 3,-1) =  0.162867503968d0
  aug_coefficient( 2,-1, 3,-2, 5,-1) = -0.117386748624d0
  aug_coefficient( 2,-1, 3,-1, 1, 0) =  0.233596680328d0
  aug_coefficient( 2,-1, 3,-1, 3, 0) =  0.059470803872d0
  aug_coefficient( 2,-1, 3,-1, 5, 0) = -0.169433177294d0
  aug_coefficient( 2,-1, 3, 0, 1, 1) =  0.143048168103d0
  aug_coefficient( 2,-1, 3, 0, 3, 1) = -0.059470803872d0
  aug_coefficient( 2,-1, 3, 0, 5, 1) = -0.214317900579d0
  aug_coefficient( 2,-1, 3, 1, 3, 2) = -0.162867503968d0
  aug_coefficient( 2,-1, 3, 1, 5, 2) = -0.245532000547d0
  aug_coefficient( 2,-1, 3, 2, 3, 3) = -0.210261043502d0
  aug_coefficient( 2,-1, 3, 2, 5, 3) = -0.253584359754d0
  aug_coefficient( 2,-1, 3, 3, 5, 4) = -0.219610497549d0
  aug_coefficient( 2, 0, 0, 0, 2, 0) =  0.282094791774d0
  aug_coefficient( 2, 0, 1,-1, 1,-1) = -0.126156626101d0
  aug_coefficient( 2, 0, 1,-1, 3,-1) =  0.202300659403d0
  aug_coefficient( 2, 0, 1, 0, 1, 0) =  0.252313252202d0
  aug_coefficient( 2, 0, 1, 0, 3, 0) =  0.247766695083d0
  aug_coefficient( 2, 0, 1, 1, 1, 1) = -0.126156626101d0
  aug_coefficient( 2, 0, 1, 1, 3, 1) =  0.202300659403d0
  aug_coefficient( 2, 0, 2,-2, 2,-2) = -0.180223751573d0
  aug_coefficient( 2, 0, 2,-2, 4,-2) =  0.156078347227d0
  aug_coefficient( 2, 0, 2,-1, 2,-1) =  0.090111875786d0
  aug_coefficient( 2, 0, 2,-1, 4,-1) =  0.220728115442d0
  aug_coefficient( 2, 0, 2, 0, 0, 0) =  0.282094791774d0
  aug_coefficient( 2, 0, 2, 0, 2, 0) =  0.180223751573d0
  aug_coefficient( 2, 0, 2, 0, 4, 0) =  0.241795535806d0
  aug_coefficient( 2, 0, 2, 1, 2, 1) =  0.090111875786d0
  aug_coefficient( 2, 0, 2, 1, 4, 1) =  0.220728115442d0
  aug_coefficient( 2, 0, 2, 2, 2, 2) = -0.180223751573d0
  aug_coefficient( 2, 0, 2, 2, 4, 2) =  0.156078347227d0
  aug_coefficient( 2, 0, 3,-3, 3,-3) = -0.210261043502d0
  aug_coefficient( 2, 0, 3,-3, 5,-3) =  0.126792179877d0
  aug_coefficient( 2, 0, 3,-2, 5,-2) =  0.190188269816d0
  aug_coefficient( 2, 0, 3,-1, 1,-1) =  0.202300659403d0
  aug_coefficient( 2, 0, 3,-1, 3,-1) =  0.126156626101d0
  aug_coefficient( 2, 0, 3,-1, 5,-1) =  0.227318461243d0
  aug_coefficient( 2, 0, 3, 0, 1, 0) =  0.247766695083d0
  aug_coefficient( 2, 0, 3, 0, 3, 0) =  0.168208834801d0
  aug_coefficient( 2, 0, 3, 0, 5, 0) =  0.239614697245d0
  aug_coefficient( 2, 0, 3, 1, 1, 1) =  0.202300659403d0
  aug_coefficient( 2, 0, 3, 1, 3, 1) =  0.126156626101d0
  aug_coefficient( 2, 0, 3, 1, 5, 1) =  0.227318461243d0
  aug_coefficient( 2, 0, 3, 2, 5, 2) =  0.190188269816d0
  aug_coefficient( 2, 0, 3, 3, 3, 3) = -0.210261043502d0
  aug_coefficient( 2, 0, 3, 3, 5, 3) =  0.126792179877d0
  aug_coefficient( 2, 1, 0, 0, 2,-1) = -0.282094791774d0
  aug_coefficient( 2, 1, 1,-1, 3,-2) = -0.261169028265d0
  aug_coefficient( 2, 1, 1, 0, 1,-1) = -0.218509686118d0
  aug_coefficient( 2, 1, 1, 0, 3,-1) = -0.233596680328d0
  aug_coefficient( 2, 1, 1, 1, 1, 0) =  0.218509686118d0
  aug_coefficient( 2, 1, 1, 1, 3, 0) = -0.143048168103d0
  aug_coefficient( 2, 1, 2,-2, 4,-3) = -0.238413613504d0
  aug_coefficient( 2, 1, 2,-1, 2,-2) = -0.220728115442d0
  aug_coefficient( 2, 1, 2,-1, 4,-2) = -0.254874873736d0
  aug_coefficient( 2, 1, 2, 0, 2,-1) = -0.090111875786d0
  aug_coefficient( 2, 1, 2, 0, 4,-1) = -0.220728115442d0
  aug_coefficient( 2, 1, 2, 1, 0, 0) =  0.282094791774d0
  aug_coefficient( 2, 1, 2, 1, 2, 0) =  0.090111875786d0
  aug_coefficient( 2, 1, 2, 1, 4, 0) = -0.161197023871d0
  aug_coefficient( 2, 1, 2, 2, 2, 1) =  0.220728115442d0
  aug_coefficient( 2, 1, 2, 2, 4, 1) = -0.090111875786d0
  aug_coefficient( 2, 1, 3,-3, 5,-4) = -0.219610497549d0
  aug_coefficient( 2, 1, 3,-2, 3,-3) = -0.210261043502d0
  aug_coefficient( 2, 1, 3,-2, 5,-3) = -0.253584359754d0
  aug_coefficient( 2, 1, 3,-1, 3,-2) = -0.162867503968d0
  aug_coefficient( 2, 1, 3,-1, 5,-2) = -0.245532000547d0
  aug_coefficient( 2, 1, 3, 0, 1,-1) =  0.143048168103d0
  aug_coefficient( 2, 1, 3, 0, 3,-1) = -0.059470803872d0
  aug_coefficient( 2, 1, 3, 0, 5,-1) = -0.214317900579d0
  aug_coefficient( 2, 1, 3, 1, 1, 0) =  0.233596680328d0
  aug_coefficient( 2, 1, 3, 1, 3, 0) =  0.059470803872d0
  aug_coefficient( 2, 1, 3, 1, 5, 0) = -0.169433177294d0
  aug_coefficient( 2, 1, 3, 2, 1, 1) =  0.261169028265d0
  aug_coefficient( 2, 1, 3, 2, 3, 1) =  0.162867503968d0
  aug_coefficient( 2, 1, 3, 2, 5, 1) = -0.117386748624d0
  aug_coefficient( 2, 1, 3, 3, 3, 2) =  0.210261043502d0
  aug_coefficient( 2, 1, 3, 3, 5, 2) = -0.063396089939d0
  aug_coefficient( 2, 2, 0, 0, 2,-2) =  0.282094791774d0
  aug_coefficient( 2, 2, 1,-1, 3,-3) =  0.319865427934d0
  aug_coefficient( 2, 2, 1, 0, 3,-2) =  0.184674390922d0
  aug_coefficient( 2, 2, 1, 1, 1,-1) = -0.309019361619d0
  aug_coefficient( 2, 2, 1, 1, 3,-1) =  0.082588898361d0
  aug_coefficient( 2, 2, 2,-2, 4,-4) =  0.337167765672d0
  aug_coefficient( 2, 2, 2,-1, 4,-3) =  0.238413613504d0
  aug_coefficient( 2, 2, 2, 0, 2,-2) = -0.180223751573d0
  aug_coefficient( 2, 2, 2, 0, 4,-2) =  0.156078347227d0
  aug_coefficient( 2, 2, 2, 1, 2,-1) = -0.220728115442d0
  aug_coefficient( 2, 2, 2, 1, 4,-1) =  0.090111875786d0
  aug_coefficient( 2, 2, 2, 2, 0, 0) =  0.282094791774d0
  aug_coefficient( 2, 2, 2, 2, 2, 0) = -0.180223751573d0
  aug_coefficient( 2, 2, 2, 2, 4, 0) =  0.040299255968d0
  aug_coefficient( 2, 2, 3,-3, 5,-5) =  0.347234685170d0
  aug_coefficient( 2, 2, 3,-2, 5,-4) =  0.268966830577d0
  aug_coefficient( 2, 2, 3,-1, 3,-3) = -0.132980760134d0
  aug_coefficient( 2, 2, 3,-1, 5,-3) =  0.200476038955d0
  aug_coefficient( 2, 2, 3, 0, 3,-2) = -0.188063194516d0
  aug_coefficient( 2, 2, 3, 0, 5,-2) =  0.141757966610d0
  aug_coefficient( 2, 2, 3, 1, 1,-1) =  0.082588898361d0
  aug_coefficient( 2, 2, 3, 1, 3,-1) = -0.206012907746d0
  aug_coefficient( 2, 2, 3, 1, 5,-1) =  0.092802373193d0
  aug_coefficient( 2, 2, 3, 2, 1, 0) =  0.184674390922d0
  aug_coefficient( 2, 2, 3, 2, 3, 0) = -0.188063194516d0
  aug_coefficient( 2, 2, 3, 2, 5, 0) =  0.053579475145d0
  aug_coefficient( 2, 2, 3, 3, 1, 1) =  0.319865427934d0
  aug_coefficient( 2, 2, 3, 3, 3, 1) = -0.132980760134d0
  aug_coefficient( 2, 2, 3, 3, 5, 1) =  0.023961469724d0
  aug_coefficient( 3,-3, 0, 0, 3, 3) = -0.282094791774d0
  aug_coefficient( 3,-3, 1,-1, 2, 2) =  0.319865427934d0
  aug_coefficient( 3,-3, 1,-1, 4, 2) = -0.061558130307d0
  aug_coefficient( 3,-3, 1, 0, 4, 3) = -0.162867503968d0
  aug_coefficient( 3,-3, 1, 1, 4, 4) = -0.325735007935d0
  aug_coefficient( 3,-3, 2,-2, 1, 1) = -0.319865427934d0
  aug_coefficient( 3,-3, 2,-2, 3, 1) =  0.132980760134d0
  aug_coefficient( 3,-3, 2,-2, 5, 1) = -0.023961469724d0
  aug_coefficient( 3,-3, 2,-1, 3, 2) =  0.210261043502d0
  aug_coefficient( 3,-3, 2,-1, 5, 2) = -0.063396089939d0
  aug_coefficient( 3,-3, 2, 0, 3, 3) =  0.210261043502d0
  aug_coefficient( 3,-3, 2, 0, 5, 3) = -0.126792179877d0
  aug_coefficient( 3,-3, 2, 1, 5, 4) = -0.219610497549d0
  aug_coefficient( 3,-3, 2, 2, 5, 5) = -0.347234685170d0
  aug_coefficient( 3,-3, 3,-3, 0, 0) =  0.282094791774d0
  aug_coefficient( 3,-3, 3,-3, 2, 0) = -0.210261043502d0
  aug_coefficient( 3,-3, 3,-3, 4, 0) =  0.076934943211d0
  aug_coefficient( 3,-3, 3,-3, 6, 0) = -0.011854396693d0
  aug_coefficient( 3,-3, 3,-2, 2, 1) = -0.210261043502d0
  aug_coefficient( 3,-3, 3,-2, 4, 1) =  0.140463346190d0
  aug_coefficient( 3,-3, 3,-2, 6, 1) = -0.031363785593d0
  aug_coefficient( 3,-3, 3,-1, 2, 2) = -0.132980760134d0
  aug_coefficient( 3,-3, 3,-1, 4, 2) =  0.188451354257d0
  aug_coefficient( 3,-3, 3,-1, 6, 2) = -0.062727571186d0
  aug_coefficient( 3,-3, 3, 0, 4, 3) =  0.203550726867d0
  aug_coefficient( 3,-3, 3, 0, 6, 3) = -0.108647340330d0
  aug_coefficient( 3,-3, 3, 1, 4, 4) =  0.166198472533d0
  aug_coefficient( 3,-3, 3, 1, 6, 4) = -0.171786528581d0
  aug_coefficient( 3,-3, 3, 2, 6, 5) = -0.254800598673d0
  aug_coefficient( 3,-3, 3, 3, 6, 6) = -0.360342462344d0
  aug_coefficient( 3,-2, 0, 0, 3, 2) =  0.282094791774d0
  aug_coefficient( 3,-2, 1,-1, 2, 1) = -0.261169028265d0
  aug_coefficient( 3,-2, 1,-1, 4, 1) =  0.106621809311d0
  aug_coefficient( 3,-2, 1, 0, 2, 2) =  0.184674390922d0
  aug_coefficient( 3,-2, 1, 0, 4, 2) =  0.213243618623d0
  aug_coefficient( 3,-2, 1, 1, 4, 3) =  0.282094791774d0
  aug_coefficient( 3,-2, 2,-2, 1, 0) =  0.184674390922d0
  aug_coefficient( 3,-2, 2,-2, 3, 0) = -0.188063194516d0
  aug_coefficient( 3,-2, 2,-2, 5, 0) =  0.053579475145d0
  aug_coefficient( 3,-2, 2,-1, 1, 1) = -0.261169028265d0
  aug_coefficient( 3,-2, 2,-1, 3, 1) = -0.162867503968d0
  aug_coefficient( 3,-2, 2,-1, 5, 1) =  0.117386748624d0
  aug_coefficient( 3,-2, 2, 0, 5, 2) =  0.190188269816d0
  aug_coefficient( 3,-2, 2, 1, 3, 3) =  0.210261043502d0
  aug_coefficient( 3,-2, 2, 1, 5, 3) =  0.253584359754d0
  aug_coefficient( 3,-2, 2, 2, 5, 4) =  0.268966830577d0
  aug_coefficient( 3,-2, 3,-3, 2,-1) =  0.210261043502d0
  aug_coefficient( 3,-2, 3,-3, 4,-1) = -0.140463346190d0
  aug_coefficient( 3,-2, 3,-3, 6,-1) =  0.031363785593d0
  aug_coefficient( 3,-2, 3,-2, 0, 0) =  0.282094791774d0
  aug_coefficient( 3,-2, 3,-2, 4, 0) = -0.179514867492d0
  aug_coefficient( 3,-2, 3,-2, 6, 0) =  0.071126380160d0
  aug_coefficient( 3,-2, 3,-1, 2, 1) = -0.162867503968d0
  aug_coefficient( 3,-2, 3,-1, 4, 1) = -0.145069920146d0
  aug_coefficient( 3,-2, 3,-1, 6, 1) =  0.121471419276d0
  aug_coefficient( 3,-2, 3, 0, 2, 2) = -0.188063194516d0
  aug_coefficient( 3,-2, 3, 0, 4, 2) = -0.044418410173d0
  aug_coefficient( 3,-2, 3, 0, 6, 2) =  0.177420363812d0
  aug_coefficient( 3,-2, 3, 1, 4, 3) =  0.095954732856d0
  aug_coefficient( 3,-2, 3, 1, 6, 3) =  0.230475813315d0
  aug_coefficient( 3,-2, 3, 2, 4, 4) =  0.214561305428d0
  aug_coefficient( 3,-2, 3, 2, 6, 4) =  0.266130545719d0
  aug_coefficient( 3,-2, 3, 3, 6, 5) =  0.254800598673d0
  aug_coefficient( 3,-1, 0, 0, 3, 1) = -0.282094791774d0
  aug_coefficient( 3,-1, 1,-1, 2, 0) =  0.202300659403d0
  aug_coefficient( 3,-1, 1,-1, 4, 0) = -0.150786008773d0
  aug_coefficient( 3,-1, 1, 0, 2, 1) = -0.233596680328d0
  aug_coefficient( 3,-1, 1, 0, 4, 1) = -0.238413613504d0
  aug_coefficient( 3,-1, 1, 1, 2, 2) =  0.082588898361d0
  aug_coefficient( 3,-1, 1, 1, 4, 2) = -0.238413613504d0
  aug_coefficient( 3,-1, 2,-2, 1,-1) = -0.082588898361d0
  aug_coefficient( 3,-1, 2,-2, 3,-1) =  0.206012907746d0
  aug_coefficient( 3,-1, 2,-2, 5,-1) = -0.092802373193d0
  aug_coefficient( 3,-1, 2,-1, 1, 0) =  0.233596680328d0
  aug_coefficient( 3,-1, 2,-1, 3, 0) =  0.059470803872d0
  aug_coefficient( 3,-1, 2,-1, 5, 0) = -0.169433177294d0
  aug_coefficient( 3,-1, 2, 0, 1, 1) = -0.202300659403d0
  aug_coefficient( 3,-1, 2, 0, 3, 1) = -0.126156626101d0
  aug_coefficient( 3,-1, 2, 0, 5, 1) = -0.227318461243d0
  aug_coefficient( 3,-1, 2, 1, 3, 2) = -0.162867503968d0
  aug_coefficient( 3,-1, 2, 1, 5, 2) = -0.245532000547d0
  aug_coefficient( 3,-1, 2, 2, 3, 3) =  0.132980760134d0
  aug_coefficient( 3,-1, 2, 2, 5, 3) = -0.200476038955d0
  aug_coefficient( 3,-1, 3,-3, 2,-2) = -0.132980760134d0
  aug_coefficient( 3,-1, 3,-3, 4,-2) =  0.188451354257d0
  aug_coefficient( 3,-1, 3,-3, 6,-2) = -0.062727571186d0
  aug_coefficient( 3,-1, 3,-2, 2,-1) =  0.162867503968d0
  aug_coefficient( 3,-1, 3,-2, 4,-1) =  0.145069920146d0
  aug_coefficient( 3,-1, 3,-2, 6,-1) = -0.121471419276d0
  aug_coefficient( 3,-1, 3,-1, 0, 0) =  0.282094791774d0
  aug_coefficient( 3,-1, 3,-1, 2, 0) =  0.126156626101d0
  aug_coefficient( 3,-1, 3,-1, 4, 0) =  0.025644981070d0
  aug_coefficient( 3,-1, 3,-1, 6, 0) = -0.177815950399d0
  aug_coefficient( 3,-1, 3, 0, 2, 1) = -0.059470803872d0
  aug_coefficient( 3,-1, 3, 0, 4, 1) = -0.099322584599d0
  aug_coefficient( 3,-1, 3, 0, 6, 1) = -0.221775454765d0
  aug_coefficient( 3,-1, 3, 1, 2, 2) = -0.206012907746d0
  aug_coefficient( 3,-1, 3, 1, 4, 2) = -0.162193101468d0
  aug_coefficient( 3,-1, 3, 1, 6, 2) = -0.242942838552d0
  aug_coefficient( 3,-1, 3, 2, 4, 3) = -0.095954732856d0
  aug_coefficient( 3,-1, 3, 2, 6, 3) = -0.230475813315d0
  aug_coefficient( 3,-1, 3, 3, 4, 4) =  0.166198472533d0
  aug_coefficient( 3,-1, 3, 3, 6, 4) = -0.171786528581d0
  aug_coefficient( 3, 0, 0, 0, 3, 0) =  0.282094791774d0
  aug_coefficient( 3, 0, 1,-1, 2,-1) = -0.143048168103d0
  aug_coefficient( 3, 0, 1,-1, 4,-1) =  0.194663900273d0
  aug_coefficient( 3, 0, 1, 0, 2, 0) =  0.247766695083d0
  aug_coefficient( 3, 0, 1, 0, 4, 0) =  0.246232521230d0
  aug_coefficient( 3, 0, 1, 1, 2, 1) = -0.143048168103d0
  aug_coefficient( 3, 0, 1, 1, 4, 1) =  0.194663900273d0
  aug_coefficient( 3, 0, 2,-2, 3,-2) = -0.188063194516d0
  aug_coefficient( 3, 0, 2,-2, 5,-2) =  0.141757966610d0
  aug_coefficient( 3, 0, 2,-1, 1,-1) = -0.143048168103d0
  aug_coefficient( 3, 0, 2,-1, 3,-1) =  0.059470803872d0
  aug_coefficient( 3, 0, 2,-1, 5,-1) =  0.214317900579d0
  aug_coefficient( 3, 0, 2, 0, 1, 0) =  0.247766695083d0
  aug_coefficient( 3, 0, 2, 0, 3, 0) =  0.168208834801d0
  aug_coefficient( 3, 0, 2, 0, 5, 0) =  0.239614697245d0
  aug_coefficient( 3, 0, 2, 1, 1, 1) = -0.143048168103d0
  aug_coefficient( 3, 0, 2, 1, 3, 1) =  0.059470803872d0
  aug_coefficient( 3, 0, 2, 1, 5, 1) =  0.214317900579d0
  aug_coefficient( 3, 0, 2, 2, 3, 2) = -0.188063194516d0
  aug_coefficient( 3, 0, 2, 2, 5, 2) =  0.141757966610d0
  aug_coefficient( 3, 0, 3,-3, 4,-3) = -0.203550726867d0
  aug_coefficient( 3, 0, 3,-3, 6,-3) =  0.108647340330d0
  aug_coefficient( 3, 0, 3,-2, 2,-2) = -0.188063194516d0
  aug_coefficient( 3, 0, 3,-2, 4,-2) = -0.044418410173d0
  aug_coefficient( 3, 0, 3,-2, 6,-2) =  0.177420363812d0
  aug_coefficient( 3, 0, 3,-1, 2,-1) =  0.059470803872d0
  aug_coefficient( 3, 0, 3,-1, 4,-1) =  0.099322584599d0
  aug_coefficient( 3, 0, 3,-1, 6,-1) =  0.221775454765d0
  aug_coefficient( 3, 0, 3, 0, 0, 0) =  0.282094791774d0
  aug_coefficient( 3, 0, 3, 0, 2, 0) =  0.168208834801d0
  aug_coefficient( 3, 0, 3, 0, 4, 0) =  0.153869886422d0
  aug_coefficient( 3, 0, 3, 0, 6, 0) =  0.237087933865d0
  aug_coefficient( 3, 0, 3, 1, 2, 1) =  0.059470803872d0
  aug_coefficient( 3, 0, 3, 1, 4, 1) =  0.099322584599d0
  aug_coefficient( 3, 0, 3, 1, 6, 1) =  0.221775454765d0
  aug_coefficient( 3, 0, 3, 2, 2, 2) = -0.188063194516d0
  aug_coefficient( 3, 0, 3, 2, 4, 2) = -0.044418410173d0
  aug_coefficient( 3, 0, 3, 2, 6, 2) =  0.177420363812d0
  aug_coefficient( 3, 0, 3, 3, 4, 3) = -0.203550726867d0
  aug_coefficient( 3, 0, 3, 3, 6, 3) =  0.108647340330d0
  aug_coefficient( 3, 1, 0, 0, 3,-1) = -0.282094791774d0
  aug_coefficient( 3, 1, 1,-1, 2,-2) =  0.082588898361d0
  aug_coefficient( 3, 1, 1,-1, 4,-2) = -0.238413613504d0
  aug_coefficient( 3, 1, 1, 0, 2,-1) = -0.233596680328d0
  aug_coefficient( 3, 1, 1, 0, 4,-1) = -0.238413613504d0
  aug_coefficient( 3, 1, 1, 1, 2, 0) =  0.202300659403d0
  aug_coefficient( 3, 1, 1, 1, 4, 0) = -0.150786008773d0
  aug_coefficient( 3, 1, 2,-2, 3,-3) =  0.132980760134d0
  aug_coefficient( 3, 1, 2,-2, 5,-3) = -0.200476038955d0
  aug_coefficient( 3, 1, 2,-1, 3,-2) = -0.162867503968d0
  aug_coefficient( 3, 1, 2,-1, 5,-2) = -0.245532000547d0
  aug_coefficient( 3, 1, 2, 0, 1,-1) = -0.202300659403d0
  aug_coefficient( 3, 1, 2, 0, 3,-1) = -0.126156626101d0
  aug_coefficient( 3, 1, 2, 0, 5,-1) = -0.227318461243d0
  aug_coefficient( 3, 1, 2, 1, 1, 0) =  0.233596680328d0
  aug_coefficient( 3, 1, 2, 1, 3, 0) =  0.059470803872d0
  aug_coefficient( 3, 1, 2, 1, 5, 0) = -0.169433177294d0
  aug_coefficient( 3, 1, 2, 2, 1, 1) = -0.082588898361d0
  aug_coefficient( 3, 1, 2, 2, 3, 1) =  0.206012907746d0
  aug_coefficient( 3, 1, 2, 2, 5, 1) = -0.092802373193d0
  aug_coefficient( 3, 1, 3,-3, 4,-4) =  0.166198472533d0
  aug_coefficient( 3, 1, 3,-3, 6,-4) = -0.171786528581d0
  aug_coefficient( 3, 1, 3,-2, 4,-3) = -0.095954732856d0
  aug_coefficient( 3, 1, 3,-2, 6,-3) = -0.230475813315d0
  aug_coefficient( 3, 1, 3,-1, 2,-2) = -0.206012907746d0
  aug_coefficient( 3, 1, 3,-1, 4,-2) = -0.162193101468d0
  aug_coefficient( 3, 1, 3,-1, 6,-2) = -0.242942838552d0
  aug_coefficient( 3, 1, 3, 0, 2,-1) = -0.059470803872d0
  aug_coefficient( 3, 1, 3, 0, 4,-1) = -0.099322584599d0
  aug_coefficient( 3, 1, 3, 0, 6,-1) = -0.221775454765d0
  aug_coefficient( 3, 1, 3, 1, 0, 0) =  0.282094791774d0
  aug_coefficient( 3, 1, 3, 1, 2, 0) =  0.126156626101d0
  aug_coefficient( 3, 1, 3, 1, 4, 0) =  0.025644981070d0
  aug_coefficient( 3, 1, 3, 1, 6, 0) = -0.177815950399d0
  aug_coefficient( 3, 1, 3, 2, 2, 1) =  0.162867503968d0
  aug_coefficient( 3, 1, 3, 2, 4, 1) =  0.145069920146d0
  aug_coefficient( 3, 1, 3, 2, 6, 1) = -0.121471419276d0
  aug_coefficient( 3, 1, 3, 3, 2, 2) = -0.132980760134d0
  aug_coefficient( 3, 1, 3, 3, 4, 2) =  0.188451354257d0
  aug_coefficient( 3, 1, 3, 3, 6, 2) = -0.062727571186d0
  aug_coefficient( 3, 2, 0, 0, 3,-2) =  0.282094791774d0
  aug_coefficient( 3, 2, 1,-1, 4,-3) =  0.282094791774d0
  aug_coefficient( 3, 2, 1, 0, 2,-2) =  0.184674390922d0
  aug_coefficient( 3, 2, 1, 0, 4,-2) =  0.213243618623d0
  aug_coefficient( 3, 2, 1, 1, 2,-1) = -0.261169028265d0
  aug_coefficient( 3, 2, 1, 1, 4,-1) =  0.106621809311d0
  aug_coefficient( 3, 2, 2,-2, 5,-4) =  0.268966830577d0
  aug_coefficient( 3, 2, 2,-1, 3,-3) =  0.210261043502d0
  aug_coefficient( 3, 2, 2,-1, 5,-3) =  0.253584359754d0
  aug_coefficient( 3, 2, 2, 0, 5,-2) =  0.190188269816d0
  aug_coefficient( 3, 2, 2, 1, 1,-1) = -0.261169028265d0
  aug_coefficient( 3, 2, 2, 1, 3,-1) = -0.162867503968d0
  aug_coefficient( 3, 2, 2, 1, 5,-1) =  0.117386748624d0
  aug_coefficient( 3, 2, 2, 2, 1, 0) =  0.184674390922d0
  aug_coefficient( 3, 2, 2, 2, 3, 0) = -0.188063194516d0
  aug_coefficient( 3, 2, 2, 2, 5, 0) =  0.053579475145d0
  aug_coefficient( 3, 2, 3,-3, 6,-5) =  0.254800598673d0
  aug_coefficient( 3, 2, 3,-2, 4,-4) =  0.214561305428d0
  aug_coefficient( 3, 2, 3,-2, 6,-4) =  0.266130545719d0
  aug_coefficient( 3, 2, 3,-1, 4,-3) =  0.095954732856d0
  aug_coefficient( 3, 2, 3,-1, 6,-3) =  0.230475813315d0
  aug_coefficient( 3, 2, 3, 0, 2,-2) = -0.188063194516d0
  aug_coefficient( 3, 2, 3, 0, 4,-2) = -0.044418410173d0
  aug_coefficient( 3, 2, 3, 0, 6,-2) =  0.177420363812d0
  aug_coefficient( 3, 2, 3, 1, 2,-1) = -0.162867503968d0
  aug_coefficient( 3, 2, 3, 1, 4,-1) = -0.145069920146d0
  aug_coefficient( 3, 2, 3, 1, 6,-1) =  0.121471419276d0
  aug_coefficient( 3, 2, 3, 2, 0, 0) =  0.282094791774d0
  aug_coefficient( 3, 2, 3, 2, 4, 0) = -0.179514867492d0
  aug_coefficient( 3, 2, 3, 2, 6, 0) =  0.071126380160d0
  aug_coefficient( 3, 2, 3, 3, 2, 1) =  0.210261043502d0
  aug_coefficient( 3, 2, 3, 3, 4, 1) = -0.140463346190d0
  aug_coefficient( 3, 2, 3, 3, 6, 1) =  0.031363785593d0
  aug_coefficient( 3, 3, 0, 0, 3,-3) = -0.282094791774d0
  aug_coefficient( 3, 3, 1,-1, 4,-4) = -0.325735007935d0
  aug_coefficient( 3, 3, 1, 0, 4,-3) = -0.162867503968d0
  aug_coefficient( 3, 3, 1, 1, 2,-2) =  0.319865427934d0
  aug_coefficient( 3, 3, 1, 1, 4,-2) = -0.061558130307d0
  aug_coefficient( 3, 3, 2,-2, 5,-5) = -0.347234685170d0
  aug_coefficient( 3, 3, 2,-1, 5,-4) = -0.219610497549d0
  aug_coefficient( 3, 3, 2, 0, 3,-3) =  0.210261043502d0
  aug_coefficient( 3, 3, 2, 0, 5,-3) = -0.126792179877d0
  aug_coefficient( 3, 3, 2, 1, 3,-2) =  0.210261043502d0
  aug_coefficient( 3, 3, 2, 1, 5,-2) = -0.063396089939d0
  aug_coefficient( 3, 3, 2, 2, 1,-1) = -0.319865427934d0
  aug_coefficient( 3, 3, 2, 2, 3,-1) =  0.132980760134d0
  aug_coefficient( 3, 3, 2, 2, 5,-1) = -0.023961469724d0
  aug_coefficient( 3, 3, 3,-3, 6,-6) = -0.360342462344d0
  aug_coefficient( 3, 3, 3,-2, 6,-5) = -0.254800598673d0
  aug_coefficient( 3, 3, 3,-1, 4,-4) =  0.166198472533d0
  aug_coefficient( 3, 3, 3,-1, 6,-4) = -0.171786528581d0
  aug_coefficient( 3, 3, 3, 0, 4,-3) =  0.203550726867d0
  aug_coefficient( 3, 3, 3, 0, 6,-3) = -0.108647340330d0
  aug_coefficient( 3, 3, 3, 1, 2,-2) = -0.132980760134d0
  aug_coefficient( 3, 3, 3, 1, 4,-2) =  0.188451354257d0
  aug_coefficient( 3, 3, 3, 1, 6,-2) = -0.062727571186d0
  aug_coefficient( 3, 3, 3, 2, 2,-1) = -0.210261043502d0
  aug_coefficient( 3, 3, 3, 2, 4,-1) =  0.140463346190d0
  aug_coefficient( 3, 3, 3, 2, 6,-1) = -0.031363785593d0
  aug_coefficient( 3, 3, 3, 3, 0, 0) =  0.282094791774d0
  aug_coefficient( 3, 3, 3, 3, 2, 0) = -0.210261043502d0
  aug_coefficient( 3, 3, 3, 3, 4, 0) =  0.076934943211d0
  aug_coefficient( 3, 3, 3, 3, 6, 0) = -0.011854396693d0

  call trace_exit('xrd_sf_gaunt_coefficients',status)
  return
end subroutine xrd_sf_gaunt_coefficients
subroutine xrd_sf_rho(wvfn,nkpts,nspins,occ,kpoint_weights,rho)
  !=========================================================================!
  ! Calculate the augmentation contribution to the valence charge           !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Arguments:                                                              !
  !                                                                         !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Parent module variables used:                                           !
  !   globalvar1, meaning                                                   !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Modules used:                                                           !
  !   model for type definitions                                            !
  !   cell for cell vectors and constraints                                 !
  !   parameters to define the job etc                                      !
  !   io for output and unit conversion routines                            !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Key Internal Variables:                                                 !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Necessary conditions:                                                   !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Modified from paw_recon_dia_aug_nmr_wv by Benjamin Shi, v1 15/09/2020   !
  !=========================================================================!
  use constants, only: dp,cmplx_1,cmplx_0,cmplx_i,pi
  use io,         only: io_abort,io_allocate_abort
  use basis,      only: current_basis
  use cell,       only: current_cell
  use parameters, only: current_params
  use wave, only: wavefunction,wave_page,wave_spin_type
  use ion, only: ion_gamma_recip_all,ion_gamma_phi_recip_all,max_bl_projectors,&
  num_bl_projectors,bl_projector_l,bl_projector_m,bl_dia_nmrshift
  use comms, only : comms_reduce_kp,comms_reduce_gv, comms_reduce_bnd
  use trace, only: trace_entry, trace_exit

  implicit none

  type(wavefunction),                  intent(in) :: wvfn
  integer,                             intent(in)    :: nkpts
  integer,                             intent(in)    :: nspins
  real(kind=dp),    dimension(:,:,:),       intent(in)  :: occ
  real(kind=dp),    dimension(:),           intent(in)  :: kpoint_weights
  complex(kind=dp), dimension(:,:,:,:),   intent(inout) :: rho
  real(kind=dp)                                     :: weighting,weighting_VCA

  complex(kind=dp), dimension(:,:,:,:), allocatable :: gamma_phi
  complex(kind=dp), dimension(:,:), allocatable     :: recip_gamma   ! **** colon problem
  complex(kind=dp), dimension(:,:), allocatable     :: temp_coeffs   ! **** colon problem


  !-------------------------------------------------------------------------!


  complex(kind=dp) :: factor
  real(kind=dp) :: A,B,C
  integer :: n,m,ni,nsp,ni1,nsp1,nb,rho_size,nk, i, j, k
  integer :: ierr, ns
  integer :: nk_i,ns_i
  integer                         :: status

  call trace_entry('xrd_sf_rho',status)
  if (wave_spin_type(wvfn)=="NON_COLLINEAR") then
    call io_abort('Error, xrd_sf_rho does not support NCM')
  end if

  ! Check in case we do something odd
  if(current_cell%num_species==0.or.wvfn%nbands_max==0) then
    call trace_exit('paw_recon_dia_aug_nmr',status)
    return
  end if


  allocate(gamma_phi(max_bl_projectors,wvfn%nbands_max,current_cell%max_ions_in_species,&
  current_cell%num_species),stat=ierr)
  if (ierr/=0) call io_allocate_abort('gamma_phi','xrd_sf_rho')


  !****************************
  !Page in the right components
  !****************************

  ! Zero this fresh memory




  allocate(recip_gamma(current_basis%max_plane_waves,max_bl_projectors),stat=ierr)             ! **** colon problem
  if (ierr/=0) call io_allocate_abort('recip_gamma','xrd_sf_rho')
  allocate(temp_coeffs(current_basis%max_plane_waves,wvfn%nbands_max),stat=ierr)
  if (ierr/=0) call io_allocate_abort('temp_coeffs','xrd_sf_rho')

  ! first build up gamma_phi


  rho = cmplx_0
  if (nspins .eq. 1) then
    do nk=1,nkpts
      ns = 1
      call zcopy(max_bl_projectors*wvfn%nbands_max*current_cell%max_ions_in_species*&
      current_cell%num_species,cmplx_0,0,gamma_phi,1)
      call zcopy(current_basis%max_plane_waves*wvfn%nbands_max,wvfn%coeffs(1,1,1,nk,ns),1,temp_coeffs(1,1),1)
      do nsp1=1,current_cell%num_species
        do ni1=1,current_cell%num_ions_in_species(nsp1)
          recip_gamma=cmplx_0
          call ion_gamma_recip_all(nsp1,ni1,nk,recip_gamma)
          call ion_gamma_phi_recip_all(num_bl_projectors(nsp1),recip_gamma,wvfn%nbands_max,&
          temp_coeffs(1,1),nk,gamma_phi(:,:,ni1,nsp1))

        end do
      end do

      do nb= 1,wvfn%nbands(nk,ns) ! * Loop over bands
        weighting=2.0_dp*occ(nb,nk,ns)*kpoint_weights(nk) ! Take spin degeneracy into account
        do nsp = 1,current_cell%num_species
          do ni = 1,current_cell%num_ions_in_species(nsp)

            weighting_VCA = weighting * current_cell%mixture_weight(ni,nsp)
            ! * Loop over ions
            do n=1,num_bl_projectors(nsp)

              do m=1,num_bl_projectors(nsp)
                rho(ni,n,m,nsp) = rho(ni,n,m,nsp) + &
                conjg(gamma_phi(n,nb,ni,nsp))*gamma_phi(m,nb,ni,nsp)*weighting_VCA
              end do ! projectors left
            end do ! projectors right

          end do ! num ion in species
        end do ! num species

      end do ! bands
    end do ! kpoints

    rho_size = current_cell%max_ions_in_species*max_bl_projectors**2*current_cell%num_species
    !call comms_reduce_gv(rho(1,1,1,1),rho_size,'SUM')
    ! JRY check next 2 calls
    call comms_reduce_bnd(rho(1,1,1,1),rho_size,'SUM')
    call comms_reduce_kp(rho(1,1,1,1),rho_size,'SUM')
  else
    do ns = 1, nspins
      do nk=1,nkpts
        call zcopy(max_bl_projectors*wvfn%nbands_max*current_cell%max_ions_in_species*&
        current_cell%num_species,cmplx_0,0,gamma_phi,1)
        call zcopy(current_basis%max_plane_waves*wvfn%nbands_max,wvfn%coeffs(1,1,1,nk,ns),1,temp_coeffs(1,1),1)
        do nsp1=1,current_cell%num_species
          do ni1=1,current_cell%num_ions_in_species(nsp1)
            recip_gamma=cmplx_0
            call ion_gamma_recip_all(nsp1,ni1,nk,recip_gamma)
            call ion_gamma_phi_recip_all(num_bl_projectors(nsp1),recip_gamma,wvfn%nbands_max,&
            temp_coeffs(1,1),nk,gamma_phi(:,:,ni1,nsp1))

          end do
        end do

        do nb= 1,wvfn%nbands(nk,ns) ! * Loop over bands
          weighting=occ(nb,nk,ns)*kpoint_weights(nk) ! Take spin degeneracy into account
          do nsp = 1,current_cell%num_species
            do ni = 1,current_cell%num_ions_in_species(nsp)

              weighting_VCA = weighting * current_cell%mixture_weight(ni,nsp)
              ! * Loop over ions
              do n=1,num_bl_projectors(nsp)
                do m=1,num_bl_projectors(nsp)
                  rho(ni,n,m,nsp) = rho(ni,n,m,nsp) + &
                  conjg(gamma_phi(n,nb,ni,nsp))*gamma_phi(m,nb,ni,nsp)*weighting_VCA
                end do ! projectors left
              end do ! projectors right

            end do ! num ion in species
          end do ! num species

        end do ! bands
      end do ! kpoints
    end do

    rho_size = current_cell%max_ions_in_species*max_bl_projectors**2*current_cell%num_species
    ! call comms_reduce_gv(rho(1,1,1,1),rho_size,'SUM')
    ! JRY check next 2 calls
    call comms_reduce_bnd(rho(1,1,1,1),rho_size,'SUM')
    call comms_reduce_kp(rho(1,1,1,1),rho_size,'SUM')
  end if

  deallocate(temp_coeffs,stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating temp_coeffs in xrd_sf_rho')

  deallocate(recip_gamma,stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating recip_gamma in xrd_sf_rho')

  deallocate(gamma_phi,stat=ierr)
  if (ierr/=0) call io_abort('Error in deallocating gamma_phi in xrd_sf_rho')

  call trace_exit('xrd_sf_rho',status)

  return
end subroutine xrd_sf_rho

function xrd_sf_apply_ylm_cmplx(theta,phi,l,m)
  !=========================================================================!
  ! This subroutine calculates the solid harmonic value at a given x,y,z    !
  ! position with angular momentum l and m.                                 !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Arguments:                                                              !
  !                                                                         !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Parent module variables used:                                           !
  !   globalvar1, meaning                                                   !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Modules used:                                                           !
  !   model for type definitions                                            !
  !   cell for cell vectors and constraints                                 !
  !   parameters to define the job etc                                      !
  !   io for output and unit conversion routines                            !
  !                                                                         !
  !-------------------------------------------------------------------------!
  ! Key Internal Variables:                                                 !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Necessary conditions:                                                   !
  ! ...                                                                     !
  !-------------------------------------------------------------------------!
  ! Written by Benjamin Shi, v1 15/09/2020                                  !
  !=========================================================================!

  use io, only : io_abort,io_allocate_abort
  use constants, only: dp, pi, cmplx_i

  integer           :: l,m
  real(kind=dp)     :: theta, phi, s, c, scale
  complex(kind=dp)  :: xrd_sf_apply_ylm_cmplx, spin_factor


  spin_factor = exp(1.0_dp*cmplx_i*real(m,dp)*phi)

  scale = 1/sqrt(pi)
  s = sin(theta)
  c = cos(theta)

  ! ** Apply and add the spherical harmonic

  select case(l)
  case(0)
    select case(m)
    case(0)
      xrd_sf_apply_ylm_cmplx = scale*sqrt(0.25_dp)
    case default
      call io_abort('Error: xrd_sf_apply_ylm_cmplx: abs(m)>l')
    end select
  case(1)
    select case(m)
    case(-1)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(3.0_dp/8.0_dp)*s)
    case(0)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(3.0_dp/4.0_dp)*c)
    case(1)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(3.0_dp/8.0_dp)*s)
    case default
      call io_abort('Error: xrd_sf_apply_ylm_cmplx: abs(m)>l')
    end select
  case(2)
    select case(m)
    case(-2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(15.0_dp/32.0_dp)*s**2)
    case(-1)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(15.0_dp/8.0_dp)*c*s)
    case(0)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(5.0_dp/16.0_dp)*(3*c**2-1.0_dp))
    case(1)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(15.0_dp/8.0_dp)*c*s)
    case(2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(15.0_dp/32.0_dp)*s**2)
    case default
      call io_abort('Error: xrd_sf_apply_ylm_cmplx: abs(m)>l')
    end select
  case(3)
    select case(m)
    case(-3)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(35.0_dp/64.0_dp)*s**3)
    case(-2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(105.0_dp/32.0_dp)*c*s**2)
    case(-1)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(21.0_dp/64.0_dp)*s*(5*c**2-1.0_dp))
    case(0)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(7.0_dp/16.0_dp)*(5*c**3 - 3*c))
    case(1)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(21.0_dp/64.0_dp)*s*(5*c**2-1.0_dp))
    case(2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(105.0_dp/32.0_dp)*c*s**2)
    case(3)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(35.0_dp/64.0_dp)*s**3)
    case default
      call io_abort('Error: xrd_sf_apply_ylm_cmplx: abs(m)>l')
    end select
  case(4)
    select case(m)
    case(-4)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(315.0_dp/512.0_dp)*s**4)
    case(-3)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(315.0_dp/64.0_dp)*c*s**3)
    case(-2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(45.0_dp/128.0_dp)*(7*c**2-1.0_dp)*s**2)
    case(-1)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(45.0_dp/64.0_dp)*(7*c**3-3*c)*s)
    case(0)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(9.0_dp/256.0_dp)*(35*c**4 - 30*c**2 + 3))
    case(1)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(45.0_dp/64.0_dp)*(7*c**3-3*c)*s)
    case(2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(45.0_dp/128.0_dp)*(7*c**2-1.0_dp)*s**2)
    case(3)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(315.0_dp/64.0_dp)*c*s**3)
    case(4)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(315.0_dp/512.0_dp)*s**4)
    case default
      call io_abort('Error: xrd_sf_apply_ylm_cmplx: abs(m)>l')
    end select
  case(5)
    select case(m)
    case(-5)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(693.0_dp/1024.0_dp)*s**5)
    case(-4)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(3465.0_dp/512.0_dp)*c*s**4)
    case(-3)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(385.0_dp/1024.0_dp)*(9*c**2 - 1.0_dp)*s**3)
    case(-2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(1155.0_dp/128.0_dp)*(3*c**3 - c)*s**2)
    case(-1)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(165.0_dp/512.0_dp)*(21*c**4 - 14*c**2 + 1.0_dp)*s)
    case(0)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(11.0_dp/256.0_dp)*(63*c**5 - 70*c**3 + 15*c))
    case(1)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(165.0_dp/512.0_dp)*(21*c**4 - 14*c**2 + 1.0_dp)*s)
    case(2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(1155.0_dp/128.0_dp)*(3*c**3 - c)*s**2)
    case(3)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(385.0_dp/1024.0_dp)*(9*c**2 - 1.0_dp)*s**3)
    case(4)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(3465.0_dp/512.0_dp)*c*s**4)
    case(5)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(693.0_dp/1024.0_dp)*s**5)
    case default
      call io_abort('Error: xrd_sf_apply_ylm_cmplx: abs(m)>l')
    end select
  case(6)
    select case(m)
    case(-6)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(3003.0_dp/4096.0_dp)*s**6)
    case(-5)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(9009.0_dp/1024.0_dp)*c*s**5)
    case(-4)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(819.0_dp/2048.0_dp)*(11*c**2 - 1.0_dp)*s**4)
    case(-3)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(1365.0_dp/1024.0_dp)*(11*c**3 - 3*c)*s**3)
    case(-2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(1365.0_dp/4096.0_dp)*(33*c**4 - 18*c**2 + 1.0_dp)*s**2)
    case(-1)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(273.0_dp/512.0_dp)*(33*c**5 - 30*c**3 + 5*c)*s)
    case(0)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(13.0_dp/1024.0_dp)*(231*c**6 - 315*c**4 + 105*c**2 - 5.0_dp))
    case(1)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(273.0_dp/512.0_dp)*(33*c**5 - 30*c**3 + 5*c)*s)
    case(2)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(1365.0_dp/4096.0_dp)*(33*c**4 - 18*c**2 + 1.0_dp)*s**2)
    case(3)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(1365.0_dp/1024.0_dp)*(11*c**3 - 3*c)*s**3)
    case(4)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(819.0_dp/2048.0_dp)*(11*c**2 - 1.0_dp)*s**4)
    case(5)
      xrd_sf_apply_ylm_cmplx = scale*(-sqrt(9009.0_dp/1024.0_dp)*c*s**5)
    case(6)
      xrd_sf_apply_ylm_cmplx = scale*(sqrt(3003.0_dp/4096.0_dp)*s**6)
    case default
      call io_abort('Error: xrd_sf_apply_ylm_cmplx: abs(m)>l')
    end select
  case default
    call io_abort('Error: xrd_sf_apply_ylm_cmplx: l>6 not programmed')
  end select

  xrd_sf_apply_ylm_cmplx = xrd_sf_apply_ylm_cmplx*spin_factor
  return
end function xrd_sf_apply_ylm_cmplx

end module xrd_sf
