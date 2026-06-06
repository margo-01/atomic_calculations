!-*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: git -*-
!=============================================================================!
!                           O T F G _ G E N E R A T E                         !
!=============================================================================!
!                                                                             !
!-----------------------------------------------------------------------------!
!               Module author: C. J. Pickard (c), 2019, Cambridge             !
!=============================================================================!

module otfg_generate

  use constants, only : dp,pi

  use otfg_basis
  use otfg_data
  use otfg_atom
  use otfg_pspot
  use otfg_scf
  use otfg_pseudise
  use otfg_io
  use otfg_algor
  use otfg_solvers

  use xc_kernels, only : xc_is_mgga ! ** CJP This is not ideal ... fix

  implicit none

  private

  !---------------------------------------------------------------------------!
  ! **                    P u b l i c   R o u t i n e s                    ** !
  !---------------------------------------------------------------------------!

  public :: otfg_inquire_psp
  public :: otfg_generate_psp

  !---------------------------------------------------------------------------!
  ! **                       P r i v a t e  f l a g s                      ** !
  !---------------------------------------------------------------------------!

  integer, save :: id_counter=0 ! * Keep track of the potential via an id number

contains

  subroutine otfg_inquire_psp(nsp,tot_num_projectors,tot_num_points,ionic_charge,num_mesh,num_kkbeta,num_nqf,&
       bl_points,tot_num_bl_projectors,num_core_projectors,spin_orbit_coupling)

    use constants, only : file_maxpath,periodic_table_name
    use io,        only : io_abort,io_allocate_abort,io_failure,io_open_file
    use cell,      only : current_cell
    use bib,       only : bib_add

    !=========================================================================!
    ! This routine finds out some important parameters for the on-the-fly     !
    ! pseudopotentials                                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) nsp                   : input  : the current species                 !
    ! 2) tot_num_projectors    : output : the total number of projectors      !
    ! 3) tot_num_points        : output : total number of radial recip points !
    ! 4) ionic_charge          : output : the ionic charge of this species    !
    ! 5) num_mesh              : output : number of mesh points               !
    ! 6) num_kkbeta            : output : number of mesh points, up to rc     !
    ! 7) num_nqf               : output : number of optimisation functions    !
    ! 8) bl_points             : output : total number of blochl recip points !
    ! 9) tot_bl_num_projectors : output : total number of blochl projectors   !
    !10) num_core_projectors   : output : number of core projectors           !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) constants                                                            !
    ! 2) io                                                                   !
    ! 3) cell                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 21st May 2002                                   !
    !=========================================================================!

    use comms, only: on_root_node, comms_gcopy
    use trace, only: trace_entry, trace_exit

    implicit none

    integer,       intent(in)  :: nsp
    integer,       intent(out) :: tot_num_projectors
    integer,       intent(out) :: tot_num_points
    real(kind=dp), intent(out) :: ionic_charge
    integer,       intent(out) :: num_mesh
    integer,       intent(out) :: num_kkbeta
    integer,       intent(out) :: num_nqf
    integer,       intent(out) :: bl_points
    integer,       intent(out) :: tot_num_bl_projectors
    integer,       intent(out) :: num_core_projectors     ! for ELNES
    logical,       intent(in)  :: spin_orbit_coupling

    !-------------------------------------------------------------------------!

    type(ae_basis) :: ab

    integer                     :: Z,i,nref,n,indx,nb,l,nl,ierr,ng,nn,nnmax  !add nn,nnmax for ELNES
    character(len=file_maxpath) :: definition,temp_def,pspfile
    character(len=80)           :: ctemp
    character(len=8)            :: cname
    character(len=7)            :: cfg
    character(len=1)            :: lchar

    real(kind=dp)      :: occ

    integer, parameter                   :: nflags = 5
    character(len=80), dimension(nflags) :: flags
    integer                              :: unit_otfglib

    ! ** Find the atomic number of this species

    integer                         :: status

    call trace_entry('otfg_inquire_psp',status)

    Z = -1 ! Set atomic number to -1 to check later
    do i=1,size(periodic_table_name)
       ! Get the atomic_number for the symbol
       indx  = index(current_cell%species_symbol(nsp),':')
       if(indx>0) then
          cname = current_cell%species_symbol(nsp)(:indx-1)
       else
          cname = current_cell%species_symbol(nsp)
       end if
       if(cname.eq.periodic_table_name(i)) Z=i
    end do
    if(Z<1) call io_abort('Error: otfg_inquire_psp: current_cell%species_symbol not known')

    ! ** Initialise the atomic basis

    call otfg_init_ae_basis(ab,Z)

    ! ** Determine the number of projectors

    ! * Set the working pseudopotential definition

    ps_set(nsp)=999

    if(trim(current_cell%species_pot(nsp)).ne."NULL.NULL") then
       select case(trim(current_cell%species_pot(nsp)))
       case("QC5")
          ps_set(nsp)=1
          current_cell%species_pot(nsp)="NULL.NULL"
       case("HARD")
          ps_set(nsp)=2
          current_cell%species_pot(nsp)="NULL.NULL"
       case("MS")
          ps_set(nsp)=3
          current_cell%species_pot(nsp)="NULL.NULL"
       case("C7")
          ps_set(nsp)=7
          current_cell%species_pot(nsp)="NULL.NULL"
       case("C8")
          ps_set(nsp)=8
          current_cell%species_pot(nsp)="NULL.NULL"
       case("C9","C16")
          ps_set(nsp)=9
          current_cell%species_pot(nsp)="NULL.NULL"
          call bib_add('OTFG_C9')
       case("NCP9")
          ps_set(nsp)=10
          current_cell%species_pot(nsp)="NULL.NULL"
       case("C17")
          ps_set(nsp)=11
          current_cell%species_pot(nsp)="NULL.NULL"
       case("NCP17")
          ps_set(nsp)=12
          current_cell%species_pot(nsp)="NULL.NULL"
       case("C18")
          ps_set(nsp)=13
          current_cell%species_pot(nsp)="NULL.NULL"
       case("NCP18")
          ps_set(nsp)=14
          current_cell%species_pot(nsp)="NULL.NULL"
       case("C19")
          ps_set(nsp)=15
          current_cell%species_pot(nsp)="NULL.NULL"
       case("NCP19","NCP")
          ps_set(nsp)=16
          current_cell%species_pot(nsp)="NULL.NULL"
       case("SOC19","SOC")
          ps_set(nsp)=17
          current_cell%species_pot(nsp)="NULL.NULL"
       end select
    else
       ! * By default, select the current recommendation
       ps_set(nsp)=15
    end if

    ! ** Read the definition string from a supplied library file (.otfglib)

    if(index(current_cell%species_pot(nsp),".otfglib").gt.0) then

       if( on_root_node ) then
          ! ** Open the OTFG library file

          flags(1) = 'READ'
          flags(2) = 'FORMATTED'
          flags(3) = 'OLD'
          flags(4) = 'REWIND'
          flags(5) = 'SEQUENTIAL'

          cname=current_cell%species_symbol(nsp)

          pspfile=current_cell%species_pot(nsp)
          if(index(pspfile,'@').gt.0) then
             cname=trim(pspfile(:index(pspfile,"@")-1))
             pspfile=trim(pspfile(index(pspfile,"@")+1:))
          end if

          call io_open_file(pspfile,nflags,flags,unit_otfglib,status)

          if(status==io_failure) &
               call io_abort('Error: otfg_define_psp: file "'//trim(current_cell%species_pot(nsp))//'" not found')

          do
             read(unit_otfglib,'(a)',end=111,err=103) ctemp
             if(trim(cname).eq.trim(ctemp(:index(ctemp," ")))) then
                current_cell%species_pot(nsp)=ctemp(index(ctemp," "):)
                goto 112
             end if
          end do
111       continue

          call io_abort('Error: otfg_inquire_psp: OTFG string for "'//trim(current_cell%species_symbol(nsp))// &
               & '" not found in '//trim(current_cell%species_pot(nsp)))

112       close(unit_otfglib)

       end if

       call comms_gcopy(current_cell%species_pot(nsp), len(current_cell%species_pot(nsp)))

    end if

    select case(trim(current_cell%species_pot(nsp)))
    case("NULL.NULL")
       definition = adjustl(ps_def(Z,ps_set(nsp)))
    case default
       definition = adjustl(current_cell%species_pot(nsp))
    end select

    if(len_trim(definition).eq.0)  call io_abort('Error: otfg_inquire_psp: no OTFG definition found')

    if(count((/(definition(i:i),i=1,len_trim(definition))/)=="|").lt.5) &
         call io_abort('Error: otfg_inquire_psp: OTFG definition or library not understood: '//trim(definition))

    ! ** Record the OTFG definition used

    current_cell%species_pot(nsp)=definition

    ! ** Set the configuration to the groundstate default

    cfg_occ = real(gs_occ(Z,:),dp)

    ! ** Check if there any test configuration flags are set

    indx = index(definition,'[')
    if(indx/=0) then

       ! * Remove the flag string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,']')
       if(indx==0) call io_abort('Error: otfg_inquire_psp: incorrect syntax for test config adjust')
       write (definition,'(a)') trim(adjustl(definition))//trim(adjustl(ctemp(indx+1:)))
    end if

    ! ** Check if there any pseudisation flags are set

    indx = index(definition,'(')
    if(indx/=0) then

       ! * Remove the flag string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,')')
       if(indx==0) call io_abort('Error: otfg_inquire_psp: incorrect syntax for pseudisation flags')
       write (definition,'(a)') trim(adjustl(definition))//trim(adjustl(ctemp(indx+1:)))
    end if

    ! ** Check if there any configuration adjustment required

    indx = index(definition,'{')
    if(indx/=0) then

       ! * Extract the configuration adjust string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,'}')
       if(indx==0) call io_abort('Error: otfg_inquire_psp: incorrect syntax for config adjust')
       ctemp = ctemp(:indx-1)

       ! * Adjust the default groundstate configuration

       do while(len_trim(ctemp)>0)
          indx = index(ctemp,',')
          if(indx==0) then
             cfg   = ctemp
             ctemp = ''
          else
             cfg   = ctemp(:indx-1)
             ctemp = ctemp(indx+1:)
          end if
          if(len_trim(cfg)<3) call io_abort('Error: otfg_inquire_psp: incorrect syntax for config adjust')

          ! * Set the principal quantum number

          read(cfg(1:1),*) nl
          nl = 10*nl

          ! * Set the angular momentum

          read(cfg(2:2),*) lchar

          select case(lchar)
          case('s') ; nl=nl+0
          case('p') ; nl=nl+1
          case('d') ; nl=nl+2
          case('f') ; nl=nl+3
          case default
             call io_abort('Error: otfg_inquire_psp: angular momentum not recognised')
          end select

          ! * Adjust the configuration

          if(count(gs_occ(0,:)==nl)/=1) call io_abort('Error: otfg_inquire_psp: nl pair not recognised or unique')

          do i=1,20
             if(gs_occ(0,i)==nl) read(cfg(3:),*) cfg_occ(i)
          end do

       end do

    end if

    ! ** Now proceed with the analysis of the string

    indx = index(definition,'|',back=.true.)
    definition = definition(indx+1:)

    if(len_trim(definition)<1) call io_abort('Error: otfg_inquire_psp: no reference states specified')

    nref=1
    do i=1,len_trim(definition)
       if(definition(i:i)==':') nref=nref+1
    end do

    ! ** Find out how many beta functions there are -- so the arrays can be allocated

    temp_def              = definition
    tot_num_projectors    = 0
    tot_num_bl_projectors = 0
    occ                   = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if

       ! * Find out the angular momentum of this reference state

       read(ctemp(2:2),*) l

       ! * Find out the reference state label

       read(ctemp(1:2),*) nl

       do i=1,20
          if(nl==gs_occ(0,i)) occ = occ+cfg_occ(i)
       end do

       temp_def = temp_def(indx+1:)
       nb = 0
       ng = 0

       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='N') nb=nb+2*l+1 ! * Count the norm conserving projectors
          if(ctemp(i:i)=='U') nb=nb+2*l+1 ! * Count the ultrasoft projectors
          if(ctemp(i:i)=='H') ng=ng+2*l+1 ! * Count the norm conserving gamma projectors
          if(ctemp(i:i)=='G') ng=ng+2*l+1 ! * Count the ultrasoft gamma projectors
       end do
       if((nb==0).and.( (index(ctemp,'P')==0) .and. index(ctemp,'L')==0)) nb=2*(2*l+1)
       if(spin_orbit_coupling.and. l>0) then ! jry double for j
          tot_num_projectors    = tot_num_projectors    + 2*nb
          tot_num_bl_projectors = tot_num_bl_projectors + 2*ng
       else
          tot_num_projectors    = tot_num_projectors    + nb
          tot_num_bl_projectors = tot_num_bl_projectors + ng
       endif
    end do

    ! * By default same number of gamma as beta projectors
    !   - In other words we only specify _extra_ gamma projectors
    !   I can't see any physical situation in which we would why we would need more betas than gamma
    !   and it will makes life so much easier if I don't code for it.

    tot_num_bl_projectors = tot_num_bl_projectors + tot_num_projectors


    ! ** For ELNES  calculation
    ! ** Find out how many core projectors, there are maybe a lot of ways to define which is core orbital
    ! depending on what is your interested initial state

    num_core_projectors = 1  ! We assume there is at least one core projectors. Even for H.

    temp_def = definition

    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if

       ! * Find out the principal quantum number label

       read(ctemp(1:1),*) nn

       if (n==1) then
          nnmax=nn
       else
          nnmax=max(nn,nnmax)
       end if
       temp_def = temp_def(indx+1:)
    end do

    ! We consider configuration of core orbitals is a noble config. .i.e. it is 1s, or 1s,2s,2p,
    ! or 1s,2s,2p,3s,3p, or etc. ---ELNES
!!!!!!! TEMPORARY FIX :::
    nnmax=min(6,nnmax)

    select case (nnmax)
    case(1,2)
       num_core_projectors = 1  ! 1s
    case(3)
       num_core_projectors = 5  ! 1s,2s,2p
    case(4)
       num_core_projectors = 9  ! 1s,2s,2p,3s,3p
    case(5)
       num_core_projectors = 18 ! 1s,2s,2p,3s,3p,3d,4s,4p
    case(6)
       num_core_projectors = 27 ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p
    case(7)
       num_core_projectors = 43 ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p
    case default
       call io_abort('Error in setting core_projectors: there is no projectors? ')
    end select

    ! ** Set the other parameters for this species

    tot_num_points     = 2001

    ionic_charge       = real(Z,dp)-sum(cfg_occ)+occ

    num_mesh           = ab%npts

    num_kkbeta         = ab%npts ! ** SORT THIS OUT ***

    num_nqf            = 10 ! ** CHECK THIS OUT ***

    bl_points          = tot_num_points ! ** same as beta proj (so one could argue that I don't need this variable,
    !    but it might come in useful later) (CJP Did it?)

    ! * Deallocate the basis arrays

    if (allocated(ab%r)) then
       deallocate(ab%r,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ab%r in otfg_inquire_psp')
    end if

    if (allocated(ab%rwgt)) then
       deallocate(ab%rwgt,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ab%rwgt in otfg_inquire_psp')
    end if

    call trace_exit('otfg_inquire_psp',status)

    return

103 call io_abort('Error: otfg_inquire_psp: problem parsing definition')

  end subroutine otfg_inquire_psp

  subroutine otfg_generate_psp(nsp,ps)

    use constants,  only : file_maxpath,periodic_table_name
    use io,         only : io_abort,stdout,io_unit_to_atomic,io_allocate_abort
    use parameters, only : current_params
    use cell,       only : current_cell
    use comms,      only : my_node_id,root_node_id

    !=========================================================================!
    ! Generate the pseudopotential                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) nsp : input  : the species for which the potential is required       !
    ! 2) ps  : output : the pseudopotential                                   !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) pspot                                                                !
    ! 2) ae_basis                                                             !
    ! 3) ae_atom                                                              !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) constants                                                            !
    ! 2) io                                                                   !
    ! 3) cell                                                                 !
    ! 4) comms                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 25th June 2002                                  !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    integer,     intent(in)    :: nsp
    type(pspot), intent(inout) :: ps

    !-------------------------------------------------------------------------!

    type(ae_basis)     :: ab
    type(ae_atom)      :: aeat
    type(pseudo_basis) :: pb
    type(pseudo_atom)  :: psat

    integer        :: Z,i,j,l,nc ! * The atomic number, and counters
    integer        :: jhnps(0:4),jhl,jhnn,ierr

    integer                     :: indx,jndx,nl
    character(len=file_maxpath) :: definition
    character(len=80)           :: ctemp
    character(len=30)           :: flags
    character(len=8)            :: cname
    character(len=7)            :: cfg
    character(len=1)            :: lchar

    real(kind=dp) :: max_err,wgt

    real(kind=dp),    dimension(100,2) :: evals
    integer,          dimension(100,2) :: nlabels


    logical, allocatable, dimension(:,:) :: ps_state

    logical :: on_root,found_ghost

    integer :: status

    call trace_entry('otfg_generate_psp',status)

    on_root = my_node_id==root_node_id

    ! ** Increment id counter

    id_counter = id_counter + 1

    ! ** Store the id for this potential

    ps%id = id_counter

    ! ** Find the atomic number of this species

    Z = -1 ! Set atomic number to -1 to check later
    do i=1,size(periodic_table_name)
       ! Get the atomic_number for the symbol
       indx  = index(current_cell%species_symbol(nsp),':')
       if(indx>0) then
          cname = current_cell%species_symbol(nsp)(:indx-1)
       else
          cname = current_cell%species_symbol(nsp)
       end if
       if(cname.eq.periodic_table_name(i)) Z=i
    end do
    if(Z<1) call io_abort('Error: otfg_generate_psp: current_cell%species_symbol not known')

    ! ** Adjust the configuration

    ! * Set the working pseudopotential definition

    select case(trim(current_cell%species_pot(nsp)))
    case("NULL.NULL")
       definition = adjustl(ps_def(Z,ps_set(nsp)))
    case default
       definition = adjustl(current_cell%species_pot(nsp))
    end select

    ! ** Set the configuration to the groundstate default

    cfg_occ = real(gs_occ(Z,:),dp)

    ! ** Check if there any configuration adjustment required

    indx = index(definition,'{')
    if(indx/=0) then

       ! * Extract the configuration adjust string

       ctemp = definition(indx+1:)
       indx = index(ctemp,'}')
       if(indx==0) call io_abort('Error: otfg_generate_psp: incorrect syntax for config adjust')
       ctemp = ctemp(:indx-1)

       ! * Adjust the default groundstate configuration

       do while(len_trim(ctemp)>0)
          indx = index(ctemp,',')
          if(indx==0) then
             cfg   = ctemp
             ctemp = ''
          else
             cfg   = ctemp(:indx-1)
             ctemp = ctemp(indx+1:)
          end if
          if(len_trim(cfg)<3) call io_abort('Error: otfg_generate_psp: incorrect syntax for config adjust')

          ! * Set the principal quantum number

          read(cfg(1:1),*) nl
          nl = 10*nl

          ! * Set the angular momentum

          read(cfg(2:2),*) lchar

          select case(lchar)
          case('s') ; nl=nl+0
          case('p') ; nl=nl+1
          case('d') ; nl=nl+2
          case('f') ; nl=nl+3
          case default
             call io_abort('Error: otfg_generate_psp: angular momentum not recognised')
          end select

          ! * Adjust the configuration

          if(count(gs_occ(0,:)==nl)/=1) call io_abort('Error: otfg_generate_psp: nl pair not recognised or unique')

          do i=1,20
             if(gs_occ(0,i)==nl) read(cfg(3:),*) cfg_occ(i)
          end do

       end do

    end if

    ! ** Store the generation configuration

    ps%gen_cfg = cfg_occ

    cfg_occ_save=cfg_occ   !  ELNES
    flags=''
    indx = index(definition,'(')
    if(indx/=0) then

       ! * Remove the flag string
       ctemp = definition(indx+1:)
       indx = index(ctemp,')')
       flags = ctemp(:indx-1)
       if(indx==0) call io_abort('Error: otfg_generate_psp: incorrect syntax for pseudisation flags')
    end if
    ! ** Compatability mode?

    compat7=.false.
    if(index(flags,"compat7")>0)  compat7=.true.

    ! * Take the values from the relativistic_treatment param, but allow string over-rides

    aeat%atomic_solver=ps%atomic_solver
    if(index(flags,"schr")>0)  aeat%atomic_solver="sh"
    if(index(flags,"dirac")>0) aeat%atomic_solver="dr"
    if(index(flags,"zora")>0)  aeat%atomic_solver="zo"
    if(index(flags,"szora")>0) aeat%atomic_solver="zs"
    if(index(flags,"kh")>0) aeat%atomic_solver="kh"
    ps%atomic_solver=aeat%atomic_solver

    ! Set spin-orbit flag for the all-electron atom
    ! could we have a case in which we want it for the all-electron atom - but then remove in the pspot?
    aeat%spin_orbit_coupling=ps%SOC


    ! ** Initialise the atomic basis

    call otfg_init_ae_basis(ab,Z)

    ! ** Initialise the atom

    call otfg_init_ae_atom(Z,ab,aeat)

    ! ** Set the all electron XC used

    aeat%theory = ps%theory

    if(index(flags,"RPBE")>0) then
       aeat%theory = "RPBE"
    else if(index(flags,"PBESOL")>0) then
       aeat%theory = "PBESOL"
    else if(index(flags,"PBE")>0) then
       aeat%theory = "PBE"
    else if(index(flags,"WC")>0) then
       aeat%theory = "WC"
    else if(index(flags,"BLYP")>0) then
       aeat%theory = "BLYP"
    else if(index(flags,"LDA_PW")>0) then
       aeat%theory = "LDA_PW"
    else if(index(flags,"LDA")>0) then
       aeat%theory = "LDA"
    else if(index(flags,"MS2")>0) then
       aeat%theory = "MS2"
    else if(index(flags,"RSCAN")>0) then
       aeat%theory = "RSCAN"
    end if


    if( xc_is_mgga(aeat%theory) ) then
       if( .not. ( aeat%atomic_solver == "kh" .or. aeat%atomic_solver == "sh" ) ) then
          call io_abort('Error: otfg_generate_psp: Only the Koelling-Harmon or Schroedinger solver &
               & may be used with meta-GGA functionals.')
       endif
    endif

    ! ** Set the finite nucleus

    ! JRY: not sure about the following line for dr
    if(aeat%atomic_solver == "kh" .or. aeat%atomic_solver == "sh")then! .or. aeat%atomic_solver == "dr") then
       aeat%fin_nuc = .true.
    else
       aeat%fin_nuc = .false.
    end if

    ! ** Do the all electron calculation for the reference atom

    call otfg_ae_scf(nsp,ab,aeat)

    !JRY: I commented out the next lines - which might break zora
    ! ** If Dirac, copy the saved large components into the orbitals
    !    if(aeat%atomic_solver == "dr") then
    !       aeat%orbr(:,:,:) = aeat%orbr_hack(:,:,:)
    !    end if

    ! ** Make the ZORA K function
    call otfg_update_zora_K(ab,aeat)

    ! ** If we're doing a ZORA atom, we've now got a Dirac AE potential to use for K
    if(aeat%atomic_solver == "zo") then

       ! Now do scalar-relativistic ZORA atom with Dirac-derived K
       call otfg_ae_scf(nsp,ab,aeat)

       ! ** If scaled ZORA, solve and rescale the solutions
    else if(aeat%atomic_solver == "zs") then

       ! Now do scalar-relativistic ZORA atom with Dirac-derived K
       call otfg_ae_scf(nsp,ab,aeat)

       call otfg_zora_scale(ab,aeat)

    end if

    ! Next lines useful for ZORA / Dirac debugging. Maybe put in on a flag? JRY
    !!    ! Dump out orbitals and potential
    !    call otfg_write_ae_potential(ab,aeat,nsp)
    !    call otfg_write_ae_orbitals(ab,aeat,nsp)

    ! ** Define the pseudopotential to be generated

    if(ps%SOC) then
       call otfg_define_psp_dir(nsp,ab,aeat,ps)
    else
       call otfg_define_psp(nsp,ab,aeat,ps)
    endif

    ! ** Construct the pseudopotential

    if(ps%SOC) then
       call otfg_construct_psp_dir(ab,aeat,ps)
    else
       call otfg_construct_psp(ab,aeat,ps)
    endif

    ! ** Descreen the pseudopotential

    call otfg_descreen_psp(ab,ps,aeat,nsp)

    ! ** Write out generated potential

    call otfg_output_psp(nsp,ab,ps,aeat)

    ! ** Convert the pseudopotential to real space form

    !call otfg_realspace_psp(ab,ps)

    ! ** Destroy the ae_basis and ae_atom

    call otfg_basis_ae_dealloc(ab,aeat)

    ! ** Test the potential in a different configuration
    !    Check if there any configuration supplied

    if(testing) then
       indx = index(definition,'[')

       ! ** By default, use the generation configuration

       cfg_occ = ps%gen_cfg

       ! * Extract the configuration adjust string
       if(indx/=0) then

          ctemp = definition(indx+1:)
          indx = index(ctemp,']')
          if(indx==0) call io_abort('Error: otfg_generate_psp: incorrect syntax for config adjust')
          ctemp = ctemp(:indx-1)

          ! * Adjust the default groundstate configuration

          do while(len_trim(ctemp)>0)
             indx = index(ctemp,',')
             if(indx==0) then
                cfg   = ctemp
                ctemp = ''
             else
                cfg   = ctemp(:indx-1)
                ctemp = ctemp(indx+1:)
             end if
             if(len_trim(cfg)<3) call io_abort('Error: otfg_generate_psp: incorrect syntax for config adjust')

             ! * Set the principal quantum number

             read(cfg(1:1),*) nl
             nl = 10*nl

             ! * Set the angular momentum

             read(cfg(2:2),*) lchar

             select case(lchar)
             case('s') ; nl=nl+0
             case('p') ; nl=nl+1
             case('d') ; nl=nl+2
             case('f') ; nl=nl+3
             case default
                call io_abort('Error: otfg_generate_psp: angular momentum not recognised')
             end select

             ! * Adjust the configuration

             if(count(gs_occ(0,:)==nl)/=1) call io_abort('Error: otfg_generate_psp: nl pair not recognised or unique')

             do i=1,20
                if(gs_occ(0,i)==nl) read(cfg(3:),*) cfg_occ(i)
             end do

          end do
       endif

       if(on_root) then
          write (stdout,*)
          write (stdout,*) "Doing a self consistent calculation for the test configuration"
       end if

       ! ** Initialise the atomic basis

       call otfg_init_ae_basis(ab,Z)  !JRY

       ! ** Initialise the atom

       call otfg_init_ae_atom(Z,ab,aeat)

       ! ** Do not fix BCs when testing in the free atom

       bc(nsp)%no(:) = 0

       ! ** Do the all electron calculation for the reference atom

       call otfg_ae_scf(nsp,ab,aeat)

       !I've turned off this which might have broken ZORA
       ! ** If dirac, copy the saved large components into the orbitals
!!$       if(aeat%atomic_solver == "dr") then
!!$          aeat%orbr(:,:,:) = aeat%orbr_hack(:,:,:)
!!$       end if

       ! ** Make the ZORA K function
       call otfg_update_zora_K(ab,aeat)

       if(aeat%atomic_solver == "zo") then

          ! Make the ZORA K function
          ! call otfg_update_zora_K(ab,aeat)

          ! Now do scalar-relativistic ZORA atom with Dirac-derived K
          call otfg_ae_scf(nsp,ab,aeat)

          ! ** If scaled ZORA, solve and rescale the solutions
       else if(aeat%atomic_solver == "zs") then

          ! Make the ZORA K function
          ! call otfg_update_zora_K(ab,aeat)

          ! Now do scalar-relativistic ZORA atom with Dirac-derived K
          call otfg_ae_scf(nsp,ab,aeat)

          call otfg_zora_scale(ab,aeat)

       end if

       if(on_root) then

          write (stdout,*) "---------------------------------------"

          indx=0
          jhnps=0
          do i=1,20
             if(cfg_occ(i)>0) then
                jhl=mod(gs_occ(0,i),10)
                jhnps(jhl)=jhnps(jhl)+1
                jhnn=(gs_occ(0,i)-mod(gs_occ(0,i),10))/10
                if(aeat%spin_orbit_coupling) then
                   write (stdout,'(1x,a,3i1,a,f18.8)') "AE eigenvalue nl2j ",jhnn,jhl,2*jhl+1," =",aeat%eval(jhnps(jhl),jhl)
                   indx=indx+1
                   evals(indx,1)=aeat%eval(jhnps(jhl),jhl)
                   nlabels(indx,1)=10*jhnn+jhl+100*(2*jhl+1)
                   if(jhl>0) then
                      write (stdout,'(1x,a,3i1,a,f18.8)') "AE eigenvalue nl2j ",jhnn,jhl,2*jhl-1," =",aeat%eval(jhnps(jhl)+1,jhl)
                      indx=indx+1
                      evals(indx,1)=aeat%eval(jhnps(jhl)+1,jhl)
                      nlabels(indx,1)=10*jhnn+jhl+100*(2*jhl-1)
                      jhnps(jhl)=jhnps(jhl)+1
                   endif
                else
                   write (stdout,'(1x,a,2i1,a,f18.8)') "AE eigenvalue nl ",jhnn,jhl," =",aeat%eval(jhnps(jhl),jhl)
                   indx=indx+1
                   evals(indx,1)=aeat%eval(jhnps(jhl),jhl)
                   nlabels(indx,1)=10*jhnn+jhl
                endif
             end if
          end do

          write (stdout,*) "---------------------------------------"

       end if

       ! ** Do the pseudo calculation for the reference atom

       call otfg_pseudo_scf_test(ps,pb,psat)

       if(on_root) then

          allocate(ps_state(psat%norbmax,0:pb%lmax))
          ps_state=.false.

          jndx=0
          jhnps=0
          do i=1,20
             if(cfg_occ(i)>0) then
                jhl=mod(gs_occ(0,i),10)
                jhnps(jhl)=jhnps(jhl)+1
                jhnn=(gs_occ(0,i)-mod(gs_occ(0,i),10))/10
                if(aeat%spin_orbit_coupling) then
                   write (stdout,'(1x,a,3i1,a,f18.8)') "PS eigenvalue nl2j ",jhnn,jhl,2*jhl+1," =",psat%eval(jhnps(jhl),jhl)
                   ps_state(jhnps(jhl),jhl)=.true.
                   jndx=jndx+1
                   evals(jndx,2)=psat%eval(jhnps(jhl),jhl)
                   nlabels(jndx,2)=10*jhnn+jhl+100*(2*jhl+1)
                   if(jhl>0) then
                      write (stdout,'(1x,a,3i1,a,f18.8)') "PS eigenvalue nl2j ",jhnn,jhl,2*jhl-1," =",psat%eval(jhnps(jhl)+1,jhl)
                      ps_state(jhnps(jhl)+1,jhl)=.true.
                      jndx=jndx+1
                      evals(jndx,2)=psat%eval(jhnps(jhl)+1,jhl)
                      nlabels(jndx,2)=10*jhnn+jhl+100*(2*jhl-1)
                   endif
                else
                   write (stdout,'(1x,a,2i1,a,f18.8)') "PS eigenvalue nl ",jhnn,jhl," =",psat%eval(jhnps(jhl),jhl)
                   ps_state(jhnps(jhl),jhl)=.true.
                   jndx=jndx+1
                   evals(jndx,2)=psat%eval(jhnps(jhl),jhl)
                   nlabels(jndx,2)=10*jhnn+jhl
                endif
             end if
          end do

          write (stdout,*) "---------------------------------------"

          max_err=0.0_dp
          do j=1,jndx
             do i=1,indx
                if(nlabels(i,1).eq.nlabels(j,2)) then
                   if(abs(evals(i,1)-evals(j,2)).gt.max_err) max_err=abs(evals(i,1)-evals(j,2))
                   exit
                end if
             end do
          end do

          if(max_err.lt.1e-2_dp) then
             write (stdout,'(a,ES10.2)') ' Maximum eigenvalue error: ',max_err
          else
             write (stdout,'(a,ES10.2,a)') ' Maximum eigenvalue error: ',max_err,' <- Warning'
          end if
          write (stdout,*)
          write (stdout,'(a,i4,a3,i4,a3,i4,a3,i4)') ' Derived cutoff energies: C=',&
               int(io_unit_to_atomic(ps%calc_c,current_params%energy_unit))+1,'M=',&
               int(io_unit_to_atomic(ps%calc_m,current_params%energy_unit))+1,'F=',&
               int(io_unit_to_atomic(ps%calc_f,current_params%energy_unit))+1,'E=',&
               int(io_unit_to_atomic(ps%calc_e,current_params%energy_unit))+1
          write (stdout,*)


          nc=nint(ps%rcmax/pb%dr)
          nc=(nc/2)*2+1

          write (stdout,'(a,f10.5,a)') ' Checking for ghost states up to',maxval(psat%eval), ' Ha'
          write (stdout,*)

          found_ghost=.false.
          do l=0,pb%lmax
             if(.not.any(psat%occ(:,l).gt.0.0_dp)) cycle
             do j=1,psat%norbs(l)

                wgt=otfg_regin(nc,pb%dr,psat%orbr(1:nc,j,l)**2*pb%r(1:nc)**2)

                if((wgt.gt.0.5_dp).and.(.not.ps_state(j,l))) then
                   write (stdout,'(a,2i5,2f10.5)') ' Ghost: ',l,j,wgt,psat%eval(j,l)
                   found_ghost=.true.
                end if

             end do
          end do

          if(found_ghost) write (stdout,*)

          if(found_ghost) then
             write (stdout,'(a)') ' Identified Ghost States'
          else
             write (stdout,'(a)') ' No Ghost States Identified'
          end if

          write (stdout,*)
          write (stdout,*) "***************************************"

          deallocate(ps_state)

       end if

       ! ** Destroy the pseudo_basis and pseudo_atom

       call otfg_basis_pseudo_dealloc(pb,psat)

       ! ** Destroy the ae_basis and ae_atom

       call otfg_basis_ae_dealloc(ab,aeat)

    end if


    ! ** Set the configuration back to the groundstate default

    cfg_occ = real(gs_occ(Z,:),dp)

    ! * The orbitals to which to apply the BCs during future regenerations

    do l=0,bc(nsp)%lmax
       if(ps%scpsp.and.(count(ps%pseudise(:,l))>1)) &
            call io_abort('Error: otfg_generate_psp: cannot generate SCPSPs with shallow cores')
       do i=1,10
          if(ps%pseudise(i,l)) then
             bc(nsp)%no(l) = i
          end if
       end do
    end do

    call trace_exit('otfg_generate_psp',status)

    return
  end subroutine otfg_generate_psp

  !************************************************************************
  !**          S T A R T   O F   P R I V A T E   R O U T I N E S         **
  !************************************************************************


  subroutine otfg_define_psp(nsp,ab,aeat,ps)

    use io,         only : io_abort,io_allocate_abort,io_failure,io_open_file,io_atomic_to_unit,&
         io_code_present
    use parameters, only : current_params
    use cell,       only : current_cell
    use constants,  only : file_maxpath,mu_0_si,elementary_charge_si,electron_mass_si,speed_light_si,&
         fine_structure_si,planck_si,avogadro_si

    !=========================================================================!
    ! Define the pseudopotential                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) nsp     : input  : The species                                       !
    ! 2) ab      : input  : The atomic basis                                  !
    ! 3) aeat    : input  : The initialised atom                              !
    ! 4) ps      : output : The pseudopotential                               !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) ae_atom                                                              !
    ! 3) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) parameters                                                           !
    ! 3) cell                                                                 !
    ! 4) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 1st July 2002                                   !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    integer,        intent(in)    :: nsp
    type(ae_basis), intent(in)    :: ab
    type(ae_atom),  intent(in)    :: aeat
    type(pspot),    intent(inout) :: ps

    !-------------------------------------------------------------------------!
!BEGIN: TEMSTEP edits
    integer       :: npts,n1,n2,n3,n4,n5,n6,n7,n8 !deleted i here, shouldn't cause issues
    integer       :: num_core_electrons=1
    real(kind=dp) :: ionic_charge_dummy
    real, dimension(:), allocatable :: rcc ! Array to hold core charge
    real(kind=dp), dimension(:), allocatable :: core_potential ! The hartree potential
    !due to core electrons, on a radial grid
    real(kind=dp), dimension(:), allocatable :: ae_potential ! The hartree potential of
    !the all-electron atom, on a radial grid
    real(kind=dp), dimension(:), allocatable :: val_potential ! The hartree potential
    !due to valence electrons of the free atom, on a radial grid
    real(kind=dp) :: ehar ! hartree energy
    real(kind=dp), dimension(:), allocatable :: rhorr_val,rhorr_core,rhorr_ae
    !core, valence, and ae charge densities, times r^2

!END: TEMSTEP edits

    integer       :: i,j,n,nb,nbb,nref,no,l,indx,jndx,ierr,nn,nnmax  ! add nn,nnmax for ELNES
    real(kind=dp) :: temp_rc,dE_use,scale,temp_rc2,dr

    real(kind=dp), dimension(:), allocatable :: rs_core_current, work, work2
    ! real(kind=dp), dimension(:,:,:), allocatable ::orbr

    character(len=file_maxpath) :: definition,temp_def
    character(len=80)           :: ctemp,ctemp2

    real(kind=dp), parameter :: dE = 0.1_dp ! * The energy step for multiple projectors

    logical :: swap,user,overide

    ! ** Define comments

    integer, parameter :: ncomm=4
    character(len=1)   :: comment(ncomm)

    integer                         :: status

    call trace_entry('otfg_define_psp',status)

    comment(1) = '#'           ! Define comment characters
    comment(2) = ';'
    comment(3) = '!'
    comment(4) = achar(13)     ! ^M

    ! ** Set the working pseudopotential definition

    select case(trim(current_cell%species_pot(nsp)))
    case("NULL.NULL")
       definition = adjustl(ps_def(aeat%Z,ps_set(nsp)))
       do j=1,ncomm
          i=index(definition,comment(j))
          if(i.gt.0) definition(i:) = ' ' ! Strip out the comment lines
       end do
    case default
       definition = adjustl(current_cell%species_pot(nsp))
    end select

    ! ** Record the pseudopotential set used

    select case(ps_set(nsp))
    case(999)
       ps%set="EXT"
    case(1)
       ps%set="QC5"
    case(2)
       ps%set="HARD"
    case(3)
       ps%set="MS"
    case(7)
       ps%set="C7"
    case(8)
       ps%set="C8"
    case(9)
       ps%set="C9"
    case(10)
       ps%set="NCP9"
    case(11)
       ps%set="C17"
    case(12)
       ps%set="NCP17"
    case(13)
       ps%set="C18"
    case(14)
       ps%set="NCP18"
    case(15)
       ps%set="C19"
    case(16)
       ps%set="NCP19"
    case(17)
       ps%set="SOC19"
    end select

    if(index(trim(definition),' ')>0) call io_abort('Error: otfg_define_psp: incorrect syntax, remove blanks')

    ! ** Output testing files?

    testing = .false.
    indx = index(definition,'[')
    if(indx/=0.or.io_code_present(current_params%devel_code,'test_pspot','OTFG')) testing = .true.

    ! ** Check if there any pseudisation flags are set

    indx = index(definition,'(')
    if(indx/=0) then

       ! * Remove the flag string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,')')
       ps%flags = ctemp(:indx-1)
       if(indx==0) call io_abort('Error: otfg_define_psp: incorrect syntax for pseudisation flags')
    end if

    ! ** Check if there any configuration adjustment required

    indx = index(definition,'{')
    if(indx/=0) then

       ! * Remove the configuration adjust string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,'}')
       if(indx==0) call io_abort('Error: otfg_define_psp: incorrect syntax for config adjust')
    end if

    indx = index(definition,'[')
    if(indx/=0) then

       ! * Remove the configuration adjust string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,']')
       if(indx==0) call io_abort('Error: otfg_define_psp: incorrect syntax for config adjust')
    end if

    ! ** Set the local channel

    ps%pv = .false. ! * Don't do the polynomial fit by default

    indx = index(definition,'|')
    ctemp = definition(:indx-1)

    indx = index(ctemp,'=')
    if(indx>0) then
       overide = .true.
       read(ctemp(indx+1:),*,err=103) ps%local_e
       read(ctemp(:indx-1),*,err=103) ps%local_l
       ps%local_l = abs(ps%local_l)
       if(index(ctemp(:indx-1),'-')>0) ps%pv=.true.
    else
       overide = .false.
       ps%local_e    = 0.0_dp                  ! * The energy
       read(ctemp,*,err=103) ps%local_l
       ps%local_l = abs(ps%local_l)
       if(index(ctemp,'-')>0) ps%pv=.true.
    end if
    ps%local_norm = 0                       ! * Do not conserve the norm

    ! ** Find out about the core radii - if only one is specified, use defaults
    !    - note, they may be specified on a projector by projector basis in addition

    ! * Count the number of "|" -- if 5, then use defaults

    n=0
    do i=1,len_trim(definition)
       if(definition(i:i)=='|') n=n+1
    end do

    ! * Read the core radii in from the definition string

    jndx = index(definition,'|')+1

    if(n==5) then      ! * Only one rc specified - use defaults
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%local_rc
       temp_rc = ps%local_rc
       ps%rinner = temp_rc*0.7_dp ! * Default setting
       definition = definition(jndx+1+indx:) ! * Remove the portion read
    else if(n==7) then ! * local_rc, beta_rc and
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%local_rc
       jndx = jndx + 1 + indx
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) temp_rc
       jndx = jndx + 1 + indx
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%rinner(0)
       ps%rinner = ps%rinner(0)
       definition = definition(jndx+indx+1:) ! * Remove the portion read
    else
       call io_abort('Error: otfg_define_psp: unable to parse pseudopotential definition - too many "|"')
    end if

    temp_rc2 = temp_rc ! Keep a copy of default beta_rc to use for the gamma projectors

    ! ** Read in the definitions of COARSE, MEDIUM and FINE
    !    - note, the should be in atomic units (Ha)

    indx = index(definition(1:),'|')-1
    read(definition(1:indx),*,err=103) ps%coarse
    jndx = 2+indx
    indx = index(definition(jndx:),'|')-1
    read(definition(jndx:jndx-1+indx),*,err=103) ps%medium
    jndx = 1+jndx+indx
    indx = index(definition(jndx:),'|')-1
    read(definition(jndx:jndx-1+indx),*,err=103) ps%fine
    definition = definition(jndx+indx+1:) ! * Remove the portion read

    ! * Set PRECISE and EXTREME (we don't really know what EXTREME is)

    ps%precise = max(1.03_dp * ps%fine + 1.57_dp,1.2_dp*ps%fine)
    ps%extreme = 1.6_dp * ps%fine

    ! ** Set the cutoff energy

    select case(current_params%basis_precision)
    case('COARSE')  ; current_params%cut_off_energy = max(ps%coarse,current_params%cut_off_energy)
    case('MEDIUM')  ; current_params%cut_off_energy = max(ps%medium,current_params%cut_off_energy)
    case('FINE')    ; current_params%cut_off_energy = max(ps%fine,current_params%cut_off_energy)
    case('PRECISE') ; current_params%cut_off_energy = max(ps%precise,current_params%cut_off_energy)
    case('EXTREME') ; current_params%cut_off_energy = max(ps%extreme,current_params%cut_off_energy)
    case('NULL')
    case default
       call io_abort('Error: otfg_define_psp: don''t recognise current_params%basis_precision')
    end select

    ! ** Now read in the projector definitions

    ! * Count the number of reference states

    nref=1
    do i=1,len_trim(definition)
       if(definition(i:i)==':') nref=nref+1
    end do

    ! ** Find out how many beta functions there are -- so the arrays can be allocated

    temp_def    = definition
    ps%num_beta = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)
       nb = 0
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='N') nb=nb+1 ! * Count the norm conserving projectors
          if(ctemp(i:i)=='U') nb=nb+1 ! * Count the ultrasoft projectors
       end do
       if(nb==0) nb=2
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='L') nb=0 ! * Don't count the local channel
          if(ctemp(i:i)=='P') nb=0 ! * Don't count the dummy channel
       end do
       ps%num_beta = ps%num_beta + nb
    end do

    ! ** Allocate the ps%beta_* arrays

    if (allocated(ps%beta_l)) then
       deallocate(ps%beta_l,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_l in otfg_define_psp')
    end if
    allocate(ps%beta_l(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_l','otfg_define_psp')
    if (allocated(ps%beta_e)) then
       deallocate(ps%beta_e,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_e in otfg_define_psp')
    end if
    allocate(ps%beta_e(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_e','otfg_define_psp')
    if (allocated(ps%beta_rc)) then
       deallocate(ps%beta_rc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_rc in otfg_define_psp')
    end if
    allocate(ps%beta_rc(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_rc','otfg_define_psp')
    if (allocated(ps%beta_delta)) then
       deallocate(ps%beta_delta,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_delta in otfg_define_psp')
    end if
    allocate(ps%beta_delta(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_delta','otfg_define_psp')
    ps%beta_delta = 0.0_dp
    if (allocated(ps%beta_nrc)) then
       deallocate(ps%beta_nrc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_nrc in otfg_define_psp')
    end if
    allocate(ps%beta_nrc(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_nrc','otfg_define_psp')
    if (allocated(ps%beta_norm)) then
       deallocate(ps%beta_norm,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_norm in otfg_define_psp')
    end if
    allocate(ps%beta_norm(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_norm','otfg_define_psp')
    if (allocated(ps%beta_scheme)) then
       deallocate(ps%beta_scheme,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_scheme in otfg_define_psp')
    end if
    allocate(ps%beta_scheme(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_scheme','otfg_define_psp')

    ! ** Set the ps%beta_* arrays
    ps%pseudise = .false.
    temp_def    = definition
    nbb         = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)

       read(ctemp(1:2),'(2i1)',err=103) no,l ! * Find out the nl for this reference state

       ps%pseudise(no-l,l) = .true.

       if(ctemp(3:3)=='P') cycle ! * Just pseudise - no projectors etc.

       if(ctemp(3:3)=='L') then

          if(l/=ps%local_l) call io_abort('Error: otfg_define_psp: local channel conflict')

          if(.not.overide) ps%local_e = aeat%eval(no-l,l)

          if(ctemp(4:4)=='G') then
             ps%local_norm = 0
          else
             ps%local_norm = 1
          end if

       else

          ps%beta_e(nbb+1) = aeat%eval(no-l,l)

          nb = 0
          do i=1,len_trim(ctemp)
             if(ctemp(i:i)=='N'.or.ctemp(i:i)=='U') then
                nb=nb+1 ! * Count the projectors

                do j=i+1,len_trim(ctemp)
                   if(ctemp(j:j)=='N'.or.ctemp(j:j)=='U'.or.ctemp(j:j)=='H'.or.ctemp(j:j)=='G') exit
                end do
                ctemp2 = ctemp(i+1:j-1)

                indx = index(ctemp2,'@')
                if(indx>0) then
                   read(ctemp2(indx+1:),*,err=103) ps%beta_delta(nbb+nb)
                   ctemp2 = ctemp2(:indx-1)
                end if

                indx = index(ctemp2,'+')
                jndx = index(ctemp2,'-')
                indx = max(indx,jndx)

                if(indx>0) then
                   read(ctemp2(indx:),*,err=103) dE_use
                   user = .true.
                   if(indx==1) then
                      ctemp2 = ' '
                   else
                      if(ctemp2(indx-1:indx-1)/='=') ctemp2 = ctemp2(:indx-1)
                   end if

                else
                   dE_use = dE*real(nb-1,dp)
                   user = .false.
                end if

                indx = index(ctemp2,'=')
                if(indx>0) then
                   if(len_trim(ctemp2(indx+1:))<1) goto 103
                   read(ctemp2(indx+1:),*,err=103) dE_use
                   ctemp2 = ctemp2(:indx-1)
                   ps%beta_e(nbb+nb) = dE_use
                else
                   if(user) then
                      ps%beta_e(nbb+nb) = aeat%eval(no-l,l) + dE_use
                   else
                      if(dE_use>0.0_dp) then
                         ps%beta_e(nbb+nb) = 0.25*dE_use/dE
                      else
                         ps%beta_e(nbb+nb) = aeat%eval(no-l,l)
                      end if
                   end if
                end if

                if(len_trim(ctemp2)>0) then
                   read(ctemp2,*,err=103) ps%beta_rc(nbb+nb)
                else
                   ps%beta_rc(nbb+nb) = temp_rc
                end if

                ps%beta_l(nbb+nb) = l

                if(ctemp(i:i)=='N') then
                   ps%beta_norm(nbb+nb) = 1 ! * Norm conserving
                else
                   ps%beta_norm(nbb+nb) = 0 ! * Ultrasoft
                end if

             end if
          end do

          ! * Didn't find any modifiers, so use defaults

          if(nb==0) then
             nb=2

             ps%beta_l(nbb+1:nbb+2) = l

             ps%beta_norm(nbb+1:nbb+2) = 0
             ps%beta_rc(nbb+1:nbb+2)   = temp_rc

             ! * Set the energies

             ps%beta_e(nbb+2) = 0.25

          end if

          nbb = nbb + nb

       end if

    end do

    ! ** Reorder the beta arrays into increasing l
    !     - ensure block diagonality

101 swap = .false.
    do nb=1,ps%num_beta-1

       if(ps%beta_l(nb)>ps%beta_l(nb+1)) then
          swap = .true.
          l = ps%beta_l(nb)
          ps%beta_l(nb)   = ps%beta_l(nb+1)
          ps%beta_l(nb+1) = l
          l = ps%beta_norm(nb)
          ps%beta_norm(nb)   = ps%beta_norm(nb+1)
          ps%beta_norm(nb+1) = l
          temp_rc = ps%beta_rc(nb)
          ps%beta_rc(nb)   = ps%beta_rc(nb+1)
          ps%beta_rc(nb+1) = temp_rc
          temp_rc = ps%beta_e(nb)
          ps%beta_e(nb)   = ps%beta_e(nb+1)
          ps%beta_e(nb+1) = temp_rc
       end if

    end do

    if(swap) goto 101

    ! ** Find the _nrc and set the _rc to exact points

    ps%local_rc = ps%local_rc
    ps%beta_rc = ps%beta_rc
    ps%rinner = ps%rinner

    ps%local_rc = otfg_round_to_grid(ps%local_rc,ps%local_nrc,ab%r)

    do nb=1,ps%num_beta
       ps%beta_rc(nb) = otfg_round_to_grid(ps%beta_rc(nb),ps%beta_nrc(nb),ab%r)
    end do

    do l=0,2*ab%lmax
       ps%rinner(l) = otfg_round_to_grid(ps%rinner(l),ps%ninner(l),ab%r)
    end do

    ! ** Set ps%mesh

    ps%mesh   = ab%npts

    ! ** Find ps%kkbeta

    if(index(ps%flags,'es')>0) then

       if(index(ps%flags,'esr=')>0) then
          read(ps%flags(index(ps%flags,'esr=')+4:),*,err=103) ps%rcmax
          temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
          ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
       else
          ps%rcmax  = otfg_round_to_grid(5.0_dp*max(ps%local_rc,maxval(ps%beta_rc),maxval(ps%rinner)),nb,ab%r)
          temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
          ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
       end if

    else
       ps%rcmax  = otfg_round_to_grid(max(ps%local_rc,maxval(ps%beta_rc),maxval(ps%rinner)),nb,ab%r)
       temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
       ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
    end if

    ! ** Set the log grid

    ps%rlog(:ps%mesh)   = ab%r(:ps%mesh)
    ps%rab(:ps%mesh)    = ab%rab(:ps%mesh)

    ! ** Set the parameters for the augmentation function optimisation

    ps%nqpts  = 401

    if(index(ps%flags,'qc=')>0) then
       read(ps%flags(index(ps%flags,'qc=')+3:),*) ps%qc
       ps%qc = ps%qc*2.0_dp
    else
       ps%qc = 10.0_dp
    end if

    if(index(ps%flags,'Qc=')>0) then
       read(ps%flags(index(ps%flags,'Qc=')+3:),*) ps%qc
    end if

    ps%nderiv = 4
    ps%nextra = 3
    ps%nqf    = ps%nderiv+ps%nextra+1

    ! ** Set the parameters for the reciprocal space representation of the pseudopotential

    ps%num_points = 2001                                ! * This needs to be the same as in inquire_psp
    ps%gmax       = 52.91772086073879_dp                ! * So as to be the same as old CASTEP

    ! ** Allocate the core current

    allocate(rs_core_current(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rs_core_current','otfg_define_psp')

    allocate(work(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('work','otfg_define_psp')

    allocate(work2(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('work2','otfg_define_psp')

    ! ** Set the ionic charge of the pseudopotential, the atomic number, and the core charge

    ps%Z = aeat%Z

    ps%ionic_charge = real(aeat%Z,dp)
    ps%rs_core_charge = 0.0_dp
    do l = 0,ab%lmax
       do n = 1,aeat%norbmax

          if(.not.ps%pseudise(n,l)) then
             ps%ionic_charge = ps%ionic_charge - aeat%occ(n,l)
             ps%rs_core_charge(:ps%mesh) = ps%rs_core_charge(:ps%mesh) + aeat%occ(n,l)*aeat%orbr(:ps%mesh,n,l)**2
          end if
       end do
    end do

    !BEGIN: TEMSTEP edits
    ! Write radial grid, radial core charge, radial ae charge, and nuclear
    ! potentials for each species
    ! We will neglect ab%r(1), since this is 0 and causes singularities
    npts=ab%npts
    if (nsp.eq.1) then
       open(unit=400,file='one.proj_pot',status='unknown')
       write(400,*) current_cell%num_species
       write(400,*) npts-1
       close(400)
    else
       open(unit=400,file='one.proj_pot',access='append',status='old')
       write(400,*) npts-1
       close(400)
    end if
    allocate(core_potential(npts))
    allocate(ae_potential(npts))
    allocate(val_potential(npts))
    allocate(rhorr_core(npts))
    allocate(rhorr_val(npts))

    if (nsp.eq.1) then
        open(unit=403,file='debug_otfg.out',status='unknown')
        ionic_charge_dummy = ps%ionic_charge
        write(403,*) "ionic_charge_dummy: ",ionic_charge_dummy
        num_core_electrons=idnint(real(aeat%Z,kind=dp)-ionic_charge_dummy)
        write(403,*) "num_core_electrons: ",num_core_electrons
        close(403)
    else
        open(unit=403,file='debug_otfg.out',access='append',status='unknown')
        ionic_charge_dummy = ps%ionic_charge
        write(403,*) "ionic_charge_dummy: ",ionic_charge_dummy
        num_core_electrons=idnint(real(aeat%Z,kind=dp)-ionic_charge_dummy)
        write(403,*) "num_core_electrons: ",num_core_electrons
        close(403)
    end if

    do i=1,npts
       rhorr_core(i)=0.0_dp
    end do

    do i=1,npts
        rhorr_core(i) = ps%rs_core_charge(i)
    end do

    do i=1,npts
       rhorr_val(i)=aeat%rhorr(i)-rhorr_core(i)
    end do

    call otfg_ae_hartree(ab,rhorr_core,core_potential,ehar)
    call otfg_ae_hartree(ab,rhorr_val,val_potential,ehar)
    call otfg_ae_hartree(ab,aeat%rhorr,ae_potential,ehar)

    if (nsp.eq.1) then
       open(unit=401,file='two.proj_pot',status='unknown')
       write(401,*) aeat%Z
       do i=2,npts
          write(401,*) ab%r(i)
       end do
       do i=2,npts
          write(401,*) core_potential(i)
       end do
       do i=2,npts
          write(401,*) ae_potential(i)
       end do
       do i=2,npts
          write(401,*) aeat%V_nuc(i)/ab%r(i)
       end do
       do i=2,npts
          write(401,*) rhorr_core(i)
       end do
       do i=2,npts
          write(401,*) rhorr_val(i)
       end do
       do i=2,npts
          write(401,*) aeat%rhorr(i)
       end do
       close(401)
    else
       open(unit=401,file='two.proj_pot',access='append',status='old')
       write(401,*) aeat%Z
       do i=2,npts
          write(401,*) ab%r(i)
       end do
       do i=2,npts
          write(401,*) core_potential(i)
       end do
       do i=2,npts
          write(401,*) ae_potential(i)
       end do
       do i=2,npts
          write(401,*) aeat%V_nuc(i)/ab%r(i)
       end do
       do i=2,npts
          write(401,*) rhorr_core(i)
       end do
       do i=2,npts
          write(401,*) rhorr_val(i)
       end do
       do i=2,npts
          write(401,*) aeat%rhorr(i)
       end do
       close(401)
    end if

    if (nsp.eq.1) then
       open(unit=410,file='rho_core.out',status='unknown')
       do i=1,npts
          write(410,*) ab%r(i),rhorr_core(i)
       end do
       close(410)
    else
       open(unit=410,file='rho_core.out',access='append',status='old')
       do i=1,npts
          write(410,*) ab%r(i),rhorr_core(i)
       end do
       close(410)
    end if


    if (nsp.eq.1) then
       open(unit=420,file='rho_val.out',status='unknown')
       do i=1,npts
          write(420,*) ab%r(i),rhorr_val(i)
       end do
       close(420)
    else
       open(unit=420,file='rho_val.out',access='append',status='old')
       do i=1,npts
          write(420,*) ab%r(i),rhorr_val(i)
       end do
       close(420)
    end if

    if (nsp.eq.1) then
       open(unit=430,file='rho_ae.out',status='unknown')
       do i=1,npts
          write(430,*) ab%r(i),aeat%rhorr(i)
       end do
       close(430)
    else
       open(unit=430,file='rho_ae.out',access='append',status='old')
       do i=1,npts
          write(430,*) ab%r(i),aeat%rhorr(i)
       end do
       close(430)
    end if

    ! ** Set the core and total charges on radial logarithmic grids for X-ray structure factor calculations

    ps%ae_rhorr_core = ps%rs_core_charge(:)
    ps%ae_rhorr_iam(:ps%mesh) = aeat%rhorr(:ps%mesh)

    ! rs_core_charge contains the real-space core charge density. We need to compute the core kinetic energy density now.
    ! KED computed from the laplacian of the charge density:
    ! \nabla \psi_{n,l}^*(\br) \nabla \psi_{n,l}(\br) = \frac{1}{4} \Biggl[
    ! \nabla^2 \rho(\br) - \sum_{n,l} f_{n,l} \Bigl( \nabla^2 \psi_{n,l}^*(\br) \psi_{n,l}(\br) +
    ! \psi_{n,l}^*(\br) \nabla^2 \psi_{n,l}(\br) \Bigr) \Biggr]


    ! Laplacian:
    ! \nabla^2 = \frac{1}{r^2} \Biggl[ \frac{\partial}{\partial r} \Bigl( r^2 \frac{\partial}{\partial r} \Bigr) - \hat{l}^2 \Biggr]

    ! the density is spherically symmetric so we only need to worry about the radial derivatives.
    ! the radial part can be re-written in terms of a second derivative which we will use
    ! throughout:
    ! \frac{1}{r^2} \frac{\partial}{\partial r} \Biggl( r^2 \frac{\partial}{\partial r} \Biggr) =
    ! \frac{1}{r} \Biggl( \frac{\partial^2}{\partial r^2} r \Biggr)

    ! the density is in fact r^2 \rho, so we need to divide by r first (see above)
    work(1) = 0.0_dp
    work(2:ps%mesh) = ps%rs_core_charge(2:ps%mesh) / ab%r(2:ps%mesh)
    call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,work,work2,ps%rs_core_kinetic_energy_density)
    ! **multiply** by r to get r^2 \nabla^2 \rho
    ps%rs_core_kinetic_energy_density(1:ps%mesh) = ps%rs_core_kinetic_energy_density(1:ps%mesh) * ab%r(1:ps%mesh)

    ! now subtract the laplacians of the orbitals. they are real, hence the factor of 2
    ! NB we work with orbital * r

    do l = 0,ab%lmax
       do n = 1,aeat%norbmax

          if(.not.ps%pseudise(n,l)) then
             call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,aeat%orbr(:,n,l),work,work2)

             work(2:ab%npts) = work2(2:ab%npts) - l*(l+1)*aeat%orbr(2:ab%npts,n,l) / ab%r(2:ab%npts)**2
             ps%rs_core_kinetic_energy_density(:ps%mesh) = ps%rs_core_kinetic_energy_density(:ps%mesh) - &
                  aeat%occ(n,l)*2.0_dp*work*aeat%orbr(:ps%mesh,n,l)

          end if
       end do
    end do

    ! extrapolate the first element of the KED
    ps%rs_core_kinetic_energy_density(1) = ps%rs_core_kinetic_energy_density(2) + &
         ( ps%rs_core_kinetic_energy_density(3)-ps%rs_core_kinetic_energy_density(2) ) / &
         ( ab%r(3) - ab%r(2) ) * ( ab%r(1) - ab%r(2) )

    ! but don't let it go below zero
    ps%rs_core_kinetic_energy_density(1) = max(0.0_dp, ps%rs_core_kinetic_energy_density(1))
    ! divide by two, see TeX formula
    ps%rs_core_kinetic_energy_density =  0.25_dp * ps%rs_core_kinetic_energy_density ! Definition with half
    ! zero out far away JIC
    ps%rs_core_kinetic_energy_density(ps%mesh+1:) = 0.0_dp
    ! ** Calculate the core current and the chemical shielding of the core electrons

    do n=2,ab%npts
       rs_core_current(n)=ps%rs_core_charge(n)/ab%r(n)
    end do
    rs_core_current(1)=0.0_dp

    ps%core_nmr_shift=(mu_0_si*speed_light_si*fine_structure_si*elementary_charge_si**2*&
         otfg_radin(ab%npts,ab%rab,rs_core_current)*1.0e6_dp)/(6.0_dp*planck_si)

    ! We now compute the magnetic suceptibility of the core electrons
    ! rather than allocate a new array we'll reuse  rs_core_current

    do n=2,ab%npts
       rs_core_current(n)=ps%rs_core_charge(n)*ab%r(n)*ab%r(n)
    end do
    rs_core_current(1)=0.0_dp

    ps%core_mag_sus=-1.0_dp*io_atomic_to_unit(1.0_dp,'cm')**3.d0*avogadro_si*1.d6*(fine_structure_si**2.d0)* &
         otfg_radin(ab%npts,ab%rab,rs_core_current)/6.0_dp

    ! ** Set the NLCC core charge

    ps%core_charge = (otfg_radin(ps%kkbeta,ab%rab,ps%rs_core_charge)>0.0_dp)

    if(index(ps%flags,'nonlcc')>0) ps%core_charge = .false.

    if(ps%core_charge) then

       ps%nnlcc = ps%ninner(0)
       ps%rnlcc = ps%rinner(0)

       if(compat7) then
          ps%rs_core_charge(:ab%npts)=otfg_core_pseudise_compat7(ps%rs_core_charge,ps,ab)
          ps%rs_core_kinetic_energy_density = otfg_core_pseudise_compat7(ps%rs_core_kinetic_energy_density,ps,ab)
       else
          ps%rs_core_charge(:ab%npts)=otfg_core_pseudise(ps%rs_core_charge,ps,ab)
          ps%rs_core_kinetic_energy_density = otfg_core_pseudise(ps%rs_core_kinetic_energy_density,ps,ab)
       end if

       call otfg_radial_transform(0,0,0,ps%mesh,ps%rlog,ps%rab,ps%num_points,ps%gmax,ps%rs_core_charge,ps%core_radial_charge)
       call otfg_derivative(ps%num_points,ps%gmax,ps%core_radial_charge,ps%dcore_radial_charge_dg)

       call otfg_radial_transform(0,0,0,ps%mesh,ps%rlog,ps%rab,ps%num_points,ps%gmax,&
            ps%rs_core_kinetic_energy_density,ps%core_radial_kinetic_energy_density)
       call otfg_derivative(ps%num_points,ps%gmax,ps%core_radial_kinetic_energy_density,&
            ps%dcore_radial_kinetic_energy_density_dg)
    else

       ps%rs_core_charge            = 0.0_dp
       ps%core_radial_charge(:)     = 0.0_dp
       ps%dcore_radial_charge_dg(:) = 0.0_dp

       ps%rs_core_kinetic_energy_density = 0.0_dp
       ps%core_radial_kinetic_energy_density = 0.0_dp
       ps%dcore_radial_kinetic_energy_density_dg = 0.0_dp
    end if

    deallocate(rs_core_current,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating rs_core_current in otfg_define_psp')

    temp_rc=temp_rc2 ! Keep a copy of default beta_rc to use for the gamma projectors

    ! ** Find out how many gamma functions there are -- so the arrays can be allocated
    !    - By default this will be the same as the number of beta projectors

    temp_def     = definition
    ps%num_gamma = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)
       nb = 0
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='H') nb=nb+1 ! * Count the norm conserving projectors
          if(ctemp(i:i)=='G') nb=nb+1 ! * Count the ultrasoft projectors
       end do
       ps%num_gamma = ps%num_gamma + nb
    end do
    ps%num_gamma= ps%num_gamma +ps%num_beta ! Total gamma = beta + extra

    ps%bl_gmax=ps%gmax ! * As beta proj
    ps%bl_num_points=ps%num_points


    ! ** Allocate the ps%gamma_* arrays

    if (allocated(ps%gamma_l)) then
       deallocate(ps%gamma_l,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_l in otfg_define_psp')
    end if
    allocate(ps%gamma_l(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_l','otfg_define_psp')
    if (allocated(ps%gamma_e)) then
       deallocate(ps%gamma_e,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_e in otfg_define_psp')
    end if
    allocate(ps%gamma_e(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_e','otfg_define_psp')
    if (allocated(ps%gamma_rc)) then
       deallocate(ps%gamma_rc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_rc in otfg_define_psp')
    end if
    allocate(ps%gamma_rc(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_rc','otfg_define_psp')
    if (allocated(ps%gamma_nrc)) then
       deallocate(ps%gamma_nrc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_nrc in otfg_define_psp')
    end if
    allocate(ps%gamma_nrc(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_nrc','otfg_define_psp')
    if (allocated(ps%gamma_norm)) then
       deallocate(ps%gamma_norm,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_norm in otfg_define_psp')
    end if
    allocate(ps%gamma_norm(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_norm','otfg_define_psp')

    if (allocated(ps%gamma_scheme)) then
       deallocate(ps%gamma_scheme,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_scheme in otfg_define_psp')
    end if
    allocate(ps%gamma_scheme(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_scheme','otfg_define_psp')

    ps%gamma_l=0

    ! ** Set the ps%gamma_* arrays

    temp_def    = definition
    nbb         = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)

       read(ctemp(1:2),'(2i1)',err=103) no,l ! * Find out the nl for this reference state

       if(((ctemp(3:3)=='L').or.(ctemp(3:3)=='P')).and..not.(ctemp(4:4)=='N'.or.ctemp(4:4)=='U'&
            .or.ctemp(4:4)=='G'.or.ctemp(4:4)=='H')) then

          !This is the local or dummy channel, there are no gammas unless explicitly requested.

       else

          ps%gamma_e(nbb+1) = aeat%eval(no-l,l)

          nb = 0
          do i=1,len_trim(ctemp)
             if(ctemp(i:i)=='N'.or.ctemp(i:i)=='U'.or.ctemp(i:i)=='G'.or.ctemp(i:i)=='H') then
                nb=nb+1 ! * Count the projectors

                do j=i+1,len_trim(ctemp)
                   if(ctemp(j:j)=='N'.or.ctemp(j:j)=='U'.or.ctemp(j:j)=='G'.or.ctemp(j:j)=='H') exit
                end do
                ctemp2 = ctemp(i+1:j-1)

                indx = index(ctemp2,'@')
                if(indx>0) then
                   ctemp2 = ctemp2(:indx-1)
                end if

                indx = index(ctemp2,'+')
                jndx = index(ctemp2,'-')
                indx = max(indx,jndx)


                if(indx>0) then
                   read(ctemp2(indx:),*,err=103) dE_use
                   user=.true.
                   if(indx==1) then
                      ctemp2 = ' '
                   else
                      if(ctemp2(indx-1:indx-1)/='=') ctemp2 = ctemp2(:indx-1)
                   end if
                else
                   dE_use = dE*real(nb-1,dp)
                   user = .false.
                end if


                indx = index(ctemp2,'=')
                if(indx>0) then
                   read(ctemp2(indx+1:),*,err=103) dE_use
                   ctemp2 = ctemp2(:indx-1)
                   ps%gamma_e(nbb+nb) = dE_use
                else
                   if(user) then
                      ps%gamma_e(nbb+nb) = aeat%eval(no-l,l) + dE_use
                   else
                      if(dE_use>0.0_dp) then
                         ps%gamma_e(nbb+nb) = 0.25*dE_use/dE
                      else
                         ps%gamma_e(nbb+nb) = aeat%eval(no-l,l)
                      end if
                   end if
                end if

                if(len_trim(ctemp2)>0) then
                   read(ctemp2,*,err=103) ps%gamma_rc(nbb+nb)
                else
                   ps%gamma_rc(nbb+nb) = temp_rc
                end if

                ps%gamma_l(nbb+nb) = l

                if(ctemp(i:i)=='N'.or.ctemp(i:i)=='H') then
                   ps%gamma_norm(nbb+nb) = 1 ! * Norm conserving
                else
                   ps%gamma_norm(nbb+nb) = 0 ! * Ultrasoft
                end if

             end if
          end do


          ! * Didn't find any modifiers, so use defaults

          if(nb==0) then
             nb=2

             ps%gamma_l(nbb+1:nbb+2) = l

             ps%gamma_norm(nbb+1:nbb+2) = 0
             ps%gamma_rc(nbb+1:nbb+2)   = temp_rc

             ! * Set the energies

             ps%gamma_e(nbb+2) = 0.25

          end if !no modifiers

          nbb = nbb + nb

       end if !if local channel with gammas or not

    end do !loop over states

102 swap = .false.
    do nb=1,ps%num_gamma-1

       if(ps%gamma_l(nb)>ps%gamma_l(nb+1)) then
          swap = .true.
          l = ps%gamma_l(nb)
          ps%gamma_l(nb)   = ps%gamma_l(nb+1)
          ps%gamma_l(nb+1) = l
          l = ps%gamma_norm(nb)
          ps%gamma_norm(nb)   = ps%gamma_norm(nb+1)
          ps%gamma_norm(nb+1) = l
          temp_rc = ps%gamma_rc(nb)
          ps%gamma_rc(nb)   = ps%gamma_rc(nb+1)
          ps%gamma_rc(nb+1) = temp_rc
          temp_rc = ps%gamma_e(nb)
          ps%gamma_e(nb)   = ps%gamma_e(nb+1)
          ps%gamma_e(nb+1) = temp_rc
       end if

    end do

    if(swap) goto 102

    ! ** Find the _nrc and set the _rc to exact points

    do nb=1,ps%num_gamma
       ps%gamma_rc(nb) = otfg_round_to_grid(ps%gamma_rc(nb),ps%gamma_nrc(nb),ab%r)
    end do

    ! ** Do any charge augmentation?

    ps%aug = .true.

    if(all(ps%beta_norm==1)) ps%aug = .false.

    if(index(ps%flags,'aug')>0) ps%aug = .true.

    ! ** Is the pseudopotential "all-electron"?

    ps%scpsp = .false.

    if(index(ps%flags,'scpsp')>0) ps%scpsp = .true.

    ! ** Sort out the initialisation of the boundary conditions

    ! * First the grids and radii

    bc(nsp)%Rc    = max(maxval(ps%beta_rc),ps%local_rc)*1.1_dp  ! * 10% more than the max core radius
    bc(nsp)%Rc    = otfg_round_to_grid(bc(nsp)%Rc,n,ab%r)
    bc(nsp)%ngc   = ngrid/6

    dr = bc(nsp)%Rc/real(bc(nsp)%ngc-1,dp)

    do i=1,ngrid
       bc(nsp)%urg(i) = real(i-1,dp)*dr
    end do
    bc(nsp)%rmax = bc(nsp)%urg(ngrid)

    ! * The angular momentum channels

    do l=4,0,-1
       if(any(ps%pseudise(:,l))) exit
    end do

    bc(nsp)%lmax = l

    ! **  define data associated with core orbitals for ELNES

    temp_def = definition

    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if

       ! * Find out the principal quantum number label

       read(ctemp(1:1),*) nn

       if (n==1) then
          nnmax=nn
       else
          nnmax=max(nn,nnmax)
       end if

       temp_def = temp_def(indx+1:)

    end do

!!!!!!! TEMPORARY FIX :::
    nnmax=min(6,nnmax)

    select case(nnmax)
    case(1,2)
       ps%num_core_projectors = 1  ! 1s
       ps%num_core_radials = 1 ! 1s
       ps%core_radials_n(1:1) = (/1/) ! 1
       ps%core_radials_l(1:1) = (/0/) ! s
    case(3)
       ps%num_core_projectors = 5  ! 1s,2s,2px,2py,2pz
       ps%num_core_radials    = 3  ! 1s,2s,2p
       ps%core_radials_n(1:3) = (/1,2,1/)
       ps%core_radials_l(1:3) = (/0,0,1/)
    case(4)
       ps%num_core_projectors = 9
       ps%num_core_radials = 5      ! 1s,2s,2p,3s,3p
       ps%core_radials_n(1:5) = (/1,2,1,3,2/)
       ps%core_radials_l(1:5) = (/0,0,1,0,1/)
    case(5)
       ps%num_core_projectors = 18
       ps%num_core_radials = 8     ! 1s,2s,2p,3s,3p,3d,4s,4p
       ps%core_radials_n(1:8) = (/1,2,1,3,2,1,4,3/)
       ps%core_radials_l(1:8) = (/0,0,1,0,1,2,0,1/)
    case(6)
       ps%num_core_projectors = 27
       ps%num_core_radials = 11    ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p
       ps%core_radials_n(1:11) = (/1,2,1,3,2,1,4,3,2,5,4/)
       ps%core_radials_l(1:11) = (/0,0,1,0,1,2,0,1,2,0,1/)
    case(7)
       ps%num_core_projectors = 43
       ps%num_core_radials = 15   ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d,6s,6p
       ps%core_radials_n(1:15) = (/1,2,1,3,2,1,4,3,2,1,5,4,3,6,5/)
       ps%core_radials_l(1:15) = (/0,0,1,0,1,2,0,1,2,3,0,1,2,0,1/)
    case default
       call io_abort('Error in otfg_define_psp: can not define core_radials')
    end select

    if(allocated(work)) then
       deallocate(work,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating work in otfg_define_psp')
    endif

    if(allocated(work2)) then
       deallocate(work2,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating work2 in otfg_define_psp')
    endif

    call trace_exit('otfg_define_psp',status)

    return

103 call io_abort('Error: otfg_define_psp: problem parsing definition')

  end subroutine otfg_define_psp



  subroutine otfg_define_psp_dir(nsp,ab,aeat,ps)

    use io,         only : io_abort,io_allocate_abort,io_failure,io_open_file,io_atomic_to_unit, &
         io_code_present
    use parameters, only : current_params
    use cell,       only : current_cell
    use constants,  only : file_maxpath,mu_0_si,elementary_charge_si,electron_mass_si,speed_light_si,&
         fine_structure_si,planck_si,avogadro_si

    !=========================================================================!
    ! Define the pseudopotential                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) nsp     : input  : The species                                       !
    ! 2) ab      : input  : The atomic basis                                  !
    ! 3) aeat    : input  : The initialised atom                              !
    ! 4) ps      : output : The pseudopotential                               !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) ae_atom                                                              !
    ! 3) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) parameters                                                           !
    ! 3) cell                                                                 !
    ! 4) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 1st July 2002                                   !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    integer,        intent(in)    :: nsp
    type(ae_basis), intent(in)    :: ab
    type(ae_atom),  intent(in)    :: aeat
    type(pspot),    intent(inout) :: ps

    !-------------------------------------------------------------------------!

    integer       :: i,j,n,nb,nbb,nref,no,l,indx,jndx,ierr,nn,nnmax  ! add nn,nnmax for ELNES
    real(kind=dp) :: temp_rc,dE_use,scale,temp_rc2,dr

    real(kind=dp), dimension(:), allocatable :: rs_core_current

    character(len=file_maxpath) :: definition,temp_def
    character(len=80)           :: ctemp,ctemp2

    real(kind=dp), parameter :: dE = 0.1_dp ! * The energy step for multiple projectors

    logical :: swap,user,overide

    ! ** Define comments

    integer, parameter :: ncomm=4
    character(len=1)   :: comment(ncomm)

    integer                         :: status

    call trace_entry('otfg_define_psp_dir',status)

    comment(1) = '#'           ! Define comment characters
    comment(2) = ';'
    comment(3) = '!'
    comment(4) = achar(13)     ! ^M

    ! ** Set the working pseudopotential definition

    select case(trim(current_cell%species_pot(nsp)))
    case("NULL.NULL")
       definition = adjustl(ps_def(aeat%Z,ps_set(nsp)))
       do j=1,ncomm
          i=index(definition,comment(j))
          if(i.gt.0) definition(i:) = ' ' ! Strip out the comment lines
       end do
    case default
       definition = adjustl(current_cell%species_pot(nsp))
    end select

    ! ** Check if we have simply been passed a scale factor

    if(index(definition,'|')==0) then
       read(definition,*,err=103) scale
       definition = adjustl(ps_def(aeat%Z,ps_set(nsp)))
    else
       scale = 1.0_dp
    end if


    ! ** Record the pseudopotential set used

    select case(ps_set(nsp))
    case(999)
       ps%set="EXT"
    case(1)
       ps%set="QC5"
    case(2)
       ps%set="HARD"
    case(3)
       ps%set="MS"
    case(7)
       ps%set="C7"
    case(8)
       ps%set="C8"
    case(9)
       ps%set="C9"
    case(10)
       ps%set="NCP9"
    case(11)
       ps%set="C17"
    case(12)
       ps%set="NCP17"
    case(13)
       ps%set="C18"
    case(14)
       ps%set="NCP18"
    case(15)
       ps%set="C19"
    case(16)
       ps%set="NCP19"
    case(17)
       ps%set="SOC19"
    end select


    if(index(trim(definition),' ')>0) call io_abort('Error: otfg_define_psp_dir: incorrect syntax, remove blanks')

    ! ** Output testing files?
    testing = .false.
    indx = index(definition,'[')
    if(indx/=0.or.io_code_present(current_params%devel_code,'test_pspot','OTFG')) testing = .true.
    ! ** Check if there any pseudisation flags are set

    indx = index(definition,'(')
    if(indx/=0) then

       ! * Remove the flag string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,')')
       ps%flags = ctemp(:indx-1)
       if(indx==0) call io_abort('Error: otfg_define_psp_dir: incorrect syntax for pseudisation flags')
    end if

    ! ** Check if there any configuration adjustment required

    indx = index(definition,'{')
    if(indx/=0) then

       ! * Remove the configuration adjust string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,'}')
       if(indx==0) call io_abort('Error: otfg_define_psp_dir: incorrect syntax for config adjust')
    end if

    indx = index(definition,'[')
    if(indx/=0) then

       ! * Remove the configuration adjust string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,']')
       if(indx==0) call io_abort('Error: otfg_define_psp_dir: incorrect syntax for config adjust')
    end if

    ! ** Set the local channel

    ps%pv = .false. ! * Don't do the polynomial fit by default

    indx = index(definition,'|')
    ctemp = definition(:indx-1)

    indx = index(ctemp,'=')
    if(indx>0) then
       overide = .true.
       read(ctemp(indx+1:),*,err=103) ps%local_e
       read(ctemp(:indx-1),*,err=103) ps%local_l
       ps%local_l = abs(ps%local_l)
       if(index(ctemp(:indx-1),'-')>0) ps%pv=.true.
    else
       overide = .false.
       ps%local_e    = 0.0_dp                  ! * The energy
       read(ctemp,*,err=103) ps%local_l
       ps%local_l = abs(ps%local_l)
       if(index(ctemp,'-')>0) ps%pv=.true.
    end if
    ps%local_norm = 0                       ! * Do not conserve the norm

    ! ** Find out about the core radii - if only one is specified, use defaults
    !    - note, they may be specified on a projector by projector basis in addition

    ! * Count the number of "|" -- if 5, then use defaults

    n=0
    do i=1,len_trim(definition)
       if(definition(i:i)=='|') n=n+1
    end do

    ! * Read the core radii in from the definition string

    jndx = index(definition,'|')+1

    if(n==5) then      ! * Only one rc specified - use defaults
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%local_rc
       temp_rc = ps%local_rc
       ps%rinner = temp_rc*0.7_dp ! * Default setting
       definition = definition(jndx+1+indx:) ! * Remove the portion read
    else if(n==7) then ! * local_rc, beta_rc and
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%local_rc
       jndx = jndx + 1 + indx
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) temp_rc
       jndx = jndx + 1 + indx
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%rinner(0)
       ps%rinner = ps%rinner(0)
       definition = definition(jndx+indx+1:) ! * Remove the portion read
    else
       call io_abort('Error: otfg_define_psp_dir: unable to parse pseudopotential definition - too many "|"')
    end if

    temp_rc2 = temp_rc ! Keep a copy of default beta_rc to use for the gamma projectors

    ! ** Read in the definitions of COARSE, MEDIUM and FINE
    !    - note, the should be in atomic units (Ha)

    indx = index(definition(1:),'|')-1
    read(definition(1:indx),*,err=103) ps%coarse
    jndx = 2+indx
    indx = index(definition(jndx:),'|')-1
    read(definition(jndx:jndx-1+indx),*,err=103) ps%medium
    jndx = 1+jndx+indx
    indx = index(definition(jndx:),'|')-1
    read(definition(jndx:jndx-1+indx),*,err=103) ps%fine
    definition = definition(jndx+indx+1:) ! * Remove the portion read

    ! * Set PRECISE and EXTREME (we don't really know what EXTREME is)

    ps%precise = max(1.03_dp * ps%fine + 1.57_dp,1.2_dp*ps%fine)
    ps%extreme = 1.6_dp * ps%fine

    ! ** Set the cutoff energy

    select case(current_params%basis_precision)
    case('COARSE')  ; current_params%cut_off_energy = max(ps%coarse,current_params%cut_off_energy)
    case('MEDIUM')  ; current_params%cut_off_energy = max(ps%medium,current_params%cut_off_energy)
    case('FINE')    ; current_params%cut_off_energy = max(ps%fine,current_params%cut_off_energy)
    case('PRECISE') ; current_params%cut_off_energy = max(ps%precise,current_params%cut_off_energy)
    case('EXTREME') ; current_params%cut_off_energy = max(ps%extreme,current_params%cut_off_energy)
    case('NULL')
    case default
       call io_abort('Error: otfg_define_psp_dir: don''t recognise current_params%basis_precision')
    end select

    ! ** Now read in the projector definitions

    ! * Count the number of reference states

    nref=1
    do i=1,len_trim(definition)
       if(definition(i:i)==':') nref=nref+1
    end do

    ! ** Find out how many beta functions there are -- so the arrays can be allocated

    temp_def    = definition
    ps%num_beta = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       read(ctemp(2:2),*) l
       temp_def = temp_def(indx+1:)
       nb = 0
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='N') nb=nb+1 ! * Count the norm conserving projectors
          if(ctemp(i:i)=='U') nb=nb+1 ! * Count the ultrasoft projectors
       end do
       if(nb==0) nb=2
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='L') nb=0 ! * Don't count the local channel
          if(ctemp(i:i)=='P') nb=0 ! * Don't count the dummy channel
       end do
       ps%num_beta = ps%num_beta + nb
       if(l>0) ps%num_beta = ps%num_beta + nb  ! twice as many for SOC
    end do

    ! ** Allocate the ps%beta_* arrays

    if (allocated(ps%beta_l)) then
       deallocate(ps%beta_l,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_l in otfg_define_psp_dir')
    end if
    allocate(ps%beta_l(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_l','otfg_define_psp_dir')
    if (allocated(ps%beta_2j)) then
       deallocate(ps%beta_2j,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_2j in otfg_define_psp_dir')
    end if
    allocate(ps%beta_2j(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_2j','otfg_define_psp_dir')
    if (allocated(ps%beta_e)) then
       deallocate(ps%beta_e,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_e in otfg_define_psp_dir')
    end if
    allocate(ps%beta_e(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_e','otfg_define_psp_dir')
    if (allocated(ps%beta_rc)) then
       deallocate(ps%beta_rc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_rc in otfg_define_psp_dir')
    end if
    allocate(ps%beta_rc(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_rc','otfg_define_psp_dir')
    if (allocated(ps%beta_delta)) then
       deallocate(ps%beta_delta,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_delta in otfg_define_psp_dir')
    end if
    allocate(ps%beta_delta(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_delta','otfg_define_psp_dir')
    ps%beta_delta = 0.0_dp
    if (allocated(ps%beta_nrc)) then
       deallocate(ps%beta_nrc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_nrc in otfg_define_psp_dir')
    end if
    allocate(ps%beta_nrc(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_nrc','otfg_define_psp_dir')
    if (allocated(ps%beta_norm)) then
       deallocate(ps%beta_norm,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_norm in otfg_define_psp_dir')
    end if
    allocate(ps%beta_norm(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_norm','otfg_define_psp_dir')
    if (allocated(ps%beta_scheme)) then
       deallocate(ps%beta_scheme,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_scheme in otfg_define_psp_dir')
    end if
    allocate(ps%beta_scheme(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_scheme','otfg_define_psp_dir')

    ! ** Set the ps%beta_* arrays
    ps%pseudise = .false.
    temp_def    = definition
    nbb         = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)

       read(ctemp(1:2),'(2i1)',err=103) no,l ! * Find out the nl for this reference state

       if(l==0) then
          ps%pseudise((no-l),l) = .true.
       else
          ps%pseudise(2*(no-l)-1,l) = .true.
          ps%pseudise(2*(no-l),l) = .true.
       endif

       if(ctemp(3:3)=='P') cycle ! * Just pseudise - no projectors etc.

       if(ctemp(3:3)=='L') then

          if(l/=ps%local_l) call io_abort('Error: otfg_define_psp_dir: local channel conflict')

          if(.not.overide) ps%local_e = aeat%eval(no-l,l)

          if(ctemp(4:4)=='G') then
             ps%local_norm = 0
          else
             ps%local_norm = 1
          end if

       else
          if(l==0) then
             ps%beta_e(nbb+1) = aeat%eval(no,l)
          else
             if(2*(no-l)-1>size(aeat%eval,1)) then
                ps%beta_e(nbb+1) = 0.0_dp
             else
                ps%beta_e(nbb+1) = aeat%eval(2*(no-l)-1,l)
             endif
             if(2*(no-l)>size(aeat%eval,1)) then
                ps%beta_e(nbb+2) = 0.0_dp
             else
                ps%beta_e(nbb+2) = aeat%eval(2*(no-l),l)
             endif
          endif

          nb = 0
          do i=1,len_trim(ctemp)
             if(ctemp(i:i)=='N'.or.ctemp(i:i)=='U') then
                nb=nb+1 ! * Count the projectors
                if(l>0.and.nb>1) nb=nb+1 !If j-split we go two at a time
                do j=i+1,len_trim(ctemp)
                   if(ctemp(j:j)=='N'.or.ctemp(j:j)=='U'.or.ctemp(j:j)=='H'.or.ctemp(j:j)=='G') exit
                end do
                ctemp2 = ctemp(i+1:j-1)

                indx = index(ctemp2,'@')
                ! '@' means that we shift the energy of the projector
                if(indx>0) then
                   read(ctemp2(indx+1:),*,err=103) ps%beta_delta(nbb+nb)
                   ! As we are going 2 at a time we set the value of the second one
                   if(l>0) ps%beta_delta(nbb+nb+1)=ps%beta_delta(nbb+nb)
                   ctemp2 = ctemp2(:indx-1)
                end if

                indx = index(ctemp2,'+')
                jndx = index(ctemp2,'-')
                indx = max(indx,jndx)

                if(indx>0) then
                   call io_abort('Error: otfg_define_psp_dir: not with SOC yet... (2) ')
                   read(ctemp2(indx:),*,err=103) dE_use
                   user = .true.
                   if(indx==1) then
                      ctemp2 = ' '
                   else
                      if(ctemp2(indx-1:indx-1)/='=') ctemp2 = ctemp2(:indx-1)
                   end if

                else
                   dE_use = dE*real(nb-1,dp) !! NO
                   user = .false.
                end if
                indx = index(ctemp2,'=')
                if(indx>0) then
                   call io_abort('Error: otfg_define_psp_dir: not with SOC yet... (3) ')
                   if(len_trim(ctemp2(indx+1:))<1) goto 103
                   read(ctemp2(indx+1:),*,err=103) dE_use
                   ctemp2 = ctemp2(:indx-1)
                   ps%beta_e(nbb+nb) = dE_use
                else
                   if(user) then
                      if(l==0) then
                         ps%beta_e(nbb+nb) = aeat%eval(no,l) + dE_use
                      else
                         if(2*(no-l)-1>size(aeat%eval,1)) then
                            ps%beta_e(nbb+nb) = 0.0_dp+ dE_use
                         else
                            ps%beta_e(nbb+nb) = aeat%eval(2*(no-l)-1,l) + dE_use
                         endif
                         if(2*(no-l)>size(aeat%eval,1)) then
                            ps%beta_e(nbb+nb+1) = 0.0_dp+ dE_use
                         else
                            ps%beta_e(nbb+nb+1) = aeat%eval(2*(no-l),l) + dE_use
                         endif
                      endif
                   else
                      if(dE_use>0.0_dp) then
                         ps%beta_e(nbb+nb) = 0.25!*dE_use/dE
                         if(l>0)ps%beta_e(nbb+nb+1) = 0.25!*dE_use/dE
                      else
                         if(l==0) then
                            ps%beta_e(nbb+nb) = aeat%eval(no,l)
                         else
                            if(2*(no-l)-1>size(aeat%eval,1)) then
                               ps%beta_e(nbb+nb) = 0.0_dp
                            else
                               ps%beta_e(nbb+nb) = aeat%eval(2*(no-l)-1,l)
                            endif
                            if(2*(no-l)>size(aeat%eval,1)) then
                               ps%beta_e(nbb+nb+1) = 0.0_dp
                            else
                               ps%beta_e(nbb+nb+1) = aeat%eval(2*(no-l),l)
                            endif
                         endif
                      end if
                   end if
                end if

                if(len_trim(ctemp2)>0) then
                   read(ctemp2,*,err=103) ps%beta_rc(nbb+nb)
                   if(l>0) ps%beta_rc(nbb+nb+1)=ps%beta_rc(nbb+nb)
                else
                   ps%beta_rc(nbb+nb) = temp_rc
                   if(l>0) ps%beta_rc(nbb+nb+1) = temp_rc
                end if

                ! The way I assign the j value is a fudge.
                ! I need a clear way to link eigenvalue with j (just copy it over from aeat??)
                ! this is the place things will break when we have more complex strings!! jry
                ps%beta_l(nbb+nb) = l
                if(l>0) ps%beta_l(nbb+nb+1) = l
                ps%beta_2j(nbb+nb) = abs(2*l+1)
                if(l>0) ps%beta_2j(nbb+nb+1) = abs(2*l-1)

                if(ctemp(i:i)=='N') then
                   ps%beta_norm(nbb+nb) = 1 ! * Norm conserving
                   if(l>0) ps%beta_norm(nbb+nb+1) = 1 ! * Norm conserving
                else
                   ps%beta_norm(nbb+nb) = 0 ! * Ultrasoft
                   if(l>0) ps%beta_norm(nbb+nb+1) = 0
                end if

             end if
          end do

          ! * Didn't find any modifiers, so use defaults

          if(nb==0) then
             !the default is 2 projectors per state
             !if s then 2 total; if p,d,f 4 in total

             nb=2

             !1st projector
             if(l>0) then
               ps%beta_l(nbb+1) = l
               ps%beta_norm(nbb+1) = 0
               ps%beta_rc(nbb+1)   = temp_rc
               ps%beta_2j(nbb+1) = abs(2*l+1)
               ps%beta_e(nbb+1) = aeat%eval(2*(no-l)-1,l)

                ps%beta_l(nbb+2) = l
                ps%beta_norm(nbb+2) = 0
                ps%beta_rc(nbb+2)   = temp_rc
                ps%beta_2j(nbb+2) = abs(2*l-1)
                ps%beta_e(nbb+2) = aeat%eval(2*(no-l),l)

                ps%beta_l(nbb+3) = l
                ps%beta_norm(nbb+3) = 0
                ps%beta_rc(nbb+3)   = temp_rc
                ps%beta_2j(nbb+3) = abs(2*l+1)
                ps%beta_e(nbb+3) = 0.25

                ps%beta_l(nbb+4) = l
                ps%beta_norm(nbb+4) = 0
                ps%beta_rc(nbb+4)   = temp_rc
                ps%beta_2j(nbb+4) = abs(2*l-1)
                ps%beta_e(nbb+4) = 0.25
              else
                ps%beta_l(nbb+1) = l
                ps%beta_norm(nbb+1) = 0
                ps%beta_rc(nbb+1)   = temp_rc
                ps%beta_2j(nbb+1) = abs(2*l+1)
                ps%beta_e(nbb+1) = aeat%eval(no,l)

                ps%beta_l(nbb+2) = l
                ps%beta_norm(nbb+2) = 0
                ps%beta_2j(nbb+2) = abs(2*l+1)
                ps%beta_rc(nbb+2)   = temp_rc
                ps%beta_e(nbb+2) = 0.25
             ! * Set the energies
           end if
          end if

          if(nb==1) then
             nbb = nbb + 1
             if(l>0) nbb = nbb + 1
          elseif(nb==3) then
             nbb = nbb + 2
             if(l>0) nbb = nbb + 2
          elseif(nb==2) then
             if(l==0) then
                nbb = nbb+nb
             else
               nbb = nbb+2*nb
             endif
          else
                call io_abort('Error: otfg_define_psp_dir: not with SOC yet... (7-2) ')
          endif

       end if

    end do

    ! ** Reorder the beta arrays into increasing l
    !     - ensure block diagonality

101 swap = .false.
    do nb=1,ps%num_beta-1

       if(ps%beta_l(nb)>ps%beta_l(nb+1)) then
          swap = .true.
          l = ps%beta_l(nb)
          ps%beta_l(nb)   = ps%beta_l(nb+1)
          ps%beta_l(nb+1) = l
          l = ps%beta_2j(nb)
          ps%beta_2j(nb)   = ps%beta_2j(nb+1)
          ps%beta_2j(nb+1) = l
          l = ps%beta_norm(nb)
          ps%beta_norm(nb)   = ps%beta_norm(nb+1)
          ps%beta_norm(nb+1) = l
          temp_rc = ps%beta_rc(nb)
          ps%beta_rc(nb)   = ps%beta_rc(nb+1)
          ps%beta_rc(nb+1) = temp_rc
          temp_rc = ps%beta_e(nb)
          ps%beta_e(nb)   = ps%beta_e(nb+1)
          ps%beta_e(nb+1) = temp_rc
       end if

    end do

    if(swap) goto 101


    ! ** Find the _nrc and set the _rc to exact points
    ps%local_rc = ps%local_rc*scale
    ps%beta_rc = ps%beta_rc*scale
    ps%rinner = ps%rinner*scale

    ps%local_rc = otfg_round_to_grid(ps%local_rc,ps%local_nrc,ab%r)

    do nb=1,ps%num_beta
       ps%beta_rc(nb) = otfg_round_to_grid(ps%beta_rc(nb),ps%beta_nrc(nb),ab%r)
    end do

    do l=0,2*ab%lmax
       ps%rinner(l) = otfg_round_to_grid(ps%rinner(l),ps%ninner(l),ab%r)
    end do


    ! ** Set ps%mesh

    ps%mesh   = ab%npts

    ! ** Find ps%kkbeta

    if(index(ps%flags,'es')>0) then

       if(index(ps%flags,'esr=')>0) then
          read(ps%flags(index(ps%flags,'esr=')+4:),*,err=103) ps%rcmax
          temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
          ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
       else
          ps%rcmax  = otfg_round_to_grid(5.0_dp*max(ps%local_rc,maxval(ps%beta_rc),maxval(ps%rinner)),nb,ab%r)
          temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
          ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
       end if

    else
       ps%rcmax  = otfg_round_to_grid(max(ps%local_rc,maxval(ps%beta_rc),maxval(ps%rinner)),nb,ab%r)
       temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
       ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
    end if

    ! ** Set the log grid

    ps%rlog(:ps%mesh)   = ab%r(:ps%mesh)
    ps%rab(:ps%mesh)    = ab%rab(:ps%mesh)

    ! ** Set the parameters for the augmentation function optimisation

    ps%nqpts  = 401

    if(index(ps%flags,'qc=')>0) then
       read(ps%flags(index(ps%flags,'qc=')+3:),*) ps%qc
       ps%qc = ps%qc*2.0_dp
    else
       ps%qc = 10.0_dp
    end if

    if(index(ps%flags,'Qc=')>0) then
       read(ps%flags(index(ps%flags,'Qc=')+3:),*) ps%qc
    end if

    ps%nderiv = 4
    ps%nextra = 3
    ps%nqf    = ps%nderiv+ps%nextra+1

    ! ** Set the parameters for the reciprocal space representation of the pseudopotential

    ps%num_points = 2001                                ! * This needs to be the same as in inquire_psp
    ps%gmax       = 52.91772086073879_dp                ! * So as to be the same as old CASTEP

    ! ** Allocate the core current

    allocate(rs_core_current(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rs_core_current','otfg_define_psp_dir')

    ! ** Set the ionic charge of the pseudopotential, the atomic number, and the core charge

    ps%Z = aeat%Z

    ps%ionic_charge = real(aeat%Z,dp)
    ps%rs_core_charge = 0.0_dp
    do l = 0,ab%lmax
       do n = 1,aeat%norbmax

          if(.not.ps%pseudise(n,l)) then
             ps%ionic_charge = ps%ionic_charge - aeat%occ(n,l)
             ps%rs_core_charge(:ps%mesh) = ps%rs_core_charge(:ps%mesh) + aeat%occ(n,l)*aeat%orbr(:ps%mesh,n,l)**2
             ps%rs_core_charge(:ps%mesh) = ps%rs_core_charge(:ps%mesh) + aeat%occ(n,l)*aeat%orbr_hack(:ps%mesh,n,l)**2
          end if
       end do
    end do

    !! Here we should expect some meta-GGA stuff

!!$
!!$    ! rs_core_charge contains the real-space core charge density. We need to compute the core kinetic energy density now.
!!$    ! KED computed from the laplacian of the charge density:
!!$    ! \nabla \psi_{n,l}^*(\br) \nabla \psi_{n,l}(\br) = \frac{1}{4} \Biggl[
!!$    ! \nabla^2 \rho(\br) - \sum_{n,l} f_{n,l} \Bigl( \nabla^2 \psi_{n,l}^*(\br) \psi_{n,l}(\br) +
!!$    ! \psi_{n,l}^*(\br) \nabla^2 \psi_{n,l}(\br) \Bigr) \Biggr]
!!$
!!$
!!$    ! Laplacian:
!!$    ! \nabla^2 = \frac{1}{r^2} \Biggl[ \frac{\partial}{\partial r} \Bigl( r^2 \frac{\partial}{\partial r} \Bigr) - \hat{l}^2 \Biggr]
!!$
!!$    ! the density is spherically symmetric so we only need to worry about the radial derivatives.
!!$    ! the radial part can be re-written in terms of a second derivative which we will use
!!$    ! throughout:
!!$    ! \frac{1}{r^2} \frac{\partial}{\partial r} \Biggl( r^2 \frac{\partial}{\partial r} \Biggr) =
!!$    ! \frac{1}{r} \Biggl( \frac{\partial^2}{\partial r^2} r \Biggr)
!!$
!!$    ! the density is in fact r^2 \rho, so we need to divide by r first (see above)
!!$    work(1) = 0.0_dp
!!$    work(2:ps%mesh) = ps%rs_core_charge(2:ps%mesh) / ab%r(2:ps%mesh)
!!$    call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,work,work2,ps%rs_core_kinetic_energy_density)
!!$    ! **multiply** by r to get r^2 \nabla^2 \rho
!!$    ps%rs_core_kinetic_energy_density(1:ps%mesh) = ps%rs_core_kinetic_energy_density(1:ps%mesh) * ab%r(1:ps%mesh)
!!$
!!$    ! now subtract the laplacians of the orbitals. they are real, hence the factor of 2
!!$    ! NB we work with orbital * r
!!$
!!$    do l = 0,ab%lmax
!!$       do n = 1,aeat%norbmax
!!$
!!$          if(.not.ps%pseudise(n,l)) then
!!$             call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,aeat%orbr(:,n,l),work,work2)
!!$
!!$             work(2:ab%npts) = work2(2:ab%npts) - l*(l+1)*aeat%orbr(2:ab%npts,n,l) / ab%r(2:ab%npts)**2
!!$             ps%rs_core_kinetic_energy_density(:ps%mesh) = ps%rs_core_kinetic_energy_density(:ps%mesh) - &
!!$                  aeat%occ(n,l)*2.0_dp*work*aeat%orbr(:ps%mesh,n,l)
!!$
!!$          end if
!!$       end do
!!$    end do
!!$
!!$    ! extrapolate the first element of the KED
!!$    ps%rs_core_kinetic_energy_density(1) = ps%rs_core_kinetic_energy_density(2) + &
!!$         ( ps%rs_core_kinetic_energy_density(3)-ps%rs_core_kinetic_energy_density(2) ) / &
!!$         ( ab%r(3) - ab%r(2) ) * ( ab%r(1) - ab%r(2) )
!!$
!!$    ! but don't let it go below zero
!!$    ps%rs_core_kinetic_energy_density(1) = max(0.0_dp, ps%rs_core_kinetic_energy_density(1))
!!$    ! divide by two, see TeX formula
!!$    ps%rs_core_kinetic_energy_density =  0.25_dp * ps%rs_core_kinetic_energy_density ! Definition with half
!!$    ! zero out far away JIC
!!$    ps%rs_core_kinetic_energy_density(ps%mesh+1:) = 0.0_dp

    ! ** Calculate the core current and the chemical shielding of the core electrons

    do n=2,ab%npts
       rs_core_current(n)=ps%rs_core_charge(n)/ab%r(n)
    end do
    rs_core_current(1)=0.0_dp

    ps%core_nmr_shift=(mu_0_si*speed_light_si*fine_structure_si*elementary_charge_si**2*&
         otfg_radin(ab%npts,ab%rab,rs_core_current)*1.0e6_dp)/(6.0_dp*planck_si)
    ! We now compute the magnetic suceptibility of the core electrons
    ! rather than allocate a new array we'll reuse  rs_core_current

    do n=2,ab%npts
       rs_core_current(n)=ps%rs_core_charge(n)*ab%r(n)*ab%r(n)
    end do
    rs_core_current(1)=0.0_dp

    ps%core_mag_sus=-1.0_dp*io_atomic_to_unit(1.0_dp,'cm')**3.d0*avogadro_si*1.d6*(fine_structure_si**2.d0)* &
         otfg_radin(ab%npts,ab%rab,rs_core_current)/6.0_dp

    ! ** Set the NLCC core charge

    ps%core_charge = (otfg_radin(ps%kkbeta,ab%rab,ps%rs_core_charge)>0.0_dp)

    if(index(ps%flags,'nonlcc')>0) ps%core_charge = .false.

    if(ps%core_charge) then

       ps%nnlcc = ps%ninner(0)
       ps%rnlcc = ps%rinner(0)

       if(compat7) then
          ps%rs_core_charge(:ab%npts)=otfg_core_pseudise_compat7(ps%rs_core_charge,ps,ab)
       else
          ps%rs_core_charge(:ab%npts)=otfg_core_pseudise(ps%rs_core_charge,ps,ab)
       end if

       call otfg_radial_transform(0,0,0,ps%mesh,ps%rlog,ps%rab,ps%num_points,ps%gmax,ps%rs_core_charge,ps%core_radial_charge)
       call otfg_derivative(ps%num_points,ps%gmax,ps%core_radial_charge,ps%dcore_radial_charge_dg)

    else

       ps%rs_core_charge            = 0.0_dp
       ps%core_radial_charge(:)     = 0.0_dp
       ps%dcore_radial_charge_dg(:) = 0.0_dp

    end if

    deallocate(rs_core_current,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating rs_core_current in otfg_define_psp_dir')

    temp_rc=temp_rc2 ! Keep a copy of default beta_rc to use for the gamma projectors

    ! ** Find out how many gamma functions there are -- so the arrays can be allocated
    !    - By default this will be the same as the number of beta projectors

    temp_def     = definition
    ps%num_gamma = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       read(ctemp(2:2),*) l
       temp_def = temp_def(indx+1:)
       nb = 0
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='H') nb=nb+1 ! * Count the norm conserving projectors
          if(ctemp(i:i)=='G') nb=nb+1 ! * Count the ultrasoft projectors
       end do
       ps%num_gamma = ps%num_gamma + nb
       if(l>0) ps%num_gamma = ps%num_gamma + nb ! twice as many for SOC
    end do
    ps%num_gamma= ps%num_gamma +ps%num_beta ! Total gamma = beta + extra

    ps%bl_gmax=ps%gmax ! * As beta proj
    ps%bl_num_points=ps%num_points


    ! ** Allocate the ps%gamma_* arrays

    if (allocated(ps%gamma_l)) then
       deallocate(ps%gamma_l,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_l in otfg_define_psp_dir')
    end if
    allocate(ps%gamma_l(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_l','otfg_define_psp_dir')
    if (allocated(ps%gamma_2j)) then
       deallocate(ps%gamma_2j,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_2j in otfg_define_psp_dir')
    end if
    allocate(ps%gamma_2j(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_2j','otfg_define_psp_dir')
    if (allocated(ps%gamma_e)) then
       deallocate(ps%gamma_e,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_e in otfg_define_psp_dir')
    end if
    allocate(ps%gamma_e(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_e','otfg_define_psp_dir')
    if (allocated(ps%gamma_rc)) then
       deallocate(ps%gamma_rc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_rc in otfg_define_psp_dir')
    end if
    allocate(ps%gamma_rc(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_rc','otfg_define_psp_dir')
    if (allocated(ps%gamma_nrc)) then
       deallocate(ps%gamma_nrc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_nrc in otfg_define_psp_dir')
    end if
    allocate(ps%gamma_nrc(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_nrc','otfg_define_psp_dir')
    if (allocated(ps%gamma_norm)) then
       deallocate(ps%gamma_norm,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_norm in otfg_define_psp_dir')
    end if
    allocate(ps%gamma_norm(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_norm','otfg_define_psp_dir')

    if (allocated(ps%gamma_scheme)) then
       deallocate(ps%gamma_scheme,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_scheme in otfg_define_psp_dir')
    end if
    allocate(ps%gamma_scheme(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_scheme','otfg_define_psp_dir')

    ps%gamma_l=0
    ps%gamma_2j=0

    ! ** Set the ps%gamma_* arrays

    temp_def    = definition
    nbb         = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)

       read(ctemp(1:2),'(2i1)',err=103) no,l ! * Find out the nl for this reference state

       if(((ctemp(3:3)=='L').or.(ctemp(3:3)=='P')).and..not.(ctemp(4:4)=='N'.or.ctemp(4:4)=='U'&
            .or.ctemp(4:4)=='G'.or.ctemp(4:4)=='H')) then

          !This is the local or dummy channel, there are no gammas unless explicitly requested.

       else


          if(l==0) then
             ps%gamma_e(nbb+1) = aeat%eval(no,l)
          else
             if(2*(no-l)-1>size(aeat%eval,1)) then
                ps%gamma_e(nbb+1) = 0.0_dp
             else
                ps%gamma_e(nbb+1) = aeat%eval(2*(no-l)-1,l)
             endif
             if(2*(no-l)>size(aeat%eval,1)) then
                ps%gamma_e(nbb+2) = 0.0_dp
             else
                ps%gamma_e(nbb+2) = aeat%eval(2*(no-l),l)
             endif
          endif

          nb = 0
          do i=1,len_trim(ctemp)
             if(ctemp(i:i)=='N'.or.ctemp(i:i)=='U'.or.ctemp(i:i)=='G'.or.ctemp(i:i)=='H') then
                nb=nb+1 ! * Count the projectors
                if(l>0.and.nb>1) nb=nb+1 !If j-split we go two at a time
                do j=i+1,len_trim(ctemp)
                   if(ctemp(j:j)=='N'.or.ctemp(j:j)=='U'.or.ctemp(j:j)=='G'.or.ctemp(j:j)=='H') exit
                end do
                ctemp2 = ctemp(i+1:j-1)

                indx = index(ctemp2,'@')
                if(indx>0) then
                   call io_abort('Error: otfg_define_psp_dir: not with SOC yet... (1G) ')
                   ! I think something is missing here
!                   ctemp2 = ctemp2(:indx-1)
                end if

                indx = index(ctemp2,'+')
                jndx = index(ctemp2,'-')
                indx = max(indx,jndx)


                if(indx>0) then
                   call io_abort('Error: otfg_define_psp_dir: not with SOC yet... (2G) ')
                   read(ctemp2(indx:),*,err=103) dE_use
                   user=.true.
                   if(indx==1) then
                      ctemp2 = ' '
                   else
                      if(ctemp2(indx-1:indx-1)/='=') ctemp2 = ctemp2(:indx-1)
                   end if
                else
                   dE_use = dE*real(nb-1,dp)
                   user = .false.
                end if


                indx = index(ctemp2,'=')
                if(indx>0) then
                   call io_abort('Error: otfg_define_psp_dir: not with SOC yet... (3G) ')
                   read(ctemp2(indx+1:),*,err=103) dE_use
                   ctemp2 = ctemp2(:indx-1)
                   ps%gamma_e(nbb+nb) = dE_use
                else
                    if(user) then
                      if(l==0) then
                         ps%gamma_e(nbb+nb) = aeat%eval(no,l) + dE_use
                      else
                         if(2*(no-l)-1>size(aeat%eval,1)) then
                            ps%gamma_e(nbb+nb) = 0.0_dp+ dE_use
                         else
                            ps%gamma_e(nbb+nb) = aeat%eval(2*(no-l)-1,l) + dE_use
                         endif
                         if(2*(no-l)>size(aeat%eval,1)) then
                            ps%gamma_e(nbb+nb+1) = 0.0_dp+ dE_use
                         else
                            ps%gamma_e(nbb+nb+1) = aeat%eval(2*(no-l),l) + dE_use
                         endif
                      endif
                   else
                      if(dE_use>0.0_dp) then
                         ps%gamma_e(nbb+nb) = 0.25!*dE_use/dE
                         if(l>0)ps%gamma_e(nbb+nb+1) = 0.25!*dE_use/dE
                      else
                         if(l==0) then
                            ps%gamma_e(nbb+nb) = aeat%eval(no,l)
                         else
                            if(2*(no-l)-1>size(aeat%eval,1)) then
                               ps%gamma_e(nbb+nb) = 0.0_dp
                            else
                               ps%gamma_e(nbb+nb) = aeat%eval(2*(no-l)-1,l)
                            endif
                            if(2*(no-l)>size(aeat%eval,1)) then
                               ps%gamma_e(nbb+nb+1) = 0.0_dp
                            else
                               ps%gamma_e(nbb+nb+1) = aeat%eval(2*(no-l),l)
                            endif
                         endif
                      end if
                   end if
                end if

                if(len_trim(ctemp2)>0) then
                   read(ctemp2,*,err=103) ps%gamma_rc(nbb+nb)
                   if(l>0) ps%gamma_rc(nbb+nb+1)=ps%gamma_rc(nbb+nb)
                else
                   ps%gamma_rc(nbb+nb) = temp_rc
                   if(l>0) ps%gamma_rc(nbb+nb+1) = temp_rc
                end if

                ! The way I assign the j value is a fudge.
                ! I need a clear way to link eigenvalue with j (just copy it over from aeat??)
                ! this is the place things will break when we have more complex strings!! jry
                ps%gamma_l(nbb+nb) = l
                if(l>0) ps%gamma_l(nbb+nb+1) = l
                ps%gamma_2j(nbb+nb) = abs(2*l+1)
                if(l>0) ps%gamma_2j(nbb+nb+1) = abs(2*l-1)


                if(ctemp(i:i)=='N'.or.ctemp(i:i)=='H') then
                   ps%gamma_norm(nbb+nb) = 1 ! * Norm conserving
                   if(l>0) ps%gamma_norm(nbb+nb+1) = 1 ! * Norm conserving
                else
                   ps%gamma_norm(nbb+nb) = 0 ! * Ultrasoft
                   if(l>0) ps%gamma_norm(nbb+nb+1) = 0
                end if

             end if
          end do
!!$
  if(nb==0) then
          ! * Didn't find any modifiers, so use defaults
          !the default is 2 projectors per state
          !if s then 2 total; if p,d,f 4 in total

          nb=2

          !1st projector

          if(l>0) then
            ps%gamma_l(nbb+1) = l
            ps%gamma_norm(nbb+1) = 0
            ps%gamma_rc(nbb+1)   = temp_rc
            ps%gamma_2j(nbb+1) = abs(2*l+1)
            ps%gamma_e(nbb+1) = aeat%eval(2*(no-l)-1,l)

            ps%gamma_l(nbb+2) = l
             ps%gamma_norm(nbb+2) = 0
             ps%gamma_rc(nbb+2)   = temp_rc
             ps%gamma_2j(nbb+2) = abs(2*l-1)
             ps%gamma_e(nbb+2) = aeat%eval(2*(no-l),l)

             ps%gamma_l(nbb+3) = l
             ps%gamma_norm(nbb+3) = 0
             ps%gamma_rc(nbb+3)   = temp_rc
             ps%gamma_2j(nbb+3) = abs(2*l+1)
             ps%gamma_e(nbb+3) = 0.25

             ps%gamma_l(nbb+4) = l
             ps%gamma_norm(nbb+4) = 0
             ps%gamma_rc(nbb+4)   = temp_rc
             ps%gamma_2j(nbb+4) = abs(2*l-1)
             ps%gamma_e(nbb+4) = 0.25
           else
             ps%gamma_l(nbb+1) = l
             ps%gamma_norm(nbb+1) = 0
             ps%gamma_rc(nbb+1)   = temp_rc
             ps%gamma_2j(nbb+1) = abs(2*l+1)
             ps%gamma_e(nbb+1) = aeat%eval(no,l)

             ps%gamma_l(nbb+2) = l
             ps%gamma_norm(nbb+2) = 0
             ps%gamma_rc(nbb+2)   = temp_rc
             ps%gamma_2j(nbb+2) = abs(2*l+1)
             ps%gamma_e(nbb+2) = 0.25
          ! * Set the energies
        end if
       end if

       if(nb==1) then
          nbb = nbb + 1
          if(l>0) nbb = nbb + 1
       elseif(nb==3) then
          nbb = nbb + 2
          if(l>0) nbb = nbb + 2
       elseif(nb==2) then
          if(l==0) then
             nbb = nbb+nb
          else
            nbb = nbb+2*nb
          endif
       else
             call io_abort('Error: otfg_define_psp_dir: not with SOC yet... (7-2) ')
       endif

       end if !if local channel with gammas or not

    end do !loop over states

102 swap = .false.
    do nb=1,ps%num_gamma-1

       if(ps%gamma_l(nb)>ps%gamma_l(nb+1)) then
          swap = .true.
          l = ps%gamma_l(nb)
          ps%gamma_l(nb)   = ps%gamma_l(nb+1)
          ps%gamma_l(nb+1) = l
          l = ps%gamma_2j(nb)
          ps%gamma_2j(nb)   = ps%gamma_2j(nb+1)
          ps%gamma_2j(nb+1) = l
          l = ps%gamma_norm(nb)
          ps%gamma_norm(nb)   = ps%gamma_norm(nb+1)
          ps%gamma_norm(nb+1) = l
          temp_rc = ps%gamma_rc(nb)
          ps%gamma_rc(nb)   = ps%gamma_rc(nb+1)
          ps%gamma_rc(nb+1) = temp_rc
          temp_rc = ps%gamma_e(nb)
          ps%gamma_e(nb)   = ps%gamma_e(nb+1)
          ps%gamma_e(nb+1) = temp_rc
       end if

    end do

    if(swap) goto 102

    ! ** Find the _nrc and set the _rc to exact points

    do nb=1,ps%num_gamma
       ps%gamma_rc(nb) = otfg_round_to_grid(ps%gamma_rc(nb),ps%gamma_nrc(nb),ab%r)
    end do

    ! ** Do any charge augmentation?

    ps%aug = .true.

    if(all(ps%beta_norm==1)) ps%aug = .false.

    if(index(ps%flags,'aug')>0) ps%aug = .true.

    ! ** Is the pseudopotential "all-electron"?

    ps%scpsp = .false.

    if(index(ps%flags,'scpsp')>0) ps%scpsp = .true.

    ! ** Sort out the initialisation of the boundary conditions

    ! * First the grids and radii

    bc(nsp)%Rc    = max(maxval(ps%beta_rc),ps%local_rc)*1.1_dp  ! * 10% more than the max core radius
    bc(nsp)%Rc    = otfg_round_to_grid(bc(nsp)%Rc,n,ab%r)
    bc(nsp)%ngc   = ngrid/6

    dr = bc(nsp)%Rc/real(bc(nsp)%ngc-1,dp)

    do i=1,ngrid
       bc(nsp)%urg(i) = real(i-1,dp)*dr
    end do
    bc(nsp)%rmax = bc(nsp)%urg(ngrid)

    ! * The angular momentum channels

    do l=4,0,-1
       if(any(ps%pseudise(:,l))) exit
    end do

    bc(nsp)%lmax = l

    ! **  define data associated with core orbitals for ELNES

    temp_def = definition

    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if

       ! * Find out the principal quantum number label

       read(ctemp(1:1),*) nn

       if (n==1) then
          nnmax=nn
       else
          nnmax=max(nn,nnmax)
       end if

       temp_def = temp_def(indx+1:)

    end do

!!!!!!! TEMPORARY FIX :::
    nnmax=min(6,nnmax)

    select case(nnmax)
      case(1,2)
        ps%num_core_projectors = 1  ! 1s
        ps%num_core_radials = 1 ! 1s
        ps%core_radials_n(1:1) = (/1/) ! 1
        ps%core_radials_l(1:1) = (/0/) ! s
      case(3)
        ps%num_core_projectors = 5  ! 1s,2s,2px,2py,2pz
        ps%num_core_radials    = 3  ! 1s,2s,2p
        ps%core_radials_n(1:3) = (/1,2,1/)
        ps%core_radials_l(1:3) = (/0,0,1/)
      case(4)
        ps%num_core_projectors = 9
        ps%num_core_radials = 5      ! 1s,2s,2p,3s,3p
        ps%core_radials_n(1:5) = (/1,2,1,3,2/)
        ps%core_radials_l(1:5) = (/0,0,1,0,1/)
      case(5)
        ps%num_core_projectors = 18
        ps%num_core_radials = 8     ! 1s,2s,2p,3s,3p,3d,4s,4p
        ps%core_radials_n(1:8) = (/1,2,1,3,2,1,4,3/)
        ps%core_radials_l(1:8) = (/0,0,1,0,1,2,0,1/)
      case(6)
        ps%num_core_projectors = 27
        ps%num_core_radials = 11    ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p
        ps%core_radials_n(1:11) = (/1,2,1,3,2,1,4,3,2,5,4/)
        ps%core_radials_l(1:11) = (/0,0,1,0,1,2,0,1,2,0,1/)
      case(7)
        ps%num_core_projectors = 43
        ps%num_core_radials = 15   ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d,6s,6p
        ps%core_radials_n(1:15) = (/1,2,1,3,2,1,4,3,2,1,5,4,3,6,5/)
        ps%core_radials_l(1:15) = (/0,0,1,0,1,2,0,1,2,3,0,1,2,0,1/)
      case default
        call io_abort('Error in otfg_define_psp_dir: can not define core_radials')
    end select
    call trace_exit('otfg_define_psp_dir',status)

    return

103 call io_abort('Error: otfg_define_psp_dir: problem parsing definition')

  end subroutine otfg_define_psp_dir


  subroutine otfg_construct_psp(ab,aeat,ps)

    use io,         only : io_abort,io_allocate_abort,stdout,io_open_file,&
         io_code_present,io_code_string,io_code_logical,io_code_real
    use comms,      only : my_node_id,root_node_id
    use parameters, only : current_params

    !=========================================================================!
    ! Construct the pseudopotential                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) ab      : input : The atomic basis                                   !
    ! 2) aeat    : input : The initialised atom                               !
    ! 3) ps      : inout : The pseudopotential                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) ae_atom                                                              !
    ! 3) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) comms                                                                !
    ! 3) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 1) psi_loc : The local partial wave                                     !
    ! 2) phi_loc : The pseudised local partial wave                           !
    ! 3) psi     : The projector partial waves                                !
    ! 4) phi     : The pseudised projector partial waves                      !
    ! 5) beta    : The beta projectors                                        !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 30th June 2002                                  !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    type(ae_basis), intent(in)    :: ab
    type(ae_atom),  intent(in)    :: aeat
    type(pspot),    intent(inout) :: ps

    !-------------------------------------------------------------------------!

    integer :: llmin,llmax,ierr
    integer :: na,nb,l,n,m,naa,nbb ! * Counters

    real(kind=dp) :: factor

    real(kind=dp), dimension(:),       allocatable :: psi_loc,phi_loc,V_loc,V_tau
    real(kind=dp), dimension(:,:),     allocatable :: psi,phi,beta,qqq,ddd,tvpsi,tvphi,psi_shift,phi_shift
    real(kind=dp), dimension(:,:,:),   allocatable :: qfunc, qfunc_bl
    real(kind=dp), dimension(:,:,:,:), allocatable :: qfcoef
    real(kind=dp), dimension(:,:,:,:), allocatable :: tfcoef

    integer,       dimension(:),       allocatable :: nnn

    real(kind=dp), dimension(:),       allocatable :: work,work2, FT_QT_func
    real(kind=dp), dimension(:),       allocatable :: laplacian_qfunc1, laplacian_qfunc2
    real(kind=dp), dimension(:,:),     allocatable :: laplacian_psi, laplacian_phi
    real(kind=dp), dimension(:,:,:,:), allocatable :: tfunc
    !real(kind=dp), dimension(:),     allocatable :: zora_K, zora_K_dr

    real(kind=dp), parameter :: ppm = 1.0e6_dp

    ! ** Gamma function arrays

    real(kind=dp), dimension(:,:),     allocatable :: psi_gamma,phi_gamma,gamma

    real(kind=dp), dimension(:,:),     allocatable :: bl_dia_nmrshield
    real(kind=dp), dimension(:,:),     allocatable :: bl_para_nmrshield
    real(kind=dp), dimension(:,:),     allocatable :: bl_para_damped
    real(kind=dp), dimension(:,:),     allocatable :: bl_efg_aug
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_fermi
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_dipolar
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_dipolar_cutoff
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_zora
    real(kind=dp), dimension(:,:),     allocatable :: bl_zke_aug_gtensor
    real(kind=dp), dimension(:,:),     allocatable :: bl_sodia_g_tensor
    real(kind=dp), dimension(:,:),     allocatable :: bl_sopara_g_tensor
    real(kind=dp), dimension(:),       allocatable :: V_ae,V_ps

    integer,       dimension(:),       allocatable :: nnn_bl

    integer                                        :: nn,nl,ng,nr,ll  ! For ELNES
    integer,       dimension(:),       allocatable :: nnn_nl        ! For ELNES
    real(kind=dp)                                  :: g             ! For ELNES
    real(kind=dp), dimension(:,:),     allocatable :: core_radials_r ! For ELNES
    real(kind=dp), dimension(:,:),     allocatable :: bl_elnes_aug,bl_ttt_temp  ! For ELNES
    real(kind=dp), dimension(:,:),     allocatable :: phi_gamma_temp,psi_gamma_temp ! For ELNES
    real(kind=dp) :: factor1,factor2   ! For ELNES

    real(kind=dp) :: betaphishift(1:ps%num_beta,1:ps%num_beta),phivphi(1:ps%num_beta,1:ps%num_beta),&
         phivphi_shift(1:ps%num_beta,1:ps%num_beta),phivphi_shift2(1:ps%num_beta,1:ps%num_beta)

    real(kind=dp), parameter :: eshift = -0.01_dp

    integer, parameter                   :: nflags = 5
    character(len=80), dimension(nflags) :: flags
    integer                              :: unit_zorafile

    character(len=1) :: cang

    logical :: on_root

    real(kind=dp) :: dtol = 1.0e-3_dp
    real(kind=dp) :: damp_scale = 0.5_dp

    integer                         :: status

    call trace_entry('otfg_construct_psp',status)

    on_root = my_node_id==root_node_id

    if(io_code_present(current_params%devel_code,'dtol','JC')) then
       dtol = io_code_real(current_params%devel_code,'dtol','JC')
    end if

    if(io_code_present(current_params%devel_code,'damp_scale','JC')) then
       damp_scale = io_code_real(current_params%devel_code,'damp_scale','JC')
    end if

    ! ** Allocate the local arrays

    allocate(psi_loc(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_loc','otfg_construct_psp')
    allocate(phi_loc(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi_loc','otfg_construct_psp')
    allocate(V_loc(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_loc','otfg_construct_psp')

    allocate(V_tau(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_tau','otfg_construct_psp')
    V_tau = 0.0_dp
    !allocate(zora_K(ab%npts),stat=ierr)
    !if (ierr/=0) call io_allocate_abort('zora_K','otfg_construct_psp')
    !allocate(zora_K_dr(ab%npts),stat=ierr)
    !if (ierr/=0) call io_allocate_abort('zora_K_dr','otfg_construct_psp')

    allocate(psi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi','otfg_construct_psp')
    allocate(phi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi','otfg_construct_psp')
    allocate(tvpsi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('tvpsi','otfg_construct_psp')
    allocate(tvphi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('tvphi','otfg_construct_psp')
    allocate(psi_shift(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_shift','otfg_construct_psp')
    allocate(phi_shift(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi_shift','otfg_construct_psp')
    allocate(beta(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('beta','otfg_construct_psp')
    allocate(qqq(ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('qqq','otfg_construct_psp')
    allocate(ddd(ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ddd','otfg_construct_psp')
    allocate(qfunc(ab%npts,ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('qfunc','otfg_construct_psp')
    allocate(qfunc_bl(ab%npts,ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('qfunc_bl','otfg_construct_psp')
    allocate(qfcoef(ps%nqf,0:2*ab%lmax,ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('qfcoef','otfg_construct_psp')
    allocate(tfcoef(ps%nqf,0:2*ab%lmax,ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('tfcoef','otfg_construct_psp')
    allocate(work(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('work','otfg_construct_psp')
    allocate(work2(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('work2','otfg_construct_psp')
    allocate(laplacian_qfunc1(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('laplacian_qfunc1','otfg_construct_psp')
    allocate(laplacian_qfunc2(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('laplacian_qfunc2','otfg_construct_psp')
    allocate(tfunc(ab%npts,0:2*ab%lmax,ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('tfunc','otfg_construct_psp')
    allocate(laplacian_psi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('laplacian_psi','otfg_construct_psp')
    allocate(laplacian_phi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('laplacian_phi','otfg_construct_psp')

    !** Gamma function arrays

    allocate(psi_gamma(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_gamma','otfg_construct_psp')
    allocate(phi_gamma(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi_gamma','otfg_construct_psp')
    allocate(gamma(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('beta_gamma','otfg_construct_psp')
    allocate(bl_dia_nmrshield(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_dia_nmrshield','otfg_construct_psp')
    allocate(bl_para_nmrshield(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_para_nmrshield','otfg_construct_psp')
    allocate(bl_para_damped(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_para_damped','otfg_construct_psp')
    allocate(bl_efg_aug(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_efg_aug','otfg_construct_psp')
    allocate(bl_hf_dipolar(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_dipolar','otfg_construct_psp')
    allocate(bl_hf_dipolar_cutoff(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_dipolar_cutoff','otfg_construct_psp')
    allocate(bl_hf_zora(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_zora','otfg_construct_psp')
    allocate(bl_hf_fermi(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_fermi','otfg_construct_psp')
    allocate(bl_zke_aug_gtensor(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_zke_aug_gtensor','otfg_construct_psp')
    allocate(bl_sodia_g_tensor(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_sodia_g_tensor','otfg_construct_psp')
    allocate(bl_sopara_g_tensor(ps%num_gamma,ps%num_gamma),stat=ierr)
    allocate(V_ae(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_ae','otfg_construct_psp')
    allocate(V_ps(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_ps','otfg_construct_psp')

    ! ** augmentation terms for ELNES calculation
    allocate(bl_elnes_aug(ps%num_core_radials,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_elnes_aug','otfg_construct_psp')
    allocate(phi_gamma_temp(ab%npts,ps%num_gamma),stat=ierr) ! ELNES
    if (ierr/=0) call io_allocate_abort('phi_gamma_temp','otfg_construct_psp')
    allocate(bl_ttt_temp(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_ttt_temp','otfg_construct_psp')
    allocate(psi_gamma_temp(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_gamma_temp','otfg_construct_psp')

    ! -------------------------------
    ! ** Pseudise the local potential
    ! -------------------------------

    ! ** Generate the partial wave

    psi_loc = otfg_partial_wave(ps%local_l,ps%local_e,ps%kkbeta,ab,aeat)

    ! ** Pseudise the partial wave, and invert the Schrodinger eqn.

    phi_loc = otfg_pseudise_partial_wave(ps%pv,ps%flags,ps%local_l,ps%local_norm,psi_loc,ps%local_nrc,ps%local_e,ab,aeat,V_loc,&
         genloc='genloc',scheme=ps%local_scheme,vtau=V_tau)

    ! --------------------------
    ! ** First the partial waves
    ! --------------------------

    beta_loop: do na=1,ps%num_beta

       ! ** Calculate the partial wave

       psi(:,na) = otfg_partial_wave(ps%beta_l(na),ps%beta_e(na),ps%kkbeta,ab,aeat)

       ! ** Pseudise the partial wave, and construct the beta-functions analytically

       phi(:,na) = otfg_pseudise_partial_wave(.false.,ps%flags,ps%beta_l(na),ps%beta_norm(na),psi(:,na),&
            ps%beta_nrc(na),ps%beta_e(na),ab,aeat,beta(:,na),vloc=V_loc,scheme=ps%beta_scheme(na),vtau=V_tau)

       ! ** Construct (T+Vae)|psi> = e|psi>

       tvpsi(:,na) = ps%beta_e(na)*psi(:,na)

       ! ** Construct (T+Vps)|psi> = -|beta> + e|phi>

       tvphi(:,na) = -beta(:,na) + ps%beta_e(na)*phi(:,na)


       if(testing) then

          ! ** Calculate the partial wave

          psi_shift(:,na) = otfg_partial_wave(ps%beta_l(na),ps%beta_e(na)+eshift,ps%kkbeta,ab,aeat)

          ! ** Pseudise the partial wave

          phi_shift(:,na) = otfg_pseudise_partial_wave(.false.,ps%flags,ps%beta_l(na),ps%beta_norm(na),psi_shift(:,na),&
               ps%beta_nrc(na),ps%beta_e(na)+eshift,ab,aeat,work,vloc=V_loc,scheme=ps%beta_scheme(na),vtau=V_tau)

       end if

    end do beta_loop

    ! ===========================
    ! ** Output the partial waves
    ! ===========================

    call otfg_write_pwave(ab,ps,psi,phi,psi_loc,phi_loc)

    ! ======================================================
    ! ** Construct the beta functions using Bloechl's scheme
    ! ** Eqns. 90 to 96 PRB 50, 17953, 1994
    ! ======================================================

    ! ** The Beta projectors

    do na = 1,ps%num_beta

       ! ------------------------
       ! ** The |beta> projectors
       ! ------------------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb)) then

             work(:) = work(:) + beta(:,nb)*otfg_radin(ps%kkbeta,ab%rab,phi(:,nb)*beta(:,na))

          end if

       end do

       beta(:,na) = beta(:,na) - work(:)

       ! ---------------
       ! ** (T+Vae)|psi>
       ! ---------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb)) then

             work(:) = work(:) + tvpsi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       tvpsi(:,na) = tvpsi(:,na) - work(:)

       ! ---------------
       ! ** (T+Vps)|phi>
       ! ---------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb)) then

             work(:) = work(:) + tvphi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       tvphi(:,na) = tvphi(:,na) - work(:)

       ! -----------------------------
       ! ** The AE partial waves |psi>
       ! -----------------------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb)) then

             work(:) = work(:) + psi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       psi(:,na) = psi(:,na) - work(:)

       ! -----------------------------
       ! ** The PS partial waves |phi>
       ! -----------------------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb)) then

             work(:) = work(:) + phi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       phi(:,na) = phi(:,na) - work(:)

       ! ------------------------------------
       ! ** Scale beta so that <beta|phi> = 1
       ! ------------------------------------

       beta(:,na) = beta(:,na)/otfg_radin(ps%kkbeta,ab%rab,beta(:,na)*phi(:,na))

       ! ------------------------------
       ! ** Balance by normalising beta
       ! ------------------------------

       factor = 1.0_dp/sqrt(otfg_radin(ps%kkbeta,ab%rab,beta(:,na)*beta(:,na)))

       beta(:,na) = factor*beta(:,na)

       psi(:,na)   = psi(:,na)/factor
       phi(:,na)   = phi(:,na)/factor
       tvpsi(:,na) = tvpsi(:,na)/factor
       tvphi(:,na) = tvphi(:,na)/factor

    end do


    if(testing) then

       phivphi        = 0.0_dp
       phivphi_shift  = 0.0_dp
       phivphi_shift2 = 0.0_dp
       betaphishift   = 0.0_dp

       do na=1,ps%num_beta
          do nb=1,ps%num_beta

             if(ps%beta_l(na).eq.ps%beta_l(nb)) betaphishift(na,nb) = otfg_radin(ps%kkbeta,ab%rab,beta(:,na)*phi_shift(:,nb))

             work(1) = 0.0_dp
             work(2:) = phi(2:,na)*(aeat%V_eff(2:)-V_loc(2:))*phi(2:,nb)/ab%r(2:)

             if(ps%beta_l(na).eq.ps%beta_l(nb)) phivphi(na,nb) = otfg_radin(ps%kkbeta,ab%rab,work)

             work(1) = 0.0_dp
             work(2:) = phi_shift(2:,na)*(aeat%V_eff(2:)-V_loc(2:))*phi_shift(2:,nb)/ab%r(2:)

             if(ps%beta_l(na).eq.ps%beta_l(nb)) phivphi_shift(na,nb) = otfg_radin(ps%kkbeta,ab%rab,work)

          end do
       end do

       do na=1,ps%num_beta
          do nb=1,ps%num_beta

             do naa=1,ps%num_beta
                do nbb=1,ps%num_beta

                   phivphi_shift2(na,nb) = phivphi_shift2(na,nb) + betaphishift(naa,na)*phivphi(naa,nbb)*betaphishift(nbb,nb)

                end do
             end do

          end do
       end do

       if(on_root) then
          write (stdout,'(a,ES8.2)') ' Identity difference: ',&
               sqrt(sum((phivphi_shift2-phivphi_shift)**2))/abs(eshift)/ps%num_beta**2
          write (stdout,*)
       end if

    end if

    ! ----------------------------------------------
    ! ** Construct the charge augmentation functions
    !     -- and the integrated coefficients, qqq
    ! ----------------------------------------------

    qfunc = 0.0_dp
    qqq   = 0.0_dp

    if(ps%aug) then
       do na=1,ps%num_beta
          ! for the kinetic energy augmentation terms, we need the laplacian of psi and phi:
          ! \nabla^2 \psi(\mathbf{r})
          ! we actually have r*psi(r) and want to end up with r*(radial part of laplacian)
          call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,psi(:,na),work2,laplacian_psi(:,na))
          laplacian_psi(2:ab%npts,na) = laplacian_psi(2:ab%npts,na) - &
               ps%beta_l(na)*(ps%beta_l(na)+1)*psi(2:ab%npts,na) / ab%r(2:ab%npts)**2
          laplacian_psi(1,na) = 0.0_dp

          call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,phi(:,na),work2,laplacian_phi(:,na))
          laplacian_phi(2:ab%npts,na) = laplacian_phi(2:ab%npts,na) - &
               ps%beta_l(na)*(ps%beta_l(na)+1)*phi(2:ab%npts,na) / ab%r(2:ab%npts)**2
          laplacian_phi(1,na) = 0.0_dp

          do nb=1,ps%num_beta
             qfunc(:,na,nb) = psi(:,na)*psi(:,nb)-phi(:,na)*phi(:,nb)
             if(ps%beta_l(na).eq.ps%beta_l(nb)) then
                qqq(na,nb) = otfg_radin(ps%kkbeta,ab%rab,qfunc(:,na,nb))
             else
                qqq(na,nb) = 0.0_dp
             end if
          end do
       end do
    end if

    ! ** Set the ddd matrix: <psi|T+Vae|psi> - <phi|T+Vps|phi>

    ddd = 0.0_dp
    do na=1,ps%num_beta
       do nb=1,ps%num_beta
          if(ps%beta_l(na).eq.ps%beta_l(nb)) ddd(na,nb) = &
               otfg_radin(ps%kkbeta,ab%rab,psi(:,na)*tvpsi(:,nb))-otfg_radin(ps%kkbeta,ab%rab,phi(:,na)*tvphi(:,nb))
!!$          if(na==nb) ddd(na,nb) = ddd(na,nb) + ps%beta_delta(na)*otfg_radin(ab%npts,ab%rab,phi(:,na)*phi(:,nb))
          if(na==nb) ps%beta_delta(na) = ps%beta_delta(na)*otfg_radin(ab%npts,ab%rab,phi(:,na)*phi(:,nb))
       end do
    end do

    ! ** Symmetrise the ddd matrix

    do na=1,ps%num_beta
       do nb=na,ps%num_beta
          ddd(na,nb) = (ddd(na,nb)+ddd(nb,na))/2.0_dp
          ddd(nb,na) = ddd(na,nb)
       end do
    end do

    ! ** Pseudise the Q-functions

    qfcoef = 0.0_dp
    tfcoef = 0.0_dp
    tfunc = 0.0_dp

    if(ps%aug) then

       do na=1,ps%num_beta
          do nb=na,ps%num_beta
             work2=qfunc(:,na,nb)
             llmin = abs(ps%beta_l(na)-ps%beta_l(nb))
             llmax = abs(ps%beta_l(na)+ps%beta_l(nb))

             ! we need \nabla^2 Q(\mathbf{r}) = \frac{1}{r^2}
             ! ( \frac{\partial}{\partial r} (r^2 \frac{\partial}{\partial r}) - \hat{l}^2) Q(\mathbf{r}) =
             ! ( \frac{1}{r} \frac{\partial^2}{\partial r^2} r - \frac{1}{r^2}\hat{l}^2)  Q(r) Y_{lm}(\hat{r})
             ! prepare second derivatives. We have r^2 Q(r), so get r Q(r) instead and differentiate it twice
             work(1) = 0.0_dp
             work(2:ab%npts) = qfunc(2:ab%npts,na,nb) / ab%r(2:ab%npts)
             call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,work,work2,laplacian_qfunc1)
             !
             ! multiply with r, this provides r^2 times the first part of the radial dependent term
             laplacian_qfunc1 = laplacian_qfunc1 * ab%r

             ! We need to divide r^2 Q(r) by r^2 to get the second part of the radial dependent term, which was
             ! otherwise left intact by the angular momentum operator
             laplacian_qfunc2(2:ab%npts) = qfunc(2:ab%npts,na,nb) / ab%r(2:ab%npts)**2
             ! extrapolate back to first point (r=0)
             laplacian_qfunc2(1) = laplacian_qfunc2(2) + ( laplacian_qfunc2(3)-laplacian_qfunc2(2) ) / &
                  ( ab%r(3) - ab%r(2) ) * ( ab%r(1) - ab%r(2) )
             !laplacian_qfunc2(1)  = laplacian_qfunc2(3) - (laplacian_qfunc2(4)-laplacian_qfunc2(2) ) / &
             !   2.0_dp / ab%rab(3)*(ab%r(3)-ab%r(1))

             ! in the first term (that doesn't depend on the angular channel L) accumulate the
             ! terms like (laplacian psi)*psi etc
             laplacian_qfunc1 = laplacian_qfunc1 - &
                  ( laplacian_psi(:,na)*psi(:,nb) + laplacian_psi(:,nb)*psi(:,na) - &
                  laplacian_phi(:,na)*phi(:,nb) - laplacian_phi(:,nb)*phi(:,na) )

             work2 = qfunc(:,na,nb)
             do l=llmin,llmax
                if(mod(l-llmin,2)/=1) then
                   ! finalise the kinetic energy density with the L-dependent term
                   tfunc(:,l,na,nb) = laplacian_qfunc1 - l*(l+1)*laplacian_qfunc2
                   tfunc(:,l,na,nb) = 0.25_dp * tfunc(:,l,na,nb) ! Definition with half
                   if(l.eq.llmin) then
                      if(compat7) then
                         qfunc(:,na,nb)=otfg_Q_pseudise_compat7(l,work2,ps,ab,qfcoef(:,l,na,nb))
                         tfunc(:,l,na,nb) = otfg_Q_pseudise_compat7(l,tfunc(:,l,na,nb),ps,ab,tfcoef(:,l,na,nb))
                      else
                         qfunc(:,na,nb)=otfg_Q_pseudise(l,work2,ps,ab,qfcoef(:,l,na,nb))
                         tfunc(:,l,na,nb) = otfg_T_pseudise(l,tfunc(:,l,na,nb),ps,ab,tfcoef(:,l,na,nb))
                      end if
                      qfunc(:,nb,na)=qfunc(:,na,nb)
                      qfcoef(:,l,nb,na)=qfcoef(:,l,na,nb)
                      tfcoef(:,l,nb,na)=tfcoef(:,l,na,nb)
                      tfunc(:,l,nb,na)=tfunc(:,l,na,nb)
                   else
                      if(compat7) then
                         work = otfg_Q_pseudise_compat7(l,work2,ps,ab,qfcoef(:,l,na,nb))
                         tfunc(:,l,na,nb) = otfg_Q_pseudise_compat7(l,tfunc(:,l,na,nb),ps,ab,tfcoef(:,l,na,nb))
                      else
                         work = otfg_Q_pseudise(l,work2,ps,ab,qfcoef(:,l,na,nb))
                         tfunc(:,l,na,nb) = otfg_T_pseudise(l,tfunc(:,l,na,nb),ps,ab,tfcoef(:,l,na,nb))
                      end if
                      qfcoef(:,l,nb,na)=qfcoef(:,l,na,nb)
                      tfcoef(:,l,nb,na)=tfcoef(:,l,na,nb)
                      tfunc(:,l,nb,na)=tfunc(:,l,na,nb)
                   end if
                end if
             end do
          end do
       end do

    end if

    ! ============================================
    ! ** Output the beta functions for examination
    ! ============================================

    call otfg_write_beta(ab,ps,beta,V_loc,aeat%V_eff)

    ! =======================================================
    ! ** Convert into the correct form: reciprocal space etc.
    ! =======================================================

    ! ** Set up indexing

    ps%num_ps_projectors = 0
    do n=1,ps%num_beta
       ps%num_ps_projectors = ps%num_ps_projectors+2*ps%beta_l(n)+1 ! 2l+1
       ps%projector_l(ps%num_ps_projectors-2*ps%beta_l(n):ps%num_ps_projectors) = ps%beta_l(n)
       ps%projector_m(ps%num_ps_projectors-2*ps%beta_l(n):ps%num_ps_projectors) = &
            (/(m,m=-ps%beta_l(n),ps%beta_l(n))/)
    end do

    allocate(nnn(ps%num_ps_projectors),stat=ierr)
    if (ierr/=0) call io_allocate_abort('nnn','otfg_construct_psp')

    m = 0
    do n=1,ps%num_beta
       m = m+2*ps%beta_l(n)+1 ! 2l+1
       nnn(m-2*ps%beta_l(n):m) = n
    end do

    ! ** The local potential

    ps%V_locr(:ps%mesh) = V_loc(:ps%mesh)
    ps%has_V_tau = xc_is_mgga(aeat%theory)
    if( ps%has_V_tau ) then
!!$       call otfg_pspot_allocate_v_tau(ps)
       ps%V_taur(:ps%mesh) = V_tau(:ps%mesh)
    endif

    ! ** The beta functions

    do n=1,ps%num_ps_projectors
       call otfg_radial_transform(ps%beta_l(nnn(n)),maxval(ps%beta_l),1,ps%kkbeta,&
            ab%r,ab%rab,ps%num_points,ps%gmax,beta(:,nnn(n)),ps%beta(:,n))

       ps%betar(1:ab%npts,n) = beta(1:ab%npts,nnn(n))
    end do
    ps%beta = ps%beta*4*pi

    ! ** The projector weights

    ps%q(:,:)  = 0.0_dp
    ps%D0(:,:) = 0.0_dp

    do n=1,ps%num_ps_projectors
       do m=1,ps%num_ps_projectors
          if((ps%projector_l(n)==ps%projector_l(m)).and.(ps%projector_m(n)==ps%projector_m(m))) then
             ps%q(n,m)  = qqq(nnn(n),nnn(m))
             ps%D0(n,m) = ddd(nnn(n),nnn(m)) !* Not D0 yet -- needs unscreening first
          else
             ps%q(n,m)  = 0.0_dp
             ps%D0(n,m) = 0.0_dp
          end if
       end do
    end do

    ! ** The augmentation functions, and grids

    do n=1,ps%num_ps_projectors
       do m=1,ps%num_ps_projectors
          ps%qfunc(1:ps%kkbeta,n,m)           = qfunc(1:ps%kkbeta,nnn(n),nnn(m))
          ps%qfcoef(1:ps%nqf,0:2*ab%lmax,n,m) = qfcoef(1:ps%nqf,0:2*ab%lmax,nnn(n),nnn(m))
          ps%tfunc(1:ps%kkbeta,0:2*ab%lmax,n,m) = tfunc(1:ps%kkbeta,0:2*ab%lmax,nnn(n),nnn(m))
          ps%tfcoef(1:ps%nqf,0:2*ab%lmax,n,m) = tfcoef(1:ps%nqf,0:2*ab%lmax,nnn(n),nnn(m))
       end do
    end do

    !************************************
    !** Now the Blochl GAMMA functions **
    !************************************

    ! Set the descreened potentials for use in the augmentations

    V_ae = aeat%V_eff
    ! reason for limit in next line is that ps%V_locr is set to maximum size across all species
    V_ps(1:ab%npts) = ps%V_locr(1:ab%npts)

    if(ps%num_beta==ps%num_gamma) then

       ! ** Fill up with the projectors already generated above

       do na=1,ps%num_beta

          psi_gamma(:,na) = psi(:,na)

          phi_gamma(:,na) = phi(:,na)

          gamma(:,na)     = beta(:,na)

       end do

    else

       ! ** Now calculate the extra ones

       gamma_loop: do na=1,ps%num_gamma

          ! ** Calculate the partial wave

          psi_gamma(:,na) = otfg_partial_wave(ps%gamma_l(na),ps%gamma_e(na),ps%kkbeta,ab,aeat)

          ! ** Pseudise the partial wave, and construct the gamma-functions analytically

          phi_gamma(:,na) = otfg_pseudise_partial_wave(.false.,ps%flags,ps%gamma_l(na),ps%gamma_norm(na),psi_gamma(:,na),&
               ps%gamma_nrc(na),ps%gamma_e(na),ab,aeat,gamma(:,na),vloc=V_loc,scheme=ps%gamma_scheme(na))

          ! ** Gamma may disappear

          if(abs(otfg_radin(ps%kkbeta,ab%rab,gamma(:,na)*gamma(:,na)))<1.0e-20_dp) then

             gamma(:,na) = ab%r**(l+1)*exp(-(ab%r/(0.7_dp*ab%r(ps%gamma_nrc(na))))**6.0_dp)

          end if

       end do gamma_loop

       ! =======================================================
       ! ** Construct the gamma functions using Bloechl's scheme
       ! ** Eqns. 90 to 96 PRB 50, 17953, 1994
       ! =======================================================

       do na = 1,ps%num_gamma

          ! ------------------------
          ! ** The |beta> projectors
          ! ------------------------

          work(:) = 0.0_dp

          do nb=1,na-1

             if(ps%gamma_l(na)==ps%gamma_l(nb)) then

                work(:) = work(:) + gamma(:,nb)*otfg_radin(ps%kkbeta,ab%rab,phi_gamma(:,nb)*gamma(:,na))

             end if

          end do

          gamma(:,na) = gamma(:,na) - work(:)

          ! -----------------------------
          ! ** The AE partial waves |psi>
          ! -----------------------------

          work(:) = 0.0_dp

          do nb=1,na-1

             if(ps%gamma_l(na)==ps%gamma_l(nb)) then

                work(:) = work(:) + psi_gamma(:,nb)*otfg_radin(ps%kkbeta,ab%rab,gamma(:,nb)*phi_gamma(:,na))

             end if

          end do

          psi_gamma(:,na) = psi_gamma(:,na) - work(:)

          ! -----------------------------
          ! ** The PS partial waves |phi>
          ! -----------------------------

          work(:) = 0.0_dp

          do nb=1,na-1

             if(ps%gamma_l(na)==ps%gamma_l(nb)) then

                work(:) = work(:) + phi_gamma(:,nb)*otfg_radin(ps%kkbeta,ab%rab,gamma(:,nb)*phi_gamma(:,na))

             end if

          end do

          phi_gamma(:,na) = phi_gamma(:,na) - work(:)

          ! ------------------------------------
          ! ** Scale beta so that <beta|phi> = 1
          ! ------------------------------------

          gamma(:,na) = gamma(:,na)/otfg_radin(ps%kkbeta,ab%rab,gamma(:,na)*phi_gamma(:,na))

          ! ------------------------------
          ! ** Balance by normalising beta
          ! ------------------------------

          factor = 1.0_dp/sqrt(otfg_radin(ps%kkbeta,ab%rab,gamma(:,na)*gamma(:,na)))

          gamma(:,na) = factor*gamma(:,na)

          psi_gamma(:,na)   = psi_gamma(:,na)/factor
          phi_gamma(:,na)   = phi_gamma(:,na)/factor

       end do

    end if

    ! ===================================
    ! ** Calculate the augmentation terms
    ! ===================================

    bl_dia_nmrshield  = 0.0_dp
    bl_para_nmrshield = 0.0_dp
    bl_para_damped    = 0.0_dp
    bl_efg_aug        = 0.0_dp
    bl_hf_dipolar     = 0.0_dp
    bl_hf_dipolar_cutoff     = 0.0_dp
    bl_hf_fermi       = 0.0_dp
    bl_hf_zora        = 0.0_dp
    bl_zke_aug_gtensor = 0.0_dp
    bl_sodia_g_tensor  = 0.0_dp
    bl_sopara_g_tensor = 0.0_dp

    do na=1,ps%num_gamma
       do nb=1,ps%num_gamma

          ! * NMR Diamagnetic Augmentation

          work=0.0_dp
          work(1)=0.0_dp
          work(2:)=( aeat%ZORA_K(2:)*(psi_gamma(2:,na)*psi_gamma(2:,nb))-(phi_gamma(2:,na)*phi_gamma(2:,nb)) ) &
               /(ab%r(2:))

          bl_dia_nmrshield(na,nb)=otfg_radin(ps%gamma_nrc(na),ab%rab,work)!*fine_structure_si**2.0_dp*ppm

          ! * Electric Field Gradient Augmentation
          ! this is the same as the hyperfine dipolar term - at least for the moment we keep them distinct

          work=0.0_dp
          work(1)=0.0_dp
          work(2:)=( (psi_gamma(2:,na)*psi_gamma(2:,nb))-(phi_gamma(2:,na)*phi_gamma(2:,nb)) ) &
               /(ab%r(2:)**3)

          bl_efg_aug(na,nb)=otfg_radin(ps%gamma_nrc(na),ab%rab,work)

          ! * Undamped spin dipolar interaction

          work=0.0_dp
          work(1)=0.0_dp
          work(2:)=( aeat%ZORA_K(2:)*(psi_gamma(2:,na)*psi_gamma(2:,nb))-(phi_gamma(2:,na)*phi_gamma(2:,nb)) ) &
               /(ab%r(2:)**3)

          bl_hf_dipolar(na,nb) = otfg_radin(ps%gamma_nrc(na),ab%rab,work)

          ! * Damped spin dipolar interaction
          work = 0.0_dp
          work(1) = 0.0_dp
          work(2:) = aeat%ZORA_K(2:)*psi_gamma(2:,na)*psi_gamma(2:,nb) / (ab%r(2:)**3) &
               - phi_gamma(2:,na)*phi_gamma(2:,nb) &
               * (1.0_dp  - exp(-(max(ab%r(2:),dtol)/damp_scale)**3)) &
               / max(ab%r(2:), dtol)**3

          bl_hf_dipolar_cutoff(na,nb) = otfg_radin(ps%gamma_nrc(na),ab%rab,work)

          ! * Extra ZORA ZSO bit
          work=0.0_dp
          work(1)=0.0_dp
          work(2:)=(aeat%ZORA_dKdr(2:)*psi_gamma(2:,na)*psi_gamma(2:,nb)) &
               /(ab%r(2:)**2)


          bl_hf_zora(na,nb) = otfg_radin(ps%gamma_nrc(na),ab%rab,work)

          ! * Hyperfine Fermi Contact Augmentation
          ! Term will only contribute for s-type partial waves

          if(ps%gamma_l(na).eq.0 .and. ps%gamma_l(nb).eq.0) then
             bl_hf_fermi(na,nb)=otfg_fermi_contact(ab%r,psi_gamma(:,na))* &
                  otfg_fermi_contact(ab%r,psi_gamma(:,nb))/(4.0_dp*pi)*aeat%ZORA_K(1)
          else
             bl_hf_fermi(na,nb)=0
          end if


          if(ps%gamma_l(na).eq.ps%gamma_l(nb)) then


             ! ** NMR Paramagnetic Augmentation
             !
             ! for l=0 the function <psi|r-3|psi> will have a divergence at the origin
             ! but as we eventually need <psi|Lr-3|psi> the contribution from l=0 is zero
             ! So it is better to set the l=0 component to zero here

             if(ps%gamma_l(na).gt.0) then
                work=0.0_dp
                work(1)=0.0_dp
                work(2:) = (aeat%ZORA_K(2:)*(psi_gamma(2:,na)*psi_gamma(2:,nb))-(phi_gamma(2:,na)*phi_gamma(2:,nb))) &
                     /(ab%r(2:)**3)

                bl_para_nmrshield(na,nb) = &
                     otfg_radin(ps%gamma_nrc(na),ab%rab,work)!*fine_structure_si**2.0_dp*ppm*2.0_dp

                work(1)=0.0_dp
                work(2:) = aeat%ZORA_K(2:)*psi_gamma(2:,na)*psi_gamma(2:,nb) / ab%r(2:)**3 &
                     - (phi_gamma(2:,na)*phi_gamma(2:,nb)) &
                     * (1.0_dp  - exp(-(max(ab%r(2:),dtol)/damp_scale)**3)) &
                     / max(ab%r(2:), dtol)**3

                bl_para_damped(na,nb) = &
                     otfg_radin(ps%gamma_nrc(na),ab%rab,work)!*fine_structure_si**2.0_dp*ppm*2.0_dp

             else

                bl_para_nmrshield(na,nb)=0.0_dp
                bl_para_damped(na,nb)=0.0_dp

             end if

             ! * EPR Zeeman kinetic-energy Correction (or PAW aug. of the kinetic energy)
             ! <Psi_na|T_l|Psi_nb>-<Phi_na|T_l|Phi_nb>, where l_na==l_nb
             ! Tl= -0.5*(d^2/dr^2 - l_n(l_n+1)/r^2)
             work=0.0_dp
             work(1) = 0.0_dp

             work(:)=psi_gamma(:,na)*otfg_apply_Tl(ps%gamma_l(nb),size(psi_gamma(:,nb)),ab%rab,ab%r,psi_gamma(:,nb))

             work(:) = work(:) - phi_gamma(:,na)* &
                  otfg_apply_Tl(ps%gamma_l(nb),size(phi_gamma(:,nb)),ab%rab,ab%r,phi_gamma(:,nb))

             bl_zke_aug_gtensor(na,nb) = otfg_radin(ps%gamma_nrc(na),ab%rab,work)

             ! * EPR spin orbit diamagnetic augmentation term
             ! new version  (bshih (16/12/2013)
             work=0.0_dp;work2=0.0_dp
             call otfg_real_derivative_kind1(ab,V_ae,work)
             call otfg_real_derivative_kind1(ab,V_ps,work2)
             work(1)=0.0_dp
             work(2:)=(psi_gamma(2:,na)*ab%r(2:)*work(2:)*psi_gamma(2:,nb)- &
                  phi_gamma(2:,na)*ab%r(2:)*work2(2:)*phi_gamma(2:,nb))

             bl_sodia_g_tensor(na,nb)=otfg_radin(ps%gamma_nrc(na),ab%rab,work)/3.0_dp
             ! bshih: why 1/3?


             ! * EPR spin orbit - paramagnetic augmentation term
             ! Add all electron part
             work=0.0_dp
             call otfg_real_derivative_kind1(ab,V_ae,work)
             work(1)=0.0_dp
             work(2:)=work(2:)*psi_gamma(2:,na)*psi_gamma(2:,nb)/ab%r(2:)
             !integrate
             bl_sopara_g_tensor(na,nb)=otfg_radin(ps%gamma_nrc(na),ab%rab,work)

             ! Subtract ps part
             work=0.0_dp
             call otfg_real_derivative_kind1(ab,V_ps,work)
             work(1)=0.0_dp
             work(2:)=work(2:)*phi_gamma(2:,na)*phi_gamma(2:,nb)/ab%r(2:)

             !integrate
             bl_sopara_g_tensor(na,nb)=bl_sopara_g_tensor(na,nb)&
                  -otfg_radin(ps%gamma_nrc(na),ab%rab,work)

          end if

       end do
    end do

    ! =============================================
    ! ** Output the gamma functions for examination
    ! =============================================

    !!if(ps%num_gamma>ps%num_beta)
    call otfg_write_gamma(ab,ps,gamma)

    ! * Copy data into pseudopotential structure

    ps%bl_dia_nmrshift  = 0.0_dp
    ps%bl_para_nmrshift = 0.0_dp
    ps%bl_para_damped   = 0.0_dp
    ps%bl_efg_aug       = 0.0_dp
    ps%bl_hf_dipolar    = 0.0_dp
    ps%bl_hf_dipolar_cutoff    = 0.0_dp
    ps%bl_hf_zora    = 0.0_dp
    ps%bl_hf_fermi      = 0.0_dp
    ps%bl_zke_aug_gtensor  = 0.0_dp
    ps%bl_sodia_g_tensor  = 0.0_dp
    ps%bl_sopara_g_tensor  = 0.0_dp

    ! ** Set up gamma indexing

    ps%num_bl_projectors = 0
    do n=1,ps%num_gamma
       ps%num_bl_projectors = ps%num_bl_projectors+2*ps%gamma_l(n)+1 ! 2l+1
       ps%bl_projector_l(ps%num_bl_projectors-2*ps%gamma_l(n):ps%num_bl_projectors) = ps%gamma_l(n)
       ps%bl_projector_m(ps%num_bl_projectors-2*ps%gamma_l(n):ps%num_bl_projectors) = &
            (/(m,m=-ps%gamma_l(n),ps%gamma_l(n))/)
    end do


    allocate(nnn_bl(ps%num_bl_projectors),stat=ierr)
    if (ierr/=0) call io_allocate_abort('nnn_bl','otfg_construct_psp')

    m = 0
    do n=1,ps%num_gamma
       m = m+2*ps%gamma_l(n)+1 ! 2l+1
       nnn_bl(m-2*ps%gamma_l(n):m) = n
    end do

    do na=1,ps%num_gamma
      do nb=1,ps%num_gamma
         qfunc_bl(:,na,nb) = psi_gamma(:,na)*psi_gamma(:,nb)-phi_gamma(:,na)*phi_gamma(:,nb)
      end do
   end do

   do n=1,ps%num_bl_projectors
      do m=1,ps%num_bl_projectors
         ps%qfunc_bl(1:ps%kkbeta,n,m) = qfunc_bl(1:ps%kkbeta,nnn_bl(n),nnn_bl(m))
      end do
   end do

    ! ** The gamma functions

    do n=1,ps%num_bl_projectors

       call otfg_radial_transform(ps%gamma_l(nnn_bl(n)),maxval(ps%gamma_l),1,maxval(ps%gamma_nrc)+2,ab%r,ab%rab,&
            ps%num_points,ps%gmax,gamma(:,nnn_bl(n)),ps%gamma(:,n))

    end do
    ps%gamma = ps%gamma*4*pi

    ! ** The projector weights

    do n=1,ps%num_bl_projectors
       do m=1,ps%num_bl_projectors
          ps%bl_dia_nmrshift(n,m)  = bl_dia_nmrshield(nnn_bl(n),nnn_bl(m))
          ps%bl_para_nmrshift(n,m) = bl_para_nmrshield(nnn_bl(n),nnn_bl(m))
          ps%bl_para_damped(n,m) = bl_para_damped(nnn_bl(n),nnn_bl(m))
          ps%bl_efg_aug(n,m)       = bl_efg_aug(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_dipolar(n,m)    = bl_hf_dipolar(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_dipolar_cutoff(n,m)    = bl_hf_dipolar_cutoff(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_zora(n,m)    = bl_hf_zora(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_fermi(n,m)      = bl_hf_fermi(nnn_bl(n),nnn_bl(m))
          ps%bl_zke_aug_gtensor(n,m)  = bl_zke_aug_gtensor(nnn_bl(n),nnn_bl(m))
          ps%bl_sodia_g_tensor(n,m)  = bl_sodia_g_tensor(nnn_bl(n),nnn_bl(m))
          ps%bl_sopara_g_tensor(n,m) = bl_sopara_g_tensor(nnn_bl(n),nnn_bl(m))
       end do
    end do

    !------------------------------------------------------------------
    ! ** Generate the terms for core excited matrix elements ---ELNES
    !------------------------------------------------------------------
    allocate(core_radials_r(ab%npts,1:15),stat=ierr)
    if (ierr/=0) call io_abort('Error in allocating core_radials_r in otfg_construct_psp')

    ! ** The core radial orbitals

    nl=0
    do nl=1,ps%num_core_radials
       nn=ps%core_radials_n(nl)
       ll=ps%core_radials_l(nl)

       ! * real space
       core_radials_r(:,nl)=aeat%orbr(:,nn,ll)

       ! * reciprocal space
       select case(ll)
       case(0)  ! for s*x/r,s*y/r,s*z/r
          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),1)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,2)=otfg_radin(ab%npts,ab%rab,work)
          end do

       case(1) ! for p*x/r, p*y/r, p*z/r
          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),0)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,1)=otfg_radin(ab%npts,ab%rab,work)
          end do

          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),2)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,3)=otfg_radin(ab%npts,ab%rab,work)
          end do

       case(2) ! for d*x/r, d*y/r, d*z/r

          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),1)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,2)=otfg_radin(ab%npts,ab%rab,work)
          end do

          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),3)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,4)=otfg_radin(ab%npts,ab%rab,work)
          end do

       case(3) ! for f*x/r, f*y/r, f*z/r
          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),0)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,1)=otfg_radin(ab%npts,ab%rab,work)
          end do

          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),1)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,2)=otfg_radin(ab%npts,ab%rab,work)
          end do

          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),2)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,3)=otfg_radin(ab%npts,ab%rab,work)
          end do

          do ng=1,ps%num_points
             g=real(ng-1,dp)/real(ps%num_points-1,dp)*ps%gmax
             do nr=1,ab%npts
                work(nr)=otfg_bessfn(g*ab%r(nr),3)*aeat%orbr(nr,nn,ll)*ab%r(nr)**2 !core_radials*r
             end do
             ps%core_radials(ng,nl,4)=otfg_radin(ab%npts,ab%rab,work)
          end do

       case default
          call io_abort('currently max ll is 3 in otfg')
       end select

    end do

    ! ** ELNES augmentation

    do na=1,ps%num_core_radials
       do nb=1,ps%num_gamma
          work=0.0_dp
          work(1)=0.0_dp
          work(2:)=((core_radials_r(2:,na)*psi_gamma(2:,nb))-(core_radials_r(2:,na)*phi_gamma(2:,nb)))*ab%r(2:)
          bl_elnes_aug(na,nb)=otfg_radin(ps%gamma_nrc(nb),ab%rab,work)
       end do
    end do

    !bl_elnes_aug        = matmul(bl_elnes_aug,transpose(bl_sss))

    ! Set up indexing
    allocate(nnn_nl(ps%num_core_projectors),stat=ierr)
    if(ierr/=0) call io_abort('Error in allocating nnn_nl in otfg_construct_psp')

    m=0
    do nl=1, ps%num_core_radials
       m=m+2*ps%core_radials_l(nl)+1
       nnn_nl(m-2*ps%core_radials_l(nl):m)=nl
    end do

    do n=1,ps%num_core_projectors
       do m=1,ps%num_bl_projectors
          ps%bl_elnes_aug(n,m)=bl_elnes_aug(nnn_nl(n),nnn_bl(m))
       end do
    end do

    ! ** Deallocate the arrays

    deallocate(psi_loc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_loc in otfg_construct_psp')
    deallocate(phi_loc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi_loc in otfg_construct_psp')
    deallocate(V_loc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_loc in otfg_construct_psp_paw')
    deallocate(V_tau,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_tau in otfg_construct_psp_paw')


    !deallocate(zora_K,stat=ierr)
    !if (ierr/=0) call io_abort('Error in deallocating zora_K in otfg_construct_psp')
    !deallocate(zora_K_dr,stat=ierr)
    !if (ierr/=0) call io_abort('Error in deallocating zora_K_dr in otfg_construct_psp')

    deallocate(psi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi in otfg_construct_psp')
    deallocate(phi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi in otfg_construct_psp')
    deallocate(tvpsi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating tvpsi in otfg_construct_psp')
    deallocate(tvphi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating tvphi in otfg_construct_psp')
    deallocate(psi_shift,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_shift in otfg_construct_psp')
    deallocate(phi_shift,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi_shift in otfg_construct_psp')
    deallocate(beta,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating beta in otfg_construct_psp')
    deallocate(qqq,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating qqq in otfg_construct_psp')
    deallocate(ddd,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating ddd in otfg_construct_psp')
    deallocate(qfunc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating qfunc in otfg_construct_psp_paw')
    deallocate(qfunc_bl,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating qfunc_bl in otfg_construct_psp_paw')
    deallocate(tfunc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating tfunc in otfg_construct_psp_paw')
    deallocate(qfcoef,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating qfcoef in otfg_construct_psp_paw')
    deallocate(tfcoef,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating tfcoef in otfg_construct_psp_paw')
    deallocate(nnn,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating nnn in otfg_construct_psp')
    deallocate(work,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating work in otfg_construct_psp')
    deallocate(work2,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating work2 in otfg_construct_psp_paw')
    deallocate(laplacian_qfunc1,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating laplacian_qfunc1 in otfg_construct_psp_paw')
    deallocate(laplacian_qfunc2,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating laplacian_qfunc2 in otfg_construct_psp_paw')
    deallocate(laplacian_psi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating laplacian_psi in otfg_construct_psp_paw')
    deallocate(laplacian_phi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating laplacian_phi in otfg_construct_psp_paw')

    ! ** The GAMMA Arrays

    deallocate(psi_gamma,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_gamma in otfg_construct_psp')
    deallocate(phi_gamma,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi_gamma in otfg_construct_psp')
    deallocate(gamma,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating beta_gamma in otfg_construct_psp')
    deallocate(nnn_bl,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating nnn_bl in otfg_construct_psp')
    deallocate(bl_dia_nmrshield,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_dia_nmrshield in otfg_construct_psp')
    deallocate(bl_para_nmrshield,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_para_nmrshield in otfg_construct_psp')
    deallocate(bl_para_damped,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_para_damped in otfg_construct_psp')
    deallocate(bl_efg_aug,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_efg_aug in otfg_construct_psp')
    deallocate(bl_hf_fermi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_hf_fermi in otfg_construct_psp')
    deallocate(bl_hf_dipolar,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_hf_dipolar in otfg_construct_psp')
    deallocate(bl_hf_dipolar_cutoff,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_hf_dipolar_cutoff in otfg_construct_psp')
    deallocate(bl_zke_aug_gtensor,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_zke_aug_gtensor in otfg_construct_psp')
    deallocate(bl_sodia_g_tensor,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_sodia_g_tensor in otfg_construct_psp')
    deallocate(bl_sopara_g_tensor,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_sopara_g_tensor in otfg_construct_psp')
    deallocate(V_ps,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_ps in otfg_construct_psp')
    deallocate(V_ae,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_ae in otfg_construct_psp')

    ! core_radials_r for ELNES
    deallocate(core_radials_r,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating core_radials_r in otfg_construct_psp')
    deallocate(phi_gamma_temp,stat=ierr) ! ELNES
    if (ierr/=0) call io_abort('Error in deallocating phi_gamma_temp in otfg_construct_psp')
    deallocate(bl_ttt_temp,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_ttt_temp in otfg_construct_psp')
    deallocate(psi_gamma_temp,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_gamma_temp in otfg_construct_psp')

    call trace_exit('otfg_construct_psp',status)

    return
  end subroutine otfg_construct_psp


  subroutine otfg_construct_psp_dir(ab,aeat,ps)

    use io,         only : io_abort,io_allocate_abort,stdout,io_open_file,&
      io_code_present,io_code_string,io_code_logical,io_code_real
    use comms,      only : my_node_id,root_node_id
    use constants,  only : periodic_table_name
    use parameters, only : current_params

    !=========================================================================!
    ! Construct the pseudopotential                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) ab      : input : The atomic basis                                   !
    ! 2) aeat    : input : The initialised atom                               !
    ! 3) ps      : inout : The pseudopotential                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) ae_atom                                                              !
    ! 3) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) comms                                                                !
    ! 3) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 1) psi_loc : The local partial wave                                     !
    ! 2) phi_loc : The pseudised local partial wave                           !
    ! 3) psi     : The projector partial waves                                !
    ! 4) phi     : The pseudised projector partial waves                      !
    ! 5) beta    : The beta projectors                                        !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 30th June 2002                                  !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    type(ae_basis), intent(in)    :: ab
    type(ae_atom),  intent(in)    :: aeat
    type(pspot),    intent(inout) :: ps

    !-------------------------------------------------------------------------!

    integer :: llmin,llmax,ierr
    integer :: na,nb,l,n,m,naa,nbb ! * Counters

    real(kind=dp) :: factor

    real(kind=dp), dimension(:),       allocatable :: psi_loc,phi_loc,V_loc,V_tau
    real(kind=dp), dimension(:,:),     allocatable :: psi,phi,beta,qqq,ddd,tvpsi,tvphi,psi_shift,phi_shift
    real(kind=dp), dimension(:,:,:),   allocatable :: qfunc
    real(kind=dp), dimension(:,:,:,:), allocatable :: qfcoef

    integer,       dimension(:),       allocatable :: nnn

    real(kind=dp), dimension(:),       allocatable :: work
    !real(kind=dp), dimension(:),     allocatable :: zora_K, zora_K_dr

    real(kind=dp), parameter :: ppm = 1.0e6_dp

    ! ** Gamma function arrays

    real(kind=dp), dimension(:,:),     allocatable :: psi_gamma,phi_gamma,gamma

    real(kind=dp), dimension(:,:),     allocatable :: bl_dia_nmrshield
    real(kind=dp), dimension(:,:),     allocatable :: bl_para_nmrshield
    real(kind=dp), dimension(:,:),     allocatable :: bl_para_damped
    real(kind=dp), dimension(:,:),     allocatable :: bl_efg_aug
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_fermi
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_dipolar
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_dipolar_cutoff
    real(kind=dp), dimension(:,:),     allocatable :: bl_hf_zora
    real(kind=dp), dimension(:,:),     allocatable :: bl_zke_aug_gtensor
    real(kind=dp), dimension(:,:),     allocatable :: bl_sodia_g_tensor
    real(kind=dp), dimension(:,:),     allocatable :: bl_sopara_g_tensor
    real(kind=dp), dimension(:),       allocatable :: V_ae,V_ps

    integer,       dimension(:),       allocatable :: nnn_bl

    integer                                        :: nn,nl,ng,nr,ll  ! For ELNES
    integer,       dimension(:),       allocatable :: nnn_nl        ! For ELNES
    real(kind=dp)                                  :: g             ! For ELNES
    real(kind=dp), dimension(:,:),     allocatable :: core_radials_r ! For ELNES
    real(kind=dp), dimension(:,:),     allocatable :: bl_elnes_aug,bl_ttt_temp  ! For ELNES
    real(kind=dp), dimension(:,:),     allocatable :: phi_gamma_temp,psi_gamma_temp ! For ELNES
    real(kind=dp) :: factor1,factor2   ! For ELNES

    real(kind=dp) :: betaphishift(1:ps%num_beta,1:ps%num_beta),phivphi(1:ps%num_beta,1:ps%num_beta),&
         phivphi_shift(1:ps%num_beta,1:ps%num_beta),phivphi_shift2(1:ps%num_beta,1:ps%num_beta)

    real(kind=dp), parameter :: eshift = -0.01_dp

    integer, parameter                   :: nflags = 5
    character(len=80), dimension(nflags) :: flags
    integer                              :: unit_zorafile

    character(len=1) :: cang

    logical :: on_root

    real(kind=dp) :: dtol = 1.0e-3_dp
    real(kind=dp) :: damp_scale = 0.5_dp

    integer                         :: status

    call trace_entry('otfg_construct_psp_dir',status)

    on_root = my_node_id==root_node_id

    if(io_code_present(current_params%devel_code,'dtol','JC')) then
      dtol = io_code_real(current_params%devel_code,'dtol','JC')
    end if

    if(io_code_present(current_params%devel_code,'damp_scale','JC')) then
      damp_scale = io_code_real(current_params%devel_code,'damp_scale','JC')
    end if

    ! ** Allocate the local arrays

    allocate(psi_loc(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_loc','otfg_construct_psp_dir')
    allocate(phi_loc(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi_loc','otfg_construct_psp_dir')
    allocate(V_loc(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_loc','otfg_construct_psp_dir')
    allocate(V_tau(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_tau','otfg_construct_psp')
    V_tau = 0.0_dp
    !allocate(zora_K(ab%npts),stat=ierr)
    !if (ierr/=0) call io_allocate_abort('zora_K','otfg_construct_psp_dir')
    !allocate(zora_K_dr(ab%npts),stat=ierr)
    !if (ierr/=0) call io_allocate_abort('zora_K_dr','otfg_construct_psp_dir')

    allocate(psi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi','otfg_construct_psp_dir')
    allocate(phi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi','otfg_construct_psp_dir')
    allocate(tvpsi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('tvpsi','otfg_construct_psp_dir')
    allocate(tvphi(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('tvphi','otfg_construct_psp_dir')
    allocate(psi_shift(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_shift','otfg_construct_psp')
    allocate(phi_shift(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi_shift','otfg_construct_psp')
    allocate(beta(ab%npts,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('beta','otfg_construct_psp_dir')
    allocate(qqq(ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('qqq','otfg_construct_psp_dir')
    allocate(ddd(ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ddd','otfg_construct_psp_dir')
    allocate(qfunc(ab%npts,ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('qfunc','otfg_construct_psp_dir')
    allocate(qfcoef(ps%nqf,0:2*ab%lmax,ps%num_beta,ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('qfcoef','otfg_construct_psp_dir')
    allocate(work(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('work','otfg_construct_psp_dir')


    !** Gamma function arrays

    allocate(psi_gamma(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_gamma','otfg_construct_psp_dir')
    allocate(phi_gamma(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('phi_gamma','otfg_construct_psp_dir')
    allocate(gamma(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('beta_gamma','otfg_construct_psp_dir')
    allocate(bl_dia_nmrshield(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_dia_nmrshield','otfg_construct_psp_dir')
    allocate(bl_para_nmrshield(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_para_nmrshield','otfg_construct_psp_dir')
    allocate(bl_para_damped(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_para_damped','otfg_construct_psp_dir')
    allocate(bl_efg_aug(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_efg_aug','otfg_construct_psp_dir')
    allocate(bl_hf_dipolar(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_dipolar','otfg_construct_psp_dir')
    allocate(bl_hf_dipolar_cutoff(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_dipolar_cutoff','otfg_construct_psp_dir')
    allocate(bl_hf_zora(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_zora','otfg_construct_psp_dir')
    allocate(bl_hf_fermi(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_hf_fermi','otfg_construct_psp_dir')
    allocate(bl_zke_aug_gtensor(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_zke_aug_gtensor','otfg_construct_psp_dir')
    allocate(bl_sodia_g_tensor(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_sodia_g_tensor','otfg_construct_psp_dir')
    allocate(bl_sopara_g_tensor(ps%num_gamma,ps%num_gamma),stat=ierr)
    allocate(V_ae(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_ae','otfg_construct_psp_dir')
    allocate(V_ps(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('V_ps','otfg_construct_psp_dir')

    ! ** augmentation terms for ELNES calculation
    allocate(bl_elnes_aug(ps%num_core_radials,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_elnes_aug','otfg_construct_psp_dir')
    allocate(phi_gamma_temp(ab%npts,ps%num_gamma),stat=ierr) ! ELNES
    if (ierr/=0) call io_allocate_abort('phi_gamma_temp','otfg_construct_psp_dir')
    allocate(bl_ttt_temp(ps%num_gamma,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bl_ttt_temp','otfg_construct_psp_dir')
    allocate(psi_gamma_temp(ab%npts,ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('psi_gamma_temp','otfg_construct_psp_dir')

    ! -------------------------------
    ! ** Pseudise the local potential
    ! -------------------------------

    ! ** Generate the partial wave

    psi_loc = otfg_partial_wave(ps%local_l,ps%local_e,ps%kkbeta,ab,aeat)

    ! ** Pseudise the partial wave, and invert the Schrodinger eqn.

    phi_loc = otfg_pseudise_partial_wave(ps%pv,ps%flags,ps%local_l,ps%local_norm,psi_loc,ps%local_nrc,&
         ps%local_e,ab,aeat,V_loc,genloc='genloc',scheme=ps%local_scheme,vtau=V_tau)

    ! --------------------------
    ! ** First the partial waves
    ! --------------------------

    beta_loop: do na=1,ps%num_beta

       ! ** Calculate the partial wave
       psi(:,na) = otfg_partial_wave(ps%beta_l(na),ps%beta_e(na),ps%kkbeta,ab,aeat,ps%beta_2j(na))

       ! ** Pseudise the partial wave, and construct the beta-functions analytically
       phi(:,na) = otfg_pseudise_partial_wave(.false.,ps%flags,ps%beta_l(na),ps%beta_norm(na),psi(:,na),&
            ps%beta_nrc(na),ps%beta_e(na),ab,aeat,beta(:,na),vloc=V_loc,scheme=ps%beta_scheme(na),vtau=V_tau)

       ! ** Construct (T+Vae)|psi> = e|psi>

       tvpsi(:,na) = ps%beta_e(na)*psi(:,na)

       ! ** Construct (T+Vps)|psi> = -|beta> + e|phi>

       tvphi(:,na) = -beta(:,na) + ps%beta_e(na)*phi(:,na)


       if(testing) then

          ! ** Calculate the partial wave

          psi_shift(:,na) = otfg_partial_wave(ps%beta_l(na),ps%beta_e(na)+eshift,ps%kkbeta,ab,aeat)

          ! ** Pseudise the partial wave

          phi_shift(:,na) = otfg_pseudise_partial_wave(.false.,ps%flags,ps%beta_l(na),ps%beta_norm(na),psi_shift(:,na),&
               ps%beta_nrc(na),ps%beta_e(na)+eshift,ab,aeat,work,vloc=V_loc,scheme=ps%beta_scheme(na),vtau=V_tau)

       end if

    end do beta_loop

    ! ===========================
    ! ** Output the partial waves
    ! ===========================

    call otfg_write_pwave(ab,ps,psi,phi,psi_loc,phi_loc)

    ! ======================================================
    ! ** Construct the beta functions using Bloechl's scheme
    ! ** Eqns. 90 to 96 PRB 50, 17953, 1994
    ! ======================================================

    ! ** The Beta projectors

    do na = 1,ps%num_beta
       ! ------------------------
       ! ** The |beta> projectors
       ! ------------------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb).and.ps%beta_2j(na)==ps%beta_2j(nb)) then

            work(:) = work(:) + beta(:,nb)*otfg_radin(ps%kkbeta,ab%rab,phi(:,nb)*beta(:,na))

          end if

       end do

       beta(:,na) = beta(:,na) - work(:)

       ! ---------------
       ! ** (T+Vae)|psi>
       ! ---------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb).and.ps%beta_2j(na)==ps%beta_2j(nb)) then

            work(:) = work(:) + tvpsi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       tvpsi(:,na) = tvpsi(:,na) - work(:)

       ! ---------------
       ! ** (T+Vps)|phi>
       ! ---------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb).and.ps%beta_2j(na)==ps%beta_2j(nb)) then

            work(:) = work(:) + tvphi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       tvphi(:,na) = tvphi(:,na) - work(:)

       ! -----------------------------
       ! ** The AE partial waves |psi>
       ! -----------------------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb).and.ps%beta_2j(na)==ps%beta_2j(nb)) then

            work(:) = work(:) + psi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       psi(:,na) = psi(:,na) - work(:)

       ! -----------------------------
       ! ** The PS partial waves |phi>
       ! -----------------------------

       work(:) = 0.0_dp

       do nb=1,na-1

          if(ps%beta_l(na)==ps%beta_l(nb).and.ps%beta_2j(na)==ps%beta_2j(nb)) then

            work(:) = work(:) + phi(:,nb)*otfg_radin(ps%kkbeta,ab%rab,beta(:,nb)*phi(:,na))

          end if

       end do

       phi(:,na) = phi(:,na) - work(:)

       ! ------------------------------------
       ! ** Scale beta so that <beta|phi> = 1
       ! ------------------------------------

       beta(:,na) = beta(:,na)/otfg_radin(ps%kkbeta,ab%rab,beta(:,na)*phi(:,na))

       ! ------------------------------
       ! ** Balance by normalising beta
       ! ------------------------------

       factor = 1.0_dp/sqrt(otfg_radin(ps%kkbeta,ab%rab,beta(:,na)*beta(:,na)))

       beta(:,na) = factor*beta(:,na)

       psi(:,na)   = psi(:,na)/factor
       phi(:,na)   = phi(:,na)/factor
       tvpsi(:,na) = tvpsi(:,na)/factor
       tvphi(:,na) = tvphi(:,na)/factor

    end do

    ! ----------------------------------------------
    ! ** Construct the charge augmentation functions
    !     -- and the integrated coefficients, qqq
    ! ----------------------------------------------

    qfunc = 0.0_dp
    qqq   = 0.0_dp

    !! CHECK for SOC
    if(ps%aug) then
       do na=1,ps%num_beta
          do nb=1,ps%num_beta
             qfunc(:,na,nb) = psi(:,na)*psi(:,nb)-phi(:,na)*phi(:,nb)
             if(ps%beta_l(na).eq.ps%beta_l(nb).and.ps%beta_2j(na)==ps%beta_2j(nb)) then
                qqq(na,nb) = otfg_radin(ps%kkbeta,ab%rab,qfunc(:,na,nb))
             else
                qqq(na,nb) = 0.0_dp
             end if
          end do
       end do
    end if
    !!

    ! ** Set the ddd matrix: <psi|T+Vae|psi> - <phi|T+Vps|phi>

    ddd = 0.0_dp
    do na=1,ps%num_beta
       do nb=1,ps%num_beta
          if(ps%beta_l(na).eq.ps%beta_l(nb).and.ps%beta_2j(na)==ps%beta_2j(nb)) ddd(na,nb) = &
               otfg_radin(ps%kkbeta,ab%rab,psi(:,na)*tvpsi(:,nb))-otfg_radin(ps%kkbeta,ab%rab,phi(:,na)*tvphi(:,nb))
!!$          if(na==nb) ddd(na,nb) = ddd(na,nb) + ps%beta_delta(na)*otfg_radin(ab%npts,ab%rab,phi(:,na)*phi(:,nb))
          if(na==nb) ps%beta_delta(na) = ps%beta_delta(na)*otfg_radin(ab%npts,ab%rab,phi(:,na)*phi(:,nb))
       end do
    end do

    ! ** Symmetrise the ddd matrix

    do na=1,ps%num_beta
       do nb=na,ps%num_beta
          ddd(na,nb) = (ddd(na,nb)+ddd(nb,na))/2.0_dp
          ddd(nb,na) = ddd(na,nb)
       end do
    end do

    ! ** Pseudise the Q-functions

    qfcoef = 0.0_dp

    if(ps%aug) then

       do na=1,ps%num_beta
          do nb=na,ps%num_beta
             llmin = abs(ps%beta_l(na)-ps%beta_l(nb))
             llmax = abs(ps%beta_l(na)+ps%beta_l(nb))
             do l=llmin,llmax
                if(mod(l-llmin,2)/=1) then
                   if(l.eq.llmin) then
                      qfunc(:,na,nb)=otfg_Q_pseudise(l,qfunc(:,na,nb),ps,ab,qfcoef(:,l,na,nb))
                      qfunc(:,nb,na)=qfunc(:,na,nb)
                      qfcoef(:,l,nb,na)=qfcoef(:,l,na,nb)
                   else
                      work = otfg_Q_pseudise(l,qfunc(:,na,nb),ps,ab,qfcoef(:,l,na,nb))
                      qfcoef(:,l,nb,na)=qfcoef(:,l,na,nb)
                   end if
                end if
             end do
          end do
       end do

    end if

    ! ============================================
    ! ** Output the beta functions for examination
    ! ============================================

    call otfg_write_beta(ab,ps,beta,V_loc,aeat%V_eff)

    ! =======================================================
    ! ** Convert into the correct form: reciprocal space etc.
    ! =======================================================

    ! ** Set up indexing

    ps%num_ps_projectors = 0
    do n=1,ps%num_beta
       ps%num_ps_projectors = ps%num_ps_projectors+2*ps%beta_l(n)+1 ! 2l+1
       ps%projector_l(ps%num_ps_projectors-2*ps%beta_l(n):ps%num_ps_projectors) = ps%beta_l(n)
       ps%projector_2j(ps%num_ps_projectors-2*ps%beta_l(n):ps%num_ps_projectors) = ps%beta_2j(n)
       ps%projector_m(ps%num_ps_projectors-2*ps%beta_l(n):ps%num_ps_projectors) = &
            (/(m,m=-ps%beta_l(n),ps%beta_l(n))/)
    end do

    allocate(nnn(ps%num_ps_projectors),stat=ierr)
    if (ierr/=0) call io_allocate_abort('nnn','otfg_construct_psp_dir')

    m = 0
    do n=1,ps%num_beta
       m = m+2*ps%beta_l(n)+1 ! 2l+1
       nnn(m-2*ps%beta_l(n):m) = n
    end do

    ! ** The local potential

    ps%V_locr(:ps%mesh) = V_loc(:ps%mesh)
    ps%has_V_tau = xc_is_mgga(aeat%theory)
    if( ps%has_V_tau ) then
!!$       call otfg_pspot_allocate_v_tau(ps)
       ps%V_taur(:ps%mesh) = V_tau(:ps%mesh)
    endif

    ! ** The beta functions

    do n=1,ps%num_ps_projectors
       call otfg_radial_transform(ps%beta_l(nnn(n)),maxval(ps%beta_l),1,ps%kkbeta,&
            ab%r,ab%rab,ps%num_points,ps%gmax,beta(:,nnn(n)),ps%beta(:,n))

       ps%betar(1:ab%npts,n) = beta(1:ab%npts,nnn(n))
    end do
    ps%beta = ps%beta*4*pi

    ! ** The projector weights

    ps%q(:,:)  = 0.0_dp
    ps%D0(:,:) = 0.0_dp

    do n=1,ps%num_ps_projectors
       do m=1,ps%num_ps_projectors
          if((ps%projector_l(n)==ps%projector_l(m))&!.and.(ps%projector_m(n)==ps%projector_m(m))&
               .and.ps%projector_2j(n)==ps%projector_2j(m)) then
             ps%q(n,m)  = qqq(nnn(n),nnn(m))
             ps%D0(n,m) = ddd(nnn(n),nnn(m)) !* Not D0 yet -- needs unscreening first
          else
             ps%q(n,m)  = 0.0_dp
             ps%D0(n,m) = 0.0_dp
          end if
       end do
    end do


    ! ** The augmentation functions, and grids

    do n=1,ps%num_ps_projectors
       do m=1,ps%num_ps_projectors
          ps%qfunc(1:ps%kkbeta,n,m)           = qfunc(1:ps%kkbeta,nnn(n),nnn(m))
          ps%qfcoef(1:ps%nqf,0:2*ab%lmax,n,m) = qfcoef(1:ps%nqf,0:2*ab%lmax,nnn(n),nnn(m))
       end do
    end do


    !************************************
    !** Now the Blochl GAMMA functions **
    !************************************

    ! Set the descreened potentials for use in the augmentations

    V_ae = aeat%V_eff
    ! reason for limit in next line is that ps%V_locr is set to maximum size across all species
    V_ps(1:ab%npts) = ps%V_locr(1:ab%npts)

!!$    if(ps%num_beta==ps%num_gamma) then
!!$
!!$       ! ** Fill up with the projectors already generated above
!!$
!!$       do na=1,ps%num_beta
!!$
!!$          psi_gamma(:,na) = psi(:,na)
!!$
!!$          phi_gamma(:,na) = phi(:,na)
!!$
!!$          gamma(:,na)     = beta(:,na)
!!$
!!$       end do
!!$
!!$    else

       ! ** Now calculate the extra ones

       gamma_loop: do na=1,ps%num_gamma

          ! ** Calculate the partial wave

          psi_gamma(:,na) = otfg_partial_wave(ps%gamma_l(na),ps%gamma_e(na),ps%kkbeta,ab,aeat)

          ! ** Pseudise the partial wave, and construct the gamma-functions analytically

          phi_gamma(:,na) = otfg_pseudise_partial_wave(.false.,ps%flags,ps%gamma_l(na),ps%gamma_norm(na),psi_gamma(:,na),&
               ps%gamma_nrc(na),ps%gamma_e(na),ab,aeat,gamma(:,na),vloc=V_loc,scheme=ps%gamma_scheme(na))

          ! ** Gamma may disappear

          if(abs(otfg_radin(ps%kkbeta,ab%rab,gamma(:,na)*gamma(:,na)))<1.0e-20_dp) then
             write(stdout,*) 'here'
             gamma(:,na) = ab%r**(l+1)*exp(-(ab%r/(0.7_dp*ab%r(ps%gamma_nrc(na))))**6.0_dp)

          end if

       end do gamma_loop

       ! =======================================================
       ! ** Construct the gamma functions using Bloechl's scheme
       ! ** Eqns. 90 to 96 PRB 50, 17953, 1994
       ! =======================================================

       do na = 1,ps%num_gamma

          ! ------------------------
          ! ** The |beta> projectors
          ! ------------------------

          work(:) = 0.0_dp

          do nb=1,na-1

             if(ps%gamma_l(na)==ps%gamma_l(nb).and.ps%gamma_2j(na)==ps%gamma_2j(nb)) then

                work(:) = work(:) + gamma(:,nb)*otfg_radin(ps%kkbeta,ab%rab,phi_gamma(:,nb)*gamma(:,na))

             end if

          end do

          gamma(:,na) = gamma(:,na) - work(:)

          ! -----------------------------
          ! ** The AE partial waves |psi>
          ! -----------------------------

          work(:) = 0.0_dp

          do nb=1,na-1

             if(ps%gamma_l(na)==ps%gamma_l(nb).and.ps%gamma_2j(na)==ps%gamma_2j(nb)) then

                work(:) = work(:) + psi_gamma(:,nb)*otfg_radin(ps%kkbeta,ab%rab,gamma(:,nb)*phi_gamma(:,na))

             end if

          end do

          psi_gamma(:,na) = psi_gamma(:,na) - work(:)

          ! -----------------------------
          ! ** The PS partial waves |phi>
          ! -----------------------------

          work(:) = 0.0_dp

          do nb=1,na-1

             if(ps%gamma_l(na)==ps%gamma_l(nb).and.ps%gamma_2j(na)==ps%gamma_2j(nb)) then

                work(:) = work(:) + phi_gamma(:,nb)*otfg_radin(ps%kkbeta,ab%rab,gamma(:,nb)*phi_gamma(:,na))

             end if

          end do

          phi_gamma(:,na) = phi_gamma(:,na) - work(:)

          ! ------------------------------------
          ! ** Scale beta so that <beta|phi> = 1
          ! ------------------------------------

          gamma(:,na) = gamma(:,na)/otfg_radin(ps%kkbeta,ab%rab,gamma(:,na)*phi_gamma(:,na))

          ! ------------------------------
          ! ** Balance by normalising beta
          ! ------------------------------

          factor = 1.0_dp/sqrt(otfg_radin(ps%kkbeta,ab%rab,gamma(:,na)*gamma(:,na)))

          gamma(:,na) = factor*gamma(:,na)

          psi_gamma(:,na)   = psi_gamma(:,na)/factor
          phi_gamma(:,na)   = phi_gamma(:,na)/factor

       end do

 !   end if

    ! ===================================
    ! ** Calculate the augmentation terms
    ! ===================================



    bl_dia_nmrshield  = 0.0_dp
    bl_para_nmrshield = 0.0_dp
    bl_para_damped    = 0.0_dp
    bl_efg_aug        = 0.0_dp
    bl_hf_dipolar     = 0.0_dp
    bl_hf_dipolar_cutoff     = 0.0_dp
    bl_hf_fermi       = 0.0_dp
    bl_hf_zora        = 0.0_dp
    bl_zke_aug_gtensor = 0.0_dp
    bl_sodia_g_tensor  = 0.0_dp
    bl_sopara_g_tensor = 0.0_dp

    do na=1,ps%num_gamma
       do nb=1,ps%num_gamma


          ! * Electric Field Gradient Augmentation
          ! this is the same as the hyperfine dipolar term - at least for the moment we keep them distinct

          work=0.0_dp
          work(1)=0.0_dp
          work(2:)=( (psi_gamma(2:,na)*psi_gamma(2:,nb))-(phi_gamma(2:,na)*phi_gamma(2:,nb)) ) &
               /(ab%r(2:)**3)

          bl_efg_aug(na,nb)=otfg_radin(ps%gamma_nrc(na),ab%rab,work)

                 end do
    end do

    ! =============================================
    ! ** Output the gamma functions for examination
    ! =============================================

    !!JRYif(ps%num_gamma>ps%num_beta)
    call otfg_write_gamma(ab,ps,gamma)

    ! * Copy data into pseudopotential structure

    ps%bl_dia_nmrshift  = 0.0_dp
    ps%bl_para_nmrshift = 0.0_dp
    ps%bl_para_damped   = 0.0_dp
    ps%bl_efg_aug       = 0.0_dp
    ps%bl_hf_dipolar    = 0.0_dp
    ps%bl_hf_dipolar_cutoff    = 0.0_dp
    ps%bl_hf_zora    = 0.0_dp
    ps%bl_hf_fermi      = 0.0_dp
    ps%bl_zke_aug_gtensor  = 0.0_dp
    ps%bl_sodia_g_tensor  = 0.0_dp
    ps%bl_sopara_g_tensor  = 0.0_dp

    ! ** Set up gamma indexing

    ps%num_bl_projectors = 0
    do n=1,ps%num_gamma
       ps%num_bl_projectors = ps%num_bl_projectors+2*ps%gamma_l(n)+1 ! 2l+1
       ps%bl_projector_l(ps%num_bl_projectors-2*ps%gamma_l(n):ps%num_bl_projectors) = ps%gamma_l(n)
       ps%bl_projector_m(ps%num_bl_projectors-2*ps%gamma_l(n):ps%num_bl_projectors) = (/(m,m=-ps%gamma_l(n),ps%gamma_l(n))/)
    end do


    allocate(nnn_bl(ps%num_bl_projectors),stat=ierr)
    if (ierr/=0) call io_allocate_abort('nnn_bl','otfg_construct_psp')

    m = 0
    do n=1,ps%num_gamma
       m = m+2*ps%gamma_l(n)+1 ! 2l+1
       nnn_bl(m-2*ps%gamma_l(n):m) = n
    end do

    ! ** The gamma functions

    do n=1,ps%num_bl_projectors

       call otfg_radial_transform(ps%gamma_l(nnn_bl(n)),maxval(ps%gamma_l),1,maxval(ps%gamma_nrc)+2,ab%r,ab%rab,&
            ps%num_points,ps%gmax,gamma(:,nnn_bl(n)),ps%gamma(:,n))

    end do
    ps%gamma = ps%gamma*4*pi

    ! ** The projector weights

    do n=1,ps%num_bl_projectors
       do m=1,ps%num_bl_projectors
          ps%bl_dia_nmrshift(n,m)  = bl_dia_nmrshield(nnn_bl(n),nnn_bl(m))
          ps%bl_para_nmrshift(n,m) = bl_para_nmrshield(nnn_bl(n),nnn_bl(m))
          ps%bl_para_damped(n,m) = bl_para_damped(nnn_bl(n),nnn_bl(m))
          ps%bl_efg_aug(n,m)       = bl_efg_aug(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_dipolar(n,m)    = bl_hf_dipolar(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_dipolar_cutoff(n,m)    = bl_hf_dipolar_cutoff(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_zora(n,m)    = bl_hf_zora(nnn_bl(n),nnn_bl(m))
          ps%bl_hf_fermi(n,m)      = bl_hf_fermi(nnn_bl(n),nnn_bl(m))
          ps%bl_zke_aug_gtensor(n,m)  = bl_zke_aug_gtensor(nnn_bl(n),nnn_bl(m))
          ps%bl_sodia_g_tensor(n,m)  = bl_sodia_g_tensor(nnn_bl(n),nnn_bl(m))
          ps%bl_sopara_g_tensor(n,m) = bl_sopara_g_tensor(nnn_bl(n),nnn_bl(m))
       end do
    end do


    !------------------------------------------------------------------
    ! ** Generate the terms for core excited matrix elements ---ELNES
    !------------------------------------------------------------------
    allocate(core_radials_r(ab%npts,1:15),stat=ierr)
    if (ierr/=0) call io_abort('Error in allocating core_radials_r in otfg_construct_psp_dir')


    ! ** Deallocate the arrays

    deallocate(psi_loc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_loc in otfg_construct_psp_dir')
    deallocate(phi_loc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi_loc in otfg_construct_psp_dir')
    deallocate(V_loc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_loc in otfg_construct_psp_dir')
    deallocate(V_tau,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_tau in otfg_construct_psp_paw')

    !deallocate(zora_K,stat=ierr)
    !if (ierr/=0) call io_abort('Error in deallocating zora_K in otfg_construct_psp_dir')
    !deallocate(zora_K_dr,stat=ierr)
    !if (ierr/=0) call io_abort('Error in deallocating zora_K_dr in otfg_construct_psp_dir')

    deallocate(psi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi in otfg_construct_psp_dir')
    deallocate(phi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi in otfg_construct_psp_dir')
    deallocate(tvpsi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating tvpsi in otfg_construct_psp_dir')
    deallocate(tvphi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating tvphi in otfg_construct_psp_dir')
    deallocate(psi_shift,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_shift in otfg_construct_psp_dir')
    deallocate(phi_shift,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi_shift in otfg_construct_psp_dir')
    deallocate(beta,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating beta in otfg_construct_psp_dir')
    deallocate(qqq,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating qqq in otfg_construct_psp_dir')
    deallocate(ddd,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating ddd in otfg_construct_psp_dir')
    deallocate(qfunc,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating qfunc in otfg_construct_psp_dir')
    deallocate(qfcoef,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating qfcoef in otfg_construct_psp_dir')
    deallocate(nnn,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating nnn in otfg_construct_psp_dir')
    deallocate(work,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating work in otfg_construct_psp_dir')


    ! ** The GAMMA Arrays

    deallocate(psi_gamma,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_gamma in otfg_construct_psp')
    deallocate(phi_gamma,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating phi_gamma in otfg_construct_psp')
    deallocate(gamma,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating beta_gamma in otfg_construct_psp')
    deallocate(nnn_bl,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating nnn_bl in otfg_construct_psp')
    deallocate(bl_dia_nmrshield,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_dia_nmrshield in otfg_construct_psp')
    deallocate(bl_para_nmrshield,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_para_nmrshield in otfg_construct_psp')
    deallocate(bl_para_damped,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_para_damped in otfg_construct_psp')
    deallocate(bl_efg_aug,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_efg_aug in otfg_construct_psp')
    deallocate(bl_hf_fermi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_hf_fermi in otfg_construct_psp')
    deallocate(bl_hf_dipolar,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_hf_dipolar in otfg_construct_psp')
    deallocate(bl_hf_dipolar_cutoff,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_hf_dipolar_cutoff in otfg_construct_psp')
    deallocate(bl_zke_aug_gtensor,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_zke_aug_gtensor in otfg_construct_psp')
    deallocate(bl_sodia_g_tensor,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_sodia_g_tensor in otfg_construct_psp')
    deallocate(bl_sopara_g_tensor,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_sopara_g_tensor in otfg_construct_psp')
    deallocate(V_ps,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_ps in otfg_construct_psp')
    deallocate(V_ae,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating V_ae in otfg_construct_psp')

    ! core_radials_r for ELNES
    deallocate(core_radials_r,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating core_radials_r in otfg_construct_psp')
    deallocate(phi_gamma_temp,stat=ierr) ! ELNES
    if (ierr/=0) call io_abort('Error in deallocating phi_gamma_temp in otfg_construct_psp')
    deallocate(bl_ttt_temp,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bl_ttt_temp in otfg_construct_psp')
    deallocate(psi_gamma_temp,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating psi_gamma_temp in otfg_construct_psp')


    call trace_exit('otfg_construct_psp_dir',status)

    return
  end subroutine otfg_construct_psp_dir



  subroutine otfg_descreen_psp(ab,ps,aeat,nsp)

    use constants,  only : periodic_table_name
    use io,         only : io_abort,io_allocate_abort,stdout
    use algor,      only : algor_erf
    use comms,      only : my_node_id,root_node_id,comms_gcopy

    !=========================================================================!
    ! Descreen the pseudopotential                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) ab      : input : The atomic basis                                   !
    ! 2) ps      : inout : The pseudopotential                                !
    ! 3) aeat    : input : The all electron atom                              !
    ! 4) nsp     : input : The species                                        !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) constants                                                            !
    ! 2) io                                                                   !
    ! 3) algor                                                                !
    ! 4) comms                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 23rd July 2002                                  !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    type(ae_basis), intent(in)    :: ab
    type(pspot),    intent(inout) :: ps
    type(ae_atom),  intent(in)    :: aeat
    integer,        intent(in)    :: nsp

    !-------------------------------------------------------------------------!

    type(ae_atom)  :: psat ! ** NOTE: This is a pseudo atom -- in the ae basis

    integer       :: i,j,n,nn,m,mm,l,no,ierr,nrc,j_loop,local_2j

    integer, dimension(:), allocatable :: nnn

    real(kind=dp) :: Ehar,Exc,g,oog,r2new,rcmax

    real(kind=dp), dimension(:), allocatable :: work,work2, v_tau

    real(kind=dp) :: proj(ps%num_ps_projectors)

    logical       :: on_root

    integer                         :: status

    call trace_entry('otfg_descreen_psp',status)

    allocate(work(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('work','otfg_descreen_psp')

    allocate(work2(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('work2','otfg_descreen_psp')

    allocate(v_tau(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('v_tau','otfg_descreen_psp')

    on_root = my_node_id==root_node_id

    ! ** Initialise the atom

    psat%spin_orbit_coupling=ps%SOC  ! need to know this for correct initialisation of atom
    call otfg_init_ae_atom(aeat%Z,ab,psat,ps)

    ! ** Calculate the occupations, orb_n, and initial eigenvalues
    if(ps%SOC) then
       call otfg_set_ae_occ_dir(ab,psat,nsp)
    else
       call otfg_set_ae_occ(ab,psat,nsp)
    endif

    psat%eval = aeat%eval

    do l=0,ab%lmax
       do no=1,psat%norbs(l)
          if(.not.ps%pseudise(no,l)) psat%occ(no,l) = 0.0_dp
          if(bc(nsp)%no(l)==no) psat%occ(no,l) = bc(nsp)%wvfn_occ(l)
       end do
    end do

    ! ** Solve the Schrodinger equation

    if( xc_is_mgga(aeat%theory) ) then
       if(ps%SOC) call io_abort('otfg_descreen_psp: no m-GGA with SOC')
       call otfg_ps_gks_solve(ab,psat,nsp)
    else
       if(ps%SOC) then
          call otfg_ps_solve_dir(ab,psat,nsp)
       else
          call otfg_ps_solve(ab,psat,nsp)
       end if
    endif

    ! ** Construct the pseudo charge density

    psat%rhorr = 0.0_dp
    do l=0,ab%lmax
       do no=1,psat%norbs(l)
          psat%rhorr = psat%rhorr + psat%occ(no,l)*psat%orbr(:,no,l)**2
       end do
    end do

    ! KED calculation
    work(1) = 0.0_dp
    work(2:) = psat%rhorr(2:) / ab%r(2:)

    call otfg_r_second_derivative(ab%npts,ab%rab,ab%r,work,work2,psat%kedrr)
    psat%kedrr = psat%kedrr * ab%r

    do l = 0, ab%lmax
       do no=1, aeat%norbs(l)
          call otfg_r_second_derivative(ab%npts, ab%rab, ab%r, &
               psat%orbr(:,no,l),work,work2)
          work(2:) = work2(2:) - l*(l+1)*psat%orbr(2:,no,l)/ab%r(2:)**2
          psat%kedrr = psat%kedrr - psat%occ(no,l)*2.0_dp*work*psat%orbr(:,no,l)
       enddo
    enddo
    psat%kedrr = 0.25_dp * psat%kedrr ! Definition with half

    ! * Augment the valence charge

!    if(psat%spin_orbit_coupling)then
if(psat%spin_orbit_coupling) then
  do l=0,ab%lmax
    do j_loop=1,2
      if(j_loop==2.and.l==0) cycle
      if(j_loop==1) then
        local_2j=2*l+1
      else
        local_2j=2*l-1
      endif

      do no=1,psat%norbs(l)

        do n=1,psat%num_proj
          if(psat%proj_l(n).eq.l.and.psat%proj_2j(n).eq.local_2j.and.psat%twoj(no,l).eq.local_2j) then
            proj(n)=otfg_radin(ab%npts,ab%rab,psat%betar(:,n)*psat%orbr(:,no,l))
          else
            proj(n)=0.0_dp
          end if
        end do

        nn=0
        do n=1,ps%num_ps_projectors
          if(ps%projector_m(n)==0) then
            nn=nn+1

            mm=0
            do m=1,ps%num_ps_projectors
              if(ps%projector_m(m)==0) then
                mm=mm+1

                psat%rhorr = psat%rhorr + psat%occ(no,l)*proj(nn)*proj(mm)*ps%qfunc(1:ab%npts,n,m)
                psat%kedrr = psat%kedrr + psat%occ(no,l)*proj(nn)*proj(mm)*ps%tfunc(1:ab%npts,0,n,m)

              end if

            end do
          end if
        end do

      end do
    end do
  end do
else
  do l=0,ab%lmax

    do no=1,psat%norbs(l)
      do n=1,psat%num_proj
        if(psat%proj_l(n).eq.l) then
          proj(n)=otfg_radin(ab%npts,ab%rab,psat%betar(:,n)*psat%orbr(:,no,l))
        else
          proj(n)=0.0_dp
        end if
      end do

      nn=0
      do n=1,ps%num_ps_projectors
        if(ps%projector_m(n)==0) then
          nn=nn+1

          mm=0
          do m=1,ps%num_ps_projectors
            if(ps%projector_m(m)==0) then
              mm=mm+1

              psat%rhorr = psat%rhorr + psat%occ(no,l)*proj(nn)*proj(mm)*ps%qfunc(1:ab%npts,n,m)
              psat%kedrr = psat%kedrr + psat%occ(no,l)*proj(nn)*proj(mm)*ps%tfunc(1:ab%npts,0,n,m)

            end if

          end do
        end if
      end do

    end do
  end do
endif


    ! ** For SCPSPs we should set rhorr above Rc

    if(any(aeat%rhorr_val>0.0_dp)) then

       do nrc=ab%npts,1,-1
          if(ab%r(nrc)<bc(nsp)%Rc) exit
       end do
       nrc=nrc+1

       psat%rhorr(nrc:) = aeat%rhorr(nrc:)

    end if
    ! ** Normalise to account for numerical inaccuracy

    if(.not.any(aeat%rhorr_val>0.0_dp)) then

       psat%rhorr = psat%rhorr/otfg_radin(ab%npts,ab%rab,psat%rhorr)*sum(psat%occ)

    end if

    ! ** Check for negative charge density

    do n=1,ab%npts
       if((psat%rhorr(n)<0.0_dp).and.on_root) then
          write (stdout,*) 'Warning: otfg_descreen_psp: negative charge density encountered'
          exit
       end if
    end do

    ! ** Unscreen the projector weights to construct ps%D0
    ! I'm doing this is a very cautious way - splitting up the SOC and scalar
    ! parts. At some piint I should safely merge them - but only once
    ! Ultrsofts are fully working ! JRY

  if(ps%SOC) then
      do n=1,ps%num_ps_projectors
        do m=1,ps%num_ps_projectors
            ! for the SOC case the D0 contains m/=m' elements as these are
            ! needed when we do the real to cmplx SP convertion inside nlpot.
            if((ps%projector_l(n)==ps%projector_l(m)).and.(ps%projector_2j(n)==ps%projector_2j(m))) then
              work(1) = 0.0_dp
              work(2:ps%kkbeta) = ps%V_locr(2:ps%kkbeta)*ps%qfunc(2:ps%kkbeta,n,m)/ab%r(2:ps%kkbeta)
              ps%D0(n,m) = ps%D0(n,m) - otfg_radin(ps%kkbeta,ab%rab,work)

              ! MGGA term
              if( ps%has_V_tau ) then
                  work(1) = 0.0_dp
                  work(2:ps%kkbeta) = ps%V_taur(2:ps%kkbeta)*ps%tfunc(2:ps%kkbeta,0,n,m) ! V_tau not actually multiplied by r!
                  ps%D0(n,m) = ps%D0(n,m) - otfg_radin(ps%kkbeta,ab%rab,work)
                endif
              else
                ps%D0(n,m) = 0.0_dp
              end if
            end do
          end do
    else
          do n=1,ps%num_ps_projectors
            do m=1,ps%num_ps_projectors
              if((ps%projector_l(n)==ps%projector_l(m)).and.(ps%projector_m(n)==ps%projector_m(m))) then
                  work(1) = 0.0_dp
                  work(2:ps%kkbeta) = ps%V_locr(2:ps%kkbeta)*ps%qfunc(2:ps%kkbeta,n,m)/ab%r(2:ps%kkbeta)
                  ps%D0(n,m) = ps%D0(n,m) - otfg_radin(ps%kkbeta,ab%rab,work)

                  ! MGGA term
                  if( ps%has_V_tau ) then
                    work(1) = 0.0_dp
                    work(2:ps%kkbeta) = ps%V_taur(2:ps%kkbeta)*ps%tfunc(2:ps%kkbeta,0,n,m) ! V_tau not actually multiplied by r!
                    ps%D0(n,m) = ps%D0(n,m) - otfg_radin(ps%kkbeta,ab%rab,work)
                  endif
                else
                  ps%D0(n,m) = 0.0_dp
                end if
              end do
            end do
          endif

    ! ** Add eigenvalue shift

    allocate(nnn(ps%num_ps_projectors),stat=ierr)
    if (ierr/=0) call io_allocate_abort('nnn','otfg_descreen_psp')

    m = 0
    do n=1,ps%num_beta
       m = m+2*ps%beta_l(n)+1 ! 2l+1
       nnn(m-2*ps%beta_l(n):m) = n
    end do

    do n=1,ps%num_ps_projectors
       ps%D0(n,n) = ps%D0(n,n) + ps%beta_delta(nnn(n))
    end do

    deallocate(nnn,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating nnn in otfg_descreen_psp')

    ! ** Unscreen the effective potential, psat%V_eff

    if(any(aeat%rhorr_val>0.0_dp)) then

       ! * If we are generating SCPSPs

       call otfg_ae_hartree(ab,psat%rhorr,work,Ehar,nrc,aeat%V_val(nrc))

       psat%V_eff = psat%V_eff - work*ab%r ! * Unscreen with the Hartree potential

    else

       call otfg_ae_hartree(ab,psat%rhorr,work,Ehar)

       psat%V_eff = psat%V_eff - work*ab%r ! * Unscreen with the Hartree potential

    end if

    call otfg_ae_xc(ab,psat%rhorr+ps%rs_core_charge(1:ab%npts),work,Exc,aeat%theory,&
         psat%kedrr+ps%rs_core_kinetic_energy_density(1:ab%npts),v_tau)

    psat%V_eff = psat%V_eff - work*ab%r ! * Unscreen with the XC potential
    psat%V_tau = psat%V_tau - v_tau

    ! ** For SCPSPs set the potential beyond Rc

    if(any(aeat%rhorr_val>0.0_dp)) then

       call otfg_ae_hartree(ab,aeat%rhorr-aeat%rhorr_val,work,Ehar)

       psat%V_eff(nrc:) = aeat%V_nuc(nrc:) + work(nrc:)*ab%r(nrc:)

       call otfg_ae_xc(ab,aeat%rhorr,work,Exc,aeat%theory)

       psat%V_eff(nrc:) = psat%V_eff(nrc:) + work(nrc:)*ab%r(nrc:)

       call otfg_ae_xc(ab,aeat%rhorr_val+ps%rs_core_charge(1:ab%npts),work,Exc,aeat%theory)

       psat%V_eff(nrc:) = psat%V_eff(nrc:) - work(nrc:)*ab%r(nrc:)

    end if

    !    --------------------------------------------------------------
    ! ** Transform to reciprocal space for use in the pseudopotential
    !    --------------------------------------------------------------

    ! ** Assume a 25 Ry cutoff, and make compatible with old database

    rcmax = 1.0_dp/sqrt(min(-100.0_dp/(4.0_dp*log(1.d-12)),0.7_dp))

    do i=1,ab%npts
       work(i)=psat%V_eff(i)+ps%ionic_charge*algor_erf(ab%r(i)/rcmax)
    enddo

    ! ** For g = 0

    r2new = 0.25_dp*(rcmax)**2
    ps%V_loc(1) = otfg_radin(ab%npts,ab%rab,work*ab%r) + ps%ionic_charge*r2new
    if( ps%has_V_tau ) &
         ps%V_tau(1) = otfg_radin(ab%npts,ab%rab,psat%V_tau*ab%r**2)

    ! ** For g > 0

    do i=2,ps%num_points

       g   = real(i-1,dp)/real(ps%num_points-1,dp)*ps%gmax
       oog = 1.0_dp/g

       do j=1,ab%npts
          work2(j)=work(j)*sin(ab%r(j)*g)
       enddo

       ps%V_loc(i)= otfg_radin(ab%npts,ab%rab,work2)*oog + ps%ionic_charge*(1.0_dp-exp(-r2new*g*g))*oog**2
       if( ps%has_V_tau ) &
            ps%V_tau(i)= otfg_radin(ab%npts,ab%rab,psat%V_tau*sin(ab%r*g)*ab%r)*oog

    end do

    deallocate(work,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating work in otfg_descreen_psp')
    deallocate(work2,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating work2 in otfg_descreen_psp')

    ! * Calculate the derivative of V_loc wrt g

    call otfg_derivative(ps%num_points,ps%gmax,ps%V_loc(:),ps%dV_loc_dg(:))
    if( ps%has_V_tau ) &
         call otfg_derivative(ps%num_points,ps%gmax,ps%V_tau(:),ps%dV_tau_dg(:))

    ! ** Destroy the pseudo_basis and pseudo_atom

    call otfg_ae_dealloc(psat)

    call trace_exit('otfg_descreen_psp',status)

    return
  end subroutine otfg_descreen_psp

  subroutine otfg_realspace_psp(ab,ps)

    use io,         only : io_abort,io_allocate_abort,stdout
    use parameters, only : current_params

    !=========================================================================!
    ! Convert the non-local part of the pseudopotential to real space         !
    ! according to the scheme of King-Smith, Payne and Lin:                   !
    ! PRB 44 13063-13066 (1991)                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) ab      : input : The atomic basis                                   !
    ! 2) ps      : inout : The pseudopotential                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 0) none                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) parameters                                                           !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 15th May 2003                                   !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    type(ae_basis), intent(in)    :: ab
    type(pspot),    intent(inout) :: ps

    !-------------------------------------------------------------------------!

    integer                                    :: nb,ierr
    integer                                    :: nGmax,ngamma
    real(kind=dp), dimension(:),   allocatable :: bb,chi
    real(kind=dp), dimension(:,:), allocatable :: aa
    real(kind=dp)                              :: wfn_gmax,pot_gmax


    integer                         :: status

    call trace_entry('otfg_realspace_psp',status)

    allocate(chi(ps%num_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('chi','otfg_realspace_psp')
    allocate(bb(ps%num_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('bb','otfg_realspace_psp')
    allocate(aa(ps%num_points,ps%num_points),stat=ierr)
    if (ierr/=0) call io_allocate_abort('aa','otfg_realspace_psp')

    ! ** Proceed according to the recipe of PRB 44 13063-13066 (1991)

    ! ** 1) Select the plane-wave cutoff Gmax and gamma

    wfn_gmax = sqrt(2.0_dp*current_params%cut_off_energy)

    if(current_params%fine_gmax<0.0_dp) then
       pot_gmax = wfn_gmax*current_params%fine_grid_scale
    else
       pot_gmax = current_params%fine_gmax
    end if

    nGmax  = 1 + int(wfn_gmax/ps%gmax*real(ps%num_points-1,dp))
    ngamma = 1 + int((2.0_dp*pot_gmax-wfn_gmax)/ps%gmax*real(ps%num_points-1,dp))

    if (nGmax>ngamma) call io_abort('Error: otfg_realspace_psp: nGmax > ngamma')
    if (ngamma>ps%num_points) call io_abort('Error: otfg_realspace_psp: ngamma>ps%num_points')

    do nb=1,ps%num_ps_projectors

       ! ** 2) Set chi_l(q) up to Gmax, zero above gamma = Gamma_1-Gmax

       chi(1:ngmax) = ps%beta(1:ngmax,nb)

       chi(ngamma+1:) = 0.0_dp

       ! ** 3) Select the real space cut-off radius (between 1.5 and 2.0 x r_cut)

       ! ** 4) Variationally set chi_l(q) for Gmax < q < gamma

       ! * Calculate the bb matrix: (lhs of eqn 8)/q**2

       ! * Calculate the aa matrix: (pi/2 q**2-A(q,q'))/q**2

       ! * Solve the simulaneous equation: bb = aa chi

    end do

    deallocate(chi,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating chi in otfg_realspace_psp')
    deallocate(bb,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating bb in otfg_realspace_psp')
    deallocate(aa,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating aa in otfg_realspace_psp')

    call trace_exit('otfg_realspace_psp',status)

    return
  end subroutine otfg_realspace_psp


  subroutine otfg_write_pwave(ab,ps,psi,phi,psi_loc,phi_loc)

    use io,         only : io_open_file,stdout
    use comms,      only : my_node_id,root_node_id
    use constants,  only : periodic_table_name

    !=========================================================================!
    ! Construct the pseudopotential                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) ab      : input : The atomic basis                                   !
    ! 2) ps      : inout : The pseudopotential                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) comms                                                                !
    ! 3) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 30th April 2004                                 !
    !=========================================================================!

    use parameters, only : current_params

    use trace, only: trace_entry, trace_exit

    implicit none

    type(ae_basis),                                intent(in) :: ab
    type(pspot),                                   intent(in) :: ps
    real(kind=dp), dimension(ab%npts,ps%num_beta), intent(in) :: psi
    real(kind=dp), dimension(ab%npts,ps%num_beta), intent(in) :: phi
    real(kind=dp), dimension(ab%npts),             intent(in) :: psi_loc
    real(kind=dp), dimension(ab%npts),             intent(in) :: phi_loc

    !-------------------------------------------------------------------------!

    integer                              :: nb,i,np,l,nm

    integer, parameter                   :: nflags = 5
    character(len=80), dimension(nflags) :: flags
    integer                              :: unit_pwave
    integer                              :: status

    real(kind=dp)                        :: ymin(0:4),ymax(0:4)

    logical :: on_root

    call trace_entry('otfg_write_pwave',status)

    do l=0,3
       if(count(ps%beta_l==l)>3) then
!          write (stdout,*) 'Warning: otfg_write_pwave not suitable for more than 3 projectors'
!          call trace_exit('otfg_write_pwave',status)
!          return
       end if
    end do

    on_root = my_node_id==root_node_id

    if(on_root.and.testing.and.current_params%write_otfg) then

       ! ** Open a .usp file

       flags(1) = 'WRITE'
       flags(2) = 'FORMATTED'
       flags(3) = 'UNKNOWN'
       flags(4) = 'REWIND'
       flags(5) = 'SEQUENTIAL'

       call io_open_file(trim(adjustl(periodic_table_name(ps%Z)))//&
            '_'//trim(adjustl(ps%set))//'_'//trim(adjustl(ps%theory))//&
            '_OTF.pwave',nflags,flags,unit_pwave,status)

       ! * Find the max-min range

       do nm=1,ab%npts
          if(ab%r(nm)>ps%rcmax*1.1_dp) exit
       end do

       ymin = 0.0_dp
       ymax = 0.0_dp

       do nb=1,ps%num_beta
          if(minval(psi(1:nm,nb))<ymin(ps%beta_l(nb))) ymin(ps%beta_l(nb)) = minval(psi(1:nm,nb))
          if(minval(phi(1:nm,nb))<ymin(ps%beta_l(nb))) ymin(ps%beta_l(nb)) = minval(phi(1:nm,nb))
          if(maxval(psi(1:nm,nb))>ymax(ps%beta_l(nb))) ymax(ps%beta_l(nb)) = maxval(psi(1:nm,nb))
          if(maxval(phi(1:nm,nb))>ymax(ps%beta_l(nb))) ymax(ps%beta_l(nb)) = maxval(phi(1:nm,nb))
       end do
       if(minval(psi_loc(1:nm))<ymin(ps%local_l)) ymin(ps%local_l) = minval(psi_loc(1:nm))
       if(minval(phi_loc(1:nm))<ymin(ps%local_l)) ymin(ps%local_l) = minval(phi_loc(1:nm))
       if(maxval(psi_loc(1:nm))>ymax(ps%local_l)) ymax(ps%local_l) = maxval(psi_loc(1:nm))
       if(maxval(phi_loc(1:nm))>ymax(ps%local_l)) ymax(ps%local_l) = maxval(phi_loc(1:nm))

       do l=0,3
          if(ymax(l)-ymin(l)==0.0_dp) ymax(l)=1.0_dp
       end do

       ymax = ymax * 1.1_dp
       ymin = ymin * 1.1_dp

       write (unit_pwave,'(a)') '@version 50109'
       write (unit_pwave,'(a)') '@default linewidth 1.5'
       write (unit_pwave,'(a)') '@g0 on'
       write (unit_pwave,'(a)') '@with g0'
       write (unit_pwave,'(a)') '@    world xmin 0'
       write (unit_pwave,'(a,f10.2)') '@    world xmax',ps%rcmax*1.1_dp
       write (unit_pwave,'(a,f10.2)') '@    world ymin',ymin(0)
       write (unit_pwave,'(a,f10.2)') '@    world ymax',ymax(0)
       write (unit_pwave,'(a)') '@    view xmin 0.050000'
       write (unit_pwave,'(a)') '@    view xmax 0.618987'
       write (unit_pwave,'(a)') '@    view ymin 0.550000'
       write (unit_pwave,'(a)') '@    view ymax 0.925000'
       if(ps%local_l==0) then
          write (unit_pwave,'(a)') '@    subtitle "s local"'
       else
          write (unit_pwave,'(a)') '@    subtitle "s"'
       end if
       write (unit_pwave,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major 1'
       write (unit_pwave,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    yaxis  tick off'
       write (unit_pwave,'(a)') '@    yaxis  ticklabel off'
       write (unit_pwave,'(a)') '@    s0 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line color 1'
       write (unit_pwave,'(a)') '@    s1 line linestyle 2'
       write (unit_pwave,'(a)') '@    s2 line color 2'
       write (unit_pwave,'(a)') '@    s2 line linestyle 3'
       write (unit_pwave,'(a)') '@    s2 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s3 line color 2'
       write (unit_pwave,'(a)') '@    s4 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s4 line color 1'
       write (unit_pwave,'(a)') '@    s4 line linestyle 2'
       write (unit_pwave,'(a)') '@    s5 line color 4'
       write (unit_pwave,'(a)') '@    s5 line linestyle 3'
       write (unit_pwave,'(a)') '@    s5 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s6 line color 4'
       write (unit_pwave,'(a)') '@    s7 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s7 line color 1'
       write (unit_pwave,'(a)') '@    s7 line linestyle 2'
       write (unit_pwave,'(a)') '@    s8 line color 3'
       write (unit_pwave,'(a)') '@    s8 line linestyle 3'
       write (unit_pwave,'(a)') '@    s8 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s9 line color 3'
       write (unit_pwave,'(a)') '@g1 on'
       write (unit_pwave,'(a)') '@with g1'
       write (unit_pwave,'(a)') '@    world xmin 0'
       write (unit_pwave,'(a,f10.2)') '@    world xmax',ps%rcmax*1.1_dp
       write (unit_pwave,'(a,f10.2)') '@    world ymin',ymin(1)
       write (unit_pwave,'(a,f10.2)') '@    world ymax',ymax(1)
       write (unit_pwave,'(a)') '@    view xmin 0.675885'
       write (unit_pwave,'(a)') '@    view xmax 1.244872'
       write (unit_pwave,'(a)') '@    view ymin 0.550000'
       write (unit_pwave,'(a)') '@    view ymax 0.925000'
       if(ps%local_l==1) then
          write (unit_pwave,'(a)') '@    subtitle "p local"'
       else
          write (unit_pwave,'(a)') '@    subtitle "p"'
       end if
       write (unit_pwave,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major 1'
       write (unit_pwave,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    yaxis  tick off'
       write (unit_pwave,'(a)') '@    yaxis  ticklabel off'
       write (unit_pwave,'(a)') '@    s0 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line color 1'
       write (unit_pwave,'(a)') '@    s1 line linestyle 2'
       write (unit_pwave,'(a)') '@    s2 line color 2'
       write (unit_pwave,'(a)') '@    s2 line linestyle 3'
       write (unit_pwave,'(a)') '@    s2 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s3 line color 2'
       write (unit_pwave,'(a)') '@    s4 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s4 line color 1'
       write (unit_pwave,'(a)') '@    s4 line linestyle 2'
       write (unit_pwave,'(a)') '@    s5 line color 4'
       write (unit_pwave,'(a)') '@    s5 line linestyle 3'
       write (unit_pwave,'(a)') '@    s5 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s6 line color 4'
       write (unit_pwave,'(a)') '@    s7 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s7 line color 1'
       write (unit_pwave,'(a)') '@    s7 line linestyle 2'
       write (unit_pwave,'(a)') '@    s8 line color 3'
       write (unit_pwave,'(a)') '@    s8 line linestyle 3'
       write (unit_pwave,'(a)') '@    s8 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s9 line color 3'
       write (unit_pwave,'(a)') '@g2 on'
       write (unit_pwave,'(a)') '@with g2'
       write (unit_pwave,'(a)') '@    world xmin 0'
       write (unit_pwave,'(a,f10.2)') '@    world xmax',ps%rcmax*1.1_dp
       write (unit_pwave,'(a,f10.2)') '@    world ymin',ymin(2)
       write (unit_pwave,'(a,f10.2)') '@    world ymax',ymax(2)
       write (unit_pwave,'(a)') '@    view xmin 0.050000'
       write (unit_pwave,'(a)') '@    view xmax 0.618987'
       write (unit_pwave,'(a)') '@    view ymin 0.100000'
       write (unit_pwave,'(a)') '@    view ymax 0.475000'
       if(ps%local_l==2) then
          write (unit_pwave,'(a)') '@    subtitle "d local"'
       else
          write (unit_pwave,'(a)') '@    subtitle "d"'
       end if
       write (unit_pwave,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major 1'
       write (unit_pwave,'(a)') '@    xaxis  label "r (Bohr)"'
       write (unit_pwave,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    yaxis  tick off'
       write (unit_pwave,'(a)') '@    yaxis  ticklabel off'
       write (unit_pwave,'(a)') '@    s0 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line color 1'
       write (unit_pwave,'(a)') '@    s1 line linestyle 2'
       write (unit_pwave,'(a)') '@    s2 line color 2'
       write (unit_pwave,'(a)') '@    s2 line linestyle 3'
       write (unit_pwave,'(a)') '@    s2 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s3 line color 2'
       write (unit_pwave,'(a)') '@    s4 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s4 line color 1'
       write (unit_pwave,'(a)') '@    s4 line linestyle 2'
       write (unit_pwave,'(a)') '@    s5 line color 4'
       write (unit_pwave,'(a)') '@    s5 line linestyle 3'
       write (unit_pwave,'(a)') '@    s5 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s6 line color 4'
       write (unit_pwave,'(a)') '@    s7 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s7 line color 1'
       write (unit_pwave,'(a)') '@    s7 line linestyle 2'
       write (unit_pwave,'(a)') '@    s8 line color 3'
       write (unit_pwave,'(a)') '@    s8 line linestyle 3'
       write (unit_pwave,'(a)') '@    s8 line linewidth 1.0'
       write (unit_pwave,'(a)') '@g3 on'
       write (unit_pwave,'(a)') '@with g3'
       write (unit_pwave,'(a)') '@    world xmin 0'
       write (unit_pwave,'(a,f10.2)') '@    world xmax',ps%rcmax*1.1_dp
       write (unit_pwave,'(a,f10.2)') '@    world ymin',ymin(3)
       write (unit_pwave,'(a,f10.2)') '@    world ymax',ymax(3)
       write (unit_pwave,'(a)') '@    view xmin 0.675885'
       write (unit_pwave,'(a)') '@    view xmax 1.244872'
       write (unit_pwave,'(a)') '@    view ymin 0.100000'
       write (unit_pwave,'(a)') '@    view ymax 0.475000'
       if(ps%local_l==3) then
          write (unit_pwave,'(a)') '@    subtitle "f local"'
       else
          write (unit_pwave,'(a)') '@    subtitle "f"'
       end if
       write (unit_pwave,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_pwave,'(a)') '@    xaxis  tick major 1'
       write (unit_pwave,'(a)') '@    xaxis  label "r (Bohr)"'
       write (unit_pwave,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_pwave,'(a)') '@    yaxis  tick off'
       write (unit_pwave,'(a)') '@    yaxis  ticklabel off'
       write (unit_pwave,'(a)') '@    s0 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s1 line color 1'
       write (unit_pwave,'(a)') '@    s1 line linestyle 2'
       write (unit_pwave,'(a)') '@    s2 line color 2'
       write (unit_pwave,'(a)') '@    s2 line linestyle 3'
       write (unit_pwave,'(a)') '@    s2 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s3 line color 2'
       write (unit_pwave,'(a)') '@    s4 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s4 line color 1'
       write (unit_pwave,'(a)') '@    s4 line linestyle 2'
       write (unit_pwave,'(a)') '@    s5 line color 4'
       write (unit_pwave,'(a)') '@    s5 line linestyle 3'
       write (unit_pwave,'(a)') '@    s5 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s6 line color 4'
       write (unit_pwave,'(a)') '@    s7 line linewidth 1.0'
       write (unit_pwave,'(a)') '@    s7 line color 1'
       write (unit_pwave,'(a)') '@    s7 line linestyle 2'
       write (unit_pwave,'(a)') '@    s8 line color 3'
       write (unit_pwave,'(a)') '@    s8 line linestyle 3'
       write (unit_pwave,'(a)') '@    s8 line linewidth 1.0'
       write (unit_pwave,'(a)') '@target G0.S0'
       write (unit_pwave,'(a)') '@type xy'
       write (unit_pwave,*) 0,0
       write (unit_pwave,*) ab%r(ps%kkbeta),0
       write (unit_pwave,'(a)')'&'
       np=0
       if(ps%local_l==0) then
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G0.S',np
          write (unit_pwave,'(a)') '@type xy'
          write (unit_pwave,*) ps%local_rc,ymin(0)
          write (unit_pwave,*) ps%local_rc,ymax(0)
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G0.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),psi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G0.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),phi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==0) then
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G0.S',np
             write (unit_pwave,'(a)') '@type xy'
             write (unit_pwave,*) ps%beta_rc(nb),ymin(0)
             write (unit_pwave,*) ps%beta_rc(nb),ymax(0)
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G0.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),psi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G0.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),phi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
          end if
       end do
       write (unit_pwave,'(a)') '@target G1.S0'
       write (unit_pwave,'(a)') '@type xy'
       write (unit_pwave,*) 0,0
       write (unit_pwave,*) ab%r(ps%kkbeta),0
       write (unit_pwave,'(a)')'&'
       np=0
       if(ps%local_l==1) then
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G1.S',np
          write (unit_pwave,'(a)') '@type xy'
          write (unit_pwave,*) ps%local_rc,ymin(1)
          write (unit_pwave,*) ps%local_rc,ymax(1)
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G1.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),psi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G1.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),phi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==1) then
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G1.S',np
             write (unit_pwave,'(a)') '@type xy'
             write (unit_pwave,*) ps%beta_rc(nb),ymin(1)
             write (unit_pwave,*) ps%beta_rc(nb),ymax(1)
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G1.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),psi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G1.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),phi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
          end if
       end do
       write (unit_pwave,'(a)') '@target G2.S0'
       write (unit_pwave,'(a)') '@type xy'
       write (unit_pwave,*) 0,0
       write (unit_pwave,*) ab%r(ps%kkbeta),0
       write (unit_pwave,'(a)')'&'
       np=0
       if(ps%local_l==2) then
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G2.S',np
          write (unit_pwave,'(a)') '@type xy'
          write (unit_pwave,*) ps%local_rc,ymin(2)
          write (unit_pwave,*) ps%local_rc,ymax(2)
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G2.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),psi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G2.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),phi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==2) then
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G2.S',np
             write (unit_pwave,'(a)') '@type xy'
             write (unit_pwave,*) ps%beta_rc(nb),ymin(2)
             write (unit_pwave,*) ps%beta_rc(nb),ymax(2)
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G2.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),psi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G2.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),phi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
          end if
       end do
       write (unit_pwave,'(a)') '@target G3.S0'
       write (unit_pwave,'(a)') '@type xy'
       write (unit_pwave,*) 0,0
       write (unit_pwave,*) ab%r(ps%kkbeta),0
       write (unit_pwave,'(a)')'&'
       np=0
       if(ps%local_l==3) then
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G3.S',np
          write (unit_pwave,'(a)') '@type xy'
          write (unit_pwave,*) ps%local_rc,ymin(3)
          write (unit_pwave,*) ps%local_rc,ymax(3)
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G3.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),psi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
          np=np+1
          write (unit_pwave,'(a,i1)') '@target G3.S',np
          write (unit_pwave,'(a)') '@type xy'
          do i=1,ps%kkbeta
             write (unit_pwave,*) ab%r(i),phi_loc(i)
          end do
          write (unit_pwave,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==3) then
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G3.S',np
             write (unit_pwave,'(a)') '@type xy'
             write (unit_pwave,*) ps%beta_rc(nb),ymin(3)
             write (unit_pwave,*) ps%beta_rc(nb),ymax(3)
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G3.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),psi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
             np=np+1
             write (unit_pwave,'(a,i1)') '@target G3.S',np
             write (unit_pwave,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_pwave,*) ab%r(i),phi(i,nb)
             end do
             write (unit_pwave,'(a)')'&'
          end if
       end do

       close(unit_pwave)

    end if

    call trace_exit('otfg_write_pwave',status)

    return
  end subroutine otfg_write_pwave

  subroutine otfg_write_beta(ab,ps,beta,ps_loc,ae_loc)

    use io,         only : io_open_file,stdout
    use comms,      only : my_node_id,root_node_id
    use constants,  only : periodic_table_name

    !=========================================================================!
    ! Write out the beta functions                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) ab      : input : The atomic basis                                   !
    ! 2) ps      : inout : The pseudopotential                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) comms                                                                !
    ! 3) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 30th April 2004                                 !
    !=========================================================================!

    use parameters, only : current_params

    use trace, only: trace_entry, trace_exit

    implicit none

    type(ae_basis),                                intent(in) :: ab
    type(pspot),                                   intent(in) :: ps
    real(kind=dp), dimension(ab%npts,ps%num_beta), intent(in) :: beta
    real(kind=dp), dimension(ab%npts),             intent(in) :: ps_loc
    real(kind=dp), dimension(ab%npts),             intent(in) :: ae_loc

    !-------------------------------------------------------------------------!

    integer                              :: nb,i,np,l

    integer, parameter                   :: nflags = 5
    character(len=80), dimension(nflags) :: flags
    integer                              :: unit_beta
    integer                              :: status

    real(kind=dp)                        :: xmax,ymin(0:4),ymax(0:4)

    logical :: on_root

    call trace_entry('otfg_write_beta',status)

    do l=0,3
       if(count(ps%beta_l==l)>3) then
!          write (stdout,*) 'Warning: otfg_write_beta not suitable for more than 3 projectors'
!          call trace_exit('otfg_write_beta',status)

!          return
       end if
    end do

    on_root = my_node_id==root_node_id
    if(on_root.and.testing.and.current_params%write_otfg) then

       ! ** Open a .usp file

       flags(1) = 'WRITE'
       flags(2) = 'FORMATTED'
       flags(3) = 'UNKNOWN'
       flags(4) = 'REWIND'
       flags(5) = 'SEQUENTIAL'

       call io_open_file(trim(adjustl(periodic_table_name(ps%Z)))//&
            '_'//trim(adjustl(ps%set))//'_'//trim(adjustl(ps%theory))//&
            '_OTF.beta',nflags,flags,unit_beta,status)

       ! * Find the max-min range

       xmax = max(maxval(ps%beta_rc),ps%local_rc)*1.2_dp

       ymin = 0.0_dp
       ymax = 0.0_dp

       do nb=1,ps%num_beta
          if(minval(beta(:,nb))<ymin(ps%beta_l(nb))) ymin(ps%beta_l(nb)) = minval(beta(:,nb))
          if(maxval(beta(:,nb))>ymax(ps%beta_l(nb))) ymax(ps%beta_l(nb)) = maxval(beta(:,nb))
       end do
       if(minval(ps_loc(2:)/ab%r(2:))<ymin(ps%local_l)) ymin(ps%local_l) = minval(ps_loc(2:)/ab%r(2:))
       if(maxval(ps_loc(2:)/ab%r(2:))>ymax(ps%local_l)) ymax(ps%local_l) = maxval(ps_loc(2:)/ab%r(2:))

       do l=0,3
          if(ymax(l)-ymin(l)==0.0_dp) ymax(l)=1.0_dp
       end do

       ymax = ymax * 1.1_dp
       ymin = ymin * 1.1_dp

       write (unit_beta,'(a)') '@version 50109'
       write (unit_beta,'(a)') '@default linewidth 1.5'
       write (unit_beta,'(a)') '@g0 on'
       write (unit_beta,'(a)') '@with g0'
       write (unit_beta,'(a)') '@    world xmin 0'
       write (unit_beta,'(a,f10.2)') '@    world xmax',xmax
       write (unit_beta,'(a,f10.2)') '@    world ymin',ymin(0)
       write (unit_beta,'(a,f10.2)') '@    world ymax',ymax(0)
       write (unit_beta,'(a)') '@    view xmin 0.050000'
       write (unit_beta,'(a)') '@    view xmax 0.618987'
       write (unit_beta,'(a)') '@    view ymin 0.550000'
       write (unit_beta,'(a)') '@    view ymax 0.925000'
       if(ps%local_l==0) then
          write (unit_beta,'(a)') '@    subtitle "s V local"'
       else
          write (unit_beta,'(a)') '@    subtitle "s"'
       end if
       write (unit_beta,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major 1'
       write (unit_beta,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    yaxis  tick off'
       write (unit_beta,'(a)') '@    yaxis  ticklabel off'
       write (unit_beta,'(a)') '@    s0 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line color 1'
       write (unit_beta,'(a)') '@    s1 line linestyle 2'
       write (unit_beta,'(a)') '@    s2 line color 2'
       write (unit_beta,'(a)') '@    s2 line linestyle 3'
       write (unit_beta,'(a)') '@    s2 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s3 line color 2'
       write (unit_beta,'(a)') '@    s4 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s4 line color 1'
       write (unit_beta,'(a)') '@    s4 line linestyle 2'
       write (unit_beta,'(a)') '@    s5 line color 4'
       write (unit_beta,'(a)') '@    s5 line linestyle 3'
       write (unit_beta,'(a)') '@    s5 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s6 line color 4'
       write (unit_beta,'(a)') '@    s7 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s7 line color 1'
       write (unit_beta,'(a)') '@    s7 line linestyle 2'
       write (unit_beta,'(a)') '@    s8 line color 3'
       write (unit_beta,'(a)') '@    s8 line linestyle 3'
       write (unit_beta,'(a)') '@    s8 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s9 line color 3'
       write (unit_beta,'(a)') '@g1 on'
       write (unit_beta,'(a)') '@with g1'
       write (unit_beta,'(a)') '@    world xmin 0'
       write (unit_beta,'(a,f10.2)') '@    world xmax',xmax
       write (unit_beta,'(a,f10.2)') '@    world ymin',ymin(1)
       write (unit_beta,'(a,f10.2)') '@    world ymax',ymax(1)
       write (unit_beta,'(a)') '@    view xmin 0.675885'
       write (unit_beta,'(a)') '@    view xmax 1.244872'
       write (unit_beta,'(a)') '@    view ymin 0.550000'
       write (unit_beta,'(a)') '@    view ymax 0.925000'
       if(ps%local_l==1) then
          write (unit_beta,'(a)') '@    subtitle "p V local"'
       else
          write (unit_beta,'(a)') '@    subtitle "p"'
       end if
       write (unit_beta,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major 1'
       write (unit_beta,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    yaxis  tick off'
       write (unit_beta,'(a)') '@    yaxis  ticklabel off'
       write (unit_beta,'(a)') '@    s0 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line color 1'
       write (unit_beta,'(a)') '@    s1 line linestyle 2'
       write (unit_beta,'(a)') '@    s2 line color 2'
       write (unit_beta,'(a)') '@    s2 line linestyle 3'
       write (unit_beta,'(a)') '@    s2 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s3 line color 2'
       write (unit_beta,'(a)') '@    s4 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s4 line color 1'
       write (unit_beta,'(a)') '@    s4 line linestyle 2'
       write (unit_beta,'(a)') '@    s5 line color 4'
       write (unit_beta,'(a)') '@    s5 line linestyle 3'
       write (unit_beta,'(a)') '@    s5 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s6 line color 4'
       write (unit_beta,'(a)') '@    s7 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s7 line color 1'
       write (unit_beta,'(a)') '@    s7 line linestyle 2'
       write (unit_beta,'(a)') '@    s8 line color 3'
       write (unit_beta,'(a)') '@    s8 line linestyle 3'
       write (unit_beta,'(a)') '@    s8 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s9 line color 3'
       write (unit_beta,'(a)') '@g2 on'
       write (unit_beta,'(a)') '@with g2'
       write (unit_beta,'(a)') '@    world xmin 0'
       write (unit_beta,'(a,f10.2)') '@    world xmax',xmax
       write (unit_beta,'(a,f10.2)') '@    world ymin',ymin(2)
       write (unit_beta,'(a,f10.2)') '@    world ymax',ymax(2)
       write (unit_beta,'(a)') '@    view xmin 0.050000'
       write (unit_beta,'(a)') '@    view xmax 0.618987'
       write (unit_beta,'(a)') '@    view ymin 0.100000'
       write (unit_beta,'(a)') '@    view ymax 0.475000'
       if(ps%local_l==2) then
          write (unit_beta,'(a)') '@    subtitle "d V local"'
       else
          write (unit_beta,'(a)') '@    subtitle "d"'
       end if
       write (unit_beta,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major 1'
       write (unit_beta,'(a)') '@    xaxis  label "r (Bohr)"'
       write (unit_beta,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    yaxis  tick off'
       write (unit_beta,'(a)') '@    yaxis  ticklabel off'
       write (unit_beta,'(a)') '@    s0 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line color 1'
       write (unit_beta,'(a)') '@    s1 line linestyle 2'
       write (unit_beta,'(a)') '@    s2 line color 2'
       write (unit_beta,'(a)') '@    s2 line linestyle 3'
       write (unit_beta,'(a)') '@    s2 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s3 line color 2'
       write (unit_beta,'(a)') '@    s4 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s4 line color 1'
       write (unit_beta,'(a)') '@    s4 line linestyle 2'
       write (unit_beta,'(a)') '@    s5 line color 4'
       write (unit_beta,'(a)') '@    s5 line linestyle 3'
       write (unit_beta,'(a)') '@    s5 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s6 line color 4'
       write (unit_beta,'(a)') '@    s7 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s7 line color 1'
       write (unit_beta,'(a)') '@    s7 line linestyle 2'
       write (unit_beta,'(a)') '@    s8 line color 3'
       write (unit_beta,'(a)') '@    s8 line linestyle 3'
       write (unit_beta,'(a)') '@    s8 line linewidth 1.0'
       write (unit_beta,'(a)') '@g3 on'
       write (unit_beta,'(a)') '@with g3'
       write (unit_beta,'(a)') '@    world xmin 0'
       write (unit_beta,'(a,f10.2)') '@    world xmax',xmax
       write (unit_beta,'(a,f10.2)') '@    world ymin',ymin(3)
       write (unit_beta,'(a,f10.2)') '@    world ymax',ymax(3)
       write (unit_beta,'(a)') '@    view xmin 0.675885'
       write (unit_beta,'(a)') '@    view xmax 1.244872'
       write (unit_beta,'(a)') '@    view ymin 0.100000'
       write (unit_beta,'(a)') '@    view ymax 0.475000'
       if(ps%local_l==3) then
          write (unit_beta,'(a)') '@    subtitle "f V local"'
       else
          write (unit_beta,'(a)') '@    subtitle "f"'
       end if
       write (unit_beta,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_beta,'(a)') '@    xaxis  tick major 1'
       write (unit_beta,'(a)') '@    xaxis  label "r (Bohr)"'
       write (unit_beta,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_beta,'(a)') '@    yaxis  tick off'
       write (unit_beta,'(a)') '@    yaxis  ticklabel off'
       write (unit_beta,'(a)') '@    s0 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s1 line color 1'
       write (unit_beta,'(a)') '@    s1 line linestyle 2'
       write (unit_beta,'(a)') '@    s2 line color 2'
       write (unit_beta,'(a)') '@    s2 line linestyle 3'
       write (unit_beta,'(a)') '@    s2 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s3 line color 2'
       write (unit_beta,'(a)') '@    s4 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s4 line color 1'
       write (unit_beta,'(a)') '@    s4 line linestyle 2'
       write (unit_beta,'(a)') '@    s5 line color 4'
       write (unit_beta,'(a)') '@    s5 line linestyle 3'
       write (unit_beta,'(a)') '@    s5 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s6 line color 4'
       write (unit_beta,'(a)') '@    s7 line linewidth 1.0'
       write (unit_beta,'(a)') '@    s7 line color 1'
       write (unit_beta,'(a)') '@    s7 line linestyle 2'
       write (unit_beta,'(a)') '@    s8 line color 3'
       write (unit_beta,'(a)') '@    s8 line linestyle 3'
       write (unit_beta,'(a)') '@    s8 line linewidth 1.0'
       write (unit_beta,'(a)') '@target G0.S0'
       write (unit_beta,'(a)') '@type xy'
       write (unit_beta,*) 0,0
       write (unit_beta,*) ab%r(ps%kkbeta),0
       write (unit_beta,'(a)')'&'
       np=0
       if(ps%local_l==0) then
          np=np+1
          write (unit_beta,'(a,i1)') '@target G0.S',np
          write (unit_beta,'(a)') '@type xy'
          write (unit_beta,*) ps%local_rc,ymin(0)
          write (unit_beta,*) ps%local_rc,ymax(0)
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G0.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ae_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G0.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ps_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==0) then
             np=np+1
             write (unit_beta,'(a,i1)') '@target G0.S',np
             write (unit_beta,'(a)') '@type xy'
             write (unit_beta,*) ps%beta_rc(nb),ymin(0)
             write (unit_beta,*) ps%beta_rc(nb),ymax(0)
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G0.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G0.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_beta,*) ab%r(i),beta(i,nb)
             end do
             write (unit_beta,'(a)')'&'
          end if
       end do
       write (unit_beta,'(a)') '@target G1.S0'
       write (unit_beta,'(a)') '@type xy'
       write (unit_beta,*) 0,0
       write (unit_beta,*) ab%r(ps%kkbeta),0
       write (unit_beta,'(a)')'&'
       np=0
       if(ps%local_l==1) then
          np=np+1
          write (unit_beta,'(a,i1)') '@target G1.S',np
          write (unit_beta,'(a)') '@type xy'
          write (unit_beta,*) ps%local_rc,ymin(1)
          write (unit_beta,*) ps%local_rc,ymax(1)
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G1.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ae_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G1.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ps_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==1) then
             np=np+1
             write (unit_beta,'(a,i1)') '@target G1.S',np
             write (unit_beta,'(a)') '@type xy'
             write (unit_beta,*) ps%beta_rc(nb),ymin(1)
             write (unit_beta,*) ps%beta_rc(nb),ymax(1)
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G1.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G1.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_beta,*) ab%r(i),beta(i,nb)
             end do
             write (unit_beta,'(a)')'&'
          end if
       end do
       write (unit_beta,'(a)') '@target G2.S0'
       write (unit_beta,'(a)') '@type xy'
       write (unit_beta,*) 0,0
       write (unit_beta,*) ab%r(ps%kkbeta),0
       write (unit_beta,'(a)')'&'
       np=0
       if(ps%local_l==2) then
          np=np+1
          write (unit_beta,'(a,i1)') '@target G2.S',np
          write (unit_beta,'(a)') '@type xy'
          write (unit_beta,*) ps%local_rc,ymin(2)
          write (unit_beta,*) ps%local_rc,ymax(2)
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G2.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ae_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G2.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ps_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==2) then
             np=np+1
             write (unit_beta,'(a,i1)') '@target G2.S',np
             write (unit_beta,'(a)') '@type xy'
             write (unit_beta,*) ps%beta_rc(nb),ymin(2)
             write (unit_beta,*) ps%beta_rc(nb),ymax(2)
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G2.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G2.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_beta,*) ab%r(i),beta(i,nb)
             end do
             write (unit_beta,'(a)')'&'
          end if
       end do
       write (unit_beta,'(a)') '@target G3.S0'
       write (unit_beta,'(a)') '@type xy'
       write (unit_beta,*) 0,0
       write (unit_beta,*) ab%r(ps%kkbeta),0
       write (unit_beta,'(a)')'&'
       np=0
       if(ps%local_l==3) then
          np=np+1
          write (unit_beta,'(a,i1)') '@target G3.S',np
          write (unit_beta,'(a)') '@type xy'
          write (unit_beta,*) ps%local_rc,ymin(3)
          write (unit_beta,*) ps%local_rc,ymax(3)
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G3.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ae_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
          np=np+1
          write (unit_beta,'(a,i1)') '@target G3.S',np
          write (unit_beta,'(a)') '@type xy'
          do i=2,ps%kkbeta
             write (unit_beta,*) ab%r(i),ps_loc(i)/ab%r(i)
          end do
          write (unit_beta,'(a)')'&'
       end if
       do nb=1,ps%num_beta
          if(ps%beta_l(nb)==3) then
             np=np+1
             write (unit_beta,'(a,i1)') '@target G3.S',np
             write (unit_beta,'(a)') '@type xy'
             write (unit_beta,*) ps%beta_rc(nb),ymin(3)
             write (unit_beta,*) ps%beta_rc(nb),ymax(3)
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G3.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_beta,'(a)')'&'
             np=np+1
             write (unit_beta,'(a,i1)') '@target G3.S',np
             write (unit_beta,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_beta,*) ab%r(i),beta(i,nb)
             end do
             write (unit_beta,'(a)')'&'
          end if
       end do

       close(unit_beta)

    end if

    call trace_exit('otfg_write_beta',status)

    return
  end subroutine otfg_write_beta

  subroutine otfg_write_gamma(ab,ps,gamma)

    use io,         only : io_open_file,stdout
    use comms,      only : my_node_id,root_node_id
    use constants,  only : periodic_table_name

    !=========================================================================!
    ! Write the gamma projectors in xmgrace format                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) ab      : input : The atomic basis                                   !
    ! 2) ps      : inout : The pseudopotential                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) comms                                                                !
    ! 3) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 30th April 2004 (cloned by JRY)                 !
    ! gamma version doesn't print the local potential                         !
    !=========================================================================!

    use parameters, only : current_params

    use trace, only: trace_entry, trace_exit

    implicit none

    type(ae_basis),                                intent(in) :: ab
    type(pspot),                                   intent(in) :: ps
    real(kind=dp), dimension(ab%npts,ps%num_gamma), intent(in) :: gamma

    !-------------------------------------------------------------------------!

    integer                              :: nb,i,np,l

    integer, parameter                   :: nflags = 5
    character(len=80), dimension(nflags) :: flags
    integer                              :: unit_gamma
    integer                              :: status

    real(kind=dp)                        :: xmax,ymin(0:3),ymax(0:3)

    logical :: on_root

    call trace_entry('otfg_write_gamma',status)

    do l=0,3
       if(count(ps%gamma_l==l)>3) then
!          write (stdout,*) 'Warning: otfg_write_gamma not suitable for more than 3 projectors'
!          call trace_exit('otfg_write_gamma',status)

!          return
       end if
    end do

    on_root = my_node_id==root_node_id

    if(on_root.and.testing.and.current_params%write_otfg) then

       ! ** Open a .usp file

       flags(1) = 'WRITE'
       flags(2) = 'FORMATTED'
       flags(3) = 'UNKNOWN'
       flags(4) = 'REWIND'
       flags(5) = 'SEQUENTIAL'

       call io_open_file(trim(adjustl(periodic_table_name(ps%Z)))//&
            '_'//trim(adjustl(ps%set))//'_'//trim(adjustl(ps%theory))//&
            '_OTF.gamma',nflags,flags,unit_gamma,status)

       ! * Find the max-min range

       xmax = max(maxval(ps%gamma_rc),ps%local_rc)*1.2_dp

       ymin = 0.0_dp
       ymax = 0.0_dp

       do nb=1,ps%num_gamma
          if(minval(gamma(:,nb))<ymin(ps%gamma_l(nb))) ymin(ps%gamma_l(nb)) = minval(gamma(:,nb))
          if(maxval(gamma(:,nb))>ymax(ps%gamma_l(nb))) ymax(ps%gamma_l(nb)) = maxval(gamma(:,nb))
       end do

       do l=0,3
          if(ymax(l)-ymin(l)==0.0_dp) ymax(l)=1.0_dp
       end do

       ymax = ymax * 1.1_dp
       ymin = ymin * 1.1_dp

       write (unit_gamma,'(a)') '@version 50109'
       write (unit_gamma,'(a)') '@default linewidth 1.5'
       write (unit_gamma,'(a)') '@g0 on'
       write (unit_gamma,'(a)') '@with g0'
       write (unit_gamma,'(a)') '@    world xmin 0'
       write (unit_gamma,'(a,f10.2)') '@    world xmax',xmax
       write (unit_gamma,'(a,f10.2)') '@    world ymin',ymin(0)
       write (unit_gamma,'(a,f10.2)') '@    world ymax',ymax(0)
       write (unit_gamma,'(a)') '@    view xmin 0.050000'
       write (unit_gamma,'(a)') '@    view xmax 0.618987'
       write (unit_gamma,'(a)') '@    view ymin 0.550000'
       write (unit_gamma,'(a)') '@    view ymax 0.925000'

       write (unit_gamma,'(a)') '@    subtitle "s"'

       write (unit_gamma,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major 1'
       write (unit_gamma,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    yaxis  tick off'
       write (unit_gamma,'(a)') '@    yaxis  ticklabel off'
       write (unit_gamma,'(a)') '@    s0 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line color 1'
       write (unit_gamma,'(a)') '@    s1 line linestyle 2'
       write (unit_gamma,'(a)') '@    s2 line color 2'
       write (unit_gamma,'(a)') '@    s2 line linestyle 3'
       write (unit_gamma,'(a)') '@    s2 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s3 line color 2'
       write (unit_gamma,'(a)') '@    s4 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s4 line color 1'
       write (unit_gamma,'(a)') '@    s4 line linestyle 2'
       write (unit_gamma,'(a)') '@    s5 line color 4'
       write (unit_gamma,'(a)') '@    s5 line linestyle 3'
       write (unit_gamma,'(a)') '@    s5 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s6 line color 4'
       write (unit_gamma,'(a)') '@    s7 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s7 line color 1'
       write (unit_gamma,'(a)') '@    s7 line linestyle 2'
       write (unit_gamma,'(a)') '@    s8 line color 3'
       write (unit_gamma,'(a)') '@    s8 line linestyle 3'
       write (unit_gamma,'(a)') '@    s8 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s9 line color 3'
       write (unit_gamma,'(a)') '@g1 on'
       write (unit_gamma,'(a)') '@with g1'
       write (unit_gamma,'(a)') '@    world xmin 0'
       write (unit_gamma,'(a,f10.2)') '@    world xmax',xmax
       write (unit_gamma,'(a,f10.2)') '@    world ymin',ymin(1)
       write (unit_gamma,'(a,f10.2)') '@    world ymax',ymax(1)
       write (unit_gamma,'(a)') '@    view xmin 0.675885'
       write (unit_gamma,'(a)') '@    view xmax 1.244872'
       write (unit_gamma,'(a)') '@    view ymin 0.550000'
       write (unit_gamma,'(a)') '@    view ymax 0.925000'

       write (unit_gamma,'(a)') '@    subtitle "p"'

       write (unit_gamma,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major 1'
       write (unit_gamma,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    yaxis  tick off'
       write (unit_gamma,'(a)') '@    yaxis  ticklabel off'
       write (unit_gamma,'(a)') '@    s0 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line color 1'
       write (unit_gamma,'(a)') '@    s1 line linestyle 2'
       write (unit_gamma,'(a)') '@    s2 line color 2'
       write (unit_gamma,'(a)') '@    s2 line linestyle 3'
       write (unit_gamma,'(a)') '@    s2 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s3 line color 2'
       write (unit_gamma,'(a)') '@    s4 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s4 line color 1'
       write (unit_gamma,'(a)') '@    s4 line linestyle 2'
       write (unit_gamma,'(a)') '@    s5 line color 4'
       write (unit_gamma,'(a)') '@    s5 line linestyle 3'
       write (unit_gamma,'(a)') '@    s5 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s6 line color 4'
       write (unit_gamma,'(a)') '@    s7 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s7 line color 1'
       write (unit_gamma,'(a)') '@    s7 line linestyle 2'
       write (unit_gamma,'(a)') '@    s8 line color 3'
       write (unit_gamma,'(a)') '@    s8 line linestyle 3'
       write (unit_gamma,'(a)') '@    s8 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s9 line color 3'
       write (unit_gamma,'(a)') '@g2 on'
       write (unit_gamma,'(a)') '@with g2'
       write (unit_gamma,'(a)') '@    world xmin 0'
       write (unit_gamma,'(a,f10.2)') '@    world xmax',xmax
       write (unit_gamma,'(a,f10.2)') '@    world ymin',ymin(2)
       write (unit_gamma,'(a,f10.2)') '@    world ymax',ymax(2)
       write (unit_gamma,'(a)') '@    view xmin 0.050000'
       write (unit_gamma,'(a)') '@    view xmax 0.618987'
       write (unit_gamma,'(a)') '@    view ymin 0.100000'
       write (unit_gamma,'(a)') '@    view ymax 0.475000'

       write (unit_gamma,'(a)') '@    subtitle "d"'

       write (unit_gamma,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major 1'
       write (unit_gamma,'(a)') '@    xaxis  label "r (Bohr)"'
       write (unit_gamma,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    yaxis  tick off'
       write (unit_gamma,'(a)') '@    yaxis  ticklabel off'
       write (unit_gamma,'(a)') '@    s0 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line color 1'
       write (unit_gamma,'(a)') '@    s1 line linestyle 2'
       write (unit_gamma,'(a)') '@    s2 line color 2'
       write (unit_gamma,'(a)') '@    s2 line linestyle 3'
       write (unit_gamma,'(a)') '@    s2 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s3 line color 2'
       write (unit_gamma,'(a)') '@    s4 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s4 line color 1'
       write (unit_gamma,'(a)') '@    s4 line linestyle 2'
       write (unit_gamma,'(a)') '@    s5 line color 4'
       write (unit_gamma,'(a)') '@    s5 line linestyle 3'
       write (unit_gamma,'(a)') '@    s5 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s6 line color 4'
       write (unit_gamma,'(a)') '@    s7 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s7 line color 1'
       write (unit_gamma,'(a)') '@    s7 line linestyle 2'
       write (unit_gamma,'(a)') '@    s8 line color 3'
       write (unit_gamma,'(a)') '@    s8 line linestyle 3'
       write (unit_gamma,'(a)') '@    s8 line linewidth 1.0'
       write (unit_gamma,'(a)') '@g3 on'
       write (unit_gamma,'(a)') '@with g3'
       write (unit_gamma,'(a)') '@    world xmin 0'
       write (unit_gamma,'(a,f10.2)') '@    world xmax',xmax
       write (unit_gamma,'(a,f10.2)') '@    world ymin',ymin(3)
       write (unit_gamma,'(a,f10.2)') '@    world ymax',ymax(3)
       write (unit_gamma,'(a)') '@    view xmin 0.675885'
       write (unit_gamma,'(a)') '@    view xmax 1.244872'
       write (unit_gamma,'(a)') '@    view ymin 0.100000'
       write (unit_gamma,'(a)') '@    view ymax 0.475000'

       write (unit_gamma,'(a)') '@    subtitle "f"'

       write (unit_gamma,'(a)') '@    xaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick minor linewidth 1.5'
       write (unit_gamma,'(a)') '@    xaxis  tick major 1'
       write (unit_gamma,'(a)') '@    xaxis  label "r (Bohr)"'
       write (unit_gamma,'(a)') '@    yaxis  bar linewidth 1.5'
       write (unit_gamma,'(a)') '@    yaxis  tick off'
       write (unit_gamma,'(a)') '@    yaxis  ticklabel off'
       write (unit_gamma,'(a)') '@    s0 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s1 line color 1'
       write (unit_gamma,'(a)') '@    s1 line linestyle 2'
       write (unit_gamma,'(a)') '@    s2 line color 2'
       write (unit_gamma,'(a)') '@    s2 line linestyle 3'
       write (unit_gamma,'(a)') '@    s2 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s3 line color 2'
       write (unit_gamma,'(a)') '@    s4 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s4 line color 1'
       write (unit_gamma,'(a)') '@    s4 line linestyle 2'
       write (unit_gamma,'(a)') '@    s5 line color 4'
       write (unit_gamma,'(a)') '@    s5 line linestyle 3'
       write (unit_gamma,'(a)') '@    s5 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s6 line color 4'
       write (unit_gamma,'(a)') '@    s7 line linewidth 1.0'
       write (unit_gamma,'(a)') '@    s7 line color 1'
       write (unit_gamma,'(a)') '@    s7 line linestyle 2'
       write (unit_gamma,'(a)') '@    s8 line color 3'
       write (unit_gamma,'(a)') '@    s8 line linestyle 3'
       write (unit_gamma,'(a)') '@    s8 line linewidth 1.0'
       write (unit_gamma,'(a)') '@target G0.S0'
       write (unit_gamma,'(a)') '@type xy'
       write (unit_gamma,*) 0,0
       write (unit_gamma,*) ab%r(ps%kkbeta),0
       write (unit_gamma,'(a)')'&'
       np=0

       do nb=1,ps%num_gamma
          if(ps%gamma_l(nb)==0) then
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G0.S',np
             write (unit_gamma,'(a)') '@type xy'
             write (unit_gamma,*) ps%gamma_rc(nb),ymin(0)
             write (unit_gamma,*) ps%gamma_rc(nb),ymax(0)
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G0.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G0.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_gamma,*) ab%r(i),gamma(i,nb)
             end do
             write (unit_gamma,'(a)')'&'
          end if
       end do
       write (unit_gamma,'(a)') '@target G1.S0'
       write (unit_gamma,'(a)') '@type xy'
       write (unit_gamma,*) 0,0
       write (unit_gamma,*) ab%r(ps%kkbeta),0
       write (unit_gamma,'(a)')'&'
       np=0

       do nb=1,ps%num_gamma
          if(ps%gamma_l(nb)==1) then
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G1.S',np
             write (unit_gamma,'(a)') '@type xy'
             write (unit_gamma,*) ps%gamma_rc(nb),ymin(1)
             write (unit_gamma,*) ps%gamma_rc(nb),ymax(1)
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G1.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G1.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_gamma,*) ab%r(i),gamma(i,nb)
             end do
             write (unit_gamma,'(a)')'&'
          end if
       end do
       write (unit_gamma,'(a)') '@target G2.S0'
       write (unit_gamma,'(a)') '@type xy'
       write (unit_gamma,*) 0,0
       write (unit_gamma,*) ab%r(ps%kkbeta),0
       write (unit_gamma,'(a)')'&'
       np=0

       do nb=1,ps%num_gamma
          if(ps%gamma_l(nb)==2) then
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G2.S',np
             write (unit_gamma,'(a)') '@type xy'
             write (unit_gamma,*) ps%gamma_rc(nb),ymin(2)
             write (unit_gamma,*) ps%gamma_rc(nb),ymax(2)
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G2.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G2.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_gamma,*) ab%r(i),gamma(i,nb)
             end do
             write (unit_gamma,'(a)')'&'
          end if
       end do
       write (unit_gamma,'(a)') '@target G3.S0'
       write (unit_gamma,'(a)') '@type xy'
       write (unit_gamma,*) 0,0
       write (unit_gamma,*) ab%r(ps%kkbeta),0
       write (unit_gamma,'(a)')'&'
       np=0

       do nb=1,ps%num_gamma
          if(ps%gamma_l(nb)==3) then
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G3.S',np
             write (unit_gamma,'(a)') '@type xy'
             write (unit_gamma,*) ps%gamma_rc(nb),ymin(3)
             write (unit_gamma,*) ps%gamma_rc(nb),ymax(3)
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G3.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
             end do
             write (unit_gamma,'(a)')'&'
             np=np+1
             write (unit_gamma,'(a,i1)') '@target G3.S',np
             write (unit_gamma,'(a)') '@type xy'
             do i=1,ps%kkbeta
                write (unit_gamma,*) ab%r(i),gamma(i,nb)
             end do
             write (unit_gamma,'(a)')'&'
          end if
       end do

       close(unit_gamma)

    end if

    call trace_exit('otfg_write_gamma',status)

    return
  end subroutine otfg_write_gamma

  subroutine otfg_core_charge(nsp,ab,aeat,ps,rcc)

    use io,         only : io_abort,io_allocate_abort,io_failure,io_open_file,io_atomic_to_unit
    use parameters, only : current_params
    use cell,       only : current_cell
    use constants,  only : file_maxpath,mu_0_si,elementary_charge_si,electron_mass_si,speed_light_si,&
         fine_structure_si,planck_si,avogadro_si

    !=========================================================================!
    ! Define the pseudopotential                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) nsp     : input  : The species                                       !
    ! 2) ab      : input  : The atomic basis                                  !
    ! 3) aeat    : input  : The initialised atom                              !
    ! 4) ps      : output : The pseudopotential                               !
    ! 5) rcc     : output : The radial core charge                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    ! 1) ae_basis                                                             !
    ! 2) ae_atom                                                              !
    ! 3) pspot                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) io                                                                   !
    ! 2) parameters                                                           !
    ! 3) cell                                                                 !
    ! 4) constants                                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    ! 0) None                                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Chris J. Pickard, V0.0, 1st July 2002                                   !
    !=========================================================================!

    use trace, only: trace_entry, trace_exit

    implicit none

    integer,        intent(in)    :: nsp
    type(ae_basis), intent(in)    :: ab
    type(ae_atom),  intent(in)    :: aeat
    type(pspot),    intent(inout) :: ps
    real,           intent(inout), dimension(:), allocatable :: rcc

    !-------------------------------------------------------------------------!

    integer       :: i,j,n,nb,nbb,nref,no,l,indx,jndx,ierr,nn,nnmax  ! add nn,nnmax for ELNES
    real(kind=dp) :: temp_rc,dE_use,scale,temp_rc2,dr

    real(kind=dp), dimension(:), allocatable :: rs_core_current

    character(len=file_maxpath) :: definition,temp_def
    character(len=80)           :: ctemp,ctemp2

    real(kind=dp), parameter :: dE = 0.1_dp ! * The energy step for multiple projectors

    logical :: swap,user,overide

    ! ** Define comments

    integer, parameter :: ncomm=4
    character(len=1)   :: comment(ncomm)

    integer                         :: status

    call trace_entry('otfg_define_psp',status)

    comment(1) = '#'           ! Define comment characters
    comment(2) = ';'
    comment(3) = '!'
    comment(4) = achar(13)     ! ^M

    ! ** Set the working pseudopotential definition

    if(trim(current_cell%species_pot(nsp)).ne."NULL.NULL") then
       select case(trim(current_cell%species_pot(nsp)))
       case("C7")
          ps_set(nsp)=1
          current_cell%species_pot(nsp)="NULL.NULL"
       case("MS")
          ps_set(nsp)=2
          current_cell%species_pot(nsp)="NULL.NULL"
       case("QC5")
          ps_set(nsp)=3
          current_cell%species_pot(nsp)="NULL.NULL"
      case("HARD")
          ps_set(nsp)=4
          current_cell%species_pot(nsp)="NULL.NULL"
      case("NCP")
          ps_set(nsp)=5
          current_cell%species_pot(nsp)="NULL.NULL"
      case("C8")
          ps_set(nsp)=6
          current_cell%species_pot(nsp)="NULL.NULL"
       end select
    end if

    select case(trim(current_cell%species_pot(nsp)))
    case("NULL.NULL")
       definition = adjustl(ps_def(aeat%Z,ps_set(nsp)))
       do j=1,ncomm
          i=index(definition,comment(j))
          if(i.gt.0) definition(i:) = ' ' ! Strip out the comment lines
       end do
    case default
       definition = adjustl(current_cell%species_pot(nsp))
    end select

    ! ** Check if we have simply been passed a scale factor

    if(index(definition,'|')==0) then
       read(definition,*,err=103) scale
       definition = adjustl(ps_def(aeat%Z,ps_set(nsp)))
    else
       scale = 1.0_dp
    end if

    if(index(trim(definition),' ')>0) call io_abort('Error: otfg_define_psp: incorrect syntax, remove blanks')

    ! ** Output testing files?

    testing = .false.
    indx = index(definition,'[')
    if(indx/=0) testing = .true.

    ! ** Check if there any pseudisation flags are set

    indx = index(definition,'(')
    if(indx/=0) then

       ! * Remove the flag string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,')')
       ps%flags = ctemp(:indx-1)
       if(indx==0) call io_abort('Error: otfg_define_psp: incorrect syntax for pseudisation flags')
    end if

    ! ** Check if there any configuration adjustment required

    indx = index(definition,'{')
    if(indx/=0) then

       ! * Remove the configuration adjust string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,'}')
       if(indx==0) call io_abort('Error: otfg_define_psp: incorrect syntax for config adjust')
    end if

    indx = index(definition,'[')
    if(indx/=0) then

       ! * Remove the configuration adjust string

       ctemp = definition(indx+1:)
       definition = definition(:indx-1)
       indx = index(ctemp,']')
       if(indx==0) call io_abort('Error: otfg_define_psp: incorrect syntax for config adjust')
    end if

    ! ** Set the local channel

    ps%pv = .false. ! * Don't do the polynomial fit by default

    indx = index(definition,'|')
    ctemp = definition(:indx-1)

    indx = index(ctemp,'=')
    if(indx>0) then
       overide = .true.
       read(ctemp(indx+1:),*,err=103) ps%local_e
       read(ctemp(:indx-1),*,err=103) ps%local_l
       ps%local_l = abs(ps%local_l)
       if(index(ctemp(:indx-1),'-')>0) ps%pv=.true.
    else
       overide = .false.
       ps%local_e    = 0.0_dp                  ! * The energy
       read(ctemp,*,err=103) ps%local_l
       ps%local_l = abs(ps%local_l)
       if(index(ctemp,'-')>0) ps%pv=.true.
    end if
    ps%local_norm = 0                       ! * Do not conserve the norm

    ! ** Find out about the core radii - if only one is specified, use defaults
    !    - note, they may be specified on a projector by projector basis in addition

    ! * Count the number of "|" -- if 5, then use defaults

    n=0
    do i=1,len_trim(definition)
       if(definition(i:i)=='|') n=n+1
    end do

    ! * Read the core radii in from the definition string

    jndx = index(definition,'|')+1

    if(n==5) then      ! * Only one rc specified - use defaults
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%local_rc
       temp_rc = ps%local_rc
       ps%rinner = temp_rc*0.7_dp ! * Default setting
       definition = definition(jndx+1+indx:) ! * Remove the portion read
    else if(n==7) then ! * local_rc, beta_rc and
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%local_rc
       jndx = jndx + 1 + indx
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) temp_rc
       jndx = jndx + 1 + indx
       indx = index(definition(jndx:),'|')-1
       read(definition(jndx:jndx-1+indx),*,err=103) ps%rinner(0)
       ps%rinner = ps%rinner(0)
       definition = definition(jndx+indx+1:) ! * Remove the portion read
    else
       call io_abort('Error: otfg_define_psp: unable to parse pseudopotential definition - too many "|"')
    end if

    temp_rc2 = temp_rc ! Keep a copy of default beta_rc to use for the gamma projectors

    ! ** Read in the definitions of COARSE, MEDIUM and FINE
    !    - note, the should be in atomic units (Ha)

    indx = index(definition(1:),'|')-1
    read(definition(1:indx),*,err=103) ps%coarse
    jndx = 2+indx
    indx = index(definition(jndx:),'|')-1
    read(definition(jndx:jndx-1+indx),*,err=103) ps%medium
    jndx = 1+jndx+indx
    indx = index(definition(jndx:),'|')-1
    read(definition(jndx:jndx-1+indx),*,err=103) ps%fine
    definition = definition(jndx+indx+1:) ! * Remove the portion read

    ! * Set PRECISE and EXTREME (we don't really know what EXTREME is)

    ps%precise = max(1.03_dp * ps%fine + 1.57_dp,1.2_dp*ps%fine)
    ps%extreme = 1.6_dp * ps%fine

    ! ** Set the cutoff energy

    select case(current_params%basis_precision)
    case('COARSE')  ; current_params%cut_off_energy = max(ps%coarse,current_params%cut_off_energy)
    case('MEDIUM')  ; current_params%cut_off_energy = max(ps%medium,current_params%cut_off_energy)
    case('FINE')    ; current_params%cut_off_energy = max(ps%fine,current_params%cut_off_energy)
    case('PRECISE') ; current_params%cut_off_energy = max(ps%precise,current_params%cut_off_energy)
    case('EXTREME') ; current_params%cut_off_energy = max(ps%extreme,current_params%cut_off_energy)
    case('NULL')
    case default
       call io_abort('Error: otfg_define_psp: don''t recognise current_params%basis_precision')
    end select

    ! ** Now read in the projector definitions

    ! * Count the number of reference states

    nref=1
    do i=1,len_trim(definition)
       if(definition(i:i)==':') nref=nref+1
    end do

    ! ** Find out how many beta functions there are -- so the arrays can be allocated

    temp_def    = definition
    ps%num_beta = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)
       nb = 0
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='N') nb=nb+1 ! * Count the norm conserving projectors
          if(ctemp(i:i)=='U') nb=nb+1 ! * Count the ultrasoft projectors
       end do
       if(nb==0) nb=2
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='L') nb=0 ! * Don't count the local channel
          if(ctemp(i:i)=='P') nb=0 ! * Don't count the dummy channel
       end do
       ps%num_beta = ps%num_beta + nb
    end do

    ! ** Allocate the ps%beta_* arrays

    if (allocated(ps%beta_l)) then
       deallocate(ps%beta_l,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_l in otfg_define_psp')
    end if
    allocate(ps%beta_l(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_l','otfg_define_psp')
    if (allocated(ps%beta_e)) then
       deallocate(ps%beta_e,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_e in otfg_define_psp')
    end if
    allocate(ps%beta_e(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_e','otfg_define_psp')
    if (allocated(ps%beta_rc)) then
       deallocate(ps%beta_rc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_rc in otfg_define_psp')
    end if
    allocate(ps%beta_rc(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_rc','otfg_define_psp')
    if (allocated(ps%beta_delta)) then
       deallocate(ps%beta_delta,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_delta in otfg_define_psp')
    end if
    allocate(ps%beta_delta(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_delta','otfg_define_psp')
    ps%beta_delta = 0.0_dp
    if (allocated(ps%beta_nrc)) then
       deallocate(ps%beta_nrc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_nrc in otfg_define_psp')
    end if
    allocate(ps%beta_nrc(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_nrc','otfg_define_psp')
    if (allocated(ps%beta_norm)) then
       deallocate(ps%beta_norm,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_norm in otfg_define_psp')
    end if
    allocate(ps%beta_norm(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_norm','otfg_define_psp')
    if (allocated(ps%beta_scheme)) then
       deallocate(ps%beta_scheme,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_scheme in otfg_define_psp')
    end if
    allocate(ps%beta_scheme(1:ps%num_beta),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%beta_scheme','otfg_define_psp')

    ! ** Set the ps%beta_* arrays
    ps%pseudise = .false.
    temp_def    = definition
    nbb         = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)

       read(ctemp(1:2),'(2i1)',err=103) no,l ! * Find out the nl for this reference state

       ps%pseudise(no-l,l) = .true.

       if(ctemp(3:3)=='P') cycle ! * Just pseudise - no projectors etc.

       if(ctemp(3:3)=='L') then

          if(l/=ps%local_l) call io_abort('Error: otfg_define_psp: local channel conflict')

          if(.not.overide) ps%local_e = aeat%eval(no-l,l)

          if(ctemp(4:4)=='G') then
             ps%local_norm = 0
          else
             ps%local_norm = 1
          end if

       else

          ps%beta_e(nbb+1) = aeat%eval(no-l,l)

          nb = 0
          do i=1,len_trim(ctemp)
             if(ctemp(i:i)=='N'.or.ctemp(i:i)=='U') then
                nb=nb+1 ! * Count the projectors

                do j=i+1,len_trim(ctemp)
                   if(ctemp(j:j)=='N'.or.ctemp(j:j)=='U'.or.ctemp(j:j)=='H'.or.ctemp(j:j)=='G') exit
                end do
                ctemp2 = ctemp(i+1:j-1)

                indx = index(ctemp2,'@')
                if(indx>0) then
                   read(ctemp2(indx+1:),*,err=103) ps%beta_delta(nbb+nb)
                   ctemp2 = ctemp2(:indx-1)
                end if

                indx = index(ctemp2,'+')
                jndx = index(ctemp2,'-')
                indx = max(indx,jndx)

                if(indx>0) then
                   read(ctemp2(indx:),*,err=103) dE_use
                   user = .true.
                   if(indx==1) then
                      ctemp2 = ' '
                   else
                      if(ctemp2(indx-1:indx-1)/='=') ctemp2 = ctemp2(:indx-1)
                   end if

                else
                   dE_use = dE*real(nb-1,dp)
                   user = .false.
                end if

                indx = index(ctemp2,'=')
                if(indx>0) then
                   if(len_trim(ctemp2(indx+1:))<1) goto 103
                   read(ctemp2(indx+1:),*,err=103) dE_use
                   ctemp2 = ctemp2(:indx-1)
                   ps%beta_e(nbb+nb) = dE_use
                else
                   if(user) then
                      ps%beta_e(nbb+nb) = aeat%eval(no-l,l) + dE_use
                   else
                      if(dE_use>0.0_dp) then
                         ps%beta_e(nbb+nb) = 0.25*dE_use/dE
                      else
                         ps%beta_e(nbb+nb) = aeat%eval(no-l,l)
                      end if
                   end if
                end if

                if(len_trim(ctemp2)>0) then
                   read(ctemp2,*,err=103) ps%beta_rc(nbb+nb)
                else
                   ps%beta_rc(nbb+nb) = temp_rc
                end if

                ps%beta_l(nbb+nb) = l

                if(ctemp(i:i)=='N') then
                   ps%beta_norm(nbb+nb) = 1 ! * Norm conserving
                else
                   ps%beta_norm(nbb+nb) = 0 ! * Ultrasoft
                end if

             end if
          end do

          ! * Didn't find any modifiers, so use defaults

          if(nb==0) then
             nb=2

             ps%beta_l(nbb+1:nbb+2) = l

             ps%beta_norm(nbb+1:nbb+2) = 0
             ps%beta_rc(nbb+1:nbb+2)   = temp_rc

             ! * Set the energies

             ps%beta_e(nbb+2) = 0.25

          end if

          nbb = nbb + nb

       end if

    end do

    ! ** Reorder the beta arrays into increasing l
    !     - ensure block diagonality

101 swap = .false.
    do nb=1,ps%num_beta-1

       if(ps%beta_l(nb)>ps%beta_l(nb+1)) then
          swap = .true.
          l = ps%beta_l(nb)
          ps%beta_l(nb)   = ps%beta_l(nb+1)
          ps%beta_l(nb+1) = l
          l = ps%beta_norm(nb)
          ps%beta_norm(nb)   = ps%beta_norm(nb+1)
          ps%beta_norm(nb+1) = l
          temp_rc = ps%beta_rc(nb)
          ps%beta_rc(nb)   = ps%beta_rc(nb+1)
          ps%beta_rc(nb+1) = temp_rc
          temp_rc = ps%beta_e(nb)
          ps%beta_e(nb)   = ps%beta_e(nb+1)
          ps%beta_e(nb+1) = temp_rc
       end if

    end do

    if(swap) goto 101

    ! ** Find the _nrc and set the _rc to exact points

    ps%local_rc = ps%local_rc*scale
    ps%beta_rc = ps%beta_rc*scale
    ps%rinner = ps%rinner*scale

    ps%local_rc = otfg_round_to_grid(ps%local_rc,ps%local_nrc,ab%r)

    do nb=1,ps%num_beta
       ps%beta_rc(nb) = otfg_round_to_grid(ps%beta_rc(nb),ps%beta_nrc(nb),ab%r)
    end do

    do l=0,2*ab%lmax
       ps%rinner(l) = otfg_round_to_grid(ps%rinner(l),ps%ninner(l),ab%r)
    end do

    ! ** Set ps%mesh

    ps%mesh   = ab%npts

    ! ** Find ps%kkbeta

    if(index(ps%flags,'es')>0) then

       if(index(ps%flags,'esr=')>0) then
          read(ps%flags(index(ps%flags,'esr=')+4:),*,err=103) ps%rcmax
          temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
          ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
       else
          ps%rcmax  = otfg_round_to_grid(5.0_dp*max(ps%local_rc,maxval(ps%beta_rc),maxval(ps%rinner)),nb,ab%r)
          temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
          ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
       end if

    else
       ps%rcmax  = otfg_round_to_grid(max(ps%local_rc,maxval(ps%beta_rc),maxval(ps%rinner)),nb,ab%r)
       temp_rc   = otfg_round_to_grid(2.0_dp*ps%rcmax,ps%kkbeta,ab%r)
       ps%kkbeta = min((ps%kkbeta/2)*2+1,ps%mesh) ! * Make odd, and less than or equal to ps%mesh
    end if

    ! ** Set the log grid

    ps%rlog(:ps%mesh)   = ab%r(:ps%mesh)
    ps%rab(:ps%mesh)    = ab%rab(:ps%mesh)

    ! ** Set the parameters for the augmentation function optimisation

    ps%nqpts  = 401

    if(index(ps%flags,'qc=')>0) then
       read(ps%flags(index(ps%flags,'qc=')+3:),*) ps%qc
       ps%qc = ps%qc*2.0_dp
    else
       ps%qc = 10.0_dp
    end if

    if(index(ps%flags,'Qc=')>0) then
       read(ps%flags(index(ps%flags,'Qc=')+3:),*) ps%qc
    end if

    ps%nderiv = 4
    ps%nextra = 3
    ps%nqf    = ps%nderiv+ps%nextra+1

    ! ** Set the parameters for the reciprocal space representation of the pseudopotential

    ps%num_points = 2001                                ! * This needs to be the same as in inquire_psp
    ps%gmax       = 52.91772086073879_dp                ! * So as to be the same as old CASTEP

    ! ** Allocate the core current

    allocate(rs_core_current(ab%npts),stat=ierr)
    if (ierr/=0) call io_allocate_abort('rs_core_current','otfg_define_psp')

    ! ** Set the ionic charge of the pseudopotential, the atomic number, and the core charge

    ps%Z = aeat%Z

    ps%ionic_charge = real(aeat%Z,dp)
    ps%rs_core_charge = 0.0_dp
    do l = 0,ab%lmax
       do n = 1,aeat%norbmax

          if(.not.ps%pseudise(n,l)) then
             ps%ionic_charge = ps%ionic_charge - aeat%occ(n,l)
             ps%rs_core_charge(:ps%mesh) = ps%rs_core_charge(:ps%mesh) + aeat%occ(n,l)*aeat%orbr(:ps%mesh,n,l)**2
          end if
       end do
    end do

    ! Store the core charge in rcc for external use

    allocate(rcc(ab%npts))
    do n=1,ab%npts
       rcc(n)=ps%rs_core_charge(n)
    end do

    ! ** Calculate the core current and the chemical shielding of the core electrons

    do n=2,ab%npts
       rs_core_current(n)=ps%rs_core_charge(n)/ab%r(n)
    end do
    rs_core_current(1)=0.0_dp

    ps%core_nmr_shift=(mu_0_si*speed_light_si*fine_structure_si*elementary_charge_si**2*&
         otfg_radin(ab%npts,ab%rab,rs_core_current)*1.0e6_dp)/(6.0_dp*planck_si)

    ! We now compute the magnetic suceptibility of the core electrons
    ! rather than allocate a new array we'll reuse  rs_core_current

    do n=2,ab%npts
       rs_core_current(n)=ps%rs_core_charge(n)*ab%r(n)*ab%r(n)
    end do
    rs_core_current(1)=0.0_dp

    ps%core_mag_sus=-1.0_dp*io_atomic_to_unit(1.0_dp,'cm')**3.d0*avogadro_si*1.d6*(fine_structure_si**2.d0)* &
         otfg_radin(ab%npts,ab%rab,rs_core_current)/6.0_dp

    ! ** Set the NLCC core charge

    ps%core_charge = (otfg_radin(ps%kkbeta,ab%rab,ps%rs_core_charge)>0.0_dp)

    if(index(ps%flags,'nonlcc')>0) ps%core_charge = .false.

    if(ps%core_charge) then

       ps%nnlcc = ps%ninner(0)
       ps%rnlcc = ps%rinner(0)

       if(compat7) then
          ps%rs_core_charge(:ab%npts)=otfg_core_pseudise_compat7(ps%rs_core_charge,ps,ab)
       else
          ps%rs_core_charge(:ab%npts)=otfg_core_pseudise(ps%rs_core_charge,ps,ab)
       end if

       call otfg_radial_transform(0,0,0,ps%mesh,ps%rlog,ps%rab,ps%num_points,ps%gmax,ps%rs_core_charge,ps%core_radial_charge)
       call otfg_derivative(ps%num_points,ps%gmax,ps%core_radial_charge,ps%dcore_radial_charge_dg)

    else

       ps%rs_core_charge            = 0.0_dp
       ps%core_radial_charge(:)     = 0.0_dp
       ps%dcore_radial_charge_dg(:) = 0.0_dp

    end if

    deallocate(rs_core_current,stat=ierr)
    if (ierr/=0) call io_abort('Error in deallocating rs_core_current in otfg_define_psp')

    temp_rc=temp_rc2 ! Keep a copy of default beta_rc to use for the gamma projectors

    ! ** Find out how many gamma functions there are -- so the arrays can be allocated
    !    - By default this will be the same as the number of beta projectors

    temp_def     = definition
    ps%num_gamma = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)
       nb = 0
       do i=1,len_trim(ctemp)
          if(ctemp(i:i)=='H') nb=nb+1 ! * Count the norm conserving projectors
          if(ctemp(i:i)=='G') nb=nb+1 ! * Count the ultrasoft projectors
       end do
       ps%num_gamma = ps%num_gamma + nb
    end do
    ps%num_gamma= ps%num_gamma +ps%num_beta ! Total gamma = beta + extra

    ps%bl_gmax=ps%gmax ! * As beta proj
    ps%bl_num_points=ps%num_points


    ! ** Allocate the ps%gamma_* arrays

    if (allocated(ps%gamma_l)) then
       deallocate(ps%gamma_l,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_l in otfg_define_psp')
    end if
    allocate(ps%gamma_l(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_l','otfg_define_psp')
    if (allocated(ps%gamma_e)) then
       deallocate(ps%gamma_e,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_e in otfg_define_psp')
    end if
    allocate(ps%gamma_e(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_e','otfg_define_psp')
    if (allocated(ps%gamma_rc)) then
       deallocate(ps%gamma_rc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_rc in otfg_define_psp')
    end if
    allocate(ps%gamma_rc(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_rc','otfg_define_psp')
    if (allocated(ps%gamma_nrc)) then
       deallocate(ps%gamma_nrc,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_nrc in otfg_define_psp')
    end if
    allocate(ps%gamma_nrc(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_nrc','otfg_define_psp')
    if (allocated(ps%gamma_norm)) then
       deallocate(ps%gamma_norm,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%gamma_norm in otfg_define_psp')
    end if
    allocate(ps%gamma_norm(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_norm','otfg_define_psp')

    if (allocated(ps%gamma_scheme)) then
       deallocate(ps%gamma_scheme,stat=ierr)
       if (ierr/=0) call io_abort('Error in deallocating ps%beta_scheme in otfg_define_psp')
    end if
    allocate(ps%gamma_scheme(1:ps%num_gamma),stat=ierr)
    if (ierr/=0) call io_allocate_abort('ps%gamma_scheme','otfg_define_psp')

    ps%gamma_l=0

    ! ** Set the ps%gamma_* arrays

    temp_def    = definition
    nbb         = 0
    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if
       temp_def = temp_def(indx+1:)

       read(ctemp(1:2),'(2i1)',err=103) no,l ! * Find out the nl for this reference state

       if(((ctemp(3:3)=='L').or.(ctemp(3:3)=='P')).and..not.(ctemp(4:4)=='N'.or.ctemp(4:4)=='U'&
            .or.ctemp(4:4)=='G'.or.ctemp(4:4)=='H')) then

          !This is the local or dummy channel, there are no gammas unless explicitly requested.

       else

          ps%gamma_e(nbb+1) = aeat%eval(no-l,l)

          nb = 0
          do i=1,len_trim(ctemp)
             if(ctemp(i:i)=='N'.or.ctemp(i:i)=='U'.or.ctemp(i:i)=='G'.or.ctemp(i:i)=='H') then
                nb=nb+1 ! * Count the projectors

                do j=i+1,len_trim(ctemp)
                   if(ctemp(j:j)=='N'.or.ctemp(j:j)=='U'.or.ctemp(j:j)=='G'.or.ctemp(j:j)=='H') exit
                end do
                ctemp2 = ctemp(i+1:j-1)

                indx = index(ctemp2,'@')
                if(indx>0) then
                   ctemp2 = ctemp2(:indx-1)
                end if

                indx = index(ctemp2,'+')
                jndx = index(ctemp2,'-')
                indx = max(indx,jndx)


                if(indx>0) then
                   read(ctemp2(indx:),*,err=103) dE_use
                   user=.true.
                   if(indx==1) then
                      ctemp2 = ' '
                   else
                      if(ctemp2(indx-1:indx-1)/='=') ctemp2 = ctemp2(:indx-1)
                   end if
                else
                   dE_use = dE*real(nb-1,dp)
                   user = .false.
                end if


                indx = index(ctemp2,'=')
                if(indx>0) then
                   read(ctemp2(indx+1:),*,err=103) dE_use
                   ctemp2 = ctemp2(:indx-1)
                   ps%gamma_e(nbb+nb) = dE_use
                else
                   if(user) then
                      ps%gamma_e(nbb+nb) = aeat%eval(no-l,l) + dE_use
                   else
                      if(dE_use>0.0_dp) then
                         ps%gamma_e(nbb+nb) = 0.25*dE_use/dE
                      else
                         ps%gamma_e(nbb+nb) = aeat%eval(no-l,l)
                      end if
                   end if
                end if

                if(len_trim(ctemp2)>0) then
                   read(ctemp2,*,err=103) ps%gamma_rc(nbb+nb)
                else
                   ps%gamma_rc(nbb+nb) = temp_rc
                end if

                ps%gamma_l(nbb+nb) = l

                if(ctemp(i:i)=='N'.or.ctemp(i:i)=='H') then
                   ps%gamma_norm(nbb+nb) = 1 ! * Norm conserving
                else
                   ps%gamma_norm(nbb+nb) = 0 ! * Ultrasoft
                end if

             end if
          end do


          ! * Didn't find any modifiers, so use defaults

          if(nb==0) then
             nb=2

             ps%gamma_l(nbb+1:nbb+2) = l

             ps%gamma_norm(nbb+1:nbb+2) = 0
             ps%gamma_rc(nbb+1:nbb+2)   = temp_rc

             ! * Set the energies

             ps%gamma_e(nbb+2) = 0.25

          end if !no modifiers

          nbb = nbb + nb

       end if !if local channel with gammas or not

    end do !loop over states

102 swap = .false.
    do nb=1,ps%num_gamma-1

       if(ps%gamma_l(nb)>ps%gamma_l(nb+1)) then
          swap = .true.
          l = ps%gamma_l(nb)
          ps%gamma_l(nb)   = ps%gamma_l(nb+1)
          ps%gamma_l(nb+1) = l
          l = ps%gamma_norm(nb)
          ps%gamma_norm(nb)   = ps%gamma_norm(nb+1)
          ps%gamma_norm(nb+1) = l
          temp_rc = ps%gamma_rc(nb)
          ps%gamma_rc(nb)   = ps%gamma_rc(nb+1)
          ps%gamma_rc(nb+1) = temp_rc
          temp_rc = ps%gamma_e(nb)
          ps%gamma_e(nb)   = ps%gamma_e(nb+1)
          ps%gamma_e(nb+1) = temp_rc
       end if

    end do

    if(swap) goto 102

    ! ** Find the _nrc and set the _rc to exact points

    do nb=1,ps%num_gamma
       ps%gamma_rc(nb) = otfg_round_to_grid(ps%gamma_rc(nb),ps%gamma_nrc(nb),ab%r)
    end do

    ! ** Do any charge augmentation?

    ps%aug = .true.

    if(all(ps%beta_norm==1)) ps%aug = .false.

    if(index(ps%flags,'aug')>0) ps%aug = .true.

    ! ** Is the pseudopotential "all-electron"?

    ps%scpsp = .false.

    if(index(ps%flags,'scpsp')>0) ps%scpsp = .true.

    ! ** Sort out the initialisation of the boundary conditions

    ! * First the grids and radii

    bc(nsp)%Rc    = max(maxval(ps%beta_rc),ps%local_rc)*1.1_dp  ! * 10% more than the max core radius
    bc(nsp)%Rc    = otfg_round_to_grid(bc(nsp)%Rc,n,ab%r)
    bc(nsp)%ngc   = ngrid/6

    dr = bc(nsp)%Rc/real(bc(nsp)%ngc-1,dp)

    do i=1,ngrid
       bc(nsp)%urg(i) = real(i-1,dp)*dr
    end do
    bc(nsp)%rmax = bc(nsp)%urg(ngrid)

    ! * The angular momentum channels

    do l=4,0,-1
       if(any(ps%pseudise(:,l))) exit
    end do

    bc(nsp)%lmax = l

    ! **  define data associated with core orbitals for ELNES

    temp_def = definition

    do n=1,nref
       indx = index(temp_def,':')
       if(indx>0) then
          ctemp = temp_def(1:indx-1)
       else
          ctemp = temp_def(1:)
       end if

       ! * Find out the principal quantum number label

       read(ctemp(1:1),*) nn

       if (n==1) then
          nnmax=nn
       else
          nnmax=max(nn,nnmax)
       end if

       temp_def = temp_def(indx+1:)

    end do

!!!!!!! TEMPORARY FIX :::
    nnmax=min(6,nnmax)

    select case(nnmax)
      case(1,2)
        ps%num_core_projectors = 1  ! 1s
        ps%num_core_radials = 1 ! 1s
        ps%core_radials_n(1:1) = (/1/) ! 1
        ps%core_radials_l(1:1) = (/0/) ! s
      case(3)
        ps%num_core_projectors = 5  ! 1s,2s,2px,2py,2pz
        ps%num_core_radials    = 3  ! 1s,2s,2p
        ps%core_radials_n(1:3) = (/1,2,1/)
        ps%core_radials_l(1:3) = (/0,0,1/)
      case(4)
        ps%num_core_projectors = 9
        ps%num_core_radials = 5      ! 1s,2s,2p,3s,3p
        ps%core_radials_n(1:5) = (/1,2,1,3,2/)
        ps%core_radials_l(1:5) = (/0,0,1,0,1/)
      case(5)
        ps%num_core_projectors = 18
        ps%num_core_radials = 8     ! 1s,2s,2p,3s,3p,3d,4s,4p
        ps%core_radials_n(1:8) = (/1,2,1,3,2,1,4,3/)
        ps%core_radials_l(1:8) = (/0,0,1,0,1,2,0,1/)
      case(6)
        ps%num_core_projectors = 27
        ps%num_core_radials = 11    ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p
        ps%core_radials_n(1:11) = (/1,2,1,3,2,1,4,3,2,5,4/)
        ps%core_radials_l(1:11) = (/0,0,1,0,1,2,0,1,2,0,1/)
      case(7)
        ps%num_core_projectors = 43
        ps%num_core_radials = 15   ! 1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d,6s,6p
        ps%core_radials_n(1:15) = (/1,2,1,3,2,1,4,3,2,1,5,4,3,6,5/)
        ps%core_radials_l(1:15) = (/0,0,1,0,1,2,0,1,2,3,0,1,2,0,1/)
      case default
        call io_abort('Error in otfg_define_psp: can not define core_radials')
    end select
    call trace_exit('otfg_define_psp',status)

    return

103 call io_abort('Error: otfg_define_psp: problem parsing definition')

  end subroutine otfg_core_charge


end module otfg_generate
