module workspace_module

  use iso_c_binding

  type :: workspace
     
     logical :: work_alloc = .false.

     integer(c_int) :: maxnx, num_ghost
     integer(c_int) :: num_eqn, num_dim, num_waves, index_capa
     integer(c_int) :: char_decomp, lim_type, multid_recon, weno_order
     logical :: fwave, tfluct_solver

     real(c_double), pointer :: ql(:,:), qr(:,:), dtdx(:)
     real(c_double), pointer :: evl(:,:,:), evr(:,:,:)
     real(c_double), pointer :: amdq(:,:), apdq(:,:), amdq2(:,:), apdq2(:,:)
     real(c_double), pointer :: wave(:,:,:), s(:,:)

     real(c_double), pointer :: xlower(:), xupper(:), dx(:)
     integer(c_int), pointer :: mthlim(:)

     ! reconstruction
     real(c_double), pointer :: dq1m(:)
     real(c_double), pointer :: uu(:,:), dq(:,:)
     real(c_double), pointer :: uh(:,:,:), gg(:,:), hh(:,:), u(:,:,:)
     real(c_double) :: epweno = 1.e-36

  end type workspace

contains
  
  subroutine alloc_workspace(cptr,maxnx,num_ghost,num_dim,num_eqn,num_waves,char_decomp,lim_type) &
       bind(c, name='alloc_workspace')

    type(c_ptr), intent(out) :: cptr
    integer(c_int), intent(in), value :: maxnx, num_ghost 
    integer(c_int), intent(in), value :: num_dim, num_eqn, num_waves
    integer(c_int), intent(in), value :: char_decomp, lim_type
    type(workspace), pointer :: wkspace

    allocate(wkspace)
    cptr = c_loc(wkspace)

    wkspace%maxnx = maxnx
    wkspace%num_ghost = num_ghost
    wkspace%num_dim = num_dim
    wkspace%num_eqn = num_eqn
    wkspace%num_waves = num_waves
    wkspace%char_decomp = char_decomp
    wkspace%lim_type = lim_type

    allocate(wkspace%ql(num_eqn,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%qr(num_eqn,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%amdq(num_eqn,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%apdq(num_eqn,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%amdq2(num_eqn,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%apdq2(num_eqn,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%wave(num_eqn,num_waves,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%s(num_waves,1-num_ghost:maxnx+num_ghost))
    allocate(wkspace%dtdx(1-num_ghost:maxnx+num_ghost))

    if (wkspace%char_decomp>1) then
       allocate(wkspace%evl(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost))
       allocate(wkspace%evr(num_eqn,num_eqn,1-num_ghost:maxnx+num_ghost))
    endif

    allocate(wkspace%xlower(num_dim))
    allocate(wkspace%xupper(num_dim))
    allocate(wkspace%dx(num_dim))
    allocate(wkspace%mthlim(num_waves))

    select case(lim_type)
    case(1)
       select case(char_decomp)
       case(1) ! Storage for tvd2_wave()
          allocate(wkspace%uu(num_waves,1-num_ghost:maxnx+num_ghost))
       case(2) ! Storage for tvd2_char()
          ! Do the array bounds here cause a bug?
          allocate(wkspace%dq(num_eqn,1-num_ghost:maxnx+num_ghost))
          allocate( wkspace%u(num_eqn,2,1-num_ghost:maxnx+num_ghost))
          allocate(wkspace%hh(-1:1,1-num_ghost:maxnx+num_ghost))
       end select
    case(2)
       select case(char_decomp)
       case(0)
          allocate(wkspace%uu(2,maxnx+2*num_ghost))
          allocate(wkspace%dq1m(maxnx+2*num_ghost))
       case(2) ! Storage for weno5_char
          allocate(wkspace%dq(num_eqn,maxnx+2*num_ghost))
          allocate(wkspace%uu(2,maxnx+2*num_ghost))
          allocate(wkspace%hh(-2:2,maxnx+2*num_ghost))
       case(3) ! Storage for weno5_trans
          allocate(wkspace%dq(num_eqn,maxnx+2*num_ghost))
          allocate(wkspace%gg(num_eqn,maxnx+2*num_ghost))
          allocate(wkspace%u(num_eqn,2,maxnx+2*num_ghost))
          allocate(wkspace%hh(-2:2,maxnx+2*num_ghost))
          allocate(wkspace%uh(num_eqn,2,maxnx+2*num_ghost))
       end select
    case(3)
       allocate(wkspace%uu(2,maxnx+2*num_ghost))
       allocate(wkspace%dq1m(maxnx+2*num_ghost))
    end select

    wkspace%work_alloc = .true.

  end subroutine alloc_workspace

  subroutine dealloc_workspace(cptr) bind(c, name='dealloc_workspace')

    type(c_ptr), intent(in) :: cptr
    type(workspace), pointer :: wkspace

    call c_f_pointer(cptr, wkspace)

    deallocate(wkspace%ql)
    deallocate(wkspace%qr)
    deallocate(wkspace%amdq)
    deallocate(wkspace%apdq)
    deallocate(wkspace%amdq2)
    deallocate(wkspace%apdq2)
    deallocate(wkspace%wave)
    deallocate(wkspace%s)
    deallocate(wkspace%dtdx)

    if (wkspace%char_decomp>1) then
       deallocate(wkspace%evl)
       deallocate(wkspace%evr)
    endif

    wkspace%work_alloc = .false.

  end subroutine dealloc_workspace

end module workspace_module
