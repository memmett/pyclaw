module workspace

  use iso_c_binding

  type, bind(c) :: problemdata
     type(c_ptr) :: name
     type(c_ptr) :: data
  end type problemdata

  ! workspace (c/python binding - see libsharpclaw.py)
  type, bind(c) :: wkspace

     ! problem
     integer(c_int) :: num_eqn, num_dim, num_waves, index_capa
     
     ! grid
     integer(c_int) :: maxmx, num_ghost

     ! reconstruction/limiters
     integer(c_int) :: char_decomp, lim_type, multid_recon, weno_order
     integer(c_int) :: fwave, tfluct_solver
     real(c_double) :: epweno = 1.e-36

     ! pointer to workspace 'blob'
     type(c_ptr) :: bptr

  end type wkspace

  ! workspace blob
  type :: wkblob
     
     logical :: allocated = .false.

     ! problem
     real(c_double), pointer :: xlower(:), xupper(:), dx(:)

     ! reconstruction/limiters
     integer(c_int), pointer :: mthlim(:)

     real(c_double), pointer :: ql(:,:), qr(:,:), dtdx(:)
     real(c_double), pointer :: evl(:,:,:), evr(:,:,:)

     real(c_double), pointer :: dq1m(:)
     real(c_double), pointer :: uu(:,:), dq(:,:)
     real(c_double), pointer :: uh(:,:,:), gg(:,:), hh(:,:), u(:,:,:)

     ! waves
     real(c_double), pointer :: amdq(:,:), apdq(:,:), amdq2(:,:), apdq2(:,:)
     real(c_double), pointer :: wave(:,:,:), s(:,:)

  end type wkblob

contains

  ! subroutine get_problem_data_scalar(pdata, num_pdata, name, scalar) 
  !   implicit none

  !   type(c_ptr), intent(in)      :: pdata
  !   integer(c_int), intent(in)   :: num_pdata
  !   character(len=*), intent(in) :: name
  !   real(c_double), intent(out)  :: scalar

  !   integer :: i
  !   type(problemdata), pointer :: pd(:)
  !   character(len=32), pointer :: pdname ! XXX
  !   real(c_double), pointer    :: data

  !   call c_f_pointer(pdata, pd, [num_pdata])

  !   do i = 1, num_pdata
  !      call c_f_pointer(pd(i)%name, pdname, [pd(i)%len])

  !      if (name == pdname(:pd(i)%len)) then
  !         call c_f_pointer(pd(i)%data, data)
  !         scalar = data
  !         return
  !      end if
  !   end do

  !   stop 'Unable to find problem data.'

  ! end subroutine get_problem_data_scalar

  subroutine create_workspace(cptr) bind(c, name='create_workspace')
    implicit none
    type(c_ptr), intent(out) :: cptr
    type(wkspace), pointer :: wks

    allocate(wks)
    cptr = c_loc(wks)

    wks%maxmx = -1

  end subroutine create_workspace

  subroutine setup_workspace(cptr) bind(c, name='setup_workspace')

    type(c_ptr), intent(inout) :: cptr
    type(wkspace), pointer :: wks
    type(wkblob), pointer :: wkb

    call c_f_pointer(cptr, wks)

    allocate(wkb)
    wks%bptr = c_loc(wkb)

    allocate(wkb%ql(wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%qr(wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%amdq(wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%apdq(wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%amdq2(wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%apdq2(wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%wave(wks%num_eqn,num_waves,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%s(num_waves,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    allocate(wkb%dtdx(1-wks%num_ghost:wks%maxmx+wks%num_ghost))

    if (wks%char_decomp>1) then
       allocate(wkb%evl(wks%num_eqn,wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
       allocate(wkb%evr(wks%num_eqn,wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
    endif

    allocate(wkb%xlower(num_dim))
    allocate(wkb%xupper(num_dim))
    allocate(wkb%dx(num_dim))
    allocate(wkb%mthlim(num_waves))

    select case(wks%lim_type)
    case(1)
       select case(wks%char_decomp)
       case(1) ! Storage for tvd2_wave()
          allocate(wkb%uu(num_waves,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
       case(2) ! Storage for tvd2_char()
          ! Do the array bounds here cause a bug?
          allocate(wkb%dq(wks%num_eqn,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
          allocate( wkb%u(wks%num_eqn,2,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
          allocate(wkb%hh(-1:1,1-wks%num_ghost:wks%maxmx+wks%num_ghost))
       end select
    case(2)
       select case(wks%char_decomp)
       case(0)
          allocate(wkb%uu(2,wks%maxmx+2*wks%num_ghost))
          allocate(wkb%dq1m(wks%maxmx+2*wks%num_ghost))
       case(2) ! Storage for weno5_char
          allocate(wkb%dq(wks%num_eqn,wks%maxmx+2*wks%num_ghost))
          allocate(wkb%uu(2,wks%maxmx+2*wks%num_ghost))
          allocate(wkb%hh(-2:2,wks%maxmx+2*wks%num_ghost))
       case(3) ! Storage for weno5_trans
          allocate(wkb%dq(wks%num_eqn,wks%maxmx+2*wks%num_ghost))
          allocate(wkb%gg(wks%num_eqn,wks%maxmx+2*wks%num_ghost))
          allocate(wkb%u(wks%num_eqn,2,wks%maxmx+2*wks%num_ghost))
          allocate(wkb%hh(-2:2,wks%maxmx+2*wks%num_ghost))
          allocate(wkb%uh(wks%num_eqn,2,wks%maxmx+2*wks%num_ghost))
       end select
    case(3)
       allocate(wkb%uu(2,wks%maxmx+2*wks%num_ghost))
       allocate(wkb%dq1m(wks%maxmx+2*wks%num_ghost))
    end select

    wkb%allocated = .true.

  end subroutine setup_workspace

  subroutine set_dx(cptr, dx, num_dim) bind(c, name='set_dx')
    implicit none

    type(c_ptr), intent(in), value :: cptr
    integer(c_int), intent(in), value :: num_dim
    real(c_double), intent(in) :: dx(num_dim)
    type(wkspace), pointer :: wks
    type(wkblob), pointer :: wkb

    call c_f_pointer(cptr, wks)
    call c_f_pointer(wks%bptr, wkb)

    wkb%dx = dx

  end subroutine set_dx

  subroutine destroy_workspace(cptr) bind(c, name='destroy_workspace')

    type(c_ptr), intent(in) :: cptr
    type(wkspace), pointer :: wks
    type(wkblob), pointer :: wkb

    call c_f_pointer(cptr, wks)
    call c_f_pointer(wks%bptr, wkb)

    deallocate(wkb%ql)
    deallocate(wkb%qr)
    deallocate(wkb%amdq)
    deallocate(wkb%apdq)
    deallocate(wkb%amdq2)
    deallocate(wkb%apdq2)
    deallocate(wkb%wave)
    deallocate(wkb%s)
    deallocate(wkb%dtdx)

    if (wks%char_decomp>1) then
       deallocate(wkb%evl)
       deallocate(wkb%evr)
    endif

    wkb%allocated = .false.

    ! XXX

  end subroutine destroy_workspace

end module workspace
