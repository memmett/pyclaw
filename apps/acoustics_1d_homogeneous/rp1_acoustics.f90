module rp_acoustics
  use iso_c_binding
  
  type :: acoustics
     real(c_double) :: z, c
  end type acoustics

contains

  ! Register this Riemann solver
  !
  ! Requests scalar problem data "zz" and "cc".
  subroutine rp_register(rpp) bind(c, name='rp_register')
    use workspace
    implicit none
    type(c_ptr), intent(in), value :: rpp
    type(acoustics), pointer :: ac

    allocate(ac)
    ac%z = 0.0d0
    ac%c = 0.0d0

    call rp_set_context(rpp, c_loc(ac))
    call rp_options_real(rpp, "zz", ac%z)
    call rp_options_real(rpp, "cc", ac%c)
  end subroutine rp_register

  subroutine rp_teardown(rpp) bind(c, name='rp_teardown')
    use workspace
    implicit none
    type(c_ptr), intent(in), value :: rpp
    type(rp), pointer :: rpi
    type(acoustics), pointer :: ac

    call c_f_pointer(rpp, rpi)
    call c_f_pointer(rpi%context, ac)

    deallocate(ac)
  end subroutine rp_teardown

  ! Riemann solver for the acoustics equations in 1d
  !
  ! On input, ql contains the state vector at the left edge of each cell
  !           qr contains the state vector at the right edge of each cell
  !
  ! On output, wave contains the waves,
  !            s the speeds,
  !
  !            amdq = A^- Delta q,
  !            apdq = A^+ Delta q,
  !                   the decomposition of the flux difference
  !                       f(qr(i-1)) - f(ql(i))
  !                   into leftgoing and rightgoing parts respectively.
  !
  !
  ! Note that the i'th Riemann problem has left state qr(i-1,:) and
  !                                       right state ql(i,:)
  !
  ! From the basic clawpack routines, this routine is called with ql = qr
  subroutine rp1(maxmx,num_eqn,num_waves,num_ghost,mx,ql,qr,&
       wave,s,amdq,apdq,context)
    use iso_c_binding
    use workspace
    implicit none

    integer, intent(in) :: maxmx, mx, num_eqn, num_waves, num_ghost
    double precision, intent(out) :: &
         s(num_waves, 1-num_ghost:mx+num_ghost), &
         wave(num_eqn, num_waves, 1-num_ghost:mx+num_ghost), &
         apdq(num_eqn, 1-num_ghost:mx+num_ghost), &
         amdq(num_eqn, 1-num_ghost:mx+num_ghost)
    double precision, intent(in) :: &
         ql(num_eqn, 1-num_ghost:mx+num_ghost), &
         qr(num_eqn, 1-num_ghost:mx+num_ghost)
    type(c_ptr), intent(in), value :: context
    type(acoustics), pointer :: ac

    integer :: i, m
    double precision :: delta(2), a1, a2

    call c_f_pointer(context, ac)

    ! split the jump in q at each interface into waves

    ! find a1 and a2, the coefficients of the 2 eigenvectors:
    do i = 2-num_ghost, mx+num_ghost
       delta(1) = ql(1,i) - qr(1,i-1)
       delta(2) = ql(2,i) - qr(2,i-1)
       a1 = (-delta(1) + ac%z*delta(2)) / (2.d0*ac%z)
       a2 =  (delta(1) + ac%z*delta(2)) / (2.d0*ac%z)

       ! compute the waves

       wave(1,1,i) = -a1*ac%z
       wave(2,1,i) = a1
       s(1,i) = -ac%c

       wave(1,2,i) = a2*ac%z
       wave(2,2,i) = a2
       s(2,i) = ac%c

    end do

    ! compute the leftgoing and rightgoing flux differences:
    ! note s(1,i) < 0   and   s(2,i) > 0

    do m=1,num_eqn
       do i = 2-num_ghost, mx+num_ghost
          amdq(m,i) = s(1,i)*wave(m,1,i)
          apdq(m,i) = s(2,i)*wave(m,2,i)
       end do
    end do

  end subroutine rp1
end module rp_acoustics
