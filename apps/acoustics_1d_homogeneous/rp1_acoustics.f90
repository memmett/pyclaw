module rp_acoustics
contains

  ! Register this Riemann solver
  !
  ! Registers scalar problem data "zz" and "cc".
  subroutine rp_register(pdata_ptr, num_pdata) bind(c, name='rp_register')
    use workspace
    implicit none
    type(c_ptr), intent(out) :: pdata_ptr
    integer(c_int), intent(out) :: num_pdata

    type(problemdata), pointer :: pdata(:)
    character(len=12,kind=c_char), pointer :: names(:)
    real(c_double), pointer :: zz, cc

    allocate(pdata(2))
    allocate(names(2))
    allocate(zz,cc)

    names(1) = "zz" // c_null_char
    pdata(1)%name = c_loc(names(1)(1:1))
    pdata(1)%data = c_loc(zz)

    names(2) = "cc" // c_null_char
    pdata(2)%name = c_loc(names(2)(1:1))
    pdata(2)%data = c_loc(cc)

    num_pdata = 2
    pdata_ptr = c_loc(pdata(1))

  end subroutine rp_register


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
       wave,s,amdq,apdq,pdata,num_pdata)
    use iso_c_binding
    use workspace
    implicit none

    integer, intent(in) :: maxmx, mx, num_eqn, num_waves, num_ghost, num_pdata
    double precision, intent(out) :: &
         s(num_waves, 1-num_ghost:mx+num_ghost), &
         wave(num_eqn, num_waves, 1-num_ghost:mx+num_ghost), &
         apdq(num_eqn, 1-num_ghost:mx+num_ghost), &
         amdq(num_eqn, 1-num_ghost:mx+num_ghost)
    double precision, intent(in) :: &
         ql(num_eqn, 1-num_ghost:mx+num_ghost), &
         qr(num_eqn, 1-num_ghost:mx+num_ghost)
    type(c_ptr), intent(in) :: pdata
    type(problemdata), pointer :: pd(:)

    integer :: i, m
    double precision :: delta(2), a1, a2
    double precision, pointer :: zz, cc

    call c_f_pointer(pdata, pd, [num_pdata])
    call c_f_pointer(pd(1)%data, zz)
    call c_f_pointer(pd(2)%data, cc)

    ! split the jump in q at each interface into waves

    ! find a1 and a2, the coefficients of the 2 eigenvectors:
    do i = 2-num_ghost, mx+num_ghost
       delta(1) = ql(1,i) - qr(1,i-1)
       delta(2) = ql(2,i) - qr(2,i-1)
       a1 = (-delta(1) + zz*delta(2)) / (2.d0*zz)
       a2 =  (delta(1) + zz*delta(2)) / (2.d0*zz)

       ! compute the waves

       wave(1,1,i) = -a1*zz
       wave(2,1,i) = a1
       s(1,i) = -cc

       wave(1,2,i) = a2*zz
       wave(2,2,i) = a2
       s(2,i) = cc

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
