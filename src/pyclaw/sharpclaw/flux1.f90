module flux1_module
  use iso_c_binding
contains

  subroutine flux1(cptr,q1d,dq1d,dt,cfl,t,ixy,num_eqn,mx,num_ghost,maxmx,rpp) &
       bind(c, name='flux1')
    ! ===================================================================
    !
    !     # Evaluate (delta t) * dq(t)/dt
    !
    !     SharpClaw
    !     Author: David Ketcheson
    !
    !     amdq, apdq, amdq2, apdq2, wave, and s are used locally:
    !
    !     amdq(num_eqn,1-num_ghost:mx+num_ghost) = left-going flux-differences
    !     apdq(num_eqn,1-num_ghost:mx+num_ghost) = right-going flux-differences
    !        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
    !                         problem (between cells i-1 and i).
    !
    !     wave(1-num_ghost:mx+num_ghost, num_eqn, num_waves) = waves from solution of
    !                                           Riemann problems,
    !            wave(m,mw,i) = mth component of jump in q across
    !                           wave in family mw in Riemann problem between
    !                           states i-1 and i.
    !
    !     s(num_waves,1-num_ghost:mx+num_ghost) = wave speeds,
    !            s(mw,i) = speed of wave in family mw in Riemann problem between
    !                      states i-1 and i.
    !
    !     Note that mx must be the size of the patch for the dimension corresponding
    !     to the value of ixy.
    !
    !     t is the time at which we want to evaluate dq/dt, which may not
    !      be the current simulation time
    ! ===================================================================
    !
    ! Modified: April 26, 2011
    ! Authors:  David Ketcheson
    !           Matteo Parsani
    !
    ! ===================================================================

    use reconstruct
    use workspace
    use rp_acoustics
    implicit none

    type(c_ptr), intent(in), value :: cptr, rpp
    integer(c_int), intent(in), value :: num_eqn, num_ghost, maxmx, mx, ixy
    real(c_double), intent(in), value :: dt, t
    real(c_double), intent(in)  ::  q1d(num_eqn,1-num_ghost:mx+num_ghost)
    real(c_double), intent(out) :: dq1d(num_eqn,1-num_ghost:maxmx+num_ghost)
    real(c_double), intent(out) :: cfl

    type(wkspace), pointer :: wks
    type(wkblob), pointer :: wkb
    type(rp), pointer :: rpi

    integer :: i, m, mw

    call c_f_pointer(cptr, wks)
    call c_f_pointer(wks%bptr, wkb)
    call c_f_pointer(rpp, rpi)

    if (wks%index_capa.gt.0) then
       ! wkb%dtdx = dt / (wkb%dx(ixy)*aux(wks%index_capa,:))
    else
       wkb%dtdx = dt/wkb%dx(ixy)
    endif
    if (wks%num_dim.gt.1) dq1d=0.d0

    select case(wks%lim_type)
       ! Non-limited reconstruction of components of q (simplest approach)
       !        case(0)
       !        select case(char_decomp)
       !            case(0)
       !                call q2qlqr_poly(q1d,ql,qr,mx)
       !            case(1)
       !                ! wave-based unlimited reconstruction
       !                call rp1(maxmx,num_eqn,num_waves,num_ghost,mx,&
       !                        q1d,q1d,aux,aux,wave,s,amdq,apdq,num_aux)
       !                call q2qlqr_poly_wave(q1d,ql,qr,wave,s,mx)
       !        end select
    case(1)
       select case(wks%char_decomp)
       case(0)
          ! TVD reconstruction w/o char. decomp.
          call tvd2(q1d,wkb%ql,wkb%qr,wkb%mthlim)
       case(1)
          ! wave-based second order reconstruction
          ! call rp1(maxmx,num_eqn,wks%num_waves,num_ghost,mx,&
          !      q1d,q1d,aux,aux,wkb%wave,wkb%s,wkb%amdq,wkb%apdq,num_aux)
          ! Need to write a tvd2_fwave routine
          call tvd2_wave(wks,q1d,wkb%ql,wkb%qr,wkb%wave,wkb%s,wkb%mthlim)
       case(2)
          ! characteristic-wise second order reconstruction
          ! call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,wkb%evl,wkb%evr)
          ! call tvd2_char(wks,q1d,wkb%ql,wkb%qr,wkb%mthlim,wkb%evl,wkb%evr)
       end select
    case(2)
       select case (wks%char_decomp)
       case (0)
          ! no characteristic decomposition
          call weno_comp(wks,q1d,wkb%ql,wkb%qr,num_eqn,maxmx,num_ghost)
       case (1)
          ! wave-based reconstruction
          ! call rp1(maxmx,num_eqn,wks%num_waves,num_ghost,mx,&
          !      q1d,q1d,aux,aux,wkb%wave,wkb%s,wkb%amdq,wkb%apdq,num_aux)
          if (wks%fwave == 1) then
             call weno5_fwave(q1d,wkb%ql,wkb%qr,wkb%wave,wkb%s)
          else
             call weno5_wave(q1d,wkb%ql,wkb%qr,wkb%wave)
          endif
       case (2)
          ! characteristic-wise reconstruction
          ! call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,wkb%evl,wkb%evr)
          ! call weno5_char(wks,q1d,wkb%ql,wkb%qr,wkb%evl,wkb%evr)
       case (3)
          ! transmission-based reconstruction
          ! call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,wkb%evl,wkb%evr)
          ! call weno5_trans(wks,q1d,wkb%ql,wkb%qr,wkb%evl,wkb%evr)
       case default
          write(*,*) 'ERROR: Unrecognized characteristic decomposition option'
          write(*,*) 'You should set 0<=char_decomp<=3'
          stop
       end select
    case(3)
       call weno5(wks,q1d,wkb%ql,wkb%qr,num_eqn,maxmx,num_ghost)
    end select


    ! solve Riemann problem at each interface 
    ! -----------------------------------------
    call rp1(maxmx,num_eqn,wks%num_waves,num_ghost,mx,wkb%ql,wkb%qr,&
         wkb%wave,wkb%s,wkb%amdq,wkb%apdq,rpi%context)

    ! compute maximum wave speed:
    cfl = 0.d0
    do mw=1,wks%num_waves
       do i=1,mx+1
          ! if s>0 use wks%dtdx(i) to compute CFL,
          ! if s<0 use wks%dtdx(i-1) to compute CFL:
          cfl = dmax1(cfl, wkb%dtdx(i)*wkb%s(mw,i), -wkb%dtdx(i-1)*wkb%s(mw,i))
       enddo
    enddo

    ! Find total fluctuation within each cell
    if (wks%tfluct_solver == 1) then
       ! tfluct should be a special solver that uses the parameters aux(i)
       ! to solve a Riemann problem with left state ql(i)
       ! and right state qr(i), and returns a total fluctuation in wks%amdq2
       ! NOTE that here wks%amdq2 is really a total fluctuation (should be
       ! called adq); we do it this way just to avoid declaring more storage
       ! call tfluct(ixy,maxmx,num_eqn,wks%num_waves,num_ghost,mx,wkb%ql,wkb%qr, &
       !      aux,aux,wkb%s,wkb%amdq2)

       ! Modify q using fluctuations:
       ! Note this may not correspond to a conservative flux-differencing
       ! for equations not in conservation form.  It is conservative if
       ! adq = f(qr(i)) - f(ql(i)).

       forall (i=1:mx, m=1:num_eqn)
          dq1d(m,i) = dq1d(m,i) - wkb%dtdx(i)*(wkb%apdq(m,i) + &
               wkb%amdq2(m,i) + wkb%amdq(m,i+1))
       end forall

    else
       ! Or we can just swap things around and use the usual Riemann solver
       ! This may be more convenient, but is less efficient. 
       ! For the moment the swapping is done working with the element of the 
       ! vectors qr, ql, auxl, auxr. 
       ! 
       ! TODO: Working with pointers!!!

       do i = 1-num_ghost+1,mx+num_ghost
          do m = 1, num_eqn
             wkb%qr(m,i-1) = wkb%ql(m,i)
             wkb%ql(m,i  ) = wkb%qr(m,i)
          enddo
       enddo

       ! if (num_aux .gt. 0) then
       !    do i = 1-num_ghost+1,mx+num_ghost
       !       do m = 1, num_aux
       !          auxr(m,i-1) = aux(m,i) !aux is not patchdat type
       !          auxl(m,i  ) = aux(m,i) !aux is not patchdat type
       !       enddo
       !    enddo
       ! endif

       call rp1(maxmx,num_eqn,wks%num_waves,num_ghost,mx,wkb%ql,wkb%qr, &
            wkb%wave,wkb%s,wkb%amdq2,wkb%apdq2,rpi%context)

       forall(i=1:mx, m=1:num_eqn)
          dq1d(m,i) = dq1d(m,i)-wkb%dtdx(i)*(wkb%amdq(m,i+1)+ &
               wkb%apdq(m,i)+wkb%amdq2(m,i)+wkb%apdq2(m,i))
       end forall
    endif

  end subroutine flux1

end module flux1_module
