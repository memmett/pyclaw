module flux1_module

  use iso_c_binding

contains

  subroutine flux1(cptr,q1d,dq1d,aux,dt,cfl,t,ixy,num_aux,num_eqn,mx,num_ghost,maxnx) &
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
    use workspace_module
    implicit none

    type(c_ptr), intent(in)     :: cptr
    integer(c_int), intent(in), value :: num_aux, num_eqn, num_ghost, maxnx, mx, ixy
    real(c_double), intent(in), value :: dt, t
    real(c_double), intent(in)  ::  q1d(num_eqn,1-num_ghost:mx+num_ghost)
    real(c_double), intent(out) :: dq1d(num_eqn,1-num_ghost:maxnx+num_ghost)
    real(c_double), intent(in)  ::   aux(num_aux,1-num_ghost:mx+num_ghost)

    real(c_double) :: auxl(num_aux,1-num_ghost:mx+num_ghost), auxr(num_aux,1-num_ghost:mx+num_ghost)
    real(c_double), intent(out) :: cfl

    type(workspace), pointer :: wks

    integer :: i, m, mw

    call c_f_pointer(cptr, wks)

    if (wks%index_capa.gt.0) then
       wks%dtdx = dt / (wks%dx(ixy)*aux(wks%index_capa,:))
    else
       wks%dtdx = dt/wks%dx(ixy)
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
       !                call rp1(maxnx,num_eqn,num_waves,num_ghost,mx,&
       !                        q1d,q1d,aux,aux,wave,s,amdq,apdq,num_aux)
       !                call q2qlqr_poly_wave(q1d,ql,qr,wave,s,mx)
       !        end select
    case(1)
       select case(wks%char_decomp)
       case(0)
          ! TVD reconstruction w/o char. decomp.
          call tvd2(q1d,wks%ql,wks%qr,wks%mthlim)
       case(1)
          ! wave-based second order reconstruction
          call rp1(maxnx,num_eqn,wks%num_waves,num_ghost,mx,&
               q1d,q1d,aux,aux,wks%wave,wks%s,wks%amdq,wks%apdq,num_aux)
          ! Need to write a tvd2_fwave routine
          call tvd2_wave(q1d,wks%ql,wks%qr,wks%wave,wks%s,wks%mthlim)
       case(2)
          ! characteristic-wise second order reconstruction
          call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,wks%evl,wks%evr)
          call tvd2_char(q1d,wks%ql,wks%qr,wks%mthlim,wks%evl,wks%evr)
       end select
    case(2)
       select case (wks%char_decomp)
       case (0)
          ! no characteristic decomposition
          call weno_comp(q1d,wks%ql,wks%qr,num_eqn,maxnx,num_ghost)
       case (1)
          ! wave-based reconstruction
          call rp1(maxnx,num_eqn,wks%num_waves,num_ghost,mx,&
               q1d,q1d,aux,aux,wks%wave,wks%s,wks%amdq,wks%apdq,num_aux)
          if (wks%fwave .eqv. .true.) then
             call weno5_fwave(q1d,wks%ql,wks%qr,wks%wave,wks%s)
          else
             call weno5_wave(q1d,wks%ql,wks%qr,wks%wave)
          endif
       case (2)
          ! characteristic-wise reconstruction
          call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,wks%evl,wks%evr)
          call weno5_char(q1d,wks%ql,wks%qr,wks%evl,wks%evr)
       case (3)
          ! transmission-based reconstruction
          call evec(mx,num_eqn,num_ghost,mx,q1d,aux,aux,wks%evl,wks%evr)
          call weno5_trans(q1d,wks%ql,wks%qr,wks%evl,wks%evr)
       case default
          write(*,*) 'ERROR: Unrecognized characteristic decomposition option'
          write(*,*) 'You should set 0<=char_decomp<=3'
          stop
       end select
    case(3)
       call weno5(q1d,wks%ql,wks%qr,num_eqn,maxnx,num_ghost)
    end select


    ! solve Riemann problem at each interface 
    ! -----------------------------------------
    call rp1(maxnx,num_eqn,wks%num_waves,num_ghost,mx,wks%ql,wks%qr,aux,aux, &
         wks%wave,wks%s,wks%amdq,wks%apdq,num_aux)

    ! compute maximum wave speed:
    cfl = 0.d0
    do mw=1,wks%num_waves
       do i=1,mx+1
          ! if s>0 use wks%dtdx(i) to compute CFL,
          ! if s<0 use wks%dtdx(i-1) to compute CFL:
          cfl = dmax1(cfl, wks%dtdx(i)*wks%s(mw,i), -wks%dtdx(i-1)*wks%s(mw,i))
       enddo
    enddo

    ! Find total fluctuation within each cell
    if (wks%tfluct_solver .eqv. .true.) then
       ! tfluct should be a special solver that uses the parameters aux(i)
       ! to solve a Riemann problem with left state ql(i)
       ! and right state qr(i), and returns a total fluctuation in wks%amdq2
       ! NOTE that here wks%amdq2 is really a total fluctuation (should be
       ! called adq); we do it this way just to avoid declaring more storage
       call tfluct(ixy,maxnx,num_eqn,wks%num_waves,num_ghost,mx,wks%ql,wks%qr, &
            aux,aux,wks%s,wks%amdq2)

       ! Modify q using fluctuations:
       ! Note this may not correspond to a conservative flux-differencing
       ! for equations not in conservation form.  It is conservative if
       ! adq = f(qr(i)) - f(ql(i)).

       forall (i=1:mx, m=1:num_eqn)
          dq1d(m,i) = dq1d(m,i) - wks%dtdx(i)*(wks%apdq(m,i) + &
               wks%amdq2(m,i) + wks%amdq(m,i+1))
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
             wks%qr(m,i-1) = wks%ql(m,i)
             wks%ql(m,i  ) = wks%qr(m,i)
          enddo
       enddo

       if (num_aux .gt. 0) then
          do i = 1-num_ghost+1,mx+num_ghost
             do m = 1, num_aux
                auxr(m,i-1) = aux(m,i) !aux is not patchdat type
                auxl(m,i  ) = aux(m,i) !aux is not patchdat type
             enddo
          enddo
       endif

       call rp1(maxnx,num_eqn,wks%num_waves,num_ghost,mx,wks%ql,wks%qr, &
            auxl,auxr,wks%wave,wks%s,wks%amdq2,wks%apdq2,num_aux)

       forall(i=1:mx, m=1:num_eqn)
          dq1d(m,i) = dq1d(m,i)-wks%dtdx(i)*(wks%amdq(m,i+1)+ &
               wks%apdq(m,i)+wks%amdq2(m,i)+wks%apdq2(m,i))
       end forall
    endif

  end subroutine flux1

end module flux1_module
