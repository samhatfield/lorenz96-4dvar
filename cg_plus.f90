subroutine cgbd(iprint, iter, nfun, gnorm, n, x, f, g, stp,finish, ndes, im, betafr, betapr, beta)
    !*********************************************************************
    !
    ! CGBD prints monitoring information.
    !
    !  Discussion:
    !
    !    The frequency and amount of output are controlled by IPRINT.
    !
    !  Modified:
    !
    !    18 December 2008
    !
    !  Author:
    !
    !    Jean Charles Gilbert, Jorge Nocedal.
    !
    !  Reference:
    !
    !    Jean Charles Gilbert, Jorge Nocedal,
    !    Global Convergence Properties of Conjugate Gradient Methods,
    !    SIAM Journal on Optimization,
    !    Volume 2, Number 1, 1992, pages 21-42.
    !
    !  Parameters
    !
    implicit none

    double precision beta
    double precision betafr
    double precision betapr
    double precision f
    logical finish
    double precision g(n)
    double precision gnorm
    integer i
    integer im
    integer iprint(2)
    integer iter
    integer n
    integer ndes
    integer nfun
    double precision stp
    double precision x(n)

    if (iter .eq. 0) then
        write (*, '(a)') ' '
        write (*, 10)
        write (*, 20) n
        write (*, 30) f, gnorm

        if (0 .lt. iprint(2)) then
            write (*, 40)
            write (*, 50) (x(i), i = 1, n)
            write (*, 60)
            write (*, 50) (g(i), i = 1, n)
        end if
        write (*, 10)
        write (*, 70)
    else
        if ((iprint(1).eq.0).and.(iter.ne.1.and..not.finish)) then
            return
        end if

        if (iprint(1).ne.0)then
            if (mod(iter-1,iprint(1)).eq.0.or.finish) then
                if (iprint(2).gt.1.and.iter.gt.1) write(*,70)
                write (*,80) iter,nfun,f,gnorm,stp,beta
            else
                return
            end if
        else
            if ( iprint(2).gt.1.and.finish) write(*,70)
            write (*,80) iter,nfun,f,gnorm,stp,beta
        end if

        if (iprint(2).eq.2.or.iprint(2).eq.3) then
            write (*,40)
            write (*,50)(x(i),i=1,n)
            if (iprint(2).eq.3) then
                write (*,60)
                write (*,50) (g(i),i=1,n)
            end if
        end if

        if (finish) then
            write (*,100)
        end if
    end if

 10 FORMAT('*************************************************')
 20 FORMAT(' N=',I5,//,'INITIAL VALUES:')
 30 FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)
 40 FORMAT(/,' VECTOR X= ')
 50 FORMAT(6(2X,1PD10.3/))
 60 FORMAT(' GRADIENT VECTOR G= ')
 70 FORMAT(/'   I  NFN',4X,'FUNC',7X,'GNORM',6X,'STEPLEN',4x,'BETA',/,&
        & ' ----------------------------------------------------')
 80 FORMAT(I4,1X,I3,2X,2(1PD10.3,2X),1PD8.1,2x,1PD8.1)
100 FORMAT(/' SUCCESSFUL CONVERGENCE (NO ERRORS).',/,' IFLAG = 0')
end subroutine  cgbd

subroutine cgfam(n, x, f, grad, srch_dir, gold, iprint, eps, w, iflag, irest, method, finish)
    ! CGFAM implements conjugate gradient methods for unconstrained nonlinear optimization.
    !
    !  Modified:
    !
    !    18 December 2008
    !
    !  Author:
    !
    !    Jean Charles Gilbert, Jorge Nocedal.
    !
    !  Reference:
    !
    !    Jean Charles Gilbert, Jorge Nocedal,
    !    Global Convergence Properties of Conjugate Gradient Methods,
    !    SIAM Journal on Optimization,
    !    Volume 2, Number 1, 1992, pages 21-42.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of variables.
    !
    !     x      =  iterate
    !     f      =  function value
    !     grad   =  gradient value
    !     gold   =  previous gradient value
    !
    !    input, integer iprint(2), controls printing.
    !    iprint(1) < 0 : no output is generated
    !    iprint(1) = 0 : output only at first and last iteration
    !    iprint(1) > 0 : output every iprint(1) iterations
    !    iprint(2)     : specifies the type of output generated;
    !                    the larger the value (between 0 and 3),
    !                    the more information
    !    iprint(2) = 0 : no additional information printed
    !    iprint(2) = 1 : initial x and gradient vectors printed
    !    iprint(2) = 2 : x vector printed every iteration
    !    iprint(2) = 3 : x vector and gradient vector printed
    !                       every iteration
    !
    !     eps    =  convergence constant
    !
    !     w      =  working array
    !
    !     iflag  =  controls termination of code, and return to main
    !               program to evaluate function and gradient
    !               iflag = -3 : improper input parameters
    !               iflag = -2 : descent was not obtained
    !               iflag = -1 : line search failure
    !               iflag =  0 : initial entry or
    !                            successful termination without error
    !               iflag =  1 : indicates a re-entry with new function values
    !               iflag =  2 : indicates a re-entry with a new iterate
    !
    !     irest  =  0 (no restarts); 1 (restart every n steps)
    !
    !     method =  1 : fletcher-reeves
    !               2 : polak-ribiere
    !               3 : positive polak-ribiere ( beta=max{beta,0} )

    implicit none

    integer n

    double precision srch_dir(n)
    double precision eps
    double precision f
    double precision grad(n)
    double precision gold(n)
    double precision w(n)
    double precision x(n)

    integer iprint(2),iflag,irest,method,im,ndes
    double precision gtol,one,zero,gnorm,ddot,stp1,ftol,xtol,stpmin,&
        & stpmax,stp,beta,betafr,betapr,dg0,gg,gg0,dgold,&
        & dgout,dg,dg1
    integer iter,nfun,maxfev,info,i,nfev,nrst,ides
    logical new,finish

    ! iter: keeps track of the number of iterations
    ! nfun: keeps track of the number of function/gradient evaluations
    common /runinf/iter,nfun
    save
    data one,zero/1.0d+0,0.0d+0/

    ! iflag = 2 indicates a re-entry with a new iterate
    if (iflag.eq.2) then
        ! call subroutine for printing output
        if (iprint(1).ge.0) then
            call cgbd(iprint,iter,nfun,gnorm,n,x,f,grad,stp,finish,ndes,im,betafr,betapr,beta)
        end if

        if (finish) then
            iflag = 0
            return
        end if
    end if

    if (iflag .eq. 0) then
        ! initialize
        !
        ! im =   number of times betapr was negative for method 2 or
        !        number of times betapr was 0 for method 3
        !eps
        ! ndes = number of line search iterations after wolfe conditions
        !        were satisfied
        iter= 0
        if (n.le.0) then
            iflag = -3
            write (*,140)
        end if
        nfun= 1
        new=.true.
        nrst= 0
        im=0
        ndes=0

        ! Set first search direction
        srch_dir = -grad

        gnorm = norm2(grad)
        stp1= one/gnorm

        ! parameters for line search routine
        !
        ! ftol and gtol are nonnegative input variables. termination
        !   occurs when the sufficient decrease condition and the
        !   directional derivative condition are satisfied.
        !
        ! xtol is a nonnegative input variable. termination occurs
        !   when the relative width of the interval of uncertainty
        !   is at most xtol.
        !
        ! stpmin and stpmax are nonnegative input variables which
        !   specify lower and upper bounds for the step.
        !
        ! maxfev is a positive integer input variable. termination
        !   occurs when the number of calls to fcn is at least
        !   maxfev by the end of an iteration.
        ftol= 1.0d-4
        gtol= 1.0d-1
        if (gtol.le.1.d-04) then
            write(*,145)
            gtol=1.d-02
        end if
        xtol= 1.0d-17
        stpmin= 1.0d-20
        stpmax= 1.0d+20
        maxfev= 40

        if (iprint(1).ge.0) then
            call cgbd(iprint,iter,nfun,gnorm,n,x,f,grad,stp,finish,ndes,im,betafr,betapr,beta)
        end if
    end if

    if (iflag .eq. 0 .or. iflag .eq. 2) then
        !     main iteration loop
        !
        iter= iter+1

        ! when nrst>n and irest=1 then restart
        !
        nrst= nrst+1
        info=0

        !  call the line search routine of mor'e and thuente
        !  (modified for our cg method)
        !
        !  Jorge More, David Thuente,
        !  Linesearch Algorithms with Guaranteed Sufficient Decrease,
        !  ACM Transactions on Mathematical Software,
        !  Volume 20, Number 3, September 1994, pages 286-307.
        !
        nfev=0

        gold = grad

        dg = dot_product(srch_dir,grad)
        dgold = dg
        stp = one

        ! shanno-phua's formula for trial step
        if (.not.new) stp = dg0/dg
        if (iter.eq.1) stp = stp1
        ides = 0
        new = .false.
    end if
72  continue

    ! call to the line search subroutine
    call cvsmod(n,x,f,grad,srch_dir,stp,ftol,gtol,xtol,stpmin,stpmax,maxfev,info,nfev,w,dg,dgout)
    ! info is an integer output variable set as follows:
    !   info = 0  improper input parameters.
    !   info =-1  a return is made to compute the function and gradient.
    !   info = 1  the sufficient decrease condition and the
    !             directional derivative condition hold.
    !   info = 2  relative width of the interval of uncertainty
    !             is at most xtol.
    !   info = 3  number of calls to fcn has reached maxfev.
    !   info = 4  the step is at the lower bound stpmin.
    !   info = 5  the step is at the upper bound stpmax.
    !   info = 6  rounding errors prevent further progress.
    !             there may not be a step which satisfies the
    !             sufficient decrease and curvature conditions.
    !             tolerances may be too small.
    if (info .eq. -1) then
        ! return to fetch function and gradient
        iflag=1
        return
    end if

    if (info .ne. 1) then
        iflag=-1
        write(*,100) info
        return
    end if

    ! test if descent direction is obtained for methods 2 and 3
    gg = dot_product(grad,grad)
    gg0 = dot_product(grad,gold)
    betapr = (gg - gg0) / gnorm**2

    if (irest.eq.1.and.nrst.gt.n) then
        nrst=0
        new=.true.
    else
        if (method.ne.1) then
            dg1 = -gg + betapr*dgout
            if (dg1.ge. 0.0d0) then
                if (iprint(1).ge.0) write(6,*) 'no descent'
                ides = ides + 1
                if (ides.gt.5) then
                    iflag=-2
                    write(*,135) i
                    return
                end if
                go to 72
            end if
        end if
    end if

    ! determine correct beta value for method chosen
    !
    ! im =   number of times betapr was negative for method 2 or
    !        number of times betapr was 0 for method 3
    !
    ! ndes = number of line search iterations after wolfe conditions
    !        were satisfied
    !
    nfun = nfun + nfev
    ndes = ndes + ides
    betafr = gg / gnorm**2

    if (nrst.eq.0) then
        beta= zero
    else
        if (method.eq.1) beta = betafr
        if (method.eq.2) beta = betapr
        if ((method.eq.2.or.method.eq.3).and.betapr.lt.0) im = im+1
        if (method.eq.3) beta = max(zero,betapr)
    end if

    !  Compute the new direction.
    srch_dir = -grad + beta*srch_dir
    dg0= dgold * stp

    ! return to driver for termination test
    gnorm = norm2(grad)
    iflag=2

    ! formats
    ! -------
100 format(/' iflag= -1 ',/' line search failed. see documentation of routine cvsmod'&
        & ,/' error return of line search: info= ',i2,/&
        & ' possible cause: function or gradient are incorrect')
135 format(/' iflag= -2',/' descent was not obtained')
140 format(/' iflag= -3',/' improper input parameters (n is not positive)')
145 format(/'  gtol is less than or equal to 1.d-04', / ' it has been reset to 1.d-02')
end subroutine cgfam

subroutine cstepm(stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, stpmin, stpmax, info)
    ! CSTEPM computes a safeguarded step for a line search.
    !
    !  Discussion:
    !
    !    The routine computes a safeguarded step for a line search, and
    !    updates an interval of uncertainty for a minimizer of the function.
    !
    !     the parameter stx contains the step with the least function
    !     value. the parameter stp contains the current step. it is
    !     assumed that the derivative at stx is negative in the
    !     direction of the step. if brackt is set true then a
    !     minimizer has been bracketed in an interval of uncertainty
    !     with endpoints stx and sty.
    !
    !     the subroutine statement is
    !
    !       subroutine cstepm(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
    !                        stpmin,stpmax,info)
    !
    !     where
    !
    !       stx, fx, and dx are variables which specify the step,
    !         the function, and the derivative at the best step obtained
    !         so far. the derivative must be negative in the direction
    !         of the step, that is, dx and stp-stx must have opposite
    !         signs. on output these parameters are updated appropriately.
    !
    !       sty, fy, and dy are variables which specify the step,
    !         the function, and the derivative at the other endpoint of
    !         the interval of uncertainty. on output these parameters are
    !         updated appropriately.
    !
    !       stp, fp, and dp are variables which specify the step,
    !         the function, and the derivative at the current step.
    !         if brackt is set true then on input stp must be
    !         between stx and sty. on output stp is set to the new step.
    !
    !       brackt is a logical variable which specifies if a minimizer
    !         has been bracketed. if the minimizer has not been bracketed
    !         then on input brackt must be set false. if the minimizer
    !         is bracketed then on output brackt is set true.
    !
    !       stpmin and stpmax are input variables which specify lower
    !         and upper bounds for the step.
    !
    !       info is an integer output variable set as follows:
    !         if info = 1,2,3,4,5, then the step has been computed
    !         according to one of the five cases below. otherwise
    !         info = 0, and this indicates improper input parameters.
    !
    !
    !     argonne national laboratory. minpack project. june 1983
    !     jorge j. more', david j. thuente

    implicit none

    integer info
    double precision stx
    double precision fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
    logical brackt,bound
    double precision gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta

    info = 0

    !  check the input parameters for errors.
    if ((brackt .and. (stp .le. min(stx,sty) .or. &
        & stp .ge. max(stx,sty))) .or. &
        & dx*(stp-stx) .ge. 0.0 .or. stpmax .lt. stpmin) return

    ! determine if the derivatives have opposite sign.
    sgnd = dp*(dx/abs(dx))

    ! first case. a higher function value.
    ! the minimum is bracketed. if the cubic step is closer
    ! to stx than the quadratic step, the cubic step is taken,
    ! else the average of the cubic and quadratic steps is taken.
    if (fp .gt. fx) then
        info = 1
        bound = .true.
        theta = 3*(fx - fp)/(stp - stx) + dx + dp
        s = max(abs(theta),abs(dx), abs(dp))
        gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
        if (stp .lt. stx) gamma = -gamma
        p = (gamma - dx) + theta
        q = ((gamma - dx) + gamma) + dp
        r = p/q
        stpc = stx + r*(stp - stx)
        stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
        if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
        else
           stpf = stpc + (stpq - stpc)/2
        end if
        brackt = .true.
    ! second case. a lower function value and derivatives of
    ! opposite sign. the minimum is bracketed. if the cubic
    ! step is closer to stx than the quadratic (secant) step,
    ! the cubic step is taken, else the quadratic step is taken.
    else if (sgnd .lt. 0.0) then
        info = 2
        bound = .false.
        theta = 3*(fx - fp)/(stp - stx) + dx + dp
        s = max(abs(theta), abs(dx), abs(dp))
        gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
        if (stp .gt. stx) gamma = -gamma
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dx
        r = p/q
        stpc = stp + r*(stx - stp)
        stpq = stp + (dp/(dp-dx))*(stx - stp)
        if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
        else
            stpf = stpq
        end if
        brackt = .true.
    ! third case. a lower function value, derivatives of the
    ! same sign, and the magnitude of the derivative decreases.
    ! the cubic step is only used if the cubic tends to infinity
    ! in the direction of the step or if the minimum of the cubic
    ! is beyond stp. otherwise the cubic step is defined to be
    ! either stpmin or stpmax. the quadratic (secant) step is also
    ! computed and if the minimum is bracketed then the the step
    ! closest to stx is taken, else the step farthest away is taken.
    else if (abs(dp) .lt. abs(dx)) then
        info = 3
        bound = .true.
        theta = 3*(fx - fp)/(stp - stx) + dx + dp
        s = max(abs(theta), abs(dx), abs(dp))
        ! the case gamma = 0 only arises if the cubic does not tend
        ! to infinity in the direction of the step.
         gamma = s*sqrt(max(0.0d0,(theta/s)**2 - (dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. 0.0 .and. gamma .ne. 0.0) then
             stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
             stpc = stpmax
         else
             stpc = stpmin
         end if
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (brackt) then
             if (abs(stp-stpc) .lt. abs(stp-stpq)) then
                 stpf = stpc
             else
                 stpf = stpq
             end if
         else
             if (abs(stp-stpc) .gt. abs(stp-stpq)) then
                 stpf = stpc
             else
                 stpf = stpq
             end if
         end if
    ! fourth case. a lower function value, derivatives of the
    ! same sign, and the magnitude of the derivative does
    ! not decrease. if the minimum is not bracketed, the step
    ! is either stpmin or stpmax, else the cubic step is taken.
    else
        info = 4
        bound = .false.
        if (brackt) then
            theta = 3*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta), abs(dy), abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
        else if (stp .gt. stx) then
            stpf = stpmax
        else
            stpf = stpmin
        end if
    end if

    ! update the interval of uncertainty. this update does not
    ! depend on the new step or the case analysis above.
    if (fp .gt. fx) then
        sty = stp
        fy = fp
        dy = dp
    else
        if (sgnd .lt. 0.0) then
            sty = stx
            fy = fx
            dy = dx
        end if
        stx = stp
        fx = fp
        dx = dp
    end if

    ! compute the new step and safeguard it.
    stpf = min(stpmax,stpf)
    stpf = max(stpmin,stpf)
    stp = stpf
    if (brackt .and. bound) then
        if (sty .gt. stx) then
            stp = min(stx+0.66*(sty-stx), stp)
        else
           stp = max(stx+0.66*(sty-stx), stp)
         end if
    end if
end subroutine cstepm

subroutine cvsmod(n, x, f, g, s, stp, ftol, gtol, xtol, stpmin, stpmax, maxfev, info, nfev, wa, dginit, dgout)
    ! CSVMOD finds a step which satisfies a sufficient decrease condition.
    !
    !  Discussion:
    !
    !    The routine finds a step which satisfies a sufficient decrease condition
    !    and a curvature condition.
    !
    !     the user must provide a subroutine which calculates the
    !     function and the gradient.
    !
    !     at each stage the subroutine updates an interval of
    !     uncertainty with endpoints stx and sty. the interval of
    !     uncertainty is initially chosen so that it contains a
    !     minimizer of the modified function
    !
    !          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
    !
    !     if a step is obtained for which the modified function
    !     has a nonpositive function value and nonnegative derivative,
    !     then the interval of uncertainty is chosen so that it
    !     contains a minimizer of f(x+stp*s).
    !
    !     the algorithm is designed to find a step which satisfies
    !     the sufficient decrease condition
    !
    !           f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),
    !
    !     and the curvature condition
    !
    !           abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).
    !
    !     if ftol is less than gtol and if, for example, the function
    !     is bounded below, then there is always a step which satisfies
    !     both conditions. if no step can be found which satisfies both
    !     conditions, then the algorithm usually stops when rounding
    !     errors prevent further progress. in this case stp only
    !     satisfies the sufficient decrease condition.
    !
    !     the subroutine statement is
    !
    !        subroutine cvsmod(n,x,f,g,s,stp,ftol,gtol,xtol,
    !                   stpmin,stpmax,maxfev,info,nfev,wa,dg,dgout)
    !     where
    !
    !       n is a positive integer input variable set to the number
    !         of variables.
    !
    !       x is an array of length n. on input it must contain the
    !         base point for the line search. on output it contains
    !         x + stp*s.
    !
    !       f is a variable. on input it must contain the value of f
    !         at x. on output it contains the value of f at x + stp*s.
    !
    !       g is an array of length n. on input it must contain the
    !         gradient of f at x. on output it contains the gradient
    !         of f at x + stp*s.
    !
    !       s is an input array of length n which specifies the
    !         search direction.
    !
    !       stp is a nonnegative variable. on input stp contains an
    !         initial estimate of a satisfactory step. on output
    !         stp contains the final estimate.
    !
    !       ftol and gtol are nonnegative input variables. termination
    !         occurs when the sufficient decrease condition and the
    !         directional derivative condition are satisfied.
    !
    !       xtol is a nonnegative input variable. termination occurs
    !         when the relative width of the interval of uncertainty
    !         is at most xtol.
    !
    !       stpmin and stpmax are nonnegative input variables which
    !         specify lower and upper bounds for the step.
    !
    !       maxfev is a positive integer input variable. termination
    !         occurs when the number of calls to fcn is at least
    !         maxfev by the end of an iteration.
    !
    !       info is an integer output variable set as follows:
    !
    !         info = 0  improper input parameters.
    !
    !         info =-1  a return is made to compute the function and gradient.
    !
    !         info = 1  the sufficient decrease condition and the
    !                   directional derivative condition hold.
    !
    !         info = 2  relative width of the interval of uncertainty
    !                   is at most xtol.
    !
    !         info = 3  number of calls to fcn has reached maxfev.
    !
    !         info = 4  the step is at the lower bound stpmin.
    !
    !         info = 5  the step is at the upper bound stpmax.
    !
    !         info = 6  rounding errors prevent further progress.
    !                   there may not be a step which satisfies the
    !                   sufficient decrease and curvature conditions.
    !                   tolerances may be too small.
    !
    !       nfev is an integer output variable set to the number of
    !         calls to fcn.
    !
    !       wa is a work array of length n.
    !
    !       *** the following two parameters are a modification to the code
    !
    !       dg is the initial directional derivative (in the original code
    !                 it was computed in this routine0
    !
    !       dgout is the value of the directional derivative when the wolfe
    !             conditions hold, and an exit is made to check descent.
    !
    !     subprograms called
    !
    !       cstepm
    !
    !       fortran-supplied...abs,max,min
    !
    !     argonne national laboratory. minpack project. june 1983
    !     jorge j. more', david j. thuente

    implicit none

    integer n

    integer info
    integer maxfev
    integer nfev
    double precision f,stp,ftol,gtol,xtol,stpmin,stpmax
    double precision x(n),g(n),s(n),wa(n)

    save

    integer infoc,j
    logical brackt,stage1
    double precision dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym, &
        & finit,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty, &
        & stmin,stmax,width,width1,xtrapf,zero,dgout

    data p5,p66,xtrapf,zero /0.5d0,0.66d0,4.0d0,0.0d0/

    if (info.eq.-1) go to 45
    if (info.eq.1) go to 321
    infoc = 1

    ! check the input parameters for errors.
    if (n .le. 0 .or. stp .le. zero .or. ftol .lt. zero .or. &
        & gtol .lt. zero .or. xtol .lt. zero .or. stpmin .lt. zero &
        & .or. stpmax .lt. stpmin .or. maxfev .le. 0) return

    ! compute the initial gradient in the search direction
    ! and check that s is a descent direction.
    if (zero .le. dginit) then
        return
    end if

    ! initialize local variables.
    brackt = .false.
    stage1 = .true.
    nfev = 0
    finit = f
    dgtest = ftol * dginit
    width = stpmax - stpmin
    width1 = width / p5
    do j = 1, n
        wa(j) = x(j)
    end do

    ! the variables stx, fx, dgx contain the values of the step,
    ! function, and directional derivative at the best step.
    ! the variables sty, fy, dgy contain the value of the step,
    ! function, and derivative at the other endpoint of
    ! the interval of uncertainty.
    ! the variables stp, f, dg contain the values of the step,
    ! function, and derivative at the current step.
    stx = zero
    fx = finit
    dgx = dginit
    sty = zero
    fy = finit
    dgy = dginit

    ! start of iteration.

   30 continue
    ! set the minimum and maximum steps to correspond
    ! to the present interval of uncertainty.
    if (brackt) then
        stmin = min(stx,sty)
        stmax = max(stx,sty)
    else
        stmin = stx
        stmax = stp + xtrapf*(stp - stx)
    end if

    ! force the step to be within the bounds stpmax and stpmin.
    stp = max(stp,stpmin)
    stp = min(stp,stpmax)

    ! if an unusual termination is to occur then let
    ! stp be the lowest point obtained so far.
    if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax)) &
        & .or. nfev .ge. maxfev-1 .or. infoc .eq. 0 &
        & .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx

    ! evaluate the function and gradient at stp
    ! and compute the directional derivative.
    do j = 1, n
        x(j) = wa(j) + stp*s(j)
    end do

    ! return to compute function value
    info=-1
    return

   45    continue

    info=0
    nfev = nfev + 1
    dg = zero
    do j = 1, n
        dg = dg + g(j)*s(j)
    end do
    ftest1 = finit + stp*dgtest

    ! test for convergence.
    if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax)) .or. infoc .eq. 0) info = 6
    if (stp .eq. stpmax .and. f .le. ftest1 .and. dg .le. dgtest) info = 5
    if (stp .eq. stpmin .and. (f .gt. ftest1 .or. dg .ge. dgtest)) info = 4
    if (nfev .ge. maxfev) info = 3
    if (brackt .and. stmax-stmin .le. xtol*stmax) info = 2

    ! more's code has been modified so that at least one new
    ! function value is computed during the line search (enforcing
    ! at least one interpolation is not easy, since the code may
    ! override an interpolation)

    if (f .le. ftest1 .and. abs(dg) .le. gtol*(-dginit).and.nfev.gt.1) info = 1
    ! check for termination.
    if (info .ne. 0)then
        dgout=dg
        return
    end if

 321     continue
    ! in the first stage we seek a step for which the modified
    ! function has a nonpositive value and nonnegative derivative.
    if (stage1 .and. f .le. ftest1 .and. dg .ge. min(ftol,gtol)*dginit) stage1 = .false.

    ! a modified function is used to predict the step only if
    ! we have not obtained a step for which the modified
    ! function has a nonpositive function value and nonnegative
    ! derivative, and if a lower function value has been
    ! obtained but the decrease is not sufficient.
    if (stage1 .and. f .le. fx .and. f .gt. ftest1) then
        ! define the modified function and derivative values.
        fm = f - stp*dgtest
        fxm = fx - stx*dgtest
        fym = fy - sty*dgtest
        dgm = dg - dgtest
        dgxm = dgx - dgtest
        dgym = dgy - dgtest
        ! call cstepm to update the interval of uncertainty
        ! and to compute the new step.
        call cstepm(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,brackt,stmin,stmax,infoc)

        ! reset the function and gradient values for f.
        fx = fxm + stx*dgtest
        fy = fym + sty*dgtest
        dgx = dgxm + dgtest
        dgy = dgym + dgtest
    else
        ! call cstepm to update the interval of uncertainty
        ! and to compute the new step.
        call cstepm(stx,fx,dgx,sty,fy,dgy,stp,f,dg,brackt,stmin,stmax,infoc)
    end if

    ! force a sufficient decrease in the size of the
    ! interval of uncertainty.
    if (brackt) then
        if (abs(sty-stx) .ge. p66*width1) then
            stp = stx + p5*(sty - stx)
        end if
        width1 = width
        width = abs(sty-stx)
    end if

    go to 30
end subroutine cvsmod
