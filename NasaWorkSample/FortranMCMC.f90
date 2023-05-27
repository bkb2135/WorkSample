module  transits
    contains
        real function integrate_simpson(f, a, b, n) result(integral)
            real, intent(in) :: a, b
            integer, intent(in) :: n
            real :: h, x, sum
            integer :: i

            h = (b-a)/n
            sum = f(a) + f(b)
            do i=1, n-1, 2
                x = a + i*h
                sum = sum + 4.0*f(x)
            end do
            do i=2, n-2, 2
                x = a + i*h
                sum = sum + 2.0*f(x)
            end do
            integral = sum*h/3.0
        end function integrate_simpson

        real function integrator(f, p, z)
        real, intent(in) :: p, z
        real :: b, a, sum
        integer :: i, n
        b = 1.0
        a = 0.0
        n = 1000
        h = (b-a)/n
        sum = f(a, p, z) + f(b, p, z)   
        do i=1, n-1, 2
            x = a + i*h
            sum = sum + 4.0*f(x, p, z)
        end do
        do i=2, n-2, 2
            x = a + i*h
            sum = sum + 2.0*f(x, p, z)
        end do
        integrator = sum*h/3.0
        end function integrator

        real function lnprob_log(time, x, y, y_err) result(log_probability)
        real,dimension(2157), intent(in) :: y, y_err, time
        real, dimension(2157) :: flux_prediction
        real, dimension(2), intent(in) :: x
        real :: p, tau

        p = EXP(x(1))
        tau = EXP(x(2))
        flux_prediction = compute_light_curve(time, p, tau)
        log_probability = -0.5*compute_chi_squared(y, flux_prediction, y_err)

        end function lnprob_log

        real function compute_chi_squared(y, pred, y_err) result(chi_squared)
        real, dimension(2157) :: y, pred, y_err
        real :: sum = 0
        integer :: i
        do i = 1, 2157
            sum  = sum+ ((y(i)-pred(i))/pred(i))**2
        end do
        chi_squared = sum
        end function compute_chi_squared

        real function delta(r, p, z)
            real, intent(in) :: p, r, z
            real, parameter :: PI=4.0*ATAN(1.0)

            IF (r >= z + p .OR. r <= z-p)THEN
                delta = 0.0
            ELSE IF (r + z <= p) THEN
                delta = 1.0
            ELSE 
                delta = (1/PI)*ACOS((z*z - p*p + r*r) / (2.0*z*r))
            END IF
        end function delta

        real function limb_darkening(r) result(ld)
            real, intent(in) :: r
            real, parameter :: u1=0.3, u2=0.2
            real :: mu
            mu = sqrt(1.0 - r*r)
            ld = 1.0 - u1*(1.0-mu) - u2*(1.0-mu)*(1.0-mu)
        end function limb_darkening

        real function denominator(r)
            real, intent(in) :: r
            denominator = limb_darkening(r)*2*r
        end function denominator

        real function numerator(r, p, z) 
        numerator = limb_darkening(r)*(1-delta(r, p, abs(z)))*2*r
        end function numerator

        real function integrate_numerator(f, p, z) result(flux_obscured)
            real, intent(in) :: p, z
            real :: b, a, sum
            integer :: i, n
            b = 1.0
            a = 0.0 
            n = 1000
            h = (b-a)/n
            sum = f(a, p, z) + f(b, p, z)
            do i=1, n-1, 2
                x = a + i*h
                sum = sum + 4.0*f(x, p, z)
            end do
            do i=2, n-2, 2
                x = a + i*h
                sum = sum + 2.0*f(x, p, z)
            end do
            integral = sum*h/3.0
            flux_obscured = integral
        end function integrate_numerator

        function compute_light_curve(time, p, tau) result(numerical_flux)
        real, dimension(2157), intent(in) :: time
        real, intent(in) :: p, tau
        real, dimension(2157) :: numerical_flux, z_arr
        z_arr = (time)/tau
        do i=1, 2157
            flux_obscured = integrator(numerator, p, z_arr(i))
            flux_unobscured = integrator(denominator, p, z_arr(i))
            numerical_flux(i) = flux_obscured/flux_unobscured
        end do
        !compute_light_curve = numerical_flux
        end function compute_light_curve

        function mcmc(starting_point, nsteps, step_size, time, flux, flux_err) result(chain)
            real, dimension(2), intent(in) :: starting_point
            integer, intent(in) :: nsteps
            real, intent(in):: step_size
            real, dimension(2157), intent(in) :: time, flux, flux_err
            real, dimension(2, nsteps) :: chain
            real, dimension(2, nsteps) :: prob_chain
            real :: log_prob, log_prob_new
            real, dimension(2) :: x_new, x_old
            integer :: i
            call srand(0)
            x_old = starting_point
            log_prob = lnprob_log(time, x_old, flux, flux_err)
            do i=1, nsteps
                x_new(1) = x_old(1) + step_size*rand() - step_size*0.5
                x_new(2) = x_old(2) + step_size*rand() - step_size*0.5
                log_prob_new = lnprob_log(time, x_new, flux, flux_err)
                if (log_prob_new > log_prob) then
                    chain(:,i) = x_new
                    x_old = x_new
                    log_prob = log_prob_new
                else if (log_prob - log_prob_new < 2) then
                    chain(:,i) = x_new
                    x_old = x_new
                    log_prob = log_prob_new
                else
                    chain(:,i) = x_old
                end if
                prob_chain(:,i) = [log_prob, log_prob_new]
            end do
            chain = prob_chain
        end function mcmc

end module transits

program transit
    use numerical
    use transits
    implicit none
    integer, parameter :: arr_size = 2157
    real, dimension(arr_size) :: time, flux, flux_err
    integer :: i
    !Read in both text files and save to arrays
    open(1, file='time.txt', status='old')
    open(2, file='flux.txt', status='old')
    open(3, file='flux_err.txt', status='old')
    do i=1, arr_size
        read(1,*) time(i)
        read(2,*) flux(i)
        read(3,*) flux_err(i)
    end do
    close(1)
    close(2)
    close(3)
    !run the mcmc function 
    print *, denominator(0.5)
    print *, mcmc([0.07, 0.2], 100, 0.1, time, flux, flux_err)


    
end program