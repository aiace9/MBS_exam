module anderson_module

	implicit none
	
	private :: kr
	public :: anderson_integration
	integer, parameter :: kr = selected_real_kind(12)


contains

	subroutine anderson_integration(vel,f,mass,beta,nu,dt,tempa)
		real(kr), dimension(:,:), intent(inout) :: vel
		real(kr), intent(out) :: tempa
		real(kr), dimension(:,:), intent(in) :: f
		real(kr), intent(in) :: beta, mass, nu, dt
		! variabili interne
		real(kr) :: sigma
		real(kr) :: rnd
		integer :: i,j

		logical :: debug = .false.
		integer :: coll = 0
		integer, save :: coll_tot, chiamate

		! tempa è in kelvin * kb
		tempa = 0
		coll = 0
		do i = 1, size(vel, dim=2), 1
	          vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/mass
	          tempa = tempa + dot_product(vel(:,i),vel(:,i))
		end do
		
		tempa = (tempa/(3.0*size(vel, dim=2)))* mass
		sigma = sqrt(1.0/(beta*mass))

		do i = 1, size(vel, dim=2), 1
			call random_number(rnd)
			if ( rnd .lt. nu * dt ) then
				if ( debug ) coll = coll +1
				!if ( debug ) print*, 'D - collisione avvenuta'
				do j = 1, size(vel, dim=1), 1
					call gasdev(vel(j,i))
					vel(j,i) = sigma * vel(j,i)
				end do
			end if

		end do
		if (debug) then
			print*, 'D - debug mode on per anderson_integration'
			print*, 'D - collissioni/particelle', coll/(size(vel, dim=2)*1.0)*100.0
			coll_tot = coll_tot + coll
			chiamate = chiamate + 1
			print*, 'D - collisioni totali', coll_tot
			print*, 'D - chiamate alla funzione', chiamate
			print*, 'D - valore medio collisioni/particelle', coll_tot/(chiamate*size(vel, dim=2)*1.0)*100.0
			print*, 'D ------------------------------------'

		endif
	end subroutine anderson_integration

    subroutine gasdev(rnd)
    ! this algorithm generate a gaussian distribution
    ! with:
    ! sigma = 1
    ! mu = 0
    real(kind=kr), INTENT(OUT) :: rnd
    real(kind=kr) :: r2,x,y
    real(kind=kr), SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.
    if (gaus_stored) then
       rnd=g
       gaus_stored=.false.
    else
       do
          call random_number(x)
          call random_number(y)
          x=(2.*x-1.)
          y=(2.*y-1.) 
          r2=x**2+y**2
          if (r2 > 0. .and. r2 < 1.) exit
       end do
       r2=sqrt(-2.*log(r2)/r2)
       rnd=x*r2
       g=y*r2
       gaus_stored=.true.
    end if
  end subroutine gasdev

end module anderson_module
