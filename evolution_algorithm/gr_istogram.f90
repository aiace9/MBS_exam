module gr_istogram_module

	implicit none
	
	private :: dp, hist, a, b, step, nbin


	public :: push_gr, init_gr, save_data_gr

	!privat declaretions
	integer, parameter :: dp = selected_real_kind(12)
	real(kind = dp), parameter:: pi = 3.1415
	
	integer, save, dimension(:), allocatable :: hist
	real(kind = dp), save :: a, b, step
	integer,save :: nbin


contains
	subroutine init_gr(side,delta)
		! this routine initialize the g(r) histogram
		real(kind = dp), intent(in) :: side,delta
		integer :: err
		! the histogram goes from 0 to half
		! of the box size
	 	a = 0
		b = side/2.0
		step = delta
		nbin= int((b-a)/step)
		allocate(hist(nbin), stat=err)
		if (err /= 0) print *, "hist: Allocation request denied"
		hist = 0		
	end subroutine init_gr

	subroutine push_gr(x, side)
		! this routine add a point in to the histogram
		! be careful x should be the distance rij between
		! the particle
		real(kind = dp), intent(in) :: x
		integer :: ibin
		logical :: debug = .true.

		ibin = floor ( x / step) + 1
		! I must do a check but I already took half of the box
		! so a point out of it is forgotten
		if (ibin <= nbin )then
      		hist(ibin) = hist(ibin) + 1
	    end if

	    if ( debug ) print*, 'D: hist(ibin)', hist(ibin)
	end subroutine push_gr

	subroutine save_data_gr(rho)
		real(kind = dp), intent(in):: rho
		integer :: myunit, ibin, points
		integer :: ios, err
		logical :: debug = .true.
		real(kind = dp) :: r,vb,nid

		open(newunit=myunit, file='ist.dat', iostat=ios, status="unknown", action="write")
		if ( ios /= 0 ) stop "Error opening file ist.dat"

		points =  sum(hist, dim=1)! compute the total points
		if ( debug ) print*, 'D: number of points', points
		if ( debug ) print*, 'D: number of bins', nbin
		if ( debug ) print*, 'D: number of hist', hist

		do ibin= 1 ,nbin
			if ( debug ) print*, 'D: saving bin', ibin, hist(ibin)
			r = (ibin-0.5)*step - a
			vb = (4.0/3.0)*pi*((nbin)**3-(nbin-1)**3)*step**3
			nid = vb * rho
    		write(unit=myunit,fmt=*) r,hist(ibin)/float(points)/step
  		end do

		if ( debug ) print*, 'D: all bin saved'
	
		if (allocated(hist)) deallocate(hist, stat=err)
		if (err /= 0) print *, "hist: Deallocation request denied"

		if ( debug ) print*, 'D: deallocation ok'

		close(unit=myunit, iostat=ios)
		if ( ios /= 0 ) stop "Error closing file ist.dat"
	end subroutine save_data_gr

end module gr_istogram_module
