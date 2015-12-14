module Diffusion

	implicit none
	
	private :: n_call, ntime, vacf, it0, t0, vx0, time0


	public :: push_D, init_D, save_D

	!privat declaretions
	integer, parameter :: dp = selected_real_kind(12)
	real(kind = dp), parameter:: pi = 3.1415
	
	integer, save, dimension(:), allocatable :: n_time, vacf, time0
	real, save, dimension(:,:,:), allocatable :: vx0
	integer,save :: n_call, it0, t0


contains
	subroutine init_D(nstep, nbody, b)
		! this routine initialize D
		! nstep: number of the time step
		! b: distance between two window stars
		integer, intent(in) :: nstep
		integer :: err
		n_call = 0
		it0 = b
		allocate(n_time(nstep), stat=err)
		if (err /= 0) print *, "n_time: Allocation request denied"
		allocate(vacf(nstep), stat=err)
		if (err /= 0) print *, "vacf: Allocation request denied"
		allocate(time0(int(nstep/it0)), stat=err)
		if (err /= 0) print *, "time0: Allocation request denied"
		allocate(vx0(3,nbody,int(nstep/it0)), stat=err)
		if (err /= 0) print *, "time0: Allocation request denied"
		n_time = 0	
		vacf = 0
		time0 = 0
		t0 = 0
	end subroutine init_D

	subroutine push_D(vel, step, nstep, nbody)
		! 
		real(kind = dp), intent(in), dimension(:,:) :: vel
		integer, intent(in) :: step, nstep, nbody
		integer :: delt
		integer :: t
		logical :: debug = .false.
		n_call = n_call + 1
		if (mod(n_call,it0) == 0) then
			t0 = t0 + 1
			tt0 = mod(t0-1) + 1
			time0(tt0) = n_call
			vx0 (:,:,tt0) = vel
		endif
		do t = 1, min(t0,nstep), 1
			! n_call = step in our case.
			delt = n_call - time0(t) +1
			if (delt < step) then
				n_time(delt) = ntime(delt) + 1
				do i = 1, nbody, 1
					! check dot product, may be not correct
					vacf(delt) = vacf(delt) + dot_product(vel(:,i),vx0(:,i,t))
				end do
		end do
	end subroutine push_D

	subroutine save_D(nstep, nbody, dt)
		! this routine normalize the data and save them.
		integer, intent(in):: nstep,nbody
		real(kind = dp), intent(in) :: dt
		real :: sum
		logical :: debug = .false.
		sum = 0
		do i = 1, nstep, 1
			sum = sum + vacf(i)/ (nbody * n_time(i))
		end do
		print *, 'the Diffusion coefficient is:', sum * dt
	end subroutine save_D

end module Diffusion