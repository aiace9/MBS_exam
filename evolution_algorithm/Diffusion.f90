module Diffusion

	implicit none
	
	private :: n_call, n_time, vacf, it0, t0, vx0, time0


	public :: push_D, init_D, save_D

	!privat declaretions
	integer, parameter :: dp = selected_real_kind(12)
	real(kind = dp), parameter:: pi = 3.1415
	
	integer, save, dimension(:), allocatable :: n_time, time0
	real, save, dimension(:), allocatable :: vacf
	real, save, dimension(:,:,:), allocatable :: vx0
	integer,save :: n_call, it0, t0, t0_max


contains
	subroutine init_D(nstep, nbody, b)
		! this routine initialize D
		! nstep: number of the time step
		! b: distance between two window stars
		integer, intent(in) :: nstep,nbody,b
		integer :: err
		n_call = 0
		it0 = b
		t0_max = int(nstep/it0)
		allocate(n_time(nstep), stat=err)
		if (err /= 0) print *, "n_time: Allocation request denied"
		allocate(vacf(nstep), stat=err)
		if (err /= 0) print *, "vacf: Allocation request denied"
		allocate(time0(t0_max), stat=err)
		if (err /= 0) print *, "time0: Allocation request denied"
		allocate(vx0(3,nbody,t0_max), stat=err)
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
		integer :: delt, tt0
		integer :: t,i
		logical :: debug = .false.
		n_call = n_call + 1
		if (debug) print*, 'n_call: ', n_call
		if (mod(n_call - 1,it0) == 0) then
			t0 = t0 + 1
			tt0 = mod(t0-1, t0_max) + 1
			time0(tt0) = n_call
			vx0 (:,:,tt0) = vel
			if (debug) then
				print*,'D - t0: ', t0
				print*,'D - time0:', time0
				print*,'D - velocity:', vx0(:,1,tt0)
			endif
		endif
		do t = 1, min(t0,t0_max), 1
			! n_call = step in our case.
			delt = n_call - time0(t) +1
			if (debug) print*,'D - delta', t ,'=',  delt
			if (delt <= step) then
				n_time(delt) = n_time(delt) + 1
				do i = 1, nbody, 1
					vacf(delt) = vacf(delt) + dot_product(vel(:,i),vx0(:,i,t))
					if (debug) then
						print*,'D - vacf(delta):', vacf(delt)
						print*,'D - dot:',dot_product(vel(:,i),vx0(:,i,t))
					end if
				end do
			end if
		end do
	end subroutine push_D

	subroutine save_D(nstep, nbody, dt)
		! this routine normalize the data and save them.
		integer, intent(in):: nstep,nbody
		real(kind = dp), intent(in) :: dt
		real :: sum
		logical :: debug = .false.
		integer ::i 
		sum = 0
		if (debug) then
			do i = 1, t0_max, 1
				print*, 'D - n_time', n_time(it0*i)
			end do
		end if
		do i = 1, nstep, 1
			if (debug) then
				print*, 'D - sum:',sum
				print*, 'D - vacf(i):',vacf(i)
			end if 
			sum = sum + vacf(i)/ (nbody * n_time(i))
		end do
		print *, 'the Diffusion coefficient is:', sum * dt
	end subroutine save_D

end module Diffusion
