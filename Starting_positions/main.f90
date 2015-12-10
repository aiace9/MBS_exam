program starting_positions
! ________________________________
!< Debug is human, de-fix divine. >
! --------------------------------
!        \   ^__^
!         \  (oo)\_______
!            (__)\       )\/\
!                ||----w |
!                ||     ||
	

	use jmol_module
	use	interaction
	
	implicit none
	integer,parameter :: dp = selected_real_kind(12)
	character(len=9),parameter :: filename = 'start.txt'
	! lattice variables
	real(dp), dimension(3) :: box_size !lx,ly,lz
	integer, dimension(3) :: quanta_box !nx,ny,nz
	real(dp) :: a !lattice_constant
	character(len=3) :: lattice !kind of lattice, only fcc supported

	! particle parameter
	integer :: n_molecule, max_molecule
	character(len=1) :: extra_molecule
	real(dp), dimension(:,:),allocatable :: pos

	! iteration index
	integer :: i,j,k,q
	integer :: m_index = 0

	!Snapshot variable
	character(len=3), dimension(:),allocatable :: atom_type
  	character(len=100) :: comment,dummy
  	character(len=10) :: snap_shot_name

  	! sample preparation variable
  	real(dp) :: e_tot_i0, e_tot_i1
  	real(dp), dimension(:,:), allocatable :: direction
  	integer :: recursion = 10
  	real(dp):: epsilon = 1.0d-8
  	real(dp) :: lambda = 0.01

  	! random variable
  	integer :: seed1, sizer
  	integer, allocatable, dimension(:) :: seed
  	double precision, dimension(3) :: rnd(3)	

	!error flag
	integer :: ios, err
	logical :: debug = .true.

	print *, "this algorithm will prepare the sample"
	print*, 'insert the box size (x,y,z)'
	read*, box_size
	print*, 'insert the lattice constant'
	read*, a
	print*, 'insert the lattice type (not implemented, FCC only)'
	!read*, lattice
	quanta_box = floor(box_size/a)
	
	!this will depend on lattice type
	print*, 'you may put up to: ',quanta_box(1)*quanta_box(2)*quanta_box(3)*4, 'molecule in the box'
	do
		print*, 'how many molecule do you want in the box?'
		read*, n_molecule
		if ( n_molecule <= quanta_box(1)*quanta_box(2)*quanta_box(3)*4 ) exit
	end do

	print*, "do you want to add an extra molecule in a random position? (y,n)"
	read*, extra_molecule
	
	if (extra_molecule == "y") then
		allocate(pos(3,n_molecule + 1), stat=err)
		if (err /= 0) print *, "pos: Allocation request denied"
		
		! random settings
		print*," seed (one integer)?"
  		read*, seed1
  		call random_seed(sizer)
  		allocate(seed(sizer), stat=err)
  		if (err /= 0) print *, "seed: Allocation request denied"
		do q = 1, sizer, 1
      	seed(q)= seed1 + 37 * (q-1)
  		end do
  		call random_seed(put=seed)
	else
		allocate(pos(3,n_molecule), stat=err)
		if (err /= 0) print *, "pos: Allocation request denied"
	end if 
	
	if ( debug ) print*, 'D: array allocation', size(pos, dim=1), size(pos, dim=2)

	allocate(atom_type(n_molecule), stat=err)
	if (err /= 0) print *, "atom_type: Allocation request denied"
	
	!sequential filling for an FCC lattice
	do k = 0, quanta_box(3)-1, 1
		do j = 0, quanta_box(2)-1, 1
			do i = 0, quanta_box(1)-1, 1
				! this is then first atom in the (i,j,k)
				! position. This position is always available
				m_index= m_index + 1
				if (m_index> n_molecule) exit
				pos(1,m_index)= i*a
				pos(2,m_index)= j*a
				pos(3,m_index)= k*a
				! form now on the position are not always available
				! because some of them could be out of the box
				! so there is an extra check
				if (j*a+a/2 < box_size(2) .and. k*a+a/2 < box_size(3) ) then
					m_index= m_index + 1
					if (m_index> n_molecule) exit
					pos(1,m_index)= i*a
					pos(2,m_index)= j*a+a/2
					pos(3,m_index)= k*a+a/2
				end if

				if (i*a+a/2 < box_size(1) .and. j*a+a/2 < box_size(2) ) then
					m_index= m_index + 1
					if (m_index> n_molecule) exit
					pos(1,m_index)= i*a+a/2
					pos(2,m_index)= j*a+a/2
					pos(3,m_index)= k*a
				end if 

				if (i*a+a/2 < box_size(1) .and. k*a+a/2 < box_size(3) ) then
					m_index= m_index + 1
					if (m_index> n_molecule) exit
					pos(1,m_index)= i*a+a/2
					pos(2,m_index)= j*a
					pos(3,m_index)= k*a+a/2
				end if 
			end do
			if (m_index> n_molecule) exit
		end do
		if (m_index> n_molecule) exit
	end do

	! fix the number  of molecules
	m_index = m_index -1
	
	if ( debug ) print*, 'D: real number of molecules', m_index

	open(unit=1, file=filename, iostat=ios, status="unknown", action="write")
	if ( ios /= 0 ) stop "Error opening file"

	write(unit=1, fmt=*, iostat=ios) pos(:,1:m_index)
	if ( ios /= 0 ) stop "Write error in file unit 1"
	
	close(unit=1, iostat=ios)
	if ( ios /= 0 ) stop "Error closing file unit 1"

  	! ATTENTION
  	! from now on the code works only for cubic box! 
  	! this is due to the potential routine

  	! potential evaluation
  	call potential(pos(:,1:m_index),m_index,e_tot_i1,box_size(1))
  	print* , "first potential evaluation: ", e_tot_i1


  	! addition of an extra molecule if required
  	if (extra_molecule == "y") then
		m_index = m_index + 1
		call random_number(rnd)
		pos(1,m_index)= rnd(1) * box_size(1)
		pos(2,m_index)= rnd(2) * box_size(2)
		pos(3,m_index)= rnd(3) * box_size(3)
		print* , "extra molecule in: ", pos(:,m_index)
	end if 

	!input for the snapshot
  	print*, 'name of the snapshot:'
  	read*, snap_shot_name
  	print*, 'comment(optional):'
  	read*, comment
  	print*, 'Kind of atom, not supported'
  	read*, dummy 
  	atom_type = 'Ar'
	call snapshot(snap_shot_name, atom_type, pos(:,1:m_index), comment, .false.)

  	call potential(pos(:,1:m_index),m_index,e_tot_i0,box_size(1))
  	print*, "second potential evaluation:",e_tot_i0
  	print*, "delta: (is expected > 0 for high density) ", e_tot_i0 - e_tot_i1

  	! start of the steepest descent
  	print*, "starting steepest descent algorithm..."

  	allocate(direction(3,m_index), stat=err)
  	if (err /= 0) print *, "direction: Allocation request denied"
  	
  	call gradient(pos(:,1:m_index),m_index,box_size(1),direction)
  	pos(:,1:m_index) = pos(:,1:m_index) + lambda * direction
  	!print*, direction
  	open(unit=7, file="energy.dat", iostat=ios, action="write")
  	if ( ios /= 0 ) stop "Error opening file energy.dat"
  	print*, '-----------------------'
  	write(unit=7, fmt=*) 0, e_tot_i0
  	
  	do i= 1, 10000000, 1
  		call potential(pos(:,1:m_index),m_index,e_tot_i1,box_size(1))
  		write(unit=7, fmt=*) i, e_tot_i1
  		if (e_tot_i0 > e_tot_i1 .and. .not. abs(e_tot_i0 - e_tot_i1)>epsilon) then
  			pos(:,1:m_index) = pos(:,1:m_index) + lambda * direction
  			e_tot_i0 = e_tot_i1
  			recursion = 10
  		else
  			recursion = recursion - 1
  			if (recursion == 0) then
  				print*, 'minimun reached'
  				exit
  			end if
  			call gradient(pos(:,1:m_index),m_index,box_size(1),direction)
  			pos(:,1:m_index) = pos(:,1:m_index) + lambda * direction
  			e_tot_i0 = e_tot_i1
  		endif
  	end do

  	close(unit=7, iostat=ios)
  	if ( ios /= 0 ) stop "Error closing file unit 7"
  	
	

	
  	if (allocated(direction)) deallocate(direction, stat=err)
  	if (err /= 0) print *, "direction: Deallocation request denied"

	if (allocated(pos)) deallocate(pos, stat=err)
	if (err /= 0) print *, "pos: Deallocation request denied"
	
	if (allocated(atom_type)) deallocate(atom_type, stat=err)
	if (err /= 0) print *, "atom_type: Deallocation request denied"
	

end program starting_positions
