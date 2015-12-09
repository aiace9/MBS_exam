program starting_positions

	use jmol_module
	
	implicit none
	integer,parameter :: dp = selected_real_kind(12)
	character(len=9),parameter :: filename = 'start.txt'
	
	!lattice variables
	real(dp), dimension(3) :: box_size !lx,ly,lz
	integer, dimension(3) :: quanta_box !nx,ny,nz
	real(dp) :: a !lattice_constant
	character(len=3) :: lattice !kind of lattice, only fcc supported

	!particle parameter
	integer :: n_molecule, max_molecule
	real(dp), dimension(:,:),allocatable :: pos

	!iteration index
	integer :: i,j,k
	integer :: m_index = 0

	!Snapshot variable
	character(len=3), dimension(:),allocatable :: atom_type
  	character(len=100) :: comment,dummy
  	character(len=10) :: snap_shot_name


	!error flag
	integer :: ios, err
	logical :: debug = .true.

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
	
	allocate(pos(3,n_molecule), stat=err)
	if (err /= 0) print *, "pos: Allocation request denied"
	
	if ( debug ) print*, 'D: array allocation', size(pos, dim=1), size(pos, dim=2)

	allocate(atom_type(n_molecule), stat=err)
	if (err /= 0) print *, "atom_type: Allocation request denied"
	
	!sequential filling for an FCC lattice
	do k = 1, quanta_box(3), 1
		do j = 1, quanta_box(2), 1
			do i = 1, quanta_box(1), 1
				
				m_index= m_index + 1
				if (m_index> n_molecule) exit
				pos(1,m_index)= i*a
				pos(2,m_index)= j*a
				pos(3,m_index)= k*a

				if (j*a+a/2 < ly .and. k*a+a/2 < lz ) then
					m_index= m_index + 1
					if (m_index> n_molecule) exit
					pos(1,m_index)= i*a
					pos(2,m_index)= j*a+a/2
					pos(3,m_index)= k*a+a/2
				end if

				if (i*a+a/2 < lx .and. j*a+a/2 < ly ) then
					m_index= m_index + 1
					if (m_index> n_molecule) exit
					pos(1,m_index)= i*a+a/2
					pos(2,m_index)= j*a+a/2
					pos(3,m_index)= k*a
				end if 

				if (i*a+a/2 < lx .and. k*a+a/2 < lz ) then
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
	
	if ( debug ) print*, 'D: molecole posizionate', m_index -1

	open(unit=1, file=filename, iostat=ios, status="unknown", action="write")
	if ( ios /= 0 ) stop "Error opening file filename"

	write(unit=1, fmt=*, iostat=ios) pos
	if ( ios /= 0 ) stop "Write error in file unit 1"
	
	close(unit=1, iostat=ios)
	if ( ios /= 0 ) stop "Error closing file unit 1"

	!input per lo snapshot
  	print*, 'nome del file dello snapshot:'
  	read*, snap_shot_name
  	print*, 'commento (opzionale)'
  	read*, comment
  	print*, 'Tipo di atomo che stiamo trattando, funzione non supportata'
  	read*, dummy 
  	atom_type = 'Ar'

  	!if I will have time I will add an energy check
	
	if ( debug )call snapshot(snap_shot_name, atom_type, pos, comment, .true.)
	

	if (allocated(pos)) deallocate(pos, stat=err)
	if (err /= 0) print *, "pos: Deallocation request denied"
	
	if (allocated(atom_type)) deallocate(atom_type, stat=err)
	if (err /= 0) print *, "atom_type: Deallocation request denied"
	

end program starting_positions
