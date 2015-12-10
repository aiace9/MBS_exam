module jmol_module

	implicit none
	
	private :: kr
	public :: snapshot

	integer, parameter :: kr = selected_real_kind(12)

contains

	subroutine snapshot(name, atom_type, pos, comment, video)
		character(len=*), intent(in) :: name
		character(len=*), intent(in), dimension(:) :: atom_type
		real(kind = kr), intent(in), dimension (:,:) :: pos 
		character(len=*), intent(in) :: comment
		logical, intent(in) :: video
		logical :: debug = .false.

		integer :: i,ios,myunit
		
		if ( debug )  then
			print*, 'debug mode on'
			print*, 'file name', name
			!print*, atom_type
		end if

		if ( .not.video ) then
			open(newunit=myunit, file=name, iostat=ios, status="unknown", action="write")
			if ( ios /= 0 ) stop "Error opening file for snapshot"
		else
			open(newunit=myunit, file=name, iostat=ios, status="unknown", action="write", position= "append")
			if ( ios /= 0 ) stop "Error opening file for video"
		end if

		if ( debug ) then
			print*, size(pos, dim=2)
			print*, comment
		endif

		write(unit=myunit, fmt=*, iostat=ios) size(pos, dim=2)
		if ( ios /= 0 ) stop "Write error in file unit myunit"

		write(unit=myunit, fmt=*, iostat=ios) comment
		if ( ios /= 0 ) stop "Write error in file unit myunit"
		
		
		do i = 1, size(pos, dim=2), 1
			
			if ( debug ) print*, atom_type(i), pos(:,i)

			write(unit=myunit, fmt=*, iostat=ios) atom_type(i), pos(:,i)
			if ( ios /= 0 ) stop "Write error in file unit myunit"
			
		end do

		close(unit=myunit, iostat=ios)
		if ( ios /= 0 ) stop "Error closing file unit myunit"

	end subroutine snapshot

end module jmol_module
