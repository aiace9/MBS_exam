module Box

	implicit none
	private :: kr
	public :: PBC
	integer,parameter::kr=selected_real_kind(12)	

  contains

	subroutine PBC (posa, posb, dist, delta, side, debug)
		! questa routine applica le pbc senza spostare le particelle
		! funziona solo per scatole cubiche
		! posa, posb, posizione delle due particelle
		! dist,lunghezza vettore fra le due particelle
		! delta, componeneti x, y, z del vettore distanza
		real(kind = kr),intent(in), dimension(:) :: posa, posb
		real(kind = kr), intent(out) :: dist
		real(kind = kr), intent(in) :: side
		real(kind = kr), dimension(3), intent(out):: delta
		logical, intent(in) :: debug
		real(kind = kr), dimension(3) :: supp
		!verificare se modulo è corretto oppure no
		supp = mod(posa(:)-posb(:),side)
		delta = supp(:) - side * int(2 * supp(:) / side)

		dist = sqrt(dot_product(delta,delta))	

		if ( debug ) then
			print*, "---------"
			print*, 'posizione prima particella:', posa
			print*, 'poszione seconda particella:', posb
			print*, 'delta:', delta
			print*, 'distanza fra le particelle:', dist
			print*, 'supp:', supp
			print*, "---------"
		end if			
	end subroutine 

	subroutine scatola(pos, side)
		!Questa routine riposiziona le particelle nella scatola 
		real(kind = kr), dimension(:,:), intent(inout) :: pos
		real(kind = kr), intent(in) :: side
		logical:: debug
		real(kind = kr), dimension(3) :: supp
		integer :: i
		debug = .false.

		do i = 1, size(pos, dim=2), 1
			if ( debug ) then
				print*, "------ particella ", i , " di ", size(pos, dim=2)
				print*, "posizione iniziale:", pos(:,i)
			end if
			pos(:,i) = modulo(pos(:,i),side)
			if ( debug ) print*, "posizione finale:", pos(:,i)
		end do

	end subroutine scatola

end module Box
