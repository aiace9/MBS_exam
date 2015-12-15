program corpi3d
  

  use Box
  
  implicit none
  integer,parameter::kr=selected_real_kind(12)
  real(kind=kr) :: dist, side
  real(kind=kr),dimension(3) :: posa ,posb, rij
  real(kind=kr),dimension(3,6) :: pos
  
  posa =(/ 1.5, 3.6, 7.5 /)
  posb =(/ 1.5, 3.6, 0.0 /)
  side = 10.


  print*, "TEST PBC"
  call PBC(posa, posb, dist, rij, side, .true.)

  side = 5.
  pos(:,1) = posa
  pos(:,2) = posb
  pos(:,3) = (/ -1.5, 9.6, -2.7 /)
  pos(:,4) = (/ -1.5, -9.6, -2.7 /)
  pos(:,5) = (/ 42.0, 3.0, 9.0 /)
  pos(:,6) = (/ -23.0, 9.6, -7.5 /)

  print*, "TEST SCATOLA"
  call scatola(pos, side)

end program corpi3d

 