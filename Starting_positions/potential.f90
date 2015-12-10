module interaction

  use Box

  implicit none

  private :: kr
  public :: potential, gradient
  
  integer,parameter::kr=selected_real_kind(12)
  real(kind=kr), parameter ::eps=1.,sigma=1.
  
  contains

  subroutine potential(pos,nbody,upot, side)
    ! questa sabroutine simula una scatola cubica di lato side
    ! ma permette alle particelle di allontanarsia all' infinito
    ! dall'origine in cui si trovano
    ! la pressione attualmente non è supportata, NON UTILIZZARE
    ! dovrebbe mancare un termine che va sommato e va divisa per dV
    real(kind=kr), intent(in), dimension(:,:) :: pos
    real(kind=kr), intent(in) :: side
    integer, intent(in) :: nbody
    real(kind=kr), intent(out) :: upot
    real(kind=kr), dimension(size(pos,1)) :: posij
    real(kind=kr) :: rij
    integer ::i,j
    upot = 0
    ! this should be improved a lot 
    do i=1,nbody
      do j=1,nbody
        if( i==j ) cycle
        call PBC(pos(:,i), pos(:,j), rij, posij, side, .false.)
        ! questa è la linea che definisce il potenziale
        upot = upot + 4*eps*(rij**(-12)-rij**(-6)) 
      end do
    end do
    upot = upot/2
  end subroutine potential

  subroutine gradient(pos,nbody,side, direction)
    real(kr), intent(in), dimension(:,:) :: pos
    integer, intent(in) :: nbody
    real(kind=kr), intent(in) :: side
    real(kr), intent(out), dimension(:,:) :: direction
    real(kind=kr), dimension(size(pos,1)) :: posij
    real(kind=kr) :: rij
    integer ::i,j
    logical :: debug = .false.
    direction = 0
    do i=1,nbody
      do j=1,nbody
        if( i==j ) cycle
        call PBC(pos(:,i), pos(:,j), rij, posij, side, .false.)
        direction(:,i) = direction(:,i) + 24*eps*(2.*rij**(-14)-rij**(-8))*posij
        if ((direction(1,i)>1 .or. direction(2,i)>1 .or. direction(3,i)>1 ).and. debug) then
          print *, i,j, direction(:,i)
          print *, rij
          print *, posij
        endif 
      end do
    end do
  end subroutine gradient

end module interaction