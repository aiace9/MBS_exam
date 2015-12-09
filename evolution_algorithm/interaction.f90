module interaction

  use Box

  implicit none

  private :: kr
  public :: interazione, interazione1
  
  !integer,parameter::kr=selected_real_kind(12)
  integer,parameter::kr=selected_real_kind(12)
  real(kind=kr), parameter ::eps=1.,sigma=1.
  
  contains

  subroutine interazione(pos,nbody,f,upot, side, Pist)
    ! questa sabroutine simula una scatola cubica di lato side
    ! ma permette alle particelle di allontanarsia all' infinito
    ! dall'origine in cui si trovano
    ! adesso aggiungo anche la pressione come parametro opzionale 
    ! la pressione attualmente non è supportata, NON UTILIZZARE
    ! dovrebbe mancare un termine che va sommato e va divisa per dV
    real(kind=kr), intent(in), dimension(:,:) :: pos
    real(kind=kr), intent(in) :: side
    integer, intent(in) :: nbody
    real(kind=kr), intent(out) :: upot
    real(kind=kr), intent(out), dimension(:,:) :: f
    real(kind=kr), intent(out), optional :: Pist
    real(kind=kr), dimension(size(pos,1)) :: posij
    real(kind=kr) :: rij
    integer ::i,j
    upot = 0
    f = 0
    Pist = 0
    do i=1,nbody
      do j=1,nbody
        if( i==j ) cycle
        call PBC(pos(:,i), pos(:,j), rij, posij, side, .false.)
        upot = upot + 4*eps*(rij**(-12)-rij**(-6)) 
        f(:,i) = f(:,i) + 24*eps*(2.*rij**(-14)-rij**(-8))*posij
        if ( j > i ) Pist = Pist + dot_product(24*eps*(2.*rij**(-14)-rij**(-8))*posij,posij)
      end do
    end do
    upot = upot/2
  end subroutine interazione

  subroutine interazione1(pos,nbody,f,upot, side)
    !questa routine simula uno spazio infinito
    real(kind=kr), intent(in), dimension(:,:) :: pos
    real(kind=kr), intent(in)::side
    integer, intent(in) :: nbody
    real(kind=kr), intent(out) :: upot
    real(kind=kr), intent(out), dimension(:,:) :: f
    real(kind=kr), dimension(size(pos,1)) :: posij
    real(kind=kr) :: rij
    integer ::i,j
    upot = 0
    f = 0
    do i=1,nbody
      do j=1,nbody
        if( i==j ) cycle
        posij = pos(:,i)-pos(:,j)
        rij=sqrt( dot_product(posij,posij) )
        upot = upot + 4*eps*(rij**(-12)-rij**(-6)) 
        f(:,i) = f(:,i) + 24*eps*(2.*rij**(-14)-rij**(-8))*posij
      end do
    end do
    upot = upot/2
  end subroutine interazione1

end module interaction