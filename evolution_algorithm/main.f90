program main
  
  use interaction
  use jmol_module
  use anderson_module
  use Box
  use gr_istogram_module

  implicit none
  integer,parameter::kr=selected_real_kind(12)
  
  !---variabili per il caricamento file--!
  character(100) :: file
  integer :: nbody
  real(kind=kr) :: massa, side
  !-variabili per l'evoluzione temporale-!
  integer :: nstep
  real(kind=kr) :: dt,vmax
  real(kind=kr) :: mvel,mepot,mekin
  real(kind=kr),dimension(:),allocatable :: velcm
  real(kind=kr),dimension(:,:),allocatable :: pos, pos_old, pos_new
  real(kind=kr),dimension(:,:),allocatable :: ekin,vel,f
  real(kind=kr),dimension(3)::D
  real(kind=kr) :: dummy_real
  !----------variabili per lo snapshot---!
  character(len=3), dimension(:),allocatable :: atom_type
  character(len=100) :: comment,dummy
  character(len=100) :: snap_shot_name
  !----------variabili per anderson------!
  logical :: anderson_evol
  real(kind=kr) :: T_bath, beta, rnd,T_ist
  real(kind=kr) :: anderson_freq
  integer :: anderson_step
  !----------flag per la dinamica--------!
  logical :: dyn
  !----------flag di errore--------------!
  integer :: err, ios, istat
  logical :: debug = .true.
  !----------pressure--------------------!
  real :: Pist
  !----------g(r)------------------------!
  real(kind=kr) :: rij
  real(kind=kr) :: delta
  real(kind=kr), dimension(3) :: posij

  !-------indici-------!
  integer :: it, i, j
  
  ! input utente
  
  write(unit=*,fmt="(a)",advance="no")"dinamica?"
  read*, dyn
  
  if (dyn) then
    write(unit=*,fmt="(a)",advance="no")"vmax?"
    read*, vmax 
  end if
  
  write(unit=*,fmt="(a)",advance="no")"n corpi?"
  read*, nbody
  
  ! set di allocazioni
  allocate(vel(3,nbody))
  allocate(ekin(3,nbody))
  allocate(pos(3,nbody))
  allocate(pos_old(3,nbody))
  allocate(pos_new(3,nbody))
  allocate(f(3,nbody))
  allocate(velcm(3))

  ! mie allocazioni
  allocate(atom_type(nbody), stat=err)
  if (err /= 0) print *, "atom_type: Allocation request denied"
  
  pos = 0
  pos_old = 0
  pos_new = 0
  atom_type = ''

  !ulteriori input
  write(unit=*,fmt="(a)",advance="no")"dt: "
  read*,dt
  write(unit=*,fmt="(a)",advance="no")"n.step: "
  read*,nstep
  print*,"body mass: "
  read*,massa
  print*, "side of the box, only cubic boxes are supported:"
  read*, side

  !-------anderson option---------!
  print*, "do you want to implement Anderson? t/f"
  read*, anderson_evol
  if (anderson_evol .eqv. .true.) then
      print*, "Bath temperature?"
      read*, T_bath
      beta = 1/T_bath
      print*, "interaction frequency:"
      read*, anderson_freq
  endif

  !--------input per lo snapshot-----------!
  print*, 'nome del file dello snapshot:'
  read*, snap_shot_name
  print*, 'commento (opzionale)'
  read*, comment
  print*, 'Tipo di atomo che stiamo trattando, funzione non supportata'
  read*, dummy 
  atom_type = 'Ar'

  !------caricamento delle coordinate------!
  print*, "nome del file con coordinate"
  read*, file
  
  open(unit=1,file=file)
  do
    read(unit=1, fmt=*, iostat=istat) pos
    if ( istat /= 0 ) exit
  end do
  close(unit=1)
  
  if ( debug ) print*, 'D - coordinate caricate con successo'


  !--------caricamento o creazione delle velocità---!
  if (dyn) then
    print*, "nome del file con velocità, (0 se il file non esiste)"
      read*, file
      if ( file == '0' ) then
        ! inizializzazine di nuove velocita
        vel=0.
        call random_number(vel)
        ! riscalazione del dominio
        vel=2*vmax*(vel-0.5)
        !rimozione della velocità del centro di massa
        do i=1,3
          velcm(i) = sum(vel(i,:))/nbody
          vel(i,:) = vel(i,:) - velcm(i)
        end do
      if ( debug ) print*, 'D - generazione velocita, successo'
      else
        !-----caricamento delle velocita-----!
        open(unit=1,file=file)
          do
            read(unit=1, fmt=*, iostat=istat) vel
            if ( istat /= 0 ) exit 
          end do
        close(unit=1)
      if ( debug ) print*, 'D - caricamento velocita, successo'
      end if
    
  end if

  !---------apertura dei file per la scrittura----!
  open(unit=1, file='pos_vel.dat', iostat=ios, status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file pos_vel.dat"

  open(unit=2, file='kin_pot_tot.dat', iostat=ios, status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file kin_pot_tot.dat"

  open(unit=3, file='T_P.dat', iostat=ios, status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file T.dat"

  !-----------algoritmo principale-----------!
  !generazione del porimo step senza verlet
  call interazione(pos,nbody,f,mepot, side)
  pos_old = pos
  do i = 1,nbody
    pos(:,i) = pos(:,i) + vel(:,i) * dt + 0.5* f(:,i)/massa * dt**2
    vel(:,i) = vel(:,i) + dt * f(:,i)/massa
  end do
  
  
  if ( debug ) print*, 'D - prima chiamata routine interazione, successo'
  do it = 1,nstep
    if (dyn) then
      !----calcolo l'interazione fra le particelle----!
      call interazione(pos,nbody,f,mepot, side)
      
      !-----aggiorno le posizioni ------!
      !-----pos_new => pos t+1
      do i = 1,nbody
        pos_new(:,i) = 2 * pos(:,i) - pos_old(:,i)+ dt**2 * f(:,i) / massa
      end do
      
      !---calcolo velocità ed energia----!
      !---al tempo t 
      T_ist = 0    
      do i=1,nbody
        call PBC(pos_new(:,i), pos_old(:,i), dummy_real, D, side, .false.)
        vel(:,i) = (D )/ (2 * dt)
        ekin(:,i) = 0.5 * massa * (vel(:,i))**2
        T_ist = T_ist + dot_product(vel(:,i),vel(:,i))
      end do
      T_ist = (T_ist/(3.0*size(vel, dim=2)))* massa

      
      mekin = sum(ekin)
      
      !----salvataggio dati-------!
      !----al tempo t
      if (mod(it,50) == 0) then
        !write(unit=1,fmt=*)it,it*dt,pos,vel
        write(unit=2,fmt=*)it,mekin,mepot,mekin+mepot
        write(unit=3,fmt=*)it, it*dt, T_ist
      endif
      
      pos_old = pos       
      pos = pos_new

      !-----riposiziono le particelle all'interno della scatola----!
      call scatola(pos,side)

      !-----percentuale----!
      if( mod(it,nstep/10) == 0) print*, floor(it/(nstep*1.0)*100)

    else
      !---- not pure dynamic evolution, thermalization of the sample-----!
      !-----primo step pos vel------!
      do i = 1,nbody
        pos(:,i) = pos(:,i) + vel(:,i) * dt + 0.5* f(:,i)/massa * dt**2
        vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/massa
      end do

      !-----riposiziono le particelle all'interno della scatola----!
      call scatola(pos,side)
      
      !----calcolo l'interazione fra le particelle----!
      call interazione(pos,nbody,f,mepot, side)
      
      !----implementazione di anderson-----!
      !----beta è 1/T attenzione alle formule ------!
      call anderson_integration(vel,f,massa,beta,anderson_freq,dt,T_ist)
      
      do i=1,nbody
        ekin(:,i) = 0.5 * massa * (vel(:,i))**2
      end do
      mekin = sum(ekin)
      
      !----salvataggio dati-------!
      if (mod(it,50) == 0) then
        write(unit=1,fmt=*)it,it*dt,pos,vel
        write(unit=2,fmt=*)it,mekin,mepot,mekin+mepot
        write(unit=3,fmt=*)it, it*dt, T_ist
      endif
      !-----percentuale----!
      if( mod(it,nstep/10) == 0) print*, floor(it/(nstep*1.0)*100)

    end if
  !-----sanapshot per il video-----!
  if (nstep < 1000) call snapshot(snap_shot_name, atom_type, pos, comment, .true.)
  end do

  ! evaluation of g(r)
  if (.true.) then
    delta = 0.01
    call init_gr(side,delta)
    do i=1,nbody
      do j=1,nbody
        if( i==j ) cycle
        call PBC(pos(:,i), pos(:,j), rij, posij, side, .false.)
        !print *, rij
        call push_gr(rij)
      end do
    end do
    call save_data_gr(nbody/side**3, nbody)
  end if

  !snapshot situazione finale
  !call snapshot(snap_shot_name, atom_type, pos, comment, .false.)

  !-----chiusura di tutti i file aperti------!
  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"

  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"

  close(unit=3, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 3"

  !-------salviamo le posizioni finali per un eventuale nuovo start------!
  open(unit=4, file='sample2_pos.dat', iostat=ios,  status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file sample2_pos.dat"
    
  write(unit=4,fmt=*) pos

  close(unit=4, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 4"

  !-------salviamo le velocità finali per un eventuale nuovo start-------!
  open(unit=4, file='sample2_vel.dat', iostat=ios,  status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file sample2_vel.dat"
  
  write(unit=4,fmt=*) vel

  close(unit=4, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 4"
  
  !mie deallocazioni
  if (allocated(atom_type)) deallocate(atom_type, stat=err)
  if (err /= 0) print *, "atom_type: Deallocation request denied"

end program main