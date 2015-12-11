program main
  
  use interaction
  use jmol_module
  use anderson_module
  use Box

  implicit none
  integer,parameter::kr=selected_real_kind(12)
  
  integer :: nstep,it,nbody,i,j,istat
  
  real(kind=kr) :: dt,mvel,mepot,massa,alfa=1,vmax,side,vol
  real(kind=kr),dimension(:),allocatable :: mekin,e_tot
  real(kind=kr),dimension(:,:),allocatable :: pos
  real(kind=kr),dimension(:),allocatable :: velcm
  real(kind=kr),dimension(:,:),allocatable :: ekin,vel,f

  ! mie variabli

  !----------variabili per lo snapshot---!
  character(len=3), dimension(:),allocatable :: atom_type
  character(len=100) :: comment,dummy
  character(len=100) :: snap_shot_name
  !----------variabili per anderson------!
  logical :: anderson_evol
  real(kind=kr) :: T_bath,beta,rnd,T_ist
  real(kind=kr) :: anderson_freq
  integer :: anderson_step
  !----------flag per la dinamica--------!
  logical :: dyn
  !----------flag di errore--------------!
  integer :: err, ios
  logical :: debug = .true.
  character(100) :: file
  !----------pressure--------------------!
  real :: Pist
  
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
  allocate(e_tot(nbody)) ! this variable seem unused
  allocate(mekin(nbody))
  allocate(vel(3,nbody))
  allocate(ekin(3,nbody))
  allocate(pos(3,nbody))
  allocate(f(3,nbody))
  allocate(velcm(3))

  ! mie allocazioni
  allocate(atom_type(nbody), stat=err)
  if (err /= 0) print *, "atom_type: Allocation request denied"
  
  pos = 0
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

  !------caricamento delle coordinate iniziali ------!
  print*, "nome del file con coordinate"
  read*, file
  
  open(unit=1,file=file)
  do
    read(unit=1, fmt=*, iostat=istat) pos
    if ( istat /= 0 ) exit
  end do
  close(unit=1)
  
  if ( debug ) print*, 'D - coordinate caricate con successo'


  !--------caricamento della configurazione iniziale delle velocità o creazione della stessa---!
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
  
  call interazione(pos,nbody,f,mepot, side,Pist)
  if ( debug ) print*, 'D - prima chiamata routine interazione, successo'
  do it = 1,nstep
    if (dyn) then
      !-----primo step pos vel------!
      do i = 1,nbody
        pos(:,i) = pos(:,i) + vel(:,i) * dt + 0.5* f(:,i)/massa * dt**2
        vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/massa
      end do

      !-----riposiziono le particelle all'interno della scatola----!
      call scatola(pos,side)
      
      !----calcolo l'interazione fra le particelle----!
      call interazione(pos,nbody,f,mepot, side,Pist)
      
      !----implementazione di anderson-----!
      if(anderson_evol .eqv. .true.) then
        !----beta è 1/T attenzione alle formule ------!
        call anderson_integration(vel,f,massa,beta,anderson_freq,dt,T_ist)
        do i=1,nbody
          ekin(:,i) = 0.5 * massa * (vel(:,i))**2
        end do

      else
        !---secondo setp pos vel senza anderson--------!
        do i=1,nbody
          vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/massa
          ekin(:,i) = 0.5 * massa * (vel(:,i))**2
        end do
      
      endif

      !-----calcolo della pressione----!
      !non rompere le balle su ro beta
      !è quasi uguale
      Pist = Pist/(3 * vol)+ro/beta
      Psum = Psum + Pist
      Psum2 = Psum2 +Pist**2
      
      mekin = sum(ekin)
      
      !----salvataggio dati-------!
      if (mod(it,50) == 0) then
        write(unit=1,fmt=*)it,it*dt,pos,vel
        write(unit=2,fmt=*)it,mekin(1),mepot,mekin(1)+mepot
        write(unit=3, fmt=*)it, it*dt, T_ist , &
                              Pist ,Psum/it, sqrt(Psum2/it - (psum /it)**2) !/sqrt(it)
      endif
      !----implementazione anneling-----!
      if (anneal) vel=alfa*vel

      !-----percentuale----!
      
      if( mod(it,nstep/10) == 0) print*, floor(it/(nstep*1.0)*100)

    else
      !---- cosa a caso che non verrà mai usata-----!
      do i = 1,nbody
        pos(:,i) = pos(:,i) + f(:,i)/massa * dt**2
      end do
      
      call interazione(pos,nbody,f,mepot, side)
    end if

  !-----sanapshot per il video-----!
  call snapshot(snap_shot_name, atom_type, pos, comment, .true.)


  end do

  !-----chiusura di tutti i file aperti------!
  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"

  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"

  close(unit=3, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 3"
  
  
  


  !salviamo le distanze fra le particelle
  !do i=1,nbody
  !  do j=i+1,nbody
  !    write(unit=3,fmt=*) i,j,sqrt(dot_product(pos(:,i)-pos(:,j),pos(:,i)-pos(:,j)))
  !  end do
  !end do

  !-------salviamo le posizioni finali per un eventuale nuovo start------!
  open(unit=4, file='pos.dat', iostat=ios,  status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file pos.dat"
    
  write(unit=4,fmt=*) pos

  close(unit=4, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 4"

  !-------salviamo le velocità finali per un eventuale nuovo start-------!
  open(unit=4, file='vel.dat', iostat=ios,  status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file vel.dat"
  
  write(unit=4,fmt=*) vel

  close(unit=4, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 4"
  
  
  !snapshot situazione finale
  !call snapshot(snap_shot_name, atom_type, pos, comment, .false.)

  !scriviamo i dati finali ottenuti interessanti
  if ( anneal ) then
    print*, 'energia finale:'
    print*, 'potenziale:', mepot
    print*, 'cinetica:', mekin(1)
    print*, 'totale:', abs(mekin(1)+mepot)
  end if

  !mie deallocazioni
  if (allocated(atom_type)) deallocate(atom_type, stat=err)
  if (err /= 0) print *, "atom_type: Deallocation request denied"

end program main