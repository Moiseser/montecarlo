Module Functions_Variables 
  implicit none 
  double precision , parameter :: epsilon_LJ=1d0, sigma_LJ=1d0
  double precision , allocatable :: Q(:,:), U(:,:)
  character(len=2) ,allocatable :: atom_type(:)
  integer :: Natoms
contains
  
  function LJ(Q1,Q2) !Calculating Energy using Lennard Jones Potential 
    implicit none 
    double precision:: Q1(3),Q2(3),LJ
    LJ=1/sum((Q1(:)-Q2(:))**2)**3
    LJ=4*(LJ**2-LJ)
  end function LJ
  
   
  function random_integer(Nmin,Nmax) !Integer between set min and max
    implicit none 
    integer::i,Nmin,Nmax,random_integer
    Double Precision :: r1
    call random_number(r1)
    random_integer=floor(r1*(Nmax-Nmin+1)) + Nmin
  end function random_integer
    
end module Functions_Variables

program Monte_Carlo
  use Functions_Variables
  implicit none
    integer:: m,N_move,i,j,nn,NMC,tt
    Double Precision :: Q_save(3),Delta_E,r1,s(3),Rmax,T
    Double Precision :: coeff=0.01
    integer :: accepts=0
    Double Precision, allocatable :: U_new(:)
    character*30 filename
    !Initialize

    !Reads # of atoms,xyzfile, and number of Monte Carlo steps
    open(unit=2,file='input')
    read(2,*) Natoms
    read(2,*) filename
    read(2,*) NMC
    !Allocates the size of matrix
    allocate (U_new(Natoms),U(Natoms,Natoms),Q(3,Natoms),atom_type(Natoms))
    !Reading XYZ file
    open(unit=15,file=filename)
    !Checks if XYZ file #atoms = input # 
    read(15,*) i
    if(i.ne.Natoms) stop 'Wrong Natoms'
    read(15,*)
    !Reads Atom type and coordinates 
    do i=1,Natoms
       read(15,*) atom_type(i), Q(:,i)
   
    enddo
    Rmax=2.5
    T=.1
    do tt,10
   

   write(*,*) 'NMC:', NMC 

    
    do nn=1,NMC
       N_move=random_integer(1,Natoms) !Random particle to move
       Q_save(:)=Q(:,N_move) !Asigns Particle with coordinates
       call random_number(s)
       Q(:,N_move)=Q(:,N_move)+coeff*(s(:)-0.5)
       s=sum(Q(:,:),DIM=2)/Natoms
       do i=1,Natoms
          if ((sum((Q(:,i)-s(:))**2))>Rmax**2) then
             Q(:,N_move)=Q_save(:)
             goto 1
          end if
       enddo
       Delta_E=0d0
       do m=1,Natoms
          if (m.ne.N_move) then
             U_new(m)=LJ(Q(:,N_move),Q(:,m))
             Delta_E=Delta_E+U_new(m)-U(m,N_move)
          else
             U_new(m)=0d0
          end if
       end do
       call random_number(r1)
       if (dexp(-Delta_E/T)>r1) then !accept
          U(:,N_move)=U_new(:)
          U(N_move,:)=U_new(:)
          accepts=accepts+1
       else
          Q(:,N_move)=Q_save(:)
       end if
1      if(mod(nn,100) == 0) then
          write (*,*) 'Testing for accept rate'
          if(accepts/100.>0.5) then
             write (*,*) 'modifying up' 
             coeff=coeff*1.01
          else
             write (*,*) 'modifying down'
             coeff=coeff*0.99
          endif
          accepts=0


       endif
       
   write (*,*) 'Newpoint:', Q_save(:) 
   write (*,*) 'Delta_E:' , Delta_E
   write (*,*) 'Atom_Moved:', N_move
enddo
T=T+.05
enddo
  end program Monte_Carlo
     
     
