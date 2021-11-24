!CREATED BY MATHEUS B. KIATAKI (kiatakimatheus@gmail.com).
MODULE baci
  IMPLICIT NONE 

  integer :: i, j, u, w, ke=0, a=5, nunit, ntrn
  integer :: numorb, lumo, nelec, stat
  integer :: numfunccart, numleitura, resto
  integer :: numorbAp=0, numorbApp=0 
  integer :: numorbAG=0, numorbAU=0, numorbBG=0, numorbBU=0
  integer :: numorbA1=0, numorbA2=0, numorbB1=0, numorbB2=0
  integer :: numorbB1G=0, numorbB2G=0, numorbB3G=0, numorbB1U=0, numorbB2U=0, numorbB3U=0
  integer :: numorbA=0 !C1
  integer :: numorbc2A=0 !C2
  integer :: numorbc2B=0 !C2
  integer, allocatable, dimension(:) :: sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8, orbx
  integer, allocatable, dimension(:) :: sinaisci, sinaisc2, sinaisc1
  integer, allocatable, dimension(:) :: sinaisproduto, sinaisyz, sinaisxy, sinaisxz

  real(16) :: ee1, ee2, ee
  real(8), allocatable, dimension(:) :: energorb, auxenergorb

  character(len=4), allocatable, dimension(:) :: orbitais, auxiliar
    
  character(len=70) :: filegmsout
  character(len=5)  :: group
  character(len=8)  :: junkgr, junkgmsout 
 
  character(len=21) :: findee1, findee2    
  
  character(len=12) :: findeigenvec="EIGENVECTORS"
  character(len=12) :: auxeigenvec
  character(len=26) :: findwavefunc="WAVEFUNCTION NORMALIZATION"
  character(len=26) :: auxwavefunc
  character(len=49) :: findsmc='SMC WARNING: NUMBER OF ELECTRONS CHANGED (PP RUN)'
  character(len=49) :: auxsmc
  character(len=39) :: findnelec='NUMBER OF ELECTRONS                   ='
  character(len=39) :: auxnelec
  character(len=46) :: findstatistics=' STATISTICS FOR GENERATION OF SYMMETRY ORBITAL'
  character(len=46) :: auxstatistics
  character(len=37) :: findnumfunccart=' NUMBER OF CARTESIAN ATOMIC ORBITALS='
  character(len=37) :: auxnumfunccart
  character(len=40) :: findnumorb=' TOTAL NUMBER OF MOS IN VARIATION SPACE='
  character(len=40) :: auxnumorb

  character(len=30) :: saidasym
END MODULE baci
!--------------------------------------------------------------------------------

PROGRAM gmssmc   
  use baci

  read(5,'(a8,1x,a5)') junkgr, group
  if(junkgr .ne. 'GROUP  =') stop'ABORT, flaggroup .ne. GROUP ='
  
  read(5,'(a8,a70)') junkgmsout, filegmsout
  if(junkgmsout .ne. 'GMSOUT =') stop'ABORT, junkgmsout .ne. GMSOUT =' 

  call readgmssmcout
   
  if (group .eq. 'Cs') then
    call CS
    open(41,file='sim-Ap')
    open(42,file='sim-App')
    ntrn=1
    do nunit=41,42
      call simetrias 
    end do
  else if (group .eq. 'Cnh 2') then
    call C2H
    open(43,file='sim-Ag')
    open(44,file='sim-Au')
    open(45,file='sim-Bg')
    open(46,file='sim-Bu')
    ntrn=2
    do nunit=43,46
      call simetrias
    end do
  else if (group .eq. 'Cnv 2') then
    call C2V
    open(47,file='sim-A1')
    open(48,file='sim-A2')
    open(49,file='sim-B1')
    open(50,file='sim-B2')
    ntrn=2
    do nunit=47,50
      call simetrias
    end do 
  else if (group .eq. 'Dnh 2') then
    call D2H
    open(51,file='sim-Ag')
    open(52,file='sim-Au')
    open(53,file='sim-B1g')
    open(54,file='sim-B2g')
    open(55,file='sim-B3g')
    open(56,file='sim-B1u')
    open(57,file='sim-B2u')
    open(58,file='sim-B3u')
    ntrn=3
    do nunit=51,58
      call simetrias
    end do 
  else if (group .eq.'C1') then
    call C1
    open(nunit,file='sim')
    ntrn=1
    nunit=59
    call simetrias
  else if (group .eq. 'Cn 2' .or. group .eq. 'cn 2') then
    call C2
    open(60,file='sim-A')
    open(61,file='sim-B')
    ntrn=1
    do nunit=60,61
      call simetrias
    end do
  else
    stop 'ABORT, YOU PUT AN UNAVAILABLE POINT GROUP'
  end if

end program gmssmc

!----------------SUBROUTINES------------------------------------------ 
subroutine readgmssmcout
  use baci 

  open (unit=99, file=filegmsout)

  do 
    read(99,'(a49)',iostat=stat) auxsmc
      if (stat .gt. 0) stop'ABORT, PROBLEM TO READ GMSOUTPUT FILE [SMC WARNING]'
      if (stat .lt. 0) stop'ABORT, END OF GMS OUTPUT FILE HAS REACHED, IT WAS NOT FOUND [SMC WARNING]'
      if (auxsmc .eq. findsmc) exit
  end do
  
  read(99,'(a39,i6)',iostat=stat) auxnelec, nelec
  if (stat .gt. 0) stop'ABORT, PROBLEM TO READ GMS OUTPUT FILE [NUMBER OF ELECTRONS]'
  if (stat .lt. 0) stop'ABORT, END OF GMS OUTPUT FILE HAS REACHED, IT WAS NOT FOUND [NUMBER OF ELECTRONS]'
  if (auxnelec .ne. findnelec) stop'ABORT, AUXNELEC .NE. FINDNELEC'

  lumo = (nelec/2) + 1

  do 
    read(99,fmt='(a46)',iostat=stat) auxstatistics
      if (stat .gt. 0) stop'ABORT, PROBLEM TO READ GMS OUTPUT FILE [STATISTICS<--->(ISPHER=1)]'
      if (stat .lt. 0) stop'ABORT, END OF GMS OUTPUT FILE HAS REACHED, IT WAS NOT FOUND [STATISTICS<--->(ISPHER=1)]'
      if (auxstatistics .eq. findstatistics) exit
  end do    

  read(99,'(a37,i11)',iostat=stat) auxnumfunccart, numfunccart
  if (stat .gt. 0) stop'ABORT, PROBLEM TO READ GMS OUTPUT FILE [NUMBER OF CARTESIAN ATOMIC ORBITALS]'
  if (stat .lt. 0) stop'ABORT, END OF GMS OUTPUT FILE HAS REACHED, IT WAS NOT FOUND [NUMBER OF CARTESIAN ATOMIC ORBITALS]'
  if (auxnumfunccart .ne. findnumfunccart) stop'ABORT, AUXNUMFUNCCART .NE. FINDNUMFUNCCART'

  do i=1,2
    read (99,*)
  end do

  read(99,'(a40,i8)',iostat=stat) auxnumorb, numorb
  if (stat .gt. 0) stop'ABORT, PROBLEM TO READ GMS OUTPUT FILE [TOTAL NUMBER OF MOS IN VARIATION SPACE]'
  if (stat .lt. 0) stop'ABORT, END OF GMS OUTPUT FILE HAS REACHED, IT WAS NOT FOUND [TOTAL NUMBER OF MOS IN VARIATION SPACE]'
  if (auxnumorb .ne. findnumorb) stop'ABORT, AUXNUMORB .NE. FINDNUMORB' 

  numleitura = numorb/5
  resto = MOD(numorb,5)    

  allocate(orbitais(numorb))
  allocate(auxiliar(a))

  allocate(sinaisxy(numfunccart)) !Cs, C2h
  allocate(sinaisc2(numorb)) !C2, C2v
  allocate(sim1(numorb)) 
  allocate(sim2(numorb)) 

  allocate(energorb(numorb))
  allocate(auxenergorb(a))  

  if (group.eq.'Cnh 2') then 
    allocate(sinaisci(numorb))
    allocate(sim3(numorb))
    allocate(sim4(numorb))
  else if (group.eq.'Cnv 2') then
    allocate(sinaisyz(numfunccart))
    allocate(sinaisproduto(numfunccart))
    allocate(sim3(numorb))
    allocate(sim4(numorb))
  else if (group.eq.'Dnh 2') then
    allocate(sinaisyz(numfunccart))
    allocate(sinaisxz(numfunccart))
    allocate(sinaisproduto(numfunccart))
    allocate(sim3(numorb))
    allocate(sim4(numorb))
    allocate(sim5(numorb))
    allocate(sim6(numorb))
    allocate(sim7(numorb))
    allocate(sim8(numorb))
  else if (group.eq.'C1') then
    allocate(sinaisc1(numorb))
  end if

  do 
    read(99,*) auxeigenvec
    if (auxeigenvec.eq.findeigenvec) exit
  end do  
  
  do i=1, 3
    read(99,*)
  end do

  do u=1, numleitura
       read(99,*) (auxenergorb(i), i=1,a)
       read(99,*) (auxiliar(i), i=1,a)
       do w=1, numfunccart+2
         read(99,*)
       end do 
       do j=1,a
         energorb(ke+j)=auxenergorb(j)
         orbitais(ke+j)=auxiliar(j)
       end do 
       ke=ke+5    
  end do

  if (resto .ne. 0) then    
       read(99,*) (auxenergorb(i), i=1,resto)
       read(99,*) (auxiliar(i), i=1,resto)    
       do j=1,resto
         energorb(ke+j)=auxenergorb(j)
         orbitais(ke+j)=auxiliar(j)
       end do
  end if
   
  do 
       read(99,fmt='(9x,a26)') auxwavefunc
       if (auxwavefunc .eq. findwavefunc) exit
  end do

  read(99,*)
  read(99,fmt='(16x,a21,F25.10)') findee1, ee1
  if(findee1 .ne. 'ONE ELECTRON ENERGY =')stop'ABORT, PROBLEM TO READ/FOUND ONE ELECTRON ENERGY'
  read(99,fmt='(16x,a21,F25.10)') findee2, ee2
  if(findee2 .ne. 'TWO ELECTRON ENERGY =')stop'ABORT, PROBLEM TO READ/FOUND TWO ELECTRON ENERGY'

  ee=ee1+ee2
 
end subroutine readgmssmcout
!================================================================================
subroutine CS ()
  use baci
 
 do i=1,numorb    
   if (orbitais(i) .eq. "A'" .or. orbitais(i) .eq. "?A'") then
      sinaisxy(i) = i      
      if (i .ge. lumo) then
        sim1(numorbAp+1)= i
        numorbAp=numorbAp+1
      end if
   elseif (orbitais(i) .eq. "A''" .or. orbitais(i) .eq. "?A''") then
      sinaisxy(i) = -i
      if (i .ge. lumo) then
        sim2(numorbApp+1)=i
        numorbApp=numorbApp+1
      end if
   else 
      stop'ABORT, THE LABEL OF ORBITALS AND POINT GROUP MISMATCH'
   end if
 end do

 do i=(numorb+1),numfunccart
  sinaisxy(i) = i
 end do

 open(1,file='orbApfmt.txt')
 write(1,*) numorbAp
 write(1, fmt='(20I4)')(sim1(i), i=1, numorbAp)
 close(1)

 open(2,file='orbAppfmt.txt')
 write(2,*) numorbApp
 write(2, fmt='(20I4)')(sim2(i), i=1, numorbApp)
 close(2)

 open(3, file="blocosinaisXY.txt")
 write(3,fmt='(15I5)')(sinaisxy(i), i=1, numfunccart)
 close(3)

end subroutine CS
!================================================================================
 subroutine C2H () 
  use baci
!vc fez com base nos hidrocarbonetos, o plano de reflexao aqui é XY. Vc nao fixou nenhum plano no gamess.

 do i=1,numorb    
   if (orbitais(i) .eq. "AG" .or. orbitais(i) .eq. "?AG") then
      sinaisxy(i) = i
      sinaisci(i) = i      
      if (i .ge. lumo) then
        sim1(numorbAG+1)= i
        numorbAG=numorbAG+1
      end if
                                                    
   elseif (orbitais(i) .eq. "AU" .or. orbitais(i) .eq. "?AU") then
      sinaisxy(i) = -i
      sinaisci(i)= -i
      if (i .ge. lumo) then
        sim2(numorbAU+1)=i
        numorbAU=numorbAU+1
      end if

   elseif (orbitais(i) .eq. "BG" .or. orbitais(i) .eq. "?BG") then
      sinaisxy(i) = -i
      sinaisci(i) = i
      if (i .ge. lumo) then
        sim3(numorbBG+1)=i
        numorbBG=numorbBG+1
      end if

   elseif (orbitais(i) .eq. "BU" .or. orbitais(i) .eq. "?BU") then
      sinaisxy(i) = i
      sinaisci(i) = -i
      if (i .ge. lumo) then
        sim4(numorbBU+1)=i
        numorbBU=numorbBU+1
      end if

    else 
       stop'ABORT, THE LABEL OF ORBITALS AND POINT GROUP MISMATCH'

   end if
 end do

 do i=(numorb+1),numfunccart
  sinaisxy(i) = i
 end do

 open(4, file='orbAGfmt.txt')
 write(4,*) numorbAG
 write(4, fmt='(20I4)')(sim1(i), i=1, numorbAG)
 close(4)

 open(7, file='orbAUfmt.txt')
 write(7,*) numorbAU
 write(7, fmt='(20I4)')(sim2(i), i=1, numorbAU)
 close(7)

 open(8, file='orbBGfmt.txt')
 write(8,*) numorbBG
 write(8, fmt='(20I4)')(sim3(i), i=1, numorbBG)
 close(8)

 open(9, file='orbBUfmt.txt')
 write(9,*) numorbBU
 write(9, fmt='(20I4)')(sim4(i), i=1, numorbBU)
 close(9)


 open(10, file="blocosinaisXY.txt")
 write(10,fmt='(15I5)')(sinaisxy(i), i=1, numfunccart)
 close(10)

 open (11, file="blocosinaiscentrodeinv.txt")
 write(11,fmt='(15I5)')(sinaisci(i), i=1, numorb)
 close(11)

 end subroutine C2H
!================================================================================
 subroutine C2V ()
 use baci
!!PLANO DE SIMETRIA YZ, XY, YZ*XY e EIXO DE SIMETRIA C2(Y)!! VC se baseou nos sinais do gms do cbr4 com o plano YZ fixado no gamess.

 do i=1,numorb    
   if (orbitais(i) .eq. "A1" .or. orbitais(i) .eq. "?A1") then
      sinaisyz(i) = i
      sinaisc2(i) = i  
      sinaisxy(i) = i 
      if (i .ge. lumo) then
        sim1(numorbA1+1)= i
        numorbA1=numorbA1+1
      end if
                                                       
   elseif (orbitais(i) .eq. "A2" .or. orbitais(i) .eq. "?A2") then
      sinaisyz(i) = -i
      sinaisc2(i)= i
      sinaisxy(i) =-i
      if (i .ge. lumo) then
        sim2(numorbA2+1)=i
        numorbA2=numorbA2+1
      end if

   elseif (orbitais(i) .eq. "B1" .or. orbitais(i) .eq. "?B1") then
      sinaisyz(i) =  i
      sinaisc2(i) = -i
      sinaisxy(i) = -i
      if (i .ge. lumo) then
        sim3(numorbB1+1)=i
        numorbB1=numorbB1+1
      end if

   elseif (orbitais(i) .eq. "B2" .or. orbitais(i) .eq. "?B2") then
      sinaisyz(i) =-i
      sinaisc2(i) = -i
      sinaisxy(i) = i
      if (i .ge. lumo) then
        sim4(numorbB2+1)=i
        numorbB2=numorbB2+1
      end if

    else 
       stop'ABORT, THE LABEL OF ORBITALS AND POINT GROUP MISMATCH'

   end if
 end do

 do i = (numorb+1), numfunccart
   sinaisyz(i) = i
   sinaisxy(i) = i
 end do

 do i=1,numfunccart
   sinaisproduto(i) = ( sinaisyz(i)*sinaisxy(i) / IABS(sinaisyz(i)) ) 
 end do

 open(12, file='orbA1fmt.txt')
 write(12,*) numorbA1
 write(12, fmt='(20I4)')(sim1(i), i=1, numorbA1)
 close(12)
  
 open(13, file='orbA2fmt.txt')
 write(13,*) numorbA2
 write(13, fmt='(20I4)')(sim2(i), i=1, numorbA2)
 close(13)

 open(14, file='orbB1fmt.txt')
 write(14,*) numorbB1
 write(14, fmt='(20I4)')(sim3(i), i=1, numorbB1)
 close(14)

 open(15, file='orbB2fmt.txt')
 write(15,*) numorbB2
 write(15, fmt='(20I4)')(sim4(i), i=1, numorbB2)
 close(15)

 open(16, file="blocosinaisYZ.txt")
 write(16,fmt='(15I5)')(sinaisyz(i), i=1, numfunccart)
 close(16)

 open(17, file="blocosinaisXY.txt")
 write(17,fmt='(15I5)')(sinaisxy(i), i=1, numfunccart)
 close(17)

 open(18, file="all3planes.txt")
 write(18,fmt='(15I5)')(sinaisyz(i), i=1, numfunccart)
 write(18,fmt='(15I5)')(sinaisxy(i), i=1, numfunccart)
 write(18,fmt='(15I5)')(sinaisproduto(i), i=1, numfunccart)
 close(18)

end subroutine C2V
!!=============================================================================== 
 subroutine D2H ()
   use baci
!vc alterou com base nos sinais dos orbitais do arquivo gamess da giseli, que fixou o plano YZ no gamess. 11/09/2018

 do i=1,numorb    
   if (orbitais(i) .eq. "AG" .or. orbitais(i) .eq. "?AG") then
      sinaisyz(i) = i
      sinaisxy(i) = i
      sinaisxz(i) = i      
      if (i .ge. lumo) then
        sim1(numorbAG+1)= i
        numorbAG=numorbAG+1
      end if
                                                   
   elseif (orbitais(i) .eq. "AU" .or. orbitais(i) .eq. "?AU") then
      sinaisyz(i) = -i
      sinaisxy(i) = -i
      sinaisxz(i) = -i 
      if (i .ge. lumo) then
        sim2(numorbAU+1)=i
        numorbAU=numorbAU+1
      end if

   elseif (orbitais(i) .eq. "B1G" .or. orbitais(i) .eq. "?B1G") then
      sinaisyz(i) = -i
      sinaisxy(i) = -i
      sinaisxz(i) =  i 
      if (i .ge. lumo) then
        sim3(numorbB1G+1)=i
        numorbB1G=numorbB1G+1
      end if

   elseif (orbitais(i) .eq. "B2G" .or. orbitais(i) .eq. "?B2G") then
      sinaisyz(i) =  i
      sinaisxy(i) = -i
      sinaisxz(i) = -i 
      if (i .ge. lumo) then
        sim4(numorbB2G+1)=i
        numorbB2G=numorbB2G+1
      end if
  
   elseif (orbitais(i) .eq. "B3G" .or. orbitais(i) .eq. "?B3G") then
      sinaisyz(i) = -i
      sinaisxy(i) =  i
      sinaisxz(i) = -i 
      if (i .ge. lumo) then
        sim5(numorbB3G+1)=i
        numorbB3G=numorbB3G+1
      end if
 
   elseif (orbitais(i) .eq. "B1U" .or. orbitais(i) .eq. "?B1U") then
      sinaisyz(i) =  i
      sinaisxy(i) =  i
      sinaisxz(i) = -i 
      if (i .ge. lumo) then
        sim6(numorbB1U+1)=i
        numorbB1U=numorbB1U+1
      end if
 
   elseif (orbitais(i) .eq. "B2U" .or. orbitais(i) .eq. "?B2U") then
      sinaisyz(i) = -i
      sinaisxy(i) =  i
      sinaisxz(i) =  i 
      if (i .ge. lumo) then
        sim7(numorbB2U+1)=i
        numorbB2U=numorbB2U+1
      end if
 
   elseif (orbitais(i) .eq. "B3U" .or. orbitais(i) .eq. "?B3U") then
      sinaisyz(i) =  i
      sinaisxy(i) = -i
      sinaisxz(i) =  i 
      if (i .ge. lumo) then
        sim8(numorbB3U+1)=i
        numorbB3U=numorbB3U+1
      end if

    else 
       stop'ABORT, THE LABEL OF ORBITALS AND POINT GROUP MISMATCH'

   end if
 end do
 
 do i=(numorb+1),numfunccart
   sinaisyz(i) = i
   sinaisxy(i) = i
   sinaisxz(i) = i 
 end do

 do i=1,numfunccart
   sinaisproduto(i) = ( sinaisyz(i)*sinaisxy(i) )/iabs(sinaisyz(i)) 
 end do

 open(19, file='orbAGfmt.txt')
 write(19,*) numorbAG
 write(19, fmt='(20I4)')(sim1(i), i=1, numorbAG)
 close(19)

 open(20, file='orbAUfmt.txt')
 write(20,*) numorbAU
 write(20, fmt='(20I4)')(sim2(i), i=1, numorbAU)
 close(20)

 open(21, file='orbB1Gfmt.txt')
 write(21,*) numorbB1G
 write(21, fmt='(20I4)')(sim3(i), i=1, numorbB1G)
 close(21)

 open(22, file='orbB2Gfmt.txt')
 write(22,*) numorbB2G
 write(22, fmt='(20I4)')(sim4(i), i=1, numorbB2G)
 close(22)

 open(23, file='orbB3Gfmt.txt')
 write(23,*) numorbB3G
 write(23, fmt='(20I4)')(sim5(i), i=1, numorbB3G)
 close(23)

 open(24, file='orbB1Ufmt.txt')
 write(24,*) numorbB1U
 write(24, fmt='(20I4)')(sim6(i), i=1, numorbB1U)
 close(24)

 open(25, file='orbB2Ufmt.txt')
 write(25,*) numorbB2U
 write(25, fmt='(20I4)')(sim7(i), i=1, numorbB2U)
 close(25)

 open(26, file='orbB3Ufmt.txt')
 write(26,*) numorbB3U
 write(26, fmt='(20I4)')(sim8(i), i=1, numorbB3U)
 close(26)

 open(27, file="blocosinaisYZ.txt")
 write(27,fmt='(15I5)')(sinaisyz(i), i=1, numfunccart)
 close(27)

 open(28, file="blocosinaisXY.txt")
 write(28,fmt='(15I5)')(sinaisxy(f), f=1, numfunccart)
 close(28)

 open(29, file="all3planes.txt")
 write(29,fmt='(15I5)')(sinaisyz(i), i=1, numfunccart)
 write(29,fmt='(15I5)')(sinaisxy(i), i=1, numfunccart)
 write(29,fmt='(15I5)')(sinaisproduto(i), i=1, numfunccart)
 close(29)

 end subroutine D2H
!================================================================================
subroutine C1 ()
  use baci
 
 do i=1,numorb    
   if (orbitais(i) .eq. "A" .or. orbitais(i) .eq. "?A") then 
      sinaisc1(i) = i
      if (i .ge. lumo) then
        numorbA=numorbA+1
        sim1(numorbA)= i
      end if
   else 
      stop'ABORT, THE LABEL OF THE ORBITALS THAT BELONGS TO C1 GROUP ARE NOT EQUAL "A"'
   end if
 end do
 
 open(30, file="orbfmt.txt")
 write(30,fmt='(15I5)')(sim1(i), i=1, numorbA)
 close(30)

 open(31, file="blocosinaisc1.txt")
 write(31,fmt='(15I5)')(sinaisc1(i), i=1, numorb)
 close(31)

end subroutine C1
!=================================================================================
subroutine C2 ()
  use baci

  do i=1,numorb
    if (orbitais(i) .eq. 'A' .or. orbitais(i) .eq. '?A') then
       sinaisc2(i) = i
       if (i .ge. lumo) then
          numorbc2A=numorbc2A+1
          sim1(numorbc2A)=i
       end if
    else if (orbitais(i) .eq. 'B' .or. orbitais(i) .eq. '?B') then
       sinaisc2(i) = -i
       if (i .ge. lumo) then
          numorbc2B=numorbc2B+1
          sim2(numorbc2B)=i
       end if
    else 
       stop'ABORT, THE LABEL OF THE ORBITALS AND POINT GROUP MISMATCH'
    end if
  end do

 open(32, file="orb-A-fmt.txt")
 write(32,fmt='(15I5)')(sim1(i), i=1, numorbc2A)
 close(32)

 open(33, file="orb-B-fmt.txt")
 write(33,fmt='(15I5)')(sim2(i), i=1, numorbc2B)
 close(33)

 open(34, file="orb-B-fmt.txt")
 write(34,fmt='(15I5)')(sinaisc2(i), i=1,numorb)
 close(34)

end subroutine C2
!===============================================================================
subroutine simetrias
  use baci 
!--------------------------------instrucoes-----------------------------------------
  write(nunit,'(a)') trim(adjustl('#!/bin/csh'))
  write(nunit,*)
  write(nunit,'(a)') trim(adjustl('set exec = ~/simetrias/build_csfs/build_csfs.x'))

  write(nunit,'(a)') trim(adjustl('# %INI_CONTROL_BLOCK'))
  write(nunit,'(a)') trim(adjustl('# PROJEC =-1 (electron), +1 (positron)'))
  write(nunit,'(a)') trim(adjustl('# NFUNCS = number of molecular orbitals'))
  write(nunit,'(a)') trim(adjustl('# NOCCUP = number of occupied orbitals'))
  write(nunit,'(a)') trim(adjustl('# CISTAR = YES (multichannel), NO (elastic only)'))
  write(nunit,'(a)') trim(adjustl('# NBLOCK = number of blocks of orbitals'))
  write(nunit,'(a)') trim(adjustl('# %END_CONTROL_BLOCK'))
  
  write(nunit,*)
 
  write(nunit,'(a)') trim(adjustl('# %INI_ORBITAL_BLOCK'))
  write(nunit,'(a)') trim(adjustl('# NOSYMM = number of symmetries'))
  write(nunit,'(a)') trim(adjustl('# For each symmetry, list of how each orbital transforms'))
  write(nunit,'(a)') trim(adjustl('# SYMM_1 = sign that define the desired symmetry (same for SYMM_2, SYMM_3, one flag per line)'))
  write(nunit,'(a)') trim(adjustl('# ENERGIES'))
  write(nunit,'(a)') trim(adjustl('# list of orbital energies'))
  write(nunit,'(a)') trim(adjustl('# %END_ORBITAL_BLOCK'))

  write(nunit,'(a)') trim(adjustl('# loop over number of blocks'))

  write(nunit,'(a)') trim(adjustl('# %INI_SET_1_BLOCK (%INI_SET_2_BLOCK, etc.)'))
  write(nunit,'(a)') trim(adjustl('# NHOLES = number of hole orbitals'))
  write(nunit,'(a)') trim(adjustl('# list of hole orbitals'))
  write(nunit,'(a)') trim(adjustl('# NPARTS = number of particle orbitals'))
  write(nunit,'(a)') trim(adjustl('# list of particle orbitals'))
  write(nunit,'(a)') trim(adjustl('# NSCATS = number of scattering orbitals'))
  write(nunit,'(a)') trim(adjustl('# list of scattering orbitals'))
  write(nunit,'(a)') trim(adjustl('# EPSCUT = energy cut'))
  write(nunit,'(a)') trim(adjustl('# %END_SET_1_BLOCK (%END_SET_2_BLOCK, etc.)'))

  write(nunit,'(a)') trim(adjustl('# This block only exists when CISTAR = YES'))

  write(nunit,'(a)') trim(adjustl('# %INI_CISTARG_BLOCK'))
  write(nunit,'(a)') trim(adjustl('# NPAIRS = number of CIS pairs'))
  write(nunit,'(a)') trim(adjustl('# HOLES'))
  write(nunit,'(a)') trim(adjustl('# list of holes'))
  write(nunit,'(a)') trim(adjustl('# PARTS'))
  write(nunit,'(a)') trim(adjustl('# list of particles'))
  write(nunit,'(a)') trim(adjustl('# NOTRIP = number of triplet states'))
  write(nunit,'(a)') trim(adjustl('# NOSING = number of singlet states'))
  write(nunit,'(a)') trim(adjustl('# %END_CISTARG_BLOCK'))
!---------------fim das instrucoes--------------------------------------------

  write(nunit,'(a)') trim(adjustl('$exec << fim > sym.lis'))

  write(nunit,'(a)') trim(adjustl('%INI_CONTROL_BLOCK'))
  write(nunit,'(a)') trim(adjustl('PROJEC = -1'))
  write(nunit,'(a,i5)') trim(adjustl('NFUNCS =')), numorb
  write(nunit,'(a,i5)') trim(adjustl('NOCCUP =')), (nelec/2)
  write(nunit,'(a)') trim(adjustl('CISTAR = NO'))
  write(nunit,'(a)') trim(adjustl('NBLOCK = 1'))
  write(nunit,'(a)') trim(adjustl('%END_CONTROL_BLOCK'))

  write(nunit,*)
      
  write(nunit,'(a)') trim(adjustl('%INI_ORBITAL_BLOCK'))
  write(nunit,'(a,i5)') trim(adjustl('NOSYMM =')), ntrn

  if(group.eq.'Cs') then
   write(nunit,'(15(i5,1x))') (sinaisxy(i), i=1,numorb)
  else if(group.eq.'Cnh 2') then
   write(nunit,'(15(i5,1x))') (sinaisxy(i), i=1,numorb)
   write(nunit,'(15(i5,1x))') (sinaisci(i), i=1,numorb)
  else if(group.eq.'Cnv 2') then
   write(nunit,'(15(i5,1x))') (sinaisyz(i),i=1,numorb)  
   write(nunit,'(15(i5,1x))') (sinaisxy(i),i=1,numorb)
  else if(group.eq.'Dnh 2') then
   write(nunit,'(15(i5,1x))') (sinaisyz(i),i=1,numorb)  
   write(nunit,'(15(i5,1x))') (sinaisxy(i),i=1,numorb)
   write(nunit,'(15(i5,1x))') (sinaisxz(i),i=1,numorb)
  else if(group.eq.'C1') then
   write(nunit,'(15(i5,1x))') (sinaisc1, i=1,numorb)
  else if(group.eq.'Cn 2') then
   write(nunit,'(15(i5,1x))') (sinaisc2(i), i=1,numorb)
  end if  

  write(nunit,'(a,1x)', advance='no') trim(adjustl('SYMM_1 ='))

  if(group.eq.'Cs') then
    if(nunit.eq.41) write(nunit,'(i2)')  +1
    if(nunit.eq.42) write(nunit,'(i2)')  -1
  elseif(group.eq.'Cnh 2') then
    if(nunit.eq.43) write(nunit,'(i2)')  +1
    if(nunit.eq.44) write(nunit,'(i2)')  -1
    if(nunit.eq.45) write(nunit,'(i2)')  -1
    if(nunit.eq.46) write(nunit,'(i2)')  +1
    write(nunit,'(a,1x)', advance='no') trim(adjustl('SYMM_2 = '))
    if(nunit.eq.43) write(nunit,'(i2)')  +1
    if(nunit.eq.44) write(nunit,'(i2)')  -1
    if(nunit.eq.45) write(nunit,'(i2)')  +1
    if(nunit.eq.46) write(nunit,'(i2)')  -1
  elseif(group.eq.'Cnv 2') then
    if(nunit.eq.47) write(nunit,'(i2)')  +1
    if(nunit.eq.48) write(nunit,'(i2)')  -1
    if(nunit.eq.49) write(nunit,'(i2)')  +1
    if(nunit.eq.50) write(nunit,'(i2)')  -1
    write(nunit,'(a,1x)', advance='no') trim(adjustl('SYMM_2 = '))
    if(nunit.eq.47) write(nunit,'(i2)')  +1
    if(nunit.eq.48) write(nunit,'(i2)')  -1
    if(nunit.eq.49) write(nunit,'(i2)')  -1
    if(nunit.eq.50) write(nunit,'(i2)')  +1
  elseif(group.eq.'Dnh 2') then
    if(nunit.eq.51) write(nunit,'(i2)')  +1
    if(nunit.eq.52) write(nunit,'(i2)')  -1
    if(nunit.eq.53) write(nunit,'(i2)')  -1
    if(nunit.eq.54) write(nunit,'(i2)')  +1
    if(nunit.eq.55) write(nunit,'(i2)')  -1
    if(nunit.eq.56) write(nunit,'(i2)')  +1
    if(nunit.eq.57) write(nunit,'(i2)')  -1
    if(nunit.eq.58) write(nunit,'(i2)')  +1
    write(nunit,'(a,1x)', advance='no') trim(adjustl('SYMM_2 = '))
    if(nunit.eq.51) write(nunit,'(i2)')  +1
    if(nunit.eq.52) write(nunit,'(i2)')  -1
    if(nunit.eq.53) write(nunit,'(i2)')  -1
    if(nunit.eq.54) write(nunit,'(i2)')  -1
    if(nunit.eq.55) write(nunit,'(i2)')  +1
    if(nunit.eq.56) write(nunit,'(i2)')  +1
    if(nunit.eq.57) write(nunit,'(i2)')  +1
    if(nunit.eq.58) write(nunit,'(i2)')  -1
    write(nunit,'(a,1x)', advance='no') trim(adjustl('SYMM_3 = '))
    if(nunit.eq.51) write(nunit,'(i2)')  +1
    if(nunit.eq.52) write(nunit,'(i2)')  -1
    if(nunit.eq.53) write(nunit,'(i2)')  +1
    if(nunit.eq.54) write(nunit,'(i2)')  -1
    if(nunit.eq.55) write(nunit,'(i2)')  -1
    if(nunit.eq.56) write(nunit,'(i2)')  -1
    if(nunit.eq.57) write(nunit,'(i2)')  +1
    if(nunit.eq.58) write(nunit,'(i2)')  +1
  elseif(group.eq.'C1') then
    if(nunit.eq.59) write(nunit,'(i2)')  +1
  elseif(group.eq.'Cn 2') then
    if(nunit.eq.60) write(nunit,'(i2)')  +1
    if(nunit.eq.61) write(nunit,'(i2)')  -1
  end if
 
  write(nunit,'(a)') trim(adjustl('ENERGIES'))
!!  do i=1, numorb
!!    write(nunit,*) energorb(i)
!!  end do
  write(nunit,'(a)') trim(adjustl('%END_ORBITAL_BLOCK'))
  
  write(nunit,*)

  write(nunit,'(a)') trim(adjustl('%INI_SET_1_BLOCK')) 
  write(nunit,'(a,i5)') trim(adjustl('NHOLES = ')), (nelec/2)
  write(nunit,'(20(i5,1x))') (i, i=1,nelec/2)
  write(nunit,'(a,i5)') trim(adjustl('NPARTS = ')), (numorb-(nelec/2))
  write(nunit,'(20(i5,1x))') (i, i=lumo, numorb)
  write(nunit,'(a,i5)') trim(adjustl('NSCATS = ')), (numorb-(nelec/2))
  write(nunit,'(20(i5,1x))') (i, i=lumo, numorb)
  write(nunit,'(a,f6.2)') trim(adjustl('EPSCUT = ')), +1.0
  write(nunit,'(a)') trim(adjustl('%END_SET_1_BLOCK')) 

  write(nunit,*)

  write(nunit,'(a)') trim(adjustl('fim')) 

end subroutine simetrias
!================================================================================
!subroutine simetriasold 
!  use baci

!  write(nunit,'(a10)') '#!/bin/csh'
!  write(nunit,*)
!  write(nunit,'(a)') trim(adjustl('#1a. linha: numero total de funcoes contraidas, numero de orbitais ocupados, numero de orbitais de buraco, numero de orbitais de particula, numero de orbitais de espalhamento, sinal (-1=eletron, +1=positron).'))
!  write(nunit,'(a)') trim(adjustl( '#2a. linha: lista dos orbitais de buraco.'))
!  write(nunit,'(a)') trim(adjustl( '#3a. linha: lista dos orbitais de particula.'))
!  write(nunit,'(a)') trim(adjustl( '#4a. linha: lista dos orbitais de espalhamento.'))
!  write(nunit,'(a)') trim(adjustl( '#5a. linha: numero de operacoes de simetria.'))
!  write(nunit,'(a)') trim(adjustl( '#6a. linha: tabela de sinais dos orbitais correspondendo a cada operacao de simetria.'))
!  write(nunit,'(a)') trim(adjustl('#7a. linha: caracteres da representacao irredutivel pra cada operacao de simetria (Ag=+1 +1 +1).'))
!  write(nunit,'(a41)') "set exec = '~/simetrias/simetrias_conf.x'"

!  write(nunit,*)'$exec << fim > sym.lis'
!  write(nunit,'(6(i5,1x))') numfunccart, (nelec/2), (nelec/2), (numorb-nelec/2), (numorb-nelec/2), -1
!  write(nunit,'(20(i5,1x))') (i, i=1,(nelec/2))
!  write(nunit,'(20(i5,1x))') (i, i=lumo,numorb)
!  write(nunit,'(20(i5,1x))') (i, i=lumo,numorb)
!  write(nunit,'(i4)') ntrn
!
!  if(group.eq.'Cs') then
!   write(nunit,'(20(i5,1x))') (sinaisxy(i), i=1,numfunccart)
!   if(nunit.eq.1) write(nunit,*) +1
!   if(nunit.eq.2) write(nunit,*) -1
!  else if(group.eq.'Cn 2') then
!   write(nunit,'(20(i5,1x))') (sinaisc2(i), i=1,numfunccart)
!   if(nunit.eq.4) write(nunit,*) +1
!   if(nunit.eq.7) write(nunit,*) -1
!  else if(group.eq.'Cnv 2') then
!   write(nunit,'(20(i5,1x))') (sinaisyz(i),i=1,numfunccart)  
!   write(nunit,'(20(i5,1x))') (sinaisxy(i),i=1,numfunccart)
!   if(nunit.eq.8)  write(nunit,'(2(i2,1x))') +1, +1
!   if(nunit.eq.9)  write(nunit,'(2(i2,1x))') -1, -1
!   if(nunit.eq.10) write(nunit,'(2(i2,1x))') +1, -1
!   if(nunit.eq.11) write(nunit,'(2(i2,1x))') -1, +1
!  else if(group.eq.'C1') then
!   write(nunit,'(20(i5,1x))') (i, i=1,numorb)
!   if(nunit.eq.3)  write(nunit,*) +1
!  end if  

!  write(nunit,'(f5.2)') -2.0
!  do i=1,numorb
!    write(nunit,'(f10.5)') energorb(i)
!  end do
!  write(nunit,'(a3)')'fim'

!end subroutine
!==================================================================================
