
module findenes

use random
use cortes

implicit none


CONTAINS

!***************************************************

!ini_seq = punto a partir de donde buscar
!ini     = inicio  de las enes - 1 (tal como se insertan los cortes)
!fin     = fin de las enes

FUNCTION findnext_n(adn,size,ini_seq,ini,fin,idum)
  integer*1  findnext_n,salida
  integer*1  adn(:)
  integer*4  i0,i,ini,fin,size,suma,ini_seq,idum


  idum=-4214231
  i0=ini_seq
  777 continue
  ini=-1
  do i=i0,size
    if (adn(i)==-1) then
      ini=i-1
      exit
    endif
  enddo

  if (ini/=-1) then
    salida=1
    do i=ini+1,size
      if (adn(i)/=-1) then
        fin=i-1
        exit
      endif
      if (i==size) fin=size
    enddo
  else
    fin=0
    salida=0
  endif

  !si el trozo de n's tiene 10 o menos se sustituyen por nucleótidos
  !aleatorios y se vuelve al principio a seguir buscando
  !if (salida/=0.and.fin-ini<=10.and.fin/=size) then
  if (salida/=0.and.fin-ini<=3000.and.fin/=size) then  !ampliacion a 3kb, el nivel habitual de coarse graining

    do i=ini+1,fin
      adn(i)=ran3(idum)*2
    enddo
    i0=fin+1
    salida=0
    goto 777
  endif
  findnext_n=salida


end FUNCTION findnext_n

!********************************************


!**************************************************


SUBROUTINE insert_enes(adn,p0,size,haygaps,idum)

integer*4    size
integer*1    adn(:)
type(corte),pointer::p,p0


integer*4    ini,ini_n,fin_n,hayenes,haygaps
integer*4    ftot(0:1),frec(0:1),idum

ini=1
p=>p0
ftot=0
haygaps=0
do
   hayenes=findnext_n(adn,size,ini,ini_n,fin_n,idum)
   if (hayenes/=1) exit
   !print*, ini_n+1,fin_n
   frec(1)=sum(1.d0*adn(ini:ini_n))
   frec(0)=ini_n-ini+1-frec(1)
   ftot=ftot+frec
   if (ini_n/=0) then
      call insert(p,ini_n,ftot)
      p=>p%next
   endif
   haygaps=1
   p%intact=1           !contain n's and shouldn't be cut

   if (fin_n/=size) then
     call insert(p,fin_n,ftot)
     p=>p%next
   endif
   ini=fin_n+1
enddo

!p=>first
!do while (associated(p%next))
!  write(*,'(5i7,i2)') p%posic+1,p%next%posic,p%next%posic-p%posic,p%next%frec-p%frec,p%intact
!  p=>p%next
!enddo

end SUBROUTINE insert_enes



!***************************************************






end MODULE findenes
