module leeadn


CONTAINS


!**************************************************


SUBROUTINE lee_adn(adn,size,enes,fen,ini,fin,tipo,f_words)
character(*)   fen
integer*1      adn(:)
integer*4      size,letra,enes,ini,fin,cont,n,i,f_words(0:)
character*83   linea
character*32   t_bases
character*2    tipo


t_bases='ATCGNRYSWKMVHDBXatcgnryswkmvhdbx'


open (1,file=fen)
if (fin==0) fin = 1000000000

linea='kkk'

do while (linea(1:2) .ne. 'SQ' .and. linea(1:6) .ne. 'ORIGIN' .and. linea(1:1) .ne. '>')
   read (1,'(A83)') linea
enddo

size=0
enes=0
cont=0
f_words=0

1  read(1,'(A83)',end=10) linea

  ! if (size>103230000) print*, trim(linea)
   n=len_trim(linea)

   do i=1,n
      letra=index(t_bases,linea(i:i))

      if (letra>0) then
        cont=cont+1
        if (cont>fin) goto 10
        if (cont>=ini) then
          size=size+1
          !if (size>103230000) print*,size

          if (letra>16) letra=letra-16
          if (letra>=5) then
             letra=5
             enes=enes+1
          endif
          if (letra>4) then
     letra=-1
          else
             if (tipo=='SW') then
        if ((letra==1).or.(letra==2)) then
                  letra=0
                else
                  letra=1
                endif
              else if (tipo=='RY') then
                if ((letra==1).or.(letra==4)) then
                  letra=0
                else
                  letra=1
                endif
              else
                letra=letra-1
              endif

              f_words(letra)=f_words(letra)+1
          endif
          adn(size)=letra
        endif
      endif
   enddo
   go to 1
10 continue
rewind(1)
close(1)

end SUBROUTINE lee_adn

!****************************************************************


end MODULE leeadn
