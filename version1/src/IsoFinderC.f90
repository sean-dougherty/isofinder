!***************************************************************************

!     Computational identification of Isochores

!     Pedro Bernaola(1), Pedro Carpena(1) & José L. Oliver(2)

!     (1) Dpto. de Física Aplicada II, Universidad de Málaga, Spain
!         rick@ctima.uma.es, pcarpena@ctima.uma.es

!     (2) Dpto. de Genética, Universidad de Granada, Spain
!         oliver@ugr.es

!     Last updated code: IsoFinderC, april 2013
!     Tested on Debian Linux 7.0 and Windows 7
!     Improvements:
!				Output table file easily to be converted in BED format
!				Better identification of chromosome-contigs 
!				Minimum segment size increased to improve the statistical reliability of t-tests

!***************************************************************************

use leeadn
use findenes
use random
use cortes
use special
use diver

implicit none

!**** Secuencia de ADN y lista de cortes *******
character*90  dna_file,fout
integer*1     adn (350000000)  !Dna sequence
integer*4     nwords,enes,haygaps
integer*4     f_words(0:1)
type(corte)  ,pointer::first,last,p_max,p
type(corte)  ,target::pp
character*80  ss
character*10  method,method_des

INTEGER*4     MINSIZE  !minimo tamaño de los segmentos
INTEGER*4     WSIZE    !tamaño de ventana para el coarse-graining
REAL*8        SIGLEVEL !nivel de signifcación (tanto por 1)
INTEGER*4     NCUTS,NSEGS    !numero de cortes, numero de segmentos
INTEGER*4     IDUM

!******* PROVISIONALES ********
real*8        tmax
integer*4     p_corte,f_corte(0:1)
real*8        mu1,mu2,nv1,nv2,t_filt,sd,n1,n2,p0
integer*4     newcuts,ini,fin,size
real*8        p01,p02,p03,p04,p00,completo

!******* LISTA DE DATOS FINALES ********
integer*4   ,dimension(100000)::l_ini,l_fin,l_size
real*8      ,dimension(100000)::l_gc
integer*4    cont

!************** SALIDA ***************
character*80  name_seq,name_seq3,name_file
character*10  name_seq2
real          xgmc,cg_inf,cg_sup,gc_prev,gc_dif,ave,var
integer       i,j,l,m,jj,iii,jjj,seq_limpia,gmc
character*1   sf
character*10  chr
real*8        SCC

!********** PARAMETROS *********
call getarg(1,dna_file)

if (trim(dna_file)=='') then
  write(*,*)  
  write(*,*) 'Sintax:'
  write(*,*) '  isofinder <DNA FASTA file> <Sig.Level> <Sig.Method> <Coarse-graining level> <File out>'
  write(*,*)
  write(*,*) 'Sig.Method:'
  write(*,*)    ' r    Randomization'
  write(*,*)    ' p1   Parametric (t-student)  '
  write(*,*)    ' p2   Parametric (t-student) with different variances'
  write(*,*)    ' p3   Parametric for maximum value'
  write(*,*)
  write(*,*) 'Output (File out): chrom,chromStart,chromEnd,name_isochore,score,strand,length_isochore,GC_isochore'
  stop
endif

call getarg(2,ss)
read(ss,*) SIGLEVEL
call getarg(3,method)
  if (method=='r') then
     method_des='   Random.'
  else if (method=='p1') then
     method_des=' t-student'
  else if (method=='p2') then
     method_des=' t dif.var'
  else if (method=='p3') then
     method_des='   Maximum'
  endif

call getarg(4,ss)
read(ss,*) WSIZE
call getarg(5,fout)

IDUM=-34241117
!MINSIZE=2*WSIZE+2
MINSIZE=10*WSIZE+2 !For increasing statistical reliability of t-tests

fin=0
write(*,*)
write(*,*) '***********************   I S O F I N D E R   **************************'
write(*,*) '* Oliver JL, Carpena P, Hackenberg M, Bernaola-Galvan P. (2004)        *'
write(*,*) '* IsoFinder: computational prediction of isochores in genome sequences *'
write(*,*) '* Nucleic Acids Research 32: W287-W292                                 *'
write(*,*) '* http://dx.doi.org/10.1093/nar/gkh399                                 *'
write(*,*) '* Last updated code: IsoFinderC, april 2013                            *'
write(*,*) '* Tested on Debian Linux 7.0 and Windows 7                             *'
write(*,*) '************************************************************************'
write(*,*) 
write(*,*) 'Reading the sequence...'
call lee_adn(adn,nwords,enes,dna_file,1,fin,'SW',f_words)
write(*,*) '=================================================================='
write(*,*)             ' DNA file:               ',trim(dna_file)
write(*,'(a26,i10)')   ' Length:                 ',nwords
write(*,'(a26,f10.3)') ' G+C (%):                ',100.d0*f_words(1)/(sum(f_words))
write(*,'(a26,i10)')   ' Undefined               ',enes
write(*,'(a26,f10.5)') ' Sig.Level:              ',SIGLEVEL
write(*,'(a26,a10)')   ' Method:                 ',method_des
write(*,'(a26,i10)')   ' Coarse-graining level:  ',wsize
write(*,*) '=================================================================='

call init_list(first,last,nwords,f_words)
call insert_enes(adn,first,nwords,haygaps,IDUM)

if (haygaps==1) then
write(*,*) 'Locating islands of Ns...'

pp=first
do while (associated(pp%next))
  if (pp%intact==1) then
     write(*,'(3i10)') pp%posic+1,pp%next%posic,pp%next%posic-pp%posic
  endif
  pp=pp%next
enddo
write(*,*) '==============================='
endif

!**********************************************
!*************** SEGMENTACION *****************
!**********************************************

write(*,*) 'Segmenting...'

completo=0.d0
NCUTS=0
print*
write(*,'(a1,a12,f7.3,a2)',advance='no') char(13),'Segmented:  ',completo,' %'
do
  p=>first
  newcuts=0
  do while (associated(p%next))
    if (p%intact==0.and.(p%next%posic-p%posic)>2*MINSIZE) then
      call tstudentmax(adn,p,tmax,p_corte,f_corte)

      call filterout(adn,p%posic+1,p_corte,WSIZE,mu1,nv1,n1)
      call filterout(adn,p_corte+1,p%next%posic,WSIZE,mu2,nv2,n2)
      if ( (n1 <= 0) .or. (n2 <= 0) ) goto 500
	    sd=dsqrt((nv1+nv2)*(1/n1+1/n2)/(n1+n2-2))	 
      t_filt=abs(mu1-mu2)/sd

      !write(*,'(i10,4f10.5,2i6)') p%next%posic-p%posic,mu1,mu2,nv1,nv2,int(n1),int(n2)
      if (method=='r') then
         p0=random_ttest(adn,p,p_corte,WSIZE,t_filt,IDUM)          !Randomization
      else if (method=='p1') then
         p0=1.d0-ttest0(t_filt,1.d0*(n1+n2-2))                     !standard t-test
      else if (method=='p2') then
         p0=1.d0-ttest_difvar(mu1,mu2,nv1/(n1-1),nv2/(n2-1),n1,n2) !standard t-test diff. variances
      else if (method=='p3') then
         p0=prob(t_filt,n1+n2)                                     !maximum test
      endif
      if (p0>=SIGLEVEL) then
      !write(*,'(i10,4f10.5,2i6,f10.6,a4)') p%next%posic-p%posic, &
	    !&                 mu1,mu2,nv1,nv2,int(n1),int(n2),p0,' ***'

        call insert(p,p_corte,f_corte)
        p=>p%next
        newcuts=newcuts+1
        NCUTS=NCUTS+1
        !write(*,'(a1,a12,i4)',advance='no') char(13),'# of cuts:  ',NCUTS
      else
        p%intact=1
        completo=completo+100.d0*(p%next%posic-p%posic)/(1.d0*nwords)
        write(*,'(a1,a12,f7.3,a2)',advance='no') char(13),'Segmented:  ',completo,' %'

        !write(*,'(i10,4f10.5,2i6,f10.6)') p%next%posic-p%posic, &
	      !&              mu1,mu2,nv1,nv2,int(n1),int(n2),p0
      endif
    endif
    500 p=>p%next
  enddo
  if (newcuts==0) exit
enddo
write(*,'(a1,a12,f7.3,a2)',advance='no') char(13),'Segmented:  ',100.d0,' %'
write(*,*)
NSEGS=NCUTS+1
!**********************************************

!*********************************
!****** Crea  la lista   *********
!*********************************

p=>first
cont=0
do while (associated(p%next))
   cont=cont+1
   l_ini(cont)=p%posic+1
   l_fin(cont)=p%next%posic
   l_size(cont)=l_fin(cont)-l_ini(cont)+1
   l_gc(cont)=100.d0*(p%next%frec(1)-p%frec(1))/l_size(cont)
   p=>p%next
enddo

nsegs=cont   !hay que rectificar el número de segmentos
             !para tener en cuenta las enes.

!*****************  AHORA TIENES LOS INICIOS, FINALES, TAMANIOS  ******
!***************** Y G+C DE CADA SEGMENTO EN LOS ARRAYS          ******
!***************** l_ini,l_fin,l_size,l_gc                       ******
!***************** que van desde 1 a NSEGS                       ******

!Imprime isocoras
  name_seq ='                    '
  name_seq2='                    '
  name_seq3='                    '
  
  l=10
  do i = 90,1,-1
   if (dna_file(i:i) == '/' .or. dna_file(i:i) == '\') exit
  enddo
  do j = i-1,1,-1
   if (dna_file(j:j) == '/' .or. dna_file(j:j) == '\') exit
   chr(l:l) = dna_file(j:j)
   l=l-1
  enddo
  l=1
  do jj = j+1,60
   if (jj > 0) then
      if(dna_file(jj:jj) == '.') exit
      name_seq3(l:l) = dna_file(jj:jj)   
      l=l+1
   endif
  enddo
  name_seq3=trim(name_seq3)

  do i = 90,1,-1
   if (dna_file(i:i) == '.') exit
  enddo
  l=10
  do j = i-1,1,-1
   if ((dna_file(j:j) == '\') .or. (dna_file(j:j) == '/')) exit
        chr(l:l) = dna_file(j:j)
        l=l-1
  enddo
  l=1
  do jj = j+1,60
   if (dna_file(jj:jj) == '.') exit
   if (jj > 0) then
       name_seq(l:l) = dna_file(jj:jj)
       l=l+1
   endif
  enddo
  name_seq2 = trim(name_seq)

open(2,file=fout)
do cont=2,NSEGS-1             
   if ( (l_gc(cont) > 0) .and. (l_gc(cont-1) > 0) .and. (l_gc(cont+1) > 0) ) then !! NO SE IMPRIMEN el primer segmento ni el ultimo de cada contig ya que son isocoras incompletas
            xgmc=l_gc(cont)/100. !gc en tanto por 1
            cg_sup=0.65d0
            cg_inf=0.3d0
            gmc=int(1000.d0*(xgmc-cg_inf)/(cg_sup-cg_inf))
                !*** Esto vale 1000 para xgmc=0.65 (65%) y 0 para xgmc=0.3 (30%)
                !*** ahora lo arreglamos para cuando tengamos valores <30% o >65%
            if (gmc<0.d0) gmc=0.d0
            if (gmc>1000.d0) gmc=1000.d0
            !output: chrom,chromStart,chromEnd,name_isochore,score,strand,length_isochore,GC_isochore
            write(2,*)name_seq2,char(9),l_ini(cont)-1,char(9),l_fin(cont),char(9),name_seq2,'_',l_ini(cont), & !0-based coordinates
            char(9),gmc,char(9),'.',char(9),l_fin(cont)-l_ini(cont)+1,char(9),xgmc*100.                        !Formato BED8, abril, 2013
   endif
enddo

!*******************************************
!*******************************************
!*******************************************

CONTAINS

!*********** PROBABILITY STUFF *********************
!***************************************************

function neff(n)
  real*8  n,neff
  neff=-11.54+4.189*dlog(n)
end function neff

!***************************************************

FUNCTION PROB(x,n)

  real*8   prob,x,n
  real*8   f,expon,nu

  expon=neff(n)

  if (abs(x)<1.0e-5) then
    prob=0
  else
    expon=neff(n)
    f=0.8
    nu=n-2
    prob=(1.0-betai(0.5*f*nu,0.5d0*f,nu/(nu+x**2)))**expon
  endif

end FUNCTION prob

!***************************************************

!*************************************************

SUBROUTINE filterout(adn,ini,fin,wsize,mu,nsigma2,rnwin)
   integer*1   adn(:)
   integer*4   ini,fin,wsize
   real*8      mu,nsigma2,rnwin

   integer*4   i,j,nwin,size
   real*8      sx,sx2,x

   size=fin-ini+1
   nwin=size/wsize
   sx=0.d0
   sx2=0.d0
   do i=1,nwin
     x=0.d0
     do j=(i-1)*wsize,i*wsize-1
       if (adn(ini+j)==1) x=x+1
     enddo
     sx=sx+x
     sx2=sx2+x**2
   enddo
   sx=sx/wsize
   sx2=sx2/(wsize**2)

   mu=sx/(1.d0*nwin)
   nsigma2=(sx2-nwin*(mu**2))
   rnwin=nwin
end subroutine filterout

!*******************************************************

!*************************************************

function random_ttest(adn,p,p_corte,wsize,t0,idum)
   integer*1   adn(:)
   integer*4   ini,fin,wsize,nexper,p_corte
   real*8      mu1,mu2,nv1,nv2,gc(1000000)
   integer*4   n1,n2,n,idum
   type(corte) ,target::p

   integer*4   i,j,size,jj
   real*8      sd,t,t0,random_ttest,prob,x,df

   nexper=50000

   ini=p%posic+1
   fin=p_corte
   size=fin-ini+1
   n1=size/wsize
   do i=1,n1
     x=0.d0
     do j=(i-1)*wsize,i*wsize-1
       if (adn(ini+j)==1) x=x+1
     enddo
     gc(i)=x/(wsize)
   enddo

   ini=p_corte+1
   fin=p%next%posic
   size=fin-ini+1
   n2=size/wsize
   n=n1+n2
   do i=1,n2
     x=0.d0
     do j=(i-1)*wsize,i*wsize-1
       if (adn(ini+j)==1) x=x+1
     enddo
     gc(i+n1)=x/(wsize)
   enddo

   df=n-2.d0
   prob=0.d0
   do i=1,nexper
      do j=n,1,-1
        jj=(ran3(idum)*j)+1
        x=gc(jj)
        gc(jj)=gc(j)
        gc(j)=x
      enddo
      mu1=0.d0
      mu2=0.d0
      nv1=0.d0
      nv2=0.d0
      do j=1,n1
        mu1=mu1+gc(j)
        nv1=nv1+gc(j)**2
      enddo
      do j=n1+1,n
        mu2=mu2+gc(j)
        nv2=nv2+gc(j)**2
      enddo
      mu1=mu1/n1
      mu2=mu2/n2
      nv1=nv1-n1*mu1**2
      nv2=nv2-n2*mu2**2
      sd=dsqrt((nv1+nv2)*(1.d0/n1+1.d0/n2)/df)
      t=abs(mu1-mu2)/sd
      !print*, t,t0
      if (t<=t0) then
         prob=prob+1.d0
      endif
   enddo
   random_ttest=prob/(1.d0*nexper)

end function random_ttest


!******************************************************************

SUBROUTINE tstudentmax(adn,p,tmax,p_corte,f_corte)
  real*8       tmax
  integer*1    adn(:)
  integer*4    p_corte,f_corte(0:1)
  type(corte) ,target::p

  integer*1    b
  integer*4    ini,fin,i
  integer*4   ,dimension(0:1)::f1,f2
  real*8       n1,n2,n,sd,t,mu1,mu2


  ini=p%posic+1
  fin=p%next%posic
  f2=p%next%frec-p%frec
  f1=0
  n=fin-ini-1.d0            !*** resto 2 para la t
  tmax=0.d0
  do i=ini,fin-1
    b=adn(i)
    f1(b)=f1(b)+1
    f2(b)=f2(b)-1
    n1=i-ini+1
    n2=fin-i
    mu1=f1(0)/n1
    mu2=f2(0)/n2
	if ( (mu1*(1-mu1)+mu2*(1-mu2))*(1/n1+1/n2)/(n1+n2-2) <=0 ) goto 600
        sd=dsqrt((mu1*(1-mu1)+mu2*(1-mu2))*(1/n1+1/n2)/(n1+n2-2))
    t=abs(mu1-mu2)/sd
    if (t>tmax.and.(i-ini)>MINSIZE.and.(fin-i-1)>MINSIZE) then
       tmax=t
       p_corte=i
       f_corte=p%frec+f1
    endif
	600 continue
    enddo
END SUBROUTINE tstudentmax

FUNCTION total_divergence(first,f_words)
  !******** PARAMETERS **********
  type(corte),pointer:: first
  integer*4             f_words(:)
  real*8                total_divergence
  !****** LOCAL VARIABLES *******
  type(corte), pointer :: p
  integer*4               size,totsize
  real*8                  totentro,inter

  totentro=entro(f_words)
  totsize=sum(f_words)


  p=>first
  inter=0
  do while (associated(p%next))
    cont=cont+1
    size=sum(p%next%frec-p%frec)

    if (size/=0) then
      inter=inter-size*entro(p%next%frec-p%frec)
    endif
    p=>p%next
  end do


  total_divergence=inter/(1.0*totsize)+totentro
END FUNCTION total_divergence
!***************************************************

end
