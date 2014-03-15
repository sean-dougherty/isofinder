module cortes




!A linked list to store all cuts
type corte
  integer*4 posic
  integer*4 frec(0:1)
  integer*1 intact
  type(corte), pointer :: next,prev
end type


CONTAINS


!***************************************************

subroutine init_list(first,last,&
           nwords,f_words)
  type(corte), pointer::first,last
  integer*4 nwords
  integer*4, dimension(:):: f_words

  allocate(first)
  allocate(last)

  first%posic=0
  first%frec=0
  first%next=>last
  first%intact=0
  nullify(first%prev)

  last%posic=nwords
  last%frec=f_words
  last%intact=0
  nullify(last%next)
  last%prev=>first
end subroutine init_list

!***********************************************

subroutine insert(p,posic,frec)
  type(corte), pointer :: p,newp
  integer*4               posic,frec(:)


  allocate(newp)

  newp%posic=posic
  newp%frec=frec
  newp%next=>p%next
  newp%prev=>p
  newp%intact=0

  p%next%prev=>newp
  p%next=>newp


END SUBROUTINE insert

!**************************************************


SUBROUTINE restore_list(first)
  type(corte),pointer :: first,p

  p=>first
  do while (associated(p%next))
    p%intact=0
    p=>p%next
  enddo
end subroutine restore_list

!***********************************************
end MODULE cortes