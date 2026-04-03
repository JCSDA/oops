! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Linked list implementation

!> Linked list subroutines

subroutine num_elements(list,num)
type(registry_t), intent(in) :: list
integer,intent(out) :: num
type(node_t), pointer :: next

num=0

next => list % head
do while(associated(next))
  num = num + 1
  next => next%next
enddo

end subroutine num_elements

integer function list_size(list) result(num)
type(registry_t), intent(in) :: list
type(node_t), pointer :: next
num=0
next => list % head
do while(associated(next))
  num = num + 1
  next => next%next
enddo
end function list_size

subroutine max_key(list,num)
type(registry_t), intent(in) :: list
integer,intent(out) :: num
type(node_t), pointer :: next

num=0

next => list % head
do while(associated(next))
  if (next % key > num) then
    num = next % key
  end if
  next => next%next
enddo

num = num + 1

end subroutine max_key

!> Initialize the linked list
subroutine init_(self)
 class(registry_t), intent(inout) :: self
end subroutine

!> Add element to the linked list
subroutine add_(self,key)
 class(registry_t), intent(inout) :: self
 integer, intent(out)             :: key

 type(node_t), pointer :: next

 call max_key(self, key)

 !allocate next element and assign key
 allocate(next)
 next%key = key

 !move the head to the front of the list
 if (associated(self % head)) then
   next%next => self%head%next
   self%head%next => next
 else
   self%head => next
 endif
end subroutine

!> Fetch element of the linked list by key
subroutine get_(self,key,ptr)
 class(registry_t), intent(in) :: self
 integer, intent(in)           :: key
 type (LISTED_TYPE), pointer   :: ptr

 type(node_t), pointer :: next
 external abor1_ftn

 !note that the list starts from self%head%next
 next => self%head
 ptr => NULL()

 !sweep the linked list to find matching key
 do while(associated(next))
   if(key.eq.next%key) then
    ptr => next%element
    exit
   else
     next => next % next
   endif
 end do
 if (.not.associated(ptr)) call abor1_ftn("registry_t%get_: key not found")
end subroutine

!> Remove element of the linked list
subroutine remove_(self,key)
 class(registry_t), intent(inout) :: self
 integer, intent(inout)           :: key

 type(node_t), pointer :: prev
 type(node_t), pointer :: next

 next => self%head
 nullify(prev)
 
 !sweep the linked list to find matching key, 
 do while(associated(next))
  if(key.eq.next%key) then
    exit
  endif
  prev => next
  next => prev%next
 enddo
 
 !reconnect the list
 if(associated(prev)) then
  if(associated(next % next)) then
    prev%next => next%next
    deallocate(next)
  else
    nullify(prev%next)
    deallocate(next)
  endif
 else
  if(associated(next % next)) then
    self%head => next%next
    deallocate(next)
  else
    nullify(self%head)
    deallocate(next)
  endif
 endif
 !remove the node and set key to 0
 key=0
end subroutine

!> linkedlist generic setup
subroutine registry_setup_(self, c_key_self, ptr)
  class(registry_t), intent(inout) :: self
  integer, intent(inout) :: c_key_self
  type (LISTED_TYPE), pointer :: ptr

  call self%init()
  call self%add(c_key_self)
  call self%get(c_key_self, ptr)
end subroutine

!> linkedlist generic delete
subroutine registry_delete_(self, c_key_self, ptr)
  class(registry_t), intent(inout) :: self
  integer, intent(inout) :: c_key_self
  type (LISTED_TYPE), pointer :: ptr


  call self%get(c_key_self, ptr)
  call self%remove(c_key_self)
end subroutine
