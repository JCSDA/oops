! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Linked list implementation

!> Linked list subroutines

!> Initialize the linked list
subroutine init_(self)
 class(registry_t), intent(inout) :: self

 !set count to zero and allocate the head of the list
 if(.not.self%l_init.or..not.associated(self%head)) then
  self%count = 0
  allocate(self%head)
  nullify(self%head%next)
  self%l_init=.true.
 endif
end subroutine

!> Add element to the linked list
subroutine add_(self,key)
 class(registry_t), intent(inout) :: self
 integer, intent(inout)           :: key

 type(node_t), pointer :: next

 !increase global counter and assign key
 self%count = self%count+1
 key = self%count

 !allocate next element and assign key
 allocate(next)
 next%key = key

 !move the head to the front of the list
 next%next => self%head%next
 self%head%next => next
end subroutine

!> Fetch element of the linked list by key
subroutine get_(self,key,ptr)
 class(registry_t), intent(in) :: self
 integer, intent(in)           :: key
 type (LISTED_TYPE), pointer   :: ptr

 type(node_t), pointer :: next

 !note that the list starts from self%head%next
 next => self%head
 ptr => NULL()

 !sweep the linked list to find matching key
 do while(associated(next))
  next=>next%next
  if(key.eq.next%key) then
   ptr => next%element
   exit
  endif
 enddo
 if (.not.associated(ptr)) call abor1_ftn("registry_t%get_: key not found")
end subroutine

!> Remove element of the linked list
subroutine remove_(self,key)
 class(registry_t), intent(inout) :: self
 integer, intent(inout)           :: key

 type(node_t), pointer :: prev
 type(node_t), pointer :: next

 next => self%head%next
 prev => NULL()

 !sweep the linked list to find matching key, 
 do while(associated(next))
  if(key.eq.next%key) exit
  prev => next
  next => next%next
 enddo
 !reconnect the list
 if(associated(next%next)) then
  if(associated(prev)) then
   prev%next => next%next
  else
   self%head%next=>next%next
  endif
 endif
 !remove the node and set key to 0
 if(associated(next)) deallocate(next)
 key=0
 return
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
