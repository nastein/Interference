program test_iso

   use dirac_matrices
   implicit none

   integer*4 :: i,j,k,l
   real*8 :: test_p4(4),w 
   complex*16 :: t1(2), t1p(2), t2(2), t2p(2), p(2), n(2), result, result2
   complex*16 :: czero = (0.0d0,0.0d0)
   complex*16 :: cone  = (1.0d0,0.0d0)
   complex*16 :: ci    = (0.0d0,1.0d0)

   w = 0.0d0
   test_p4 = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)


   call dirac_matrices_in(938.0d0,938.0d0,100.0d0)
   call current_init(w,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,2)
   call define_spinors()

   p(1) = cone
   p(2) = czero

   n(1) = czero
   n(2) = cone

   t1 = n
   t2 = p
   t1p = p 
   t2p = p

   !Ivminus is called Ivminus(it1,it2,itp1,itp2)

   result = -Ivminus(p,p,n,p)
   write(6,*) '<np | Iv^dag | pp> = ', result

   result = -Ivminus(n,p,n,n)
   write(6,*) '<nn | Iv^dag | np> = ', result

   Write(6,*) '12 -> 21 exchange:'
   result2 = IDeltaB(p,p,n,p)
   write(6,*) '<pp | DeltaB | np> = ', result2

   result2 = IDeltaB(n,p,n,n)
   write(6,*) '<np | DeltaB | nn> = ', result2

   Write(6,*) '21 -> 12 exchange:'
   result2 = IDeltaB(p,p,p,n)
   write(6,*) '<pp | DeltaB | pn> = ', result2

   result2 = IDeltaB(p,n,n,n)
   write(6,*) '<pn | DeltaB | nn> = ', result2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Write(6,*) '12 -> 21 exchange:'
   result2 = IDeltaC(p,p,n,p)
   write(6,*) '<pp | DeltaC | np> = ', result2

   result2 = IDeltaC(n,p,n,n)
   write(6,*) '<np | DeltaC | nn> = ', result2

   Write(6,*) '21 -> 12 exchange:'
   result2 = IDeltaC(p,p,p,n)
   write(6,*) '<pp | DeltaC | pn> = ', result2

   result2 = IDeltaC(p,n,n,n)
   write(6,*) '<pn | DeltaC | nn> = ', result2
end program test_iso


