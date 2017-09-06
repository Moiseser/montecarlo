program test
  double precision :: T,x
   
T=.1
do i=1,5
   write (*,*)
  x = 10*T
write (*,*) 'x' , x
T=T+.1
write (*,*) 'i:', i
end do

end program test 
