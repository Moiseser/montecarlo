program timetest
  implicit none
  double precision :: totaltime,START,END
  CALL CPU_TIME(START)  

  CALL CPU_TIME(END) 

  totaltime = START - END
  write (*,*) totaltime 
  end program 
