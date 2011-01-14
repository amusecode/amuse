*
*
*
*         Program Monte-Carlo  MONTCARL
*         -----------------------------
*
*
*
      include 'common.h'
*
      integer iphase
*
*
*
*
*
*             read initial parameters
*
*
      print*,'calling input'
      call input
*     
      if(istart.eq.1) then
*
*
*             compute initial model and save it
*
*
         call start
*
         iphase = 1
*
      else
*
*
*             read restart model
*
*
         call mydump(2)
*            
*            read parameters which will be changed
*
         call input
*
         go to 40
*
      endif
*
*
 10   continue
*
      go to (20, 30) iphase
*
*
*             devision on zones and super-zones and output
*
*
 20   call output
 40   call zone
      iphase = 2
      go to 10
*
*
*            main flow rutine: calculation of cluster evolution
*
*
 30   call relaxt
      iphase = 1
      go to 10
*
      end
*
*
*
*
 
