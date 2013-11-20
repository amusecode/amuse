!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert1d(x,order,npoint)
  use amr_parameters, ONLY: qdp
  implicit none
  integer     ,INTENT(IN)                     ::npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x
  real(qdp),INTENT(OUT),dimension(1:npoint)::order

  integer::ip

  do ip=1,npoint
#ifdef QUADHILBERT
     order(ip)=real(x(ip),kind=16)
#else
     order(ip)=real(x(ip),kind=8)
#endif
  end do

end subroutine hilbert1d
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert2d(x,y,order,bit_length,npoint)

  use amr_parameters, ONLY: qdp
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y
  real(qdp),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:2*bit_length-1)::i_bit_mask 
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask
  integer,dimension(0:3,0:1,0:3)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert2d'
     call clean_stop
  endif

  state_diagram = RESHAPE( (/ 1, 0, 2, 0, &
                            & 0, 1, 3, 2, &
                            & 0, 3, 1, 1, &
                            & 0, 3, 1, 2, &
                            & 2, 2, 0, 3, &
                            & 2, 1, 3, 0, &
                            & 3, 1, 3, 2, &
                            & 2, 3, 1, 0  /), &
                         & (/ 4, 2, 4 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
     enddo
     
     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(2*i+1)=x_bit_mask(i)
        i_bit_mask(2*i  )=y_bit_mask(i)
     end do
     
     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b1=0 ; if(i_bit_mask(2*i+1))b1=1
        b0=0 ; if(i_bit_mask(2*i)  )b0=1
        sdigit=b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(2*i+1)=btest(hdigit,1)
        i_bit_mask(2*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,2*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
#ifdef QUADHILBERT
        order(ip)=order(ip)+real(b0,kind=16)*real(2,kind=16)**i
#else
        order(ip)=order(ip)+real(b0,kind=8)*real(2,kind=8)**i
#endif
     end do

  end do

end subroutine hilbert2d
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  use amr_parameters, ONLY: qdp
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(qdp),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     call clean_stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo
     
     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do
     
     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
#ifdef QUADHILBERT
        order(ip)=order(ip)+real(b0,kind=16)*real(2,kind=16)**i
#else
        order(ip)=order(ip)+real(b0,kind=8)*real(2,kind=8)**i
#endif
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================

