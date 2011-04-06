FUNCTION call_aarseth_zare(TIME, BODY, x, y, z, vx, vy, vz, tend, nl)
    implicit none
    INTEGER call_aarseth_zare, AarsethZare
    REAL*8 TIME(3)
    REAL*8 BODY(3)
    REAL*8 x(3), y(3), z(3)
    REAL*8 POS(3,3)
    REAL*8 vx(3), vy(3), vz(3)
    REAL*8 VEL(3,3)
    REAL*8 tend(3)
    REAL*8 ETot, TCRIT
    INTEGER NREG, i, nl
    write(*,*) "BODY;", BODY
    TCRIT = tend(1)
    NREG = 0
    ETot = 0
    do i=1,3 ! first index is Ccartesian coordinate, second the particle number
       POS(1,i) = x(i)
       POS(2,i) = y(i)
       POS(3,i) = z(i)
       VEL(1,i) = vx(i)
       VEL(2,i) = vy(i)
       VEL(3,i) = vz(i)
    enddo
    call_aarseth_zare = AarsethZare(TIME,BODY,POS,VEL,TCRIT,ETot,NREG)
    write (*,*) 'call_aarseth_zare=', call_aarseth_zare, TIME, pos
    do i=1,3 ! first index is Ccartesian coordinate, second the particle number
       x(i)  = POS(1,i)
       y(i)  = POS(2,i)
       z(i)  = POS(3,i)
       vx(i) = VEL(1,i)
       vy(i) = VEL(2,i)
       vz(i) = VEL(3,i)
    enddo
    write (*,*) 'x=', x
END FUNCTION

FUNCTION construct_orbital_elements(m, r, v, e1, e2, nl)
    implicit none
    INTEGER construct_orbital_elements, transform1new 
    REAL*8 m(3), r(3), v(3), e1(3), e2(3)
    REAL*8 m1, m2, element(6)
    INTEGER nl
    m1 = m(1)
    m2 = m(2)
    construct_orbital_elements = transform1new(m1, m2, r, v, element)
    write (*,*) 'call_transform1new=', construct_orbital_elements, element
    e1(1) = element(1)
    e1(2) = element(2)
    e1(3) = element(3)
    e2(1) = element(4)
    e2(2) = element(5)
    e2(3) = element(6)
END FUNCTION

FUNCTION construct_orbital_coordinates(m, e1, e2, r, v, nl)
    implicit none
    INTEGER construct_orbital_coordinates, transform2
    REAL*8 m(3), r(3), v(3), e1(3), e2(3)
    REAL*8 mtot, element(6)
    INTEGER nl
    mtot = m(1)+m(2)
    element(1) = e1(1)
    element(2) = e1(2)
    element(3) = e1(3)
    element(4) = e2(1)
    element(5) = e2(2)
    element(6) = e2(3)
    construct_orbital_coordinates = transform2(mtot,element,r,v)
    write (*,*) 'call_transform2=', construct_orbital_coordinates, r, v
END FUNCTION

