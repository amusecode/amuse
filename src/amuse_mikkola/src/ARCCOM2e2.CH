c       ARCCOM.CH =filename                                                    
        IMPLICIT REAL*8 (A-H,M,O-Z)                                            
        PARAMETER (NMX=200,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &   NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
        COMMON/DataForRoutines1/X(NMX3),V(NMX3),WTTL,M(NMX),
     &   XC(NMX3),WC(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),N
        COMMON/DataForChainRoutinesTwo/MMIJ,CMX(3),CMV(3),
     &   ENERGY,Energr,CHTIME
        common/softening/ee,cmethod(3),Clight,NofBH
        common/TIMECOMMON/Taika,timecomparison
        common/spincommon/spin(3)! the relative spin of M(1) !Spin=spin*G*M^2/c
        common/tolerancecommon/EPS

