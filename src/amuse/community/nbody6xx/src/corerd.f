***
      REAL*8 FUNCTION CORERD(KW,MC,M0,MFLASH)
*     
* A function to determine the radius of the core of a giant-like star.
* NOTE: this is out of date so rc should be obtained using HRDIAG!
* It is still OK to use but bear in mind that the core radius calculated
* for non-degenerate giant cores is only a rough estimate.
*
*     Author : C. A. Tout
*     Date :   26th February 1997
*     Updated 6/1/98 by J. Hurley
*
      INTEGER KW
      REAL*8 MC,MCH,M0,MFLASH
      PARAMETER (MCH = 1.44d0)
*
* First do the black holes and neutron stars.
*
      IF(KW.EQ.14)THEN
         CORERD = 4.24d-06*MC
      ELSEIF(KW.EQ.13)THEN
         CORERD = 1.4d-05
*
* Main sequence stars.
*
      ELSEIF(KW.LE.1.OR.KW.EQ.7)THEN
         CORERD = 0.d0
*
* Core-helium-burning stars, FAGB stars and non-degenerate giant cores.
*
      ELSEIF(KW.EQ.4.OR.KW.EQ.5.OR.(KW.LE.3.AND.M0.GT.MFLASH))THEN
         CORERD = 0.2239d0*MC**0.62d0
*
* The degenerate giants and white dwarfs.
*
      ELSE
         CORERD = 0.0115d0*SQRT(MAX(1.48204d-06,(MCH/MC)**(2.d0/3.d0)
     &                                        - (MC/MCH)**(2.d0/3.d0)))
* 
* Degenerate giants have hot subdwarf cores.
*
         IF(KW.LE.9) CORERD = 5.d0*CORERD
      ENDIF
*
      RETURN
      END
***
