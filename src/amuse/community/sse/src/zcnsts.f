***
      SUBROUTINE zcnsts(z,zpars)
* 
      implicit none
      integer kw
*
      real*8 z,zpars(20)
      real*8 tm,tn,tscls(20),lums(10),GB(10)
      real*8 lzs,dlzs,lz,lzd,dum1,m1,m2,rr,rb,mhefl,lhefl,thefl,lx
      real*8 tbgbf,thef,lbagbf,lheif,lhef,lzahbf
      real*8 rgbf,ragbf,rminf,mcgbf
      external tbgbf,thef,lbagbf,lheif,lhef,lzahbf
      external rgbf,ragbf,rminf,mcgbf
*
      include 'zdata.h'
      real*8 msp(200),gbp(200),c(5)
      common /MSCFF/ msp
      common /GBCFF/ gbp
      data c /3.040581d-01, 8.049509d-02, 8.967485d-02,
     &        8.780198d-02, 2.219170d-02/
*
*       ------------------------------------------------------------
*
*      zpars:  1; M below which hook doesn't appear on MS, Mhook.
*              2; M above which He ignition occurs non-degenerately, Mhef.
*              3; M above which He ignition occurs on the HG, Mfgb. 
*              4; M below which C/O ignition doesn't occur, Mup.
*              5; M above which C ignites in the centre, Mec.
*              6; value of log D for M<= zpars(3)
*              7; value of x for Rgb propto M^(-x)
*              8; value of x for tMS = MAX(tHOOK,x*tBGB)
*              9; constant for McHeIf when computing Mc,BGB, mchefl.
*             10; constant for McHeIf when computing Mc,HeI, mchefl.
*             11; hydrogen abundance.
*             12; helium abundance.
*             13; constant x in rmin = rgb*x**y used by LM CHeB.
*             14; z**0.4 to be used for WD L formula.
*
*       ------------------------------------------------------------
*
      lzs = log10(z/0.02d0)
      dlzs = 1.d0/(z*log(10.d0))
      lz = log10(z)
      lzd = lzs + 1.d0
*
      zpars(1) = 1.0185d0 + lzs*(0.16015d0 + lzs*0.0892d0)
      zpars(2) = 1.995d0 + lzs*(0.25d0 + lzs*0.087d0)
      zpars(3) = 16.5d0*z**0.06d0/(1.d0 + (1.0d-04/z)**1.27d0)
      zpars(4) = MAX(6.11044d0 + 1.02167d0*lzs, 5.d0)
      zpars(5) = zpars(4) + 1.8d0
      zpars(6) = 5.37d0 + lzs*0.135d0
      zpars(7) = c(1) + lzs*(c(2) + lzs*(c(3) + lzs*(c(4) + lzs*c(5))))
      zpars(8) = MAX(0.95d0,MAX(0.95d0-(10.d0/3.d0)*(z-0.01d0),
     &           MIN(0.99d0,0.98d0-(100.d0/7.d0)*(z-0.001d0))))
***
* Lzams
      msp(1) = xz(1)+lzs*(xz(2)+lzs*(xz(3)+lzs*(xz(4)+lzs*xz(5))))
      msp(2) = xz(6)+lzs*(xz(7)+lzs*(xz(8)+lzs*(xz(9)+lzs*xz(10))))
      msp(3) = xz(11)+lzs*(xz(12)+lzs*(xz(13)+lzs*(xz(14)+lzs*xz(15))))
      msp(4) = xz(16)+lzs*(xz(17)+lzs*(xz(18)+lzs*(xz(19)+lzs*xz(20))))
      msp(5) = xz(21)+lzs*(xz(22)+lzs*(xz(23)+lzs*(xz(24)+lzs*xz(25))))
      msp(6) = xz(26)+lzs*(xz(27)+lzs*(xz(28)+lzs*(xz(29)+lzs*xz(30))))
      msp(7) = xz(31)+lzs*(xz(32)+lzs*(xz(33)+lzs*(xz(34)+lzs*xz(35))))
* Rzams
      msp(8) = xz(36)+lzs*(xz(37)+lzs*(xz(38)+lzs*(xz(39)+lzs*xz(40))))
      msp(9) = xz(41)+lzs*(xz(42)+lzs*(xz(43)+lzs*(xz(44)+lzs*xz(45))))
      msp(10) = xz(46)+lzs*(xz(47)+lzs*(xz(48)+lzs*(xz(49)+lzs*xz(50))))
      msp(11) = xz(51)+lzs*(xz(52)+lzs*(xz(53)+lzs*(xz(54)+lzs*xz(55))))
      msp(12) = xz(56)+lzs*(xz(57)+lzs*(xz(58)+lzs*(xz(59)+lzs*xz(60))))
      msp(13) = xz(61)
      msp(14) = xz(62)+lzs*(xz(63)+lzs*(xz(64)+lzs*(xz(65)+lzs*xz(66))))
      msp(15) = xz(67)+lzs*(xz(68)+lzs*(xz(69)+lzs*(xz(70)+lzs*xz(71))))
      msp(16) = xz(72)+lzs*(xz(73)+lzs*(xz(74)+lzs*(xz(75)+lzs*xz(76))))
* Tbgb 
      msp(17) = xt(1)+lzs*(xt(2)+lzs*(xt(3)+lzs*xt(4)))
      msp(18) = xt(5)+lzs*(xt(6)+lzs*(xt(7)+lzs*xt(8)))
      msp(19) = xt(9)+lzs*(xt(10)+lzs*(xt(11)+lzs*xt(12)))
      msp(20) = xt(13)+lzs*(xt(14)+lzs*(xt(15)+lzs*xt(16)))
      msp(21) = xt(17)
* dTbgb/dz
      msp(117) = dlzs*(xt(2)+lzs*(2.d0*xt(3)+3.d0*lzs*xt(4)))
      msp(118) = dlzs*(xt(6)+lzs*(2.d0*xt(7)+3.d0*lzs*xt(8)))
      msp(119) = dlzs*(xt(10)+lzs*(2.d0*xt(11)+3.d0*lzs*xt(12)))
      msp(120) = dlzs*(xt(14)+lzs*(2.d0*xt(15)+3.d0*lzs*xt(16)))
* Thook
      msp(22) = xt(18)+lzs*(xt(19)+lzs*(xt(20)+lzs*xt(21)))
      msp(23) = xt(22)
      msp(24) = xt(23)+lzs*(xt(24)+lzs*(xt(25)+lzs*xt(26)))
      msp(25) = xt(27)+lzs*(xt(28)+lzs*(xt(29)+lzs*xt(30)))
      msp(26) = xt(31)
* Ltms 
      msp(27) = xl(1)+lzs*(xl(2)+lzs*(xl(3)+lzs*(xl(4)+lzs*xl(5))))
      msp(28) = xl(6)+lzs*(xl(7)+lzs*(xl(8)+lzs*(xl(9)+lzs*xl(10))))
      msp(29) = xl(11)+lzs*(xl(12)+lzs*(xl(13)+lzs*xl(14)))
      msp(30) = xl(15)+lzs*(xl(16)+lzs*(xl(17)+lzs*(xl(18)+lzs*xl(19))))
      msp(27) = msp(27)*msp(30)
      msp(28) = msp(28)*msp(30)
      msp(31) = xl(20)+lzs*(xl(21)+lzs*(xl(22)+lzs*xl(23)))
      msp(32) = xl(24)+lzs*(xl(25)+lzs*(xl(26)+lzs*xl(27)))
* Lalpha
      m2 = 2.d0
      msp(33) = xl(28)+lzs*(xl(29)+lzs*(xl(30)+lzs*xl(31)))
      msp(34) = xl(32)+lzs*(xl(33)+lzs*(xl(34)+lzs*xl(35)))
      msp(35) = xl(36)+lzs*(xl(37)+lzs*(xl(38)+lzs*xl(39)))
      msp(36) = xl(40)+lzs*(xl(41)+lzs*(xl(42)+lzs*xl(43)))
      msp(37) = MAX(0.9d0,1.1064d0+lzs*(0.415d0+0.18d0*lzs))
      msp(38) = MAX(1.d0,1.19d0+lzs*(0.377d0+0.176d0*lzs))
      if(z.gt.0.01d0)then
         msp(37) = MIN(msp(37),1.d0)
         msp(38) = MIN(msp(38),1.1d0)
      endif
      msp(39) = MAX(0.145d0,0.0977d0-lzs*(0.231d0+0.0753d0*lzs))
      msp(40) = MIN(0.24d0+lzs*(0.18d0+0.595d0*lzs),0.306d0+0.053d0*lzs)
      msp(41) = MIN(0.33d0+lzs*(0.132d0+0.218d0*lzs),
     &              0.3625d0+0.062d0*lzs)
      msp(42) = (msp(33)+msp(34)*m2**msp(36))/
     &          (m2**0.4d0+msp(35)*m2**1.9d0)
* Lbeta
      msp(43) = xl(44)+lzs*(xl(45)+lzs*(xl(46)+lzs*(xl(47)+lzs*xl(48))))
      msp(44) = xl(49)+lzs*(xl(50)+lzs*(xl(51)+lzs*(xl(52)+lzs*xl(53))))
      msp(45) = xl(54)+lzs*(xl(55)+lzs*xl(56))
      msp(46) = MIN(1.4d0,1.5135d0+0.3769d0*lzs)
      msp(46) = MAX(0.6355d0-0.4192d0*lzs,MAX(1.25d0,msp(46)))
* Lhook
      msp(47) = xl(57)+lzs*(xl(58)+lzs*(xl(59)+lzs*xl(60)))
      msp(48) = xl(61)+lzs*(xl(62)+lzs*(xl(63)+lzs*xl(64)))
      msp(49) = xl(65)+lzs*(xl(66)+lzs*(xl(67)+lzs*xl(68)))
      msp(50) = xl(69)+lzs*(xl(70)+lzs*(xl(71)+lzs*xl(72)))
      msp(51) = MIN(1.4d0,1.5135d0+0.3769d0*lzs)
      msp(51) = MAX(0.6355d0-0.4192d0*lzs,MAX(1.25d0,msp(51)))
* Rtms
      msp(52) = xr(1)+lzs*(xr(2)+lzs*(xr(3)+lzs*(xr(4)+lzs*xr(5))))
      msp(53) = xr(6)+lzs*(xr(7)+lzs*(xr(8)+lzs*(xr(9)+lzs*xr(10))))
      msp(54) = xr(11)+lzs*(xr(12)+lzs*(xr(13)+lzs*(xr(14)+lzs*xr(15))))
      msp(55) = xr(16)+lzs*(xr(17)+lzs*(xr(18)+lzs*xr(19)))
      msp(56) = xr(20)+lzs*(xr(21)+lzs*(xr(22)+lzs*xr(23)))
      msp(52) = msp(52)*msp(54)
      msp(53) = msp(53)*msp(54)
      msp(57) = xr(24)
      msp(58) = xr(25)+lzs*(xr(26)+lzs*(xr(27)+lzs*xr(28)))
      msp(59) = xr(29)+lzs*(xr(30)+lzs*(xr(31)+lzs*xr(32)))
      msp(60) = xr(33)+lzs*(xr(34)+lzs*(xr(35)+lzs*xr(36)))
      msp(61) = xr(37)+lzs*(xr(38)+lzs*(xr(39)+lzs*xr(40)))
*
      msp(62) = MAX(0.097d0-0.1072d0*(lz+3.d0),MAX(0.097d0,MIN(0.1461d0,
     &              0.1461d0+0.1237d0*(lz+2.d0))))
      msp(62) = 10.d0**msp(62)
      m2 = msp(62) + 0.1d0
      msp(63) = (msp(52)+msp(53)*msp(62)**msp(55))/
     &          (msp(54)+msp(62)**msp(56))
      msp(64) = (msp(57)*m2**3+msp(58)*m2**msp(61)+
     &           msp(59)*m2**(msp(61)+1.5d0))/(msp(60)+m2**5)
* Ralpha
      msp(65) = xr(41)+lzs*(xr(42)+lzs*(xr(43)+lzs*xr(44)))
      msp(66) = xr(45)+lzs*(xr(46)+lzs*(xr(47)+lzs*xr(48)))
      msp(67) = xr(49)+lzs*(xr(50)+lzs*(xr(51)+lzs*xr(52)))
      msp(68) = xr(53)+lzs*(xr(54)+lzs*(xr(55)+lzs*xr(56)))
      msp(69) = xr(57)+lzs*(xr(58)+lzs*(xr(59)+lzs*(xr(60)+lzs*xr(61))))
      msp(70) = MAX(0.9d0,MIN(1.d0,1.116d0+0.166d0*lzs))
      msp(71) = MAX(1.477d0+0.296d0*lzs,MIN(1.6d0,-0.308d0-1.046d0*lzs))
      msp(71) = MAX(0.8d0,MIN(0.8d0-2.d0*lzs,msp(71)))
      msp(72) = xr(62)+lzs*(xr(63)+lzs*xr(64))
      msp(73) = MAX(0.065d0,0.0843d0-lzs*(0.0475d0+0.0352d0*lzs))
      msp(74) = 0.0736d0+lzs*(0.0749d0+0.04426d0*lzs)
      if(z.lt.0.004d0) msp(74) = MIN(0.055d0,msp(74))
      msp(75) = MAX(0.091d0,MIN(0.121d0,0.136d0+0.0352d0*lzs))
      msp(76) = (msp(65)*msp(71)**msp(67))/(msp(66) + msp(71)**msp(68))
      if(msp(70).gt.msp(71))then
         msp(70) = msp(71)
         msp(75) = msp(76)
      endif
* Rbeta
      msp(77) = xr(65)+lzs*(xr(66)+lzs*(xr(67)+lzs*xr(68)))
      msp(78) = xr(69)+lzs*(xr(70)+lzs*(xr(71)+lzs*xr(72)))
      msp(79) = xr(73)+lzs*(xr(74)+lzs*(xr(75)+lzs*xr(76)))
      msp(80) = xr(77)+lzs*(xr(78)+lzs*(xr(79)+lzs*xr(80)))
      msp(81) = xr(81)+lzs*(xr(82)+lzs*lzs*xr(83))
      if(z.gt.0.01d0) msp(81) = MAX(msp(81),0.95d0)
      msp(82) = MAX(1.4d0,MIN(1.6d0,1.6d0+lzs*(0.764d0+0.3322d0*lzs)))
* Rgamma
      msp(83) = MAX(xr(84)+lzs*(xr(85)+lzs*(xr(86)+lzs*xr(87))),
     &              xr(96)+lzs*(xr(97)+lzs*xr(98)))
      msp(84) = MIN(0.d0,xr(88)+lzs*(xr(89)+lzs*(xr(90)+lzs*xr(91))))
      msp(84) = MAX(msp(84),xr(99)+lzs*(xr(100)+lzs*xr(101)))
      msp(85) = xr(92)+lzs*(xr(93)+lzs*(xr(94)+lzs*xr(95)))
      msp(85) = MAX(0.d0,MIN(msp(85),7.454d0+9.046d0*lzs))
      msp(86) = MIN(xr(102)+lzs*xr(103),MAX(2.d0,-13.3d0-18.6d0*lzs))
      msp(87) = MIN(1.5d0,MAX(0.4d0,2.493d0+1.1475d0*lzs))
      msp(88) = MAX(1.d0,MIN(1.27d0,0.8109d0-0.6282d0*lzs))
      msp(88) = MAX(msp(88),0.6355d0-0.4192d0*lzs)
      msp(89) = MAX(5.855420d-02,-0.2711d0-lzs*(0.5756d0+0.0838d0*lzs))
* Rhook
      msp(90) = xr(104)+lzs*(xr(105)+lzs*(xr(106)+lzs*xr(107)))
      msp(91) = xr(108)+lzs*(xr(109)+lzs*(xr(110)+lzs*xr(111)))
      msp(92) = xr(112)+lzs*(xr(113)+lzs*(xr(114)+lzs*xr(115)))
      msp(93) = xr(116)+lzs*(xr(117)+lzs*(xr(118)+lzs*xr(119)))
      msp(94) = MIN(1.25d0,
     &          MAX(1.1d0,1.9848d0+lzs*(1.1386d0+0.3564d0*lzs)))
      msp(95) = 0.063d0 + lzs*(0.0481d0 + 0.00984d0*lzs)
      msp(96) = MIN(1.3d0,MAX(0.45d0,1.2d0+2.45d0*lzs))
* Lneta
      if(z.gt.0.0009d0)then
         msp(97) = 10.d0
      else
         msp(97) = 20.d0
      endif
* Lbgb 
      gbp(1) = xg(1)+lzs*(xg(2)+lzs*(xg(3)+lzs*xg(4)))
      gbp(2) = xg(5)+lzs*(xg(6)+lzs*(xg(7)+lzs*xg(8)))
      gbp(3) = xg(9)+lzs*(xg(10)+lzs*(xg(11)+lzs*xg(12)))
      gbp(4) = xg(13)+lzs*(xg(14)+lzs*(xg(15)+lzs*xg(16)))
      gbp(5) = xg(17)+lzs*(xg(18)+lzs*xg(19))
      gbp(6) = xg(20)+lzs*(xg(21)+lzs*xg(22))
      gbp(3) = gbp(3)**gbp(6)
      gbp(7) = xg(23)
      gbp(8) = xg(24)
* Lbagb
* set gbp(16) = 1.d0 until it is reset later with an initial
* call to Lbagbf using mass = zpars(2) and mhefl = 0.0
      gbp(9) = xg(25) + lzs*(xg(26) + lzs*xg(27))
      gbp(10) = xg(28) + lzs*(xg(29) + lzs*xg(30))
      gbp(11) = 15.d0
      gbp(12) = xg(31)+lzs*(xg(32)+lzs*(xg(33)+lzs*xg(34)))
      gbp(13) = xg(35)+lzs*(xg(36)+lzs*(xg(37)+lzs*xg(38)))
      gbp(14) = xg(39)+lzs*(xg(40)+lzs*(xg(41)+lzs*xg(42)))
      gbp(15) = xg(43)+lzs*xg(44)
      gbp(12) = gbp(12)**gbp(15)
      gbp(14) = gbp(14)**gbp(15)
      gbp(16) = 1.d0
* Rgb
      gbp(17) = -4.6739d0-0.9394d0*lz
      gbp(17) = 10.d0**gbp(17)
      gbp(17) = MAX(gbp(17),-0.04167d0+55.67d0*z)
      gbp(17) = MIN(gbp(17),0.4771d0-9329.21d0*z**2.94d0)
      gbp(18) = MIN(0.54d0,0.397d0+lzs*(0.28826d0+0.5293d0*lzs))
      gbp(19) = MAX(-0.1451d0,-2.2794d0-lz*(1.5175d0+0.254d0*lz))
      gbp(19) = 10.d0**gbp(19)
      if(z.gt.0.004d0)then
         gbp(19) = MAX(gbp(19),0.7307d0+14265.1d0*z**3.395d0)
      endif
      gbp(20) = xg(45)+lzs*(xg(46)+lzs*(xg(47)+lzs*(xg(48)+
     &          lzs*(xg(49)+lzs*xg(50)))))
      gbp(21) = xg(51)+lzs*(xg(52)+lzs*(xg(53)+lzs*(xg(54)+lzs*xg(55))))
      gbp(22) = xg(56)+lzs*(xg(57)+lzs*(xg(58)+lzs*(xg(59)+
     &          lzs*(xg(60)+lzs*xg(61)))))
      gbp(23) = xg(62)+lzs*(xg(63)+lzs*(xg(64)+lzs*(xg(65)+lzs*xg(66))))
* Ragb
      gbp(24) = MIN(0.99164d0-743.123d0*z**2.83d0,
     &              1.0422d0+lzs*(0.13156d0+0.045d0*lzs))
      gbp(25) = xg(67)+lzs*(xg(68)+lzs*(xg(69)+lzs*(xg(70)+
     &          lzs*(xg(71)+lzs*xg(72)))))
      gbp(26) = xg(73)+lzs*(xg(74)+lzs*(xg(75)+lzs*(xg(76)+lzs*xg(77))))
      gbp(27) = xg(78)+lzs*(xg(79)+lzs*(xg(80)+lzs*(xg(81)+
     &          lzs*(xg(82)+lzs*xg(83)))))
      gbp(28) = xg(84)+lzs*(xg(85)+lzs*(xg(86)+lzs*(xg(87)+lzs*xg(88))))
      gbp(29) = xg(89)+lzs*(xg(90)+lzs*(xg(91)+lzs*(xg(92)+
     &          lzs*(xg(93)+lzs*xg(94)))))
      gbp(30) = xg(95)+lzs*(xg(96)+lzs*(xg(97)+lzs*(xg(98)+
     &          lzs*(xg(99)+lzs*xg(100)))))
      m1 = zpars(2) - 0.2d0
      gbp(31) = gbp(29) + gbp(30)*m1
      gbp(32) = MIN(gbp(25)/zpars(2)**gbp(26),gbp(27)/zpars(2)**gbp(28))
* Mchei
      gbp(33) = xg(101)**4
      gbp(34) = xg(102)*4.d0
* Mcagb
      gbp(35) = xg(103)+lzs*(xg(104)+lzs*(xg(105)+lzs*xg(106)))
      gbp(36) = xg(107)+lzs*(xg(108)+lzs*(xg(109)+lzs*xg(110)))
      gbp(37) = xg(111)+lzs*xg(112)
      gbp(35) = gbp(35)**4
      gbp(36) = gbp(36)*4.d0
      gbp(37) = gbp(37)**4
* Lhei
* set gbp(41) = -1.d0 until it is reset later with an initial
* call to Lheif using mass = zpars(2) and mhefl = 0.0
      gbp(38) = xh(1)+lzs*xh(2)
      gbp(39) = xh(3)+lzs*xh(4)
      gbp(40) = xh(5)
      gbp(41) = -1.d0
      gbp(42) = xh(6)+lzs*(xh(7)+lzs*xh(8))
      gbp(43) = xh(9)+lzs*(xh(10)+lzs*xh(11))
      gbp(44) = xh(12)+lzs*(xh(13)+lzs*xh(14))
      gbp(42) = gbp(42)**2
      gbp(44) = gbp(44)**2
* Lhe
      gbp(45) = xh(15)+lzs*(xh(16)+lzs*xh(17))
      if(lzs.gt.-1.d0)then
         gbp(46) = 1.d0 - xh(19)*(lzs+1.d0)**xh(18)
      else
         gbp(46) = 1.d0
      endif
      gbp(47) = xh(20)+lzs*(xh(21)+lzs*xh(22))
      gbp(48) = xh(23)+lzs*(xh(24)+lzs*xh(25))
      gbp(45) = gbp(45)**gbp(48)
      gbp(47) = gbp(47)**gbp(48)
      gbp(46) = gbp(46)/zpars(3)**0.1d0+(gbp(46)*gbp(47)-gbp(45))/
     &          zpars(3)**(gbp(48)+0.1d0)
* Rmin
      gbp(49) = xh(26)+lzs*(xh(27)+lzs*(xh(28)+lzs*xh(29)))
      gbp(50) = xh(30)+lzs*(xh(31)+lzs*(xh(32)+lzs*xh(33)))
      gbp(51) = xh(34)+lzs*(xh(35)+lzs*(xh(36)+lzs*xh(37)))
      gbp(52) = 5.d0+xh(38)*z**xh(39)
      gbp(53) = xh(40)+lzs*(xh(41)+lzs*(xh(42)+lzs*xh(43)))
      gbp(49) = gbp(49)**gbp(53)
      gbp(51) = gbp(51)**(2.d0*gbp(53))
* The
* set gbp(57) = -1.d0 until it is reset later with an initial
* call to Thef using mass = zpars(2), mc = 0.0  and mhefl = 0.0
      gbp(54) = xh(44)+lzs*(xh(45)+lzs*(xh(46)+lzs*xh(47)))
      gbp(55) = xh(48)+lzs*(xh(49)+lzs*xh(50))
      gbp(55) = MAX(gbp(55),1.d0)
      gbp(56) = xh(51)
      gbp(57) = -1.d0
      gbp(58) = xh(52)+lzs*(xh(53)+lzs*(xh(54)+lzs*xh(55)))
      gbp(59) = xh(56)+lzs*(xh(57)+lzs*(xh(58)+lzs*xh(59)))
      gbp(60) = xh(60)+lzs*(xh(61)+lzs*(xh(62)+lzs*xh(63)))
      gbp(61) = xh(64)+lzs*xh(65)
      gbp(58) = gbp(58)**gbp(61)
      gbp(60) = gbp(60)**5
* Tbl
      dum1 = zpars(2)/zpars(3)
      gbp(62) = xh(66)+lzs*xh(67)
      gbp(62) = -gbp(62)*log10(dum1)
      gbp(63) = xh(68)
      if(lzd.gt.0.d0) then
         gbp(64) = 1.d0-lzd*(xh(69)+lzd*(xh(70)+lzd*xh(71)))
      else
         gbp(64) = 1.d0
      end if
      gbp(65) = 1.d0-gbp(64)*dum1**gbp(63)
      gbp(66) = 1.d0 - lzd*(xh(77) + lzd*(xh(78) + lzd*xh(79)))
      gbp(67) = xh(72) + lzs*(xh(73) + lzs*(xh(74) + lzs*xh(75)))
      gbp(68) = xh(76)
* Lzahb
      gbp(69) = xh(80) + lzs*(xh(81) + lzs*xh(82))
      gbp(70) = xh(83) + lzs*(xh(84) + lzs*xh(85))
      gbp(71) = 15.d0
      gbp(72) = xh(86)
      gbp(73) = xh(87)
* Rzahb
      gbp(75) = xh(88) + lzs*(xh(89) + lzs*(xh(90) + lzs*xh(91)))
      gbp(76) = xh(92) + lzs*(xh(93) + lzs*(xh(94) + lzs*xh(95)))
      gbp(77) = xh(96) + lzs*(xh(97) + lzs*(xh(98) + lzs*xh(99)))
***
* finish Lbagb
      mhefl = 0.d0
      lx = lbagbf(zpars(2),mhefl)
      gbp(16) = lx
* finish LHeI
      dum1 = 0.d0
      lhefl = lheif(zpars(2),mhefl)
      gbp(41) = (gbp(38)*zpars(2)**gbp(39)-lhefl)/
     &          (EXP(zpars(2)*gbp(40))*lhefl)
* finish THe
      thefl = thef(zpars(2),dum1,mhefl)*tbgbf(zpars(2))
      gbp(57) = (thefl-gbp(54))/(gbp(54)*EXP(gbp(56)*zpars(2)))
* finish Tblf
      rb = ragbf(zpars(3),lheif(zpars(3),zpars(2)),mhefl)
      rr = 1.d0 - rminf(zpars(3))/rb
      rr = MAX(rr,1.0d-12)
      gbp(66) = gbp(66)/(zpars(3)**gbp(67)*rr**gbp(68))
* finish Lzahb
      gbp(74) = lhefl*lHef(zpars(2))
***
      kw = 0
      tm = 0.d0
      tn = 0.d0
      CALL star(kw,zpars(2),zpars(2),tm,tn,tscls,lums,GB,zpars)
      zpars(9) = mcgbf(lums(3),GB,lums(6))
      zpars(10) = mcgbf(lums(4),GB,lums(6))
* set the hydrogen and helium abundances
      zpars(11) = 0.76d0 - 3.d0*z
      zpars(12) = 0.24d0 + 2.d0*z
* set constant for low-mass CHeB stars
      zpars(13) = rminf(zpars(2))/
     &            rgbf(zpars(2),lzahbf(zpars(2),zpars(9),zpars(2)))
* 
      zpars(14) = z**0.4d0
*
      return
      end
***
