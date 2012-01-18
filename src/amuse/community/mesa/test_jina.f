! ***********************************************************************
!
!   Copyright (C) 2010  Ed Brown, Bill Paxton
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful, 
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, writeto the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
 
      module test_jina_support
      
      use chem_def, only: nuclide_set, free_nuclide_data
      use chem_lib
   	use jina_def
   	use jina_lib
   	use alert_lib
   	use const_lib
   	
   	implicit none

   	
      
      contains
      

      subroutine set_up_430(Z, A)
      	integer, parameter :: N = 430
      	integer, dimension(N) :: Z, A
      	integer :: i
      	!n
      	Z(1) = 0; A(1) = 1
      	!H
      	Z(2:4) = 1; A(2:4) = [1, 2, 3]
      	!He
      	Z(5:6) = 2; A(5:6) = [3, 4]
      	!Li
      	Z(7:9) = 3; A(7:9) = [6, 7, 8]
      	!Be
      	Z(10:13) = 4; A(10:13) = [7, 9, 10, 11]
      	!B
      	Z(14:19) = 5; A(14:19) = [8, 10, 11, 12, 13, 14]
      	!C
      	Z(20:27) = 6; A(20:27) = [(i, i=9, 16)]
      	!N
      	Z(28:36) = 7; A(28:36) = [(i, i=12, 20)]
      	!O
      	Z(37:44) = 8; A(37:44) = [(i, i=13, 20)]
      	!F
      	Z(45:54) = 9; A(45:54) = [(i, i=15, 24)]
      	!Ne
      	Z(55:66) = 10; A(55:66) = [(i, i=17, 28)]
      	!Na
      	Z(67:78) = 11; A(67:78) = [(i, i=20, 31)]
      	!Mg
      	Z(79:92) = 12; A(79:92) = [(i, i=20, 33)]
      	!Al
      	Z(93:106) = 13; A(93:106) = [(i, i=22, 35)]
      	!Si
      	Z(107:123) = 14; A(107:123) = [(i, i=22, 38)]
      	!P
      	Z(124:138) = 15; A(124:138) = [(i, i=26, 40)]
      	!S
      	Z(139:154) = 16; A(139:154) = [(i, i=27, 42)]
      	!Cl
      	Z(155:168) = 17; A(155:168) = [(i, i=31, 44)]
      	!Ar
      	Z(169:185) = 18; A(169:185) = [(i, i=31, 47)]
      	!K
      	Z(186:197) = 19; A(186:197) = [(i, i=35, 46)]
      	!Ca
      	Z(198:216) = 20; A(198:216) = [(i, i=35, 53)]
      	!Sc
      	Z(217:230) = 21; A(217:230) = [(i, i=40, 53)]
      	!Ti
      	Z(231:247) = 22; A(231:247) = [(i, i=39, 55)]
      	!V
      	Z(248:262) = 23; A(248:262) = [(i, i=43, 57)]
      	!Cr
      	Z(263:280) = 24; A(263:280) = [(i, i=43, 60)]
      	!Mn
      	Z(281:298) = 25; A(281:298) = [(i, i=46, 63)]
      	!Fe
      	Z(299:319) = 26; A(299:319) = [(i, i=46, 66)]
      	!Co
      	Z(320:337) = 27; A(320:337) = [(i, i=50, 67)]
      	!Ni
      	Z(338:361) = 28; A(338:361) = [(i, i=50, 73)]
      	!Cu
      	Z(362:378) = 29; A(362:378) = [(i, i=56, 72)]
      	!Zn
      	Z(379:396) = 30; A(379:396) = [(i, i=55, 72)]
      	!Ga
      	Z(397:412) = 31; A(397:412) = [(i, i=60, 75)]
      	!Ge
      	Z(413:430) = 32; A(413:430) = [(i, i=59, 76)]
      end subroutine set_up_430


      subroutine init_430(set)
      	integer, parameter :: N = 430
      	integer, dimension(N) :: Z, A
      	character(len=iso_name_length) :: names(N), npa(3)
      	type(nuclide_set), dimension(N), intent(out) :: set
      	integer :: i      	
      	call set_up_430(Z, A)	
      	call generate_nuclide_names(Z,A,names)	
      	! move n, p, he4 to the end
      	npa = names([1,2,6])
      	! pull off n,p
      	names = cshift(names,2)
      	names(4:427) = names(5:428)
      	names(428:430) = npa
      	call generate_nuclide_set(names,set)	
      end subroutine init_430


      subroutine init_25(set)
         use chem_lib, only: chem_create_set_from_file
      	type(nuclide_set), dimension(:), pointer :: set
      	integer :: ierr
      	call chem_create_set_from_file('test_ZA.txt', '', set, ierr)
      	if (ierr /= 0) then
      	   write(*,*) 'failed in chem_create_set_from_file'
      	   stop 1
      	end if
      end subroutine init_25


      subroutine init_100(set)
         use chem_lib, only: chem_create_set_from_file
      	type(nuclide_set), dimension(:), pointer :: set
      	integer :: ierr
      	call chem_create_set_from_file('net100.jina', '', set, ierr)
      	if (ierr /= 0) then
      	   write(*,*) 'failed in chem_create_set_from_file'
      	   stop 1
      	end if
      end subroutine init_100


      subroutine init_200(set)
         use chem_lib, only: chem_create_set_from_file
      	type(nuclide_set), dimension(:), pointer :: set
      	integer :: ierr
      	call chem_create_set_from_file('net200.data', '', set, ierr)
      	if (ierr /= 0) then
      	   write(*,*) 'failed in chem_create_set_from_file'
      	   stop 1
      	end if
      end subroutine init_200


      subroutine init_Limongi_et_al_2000(set)
      	type(nuclide_set), dimension(:), pointer :: set
      	integer :: ierr
      	call chem_create_set_from_file('net_Limongi_et_al_2000.data', '', set, ierr)
      	if (ierr /= 0) then
      	   write(*,*) 'failed in chem_create_set_from_file'
      	   stop 1
      	end if
      end subroutine init_Limongi_et_al_2000
      
      

      subroutine set_up_690(Z, A)
      	integer, parameter :: N = 690
      	integer, dimension(N) :: Z, A
      	integer :: ilo, ihi
      	ilo=0; ihi=0
      	write(*,*) 'set_up_690'
      	!n
      	call setZA(0,1,1)
      	!H
      	call setZA(1,1,3)
      	!He
      	call setZA(2,3,4)
      	!Li
      	call setZA(3,6,7)
      	!Be
      	call setZA(4,7,7)
      	call setZA(4,9,10)
      	!B
      	call setZA(5,8,8)
      	call setZA(5,10,11)
      	!C
      	call setZA(6,11,14)
      	!N
      	call setZA(7,13,15)
      	!O
      	call setZA(8,14,18)
      	!F
      	call setZA(9,17,19)
      	!Ne
      	call setZA(10,18,23)
      	!Na
      	call setZA(11,20,24)
      	!Mg
      	call setZA(12,21,27)
      	!Al
      	call setZA(13,23,28)
      	!Si
      	call setZA(14,24,32)
      	!P
      	call setZA(15,27,34)
      	!S
      	call setZA(16,29,37)
      	!Cl
      	call setZA(17,31,38)
      	!Ar
      	call setZA(18,33,41)
      	!K
      	call setZA(19,36,42)
      	!Ca
      	call setZA(20,37,49)
      	!Sc
      	call setZA(21,40,50)
      	!Ti
      	call setZA(22,41,51)
      	!V
      	call setZA(23,44,52)
      	!Cr
      	call setZA(24,45,55)
      	!Mn
      	call setZA(25,48,57)
      	!Fe
      	call setZA(26,49,61)
      	!Co
      	call setZA(27,50,62)
      	!Ni
      	call setZA(28,53,66)
      	!Cu
      	call setZA(29,56,66)
      	!Zn
      	call setZA(30,55,70)
      	!Ga
      	call setZA(31,60,72)
      	!Ge
      	call setZA(32,59,76)
      	!arsenic
      	call setZA(33,64,77)
      	!selenium 
      	call setZA(34,65,82)
      	!bromine 
      	call setZA(35,70,83)
         ! krypton 
         call setZA(36,68,86)
         ! rubidium 
         call setZA(37,74,89)
         ! strontium 
         call setZA(38,72,91)
         ! yttrium 
         call setZA(39,77,94)
         ! zirconium 
         call setZA(40,76,95)
         ! niobium 
         call setZA(41,82,97)
         ! molybdenum 
         call setZA(42,80,98)
         ! technetium 
         call setZA(43,86,99)
         ! ruthenium 
         call setZA(44,86,101)
         ! rhodium 
         call setZA(45,89,103)
         ! palladium 
         call setZA(46,88,105)
         ! silver 
         call setZA(47,93,109)
         ! cadmium 
         call setZA(48,92,110)
         ! indium 
         call setZA(49,98,117)
         ! tin 
         call setZA(50,96,124)
         ! antimony 
         call setZA(51,104,125)
         ! tellurium 
         call setZA(52,104,130)
         ! iodine 
         call setZA(53,108,131)
         ! xenon 
         call setZA(54,108,136)

      	if (ihi /= N) then
         	write(*,*) 'set_up_700: ihi /= N', ihi, N
         	stop
      	end if
      	
      	contains
      
         subroutine setZA(zee,alo,ahi)
            integer, intent(in) :: zee, alo, ahi
            integer :: i
         	ilo = ihi+1; ihi=ilo+ahi-alo
         	if (ihi > N) return
         	Z(ilo:ihi) = zee
         	A(ilo:ihi) = [(i, i=alo,ahi)]
         	!write(*,'(a,i5,5x,3i5,5x,99i5)') 'ihi, zee, alo, ahi, A', ihi, zee, alo, ahi, A(ilo:ihi)
         end subroutine setZA
         
      end subroutine set_up_690


      subroutine init_690(set)
      	integer, dimension(690) :: Z, A
      	character(len=iso_name_length) :: names(690), npa(3)
      	type(nuclide_set), dimension(690), intent(out) :: set
      	integer :: i      	
      	call set_up_690(Z, A)	
      	call generate_nuclide_names(Z,A,names)	
      	! move n, p, he4 to the end
      	npa = names([1,2,6])
      	! pull off n,p
      	names = cshift(names,2)
      	names(4:687) = names(5:688)
      	names(688:690) = npa
      	call generate_nuclide_set(names,set)	
      end subroutine init_690          
      

      subroutine set_up_cno(Z, A) ! set up a simple CNO network
      	integer, parameter :: N = 10
      	integer, dimension(N) :: Z, A
      	integer :: i
      	!n
      	Z(1) = 0; A(1) = 1
      	!H
      	Z(2:3) = 1; A(2:3) = [1,2]
      	!He
      	Z(4) = 2; A(4) = 4
      	! C
      	Z(5:6) = 6; A(5:6) = [12,13]
      	!N
      	Z(7:9) = 7; A(7:9) = [13,14,15]
      	!O
      	Z(10) = 8; A(10) = 15
      end subroutine set_up_cno


      subroutine init_cno(set)
      	integer, dimension(10) :: Z, A
      	character(len=iso_name_length) :: names(10), npa(3)
      	type(nuclide_set), dimension(10), intent(out) :: set
      	integer :: i
      	call set_up_cno(Z, A)
      	call generate_nuclide_names(Z,A,names)
      	! move n, p, he4 to the end
      	npa = names([1,2,4])
      	! pull off n,p
      	names = cshift(names,2)
      	names(2:7) = names(3:8)
      	names(8:10) = npa
      	call generate_nuclide_set(names,set)
      end subroutine init_cno


      subroutine makenet_686(nuclei_names) ! used in X-ray burst calculations
         character(len=5), dimension(686), intent(out) :: nuclei_names
         nuclei_names(1) = '    n'
         nuclei_names(2) = '    p'
         nuclei_names(3) = '    d'
         nuclei_names(4) = '    t'
         nuclei_names(5) = '  he3'
         nuclei_names(6) = '  he4'
         nuclei_names(7) = '  he6'
         nuclei_names(8) = '  li6'
         nuclei_names(9) = '  li7'
         nuclei_names(10) = '  li8'
         nuclei_names(11) = '  li9'
         nuclei_names(12) = '  be7'
         nuclei_names(13) = '  be9'
         nuclei_names(14) = '   b8'
         nuclei_names(15) = '  b10'
         nuclei_names(16) = '  b11'
         nuclei_names(17) = '   c9'
         nuclei_names(18) = '  c10'
         nuclei_names(19) = '  c11'
         nuclei_names(20) = '  c12'
         nuclei_names(21) = '  c13'
         nuclei_names(22) = '  c14'
         nuclei_names(23) = '  n12'
         nuclei_names(24) = '  n13'
         nuclei_names(25) = '  n14'
         nuclei_names(26) = '  n15'
         nuclei_names(27) = '  o13'
         nuclei_names(28) = '  o14'
         nuclei_names(29) = '  o15'
         nuclei_names(30) = '  o16'
         nuclei_names(31) = '  o17'
         nuclei_names(32) = '  o18'
         nuclei_names(33) = '  f17'
         nuclei_names(34) = '  f18'
         nuclei_names(35) = '  f19'
         nuclei_names(36) = ' ne17'
         nuclei_names(37) = ' ne18'
         nuclei_names(38) = ' ne19'
         nuclei_names(39) = ' ne20'
         nuclei_names(40) = ' ne21'
         nuclei_names(41) = ' ne22'
         nuclei_names(42) = ' ne23'
         nuclei_names(43) = ' na20'
         nuclei_names(44) = ' na21'
         nuclei_names(45) = ' na22'
         nuclei_names(46) = ' na23'
         nuclei_names(47) = ' mg20'
         nuclei_names(48) = ' mg21'
         nuclei_names(49) = ' mg22'
         nuclei_names(50) = ' mg23'
         nuclei_names(51) = ' mg24'
         nuclei_names(52) = ' mg25'
         nuclei_names(53) = ' mg26'
         nuclei_names(54) = ' al22'
         nuclei_names(55) = ' al23'
         nuclei_names(56) = ' al24'
         nuclei_names(57) = ' al25'
         nuclei_names(58) = ' al26'
         nuclei_names(59) = ' al27'
         nuclei_names(60) = ' si22'
         nuclei_names(61) = ' si23'
         nuclei_names(62) = ' si24'
         nuclei_names(63) = ' si25'
         nuclei_names(64) = ' si26'
         nuclei_names(65) = ' si27'
         nuclei_names(66) = ' si28'
         nuclei_names(67) = ' si29'
         nuclei_names(68) = ' si30'
         nuclei_names(69) = '  p26'
         nuclei_names(70) = '  p27'
         nuclei_names(71) = '  p28'
         nuclei_names(72) = '  p29'
         nuclei_names(73) = '  p30'
         nuclei_names(74) = '  p31'
         nuclei_names(75) = '  s27'
         nuclei_names(76) = '  s28'
         nuclei_names(77) = '  s29'
         nuclei_names(78) = '  s30'
         nuclei_names(79) = '  s31'
         nuclei_names(80) = '  s32'
         nuclei_names(81) = '  s33'
         nuclei_names(82) = '  s34'
         nuclei_names(83) = '  s35'
         nuclei_names(84) = '  s36'
         nuclei_names(85) = ' cl31'
         nuclei_names(86) = ' cl32'
         nuclei_names(87) = ' cl33'
         nuclei_names(88) = ' cl34'
         nuclei_names(89) = ' cl35'
         nuclei_names(90) = ' cl36'
         nuclei_names(91) = ' cl37'
         nuclei_names(92) = ' ar31'
         nuclei_names(93) = ' ar32'
         nuclei_names(94) = ' ar33'
         nuclei_names(95) = ' ar34'
         nuclei_names(96) = ' ar35'
         nuclei_names(97) = ' ar36'
         nuclei_names(98) = ' ar37'
         nuclei_names(99) = ' ar38'
         nuclei_names(100) = ' ar39'
         nuclei_names(101) = ' ar40'
         nuclei_names(102) = '  k35'
         nuclei_names(103) = '  k36'
         nuclei_names(104) = '  k37'
         nuclei_names(105) = '  k38'
         nuclei_names(106) = '  k39'
         nuclei_names(107) = '  k40'
         nuclei_names(108) = '  k41'
         nuclei_names(109) = ' ca35'
         nuclei_names(110) = ' ca36'
         nuclei_names(111) = ' ca37'
         nuclei_names(112) = ' ca38'
         nuclei_names(113) = ' ca39'
         nuclei_names(114) = ' ca40'
         nuclei_names(115) = ' ca41'
         nuclei_names(116) = ' ca42'
         nuclei_names(117) = ' ca43'
         nuclei_names(118) = ' ca44'
         nuclei_names(119) = ' sc40'
         nuclei_names(120) = ' sc41'
         nuclei_names(121) = ' sc42'
         nuclei_names(122) = ' sc43'
         nuclei_names(123) = ' sc44'
         nuclei_names(124) = ' sc45'
         nuclei_names(125) = ' ti39'
         nuclei_names(126) = ' ti40'
         nuclei_names(127) = ' ti41'
         nuclei_names(128) = ' ti42'
         nuclei_names(129) = ' ti43'
         nuclei_names(130) = ' ti44'
         nuclei_names(131) = ' ti45'
         nuclei_names(132) = ' ti46'
         nuclei_names(133) = ' ti47'
         nuclei_names(134) = ' ti48'
         nuclei_names(135) = ' ti49'
         nuclei_names(136) = ' ti50'
         nuclei_names(137) = '  v43'
         nuclei_names(138) = '  v44'
         nuclei_names(139) = '  v45'
         nuclei_names(140) = '  v46'
         nuclei_names(141) = '  v47'
         nuclei_names(142) = '  v48'
         nuclei_names(143) = '  v49'
         nuclei_names(144) = '  v50'
         nuclei_names(145) = '  v51'
         nuclei_names(146) = ' cr43'
         nuclei_names(147) = ' cr44'
         nuclei_names(148) = ' cr45'
         nuclei_names(149) = ' cr46'
         nuclei_names(150) = ' cr47'
         nuclei_names(151) = ' cr48'
         nuclei_names(152) = ' cr49'
         nuclei_names(153) = ' cr50'
         nuclei_names(154) = ' cr51'
         nuclei_names(155) = ' cr52'
         nuclei_names(156) = ' cr53'
         nuclei_names(157) = ' cr54'
         nuclei_names(158) = ' mn46'
         nuclei_names(159) = ' mn47'
         nuclei_names(160) = ' mn48'
         nuclei_names(161) = ' mn49'
         nuclei_names(162) = ' mn50'
         nuclei_names(163) = ' mn51'
         nuclei_names(164) = ' mn52'
         nuclei_names(165) = ' mn53'
         nuclei_names(166) = ' mn54'
         nuclei_names(167) = ' mn55'
         nuclei_names(168) = ' fe46'
         nuclei_names(169) = ' fe47'
         nuclei_names(170) = ' fe48'
         nuclei_names(171) = ' fe49'
         nuclei_names(172) = ' fe50'
         nuclei_names(173) = ' fe51'
         nuclei_names(174) = ' fe52'
         nuclei_names(175) = ' fe53'
         nuclei_names(176) = ' fe54'
         nuclei_names(177) = ' fe55'
         nuclei_names(178) = ' fe56'
         nuclei_names(179) = ' fe57'
         nuclei_names(180) = ' fe58'
         nuclei_names(181) = ' co50'
         nuclei_names(182) = ' co51'
         nuclei_names(183) = ' co52'
         nuclei_names(184) = ' co53'
         nuclei_names(185) = ' co54'
         nuclei_names(186) = ' co55'
         nuclei_names(187) = ' co56'
         nuclei_names(188) = ' co57'
         nuclei_names(189) = ' co58'
         nuclei_names(190) = ' co59'
         nuclei_names(191) = ' ni50'
         nuclei_names(192) = ' ni51'
         nuclei_names(193) = ' ni52'
         nuclei_names(194) = ' ni53'
         nuclei_names(195) = ' ni54'
         nuclei_names(196) = ' ni55'
         nuclei_names(197) = ' ni56'
         nuclei_names(198) = ' ni57'
         nuclei_names(199) = ' ni58'
         nuclei_names(200) = ' ni59'
         nuclei_names(201) = ' ni60'
         nuclei_names(202) = ' ni61'
         nuclei_names(203) = ' ni62'
         nuclei_names(204) = ' ni63'
         nuclei_names(205) = ' ni64'
         nuclei_names(206) = ' cu56'
         nuclei_names(207) = ' cu57'
         nuclei_names(208) = ' cu58'
         nuclei_names(209) = ' cu59'
         nuclei_names(210) = ' cu60'
         nuclei_names(211) = ' cu61'
         nuclei_names(212) = ' cu62'
         nuclei_names(213) = ' cu63'
         nuclei_names(214) = ' cu64'
         nuclei_names(215) = ' cu65'
         nuclei_names(216) = ' zn55'
         nuclei_names(217) = ' zn56'
         nuclei_names(218) = ' zn57'
         nuclei_names(219) = ' zn58'
         nuclei_names(220) = ' zn59'
         nuclei_names(221) = ' zn60'
         nuclei_names(222) = ' zn61'
         nuclei_names(223) = ' zn62'
         nuclei_names(224) = ' zn63'
         nuclei_names(225) = ' zn64'
         nuclei_names(226) = ' zn65'
         nuclei_names(227) = ' zn66'
         nuclei_names(228) = ' zn67'
         nuclei_names(229) = ' zn68'
         nuclei_names(230) = ' zn69'
         nuclei_names(231) = ' zn70'
         nuclei_names(232) = ' ga60'
         nuclei_names(233) = ' ga61'
         nuclei_names(234) = ' ga62'
         nuclei_names(235) = ' ga63'
         nuclei_names(236) = ' ga64'
         nuclei_names(237) = ' ga65'
         nuclei_names(238) = ' ga66'
         nuclei_names(239) = ' ga67'
         nuclei_names(240) = ' ga68'
         nuclei_names(241) = ' ga69'
         nuclei_names(242) = ' ga70'
         nuclei_names(243) = ' ga71'
         nuclei_names(244) = ' ga72'
         nuclei_names(245) = ' ge59'
         nuclei_names(246) = ' ge60'
         nuclei_names(247) = ' ge61'
         nuclei_names(248) = ' ge62'
         nuclei_names(249) = ' ge63'
         nuclei_names(250) = ' ge64'
         nuclei_names(251) = ' ge65'
         nuclei_names(252) = ' ge66'
         nuclei_names(253) = ' ge67'
         nuclei_names(254) = ' ge68'
         nuclei_names(255) = ' ge69'
         nuclei_names(256) = ' ge70'
         nuclei_names(257) = ' ge71'
         nuclei_names(258) = ' ge72'
         nuclei_names(259) = ' ge73'
         nuclei_names(260) = ' ge74'
         nuclei_names(261) = ' ge75'
         nuclei_names(262) = ' ge76'
         nuclei_names(263) = ' as64'
         nuclei_names(264) = ' as65'
         nuclei_names(265) = ' as66'
         nuclei_names(266) = ' as67'
         nuclei_names(267) = ' as68'
         nuclei_names(268) = ' as69'
         nuclei_names(269) = ' as70'
         nuclei_names(270) = ' as71'
         nuclei_names(271) = ' as72'
         nuclei_names(272) = ' as73'
         nuclei_names(273) = ' as74'
         nuclei_names(274) = ' as75'
         nuclei_names(275) = ' as76'
         nuclei_names(276) = ' as77'
         nuclei_names(277) = ' se65'
         nuclei_names(278) = ' se66'
         nuclei_names(279) = ' se67'
         nuclei_names(280) = ' se68'
         nuclei_names(281) = ' se69'
         nuclei_names(282) = ' se70'
         nuclei_names(283) = ' se71'
         nuclei_names(284) = ' se72'
         nuclei_names(285) = ' se73'
         nuclei_names(286) = ' se74'
         nuclei_names(287) = ' se75'
         nuclei_names(288) = ' se76'
         nuclei_names(289) = ' se77'
         nuclei_names(290) = ' se78'
         nuclei_names(291) = ' se79'
         nuclei_names(292) = ' se80'
         nuclei_names(293) = ' se81'
         nuclei_names(294) = ' se82'
         nuclei_names(295) = ' br70'
         nuclei_names(296) = ' br71'
         nuclei_names(297) = ' br72'
         nuclei_names(298) = ' br73'
         nuclei_names(299) = ' br74'
         nuclei_names(300) = ' br75'
         nuclei_names(301) = ' br76'
         nuclei_names(302) = ' br77'
         nuclei_names(303) = ' br78'
         nuclei_names(304) = ' br79'
         nuclei_names(305) = ' br80'
         nuclei_names(306) = ' br81'
         nuclei_names(307) = ' br82'
         nuclei_names(308) = ' br83'
         nuclei_names(309) = ' kr67'
         nuclei_names(310) = ' kr68'
         nuclei_names(311) = ' kr69'
         nuclei_names(312) = ' kr70'
         nuclei_names(313) = ' kr71'
         nuclei_names(314) = ' kr72'
         nuclei_names(315) = ' kr73'
         nuclei_names(316) = ' kr74'
         nuclei_names(317) = ' kr75'
         nuclei_names(318) = ' kr76'
         nuclei_names(319) = ' kr77'
         nuclei_names(320) = ' kr78'
         nuclei_names(321) = ' kr79'
         nuclei_names(322) = ' kr80'
         nuclei_names(323) = ' kr81'
         nuclei_names(324) = ' kr82'
         nuclei_names(325) = ' kr83'
         nuclei_names(326) = ' kr84'
         nuclei_names(327) = ' kr85'
         nuclei_names(328) = ' kr86'
         nuclei_names(329) = ' rb74'
         nuclei_names(330) = ' rb75'
         nuclei_names(331) = ' rb76'
         nuclei_names(332) = ' rb77'
         nuclei_names(333) = ' rb78'
         nuclei_names(334) = ' rb79'
         nuclei_names(335) = ' rb80'
         nuclei_names(336) = ' rb81'
         nuclei_names(337) = ' rb82'
         nuclei_names(338) = ' rb83'
         nuclei_names(339) = ' rb84'
         nuclei_names(340) = ' rb85'
         nuclei_names(341) = ' rb86'
         nuclei_names(342) = ' rb87'
         nuclei_names(343) = ' rb88'
         nuclei_names(344) = ' rb89'
         nuclei_names(345) = ' sr72'
         nuclei_names(346) = ' sr73'
         nuclei_names(347) = ' sr74'
         nuclei_names(348) = ' sr75'
         nuclei_names(349) = ' sr76'
         nuclei_names(350) = ' sr77'
         nuclei_names(351) = ' sr78'
         nuclei_names(352) = ' sr79'
         nuclei_names(353) = ' sr80'
         nuclei_names(354) = ' sr81'
         nuclei_names(355) = ' sr82'
         nuclei_names(356) = ' sr83'
         nuclei_names(357) = ' sr84'
         nuclei_names(358) = ' sr85'
         nuclei_names(359) = ' sr86'
         nuclei_names(360) = ' sr87'
         nuclei_names(361) = ' sr88'
         nuclei_names(362) = ' sr89'
         nuclei_names(363) = ' sr90'
         nuclei_names(364) = ' sr91'
         nuclei_names(365) = '  y77'
         nuclei_names(366) = '  y78'
         nuclei_names(367) = '  y79'
         nuclei_names(368) = '  y80'
         nuclei_names(369) = '  y81'
         nuclei_names(370) = '  y82'
         nuclei_names(371) = '  y83'
         nuclei_names(372) = '  y84'
         nuclei_names(373) = '  y85'
         nuclei_names(374) = '  y86'
         nuclei_names(375) = '  y87'
         nuclei_names(376) = '  y88'
         nuclei_names(377) = '  y89'
         nuclei_names(378) = '  y90'
         nuclei_names(379) = '  y91'
         nuclei_names(380) = '  y92'
         nuclei_names(381) = '  y93'
         nuclei_names(382) = '  y94'
         nuclei_names(383) = ' zr76'
         nuclei_names(384) = ' zr77'
         nuclei_names(385) = ' zr78'
         nuclei_names(386) = ' zr79'
         nuclei_names(387) = ' zr80'
         nuclei_names(388) = ' zr81'
         nuclei_names(389) = ' zr82'
         nuclei_names(390) = ' zr83'
         nuclei_names(391) = ' zr84'
         nuclei_names(392) = ' zr85'
         nuclei_names(393) = ' zr86'
         nuclei_names(394) = ' zr87'
         nuclei_names(395) = ' zr88'
         nuclei_names(396) = ' zr89'
         nuclei_names(397) = ' zr90'
         nuclei_names(398) = ' zr91'
         nuclei_names(399) = ' zr92'
         nuclei_names(400) = ' zr93'
         nuclei_names(401) = ' zr94'
         nuclei_names(402) = ' zr95'
         nuclei_names(403) = ' nb82'
         nuclei_names(404) = ' nb83'
         nuclei_names(405) = ' nb84'
         nuclei_names(406) = ' nb85'
         nuclei_names(407) = ' nb86'
         nuclei_names(408) = ' nb87'
         nuclei_names(409) = ' nb88'
         nuclei_names(410) = ' nb89'
         nuclei_names(411) = ' nb90'
         nuclei_names(412) = ' nb91'
         nuclei_names(413) = ' nb92'
         nuclei_names(414) = ' nb93'
         nuclei_names(415) = ' nb94'
         nuclei_names(416) = ' nb95'
         nuclei_names(417) = ' nb96'
         nuclei_names(418) = ' nb97'
         nuclei_names(419) = ' mo80'
         nuclei_names(420) = ' mo81'
         nuclei_names(421) = ' mo82'
         nuclei_names(422) = ' mo83'
         nuclei_names(423) = ' mo84'
         nuclei_names(424) = ' mo85'
         nuclei_names(425) = ' mo86'
         nuclei_names(426) = ' mo87'
         nuclei_names(427) = ' mo88'
         nuclei_names(428) = ' mo89'
         nuclei_names(429) = ' mo90'
         nuclei_names(430) = ' mo91'
         nuclei_names(431) = ' mo92'
         nuclei_names(432) = ' mo93'
         nuclei_names(433) = ' mo94'
         nuclei_names(434) = ' mo95'
         nuclei_names(435) = ' mo96'
         nuclei_names(436) = ' mo97'
         nuclei_names(437) = ' mo98'
         nuclei_names(438) = ' tc86'
         nuclei_names(439) = ' tc87'
         nuclei_names(440) = ' tc88'
         nuclei_names(441) = ' tc89'
         nuclei_names(442) = ' tc90'
         nuclei_names(443) = ' tc91'
         nuclei_names(444) = ' tc92'
         nuclei_names(445) = ' tc93'
         nuclei_names(446) = ' tc94'
         nuclei_names(447) = ' tc95'
         nuclei_names(448) = ' tc96'
         nuclei_names(449) = ' tc97'
         nuclei_names(450) = ' tc98'
         nuclei_names(451) = ' tc99'
         nuclei_names(452) = ' ru86'
         nuclei_names(453) = ' ru87'
         nuclei_names(454) = ' ru88'
         nuclei_names(455) = ' ru89'
         nuclei_names(456) = ' ru90'
         nuclei_names(457) = ' ru91'
         nuclei_names(458) = ' ru92'
         nuclei_names(459) = ' ru93'
         nuclei_names(460) = ' ru94'
         nuclei_names(461) = ' ru95'
         nuclei_names(462) = ' ru96'
         nuclei_names(463) = ' ru97'
         nuclei_names(464) = ' ru98'
         nuclei_names(465) = ' ru99'
         nuclei_names(466) = 'ru100'
         nuclei_names(467) = 'ru101'
         nuclei_names(468) = ' rh89'
         nuclei_names(469) = ' rh90'
         nuclei_names(470) = ' rh91'
         nuclei_names(471) = ' rh92'
         nuclei_names(472) = ' rh93'
         nuclei_names(473) = ' rh94'
         nuclei_names(474) = ' rh95'
         nuclei_names(475) = ' rh96'
         nuclei_names(476) = ' rh97'
         nuclei_names(477) = ' rh98'
         nuclei_names(478) = ' rh99'
         nuclei_names(479) = 'rh100'
         nuclei_names(480) = 'rh101'
         nuclei_names(481) = 'rh102'
         nuclei_names(482) = 'rh103'
         nuclei_names(483) = ' pd88'
         nuclei_names(484) = ' pd89'
         nuclei_names(485) = ' pd90'
         nuclei_names(486) = ' pd91'
         nuclei_names(487) = ' pd92'
         nuclei_names(488) = ' pd93'
         nuclei_names(489) = ' pd94'
         nuclei_names(490) = ' pd95'
         nuclei_names(491) = ' pd96'
         nuclei_names(492) = ' pd97'
         nuclei_names(493) = ' pd98'
         nuclei_names(494) = ' pd99'
         nuclei_names(495) = 'pd100'
         nuclei_names(496) = 'pd101'
         nuclei_names(497) = 'pd102'
         nuclei_names(498) = 'pd103'
         nuclei_names(499) = 'pd104'
         nuclei_names(500) = 'pd105'
         nuclei_names(501) = ' ag93'
         nuclei_names(502) = ' ag94'
         nuclei_names(503) = ' ag95'
         nuclei_names(504) = ' ag96'
         nuclei_names(505) = ' ag97'
         nuclei_names(506) = ' ag98'
         nuclei_names(507) = ' ag99'
         nuclei_names(508) = 'ag100'
         nuclei_names(509) = 'ag101'
         nuclei_names(510) = 'ag102'
         nuclei_names(511) = 'ag103'
         nuclei_names(512) = 'ag104'
         nuclei_names(513) = 'ag105'
         nuclei_names(514) = 'ag106'
         nuclei_names(515) = 'ag107'
         nuclei_names(516) = 'ag108'
         nuclei_names(517) = 'ag109'
         nuclei_names(518) = ' cd92'
         nuclei_names(519) = ' cd93'
         nuclei_names(520) = ' cd94'
         nuclei_names(521) = ' cd95'
         nuclei_names(522) = ' cd96'
         nuclei_names(523) = ' cd97'
         nuclei_names(524) = ' cd98'
         nuclei_names(525) = ' cd99'
         nuclei_names(526) = 'cd100'
         nuclei_names(527) = 'cd101'
         nuclei_names(528) = 'cd102'
         nuclei_names(529) = 'cd103'
         nuclei_names(530) = 'cd104'
         nuclei_names(531) = 'cd105'
         nuclei_names(532) = 'cd106'
         nuclei_names(533) = 'cd107'
         nuclei_names(534) = 'cd108'
         nuclei_names(535) = 'cd109'
         nuclei_names(536) = 'cd110'
         nuclei_names(537) = ' in98'
         nuclei_names(538) = ' in99'
         nuclei_names(539) = 'in100'
         nuclei_names(540) = 'in101'
         nuclei_names(541) = 'in102'
         nuclei_names(542) = 'in103'
         nuclei_names(543) = 'in104'
         nuclei_names(544) = 'in105'
         nuclei_names(545) = 'in106'
         nuclei_names(546) = 'in107'
         nuclei_names(547) = 'in108'
         nuclei_names(548) = 'in109'
         nuclei_names(549) = 'in110'
         nuclei_names(550) = 'in111'
         nuclei_names(551) = 'in112'
         nuclei_names(552) = 'in113'
         nuclei_names(553) = 'in114'
         nuclei_names(554) = 'in115'
         nuclei_names(555) = 'in116'
         nuclei_names(556) = 'in117'
         nuclei_names(557) = ' sn96'
         nuclei_names(558) = ' sn97'
         nuclei_names(559) = ' sn98'
         nuclei_names(560) = ' sn99'
         nuclei_names(561) = 'sn100'
         nuclei_names(562) = 'sn101'
         nuclei_names(563) = 'sn102'
         nuclei_names(564) = 'sn103'
         nuclei_names(565) = 'sn104'
         nuclei_names(566) = 'sn105'
         nuclei_names(567) = 'sn106'
         nuclei_names(568) = 'sn107'
         nuclei_names(569) = 'sn108'
         nuclei_names(570) = 'sn109'
         nuclei_names(571) = 'sn110'
         nuclei_names(572) = 'sn111'
         nuclei_names(573) = 'sn112'
         nuclei_names(574) = 'sn113'
         nuclei_names(575) = 'sn114'
         nuclei_names(576) = 'sn115'
         nuclei_names(577) = 'sn116'
         nuclei_names(578) = 'sn117'
         nuclei_names(579) = 'sn118'
         nuclei_names(580) = 'sn119'
         nuclei_names(581) = 'sn120'
         nuclei_names(582) = 'sn121'
         nuclei_names(583) = 'sn122'
         nuclei_names(584) = 'sn124'
         nuclei_names(585) = 'sb104'
         nuclei_names(586) = 'sb105'
         nuclei_names(587) = 'sb106'
         nuclei_names(588) = 'sb107'
         nuclei_names(589) = 'sb108'
         nuclei_names(590) = 'sb109'
         nuclei_names(591) = 'sb110'
         nuclei_names(592) = 'sb111'
         nuclei_names(593) = 'sb112'
         nuclei_names(594) = 'sb113'
         nuclei_names(595) = 'sb114'
         nuclei_names(596) = 'sb115'
         nuclei_names(597) = 'sb116'
         nuclei_names(598) = 'sb117'
         nuclei_names(599) = 'sb118'
         nuclei_names(600) = 'sb119'
         nuclei_names(601) = 'sb120'
         nuclei_names(602) = 'sb121'
         nuclei_names(603) = 'sb122'
         nuclei_names(604) = 'sb123'
         nuclei_names(605) = 'sb124'
         nuclei_names(606) = 'sb125'
         nuclei_names(607) = 'te104'
         nuclei_names(608) = 'te105'
         nuclei_names(609) = 'te106'
         nuclei_names(610) = 'te107'
         nuclei_names(611) = 'te108'
         nuclei_names(612) = 'te109'
         nuclei_names(613) = 'te110'
         nuclei_names(614) = 'te111'
         nuclei_names(615) = 'te112'
         nuclei_names(616) = 'te113'
         nuclei_names(617) = 'te114'
         nuclei_names(618) = 'te115'
         nuclei_names(619) = 'te116'
         nuclei_names(620) = 'te117'
         nuclei_names(621) = 'te118'
         nuclei_names(622) = 'te119'
         nuclei_names(623) = 'te120'
         nuclei_names(624) = 'te121'
         nuclei_names(625) = 'te122'
         nuclei_names(626) = 'te123'
         nuclei_names(627) = 'te124'
         nuclei_names(628) = 'te125'
         nuclei_names(629) = 'te126'
         nuclei_names(630) = 'te127'
         nuclei_names(631) = 'te128'
         nuclei_names(632) = 'te129'
         nuclei_names(633) = 'te130'
         nuclei_names(634) = ' i108'
         nuclei_names(635) = ' i109'
         nuclei_names(636) = ' i110'
         nuclei_names(637) = ' i111'
         nuclei_names(638) = ' i112'
         nuclei_names(639) = ' i113'
         nuclei_names(640) = ' i114'
         nuclei_names(641) = ' i115'
         nuclei_names(642) = ' i116'
         nuclei_names(643) = ' i117'
         nuclei_names(644) = ' i118'
         nuclei_names(645) = ' i119'
         nuclei_names(646) = ' i120'
         nuclei_names(647) = ' i121'
         nuclei_names(648) = ' i122'
         nuclei_names(649) = ' i123'
         nuclei_names(650) = ' i124'
         nuclei_names(651) = ' i125'
         nuclei_names(652) = ' i126'
         nuclei_names(653) = ' i127'
         nuclei_names(654) = ' i128'
         nuclei_names(655) = ' i129'
         nuclei_names(656) = ' i130'
         nuclei_names(657) = ' i131'
         nuclei_names(658) = 'xe108'
         nuclei_names(659) = 'xe109'
         nuclei_names(660) = 'xe110'
         nuclei_names(661) = 'xe111'
         nuclei_names(662) = 'xe112'
         nuclei_names(663) = 'xe113'
         nuclei_names(664) = 'xe114'
         nuclei_names(665) = 'xe115'
         nuclei_names(666) = 'xe116'
         nuclei_names(667) = 'xe117'
         nuclei_names(668) = 'xe118'
         nuclei_names(669) = 'xe119'
         nuclei_names(670) = 'xe120'
         nuclei_names(671) = 'xe121'
         nuclei_names(672) = 'xe122'
         nuclei_names(673) = 'xe123'
         nuclei_names(674) = 'xe124'
         nuclei_names(675) = 'xe125'
         nuclei_names(676) = 'xe126'
         nuclei_names(677) = 'xe127'
         nuclei_names(678) = 'xe128'
         nuclei_names(679) = 'xe129'
         nuclei_names(680) = 'xe130'
         nuclei_names(681) = 'xe131'
         nuclei_names(682) = 'xe132'
         nuclei_names(683) = 'xe133'
         nuclei_names(684) = 'xe134'
         nuclei_names(685) = 'xe135'
         nuclei_names(686) = 'xe136'
      end subroutine makenet_686


      subroutine init_686(set)
      	character(len=iso_name_length) :: names(686), npa(3)
      	type(nuclide_set), dimension(686), intent(out) :: set
      	integer :: i   
      	call makenet_686(names)   	
      	! move n, p, he4 to the end
      	npa = names([1,2,6])
      	! pull off n,p
      	names = cshift(names,2)
      	names(4:683) = names(5:684)
      	names(684:686) = npa
      	call generate_nuclide_set(names,set)	
      end subroutine init_686
      
      
      subroutine test_burn(nnuc)
         use num_def
         use num_lib 
         use mtx_lib
         use rates_def, only: num_categories
      	integer, intent(in) :: nnuc ! 430 !200 !10 ! 10
      	character(len=80):: message
      	type(nuclide_data) :: nuclides
      	type(reaction_data) :: rates
      	type(nuclide_set), pointer :: set(:) ! (nnuc)
      	integer :: i, ierr, ios, length, k, j
      	double precision, dimension(:), pointer :: X, Y, Ynew, Xnew ! (nnuc)
      	double precision :: rho, temp, time, tend, deltaB, abar, zbar, z2bar, ye
      	logical :: use_graboske_et_al_screening, use_weaklib
      	integer :: caller_id, handle, jh1, jn, jh2, jhe4, jc12, jc13, jn13, jn14, jn15, jo15, jo16
      	integer, pointer :: chem_id(:)
         double precision :: category_factors(num_categories)
      	
      	! for jina_1_zone_burn
      	integer, parameter :: num_times = 1
      	integer :: which_solver, which_decsol, max_steps, iout, &
      	   nfcn, njac, nstep, naccpt, nrejct
         double precision, dimension(num_times) :: &
            times, temperatures, densities, thetas, etas, d_eta_dlnTs, d_eta_dlnRhos
         double precision, target :: y_history_array(nnuc,num_times)
         double precision, pointer :: y_history(:,:)
         double precision :: rtol, atol, max_step_size, h
      	
      	include 'formats.dek'
      	
      	allocate(set(nnuc), X(nnuc), Y(nnuc), Ynew(nnuc), Xnew(nnuc))
      	
      	select case(nnuc)
      	case(10)
         	call init_cno(set)
      	case(25)
         	call init_25(set)
      	case(100)
         	call init_100(set)
      	case(200)
         	call init_200(set)
      	case(430)
         	call init_430(set)
      	case(686)
         	call init_686(set)
      	case(690)
         	call init_690(set)
      	case default
      	   write(*,*) 'bad nnuc for test_burn'
      	   stop 1
      	end select
      	
      	jh1 = get_nuclide_index_in_set('h1', set)
      	jn = get_nuclide_index_in_set('h2', set)
      	jh2 = get_nuclide_index_in_set('h3', set)
      	jhe4 = get_nuclide_index_in_set('he4', set)
      	jc12 = get_nuclide_index_in_set('c12', set)
      	jc13 = get_nuclide_index_in_set('c13', set)
      	jn13 = get_nuclide_index_in_set('n13', set)
      	jn14 = get_nuclide_index_in_set('n14', set)
      	jn15 = get_nuclide_index_in_set('n15', set)
      	jo15 = get_nuclide_index_in_set('o15', set)
      	jo16 = get_nuclide_index_in_set('o16', set)
	
      	call extract_nuclides_from_chem_isos(set,nuclides,chem_id,ierr)
      	if (ierr /= 0) call alert(ierr,'problem extracting nuclides')

         use_weaklib = .true.
      	call extract_reaclib_rates(set, nuclides, rates, use_weaklib, ierr)
      	if (ierr /= 0) call alert(ierr,'problem extracting rates')
      	!call output_rates(6,rates,nuclides,pretty_print_format)
      	
      	category_factors = 1
      	time = 0
         rtol = 1d-7
         atol = 1d-8
      	
      	X = 0
      	if (nnuc > 10) then
         	X(jhe4) = 0.95d0
         	X(jc12) = 0.025d0
         	X(jo16) = 0.025d0
         	tend = 100 ! seconds
         	!rho = 1.5d5
         	!temp = 2d9
         	rho = 1d7
         	temp = 4d9
      	else
         	X(jh1) = 0.98d0
         	X(jn14) = 0.02d0
         	tend = 1d2*secyer
         	rho = 1d4
         	temp = 4d7
      	end if
      	
      	Y = X/nuclides% A
      	
      	call nuclides_composition_info(Y, nuclides, abar, zbar, z2bar, ye)
	      
	      write(*,*)
	      write(*,*)
	      write(*,2) 'num nuclides', nnuc
	      write(*,*)
	      write(*,1) 'rho', rho
	      write(*,1) 'temp', temp
	      write(*,1) 'tend', tend
	      write(*,1) 'tend/secyer', tend/secyer
	      write(*,*)
	      write(*,2) 'initial abundances'
	      call show_X(X)
	      write(*,*)
	      write(*,1) 'ye', ye
	      write(*,1) 'abar', abar
	      write(*,1) 'zbar', zbar
	      write(*,1) 'z2bar', z2bar
	      write(*,*)
      	
      	write(*,*) 'call alloc_jina_burn'
      	caller_id = 0
      	which_decsol = sparskit
      	use_weaklib = .true.
	      handle = alloc_jina_burn_handle(which_decsol, rtol, atol, set, use_weaklib, ierr)
	      if (ierr /= 0) stop 'failed in alloc_jina_burn'
      	
      	y_history => y_history_array
      	!which_solver = 9 ! sodex_solver
      	!which_solver = 3 ! ros3p_solver
      	which_solver = 8 ! seulex_solver
      	h = 0
      	max_step_size = 0
      	max_steps = 10000
         iout = 1
         
         times(1) = tend
         temperatures(1) = temp
         densities(1) = rho
         thetas(1) = 0
         etas(1) = 0
         d_eta_dlnTs(1) = 0
         d_eta_dlnRhos(1) = 0
         use_graboske_et_al_screening = .true.

      	write(*,*) 'call jina_1_zone_burn'
	      call jina_1_zone_burn( &
            handle, caller_id, Y, which_solver, &
            num_times, times, temperatures, densities, &
            thetas, etas, d_eta_dlnTs, d_eta_dlnRhos, &
            use_graboske_et_al_screening, category_factors, &
            h, max_step_size, max_steps, & 
            jina_burn_solout, iout, y_history, nfcn, njac, nstep, naccpt, nrejct, ierr)
	      if (ierr /= 0) stop 'failed in jina_1_zone_burn'
	      Ynew = y_history(:,num_times)
      	
      	write(*,*) 'call free_jina_burn'
	      call free_jina_burn_handle(handle,ierr)
	      if (ierr /= 0) stop 'failed in free_jina_burn'
	      
      	Xnew= Ynew*nuclides% A
	      
	      if (nnuc > 10) then
   	      write(*,*)
   	      write(*,2) 'nfcn', nfcn
   	      write(*,2) 'njac', njac
   	      write(*,2) 'nstep', nstep
   	      write(*,2) 'naccpt', naccpt
   	      write(*,2) 'nrejct', nrejct
	      end if
	      
	      write(*,*)
	      write(*,2) 'new abundances'
	      call show_X(Xnew)
	      
      	
      	call nuclides_composition_info(Ynew, nuclides, abar, zbar, z2bar, ye)
	      write(*,*)
	      write(*,1) 'ye', ye
	      write(*,1) 'abar', abar
	      write(*,1) 'zbar', zbar
	      write(*,1) 'z2bar', z2bar
	      write(*,*)
	      
	      
      	call free_reaction_data(rates)
      	call free_nuclide_data(nuclides)

      	deallocate(set, X, Y, Ynew, Xnew)
      	
      	contains
      	
	
      	subroutine show_X(X)
      	   double precision :: X(:)
      	   include 'formats.dek'
      	   integer :: j
      	   double precision :: xsum
      	   xsum = 0
         	do j=1, nuclides% nnuclides
         	   if (x(j) > 1d-3) &
         	      write(*,1) trim(nuclides% name(j)), x(j)
         	   if (x(j) > 1.1 .or. x(j) < -0.1) then
         	      write(*,1) 'bad x for ' // trim(nuclides% name(j)), x(j)
         	      stop 1
         	   end if
         	   xsum = xsum + x(j)
         	end do
         	write(*,1) 'xsum', xsum
         	write(*,*)
   	   end subroutine show_X

      	
      end subroutine test_burn


      subroutine jina_burn_solout(nr, told, t, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
         ! nr is the step number.
         ! t is the current time value; told is the previous t value.
         ! y is the current y vector.
         ! irtrn negative means terminate integration.
         ! rwork_y and iwork_y hold info for interp_y
         ! note that these are not the same as the rwork and iwork arrays for the solver.
         integer, intent(in) :: nr, n, lrpar, lipar
         double precision, intent(in) :: told, t
         double precision, intent(inout) :: y(n)
         ! y can be modified if necessary to keep it in valid range of possible solutions.
         double precision, intent(inout), target :: rpar(lrpar), rwork_y(*)
         integer, intent(inout), target :: ipar(lipar), iwork_y(*)
         interface
            include 'num_interp_y.dek'
         end interface
         integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
         type (Jina_Info), pointer :: g
         integer :: handle, j, ierr
         double precision :: x, xsum
   	   include 'formats.dek'
         handle = ipar(i_handle)
         ierr = 0
      	call get_jina_burner_ptr(handle, g, ierr)
      	if (ierr /= 0) then
      	   irtrn = -1; return
      	end if
      	irtrn = 0
      	if (n > 10 .and. mod(nr,10) == 0) write(*,2) 'step t', nr, t
      	return
      	
      	
      	xsum = 0
      	do j=1, g% nuclides% nnuclides
      	   x = y(j)*g% nuclides% A(j)
      	   if (x > 1d-4) write(*,1) trim(g% nuclides% name(j)), x
      	   if (x > 1.1 .or. x < -0.1) then
      	      write(*,1) 'bad x: stop in jina_burn_solout ' // trim(g% nuclides% name(j)), x
      	      stop 1
      	   end if
      	   xsum = xsum + y(j)*g% nuclides% A(j)
      	end do
      	!write(*,2) 'nnuclides', g% nuclides% nnuclides
      	!write(*,1) 'xsum', xsum
      	if (xsum > 1.1) then
      	   write(*,*) 'bad xsum: stop in jina_burn_solout'
      	   stop 1
      	end if
         write(*,*)
      end subroutine jina_burn_solout
      
      
      subroutine test_jina_net_get
         use rates_def, only: i_rate, num_rvs, num_categories, category_name
      	type(nuclide_set), pointer :: set(:) ! (nnuc)
         type (Jina_Info), pointer :: g
         integer, dimension(:), pointer :: Z, A ! will be allocated
      	double precision, dimension(:), pointer :: X, Y ! (nnuc)
      	double precision :: &
      	   abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
      	   T, logT, rho, logRho, theta_e_for_graboske_et_al, &
            eps_nuc, d_eps_nuc_dlnRho, d_eps_nuc_dlnT, &
            eps_neu, d_eps_neu_dlnT, d_eps_neu_dlnRho
      	double precision, dimension(:), pointer :: &
      	   dxdt, d_eps_nuc_dx, d_dxdt_dlnRho, d_dxdt_dlnT ! (nnuc)
      	double precision, dimension(:,:), pointer :: d_dxdt_dx, reaction_eps_nuc, eps_nuc_categories
      	logical :: use_graboske_et_al_screening, use_weaklib
      	integer :: i, nnuc, ierr, handle, nr
      	integer :: &
            jneut, jh1, jh2, jh3, jhe3, jhe4, jli6, jli7, jbe7, jbe8, jbe9, jbe10, jb8, jb10, jb11, &
            jc12, jc13, jc14, jn13, jn14, jn15, jn16, jo14, jo15, jo16, jo17, jo18, jo19, &
            jf17, jf18, jf19, jf20, jne18, jne19, jne20, jne21, jne22, jne23, jna21, jna22, jna23, jna24, &
            jmg22, jmg23, jmg24, jmg25, jmg26, jmg27, jal25, jal26, jal27, jal28, &
            jsi27, jsi28, jsi29, jsi30, jsi31, jsi32, jp29, jp30, jp31, jp32, js31, js32
         double precision :: category_factors(num_categories)
      	character(len=iso_name_length), pointer :: names(:)
      	character (len=4096) :: fname
      
      	include 'formats.dek'
      	
      	!write(*,*) 'test_jina_net_get'
      	
      	ierr = 0
      	
      	category_factors = 1
      	
      	!fname = 'test_ZA.txt'
      	fname = '../../data/jina_data/net_A32.jina'
      	call chem_read_ZA(fname, Z, A, nnuc, ierr)
      	if (ierr /= 0) then
      	   write(*,*) 'failed in chem_read_ZA ' // trim(fname)
      	   stop 1
      	end if
      	
      	allocate(names(nnuc), set(nnuc), X(nnuc), Y(nnuc), dxdt(nnuc), &
      	   d_eps_nuc_dx(nnuc), d_dxdt_dlnRho(nnuc), d_dxdt_dlnT(nnuc), d_dxdt_dx(nnuc,nnuc))
      	
      	call generate_nuclide_names(Z,A,names)	
      	call generate_nuclide_set(names,set)
      	
      	!call write_nuclides
      	
      	call set_index_in_set
                  
         use_weaklib = .true.
         !use_weaklib = .false. ! for comparison with net
         
         write(*,*) 'use_weaklib', use_weaklib
         
         handle = alloc_jina_net_handle(set, use_weaklib, ierr)
      	if (ierr /= 0) return
      	
      	call get_jina_net_ptr(handle, g, ierr)
      	if (ierr /= 0) return
      	
      	nr = g% rates% nreactions
      	allocate(reaction_eps_nuc(num_rvs, nr), eps_nuc_categories(num_rvs, num_categories), &
      	   g% rates% forward_rate_num(nr), g% rates% reverse_rate_num(nr), &
      	   g% rates% forward_rate_factor(nr), g% rates% reverse_rate_factor(nr))
      	   
      	if (.false.) then ! only use selected rates
         	g% rates% num_forward_rate_factors = nr
         	g% rates% num_reverse_rate_factors = nr
         	forall (i=1:nr)
         	   g% rates% forward_rate_num(i) = i
         	   g% rates% reverse_rate_num(i) = i
         	end forall
         	if (.false.) then
            	g% rates% forward_rate_factor = 0
            	g% rates% reverse_rate_factor = 0
            	! now selectively turn some back on
            	do i=1,nr
            	   if (g% rates% reaction_handle(i) == 'r_h1_h2_to_he3') then
                  	g% rates% forward_rate_factor(i) = 1
                  	g% rates% reverse_rate_factor(i) = 1
            	   end if
            	end do
         	else
            	g% rates% forward_rate_factor = 1
            	g% rates% reverse_rate_factor = 1
            	! now selectively turn some off
            	do i=1,nr
            	   if (g% rates% reaction_handle(i) == 'r_he4_n15_to_f19') then
                  	g% rates% forward_rate_factor(i) = 0
                  	g% rates% reverse_rate_factor(i) = 0
            	   end if
            	   if (g% rates% reaction_handle(i) == 'r_he4_o16_to_ne20') then
                  	g% rates% forward_rate_factor(i) = 0
                  	g% rates% reverse_rate_factor(i) = 0
            	   end if
            	   if (g% rates% reaction_handle(i) == 'r_c12_c12_to_he4_ne20') then
                  	g% rates% forward_rate_factor(i) = 0
                  	g% rates% reverse_rate_factor(i) = 0
            	   end if
            	   if (g% rates% reaction_handle(i) == 'r_he4_f17_to_h1_ne20') then
                  	g% rates% forward_rate_factor(i) = 0
                  	g% rates% reverse_rate_factor(i) = 0
            	   end if
            	   if (g% rates% reaction_handle(i) == 'r_h2_be7_to_h1_he4_he4') then
                  	g% rates% forward_rate_factor(i) = 0
                  	g% rates% reverse_rate_factor(i) = 0
            	   end if
            	   if (g% rates% reaction_handle(i) == 'r_h2_he3_to_h1_he4') then
                  	g% rates% forward_rate_factor(i) = 0
                  	g% rates% reverse_rate_factor(i) = 0
            	   end if
            	   if (g% rates% reaction_handle(i) == 'r_he3_be7_to_h1_h1_he4_he4') then
                  	g% rates% forward_rate_factor(i) = 0
                  	g% rates% reverse_rate_factor(i) = 0
            	   end if
            	end do
         	end if
         end if
      	
      	if (.false.) then
         	write(*,*)
         	write(*,*) nnuc, 'nuclides'
         	do i = 1, nnuc
         	   write(*,*)  trim(names(i))
         	end do
         	write(*,*)
         	write(*,*) nr, 'rates'
         	do i=1,nr
         	   !if (g% rates% weak_mask(i) /= 0d0) cycle
         	   write(*,*) g% rates% reaction_handle(i)
         	end do
         	write(*,*)
            call output_rates(6,g% rates,g% nuclides,pretty_print_format)
         	stop 'test_jina_net_get'
      	end if


      	Y = X/g% nuclides% A

         use_graboske_et_al_screening = .true.
      	
      	call nuclides_composition_info(Y, g% nuclides, abar, zbar, z2bar, ye)



                              x(jneut)=     1.0000000002491754D-99
                                x(jh1)=     5.7065323898114011D-22
                                x(jh2)=     3.9990730133971550D-39
                                x(jh3)=     1.0000000002491754D-99
                               x(jhe3)=     3.5066319623102392D-42
                               x(jhe4)=     2.5547165429352720D-14
                               x(jli6)=     4.2281198243494472D-36
                               x(jli7)=     9.2303859640055423D-40
                               x(jbe7)=     6.5797980137297868D-37
                               x(jbe8)=     1.0000000002491754D-99
                               x(jbe9)=     1.1224105005208780D-40
                              x(jbe10)=     1.0000000002491754D-99
                                x(jb8)=     1.4047354460555623D-51
                               x(jb10)=     1.8369653669845318D-34
                               x(jb11)=     1.8210704354282600D-35
                               x(jc12)=     7.0590426262982196D-03
                               x(jc13)=     5.6308797885553970D-05
                               x(jc14)=     2.5175565222221229D-34
                               x(jn13)=     1.5727947851062827D-37
                               x(jn14)=     1.2767383346650613D-02
                               x(jn15)=     1.9717785040260599D-07
                               x(jn16)=     1.0000000000432924D-99
                               x(jo15)=     6.9496354703012468D-22
                               x(jo16)=     5.3449870801142374D-01
                               x(jo17)=     2.4124657464050893D-07
                               x(jo18)=     7.1241599420610863D-11
                               x(jo19)=     2.1013250663180505D-37
                               x(jf17)=     3.5559030954878191D-40
                               x(jf18)=     3.9110788433943048D-18
                               x(jf19)=     8.2845547286891774D-11
                               x(jf20)=     1.8937676158246688D-34
                              x(jne20)=     4.1817088702924177D-01
                              x(jne21)=     4.6346310346123262D-06
                              x(jne22)=     9.1730427356887773D-06
                              x(jne23)=     4.5579316339565057D-24
                              x(jna21)=     5.1181961327968276D-15
                              x(jna22)=     5.7944107761720402D-13
                              x(jna23)=     1.8691288446152685D-04
                              x(jna24)=     0.0000000000000000D+00
                              x(jmg23)=     0.0000000000000000D+00
                              x(jmg24)=     2.5539942127593723D-02
                              x(jmg25)=     7.8953614647545440D-05
                              x(jmg26)=     1.1368133968950456D-04
                              x(jmg27)=     5.7742923713829461D-22
                              x(jal25)=     6.1574039585880805D-15
                              x(jal26)=     2.9024202189705144D-14
                              x(jal27)=     6.5112005006150452D-05
                              x(jal28)=     0.0000000000000000D+00
                              x(jsi27)=     0.0000000000000000D+00
                              x(jsi28)=     1.0311511716025569D-03
                              x(jsi29)=     3.0842073773427939D-05
                              x(jsi30)=     2.1576086089527883D-05
                              x(jsi31)=     8.5397581410466909D-16
                              x(jsi32)=     1.3237382302609506D-16
                               x(jp29)=     9.8104071310305693D-18
                               x(jp30)=     7.7338174133713327D-16
                               x(jp31)=     7.4519168421373791D-06
                               x(jp32)=     1.6575577535197561D-14
                               x(js31)=     0.0000000000000000D+00
                               x(js32)=     3.5780071584807751D-04


                                  logT =    8.9649159615935936D+00
                                logRho =    5.7293487631211777D+00
                                    ye =    4.9998470161001080D-01
                                   eta =    1.7608098153072727D-01
                          d_eta_dlnRho =   -1.9141839505739286D-06
                            d_eta_dlnT =   -8.4837789790765888D-07
 use_graboske_et_al_screening = .false.
               theta_e_for_graboske_et_al =  0.0000000000000000D+00


             T =    10**logT
           rho =    10**logRho
         
         write(*,1) 'logT', logT
         write(*,1) 'logRho', logRho
         
         call jina_net_get(handle, &
               x, T, logT, rho, logRho, &
               abar, zbar, z2bar, eta, d_eta_dlnRho, d_eta_dlnT, category_factors, &
               eps_nuc, d_eps_nuc_dlnRho, d_eps_nuc_dlnT, d_eps_nuc_dx, &
               reaction_eps_nuc, eps_nuc_categories, &
               eps_neu, d_eps_neu_dlnT, d_eps_neu_dlnRho, &
               dxdt, d_dxdt_dlnRho, d_dxdt_dlnT, d_dxdt_dx, &
               use_graboske_et_al_screening, theta_e_for_graboske_et_al, &
               ierr)
         
         if (ierr /= 0) then
            write(*,*) 'jina_net_get ierr', ierr
            stop
         end if
      	
      	if (.true.) then
         	write(*,*)
      	   write(*,1) 'eps_neu', eps_neu
      	   write(*,1) 'd_eps_neu_dlnT', d_eps_neu_dlnT
      	   write(*,1) 'd_eps_neu_dlnRho', d_eps_neu_dlnRho
         	write(*,*)
      	   write(*,1) 'eps_nuc', eps_nuc
         	write(*,*)
         	do i = 1, num_categories
         	   if (abs(eps_nuc_categories(i_rate,i)) < 1d-20) cycle
         	   write(*,1)  'eps_nuc_cat ' // trim(category_name(i)), eps_nuc_categories(i_rate,i)
         	end do
         	write(*,*)
      	   do i = 1, nr
      	      write(*,1) 'eps nuc    ' // g% rates% reaction_handle(i), reaction_eps_nuc(i_rate,i)
      	   end do
         	write(*,*)
         	do i = 1, nnuc
         	   write(*,1)  'x ' // trim(names(i)), x(i)
         	end do
         	write(*,*)
         	do i = 1, nnuc
         	   write(*,1)  'dxdt ' // trim(names(i)), dxdt(i)
         	end do
         	write(*,*)
         	do i = 1, nnuc
         	   write(*,1)  'd_dxdt_dlnRho ' // trim(names(i)), d_dxdt_dlnRho(i)
         	end do
         	write(*,*)
         	do i = 1, nnuc
         	   write(*,1)  'd_dxdt_dlnT ' // trim(names(i)), d_dxdt_dlnT(i)
         	end do
         	write(*,*)
         	do i = 1, nnuc
         	   write(*,1)  'd_dxdt_dx(1,:) ' // trim(names(i)), d_dxdt_dx(1,i)
         	end do
         	write(*,*)
         	do i = 1, num_categories
         	   if (abs(eps_nuc_categories(i_rate,i)) < 1d-20) cycle
         	   write(*,1)  'eps_nuc_cat ' // trim(category_name(i)), eps_nuc_categories(i_rate,i)
         	end do
         	write(*,*)
      	   write(*,1) 'eps_neu', eps_neu
      	   write(*,1) 'd_eps_neu_dlnT', d_eps_neu_dlnT
      	   write(*,1) 'd_eps_neu_dlnRho', d_eps_neu_dlnRho
         	write(*,*)
      	   write(*,1) 'eps_nuc', eps_nuc
      	   write(*,1) 'd_epsnuc_dlnd', d_eps_nuc_dlnRho
      	   write(*,1) 'd_epsnuc_dlnT', d_eps_nuc_dlnT
         	do i = 1, nnuc
         	   write(*,1)  'd_eps_nuc_dx ' // trim(names(i)), d_eps_nuc_dx(i)
         	end do
         	write(*,*)
      	   stop
      	end if
      
	      call free_jina_net_handle(handle,ierr)

      	
      	deallocate(names, set, X, Y, dxdt, &
      	   d_eps_nuc_dx, d_dxdt_dlnRho, d_dxdt_dlnT, d_dxdt_dx, reaction_eps_nuc)
      	   
      	   
      	contains
      	
      	
      	subroutine write_nuclides
      	   integer :: j
      	   do j=1,size(names)
      	      write(6,fmt='(a)',advance='no') 'j' // trim(adjustl(names(j))) // ', '
      	   end do
      	   write(*,*)
      	   do j=1,size(names)
      	      write(6,fmt='(a)') 'j' // trim(adjustl(names(j))) // &
      	         ' =  get_nuclide_index_in_set("' // trim(adjustl(names(j))) // '", set)'
      	   end do
      	   write(*,*)
      	   stop 'write_nuclides'
      	end subroutine write_nuclides
      	
      	
      	subroutine set_index_in_set
            jneut =  get_nuclide_index_in_set("neut", set)
            jh1 =  get_nuclide_index_in_set("h1", set)
            jh2 =  get_nuclide_index_in_set("h2", set)
            jh3 =  get_nuclide_index_in_set("h3", set)
            jhe3 =  get_nuclide_index_in_set("he3", set)
            jhe4 =  get_nuclide_index_in_set("he4", set)
            jli6 =  get_nuclide_index_in_set("li6", set)
            jli7 =  get_nuclide_index_in_set("li7", set)
            jbe7 =  get_nuclide_index_in_set("be7", set)
            jbe8 =  get_nuclide_index_in_set("be8", set)
            jbe9 =  get_nuclide_index_in_set("be9", set)
            jbe10 =  get_nuclide_index_in_set("be10", set)
            jb8 =  get_nuclide_index_in_set("b8", set)
            jb10 =  get_nuclide_index_in_set("b10", set)
            jb11 =  get_nuclide_index_in_set("b11", set)
            jc12 =  get_nuclide_index_in_set("c12", set)
            jc13 =  get_nuclide_index_in_set("c13", set)
            jc14 =  get_nuclide_index_in_set("c14", set)
            jn13 =  get_nuclide_index_in_set("n13", set)
            jn14 =  get_nuclide_index_in_set("n14", set)
            jn15 =  get_nuclide_index_in_set("n15", set)
            jn16 =  get_nuclide_index_in_set("n16", set)
            jo14 =  get_nuclide_index_in_set("o14", set)
            jo15 =  get_nuclide_index_in_set("o15", set)
            jo16 =  get_nuclide_index_in_set("o16", set)
            jo17 =  get_nuclide_index_in_set("o17", set)
            jo18 =  get_nuclide_index_in_set("o18", set)
            jo19 =  get_nuclide_index_in_set("o19", set)
            jf17 =  get_nuclide_index_in_set("f17", set)
            jf18 =  get_nuclide_index_in_set("f18", set)
            jf19 =  get_nuclide_index_in_set("f19", set)
            jf20 =  get_nuclide_index_in_set("f20", set)
            jne18 =  get_nuclide_index_in_set("ne18", set)
            jne19 =  get_nuclide_index_in_set("ne19", set)
            jne20 =  get_nuclide_index_in_set("ne20", set)
            jne21 =  get_nuclide_index_in_set("ne21", set)
            jne22 =  get_nuclide_index_in_set("ne22", set)
            jne23 =  get_nuclide_index_in_set("ne23", set)
            jna21 =  get_nuclide_index_in_set("na21", set)
            jna22 =  get_nuclide_index_in_set("na22", set)
            jna23 =  get_nuclide_index_in_set("na23", set)
            jna24 =  get_nuclide_index_in_set("na24", set)
            jmg22 =  get_nuclide_index_in_set("mg22", set)
            jmg23 =  get_nuclide_index_in_set("mg23", set)
            jmg24 =  get_nuclide_index_in_set("mg24", set)
            jmg25 =  get_nuclide_index_in_set("mg25", set)
            jmg26 =  get_nuclide_index_in_set("mg26", set)
            jmg27 =  get_nuclide_index_in_set("mg27", set)
            jal25 =  get_nuclide_index_in_set("al25", set)
            jal26 =  get_nuclide_index_in_set("al26", set)
            jal27 =  get_nuclide_index_in_set("al27", set)
            jal28 =  get_nuclide_index_in_set("al28", set)
            jsi27 =  get_nuclide_index_in_set("si27", set)
            jsi28 =  get_nuclide_index_in_set("si28", set)
            jsi29 =  get_nuclide_index_in_set("si29", set)
            jsi30 =  get_nuclide_index_in_set("si30", set)
            jsi31 =  get_nuclide_index_in_set("si31", set)
            jsi32 =  get_nuclide_index_in_set("si32", set)
            jp29 =  get_nuclide_index_in_set("p29", set)
            jp30 =  get_nuclide_index_in_set("p30", set)
            jp31 =  get_nuclide_index_in_set("p31", set)
            jp32 =  get_nuclide_index_in_set("p32", set)
            js31 =  get_nuclide_index_in_set("s31", set)
            js32 =  get_nuclide_index_in_set("s32", set)
       	end subroutine set_index_in_set
      	
      
      end subroutine test_jina_net_get
      
      
      
      end module test_jina_support


      
      program test_jina
      use jina_lib
      use chem_lib
      use rates_lib
      use weak_lib
      use test_jina_support
   	
   	integer :: ierr
   	character(len=4096) :: data_dir
      data_dir = '../../data'   	

   	call const_init

   	ierr = 0
   	call chem_init(data_dir, ierr)
   	if (ierr /= 0) then
   	   write(*,*) 'chem_init failed'
   	   stop 1
   	end if
   	
   	call rates_init(data_dir, 'rates', ierr)
   	if (ierr /= 0) then
   	   write(*,*) 'rates_init failed'
   	   stop 1
   	end if
   	
      call weak_init(data_dir, ierr)
      if (ierr /= 0) then
         write(*,*) 'weak_init failed'
         stop 1
      end if
   	
   	call jina_init(data_dir, ierr)
   	if (ierr /= 0) then
   	   write(*,*) 'jina_init failed'
   	   stop 1
   	end if
   	
      !call test_jina_net_get
   	
   	!call test_burn(10)
   	!call test_burn(25)
!~   	call test_burn(100)  ! use this for testing
   	!call test_burn(200)
   	!call test_burn(430)
   	!call test_burn(686)
   	!call test_burn(690)
   	stop

      end program test_jina
      
