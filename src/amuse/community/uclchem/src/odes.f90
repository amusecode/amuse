totalSwap=RATE(1)*Y(19)*bulkLayersReciprocal+RATE(2)*Y(73)&
    &*bulkLayersReciprocal+RATE(3)*Y(78)*bulkLayersReciprocal+RATE(4)*Y(85)&
    &*bulkLayersReciprocal+RATE(5)*Y(96)*bulkLayersReciprocal+RATE(6)*Y(114)&
    &*bulkLayersReciprocal+RATE(7)*Y(127)*bulkLayersReciprocal+RATE(8)*Y(202)&
    &*bulkLayersReciprocal+RATE(9)*Y(203)*bulkLayersReciprocal+RATE(10)*Y(286&
    &)*bulkLayersReciprocal+RATE(11)*Y(283)*bulkLayersReciprocal+RATE(12)&
    &*Y(311)*bulkLayersReciprocal+RATE(13)*Y(23)*bulkLayersReciprocal+RATE(14&
    &)*Y(30)*bulkLayersReciprocal+RATE(15)*Y(223)*bulkLayersReciprocal&
    &+RATE(16)*Y(154)*bulkLayersReciprocal+RATE(17)*Y(38)&
    &*bulkLayersReciprocal+RATE(18)*Y(210)*bulkLayersReciprocal+RATE(19)&
    &*Y(242)*bulkLayersReciprocal+RATE(20)*Y(215)*bulkLayersReciprocal&
    &+RATE(21)*Y(224)*bulkLayersReciprocal+RATE(22)*Y(169)&
    &*bulkLayersReciprocal+RATE(23)*Y(49)*bulkLayersReciprocal+RATE(24)*Y(190&
    &)*bulkLayersReciprocal+RATE(25)*Y(86)*bulkLayersReciprocal+RATE(26)&
    &*Y(111)*bulkLayersReciprocal+RATE(27)*Y(241)*bulkLayersReciprocal&
    &+RATE(28)*Y(240)*bulkLayersReciprocal+RATE(29)*Y(5)*bulkLayersReciprocal&
    &+RATE(30)*Y(9)*bulkLayersReciprocal+RATE(31)*Y(113)*bulkLayersReciprocal&
    &+RATE(32)*Y(139)*bulkLayersReciprocal+RATE(33)*Y(269)&
    &*bulkLayersReciprocal+RATE(34)*Y(65)*bulkLayersReciprocal+RATE(35)*Y(185&
    &)*bulkLayersReciprocal+RATE(36)*Y(332)*bulkLayersReciprocal+RATE(37)&
    &*Y(272)*bulkLayersReciprocal+RATE(38)*Y(152)*bulkLayersReciprocal&
    &+RATE(39)*Y(289)*bulkLayersReciprocal+RATE(40)*Y(195)&
    &*bulkLayersReciprocal+RATE(41)*Y(94)*bulkLayersReciprocal+RATE(42)*Y(126&
    &)*bulkLayersReciprocal+RATE(43)*Y(256)*bulkLayersReciprocal+RATE(44)&
    &*Y(273)*bulkLayersReciprocal+RATE(45)*Y(253)*bulkLayersReciprocal&
    &+RATE(46)*Y(14)*bulkLayersReciprocal+RATE(47)*Y(95)*bulkLayersReciprocal&
    &+RATE(48)*Y(228)*bulkLayersReciprocal+RATE(49)*Y(153)&
    &*bulkLayersReciprocal+RATE(50)*Y(181)*bulkLayersReciprocal+RATE(51)&
    &*Y(328)*bulkLayersReciprocal+RATE(52)*Y(254)*bulkLayersReciprocal&
    &+RATE(53)*Y(74)*bulkLayersReciprocal+RATE(54)*Y(31)*bulkLayersReciprocal&
    &+RATE(55)*Y(112)*bulkLayersReciprocal+RATE(56)*Y(296)&
    &*bulkLayersReciprocal+RATE(57)*Y(39)*bulkLayersReciprocal+RATE(58)*Y(50)&
    &*bulkLayersReciprocal+RATE(59)*Y(255)*bulkLayersReciprocal+RATE(60)*Y(59&
    &)*bulkLayersReciprocal+RATE(61)*Y(138)*bulkLayersReciprocal+RATE(62)&
    &*Y(270)*bulkLayersReciprocal+RATE(63)*Y(271)*bulkLayersReciprocal&
    &+RATE(64)*Y(51)*bulkLayersReciprocal+RATE(65)*Y(168)&
    &*bulkLayersReciprocal+RATE(66)*Y(180)*bulkLayersReciprocal+RATE(67)&
    &*Y(225)*bulkLayersReciprocal+RATE(68)*Y(305)*bulkLayersReciprocal&
    &+RATE(69)*Y(60)*bulkLayersReciprocal+RATE(70)*Y(171)&
    &*bulkLayersReciprocal+RATE(71)*Y(323)*bulkLayersReciprocal+RATE(72)&
    &*Y(115)*bulkLayersReciprocal+RATE(73)*Y(211)*bulkLayersReciprocal&
    &+RATE(74)*Y(297)*bulkLayersReciprocal+RATE(75)*Y(321)&
    &*bulkLayersReciprocal+RATE(76)*Y(128)*bulkLayersReciprocal+RATE(77)&
    &*Y(140)*bulkLayersReciprocal+RATE(78)*Y(155)*bulkLayersReciprocal&
    &+RATE(79)*Y(170)*bulkLayersReciprocal+RATE(80)*Y(239)&
    &*bulkLayersReciprocal+RATE(81)*Y(304)*bulkLayersReciprocal+RATE(82)&
    &*Y(279)*bulkLayersReciprocal+RATE(83)*Y(322)*bulkLayersReciprocal

    PROD = RATE(776)*Y(229)*Y(3)+RATE(828)*Y(239)*Y(5)&
    &*bulkLayersReciprocal
    YDOT(1) = PROD
    LOSS = RATE(102)*Y(2)+RATE(143)*Y(2)+RATE(453)*D*Y(2)*Y(3&
    &)/safeMantle+RATE(504)*D*Y(2)+RATE(620)*D*Y(2)+RATE(620)*D*Y(2)&
    &+RATE(1742)*D*Y(187)*Y(2)+RATE(1798)*D*Y(2)*Y(20)+RATE(1799)*D*Y(2)*Y(6)&
    &+RATE(1800)*D*Y(2)*Y(62)+RATE(1801)*D*Y(2)*Y(161)+RATE(1802)*D*Y(2)*Y(56&
    &)+RATE(1803)*D*Y(2)*Y(83)+RATE(1804)*D*Y(2)*Y(100)+RATE(1805)*D*Y(2)*Y(8&
    &)+RATE(1806)*D*Y(2)*Y(91)+RATE(1807)*D*Y(2)*Y(13)+RATE(1808)*D*Y(2)*Y(48&
    &)+RATE(1809)*D*Y(2)*Y(290)+RATE(1810)*D*Y(2)*Y(21)+RATE(1811)*D*Y(2)&
    &*Y(25)+RATE(1812)*D*Y(2)*Y(33)+RATE(1813)*D*Y(2)*Y(42)+RATE(1814)*D*Y(2)&
    &*Y(53)+RATE(1815)*D*Y(2)*Y(184)+RATE(1816)*D*Y(2)*Y(189)+RATE(1817)*D&
    &*Y(2)*Y(174)+RATE(1818)*D*Y(2)*Y(15)+RATE(1819)*D*Y(2)*Y(319)+RATE(1820)&
    &*D*Y(2)*Y(122)+RATE(1821)*D*Y(2)*Y(303)+RATE(1822)*D*Y(2)*Y(67)&
    &+RATE(1823)*D*Y(2)*Y(80)+RATE(1824)*D*Y(2)*Y(89)+RATE(1825)*D*Y(2)*Y(24)&
    &+RATE(1826)*D*Y(2)*Y(218)+RATE(1827)*D*Y(2)*Y(32)+RATE(1828)*D*Y(2)*Y(41&
    &)+RATE(1829)*D*Y(2)*Y(20)+RATE(1830)*D*Y(2)*Y(232)+RATE(1831)*D*Y(2)&
    &*Y(99)+RATE(1832)*D*Y(2)*Y(101)+RATE(1833)*D*Y(2)*Y(131)+RATE(1834)*D&
    &*Y(2)*Y(62)+RATE(1835)*D*Y(2)*Y(183)+RATE(1836)*D*Y(2)*Y(90)+RATE(1837)&
    &*D*Y(2)*Y(116)+RATE(1838)*D*Y(2)*Y(116)+RATE(1839)*D*Y(2)*Y(244)&
    &+RATE(1840)*D*Y(2)*Y(92)+RATE(1841)*D*Y(2)*Y(144)+RATE(1842)*D*Y(2)&
    &*Y(144)+RATE(1843)*D*Y(2)*Y(144)+RATE(1844)*D*Y(2)*Y(173)+RATE(1845)*D&
    &*Y(2)*Y(291)+RATE(1846)*D*Y(2)*Y(43)+RATE(1847)*D*Y(2)*Y(54)+RATE(1848)&
    &*D*Y(2)*Y(35)+RATE(1849)*D*Y(2)*Y(264)+RATE(1850)*D*Y(2)*Y(133)&
    &+RATE(1851)*D*Y(2)*Y(133)+RATE(1852)*D*Y(2)*Y(265)+RATE(1853)*D*Y(2)&
    &*Y(265)+RATE(1854)*D*Y(2)*Y(161)+RATE(1855)*D*Y(2)*Y(175)+RATE(1856)*D&
    &*Y(2)*Y(175)+RATE(1857)*D*Y(2)*Y(175)+RATE(1858)*D*Y(2)*Y(221)+RATE(1859&
    &)*D*Y(2)*Y(221)+RATE(1860)*D*Y(2)*Y(221)+RATE(1861)*D*Y(2)*Y(300)&
    &+RATE(1862)*D*Y(2)*Y(56)+RATE(1863)*D*Y(2)*Y(316)+RATE(1864)*D*Y(2)&
    &*Y(277)+RATE(1865)*D*Y(2)*Y(277)+RATE(1866)*D*Y(2)*Y(18)+RATE(1867)*D&
    &*Y(2)*Y(16)+RATE(1868)*D*Y(2)*Y(46)+RATE(1869)*D*Y(2)*Y(56)+RATE(1870)*D&
    &*Y(2)*Y(106)+RATE(1942)*D*Y(4)*Y(2)+RATE(2239)*D*Y(193)*Y(2)
    PROD = RATE(92)*Y(182)+RATE(94)*Y(268)+RATE(103)*Y(6)+RATE(105)*Y(6)&
    &+RATE(105)*Y(6)+RATE(112)*Y(75)+RATE(115)*Y(80)+RATE(116)*Y(89)+RATE(123&
    &)*Y(20)+RATE(124)*Y(21)+RATE(126)*Y(24)+RATE(128)*Y(32)+RATE(144)*Y(101)&
    &+RATE(147)*Y(62)+RATE(153)*Y(193)+RATE(154)*Y(90)+RATE(155)*Y(116)&
    &+RATE(160)*Y(92)+RATE(162)*Y(144)+RATE(163)*Y(173)+RATE(169)*Y(35)&
    &+RATE(172)*Y(43)+RATE(173)*Y(54)+RATE(183)*Y(175)+RATE(187)*Y(56)&
    &+RATE(194)*Y(121)+RATE(195)*Y(135)+RATE(196)*Y(146)+RATE(231)*Y(3&
    &)/safeMantle+RATE(314)*D*Y(3)/safeMantle*Y(2)+RATE(397)*Y(3)/safeMantle&
    &+RATE(471)*D*Y(213)+RATE(490)*D*Y(172)+RATE(493)*D*Y(53)+RATE(508)*D&
    &*Y(196)+RATE(514)*D*Y(160)+RATE(525)*D*Y(274)+RATE(526)*D*Y(66)+RATE(527&
    &)*D*Y(189)+RATE(533)*D*Y(102)+RATE(536)*D*Y(243)+RATE(548)*D*Y(275)&
    &+RATE(550)*D*Y(306)+RATE(555)*D*Y(307)+RATE(556)*D*Y(281)+RATE(557)*D&
    &*Y(326)+RATE(564)*D*Y(120)+RATE(573)*D*Y(64)+RATE(610)*D*Y(177)+RATE(661&
    &)*Y(45)*Y(130)+RATE(713)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(765)&
    &*Y(45)*Y(130)+RATE(817)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(834)*Y(75&
    &)+RATE(836)*Y(76)+RATE(838)*Y(80)+RATE(839)*Y(89)+RATE(846)*Y(20)&
    &+RATE(850)*Y(24)+RATE(852)*Y(25)+RATE(855)*Y(32)+RATE(859)*Y(33)&
    &+RATE(864)*Y(159)+RATE(867)*Y(41)+RATE(869)*Y(41)+RATE(871)*Y(42)&
    &+RATE(880)*Y(6)+RATE(880)*Y(6)+RATE(881)*Y(8)+RATE(882)*Y(101)+RATE(884)&
    &*Y(131)+RATE(884)*Y(131)+RATE(886)*Y(131)+RATE(889)*Y(62)+RATE(890)*Y(63&
    &)+RATE(892)*Y(183)+RATE(896)*Y(263)+RATE(896)*Y(263)+RATE(897)*Y(10)&
    &+RATE(900)*Y(193)+RATE(902)*Y(90)+RATE(903)*Y(116)+RATE(905)*Y(117)&
    &+RATE(908)*Y(92)+RATE(910)*Y(144)+RATE(911)*Y(173)+RATE(912)*Y(174)&
    &+RATE(920)*Y(35)+RATE(924)*Y(43)+RATE(925)*Y(54)+RATE(935)*Y(175)&
    &+RATE(940)*Y(56)+RATE(942)*Y(57)+RATE(951)*Y(121)+RATE(952)*Y(122)&
    &+RATE(954)*Y(135)+RATE(955)*Y(146)+RATE(959)*Y(166)+RATE(960)*Y(166)&
    &+RATE(1079)*Y(3)+RATE(1162)*Y(5)+RATE(1222)*D*Y(16)*Y(76)+RATE(1223)*D&
    &*Y(16)*Y(21)+RATE(1224)*D*Y(16)*Y(25)+RATE(1228)*D*Y(16)*Y(184)&
    &+RATE(1234)*D*Y(16)*Y(174)+RATE(1240)*D*Y(16)*Y(122)+RATE(1242)*D*Y(16)&
    &*Y(89)+RATE(1243)*D*Y(16)*Y(123)+RATE(1245)*D*Y(16)*Y(200)+RATE(1246)*D&
    &*Y(16)*Y(24)+RATE(1248)*D*Y(16)*Y(32)+RATE(1249)*D*Y(16)*Y(20)+RATE(1255&
    &)*D*Y(16)*Y(173)+RATE(1259)*D*Y(16)*Y(43)+RATE(1260)*D*Y(16)*Y(43)&
    &+RATE(1262)*D*Y(16)*Y(35)+RATE(1271)*D*Y(16)*Y(56)+RATE(1277)*D*Y(16)&
    &*Y(121)+RATE(1303)*D*Y(18)*Y(75)+RATE(1304)*D*Y(18)*Y(24)+RATE(1306)*D&
    &*Y(18)*Y(32)+RATE(1311)*D*Y(18)*Y(20)+RATE(1316)*D*Y(18)*Y(62)+RATE(1317&
    &)*D*Y(18)*Y(62)+RATE(1318)*D*Y(18)*Y(183)+RATE(1321)*D*Y(18)*Y(288)&
    &+RATE(1323)*D*Y(18)*Y(92)+RATE(1324)*D*Y(18)*Y(173)+RATE(1325)*D*Y(18)&
    &*Y(43)+RATE(1327)*D*Y(18)*Y(35)+RATE(1333)*D*Y(18)*Y(56)+RATE(1340)*D&
    &*Y(18)*Y(121)+RATE(1361)*D*Y(67)*Y(80)+RATE(1362)*D*Y(67)*Y(90)&
    &+RATE(1383)*D*Y(75)*Y(106)+RATE(1384)*D*Y(75)*Y(90)+RATE(1385)*D*Y(75)&
    &*Y(92)+RATE(1391)*D*Y(76)*Y(333)+RATE(1400)*D*Y(80)*Y(133)+RATE(1406)*D&
    &*Y(81)*Y(333)+RATE(1406)*D*Y(81)*Y(333)+RATE(1407)*D*Y(81)*Y(333)&
    &+RATE(1426)*D*Y(204)*Y(333)+RATE(1428)*D*Y(213)*Y(333)+RATE(1443)*D*Y(20&
    &)*Y(68)+RATE(1458)*D*Y(20)*Y(28)+RATE(1463)*D*Y(20)*Y(48)+RATE(1467)*D&
    &*Y(20)*Y(165)+RATE(1468)*D*Y(20)*Y(106)+RATE(1471)*D*Y(20)*Y(80)&
    &+RATE(1472)*D*Y(20)*Y(109)+RATE(1473)*D*Y(20)*Y(159)+RATE(1474)*D*Y(20)&
    &*Y(41)+RATE(1480)*D*Y(20)*Y(27)+RATE(1484)*D*Y(20)*Y(133)+RATE(1485)*D&
    &*Y(20)*Y(161)+RATE(1486)*D*Y(20)*Y(161)+RATE(1491)*D*Y(20)*Y(46)&
    &+RATE(1493)*D*Y(20)*Y(300)+RATE(1494)*D*Y(20)*Y(56)+RATE(1495)*D*Y(20)&
    &*Y(163)+RATE(1498)*D*Y(20)*Y(277)+RATE(1505)*D*Y(21)*Y(333)+RATE(1506)*D&
    &*Y(21)*Y(67)+RATE(1512)*D*Y(21)*Y(41)+RATE(1514)*D*Y(21)*Y(82)+RATE(1519&
    &)*D*Y(21)*Y(62)+RATE(1525)*D*Y(21)*Y(90)+RATE(1529)*D*Y(21)*Y(27)&
    &+RATE(1536)*D*Y(21)*Y(46)+RATE(1540)*D*Y(21)*Y(163)+RATE(1570)*D*Y(24)&
    &*Y(165)+RATE(1573)*D*Y(24)*Y(24)+RATE(1573)*D*Y(24)*Y(24)+RATE(1574)*D&
    &*Y(24)*Y(24)+RATE(1585)*D*Y(24)*Y(133)+RATE(1587)*D*Y(24)*Y(161)&
    &+RATE(1587)*D*Y(24)*Y(161)+RATE(1592)*D*Y(24)*Y(46)+RATE(1592)*D*Y(24)&
    &*Y(46)+RATE(1593)*D*Y(24)*Y(46)+RATE(1595)*D*Y(24)*Y(56)+RATE(1599)*D&
    &*Y(24)*Y(163)+RATE(1602)*D*Y(25)*Y(333)+RATE(1602)*D*Y(25)*Y(333)&
    &+RATE(1603)*D*Y(25)*Y(333)+RATE(1606)*D*Y(25)*Y(62)+RATE(1607)*D*Y(25)&
    &*Y(183)+RATE(1609)*D*Y(25)*Y(183)+RATE(1613)*D*Y(25)*Y(46)+RATE(1616)*D&
    &*Y(25)*Y(163)+RATE(1617)*D*Y(32)*Y(165)+RATE(1620)*D*Y(32)*Y(32)&
    &+RATE(1636)*D*Y(32)*Y(46)+RATE(1637)*D*Y(32)*Y(46)+RATE(1641)*D*Y(32)&
    &*Y(163)+RATE(1646)*D*Y(33)*Y(333)+RATE(1648)*D*Y(33)*Y(333)+RATE(1648)*D&
    &*Y(33)*Y(333)+RATE(1659)*D*Y(33)*Y(46)+RATE(1670)*D*Y(220)*Y(333)&
    &+RATE(1673)*D*Y(172)*Y(333)+RATE(1675)*D*Y(172)*Y(333)+RATE(1676)*D&
    &*Y(172)*Y(333)+RATE(1677)*D*Y(172)*Y(333)+RATE(1683)*D*Y(41)*Y(81)&
    &+RATE(1693)*D*Y(41)*Y(104)+RATE(1698)*D*Y(41)*Y(165)+RATE(1699)*D*Y(41)&
    &*Y(165)+RATE(1710)*D*Y(42)*Y(333)+RATE(1710)*D*Y(42)*Y(333)+RATE(1711)*D&
    &*Y(42)*Y(333)+RATE(1721)*D*Y(53)*Y(333)+RATE(1723)*D*Y(53)*Y(333)&
    &+RATE(1723)*D*Y(53)*Y(333)+RATE(1724)*D*Y(53)*Y(333)+RATE(1737)*D*Y(53)&
    &*Y(69)+RATE(1740)*D*Y(186)*Y(4)+RATE(1741)*D*Y(186)*Y(6)+RATE(1747)*D&
    &*Y(82)*Y(80)+RATE(1751)*D*Y(82)*Y(90)+RATE(1753)*D*Y(82)*Y(92)+RATE(1771&
    &)*D*Y(83)*Y(90)+RATE(1798)*D*Y(2)*Y(20)+RATE(1798)*D*Y(2)*Y(20)&
    &+RATE(1799)*D*Y(2)*Y(6)+RATE(1799)*D*Y(2)*Y(6)+RATE(1799)*D*Y(2)*Y(6)&
    &+RATE(1800)*D*Y(2)*Y(62)+RATE(1800)*D*Y(2)*Y(62)+RATE(1801)*D*Y(2)*Y(161&
    &)+RATE(1802)*D*Y(2)*Y(56)+RATE(1802)*D*Y(2)*Y(56)+RATE(1840)*D*Y(2)*Y(92&
    &)+RATE(1872)*D*Y(4)*Y(67)+RATE(1873)*D*Y(4)*Y(80)+RATE(1874)*D*Y(4)*Y(75&
    &)+RATE(1875)*D*Y(4)*Y(197)+RATE(1876)*D*Y(4)*Y(24)+RATE(1877)*D*Y(4)&
    &*Y(32)+RATE(1878)*D*Y(4)*Y(41)+RATE(1879)*D*Y(4)*Y(20)+RATE(1880)*D*Y(4)&
    &*Y(233)+RATE(1881)*D*Y(4)*Y(131)+RATE(1882)*D*Y(4)*Y(261)+RATE(1883)*D&
    &*Y(4)*Y(62)+RATE(1884)*D*Y(4)*Y(330)+RATE(1885)*D*Y(4)*Y(183)+RATE(1886)&
    &*D*Y(4)*Y(90)+RATE(1887)*D*Y(4)*Y(116)+RATE(1888)*D*Y(4)*Y(193)&
    &+RATE(1889)*D*Y(4)*Y(324)+RATE(1890)*D*Y(4)*Y(173)+RATE(1891)*D*Y(4)&
    &*Y(69)+RATE(1892)*D*Y(4)*Y(43)+RATE(1893)*D*Y(4)*Y(54)+RATE(1894)*D*Y(4)&
    &*Y(35)+RATE(1895)*D*Y(4)*Y(133)+RATE(1896)*D*Y(4)*Y(265)+RATE(1897)*D&
    &*Y(4)*Y(161)+RATE(1898)*D*Y(4)*Y(46)+RATE(1899)*D*Y(4)*Y(300)+RATE(1900)&
    &*D*Y(4)*Y(56)+RATE(1901)*D*Y(4)*Y(316)+RATE(1902)*D*Y(4)*Y(163)&
    &+RATE(1903)*D*Y(4)*Y(320)+RATE(1904)*D*Y(4)*Y(277)+RATE(1905)*D*Y(4)&
    &*Y(105)+RATE(1906)*D*Y(4)*Y(292)+RATE(1907)*D*Y(4)*Y(314)+RATE(1908)*D&
    &*Y(4)*Y(207)+RATE(1909)*D*Y(4)*Y(135)+RATE(1910)*D*Y(4)*Y(146)+RATE(1911&
    &)*D*Y(4)*Y(166)+RATE(1912)*D*Y(4)*Y(121)+RATE(1913)*D*Y(4)*Y(235)&
    &+RATE(1914)*D*Y(4)*Y(302)+RATE(1916)*D*Y(4)*Y(109)+RATE(1925)*D*Y(4)&
    &*Y(131)+RATE(1928)*D*Y(4)*Y(183)+RATE(1944)*D*Y(4)*Y(333)+RATE(1945)*D&
    &*Y(6)*Y(20)+RATE(1946)*D*Y(6)*Y(6)+RATE(1946)*D*Y(6)*Y(6)+RATE(1947)*D&
    &*Y(6)*Y(62)+RATE(1950)*D*Y(6)*Y(56)+RATE(1951)*D*Y(6)*Y(333)+RATE(1951)&
    &*D*Y(6)*Y(333)+RATE(1953)*D*Y(6)*Y(18)+RATE(1954)*D*Y(6)*Y(68)+RATE(1955&
    &)*D*Y(6)*Y(76)+RATE(1956)*D*Y(6)*Y(21)+RATE(1957)*D*Y(6)*Y(25)+RATE(1958&
    &)*D*Y(6)*Y(42)+RATE(1959)*D*Y(6)*Y(83)+RATE(1960)*D*Y(6)*Y(100)&
    &+RATE(1961)*D*Y(6)*Y(100)+RATE(1962)*D*Y(6)*Y(234)+RATE(1963)*D*Y(6)&
    &*Y(187)+RATE(1964)*D*Y(6)*Y(63)+RATE(1965)*D*Y(6)*Y(184)+RATE(1966)*D&
    &*Y(6)*Y(91)+RATE(1967)*D*Y(6)*Y(194)+RATE(1968)*D*Y(6)*Y(174)+RATE(1969)&
    &*D*Y(6)*Y(13)+RATE(1971)*D*Y(6)*Y(28)+RATE(1972)*D*Y(6)*Y(104)+RATE(1974&
    &)*D*Y(6)*Y(36)+RATE(1975)*D*Y(6)*Y(44)+RATE(1976)*D*Y(6)*Y(55)+RATE(1977&
    &)*D*Y(6)*Y(48)+RATE(1979)*D*Y(6)*Y(57)+RATE(1980)*D*Y(6)*Y(165)&
    &+RATE(1981)*D*Y(6)*Y(319)+RATE(1982)*D*Y(6)*Y(167)+RATE(1983)*D*Y(6)&
    &*Y(236)+RATE(1984)*D*Y(6)*Y(75)+RATE(1985)*D*Y(6)*Y(16)+RATE(1986)*D*Y(6&
    &)*Y(24)+RATE(1987)*D*Y(6)*Y(32)+RATE(1988)*D*Y(6)*Y(20)+RATE(1989)*D*Y(6&
    &)*Y(82)+RATE(1990)*D*Y(6)*Y(173)+RATE(1991)*D*Y(6)*Y(27)+RATE(1992)*D&
    &*Y(6)*Y(43)+RATE(1993)*D*Y(6)*Y(35)+RATE(1994)*D*Y(6)*Y(161)+RATE(1996)&
    &*D*Y(6)*Y(46)+RATE(1997)*D*Y(6)*Y(56)+RATE(1998)*D*Y(6)*Y(163)+RATE(2027&
    &)*D*Y(8)*Y(333)+RATE(2027)*D*Y(8)*Y(333)+RATE(2028)*D*Y(8)*Y(67)&
    &+RATE(2030)*D*Y(8)*Y(75)+RATE(2031)*D*Y(8)*Y(16)+RATE(2032)*D*Y(8)*Y(24)&
    &+RATE(2033)*D*Y(8)*Y(41)+RATE(2034)*D*Y(8)*Y(41)+RATE(2035)*D*Y(8)*Y(20)&
    &+RATE(2036)*D*Y(8)*Y(82)+RATE(2037)*D*Y(8)*Y(232)+RATE(2038)*D*Y(8)*Y(99&
    &)+RATE(2039)*D*Y(8)*Y(6)+RATE(2040)*D*Y(8)*Y(131)+RATE(2041)*D*Y(8)*Y(62&
    &)+RATE(2042)*D*Y(8)*Y(183)+RATE(2045)*D*Y(8)*Y(11)+RATE(2046)*D*Y(8)&
    &*Y(103)+RATE(2047)*D*Y(8)*Y(27)+RATE(2048)*D*Y(8)*Y(35)+RATE(2049)*D*Y(8&
    &)*Y(133)+RATE(2050)*D*Y(8)*Y(161)+RATE(2051)*D*Y(8)*Y(46)+RATE(2052)*D&
    &*Y(8)*Y(56)+RATE(2053)*D*Y(196)*Y(333)+RATE(2053)*D*Y(196)*Y(333)&
    &+RATE(2054)*D*Y(196)*Y(333)+RATE(2058)*D*Y(131)*Y(162)+RATE(2065)*D&
    &*Y(132)*Y(333)+RATE(2065)*D*Y(132)*Y(333)+RATE(2066)*D*Y(132)*Y(333)&
    &+RATE(2072)*D*Y(262)*Y(333)+RATE(2072)*D*Y(262)*Y(333)+RATE(2073)*D&
    &*Y(262)*Y(333)+RATE(2075)*D*Y(160)*Y(333)+RATE(2103)*D*Y(62)*Y(106)&
    &+RATE(2120)*D*Y(63)*Y(333)+RATE(2120)*D*Y(63)*Y(333)+RATE(2121)*D*Y(63)&
    &*Y(333)+RATE(2134)*D*Y(63)*Y(163)+RATE(2141)*D*Y(184)*Y(333)+RATE(2142)&
    &*D*Y(184)*Y(333)+RATE(2142)*D*Y(184)*Y(333)+RATE(2145)*D*Y(329)*Y(333)&
    &+RATE(2147)*D*Y(10)*Y(333)+RATE(2148)*D*Y(10)*Y(333)+RATE(2148)*D*Y(10)&
    &*Y(333)+RATE(2148)*D*Y(10)*Y(333)+RATE(2183)*D*Y(10)*Y(69)+RATE(2192)*D&
    &*Y(10)*Y(46)+RATE(2209)*D*Y(142)*Y(333)+RATE(2210)*D*Y(142)*Y(333)&
    &+RATE(2211)*D*Y(142)*Y(333)+RATE(2211)*D*Y(142)*Y(333)+RATE(2214)*D&
    &*Y(274)*Y(333)+RATE(2215)*D*Y(274)*Y(333)+RATE(2216)*D*Y(66)*Y(333)&
    &+RATE(2217)*D*Y(66)*Y(333)+RATE(2219)*D*Y(66)*Y(333)+RATE(2219)*D*Y(66)&
    &*Y(333)+RATE(2235)*D*Y(189)*Y(333)+RATE(2237)*D*Y(189)*Y(333)+RATE(2237)&
    &*D*Y(189)*Y(333)+RATE(2238)*D*Y(189)*Y(333)+RATE(2240)*D*Y(194)*Y(333)&
    &+RATE(2244)*D*Y(90)*Y(192)+RATE(2258)*D*Y(91)*Y(333)+RATE(2267)*D*Y(102)&
    &*Y(333)+RATE(2267)*D*Y(102)*Y(333)+RATE(2268)*D*Y(102)*Y(333)+RATE(2269)&
    &*D*Y(102)*Y(333)+RATE(2294)*D*Y(117)*Y(333)+RATE(2315)*D*Y(243)*Y(333)&
    &+RATE(2316)*D*Y(243)*Y(333)+RATE(2319)*D*Y(245)*Y(333)+RATE(2335)*D*Y(13&
    &)*Y(80)+RATE(2338)*D*Y(13)*Y(89)+RATE(2339)*D*Y(13)*Y(109)+RATE(2342)*D&
    &*Y(13)*Y(75)+RATE(2349)*D*Y(13)*Y(24)+RATE(2359)*D*Y(13)*Y(41)+RATE(2361&
    &)*D*Y(13)*Y(41)+RATE(2363)*D*Y(13)*Y(20)+RATE(2374)*D*Y(13)*Y(131)&
    &+RATE(2379)*D*Y(13)*Y(62)+RATE(2382)*D*Y(13)*Y(330)+RATE(2383)*D*Y(13)&
    &*Y(183)+RATE(2385)*D*Y(13)*Y(263)+RATE(2388)*D*Y(13)*Y(90)+RATE(2390)*D&
    &*Y(13)*Y(90)+RATE(2392)*D*Y(13)*Y(116)+RATE(2395)*D*Y(13)*Y(244)&
    &+RATE(2397)*D*Y(13)*Y(193)+RATE(2398)*D*Y(13)*Y(92)+RATE(2399)*D*Y(13)&
    &*Y(92)+RATE(2401)*D*Y(13)*Y(144)+RATE(2404)*D*Y(13)*Y(324)+RATE(2405)*D&
    &*Y(13)*Y(173)+RATE(2409)*D*Y(13)*Y(43)+RATE(2411)*D*Y(13)*Y(54)&
    &+RATE(2412)*D*Y(13)*Y(35)+RATE(2424)*D*Y(13)*Y(56)+RATE(2435)*D*Y(13)&
    &*Y(135)+RATE(2437)*D*Y(13)*Y(146)+RATE(2439)*D*Y(13)*Y(166)+RATE(2440)*D&
    &*Y(13)*Y(121)+RATE(2446)*D*Y(15)*Y(333)+RATE(2456)*D*Y(145)*Y(333)&
    &+RATE(2459)*D*Y(275)*Y(333)+RATE(2460)*D*Y(119)*Y(333)+RATE(2462)*D&
    &*Y(306)*Y(333)+RATE(2464)*D*Y(174)*Y(333)+RATE(2468)*D*Y(325)*Y(333)&
    &+RATE(2470)*D*Y(307)*Y(333)+RATE(2471)*D*Y(281)*Y(333)+RATE(2472)*D&
    &*Y(326)*Y(333)+RATE(2473)*D*Y(326)*Y(333)+RATE(2491)*D*Y(27)*Y(76)&
    &+RATE(2494)*D*Y(27)*Y(81)+RATE(2496)*D*Y(27)*Y(25)+RATE(2498)*D*Y(27)&
    &*Y(63)+RATE(2501)*D*Y(27)*Y(174)+RATE(2502)*D*Y(27)*Y(36)+RATE(2503)*D&
    &*Y(27)*Y(44)+RATE(2505)*D*Y(27)*Y(57)+RATE(2515)*D*Y(27)*Y(75)+RATE(2517&
    &)*D*Y(27)*Y(200)+RATE(2519)*D*Y(27)*Y(280)+RATE(2521)*D*Y(27)*Y(24)&
    &+RATE(2522)*D*Y(27)*Y(24)+RATE(2524)*D*Y(27)*Y(32)+RATE(2526)*D*Y(27)&
    &*Y(32)+RATE(2526)*D*Y(27)*Y(32)+RATE(2533)*D*Y(27)*Y(116)+RATE(2536)*D&
    &*Y(27)*Y(173)+RATE(2539)*D*Y(27)*Y(35)+RATE(2547)*D*Y(27)*Y(56)&
    &+RATE(2572)*D*Y(28)*Y(159)+RATE(2574)*D*Y(28)*Y(159)+RATE(2575)*D*Y(28)&
    &*Y(159)+RATE(2576)*D*Y(28)*Y(41)+RATE(2577)*D*Y(28)*Y(41)+RATE(2578)*D&
    &*Y(28)*Y(41)+RATE(2578)*D*Y(28)*Y(41)+RATE(2585)*D*Y(28)*Y(183)&
    &+RATE(2591)*D*Y(28)*Y(35)+RATE(2609)*D*Y(104)*Y(131)+RATE(2610)*D*Y(104)&
    &*Y(183)+RATE(2614)*D*Y(120)*Y(333)+RATE(2625)*D*Y(35)*Y(68)+RATE(2637)*D&
    &*Y(35)*Y(48)+RATE(2641)*D*Y(35)*Y(165)+RATE(2647)*D*Y(35)*Y(35)&
    &+RATE(2647)*D*Y(35)*Y(35)+RATE(2650)*D*Y(35)*Y(133)+RATE(2654)*D*Y(35)&
    &*Y(46)+RATE(2657)*D*Y(35)*Y(56)+RATE(2660)*D*Y(35)*Y(163)+RATE(2667)*D&
    &*Y(36)*Y(333)+RATE(2669)*D*Y(36)*Y(67)+RATE(2696)*D*Y(36)*Y(163)&
    &+RATE(2704)*D*Y(43)*Y(68)+RATE(2724)*D*Y(43)*Y(133)+RATE(2727)*D*Y(43)&
    &*Y(131)+RATE(2733)*D*Y(44)*Y(333)+RATE(2733)*D*Y(44)*Y(333)+RATE(2734)*D&
    &*Y(44)*Y(333)+RATE(2753)*D*Y(44)*Y(163)+RATE(2797)*D*Y(55)*Y(333)&
    &+RATE(2798)*D*Y(55)*Y(333)+RATE(2798)*D*Y(55)*Y(333)+RATE(2806)*D*Y(64)&
    &*Y(333)+RATE(2806)*D*Y(64)*Y(333)+RATE(2807)*D*Y(64)*Y(333)+RATE(2838)*D&
    &*Y(46)*Y(245)+RATE(2841)*D*Y(46)*Y(174)+RATE(2844)*D*Y(46)*Y(44)&
    &+RATE(2848)*D*Y(46)*Y(57)+RATE(2850)*D*Y(46)*Y(122)+RATE(2851)*D*Y(46)&
    &*Y(136)+RATE(2856)*D*Y(46)*Y(89)+RATE(2861)*D*Y(46)*Y(123)+RATE(2879)*D&
    &*Y(46)*Y(90)+RATE(2880)*D*Y(46)*Y(116)+RATE(2883)*D*Y(46)*Y(244)&
    &+RATE(2884)*D*Y(46)*Y(144)+RATE(2888)*D*Y(46)*Y(173)+RATE(2890)*D*Y(46)&
    &*Y(43)+RATE(2902)*D*Y(46)*Y(56)+RATE(2911)*D*Y(46)*Y(135)+RATE(2911)*D&
    &*Y(46)*Y(135)+RATE(2912)*D*Y(46)*Y(146)+RATE(2914)*D*Y(46)*Y(121)&
    &+RATE(2951)*D*Y(48)*Y(56)+RATE(2968)*D*Y(162)*Y(80)+RATE(2969)*D*Y(162)&
    &*Y(159)+RATE(2972)*D*Y(176)*Y(333)+RATE(2986)*D*Y(56)*Y(117)+RATE(2991)&
    &*D*Y(56)*Y(165)+RATE(2992)*D*Y(56)*Y(106)+RATE(2994)*D*Y(56)*Y(80)&
    &+RATE(2999)*D*Y(56)*Y(82)+RATE(3000)*D*Y(56)*Y(99)+RATE(3002)*D*Y(56)&
    &*Y(233)+RATE(3004)*D*Y(56)*Y(131)+RATE(3012)*D*Y(56)*Y(133)+RATE(3015)*D&
    &*Y(56)*Y(163)+RATE(3016)*D*Y(56)*Y(277)+RATE(3017)*D*Y(56)*Y(105)&
    &+RATE(3029)*D*Y(57)*Y(333)+RATE(3047)*D*Y(57)*Y(163)+RATE(3054)*D*Y(163)&
    &*Y(189)+RATE(3057)*D*Y(163)*Y(116)+RATE(3059)*D*Y(163)*Y(173)+RATE(3064)&
    &*D*Y(165)*Y(183)+RATE(3090)*D*Y(121)*Y(165)+RATE(3091)*D*Y(122)*Y(333)&
    &+RATE(3093)*D*Y(136)*Y(333)+RATE(3093)*D*Y(136)*Y(333)+RATE(3094)*D&
    &*Y(136)*Y(333)+RATE(3096)*D*Y(136)*Y(163)+RATE(3097)*D*Y(147)*Y(333)&
    &+RATE(3100)*D*Y(167)*Y(333)+RATE(3102)*D*Y(177)*Y(333)+RATE(3105)*D&
    &*Y(247)*Y(333)
    YDOT(2) = PROD-LOSS
    LOSS = RATE(231)*Y(3)/safeMantle+RATE(314)*D*Y(3)/safeMantle*Y(2)&
    &+RATE(397)*Y(3)/safeMantle+RATE(452)*D*Y(2)*Y(3)/safeMantle+RATE(452)*D&
    &*Y(2)*Y(3)/safeMantle+RATE(453)*D*Y(2)*Y(3)/safeMantle+RATE(624)*Y(188)&
    &*Y(3)+RATE(625)*Y(84)*Y(3)+RATE(627)*Y(230)*Y(3)+RATE(628)*Y(3)*Y(3)&
    &+RATE(628)*Y(3)*Y(3)+RATE(629)*Y(3)*Y(17)+RATE(630)*Y(3)*Y(22)+RATE(631)&
    &*Y(3)*Y(26)+RATE(632)*Y(3)*Y(34)+RATE(633)*Y(3)*Y(47)+RATE(634)*Y(3)&
    &*Y(58)+RATE(635)*Y(3)*Y(29)+RATE(636)*Y(3)*Y(37)+RATE(637)*Y(3)*Y(45)&
    &+RATE(638)*Y(3)*Y(97)+RATE(639)*Y(3)*Y(118)+RATE(640)*Y(3)*Y(118)&
    &+RATE(641)*Y(3)*Y(130)+RATE(642)*Y(3)*Y(130)+RATE(643)*Y(3)*Y(130)&
    &+RATE(644)*Y(3)*Y(150)+RATE(645)*Y(3)*Y(143)+RATE(646)*Y(3)*Y(222)&
    &+RATE(647)*Y(3)*Y(226)+RATE(651)*Y(246)*Y(3)+RATE(652)*Y(179)*Y(3)&
    &+RATE(653)*Y(327)*Y(3)+RATE(654)*Y(248)*Y(3)+RATE(666)*Y(164)*Y(3)&
    &+RATE(667)*Y(317)*Y(3)+RATE(668)*Y(110)*Y(3)+RATE(669)*Y(125)*Y(3)&
    &+RATE(670)*Y(137)*Y(3)+RATE(671)*Y(151)*Y(3)+RATE(672)*Y(229)*Y(3)&
    &+RATE(728)*Y(188)*Y(3)+RATE(729)*Y(84)*Y(3)+RATE(731)*Y(230)*Y(3)&
    &+RATE(732)*Y(3)*Y(3)+RATE(732)*Y(3)*Y(3)+RATE(733)*Y(3)*Y(17)+RATE(734)&
    &*Y(3)*Y(22)+RATE(735)*Y(3)*Y(26)+RATE(736)*Y(3)*Y(34)+RATE(737)*Y(3)&
    &*Y(47)+RATE(738)*Y(3)*Y(58)+RATE(739)*Y(3)*Y(29)+RATE(740)*Y(3)*Y(37)&
    &+RATE(741)*Y(3)*Y(45)+RATE(742)*Y(3)*Y(97)+RATE(743)*Y(3)*Y(118)&
    &+RATE(744)*Y(3)*Y(118)+RATE(745)*Y(3)*Y(130)+RATE(746)*Y(3)*Y(130)&
    &+RATE(747)*Y(3)*Y(130)+RATE(748)*Y(3)*Y(150)+RATE(749)*Y(3)*Y(143)&
    &+RATE(750)*Y(3)*Y(222)+RATE(751)*Y(3)*Y(226)+RATE(755)*Y(246)*Y(3)&
    &+RATE(756)*Y(179)*Y(3)+RATE(757)*Y(327)*Y(3)+RATE(758)*Y(248)*Y(3)&
    &+RATE(770)*Y(164)*Y(3)+RATE(771)*Y(317)*Y(3)+RATE(772)*Y(110)*Y(3)&
    &+RATE(773)*Y(125)*Y(3)+RATE(774)*Y(137)*Y(3)+RATE(775)*Y(151)*Y(3)&
    &+RATE(776)*Y(229)*Y(3)+RATE(996)*Y(3)*totalSwap/safeMantle+RATE(1079)&
    &*Y(3)
    PROD = RATE(29)*Y(5)*bulkLayersReciprocal+RATE(84)*Y(34)+RATE(85)&
    &*Y(157)+RATE(86)*Y(157)+RATE(89)*Y(130)+RATE(91)*Y(61)+RATE(93)*Y(118)&
    &+RATE(95)*Y(246)+RATE(96)*Y(148)+RATE(98)*Y(52)+RATE(504)*D*Y(2)&
    &+RATE(505)*D*Y(4)+RATE(522)*D*Y(10)+RATE(522)*D*Y(10)+RATE(522)*D*Y(10)&
    &+RATE(543)*D*Y(15)
    YDOT(3) = PROD-LOSS
    LOSS = RATE(505)*D*Y(4)+RATE(1740)*D*Y(186)*Y(4)+RATE(1871)*D*Y(4)&
    &*Y(92)+RATE(1872)*D*Y(4)*Y(67)+RATE(1873)*D*Y(4)*Y(80)+RATE(1874)*D*Y(4)&
    &*Y(75)+RATE(1875)*D*Y(4)*Y(197)+RATE(1876)*D*Y(4)*Y(24)+RATE(1877)*D*Y(4&
    &)*Y(32)+RATE(1878)*D*Y(4)*Y(41)+RATE(1879)*D*Y(4)*Y(20)+RATE(1880)*D*Y(4&
    &)*Y(233)+RATE(1881)*D*Y(4)*Y(131)+RATE(1882)*D*Y(4)*Y(261)+RATE(1883)*D&
    &*Y(4)*Y(62)+RATE(1884)*D*Y(4)*Y(330)+RATE(1885)*D*Y(4)*Y(183)+RATE(1886)&
    &*D*Y(4)*Y(90)+RATE(1887)*D*Y(4)*Y(116)+RATE(1888)*D*Y(4)*Y(193)&
    &+RATE(1889)*D*Y(4)*Y(324)+RATE(1890)*D*Y(4)*Y(173)+RATE(1891)*D*Y(4)&
    &*Y(69)+RATE(1892)*D*Y(4)*Y(43)+RATE(1893)*D*Y(4)*Y(54)+RATE(1894)*D*Y(4)&
    &*Y(35)+RATE(1895)*D*Y(4)*Y(133)+RATE(1896)*D*Y(4)*Y(265)+RATE(1897)*D&
    &*Y(4)*Y(161)+RATE(1898)*D*Y(4)*Y(46)+RATE(1899)*D*Y(4)*Y(300)+RATE(1900)&
    &*D*Y(4)*Y(56)+RATE(1901)*D*Y(4)*Y(316)+RATE(1902)*D*Y(4)*Y(163)&
    &+RATE(1903)*D*Y(4)*Y(320)+RATE(1904)*D*Y(4)*Y(277)+RATE(1905)*D*Y(4)&
    &*Y(105)+RATE(1906)*D*Y(4)*Y(292)+RATE(1907)*D*Y(4)*Y(314)+RATE(1908)*D&
    &*Y(4)*Y(207)+RATE(1909)*D*Y(4)*Y(135)+RATE(1910)*D*Y(4)*Y(146)+RATE(1911&
    &)*D*Y(4)*Y(166)+RATE(1912)*D*Y(4)*Y(121)+RATE(1913)*D*Y(4)*Y(235)&
    &+RATE(1914)*D*Y(4)*Y(302)+RATE(1915)*D*Y(4)*Y(89)+RATE(1916)*D*Y(4)&
    &*Y(109)+RATE(1917)*D*Y(4)*Y(75)+RATE(1918)*D*Y(4)*Y(24)+RATE(1919)*D*Y(4&
    &)*Y(214)+RATE(1920)*D*Y(4)*Y(159)+RATE(1921)*D*Y(4)*Y(159)+RATE(1922)*D&
    &*Y(4)*Y(159)+RATE(1923)*D*Y(4)*Y(41)+RATE(1924)*D*Y(4)*Y(232)+RATE(1925)&
    &*D*Y(4)*Y(131)+RATE(1926)*D*Y(4)*Y(131)+RATE(1927)*D*Y(4)*Y(183)&
    &+RATE(1928)*D*Y(4)*Y(183)+RATE(1929)*D*Y(4)*Y(263)+RATE(1930)*D*Y(4)&
    &*Y(116)+RATE(1931)*D*Y(4)*Y(116)+RATE(1932)*D*Y(4)*Y(244)+RATE(1933)*D&
    &*Y(4)*Y(227)+RATE(1934)*D*Y(4)*Y(144)+RATE(1935)*D*Y(4)*Y(173)+RATE(1936&
    &)*D*Y(4)*Y(264)+RATE(1937)*D*Y(4)*Y(300)+RATE(1938)*D*Y(4)*Y(135)&
    &+RATE(1939)*D*Y(4)*Y(146)+RATE(1940)*D*Y(4)*Y(166)+RATE(1941)*D*Y(4)&
    &*Y(121)+RATE(1942)*D*Y(4)*Y(2)+RATE(1943)*D*Y(4)*Y(11)+RATE(1944)*D*Y(4)&
    &*Y(333)
    PROD = RATE(102)*Y(2)+RATE(103)*Y(6)+RATE(143)*Y(2)+RATE(848)*Y(21)&
    &+RATE(853)*Y(25)+RATE(881)*Y(8)+RATE(898)*Y(10)+RATE(913)*Y(174)&
    &+RATE(922)*Y(36)+RATE(1742)*D*Y(187)*Y(2)+RATE(1803)*D*Y(2)*Y(83)&
    &+RATE(1804)*D*Y(2)*Y(100)+RATE(1805)*D*Y(2)*Y(8)+RATE(1806)*D*Y(2)*Y(91)&
    &+RATE(1807)*D*Y(2)*Y(13)+RATE(1808)*D*Y(2)*Y(48)+RATE(1871)*D*Y(4)*Y(92)&
    &+RATE(1969)*D*Y(6)*Y(13)+RATE(2362)*D*Y(13)*Y(41)+RATE(2380)*D*Y(13)&
    &*Y(62)+RATE(2396)*D*Y(13)*Y(244)+RATE(2402)*D*Y(13)*Y(144)
    YDOT(4) = PROD-LOSS
    LOSS = RATE(29)*Y(5)*bulkLayersReciprocal+RATE(676)*Y(190)*Y(5)&
    &*bulkLayersReciprocal+RATE(677)*Y(86)*Y(5)*bulkLayersReciprocal+RATE(679&
    &)*Y(240)*Y(5)*bulkLayersReciprocal+RATE(680)*Y(5)*Y(5)&
    &*bulkLayersReciprocal+RATE(680)*Y(5)*Y(5)*bulkLayersReciprocal+RATE(681)&
    &*Y(5)*Y(19)*bulkLayersReciprocal+RATE(682)*Y(5)*Y(23)&
    &*bulkLayersReciprocal+RATE(683)*Y(5)*Y(30)*bulkLayersReciprocal+RATE(684&
    &)*Y(5)*Y(38)*bulkLayersReciprocal+RATE(685)*Y(5)*Y(51)&
    &*bulkLayersReciprocal+RATE(686)*Y(5)*Y(60)*bulkLayersReciprocal+RATE(687&
    &)*Y(5)*Y(31)*bulkLayersReciprocal+RATE(688)*Y(5)*Y(39)&
    &*bulkLayersReciprocal+RATE(689)*Y(5)*Y(50)*bulkLayersReciprocal+RATE(690&
    &)*Y(5)*Y(111)*bulkLayersReciprocal+RATE(691)*Y(5)*Y(126)&
    &*bulkLayersReciprocal+RATE(692)*Y(5)*Y(126)*bulkLayersReciprocal&
    &+RATE(693)*Y(5)*Y(139)*bulkLayersReciprocal+RATE(694)*Y(5)*Y(139)&
    &*bulkLayersReciprocal+RATE(695)*Y(5)*Y(139)*bulkLayersReciprocal&
    &+RATE(696)*Y(5)*Y(154)*bulkLayersReciprocal+RATE(697)*Y(5)*Y(152)&
    &*bulkLayersReciprocal+RATE(698)*Y(5)*Y(225)*bulkLayersReciprocal&
    &+RATE(699)*Y(5)*Y(228)*bulkLayersReciprocal+RATE(703)*Y(253)*Y(5)&
    &*bulkLayersReciprocal+RATE(704)*Y(181)*Y(5)*bulkLayersReciprocal&
    &+RATE(705)*Y(328)*Y(5)*bulkLayersReciprocal+RATE(706)*Y(254)*Y(5)&
    &*bulkLayersReciprocal+RATE(718)*Y(171)*Y(5)*bulkLayersReciprocal&
    &+RATE(719)*Y(323)*Y(5)*bulkLayersReciprocal+RATE(720)*Y(115)*Y(5)&
    &*bulkLayersReciprocal+RATE(721)*Y(128)*Y(5)*bulkLayersReciprocal&
    &+RATE(722)*Y(140)*Y(5)*bulkLayersReciprocal+RATE(723)*Y(155)*Y(5)&
    &*bulkLayersReciprocal+RATE(724)*Y(239)*Y(5)*bulkLayersReciprocal&
    &+RATE(780)*Y(190)*Y(5)*bulkLayersReciprocal+RATE(781)*Y(86)*Y(5)&
    &*bulkLayersReciprocal+RATE(783)*Y(240)*Y(5)*bulkLayersReciprocal&
    &+RATE(784)*Y(5)*Y(5)*bulkLayersReciprocal+RATE(784)*Y(5)*Y(5)&
    &*bulkLayersReciprocal+RATE(785)*Y(5)*Y(19)*bulkLayersReciprocal+RATE(786&
    &)*Y(5)*Y(23)*bulkLayersReciprocal+RATE(787)*Y(5)*Y(30)&
    &*bulkLayersReciprocal+RATE(788)*Y(5)*Y(38)*bulkLayersReciprocal+RATE(789&
    &)*Y(5)*Y(51)*bulkLayersReciprocal+RATE(790)*Y(5)*Y(60)&
    &*bulkLayersReciprocal+RATE(791)*Y(5)*Y(31)*bulkLayersReciprocal+RATE(792&
    &)*Y(5)*Y(39)*bulkLayersReciprocal+RATE(793)*Y(5)*Y(50)&
    &*bulkLayersReciprocal+RATE(794)*Y(5)*Y(111)*bulkLayersReciprocal&
    &+RATE(795)*Y(5)*Y(126)*bulkLayersReciprocal+RATE(796)*Y(5)*Y(126)&
    &*bulkLayersReciprocal+RATE(797)*Y(5)*Y(139)*bulkLayersReciprocal&
    &+RATE(798)*Y(5)*Y(139)*bulkLayersReciprocal+RATE(799)*Y(5)*Y(139)&
    &*bulkLayersReciprocal+RATE(800)*Y(5)*Y(154)*bulkLayersReciprocal&
    &+RATE(801)*Y(5)*Y(152)*bulkLayersReciprocal+RATE(802)*Y(5)*Y(225)&
    &*bulkLayersReciprocal+RATE(803)*Y(5)*Y(228)*bulkLayersReciprocal&
    &+RATE(807)*Y(253)*Y(5)*bulkLayersReciprocal+RATE(808)*Y(181)*Y(5)&
    &*bulkLayersReciprocal+RATE(809)*Y(328)*Y(5)*bulkLayersReciprocal&
    &+RATE(810)*Y(254)*Y(5)*bulkLayersReciprocal+RATE(822)*Y(171)*Y(5)&
    &*bulkLayersReciprocal+RATE(823)*Y(323)*Y(5)*bulkLayersReciprocal&
    &+RATE(824)*Y(115)*Y(5)*bulkLayersReciprocal+RATE(825)*Y(128)*Y(5)&
    &*bulkLayersReciprocal+RATE(826)*Y(140)*Y(5)*bulkLayersReciprocal&
    &+RATE(827)*Y(155)*Y(5)*bulkLayersReciprocal+RATE(828)*Y(239)*Y(5)&
    &*bulkLayersReciprocal+RATE(1162)*Y(5)
    PROD = RATE(996)*Y(3)*totalSwap/safeMantle
    YDOT(5) = PROD-LOSS
    LOSS = RATE(103)*Y(6)+RATE(104)*Y(6)+RATE(105)*Y(6)+RATE(506)*D*Y(6)&
    &+RATE(880)*Y(6)+RATE(1741)*D*Y(186)*Y(6)+RATE(1799)*D*Y(2)*Y(6)&
    &+RATE(1945)*D*Y(6)*Y(20)+RATE(1946)*D*Y(6)*Y(6)+RATE(1946)*D*Y(6)*Y(6)&
    &+RATE(1947)*D*Y(6)*Y(62)+RATE(1948)*D*Y(6)*Y(119)+RATE(1949)*D*Y(6)&
    &*Y(161)+RATE(1950)*D*Y(6)*Y(56)+RATE(1951)*D*Y(6)*Y(333)+RATE(1952)*D&
    &*Y(6)*Y(13)+RATE(1953)*D*Y(6)*Y(18)+RATE(1954)*D*Y(6)*Y(68)+RATE(1955)*D&
    &*Y(6)*Y(76)+RATE(1956)*D*Y(6)*Y(21)+RATE(1957)*D*Y(6)*Y(25)+RATE(1958)*D&
    &*Y(6)*Y(42)+RATE(1959)*D*Y(6)*Y(83)+RATE(1960)*D*Y(6)*Y(100)+RATE(1961)&
    &*D*Y(6)*Y(100)+RATE(1962)*D*Y(6)*Y(234)+RATE(1963)*D*Y(6)*Y(187)&
    &+RATE(1964)*D*Y(6)*Y(63)+RATE(1965)*D*Y(6)*Y(184)+RATE(1966)*D*Y(6)*Y(91&
    &)+RATE(1967)*D*Y(6)*Y(194)+RATE(1968)*D*Y(6)*Y(174)+RATE(1969)*D*Y(6)&
    &*Y(13)+RATE(1970)*D*Y(6)*Y(15)+RATE(1971)*D*Y(6)*Y(28)+RATE(1972)*D*Y(6)&
    &*Y(104)+RATE(1973)*D*Y(6)*Y(36)+RATE(1974)*D*Y(6)*Y(36)+RATE(1975)*D*Y(6&
    &)*Y(44)+RATE(1976)*D*Y(6)*Y(55)+RATE(1977)*D*Y(6)*Y(48)+RATE(1978)*D*Y(6&
    &)*Y(176)+RATE(1979)*D*Y(6)*Y(57)+RATE(1980)*D*Y(6)*Y(165)+RATE(1981)*D&
    &*Y(6)*Y(319)+RATE(1982)*D*Y(6)*Y(167)+RATE(1983)*D*Y(6)*Y(236)+RATE(1984&
    &)*D*Y(6)*Y(75)+RATE(1985)*D*Y(6)*Y(16)+RATE(1986)*D*Y(6)*Y(24)+RATE(1987&
    &)*D*Y(6)*Y(32)+RATE(1988)*D*Y(6)*Y(20)+RATE(1989)*D*Y(6)*Y(82)+RATE(1990&
    &)*D*Y(6)*Y(173)+RATE(1991)*D*Y(6)*Y(27)+RATE(1992)*D*Y(6)*Y(43)&
    &+RATE(1993)*D*Y(6)*Y(35)+RATE(1994)*D*Y(6)*Y(161)+RATE(1995)*D*Y(6)&
    &*Y(161)+RATE(1996)*D*Y(6)*Y(46)+RATE(1997)*D*Y(6)*Y(56)+RATE(1998)*D*Y(6&
    &)*Y(163)+RATE(1999)*D*Y(6)*Y(18)+RATE(2000)*D*Y(6)*Y(16)+RATE(2001)*D&
    &*Y(6)*Y(33)+RATE(2002)*D*Y(6)*Y(20)+RATE(2003)*D*Y(6)*Y(174)+RATE(2004)&
    &*D*Y(6)*Y(165)+RATE(2005)*D*Y(6)*Y(106)+RATE(2006)*D*Y(6)*Y(122)&
    &+RATE(2007)*D*Y(6)*Y(147)+RATE(2039)*D*Y(8)*Y(6)
    PROD = RATE(90)*Y(130)+RATE(117)*Y(109)+RATE(118)*Y(123)+RATE(130)&
    &*Y(32)+RATE(134)*Y(159)+RATE(136)*Y(41)+RATE(145)*Y(131)+RATE(146)*Y(261&
    &)+RATE(149)*Y(183)+RATE(151)*Y(263)+RATE(175)*Y(54)+RATE(197)*Y(166)&
    &+RATE(232)*Y(7)/safeMantle+RATE(315)*D*Y(7)/safeMantle*Y(2)+RATE(398)&
    &*Y(7)/safeMantle+RATE(453)*D*Y(2)*Y(3)/safeMantle+RATE(620)*D*Y(2)&
    &+RATE(648)*Y(130)*Y(47)+RATE(700)*Y(139)*Y(51)*bulkLayersReciprocal&
    &+RATE(732)*Y(3)*Y(3)+RATE(744)*Y(3)*Y(118)+RATE(747)*Y(3)*Y(130)&
    &+RATE(752)*Y(130)*Y(47)+RATE(784)*Y(5)*Y(5)*bulkLayersReciprocal&
    &+RATE(796)*Y(5)*Y(126)*bulkLayersReciprocal+RATE(799)*Y(5)*Y(139)&
    &*bulkLayersReciprocal+RATE(804)*Y(139)*Y(51)*bulkLayersReciprocal&
    &+RATE(840)*Y(109)+RATE(841)*Y(123)+RATE(851)*Y(25)+RATE(857)*Y(32)&
    &+RATE(858)*Y(33)+RATE(863)*Y(159)+RATE(866)*Y(41)+RATE(869)*Y(41)&
    &+RATE(870)*Y(42)+RATE(883)*Y(131)+RATE(887)*Y(261)+RATE(893)*Y(183)&
    &+RATE(895)*Y(263)+RATE(898)*Y(10)+RATE(927)*Y(54)+RATE(957)*Y(146)&
    &+RATE(958)*Y(166)+RATE(960)*Y(166)+RATE(1080)*Y(7)+RATE(1163)*Y(9)&
    &+RATE(1225)*D*Y(16)*Y(33)+RATE(1229)*D*Y(16)*Y(66)+RATE(1253)*D*Y(16)&
    &*Y(101)+RATE(1305)*D*Y(18)*Y(32)+RATE(1310)*D*Y(18)*Y(41)+RATE(1326)*D&
    &*Y(18)*Y(54)+RATE(1339)*D*Y(18)*Y(135)+RATE(1401)*D*Y(80)*Y(105)&
    &+RATE(1414)*D*Y(81)*Y(105)+RATE(1415)*D*Y(81)*Y(166)+RATE(1445)*D*Y(20)&
    &*Y(33)+RATE(1507)*D*Y(21)*Y(75)+RATE(1508)*D*Y(21)*Y(24)+RATE(1512)*D&
    &*Y(21)*Y(41)+RATE(1513)*D*Y(21)*Y(20)+RATE(1521)*D*Y(21)*Y(62)+RATE(1523&
    &)*D*Y(21)*Y(183)+RATE(1524)*D*Y(21)*Y(90)+RATE(1530)*D*Y(21)*Y(43)&
    &+RATE(1532)*D*Y(21)*Y(35)+RATE(1539)*D*Y(21)*Y(56)+RATE(1572)*D*Y(24)&
    &*Y(24)+RATE(1586)*D*Y(24)*Y(161)+RATE(1591)*D*Y(24)*Y(46)+RATE(1598)*D&
    &*Y(24)*Y(163)+RATE(1601)*D*Y(25)*Y(333)+RATE(1609)*D*Y(25)*Y(183)&
    &+RATE(1619)*D*Y(32)*Y(32)+RATE(1636)*D*Y(32)*Y(46)+RATE(1639)*D*Y(32)&
    &*Y(56)+RATE(1647)*D*Y(33)*Y(333)+RATE(1649)*D*Y(33)*Y(109)+RATE(1654)*D&
    &*Y(33)*Y(183)+RATE(1656)*D*Y(33)*Y(173)+RATE(1660)*D*Y(33)*Y(46)&
    &+RATE(1662)*D*Y(33)*Y(56)+RATE(1663)*D*Y(33)*Y(163)+RATE(1664)*D*Y(33)&
    &*Y(277)+RATE(1677)*D*Y(172)*Y(333)+RATE(1691)*D*Y(41)*Y(174)+RATE(1692)&
    &*D*Y(41)*Y(104)+RATE(1699)*D*Y(41)*Y(165)+RATE(1721)*D*Y(53)*Y(333)&
    &+RATE(1722)*D*Y(53)*Y(333)+RATE(1725)*D*Y(53)*Y(333)+RATE(1725)*D*Y(53)&
    &*Y(333)+RATE(1739)*D*Y(53)*Y(166)+RATE(1805)*D*Y(2)*Y(8)+RATE(1810)*D&
    &*Y(2)*Y(21)+RATE(1811)*D*Y(2)*Y(25)+RATE(1812)*D*Y(2)*Y(33)+RATE(1813)*D&
    &*Y(2)*Y(42)+RATE(1814)*D*Y(2)*Y(53)+RATE(1815)*D*Y(2)*Y(184)+RATE(1816)&
    &*D*Y(2)*Y(189)+RATE(1817)*D*Y(2)*Y(174)+RATE(1820)*D*Y(2)*Y(122)&
    &+RATE(1823)*D*Y(2)*Y(80)+RATE(1824)*D*Y(2)*Y(89)+RATE(1825)*D*Y(2)*Y(24)&
    &+RATE(1827)*D*Y(2)*Y(32)+RATE(1828)*D*Y(2)*Y(41)+RATE(1829)*D*Y(2)*Y(20)&
    &+RATE(1832)*D*Y(2)*Y(101)+RATE(1833)*D*Y(2)*Y(131)+RATE(1834)*D*Y(2)&
    &*Y(62)+RATE(1835)*D*Y(2)*Y(183)+RATE(1836)*D*Y(2)*Y(90)+RATE(1837)*D*Y(2&
    &)*Y(116)+RATE(1839)*D*Y(2)*Y(244)+RATE(1842)*D*Y(2)*Y(144)+RATE(1844)*D&
    &*Y(2)*Y(173)+RATE(1846)*D*Y(2)*Y(43)+RATE(1847)*D*Y(2)*Y(54)+RATE(1848)&
    &*D*Y(2)*Y(35)+RATE(1856)*D*Y(2)*Y(175)+RATE(1862)*D*Y(2)*Y(56)+RATE(1915&
    &)*D*Y(4)*Y(89)+RATE(1916)*D*Y(4)*Y(109)+RATE(1917)*D*Y(4)*Y(75)&
    &+RATE(1918)*D*Y(4)*Y(24)+RATE(1921)*D*Y(4)*Y(159)+RATE(1922)*D*Y(4)&
    &*Y(159)+RATE(1922)*D*Y(4)*Y(159)+RATE(1923)*D*Y(4)*Y(41)+RATE(1925)*D&
    &*Y(4)*Y(131)+RATE(1926)*D*Y(4)*Y(131)+RATE(1927)*D*Y(4)*Y(183)+RATE(1928&
    &)*D*Y(4)*Y(183)+RATE(1929)*D*Y(4)*Y(263)+RATE(1930)*D*Y(4)*Y(116)&
    &+RATE(1932)*D*Y(4)*Y(244)+RATE(1934)*D*Y(4)*Y(144)+RATE(1935)*D*Y(4)&
    &*Y(173)+RATE(1938)*D*Y(4)*Y(135)+RATE(1939)*D*Y(4)*Y(146)+RATE(1940)*D&
    &*Y(4)*Y(166)+RATE(1941)*D*Y(4)*Y(121)+RATE(1945)*D*Y(6)*Y(20)+RATE(1946)&
    &*D*Y(6)*Y(6)+RATE(1947)*D*Y(6)*Y(62)+RATE(1948)*D*Y(6)*Y(119)+RATE(1949)&
    &*D*Y(6)*Y(161)+RATE(1950)*D*Y(6)*Y(56)+RATE(2008)*D*Y(8)*Y(67)+RATE(2009&
    &)*D*Y(8)*Y(80)+RATE(2010)*D*Y(8)*Y(75)+RATE(2011)*D*Y(8)*Y(24)+RATE(2012&
    &)*D*Y(8)*Y(41)+RATE(2013)*D*Y(8)*Y(20)+RATE(2014)*D*Y(8)*Y(82)+RATE(2015&
    &)*D*Y(8)*Y(99)+RATE(2016)*D*Y(8)*Y(131)+RATE(2017)*D*Y(8)*Y(62)&
    &+RATE(2018)*D*Y(8)*Y(183)+RATE(2019)*D*Y(8)*Y(90)+RATE(2020)*D*Y(8)&
    &*Y(116)+RATE(2021)*D*Y(8)*Y(43)+RATE(2022)*D*Y(8)*Y(54)+RATE(2023)*D*Y(8&
    &)*Y(35)+RATE(2024)*D*Y(8)*Y(133)+RATE(2025)*D*Y(8)*Y(161)+RATE(2026)*D&
    &*Y(8)*Y(56)+RATE(2029)*D*Y(8)*Y(109)+RATE(2029)*D*Y(8)*Y(109)+RATE(2033)&
    &*D*Y(8)*Y(41)+RATE(2040)*D*Y(8)*Y(131)+RATE(2042)*D*Y(8)*Y(183)&
    &+RATE(2043)*D*Y(8)*Y(183)+RATE(2043)*D*Y(8)*Y(183)+RATE(2064)*D*Y(132)&
    &*Y(333)+RATE(2076)*D*Y(160)*Y(333)+RATE(2119)*D*Y(63)*Y(333)+RATE(2147)&
    &*D*Y(10)*Y(333)+RATE(2149)*D*Y(10)*Y(67)+RATE(2150)*D*Y(10)*Y(75)&
    &+RATE(2151)*D*Y(10)*Y(197)+RATE(2152)*D*Y(10)*Y(16)+RATE(2153)*D*Y(10)&
    &*Y(24)+RATE(2154)*D*Y(10)*Y(32)+RATE(2158)*D*Y(10)*Y(214)+RATE(2159)*D&
    &*Y(10)*Y(159)+RATE(2160)*D*Y(10)*Y(159)+RATE(2161)*D*Y(10)*Y(41)&
    &+RATE(2162)*D*Y(10)*Y(20)+RATE(2163)*D*Y(10)*Y(82)+RATE(2164)*D*Y(10)&
    &*Y(232)+RATE(2165)*D*Y(10)*Y(99)+RATE(2166)*D*Y(10)*Y(99)+RATE(2167)*D&
    &*Y(10)*Y(233)+RATE(2168)*D*Y(10)*Y(186)+RATE(2169)*D*Y(10)*Y(131)&
    &+RATE(2170)*D*Y(10)*Y(261)+RATE(2171)*D*Y(10)*Y(62)+RATE(2172)*D*Y(10)&
    &*Y(183)+RATE(2173)*D*Y(10)*Y(90)+RATE(2174)*D*Y(10)*Y(116)+RATE(2175)*D&
    &*Y(10)*Y(267)+RATE(2176)*D*Y(10)*Y(267)+RATE(2177)*D*Y(10)*Y(244)&
    &+RATE(2178)*D*Y(10)*Y(193)+RATE(2179)*D*Y(10)*Y(92)+RATE(2180)*D*Y(10)&
    &*Y(144)+RATE(2181)*D*Y(10)*Y(324)+RATE(2182)*D*Y(10)*Y(173)+RATE(2183)*D&
    &*Y(10)*Y(69)+RATE(2184)*D*Y(10)*Y(103)+RATE(2185)*D*Y(10)*Y(43)&
    &+RATE(2186)*D*Y(10)*Y(54)+RATE(2187)*D*Y(10)*Y(35)+RATE(2188)*D*Y(10)&
    &*Y(264)+RATE(2189)*D*Y(10)*Y(133)+RATE(2190)*D*Y(10)*Y(265)+RATE(2191)*D&
    &*Y(10)*Y(161)+RATE(2193)*D*Y(10)*Y(46)+RATE(2194)*D*Y(10)*Y(300)&
    &+RATE(2195)*D*Y(10)*Y(56)+RATE(2196)*D*Y(10)*Y(316)+RATE(2197)*D*Y(10)&
    &*Y(163)+RATE(2198)*D*Y(10)*Y(320)+RATE(2199)*D*Y(10)*Y(277)+RATE(2200)*D&
    &*Y(10)*Y(105)+RATE(2201)*D*Y(10)*Y(135)+RATE(2202)*D*Y(10)*Y(146)&
    &+RATE(2203)*D*Y(10)*Y(166)+RATE(2204)*D*Y(10)*Y(121)+RATE(2205)*D*Y(10)&
    &*Y(235)+RATE(2206)*D*Y(10)*Y(302)+RATE(2209)*D*Y(142)*Y(333)+RATE(2214)&
    &*D*Y(274)*Y(333)+RATE(2217)*D*Y(66)*Y(333)+RATE(2218)*D*Y(66)*Y(333)&
    &+RATE(2236)*D*Y(189)*Y(333)+RATE(2238)*D*Y(189)*Y(333)+RATE(2239)*D&
    &*Y(193)*Y(2)+RATE(2287)*D*Y(116)*Y(116)+RATE(2334)*D*Y(13)*Y(80)&
    &+RATE(2337)*D*Y(13)*Y(89)+RATE(2339)*D*Y(13)*Y(109)+RATE(2340)*D*Y(13)&
    &*Y(109)+RATE(2348)*D*Y(13)*Y(24)+RATE(2352)*D*Y(13)*Y(32)+RATE(2353)*D&
    &*Y(13)*Y(206)+RATE(2353)*D*Y(13)*Y(206)+RATE(2359)*D*Y(13)*Y(41)&
    &+RATE(2360)*D*Y(13)*Y(41)+RATE(2373)*D*Y(13)*Y(131)+RATE(2376)*D*Y(13)&
    &*Y(261)+RATE(2384)*D*Y(13)*Y(183)+RATE(2408)*D*Y(13)*Y(43)+RATE(2410)*D&
    &*Y(13)*Y(54)+RATE(2434)*D*Y(13)*Y(135)+RATE(2436)*D*Y(13)*Y(146)&
    &+RATE(2438)*D*Y(13)*Y(166)+RATE(2438)*D*Y(13)*Y(166)+RATE(2439)*D*Y(13)&
    &*Y(166)+RATE(2466)*D*Y(174)*Y(183)+RATE(2493)*D*Y(27)*Y(81)+RATE(2499)*D&
    &*Y(27)*Y(63)+RATE(2500)*D*Y(27)*Y(184)+RATE(2525)*D*Y(27)*Y(32)&
    &+RATE(2577)*D*Y(28)*Y(41)+RATE(2589)*D*Y(28)*Y(54)+RATE(2611)*D*Y(104)&
    &*Y(183)+RATE(2626)*D*Y(35)*Y(33)+RATE(2646)*D*Y(35)*Y(35)+RATE(2680)*D&
    &*Y(36)*Y(62)+RATE(2805)*D*Y(64)*Y(333)+RATE(2831)*D*Y(46)*Y(53)&
    &+RATE(2834)*D*Y(46)*Y(63)+RATE(2836)*D*Y(46)*Y(184)+RATE(2845)*D*Y(46)&
    &*Y(55)+RATE(2852)*D*Y(46)*Y(147)+RATE(2858)*D*Y(46)*Y(109)+RATE(2873)*D&
    &*Y(46)*Y(101)+RATE(2910)*D*Y(46)*Y(135)+RATE(3065)*D*Y(165)*Y(183)&
    &+RATE(3092)*D*Y(136)*Y(333)+RATE(3098)*D*Y(147)*Y(333)+RATE(3099)*D&
    &*Y(167)*Y(333)+RATE(3101)*D*Y(177)*Y(333)
    YDOT(6) = PROD-LOSS
    LOSS = RATE(232)*Y(7)/safeMantle+RATE(315)*D*Y(7)/safeMantle*Y(2)&
    &+RATE(398)*Y(7)/safeMantle+RATE(997)*Y(7)*totalSwap/safeMantle+RATE(1080&
    &)*Y(7)
    PROD = RATE(30)*Y(9)*bulkLayersReciprocal+RATE(88)*Y(40)+RATE(97)&
    &*Y(52)+RATE(452)*D*Y(2)*Y(3)/safeMantle+RATE(506)*D*Y(6)+RATE(507)*D*Y(8&
    &)+RATE(628)*Y(3)*Y(3)+RATE(640)*Y(3)*Y(118)+RATE(643)*Y(3)*Y(130)
    YDOT(7) = PROD-LOSS
    LOSS = RATE(507)*D*Y(8)+RATE(881)*Y(8)+RATE(1805)*D*Y(2)*Y(8)&
    &+RATE(2008)*D*Y(8)*Y(67)+RATE(2009)*D*Y(8)*Y(80)+RATE(2010)*D*Y(8)*Y(75)&
    &+RATE(2011)*D*Y(8)*Y(24)+RATE(2012)*D*Y(8)*Y(41)+RATE(2013)*D*Y(8)*Y(20)&
    &+RATE(2014)*D*Y(8)*Y(82)+RATE(2015)*D*Y(8)*Y(99)+RATE(2016)*D*Y(8)*Y(131&
    &)+RATE(2017)*D*Y(8)*Y(62)+RATE(2018)*D*Y(8)*Y(183)+RATE(2019)*D*Y(8)&
    &*Y(90)+RATE(2020)*D*Y(8)*Y(116)+RATE(2021)*D*Y(8)*Y(43)+RATE(2022)*D*Y(8&
    &)*Y(54)+RATE(2023)*D*Y(8)*Y(35)+RATE(2024)*D*Y(8)*Y(133)+RATE(2025)*D&
    &*Y(8)*Y(161)+RATE(2026)*D*Y(8)*Y(56)+RATE(2027)*D*Y(8)*Y(333)+RATE(2028)&
    &*D*Y(8)*Y(67)+RATE(2029)*D*Y(8)*Y(109)+RATE(2030)*D*Y(8)*Y(75)+RATE(2031&
    &)*D*Y(8)*Y(16)+RATE(2032)*D*Y(8)*Y(24)+RATE(2033)*D*Y(8)*Y(41)+RATE(2034&
    &)*D*Y(8)*Y(41)+RATE(2035)*D*Y(8)*Y(20)+RATE(2036)*D*Y(8)*Y(82)+RATE(2037&
    &)*D*Y(8)*Y(232)+RATE(2038)*D*Y(8)*Y(99)+RATE(2039)*D*Y(8)*Y(6)+RATE(2040&
    &)*D*Y(8)*Y(131)+RATE(2041)*D*Y(8)*Y(62)+RATE(2042)*D*Y(8)*Y(183)&
    &+RATE(2043)*D*Y(8)*Y(183)+RATE(2044)*D*Y(8)*Y(116)+RATE(2045)*D*Y(8)&
    &*Y(11)+RATE(2046)*D*Y(8)*Y(103)+RATE(2047)*D*Y(8)*Y(27)+RATE(2048)*D*Y(8&
    &)*Y(35)+RATE(2049)*D*Y(8)*Y(133)+RATE(2050)*D*Y(8)*Y(161)+RATE(2051)*D&
    &*Y(8)*Y(46)+RATE(2052)*D*Y(8)*Y(56)
    PROD = RATE(104)*Y(6)+RATE(897)*Y(10)+RATE(1818)*D*Y(2)*Y(15)&
    &+RATE(1931)*D*Y(4)*Y(116)+RATE(1942)*D*Y(4)*Y(2)+RATE(1952)*D*Y(6)*Y(13&
    &)
    YDOT(8) = PROD-LOSS
    LOSS = RATE(30)*Y(9)*bulkLayersReciprocal+RATE(1163)*Y(9)
    PROD = RATE(680)*Y(5)*Y(5)*bulkLayersReciprocal+RATE(692)*Y(5)*Y(126&
    &)*bulkLayersReciprocal+RATE(695)*Y(5)*Y(139)*bulkLayersReciprocal&
    &+RATE(997)*Y(7)*totalSwap/safeMantle
    YDOT(9) = PROD-LOSS
    LOSS = RATE(522)*D*Y(10)+RATE(897)*Y(10)+RATE(898)*Y(10)+RATE(2147)&
    &*D*Y(10)*Y(333)+RATE(2148)*D*Y(10)*Y(333)+RATE(2149)*D*Y(10)*Y(67)&
    &+RATE(2150)*D*Y(10)*Y(75)+RATE(2151)*D*Y(10)*Y(197)+RATE(2152)*D*Y(10)&
    &*Y(16)+RATE(2153)*D*Y(10)*Y(24)+RATE(2154)*D*Y(10)*Y(32)+RATE(2155)*D&
    &*Y(10)*Y(237)+RATE(2156)*D*Y(10)*Y(237)+RATE(2157)*D*Y(10)*Y(237)&
    &+RATE(2158)*D*Y(10)*Y(214)+RATE(2159)*D*Y(10)*Y(159)+RATE(2160)*D*Y(10)&
    &*Y(159)+RATE(2161)*D*Y(10)*Y(41)+RATE(2162)*D*Y(10)*Y(20)+RATE(2163)*D&
    &*Y(10)*Y(82)+RATE(2164)*D*Y(10)*Y(232)+RATE(2165)*D*Y(10)*Y(99)&
    &+RATE(2166)*D*Y(10)*Y(99)+RATE(2167)*D*Y(10)*Y(233)+RATE(2168)*D*Y(10)&
    &*Y(186)+RATE(2169)*D*Y(10)*Y(131)+RATE(2170)*D*Y(10)*Y(261)+RATE(2171)*D&
    &*Y(10)*Y(62)+RATE(2172)*D*Y(10)*Y(183)+RATE(2173)*D*Y(10)*Y(90)&
    &+RATE(2174)*D*Y(10)*Y(116)+RATE(2175)*D*Y(10)*Y(267)+RATE(2176)*D*Y(10)&
    &*Y(267)+RATE(2177)*D*Y(10)*Y(244)+RATE(2178)*D*Y(10)*Y(193)+RATE(2179)*D&
    &*Y(10)*Y(92)+RATE(2180)*D*Y(10)*Y(144)+RATE(2181)*D*Y(10)*Y(324)&
    &+RATE(2182)*D*Y(10)*Y(173)+RATE(2183)*D*Y(10)*Y(69)+RATE(2184)*D*Y(10)&
    &*Y(103)+RATE(2185)*D*Y(10)*Y(43)+RATE(2186)*D*Y(10)*Y(54)+RATE(2187)*D&
    &*Y(10)*Y(35)+RATE(2188)*D*Y(10)*Y(264)+RATE(2189)*D*Y(10)*Y(133)&
    &+RATE(2190)*D*Y(10)*Y(265)+RATE(2191)*D*Y(10)*Y(161)+RATE(2192)*D*Y(10)&
    &*Y(46)+RATE(2193)*D*Y(10)*Y(46)+RATE(2194)*D*Y(10)*Y(300)+RATE(2195)*D&
    &*Y(10)*Y(56)+RATE(2196)*D*Y(10)*Y(316)+RATE(2197)*D*Y(10)*Y(163)&
    &+RATE(2198)*D*Y(10)*Y(320)+RATE(2199)*D*Y(10)*Y(277)+RATE(2200)*D*Y(10)&
    &*Y(105)+RATE(2201)*D*Y(10)*Y(135)+RATE(2202)*D*Y(10)*Y(146)+RATE(2203)*D&
    &*Y(10)*Y(166)+RATE(2204)*D*Y(10)*Y(121)+RATE(2205)*D*Y(10)*Y(235)&
    &+RATE(2206)*D*Y(10)*Y(302)
    PROD = RATE(1970)*D*Y(6)*Y(15)+RATE(1973)*D*Y(6)*Y(36)+RATE(1978)*D&
    &*Y(6)*Y(176)+RATE(2039)*D*Y(8)*Y(6)+RATE(2044)*D*Y(8)*Y(116)
    YDOT(10) = PROD-LOSS
    LOSS = RATE(106)*Y(11)+RATE(159)*Y(11)+RATE(541)*D*Y(11)+RATE(1943)&
    &*D*Y(4)*Y(11)+RATE(2045)*D*Y(8)*Y(11)
    PROD = RATE(248)*Y(12)/safeMantle+RATE(331)*D*Y(12)/safeMantle*Y(2)&
    &+RATE(414)*Y(12)/safeMantle+RATE(1096)*Y(12)+RATE(1179)*Y(14)+RATE(1807)&
    &*D*Y(2)*Y(13)+RATE(1818)*D*Y(2)*Y(15)+RATE(1952)*D*Y(6)*Y(13)+RATE(1969)&
    &*D*Y(6)*Y(13)+RATE(1970)*D*Y(6)*Y(15)+RATE(2320)*D*Y(13)*Y(67)+RATE(2321&
    &)*D*Y(13)*Y(80)+RATE(2322)*D*Y(13)*Y(16)+RATE(2323)*D*Y(13)*Y(41)&
    &+RATE(2324)*D*Y(13)*Y(20)+RATE(2325)*D*Y(13)*Y(131)+RATE(2326)*D*Y(13)&
    &*Y(62)+RATE(2327)*D*Y(13)*Y(183)+RATE(2328)*D*Y(13)*Y(103)+RATE(2329)*D&
    &*Y(13)*Y(54)+RATE(2330)*D*Y(13)*Y(161)+RATE(2331)*D*Y(13)*Y(320)&
    &+RATE(2332)*D*Y(13)*Y(105)+RATE(2333)*D*Y(13)*Y(67)+RATE(2334)*D*Y(13)&
    &*Y(80)+RATE(2335)*D*Y(13)*Y(80)+RATE(2336)*D*Y(13)*Y(80)+RATE(2337)*D&
    &*Y(13)*Y(89)+RATE(2338)*D*Y(13)*Y(89)+RATE(2339)*D*Y(13)*Y(109)&
    &+RATE(2340)*D*Y(13)*Y(109)+RATE(2341)*D*Y(13)*Y(109)+RATE(2342)*D*Y(13)&
    &*Y(75)+RATE(2343)*D*Y(13)*Y(75)+RATE(2344)*D*Y(13)*Y(75)+RATE(2345)*D&
    &*Y(13)*Y(197)+RATE(2346)*D*Y(13)*Y(284)+RATE(2347)*D*Y(13)*Y(280)&
    &+RATE(2348)*D*Y(13)*Y(24)+RATE(2349)*D*Y(13)*Y(24)+RATE(2350)*D*Y(13)&
    &*Y(218)+RATE(2351)*D*Y(13)*Y(218)+RATE(2352)*D*Y(13)*Y(32)+RATE(2353)*D&
    &*Y(13)*Y(206)+RATE(2354)*D*Y(13)*Y(237)+RATE(2355)*D*Y(13)*Y(214)&
    &+RATE(2356)*D*Y(13)*Y(214)+RATE(2357)*D*Y(13)*Y(159)+RATE(2358)*D*Y(13)&
    &*Y(159)+RATE(2359)*D*Y(13)*Y(41)+RATE(2360)*D*Y(13)*Y(41)+RATE(2361)*D&
    &*Y(13)*Y(41)+RATE(2362)*D*Y(13)*Y(41)+RATE(2363)*D*Y(13)*Y(20)+RATE(2364&
    &)*D*Y(13)*Y(82)+RATE(2365)*D*Y(13)*Y(82)+RATE(2366)*D*Y(13)*Y(232)&
    &+RATE(2367)*D*Y(13)*Y(232)+RATE(2368)*D*Y(13)*Y(232)+RATE(2369)*D*Y(13)&
    &*Y(232)+RATE(2370)*D*Y(13)*Y(99)+RATE(2371)*D*Y(13)*Y(233)+RATE(2372)*D&
    &*Y(13)*Y(233)+RATE(2373)*D*Y(13)*Y(131)+RATE(2374)*D*Y(13)*Y(131)&
    &+RATE(2375)*D*Y(13)*Y(131)+RATE(2376)*D*Y(13)*Y(261)+RATE(2377)*D*Y(13)&
    &*Y(261)+RATE(2378)*D*Y(13)*Y(261)+RATE(2379)*D*Y(13)*Y(62)+RATE(2380)*D&
    &*Y(13)*Y(62)+RATE(2381)*D*Y(13)*Y(330)+RATE(2382)*D*Y(13)*Y(330)&
    &+RATE(2383)*D*Y(13)*Y(183)+RATE(2384)*D*Y(13)*Y(183)+RATE(2385)*D*Y(13)&
    &*Y(263)+RATE(2386)*D*Y(13)*Y(288)+RATE(2387)*D*Y(13)*Y(288)+RATE(2388)*D&
    &*Y(13)*Y(90)+RATE(2389)*D*Y(13)*Y(90)+RATE(2390)*D*Y(13)*Y(90)+RATE(2391&
    &)*D*Y(13)*Y(90)+RATE(2392)*D*Y(13)*Y(116)+RATE(2394)*D*Y(13)*Y(116)&
    &+RATE(2395)*D*Y(13)*Y(244)+RATE(2396)*D*Y(13)*Y(244)+RATE(2397)*D*Y(13)&
    &*Y(193)+RATE(2398)*D*Y(13)*Y(92)+RATE(2399)*D*Y(13)*Y(92)+RATE(2400)*D&
    &*Y(13)*Y(92)+RATE(2401)*D*Y(13)*Y(144)+RATE(2402)*D*Y(13)*Y(144)&
    &+RATE(2403)*D*Y(13)*Y(324)+RATE(2404)*D*Y(13)*Y(324)+RATE(2405)*D*Y(13)&
    &*Y(173)+RATE(2406)*D*Y(13)*Y(103)+RATE(2407)*D*Y(13)*Y(291)+RATE(2408)*D&
    &*Y(13)*Y(43)+RATE(2409)*D*Y(13)*Y(43)+RATE(2410)*D*Y(13)*Y(54)+RATE(2411&
    &)*D*Y(13)*Y(54)+RATE(2412)*D*Y(13)*Y(35)+RATE(2413)*D*Y(13)*Y(133)&
    &+RATE(2414)*D*Y(13)*Y(133)+RATE(2415)*D*Y(13)*Y(265)+RATE(2416)*D*Y(13)&
    &*Y(265)+RATE(2417)*D*Y(13)*Y(161)+RATE(2418)*D*Y(13)*Y(221)+RATE(2419)*D&
    &*Y(13)*Y(221)+RATE(2420)*D*Y(13)*Y(300)+RATE(2421)*D*Y(13)*Y(300)&
    &+RATE(2422)*D*Y(13)*Y(300)+RATE(2423)*D*Y(13)*Y(300)+RATE(2424)*D*Y(13)&
    &*Y(56)+RATE(2425)*D*Y(13)*Y(316)+RATE(2426)*D*Y(13)*Y(320)+RATE(2427)*D&
    &*Y(13)*Y(320)+RATE(2428)*D*Y(13)*Y(277)+RATE(2429)*D*Y(13)*Y(277)&
    &+RATE(2430)*D*Y(13)*Y(292)+RATE(2431)*D*Y(13)*Y(314)+RATE(2432)*D*Y(13)&
    &*Y(207)+RATE(2433)*D*Y(13)*Y(207)+RATE(2434)*D*Y(13)*Y(135)+RATE(2435)*D&
    &*Y(13)*Y(135)+RATE(2436)*D*Y(13)*Y(146)+RATE(2437)*D*Y(13)*Y(146)&
    &+RATE(2438)*D*Y(13)*Y(166)+RATE(2439)*D*Y(13)*Y(166)+RATE(2440)*D*Y(13)&
    &*Y(121)+RATE(2441)*D*Y(13)*Y(235)+RATE(2442)*D*Y(13)*Y(235)+RATE(2443)*D&
    &*Y(13)*Y(302)+RATE(2444)*D*Y(13)*Y(302)+RATE(2445)*D*Y(13)*Y(333)&
    &+RATE(2446)*D*Y(15)*Y(333)
    YDOT(11) = PROD-LOSS
    LOSS = RATE(248)*Y(12)/safeMantle+RATE(331)*D*Y(12)/safeMantle*Y(2)&
    &+RATE(414)*Y(12)/safeMantle+RATE(1013)*Y(12)*totalSwap/safeMantle&
    &+RATE(1096)*Y(12)
    PROD = RATE(46)*Y(14)*bulkLayersReciprocal+RATE(541)*D*Y(11)&
    &+RATE(542)*D*Y(13)+RATE(543)*D*Y(15)
    YDOT(12) = PROD-LOSS
    LOSS = RATE(542)*D*Y(13)+RATE(1807)*D*Y(2)*Y(13)+RATE(1952)*D*Y(6)&
    &*Y(13)+RATE(1969)*D*Y(6)*Y(13)+RATE(2320)*D*Y(13)*Y(67)+RATE(2321)*D&
    &*Y(13)*Y(80)+RATE(2322)*D*Y(13)*Y(16)+RATE(2323)*D*Y(13)*Y(41)+RATE(2324&
    &)*D*Y(13)*Y(20)+RATE(2325)*D*Y(13)*Y(131)+RATE(2326)*D*Y(13)*Y(62)&
    &+RATE(2327)*D*Y(13)*Y(183)+RATE(2328)*D*Y(13)*Y(103)+RATE(2329)*D*Y(13)&
    &*Y(54)+RATE(2330)*D*Y(13)*Y(161)+RATE(2331)*D*Y(13)*Y(320)+RATE(2332)*D&
    &*Y(13)*Y(105)+RATE(2333)*D*Y(13)*Y(67)+RATE(2334)*D*Y(13)*Y(80)&
    &+RATE(2335)*D*Y(13)*Y(80)+RATE(2336)*D*Y(13)*Y(80)+RATE(2337)*D*Y(13)&
    &*Y(89)+RATE(2338)*D*Y(13)*Y(89)+RATE(2339)*D*Y(13)*Y(109)+RATE(2340)*D&
    &*Y(13)*Y(109)+RATE(2341)*D*Y(13)*Y(109)+RATE(2342)*D*Y(13)*Y(75)&
    &+RATE(2343)*D*Y(13)*Y(75)+RATE(2344)*D*Y(13)*Y(75)+RATE(2345)*D*Y(13)&
    &*Y(197)+RATE(2346)*D*Y(13)*Y(284)+RATE(2347)*D*Y(13)*Y(280)+RATE(2348)*D&
    &*Y(13)*Y(24)+RATE(2349)*D*Y(13)*Y(24)+RATE(2350)*D*Y(13)*Y(218)&
    &+RATE(2351)*D*Y(13)*Y(218)+RATE(2352)*D*Y(13)*Y(32)+RATE(2353)*D*Y(13)&
    &*Y(206)+RATE(2354)*D*Y(13)*Y(237)+RATE(2355)*D*Y(13)*Y(214)+RATE(2356)*D&
    &*Y(13)*Y(214)+RATE(2357)*D*Y(13)*Y(159)+RATE(2358)*D*Y(13)*Y(159)&
    &+RATE(2359)*D*Y(13)*Y(41)+RATE(2360)*D*Y(13)*Y(41)+RATE(2361)*D*Y(13)&
    &*Y(41)+RATE(2362)*D*Y(13)*Y(41)+RATE(2363)*D*Y(13)*Y(20)+RATE(2364)*D&
    &*Y(13)*Y(82)+RATE(2365)*D*Y(13)*Y(82)+RATE(2366)*D*Y(13)*Y(232)&
    &+RATE(2367)*D*Y(13)*Y(232)+RATE(2368)*D*Y(13)*Y(232)+RATE(2369)*D*Y(13)&
    &*Y(232)+RATE(2370)*D*Y(13)*Y(99)+RATE(2371)*D*Y(13)*Y(233)+RATE(2372)*D&
    &*Y(13)*Y(233)+RATE(2373)*D*Y(13)*Y(131)+RATE(2374)*D*Y(13)*Y(131)&
    &+RATE(2375)*D*Y(13)*Y(131)+RATE(2376)*D*Y(13)*Y(261)+RATE(2377)*D*Y(13)&
    &*Y(261)+RATE(2378)*D*Y(13)*Y(261)+RATE(2379)*D*Y(13)*Y(62)+RATE(2380)*D&
    &*Y(13)*Y(62)+RATE(2381)*D*Y(13)*Y(330)+RATE(2382)*D*Y(13)*Y(330)&
    &+RATE(2383)*D*Y(13)*Y(183)+RATE(2384)*D*Y(13)*Y(183)+RATE(2385)*D*Y(13)&
    &*Y(263)+RATE(2386)*D*Y(13)*Y(288)+RATE(2387)*D*Y(13)*Y(288)+RATE(2388)*D&
    &*Y(13)*Y(90)+RATE(2389)*D*Y(13)*Y(90)+RATE(2390)*D*Y(13)*Y(90)+RATE(2391&
    &)*D*Y(13)*Y(90)+RATE(2392)*D*Y(13)*Y(116)+RATE(2393)*D*Y(13)*Y(116)&
    &+RATE(2394)*D*Y(13)*Y(116)+RATE(2395)*D*Y(13)*Y(244)+RATE(2396)*D*Y(13)&
    &*Y(244)+RATE(2397)*D*Y(13)*Y(193)+RATE(2398)*D*Y(13)*Y(92)+RATE(2399)*D&
    &*Y(13)*Y(92)+RATE(2400)*D*Y(13)*Y(92)+RATE(2401)*D*Y(13)*Y(144)&
    &+RATE(2402)*D*Y(13)*Y(144)+RATE(2403)*D*Y(13)*Y(324)+RATE(2404)*D*Y(13)&
    &*Y(324)+RATE(2405)*D*Y(13)*Y(173)+RATE(2406)*D*Y(13)*Y(103)+RATE(2407)*D&
    &*Y(13)*Y(291)+RATE(2408)*D*Y(13)*Y(43)+RATE(2409)*D*Y(13)*Y(43)&
    &+RATE(2410)*D*Y(13)*Y(54)+RATE(2411)*D*Y(13)*Y(54)+RATE(2412)*D*Y(13)&
    &*Y(35)+RATE(2413)*D*Y(13)*Y(133)+RATE(2414)*D*Y(13)*Y(133)+RATE(2415)*D&
    &*Y(13)*Y(265)+RATE(2416)*D*Y(13)*Y(265)+RATE(2417)*D*Y(13)*Y(161)&
    &+RATE(2418)*D*Y(13)*Y(221)+RATE(2419)*D*Y(13)*Y(221)+RATE(2420)*D*Y(13)&
    &*Y(300)+RATE(2421)*D*Y(13)*Y(300)+RATE(2422)*D*Y(13)*Y(300)+RATE(2423)*D&
    &*Y(13)*Y(300)+RATE(2424)*D*Y(13)*Y(56)+RATE(2425)*D*Y(13)*Y(316)&
    &+RATE(2426)*D*Y(13)*Y(320)+RATE(2427)*D*Y(13)*Y(320)+RATE(2428)*D*Y(13)&
    &*Y(277)+RATE(2429)*D*Y(13)*Y(277)+RATE(2430)*D*Y(13)*Y(292)+RATE(2431)*D&
    &*Y(13)*Y(314)+RATE(2432)*D*Y(13)*Y(207)+RATE(2433)*D*Y(13)*Y(207)&
    &+RATE(2434)*D*Y(13)*Y(135)+RATE(2435)*D*Y(13)*Y(135)+RATE(2436)*D*Y(13)&
    &*Y(146)+RATE(2437)*D*Y(13)*Y(146)+RATE(2438)*D*Y(13)*Y(166)+RATE(2439)*D&
    &*Y(13)*Y(166)+RATE(2440)*D*Y(13)*Y(121)+RATE(2441)*D*Y(13)*Y(235)&
    &+RATE(2442)*D*Y(13)*Y(235)+RATE(2443)*D*Y(13)*Y(302)+RATE(2444)*D*Y(13)&
    &*Y(302)+RATE(2445)*D*Y(13)*Y(333)
    PROD = RATE(106)*Y(11)+RATE(159)*Y(11)
    YDOT(13) = PROD-LOSS
    LOSS = RATE(46)*Y(14)*bulkLayersReciprocal+RATE(1179)*Y(14)
    PROD = RATE(1013)*Y(12)*totalSwap/safeMantle
    YDOT(14) = PROD-LOSS
    LOSS = RATE(543)*D*Y(15)+RATE(1818)*D*Y(2)*Y(15)+RATE(1970)*D*Y(6)&
    &*Y(15)+RATE(2446)*D*Y(15)*Y(333)
    PROD = RATE(1943)*D*Y(4)*Y(11)+RATE(2045)*D*Y(8)*Y(11)+RATE(2393)*D&
    &*Y(13)*Y(116)
    YDOT(15) = PROD-LOSS
    LOSS = RATE(99)*Y(16)+RATE(110)*Y(16)+RATE(454)*D*Y(16)+RATE(830)&
    &*Y(16)+RATE(1217)*D*Y(16)*Y(68)+RATE(1218)*D*Y(16)*Y(83)+RATE(1219)*D&
    &*Y(16)*Y(100)+RATE(1220)*D*Y(16)*Y(104)+RATE(1221)*D*Y(16)*Y(162)&
    &+RATE(1222)*D*Y(16)*Y(76)+RATE(1223)*D*Y(16)*Y(21)+RATE(1224)*D*Y(16)&
    &*Y(25)+RATE(1225)*D*Y(16)*Y(33)+RATE(1226)*D*Y(16)*Y(53)+RATE(1227)*D&
    &*Y(16)*Y(63)+RATE(1228)*D*Y(16)*Y(184)+RATE(1229)*D*Y(16)*Y(66)&
    &+RATE(1230)*D*Y(16)*Y(91)+RATE(1231)*D*Y(16)*Y(117)+RATE(1232)*D*Y(16)&
    &*Y(243)+RATE(1233)*D*Y(16)*Y(145)+RATE(1234)*D*Y(16)*Y(174)+RATE(1235)*D&
    &*Y(16)*Y(120)+RATE(1236)*D*Y(16)*Y(36)+RATE(1237)*D*Y(16)*Y(162)&
    &+RATE(1238)*D*Y(16)*Y(176)+RATE(1239)*D*Y(16)*Y(57)+RATE(1240)*D*Y(16)&
    &*Y(122)+RATE(1241)*D*Y(16)*Y(236)+RATE(1242)*D*Y(16)*Y(89)+RATE(1243)*D&
    &*Y(16)*Y(123)+RATE(1244)*D*Y(16)*Y(197)+RATE(1245)*D*Y(16)*Y(200)&
    &+RATE(1246)*D*Y(16)*Y(24)+RATE(1247)*D*Y(16)*Y(24)+RATE(1248)*D*Y(16)&
    &*Y(32)+RATE(1249)*D*Y(16)*Y(20)+RATE(1250)*D*Y(16)*Y(82)+RATE(1251)*D&
    &*Y(16)*Y(99)+RATE(1252)*D*Y(16)*Y(233)+RATE(1253)*D*Y(16)*Y(101)&
    &+RATE(1254)*D*Y(16)*Y(116)+RATE(1255)*D*Y(16)*Y(173)+RATE(1256)*D*Y(16)&
    &*Y(173)+RATE(1257)*D*Y(16)*Y(103)+RATE(1258)*D*Y(16)*Y(291)+RATE(1259)*D&
    &*Y(16)*Y(43)+RATE(1260)*D*Y(16)*Y(43)+RATE(1261)*D*Y(16)*Y(43)+RATE(1262&
    &)*D*Y(16)*Y(35)+RATE(1263)*D*Y(16)*Y(35)+RATE(1264)*D*Y(16)*Y(133)&
    &+RATE(1265)*D*Y(16)*Y(133)+RATE(1266)*D*Y(16)*Y(265)+RATE(1267)*D*Y(16)&
    &*Y(265)+RATE(1268)*D*Y(16)*Y(161)+RATE(1269)*D*Y(16)*Y(221)+RATE(1270)*D&
    &*Y(16)*Y(300)+RATE(1271)*D*Y(16)*Y(56)+RATE(1272)*D*Y(16)*Y(56)&
    &+RATE(1273)*D*Y(16)*Y(316)+RATE(1274)*D*Y(16)*Y(320)+RATE(1275)*D*Y(16)&
    &*Y(277)+RATE(1276)*D*Y(16)*Y(277)+RATE(1277)*D*Y(16)*Y(121)+RATE(1278)*D&
    &*Y(16)*Y(16)+RATE(1278)*D*Y(16)*Y(16)+RATE(1279)*D*Y(16)*Y(27)+RATE(1280&
    &)*D*Y(16)*Y(48)+RATE(1281)*D*Y(16)*Y(46)+RATE(1282)*D*Y(16)*Y(165)&
    &+RATE(1283)*D*Y(16)*Y(163)+RATE(1343)*D*Y(18)*Y(16)+RATE(1867)*D*Y(2)&
    &*Y(16)+RATE(1985)*D*Y(6)*Y(16)+RATE(2000)*D*Y(6)*Y(16)+RATE(2031)*D*Y(8)&
    &*Y(16)+RATE(2152)*D*Y(10)*Y(16)+RATE(2322)*D*Y(13)*Y(16)+RATE(2455)*D&
    &*Y(227)*Y(16)
    PROD = RATE(111)*Y(67)+RATE(111)*Y(67)+RATE(120)*Y(197)+RATE(123)&
    &*Y(20)+RATE(138)*Y(82)+RATE(139)*Y(99)+RATE(142)*Y(233)+RATE(191)*Y(207)&
    &+RATE(192)*Y(292)+RATE(193)*Y(314)+RATE(203)*Y(17)/safeMantle+RATE(286)&
    &*D*Y(17)/safeMantle*Y(2)+RATE(369)*Y(17)/safeMantle+RATE(469)*D*Y(192)&
    &+RATE(832)*Y(67)+RATE(832)*Y(67)+RATE(833)*Y(68)+RATE(843)*Y(197)&
    &+RATE(846)*Y(20)+RATE(848)*Y(21)+RATE(873)*Y(82)+RATE(874)*Y(99)&
    &+RATE(878)*Y(233)+RATE(879)*Y(234)+RATE(947)*Y(207)+RATE(950)*Y(314)&
    &+RATE(1051)*Y(17)+RATE(1134)*Y(19)+RATE(1284)*D*Y(18)*Y(24)+RATE(1285)*D&
    &*Y(18)*Y(20)+RATE(1286)*D*Y(18)*Y(131)+RATE(1287)*D*Y(18)*Y(183)&
    &+RATE(1288)*D*Y(18)*Y(116)+RATE(1289)*D*Y(18)*Y(69)+RATE(1291)*D*Y(18)&
    &*Y(54)+RATE(1292)*D*Y(18)*Y(133)+RATE(1293)*D*Y(18)*Y(265)+RATE(1294)*D&
    &*Y(18)*Y(300)+RATE(1295)*D*Y(18)*Y(277)+RATE(1296)*D*Y(18)*Y(105)&
    &+RATE(1297)*D*Y(18)*Y(292)+RATE(1298)*D*Y(18)*Y(314)+RATE(1299)*D*Y(18)&
    &*Y(207)+RATE(1300)*D*Y(18)*Y(135)+RATE(1301)*D*Y(18)*Y(146)+RATE(1302)*D&
    &*Y(18)*Y(302)+RATE(1347)*D*Y(18)*Y(333)+RATE(1359)*D*Y(67)*Y(165)&
    &+RATE(1364)*D*Y(67)*Y(163)+RATE(1368)*D*Y(68)*Y(333)+RATE(1368)*D*Y(68)&
    &*Y(333)+RATE(1369)*D*Y(68)*Y(67)+RATE(1372)*D*Y(68)*Y(163)+RATE(1392)*D&
    &*Y(76)*Y(333)+RATE(1423)*D*Y(199)*Y(333)+RATE(1427)*D*Y(192)*Y(333)&
    &+RATE(1430)*D*Y(309)*Y(333)+RATE(1447)*D*Y(20)*Y(100)+RATE(1462)*D*Y(20)&
    &*Y(55)+RATE(1481)*D*Y(20)*Y(27)+RATE(1492)*D*Y(20)*Y(46)+RATE(1496)*D&
    &*Y(20)*Y(163)+RATE(1505)*D*Y(21)*Y(333)+RATE(1509)*D*Y(21)*Y(159)&
    &+RATE(1517)*D*Y(21)*Y(131)+RATE(1520)*D*Y(21)*Y(62)+RATE(1522)*D*Y(21)&
    &*Y(183)+RATE(1526)*D*Y(21)*Y(90)+RATE(1528)*D*Y(21)*Y(92)+RATE(1531)*D&
    &*Y(21)*Y(54)+RATE(1538)*D*Y(21)*Y(300)+RATE(1541)*D*Y(21)*Y(163)&
    &+RATE(1601)*D*Y(25)*Y(333)+RATE(1602)*D*Y(25)*Y(333)+RATE(1760)*D*Y(82)&
    &*Y(163)+RATE(1769)*D*Y(83)*Y(333)+RATE(1793)*D*Y(100)*Y(333)+RATE(1797)&
    &*D*Y(234)*Y(333)+RATE(1798)*D*Y(2)*Y(20)+RATE(1822)*D*Y(2)*Y(67)&
    &+RATE(1829)*D*Y(2)*Y(20)+RATE(1831)*D*Y(2)*Y(99)+RATE(1945)*D*Y(6)*Y(20)&
    &+RATE(2333)*D*Y(13)*Y(67)+RATE(2343)*D*Y(13)*Y(75)+RATE(2364)*D*Y(13)&
    &*Y(82)+RATE(2368)*D*Y(13)*Y(232)+RATE(2371)*D*Y(13)*Y(233)+RATE(2400)*D&
    &*Y(13)*Y(92)+RATE(2431)*D*Y(13)*Y(314)+RATE(2432)*D*Y(13)*Y(207)&
    &+RATE(2497)*D*Y(27)*Y(83)+RATE(2510)*D*Y(27)*Y(67)+RATE(2527)*D*Y(27)&
    &*Y(82)+RATE(2580)*D*Y(28)*Y(99)+RATE(2670)*D*Y(36)*Y(67)+RATE(2827)*D&
    &*Y(46)*Y(68)+RATE(2828)*D*Y(46)*Y(76)+RATE(2849)*D*Y(46)*Y(208)&
    &+RATE(2854)*D*Y(46)*Y(67)+RATE(2869)*D*Y(46)*Y(82)+RATE(2872)*D*Y(46)&
    &*Y(233)+RATE(2909)*D*Y(46)*Y(207)+RATE(2934)*D*Y(48)*Y(67)+RATE(2940)*D&
    &*Y(48)*Y(82)+RATE(2974)*D*Y(301)*Y(333)+RATE(3051)*D*Y(163)*Y(18)&
    &+RATE(3078)*D*Y(105)*Y(99)+RATE(3084)*D*Y(208)*Y(333)+RATE(3086)*D*Y(293&
    &)*Y(333)+RATE(3087)*D*Y(313)*Y(333)
    YDOT(16) = PROD-LOSS
    LOSS = RATE(203)*Y(17)/safeMantle+RATE(286)*D*Y(17)/safeMantle*Y(2)&
    &+RATE(369)*Y(17)/safeMantle+RATE(621)*Y(17)*Y(222)+RATE(629)*Y(3)*Y(17)&
    &+RATE(725)*Y(17)*Y(222)+RATE(733)*Y(3)*Y(17)+RATE(968)*Y(17)&
    &*totalSwap/safeMantle+RATE(1051)*Y(17)
    PROD = RATE(1)*Y(19)*bulkLayersReciprocal+RATE(454)*D*Y(16)+RATE(455&
    &)*D*Y(18)
    YDOT(17) = PROD-LOSS
    LOSS = RATE(455)*D*Y(18)+RATE(1284)*D*Y(18)*Y(24)+RATE(1285)*D*Y(18)&
    &*Y(20)+RATE(1286)*D*Y(18)*Y(131)+RATE(1287)*D*Y(18)*Y(183)+RATE(1288)*D&
    &*Y(18)*Y(116)+RATE(1289)*D*Y(18)*Y(69)+RATE(1290)*D*Y(18)*Y(291)&
    &+RATE(1291)*D*Y(18)*Y(54)+RATE(1292)*D*Y(18)*Y(133)+RATE(1293)*D*Y(18)&
    &*Y(265)+RATE(1294)*D*Y(18)*Y(300)+RATE(1295)*D*Y(18)*Y(277)+RATE(1296)*D&
    &*Y(18)*Y(105)+RATE(1297)*D*Y(18)*Y(292)+RATE(1298)*D*Y(18)*Y(314)&
    &+RATE(1299)*D*Y(18)*Y(207)+RATE(1300)*D*Y(18)*Y(135)+RATE(1301)*D*Y(18)&
    &*Y(146)+RATE(1302)*D*Y(18)*Y(302)+RATE(1303)*D*Y(18)*Y(75)+RATE(1304)*D&
    &*Y(18)*Y(24)+RATE(1305)*D*Y(18)*Y(32)+RATE(1306)*D*Y(18)*Y(32)+RATE(1307&
    &)*D*Y(18)*Y(206)+RATE(1308)*D*Y(18)*Y(159)+RATE(1309)*D*Y(18)*Y(159)&
    &+RATE(1310)*D*Y(18)*Y(41)+RATE(1311)*D*Y(18)*Y(20)+RATE(1312)*D*Y(18)&
    &*Y(232)+RATE(1313)*D*Y(18)*Y(131)+RATE(1314)*D*Y(18)*Y(131)+RATE(1315)*D&
    &*Y(18)*Y(261)+RATE(1316)*D*Y(18)*Y(62)+RATE(1317)*D*Y(18)*Y(62)&
    &+RATE(1318)*D*Y(18)*Y(183)+RATE(1319)*D*Y(18)*Y(288)+RATE(1320)*D*Y(18)&
    &*Y(288)+RATE(1321)*D*Y(18)*Y(288)+RATE(1322)*D*Y(18)*Y(116)+RATE(1323)*D&
    &*Y(18)*Y(92)+RATE(1324)*D*Y(18)*Y(173)+RATE(1325)*D*Y(18)*Y(43)&
    &+RATE(1326)*D*Y(18)*Y(54)+RATE(1327)*D*Y(18)*Y(35)+RATE(1328)*D*Y(18)&
    &*Y(265)+RATE(1329)*D*Y(18)*Y(161)+RATE(1330)*D*Y(18)*Y(161)+RATE(1331)*D&
    &*Y(18)*Y(221)+RATE(1332)*D*Y(18)*Y(300)+RATE(1333)*D*Y(18)*Y(56)&
    &+RATE(1334)*D*Y(18)*Y(320)+RATE(1335)*D*Y(18)*Y(277)+RATE(1336)*D*Y(18)&
    &*Y(277)+RATE(1337)*D*Y(18)*Y(277)+RATE(1338)*D*Y(18)*Y(207)+RATE(1339)*D&
    &*Y(18)*Y(135)+RATE(1340)*D*Y(18)*Y(121)+RATE(1341)*D*Y(18)*Y(235)&
    &+RATE(1342)*D*Y(18)*Y(302)+RATE(1343)*D*Y(18)*Y(16)+RATE(1344)*D*Y(18)&
    &*Y(27)+RATE(1345)*D*Y(18)*Y(46)+RATE(1346)*D*Y(18)*Y(163)+RATE(1347)*D&
    &*Y(18)*Y(333)+RATE(1866)*D*Y(2)*Y(18)+RATE(1953)*D*Y(6)*Y(18)+RATE(1999)&
    &*D*Y(6)*Y(18)+RATE(3051)*D*Y(163)*Y(18)
    PROD = RATE(99)*Y(16)+RATE(110)*Y(16)+RATE(124)*Y(21)+RATE(830)*Y(16&
    &)+RATE(833)*Y(68)+RATE(851)*Y(25)+RATE(875)*Y(100)+RATE(1217)*D*Y(16)&
    &*Y(68)+RATE(1218)*D*Y(16)*Y(83)+RATE(1219)*D*Y(16)*Y(100)+RATE(1220)*D&
    &*Y(16)*Y(104)+RATE(1221)*D*Y(16)*Y(162)+RATE(1810)*D*Y(2)*Y(21)&
    &+RATE(2322)*D*Y(13)*Y(16)+RATE(2333)*D*Y(13)*Y(67)+RATE(2344)*D*Y(13)&
    &*Y(75)+RATE(2345)*D*Y(13)*Y(197)+RATE(2348)*D*Y(13)*Y(24)+RATE(2363)*D&
    &*Y(13)*Y(20)+RATE(2365)*D*Y(13)*Y(82)+RATE(2369)*D*Y(13)*Y(232)&
    &+RATE(2370)*D*Y(13)*Y(99)+RATE(2372)*D*Y(13)*Y(233)+RATE(2390)*D*Y(13)&
    &*Y(90)+RATE(2399)*D*Y(13)*Y(92)+RATE(2433)*D*Y(13)*Y(207)+RATE(2490)*D&
    &*Y(27)*Y(68)
    YDOT(18) = PROD-LOSS
    LOSS = RATE(1)*Y(19)*bulkLayersReciprocal+RATE(673)*Y(19)*Y(225)&
    &*bulkLayersReciprocal+RATE(681)*Y(5)*Y(19)*bulkLayersReciprocal+RATE(777&
    &)*Y(19)*Y(225)*bulkLayersReciprocal+RATE(785)*Y(5)*Y(19)&
    &*bulkLayersReciprocal+RATE(1134)*Y(19)
    PROD = RATE(968)*Y(17)*totalSwap/safeMantle
    YDOT(19) = PROD-LOSS
    LOSS = RATE(123)*Y(20)+RATE(476)*D*Y(20)+RATE(846)*Y(20)+RATE(847)&
    &*Y(20)+RATE(1249)*D*Y(16)*Y(20)+RATE(1285)*D*Y(18)*Y(20)+RATE(1311)*D&
    &*Y(18)*Y(20)+RATE(1431)*D*Y(20)*Y(46)+RATE(1432)*D*Y(20)*Y(68)+RATE(1433&
    &)*D*Y(20)*Y(83)+RATE(1434)*D*Y(20)*Y(100)+RATE(1435)*D*Y(20)*Y(132)&
    &+RATE(1436)*D*Y(20)*Y(63)+RATE(1437)*D*Y(20)*Y(28)+RATE(1438)*D*Y(20)&
    &*Y(104)+RATE(1439)*D*Y(20)*Y(44)+RATE(1440)*D*Y(20)*Y(48)+RATE(1441)*D&
    &*Y(20)*Y(162)+RATE(1442)*D*Y(20)*Y(57)+RATE(1443)*D*Y(20)*Y(68)&
    &+RATE(1444)*D*Y(20)*Y(76)+RATE(1445)*D*Y(20)*Y(33)+RATE(1446)*D*Y(20)&
    &*Y(53)+RATE(1447)*D*Y(20)*Y(100)+RATE(1448)*D*Y(20)*Y(132)+RATE(1449)*D&
    &*Y(20)*Y(63)+RATE(1450)*D*Y(20)*Y(142)+RATE(1451)*D*Y(20)*Y(66)&
    &+RATE(1452)*D*Y(20)*Y(91)+RATE(1453)*D*Y(20)*Y(102)+RATE(1454)*D*Y(20)&
    &*Y(102)+RATE(1455)*D*Y(20)*Y(117)+RATE(1456)*D*Y(20)*Y(145)+RATE(1457)*D&
    &*Y(20)*Y(174)+RATE(1458)*D*Y(20)*Y(28)+RATE(1459)*D*Y(20)*Y(120)&
    &+RATE(1460)*D*Y(20)*Y(36)+RATE(1461)*D*Y(20)*Y(44)+RATE(1462)*D*Y(20)&
    &*Y(55)+RATE(1463)*D*Y(20)*Y(48)+RATE(1464)*D*Y(20)*Y(162)+RATE(1465)*D&
    &*Y(20)*Y(176)+RATE(1466)*D*Y(20)*Y(57)+RATE(1467)*D*Y(20)*Y(165)&
    &+RATE(1468)*D*Y(20)*Y(106)+RATE(1469)*D*Y(20)*Y(122)+RATE(1470)*D*Y(20)&
    &*Y(236)+RATE(1471)*D*Y(20)*Y(80)+RATE(1472)*D*Y(20)*Y(109)+RATE(1473)*D&
    &*Y(20)*Y(159)+RATE(1474)*D*Y(20)*Y(41)+RATE(1475)*D*Y(20)*Y(232)&
    &+RATE(1476)*D*Y(20)*Y(131)+RATE(1477)*D*Y(20)*Y(116)+RATE(1478)*D*Y(20)&
    &*Y(144)+RATE(1479)*D*Y(20)*Y(103)+RATE(1480)*D*Y(20)*Y(27)+RATE(1481)*D&
    &*Y(20)*Y(27)+RATE(1482)*D*Y(20)*Y(133)+RATE(1483)*D*Y(20)*Y(133)&
    &+RATE(1484)*D*Y(20)*Y(133)+RATE(1485)*D*Y(20)*Y(161)+RATE(1486)*D*Y(20)&
    &*Y(161)+RATE(1487)*D*Y(20)*Y(161)+RATE(1488)*D*Y(20)*Y(161)+RATE(1489)*D&
    &*Y(20)*Y(175)+RATE(1490)*D*Y(20)*Y(175)+RATE(1491)*D*Y(20)*Y(46)&
    &+RATE(1492)*D*Y(20)*Y(46)+RATE(1493)*D*Y(20)*Y(300)+RATE(1494)*D*Y(20)&
    &*Y(56)+RATE(1495)*D*Y(20)*Y(163)+RATE(1496)*D*Y(20)*Y(163)+RATE(1497)*D&
    &*Y(20)*Y(277)+RATE(1498)*D*Y(20)*Y(277)+RATE(1513)*D*Y(21)*Y(20)&
    &+RATE(1798)*D*Y(2)*Y(20)+RATE(1829)*D*Y(2)*Y(20)+RATE(1879)*D*Y(4)*Y(20)&
    &+RATE(1945)*D*Y(6)*Y(20)+RATE(1988)*D*Y(6)*Y(20)+RATE(2002)*D*Y(6)*Y(20)&
    &+RATE(2013)*D*Y(8)*Y(20)+RATE(2035)*D*Y(8)*Y(20)+RATE(2162)*D*Y(10)*Y(20&
    &)+RATE(2324)*D*Y(13)*Y(20)+RATE(2363)*D*Y(13)*Y(20)
    PROD = RATE(126)*Y(24)+RATE(130)*Y(32)+RATE(215)*Y(22)/safeMantle&
    &+RATE(298)*D*Y(22)/safeMantle*Y(2)+RATE(381)*Y(22)/safeMantle+RATE(733)&
    &*Y(3)*Y(17)+RATE(785)*Y(5)*Y(19)*bulkLayersReciprocal+RATE(850)*Y(24)&
    &+RATE(853)*Y(25)+RATE(857)*Y(32)+RATE(869)*Y(41)+RATE(1063)*Y(22)&
    &+RATE(1146)*Y(23)+RATE(1247)*D*Y(16)*Y(24)+RATE(1247)*D*Y(16)*Y(24)&
    &+RATE(1254)*D*Y(16)*Y(116)+RATE(1256)*D*Y(16)*Y(173)+RATE(1261)*D*Y(16)&
    &*Y(43)+RATE(1263)*D*Y(16)*Y(35)+RATE(1272)*D*Y(16)*Y(56)+RATE(1308)*D&
    &*Y(18)*Y(159)+RATE(1314)*D*Y(18)*Y(131)+RATE(1392)*D*Y(76)*Y(333)&
    &+RATE(1408)*D*Y(81)*Y(333)+RATE(1408)*D*Y(81)*Y(333)+RATE(1499)*D*Y(21)&
    &*Y(116)+RATE(1500)*D*Y(21)*Y(69)+RATE(1501)*D*Y(21)*Y(54)+RATE(1502)*D&
    &*Y(21)*Y(133)+RATE(1503)*D*Y(21)*Y(163)+RATE(1504)*D*Y(21)*Y(105)&
    &+RATE(1554)*D*Y(24)*Y(100)+RATE(1575)*D*Y(24)*Y(24)+RATE(1577)*D*Y(24)&
    &*Y(82)+RATE(1594)*D*Y(24)*Y(46)+RATE(1596)*D*Y(24)*Y(56)+RATE(1603)*D&
    &*Y(25)*Y(333)+RATE(1608)*D*Y(25)*Y(183)+RATE(1611)*D*Y(25)*Y(54)&
    &+RATE(1647)*D*Y(33)*Y(333)+RATE(1648)*D*Y(33)*Y(333)+RATE(1725)*D*Y(53)&
    &*Y(333)+RATE(1822)*D*Y(2)*Y(67)+RATE(1825)*D*Y(2)*Y(24)+RATE(1867)*D*Y(2&
    &)*Y(16)+RATE(1985)*D*Y(6)*Y(16)+RATE(2208)*D*Y(142)*Y(333)+RATE(2318)*D&
    &*Y(245)*Y(333)+RATE(2336)*D*Y(13)*Y(80)+RATE(2344)*D*Y(13)*Y(75)&
    &+RATE(2386)*D*Y(13)*Y(288)+RATE(2389)*D*Y(13)*Y(90)+RATE(2523)*D*Y(27)&
    &*Y(24)+RATE(2829)*D*Y(46)*Y(81)+RATE(2863)*D*Y(46)*Y(75)+RATE(2936)*D&
    &*Y(48)*Y(75)+RATE(2946)*D*Y(48)*Y(90)
    YDOT(20) = PROD-LOSS
    LOSS = RATE(124)*Y(21)+RATE(477)*D*Y(21)+RATE(848)*Y(21)+RATE(1223)&
    &*D*Y(16)*Y(21)+RATE(1499)*D*Y(21)*Y(116)+RATE(1500)*D*Y(21)*Y(69)&
    &+RATE(1501)*D*Y(21)*Y(54)+RATE(1502)*D*Y(21)*Y(133)+RATE(1503)*D*Y(21)&
    &*Y(163)+RATE(1504)*D*Y(21)*Y(105)+RATE(1505)*D*Y(21)*Y(333)+RATE(1506)*D&
    &*Y(21)*Y(67)+RATE(1507)*D*Y(21)*Y(75)+RATE(1508)*D*Y(21)*Y(24)+RATE(1509&
    &)*D*Y(21)*Y(159)+RATE(1510)*D*Y(21)*Y(159)+RATE(1511)*D*Y(21)*Y(159)&
    &+RATE(1512)*D*Y(21)*Y(41)+RATE(1513)*D*Y(21)*Y(20)+RATE(1514)*D*Y(21)&
    &*Y(82)+RATE(1515)*D*Y(21)*Y(232)+RATE(1516)*D*Y(21)*Y(131)+RATE(1517)*D&
    &*Y(21)*Y(131)+RATE(1518)*D*Y(21)*Y(131)+RATE(1519)*D*Y(21)*Y(62)&
    &+RATE(1520)*D*Y(21)*Y(62)+RATE(1521)*D*Y(21)*Y(62)+RATE(1522)*D*Y(21)&
    &*Y(183)+RATE(1523)*D*Y(21)*Y(183)+RATE(1524)*D*Y(21)*Y(90)+RATE(1525)*D&
    &*Y(21)*Y(90)+RATE(1526)*D*Y(21)*Y(90)+RATE(1527)*D*Y(21)*Y(116)&
    &+RATE(1528)*D*Y(21)*Y(92)+RATE(1529)*D*Y(21)*Y(27)+RATE(1530)*D*Y(21)&
    &*Y(43)+RATE(1531)*D*Y(21)*Y(54)+RATE(1532)*D*Y(21)*Y(35)+RATE(1533)*D&
    &*Y(21)*Y(161)+RATE(1534)*D*Y(21)*Y(161)+RATE(1535)*D*Y(21)*Y(161)&
    &+RATE(1536)*D*Y(21)*Y(46)+RATE(1537)*D*Y(21)*Y(300)+RATE(1538)*D*Y(21)&
    &*Y(300)+RATE(1539)*D*Y(21)*Y(56)+RATE(1540)*D*Y(21)*Y(163)+RATE(1541)*D&
    &*Y(21)*Y(163)+RATE(1810)*D*Y(2)*Y(21)+RATE(1956)*D*Y(6)*Y(21)
    PROD = RATE(847)*Y(20)+RATE(852)*Y(25)+RATE(858)*Y(33)+RATE(1226)*D&
    &*Y(16)*Y(53)+RATE(1227)*D*Y(16)*Y(63)+RATE(1230)*D*Y(16)*Y(91)+RATE(1231&
    &)*D*Y(16)*Y(117)+RATE(1232)*D*Y(16)*Y(243)+RATE(1233)*D*Y(16)*Y(145)&
    &+RATE(1235)*D*Y(16)*Y(120)+RATE(1236)*D*Y(16)*Y(36)+RATE(1238)*D*Y(16)&
    &*Y(176)+RATE(1239)*D*Y(16)*Y(57)+RATE(1285)*D*Y(18)*Y(20)+RATE(1322)*D&
    &*Y(18)*Y(116)+RATE(1432)*D*Y(20)*Y(68)+RATE(1433)*D*Y(20)*Y(83)&
    &+RATE(1434)*D*Y(20)*Y(100)+RATE(1435)*D*Y(20)*Y(132)+RATE(1436)*D*Y(20)&
    &*Y(63)+RATE(1437)*D*Y(20)*Y(28)+RATE(1438)*D*Y(20)*Y(104)+RATE(1439)*D&
    &*Y(20)*Y(44)+RATE(1440)*D*Y(20)*Y(48)+RATE(1441)*D*Y(20)*Y(162)&
    &+RATE(1442)*D*Y(20)*Y(57)+RATE(1811)*D*Y(2)*Y(25)+RATE(1866)*D*Y(2)*Y(18&
    &)+RATE(1879)*D*Y(4)*Y(20)+RATE(1918)*D*Y(4)*Y(24)+RATE(1953)*D*Y(6)*Y(18&
    &)+RATE(2013)*D*Y(8)*Y(20)+RATE(2031)*D*Y(8)*Y(16)+RATE(2152)*D*Y(10)&
    &*Y(16)+RATE(2324)*D*Y(13)*Y(20)+RATE(2336)*D*Y(13)*Y(80)+RATE(2343)*D&
    &*Y(13)*Y(75)+RATE(2349)*D*Y(13)*Y(24)+RATE(2352)*D*Y(13)*Y(32)+RATE(2359&
    &)*D*Y(13)*Y(41)+RATE(2391)*D*Y(13)*Y(90)+RATE(2394)*D*Y(13)*Y(116)&
    &+RATE(2492)*D*Y(27)*Y(76)+RATE(2495)*D*Y(27)*Y(81)
    YDOT(21) = PROD-LOSS
    LOSS = RATE(215)*Y(22)/safeMantle+RATE(298)*D*Y(22)/safeMantle*Y(2)&
    &+RATE(381)*Y(22)/safeMantle+RATE(630)*Y(3)*Y(22)+RATE(734)*Y(3)*Y(22)&
    &+RATE(980)*Y(22)*totalSwap/safeMantle+RATE(1063)*Y(22)
    PROD = RATE(13)*Y(23)*bulkLayersReciprocal+RATE(476)*D*Y(20)&
    &+RATE(477)*D*Y(21)+RATE(629)*Y(3)*Y(17)
    YDOT(22) = PROD-LOSS
    LOSS = RATE(13)*Y(23)*bulkLayersReciprocal+RATE(682)*Y(5)*Y(23)&
    &*bulkLayersReciprocal+RATE(786)*Y(5)*Y(23)*bulkLayersReciprocal&
    &+RATE(1146)*Y(23)
    PROD = RATE(681)*Y(5)*Y(19)*bulkLayersReciprocal+RATE(980)*Y(22)&
    &*totalSwap/safeMantle
    YDOT(23) = PROD-LOSS
    LOSS = RATE(125)*Y(24)+RATE(126)*Y(24)+RATE(478)*D*Y(24)+RATE(849)&
    &*Y(24)+RATE(850)*Y(24)+RATE(1246)*D*Y(16)*Y(24)+RATE(1247)*D*Y(16)*Y(24)&
    &+RATE(1284)*D*Y(18)*Y(24)+RATE(1304)*D*Y(18)*Y(24)+RATE(1508)*D*Y(21)&
    &*Y(24)+RATE(1542)*D*Y(24)*Y(68)+RATE(1543)*D*Y(24)*Y(83)+RATE(1544)*D&
    &*Y(24)*Y(100)+RATE(1545)*D*Y(24)*Y(132)+RATE(1546)*D*Y(24)*Y(63)&
    &+RATE(1547)*D*Y(24)*Y(104)+RATE(1548)*D*Y(24)*Y(44)+RATE(1549)*D*Y(24)&
    &*Y(48)+RATE(1550)*D*Y(24)*Y(162)+RATE(1551)*D*Y(24)*Y(57)+RATE(1552)*D&
    &*Y(24)*Y(76)+RATE(1553)*D*Y(24)*Y(53)+RATE(1554)*D*Y(24)*Y(100)&
    &+RATE(1555)*D*Y(24)*Y(132)+RATE(1556)*D*Y(24)*Y(63)+RATE(1557)*D*Y(24)&
    &*Y(66)+RATE(1558)*D*Y(24)*Y(91)+RATE(1559)*D*Y(24)*Y(102)+RATE(1560)*D&
    &*Y(24)*Y(102)+RATE(1561)*D*Y(24)*Y(117)+RATE(1562)*D*Y(24)*Y(145)&
    &+RATE(1563)*D*Y(24)*Y(120)+RATE(1564)*D*Y(24)*Y(36)+RATE(1565)*D*Y(24)&
    &*Y(44)+RATE(1566)*D*Y(24)*Y(55)+RATE(1567)*D*Y(24)*Y(162)+RATE(1568)*D&
    &*Y(24)*Y(176)+RATE(1569)*D*Y(24)*Y(57)+RATE(1570)*D*Y(24)*Y(165)&
    &+RATE(1571)*D*Y(24)*Y(236)+RATE(1572)*D*Y(24)*Y(24)+RATE(1572)*D*Y(24)&
    &*Y(24)+RATE(1573)*D*Y(24)*Y(24)+RATE(1573)*D*Y(24)*Y(24)+RATE(1574)*D&
    &*Y(24)*Y(24)+RATE(1574)*D*Y(24)*Y(24)+RATE(1575)*D*Y(24)*Y(24)+RATE(1575&
    &)*D*Y(24)*Y(24)+RATE(1576)*D*Y(24)*Y(41)+RATE(1577)*D*Y(24)*Y(82)&
    &+RATE(1578)*D*Y(24)*Y(131)+RATE(1579)*D*Y(24)*Y(116)+RATE(1580)*D*Y(24)&
    &*Y(144)+RATE(1581)*D*Y(24)*Y(103)+RATE(1582)*D*Y(24)*Y(264)+RATE(1583)*D&
    &*Y(24)*Y(133)+RATE(1584)*D*Y(24)*Y(133)+RATE(1585)*D*Y(24)*Y(133)&
    &+RATE(1586)*D*Y(24)*Y(161)+RATE(1587)*D*Y(24)*Y(161)+RATE(1588)*D*Y(24)&
    &*Y(161)+RATE(1589)*D*Y(24)*Y(161)+RATE(1590)*D*Y(24)*Y(161)+RATE(1591)*D&
    &*Y(24)*Y(46)+RATE(1592)*D*Y(24)*Y(46)+RATE(1593)*D*Y(24)*Y(46)+RATE(1594&
    &)*D*Y(24)*Y(46)+RATE(1595)*D*Y(24)*Y(56)+RATE(1596)*D*Y(24)*Y(56)&
    &+RATE(1597)*D*Y(24)*Y(56)+RATE(1598)*D*Y(24)*Y(163)+RATE(1599)*D*Y(24)&
    &*Y(163)+RATE(1825)*D*Y(2)*Y(24)+RATE(1876)*D*Y(4)*Y(24)+RATE(1918)*D*Y(4&
    &)*Y(24)+RATE(1986)*D*Y(6)*Y(24)+RATE(2011)*D*Y(8)*Y(24)+RATE(2032)*D*Y(8&
    &)*Y(24)+RATE(2153)*D*Y(10)*Y(24)+RATE(2348)*D*Y(13)*Y(24)+RATE(2349)*D&
    &*Y(13)*Y(24)+RATE(2521)*D*Y(27)*Y(24)+RATE(2522)*D*Y(27)*Y(24)+RATE(2523&
    &)*D*Y(27)*Y(24)+RATE(2555)*D*Y(28)*Y(24)
    PROD = RATE(127)*Y(218)+RATE(128)*Y(32)+RATE(136)*Y(41)+RATE(216)&
    &*Y(26)/safeMantle+RATE(299)*D*Y(26)/safeMantle*Y(2)+RATE(382)*Y(26&
    &)/safeMantle+RATE(734)*Y(3)*Y(22)+RATE(786)*Y(5)*Y(23)&
    &*bulkLayersReciprocal+RATE(854)*Y(218)+RATE(855)*Y(32)+RATE(866)*Y(41)&
    &+RATE(1064)*Y(26)+RATE(1147)*Y(30)+RATE(1476)*D*Y(20)*Y(131)+RATE(1477)&
    &*D*Y(20)*Y(116)+RATE(1478)*D*Y(20)*Y(144)+RATE(1490)*D*Y(20)*Y(175)&
    &+RATE(1511)*D*Y(21)*Y(159)+RATE(1518)*D*Y(21)*Y(131)+RATE(1600)*D*Y(25)&
    &*Y(133)+RATE(1621)*D*Y(32)*Y(32)+RATE(1622)*D*Y(32)*Y(82)+RATE(1634)*D&
    &*Y(32)*Y(161)+RATE(1640)*D*Y(32)*Y(56)+RATE(1646)*D*Y(33)*Y(333)&
    &+RATE(1657)*D*Y(33)*Y(54)+RATE(1673)*D*Y(172)*Y(333)+RATE(1681)*D*Y(41)&
    &*Y(68)+RATE(1697)*D*Y(41)*Y(57)+RATE(1710)*D*Y(42)*Y(333)+RATE(1721)*D&
    &*Y(53)*Y(333)+RATE(1827)*D*Y(2)*Y(32)+RATE(1838)*D*Y(2)*Y(116)+RATE(1988&
    &)*D*Y(6)*Y(20)+RATE(2000)*D*Y(6)*Y(16)+RATE(2063)*D*Y(132)*Y(333)&
    &+RATE(2207)*D*Y(142)*Y(333)+RATE(2341)*D*Y(13)*Y(109)+RATE(2350)*D*Y(13)&
    &*Y(218)+RATE(2377)*D*Y(13)*Y(261)+RATE(2582)*D*Y(28)*Y(131)+RATE(2832)*D&
    &*Y(46)*Y(53)+RATE(2855)*D*Y(46)*Y(80)+RATE(2859)*D*Y(46)*Y(109)
    YDOT(24) = PROD-LOSS
    LOSS = RATE(479)*D*Y(25)+RATE(851)*Y(25)+RATE(852)*Y(25)+RATE(853)&
    &*Y(25)+RATE(1224)*D*Y(16)*Y(25)+RATE(1600)*D*Y(25)*Y(133)+RATE(1601)*D&
    &*Y(25)*Y(333)+RATE(1602)*D*Y(25)*Y(333)+RATE(1603)*D*Y(25)*Y(333)&
    &+RATE(1604)*D*Y(25)*Y(232)+RATE(1605)*D*Y(25)*Y(131)+RATE(1606)*D*Y(25)&
    &*Y(62)+RATE(1607)*D*Y(25)*Y(183)+RATE(1608)*D*Y(25)*Y(183)+RATE(1609)*D&
    &*Y(25)*Y(183)+RATE(1610)*D*Y(25)*Y(116)+RATE(1611)*D*Y(25)*Y(54)&
    &+RATE(1612)*D*Y(25)*Y(161)+RATE(1613)*D*Y(25)*Y(46)+RATE(1614)*D*Y(25)&
    &*Y(300)+RATE(1615)*D*Y(25)*Y(300)+RATE(1616)*D*Y(25)*Y(163)+RATE(1811)*D&
    &*Y(2)*Y(25)+RATE(1957)*D*Y(6)*Y(25)+RATE(2496)*D*Y(27)*Y(25)
    PROD = RATE(125)*Y(24)+RATE(849)*Y(24)+RATE(859)*Y(33)+RATE(870)&
    &*Y(42)+RATE(1284)*D*Y(18)*Y(24)+RATE(1313)*D*Y(18)*Y(131)+RATE(1315)*D&
    &*Y(18)*Y(261)+RATE(1444)*D*Y(20)*Y(76)+RATE(1446)*D*Y(20)*Y(53)&
    &+RATE(1448)*D*Y(20)*Y(132)+RATE(1449)*D*Y(20)*Y(63)+RATE(1450)*D*Y(20)&
    &*Y(142)+RATE(1451)*D*Y(20)*Y(66)+RATE(1452)*D*Y(20)*Y(91)+RATE(1453)*D&
    &*Y(20)*Y(102)+RATE(1454)*D*Y(20)*Y(102)+RATE(1455)*D*Y(20)*Y(117)&
    &+RATE(1456)*D*Y(20)*Y(145)+RATE(1457)*D*Y(20)*Y(174)+RATE(1459)*D*Y(20)&
    &*Y(120)+RATE(1460)*D*Y(20)*Y(36)+RATE(1461)*D*Y(20)*Y(44)+RATE(1465)*D&
    &*Y(20)*Y(176)+RATE(1466)*D*Y(20)*Y(57)+RATE(1469)*D*Y(20)*Y(122)&
    &+RATE(1527)*D*Y(21)*Y(116)+RATE(1542)*D*Y(24)*Y(68)+RATE(1543)*D*Y(24)&
    &*Y(83)+RATE(1544)*D*Y(24)*Y(100)+RATE(1545)*D*Y(24)*Y(132)+RATE(1546)*D&
    &*Y(24)*Y(63)+RATE(1547)*D*Y(24)*Y(104)+RATE(1548)*D*Y(24)*Y(44)&
    &+RATE(1549)*D*Y(24)*Y(48)+RATE(1550)*D*Y(24)*Y(162)+RATE(1551)*D*Y(24)&
    &*Y(57)+RATE(1692)*D*Y(41)*Y(104)+RATE(1812)*D*Y(2)*Y(33)+RATE(1876)*D&
    &*Y(4)*Y(24)+RATE(1956)*D*Y(6)*Y(21)+RATE(1999)*D*Y(6)*Y(18)+RATE(2011)*D&
    &*Y(8)*Y(24)+RATE(2035)*D*Y(8)*Y(20)+RATE(2162)*D*Y(10)*Y(20)+RATE(2341)&
    &*D*Y(13)*Y(109)+RATE(2351)*D*Y(13)*Y(218)+RATE(2360)*D*Y(13)*Y(41)&
    &+RATE(2375)*D*Y(13)*Y(131)+RATE(2378)*D*Y(13)*Y(261)+RATE(2555)*D*Y(28)&
    &*Y(24)
    YDOT(25) = PROD-LOSS
    LOSS = RATE(216)*Y(26)/safeMantle+RATE(299)*D*Y(26)/safeMantle*Y(2)&
    &+RATE(382)*Y(26)/safeMantle+RATE(631)*Y(3)*Y(26)+RATE(735)*Y(3)*Y(26)&
    &+RATE(981)*Y(26)*totalSwap/safeMantle+RATE(1064)*Y(26)
    PROD = RATE(14)*Y(30)*bulkLayersReciprocal+RATE(84)*Y(34)+RATE(88)&
    &*Y(40)+RATE(478)*D*Y(24)+RATE(479)*D*Y(25)+RATE(630)*Y(3)*Y(22)
    YDOT(26) = PROD-LOSS
    LOSS = RATE(107)*Y(27)+RATE(166)*Y(27)+RATE(560)*D*Y(27)+RATE(1279)&
    &*D*Y(16)*Y(27)+RATE(1344)*D*Y(18)*Y(27)+RATE(1480)*D*Y(20)*Y(27)&
    &+RATE(1481)*D*Y(20)*Y(27)+RATE(1529)*D*Y(21)*Y(27)+RATE(1991)*D*Y(6)&
    &*Y(27)+RATE(2047)*D*Y(8)*Y(27)+RATE(2489)*D*Y(27)*Y(104)+RATE(2490)*D&
    &*Y(27)*Y(68)+RATE(2491)*D*Y(27)*Y(76)+RATE(2492)*D*Y(27)*Y(76)+RATE(2493&
    &)*D*Y(27)*Y(81)+RATE(2494)*D*Y(27)*Y(81)+RATE(2495)*D*Y(27)*Y(81)&
    &+RATE(2496)*D*Y(27)*Y(25)+RATE(2497)*D*Y(27)*Y(83)+RATE(2498)*D*Y(27)&
    &*Y(63)+RATE(2499)*D*Y(27)*Y(63)+RATE(2500)*D*Y(27)*Y(184)+RATE(2501)*D&
    &*Y(27)*Y(174)+RATE(2502)*D*Y(27)*Y(36)+RATE(2503)*D*Y(27)*Y(44)&
    &+RATE(2504)*D*Y(27)*Y(162)+RATE(2505)*D*Y(27)*Y(57)+RATE(2506)*D*Y(27)&
    &*Y(278)+RATE(2507)*D*Y(27)*Y(208)+RATE(2508)*D*Y(27)*Y(236)+RATE(2509)*D&
    &*Y(27)*Y(236)+RATE(2510)*D*Y(27)*Y(67)+RATE(2511)*D*Y(27)*Y(89)&
    &+RATE(2512)*D*Y(27)*Y(109)+RATE(2513)*D*Y(27)*Y(123)+RATE(2514)*D*Y(27)&
    &*Y(123)+RATE(2515)*D*Y(27)*Y(75)+RATE(2516)*D*Y(27)*Y(197)+RATE(2517)*D&
    &*Y(27)*Y(200)+RATE(2518)*D*Y(27)*Y(284)+RATE(2519)*D*Y(27)*Y(280)&
    &+RATE(2520)*D*Y(27)*Y(308)+RATE(2521)*D*Y(27)*Y(24)+RATE(2522)*D*Y(27)&
    &*Y(24)+RATE(2523)*D*Y(27)*Y(24)+RATE(2524)*D*Y(27)*Y(32)+RATE(2525)*D&
    &*Y(27)*Y(32)+RATE(2526)*D*Y(27)*Y(32)+RATE(2527)*D*Y(27)*Y(82)+RATE(2528&
    &)*D*Y(27)*Y(232)+RATE(2529)*D*Y(27)*Y(233)+RATE(2530)*D*Y(27)*Y(101)&
    &+RATE(2531)*D*Y(27)*Y(116)+RATE(2532)*D*Y(27)*Y(116)+RATE(2533)*D*Y(27)&
    &*Y(116)+RATE(2534)*D*Y(27)*Y(244)+RATE(2535)*D*Y(27)*Y(144)+RATE(2536)*D&
    &*Y(27)*Y(173)+RATE(2537)*D*Y(27)*Y(173)+RATE(2538)*D*Y(27)*Y(291)&
    &+RATE(2539)*D*Y(27)*Y(35)+RATE(2540)*D*Y(27)*Y(264)+RATE(2541)*D*Y(27)&
    &*Y(264)+RATE(2542)*D*Y(27)*Y(264)+RATE(2543)*D*Y(27)*Y(133)+RATE(2544)*D&
    &*Y(27)*Y(265)+RATE(2545)*D*Y(27)*Y(161)+RATE(2546)*D*Y(27)*Y(175)&
    &+RATE(2547)*D*Y(27)*Y(56)+RATE(2548)*D*Y(27)*Y(56)+RATE(2549)*D*Y(27)&
    &*Y(316)+RATE(2550)*D*Y(27)*Y(277)+RATE(2551)*D*Y(27)*Y(277)+RATE(2552)*D&
    &*Y(27)*Y(207)+RATE(2597)*D*Y(28)*Y(27)
    PROD = RATE(119)*Y(197)+RATE(138)*Y(82)+RATE(167)*Y(103)+RATE(167)&
    &*Y(103)+RATE(169)*Y(35)+RATE(177)*Y(133)+RATE(179)*Y(265)+RATE(256)*Y(29&
    &)/safeMantle+RATE(339)*D*Y(29)/safeMantle*Y(2)+RATE(422)*Y(29&
    &)/safeMantle+RATE(842)*Y(197)+RATE(873)*Y(82)+RATE(917)*Y(103)+RATE(917)&
    &*Y(103)+RATE(920)*Y(35)+RATE(922)*Y(36)+RATE(929)*Y(133)+RATE(931)*Y(265&
    &)+RATE(1104)*Y(29)+RATE(1187)*Y(31)+RATE(1236)*D*Y(16)*Y(36)+RATE(1250)&
    &*D*Y(16)*Y(82)+RATE(1257)*D*Y(16)*Y(103)+RATE(1263)*D*Y(16)*Y(35)&
    &+RATE(1265)*D*Y(16)*Y(133)+RATE(1266)*D*Y(16)*Y(265)+RATE(1328)*D*Y(18)&
    &*Y(265)+RATE(1422)*D*Y(199)*Y(333)+RATE(1424)*D*Y(290)*Y(333)+RATE(1437)&
    &*D*Y(20)*Y(28)+RATE(1460)*D*Y(20)*Y(36)+RATE(1479)*D*Y(20)*Y(103)&
    &+RATE(1483)*D*Y(20)*Y(133)+RATE(1564)*D*Y(24)*Y(36)+RATE(1583)*D*Y(24)&
    &*Y(133)+RATE(1757)*D*Y(82)*Y(133)+RATE(1769)*D*Y(83)*Y(333)+RATE(1848)*D&
    &*Y(2)*Y(35)+RATE(1851)*D*Y(2)*Y(133)+RATE(1852)*D*Y(2)*Y(265)+RATE(1973)&
    &*D*Y(6)*Y(36)+RATE(2365)*D*Y(13)*Y(82)+RATE(2390)*D*Y(13)*Y(90)&
    &+RATE(2391)*D*Y(13)*Y(90)+RATE(2399)*D*Y(13)*Y(92)+RATE(2406)*D*Y(13)&
    &*Y(103)+RATE(2413)*D*Y(13)*Y(133)+RATE(2415)*D*Y(13)*Y(265)+RATE(2553)*D&
    &*Y(28)*Y(67)+RATE(2554)*D*Y(28)*Y(75)+RATE(2555)*D*Y(28)*Y(24)+RATE(2556&
    &)*D*Y(28)*Y(41)+RATE(2557)*D*Y(28)*Y(82)+RATE(2558)*D*Y(28)*Y(99)&
    &+RATE(2559)*D*Y(28)*Y(131)+RATE(2560)*D*Y(28)*Y(62)+RATE(2561)*D*Y(28)&
    &*Y(183)+RATE(2562)*D*Y(28)*Y(90)+RATE(2563)*D*Y(28)*Y(116)+RATE(2564)*D&
    &*Y(28)*Y(69)+RATE(2565)*D*Y(28)*Y(43)+RATE(2566)*D*Y(28)*Y(54)+RATE(2567&
    &)*D*Y(28)*Y(35)+RATE(2568)*D*Y(28)*Y(133)+RATE(2569)*D*Y(28)*Y(161)&
    &+RATE(2570)*D*Y(28)*Y(300)+RATE(2571)*D*Y(28)*Y(56)+RATE(2576)*D*Y(28)&
    &*Y(41)+RATE(2588)*D*Y(28)*Y(291)+RATE(2596)*D*Y(28)*Y(300)+RATE(2598)*D&
    &*Y(28)*Y(333)+RATE(2608)*D*Y(104)*Y(333)+RATE(2608)*D*Y(104)*Y(333)&
    &+RATE(2615)*D*Y(120)*Y(333)+RATE(2624)*D*Y(35)*Y(68)+RATE(2628)*D*Y(35)&
    &*Y(100)+RATE(2629)*D*Y(35)*Y(132)+RATE(2630)*D*Y(35)*Y(63)+RATE(2635)*D&
    &*Y(35)*Y(44)+RATE(2636)*D*Y(35)*Y(55)+RATE(2643)*D*Y(35)*Y(82)+RATE(2648&
    &)*D*Y(35)*Y(35)+RATE(2655)*D*Y(35)*Y(46)+RATE(2656)*D*Y(35)*Y(56)&
    &+RATE(2659)*D*Y(35)*Y(163)+RATE(2667)*D*Y(36)*Y(333)+RATE(2668)*D*Y(36)&
    &*Y(67)+RATE(2671)*D*Y(36)*Y(75)+RATE(2672)*D*Y(36)*Y(82)+RATE(2673)*D&
    &*Y(36)*Y(232)+RATE(2676)*D*Y(36)*Y(99)+RATE(2677)*D*Y(36)*Y(131)&
    &+RATE(2679)*D*Y(36)*Y(62)+RATE(2683)*D*Y(36)*Y(90)+RATE(2684)*D*Y(36)&
    &*Y(116)+RATE(2685)*D*Y(36)*Y(92)+RATE(2686)*D*Y(36)*Y(103)+RATE(2687)*D&
    &*Y(36)*Y(43)+RATE(2688)*D*Y(36)*Y(54)+RATE(2689)*D*Y(36)*Y(35)+RATE(2692&
    &)*D*Y(36)*Y(161)+RATE(2693)*D*Y(36)*Y(46)+RATE(2694)*D*Y(36)*Y(56)&
    &+RATE(2695)*D*Y(36)*Y(163)+RATE(2733)*D*Y(44)*Y(333)+RATE(2821)*D*Y(133)&
    &*Y(163)+RATE(2822)*D*Y(134)*Y(333)+RATE(2823)*D*Y(266)*Y(333)+RATE(2842)&
    &*D*Y(46)*Y(104)+RATE(2868)*D*Y(46)*Y(82)+RATE(2889)*D*Y(46)*Y(103)&
    &+RATE(2894)*D*Y(46)*Y(133)+RATE(2896)*D*Y(46)*Y(265)+RATE(2945)*D*Y(48)&
    &*Y(90)+RATE(2948)*D*Y(48)*Y(103)+RATE(3079)*D*Y(105)*Y(133)
    YDOT(27) = PROD-LOSS
    LOSS = RATE(561)*D*Y(28)+RATE(1437)*D*Y(20)*Y(28)+RATE(1458)*D*Y(20)&
    &*Y(28)+RATE(1971)*D*Y(6)*Y(28)+RATE(2553)*D*Y(28)*Y(67)+RATE(2554)*D&
    &*Y(28)*Y(75)+RATE(2555)*D*Y(28)*Y(24)+RATE(2556)*D*Y(28)*Y(41)+RATE(2557&
    &)*D*Y(28)*Y(82)+RATE(2558)*D*Y(28)*Y(99)+RATE(2559)*D*Y(28)*Y(131)&
    &+RATE(2560)*D*Y(28)*Y(62)+RATE(2561)*D*Y(28)*Y(183)+RATE(2562)*D*Y(28)&
    &*Y(90)+RATE(2563)*D*Y(28)*Y(116)+RATE(2564)*D*Y(28)*Y(69)+RATE(2565)*D&
    &*Y(28)*Y(43)+RATE(2566)*D*Y(28)*Y(54)+RATE(2567)*D*Y(28)*Y(35)+RATE(2568&
    &)*D*Y(28)*Y(133)+RATE(2569)*D*Y(28)*Y(161)+RATE(2570)*D*Y(28)*Y(300)&
    &+RATE(2571)*D*Y(28)*Y(56)+RATE(2572)*D*Y(28)*Y(159)+RATE(2573)*D*Y(28)&
    &*Y(159)+RATE(2574)*D*Y(28)*Y(159)+RATE(2575)*D*Y(28)*Y(159)+RATE(2576)*D&
    &*Y(28)*Y(41)+RATE(2577)*D*Y(28)*Y(41)+RATE(2578)*D*Y(28)*Y(41)+RATE(2579&
    &)*D*Y(28)*Y(232)+RATE(2580)*D*Y(28)*Y(99)+RATE(2581)*D*Y(28)*Y(131)&
    &+RATE(2582)*D*Y(28)*Y(131)+RATE(2583)*D*Y(28)*Y(183)+RATE(2584)*D*Y(28)&
    &*Y(183)+RATE(2585)*D*Y(28)*Y(183)+RATE(2586)*D*Y(28)*Y(116)+RATE(2587)*D&
    &*Y(28)*Y(291)+RATE(2588)*D*Y(28)*Y(291)+RATE(2589)*D*Y(28)*Y(54)&
    &+RATE(2590)*D*Y(28)*Y(54)+RATE(2591)*D*Y(28)*Y(35)+RATE(2592)*D*Y(28)&
    &*Y(133)+RATE(2593)*D*Y(28)*Y(161)+RATE(2594)*D*Y(28)*Y(161)+RATE(2595)*D&
    &*Y(28)*Y(300)+RATE(2596)*D*Y(28)*Y(300)+RATE(2597)*D*Y(28)*Y(27)&
    &+RATE(2598)*D*Y(28)*Y(333)
    PROD = RATE(107)*Y(27)+RATE(166)*Y(27)+RATE(2364)*D*Y(13)*Y(82)&
    &+RATE(2389)*D*Y(13)*Y(90)+RATE(2406)*D*Y(13)*Y(103)+RATE(2408)*D*Y(13)&
    &*Y(43)+RATE(2412)*D*Y(13)*Y(35)+RATE(2414)*D*Y(13)*Y(133)+RATE(2416)*D&
    &*Y(13)*Y(265)+RATE(2489)*D*Y(27)*Y(104)
    YDOT(28) = PROD-LOSS
    LOSS = RATE(256)*Y(29)/safeMantle+RATE(339)*D*Y(29)/safeMantle*Y(2)&
    &+RATE(422)*Y(29)/safeMantle+RATE(635)*Y(3)*Y(29)+RATE(655)*Y(29)*Y(29)&
    &+RATE(655)*Y(29)*Y(29)+RATE(739)*Y(3)*Y(29)+RATE(759)*Y(29)*Y(29)&
    &+RATE(759)*Y(29)*Y(29)+RATE(1021)*Y(29)*totalSwap/safeMantle+RATE(1104)&
    &*Y(29)
    PROD = RATE(54)*Y(31)*bulkLayersReciprocal+RATE(560)*D*Y(27)&
    &+RATE(561)*D*Y(28)
    YDOT(29) = PROD-LOSS
    LOSS = RATE(14)*Y(30)*bulkLayersReciprocal+RATE(683)*Y(5)*Y(30)&
    &*bulkLayersReciprocal+RATE(787)*Y(5)*Y(30)*bulkLayersReciprocal&
    &+RATE(1147)*Y(30)
    PROD = RATE(682)*Y(5)*Y(23)*bulkLayersReciprocal+RATE(981)*Y(26)&
    &*totalSwap/safeMantle
    YDOT(30) = PROD-LOSS
    LOSS = RATE(54)*Y(31)*bulkLayersReciprocal+RATE(687)*Y(5)*Y(31)&
    &*bulkLayersReciprocal+RATE(707)*Y(31)*Y(31)*bulkLayersReciprocal&
    &+RATE(707)*Y(31)*Y(31)*bulkLayersReciprocal+RATE(791)*Y(5)*Y(31)&
    &*bulkLayersReciprocal+RATE(811)*Y(31)*Y(31)*bulkLayersReciprocal&
    &+RATE(811)*Y(31)*Y(31)*bulkLayersReciprocal+RATE(1187)*Y(31)
    PROD = RATE(1021)*Y(29)*totalSwap/safeMantle
    YDOT(31) = PROD-LOSS
    LOSS = RATE(128)*Y(32)+RATE(129)*Y(32)+RATE(130)*Y(32)+RATE(482)*D&
    &*Y(32)+RATE(855)*Y(32)+RATE(856)*Y(32)+RATE(857)*Y(32)+RATE(1248)*D*Y(16&
    &)*Y(32)+RATE(1305)*D*Y(18)*Y(32)+RATE(1306)*D*Y(18)*Y(32)+RATE(1617)*D&
    &*Y(32)*Y(165)+RATE(1618)*D*Y(32)*Y(89)+RATE(1619)*D*Y(32)*Y(32)&
    &+RATE(1619)*D*Y(32)*Y(32)+RATE(1620)*D*Y(32)*Y(32)+RATE(1620)*D*Y(32)&
    &*Y(32)+RATE(1621)*D*Y(32)*Y(32)+RATE(1621)*D*Y(32)*Y(32)+RATE(1622)*D&
    &*Y(32)*Y(82)+RATE(1623)*D*Y(32)*Y(131)+RATE(1624)*D*Y(32)*Y(62)&
    &+RATE(1625)*D*Y(32)*Y(183)+RATE(1626)*D*Y(32)*Y(116)+RATE(1627)*D*Y(32)&
    &*Y(144)+RATE(1628)*D*Y(32)*Y(43)+RATE(1629)*D*Y(32)*Y(54)+RATE(1630)*D&
    &*Y(32)*Y(264)+RATE(1631)*D*Y(32)*Y(133)+RATE(1632)*D*Y(32)*Y(161)&
    &+RATE(1633)*D*Y(32)*Y(161)+RATE(1634)*D*Y(32)*Y(161)+RATE(1635)*D*Y(32)&
    &*Y(175)+RATE(1636)*D*Y(32)*Y(46)+RATE(1637)*D*Y(32)*Y(46)+RATE(1638)*D&
    &*Y(32)*Y(56)+RATE(1639)*D*Y(32)*Y(56)+RATE(1640)*D*Y(32)*Y(56)+RATE(1641&
    &)*D*Y(32)*Y(163)+RATE(1642)*D*Y(32)*Y(82)+RATE(1827)*D*Y(2)*Y(32)&
    &+RATE(1877)*D*Y(4)*Y(32)+RATE(1987)*D*Y(6)*Y(32)+RATE(2154)*D*Y(10)*Y(32&
    &)+RATE(2352)*D*Y(13)*Y(32)+RATE(2524)*D*Y(27)*Y(32)+RATE(2525)*D*Y(27)&
    &*Y(32)+RATE(2526)*D*Y(27)*Y(32)
    PROD = RATE(132)*Y(237)+RATE(133)*Y(214)+RATE(135)*Y(159)+RATE(219)&
    &*Y(34)/safeMantle+RATE(302)*D*Y(34)/safeMantle*Y(2)+RATE(385)*Y(34&
    &)/safeMantle+RATE(735)*Y(3)*Y(26)+RATE(787)*Y(5)*Y(30)&
    &*bulkLayersReciprocal+RATE(861)*Y(237)+RATE(862)*Y(214)+RATE(865)*Y(159)&
    &+RATE(867)*Y(41)+RATE(1067)*Y(34)+RATE(1150)*Y(38)+RATE(1421)*D*Y(109)&
    &*Y(165)+RATE(1575)*D*Y(24)*Y(24)+RATE(1576)*D*Y(24)*Y(41)+RATE(1576)*D&
    &*Y(24)*Y(41)+RATE(1578)*D*Y(24)*Y(131)+RATE(1579)*D*Y(24)*Y(116)&
    &+RATE(1580)*D*Y(24)*Y(144)+RATE(1597)*D*Y(24)*Y(56)+RATE(1605)*D*Y(25)&
    &*Y(131)+RATE(1643)*D*Y(33)*Y(116)+RATE(1644)*D*Y(33)*Y(69)+RATE(1645)*D&
    &*Y(33)*Y(133)+RATE(1668)*D*Y(33)*Y(333)+RATE(1671)*D*Y(220)*Y(333)&
    &+RATE(1674)*D*Y(172)*Y(333)+RATE(1675)*D*Y(172)*Y(333)+RATE(1680)*D*Y(41&
    &)*Y(68)+RATE(1682)*D*Y(41)*Y(76)+RATE(1684)*D*Y(41)*Y(100)+RATE(1685)*D&
    &*Y(41)*Y(234)+RATE(1686)*D*Y(41)*Y(132)+RATE(1687)*D*Y(41)*Y(63)&
    &+RATE(1688)*D*Y(41)*Y(91)+RATE(1695)*D*Y(41)*Y(55)+RATE(1700)*D*Y(41)&
    &*Y(82)+RATE(1701)*D*Y(41)*Y(161)+RATE(1702)*D*Y(41)*Y(56)+RATE(1703)*D&
    &*Y(41)*Y(163)+RATE(1711)*D*Y(42)*Y(333)+RATE(1712)*D*Y(42)*Y(159)&
    &+RATE(1713)*D*Y(42)*Y(41)+RATE(1714)*D*Y(42)*Y(232)+RATE(1715)*D*Y(42)&
    &*Y(99)+RATE(1716)*D*Y(42)*Y(131)+RATE(1717)*D*Y(42)*Y(62)+RATE(1718)*D&
    &*Y(42)*Y(183)+RATE(1719)*D*Y(42)*Y(54)+RATE(1720)*D*Y(42)*Y(300)&
    &+RATE(1722)*D*Y(53)*Y(333)+RATE(1723)*D*Y(53)*Y(333)+RATE(1826)*D*Y(2)&
    &*Y(218)+RATE(1828)*D*Y(2)*Y(41)+RATE(1986)*D*Y(6)*Y(24)+RATE(2002)*D*Y(6&
    &)*Y(20)+RATE(2355)*D*Y(13)*Y(214)+RATE(2357)*D*Y(13)*Y(159)+RATE(2362)*D&
    &*Y(13)*Y(41)+RATE(2512)*D*Y(27)*Y(109)+RATE(2514)*D*Y(27)*Y(123)&
    &+RATE(2574)*D*Y(28)*Y(159)+RATE(2642)*D*Y(35)*Y(41)+RATE(2722)*D*Y(43)&
    &*Y(41)+RATE(2860)*D*Y(46)*Y(109)+RATE(2862)*D*Y(46)*Y(123)+RATE(2867)*D&
    &*Y(46)*Y(41)+RATE(2995)*D*Y(56)*Y(80)+RATE(3081)*D*Y(106)*Y(159)
    YDOT(32) = PROD-LOSS
    LOSS = RATE(483)*D*Y(33)+RATE(858)*Y(33)+RATE(859)*Y(33)+RATE(1225)&
    &*D*Y(16)*Y(33)+RATE(1445)*D*Y(20)*Y(33)+RATE(1643)*D*Y(33)*Y(116)&
    &+RATE(1644)*D*Y(33)*Y(69)+RATE(1645)*D*Y(33)*Y(133)+RATE(1646)*D*Y(33)&
    &*Y(333)+RATE(1647)*D*Y(33)*Y(333)+RATE(1648)*D*Y(33)*Y(333)+RATE(1649)*D&
    &*Y(33)*Y(109)+RATE(1650)*D*Y(33)*Y(237)+RATE(1651)*D*Y(33)*Y(214)&
    &+RATE(1652)*D*Y(33)*Y(159)+RATE(1653)*D*Y(33)*Y(131)+RATE(1654)*D*Y(33)&
    &*Y(183)+RATE(1655)*D*Y(33)*Y(116)+RATE(1656)*D*Y(33)*Y(173)+RATE(1657)*D&
    &*Y(33)*Y(54)+RATE(1658)*D*Y(33)*Y(161)+RATE(1659)*D*Y(33)*Y(46)&
    &+RATE(1660)*D*Y(33)*Y(46)+RATE(1661)*D*Y(33)*Y(300)+RATE(1662)*D*Y(33)&
    &*Y(56)+RATE(1663)*D*Y(33)*Y(163)+RATE(1664)*D*Y(33)*Y(277)+RATE(1665)*D&
    &*Y(33)*Y(166)+RATE(1666)*D*Y(33)*Y(62)+RATE(1667)*D*Y(33)*Y(90)&
    &+RATE(1668)*D*Y(33)*Y(333)+RATE(1812)*D*Y(2)*Y(33)+RATE(2001)*D*Y(6)&
    &*Y(33)+RATE(2626)*D*Y(35)*Y(33)
    PROD = RATE(129)*Y(32)+RATE(856)*Y(32)+RATE(871)*Y(42)+RATE(1309)*D&
    &*Y(18)*Y(159)+RATE(1510)*D*Y(21)*Y(159)+RATE(1516)*D*Y(21)*Y(131)&
    &+RATE(1552)*D*Y(24)*Y(76)+RATE(1553)*D*Y(24)*Y(53)+RATE(1555)*D*Y(24)&
    &*Y(132)+RATE(1556)*D*Y(24)*Y(63)+RATE(1557)*D*Y(24)*Y(66)+RATE(1558)*D&
    &*Y(24)*Y(91)+RATE(1559)*D*Y(24)*Y(102)+RATE(1560)*D*Y(24)*Y(102)&
    &+RATE(1561)*D*Y(24)*Y(117)+RATE(1562)*D*Y(24)*Y(145)+RATE(1563)*D*Y(24)&
    &*Y(120)+RATE(1564)*D*Y(24)*Y(36)+RATE(1565)*D*Y(24)*Y(44)+RATE(1566)*D&
    &*Y(24)*Y(55)+RATE(1568)*D*Y(24)*Y(176)+RATE(1569)*D*Y(24)*Y(57)&
    &+RATE(1610)*D*Y(25)*Y(116)+RATE(1693)*D*Y(41)*Y(104)+RATE(1813)*D*Y(2)&
    &*Y(42)+RATE(1877)*D*Y(4)*Y(32)+RATE(1919)*D*Y(4)*Y(214)+RATE(1920)*D*Y(4&
    &)*Y(159)+RATE(1923)*D*Y(4)*Y(41)+RATE(1957)*D*Y(6)*Y(25)+RATE(2032)*D&
    &*Y(8)*Y(24)+RATE(2033)*D*Y(8)*Y(41)+RATE(2153)*D*Y(10)*Y(24)+RATE(2155)&
    &*D*Y(10)*Y(237)+RATE(2159)*D*Y(10)*Y(159)+RATE(2354)*D*Y(13)*Y(237)&
    &+RATE(2356)*D*Y(13)*Y(214)+RATE(2358)*D*Y(13)*Y(159)+RATE(2361)*D*Y(13)&
    &*Y(41)+RATE(2575)*D*Y(28)*Y(159)+RATE(2576)*D*Y(28)*Y(41)+RATE(2830)*D&
    &*Y(46)*Y(42)+RATE(2939)*D*Y(48)*Y(41)
    YDOT(33) = PROD-LOSS
    LOSS = RATE(84)*Y(34)+RATE(219)*Y(34)/safeMantle+RATE(302)*D*Y(34&
    &)/safeMantle*Y(2)+RATE(385)*Y(34)/safeMantle+RATE(622)*Y(34)*Y(226)&
    &+RATE(623)*Y(34)*Y(118)+RATE(632)*Y(3)*Y(34)+RATE(726)*Y(34)*Y(226)&
    &+RATE(727)*Y(34)*Y(118)+RATE(736)*Y(3)*Y(34)+RATE(984)*Y(34)&
    &*totalSwap/safeMantle+RATE(1067)*Y(34)
    PROD = RATE(17)*Y(38)*bulkLayersReciprocal+RATE(87)*Y(157)+RATE(482)&
    &*D*Y(32)+RATE(483)*D*Y(33)+RATE(631)*Y(3)*Y(26)
    YDOT(34) = PROD-LOSS
    LOSS = RATE(169)*Y(35)+RATE(170)*Y(35)+RATE(566)*D*Y(35)+RATE(920)&
    &*Y(35)+RATE(921)*Y(35)+RATE(1262)*D*Y(16)*Y(35)+RATE(1263)*D*Y(16)*Y(35)&
    &+RATE(1327)*D*Y(18)*Y(35)+RATE(1532)*D*Y(21)*Y(35)+RATE(1848)*D*Y(2)&
    &*Y(35)+RATE(1894)*D*Y(4)*Y(35)+RATE(1993)*D*Y(6)*Y(35)+RATE(2023)*D*Y(8)&
    &*Y(35)+RATE(2048)*D*Y(8)*Y(35)+RATE(2187)*D*Y(10)*Y(35)+RATE(2412)*D&
    &*Y(13)*Y(35)+RATE(2539)*D*Y(27)*Y(35)+RATE(2567)*D*Y(28)*Y(35)+RATE(2591&
    &)*D*Y(28)*Y(35)+RATE(2620)*D*Y(35)*Y(83)+RATE(2621)*D*Y(35)*Y(100)&
    &+RATE(2622)*D*Y(35)*Y(104)+RATE(2623)*D*Y(35)*Y(48)+RATE(2624)*D*Y(35)&
    &*Y(68)+RATE(2625)*D*Y(35)*Y(68)+RATE(2626)*D*Y(35)*Y(33)+RATE(2627)*D&
    &*Y(35)*Y(53)+RATE(2628)*D*Y(35)*Y(100)+RATE(2629)*D*Y(35)*Y(132)&
    &+RATE(2630)*D*Y(35)*Y(63)+RATE(2631)*D*Y(35)*Y(91)+RATE(2632)*D*Y(35)&
    &*Y(117)+RATE(2633)*D*Y(35)*Y(145)+RATE(2634)*D*Y(35)*Y(120)+RATE(2635)*D&
    &*Y(35)*Y(44)+RATE(2636)*D*Y(35)*Y(55)+RATE(2637)*D*Y(35)*Y(48)+RATE(2638&
    &)*D*Y(35)*Y(162)+RATE(2639)*D*Y(35)*Y(176)+RATE(2640)*D*Y(35)*Y(57)&
    &+RATE(2641)*D*Y(35)*Y(165)+RATE(2642)*D*Y(35)*Y(41)+RATE(2643)*D*Y(35)&
    &*Y(82)+RATE(2644)*D*Y(35)*Y(62)+RATE(2645)*D*Y(35)*Y(54)+RATE(2646)*D&
    &*Y(35)*Y(35)+RATE(2646)*D*Y(35)*Y(35)+RATE(2647)*D*Y(35)*Y(35)+RATE(2647&
    &)*D*Y(35)*Y(35)+RATE(2648)*D*Y(35)*Y(35)+RATE(2648)*D*Y(35)*Y(35)&
    &+RATE(2649)*D*Y(35)*Y(264)+RATE(2650)*D*Y(35)*Y(133)+RATE(2651)*D*Y(35)&
    &*Y(133)+RATE(2652)*D*Y(35)*Y(161)+RATE(2653)*D*Y(35)*Y(161)+RATE(2654)*D&
    &*Y(35)*Y(46)+RATE(2655)*D*Y(35)*Y(46)+RATE(2656)*D*Y(35)*Y(56)+RATE(2657&
    &)*D*Y(35)*Y(56)+RATE(2658)*D*Y(35)*Y(56)+RATE(2659)*D*Y(35)*Y(163)&
    &+RATE(2660)*D*Y(35)*Y(163)+RATE(2689)*D*Y(36)*Y(35)
    PROD = RATE(161)*Y(227)+RATE(172)*Y(43)+RATE(175)*Y(54)+RATE(259)&
    &*Y(37)/safeMantle+RATE(342)*D*Y(37)/safeMantle*Y(2)+RATE(425)*Y(37&
    &)/safeMantle+RATE(739)*Y(3)*Y(29)+RATE(791)*Y(5)*Y(31)&
    &*bulkLayersReciprocal+RATE(909)*Y(227)+RATE(924)*Y(43)+RATE(927)*Y(54)&
    &+RATE(1107)*Y(37)+RATE(1190)*Y(39)+RATE(1261)*D*Y(16)*Y(43)+RATE(1461)*D&
    &*Y(20)*Y(44)+RATE(1481)*D*Y(20)*Y(27)+RATE(1565)*D*Y(24)*Y(44)+RATE(1581&
    &)*D*Y(24)*Y(103)+RATE(1628)*D*Y(32)*Y(43)+RATE(1783)*D*Y(99)*Y(144)&
    &+RATE(1843)*D*Y(2)*Y(144)+RATE(1846)*D*Y(2)*Y(43)+RATE(1850)*D*Y(2)&
    &*Y(133)+RATE(1853)*D*Y(2)*Y(265)+RATE(1859)*D*Y(2)*Y(221)+RATE(1991)*D&
    &*Y(6)*Y(27)+RATE(2086)*D*Y(62)*Y(83)+RATE(2511)*D*Y(27)*Y(89)+RATE(2513)&
    &*D*Y(27)*Y(123)+RATE(2523)*D*Y(27)*Y(24)+RATE(2530)*D*Y(27)*Y(101)&
    &+RATE(2531)*D*Y(27)*Y(116)+RATE(2535)*D*Y(27)*Y(144)+RATE(2537)*D*Y(27)&
    &*Y(173)+RATE(2546)*D*Y(27)*Y(175)+RATE(2548)*D*Y(27)*Y(56)+RATE(2572)*D&
    &*Y(28)*Y(159)+RATE(2573)*D*Y(28)*Y(159)+RATE(2581)*D*Y(28)*Y(131)&
    &+RATE(2583)*D*Y(28)*Y(183)+RATE(2585)*D*Y(28)*Y(183)+RATE(2590)*D*Y(28)&
    &*Y(54)+RATE(2615)*D*Y(120)*Y(333)+RATE(2661)*D*Y(36)*Y(131)+RATE(2662)*D&
    &*Y(36)*Y(62)+RATE(2663)*D*Y(36)*Y(54)+RATE(2664)*D*Y(36)*Y(133)&
    &+RATE(2665)*D*Y(36)*Y(161)+RATE(2666)*D*Y(36)*Y(163)+RATE(2708)*D*Y(43)&
    &*Y(100)+RATE(2719)*D*Y(43)*Y(55)+RATE(2725)*D*Y(43)*Y(56)+RATE(2734)*D&
    &*Y(44)*Y(333)+RATE(2735)*D*Y(44)*Y(67)+RATE(2736)*D*Y(44)*Y(75)&
    &+RATE(2737)*D*Y(44)*Y(131)+RATE(2739)*D*Y(44)*Y(62)+RATE(2742)*D*Y(44)&
    &*Y(183)+RATE(2746)*D*Y(44)*Y(90)+RATE(2747)*D*Y(44)*Y(116)+RATE(2748)*D&
    &*Y(44)*Y(92)+RATE(2749)*D*Y(44)*Y(43)+RATE(2750)*D*Y(44)*Y(54)+RATE(2754&
    &)*D*Y(44)*Y(163)+RATE(2798)*D*Y(55)*Y(333)+RATE(2799)*D*Y(55)*Y(67)&
    &+RATE(2878)*D*Y(46)*Y(90)+RATE(2886)*D*Y(46)*Y(144)+RATE(2891)*D*Y(46)&
    &*Y(43)
    YDOT(35) = PROD-LOSS
    LOSS = RATE(567)*D*Y(36)+RATE(922)*Y(36)+RATE(1236)*D*Y(16)*Y(36)&
    &+RATE(1460)*D*Y(20)*Y(36)+RATE(1564)*D*Y(24)*Y(36)+RATE(1973)*D*Y(6)&
    &*Y(36)+RATE(1974)*D*Y(6)*Y(36)+RATE(2502)*D*Y(27)*Y(36)+RATE(2661)*D&
    &*Y(36)*Y(131)+RATE(2662)*D*Y(36)*Y(62)+RATE(2663)*D*Y(36)*Y(54)&
    &+RATE(2664)*D*Y(36)*Y(133)+RATE(2665)*D*Y(36)*Y(161)+RATE(2666)*D*Y(36)&
    &*Y(163)+RATE(2667)*D*Y(36)*Y(333)+RATE(2668)*D*Y(36)*Y(67)+RATE(2669)*D&
    &*Y(36)*Y(67)+RATE(2670)*D*Y(36)*Y(67)+RATE(2671)*D*Y(36)*Y(75)+RATE(2672&
    &)*D*Y(36)*Y(82)+RATE(2673)*D*Y(36)*Y(232)+RATE(2674)*D*Y(36)*Y(232)&
    &+RATE(2675)*D*Y(36)*Y(232)+RATE(2676)*D*Y(36)*Y(99)+RATE(2677)*D*Y(36)&
    &*Y(131)+RATE(2678)*D*Y(36)*Y(131)+RATE(2679)*D*Y(36)*Y(62)+RATE(2680)*D&
    &*Y(36)*Y(62)+RATE(2681)*D*Y(36)*Y(62)+RATE(2682)*D*Y(36)*Y(62)+RATE(2683&
    &)*D*Y(36)*Y(90)+RATE(2684)*D*Y(36)*Y(116)+RATE(2685)*D*Y(36)*Y(92)&
    &+RATE(2686)*D*Y(36)*Y(103)+RATE(2687)*D*Y(36)*Y(43)+RATE(2688)*D*Y(36)&
    &*Y(54)+RATE(2689)*D*Y(36)*Y(35)+RATE(2690)*D*Y(36)*Y(133)+RATE(2691)*D&
    &*Y(36)*Y(161)+RATE(2692)*D*Y(36)*Y(161)+RATE(2693)*D*Y(36)*Y(46)&
    &+RATE(2694)*D*Y(36)*Y(56)+RATE(2695)*D*Y(36)*Y(163)+RATE(2696)*D*Y(36)&
    &*Y(163)
    PROD = RATE(170)*Y(35)+RATE(921)*Y(35)+RATE(1894)*D*Y(4)*Y(35)&
    &+RATE(1971)*D*Y(6)*Y(28)+RATE(2023)*D*Y(8)*Y(35)+RATE(2047)*D*Y(8)*Y(27)&
    &+RATE(2400)*D*Y(13)*Y(92)+RATE(2409)*D*Y(13)*Y(43)+RATE(2410)*D*Y(13)&
    &*Y(54)+RATE(2567)*D*Y(28)*Y(35)+RATE(2584)*D*Y(28)*Y(183)+RATE(2586)*D&
    &*Y(28)*Y(116)+RATE(2620)*D*Y(35)*Y(83)+RATE(2621)*D*Y(35)*Y(100)&
    &+RATE(2622)*D*Y(35)*Y(104)+RATE(2623)*D*Y(35)*Y(48)
    YDOT(36) = PROD-LOSS
    LOSS = RATE(259)*Y(37)/safeMantle+RATE(342)*D*Y(37)/safeMantle*Y(2)&
    &+RATE(425)*Y(37)/safeMantle+RATE(636)*Y(3)*Y(37)+RATE(656)*Y(37)*Y(97)&
    &+RATE(657)*Y(37)*Y(226)+RATE(740)*Y(3)*Y(37)+RATE(760)*Y(37)*Y(97)&
    &+RATE(761)*Y(37)*Y(226)+RATE(1024)*Y(37)*totalSwap/safeMantle+RATE(1107)&
    &*Y(37)
    PROD = RATE(57)*Y(39)*bulkLayersReciprocal+RATE(97)*Y(52)+RATE(109)&
    &*Y(226)+RATE(566)*D*Y(35)+RATE(567)*D*Y(36)+RATE(635)*Y(3)*Y(29)&
    &+RATE(829)*Y(226)
    YDOT(37) = PROD-LOSS
    LOSS = RATE(17)*Y(38)*bulkLayersReciprocal+RATE(674)*Y(38)*Y(228)&
    &*bulkLayersReciprocal+RATE(675)*Y(38)*Y(126)*bulkLayersReciprocal&
    &+RATE(684)*Y(5)*Y(38)*bulkLayersReciprocal+RATE(778)*Y(38)*Y(228)&
    &*bulkLayersReciprocal+RATE(779)*Y(38)*Y(126)*bulkLayersReciprocal&
    &+RATE(788)*Y(5)*Y(38)*bulkLayersReciprocal+RATE(1150)*Y(38)
    PROD = RATE(683)*Y(5)*Y(30)*bulkLayersReciprocal+RATE(984)*Y(34)&
    &*totalSwap/safeMantle
    YDOT(38) = PROD-LOSS
    LOSS = RATE(57)*Y(39)*bulkLayersReciprocal+RATE(688)*Y(5)*Y(39)&
    &*bulkLayersReciprocal+RATE(708)*Y(39)*Y(111)*bulkLayersReciprocal&
    &+RATE(709)*Y(39)*Y(228)*bulkLayersReciprocal+RATE(792)*Y(5)*Y(39)&
    &*bulkLayersReciprocal+RATE(812)*Y(39)*Y(111)*bulkLayersReciprocal&
    &+RATE(813)*Y(39)*Y(228)*bulkLayersReciprocal+RATE(1190)*Y(39)
    PROD = RATE(687)*Y(5)*Y(31)*bulkLayersReciprocal+RATE(1024)*Y(37)&
    &*totalSwap/safeMantle
    YDOT(39) = PROD-LOSS
    LOSS = RATE(88)*Y(40)+RATE(225)*Y(40)/safeMantle+RATE(308)*D*Y(40&
    &)/safeMantle*Y(2)+RATE(391)*Y(40)/safeMantle+RATE(990)*Y(40)&
    &*totalSwap/safeMantle+RATE(1073)*Y(40)
    PROD = RATE(23)*Y(49)*bulkLayersReciprocal+RATE(491)*D*Y(41)&
    &+RATE(492)*D*Y(42)+RATE(493)*D*Y(53)+RATE(622)*Y(34)*Y(226)+RATE(632)&
    &*Y(3)*Y(34)
    YDOT(40) = PROD-LOSS
    LOSS = RATE(136)*Y(41)+RATE(491)*D*Y(41)+RATE(866)*Y(41)+RATE(867)&
    &*Y(41)+RATE(868)*Y(41)+RATE(869)*Y(41)+RATE(1310)*D*Y(18)*Y(41)&
    &+RATE(1474)*D*Y(20)*Y(41)+RATE(1512)*D*Y(21)*Y(41)+RATE(1576)*D*Y(24)&
    &*Y(41)+RATE(1679)*D*Y(41)*Y(100)+RATE(1680)*D*Y(41)*Y(68)+RATE(1681)*D&
    &*Y(41)*Y(68)+RATE(1682)*D*Y(41)*Y(76)+RATE(1683)*D*Y(41)*Y(81)+RATE(1684&
    &)*D*Y(41)*Y(100)+RATE(1685)*D*Y(41)*Y(234)+RATE(1686)*D*Y(41)*Y(132)&
    &+RATE(1687)*D*Y(41)*Y(63)+RATE(1688)*D*Y(41)*Y(91)+RATE(1689)*D*Y(41)&
    &*Y(243)+RATE(1690)*D*Y(41)*Y(145)+RATE(1691)*D*Y(41)*Y(174)+RATE(1692)*D&
    &*Y(41)*Y(104)+RATE(1693)*D*Y(41)*Y(104)+RATE(1694)*D*Y(41)*Y(120)&
    &+RATE(1695)*D*Y(41)*Y(55)+RATE(1696)*D*Y(41)*Y(57)+RATE(1697)*D*Y(41)&
    &*Y(57)+RATE(1698)*D*Y(41)*Y(165)+RATE(1699)*D*Y(41)*Y(165)+RATE(1700)*D&
    &*Y(41)*Y(82)+RATE(1701)*D*Y(41)*Y(161)+RATE(1702)*D*Y(41)*Y(56)&
    &+RATE(1703)*D*Y(41)*Y(163)+RATE(1713)*D*Y(42)*Y(41)+RATE(1828)*D*Y(2)&
    &*Y(41)+RATE(1878)*D*Y(4)*Y(41)+RATE(1923)*D*Y(4)*Y(41)+RATE(2012)*D*Y(8)&
    &*Y(41)+RATE(2033)*D*Y(8)*Y(41)+RATE(2034)*D*Y(8)*Y(41)+RATE(2161)*D*Y(10&
    &)*Y(41)+RATE(2323)*D*Y(13)*Y(41)+RATE(2359)*D*Y(13)*Y(41)+RATE(2360)*D&
    &*Y(13)*Y(41)+RATE(2361)*D*Y(13)*Y(41)+RATE(2362)*D*Y(13)*Y(41)+RATE(2556&
    &)*D*Y(28)*Y(41)+RATE(2576)*D*Y(28)*Y(41)+RATE(2577)*D*Y(28)*Y(41)&
    &+RATE(2578)*D*Y(28)*Y(41)+RATE(2642)*D*Y(35)*Y(41)+RATE(2722)*D*Y(43)&
    &*Y(41)+RATE(2867)*D*Y(46)*Y(41)+RATE(2922)*D*Y(48)*Y(41)+RATE(2939)*D&
    &*Y(48)*Y(41)
    PROD = RATE(131)*Y(237)+RATE(225)*Y(40)/safeMantle+RATE(308)*D*Y(40&
    &)/safeMantle*Y(2)+RATE(391)*Y(40)/safeMantle+RATE(726)*Y(34)*Y(226)&
    &+RATE(736)*Y(3)*Y(34)+RATE(778)*Y(38)*Y(228)*bulkLayersReciprocal&
    &+RATE(788)*Y(5)*Y(38)*bulkLayersReciprocal+RATE(860)*Y(237)+RATE(1073)&
    &*Y(40)+RATE(1156)*Y(49)+RATE(1226)*D*Y(16)*Y(53)+RATE(1446)*D*Y(20)*Y(53&
    &)+RATE(1553)*D*Y(24)*Y(53)+RATE(1618)*D*Y(32)*Y(89)+RATE(1621)*D*Y(32)&
    &*Y(32)+RATE(1623)*D*Y(32)*Y(131)+RATE(1624)*D*Y(32)*Y(62)+RATE(1625)*D&
    &*Y(32)*Y(183)+RATE(1626)*D*Y(32)*Y(116)+RATE(1627)*D*Y(32)*Y(144)&
    &+RATE(1628)*D*Y(32)*Y(43)+RATE(1629)*D*Y(32)*Y(54)+RATE(1635)*D*Y(32)&
    &*Y(175)+RATE(1638)*D*Y(32)*Y(56)+RATE(1652)*D*Y(33)*Y(159)+RATE(1653)*D&
    &*Y(33)*Y(131)+RATE(1665)*D*Y(33)*Y(166)+RATE(1704)*D*Y(42)*Y(80)&
    &+RATE(1705)*D*Y(42)*Y(131)+RATE(1706)*D*Y(42)*Y(183)+RATE(1707)*D*Y(42)&
    &*Y(54)+RATE(1708)*D*Y(42)*Y(161)+RATE(1709)*D*Y(42)*Y(300)+RATE(1724)*D&
    &*Y(53)*Y(333)+RATE(1726)*D*Y(53)*Y(67)+RATE(1727)*D*Y(53)*Y(75)&
    &+RATE(1728)*D*Y(53)*Y(232)+RATE(1729)*D*Y(53)*Y(99)+RATE(1730)*D*Y(53)&
    &*Y(131)+RATE(1731)*D*Y(53)*Y(62)+RATE(1732)*D*Y(53)*Y(183)+RATE(1733)*D&
    &*Y(53)*Y(90)+RATE(1734)*D*Y(53)*Y(116)+RATE(1735)*D*Y(53)*Y(193)&
    &+RATE(1736)*D*Y(53)*Y(92)+RATE(1737)*D*Y(53)*Y(69)+RATE(1738)*D*Y(53)&
    &*Y(163)+RATE(1739)*D*Y(53)*Y(166)+RATE(1987)*D*Y(6)*Y(32)+RATE(2627)*D&
    &*Y(35)*Y(53)+RATE(2707)*D*Y(43)*Y(53)+RATE(2770)*D*Y(54)*Y(53)+RATE(2981&
    &)*D*Y(56)*Y(53)
    YDOT(41) = PROD-LOSS
    LOSS = RATE(492)*D*Y(42)+RATE(870)*Y(42)+RATE(871)*Y(42)+RATE(1704)&
    &*D*Y(42)*Y(80)+RATE(1705)*D*Y(42)*Y(131)+RATE(1706)*D*Y(42)*Y(183)&
    &+RATE(1707)*D*Y(42)*Y(54)+RATE(1708)*D*Y(42)*Y(161)+RATE(1709)*D*Y(42)&
    &*Y(300)+RATE(1710)*D*Y(42)*Y(333)+RATE(1711)*D*Y(42)*Y(333)+RATE(1712)*D&
    &*Y(42)*Y(159)+RATE(1713)*D*Y(42)*Y(41)+RATE(1714)*D*Y(42)*Y(232)&
    &+RATE(1715)*D*Y(42)*Y(99)+RATE(1716)*D*Y(42)*Y(131)+RATE(1717)*D*Y(42)&
    &*Y(62)+RATE(1718)*D*Y(42)*Y(183)+RATE(1719)*D*Y(42)*Y(54)+RATE(1720)*D&
    &*Y(42)*Y(300)+RATE(1813)*D*Y(2)*Y(42)+RATE(1958)*D*Y(6)*Y(42)+RATE(2830)&
    &*D*Y(46)*Y(42)
    PROD = RATE(868)*Y(41)+RATE(1655)*D*Y(33)*Y(116)+RATE(1679)*D*Y(41)&
    &*Y(100)+RATE(1814)*D*Y(2)*Y(53)+RATE(1878)*D*Y(4)*Y(41)+RATE(2012)*D*Y(8&
    &)*Y(41)+RATE(2154)*D*Y(10)*Y(32)+RATE(2323)*D*Y(13)*Y(41)+RATE(2556)*D&
    &*Y(28)*Y(41)+RATE(2922)*D*Y(48)*Y(41)
    YDOT(42) = PROD-LOSS
    LOSS = RATE(171)*Y(43)+RATE(172)*Y(43)+RATE(568)*D*Y(43)+RATE(923)&
    &*Y(43)+RATE(924)*Y(43)+RATE(1259)*D*Y(16)*Y(43)+RATE(1260)*D*Y(16)*Y(43)&
    &+RATE(1261)*D*Y(16)*Y(43)+RATE(1325)*D*Y(18)*Y(43)+RATE(1530)*D*Y(21)&
    &*Y(43)+RATE(1628)*D*Y(32)*Y(43)+RATE(1846)*D*Y(2)*Y(43)+RATE(1892)*D*Y(4&
    &)*Y(43)+RATE(1992)*D*Y(6)*Y(43)+RATE(2021)*D*Y(8)*Y(43)+RATE(2185)*D&
    &*Y(10)*Y(43)+RATE(2408)*D*Y(13)*Y(43)+RATE(2409)*D*Y(13)*Y(43)+RATE(2565&
    &)*D*Y(28)*Y(43)+RATE(2687)*D*Y(36)*Y(43)+RATE(2697)*D*Y(43)*Y(68)&
    &+RATE(2698)*D*Y(43)*Y(83)+RATE(2699)*D*Y(43)*Y(100)+RATE(2700)*D*Y(43)&
    &*Y(63)+RATE(2701)*D*Y(43)*Y(104)+RATE(2702)*D*Y(43)*Y(162)+RATE(2703)*D&
    &*Y(43)*Y(57)+RATE(2704)*D*Y(43)*Y(68)+RATE(2705)*D*Y(43)*Y(76)+RATE(2706&
    &)*D*Y(43)*Y(81)+RATE(2707)*D*Y(43)*Y(53)+RATE(2708)*D*Y(43)*Y(100)&
    &+RATE(2709)*D*Y(43)*Y(132)+RATE(2710)*D*Y(43)*Y(63)+RATE(2711)*D*Y(43)&
    &*Y(142)+RATE(2712)*D*Y(43)*Y(66)+RATE(2713)*D*Y(43)*Y(91)+RATE(2714)*D&
    &*Y(43)*Y(102)+RATE(2715)*D*Y(43)*Y(102)+RATE(2716)*D*Y(43)*Y(117)&
    &+RATE(2717)*D*Y(43)*Y(145)+RATE(2718)*D*Y(43)*Y(120)+RATE(2719)*D*Y(43)&
    &*Y(55)+RATE(2720)*D*Y(43)*Y(176)+RATE(2721)*D*Y(43)*Y(57)+RATE(2722)*D&
    &*Y(43)*Y(41)+RATE(2723)*D*Y(43)*Y(133)+RATE(2724)*D*Y(43)*Y(133)&
    &+RATE(2725)*D*Y(43)*Y(56)+RATE(2726)*D*Y(43)*Y(56)+RATE(2727)*D*Y(43)&
    &*Y(131)+RATE(2749)*D*Y(44)*Y(43)+RATE(2890)*D*Y(46)*Y(43)+RATE(2891)*D&
    &*Y(46)*Y(43)+RATE(2928)*D*Y(48)*Y(43)
    PROD = RATE(173)*Y(54)+RATE(260)*Y(45)/safeMantle+RATE(343)*D*Y(45&
    &)/safeMantle*Y(2)+RATE(426)*Y(45)/safeMantle+RATE(740)*Y(3)*Y(37)&
    &+RATE(751)*Y(3)*Y(226)+RATE(761)*Y(37)*Y(226)+RATE(792)*Y(5)*Y(39)&
    &*bulkLayersReciprocal+RATE(803)*Y(5)*Y(228)*bulkLayersReciprocal&
    &+RATE(813)*Y(39)*Y(228)*bulkLayersReciprocal+RATE(925)*Y(54)+RATE(1108)&
    &*Y(45)+RATE(1191)*Y(50)+RATE(1439)*D*Y(20)*Y(44)+RATE(1548)*D*Y(24)*Y(44&
    &)+RATE(1566)*D*Y(24)*Y(55)+RATE(1629)*D*Y(32)*Y(54)+RATE(1841)*D*Y(2)&
    &*Y(144)+RATE(1847)*D*Y(2)*Y(54)+RATE(1993)*D*Y(6)*Y(35)+RATE(2642)*D&
    &*Y(35)*Y(41)+RATE(2644)*D*Y(35)*Y(62)+RATE(2645)*D*Y(35)*Y(54)+RATE(2645&
    &)*D*Y(35)*Y(54)+RATE(2648)*D*Y(35)*Y(35)+RATE(2658)*D*Y(35)*Y(56)&
    &+RATE(2678)*D*Y(36)*Y(131)+RATE(2728)*D*Y(44)*Y(183)+RATE(2729)*D*Y(44)&
    &*Y(116)+RATE(2730)*D*Y(44)*Y(54)+RATE(2731)*D*Y(44)*Y(133)+RATE(2732)*D&
    &*Y(44)*Y(163)+RATE(2771)*D*Y(54)*Y(100)+RATE(2778)*D*Y(54)*Y(91)&
    &+RATE(2792)*D*Y(54)*Y(82)+RATE(2797)*D*Y(55)*Y(333)+RATE(2804)*D*Y(55)&
    &*Y(54)+RATE(2805)*D*Y(64)*Y(333)+RATE(2806)*D*Y(64)*Y(333)+RATE(2892)*D&
    &*Y(46)*Y(54)+RATE(3007)*D*Y(56)*Y(90)+RATE(3011)*D*Y(56)*Y(54)
    YDOT(43) = PROD-LOSS
    LOSS = RATE(569)*D*Y(44)+RATE(1439)*D*Y(20)*Y(44)+RATE(1461)*D*Y(20)&
    &*Y(44)+RATE(1548)*D*Y(24)*Y(44)+RATE(1565)*D*Y(24)*Y(44)+RATE(1975)*D&
    &*Y(6)*Y(44)+RATE(2503)*D*Y(27)*Y(44)+RATE(2635)*D*Y(35)*Y(44)+RATE(2728)&
    &*D*Y(44)*Y(183)+RATE(2729)*D*Y(44)*Y(116)+RATE(2730)*D*Y(44)*Y(54)&
    &+RATE(2731)*D*Y(44)*Y(133)+RATE(2732)*D*Y(44)*Y(163)+RATE(2733)*D*Y(44)&
    &*Y(333)+RATE(2734)*D*Y(44)*Y(333)+RATE(2735)*D*Y(44)*Y(67)+RATE(2736)*D&
    &*Y(44)*Y(75)+RATE(2737)*D*Y(44)*Y(131)+RATE(2738)*D*Y(44)*Y(131)&
    &+RATE(2739)*D*Y(44)*Y(62)+RATE(2740)*D*Y(44)*Y(62)+RATE(2741)*D*Y(44)&
    &*Y(62)+RATE(2742)*D*Y(44)*Y(183)+RATE(2743)*D*Y(44)*Y(183)+RATE(2744)*D&
    &*Y(44)*Y(183)+RATE(2745)*D*Y(44)*Y(183)+RATE(2746)*D*Y(44)*Y(90)&
    &+RATE(2747)*D*Y(44)*Y(116)+RATE(2748)*D*Y(44)*Y(92)+RATE(2749)*D*Y(44)&
    &*Y(43)+RATE(2750)*D*Y(44)*Y(54)+RATE(2751)*D*Y(44)*Y(161)+RATE(2752)*D&
    &*Y(44)*Y(161)+RATE(2753)*D*Y(44)*Y(163)+RATE(2754)*D*Y(44)*Y(163)&
    &+RATE(2844)*D*Y(46)*Y(44)
    PROD = RATE(171)*Y(43)+RATE(923)*Y(43)+RATE(1892)*D*Y(4)*Y(43)&
    &+RATE(1933)*D*Y(4)*Y(227)+RATE(1974)*D*Y(6)*Y(36)+RATE(2021)*D*Y(8)*Y(43&
    &)+RATE(2048)*D*Y(8)*Y(35)+RATE(2187)*D*Y(10)*Y(35)+RATE(2411)*D*Y(13)&
    &*Y(54)+RATE(2565)*D*Y(28)*Y(43)+RATE(2590)*D*Y(28)*Y(54)+RATE(2627)*D&
    &*Y(35)*Y(53)+RATE(2631)*D*Y(35)*Y(91)+RATE(2632)*D*Y(35)*Y(117)&
    &+RATE(2633)*D*Y(35)*Y(145)+RATE(2634)*D*Y(35)*Y(120)+RATE(2639)*D*Y(35)&
    &*Y(176)+RATE(2640)*D*Y(35)*Y(57)+RATE(2682)*D*Y(36)*Y(62)+RATE(2689)*D&
    &*Y(36)*Y(35)+RATE(2697)*D*Y(43)*Y(68)+RATE(2698)*D*Y(43)*Y(83)+RATE(2699&
    &)*D*Y(43)*Y(100)+RATE(2700)*D*Y(43)*Y(63)+RATE(2701)*D*Y(43)*Y(104)&
    &+RATE(2702)*D*Y(43)*Y(162)+RATE(2703)*D*Y(43)*Y(57)+RATE(2928)*D*Y(48)&
    &*Y(43)
    YDOT(44) = PROD-LOSS
    LOSS = RATE(260)*Y(45)/safeMantle+RATE(343)*D*Y(45)/safeMantle*Y(2)&
    &+RATE(426)*Y(45)/safeMantle+RATE(637)*Y(3)*Y(45)+RATE(658)*Y(45)*Y(226)&
    &+RATE(659)*Y(45)*Y(118)+RATE(660)*Y(45)*Y(118)+RATE(661)*Y(45)*Y(130)&
    &+RATE(662)*Y(45)*Y(130)+RATE(741)*Y(3)*Y(45)+RATE(762)*Y(45)*Y(226)&
    &+RATE(763)*Y(45)*Y(118)+RATE(764)*Y(45)*Y(118)+RATE(765)*Y(45)*Y(130)&
    &+RATE(766)*Y(45)*Y(130)+RATE(1025)*Y(45)*totalSwap/safeMantle+RATE(1108)&
    &*Y(45)
    PROD = RATE(58)*Y(50)*bulkLayersReciprocal+RATE(98)*Y(52)+RATE(568)&
    &*D*Y(43)+RATE(569)*D*Y(44)+RATE(636)*Y(3)*Y(37)+RATE(647)*Y(3)*Y(226)&
    &+RATE(657)*Y(37)*Y(226)
    YDOT(45) = PROD-LOSS
    LOSS = RATE(108)*Y(46)+RATE(180)*Y(46)+RATE(579)*D*Y(46)+RATE(1281)&
    &*D*Y(16)*Y(46)+RATE(1345)*D*Y(18)*Y(46)+RATE(1431)*D*Y(20)*Y(46)&
    &+RATE(1491)*D*Y(20)*Y(46)+RATE(1492)*D*Y(20)*Y(46)+RATE(1536)*D*Y(21)&
    &*Y(46)+RATE(1591)*D*Y(24)*Y(46)+RATE(1592)*D*Y(24)*Y(46)+RATE(1593)*D&
    &*Y(24)*Y(46)+RATE(1594)*D*Y(24)*Y(46)+RATE(1613)*D*Y(25)*Y(46)+RATE(1636&
    &)*D*Y(32)*Y(46)+RATE(1637)*D*Y(32)*Y(46)+RATE(1659)*D*Y(33)*Y(46)&
    &+RATE(1660)*D*Y(33)*Y(46)+RATE(1868)*D*Y(2)*Y(46)+RATE(1898)*D*Y(4)*Y(46&
    &)+RATE(1996)*D*Y(6)*Y(46)+RATE(2051)*D*Y(8)*Y(46)+RATE(2192)*D*Y(10)&
    &*Y(46)+RATE(2193)*D*Y(10)*Y(46)+RATE(2654)*D*Y(35)*Y(46)+RATE(2655)*D&
    &*Y(35)*Y(46)+RATE(2693)*D*Y(36)*Y(46)+RATE(2824)*D*Y(46)*Y(83)+RATE(2825&
    &)*D*Y(46)*Y(100)+RATE(2826)*D*Y(46)*Y(104)+RATE(2827)*D*Y(46)*Y(68)&
    &+RATE(2828)*D*Y(46)*Y(76)+RATE(2829)*D*Y(46)*Y(81)+RATE(2830)*D*Y(46)&
    &*Y(42)+RATE(2831)*D*Y(46)*Y(53)+RATE(2832)*D*Y(46)*Y(53)+RATE(2833)*D&
    &*Y(46)*Y(234)+RATE(2834)*D*Y(46)*Y(63)+RATE(2835)*D*Y(46)*Y(184)&
    &+RATE(2836)*D*Y(46)*Y(184)+RATE(2837)*D*Y(46)*Y(243)+RATE(2838)*D*Y(46)&
    &*Y(245)+RATE(2839)*D*Y(46)*Y(245)+RATE(2840)*D*Y(46)*Y(174)+RATE(2841)*D&
    &*Y(46)*Y(174)+RATE(2842)*D*Y(46)*Y(104)+RATE(2843)*D*Y(46)*Y(120)&
    &+RATE(2844)*D*Y(46)*Y(44)+RATE(2845)*D*Y(46)*Y(55)+RATE(2846)*D*Y(46)&
    &*Y(266)+RATE(2847)*D*Y(46)*Y(176)+RATE(2848)*D*Y(46)*Y(57)+RATE(2849)*D&
    &*Y(46)*Y(208)+RATE(2850)*D*Y(46)*Y(122)+RATE(2851)*D*Y(46)*Y(136)&
    &+RATE(2852)*D*Y(46)*Y(147)+RATE(2853)*D*Y(46)*Y(236)+RATE(2854)*D*Y(46)&
    &*Y(67)+RATE(2855)*D*Y(46)*Y(80)+RATE(2856)*D*Y(46)*Y(89)+RATE(2857)*D&
    &*Y(46)*Y(109)+RATE(2858)*D*Y(46)*Y(109)+RATE(2859)*D*Y(46)*Y(109)&
    &+RATE(2860)*D*Y(46)*Y(109)+RATE(2861)*D*Y(46)*Y(123)+RATE(2862)*D*Y(46)&
    &*Y(123)+RATE(2863)*D*Y(46)*Y(75)+RATE(2864)*D*Y(46)*Y(197)+RATE(2865)*D&
    &*Y(46)*Y(284)+RATE(2866)*D*Y(46)*Y(308)+RATE(2867)*D*Y(46)*Y(41)&
    &+RATE(2868)*D*Y(46)*Y(82)+RATE(2869)*D*Y(46)*Y(82)+RATE(2870)*D*Y(46)&
    &*Y(232)+RATE(2871)*D*Y(46)*Y(233)+RATE(2872)*D*Y(46)*Y(233)+RATE(2873)*D&
    &*Y(46)*Y(101)+RATE(2874)*D*Y(46)*Y(131)+RATE(2875)*D*Y(46)*Y(62)&
    &+RATE(2876)*D*Y(46)*Y(183)+RATE(2877)*D*Y(46)*Y(90)+RATE(2878)*D*Y(46)&
    &*Y(90)+RATE(2879)*D*Y(46)*Y(90)+RATE(2880)*D*Y(46)*Y(116)+RATE(2881)*D&
    &*Y(46)*Y(116)+RATE(2882)*D*Y(46)*Y(244)+RATE(2883)*D*Y(46)*Y(244)&
    &+RATE(2884)*D*Y(46)*Y(144)+RATE(2885)*D*Y(46)*Y(144)+RATE(2886)*D*Y(46)&
    &*Y(144)+RATE(2887)*D*Y(46)*Y(173)+RATE(2888)*D*Y(46)*Y(173)+RATE(2889)*D&
    &*Y(46)*Y(103)+RATE(2890)*D*Y(46)*Y(43)+RATE(2891)*D*Y(46)*Y(43)&
    &+RATE(2892)*D*Y(46)*Y(54)+RATE(2893)*D*Y(46)*Y(264)+RATE(2894)*D*Y(46)&
    &*Y(133)+RATE(2895)*D*Y(46)*Y(265)+RATE(2896)*D*Y(46)*Y(265)+RATE(2897)*D&
    &*Y(46)*Y(175)+RATE(2898)*D*Y(46)*Y(221)+RATE(2899)*D*Y(46)*Y(221)&
    &+RATE(2900)*D*Y(46)*Y(300)+RATE(2901)*D*Y(46)*Y(300)+RATE(2902)*D*Y(46)&
    &*Y(56)+RATE(2903)*D*Y(46)*Y(316)+RATE(2904)*D*Y(46)*Y(320)+RATE(2905)*D&
    &*Y(46)*Y(277)+RATE(2906)*D*Y(46)*Y(292)+RATE(2907)*D*Y(46)*Y(314)&
    &+RATE(2908)*D*Y(46)*Y(207)+RATE(2909)*D*Y(46)*Y(207)+RATE(2910)*D*Y(46)&
    &*Y(135)+RATE(2911)*D*Y(46)*Y(135)+RATE(2912)*D*Y(46)*Y(146)+RATE(2913)*D&
    &*Y(46)*Y(166)+RATE(2914)*D*Y(46)*Y(121)+RATE(2915)*D*Y(46)*Y(46)&
    &+RATE(2915)*D*Y(46)*Y(46)+RATE(2916)*D*Y(46)*Y(277)+RATE(2917)*D*Y(46)&
    &*Y(106)+RATE(2918)*D*Y(46)*Y(105)
    PROD = RATE(139)*Y(99)+RATE(140)*Y(232)+RATE(177)*Y(133)+RATE(178)&
    &*Y(264)+RATE(182)*Y(161)+RATE(182)*Y(161)+RATE(184)*Y(221)+RATE(187)&
    &*Y(56)+RATE(198)*Y(235)+RATE(200)*Y(277)+RATE(202)*Y(320)+RATE(266)*Y(47&
    &)/safeMantle+RATE(349)*D*Y(47)/safeMantle*Y(2)+RATE(432)*Y(47&
    &)/safeMantle+RATE(874)*Y(99)+RATE(875)*Y(100)+RATE(876)*Y(232)+RATE(929)&
    &*Y(133)+RATE(930)*Y(264)+RATE(933)*Y(161)+RATE(933)*Y(161)+RATE(934)&
    &*Y(162)+RATE(936)*Y(175)+RATE(937)*Y(221)+RATE(940)*Y(56)+RATE(961)&
    &*Y(235)+RATE(963)*Y(236)+RATE(965)*Y(277)+RATE(967)*Y(320)+RATE(1114)&
    &*Y(47)+RATE(1197)*Y(51)+RATE(1237)*D*Y(16)*Y(162)+RATE(1239)*D*Y(16)&
    &*Y(57)+RATE(1251)*D*Y(16)*Y(99)+RATE(1264)*D*Y(16)*Y(133)+RATE(1268)*D&
    &*Y(16)*Y(161)+RATE(1272)*D*Y(16)*Y(56)+RATE(1275)*D*Y(16)*Y(277)&
    &+RATE(1329)*D*Y(18)*Y(161)+RATE(1335)*D*Y(18)*Y(277)+RATE(1440)*D*Y(20)&
    &*Y(48)+RATE(1464)*D*Y(20)*Y(162)+RATE(1466)*D*Y(20)*Y(57)+RATE(1482)*D&
    &*Y(20)*Y(133)+RATE(1486)*D*Y(20)*Y(161)+RATE(1488)*D*Y(20)*Y(161)&
    &+RATE(1534)*D*Y(21)*Y(161)+RATE(1549)*D*Y(24)*Y(48)+RATE(1567)*D*Y(24)&
    &*Y(162)+RATE(1569)*D*Y(24)*Y(57)+RATE(1589)*D*Y(24)*Y(161)+RATE(1597)*D&
    &*Y(24)*Y(56)+RATE(1638)*D*Y(32)*Y(56)+RATE(1658)*D*Y(33)*Y(161)&
    &+RATE(1696)*D*Y(41)*Y(57)+RATE(1759)*D*Y(82)*Y(161)+RATE(1785)*D*Y(99)&
    &*Y(161)+RATE(1793)*D*Y(100)*Y(333)+RATE(1801)*D*Y(2)*Y(161)+RATE(1801)*D&
    &*Y(2)*Y(161)+RATE(1802)*D*Y(2)*Y(56)+RATE(1808)*D*Y(2)*Y(48)+RATE(1838)&
    &*D*Y(2)*Y(116)+RATE(1841)*D*Y(2)*Y(144)+RATE(1850)*D*Y(2)*Y(133)&
    &+RATE(1854)*D*Y(2)*Y(161)+RATE(1855)*D*Y(2)*Y(175)+RATE(1858)*D*Y(2)&
    &*Y(221)+RATE(1862)*D*Y(2)*Y(56)+RATE(1864)*D*Y(2)*Y(277)+RATE(1924)*D&
    &*Y(4)*Y(232)+RATE(1949)*D*Y(6)*Y(161)+RATE(1949)*D*Y(6)*Y(161)+RATE(1950&
    &)*D*Y(6)*Y(56)+RATE(2063)*D*Y(132)*Y(333)+RATE(2119)*D*Y(63)*Y(333)&
    &+RATE(2120)*D*Y(63)*Y(333)+RATE(2217)*D*Y(66)*Y(333)+RATE(2316)*D*Y(243)&
    &*Y(333)+RATE(2366)*D*Y(13)*Y(232)+RATE(2370)*D*Y(13)*Y(99)+RATE(2375)*D&
    &*Y(13)*Y(131)+RATE(2394)*D*Y(13)*Y(116)+RATE(2414)*D*Y(13)*Y(133)&
    &+RATE(2417)*D*Y(13)*Y(161)+RATE(2418)*D*Y(13)*Y(221)+RATE(2420)*D*Y(13)&
    &*Y(300)+RATE(2427)*D*Y(13)*Y(320)+RATE(2428)*D*Y(13)*Y(277)+RATE(2441)*D&
    &*Y(13)*Y(235)+RATE(2473)*D*Y(326)*Y(333)+RATE(2504)*D*Y(27)*Y(162)&
    &+RATE(2506)*D*Y(27)*Y(278)+RATE(2532)*D*Y(27)*Y(116)+RATE(2540)*D*Y(27)&
    &*Y(264)+RATE(2540)*D*Y(27)*Y(264)+RATE(2543)*D*Y(27)*Y(133)+RATE(2545)*D&
    &*Y(27)*Y(161)+RATE(2548)*D*Y(27)*Y(56)+RATE(2550)*D*Y(27)*Y(277)&
    &+RATE(2592)*D*Y(28)*Y(133)+RATE(2593)*D*Y(28)*Y(161)+RATE(2623)*D*Y(35)&
    &*Y(48)+RATE(2638)*D*Y(35)*Y(162)+RATE(2640)*D*Y(35)*Y(57)+RATE(2650)*D&
    &*Y(35)*Y(133)+RATE(2652)*D*Y(35)*Y(161)+RATE(2658)*D*Y(35)*Y(56)&
    &+RATE(2681)*D*Y(36)*Y(62)+RATE(2690)*D*Y(36)*Y(133)+RATE(2721)*D*Y(43)&
    &*Y(57)+RATE(2726)*D*Y(43)*Y(56)+RATE(2741)*D*Y(44)*Y(62)+RATE(2751)*D&
    &*Y(44)*Y(161)+RATE(2818)*D*Y(133)*Y(161)+RATE(2820)*D*Y(133)*Y(163)&
    &+RATE(2822)*D*Y(134)*Y(333)+RATE(2919)*D*Y(48)*Y(67)+RATE(2920)*D*Y(48)&
    &*Y(80)+RATE(2921)*D*Y(48)*Y(75)+RATE(2922)*D*Y(48)*Y(41)+RATE(2923)*D&
    &*Y(48)*Y(99)+RATE(2924)*D*Y(48)*Y(131)+RATE(2925)*D*Y(48)*Y(62)&
    &+RATE(2926)*D*Y(48)*Y(183)+RATE(2927)*D*Y(48)*Y(116)+RATE(2928)*D*Y(48)&
    &*Y(43)+RATE(2929)*D*Y(48)*Y(54)+RATE(2930)*D*Y(48)*Y(161)+RATE(2931)*D&
    &*Y(48)*Y(300)+RATE(2932)*D*Y(48)*Y(56)+RATE(2933)*D*Y(48)*Y(320)&
    &+RATE(2953)*D*Y(48)*Y(333)+RATE(2956)*D*Y(161)*Y(234)+RATE(2957)*D*Y(161&
    &)*Y(165)+RATE(2962)*D*Y(161)*Y(163)+RATE(2963)*D*Y(161)*Y(277)+RATE(2967&
    &)*D*Y(162)*Y(333)+RATE(2967)*D*Y(162)*Y(333)+RATE(2971)*D*Y(162)*Y(163)&
    &+RATE(2975)*D*Y(301)*Y(333)+RATE(2982)*D*Y(56)*Y(100)+RATE(2983)*D*Y(56)&
    &*Y(63)+RATE(2989)*D*Y(56)*Y(55)+RATE(2998)*D*Y(56)*Y(82)+RATE(3014)*D&
    &*Y(56)*Y(56)+RATE(3029)*D*Y(57)*Y(333)+RATE(3030)*D*Y(57)*Y(67)&
    &+RATE(3031)*D*Y(57)*Y(75)+RATE(3032)*D*Y(57)*Y(82)+RATE(3033)*D*Y(57)&
    &*Y(232)+RATE(3034)*D*Y(57)*Y(99)+RATE(3035)*D*Y(57)*Y(131)+RATE(3036)*D&
    &*Y(57)*Y(62)+RATE(3037)*D*Y(57)*Y(183)+RATE(3038)*D*Y(57)*Y(90)&
    &+RATE(3040)*D*Y(57)*Y(116)+RATE(3041)*D*Y(57)*Y(92)+RATE(3042)*D*Y(57)&
    &*Y(103)+RATE(3043)*D*Y(57)*Y(54)+RATE(3044)*D*Y(57)*Y(133)+RATE(3045)*D&
    &*Y(57)*Y(56)+RATE(3046)*D*Y(57)*Y(163)+RATE(3048)*D*Y(57)*Y(105)&
    &+RATE(3049)*D*Y(57)*Y(121)+RATE(3050)*D*Y(57)*Y(235)+RATE(3061)*D*Y(163)&
    &*Y(277)+RATE(3080)*D*Y(105)*Y(161)+RATE(3103)*D*Y(236)*Y(333)+RATE(3107)&
    &*D*Y(278)*Y(333)+RATE(3115)*D*Y(319)*Y(333)+RATE(3115)*D*Y(319)*Y(333)&
    &+RATE(3116)*D*Y(319)*Y(333)
    YDOT(46) = PROD-LOSS
    LOSS = RATE(266)*Y(47)/safeMantle+RATE(349)*D*Y(47)/safeMantle*Y(2)&
    &+RATE(432)*Y(47)/safeMantle+RATE(633)*Y(3)*Y(47)+RATE(648)*Y(130)*Y(47)&
    &+RATE(663)*Y(47)*Y(47)+RATE(663)*Y(47)*Y(47)+RATE(664)*Y(47)*Y(222)&
    &+RATE(737)*Y(3)*Y(47)+RATE(752)*Y(130)*Y(47)+RATE(767)*Y(47)*Y(47)&
    &+RATE(767)*Y(47)*Y(47)+RATE(768)*Y(47)*Y(222)+RATE(1031)*Y(47)&
    &*totalSwap/safeMantle+RATE(1114)*Y(47)
    PROD = RATE(64)*Y(51)*bulkLayersReciprocal+RATE(579)*D*Y(46)&
    &+RATE(580)*D*Y(48)
    YDOT(47) = PROD-LOSS
    LOSS = RATE(580)*D*Y(48)+RATE(1280)*D*Y(16)*Y(48)+RATE(1440)*D*Y(20)&
    &*Y(48)+RATE(1463)*D*Y(20)*Y(48)+RATE(1549)*D*Y(24)*Y(48)+RATE(1808)*D&
    &*Y(2)*Y(48)+RATE(1977)*D*Y(6)*Y(48)+RATE(2623)*D*Y(35)*Y(48)+RATE(2637)&
    &*D*Y(35)*Y(48)+RATE(2919)*D*Y(48)*Y(67)+RATE(2920)*D*Y(48)*Y(80)&
    &+RATE(2921)*D*Y(48)*Y(75)+RATE(2922)*D*Y(48)*Y(41)+RATE(2923)*D*Y(48)&
    &*Y(99)+RATE(2924)*D*Y(48)*Y(131)+RATE(2925)*D*Y(48)*Y(62)+RATE(2926)*D&
    &*Y(48)*Y(183)+RATE(2927)*D*Y(48)*Y(116)+RATE(2928)*D*Y(48)*Y(43)&
    &+RATE(2929)*D*Y(48)*Y(54)+RATE(2930)*D*Y(48)*Y(161)+RATE(2931)*D*Y(48)&
    &*Y(300)+RATE(2932)*D*Y(48)*Y(56)+RATE(2933)*D*Y(48)*Y(320)+RATE(2934)*D&
    &*Y(48)*Y(67)+RATE(2935)*D*Y(48)*Y(109)+RATE(2936)*D*Y(48)*Y(75)&
    &+RATE(2937)*D*Y(48)*Y(159)+RATE(2938)*D*Y(48)*Y(159)+RATE(2939)*D*Y(48)&
    &*Y(41)+RATE(2940)*D*Y(48)*Y(82)+RATE(2941)*D*Y(48)*Y(232)+RATE(2942)*D&
    &*Y(48)*Y(131)+RATE(2943)*D*Y(48)*Y(183)+RATE(2944)*D*Y(48)*Y(183)&
    &+RATE(2945)*D*Y(48)*Y(90)+RATE(2946)*D*Y(48)*Y(90)+RATE(2947)*D*Y(48)&
    &*Y(116)+RATE(2948)*D*Y(48)*Y(103)+RATE(2949)*D*Y(48)*Y(264)+RATE(2950)*D&
    &*Y(48)*Y(300)+RATE(2951)*D*Y(48)*Y(56)+RATE(2952)*D*Y(48)*Y(320)&
    &+RATE(2953)*D*Y(48)*Y(333)
    PROD = RATE(108)*Y(46)+RATE(180)*Y(46)+RATE(934)*Y(162)+RATE(942)&
    &*Y(57)+RATE(1330)*D*Y(18)*Y(161)+RATE(1535)*D*Y(21)*Y(161)+RATE(1898)*D&
    &*Y(4)*Y(46)+RATE(2367)*D*Y(13)*Y(232)+RATE(2413)*D*Y(13)*Y(133)&
    &+RATE(2417)*D*Y(13)*Y(161)+RATE(2419)*D*Y(13)*Y(221)+RATE(2421)*D*Y(13)&
    &*Y(300)+RATE(2424)*D*Y(13)*Y(56)+RATE(2429)*D*Y(13)*Y(277)+RATE(2442)*D&
    &*Y(13)*Y(235)+RATE(2594)*D*Y(28)*Y(161)+RATE(2824)*D*Y(46)*Y(83)&
    &+RATE(2825)*D*Y(46)*Y(100)+RATE(2826)*D*Y(46)*Y(104)
    YDOT(48) = PROD-LOSS
    LOSS = RATE(23)*Y(49)*bulkLayersReciprocal+RATE(1156)*Y(49)
    PROD = RATE(674)*Y(38)*Y(228)*bulkLayersReciprocal+RATE(684)*Y(5)&
    &*Y(38)*bulkLayersReciprocal+RATE(990)*Y(40)*totalSwap/safeMantle
    YDOT(49) = PROD-LOSS
    LOSS = RATE(58)*Y(50)*bulkLayersReciprocal+RATE(689)*Y(5)*Y(50)&
    &*bulkLayersReciprocal+RATE(710)*Y(50)*Y(228)*bulkLayersReciprocal&
    &+RATE(711)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(712)*Y(50)*Y(126)&
    &*bulkLayersReciprocal+RATE(713)*Y(50)*Y(139)*bulkLayersReciprocal&
    &+RATE(714)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(793)*Y(5)*Y(50)&
    &*bulkLayersReciprocal+RATE(814)*Y(50)*Y(228)*bulkLayersReciprocal&
    &+RATE(815)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(816)*Y(50)*Y(126)&
    &*bulkLayersReciprocal+RATE(817)*Y(50)*Y(139)*bulkLayersReciprocal&
    &+RATE(818)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(1191)*Y(50)
    PROD = RATE(688)*Y(5)*Y(39)*bulkLayersReciprocal+RATE(699)*Y(5)&
    &*Y(228)*bulkLayersReciprocal+RATE(709)*Y(39)*Y(228)*bulkLayersReciprocal&
    &+RATE(1025)*Y(45)*totalSwap/safeMantle
    YDOT(50) = PROD-LOSS
    LOSS = RATE(64)*Y(51)*bulkLayersReciprocal+RATE(685)*Y(5)*Y(51)&
    &*bulkLayersReciprocal+RATE(700)*Y(139)*Y(51)*bulkLayersReciprocal&
    &+RATE(715)*Y(51)*Y(51)*bulkLayersReciprocal+RATE(715)*Y(51)*Y(51)&
    &*bulkLayersReciprocal+RATE(716)*Y(51)*Y(225)*bulkLayersReciprocal&
    &+RATE(789)*Y(5)*Y(51)*bulkLayersReciprocal+RATE(804)*Y(139)*Y(51)&
    &*bulkLayersReciprocal+RATE(819)*Y(51)*Y(51)*bulkLayersReciprocal&
    &+RATE(819)*Y(51)*Y(51)*bulkLayersReciprocal+RATE(820)*Y(51)*Y(225)&
    &*bulkLayersReciprocal+RATE(1197)*Y(51)
    PROD = RATE(1031)*Y(47)*totalSwap/safeMantle
    YDOT(51) = PROD-LOSS
    LOSS = RATE(97)*Y(52)+RATE(98)*Y(52)+RATE(262)*Y(52)/safeMantle&
    &+RATE(345)*D*Y(52)/safeMantle*Y(2)+RATE(428)*Y(52)/safeMantle+RATE(1027)&
    &*Y(52)*totalSwap/safeMantle+RATE(1110)*Y(52)
    PROD = RATE(60)*Y(59)*bulkLayersReciprocal+RATE(571)*D*Y(54)&
    &+RATE(572)*D*Y(55)+RATE(573)*D*Y(64)+RATE(637)*Y(3)*Y(45)+RATE(658)*Y(45&
    &)*Y(226)+RATE(660)*Y(45)*Y(118)+RATE(662)*Y(45)*Y(130)
    YDOT(52) = PROD-LOSS
    LOSS = RATE(493)*D*Y(53)+RATE(1226)*D*Y(16)*Y(53)+RATE(1446)*D*Y(20)&
    &*Y(53)+RATE(1553)*D*Y(24)*Y(53)+RATE(1721)*D*Y(53)*Y(333)+RATE(1722)*D&
    &*Y(53)*Y(333)+RATE(1723)*D*Y(53)*Y(333)+RATE(1724)*D*Y(53)*Y(333)&
    &+RATE(1725)*D*Y(53)*Y(333)+RATE(1726)*D*Y(53)*Y(67)+RATE(1727)*D*Y(53)&
    &*Y(75)+RATE(1728)*D*Y(53)*Y(232)+RATE(1729)*D*Y(53)*Y(99)+RATE(1730)*D&
    &*Y(53)*Y(131)+RATE(1731)*D*Y(53)*Y(62)+RATE(1732)*D*Y(53)*Y(183)&
    &+RATE(1733)*D*Y(53)*Y(90)+RATE(1734)*D*Y(53)*Y(116)+RATE(1735)*D*Y(53)&
    &*Y(193)+RATE(1736)*D*Y(53)*Y(92)+RATE(1737)*D*Y(53)*Y(69)+RATE(1738)*D&
    &*Y(53)*Y(163)+RATE(1739)*D*Y(53)*Y(166)+RATE(1814)*D*Y(2)*Y(53)&
    &+RATE(2627)*D*Y(35)*Y(53)+RATE(2707)*D*Y(43)*Y(53)+RATE(2770)*D*Y(54)&
    &*Y(53)+RATE(2831)*D*Y(46)*Y(53)+RATE(2832)*D*Y(46)*Y(53)+RATE(2981)*D&
    &*Y(56)*Y(53)
    PROD = RATE(1689)*D*Y(41)*Y(243)+RATE(1690)*D*Y(41)*Y(145)+RATE(1694&
    &)*D*Y(41)*Y(120)+RATE(1696)*D*Y(41)*Y(57)+RATE(1713)*D*Y(42)*Y(41)&
    &+RATE(1958)*D*Y(6)*Y(42)+RATE(2001)*D*Y(6)*Y(33)+RATE(2034)*D*Y(8)*Y(41)&
    &+RATE(2156)*D*Y(10)*Y(237)+RATE(2161)*D*Y(10)*Y(41)
    YDOT(53) = PROD-LOSS
    LOSS = RATE(173)*Y(54)+RATE(174)*Y(54)+RATE(175)*Y(54)+RATE(571)*D&
    &*Y(54)+RATE(925)*Y(54)+RATE(926)*Y(54)+RATE(927)*Y(54)+RATE(1291)*D*Y(18&
    &)*Y(54)+RATE(1326)*D*Y(18)*Y(54)+RATE(1501)*D*Y(21)*Y(54)+RATE(1531)*D&
    &*Y(21)*Y(54)+RATE(1611)*D*Y(25)*Y(54)+RATE(1629)*D*Y(32)*Y(54)+RATE(1657&
    &)*D*Y(33)*Y(54)+RATE(1678)*D*Y(172)*Y(54)+RATE(1707)*D*Y(42)*Y(54)&
    &+RATE(1719)*D*Y(42)*Y(54)+RATE(1847)*D*Y(2)*Y(54)+RATE(1893)*D*Y(4)*Y(54&
    &)+RATE(2022)*D*Y(8)*Y(54)+RATE(2186)*D*Y(10)*Y(54)+RATE(2329)*D*Y(13)&
    &*Y(54)+RATE(2410)*D*Y(13)*Y(54)+RATE(2411)*D*Y(13)*Y(54)+RATE(2566)*D&
    &*Y(28)*Y(54)+RATE(2589)*D*Y(28)*Y(54)+RATE(2590)*D*Y(28)*Y(54)+RATE(2645&
    &)*D*Y(35)*Y(54)+RATE(2663)*D*Y(36)*Y(54)+RATE(2688)*D*Y(36)*Y(54)&
    &+RATE(2730)*D*Y(44)*Y(54)+RATE(2750)*D*Y(44)*Y(54)+RATE(2755)*D*Y(54)&
    &*Y(81)+RATE(2756)*D*Y(54)*Y(100)+RATE(2757)*D*Y(54)*Y(132)+RATE(2758)*D&
    &*Y(54)*Y(63)+RATE(2759)*D*Y(54)*Y(184)+RATE(2760)*D*Y(54)*Y(91)&
    &+RATE(2761)*D*Y(54)*Y(174)+RATE(2762)*D*Y(54)*Y(104)+RATE(2763)*D*Y(54)&
    &*Y(162)+RATE(2764)*D*Y(54)*Y(301)+RATE(2765)*D*Y(54)*Y(165)+RATE(2766)*D&
    &*Y(54)*Y(278)+RATE(2767)*D*Y(54)*Y(76)+RATE(2768)*D*Y(54)*Y(81)&
    &+RATE(2769)*D*Y(54)*Y(199)+RATE(2770)*D*Y(54)*Y(53)+RATE(2771)*D*Y(54)&
    &*Y(100)+RATE(2772)*D*Y(54)*Y(132)+RATE(2773)*D*Y(54)*Y(63)+RATE(2774)*D&
    &*Y(54)*Y(184)+RATE(2775)*D*Y(54)*Y(142)+RATE(2776)*D*Y(54)*Y(66)&
    &+RATE(2777)*D*Y(54)*Y(189)+RATE(2778)*D*Y(54)*Y(91)+RATE(2779)*D*Y(54)&
    &*Y(102)+RATE(2780)*D*Y(54)*Y(102)+RATE(2781)*D*Y(54)*Y(117)+RATE(2782)*D&
    &*Y(54)*Y(243)+RATE(2783)*D*Y(54)*Y(245)+RATE(2784)*D*Y(54)*Y(145)&
    &+RATE(2785)*D*Y(54)*Y(174)+RATE(2786)*D*Y(54)*Y(326)+RATE(2787)*D*Y(54)&
    &*Y(307)+RATE(2788)*D*Y(54)*Y(120)+RATE(2789)*D*Y(54)*Y(176)+RATE(2790)*D&
    &*Y(54)*Y(122)+RATE(2791)*D*Y(54)*Y(247)+RATE(2792)*D*Y(54)*Y(82)&
    &+RATE(2804)*D*Y(55)*Y(54)+RATE(2892)*D*Y(46)*Y(54)+RATE(2929)*D*Y(48)&
    &*Y(54)+RATE(3011)*D*Y(56)*Y(54)+RATE(3025)*D*Y(57)*Y(54)+RATE(3043)*D&
    &*Y(57)*Y(54)
    PROD = RATE(262)*Y(52)/safeMantle+RATE(345)*D*Y(52)/safeMantle*Y(2)&
    &+RATE(428)*Y(52)/safeMantle+RATE(741)*Y(3)*Y(45)+RATE(762)*Y(45)*Y(226)&
    &+RATE(764)*Y(45)*Y(118)+RATE(766)*Y(45)*Y(130)+RATE(793)*Y(5)*Y(50)&
    &*bulkLayersReciprocal+RATE(814)*Y(50)*Y(228)*bulkLayersReciprocal&
    &+RATE(816)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(818)*Y(50)*Y(139)&
    &*bulkLayersReciprocal+RATE(1110)*Y(52)+RATE(1193)*Y(59)+RATE(1992)*D*Y(6&
    &)*Y(43)+RATE(2722)*D*Y(43)*Y(41)+RATE(2726)*D*Y(43)*Y(56)+RATE(2743)*D&
    &*Y(44)*Y(183)+RATE(2793)*D*Y(55)*Y(116)+RATE(2794)*D*Y(55)*Y(69)&
    &+RATE(2795)*D*Y(55)*Y(133)+RATE(2796)*D*Y(55)*Y(105)+RATE(2807)*D*Y(64)&
    &*Y(333)
    YDOT(54) = PROD-LOSS
    LOSS = RATE(572)*D*Y(55)+RATE(1462)*D*Y(20)*Y(55)+RATE(1566)*D*Y(24)&
    &*Y(55)+RATE(1695)*D*Y(41)*Y(55)+RATE(1976)*D*Y(6)*Y(55)+RATE(2636)*D&
    &*Y(35)*Y(55)+RATE(2719)*D*Y(43)*Y(55)+RATE(2793)*D*Y(55)*Y(116)&
    &+RATE(2794)*D*Y(55)*Y(69)+RATE(2795)*D*Y(55)*Y(133)+RATE(2796)*D*Y(55)&
    &*Y(105)+RATE(2797)*D*Y(55)*Y(333)+RATE(2798)*D*Y(55)*Y(333)+RATE(2799)*D&
    &*Y(55)*Y(67)+RATE(2800)*D*Y(55)*Y(131)+RATE(2801)*D*Y(55)*Y(62)&
    &+RATE(2802)*D*Y(55)*Y(183)+RATE(2803)*D*Y(55)*Y(116)+RATE(2804)*D*Y(55)&
    &*Y(54)+RATE(2845)*D*Y(46)*Y(55)+RATE(2989)*D*Y(56)*Y(55)
    PROD = RATE(174)*Y(54)+RATE(926)*Y(54)+RATE(1291)*D*Y(18)*Y(54)&
    &+RATE(1501)*D*Y(21)*Y(54)+RATE(1707)*D*Y(42)*Y(54)+RATE(1893)*D*Y(4)&
    &*Y(54)+RATE(1975)*D*Y(6)*Y(44)+RATE(2022)*D*Y(8)*Y(54)+RATE(2185)*D*Y(10&
    &)*Y(43)+RATE(2329)*D*Y(13)*Y(54)+RATE(2566)*D*Y(28)*Y(54)+RATE(2635)*D&
    &*Y(35)*Y(44)+RATE(2663)*D*Y(36)*Y(54)+RATE(2681)*D*Y(36)*Y(62)+RATE(2687&
    &)*D*Y(36)*Y(43)+RATE(2705)*D*Y(43)*Y(76)+RATE(2706)*D*Y(43)*Y(81)&
    &+RATE(2707)*D*Y(43)*Y(53)+RATE(2709)*D*Y(43)*Y(132)+RATE(2710)*D*Y(43)&
    &*Y(63)+RATE(2711)*D*Y(43)*Y(142)+RATE(2712)*D*Y(43)*Y(66)+RATE(2713)*D&
    &*Y(43)*Y(91)+RATE(2714)*D*Y(43)*Y(102)+RATE(2715)*D*Y(43)*Y(102)&
    &+RATE(2716)*D*Y(43)*Y(117)+RATE(2717)*D*Y(43)*Y(145)+RATE(2718)*D*Y(43)&
    &*Y(120)+RATE(2720)*D*Y(43)*Y(176)+RATE(2721)*D*Y(43)*Y(57)+RATE(2730)*D&
    &*Y(44)*Y(54)+RATE(2738)*D*Y(44)*Y(131)+RATE(2740)*D*Y(44)*Y(62)&
    &+RATE(2744)*D*Y(44)*Y(183)+RATE(2749)*D*Y(44)*Y(43)+RATE(2755)*D*Y(54)&
    &*Y(81)+RATE(2756)*D*Y(54)*Y(100)+RATE(2757)*D*Y(54)*Y(132)+RATE(2758)*D&
    &*Y(54)*Y(63)+RATE(2759)*D*Y(54)*Y(184)+RATE(2760)*D*Y(54)*Y(91)&
    &+RATE(2761)*D*Y(54)*Y(174)+RATE(2762)*D*Y(54)*Y(104)+RATE(2763)*D*Y(54)&
    &*Y(162)+RATE(2764)*D*Y(54)*Y(301)+RATE(2765)*D*Y(54)*Y(165)+RATE(2766)*D&
    &*Y(54)*Y(278)+RATE(2929)*D*Y(48)*Y(54)+RATE(3025)*D*Y(57)*Y(54)
    YDOT(55) = PROD-LOSS
    LOSS = RATE(187)*Y(56)+RATE(588)*D*Y(56)+RATE(940)*Y(56)+RATE(941)&
    &*Y(56)+RATE(1271)*D*Y(16)*Y(56)+RATE(1272)*D*Y(16)*Y(56)+RATE(1333)*D&
    &*Y(18)*Y(56)+RATE(1494)*D*Y(20)*Y(56)+RATE(1539)*D*Y(21)*Y(56)+RATE(1595&
    &)*D*Y(24)*Y(56)+RATE(1596)*D*Y(24)*Y(56)+RATE(1597)*D*Y(24)*Y(56)&
    &+RATE(1638)*D*Y(32)*Y(56)+RATE(1639)*D*Y(32)*Y(56)+RATE(1640)*D*Y(32)&
    &*Y(56)+RATE(1662)*D*Y(33)*Y(56)+RATE(1702)*D*Y(41)*Y(56)+RATE(1802)*D&
    &*Y(2)*Y(56)+RATE(1862)*D*Y(2)*Y(56)+RATE(1869)*D*Y(2)*Y(56)+RATE(1900)*D&
    &*Y(4)*Y(56)+RATE(1950)*D*Y(6)*Y(56)+RATE(1997)*D*Y(6)*Y(56)+RATE(2026)*D&
    &*Y(8)*Y(56)+RATE(2052)*D*Y(8)*Y(56)+RATE(2195)*D*Y(10)*Y(56)+RATE(2424)&
    &*D*Y(13)*Y(56)+RATE(2547)*D*Y(27)*Y(56)+RATE(2548)*D*Y(27)*Y(56)&
    &+RATE(2571)*D*Y(28)*Y(56)+RATE(2656)*D*Y(35)*Y(56)+RATE(2657)*D*Y(35)&
    &*Y(56)+RATE(2658)*D*Y(35)*Y(56)+RATE(2694)*D*Y(36)*Y(56)+RATE(2725)*D&
    &*Y(43)*Y(56)+RATE(2726)*D*Y(43)*Y(56)+RATE(2902)*D*Y(46)*Y(56)+RATE(2932&
    &)*D*Y(48)*Y(56)+RATE(2951)*D*Y(48)*Y(56)+RATE(2977)*D*Y(56)*Y(68)&
    &+RATE(2978)*D*Y(56)*Y(83)+RATE(2979)*D*Y(56)*Y(100)+RATE(2980)*D*Y(56)&
    &*Y(104)+RATE(2981)*D*Y(56)*Y(53)+RATE(2982)*D*Y(56)*Y(100)+RATE(2983)*D&
    &*Y(56)*Y(63)+RATE(2984)*D*Y(56)*Y(91)+RATE(2985)*D*Y(56)*Y(117)&
    &+RATE(2986)*D*Y(56)*Y(117)+RATE(2987)*D*Y(56)*Y(145)+RATE(2988)*D*Y(56)&
    &*Y(120)+RATE(2989)*D*Y(56)*Y(55)+RATE(2990)*D*Y(56)*Y(176)+RATE(2991)*D&
    &*Y(56)*Y(165)+RATE(2992)*D*Y(56)*Y(106)+RATE(2993)*D*Y(56)*Y(80)&
    &+RATE(2994)*D*Y(56)*Y(80)+RATE(2995)*D*Y(56)*Y(80)+RATE(2996)*D*Y(56)&
    &*Y(89)+RATE(2997)*D*Y(56)*Y(123)+RATE(2998)*D*Y(56)*Y(82)+RATE(2999)*D&
    &*Y(56)*Y(82)+RATE(3000)*D*Y(56)*Y(99)+RATE(3001)*D*Y(56)*Y(233)&
    &+RATE(3002)*D*Y(56)*Y(233)+RATE(3003)*D*Y(56)*Y(131)+RATE(3004)*D*Y(56)&
    &*Y(131)+RATE(3005)*D*Y(56)*Y(183)+RATE(3006)*D*Y(56)*Y(90)+RATE(3007)*D&
    &*Y(56)*Y(90)+RATE(3008)*D*Y(56)*Y(116)+RATE(3009)*D*Y(56)*Y(144)&
    &+RATE(3010)*D*Y(56)*Y(291)+RATE(3011)*D*Y(56)*Y(54)+RATE(3012)*D*Y(56)&
    &*Y(133)+RATE(3013)*D*Y(56)*Y(175)+RATE(3014)*D*Y(56)*Y(56)+RATE(3014)*D&
    &*Y(56)*Y(56)+RATE(3015)*D*Y(56)*Y(163)+RATE(3016)*D*Y(56)*Y(277)&
    &+RATE(3017)*D*Y(56)*Y(105)+RATE(3018)*D*Y(56)*Y(89)+RATE(3045)*D*Y(57)&
    &*Y(56)
    PROD = RATE(135)*Y(159)+RATE(147)*Y(62)+RATE(157)*Y(267)+RATE(271)&
    &*Y(58)/safeMantle+RATE(354)*D*Y(58)/safeMantle*Y(2)+RATE(437)*Y(58&
    &)/safeMantle+RATE(737)*Y(3)*Y(47)+RATE(789)*Y(5)*Y(51)&
    &*bulkLayersReciprocal+RATE(865)*Y(159)+RATE(889)*Y(62)+RATE(906)*Y(267)&
    &+RATE(936)*Y(175)+RATE(1119)*Y(58)+RATE(1202)*Y(60)+RATE(1227)*D*Y(16)&
    &*Y(63)+RATE(1442)*D*Y(20)*Y(57)+RATE(1449)*D*Y(20)*Y(63)+RATE(1487)*D&
    &*Y(20)*Y(161)+RATE(1489)*D*Y(20)*Y(175)+RATE(1492)*D*Y(20)*Y(46)&
    &+RATE(1533)*D*Y(21)*Y(161)+RATE(1551)*D*Y(24)*Y(57)+RATE(1556)*D*Y(24)&
    &*Y(63)+RATE(1584)*D*Y(24)*Y(133)+RATE(1590)*D*Y(24)*Y(161)+RATE(1594)*D&
    &*Y(24)*Y(46)+RATE(1612)*D*Y(25)*Y(161)+RATE(1624)*D*Y(32)*Y(62)&
    &+RATE(1632)*D*Y(32)*Y(161)+RATE(1675)*D*Y(172)*Y(333)+RATE(1786)*D*Y(99)&
    &*Y(175)+RATE(1800)*D*Y(2)*Y(62)+RATE(1819)*D*Y(2)*Y(319)+RATE(1830)*D&
    &*Y(2)*Y(232)+RATE(1831)*D*Y(2)*Y(99)+RATE(1834)*D*Y(2)*Y(62)+RATE(1843)&
    &*D*Y(2)*Y(144)+RATE(1849)*D*Y(2)*Y(264)+RATE(1851)*D*Y(2)*Y(133)&
    &+RATE(1854)*D*Y(2)*Y(161)+RATE(1857)*D*Y(2)*Y(175)+RATE(1857)*D*Y(2)&
    &*Y(175)+RATE(1860)*D*Y(2)*Y(221)+RATE(1865)*D*Y(2)*Y(277)+RATE(1868)*D&
    &*Y(2)*Y(46)+RATE(1936)*D*Y(4)*Y(264)+RATE(1947)*D*Y(6)*Y(62)+RATE(1995)&
    &*D*Y(6)*Y(161)+RATE(1995)*D*Y(6)*Y(161)+RATE(1996)*D*Y(6)*Y(46)&
    &+RATE(2080)*D*Y(62)*Y(68)+RATE(2082)*D*Y(62)*Y(199)+RATE(2085)*D*Y(62)&
    &*Y(83)+RATE(2087)*D*Y(62)*Y(100)+RATE(2100)*D*Y(62)*Y(104)+RATE(2121)*D&
    &*Y(63)*Y(333)+RATE(2122)*D*Y(63)*Y(67)+RATE(2123)*D*Y(63)*Y(75)&
    &+RATE(2124)*D*Y(63)*Y(99)+RATE(2125)*D*Y(63)*Y(131)+RATE(2126)*D*Y(63)&
    &*Y(62)+RATE(2127)*D*Y(63)*Y(183)+RATE(2129)*D*Y(63)*Y(90)+RATE(2131)*D&
    &*Y(63)*Y(116)+RATE(2132)*D*Y(63)*Y(92)+RATE(2133)*D*Y(63)*Y(163)&
    &+RATE(2135)*D*Y(63)*Y(320)+RATE(2188)*D*Y(10)*Y(264)+RATE(2207)*D*Y(142)&
    &*Y(333)+RATE(2218)*D*Y(66)*Y(333)+RATE(2219)*D*Y(66)*Y(333)+RATE(2291)*D&
    &*Y(116)*Y(161)+RATE(2317)*D*Y(243)*Y(333)+RATE(2358)*D*Y(13)*Y(159)&
    &+RATE(2380)*D*Y(13)*Y(62)+RATE(2461)*D*Y(306)*Y(333)+RATE(2474)*D*Y(326)&
    &*Y(333)+RATE(2644)*D*Y(35)*Y(62)+RATE(2651)*D*Y(35)*Y(133)+RATE(2653)*D&
    &*Y(35)*Y(161)+RATE(2655)*D*Y(35)*Y(46)+RATE(2682)*D*Y(36)*Y(62)&
    &+RATE(2691)*D*Y(36)*Y(161)+RATE(2703)*D*Y(43)*Y(57)+RATE(2710)*D*Y(43)&
    &*Y(63)+RATE(2724)*D*Y(43)*Y(133)+RATE(2740)*D*Y(44)*Y(62)+RATE(2752)*D&
    &*Y(44)*Y(161)+RATE(2773)*D*Y(54)*Y(63)+RATE(2801)*D*Y(55)*Y(62)&
    &+RATE(2830)*D*Y(46)*Y(42)+RATE(2835)*D*Y(46)*Y(184)+RATE(2840)*D*Y(46)&
    &*Y(174)+RATE(2857)*D*Y(46)*Y(109)+RATE(2867)*D*Y(46)*Y(41)+RATE(2874)*D&
    &*Y(46)*Y(131)+RATE(2875)*D*Y(46)*Y(62)+RATE(2875)*D*Y(46)*Y(62)&
    &+RATE(2876)*D*Y(46)*Y(183)+RATE(2877)*D*Y(46)*Y(90)+RATE(2881)*D*Y(46)&
    &*Y(116)+RATE(2885)*D*Y(46)*Y(144)+RATE(2887)*D*Y(46)*Y(173)+RATE(2891)*D&
    &*Y(46)*Y(43)+RATE(2892)*D*Y(46)*Y(54)+RATE(2897)*D*Y(46)*Y(175)&
    &+RATE(2913)*D*Y(46)*Y(166)+RATE(2938)*D*Y(48)*Y(159)+RATE(2939)*D*Y(48)&
    &*Y(41)+RATE(2942)*D*Y(48)*Y(131)+RATE(2943)*D*Y(48)*Y(183)+RATE(3019)*D&
    &*Y(57)*Y(67)+RATE(3020)*D*Y(57)*Y(75)+RATE(3021)*D*Y(57)*Y(131)&
    &+RATE(3022)*D*Y(57)*Y(62)+RATE(3023)*D*Y(57)*Y(183)+RATE(3024)*D*Y(57)&
    &*Y(116)+RATE(3025)*D*Y(57)*Y(54)+RATE(3026)*D*Y(57)*Y(133)+RATE(3027)*D&
    &*Y(57)*Y(161)+RATE(3028)*D*Y(57)*Y(163)+RATE(3095)*D*Y(136)*Y(161)&
    &+RATE(3104)*D*Y(247)*Y(333)
    YDOT(56) = PROD-LOSS
    LOSS = RATE(589)*D*Y(57)+RATE(942)*Y(57)+RATE(1239)*D*Y(16)*Y(57)&
    &+RATE(1442)*D*Y(20)*Y(57)+RATE(1466)*D*Y(20)*Y(57)+RATE(1551)*D*Y(24)&
    &*Y(57)+RATE(1569)*D*Y(24)*Y(57)+RATE(1696)*D*Y(41)*Y(57)+RATE(1697)*D&
    &*Y(41)*Y(57)+RATE(1979)*D*Y(6)*Y(57)+RATE(2505)*D*Y(27)*Y(57)+RATE(2640)&
    &*D*Y(35)*Y(57)+RATE(2703)*D*Y(43)*Y(57)+RATE(2721)*D*Y(43)*Y(57)&
    &+RATE(2848)*D*Y(46)*Y(57)+RATE(3019)*D*Y(57)*Y(67)+RATE(3020)*D*Y(57)&
    &*Y(75)+RATE(3021)*D*Y(57)*Y(131)+RATE(3022)*D*Y(57)*Y(62)+RATE(3023)*D&
    &*Y(57)*Y(183)+RATE(3024)*D*Y(57)*Y(116)+RATE(3025)*D*Y(57)*Y(54)&
    &+RATE(3026)*D*Y(57)*Y(133)+RATE(3027)*D*Y(57)*Y(161)+RATE(3028)*D*Y(57)&
    &*Y(163)+RATE(3029)*D*Y(57)*Y(333)+RATE(3030)*D*Y(57)*Y(67)+RATE(3031)*D&
    &*Y(57)*Y(75)+RATE(3032)*D*Y(57)*Y(82)+RATE(3033)*D*Y(57)*Y(232)&
    &+RATE(3034)*D*Y(57)*Y(99)+RATE(3035)*D*Y(57)*Y(131)+RATE(3036)*D*Y(57)&
    &*Y(62)+RATE(3037)*D*Y(57)*Y(183)+RATE(3038)*D*Y(57)*Y(90)+RATE(3039)*D&
    &*Y(57)*Y(116)+RATE(3040)*D*Y(57)*Y(116)+RATE(3041)*D*Y(57)*Y(92)&
    &+RATE(3042)*D*Y(57)*Y(103)+RATE(3043)*D*Y(57)*Y(54)+RATE(3044)*D*Y(57)&
    &*Y(133)+RATE(3045)*D*Y(57)*Y(56)+RATE(3046)*D*Y(57)*Y(163)+RATE(3047)*D&
    &*Y(57)*Y(163)+RATE(3048)*D*Y(57)*Y(105)+RATE(3049)*D*Y(57)*Y(121)&
    &+RATE(3050)*D*Y(57)*Y(235)
    PROD = RATE(890)*Y(63)+RATE(941)*Y(56)+RATE(1900)*D*Y(4)*Y(56)&
    &+RATE(1977)*D*Y(6)*Y(48)+RATE(2026)*D*Y(8)*Y(56)+RATE(2051)*D*Y(8)*Y(46)&
    &+RATE(2193)*D*Y(10)*Y(46)+RATE(2357)*D*Y(13)*Y(159)+RATE(2379)*D*Y(13)&
    &*Y(62)+RATE(2571)*D*Y(28)*Y(56)+RATE(2693)*D*Y(36)*Y(46)+RATE(2843)*D&
    &*Y(46)*Y(120)+RATE(2847)*D*Y(46)*Y(176)+RATE(2932)*D*Y(48)*Y(56)&
    &+RATE(2947)*D*Y(48)*Y(116)+RATE(2977)*D*Y(56)*Y(68)+RATE(2978)*D*Y(56)&
    &*Y(83)+RATE(2979)*D*Y(56)*Y(100)+RATE(2980)*D*Y(56)*Y(104)
    YDOT(57) = PROD-LOSS
    LOSS = RATE(271)*Y(58)/safeMantle+RATE(354)*D*Y(58)/safeMantle*Y(2)&
    &+RATE(437)*Y(58)/safeMantle+RATE(626)*Y(97)*Y(58)+RATE(634)*Y(3)*Y(58)&
    &+RATE(650)*Y(118)*Y(58)+RATE(665)*Y(58)*Y(226)+RATE(730)*Y(97)*Y(58)&
    &+RATE(738)*Y(3)*Y(58)+RATE(754)*Y(118)*Y(58)+RATE(769)*Y(58)*Y(226)&
    &+RATE(1036)*Y(58)*totalSwap/safeMantle+RATE(1119)*Y(58)
    PROD = RATE(69)*Y(60)*bulkLayersReciprocal+RATE(87)*Y(157)+RATE(91)&
    &*Y(61)+RATE(588)*D*Y(56)+RATE(589)*D*Y(57)+RATE(633)*Y(3)*Y(47)
    YDOT(58) = PROD-LOSS
    LOSS = RATE(60)*Y(59)*bulkLayersReciprocal+RATE(1193)*Y(59)
    PROD = RATE(689)*Y(5)*Y(50)*bulkLayersReciprocal+RATE(710)*Y(50)&
    &*Y(228)*bulkLayersReciprocal+RATE(712)*Y(50)*Y(126)*bulkLayersReciprocal&
    &+RATE(714)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(1027)*Y(52)&
    &*totalSwap/safeMantle
    YDOT(59) = PROD-LOSS
    LOSS = RATE(69)*Y(60)*bulkLayersReciprocal+RATE(678)*Y(111)*Y(60)&
    &*bulkLayersReciprocal+RATE(686)*Y(5)*Y(60)*bulkLayersReciprocal+RATE(702&
    &)*Y(126)*Y(60)*bulkLayersReciprocal+RATE(717)*Y(60)*Y(228)&
    &*bulkLayersReciprocal+RATE(782)*Y(111)*Y(60)*bulkLayersReciprocal&
    &+RATE(790)*Y(5)*Y(60)*bulkLayersReciprocal+RATE(806)*Y(126)*Y(60)&
    &*bulkLayersReciprocal+RATE(821)*Y(60)*Y(228)*bulkLayersReciprocal&
    &+RATE(1202)*Y(60)
    PROD = RATE(685)*Y(5)*Y(51)*bulkLayersReciprocal+RATE(1036)*Y(58)&
    &*totalSwap/safeMantle
    YDOT(60) = PROD-LOSS
    LOSS = RATE(91)*Y(61)+RATE(236)*Y(61)/safeMantle+RATE(319)*D*Y(61&
    &)/safeMantle*Y(2)+RATE(402)*Y(61)/safeMantle+RATE(1001)*Y(61)&
    &*totalSwap/safeMantle+RATE(1084)*Y(61)
    PROD = RATE(34)*Y(65)*bulkLayersReciprocal+RATE(515)*D*Y(62)&
    &+RATE(516)*D*Y(63)+RATE(526)*D*Y(66)+RATE(634)*Y(3)*Y(58)+RATE(665)*Y(58&
    &)*Y(226)
    YDOT(61) = PROD-LOSS
    LOSS = RATE(147)*Y(62)+RATE(515)*D*Y(62)+RATE(888)*Y(62)+RATE(889)&
    &*Y(62)+RATE(1316)*D*Y(18)*Y(62)+RATE(1317)*D*Y(18)*Y(62)+RATE(1519)*D&
    &*Y(21)*Y(62)+RATE(1520)*D*Y(21)*Y(62)+RATE(1521)*D*Y(21)*Y(62)+RATE(1606&
    &)*D*Y(25)*Y(62)+RATE(1624)*D*Y(32)*Y(62)+RATE(1666)*D*Y(33)*Y(62)&
    &+RATE(1717)*D*Y(42)*Y(62)+RATE(1731)*D*Y(53)*Y(62)+RATE(1800)*D*Y(2)&
    &*Y(62)+RATE(1834)*D*Y(2)*Y(62)+RATE(1883)*D*Y(4)*Y(62)+RATE(1947)*D*Y(6)&
    &*Y(62)+RATE(2017)*D*Y(8)*Y(62)+RATE(2041)*D*Y(8)*Y(62)+RATE(2077)*D*Y(62&
    &)*Y(100)+RATE(2078)*D*Y(62)*Y(91)+RATE(2079)*D*Y(62)*Y(104)+RATE(2080)*D&
    &*Y(62)*Y(68)+RATE(2081)*D*Y(62)*Y(81)+RATE(2082)*D*Y(62)*Y(199)&
    &+RATE(2083)*D*Y(62)*Y(199)+RATE(2084)*D*Y(62)*Y(309)+RATE(2085)*D*Y(62)&
    &*Y(83)+RATE(2086)*D*Y(62)*Y(83)+RATE(2087)*D*Y(62)*Y(100)+RATE(2088)*D&
    &*Y(62)*Y(132)+RATE(2089)*D*Y(62)*Y(196)+RATE(2090)*D*Y(62)*Y(184)&
    &+RATE(2091)*D*Y(62)*Y(142)+RATE(2092)*D*Y(62)*Y(91)+RATE(2093)*D*Y(62)&
    &*Y(117)+RATE(2094)*D*Y(62)*Y(243)+RATE(2095)*D*Y(62)*Y(145)+RATE(2096)*D&
    &*Y(62)*Y(306)+RATE(2097)*D*Y(62)*Y(174)+RATE(2098)*D*Y(62)*Y(326)&
    &+RATE(2099)*D*Y(62)*Y(307)+RATE(2100)*D*Y(62)*Y(104)+RATE(2101)*D*Y(62)&
    &*Y(120)+RATE(2102)*D*Y(62)*Y(176)+RATE(2103)*D*Y(62)*Y(106)+RATE(2104)*D&
    &*Y(62)*Y(122)+RATE(2105)*D*Y(62)*Y(167)+RATE(2106)*D*Y(62)*Y(177)&
    &+RATE(2126)*D*Y(63)*Y(62)+RATE(2171)*D*Y(10)*Y(62)+RATE(2326)*D*Y(13)&
    &*Y(62)+RATE(2379)*D*Y(13)*Y(62)+RATE(2380)*D*Y(13)*Y(62)+RATE(2560)*D&
    &*Y(28)*Y(62)+RATE(2644)*D*Y(35)*Y(62)+RATE(2662)*D*Y(36)*Y(62)+RATE(2679&
    &)*D*Y(36)*Y(62)+RATE(2680)*D*Y(36)*Y(62)+RATE(2681)*D*Y(36)*Y(62)&
    &+RATE(2682)*D*Y(36)*Y(62)+RATE(2739)*D*Y(44)*Y(62)+RATE(2740)*D*Y(44)&
    &*Y(62)+RATE(2741)*D*Y(44)*Y(62)+RATE(2801)*D*Y(55)*Y(62)+RATE(2875)*D&
    &*Y(46)*Y(62)+RATE(2925)*D*Y(48)*Y(62)+RATE(3022)*D*Y(57)*Y(62)+RATE(3036&
    &)*D*Y(57)*Y(62)
    PROD = RATE(236)*Y(61)/safeMantle+RATE(319)*D*Y(61)/safeMantle*Y(2)&
    &+RATE(402)*Y(61)/safeMantle+RATE(738)*Y(3)*Y(58)+RATE(769)*Y(58)*Y(226)&
    &+RATE(790)*Y(5)*Y(60)*bulkLayersReciprocal+RATE(821)*Y(60)*Y(228)&
    &*bulkLayersReciprocal+RATE(1084)*Y(61)+RATE(1167)*Y(65)+RATE(1436)*D&
    &*Y(20)*Y(63)+RATE(1451)*D*Y(20)*Y(66)+RATE(1546)*D*Y(24)*Y(63)+RATE(1557&
    &)*D*Y(24)*Y(66)+RATE(1588)*D*Y(24)*Y(161)+RATE(1596)*D*Y(24)*Y(56)&
    &+RATE(1631)*D*Y(32)*Y(133)+RATE(1633)*D*Y(32)*Y(161)+RATE(1640)*D*Y(32)&
    &*Y(56)+RATE(1673)*D*Y(172)*Y(333)+RATE(1674)*D*Y(172)*Y(333)+RATE(1702)&
    &*D*Y(41)*Y(56)+RATE(1855)*D*Y(2)*Y(175)+RATE(1869)*D*Y(2)*Y(56)&
    &+RATE(1920)*D*Y(4)*Y(159)+RATE(1997)*D*Y(6)*Y(56)+RATE(2107)*D*Y(63)&
    &*Y(67)+RATE(2108)*D*Y(63)*Y(80)+RATE(2109)*D*Y(63)*Y(75)+RATE(2110)*D&
    &*Y(63)*Y(131)+RATE(2111)*D*Y(63)*Y(183)+RATE(2112)*D*Y(63)*Y(116)&
    &+RATE(2113)*D*Y(63)*Y(69)+RATE(2114)*D*Y(63)*Y(133)+RATE(2115)*D*Y(63)&
    &*Y(161)+RATE(2116)*D*Y(63)*Y(300)+RATE(2117)*D*Y(63)*Y(163)+RATE(2118)*D&
    &*Y(63)*Y(105)+RATE(2140)*D*Y(183)*Y(278)+RATE(2159)*D*Y(10)*Y(159)&
    &+RATE(2176)*D*Y(10)*Y(267)+RATE(2208)*D*Y(142)*Y(333)+RATE(2216)*D*Y(66)&
    &*Y(333)+RATE(2220)*D*Y(66)*Y(67)+RATE(2221)*D*Y(66)*Y(206)+RATE(2222)*D&
    &*Y(66)*Y(214)+RATE(2223)*D*Y(66)*Y(159)+RATE(2224)*D*Y(66)*Y(233)&
    &+RATE(2225)*D*Y(66)*Y(131)+RATE(2226)*D*Y(66)*Y(183)+RATE(2227)*D*Y(66)&
    &*Y(90)+RATE(2228)*D*Y(66)*Y(92)+RATE(2229)*D*Y(66)*Y(324)+RATE(2230)*D&
    &*Y(66)*Y(316)+RATE(2231)*D*Y(66)*Y(105)+RATE(2232)*D*Y(66)*Y(135)&
    &+RATE(2233)*D*Y(66)*Y(121)+RATE(2234)*D*Y(66)*Y(235)+RATE(2245)*D*Y(90)&
    &*Y(172)+RATE(2656)*D*Y(35)*Y(56)+RATE(2700)*D*Y(43)*Y(63)+RATE(2712)*D&
    &*Y(43)*Y(66)+RATE(2723)*D*Y(43)*Y(133)+RATE(2725)*D*Y(43)*Y(56)&
    &+RATE(2758)*D*Y(54)*Y(63)+RATE(2776)*D*Y(54)*Y(66)+RATE(2935)*D*Y(48)&
    &*Y(109)+RATE(2937)*D*Y(48)*Y(159)+RATE(2944)*D*Y(48)*Y(183)+RATE(2993)*D&
    &*Y(56)*Y(80)+RATE(2996)*D*Y(56)*Y(89)+RATE(2997)*D*Y(56)*Y(123)&
    &+RATE(3003)*D*Y(56)*Y(131)+RATE(3005)*D*Y(56)*Y(183)+RATE(3006)*D*Y(56)&
    &*Y(90)+RATE(3008)*D*Y(56)*Y(116)+RATE(3009)*D*Y(56)*Y(144)+RATE(3011)*D&
    &*Y(56)*Y(54)+RATE(3013)*D*Y(56)*Y(175)+RATE(3014)*D*Y(56)*Y(56)
    YDOT(62) = PROD-LOSS
    LOSS = RATE(516)*D*Y(63)+RATE(890)*Y(63)+RATE(1227)*D*Y(16)*Y(63)&
    &+RATE(1436)*D*Y(20)*Y(63)+RATE(1449)*D*Y(20)*Y(63)+RATE(1546)*D*Y(24)&
    &*Y(63)+RATE(1556)*D*Y(24)*Y(63)+RATE(1687)*D*Y(41)*Y(63)+RATE(1964)*D&
    &*Y(6)*Y(63)+RATE(2107)*D*Y(63)*Y(67)+RATE(2108)*D*Y(63)*Y(80)+RATE(2109)&
    &*D*Y(63)*Y(75)+RATE(2110)*D*Y(63)*Y(131)+RATE(2111)*D*Y(63)*Y(183)&
    &+RATE(2112)*D*Y(63)*Y(116)+RATE(2113)*D*Y(63)*Y(69)+RATE(2114)*D*Y(63)&
    &*Y(133)+RATE(2115)*D*Y(63)*Y(161)+RATE(2116)*D*Y(63)*Y(300)+RATE(2117)*D&
    &*Y(63)*Y(163)+RATE(2118)*D*Y(63)*Y(105)+RATE(2119)*D*Y(63)*Y(333)&
    &+RATE(2120)*D*Y(63)*Y(333)+RATE(2121)*D*Y(63)*Y(333)+RATE(2122)*D*Y(63)&
    &*Y(67)+RATE(2123)*D*Y(63)*Y(75)+RATE(2124)*D*Y(63)*Y(99)+RATE(2125)*D&
    &*Y(63)*Y(131)+RATE(2126)*D*Y(63)*Y(62)+RATE(2127)*D*Y(63)*Y(183)&
    &+RATE(2128)*D*Y(63)*Y(183)+RATE(2129)*D*Y(63)*Y(90)+RATE(2130)*D*Y(63)&
    &*Y(116)+RATE(2131)*D*Y(63)*Y(116)+RATE(2132)*D*Y(63)*Y(92)+RATE(2133)*D&
    &*Y(63)*Y(163)+RATE(2134)*D*Y(63)*Y(163)+RATE(2135)*D*Y(63)*Y(320)&
    &+RATE(2498)*D*Y(27)*Y(63)+RATE(2499)*D*Y(27)*Y(63)+RATE(2630)*D*Y(35)&
    &*Y(63)+RATE(2700)*D*Y(43)*Y(63)+RATE(2710)*D*Y(43)*Y(63)+RATE(2758)*D&
    &*Y(54)*Y(63)+RATE(2773)*D*Y(54)*Y(63)+RATE(2834)*D*Y(46)*Y(63)+RATE(2983&
    &)*D*Y(56)*Y(63)
    PROD = RATE(888)*Y(62)+RATE(1883)*D*Y(4)*Y(62)+RATE(1979)*D*Y(6)&
    &*Y(57)+RATE(2017)*D*Y(8)*Y(62)+RATE(2052)*D*Y(8)*Y(56)+RATE(2077)*D*Y(62&
    &)*Y(100)+RATE(2078)*D*Y(62)*Y(91)+RATE(2079)*D*Y(62)*Y(104)+RATE(2192)*D&
    &*Y(10)*Y(46)+RATE(2195)*D*Y(10)*Y(56)+RATE(2326)*D*Y(13)*Y(62)+RATE(2560&
    &)*D*Y(28)*Y(62)+RATE(2662)*D*Y(36)*Y(62)+RATE(2694)*D*Y(36)*Y(56)&
    &+RATE(2925)*D*Y(48)*Y(62)+RATE(2981)*D*Y(56)*Y(53)+RATE(2984)*D*Y(56)&
    &*Y(91)+RATE(2985)*D*Y(56)*Y(117)+RATE(2987)*D*Y(56)*Y(145)+RATE(2988)*D&
    &*Y(56)*Y(120)+RATE(2990)*D*Y(56)*Y(176)+RATE(3022)*D*Y(57)*Y(62)&
    &+RATE(3039)*D*Y(57)*Y(116)+RATE(3045)*D*Y(57)*Y(56)
    YDOT(63) = PROD-LOSS
    LOSS = RATE(573)*D*Y(64)+RATE(2805)*D*Y(64)*Y(333)+RATE(2806)*D*Y(64&
    &)*Y(333)+RATE(2807)*D*Y(64)*Y(333)
    PROD = RATE(1462)*D*Y(20)*Y(55)+RATE(1531)*D*Y(21)*Y(54)+RATE(1611)&
    &*D*Y(25)*Y(54)+RATE(1657)*D*Y(33)*Y(54)+RATE(1678)*D*Y(172)*Y(54)&
    &+RATE(1695)*D*Y(41)*Y(55)+RATE(1719)*D*Y(42)*Y(54)+RATE(1976)*D*Y(6)&
    &*Y(55)+RATE(2186)*D*Y(10)*Y(54)+RATE(2636)*D*Y(35)*Y(55)+RATE(2688)*D&
    &*Y(36)*Y(54)+RATE(2719)*D*Y(43)*Y(55)+RATE(2741)*D*Y(44)*Y(62)+RATE(2745&
    &)*D*Y(44)*Y(183)+RATE(2750)*D*Y(44)*Y(54)+RATE(2767)*D*Y(54)*Y(76)&
    &+RATE(2768)*D*Y(54)*Y(81)+RATE(2770)*D*Y(54)*Y(53)+RATE(2772)*D*Y(54)&
    &*Y(132)+RATE(2773)*D*Y(54)*Y(63)+RATE(2774)*D*Y(54)*Y(184)+RATE(2775)*D&
    &*Y(54)*Y(142)+RATE(2776)*D*Y(54)*Y(66)+RATE(2777)*D*Y(54)*Y(189)&
    &+RATE(2779)*D*Y(54)*Y(102)+RATE(2780)*D*Y(54)*Y(102)+RATE(2781)*D*Y(54)&
    &*Y(117)+RATE(2782)*D*Y(54)*Y(243)+RATE(2783)*D*Y(54)*Y(245)+RATE(2784)*D&
    &*Y(54)*Y(145)+RATE(2785)*D*Y(54)*Y(174)+RATE(2786)*D*Y(54)*Y(326)&
    &+RATE(2787)*D*Y(54)*Y(307)+RATE(2788)*D*Y(54)*Y(120)+RATE(2789)*D*Y(54)&
    &*Y(176)+RATE(2790)*D*Y(54)*Y(122)+RATE(2791)*D*Y(54)*Y(247)+RATE(2800)*D&
    &*Y(55)*Y(131)+RATE(2801)*D*Y(55)*Y(62)+RATE(2802)*D*Y(55)*Y(183)&
    &+RATE(2803)*D*Y(55)*Y(116)+RATE(2804)*D*Y(55)*Y(54)+RATE(2989)*D*Y(56)&
    &*Y(55)+RATE(3043)*D*Y(57)*Y(54)
    YDOT(64) = PROD-LOSS
    LOSS = RATE(34)*Y(65)*bulkLayersReciprocal+RATE(1167)*Y(65)
    PROD = RATE(686)*Y(5)*Y(60)*bulkLayersReciprocal+RATE(717)*Y(60)&
    &*Y(228)*bulkLayersReciprocal+RATE(1001)*Y(61)*totalSwap/safeMantle
    YDOT(65) = PROD-LOSS
    LOSS = RATE(526)*D*Y(66)+RATE(1229)*D*Y(16)*Y(66)+RATE(1451)*D*Y(20)&
    &*Y(66)+RATE(1557)*D*Y(24)*Y(66)+RATE(2216)*D*Y(66)*Y(333)+RATE(2217)*D&
    &*Y(66)*Y(333)+RATE(2218)*D*Y(66)*Y(333)+RATE(2219)*D*Y(66)*Y(333)&
    &+RATE(2220)*D*Y(66)*Y(67)+RATE(2221)*D*Y(66)*Y(206)+RATE(2222)*D*Y(66)&
    &*Y(214)+RATE(2223)*D*Y(66)*Y(159)+RATE(2224)*D*Y(66)*Y(233)+RATE(2225)*D&
    &*Y(66)*Y(131)+RATE(2226)*D*Y(66)*Y(183)+RATE(2227)*D*Y(66)*Y(90)&
    &+RATE(2228)*D*Y(66)*Y(92)+RATE(2229)*D*Y(66)*Y(324)+RATE(2230)*D*Y(66)&
    &*Y(316)+RATE(2231)*D*Y(66)*Y(105)+RATE(2232)*D*Y(66)*Y(135)+RATE(2233)*D&
    &*Y(66)*Y(121)+RATE(2234)*D*Y(66)*Y(235)+RATE(2712)*D*Y(43)*Y(66)&
    &+RATE(2776)*D*Y(54)*Y(66)
    PROD = RATE(1520)*D*Y(21)*Y(62)+RATE(1687)*D*Y(41)*Y(63)+RATE(1697)&
    &*D*Y(41)*Y(57)+RATE(1717)*D*Y(42)*Y(62)+RATE(1731)*D*Y(53)*Y(62)&
    &+RATE(1964)*D*Y(6)*Y(63)+RATE(2041)*D*Y(8)*Y(62)+RATE(2081)*D*Y(62)*Y(81&
    &)+RATE(2088)*D*Y(62)*Y(132)+RATE(2089)*D*Y(62)*Y(196)+RATE(2090)*D*Y(62)&
    &*Y(184)+RATE(2091)*D*Y(62)*Y(142)+RATE(2092)*D*Y(62)*Y(91)+RATE(2093)*D&
    &*Y(62)*Y(117)+RATE(2094)*D*Y(62)*Y(243)+RATE(2095)*D*Y(62)*Y(145)&
    &+RATE(2096)*D*Y(62)*Y(306)+RATE(2097)*D*Y(62)*Y(174)+RATE(2098)*D*Y(62)&
    &*Y(326)+RATE(2101)*D*Y(62)*Y(120)+RATE(2102)*D*Y(62)*Y(176)+RATE(2104)*D&
    &*Y(62)*Y(122)+RATE(2105)*D*Y(62)*Y(167)+RATE(2106)*D*Y(62)*Y(177)&
    &+RATE(2126)*D*Y(63)*Y(62)+RATE(2128)*D*Y(63)*Y(183)+RATE(2130)*D*Y(63)&
    &*Y(116)+RATE(2157)*D*Y(10)*Y(237)+RATE(2171)*D*Y(10)*Y(62)+RATE(2175)*D&
    &*Y(10)*Y(267)+RATE(2630)*D*Y(35)*Y(63)+RATE(2679)*D*Y(36)*Y(62)&
    &+RATE(2739)*D*Y(44)*Y(62)+RATE(2832)*D*Y(46)*Y(53)+RATE(2983)*D*Y(56)&
    &*Y(63)+RATE(3036)*D*Y(57)*Y(62)
    YDOT(66) = PROD-LOSS
    LOSS = RATE(111)*Y(67)+RATE(456)*D*Y(67)+RATE(831)*Y(67)+RATE(832)&
    &*Y(67)+RATE(1348)*D*Y(67)*Y(83)+RATE(1349)*D*Y(67)*Y(100)+RATE(1350)*D&
    &*Y(67)*Y(104)+RATE(1351)*D*Y(67)*Y(162)+RATE(1352)*D*Y(67)*Y(132)&
    &+RATE(1353)*D*Y(67)*Y(91)+RATE(1354)*D*Y(67)*Y(117)+RATE(1355)*D*Y(67)&
    &*Y(145)+RATE(1356)*D*Y(67)*Y(120)+RATE(1357)*D*Y(67)*Y(162)+RATE(1358)*D&
    &*Y(67)*Y(176)+RATE(1359)*D*Y(67)*Y(165)+RATE(1360)*D*Y(67)*Y(236)&
    &+RATE(1361)*D*Y(67)*Y(80)+RATE(1362)*D*Y(67)*Y(90)+RATE(1363)*D*Y(67)&
    &*Y(161)+RATE(1364)*D*Y(67)*Y(163)+RATE(1369)*D*Y(68)*Y(67)+RATE(1506)*D&
    &*Y(21)*Y(67)+RATE(1726)*D*Y(53)*Y(67)+RATE(1822)*D*Y(2)*Y(67)+RATE(1872)&
    &*D*Y(4)*Y(67)+RATE(2008)*D*Y(8)*Y(67)+RATE(2028)*D*Y(8)*Y(67)+RATE(2107)&
    &*D*Y(63)*Y(67)+RATE(2122)*D*Y(63)*Y(67)+RATE(2149)*D*Y(10)*Y(67)&
    &+RATE(2220)*D*Y(66)*Y(67)+RATE(2320)*D*Y(13)*Y(67)+RATE(2333)*D*Y(13)&
    &*Y(67)+RATE(2510)*D*Y(27)*Y(67)+RATE(2553)*D*Y(28)*Y(67)+RATE(2668)*D&
    &*Y(36)*Y(67)+RATE(2669)*D*Y(36)*Y(67)+RATE(2670)*D*Y(36)*Y(67)+RATE(2735&
    &)*D*Y(44)*Y(67)+RATE(2799)*D*Y(55)*Y(67)+RATE(2854)*D*Y(46)*Y(67)&
    &+RATE(2919)*D*Y(48)*Y(67)+RATE(2934)*D*Y(48)*Y(67)+RATE(3019)*D*Y(57)&
    &*Y(67)+RATE(3030)*D*Y(57)*Y(67)
    PROD = RATE(112)*Y(75)+RATE(119)*Y(197)+RATE(121)*Y(284)+RATE(122)&
    &*Y(280)+RATE(204)*Y(71)/safeMantle+RATE(287)*D*Y(71)/safeMantle*Y(2)&
    &+RATE(370)*Y(71)/safeMantle+RATE(834)*Y(75)+RATE(842)*Y(197)+RATE(844)&
    &*Y(284)+RATE(845)*Y(280)+RATE(948)*Y(292)+RATE(949)*Y(314)+RATE(1052)&
    &*Y(71)+RATE(1135)*Y(73)+RATE(1217)*D*Y(16)*Y(68)+RATE(1244)*D*Y(16)&
    &*Y(197)+RATE(1249)*D*Y(16)*Y(20)+RATE(1250)*D*Y(16)*Y(82)+RATE(1251)*D&
    &*Y(16)*Y(99)+RATE(1252)*D*Y(16)*Y(233)+RATE(1278)*D*Y(16)*Y(16)&
    &+RATE(1338)*D*Y(18)*Y(207)+RATE(1365)*D*Y(68)*Y(116)+RATE(1366)*D*Y(68)&
    &*Y(133)+RATE(1367)*D*Y(68)*Y(163)+RATE(1376)*D*Y(75)*Y(100)+RATE(1391)*D&
    &*Y(76)*Y(333)+RATE(1393)*D*Y(76)*Y(131)+RATE(1395)*D*Y(76)*Y(90)&
    &+RATE(1397)*D*Y(76)*Y(92)+RATE(1406)*D*Y(81)*Y(333)+RATE(1422)*D*Y(199)&
    &*Y(333)+RATE(1427)*D*Y(192)*Y(333)+RATE(1429)*D*Y(309)*Y(333)+RATE(1432)&
    &*D*Y(20)*Y(68)+RATE(1444)*D*Y(20)*Y(76)+RATE(1542)*D*Y(24)*Y(68)&
    &+RATE(1552)*D*Y(24)*Y(76)+RATE(1749)*D*Y(82)*Y(82)+RATE(2347)*D*Y(13)&
    &*Y(280)+RATE(2430)*D*Y(13)*Y(292)+RATE(2697)*D*Y(43)*Y(68)+RATE(2705)*D&
    &*Y(43)*Y(76)+RATE(2767)*D*Y(54)*Y(76)+RATE(2977)*D*Y(56)*Y(68)+RATE(3085&
    &)*D*Y(293)*Y(333)+RATE(3088)*D*Y(313)*Y(333)
    YDOT(67) = PROD-LOSS
    LOSS = RATE(457)*D*Y(68)+RATE(833)*Y(68)+RATE(1217)*D*Y(16)*Y(68)&
    &+RATE(1365)*D*Y(68)*Y(116)+RATE(1366)*D*Y(68)*Y(133)+RATE(1367)*D*Y(68)&
    &*Y(163)+RATE(1368)*D*Y(68)*Y(333)+RATE(1369)*D*Y(68)*Y(67)+RATE(1370)*D&
    &*Y(68)*Y(116)+RATE(1371)*D*Y(68)*Y(161)+RATE(1372)*D*Y(68)*Y(163)&
    &+RATE(1432)*D*Y(20)*Y(68)+RATE(1443)*D*Y(20)*Y(68)+RATE(1542)*D*Y(24)&
    &*Y(68)+RATE(1680)*D*Y(41)*Y(68)+RATE(1681)*D*Y(41)*Y(68)+RATE(1954)*D&
    &*Y(6)*Y(68)+RATE(2080)*D*Y(62)*Y(68)+RATE(2490)*D*Y(27)*Y(68)+RATE(2624)&
    &*D*Y(35)*Y(68)+RATE(2625)*D*Y(35)*Y(68)+RATE(2697)*D*Y(43)*Y(68)&
    &+RATE(2704)*D*Y(43)*Y(68)+RATE(2827)*D*Y(46)*Y(68)+RATE(2977)*D*Y(56)&
    &*Y(68)
    PROD = RATE(831)*Y(67)+RATE(836)*Y(76)+RATE(1223)*D*Y(16)*Y(21)&
    &+RATE(1311)*D*Y(18)*Y(20)+RATE(1343)*D*Y(18)*Y(16)+RATE(1348)*D*Y(67)&
    &*Y(83)+RATE(1349)*D*Y(67)*Y(100)+RATE(1350)*D*Y(67)*Y(104)+RATE(1351)*D&
    &*Y(67)*Y(162)+RATE(1513)*D*Y(21)*Y(20)+RATE(1872)*D*Y(4)*Y(67)+RATE(1917&
    &)*D*Y(4)*Y(75)+RATE(2008)*D*Y(8)*Y(67)+RATE(2107)*D*Y(63)*Y(67)&
    &+RATE(2320)*D*Y(13)*Y(67)+RATE(2334)*D*Y(13)*Y(80)+RATE(2342)*D*Y(13)&
    &*Y(75)+RATE(2346)*D*Y(13)*Y(284)+RATE(2553)*D*Y(28)*Y(67)+RATE(2919)*D&
    &*Y(48)*Y(67)+RATE(3019)*D*Y(57)*Y(67)
    YDOT(68) = PROD-LOSS
    LOSS = RATE(165)*Y(69)+RATE(558)*D*Y(69)+RATE(916)*Y(69)+RATE(1289)&
    &*D*Y(18)*Y(69)+RATE(1500)*D*Y(21)*Y(69)+RATE(1644)*D*Y(33)*Y(69)&
    &+RATE(1737)*D*Y(53)*Y(69)+RATE(1891)*D*Y(4)*Y(69)+RATE(2113)*D*Y(63)&
    &*Y(69)+RATE(2183)*D*Y(10)*Y(69)+RATE(2475)*D*Y(69)*Y(81)+RATE(2476)*D&
    &*Y(69)*Y(234)+RATE(2477)*D*Y(69)*Y(132)+RATE(2478)*D*Y(69)*Y(184)&
    &+RATE(2479)*D*Y(69)*Y(117)+RATE(2480)*D*Y(69)*Y(174)+RATE(2481)*D*Y(69)&
    &*Y(104)+RATE(2482)*D*Y(69)*Y(134)+RATE(2483)*D*Y(69)*Y(162)+RATE(2484)*D&
    &*Y(69)*Y(165)+RATE(2485)*D*Y(69)*Y(278)+RATE(2486)*D*Y(69)*Y(106)&
    &+RATE(2487)*D*Y(69)*Y(236)+RATE(2564)*D*Y(28)*Y(69)+RATE(2794)*D*Y(55)&
    &*Y(69)
    PROD = RATE(255)*Y(72)/safeMantle+RATE(338)*D*Y(72)/safeMantle*Y(2)&
    &+RATE(421)*Y(72)/safeMantle+RATE(1103)*Y(72)+RATE(1186)*Y(74)+RATE(2488)&
    &*D*Y(70)*Y(333)
    YDOT(69) = PROD-LOSS
    LOSS = RATE(559)*D*Y(70)+RATE(2488)*D*Y(70)*Y(333)
    PROD = RATE(165)*Y(69)+RATE(916)*Y(69)+RATE(1289)*D*Y(18)*Y(69)&
    &+RATE(1500)*D*Y(21)*Y(69)+RATE(1644)*D*Y(33)*Y(69)+RATE(1737)*D*Y(53)&
    &*Y(69)+RATE(1891)*D*Y(4)*Y(69)+RATE(2113)*D*Y(63)*Y(69)+RATE(2183)*D&
    &*Y(10)*Y(69)+RATE(2475)*D*Y(69)*Y(81)+RATE(2476)*D*Y(69)*Y(234)&
    &+RATE(2477)*D*Y(69)*Y(132)+RATE(2478)*D*Y(69)*Y(184)+RATE(2479)*D*Y(69)&
    &*Y(117)+RATE(2480)*D*Y(69)*Y(174)+RATE(2481)*D*Y(69)*Y(104)+RATE(2482)*D&
    &*Y(69)*Y(134)+RATE(2483)*D*Y(69)*Y(162)+RATE(2484)*D*Y(69)*Y(165)&
    &+RATE(2485)*D*Y(69)*Y(278)+RATE(2486)*D*Y(69)*Y(106)+RATE(2487)*D*Y(69)&
    &*Y(236)+RATE(2564)*D*Y(28)*Y(69)+RATE(2794)*D*Y(55)*Y(69)
    YDOT(70) = PROD-LOSS
    LOSS = RATE(204)*Y(71)/safeMantle+RATE(287)*D*Y(71)/safeMantle*Y(2)&
    &+RATE(370)*Y(71)/safeMantle+RATE(969)*Y(71)*totalSwap/safeMantle&
    &+RATE(1052)*Y(71)
    PROD = RATE(2)*Y(73)*bulkLayersReciprocal+RATE(456)*D*Y(67)+RATE(457&
    &)*D*Y(68)+RATE(469)*D*Y(192)
    YDOT(71) = PROD-LOSS
    LOSS = RATE(255)*Y(72)/safeMantle+RATE(338)*D*Y(72)/safeMantle*Y(2)&
    &+RATE(421)*Y(72)/safeMantle+RATE(1020)*Y(72)*totalSwap/safeMantle&
    &+RATE(1103)*Y(72)
    PROD = RATE(53)*Y(74)*bulkLayersReciprocal+RATE(558)*D*Y(69)&
    &+RATE(559)*D*Y(70)
    YDOT(72) = PROD-LOSS
    LOSS = RATE(2)*Y(73)*bulkLayersReciprocal+RATE(1135)*Y(73)
    PROD = RATE(969)*Y(71)*totalSwap/safeMantle
    YDOT(73) = PROD-LOSS
    LOSS = RATE(53)*Y(74)*bulkLayersReciprocal+RATE(1186)*Y(74)
    PROD = RATE(1020)*Y(72)*totalSwap/safeMantle
    YDOT(74) = PROD-LOSS
    LOSS = RATE(112)*Y(75)+RATE(113)*Y(75)+RATE(458)*D*Y(75)+RATE(834)&
    &*Y(75)+RATE(835)*Y(75)+RATE(1303)*D*Y(18)*Y(75)+RATE(1373)*D*Y(75)*Y(83)&
    &+RATE(1374)*D*Y(75)*Y(100)+RATE(1375)*D*Y(75)*Y(104)+RATE(1376)*D*Y(75)&
    &*Y(100)+RATE(1377)*D*Y(75)*Y(132)+RATE(1378)*D*Y(75)*Y(91)+RATE(1379)*D&
    &*Y(75)*Y(117)+RATE(1380)*D*Y(75)*Y(145)+RATE(1381)*D*Y(75)*Y(120)&
    &+RATE(1382)*D*Y(75)*Y(176)+RATE(1383)*D*Y(75)*Y(106)+RATE(1384)*D*Y(75)&
    &*Y(90)+RATE(1385)*D*Y(75)*Y(92)+RATE(1386)*D*Y(75)*Y(291)+RATE(1387)*D&
    &*Y(75)*Y(161)+RATE(1388)*D*Y(75)*Y(82)+RATE(1507)*D*Y(21)*Y(75)&
    &+RATE(1727)*D*Y(53)*Y(75)+RATE(1874)*D*Y(4)*Y(75)+RATE(1917)*D*Y(4)*Y(75&
    &)+RATE(1984)*D*Y(6)*Y(75)+RATE(2010)*D*Y(8)*Y(75)+RATE(2030)*D*Y(8)*Y(75&
    &)+RATE(2109)*D*Y(63)*Y(75)+RATE(2123)*D*Y(63)*Y(75)+RATE(2150)*D*Y(10)&
    &*Y(75)+RATE(2342)*D*Y(13)*Y(75)+RATE(2343)*D*Y(13)*Y(75)+RATE(2344)*D&
    &*Y(13)*Y(75)+RATE(2515)*D*Y(27)*Y(75)+RATE(2554)*D*Y(28)*Y(75)+RATE(2671&
    &)*D*Y(36)*Y(75)+RATE(2736)*D*Y(44)*Y(75)+RATE(2863)*D*Y(46)*Y(75)&
    &+RATE(2921)*D*Y(48)*Y(75)+RATE(2936)*D*Y(48)*Y(75)+RATE(3020)*D*Y(57)&
    &*Y(75)+RATE(3031)*D*Y(57)*Y(75)
    PROD = RATE(115)*Y(80)+RATE(122)*Y(280)+RATE(152)*Y(288)+RATE(205)&
    &*Y(77)/safeMantle+RATE(288)*D*Y(77)/safeMantle*Y(2)+RATE(371)*Y(77&
    &)/safeMantle+RATE(838)*Y(80)+RATE(845)*Y(280)+RATE(899)*Y(288)+RATE(1053&
    &)*Y(77)+RATE(1136)*Y(78)+RATE(1246)*D*Y(16)*Y(24)+RATE(1319)*D*Y(18)&
    &*Y(288)+RATE(1389)*D*Y(76)*Y(133)+RATE(1390)*D*Y(76)*Y(163)+RATE(1407)*D&
    &*Y(81)*Y(333)+RATE(1410)*D*Y(81)*Y(214)+RATE(1411)*D*Y(81)*Y(183)&
    &+RATE(1412)*D*Y(81)*Y(90)+RATE(1413)*D*Y(81)*Y(92)+RATE(1823)*D*Y(2)&
    &*Y(80)+RATE(2081)*D*Y(62)*Y(81)+RATE(2706)*D*Y(43)*Y(81)+RATE(2768)*D&
    &*Y(54)*Y(81)+RATE(2993)*D*Y(56)*Y(80)
    YDOT(75) = PROD-LOSS
    LOSS = RATE(459)*D*Y(76)+RATE(836)*Y(76)+RATE(1222)*D*Y(16)*Y(76)&
    &+RATE(1389)*D*Y(76)*Y(133)+RATE(1390)*D*Y(76)*Y(163)+RATE(1391)*D*Y(76)&
    &*Y(333)+RATE(1392)*D*Y(76)*Y(333)+RATE(1393)*D*Y(76)*Y(131)+RATE(1394)*D&
    &*Y(76)*Y(90)+RATE(1395)*D*Y(76)*Y(90)+RATE(1396)*D*Y(76)*Y(116)&
    &+RATE(1397)*D*Y(76)*Y(92)+RATE(1444)*D*Y(20)*Y(76)+RATE(1552)*D*Y(24)&
    &*Y(76)+RATE(1682)*D*Y(41)*Y(76)+RATE(1955)*D*Y(6)*Y(76)+RATE(2491)*D&
    &*Y(27)*Y(76)+RATE(2492)*D*Y(27)*Y(76)+RATE(2705)*D*Y(43)*Y(76)+RATE(2767&
    &)*D*Y(54)*Y(76)+RATE(2828)*D*Y(46)*Y(76)
    PROD = RATE(113)*Y(75)+RATE(835)*Y(75)+RATE(1224)*D*Y(16)*Y(25)&
    &+RATE(1225)*D*Y(16)*Y(33)+RATE(1304)*D*Y(18)*Y(24)+RATE(1305)*D*Y(18)&
    &*Y(32)+RATE(1352)*D*Y(67)*Y(132)+RATE(1353)*D*Y(67)*Y(91)+RATE(1354)*D&
    &*Y(67)*Y(117)+RATE(1355)*D*Y(67)*Y(145)+RATE(1356)*D*Y(67)*Y(120)&
    &+RATE(1358)*D*Y(67)*Y(176)+RATE(1370)*D*Y(68)*Y(116)+RATE(1373)*D*Y(75)&
    &*Y(83)+RATE(1374)*D*Y(75)*Y(100)+RATE(1375)*D*Y(75)*Y(104)+RATE(1508)*D&
    &*Y(21)*Y(24)+RATE(1680)*D*Y(41)*Y(68)+RATE(1726)*D*Y(53)*Y(67)+RATE(1874&
    &)*D*Y(4)*Y(75)+RATE(1954)*D*Y(6)*Y(68)+RATE(2010)*D*Y(8)*Y(75)+RATE(2028&
    &)*D*Y(8)*Y(67)+RATE(2080)*D*Y(62)*Y(68)+RATE(2109)*D*Y(63)*Y(75)&
    &+RATE(2122)*D*Y(63)*Y(67)+RATE(2149)*D*Y(10)*Y(67)+RATE(2220)*D*Y(66)&
    &*Y(67)+RATE(2335)*D*Y(13)*Y(80)+RATE(2337)*D*Y(13)*Y(89)+RATE(2339)*D&
    &*Y(13)*Y(109)+RATE(2347)*D*Y(13)*Y(280)+RATE(2387)*D*Y(13)*Y(288)&
    &+RATE(2554)*D*Y(28)*Y(75)+RATE(2624)*D*Y(35)*Y(68)+RATE(2668)*D*Y(36)&
    &*Y(67)+RATE(2735)*D*Y(44)*Y(67)+RATE(2921)*D*Y(48)*Y(75)+RATE(3020)*D&
    &*Y(57)*Y(75)+RATE(3030)*D*Y(57)*Y(67)
    YDOT(76) = PROD-LOSS
    LOSS = RATE(205)*Y(77)/safeMantle+RATE(288)*D*Y(77)/safeMantle*Y(2)&
    &+RATE(371)*Y(77)/safeMantle+RATE(970)*Y(77)*totalSwap/safeMantle&
    &+RATE(1053)*Y(77)
    PROD = RATE(3)*Y(78)*bulkLayersReciprocal+RATE(458)*D*Y(75)+RATE(459&
    &)*D*Y(76)
    YDOT(77) = PROD-LOSS
    LOSS = RATE(3)*Y(78)*bulkLayersReciprocal+RATE(1136)*Y(78)
    PROD = RATE(970)*Y(77)*totalSwap/safeMantle
    YDOT(78) = PROD-LOSS
    LOSS = RATE(206)*Y(79)/safeMantle+RATE(289)*D*Y(79)/safeMantle*Y(2)&
    &+RATE(372)*Y(79)/safeMantle+RATE(971)*Y(79)*totalSwap/safeMantle&
    &+RATE(1054)*Y(79)
    PROD = RATE(4)*Y(85)*bulkLayersReciprocal+RATE(460)*D*Y(80)+RATE(461&
    &)*D*Y(81)
    YDOT(79) = PROD-LOSS
    LOSS = RATE(114)*Y(80)+RATE(115)*Y(80)+RATE(460)*D*Y(80)+RATE(837)&
    &*Y(80)+RATE(838)*Y(80)+RATE(1361)*D*Y(67)*Y(80)+RATE(1398)*D*Y(80)*Y(91)&
    &+RATE(1399)*D*Y(80)*Y(290)+RATE(1400)*D*Y(80)*Y(133)+RATE(1401)*D*Y(80)&
    &*Y(105)+RATE(1471)*D*Y(20)*Y(80)+RATE(1704)*D*Y(42)*Y(80)+RATE(1747)*D&
    &*Y(82)*Y(80)+RATE(1823)*D*Y(2)*Y(80)+RATE(1873)*D*Y(4)*Y(80)+RATE(2009)&
    &*D*Y(8)*Y(80)+RATE(2108)*D*Y(63)*Y(80)+RATE(2321)*D*Y(13)*Y(80)&
    &+RATE(2334)*D*Y(13)*Y(80)+RATE(2335)*D*Y(13)*Y(80)+RATE(2336)*D*Y(13)&
    &*Y(80)+RATE(2855)*D*Y(46)*Y(80)+RATE(2920)*D*Y(48)*Y(80)+RATE(2964)*D&
    &*Y(162)*Y(80)+RATE(2968)*D*Y(162)*Y(80)+RATE(2993)*D*Y(56)*Y(80)&
    &+RATE(2994)*D*Y(56)*Y(80)+RATE(2995)*D*Y(56)*Y(80)+RATE(3108)*D*Y(278)&
    &*Y(80)+RATE(3109)*D*Y(278)*Y(80)+RATE(3110)*D*Y(278)*Y(80)
    PROD = RATE(116)*Y(89)+RATE(117)*Y(109)+RATE(206)*Y(79)/safeMantle&
    &+RATE(289)*D*Y(79)/safeMantle*Y(2)+RATE(372)*Y(79)/safeMantle+RATE(839)&
    &*Y(89)+RATE(840)*Y(109)+RATE(1054)*Y(79)+RATE(1137)*Y(85)+RATE(1248)*D&
    &*Y(16)*Y(32)+RATE(1307)*D*Y(18)*Y(206)+RATE(1402)*D*Y(81)*Y(131)&
    &+RATE(1403)*D*Y(81)*Y(183)+RATE(1404)*D*Y(81)*Y(116)+RATE(1405)*D*Y(81)&
    &*Y(133)+RATE(1420)*D*Y(89)*Y(161)+RATE(1572)*D*Y(24)*Y(24)+RATE(1573)*D&
    &*Y(24)*Y(24)+RATE(1618)*D*Y(32)*Y(89)+RATE(1824)*D*Y(2)*Y(89)+RATE(1984)&
    &*D*Y(6)*Y(75)+RATE(2475)*D*Y(69)*Y(81)+RATE(2511)*D*Y(27)*Y(89)&
    &+RATE(2755)*D*Y(54)*Y(81)+RATE(2996)*D*Y(56)*Y(89)
    YDOT(80) = PROD-LOSS
    LOSS = RATE(461)*D*Y(81)+RATE(1402)*D*Y(81)*Y(131)+RATE(1403)*D*Y(81&
    &)*Y(183)+RATE(1404)*D*Y(81)*Y(116)+RATE(1405)*D*Y(81)*Y(133)+RATE(1406)&
    &*D*Y(81)*Y(333)+RATE(1407)*D*Y(81)*Y(333)+RATE(1408)*D*Y(81)*Y(333)&
    &+RATE(1409)*D*Y(81)*Y(214)+RATE(1410)*D*Y(81)*Y(214)+RATE(1411)*D*Y(81)&
    &*Y(183)+RATE(1412)*D*Y(81)*Y(90)+RATE(1413)*D*Y(81)*Y(92)+RATE(1414)*D&
    &*Y(81)*Y(105)+RATE(1415)*D*Y(81)*Y(166)+RATE(1416)*D*Y(81)*Y(166)&
    &+RATE(1417)*D*Y(81)*Y(166)+RATE(1418)*D*Y(81)*Y(166)+RATE(1683)*D*Y(41)&
    &*Y(81)+RATE(2081)*D*Y(62)*Y(81)+RATE(2475)*D*Y(69)*Y(81)+RATE(2493)*D&
    &*Y(27)*Y(81)+RATE(2494)*D*Y(27)*Y(81)+RATE(2495)*D*Y(27)*Y(81)+RATE(2706&
    &)*D*Y(43)*Y(81)+RATE(2755)*D*Y(54)*Y(81)+RATE(2768)*D*Y(54)*Y(81)&
    &+RATE(2829)*D*Y(46)*Y(81)
    PROD = RATE(114)*Y(80)+RATE(837)*Y(80)+RATE(1306)*D*Y(18)*Y(32)&
    &+RATE(1307)*D*Y(18)*Y(206)+RATE(1310)*D*Y(18)*Y(41)+RATE(1377)*D*Y(75)&
    &*Y(132)+RATE(1378)*D*Y(75)*Y(91)+RATE(1379)*D*Y(75)*Y(117)+RATE(1380)*D&
    &*Y(75)*Y(145)+RATE(1381)*D*Y(75)*Y(120)+RATE(1382)*D*Y(75)*Y(176)&
    &+RATE(1394)*D*Y(76)*Y(90)+RATE(1396)*D*Y(76)*Y(116)+RATE(1398)*D*Y(80)&
    &*Y(91)+RATE(1399)*D*Y(80)*Y(290)+RATE(1445)*D*Y(20)*Y(33)+RATE(1512)*D&
    &*Y(21)*Y(41)+RATE(1681)*D*Y(41)*Y(68)+RATE(1682)*D*Y(41)*Y(76)+RATE(1704&
    &)*D*Y(42)*Y(80)+RATE(1727)*D*Y(53)*Y(75)+RATE(1873)*D*Y(4)*Y(80)&
    &+RATE(1915)*D*Y(4)*Y(89)+RATE(1916)*D*Y(4)*Y(109)+RATE(1955)*D*Y(6)*Y(76&
    &)+RATE(2009)*D*Y(8)*Y(80)+RATE(2029)*D*Y(8)*Y(109)+RATE(2030)*D*Y(8)&
    &*Y(75)+RATE(2108)*D*Y(63)*Y(80)+RATE(2123)*D*Y(63)*Y(75)+RATE(2150)*D&
    &*Y(10)*Y(75)+RATE(2321)*D*Y(13)*Y(80)+RATE(2338)*D*Y(13)*Y(89)+RATE(2340&
    &)*D*Y(13)*Y(109)+RATE(2671)*D*Y(36)*Y(75)+RATE(2736)*D*Y(44)*Y(75)&
    &+RATE(2799)*D*Y(55)*Y(67)+RATE(2920)*D*Y(48)*Y(80)+RATE(2935)*D*Y(48)&
    &*Y(109)+RATE(2964)*D*Y(162)*Y(80)+RATE(3031)*D*Y(57)*Y(75)
    YDOT(81) = PROD-LOSS
    LOSS = RATE(138)*Y(82)+RATE(496)*D*Y(82)+RATE(873)*Y(82)+RATE(1250)&
    &*D*Y(16)*Y(82)+RATE(1388)*D*Y(75)*Y(82)+RATE(1514)*D*Y(21)*Y(82)&
    &+RATE(1577)*D*Y(24)*Y(82)+RATE(1622)*D*Y(32)*Y(82)+RATE(1642)*D*Y(32)&
    &*Y(82)+RATE(1700)*D*Y(41)*Y(82)+RATE(1744)*D*Y(82)*Y(104)+RATE(1745)*D&
    &*Y(82)*Y(145)+RATE(1746)*D*Y(82)*Y(176)+RATE(1747)*D*Y(82)*Y(80)&
    &+RATE(1748)*D*Y(82)*Y(109)+RATE(1749)*D*Y(82)*Y(82)+RATE(1749)*D*Y(82)&
    &*Y(82)+RATE(1750)*D*Y(82)*Y(131)+RATE(1751)*D*Y(82)*Y(90)+RATE(1752)*D&
    &*Y(82)*Y(116)+RATE(1753)*D*Y(82)*Y(92)+RATE(1754)*D*Y(82)*Y(144)&
    &+RATE(1755)*D*Y(82)*Y(264)+RATE(1756)*D*Y(82)*Y(133)+RATE(1757)*D*Y(82)&
    &*Y(133)+RATE(1758)*D*Y(82)*Y(161)+RATE(1759)*D*Y(82)*Y(161)+RATE(1760)*D&
    &*Y(82)*Y(163)+RATE(1761)*D*Y(82)*Y(166)+RATE(1989)*D*Y(6)*Y(82)&
    &+RATE(2014)*D*Y(8)*Y(82)+RATE(2036)*D*Y(8)*Y(82)+RATE(2163)*D*Y(10)*Y(82&
    &)+RATE(2364)*D*Y(13)*Y(82)+RATE(2365)*D*Y(13)*Y(82)+RATE(2527)*D*Y(27)&
    &*Y(82)+RATE(2557)*D*Y(28)*Y(82)+RATE(2643)*D*Y(35)*Y(82)+RATE(2672)*D&
    &*Y(36)*Y(82)+RATE(2792)*D*Y(54)*Y(82)+RATE(2868)*D*Y(46)*Y(82)+RATE(2869&
    &)*D*Y(46)*Y(82)+RATE(2940)*D*Y(48)*Y(82)+RATE(2998)*D*Y(56)*Y(82)&
    &+RATE(2999)*D*Y(56)*Y(82)+RATE(3032)*D*Y(57)*Y(82)
    PROD = RATE(120)*Y(197)+RATE(121)*Y(284)+RATE(133)*Y(214)+RATE(152)&
    &*Y(288)+RATE(154)*Y(90)+RATE(160)*Y(92)+RATE(168)*Y(291)+RATE(168)*Y(291&
    &)+RATE(184)*Y(221)+RATE(227)*Y(84)/safeMantle+RATE(310)*D*Y(84&
    &)/safeMantle*Y(2)+RATE(393)*Y(84)/safeMantle+RATE(725)*Y(17)*Y(222)&
    &+RATE(777)*Y(19)*Y(225)*bulkLayersReciprocal+RATE(843)*Y(197)+RATE(844)&
    &*Y(284)+RATE(862)*Y(214)+RATE(899)*Y(288)+RATE(902)*Y(90)+RATE(908)*Y(92&
    &)+RATE(919)*Y(291)+RATE(919)*Y(291)+RATE(937)*Y(221)+RATE(1075)*Y(84)&
    &+RATE(1158)*Y(86)+RATE(1218)*D*Y(16)*Y(83)+RATE(1230)*D*Y(16)*Y(91)&
    &+RATE(1244)*D*Y(16)*Y(197)+RATE(1257)*D*Y(16)*Y(103)+RATE(1258)*D*Y(16)&
    &*Y(291)+RATE(1262)*D*Y(16)*Y(35)+RATE(1264)*D*Y(16)*Y(133)+RATE(1267)*D&
    &*Y(16)*Y(265)+RATE(1269)*D*Y(16)*Y(221)+RATE(1279)*D*Y(16)*Y(27)&
    &+RATE(1290)*D*Y(18)*Y(291)+RATE(1331)*D*Y(18)*Y(221)+RATE(1348)*D*Y(67)&
    &*Y(83)+RATE(1353)*D*Y(67)*Y(91)+RATE(1373)*D*Y(75)*Y(83)+RATE(1378)*D&
    &*Y(75)*Y(91)+RATE(1386)*D*Y(75)*Y(291)+RATE(1394)*D*Y(76)*Y(90)&
    &+RATE(1409)*D*Y(81)*Y(214)+RATE(1423)*D*Y(199)*Y(333)+RATE(1425)*D*Y(290&
    &)*Y(333)+RATE(1425)*D*Y(290)*Y(333)+RATE(1433)*D*Y(20)*Y(83)+RATE(1452)&
    &*D*Y(20)*Y(91)+RATE(1480)*D*Y(20)*Y(27)+RATE(1543)*D*Y(24)*Y(83)&
    &+RATE(1558)*D*Y(24)*Y(91)+RATE(1762)*D*Y(83)*Y(99)+RATE(1763)*D*Y(83)&
    &*Y(131)+RATE(1764)*D*Y(83)*Y(90)+RATE(1765)*D*Y(83)*Y(116)+RATE(1766)*D&
    &*Y(83)*Y(133)+RATE(1767)*D*Y(83)*Y(161)+RATE(1768)*D*Y(83)*Y(163)&
    &+RATE(1803)*D*Y(2)*Y(83)+RATE(1809)*D*Y(2)*Y(290)+RATE(1836)*D*Y(2)*Y(90&
    &)+RATE(1845)*D*Y(2)*Y(291)+RATE(1860)*D*Y(2)*Y(221)+RATE(2092)*D*Y(62)&
    &*Y(91)+RATE(2258)*D*Y(91)*Y(333)+RATE(2259)*D*Y(91)*Y(232)+RATE(2260)*D&
    &*Y(91)*Y(99)+RATE(2261)*D*Y(91)*Y(131)+RATE(2262)*D*Y(91)*Y(90)&
    &+RATE(2263)*D*Y(91)*Y(116)+RATE(2265)*D*Y(91)*Y(92)+RATE(2266)*D*Y(91)&
    &*Y(163)+RATE(2267)*D*Y(102)*Y(333)+RATE(2345)*D*Y(13)*Y(197)+RATE(2346)&
    &*D*Y(13)*Y(284)+RATE(2356)*D*Y(13)*Y(214)+RATE(2387)*D*Y(13)*Y(288)&
    &+RATE(2407)*D*Y(13)*Y(291)+RATE(2419)*D*Y(13)*Y(221)+RATE(2490)*D*Y(27)&
    &*Y(68)+RATE(2492)*D*Y(27)*Y(76)+RATE(2507)*D*Y(27)*Y(208)+RATE(2510)*D&
    &*Y(27)*Y(67)+RATE(2516)*D*Y(27)*Y(197)+RATE(2516)*D*Y(27)*Y(197)&
    &+RATE(2518)*D*Y(27)*Y(284)+RATE(2520)*D*Y(27)*Y(308)+RATE(2529)*D*Y(27)&
    &*Y(233)+RATE(2552)*D*Y(27)*Y(207)+RATE(2620)*D*Y(35)*Y(83)+RATE(2631)*D&
    &*Y(35)*Y(91)+RATE(2698)*D*Y(43)*Y(83)+RATE(2713)*D*Y(43)*Y(91)+RATE(2824&
    &)*D*Y(46)*Y(83)+RATE(2864)*D*Y(46)*Y(197)+RATE(2877)*D*Y(46)*Y(90)&
    &+RATE(2899)*D*Y(46)*Y(221)+RATE(2978)*D*Y(56)*Y(83)+RATE(2984)*D*Y(56)&
    &*Y(91)+RATE(3006)*D*Y(56)*Y(90)
    YDOT(82) = PROD-LOSS
    LOSS = RATE(497)*D*Y(83)+RATE(1218)*D*Y(16)*Y(83)+RATE(1348)*D*Y(67)&
    &*Y(83)+RATE(1373)*D*Y(75)*Y(83)+RATE(1433)*D*Y(20)*Y(83)+RATE(1543)*D&
    &*Y(24)*Y(83)+RATE(1762)*D*Y(83)*Y(99)+RATE(1763)*D*Y(83)*Y(131)&
    &+RATE(1764)*D*Y(83)*Y(90)+RATE(1765)*D*Y(83)*Y(116)+RATE(1766)*D*Y(83)&
    &*Y(133)+RATE(1767)*D*Y(83)*Y(161)+RATE(1768)*D*Y(83)*Y(163)+RATE(1769)*D&
    &*Y(83)*Y(333)+RATE(1770)*D*Y(83)*Y(131)+RATE(1771)*D*Y(83)*Y(90)&
    &+RATE(1772)*D*Y(83)*Y(116)+RATE(1773)*D*Y(83)*Y(161)+RATE(1803)*D*Y(2)&
    &*Y(83)+RATE(1959)*D*Y(6)*Y(83)+RATE(2085)*D*Y(62)*Y(83)+RATE(2086)*D&
    &*Y(62)*Y(83)+RATE(2497)*D*Y(27)*Y(83)+RATE(2620)*D*Y(35)*Y(83)+RATE(2698&
    &)*D*Y(43)*Y(83)+RATE(2824)*D*Y(46)*Y(83)+RATE(2978)*D*Y(56)*Y(83)
    PROD = RATE(1327)*D*Y(18)*Y(35)+RATE(1344)*D*Y(18)*Y(27)+RATE(1458)&
    &*D*Y(20)*Y(28)+RATE(1529)*D*Y(21)*Y(27)+RATE(1532)*D*Y(21)*Y(35)&
    &+RATE(1744)*D*Y(82)*Y(104)+RATE(2014)*D*Y(8)*Y(82)+RATE(2355)*D*Y(13)&
    &*Y(214)+RATE(2388)*D*Y(13)*Y(90)+RATE(2398)*D*Y(13)*Y(92)+RATE(2407)*D&
    &*Y(13)*Y(291)+RATE(2418)*D*Y(13)*Y(221)+RATE(2557)*D*Y(28)*Y(82)
    YDOT(83) = PROD-LOSS
    LOSS = RATE(227)*Y(84)/safeMantle+RATE(310)*D*Y(84)/safeMantle*Y(2)&
    &+RATE(393)*Y(84)/safeMantle+RATE(625)*Y(84)*Y(3)+RATE(729)*Y(84)*Y(3)&
    &+RATE(992)*Y(84)*totalSwap/safeMantle+RATE(1075)*Y(84)
    PROD = RATE(25)*Y(86)*bulkLayersReciprocal+RATE(496)*D*Y(82)&
    &+RATE(497)*D*Y(83)+RATE(621)*Y(17)*Y(222)
    YDOT(84) = PROD-LOSS
    LOSS = RATE(4)*Y(85)*bulkLayersReciprocal+RATE(1137)*Y(85)
    PROD = RATE(971)*Y(79)*totalSwap/safeMantle
    YDOT(85) = PROD-LOSS
    LOSS = RATE(25)*Y(86)*bulkLayersReciprocal+RATE(677)*Y(86)*Y(5)&
    &*bulkLayersReciprocal+RATE(781)*Y(86)*Y(5)*bulkLayersReciprocal&
    &+RATE(1158)*Y(86)
    PROD = RATE(673)*Y(19)*Y(225)*bulkLayersReciprocal+RATE(992)*Y(84)&
    &*totalSwap/safeMantle
    YDOT(86) = PROD-LOSS
    LOSS = RATE(243)*Y(87)/safeMantle+RATE(326)*D*Y(87)/safeMantle*Y(2)&
    &+RATE(409)*Y(87)/safeMantle+RATE(1008)*Y(87)*totalSwap/safeMantle&
    &+RATE(1091)*Y(87)
    PROD = RATE(41)*Y(94)*bulkLayersReciprocal+RATE(531)*D*Y(90)&
    &+RATE(532)*D*Y(91)+RATE(533)*D*Y(102)+RATE(625)*Y(84)*Y(3)
    YDOT(87) = PROD-LOSS
    LOSS = RATE(249)*Y(88)/safeMantle+RATE(332)*D*Y(88)/safeMantle*Y(2)&
    &+RATE(415)*Y(88)/safeMantle+RATE(1014)*Y(88)*totalSwap/safeMantle&
    &+RATE(1097)*Y(88)
    PROD = RATE(47)*Y(95)*bulkLayersReciprocal+RATE(544)*D*Y(92)
    YDOT(88) = PROD-LOSS
    LOSS = RATE(116)*Y(89)+RATE(462)*D*Y(89)+RATE(839)*Y(89)+RATE(1242)&
    &*D*Y(16)*Y(89)+RATE(1419)*D*Y(89)*Y(161)+RATE(1420)*D*Y(89)*Y(161)&
    &+RATE(1618)*D*Y(32)*Y(89)+RATE(1824)*D*Y(2)*Y(89)+RATE(1915)*D*Y(4)*Y(89&
    &)+RATE(2337)*D*Y(13)*Y(89)+RATE(2338)*D*Y(13)*Y(89)+RATE(2511)*D*Y(27)&
    &*Y(89)+RATE(2856)*D*Y(46)*Y(89)+RATE(2996)*D*Y(56)*Y(89)+RATE(3018)*D&
    &*Y(56)*Y(89)
    PROD = RATE(118)*Y(123)+RATE(207)*Y(93)/safeMantle+RATE(290)*D*Y(93&
    &)/safeMantle*Y(2)+RATE(373)*Y(93)/safeMantle+RATE(841)*Y(123)+RATE(1055)&
    &*Y(93)+RATE(1138)*Y(96)+RATE(1418)*D*Y(81)*Y(166)+RATE(1574)*D*Y(24)&
    &*Y(24)+RATE(1748)*D*Y(82)*Y(109)+RATE(2857)*D*Y(46)*Y(109)
    YDOT(89) = PROD-LOSS
    LOSS = RATE(154)*Y(90)+RATE(531)*D*Y(90)+RATE(902)*Y(90)+RATE(1362)&
    &*D*Y(67)*Y(90)+RATE(1384)*D*Y(75)*Y(90)+RATE(1394)*D*Y(76)*Y(90)&
    &+RATE(1395)*D*Y(76)*Y(90)+RATE(1412)*D*Y(81)*Y(90)+RATE(1524)*D*Y(21)&
    &*Y(90)+RATE(1525)*D*Y(21)*Y(90)+RATE(1526)*D*Y(21)*Y(90)+RATE(1667)*D&
    &*Y(33)*Y(90)+RATE(1733)*D*Y(53)*Y(90)+RATE(1751)*D*Y(82)*Y(90)+RATE(1764&
    &)*D*Y(83)*Y(90)+RATE(1771)*D*Y(83)*Y(90)+RATE(1836)*D*Y(2)*Y(90)&
    &+RATE(1886)*D*Y(4)*Y(90)+RATE(2019)*D*Y(8)*Y(90)+RATE(2129)*D*Y(63)*Y(90&
    &)+RATE(2173)*D*Y(10)*Y(90)+RATE(2227)*D*Y(66)*Y(90)+RATE(2241)*D*Y(90)&
    &*Y(100)+RATE(2242)*D*Y(90)*Y(104)+RATE(2243)*D*Y(90)*Y(290)+RATE(2244)*D&
    &*Y(90)*Y(192)+RATE(2245)*D*Y(90)*Y(172)+RATE(2246)*D*Y(90)*Y(132)&
    &+RATE(2247)*D*Y(90)*Y(142)+RATE(2248)*D*Y(90)*Y(189)+RATE(2249)*D*Y(90)&
    &*Y(117)+RATE(2250)*D*Y(90)*Y(145)+RATE(2251)*D*Y(90)*Y(174)+RATE(2252)*D&
    &*Y(90)*Y(307)+RATE(2253)*D*Y(90)*Y(120)+RATE(2254)*D*Y(90)*Y(176)&
    &+RATE(2262)*D*Y(91)*Y(90)+RATE(2388)*D*Y(13)*Y(90)+RATE(2389)*D*Y(13)&
    &*Y(90)+RATE(2390)*D*Y(13)*Y(90)+RATE(2391)*D*Y(13)*Y(90)+RATE(2562)*D&
    &*Y(28)*Y(90)+RATE(2683)*D*Y(36)*Y(90)+RATE(2746)*D*Y(44)*Y(90)+RATE(2877&
    &)*D*Y(46)*Y(90)+RATE(2878)*D*Y(46)*Y(90)+RATE(2879)*D*Y(46)*Y(90)&
    &+RATE(2945)*D*Y(48)*Y(90)+RATE(2946)*D*Y(48)*Y(90)+RATE(3006)*D*Y(56)&
    &*Y(90)+RATE(3007)*D*Y(56)*Y(90)+RATE(3038)*D*Y(57)*Y(90)
    PROD = RATE(144)*Y(101)+RATE(243)*Y(87)/safeMantle+RATE(326)*D*Y(87&
    &)/safeMantle*Y(2)+RATE(409)*Y(87)/safeMantle+RATE(729)*Y(84)*Y(3)&
    &+RATE(781)*Y(86)*Y(5)*bulkLayersReciprocal+RATE(882)*Y(101)+RATE(1091)&
    &*Y(87)+RATE(1174)*Y(94)+RATE(1259)*D*Y(16)*Y(43)+RATE(1320)*D*Y(18)&
    &*Y(288)+RATE(1398)*D*Y(80)*Y(91)+RATE(1400)*D*Y(80)*Y(133)+RATE(1453)*D&
    &*Y(20)*Y(102)+RATE(1479)*D*Y(20)*Y(103)+RATE(1482)*D*Y(20)*Y(133)&
    &+RATE(1559)*D*Y(24)*Y(102)+RATE(1577)*D*Y(24)*Y(82)+RATE(1581)*D*Y(24)&
    &*Y(103)+RATE(1584)*D*Y(24)*Y(133)+RATE(1622)*D*Y(32)*Y(82)+RATE(1631)*D&
    &*Y(32)*Y(133)+RATE(1700)*D*Y(41)*Y(82)+RATE(1748)*D*Y(82)*Y(109)&
    &+RATE(1750)*D*Y(82)*Y(131)+RATE(1752)*D*Y(82)*Y(116)+RATE(1754)*D*Y(82)&
    &*Y(144)+RATE(1761)*D*Y(82)*Y(166)+RATE(1770)*D*Y(83)*Y(131)+RATE(1806)*D&
    &*Y(2)*Y(91)+RATE(1832)*D*Y(2)*Y(101)+RATE(1840)*D*Y(2)*Y(92)+RATE(1845)&
    &*D*Y(2)*Y(291)+RATE(1858)*D*Y(2)*Y(221)+RATE(1871)*D*Y(4)*Y(92)&
    &+RATE(1919)*D*Y(4)*Y(214)+RATE(1989)*D*Y(6)*Y(82)+RATE(2078)*D*Y(62)&
    &*Y(91)+RATE(2083)*D*Y(62)*Y(199)+RATE(2137)*D*Y(183)*Y(199)+RATE(2255)*D&
    &*Y(91)*Y(133)+RATE(2256)*D*Y(91)*Y(161)+RATE(2257)*D*Y(91)*Y(163)&
    &+RATE(2268)*D*Y(102)*Y(333)+RATE(2270)*D*Y(102)*Y(214)+RATE(2272)*D&
    &*Y(102)*Y(131)+RATE(2274)*D*Y(102)*Y(183)+RATE(2495)*D*Y(27)*Y(81)&
    &+RATE(2512)*D*Y(27)*Y(109)+RATE(2521)*D*Y(27)*Y(24)+RATE(2525)*D*Y(27)&
    &*Y(32)+RATE(2526)*D*Y(27)*Y(32)+RATE(2530)*D*Y(27)*Y(101)+RATE(2532)*D&
    &*Y(27)*Y(116)+RATE(2534)*D*Y(27)*Y(244)+RATE(2643)*D*Y(35)*Y(82)&
    &+RATE(2714)*D*Y(43)*Y(102)+RATE(2760)*D*Y(54)*Y(91)+RATE(2769)*D*Y(54)&
    &*Y(199)+RATE(2779)*D*Y(54)*Y(102)+RATE(2792)*D*Y(54)*Y(82)+RATE(2998)*D&
    &*Y(56)*Y(82)+RATE(3010)*D*Y(56)*Y(291)
    YDOT(90) = PROD-LOSS
    LOSS = RATE(532)*D*Y(91)+RATE(1230)*D*Y(16)*Y(91)+RATE(1353)*D*Y(67)&
    &*Y(91)+RATE(1378)*D*Y(75)*Y(91)+RATE(1398)*D*Y(80)*Y(91)+RATE(1452)*D&
    &*Y(20)*Y(91)+RATE(1558)*D*Y(24)*Y(91)+RATE(1688)*D*Y(41)*Y(91)+RATE(1806&
    &)*D*Y(2)*Y(91)+RATE(1966)*D*Y(6)*Y(91)+RATE(2078)*D*Y(62)*Y(91)&
    &+RATE(2092)*D*Y(62)*Y(91)+RATE(2255)*D*Y(91)*Y(133)+RATE(2256)*D*Y(91)&
    &*Y(161)+RATE(2257)*D*Y(91)*Y(163)+RATE(2258)*D*Y(91)*Y(333)+RATE(2259)*D&
    &*Y(91)*Y(232)+RATE(2260)*D*Y(91)*Y(99)+RATE(2261)*D*Y(91)*Y(131)&
    &+RATE(2262)*D*Y(91)*Y(90)+RATE(2263)*D*Y(91)*Y(116)+RATE(2264)*D*Y(91)&
    &*Y(116)+RATE(2265)*D*Y(91)*Y(92)+RATE(2266)*D*Y(91)*Y(163)+RATE(2631)*D&
    &*Y(35)*Y(91)+RATE(2713)*D*Y(43)*Y(91)+RATE(2760)*D*Y(54)*Y(91)+RATE(2778&
    &)*D*Y(54)*Y(91)+RATE(2984)*D*Y(56)*Y(91)
    PROD = RATE(1325)*D*Y(18)*Y(43)+RATE(1326)*D*Y(18)*Y(54)+RATE(1530)&
    &*D*Y(21)*Y(43)+RATE(1745)*D*Y(82)*Y(145)+RATE(1746)*D*Y(82)*Y(176)&
    &+RATE(1764)*D*Y(83)*Y(90)+RATE(1772)*D*Y(83)*Y(116)+RATE(1809)*D*Y(2)&
    &*Y(290)+RATE(1886)*D*Y(4)*Y(90)+RATE(1959)*D*Y(6)*Y(83)+RATE(2019)*D*Y(8&
    &)*Y(90)+RATE(2036)*D*Y(8)*Y(82)+RATE(2085)*D*Y(62)*Y(83)+RATE(2163)*D&
    &*Y(10)*Y(82)+RATE(2241)*D*Y(90)*Y(100)+RATE(2242)*D*Y(90)*Y(104)&
    &+RATE(2243)*D*Y(90)*Y(290)+RATE(2496)*D*Y(27)*Y(25)+RATE(2562)*D*Y(28)&
    &*Y(90)+RATE(2577)*D*Y(28)*Y(41)+RATE(2670)*D*Y(36)*Y(67)+RATE(2672)*D&
    &*Y(36)*Y(82)+RATE(3032)*D*Y(57)*Y(82)
    YDOT(91) = PROD-LOSS
    LOSS = RATE(160)*Y(92)+RATE(544)*D*Y(92)+RATE(908)*Y(92)+RATE(1323)&
    &*D*Y(18)*Y(92)+RATE(1385)*D*Y(75)*Y(92)+RATE(1397)*D*Y(76)*Y(92)&
    &+RATE(1413)*D*Y(81)*Y(92)+RATE(1528)*D*Y(21)*Y(92)+RATE(1736)*D*Y(53)&
    &*Y(92)+RATE(1753)*D*Y(82)*Y(92)+RATE(1840)*D*Y(2)*Y(92)+RATE(1871)*D*Y(4&
    &)*Y(92)+RATE(2132)*D*Y(63)*Y(92)+RATE(2179)*D*Y(10)*Y(92)+RATE(2228)*D&
    &*Y(66)*Y(92)+RATE(2265)*D*Y(91)*Y(92)+RATE(2398)*D*Y(13)*Y(92)+RATE(2399&
    &)*D*Y(13)*Y(92)+RATE(2400)*D*Y(13)*Y(92)+RATE(2447)*D*Y(92)*Y(132)&
    &+RATE(2448)*D*Y(92)*Y(142)+RATE(2449)*D*Y(92)*Y(189)+RATE(2450)*D*Y(92)&
    &*Y(117)+RATE(2451)*D*Y(92)*Y(145)+RATE(2452)*D*Y(92)*Y(174)+RATE(2453)*D&
    &*Y(92)*Y(120)+RATE(2454)*D*Y(92)*Y(176)+RATE(2685)*D*Y(36)*Y(92)&
    &+RATE(2748)*D*Y(44)*Y(92)+RATE(3041)*D*Y(57)*Y(92)
    PROD = RATE(249)*Y(88)/safeMantle+RATE(332)*D*Y(88)/safeMantle*Y(2)&
    &+RATE(415)*Y(88)/safeMantle+RATE(1097)*Y(88)+RATE(1180)*Y(95)+RATE(1260)&
    &*D*Y(16)*Y(43)+RATE(1454)*D*Y(20)*Y(102)+RATE(1560)*D*Y(24)*Y(102)&
    &+RATE(1671)*D*Y(220)*Y(333)+RATE(2269)*D*Y(102)*Y(333)+RATE(2271)*D&
    &*Y(102)*Y(214)+RATE(2273)*D*Y(102)*Y(131)+RATE(2275)*D*Y(102)*Y(183)&
    &+RATE(2455)*D*Y(227)*Y(16)+RATE(2522)*D*Y(27)*Y(24)+RATE(2715)*D*Y(43)&
    &*Y(102)+RATE(2780)*D*Y(54)*Y(102)
    YDOT(92) = PROD-LOSS
    LOSS = RATE(207)*Y(93)/safeMantle+RATE(290)*D*Y(93)/safeMantle*Y(2)&
    &+RATE(373)*Y(93)/safeMantle+RATE(972)*Y(93)*totalSwap/safeMantle&
    &+RATE(1055)*Y(93)
    PROD = RATE(5)*Y(96)*bulkLayersReciprocal+RATE(462)*D*Y(89)
    YDOT(93) = PROD-LOSS
    LOSS = RATE(41)*Y(94)*bulkLayersReciprocal+RATE(1174)*Y(94)
    PROD = RATE(677)*Y(86)*Y(5)*bulkLayersReciprocal+RATE(1008)*Y(87)&
    &*totalSwap/safeMantle
    YDOT(94) = PROD-LOSS
    LOSS = RATE(47)*Y(95)*bulkLayersReciprocal+RATE(1180)*Y(95)
    PROD = RATE(1014)*Y(88)*totalSwap/safeMantle
    YDOT(95) = PROD-LOSS
    LOSS = RATE(5)*Y(96)*bulkLayersReciprocal+RATE(1138)*Y(96)
    PROD = RATE(972)*Y(93)*totalSwap/safeMantle
    YDOT(96) = PROD-LOSS
    LOSS = RATE(228)*Y(97)/safeMantle+RATE(311)*D*Y(97)/safeMantle*Y(2)&
    &+RATE(394)*Y(97)/safeMantle+RATE(626)*Y(97)*Y(58)+RATE(638)*Y(3)*Y(97)&
    &+RATE(656)*Y(37)*Y(97)+RATE(730)*Y(97)*Y(58)+RATE(742)*Y(3)*Y(97)&
    &+RATE(760)*Y(37)*Y(97)+RATE(993)*Y(97)*totalSwap/safeMantle+RATE(1076)&
    &*Y(97)
    PROD = RATE(26)*Y(111)*bulkLayersReciprocal+RATE(90)*Y(130)+RATE(93)&
    &*Y(118)+RATE(109)*Y(226)+RATE(498)*D*Y(99)+RATE(499)*D*Y(100)+RATE(621)&
    &*Y(17)*Y(222)+RATE(640)*Y(3)*Y(118)+RATE(647)*Y(3)*Y(226)+RATE(660)*Y(45&
    &)*Y(118)+RATE(664)*Y(47)*Y(222)+RATE(829)*Y(226)
    YDOT(97) = PROD-LOSS
    LOSS = RATE(257)*Y(98)/safeMantle+RATE(340)*D*Y(98)/safeMantle*Y(2)&
    &+RATE(423)*Y(98)/safeMantle+RATE(1022)*Y(98)*totalSwap/safeMantle&
    &+RATE(1105)*Y(98)
    PROD = RATE(55)*Y(112)*bulkLayersReciprocal+RATE(562)*D*Y(103)&
    &+RATE(563)*D*Y(104)+RATE(564)*D*Y(120)+RATE(655)*Y(29)*Y(29)
    YDOT(98) = PROD-LOSS
    LOSS = RATE(101)*Y(99)+RATE(139)*Y(99)+RATE(498)*D*Y(99)+RATE(874)&
    &*Y(99)+RATE(1251)*D*Y(16)*Y(99)+RATE(1715)*D*Y(42)*Y(99)+RATE(1729)*D&
    &*Y(53)*Y(99)+RATE(1762)*D*Y(83)*Y(99)+RATE(1774)*D*Y(99)*Y(104)&
    &+RATE(1775)*D*Y(99)*Y(196)+RATE(1776)*D*Y(99)*Y(243)+RATE(1777)*D*Y(99)&
    &*Y(145)+RATE(1778)*D*Y(99)*Y(120)+RATE(1779)*D*Y(99)*Y(176)+RATE(1780)*D&
    &*Y(99)*Y(319)+RATE(1781)*D*Y(99)*Y(167)+RATE(1782)*D*Y(99)*Y(236)&
    &+RATE(1783)*D*Y(99)*Y(144)+RATE(1784)*D*Y(99)*Y(264)+RATE(1785)*D*Y(99)&
    &*Y(161)+RATE(1786)*D*Y(99)*Y(175)+RATE(1831)*D*Y(2)*Y(99)+RATE(2015)*D&
    &*Y(8)*Y(99)+RATE(2038)*D*Y(8)*Y(99)+RATE(2124)*D*Y(63)*Y(99)+RATE(2165)&
    &*D*Y(10)*Y(99)+RATE(2166)*D*Y(10)*Y(99)+RATE(2260)*D*Y(91)*Y(99)&
    &+RATE(2370)*D*Y(13)*Y(99)+RATE(2558)*D*Y(28)*Y(99)+RATE(2580)*D*Y(28)&
    &*Y(99)+RATE(2676)*D*Y(36)*Y(99)+RATE(2923)*D*Y(48)*Y(99)+RATE(3000)*D&
    &*Y(56)*Y(99)+RATE(3034)*D*Y(57)*Y(99)+RATE(3078)*D*Y(105)*Y(99)
    PROD = RATE(127)*Y(218)+RATE(131)*Y(237)+RATE(140)*Y(232)+RATE(145)&
    &*Y(131)+RATE(155)*Y(116)+RATE(161)*Y(227)+RATE(186)*Y(300)+RATE(228)&
    &*Y(97)/safeMantle+RATE(311)*D*Y(97)/safeMantle*Y(2)+RATE(394)*Y(97&
    &)/safeMantle+RATE(725)*Y(17)*Y(222)+RATE(744)*Y(3)*Y(118)+RATE(751)*Y(3)&
    &*Y(226)+RATE(764)*Y(45)*Y(118)+RATE(768)*Y(47)*Y(222)+RATE(777)*Y(19)&
    &*Y(225)*bulkLayersReciprocal+RATE(796)*Y(5)*Y(126)*bulkLayersReciprocal&
    &+RATE(803)*Y(5)*Y(228)*bulkLayersReciprocal+RATE(816)*Y(50)*Y(126)&
    &*bulkLayersReciprocal+RATE(820)*Y(51)*Y(225)*bulkLayersReciprocal&
    &+RATE(854)*Y(218)+RATE(860)*Y(237)+RATE(876)*Y(232)+RATE(883)*Y(131)&
    &+RATE(884)*Y(131)+RATE(903)*Y(116)+RATE(909)*Y(227)+RATE(939)*Y(300)&
    &+RATE(1076)*Y(97)+RATE(1159)*Y(111)+RATE(1219)*D*Y(16)*Y(100)+RATE(1231)&
    &*D*Y(16)*Y(117)+RATE(1241)*D*Y(16)*Y(236)+RATE(1254)*D*Y(16)*Y(116)&
    &+RATE(1265)*D*Y(16)*Y(133)+RATE(1268)*D*Y(16)*Y(161)+RATE(1269)*D*Y(16)&
    &*Y(221)+RATE(1270)*D*Y(16)*Y(300)+RATE(1271)*D*Y(16)*Y(56)+RATE(1274)*D&
    &*Y(16)*Y(320)+RATE(1276)*D*Y(16)*Y(277)+RATE(1281)*D*Y(16)*Y(46)&
    &+RATE(1312)*D*Y(18)*Y(232)+RATE(1313)*D*Y(18)*Y(131)+RATE(1322)*D*Y(18)&
    &*Y(116)+RATE(1330)*D*Y(18)*Y(161)+RATE(1332)*D*Y(18)*Y(300)+RATE(1334)*D&
    &*Y(18)*Y(320)+RATE(1336)*D*Y(18)*Y(277)+RATE(1341)*D*Y(18)*Y(235)&
    &+RATE(1349)*D*Y(67)*Y(100)+RATE(1354)*D*Y(67)*Y(117)+RATE(1357)*D*Y(67)&
    &*Y(162)+RATE(1360)*D*Y(67)*Y(236)+RATE(1363)*D*Y(67)*Y(161)+RATE(1363)*D&
    &*Y(67)*Y(161)+RATE(1370)*D*Y(68)*Y(116)+RATE(1371)*D*Y(68)*Y(161)&
    &+RATE(1374)*D*Y(75)*Y(100)+RATE(1379)*D*Y(75)*Y(117)+RATE(1387)*D*Y(75)&
    &*Y(161)+RATE(1396)*D*Y(76)*Y(116)+RATE(1400)*D*Y(80)*Y(133)+RATE(1434)*D&
    &*Y(20)*Y(100)+RATE(1455)*D*Y(20)*Y(117)+RATE(1475)*D*Y(20)*Y(232)&
    &+RATE(1477)*D*Y(20)*Y(116)+RATE(1486)*D*Y(20)*Y(161)+RATE(1487)*D*Y(20)&
    &*Y(161)+RATE(1491)*D*Y(20)*Y(46)+RATE(1493)*D*Y(20)*Y(300)+RATE(1497)*D&
    &*Y(20)*Y(277)+RATE(1515)*D*Y(21)*Y(232)+RATE(1516)*D*Y(21)*Y(131)&
    &+RATE(1527)*D*Y(21)*Y(116)+RATE(1537)*D*Y(21)*Y(300)+RATE(1544)*D*Y(24)&
    &*Y(100)+RATE(1561)*D*Y(24)*Y(117)+RATE(1579)*D*Y(24)*Y(116)+RATE(1588)*D&
    &*Y(24)*Y(161)+RATE(1591)*D*Y(24)*Y(46)+RATE(1592)*D*Y(24)*Y(46)&
    &+RATE(1604)*D*Y(25)*Y(232)+RATE(1610)*D*Y(25)*Y(116)+RATE(1614)*D*Y(25)&
    &*Y(300)+RATE(1626)*D*Y(32)*Y(116)+RATE(1636)*D*Y(32)*Y(46)+RATE(1655)*D&
    &*Y(33)*Y(116)+RATE(1661)*D*Y(33)*Y(300)+RATE(1679)*D*Y(41)*Y(100)&
    &+RATE(1752)*D*Y(82)*Y(116)+RATE(1756)*D*Y(82)*Y(133)+RATE(1758)*D*Y(82)&
    &*Y(161)+RATE(1772)*D*Y(83)*Y(116)+RATE(1773)*D*Y(83)*Y(161)+RATE(1787)*D&
    &*Y(100)*Y(131)+RATE(1788)*D*Y(100)*Y(183)+RATE(1789)*D*Y(100)*Y(116)&
    &+RATE(1790)*D*Y(100)*Y(133)+RATE(1791)*D*Y(100)*Y(161)+RATE(1792)*D&
    &*Y(100)*Y(163)+RATE(1804)*D*Y(2)*Y(100)+RATE(1826)*D*Y(2)*Y(218)&
    &+RATE(1830)*D*Y(2)*Y(232)+RATE(1837)*D*Y(2)*Y(116)+RATE(1859)*D*Y(2)&
    &*Y(221)+RATE(1861)*D*Y(2)*Y(300)+RATE(1931)*D*Y(4)*Y(116)+RATE(1933)*D&
    &*Y(4)*Y(227)+RATE(1937)*D*Y(4)*Y(300)+RATE(2044)*D*Y(8)*Y(116)+RATE(2060&
    &)*D*Y(131)*Y(165)+RATE(2064)*D*Y(132)*Y(333)+RATE(2065)*D*Y(132)*Y(333)&
    &+RATE(2077)*D*Y(62)*Y(100)+RATE(2093)*D*Y(62)*Y(117)+RATE(2130)*D*Y(63)&
    &*Y(116)+RATE(2175)*D*Y(10)*Y(267)+RATE(2209)*D*Y(142)*Y(333)+RATE(2241)&
    &*D*Y(90)*Y(100)+RATE(2249)*D*Y(90)*Y(117)+RATE(2264)*D*Y(91)*Y(116)&
    &+RATE(2281)*D*Y(116)*Y(132)+RATE(2284)*D*Y(116)*Y(162)+RATE(2286)*D&
    &*Y(116)*Y(165)+RATE(2287)*D*Y(116)*Y(116)+RATE(2287)*D*Y(116)*Y(116)&
    &+RATE(2288)*D*Y(116)*Y(116)+RATE(2290)*D*Y(116)*Y(133)+RATE(2292)*D&
    &*Y(116)*Y(161)+RATE(2294)*D*Y(117)*Y(333)+RATE(2295)*D*Y(117)*Y(206)&
    &+RATE(2296)*D*Y(117)*Y(214)+RATE(2297)*D*Y(117)*Y(159)+RATE(2298)*D&
    &*Y(117)*Y(233)+RATE(2299)*D*Y(117)*Y(131)+RATE(2300)*D*Y(117)*Y(261)&
    &+RATE(2301)*D*Y(117)*Y(183)+RATE(2302)*D*Y(117)*Y(116)+RATE(2303)*D&
    &*Y(117)*Y(324)+RATE(2304)*D*Y(117)*Y(173)+RATE(2305)*D*Y(117)*Y(265)&
    &+RATE(2306)*D*Y(117)*Y(300)+RATE(2307)*D*Y(117)*Y(316)+RATE(2308)*D&
    &*Y(117)*Y(163)+RATE(2309)*D*Y(117)*Y(277)+RATE(2310)*D*Y(117)*Y(135)&
    &+RATE(2311)*D*Y(117)*Y(166)+RATE(2312)*D*Y(117)*Y(121)+RATE(2313)*D&
    &*Y(117)*Y(235)+RATE(2314)*D*Y(117)*Y(302)+RATE(2316)*D*Y(243)*Y(333)&
    &+RATE(2317)*D*Y(243)*Y(333)+RATE(2351)*D*Y(13)*Y(218)+RATE(2367)*D*Y(13)&
    &*Y(232)+RATE(2393)*D*Y(13)*Y(116)+RATE(2422)*D*Y(13)*Y(300)+RATE(2450)*D&
    &*Y(92)*Y(117)+RATE(2455)*D*Y(227)*Y(16)+RATE(2460)*D*Y(119)*Y(333)&
    &+RATE(2528)*D*Y(27)*Y(232)+RATE(2531)*D*Y(27)*Y(116)+RATE(2586)*D*Y(28)&
    &*Y(116)+RATE(2596)*D*Y(28)*Y(300)+RATE(2612)*D*Y(104)*Y(116)+RATE(2613)&
    &*D*Y(104)*Y(300)+RATE(2621)*D*Y(35)*Y(100)+RATE(2632)*D*Y(35)*Y(117)&
    &+RATE(2674)*D*Y(36)*Y(232)+RATE(2699)*D*Y(43)*Y(100)+RATE(2716)*D*Y(43)&
    &*Y(117)+RATE(2756)*D*Y(54)*Y(100)+RATE(2781)*D*Y(54)*Y(117)+RATE(2803)*D&
    &*Y(55)*Y(116)+RATE(2825)*D*Y(46)*Y(100)+RATE(2854)*D*Y(46)*Y(67)&
    &+RATE(2855)*D*Y(46)*Y(80)+RATE(2863)*D*Y(46)*Y(75)+RATE(2864)*D*Y(46)&
    &*Y(197)+RATE(2865)*D*Y(46)*Y(284)+RATE(2866)*D*Y(46)*Y(308)+RATE(2868)*D&
    &*Y(46)*Y(82)+RATE(2870)*D*Y(46)*Y(232)+RATE(2871)*D*Y(46)*Y(233)&
    &+RATE(2878)*D*Y(46)*Y(90)+RATE(2881)*D*Y(46)*Y(116)+RATE(2882)*D*Y(46)&
    &*Y(244)+RATE(2898)*D*Y(46)*Y(221)+RATE(2901)*D*Y(46)*Y(300)+RATE(2906)*D&
    &*Y(46)*Y(292)+RATE(2907)*D*Y(46)*Y(314)+RATE(2908)*D*Y(46)*Y(207)&
    &+RATE(2941)*D*Y(48)*Y(232)+RATE(2947)*D*Y(48)*Y(116)+RATE(2961)*D*Y(161)&
    &*Y(221)+RATE(2968)*D*Y(162)*Y(80)+RATE(2976)*D*Y(301)*Y(333)+RATE(2979)&
    &*D*Y(56)*Y(100)+RATE(2985)*D*Y(56)*Y(117)+RATE(2995)*D*Y(56)*Y(80)&
    &+RATE(3001)*D*Y(56)*Y(233)+RATE(3007)*D*Y(56)*Y(90)+RATE(3008)*D*Y(56)&
    &*Y(116)+RATE(3039)*D*Y(57)*Y(116)+RATE(3058)*D*Y(163)*Y(116)+RATE(3066)&
    &*D*Y(165)*Y(300)+RATE(3076)*D*Y(105)*Y(117)+RATE(3077)*D*Y(105)*Y(232)&
    &+RATE(3082)*D*Y(106)*Y(300)+RATE(3108)*D*Y(278)*Y(80)
    YDOT(99) = PROD-LOSS
    LOSS = RATE(499)*D*Y(100)+RATE(875)*Y(100)+RATE(1219)*D*Y(16)*Y(100)&
    &+RATE(1349)*D*Y(67)*Y(100)+RATE(1374)*D*Y(75)*Y(100)+RATE(1376)*D*Y(75)&
    &*Y(100)+RATE(1434)*D*Y(20)*Y(100)+RATE(1447)*D*Y(20)*Y(100)+RATE(1544)*D&
    &*Y(24)*Y(100)+RATE(1554)*D*Y(24)*Y(100)+RATE(1679)*D*Y(41)*Y(100)&
    &+RATE(1684)*D*Y(41)*Y(100)+RATE(1787)*D*Y(100)*Y(131)+RATE(1788)*D*Y(100&
    &)*Y(183)+RATE(1789)*D*Y(100)*Y(116)+RATE(1790)*D*Y(100)*Y(133)+RATE(1791&
    &)*D*Y(100)*Y(161)+RATE(1792)*D*Y(100)*Y(163)+RATE(1793)*D*Y(100)*Y(333)&
    &+RATE(1794)*D*Y(100)*Y(131)+RATE(1795)*D*Y(100)*Y(183)+RATE(1796)*D&
    &*Y(100)*Y(320)+RATE(1804)*D*Y(2)*Y(100)+RATE(1960)*D*Y(6)*Y(100)&
    &+RATE(1961)*D*Y(6)*Y(100)+RATE(2077)*D*Y(62)*Y(100)+RATE(2087)*D*Y(62)&
    &*Y(100)+RATE(2241)*D*Y(90)*Y(100)+RATE(2621)*D*Y(35)*Y(100)+RATE(2628)*D&
    &*Y(35)*Y(100)+RATE(2699)*D*Y(43)*Y(100)+RATE(2708)*D*Y(43)*Y(100)&
    &+RATE(2756)*D*Y(54)*Y(100)+RATE(2771)*D*Y(54)*Y(100)+RATE(2825)*D*Y(46)&
    &*Y(100)+RATE(2979)*D*Y(56)*Y(100)+RATE(2982)*D*Y(56)*Y(100)
    PROD = RATE(101)*Y(99)+RATE(905)*Y(117)+RATE(1237)*D*Y(16)*Y(162)&
    &+RATE(1280)*D*Y(16)*Y(48)+RATE(1312)*D*Y(18)*Y(232)+RATE(1329)*D*Y(18)&
    &*Y(161)+RATE(1331)*D*Y(18)*Y(221)+RATE(1333)*D*Y(18)*Y(56)+RATE(1337)*D&
    &*Y(18)*Y(277)+RATE(1345)*D*Y(18)*Y(46)+RATE(1357)*D*Y(67)*Y(162)&
    &+RATE(1371)*D*Y(68)*Y(161)+RATE(1463)*D*Y(20)*Y(48)+RATE(1533)*D*Y(21)&
    &*Y(161)+RATE(1536)*D*Y(21)*Y(46)+RATE(1539)*D*Y(21)*Y(56)+RATE(1762)*D&
    &*Y(83)*Y(99)+RATE(1774)*D*Y(99)*Y(104)+RATE(1925)*D*Y(4)*Y(131)&
    &+RATE(1930)*D*Y(4)*Y(116)+RATE(2015)*D*Y(8)*Y(99)+RATE(2350)*D*Y(13)&
    &*Y(218)+RATE(2366)*D*Y(13)*Y(232)+RATE(2373)*D*Y(13)*Y(131)+RATE(2392)*D&
    &*Y(13)*Y(116)+RATE(2423)*D*Y(13)*Y(300)+RATE(2558)*D*Y(28)*Y(99)&
    &+RATE(2579)*D*Y(28)*Y(232)+RATE(2827)*D*Y(46)*Y(68)+RATE(2833)*D*Y(46)&
    &*Y(234)+RATE(2923)*D*Y(48)*Y(99)+RATE(2934)*D*Y(48)*Y(67)+RATE(2936)*D&
    &*Y(48)*Y(75)
    YDOT(100) = PROD-LOSS
    LOSS = RATE(144)*Y(101)+RATE(509)*D*Y(101)+RATE(882)*Y(101)&
    &+RATE(1253)*D*Y(16)*Y(101)+RATE(1832)*D*Y(2)*Y(101)+RATE(2530)*D*Y(27)&
    &*Y(101)+RATE(2873)*D*Y(46)*Y(101)
    PROD = RATE(233)*Y(107)/safeMantle+RATE(316)*D*Y(107)/safeMantle*Y(2&
    &)+RATE(399)*Y(107)/safeMantle+RATE(1081)*Y(107)+RATE(1164)*Y(113)&
    &+RATE(2514)*D*Y(27)*Y(123)+RATE(2524)*D*Y(27)*Y(32)
    YDOT(101) = PROD-LOSS
    LOSS = RATE(533)*D*Y(102)+RATE(1453)*D*Y(20)*Y(102)+RATE(1454)*D&
    &*Y(20)*Y(102)+RATE(1559)*D*Y(24)*Y(102)+RATE(1560)*D*Y(24)*Y(102)&
    &+RATE(2267)*D*Y(102)*Y(333)+RATE(2268)*D*Y(102)*Y(333)+RATE(2269)*D&
    &*Y(102)*Y(333)+RATE(2270)*D*Y(102)*Y(214)+RATE(2271)*D*Y(102)*Y(214)&
    &+RATE(2272)*D*Y(102)*Y(131)+RATE(2273)*D*Y(102)*Y(131)+RATE(2274)*D&
    &*Y(102)*Y(183)+RATE(2275)*D*Y(102)*Y(183)+RATE(2714)*D*Y(43)*Y(102)&
    &+RATE(2715)*D*Y(43)*Y(102)+RATE(2779)*D*Y(54)*Y(102)+RATE(2780)*D*Y(54)&
    &*Y(102)
    PROD = RATE(1395)*D*Y(76)*Y(90)+RATE(1397)*D*Y(76)*Y(92)+RATE(1412)&
    &*D*Y(81)*Y(90)+RATE(1413)*D*Y(81)*Y(92)+RATE(1526)*D*Y(21)*Y(90)&
    &+RATE(1528)*D*Y(21)*Y(92)+RATE(1651)*D*Y(33)*Y(214)+RATE(1688)*D*Y(41)&
    &*Y(91)+RATE(1733)*D*Y(53)*Y(90)+RATE(1736)*D*Y(53)*Y(92)+RATE(1966)*D&
    &*Y(6)*Y(91)+RATE(2129)*D*Y(63)*Y(90)+RATE(2132)*D*Y(63)*Y(92)+RATE(2173)&
    &*D*Y(10)*Y(90)+RATE(2179)*D*Y(10)*Y(92)+RATE(2227)*D*Y(66)*Y(90)&
    &+RATE(2228)*D*Y(66)*Y(92)+RATE(2246)*D*Y(90)*Y(132)+RATE(2247)*D*Y(90)&
    &*Y(142)+RATE(2248)*D*Y(90)*Y(189)+RATE(2249)*D*Y(90)*Y(117)+RATE(2250)*D&
    &*Y(90)*Y(145)+RATE(2251)*D*Y(90)*Y(174)+RATE(2252)*D*Y(90)*Y(307)&
    &+RATE(2253)*D*Y(90)*Y(120)+RATE(2254)*D*Y(90)*Y(176)+RATE(2262)*D*Y(91)&
    &*Y(90)+RATE(2264)*D*Y(91)*Y(116)+RATE(2265)*D*Y(91)*Y(92)+RATE(2447)*D&
    &*Y(92)*Y(132)+RATE(2448)*D*Y(92)*Y(142)+RATE(2449)*D*Y(92)*Y(189)&
    &+RATE(2450)*D*Y(92)*Y(117)+RATE(2451)*D*Y(92)*Y(145)+RATE(2452)*D*Y(92)&
    &*Y(174)+RATE(2453)*D*Y(92)*Y(120)+RATE(2454)*D*Y(92)*Y(176)+RATE(2578)*D&
    &*Y(28)*Y(41)+RATE(2626)*D*Y(35)*Y(33)+RATE(2683)*D*Y(36)*Y(90)+RATE(2685&
    &)*D*Y(36)*Y(92)+RATE(2746)*D*Y(44)*Y(90)+RATE(2748)*D*Y(44)*Y(92)&
    &+RATE(2769)*D*Y(54)*Y(199)+RATE(2778)*D*Y(54)*Y(91)+RATE(3038)*D*Y(57)&
    &*Y(90)+RATE(3041)*D*Y(57)*Y(92)
    YDOT(102) = PROD-LOSS
    LOSS = RATE(167)*Y(103)+RATE(562)*D*Y(103)+RATE(917)*Y(103)&
    &+RATE(1257)*D*Y(16)*Y(103)+RATE(1479)*D*Y(20)*Y(103)+RATE(1581)*D*Y(24)&
    &*Y(103)+RATE(2046)*D*Y(8)*Y(103)+RATE(2184)*D*Y(10)*Y(103)+RATE(2328)*D&
    &*Y(13)*Y(103)+RATE(2406)*D*Y(13)*Y(103)+RATE(2599)*D*Y(103)*Y(145)&
    &+RATE(2600)*D*Y(103)*Y(176)+RATE(2686)*D*Y(36)*Y(103)+RATE(2889)*D*Y(46)&
    &*Y(103)+RATE(2948)*D*Y(48)*Y(103)+RATE(3042)*D*Y(57)*Y(103)
    PROD = RATE(257)*Y(98)/safeMantle+RATE(340)*D*Y(98)/safeMantle*Y(2)&
    &+RATE(423)*Y(98)/safeMantle+RATE(759)*Y(29)*Y(29)+RATE(811)*Y(31)*Y(31)&
    &*bulkLayersReciprocal+RATE(1105)*Y(98)+RATE(1188)*Y(112)+RATE(1220)*D&
    &*Y(16)*Y(104)+RATE(1235)*D*Y(16)*Y(120)+RATE(1350)*D*Y(67)*Y(104)&
    &+RATE(1356)*D*Y(67)*Y(120)+RATE(1375)*D*Y(75)*Y(104)+RATE(1381)*D*Y(75)&
    &*Y(120)+RATE(1438)*D*Y(20)*Y(104)+RATE(1459)*D*Y(20)*Y(120)+RATE(1547)*D&
    &*Y(24)*Y(104)+RATE(1563)*D*Y(24)*Y(120)+RATE(1692)*D*Y(41)*Y(104)&
    &+RATE(1693)*D*Y(41)*Y(104)+RATE(1694)*D*Y(41)*Y(120)+RATE(1744)*D*Y(82)&
    &*Y(104)+RATE(1749)*D*Y(82)*Y(82)+RATE(1756)*D*Y(82)*Y(133)+RATE(1774)*D&
    &*Y(99)*Y(104)+RATE(1778)*D*Y(99)*Y(120)+RATE(2079)*D*Y(62)*Y(104)&
    &+RATE(2101)*D*Y(62)*Y(120)+RATE(2242)*D*Y(90)*Y(104)+RATE(2253)*D*Y(90)&
    &*Y(120)+RATE(2283)*D*Y(116)*Y(120)+RATE(2453)*D*Y(92)*Y(120)+RATE(2481)&
    &*D*Y(69)*Y(104)+RATE(2489)*D*Y(27)*Y(104)+RATE(2527)*D*Y(27)*Y(82)&
    &+RATE(2538)*D*Y(27)*Y(291)+RATE(2539)*D*Y(27)*Y(35)+RATE(2540)*D*Y(27)&
    &*Y(264)+RATE(2542)*D*Y(27)*Y(264)+RATE(2543)*D*Y(27)*Y(133)+RATE(2544)*D&
    &*Y(27)*Y(265)+RATE(2587)*D*Y(28)*Y(291)+RATE(2601)*D*Y(104)*Y(131)&
    &+RATE(2602)*D*Y(104)*Y(183)+RATE(2603)*D*Y(104)*Y(116)+RATE(2604)*D&
    &*Y(104)*Y(133)+RATE(2605)*D*Y(104)*Y(161)+RATE(2606)*D*Y(104)*Y(300)&
    &+RATE(2607)*D*Y(104)*Y(163)+RATE(2609)*D*Y(104)*Y(131)+RATE(2610)*D&
    &*Y(104)*Y(183)+RATE(2611)*D*Y(104)*Y(183)+RATE(2613)*D*Y(104)*Y(300)&
    &+RATE(2614)*D*Y(120)*Y(333)+RATE(2616)*D*Y(120)*Y(214)+RATE(2617)*D&
    &*Y(120)*Y(232)+RATE(2618)*D*Y(120)*Y(131)+RATE(2619)*D*Y(120)*Y(163)&
    &+RATE(2622)*D*Y(35)*Y(104)+RATE(2634)*D*Y(35)*Y(120)+RATE(2646)*D*Y(35)&
    &*Y(35)+RATE(2647)*D*Y(35)*Y(35)+RATE(2650)*D*Y(35)*Y(133)+RATE(2651)*D&
    &*Y(35)*Y(133)+RATE(2701)*D*Y(43)*Y(104)+RATE(2718)*D*Y(43)*Y(120)&
    &+RATE(2723)*D*Y(43)*Y(133)+RATE(2724)*D*Y(43)*Y(133)+RATE(2762)*D*Y(54)&
    &*Y(104)+RATE(2788)*D*Y(54)*Y(120)+RATE(2817)*D*Y(133)*Y(133)+RATE(2819)&
    &*D*Y(133)*Y(221)+RATE(2826)*D*Y(46)*Y(104)+RATE(2843)*D*Y(46)*Y(120)&
    &+RATE(2980)*D*Y(56)*Y(104)+RATE(2988)*D*Y(56)*Y(120)
    YDOT(103) = PROD-LOSS
    LOSS = RATE(563)*D*Y(104)+RATE(1220)*D*Y(16)*Y(104)+RATE(1350)*D&
    &*Y(67)*Y(104)+RATE(1375)*D*Y(75)*Y(104)+RATE(1438)*D*Y(20)*Y(104)&
    &+RATE(1547)*D*Y(24)*Y(104)+RATE(1692)*D*Y(41)*Y(104)+RATE(1693)*D*Y(41)&
    &*Y(104)+RATE(1744)*D*Y(82)*Y(104)+RATE(1774)*D*Y(99)*Y(104)+RATE(1972)*D&
    &*Y(6)*Y(104)+RATE(2079)*D*Y(62)*Y(104)+RATE(2100)*D*Y(62)*Y(104)&
    &+RATE(2242)*D*Y(90)*Y(104)+RATE(2481)*D*Y(69)*Y(104)+RATE(2489)*D*Y(27)&
    &*Y(104)+RATE(2601)*D*Y(104)*Y(131)+RATE(2602)*D*Y(104)*Y(183)+RATE(2603)&
    &*D*Y(104)*Y(116)+RATE(2604)*D*Y(104)*Y(133)+RATE(2605)*D*Y(104)*Y(161)&
    &+RATE(2606)*D*Y(104)*Y(300)+RATE(2607)*D*Y(104)*Y(163)+RATE(2608)*D&
    &*Y(104)*Y(333)+RATE(2609)*D*Y(104)*Y(131)+RATE(2610)*D*Y(104)*Y(183)&
    &+RATE(2611)*D*Y(104)*Y(183)+RATE(2612)*D*Y(104)*Y(116)+RATE(2613)*D&
    &*Y(104)*Y(300)+RATE(2622)*D*Y(35)*Y(104)+RATE(2701)*D*Y(43)*Y(104)&
    &+RATE(2762)*D*Y(54)*Y(104)+RATE(2826)*D*Y(46)*Y(104)+RATE(2842)*D*Y(46)&
    &*Y(104)+RATE(2980)*D*Y(56)*Y(104)
    PROD = RATE(2328)*D*Y(13)*Y(103)+RATE(2497)*D*Y(27)*Y(83)+RATE(2502)&
    &*D*Y(27)*Y(36)+RATE(2591)*D*Y(28)*Y(35)+RATE(2592)*D*Y(28)*Y(133)&
    &+RATE(2597)*D*Y(28)*Y(27)
    YDOT(104) = PROD-LOSS
    LOSS = RATE(190)*Y(105)+RATE(594)*D*Y(105)+RATE(946)*Y(105)&
    &+RATE(1296)*D*Y(18)*Y(105)+RATE(1401)*D*Y(80)*Y(105)+RATE(1414)*D*Y(81)&
    &*Y(105)+RATE(1504)*D*Y(21)*Y(105)+RATE(1905)*D*Y(4)*Y(105)+RATE(2118)*D&
    &*Y(63)*Y(105)+RATE(2200)*D*Y(10)*Y(105)+RATE(2231)*D*Y(66)*Y(105)&
    &+RATE(2332)*D*Y(13)*Y(105)+RATE(2796)*D*Y(55)*Y(105)+RATE(2918)*D*Y(46)&
    &*Y(105)+RATE(3017)*D*Y(56)*Y(105)+RATE(3048)*D*Y(57)*Y(105)+RATE(3069)*D&
    &*Y(105)*Y(234)+RATE(3070)*D*Y(105)*Y(132)+RATE(3071)*D*Y(105)*Y(184)&
    &+RATE(3072)*D*Y(105)*Y(174)+RATE(3073)*D*Y(105)*Y(134)+RATE(3074)*D&
    &*Y(105)*Y(162)+RATE(3075)*D*Y(105)*Y(165)+RATE(3076)*D*Y(105)*Y(117)&
    &+RATE(3077)*D*Y(105)*Y(232)+RATE(3078)*D*Y(105)*Y(99)+RATE(3079)*D*Y(105&
    &)*Y(133)+RATE(3080)*D*Y(105)*Y(161)
    PROD = RATE(191)*Y(207)+RATE(194)*Y(121)+RATE(198)*Y(235)+RATE(199)&
    &*Y(302)+RATE(274)*Y(110)/safeMantle+RATE(357)*D*Y(110)/safeMantle*Y(2)&
    &+RATE(440)*Y(110)/safeMantle+RATE(947)*Y(207)+RATE(948)*Y(292)+RATE(951)&
    &*Y(121)+RATE(961)*Y(235)+RATE(964)*Y(302)+RATE(1122)*Y(110)+RATE(1205)&
    &*Y(115)+RATE(1469)*D*Y(20)*Y(122)+RATE(1470)*D*Y(20)*Y(236)+RATE(2104)*D&
    &*Y(62)*Y(122)+RATE(2433)*D*Y(13)*Y(207)+RATE(2442)*D*Y(13)*Y(235)&
    &+RATE(2443)*D*Y(13)*Y(302)+RATE(2469)*D*Y(307)*Y(333)+RATE(2486)*D*Y(69)&
    &*Y(106)+RATE(2508)*D*Y(27)*Y(236)+RATE(2552)*D*Y(27)*Y(207)+RATE(2790)*D&
    &*Y(54)*Y(122)+RATE(2908)*D*Y(46)*Y(207)+RATE(3083)*D*Y(106)*Y(333)&
    &+RATE(3084)*D*Y(208)*Y(333)+RATE(3085)*D*Y(293)*Y(333)+RATE(3091)*D&
    &*Y(122)*Y(333)+RATE(3092)*D*Y(136)*Y(333)+RATE(3093)*D*Y(136)*Y(333)&
    &+RATE(3103)*D*Y(236)*Y(333)+RATE(3104)*D*Y(247)*Y(333)+RATE(3106)*D&
    &*Y(303)*Y(333)
    YDOT(105) = PROD-LOSS
    LOSS = RATE(595)*D*Y(106)+RATE(1383)*D*Y(75)*Y(106)+RATE(1468)*D&
    &*Y(20)*Y(106)+RATE(1870)*D*Y(2)*Y(106)+RATE(2005)*D*Y(6)*Y(106)&
    &+RATE(2103)*D*Y(62)*Y(106)+RATE(2486)*D*Y(69)*Y(106)+RATE(2917)*D*Y(46)&
    &*Y(106)+RATE(2992)*D*Y(56)*Y(106)+RATE(3081)*D*Y(106)*Y(159)+RATE(3082)&
    &*D*Y(106)*Y(300)+RATE(3083)*D*Y(106)*Y(333)
    PROD = RATE(190)*Y(105)+RATE(946)*Y(105)+RATE(952)*Y(122)+RATE(963)&
    &*Y(236)+RATE(1241)*D*Y(16)*Y(236)+RATE(1296)*D*Y(18)*Y(105)+RATE(1338)*D&
    &*Y(18)*Y(207)+RATE(1341)*D*Y(18)*Y(235)+RATE(1415)*D*Y(81)*Y(166)&
    &+RATE(1504)*D*Y(21)*Y(105)+RATE(1571)*D*Y(24)*Y(236)+RATE(1782)*D*Y(99)&
    &*Y(236)+RATE(1820)*D*Y(2)*Y(122)+RATE(1821)*D*Y(2)*Y(303)+RATE(1905)*D&
    &*Y(4)*Y(105)+RATE(1941)*D*Y(4)*Y(121)+RATE(2118)*D*Y(63)*Y(105)&
    &+RATE(2332)*D*Y(13)*Y(105)+RATE(2430)*D*Y(13)*Y(292)+RATE(2432)*D*Y(13)&
    &*Y(207)+RATE(2434)*D*Y(13)*Y(135)+RATE(2438)*D*Y(13)*Y(166)+RATE(2440)*D&
    &*Y(13)*Y(121)+RATE(2441)*D*Y(13)*Y(235)+RATE(2444)*D*Y(13)*Y(302)&
    &+RATE(2507)*D*Y(27)*Y(208)+RATE(2509)*D*Y(27)*Y(236)+RATE(2796)*D*Y(55)&
    &*Y(105)+RATE(2853)*D*Y(46)*Y(236)+RATE(3056)*D*Y(163)*Y(236)+RATE(3069)&
    &*D*Y(105)*Y(234)+RATE(3070)*D*Y(105)*Y(132)+RATE(3071)*D*Y(105)*Y(184)&
    &+RATE(3072)*D*Y(105)*Y(174)+RATE(3073)*D*Y(105)*Y(134)+RATE(3074)*D&
    &*Y(105)*Y(162)+RATE(3075)*D*Y(105)*Y(165)
    YDOT(106) = PROD-LOSS
    LOSS = RATE(233)*Y(107)/safeMantle+RATE(316)*D*Y(107)/safeMantle*Y(2&
    &)+RATE(399)*Y(107)/safeMantle+RATE(998)*Y(107)*totalSwap/safeMantle&
    &+RATE(1081)*Y(107)
    PROD = RATE(31)*Y(113)*bulkLayersReciprocal+RATE(509)*D*Y(101)
    YDOT(107) = PROD-LOSS
    LOSS = RATE(208)*Y(108)/safeMantle+RATE(291)*D*Y(108)/safeMantle*Y(2&
    &)+RATE(374)*Y(108)/safeMantle+RATE(973)*Y(108)*totalSwap/safeMantle&
    &+RATE(1056)*Y(108)
    PROD = RATE(6)*Y(114)*bulkLayersReciprocal+RATE(463)*D*Y(109)
    YDOT(108) = PROD-LOSS
    LOSS = RATE(117)*Y(109)+RATE(463)*D*Y(109)+RATE(840)*Y(109)&
    &+RATE(1421)*D*Y(109)*Y(165)+RATE(1472)*D*Y(20)*Y(109)+RATE(1649)*D*Y(33)&
    &*Y(109)+RATE(1748)*D*Y(82)*Y(109)+RATE(1916)*D*Y(4)*Y(109)+RATE(2029)*D&
    &*Y(8)*Y(109)+RATE(2339)*D*Y(13)*Y(109)+RATE(2340)*D*Y(13)*Y(109)&
    &+RATE(2341)*D*Y(13)*Y(109)+RATE(2512)*D*Y(27)*Y(109)+RATE(2857)*D*Y(46)&
    &*Y(109)+RATE(2858)*D*Y(46)*Y(109)+RATE(2859)*D*Y(46)*Y(109)+RATE(2860)*D&
    &*Y(46)*Y(109)+RATE(2935)*D*Y(48)*Y(109)+RATE(3111)*D*Y(278)*Y(109)&
    &+RATE(3112)*D*Y(278)*Y(109)+RATE(3113)*D*Y(278)*Y(109)
    PROD = RATE(208)*Y(108)/safeMantle+RATE(291)*D*Y(108)/safeMantle*Y(2&
    &)+RATE(374)*Y(108)/safeMantle+RATE(1056)*Y(108)+RATE(1139)*Y(114)&
    &+RATE(1415)*D*Y(81)*Y(166)+RATE(1417)*D*Y(81)*Y(166)+RATE(1474)*D*Y(20)&
    &*Y(41)+RATE(1619)*D*Y(32)*Y(32)+RATE(1650)*D*Y(33)*Y(237)+RATE(1651)*D&
    &*Y(33)*Y(214)+RATE(2157)*D*Y(10)*Y(237)+RATE(2513)*D*Y(27)*Y(123)&
    &+RATE(2997)*D*Y(56)*Y(123)
    YDOT(109) = PROD-LOSS
    LOSS = RATE(274)*Y(110)/safeMantle+RATE(357)*D*Y(110)/safeMantle*Y(2&
    &)+RATE(440)*Y(110)/safeMantle+RATE(668)*Y(110)*Y(3)+RATE(772)*Y(110)*Y(3&
    &)+RATE(1039)*Y(110)*totalSwap/safeMantle+RATE(1122)*Y(110)
    PROD = RATE(72)*Y(115)*bulkLayersReciprocal+RATE(594)*D*Y(105)&
    &+RATE(595)*D*Y(106)
    YDOT(110) = PROD-LOSS
    LOSS = RATE(26)*Y(111)*bulkLayersReciprocal+RATE(678)*Y(111)*Y(60)&
    &*bulkLayersReciprocal+RATE(690)*Y(5)*Y(111)*bulkLayersReciprocal&
    &+RATE(708)*Y(39)*Y(111)*bulkLayersReciprocal+RATE(782)*Y(111)*Y(60)&
    &*bulkLayersReciprocal+RATE(794)*Y(5)*Y(111)*bulkLayersReciprocal&
    &+RATE(812)*Y(39)*Y(111)*bulkLayersReciprocal+RATE(1159)*Y(111)
    PROD = RATE(673)*Y(19)*Y(225)*bulkLayersReciprocal+RATE(692)*Y(5)&
    &*Y(126)*bulkLayersReciprocal+RATE(699)*Y(5)*Y(228)*bulkLayersReciprocal&
    &+RATE(712)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(716)*Y(51)*Y(225)&
    &*bulkLayersReciprocal+RATE(993)*Y(97)*totalSwap/safeMantle
    YDOT(111) = PROD-LOSS
    LOSS = RATE(55)*Y(112)*bulkLayersReciprocal+RATE(1188)*Y(112)
    PROD = RATE(707)*Y(31)*Y(31)*bulkLayersReciprocal+RATE(1022)*Y(98)&
    &*totalSwap/safeMantle
    YDOT(112) = PROD-LOSS
    LOSS = RATE(31)*Y(113)*bulkLayersReciprocal+RATE(1164)*Y(113)
    PROD = RATE(998)*Y(107)*totalSwap/safeMantle
    YDOT(113) = PROD-LOSS
    LOSS = RATE(6)*Y(114)*bulkLayersReciprocal+RATE(1139)*Y(114)
    PROD = RATE(973)*Y(108)*totalSwap/safeMantle
    YDOT(114) = PROD-LOSS
    LOSS = RATE(72)*Y(115)*bulkLayersReciprocal+RATE(720)*Y(115)*Y(5)&
    &*bulkLayersReciprocal+RATE(824)*Y(115)*Y(5)*bulkLayersReciprocal&
    &+RATE(1205)*Y(115)
    PROD = RATE(1039)*Y(110)*totalSwap/safeMantle
    YDOT(115) = PROD-LOSS
    LOSS = RATE(155)*Y(116)+RATE(156)*Y(116)+RATE(534)*D*Y(116)+RATE(903&
    &)*Y(116)+RATE(904)*Y(116)+RATE(1254)*D*Y(16)*Y(116)+RATE(1288)*D*Y(18)&
    &*Y(116)+RATE(1322)*D*Y(18)*Y(116)+RATE(1365)*D*Y(68)*Y(116)+RATE(1370)*D&
    &*Y(68)*Y(116)+RATE(1396)*D*Y(76)*Y(116)+RATE(1404)*D*Y(81)*Y(116)&
    &+RATE(1477)*D*Y(20)*Y(116)+RATE(1499)*D*Y(21)*Y(116)+RATE(1527)*D*Y(21)&
    &*Y(116)+RATE(1579)*D*Y(24)*Y(116)+RATE(1610)*D*Y(25)*Y(116)+RATE(1626)*D&
    &*Y(32)*Y(116)+RATE(1643)*D*Y(33)*Y(116)+RATE(1655)*D*Y(33)*Y(116)&
    &+RATE(1734)*D*Y(53)*Y(116)+RATE(1752)*D*Y(82)*Y(116)+RATE(1765)*D*Y(83)&
    &*Y(116)+RATE(1772)*D*Y(83)*Y(116)+RATE(1789)*D*Y(100)*Y(116)+RATE(1837)&
    &*D*Y(2)*Y(116)+RATE(1838)*D*Y(2)*Y(116)+RATE(1887)*D*Y(4)*Y(116)&
    &+RATE(1930)*D*Y(4)*Y(116)+RATE(1931)*D*Y(4)*Y(116)+RATE(2020)*D*Y(8)&
    &*Y(116)+RATE(2044)*D*Y(8)*Y(116)+RATE(2112)*D*Y(63)*Y(116)+RATE(2130)*D&
    &*Y(63)*Y(116)+RATE(2131)*D*Y(63)*Y(116)+RATE(2174)*D*Y(10)*Y(116)&
    &+RATE(2263)*D*Y(91)*Y(116)+RATE(2264)*D*Y(91)*Y(116)+RATE(2276)*D*Y(116)&
    &*Y(132)+RATE(2277)*D*Y(116)*Y(184)+RATE(2278)*D*Y(116)*Y(162)+RATE(2279)&
    &*D*Y(116)*Y(165)+RATE(2280)*D*Y(116)*Y(236)+RATE(2281)*D*Y(116)*Y(132)&
    &+RATE(2282)*D*Y(116)*Y(145)+RATE(2283)*D*Y(116)*Y(120)+RATE(2284)*D&
    &*Y(116)*Y(162)+RATE(2285)*D*Y(116)*Y(176)+RATE(2286)*D*Y(116)*Y(165)&
    &+RATE(2287)*D*Y(116)*Y(116)+RATE(2287)*D*Y(116)*Y(116)+RATE(2288)*D&
    &*Y(116)*Y(116)+RATE(2288)*D*Y(116)*Y(116)+RATE(2289)*D*Y(116)*Y(144)&
    &+RATE(2290)*D*Y(116)*Y(133)+RATE(2291)*D*Y(116)*Y(161)+RATE(2292)*D&
    &*Y(116)*Y(161)+RATE(2293)*D*Y(116)*Y(175)+RATE(2302)*D*Y(117)*Y(116)&
    &+RATE(2392)*D*Y(13)*Y(116)+RATE(2393)*D*Y(13)*Y(116)+RATE(2394)*D*Y(13)&
    &*Y(116)+RATE(2531)*D*Y(27)*Y(116)+RATE(2532)*D*Y(27)*Y(116)+RATE(2533)*D&
    &*Y(27)*Y(116)+RATE(2563)*D*Y(28)*Y(116)+RATE(2586)*D*Y(28)*Y(116)&
    &+RATE(2603)*D*Y(104)*Y(116)+RATE(2612)*D*Y(104)*Y(116)+RATE(2684)*D*Y(36&
    &)*Y(116)+RATE(2729)*D*Y(44)*Y(116)+RATE(2747)*D*Y(44)*Y(116)+RATE(2793)&
    &*D*Y(55)*Y(116)+RATE(2803)*D*Y(55)*Y(116)+RATE(2880)*D*Y(46)*Y(116)&
    &+RATE(2881)*D*Y(46)*Y(116)+RATE(2927)*D*Y(48)*Y(116)+RATE(2947)*D*Y(48)&
    &*Y(116)+RATE(3008)*D*Y(56)*Y(116)+RATE(3024)*D*Y(57)*Y(116)+RATE(3039)*D&
    &*Y(57)*Y(116)+RATE(3040)*D*Y(57)*Y(116)+RATE(3057)*D*Y(163)*Y(116)&
    &+RATE(3058)*D*Y(163)*Y(116)
    PROD = RATE(132)*Y(237)+RATE(157)*Y(267)+RATE(244)*Y(118)/safeMantle&
    &+RATE(327)*D*Y(118)/safeMantle*Y(2)+RATE(410)*Y(118)/safeMantle+RATE(742&
    &)*Y(3)*Y(97)+RATE(747)*Y(3)*Y(130)+RATE(766)*Y(45)*Y(130)+RATE(794)*Y(5)&
    &*Y(111)*bulkLayersReciprocal+RATE(799)*Y(5)*Y(139)*bulkLayersReciprocal&
    &+RATE(818)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(861)*Y(237)+RATE(906)&
    &*Y(267)+RATE(1092)*Y(118)+RATE(1175)*Y(126)+RATE(1309)*D*Y(18)*Y(159)&
    &+RATE(1352)*D*Y(67)*Y(132)+RATE(1377)*D*Y(75)*Y(132)+RATE(1387)*D*Y(75)&
    &*Y(161)+RATE(1419)*D*Y(89)*Y(161)+RATE(1448)*D*Y(20)*Y(132)+RATE(1475)*D&
    &*Y(20)*Y(232)+RATE(1476)*D*Y(20)*Y(131)+RATE(1483)*D*Y(20)*Y(133)&
    &+RATE(1488)*D*Y(20)*Y(161)+RATE(1489)*D*Y(20)*Y(175)+RATE(1494)*D*Y(20)&
    &*Y(56)+RATE(1535)*D*Y(21)*Y(161)+RATE(1555)*D*Y(24)*Y(132)+RATE(1578)*D&
    &*Y(24)*Y(131)+RATE(1590)*D*Y(24)*Y(161)+RATE(1593)*D*Y(24)*Y(46)&
    &+RATE(1615)*D*Y(25)*Y(300)+RATE(1623)*D*Y(32)*Y(131)+RATE(1633)*D*Y(32)&
    &*Y(161)+RATE(1750)*D*Y(82)*Y(131)+RATE(1794)*D*Y(100)*Y(131)+RATE(1833)&
    &*D*Y(2)*Y(131)+RATE(2066)*D*Y(132)*Y(333)+RATE(2067)*D*Y(132)*Y(159)&
    &+RATE(2068)*D*Y(132)*Y(131)+RATE(2070)*D*Y(132)*Y(163)+RATE(2088)*D*Y(62&
    &)*Y(132)+RATE(2211)*D*Y(142)*Y(333)+RATE(2246)*D*Y(90)*Y(132)+RATE(2354)&
    &*D*Y(13)*Y(237)+RATE(2447)*D*Y(92)*Y(132)+RATE(2479)*D*Y(69)*Y(117)&
    &+RATE(2675)*D*Y(36)*Y(232)+RATE(2709)*D*Y(43)*Y(132)+RATE(2738)*D*Y(44)&
    &*Y(131)+RATE(2772)*D*Y(54)*Y(132)+RATE(2800)*D*Y(55)*Y(131)+RATE(2860)*D&
    &*Y(46)*Y(109)+RATE(2874)*D*Y(46)*Y(131)+RATE(3003)*D*Y(56)*Y(131)&
    &+RATE(3110)*D*Y(278)*Y(80)+RATE(3113)*D*Y(278)*Y(109)
    YDOT(116) = PROD-LOSS
    LOSS = RATE(535)*D*Y(117)+RATE(905)*Y(117)+RATE(1231)*D*Y(16)*Y(117)&
    &+RATE(1354)*D*Y(67)*Y(117)+RATE(1379)*D*Y(75)*Y(117)+RATE(1455)*D*Y(20)&
    &*Y(117)+RATE(1561)*D*Y(24)*Y(117)+RATE(2093)*D*Y(62)*Y(117)+RATE(2249)*D&
    &*Y(90)*Y(117)+RATE(2294)*D*Y(117)*Y(333)+RATE(2295)*D*Y(117)*Y(206)&
    &+RATE(2296)*D*Y(117)*Y(214)+RATE(2297)*D*Y(117)*Y(159)+RATE(2298)*D&
    &*Y(117)*Y(233)+RATE(2299)*D*Y(117)*Y(131)+RATE(2300)*D*Y(117)*Y(261)&
    &+RATE(2301)*D*Y(117)*Y(183)+RATE(2302)*D*Y(117)*Y(116)+RATE(2303)*D&
    &*Y(117)*Y(324)+RATE(2304)*D*Y(117)*Y(173)+RATE(2305)*D*Y(117)*Y(265)&
    &+RATE(2306)*D*Y(117)*Y(300)+RATE(2307)*D*Y(117)*Y(316)+RATE(2308)*D&
    &*Y(117)*Y(163)+RATE(2309)*D*Y(117)*Y(277)+RATE(2310)*D*Y(117)*Y(135)&
    &+RATE(2311)*D*Y(117)*Y(166)+RATE(2312)*D*Y(117)*Y(121)+RATE(2313)*D&
    &*Y(117)*Y(235)+RATE(2314)*D*Y(117)*Y(302)+RATE(2450)*D*Y(92)*Y(117)&
    &+RATE(2479)*D*Y(69)*Y(117)+RATE(2632)*D*Y(35)*Y(117)+RATE(2716)*D*Y(43)&
    &*Y(117)+RATE(2781)*D*Y(54)*Y(117)+RATE(2985)*D*Y(56)*Y(117)+RATE(2986)*D&
    &*Y(56)*Y(117)+RATE(3076)*D*Y(105)*Y(117)
    PROD = RATE(156)*Y(116)+RATE(886)*Y(131)+RATE(904)*Y(116)+RATE(1229)&
    &*D*Y(16)*Y(66)+RATE(1288)*D*Y(18)*Y(116)+RATE(1314)*D*Y(18)*Y(131)&
    &+RATE(1316)*D*Y(18)*Y(62)+RATE(1365)*D*Y(68)*Y(116)+RATE(1376)*D*Y(75)&
    &*Y(100)+RATE(1404)*D*Y(81)*Y(116)+RATE(1431)*D*Y(20)*Y(46)+RATE(1447)*D&
    &*Y(20)*Y(100)+RATE(1464)*D*Y(20)*Y(162)+RATE(1470)*D*Y(20)*Y(236)&
    &+RATE(1499)*D*Y(21)*Y(116)+RATE(1515)*D*Y(21)*Y(232)+RATE(1518)*D*Y(21)&
    &*Y(131)+RATE(1521)*D*Y(21)*Y(62)+RATE(1534)*D*Y(21)*Y(161)+RATE(1554)*D&
    &*Y(24)*Y(100)+RATE(1605)*D*Y(25)*Y(131)+RATE(1612)*D*Y(25)*Y(161)&
    &+RATE(1613)*D*Y(25)*Y(46)+RATE(1643)*D*Y(33)*Y(116)+RATE(1653)*D*Y(33)&
    &*Y(131)+RATE(1660)*D*Y(33)*Y(46)+RATE(1684)*D*Y(41)*Y(100)+RATE(1715)*D&
    &*Y(42)*Y(99)+RATE(1729)*D*Y(53)*Y(99)+RATE(1765)*D*Y(83)*Y(116)&
    &+RATE(1770)*D*Y(83)*Y(131)+RATE(1775)*D*Y(99)*Y(196)+RATE(1776)*D*Y(99)&
    &*Y(243)+RATE(1777)*D*Y(99)*Y(145)+RATE(1778)*D*Y(99)*Y(120)+RATE(1779)*D&
    &*Y(99)*Y(176)+RATE(1781)*D*Y(99)*Y(167)+RATE(1789)*D*Y(100)*Y(116)&
    &+RATE(1794)*D*Y(100)*Y(131)+RATE(1795)*D*Y(100)*Y(183)+RATE(1887)*D*Y(4)&
    &*Y(116)+RATE(1922)*D*Y(4)*Y(159)+RATE(1924)*D*Y(4)*Y(232)+RATE(1926)*D&
    &*Y(4)*Y(131)+RATE(1948)*D*Y(6)*Y(119)+RATE(1960)*D*Y(6)*Y(100)+RATE(2020&
    &)*D*Y(8)*Y(116)+RATE(2038)*D*Y(8)*Y(99)+RATE(2040)*D*Y(8)*Y(131)&
    &+RATE(2058)*D*Y(131)*Y(162)+RATE(2061)*D*Y(131)*Y(165)+RATE(2069)*D&
    &*Y(132)*Y(161)+RATE(2083)*D*Y(62)*Y(199)+RATE(2084)*D*Y(62)*Y(309)&
    &+RATE(2086)*D*Y(62)*Y(83)+RATE(2087)*D*Y(62)*Y(100)+RATE(2112)*D*Y(63)&
    &*Y(116)+RATE(2124)*D*Y(63)*Y(99)+RATE(2165)*D*Y(10)*Y(99)+RATE(2176)*D&
    &*Y(10)*Y(267)+RATE(2260)*D*Y(91)*Y(99)+RATE(2276)*D*Y(116)*Y(132)&
    &+RATE(2277)*D*Y(116)*Y(184)+RATE(2278)*D*Y(116)*Y(162)+RATE(2279)*D&
    &*Y(116)*Y(165)+RATE(2280)*D*Y(116)*Y(236)+RATE(2374)*D*Y(13)*Y(131)&
    &+RATE(2563)*D*Y(28)*Y(116)+RATE(2581)*D*Y(28)*Y(131)+RATE(2603)*D*Y(104)&
    &*Y(116)+RATE(2609)*D*Y(104)*Y(131)+RATE(2628)*D*Y(35)*Y(100)+RATE(2676)&
    &*D*Y(36)*Y(99)+RATE(2678)*D*Y(36)*Y(131)+RATE(2708)*D*Y(43)*Y(100)&
    &+RATE(2729)*D*Y(44)*Y(116)+RATE(2771)*D*Y(54)*Y(100)+RATE(2793)*D*Y(55)&
    &*Y(116)+RATE(2828)*D*Y(46)*Y(76)+RATE(2829)*D*Y(46)*Y(81)+RATE(2837)*D&
    &*Y(46)*Y(243)+RATE(2839)*D*Y(46)*Y(245)+RATE(2927)*D*Y(48)*Y(116)&
    &+RATE(2942)*D*Y(48)*Y(131)+RATE(2945)*D*Y(48)*Y(90)+RATE(2968)*D*Y(162)&
    &*Y(80)+RATE(2982)*D*Y(56)*Y(100)+RATE(3024)*D*Y(57)*Y(116)+RATE(3034)*D&
    &*Y(57)*Y(99)+RATE(3109)*D*Y(278)*Y(80)
    YDOT(117) = PROD-LOSS
    LOSS = RATE(93)*Y(118)+RATE(244)*Y(118)/safeMantle+RATE(327)*D*Y(118&
    &)/safeMantle*Y(2)+RATE(410)*Y(118)/safeMantle+RATE(623)*Y(34)*Y(118)&
    &+RATE(639)*Y(3)*Y(118)+RATE(640)*Y(3)*Y(118)+RATE(650)*Y(118)*Y(58)&
    &+RATE(659)*Y(45)*Y(118)+RATE(660)*Y(45)*Y(118)+RATE(727)*Y(34)*Y(118)&
    &+RATE(743)*Y(3)*Y(118)+RATE(744)*Y(3)*Y(118)+RATE(754)*Y(118)*Y(58)&
    &+RATE(763)*Y(45)*Y(118)+RATE(764)*Y(45)*Y(118)+RATE(1009)*Y(118)&
    &*totalSwap/safeMantle+RATE(1092)*Y(118)
    PROD = RATE(42)*Y(126)*bulkLayersReciprocal+RATE(89)*Y(130)+RATE(534&
    &)*D*Y(116)+RATE(535)*D*Y(117)+RATE(549)*D*Y(119)+RATE(638)*Y(3)*Y(97)&
    &+RATE(643)*Y(3)*Y(130)+RATE(662)*Y(45)*Y(130)
    YDOT(118) = PROD-LOSS
    LOSS = RATE(549)*D*Y(119)+RATE(1948)*D*Y(6)*Y(119)+RATE(2460)*D&
    &*Y(119)*Y(333)
    PROD = RATE(1317)*D*Y(18)*Y(62)+RATE(1961)*D*Y(6)*Y(100)+RATE(2166)&
    &*D*Y(10)*Y(99)
    YDOT(119) = PROD-LOSS
    LOSS = RATE(564)*D*Y(120)+RATE(1235)*D*Y(16)*Y(120)+RATE(1356)*D&
    &*Y(67)*Y(120)+RATE(1381)*D*Y(75)*Y(120)+RATE(1459)*D*Y(20)*Y(120)&
    &+RATE(1563)*D*Y(24)*Y(120)+RATE(1694)*D*Y(41)*Y(120)+RATE(1778)*D*Y(99)&
    &*Y(120)+RATE(2101)*D*Y(62)*Y(120)+RATE(2253)*D*Y(90)*Y(120)+RATE(2283)*D&
    &*Y(116)*Y(120)+RATE(2453)*D*Y(92)*Y(120)+RATE(2614)*D*Y(120)*Y(333)&
    &+RATE(2615)*D*Y(120)*Y(333)+RATE(2616)*D*Y(120)*Y(214)+RATE(2617)*D&
    &*Y(120)*Y(232)+RATE(2618)*D*Y(120)*Y(131)+RATE(2619)*D*Y(120)*Y(163)&
    &+RATE(2634)*D*Y(35)*Y(120)+RATE(2718)*D*Y(43)*Y(120)+RATE(2788)*D*Y(54)&
    &*Y(120)+RATE(2843)*D*Y(46)*Y(120)+RATE(2988)*D*Y(56)*Y(120)
    PROD = RATE(1972)*D*Y(6)*Y(104)+RATE(2046)*D*Y(8)*Y(103)+RATE(2100)&
    &*D*Y(62)*Y(104)+RATE(2184)*D*Y(10)*Y(103)+RATE(2503)*D*Y(27)*Y(44)&
    &+RATE(2589)*D*Y(28)*Y(54)+RATE(2599)*D*Y(103)*Y(145)+RATE(2600)*D*Y(103)&
    &*Y(176)+RATE(2612)*D*Y(104)*Y(116)+RATE(2686)*D*Y(36)*Y(103)+RATE(2690)&
    &*D*Y(36)*Y(133)+RATE(3042)*D*Y(57)*Y(103)
    YDOT(120) = PROD-LOSS
    LOSS = RATE(194)*Y(121)+RATE(602)*D*Y(121)+RATE(951)*Y(121)&
    &+RATE(1277)*D*Y(16)*Y(121)+RATE(1340)*D*Y(18)*Y(121)+RATE(1912)*D*Y(4)&
    &*Y(121)+RATE(1941)*D*Y(4)*Y(121)+RATE(2204)*D*Y(10)*Y(121)+RATE(2233)*D&
    &*Y(66)*Y(121)+RATE(2312)*D*Y(117)*Y(121)+RATE(2440)*D*Y(13)*Y(121)&
    &+RATE(2914)*D*Y(46)*Y(121)+RATE(3049)*D*Y(57)*Y(121)+RATE(3089)*D*Y(121)&
    &*Y(165)+RATE(3090)*D*Y(121)*Y(165)
    PROD = RATE(195)*Y(135)+RATE(278)*Y(125)/safeMantle+RATE(361)*D&
    &*Y(125)/safeMantle*Y(2)+RATE(444)*Y(125)/safeMantle+RATE(772)*Y(110)*Y(3&
    &)+RATE(824)*Y(115)*Y(5)*bulkLayersReciprocal+RATE(954)*Y(135)+RATE(957)&
    &*Y(146)+RATE(960)*Y(166)+RATE(1126)*Y(125)+RATE(1209)*Y(128)+RATE(3094)&
    &*D*Y(136)*Y(333)+RATE(3098)*D*Y(147)*Y(333)
    YDOT(121) = PROD-LOSS
    LOSS = RATE(603)*D*Y(122)+RATE(952)*Y(122)+RATE(1240)*D*Y(16)*Y(122)&
    &+RATE(1469)*D*Y(20)*Y(122)+RATE(1820)*D*Y(2)*Y(122)+RATE(2006)*D*Y(6)&
    &*Y(122)+RATE(2104)*D*Y(62)*Y(122)+RATE(2790)*D*Y(54)*Y(122)+RATE(2850)*D&
    &*Y(46)*Y(122)+RATE(3091)*D*Y(122)*Y(333)
    PROD = RATE(1416)*D*Y(81)*Y(166)+RATE(1870)*D*Y(2)*Y(106)+RATE(1912)&
    &*D*Y(4)*Y(121)+RATE(1938)*D*Y(4)*Y(135)+RATE(2200)*D*Y(10)*Y(105)&
    &+RATE(2231)*D*Y(66)*Y(105)+RATE(2435)*D*Y(13)*Y(135)+RATE(2436)*D*Y(13)&
    &*Y(146)+RATE(2439)*D*Y(13)*Y(166)+RATE(3048)*D*Y(57)*Y(105)+RATE(3076)*D&
    &*Y(105)*Y(117)+RATE(3089)*D*Y(121)*Y(165)
    YDOT(122) = PROD-LOSS
    LOSS = RATE(118)*Y(123)+RATE(464)*D*Y(123)+RATE(841)*Y(123)&
    &+RATE(1243)*D*Y(16)*Y(123)+RATE(2513)*D*Y(27)*Y(123)+RATE(2514)*D*Y(27)&
    &*Y(123)+RATE(2861)*D*Y(46)*Y(123)+RATE(2862)*D*Y(46)*Y(123)+RATE(2997)*D&
    &*Y(56)*Y(123)
    PROD = RATE(209)*Y(124)/safeMantle+RATE(292)*D*Y(124)/safeMantle*Y(2&
    &)+RATE(375)*Y(124)/safeMantle+RATE(1057)*Y(124)+RATE(1140)*Y(127)&
    &+RATE(1416)*D*Y(81)*Y(166)+RATE(1620)*D*Y(32)*Y(32)
    YDOT(123) = PROD-LOSS
    LOSS = RATE(209)*Y(124)/safeMantle+RATE(292)*D*Y(124)/safeMantle*Y(2&
    &)+RATE(375)*Y(124)/safeMantle+RATE(974)*Y(124)*totalSwap/safeMantle&
    &+RATE(1057)*Y(124)
    PROD = RATE(7)*Y(127)*bulkLayersReciprocal+RATE(464)*D*Y(123)
    YDOT(124) = PROD-LOSS
    LOSS = RATE(278)*Y(125)/safeMantle+RATE(361)*D*Y(125)/safeMantle*Y(2&
    &)+RATE(444)*Y(125)/safeMantle+RATE(669)*Y(125)*Y(3)+RATE(773)*Y(125)*Y(3&
    &)+RATE(1043)*Y(125)*totalSwap/safeMantle+RATE(1126)*Y(125)
    PROD = RATE(76)*Y(128)*bulkLayersReciprocal+RATE(602)*D*Y(121)&
    &+RATE(603)*D*Y(122)+RATE(668)*Y(110)*Y(3)
    YDOT(125) = PROD-LOSS
    LOSS = RATE(42)*Y(126)*bulkLayersReciprocal+RATE(675)*Y(38)*Y(126)&
    &*bulkLayersReciprocal+RATE(691)*Y(5)*Y(126)*bulkLayersReciprocal&
    &+RATE(692)*Y(5)*Y(126)*bulkLayersReciprocal+RATE(702)*Y(126)*Y(60)&
    &*bulkLayersReciprocal+RATE(711)*Y(50)*Y(126)*bulkLayersReciprocal&
    &+RATE(712)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(779)*Y(38)*Y(126)&
    &*bulkLayersReciprocal+RATE(795)*Y(5)*Y(126)*bulkLayersReciprocal&
    &+RATE(796)*Y(5)*Y(126)*bulkLayersReciprocal+RATE(806)*Y(126)*Y(60)&
    &*bulkLayersReciprocal+RATE(815)*Y(50)*Y(126)*bulkLayersReciprocal&
    &+RATE(816)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(1175)*Y(126)
    PROD = RATE(690)*Y(5)*Y(111)*bulkLayersReciprocal+RATE(695)*Y(5)&
    &*Y(139)*bulkLayersReciprocal+RATE(714)*Y(50)*Y(139)*bulkLayersReciprocal&
    &+RATE(1009)*Y(118)*totalSwap/safeMantle
    YDOT(126) = PROD-LOSS
    LOSS = RATE(7)*Y(127)*bulkLayersReciprocal+RATE(1140)*Y(127)
    PROD = RATE(974)*Y(124)*totalSwap/safeMantle
    YDOT(127) = PROD-LOSS
    LOSS = RATE(76)*Y(128)*bulkLayersReciprocal+RATE(721)*Y(128)*Y(5)&
    &*bulkLayersReciprocal+RATE(825)*Y(128)*Y(5)*bulkLayersReciprocal&
    &+RATE(1209)*Y(128)
    PROD = RATE(720)*Y(115)*Y(5)*bulkLayersReciprocal+RATE(1043)*Y(125)&
    &*totalSwap/safeMantle
    YDOT(128) = PROD-LOSS
    LOSS = RATE(263)*Y(129)/safeMantle+RATE(346)*D*Y(129)/safeMantle*Y(2&
    &)+RATE(429)*Y(129)/safeMantle+RATE(1028)*Y(129)*totalSwap/safeMantle&
    &+RATE(1111)*Y(129)
    PROD = RATE(61)*Y(138)*bulkLayersReciprocal+RATE(96)*Y(148)+RATE(574&
    &)*D*Y(133)+RATE(575)*D*Y(134)+RATE(664)*Y(47)*Y(222)
    YDOT(129) = PROD-LOSS
    LOSS = RATE(89)*Y(130)+RATE(90)*Y(130)+RATE(234)*Y(130)/safeMantle&
    &+RATE(317)*D*Y(130)/safeMantle*Y(2)+RATE(400)*Y(130)/safeMantle+RATE(641&
    &)*Y(3)*Y(130)+RATE(642)*Y(3)*Y(130)+RATE(643)*Y(3)*Y(130)+RATE(648)&
    &*Y(130)*Y(47)+RATE(661)*Y(45)*Y(130)+RATE(662)*Y(45)*Y(130)+RATE(745)&
    &*Y(3)*Y(130)+RATE(746)*Y(3)*Y(130)+RATE(747)*Y(3)*Y(130)+RATE(752)*Y(130&
    &)*Y(47)+RATE(765)*Y(45)*Y(130)+RATE(766)*Y(45)*Y(130)+RATE(999)*Y(130)&
    &*totalSwap/safeMantle+RATE(1082)*Y(130)
    PROD = RATE(32)*Y(139)*bulkLayersReciprocal+RATE(510)*D*Y(131)&
    &+RATE(511)*D*Y(132)+RATE(639)*Y(3)*Y(118)
    YDOT(130) = PROD-LOSS
    LOSS = RATE(145)*Y(131)+RATE(510)*D*Y(131)+RATE(883)*Y(131)+RATE(884&
    &)*Y(131)+RATE(885)*Y(131)+RATE(886)*Y(131)+RATE(1286)*D*Y(18)*Y(131)&
    &+RATE(1313)*D*Y(18)*Y(131)+RATE(1314)*D*Y(18)*Y(131)+RATE(1393)*D*Y(76)&
    &*Y(131)+RATE(1402)*D*Y(81)*Y(131)+RATE(1476)*D*Y(20)*Y(131)+RATE(1516)*D&
    &*Y(21)*Y(131)+RATE(1517)*D*Y(21)*Y(131)+RATE(1518)*D*Y(21)*Y(131)&
    &+RATE(1578)*D*Y(24)*Y(131)+RATE(1605)*D*Y(25)*Y(131)+RATE(1623)*D*Y(32)&
    &*Y(131)+RATE(1653)*D*Y(33)*Y(131)+RATE(1705)*D*Y(42)*Y(131)+RATE(1716)*D&
    &*Y(42)*Y(131)+RATE(1730)*D*Y(53)*Y(131)+RATE(1750)*D*Y(82)*Y(131)&
    &+RATE(1763)*D*Y(83)*Y(131)+RATE(1770)*D*Y(83)*Y(131)+RATE(1787)*D*Y(100)&
    &*Y(131)+RATE(1794)*D*Y(100)*Y(131)+RATE(1833)*D*Y(2)*Y(131)+RATE(1881)*D&
    &*Y(4)*Y(131)+RATE(1925)*D*Y(4)*Y(131)+RATE(1926)*D*Y(4)*Y(131)+RATE(2016&
    &)*D*Y(8)*Y(131)+RATE(2040)*D*Y(8)*Y(131)+RATE(2055)*D*Y(131)*Y(162)&
    &+RATE(2056)*D*Y(131)*Y(189)+RATE(2057)*D*Y(131)*Y(145)+RATE(2058)*D&
    &*Y(131)*Y(162)+RATE(2059)*D*Y(131)*Y(176)+RATE(2060)*D*Y(131)*Y(165)&
    &+RATE(2061)*D*Y(131)*Y(165)+RATE(2068)*D*Y(132)*Y(131)+RATE(2110)*D*Y(63&
    &)*Y(131)+RATE(2125)*D*Y(63)*Y(131)+RATE(2169)*D*Y(10)*Y(131)+RATE(2225)&
    &*D*Y(66)*Y(131)+RATE(2261)*D*Y(91)*Y(131)+RATE(2272)*D*Y(102)*Y(131)&
    &+RATE(2273)*D*Y(102)*Y(131)+RATE(2299)*D*Y(117)*Y(131)+RATE(2325)*D*Y(13&
    &)*Y(131)+RATE(2373)*D*Y(13)*Y(131)+RATE(2374)*D*Y(13)*Y(131)+RATE(2375)&
    &*D*Y(13)*Y(131)+RATE(2559)*D*Y(28)*Y(131)+RATE(2581)*D*Y(28)*Y(131)&
    &+RATE(2582)*D*Y(28)*Y(131)+RATE(2601)*D*Y(104)*Y(131)+RATE(2609)*D*Y(104&
    &)*Y(131)+RATE(2618)*D*Y(120)*Y(131)+RATE(2661)*D*Y(36)*Y(131)+RATE(2677)&
    &*D*Y(36)*Y(131)+RATE(2678)*D*Y(36)*Y(131)+RATE(2727)*D*Y(43)*Y(131)&
    &+RATE(2737)*D*Y(44)*Y(131)+RATE(2738)*D*Y(44)*Y(131)+RATE(2800)*D*Y(55)&
    &*Y(131)+RATE(2874)*D*Y(46)*Y(131)+RATE(2924)*D*Y(48)*Y(131)+RATE(2942)*D&
    &*Y(48)*Y(131)+RATE(3003)*D*Y(56)*Y(131)+RATE(3004)*D*Y(56)*Y(131)&
    &+RATE(3021)*D*Y(57)*Y(131)+RATE(3035)*D*Y(57)*Y(131)
    PROD = RATE(134)*Y(159)+RATE(234)*Y(130)/safeMantle+RATE(317)*D&
    &*Y(130)/safeMantle*Y(2)+RATE(400)*Y(130)/safeMantle+RATE(743)*Y(3)*Y(118&
    &)+RATE(795)*Y(5)*Y(126)*bulkLayersReciprocal+RATE(863)*Y(159)+RATE(1082)&
    &*Y(130)+RATE(1165)*Y(139)+RATE(1419)*D*Y(89)*Y(161)+RATE(1435)*D*Y(20)&
    &*Y(132)+RATE(1450)*D*Y(20)*Y(142)+RATE(1510)*D*Y(21)*Y(159)+RATE(1545)*D&
    &*Y(24)*Y(132)+RATE(1571)*D*Y(24)*Y(236)+RATE(1582)*D*Y(24)*Y(264)&
    &+RATE(1583)*D*Y(24)*Y(133)+RATE(1589)*D*Y(24)*Y(161)+RATE(1595)*D*Y(24)&
    &*Y(56)+RATE(1630)*D*Y(32)*Y(264)+RATE(1632)*D*Y(32)*Y(161)+RATE(1637)*D&
    &*Y(32)*Y(46)+RATE(1639)*D*Y(32)*Y(56)+RATE(1672)*D*Y(159)*Y(318)&
    &+RATE(1677)*D*Y(172)*Y(333)+RATE(2062)*D*Y(132)*Y(163)+RATE(2071)*D&
    &*Y(132)*Y(333)+RATE(2091)*D*Y(62)*Y(142)+RATE(2156)*D*Y(10)*Y(237)&
    &+RATE(2210)*D*Y(142)*Y(333)+RATE(2212)*D*Y(142)*Y(159)+RATE(2213)*D&
    &*Y(142)*Y(183)+RATE(2247)*D*Y(90)*Y(142)+RATE(2276)*D*Y(116)*Y(132)&
    &+RATE(2288)*D*Y(116)*Y(116)+RATE(2289)*D*Y(116)*Y(144)+RATE(2293)*D&
    &*Y(116)*Y(175)+RATE(2448)*D*Y(92)*Y(142)+RATE(2477)*D*Y(69)*Y(132)&
    &+RATE(2711)*D*Y(43)*Y(142)+RATE(2757)*D*Y(54)*Y(132)+RATE(2775)*D*Y(54)&
    &*Y(142)+RATE(2808)*D*Y(133)*Y(132)+RATE(2859)*D*Y(46)*Y(109)+RATE(2862)&
    &*D*Y(46)*Y(123)+RATE(3070)*D*Y(105)*Y(132)+RATE(3111)*D*Y(278)*Y(109)
    YDOT(131) = PROD-LOSS
    LOSS = RATE(511)*D*Y(132)+RATE(1352)*D*Y(67)*Y(132)+RATE(1377)*D&
    &*Y(75)*Y(132)+RATE(1435)*D*Y(20)*Y(132)+RATE(1448)*D*Y(20)*Y(132)&
    &+RATE(1545)*D*Y(24)*Y(132)+RATE(1555)*D*Y(24)*Y(132)+RATE(1686)*D*Y(41)&
    &*Y(132)+RATE(2062)*D*Y(132)*Y(163)+RATE(2063)*D*Y(132)*Y(333)+RATE(2064)&
    &*D*Y(132)*Y(333)+RATE(2065)*D*Y(132)*Y(333)+RATE(2066)*D*Y(132)*Y(333)&
    &+RATE(2067)*D*Y(132)*Y(159)+RATE(2068)*D*Y(132)*Y(131)+RATE(2069)*D&
    &*Y(132)*Y(161)+RATE(2070)*D*Y(132)*Y(163)+RATE(2071)*D*Y(132)*Y(333)&
    &+RATE(2088)*D*Y(62)*Y(132)+RATE(2246)*D*Y(90)*Y(132)+RATE(2276)*D*Y(116)&
    &*Y(132)+RATE(2281)*D*Y(116)*Y(132)+RATE(2447)*D*Y(92)*Y(132)+RATE(2477)&
    &*D*Y(69)*Y(132)+RATE(2629)*D*Y(35)*Y(132)+RATE(2709)*D*Y(43)*Y(132)&
    &+RATE(2757)*D*Y(54)*Y(132)+RATE(2772)*D*Y(54)*Y(132)+RATE(2808)*D*Y(133)&
    &*Y(132)+RATE(3070)*D*Y(105)*Y(132)
    PROD = RATE(885)*Y(131)+RATE(1286)*D*Y(18)*Y(131)+RATE(1402)*D*Y(81)&
    &*Y(131)+RATE(1519)*D*Y(21)*Y(62)+RATE(1567)*D*Y(24)*Y(162)+RATE(1604)*D&
    &*Y(25)*Y(232)+RATE(1659)*D*Y(33)*Y(46)+RATE(1662)*D*Y(33)*Y(56)&
    &+RATE(1705)*D*Y(42)*Y(131)+RATE(1734)*D*Y(53)*Y(116)+RATE(1763)*D*Y(83)&
    &*Y(131)+RATE(1787)*D*Y(100)*Y(131)+RATE(1881)*D*Y(4)*Y(131)+RATE(2016)*D&
    &*Y(8)*Y(131)+RATE(2055)*D*Y(131)*Y(162)+RATE(2110)*D*Y(63)*Y(131)&
    &+RATE(2131)*D*Y(63)*Y(116)+RATE(2174)*D*Y(10)*Y(116)+RATE(2263)*D*Y(91)&
    &*Y(116)+RATE(2282)*D*Y(116)*Y(145)+RATE(2283)*D*Y(116)*Y(120)+RATE(2285)&
    &*D*Y(116)*Y(176)+RATE(2302)*D*Y(117)*Y(116)+RATE(2325)*D*Y(13)*Y(131)&
    &+RATE(2559)*D*Y(28)*Y(131)+RATE(2572)*D*Y(28)*Y(159)+RATE(2601)*D*Y(104)&
    &*Y(131)+RATE(2661)*D*Y(36)*Y(131)+RATE(2684)*D*Y(36)*Y(116)+RATE(2747)*D&
    &*Y(44)*Y(116)+RATE(2924)*D*Y(48)*Y(131)+RATE(2937)*D*Y(48)*Y(159)&
    &+RATE(3021)*D*Y(57)*Y(131)+RATE(3040)*D*Y(57)*Y(116)
    YDOT(132) = PROD-LOSS
    LOSS = RATE(176)*Y(133)+RATE(177)*Y(133)+RATE(574)*D*Y(133)+RATE(928&
    &)*Y(133)+RATE(929)*Y(133)+RATE(1264)*D*Y(16)*Y(133)+RATE(1265)*D*Y(16)&
    &*Y(133)+RATE(1292)*D*Y(18)*Y(133)+RATE(1366)*D*Y(68)*Y(133)+RATE(1389)*D&
    &*Y(76)*Y(133)+RATE(1400)*D*Y(80)*Y(133)+RATE(1405)*D*Y(81)*Y(133)&
    &+RATE(1482)*D*Y(20)*Y(133)+RATE(1483)*D*Y(20)*Y(133)+RATE(1484)*D*Y(20)&
    &*Y(133)+RATE(1502)*D*Y(21)*Y(133)+RATE(1583)*D*Y(24)*Y(133)+RATE(1584)*D&
    &*Y(24)*Y(133)+RATE(1585)*D*Y(24)*Y(133)+RATE(1600)*D*Y(25)*Y(133)&
    &+RATE(1631)*D*Y(32)*Y(133)+RATE(1645)*D*Y(33)*Y(133)+RATE(1756)*D*Y(82)&
    &*Y(133)+RATE(1757)*D*Y(82)*Y(133)+RATE(1766)*D*Y(83)*Y(133)+RATE(1790)*D&
    &*Y(100)*Y(133)+RATE(1850)*D*Y(2)*Y(133)+RATE(1851)*D*Y(2)*Y(133)&
    &+RATE(1895)*D*Y(4)*Y(133)+RATE(2024)*D*Y(8)*Y(133)+RATE(2049)*D*Y(8)&
    &*Y(133)+RATE(2114)*D*Y(63)*Y(133)+RATE(2189)*D*Y(10)*Y(133)+RATE(2255)*D&
    &*Y(91)*Y(133)+RATE(2290)*D*Y(116)*Y(133)+RATE(2413)*D*Y(13)*Y(133)&
    &+RATE(2414)*D*Y(13)*Y(133)+RATE(2543)*D*Y(27)*Y(133)+RATE(2568)*D*Y(28)&
    &*Y(133)+RATE(2592)*D*Y(28)*Y(133)+RATE(2604)*D*Y(104)*Y(133)+RATE(2650)&
    &*D*Y(35)*Y(133)+RATE(2651)*D*Y(35)*Y(133)+RATE(2664)*D*Y(36)*Y(133)&
    &+RATE(2690)*D*Y(36)*Y(133)+RATE(2723)*D*Y(43)*Y(133)+RATE(2724)*D*Y(43)&
    &*Y(133)+RATE(2731)*D*Y(44)*Y(133)+RATE(2795)*D*Y(55)*Y(133)+RATE(2808)*D&
    &*Y(133)*Y(132)+RATE(2809)*D*Y(133)*Y(184)+RATE(2810)*D*Y(133)*Y(145)&
    &+RATE(2811)*D*Y(133)*Y(174)+RATE(2812)*D*Y(133)*Y(162)+RATE(2813)*D&
    &*Y(133)*Y(165)+RATE(2814)*D*Y(133)*Y(318)+RATE(2815)*D*Y(133)*Y(236)&
    &+RATE(2816)*D*Y(133)*Y(176)+RATE(2817)*D*Y(133)*Y(133)+RATE(2817)*D&
    &*Y(133)*Y(133)+RATE(2818)*D*Y(133)*Y(161)+RATE(2819)*D*Y(133)*Y(221)&
    &+RATE(2820)*D*Y(133)*Y(163)+RATE(2821)*D*Y(133)*Y(163)+RATE(2894)*D*Y(46&
    &)*Y(133)+RATE(3012)*D*Y(56)*Y(133)+RATE(3026)*D*Y(57)*Y(133)+RATE(3044)&
    &*D*Y(57)*Y(133)+RATE(3079)*D*Y(105)*Y(133)
    PROD = RATE(162)*Y(144)+RATE(178)*Y(264)+RATE(263)*Y(129)/safeMantle&
    &+RATE(346)*D*Y(129)/safeMantle*Y(2)+RATE(429)*Y(129)/safeMantle+RATE(768&
    &)*Y(47)*Y(222)+RATE(820)*Y(51)*Y(225)*bulkLayersReciprocal+RATE(910)&
    &*Y(144)+RATE(930)*Y(264)+RATE(1111)*Y(129)+RATE(1194)*Y(138)+RATE(1233)&
    &*D*Y(16)*Y(145)+RATE(1355)*D*Y(67)*Y(145)+RATE(1380)*D*Y(75)*Y(145)&
    &+RATE(1456)*D*Y(20)*Y(145)+RATE(1478)*D*Y(20)*Y(144)+RATE(1562)*D*Y(24)&
    &*Y(145)+RATE(1580)*D*Y(24)*Y(144)+RATE(1582)*D*Y(24)*Y(264)+RATE(1627)*D&
    &*Y(32)*Y(144)+RATE(1690)*D*Y(41)*Y(145)+RATE(1745)*D*Y(82)*Y(145)&
    &+RATE(1754)*D*Y(82)*Y(144)+RATE(1755)*D*Y(82)*Y(264)+RATE(1758)*D*Y(82)&
    &*Y(161)+RATE(1777)*D*Y(99)*Y(145)+RATE(1784)*D*Y(99)*Y(264)+RATE(1842)*D&
    &*Y(2)*Y(144)+RATE(1849)*D*Y(2)*Y(264)+RATE(2057)*D*Y(131)*Y(145)&
    &+RATE(2076)*D*Y(160)*Y(333)+RATE(2095)*D*Y(62)*Y(145)+RATE(2250)*D*Y(90)&
    &*Y(145)+RATE(2282)*D*Y(116)*Y(145)+RATE(2289)*D*Y(116)*Y(144)+RATE(2402)&
    &*D*Y(13)*Y(144)+RATE(2451)*D*Y(92)*Y(145)+RATE(2456)*D*Y(145)*Y(333)&
    &+RATE(2457)*D*Y(145)*Y(232)+RATE(2458)*D*Y(145)*Y(163)+RATE(2482)*D*Y(69&
    &)*Y(134)+RATE(2509)*D*Y(27)*Y(236)+RATE(2528)*D*Y(27)*Y(232)+RATE(2535)&
    &*D*Y(27)*Y(144)+RATE(2541)*D*Y(27)*Y(264)+RATE(2541)*D*Y(27)*Y(264)&
    &+RATE(2545)*D*Y(27)*Y(161)+RATE(2547)*D*Y(27)*Y(56)+RATE(2551)*D*Y(27)&
    &*Y(277)+RATE(2575)*D*Y(28)*Y(159)+RATE(2579)*D*Y(28)*Y(232)+RATE(2594)*D&
    &*Y(28)*Y(161)+RATE(2595)*D*Y(28)*Y(300)+RATE(2599)*D*Y(103)*Y(145)&
    &+RATE(2633)*D*Y(35)*Y(145)+RATE(2649)*D*Y(35)*Y(264)+RATE(2653)*D*Y(35)&
    &*Y(161)+RATE(2654)*D*Y(35)*Y(46)+RATE(2717)*D*Y(43)*Y(145)+RATE(2784)*D&
    &*Y(54)*Y(145)+RATE(2869)*D*Y(46)*Y(82)+RATE(2885)*D*Y(46)*Y(144)&
    &+RATE(2889)*D*Y(46)*Y(103)+RATE(2893)*D*Y(46)*Y(264)+RATE(2895)*D*Y(46)&
    &*Y(265)+RATE(2898)*D*Y(46)*Y(221)+RATE(2960)*D*Y(161)*Y(221)+RATE(2987)&
    &*D*Y(56)*Y(145)+RATE(3009)*D*Y(56)*Y(144)+RATE(3073)*D*Y(105)*Y(134)
    YDOT(133) = PROD-LOSS
    LOSS = RATE(575)*D*Y(134)+RATE(2482)*D*Y(69)*Y(134)+RATE(2822)*D&
    &*Y(134)*Y(333)+RATE(3073)*D*Y(105)*Y(134)
    PROD = RATE(176)*Y(133)+RATE(928)*Y(133)+RATE(1292)*D*Y(18)*Y(133)&
    &+RATE(1366)*D*Y(68)*Y(133)+RATE(1389)*D*Y(76)*Y(133)+RATE(1405)*D*Y(81)&
    &*Y(133)+RATE(1502)*D*Y(21)*Y(133)+RATE(1600)*D*Y(25)*Y(133)+RATE(1645)*D&
    &*Y(33)*Y(133)+RATE(1766)*D*Y(83)*Y(133)+RATE(1773)*D*Y(83)*Y(161)&
    &+RATE(1790)*D*Y(100)*Y(133)+RATE(1895)*D*Y(4)*Y(133)+RATE(1934)*D*Y(4)&
    &*Y(144)+RATE(1936)*D*Y(4)*Y(264)+RATE(2024)*D*Y(8)*Y(133)+RATE(2114)*D&
    &*Y(63)*Y(133)+RATE(2188)*D*Y(10)*Y(264)+RATE(2255)*D*Y(91)*Y(133)&
    &+RATE(2401)*D*Y(13)*Y(144)+RATE(2499)*D*Y(27)*Y(63)+RATE(2504)*D*Y(27)&
    &*Y(162)+RATE(2505)*D*Y(27)*Y(57)+RATE(2508)*D*Y(27)*Y(236)+RATE(2568)*D&
    &*Y(28)*Y(133)+RATE(2574)*D*Y(28)*Y(159)+RATE(2580)*D*Y(28)*Y(99)&
    &+RATE(2582)*D*Y(28)*Y(131)+RATE(2593)*D*Y(28)*Y(161)+RATE(2604)*D*Y(104)&
    &*Y(133)+RATE(2637)*D*Y(35)*Y(48)+RATE(2664)*D*Y(36)*Y(133)+RATE(2675)*D&
    &*Y(36)*Y(232)+RATE(2691)*D*Y(36)*Y(161)+RATE(2731)*D*Y(44)*Y(133)&
    &+RATE(2795)*D*Y(55)*Y(133)+RATE(2808)*D*Y(133)*Y(132)+RATE(2809)*D*Y(133&
    &)*Y(184)+RATE(2810)*D*Y(133)*Y(145)+RATE(2811)*D*Y(133)*Y(174)+RATE(2812&
    &)*D*Y(133)*Y(162)+RATE(2813)*D*Y(133)*Y(165)+RATE(2814)*D*Y(133)*Y(318)&
    &+RATE(2815)*D*Y(133)*Y(236)+RATE(2842)*D*Y(46)*Y(104)+RATE(2846)*D*Y(46)&
    &*Y(266)+RATE(2940)*D*Y(48)*Y(82)+RATE(2946)*D*Y(48)*Y(90)+RATE(2948)*D&
    &*Y(48)*Y(103)+RATE(2949)*D*Y(48)*Y(264)+RATE(3026)*D*Y(57)*Y(133)
    YDOT(134) = PROD-LOSS
    LOSS = RATE(195)*Y(135)+RATE(604)*D*Y(135)+RATE(953)*Y(135)+RATE(954&
    &)*Y(135)+RATE(1300)*D*Y(18)*Y(135)+RATE(1339)*D*Y(18)*Y(135)+RATE(1909)&
    &*D*Y(4)*Y(135)+RATE(1938)*D*Y(4)*Y(135)+RATE(2201)*D*Y(10)*Y(135)&
    &+RATE(2232)*D*Y(66)*Y(135)+RATE(2310)*D*Y(117)*Y(135)+RATE(2434)*D*Y(13)&
    &*Y(135)+RATE(2435)*D*Y(13)*Y(135)+RATE(2910)*D*Y(46)*Y(135)+RATE(2911)*D&
    &*Y(46)*Y(135)
    PROD = RATE(196)*Y(146)+RATE(197)*Y(166)+RATE(279)*Y(137)/safeMantle&
    &+RATE(362)*D*Y(137)/safeMantle*Y(2)+RATE(445)*Y(137)/safeMantle+RATE(773&
    &)*Y(125)*Y(3)+RATE(825)*Y(128)*Y(5)*bulkLayersReciprocal+RATE(955)*Y(146&
    &)+RATE(958)*Y(166)+RATE(1127)*Y(137)+RATE(1210)*Y(140)+RATE(3097)*D&
    &*Y(147)*Y(333)+RATE(3099)*D*Y(167)*Y(333)
    YDOT(135) = PROD-LOSS
    LOSS = RATE(605)*D*Y(136)+RATE(2851)*D*Y(46)*Y(136)+RATE(3092)*D&
    &*Y(136)*Y(333)+RATE(3093)*D*Y(136)*Y(333)+RATE(3094)*D*Y(136)*Y(333)&
    &+RATE(3095)*D*Y(136)*Y(161)+RATE(3096)*D*Y(136)*Y(163)
    PROD = RATE(953)*Y(135)+RATE(1300)*D*Y(18)*Y(135)+RATE(1417)*D*Y(81)&
    &*Y(166)+RATE(1909)*D*Y(4)*Y(135)+RATE(1939)*D*Y(4)*Y(146)+RATE(2005)*D&
    &*Y(6)*Y(106)+RATE(2204)*D*Y(10)*Y(121)+RATE(2233)*D*Y(66)*Y(121)&
    &+RATE(2312)*D*Y(117)*Y(121)+RATE(2437)*D*Y(13)*Y(146)+RATE(3049)*D*Y(57)&
    &*Y(121)
    YDOT(136) = PROD-LOSS
    LOSS = RATE(279)*Y(137)/safeMantle+RATE(362)*D*Y(137)/safeMantle*Y(2&
    &)+RATE(445)*Y(137)/safeMantle+RATE(670)*Y(137)*Y(3)+RATE(774)*Y(137)*Y(3&
    &)+RATE(1044)*Y(137)*totalSwap/safeMantle+RATE(1127)*Y(137)
    PROD = RATE(77)*Y(140)*bulkLayersReciprocal+RATE(604)*D*Y(135)&
    &+RATE(605)*D*Y(136)+RATE(669)*Y(125)*Y(3)
    YDOT(137) = PROD-LOSS
    LOSS = RATE(61)*Y(138)*bulkLayersReciprocal+RATE(1194)*Y(138)
    PROD = RATE(716)*Y(51)*Y(225)*bulkLayersReciprocal+RATE(1028)*Y(129)&
    &*totalSwap/safeMantle
    YDOT(138) = PROD-LOSS
    LOSS = RATE(32)*Y(139)*bulkLayersReciprocal+RATE(693)*Y(5)*Y(139)&
    &*bulkLayersReciprocal+RATE(694)*Y(5)*Y(139)*bulkLayersReciprocal&
    &+RATE(695)*Y(5)*Y(139)*bulkLayersReciprocal+RATE(700)*Y(139)*Y(51)&
    &*bulkLayersReciprocal+RATE(713)*Y(50)*Y(139)*bulkLayersReciprocal&
    &+RATE(714)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(797)*Y(5)*Y(139)&
    &*bulkLayersReciprocal+RATE(798)*Y(5)*Y(139)*bulkLayersReciprocal&
    &+RATE(799)*Y(5)*Y(139)*bulkLayersReciprocal+RATE(804)*Y(139)*Y(51)&
    &*bulkLayersReciprocal+RATE(817)*Y(50)*Y(139)*bulkLayersReciprocal&
    &+RATE(818)*Y(50)*Y(139)*bulkLayersReciprocal+RATE(1165)*Y(139)
    PROD = RATE(691)*Y(5)*Y(126)*bulkLayersReciprocal+RATE(999)*Y(130)&
    &*totalSwap/safeMantle
    YDOT(139) = PROD-LOSS
    LOSS = RATE(77)*Y(140)*bulkLayersReciprocal+RATE(722)*Y(140)*Y(5)&
    &*bulkLayersReciprocal+RATE(826)*Y(140)*Y(5)*bulkLayersReciprocal&
    &+RATE(1210)*Y(140)
    PROD = RATE(721)*Y(128)*Y(5)*bulkLayersReciprocal+RATE(1044)*Y(137)&
    &*totalSwap/safeMantle
    YDOT(140) = PROD-LOSS
    LOSS = RATE(523)*D*Y(141)
    PROD = RATE(240)*Y(143)/safeMantle+RATE(323)*D*Y(143)/safeMantle*Y(2&
    &)+RATE(406)*Y(143)/safeMantle+RATE(746)*Y(3)*Y(130)+RATE(798)*Y(5)*Y(139&
    &)*bulkLayersReciprocal+RATE(1088)*Y(143)+RATE(1171)*Y(152)
    YDOT(141) = PROD-LOSS
    LOSS = RATE(524)*D*Y(142)+RATE(1450)*D*Y(20)*Y(142)+RATE(2091)*D&
    &*Y(62)*Y(142)+RATE(2207)*D*Y(142)*Y(333)+RATE(2208)*D*Y(142)*Y(333)&
    &+RATE(2209)*D*Y(142)*Y(333)+RATE(2210)*D*Y(142)*Y(333)+RATE(2211)*D&
    &*Y(142)*Y(333)+RATE(2212)*D*Y(142)*Y(159)+RATE(2213)*D*Y(142)*Y(183)&
    &+RATE(2247)*D*Y(90)*Y(142)+RATE(2448)*D*Y(92)*Y(142)+RATE(2711)*D*Y(43)&
    &*Y(142)+RATE(2775)*D*Y(54)*Y(142)
    PROD = RATE(864)*Y(159)+RATE(1308)*D*Y(18)*Y(159)+RATE(1393)*D*Y(76)&
    &*Y(131)+RATE(1511)*D*Y(21)*Y(159)+RATE(1517)*D*Y(21)*Y(131)+RATE(1606)*D&
    &*Y(25)*Y(62)+RATE(1650)*D*Y(33)*Y(237)+RATE(1652)*D*Y(33)*Y(159)&
    &+RATE(1658)*D*Y(33)*Y(161)+RATE(1686)*D*Y(41)*Y(132)+RATE(1716)*D*Y(42)&
    &*Y(131)+RATE(1730)*D*Y(53)*Y(131)+RATE(1921)*D*Y(4)*Y(159)+RATE(2056)*D&
    &*Y(131)*Y(189)+RATE(2057)*D*Y(131)*Y(145)+RATE(2059)*D*Y(131)*Y(176)&
    &+RATE(2068)*D*Y(132)*Y(131)+RATE(2125)*D*Y(63)*Y(131)+RATE(2169)*D*Y(10)&
    &*Y(131)+RATE(2225)*D*Y(66)*Y(131)+RATE(2261)*D*Y(91)*Y(131)+RATE(2272)*D&
    &*Y(102)*Y(131)+RATE(2273)*D*Y(102)*Y(131)+RATE(2281)*D*Y(116)*Y(132)&
    &+RATE(2299)*D*Y(117)*Y(131)+RATE(2573)*D*Y(28)*Y(159)+RATE(2618)*D*Y(120&
    &)*Y(131)+RATE(2629)*D*Y(35)*Y(132)+RATE(2677)*D*Y(36)*Y(131)+RATE(2737)&
    &*D*Y(44)*Y(131)+RATE(2831)*D*Y(46)*Y(53)+RATE(2938)*D*Y(48)*Y(159)&
    &+RATE(2969)*D*Y(162)*Y(159)+RATE(3035)*D*Y(57)*Y(131)+RATE(3112)*D*Y(278&
    &)*Y(109)
    YDOT(142) = PROD-LOSS
    LOSS = RATE(240)*Y(143)/safeMantle+RATE(323)*D*Y(143)/safeMantle*Y(2&
    &)+RATE(406)*Y(143)/safeMantle+RATE(645)*Y(3)*Y(143)+RATE(649)*Y(143)&
    &*Y(226)+RATE(749)*Y(3)*Y(143)+RATE(753)*Y(143)*Y(226)+RATE(1005)*Y(143)&
    &*totalSwap/safeMantle+RATE(1088)*Y(143)
    PROD = RATE(38)*Y(152)*bulkLayersReciprocal+RATE(86)*Y(157)+RATE(523&
    &)*D*Y(141)+RATE(524)*D*Y(142)+RATE(642)*Y(3)*Y(130)
    YDOT(143) = PROD-LOSS
    LOSS = RATE(162)*Y(144)+RATE(546)*D*Y(144)+RATE(910)*Y(144)&
    &+RATE(1478)*D*Y(20)*Y(144)+RATE(1580)*D*Y(24)*Y(144)+RATE(1627)*D*Y(32)&
    &*Y(144)+RATE(1754)*D*Y(82)*Y(144)+RATE(1783)*D*Y(99)*Y(144)+RATE(1841)*D&
    &*Y(2)*Y(144)+RATE(1842)*D*Y(2)*Y(144)+RATE(1843)*D*Y(2)*Y(144)+RATE(1934&
    &)*D*Y(4)*Y(144)+RATE(2180)*D*Y(10)*Y(144)+RATE(2289)*D*Y(116)*Y(144)&
    &+RATE(2401)*D*Y(13)*Y(144)+RATE(2402)*D*Y(13)*Y(144)+RATE(2535)*D*Y(27)&
    &*Y(144)+RATE(2884)*D*Y(46)*Y(144)+RATE(2885)*D*Y(46)*Y(144)+RATE(2886)*D&
    &*Y(46)*Y(144)+RATE(3009)*D*Y(56)*Y(144)
    PROD = RATE(251)*Y(148)/safeMantle+RATE(334)*D*Y(148)/safeMantle*Y(2&
    &)+RATE(417)*Y(148)/safeMantle+RATE(1099)*Y(148)+RATE(1182)*Y(153)&
    &+RATE(1630)*D*Y(32)*Y(264)+RATE(2075)*D*Y(160)*Y(333)+RATE(2290)*D*Y(116&
    &)*Y(133)+RATE(2649)*D*Y(35)*Y(264)+RATE(2652)*D*Y(35)*Y(161)+RATE(2657)&
    &*D*Y(35)*Y(56)+RATE(2810)*D*Y(133)*Y(145)+RATE(2890)*D*Y(46)*Y(43)
    YDOT(144) = PROD-LOSS
    LOSS = RATE(547)*D*Y(145)+RATE(1233)*D*Y(16)*Y(145)+RATE(1355)*D&
    &*Y(67)*Y(145)+RATE(1380)*D*Y(75)*Y(145)+RATE(1456)*D*Y(20)*Y(145)&
    &+RATE(1562)*D*Y(24)*Y(145)+RATE(1690)*D*Y(41)*Y(145)+RATE(1745)*D*Y(82)&
    &*Y(145)+RATE(1777)*D*Y(99)*Y(145)+RATE(2057)*D*Y(131)*Y(145)+RATE(2095)&
    &*D*Y(62)*Y(145)+RATE(2250)*D*Y(90)*Y(145)+RATE(2282)*D*Y(116)*Y(145)&
    &+RATE(2451)*D*Y(92)*Y(145)+RATE(2456)*D*Y(145)*Y(333)+RATE(2457)*D*Y(145&
    &)*Y(232)+RATE(2458)*D*Y(145)*Y(163)+RATE(2599)*D*Y(103)*Y(145)+RATE(2633&
    &)*D*Y(35)*Y(145)+RATE(2717)*D*Y(43)*Y(145)+RATE(2784)*D*Y(54)*Y(145)&
    &+RATE(2810)*D*Y(133)*Y(145)+RATE(2987)*D*Y(56)*Y(145)
    PROD = RATE(2049)*D*Y(8)*Y(133)+RATE(2189)*D*Y(10)*Y(133)+RATE(2498)&
    &*D*Y(27)*Y(63)+RATE(2638)*D*Y(35)*Y(162)+RATE(2674)*D*Y(36)*Y(232)&
    &+RATE(2680)*D*Y(36)*Y(62)+RATE(2752)*D*Y(44)*Y(161)+RATE(2816)*D*Y(133)&
    &*Y(176)+RATE(2844)*D*Y(46)*Y(44)+RATE(2845)*D*Y(46)*Y(55)+RATE(3044)*D&
    &*Y(57)*Y(133)
    YDOT(145) = PROD-LOSS
    LOSS = RATE(196)*Y(146)+RATE(606)*D*Y(146)+RATE(955)*Y(146)+RATE(956&
    &)*Y(146)+RATE(957)*Y(146)+RATE(1301)*D*Y(18)*Y(146)+RATE(1910)*D*Y(4)&
    &*Y(146)+RATE(1939)*D*Y(4)*Y(146)+RATE(2202)*D*Y(10)*Y(146)+RATE(2436)*D&
    &*Y(13)*Y(146)+RATE(2437)*D*Y(13)*Y(146)+RATE(2912)*D*Y(46)*Y(146)
    PROD = RATE(280)*Y(151)/safeMantle+RATE(363)*D*Y(151)/safeMantle*Y(2&
    &)+RATE(446)*Y(151)/safeMantle+RATE(774)*Y(137)*Y(3)+RATE(826)*Y(140)*Y(5&
    &)*bulkLayersReciprocal+RATE(959)*Y(166)+RATE(1128)*Y(151)+RATE(1211)&
    &*Y(155)+RATE(1761)*D*Y(82)*Y(166)+RATE(1781)*D*Y(99)*Y(167)+RATE(2105)*D&
    &*Y(62)*Y(167)+RATE(2913)*D*Y(46)*Y(166)+RATE(3100)*D*Y(167)*Y(333)&
    &+RATE(3101)*D*Y(177)*Y(333)
    YDOT(146) = PROD-LOSS
    LOSS = RATE(607)*D*Y(147)+RATE(2007)*D*Y(6)*Y(147)+RATE(2852)*D*Y(46&
    &)*Y(147)+RATE(3097)*D*Y(147)*Y(333)+RATE(3098)*D*Y(147)*Y(333)
    PROD = RATE(956)*Y(146)+RATE(1301)*D*Y(18)*Y(146)+RATE(1418)*D*Y(81)&
    &*Y(166)+RATE(1665)*D*Y(33)*Y(166)+RATE(1739)*D*Y(53)*Y(166)+RATE(1910)*D&
    &*Y(4)*Y(146)+RATE(1940)*D*Y(4)*Y(166)+RATE(2006)*D*Y(6)*Y(122)+RATE(2201&
    &)*D*Y(10)*Y(135)+RATE(2232)*D*Y(66)*Y(135)+RATE(2310)*D*Y(117)*Y(135)
    YDOT(147) = PROD-LOSS
    LOSS = RATE(96)*Y(148)+RATE(251)*Y(148)/safeMantle+RATE(334)*D*Y(148&
    &)/safeMantle*Y(2)+RATE(417)*Y(148)/safeMantle+RATE(1016)*Y(148)&
    &*totalSwap/safeMantle+RATE(1099)*Y(148)
    PROD = RATE(49)*Y(153)*bulkLayersReciprocal+RATE(514)*D*Y(160)&
    &+RATE(546)*D*Y(144)+RATE(547)*D*Y(145)
    YDOT(148) = PROD-LOSS
    LOSS = RATE(481)*D*Y(149)
    PROD = RATE(218)*Y(150)/safeMantle+RATE(301)*D*Y(150)/safeMantle*Y(2&
    &)+RATE(384)*Y(150)/safeMantle+RATE(745)*Y(3)*Y(130)+RATE(797)*Y(5)*Y(139&
    &)*bulkLayersReciprocal+RATE(1066)*Y(150)+RATE(1149)*Y(154)
    YDOT(149) = PROD-LOSS
    LOSS = RATE(218)*Y(150)/safeMantle+RATE(301)*D*Y(150)/safeMantle*Y(2&
    &)+RATE(384)*Y(150)/safeMantle+RATE(644)*Y(3)*Y(150)+RATE(748)*Y(3)*Y(150&
    &)+RATE(983)*Y(150)*totalSwap/safeMantle+RATE(1066)*Y(150)
    PROD = RATE(16)*Y(154)*bulkLayersReciprocal+RATE(85)*Y(157)+RATE(481&
    &)*D*Y(149)+RATE(641)*Y(3)*Y(130)
    YDOT(150) = PROD-LOSS
    LOSS = RATE(280)*Y(151)/safeMantle+RATE(363)*D*Y(151)/safeMantle*Y(2&
    &)+RATE(446)*Y(151)/safeMantle+RATE(671)*Y(151)*Y(3)+RATE(775)*Y(151)*Y(3&
    &)+RATE(1045)*Y(151)*totalSwap/safeMantle+RATE(1128)*Y(151)
    PROD = RATE(78)*Y(155)*bulkLayersReciprocal+RATE(606)*D*Y(146)&
    &+RATE(607)*D*Y(147)+RATE(670)*Y(137)*Y(3)
    YDOT(151) = PROD-LOSS
    LOSS = RATE(38)*Y(152)*bulkLayersReciprocal+RATE(697)*Y(5)*Y(152)&
    &*bulkLayersReciprocal+RATE(701)*Y(152)*Y(228)*bulkLayersReciprocal&
    &+RATE(801)*Y(5)*Y(152)*bulkLayersReciprocal+RATE(805)*Y(152)*Y(228)&
    &*bulkLayersReciprocal+RATE(1171)*Y(152)
    PROD = RATE(694)*Y(5)*Y(139)*bulkLayersReciprocal+RATE(1005)*Y(143)&
    &*totalSwap/safeMantle
    YDOT(152) = PROD-LOSS
    LOSS = RATE(49)*Y(153)*bulkLayersReciprocal+RATE(1182)*Y(153)
    PROD = RATE(1016)*Y(148)*totalSwap/safeMantle
    YDOT(153) = PROD-LOSS
    LOSS = RATE(16)*Y(154)*bulkLayersReciprocal+RATE(696)*Y(5)*Y(154)&
    &*bulkLayersReciprocal+RATE(800)*Y(5)*Y(154)*bulkLayersReciprocal&
    &+RATE(1149)*Y(154)
    PROD = RATE(693)*Y(5)*Y(139)*bulkLayersReciprocal+RATE(983)*Y(150)&
    &*totalSwap/safeMantle
    YDOT(154) = PROD-LOSS
    LOSS = RATE(78)*Y(155)*bulkLayersReciprocal+RATE(723)*Y(155)*Y(5)&
    &*bulkLayersReciprocal+RATE(827)*Y(155)*Y(5)*bulkLayersReciprocal&
    &+RATE(1211)*Y(155)
    PROD = RATE(722)*Y(140)*Y(5)*bulkLayersReciprocal+RATE(1045)*Y(151)&
    &*totalSwap/safeMantle
    YDOT(155) = PROD-LOSS
    LOSS = RATE(267)*Y(156)/safeMantle+RATE(350)*D*Y(156)/safeMantle*Y(2&
    &)+RATE(433)*Y(156)/safeMantle+RATE(1032)*Y(156)*totalSwap/safeMantle&
    &+RATE(1115)*Y(156)
    PROD = RATE(65)*Y(168)*bulkLayersReciprocal+RATE(581)*D*Y(161)&
    &+RATE(582)*D*Y(162)+RATE(663)*Y(47)*Y(47)
    YDOT(156) = PROD-LOSS
    LOSS = RATE(85)*Y(157)+RATE(86)*Y(157)+RATE(87)*Y(157)+RATE(224)&
    &*Y(157)/safeMantle+RATE(307)*D*Y(157)/safeMantle*Y(2)+RATE(390)*Y(157&
    &)/safeMantle+RATE(989)*Y(157)*totalSwap/safeMantle+RATE(1072)*Y(157)
    PROD = RATE(22)*Y(169)*bulkLayersReciprocal+RATE(489)*D*Y(159)&
    &+RATE(490)*D*Y(172)+RATE(644)*Y(3)*Y(150)+RATE(645)*Y(3)*Y(143)+RATE(649&
    &)*Y(143)*Y(226)
    YDOT(157) = PROD-LOSS
    LOSS = RATE(281)*Y(158)/safeMantle+RATE(364)*D*Y(158)/safeMantle*Y(2&
    &)+RATE(447)*Y(158)/safeMantle+RATE(1046)*Y(158)*totalSwap/safeMantle&
    &+RATE(1129)*Y(158)
    PROD = RATE(79)*Y(170)*bulkLayersReciprocal+RATE(608)*D*Y(166)&
    &+RATE(609)*D*Y(167)+RATE(610)*D*Y(177)+RATE(671)*Y(151)*Y(3)
    YDOT(158) = PROD-LOSS
    LOSS = RATE(134)*Y(159)+RATE(135)*Y(159)+RATE(489)*D*Y(159)+RATE(863&
    &)*Y(159)+RATE(864)*Y(159)+RATE(865)*Y(159)+RATE(1308)*D*Y(18)*Y(159)&
    &+RATE(1309)*D*Y(18)*Y(159)+RATE(1473)*D*Y(20)*Y(159)+RATE(1509)*D*Y(21)&
    &*Y(159)+RATE(1510)*D*Y(21)*Y(159)+RATE(1511)*D*Y(21)*Y(159)+RATE(1652)*D&
    &*Y(33)*Y(159)+RATE(1672)*D*Y(159)*Y(318)+RATE(1712)*D*Y(42)*Y(159)&
    &+RATE(1920)*D*Y(4)*Y(159)+RATE(1921)*D*Y(4)*Y(159)+RATE(1922)*D*Y(4)&
    &*Y(159)+RATE(2067)*D*Y(132)*Y(159)+RATE(2159)*D*Y(10)*Y(159)+RATE(2160)&
    &*D*Y(10)*Y(159)+RATE(2212)*D*Y(142)*Y(159)+RATE(2223)*D*Y(66)*Y(159)&
    &+RATE(2297)*D*Y(117)*Y(159)+RATE(2357)*D*Y(13)*Y(159)+RATE(2358)*D*Y(13)&
    &*Y(159)+RATE(2572)*D*Y(28)*Y(159)+RATE(2573)*D*Y(28)*Y(159)+RATE(2574)*D&
    &*Y(28)*Y(159)+RATE(2575)*D*Y(28)*Y(159)+RATE(2937)*D*Y(48)*Y(159)&
    &+RATE(2938)*D*Y(48)*Y(159)+RATE(2969)*D*Y(162)*Y(159)+RATE(3081)*D*Y(106&
    &)*Y(159)
    PROD = RATE(224)*Y(157)/safeMantle+RATE(307)*D*Y(157)/safeMantle*Y(2&
    &)+RATE(390)*Y(157)/safeMantle+RATE(748)*Y(3)*Y(150)+RATE(749)*Y(3)*Y(143&
    &)+RATE(753)*Y(143)*Y(226)+RATE(800)*Y(5)*Y(154)*bulkLayersReciprocal&
    &+RATE(801)*Y(5)*Y(152)*bulkLayersReciprocal+RATE(805)*Y(152)*Y(228)&
    &*bulkLayersReciprocal+RATE(1072)*Y(157)+RATE(1155)*Y(169)+RATE(1676)*D&
    &*Y(172)*Y(333)+RATE(1678)*D*Y(172)*Y(54)+RATE(2155)*D*Y(10)*Y(237)
    YDOT(159) = PROD-LOSS
    LOSS = RATE(514)*D*Y(160)+RATE(2075)*D*Y(160)*Y(333)+RATE(2076)*D&
    &*Y(160)*Y(333)
    PROD = RATE(2180)*D*Y(10)*Y(144)+RATE(2751)*D*Y(44)*Y(161)
    YDOT(160) = PROD-LOSS
    LOSS = RATE(181)*Y(161)+RATE(182)*Y(161)+RATE(581)*D*Y(161)+RATE(932&
    &)*Y(161)+RATE(933)*Y(161)+RATE(1268)*D*Y(16)*Y(161)+RATE(1329)*D*Y(18)&
    &*Y(161)+RATE(1330)*D*Y(18)*Y(161)+RATE(1363)*D*Y(67)*Y(161)+RATE(1371)*D&
    &*Y(68)*Y(161)+RATE(1387)*D*Y(75)*Y(161)+RATE(1419)*D*Y(89)*Y(161)&
    &+RATE(1420)*D*Y(89)*Y(161)+RATE(1485)*D*Y(20)*Y(161)+RATE(1486)*D*Y(20)&
    &*Y(161)+RATE(1487)*D*Y(20)*Y(161)+RATE(1488)*D*Y(20)*Y(161)+RATE(1533)*D&
    &*Y(21)*Y(161)+RATE(1534)*D*Y(21)*Y(161)+RATE(1535)*D*Y(21)*Y(161)&
    &+RATE(1586)*D*Y(24)*Y(161)+RATE(1587)*D*Y(24)*Y(161)+RATE(1588)*D*Y(24)&
    &*Y(161)+RATE(1589)*D*Y(24)*Y(161)+RATE(1590)*D*Y(24)*Y(161)+RATE(1612)*D&
    &*Y(25)*Y(161)+RATE(1632)*D*Y(32)*Y(161)+RATE(1633)*D*Y(32)*Y(161)&
    &+RATE(1634)*D*Y(32)*Y(161)+RATE(1658)*D*Y(33)*Y(161)+RATE(1701)*D*Y(41)&
    &*Y(161)+RATE(1708)*D*Y(42)*Y(161)+RATE(1758)*D*Y(82)*Y(161)+RATE(1759)*D&
    &*Y(82)*Y(161)+RATE(1767)*D*Y(83)*Y(161)+RATE(1773)*D*Y(83)*Y(161)&
    &+RATE(1785)*D*Y(99)*Y(161)+RATE(1791)*D*Y(100)*Y(161)+RATE(1801)*D*Y(2)&
    &*Y(161)+RATE(1854)*D*Y(2)*Y(161)+RATE(1897)*D*Y(4)*Y(161)+RATE(1949)*D&
    &*Y(6)*Y(161)+RATE(1994)*D*Y(6)*Y(161)+RATE(1995)*D*Y(6)*Y(161)+RATE(2025&
    &)*D*Y(8)*Y(161)+RATE(2050)*D*Y(8)*Y(161)+RATE(2069)*D*Y(132)*Y(161)&
    &+RATE(2115)*D*Y(63)*Y(161)+RATE(2191)*D*Y(10)*Y(161)+RATE(2256)*D*Y(91)&
    &*Y(161)+RATE(2291)*D*Y(116)*Y(161)+RATE(2292)*D*Y(116)*Y(161)+RATE(2330)&
    &*D*Y(13)*Y(161)+RATE(2417)*D*Y(13)*Y(161)+RATE(2545)*D*Y(27)*Y(161)&
    &+RATE(2569)*D*Y(28)*Y(161)+RATE(2593)*D*Y(28)*Y(161)+RATE(2594)*D*Y(28)&
    &*Y(161)+RATE(2605)*D*Y(104)*Y(161)+RATE(2652)*D*Y(35)*Y(161)+RATE(2653)&
    &*D*Y(35)*Y(161)+RATE(2665)*D*Y(36)*Y(161)+RATE(2691)*D*Y(36)*Y(161)&
    &+RATE(2692)*D*Y(36)*Y(161)+RATE(2751)*D*Y(44)*Y(161)+RATE(2752)*D*Y(44)&
    &*Y(161)+RATE(2818)*D*Y(133)*Y(161)+RATE(2930)*D*Y(48)*Y(161)+RATE(2954)&
    &*D*Y(161)*Y(187)+RATE(2955)*D*Y(161)*Y(319)+RATE(2956)*D*Y(161)*Y(234)&
    &+RATE(2957)*D*Y(161)*Y(165)+RATE(2958)*D*Y(161)*Y(303)+RATE(2959)*D&
    &*Y(161)*Y(303)+RATE(2960)*D*Y(161)*Y(221)+RATE(2961)*D*Y(161)*Y(221)&
    &+RATE(2962)*D*Y(161)*Y(163)+RATE(2963)*D*Y(161)*Y(277)+RATE(3027)*D*Y(57&
    &)*Y(161)+RATE(3080)*D*Y(105)*Y(161)+RATE(3095)*D*Y(136)*Y(161)
    PROD = RATE(183)*Y(175)+RATE(267)*Y(156)/safeMantle+RATE(350)*D&
    &*Y(156)/safeMantle*Y(2)+RATE(433)*Y(156)/safeMantle+RATE(767)*Y(47)*Y(47&
    &)+RATE(819)*Y(51)*Y(51)*bulkLayersReciprocal+RATE(935)*Y(175)+RATE(1115)&
    &*Y(156)+RATE(1198)*Y(168)+RATE(1221)*D*Y(16)*Y(162)+RATE(1238)*D*Y(16)&
    &*Y(176)+RATE(1351)*D*Y(67)*Y(162)+RATE(1358)*D*Y(67)*Y(176)+RATE(1382)*D&
    &*Y(75)*Y(176)+RATE(1441)*D*Y(20)*Y(162)+RATE(1465)*D*Y(20)*Y(176)&
    &+RATE(1490)*D*Y(20)*Y(175)+RATE(1550)*D*Y(24)*Y(162)+RATE(1568)*D*Y(24)&
    &*Y(176)+RATE(1635)*D*Y(32)*Y(175)+RATE(1746)*D*Y(82)*Y(176)+RATE(1779)*D&
    &*Y(99)*Y(176)+RATE(1856)*D*Y(2)*Y(175)+RATE(1978)*D*Y(6)*Y(176)&
    &+RATE(2055)*D*Y(131)*Y(162)+RATE(2058)*D*Y(131)*Y(162)+RATE(2059)*D&
    &*Y(131)*Y(176)+RATE(2102)*D*Y(62)*Y(176)+RATE(2254)*D*Y(90)*Y(176)&
    &+RATE(2278)*D*Y(116)*Y(162)+RATE(2285)*D*Y(116)*Y(176)+RATE(2293)*D&
    &*Y(116)*Y(175)+RATE(2369)*D*Y(13)*Y(232)+RATE(2426)*D*Y(13)*Y(320)&
    &+RATE(2454)*D*Y(92)*Y(176)+RATE(2483)*D*Y(69)*Y(162)+RATE(2542)*D*Y(27)&
    &*Y(264)+RATE(2546)*D*Y(27)*Y(175)+RATE(2600)*D*Y(103)*Y(176)+RATE(2639)&
    &*D*Y(35)*Y(176)+RATE(2702)*D*Y(43)*Y(162)+RATE(2720)*D*Y(43)*Y(176)&
    &+RATE(2763)*D*Y(54)*Y(162)+RATE(2789)*D*Y(54)*Y(176)+RATE(2812)*D*Y(133)&
    &*Y(162)+RATE(2816)*D*Y(133)*Y(176)+RATE(2817)*D*Y(133)*Y(133)+RATE(2837)&
    &*D*Y(46)*Y(243)+RATE(2847)*D*Y(46)*Y(176)+RATE(2853)*D*Y(46)*Y(236)&
    &+RATE(2870)*D*Y(46)*Y(232)+RATE(2886)*D*Y(46)*Y(144)+RATE(2893)*D*Y(46)&
    &*Y(264)+RATE(2894)*D*Y(46)*Y(133)+RATE(2897)*D*Y(46)*Y(175)+RATE(2899)*D&
    &*Y(46)*Y(221)+RATE(2902)*D*Y(46)*Y(56)+RATE(2904)*D*Y(46)*Y(320)&
    &+RATE(2905)*D*Y(46)*Y(277)+RATE(2915)*D*Y(46)*Y(46)+RATE(2949)*D*Y(48)&
    &*Y(264)+RATE(2952)*D*Y(48)*Y(320)+RATE(2964)*D*Y(162)*Y(80)+RATE(2965)*D&
    &*Y(162)*Y(183)+RATE(2966)*D*Y(162)*Y(163)+RATE(2969)*D*Y(162)*Y(159)&
    &+RATE(2972)*D*Y(176)*Y(333)+RATE(2973)*D*Y(176)*Y(232)+RATE(2990)*D*Y(56&
    &)*Y(176)+RATE(3013)*D*Y(56)*Y(175)+RATE(3055)*D*Y(163)*Y(176)+RATE(3074)&
    &*D*Y(105)*Y(162)
    YDOT(161) = PROD-LOSS
    LOSS = RATE(582)*D*Y(162)+RATE(934)*Y(162)+RATE(1221)*D*Y(16)*Y(162)&
    &+RATE(1237)*D*Y(16)*Y(162)+RATE(1351)*D*Y(67)*Y(162)+RATE(1357)*D*Y(67)&
    &*Y(162)+RATE(1441)*D*Y(20)*Y(162)+RATE(1464)*D*Y(20)*Y(162)+RATE(1550)*D&
    &*Y(24)*Y(162)+RATE(1567)*D*Y(24)*Y(162)+RATE(2055)*D*Y(131)*Y(162)&
    &+RATE(2058)*D*Y(131)*Y(162)+RATE(2278)*D*Y(116)*Y(162)+RATE(2284)*D&
    &*Y(116)*Y(162)+RATE(2483)*D*Y(69)*Y(162)+RATE(2504)*D*Y(27)*Y(162)&
    &+RATE(2638)*D*Y(35)*Y(162)+RATE(2702)*D*Y(43)*Y(162)+RATE(2763)*D*Y(54)&
    &*Y(162)+RATE(2812)*D*Y(133)*Y(162)+RATE(2964)*D*Y(162)*Y(80)+RATE(2965)&
    &*D*Y(162)*Y(183)+RATE(2966)*D*Y(162)*Y(163)+RATE(2967)*D*Y(162)*Y(333)&
    &+RATE(2968)*D*Y(162)*Y(80)+RATE(2969)*D*Y(162)*Y(159)+RATE(2970)*D*Y(162&
    &)*Y(267)+RATE(2971)*D*Y(162)*Y(163)+RATE(3074)*D*Y(105)*Y(162)
    PROD = RATE(181)*Y(161)+RATE(932)*Y(161)+RATE(1708)*D*Y(42)*Y(161)&
    &+RATE(1767)*D*Y(83)*Y(161)+RATE(1791)*D*Y(100)*Y(161)+RATE(1897)*D*Y(4)&
    &*Y(161)+RATE(2025)*D*Y(8)*Y(161)+RATE(2115)*D*Y(63)*Y(161)+RATE(2256)*D&
    &*Y(91)*Y(161)+RATE(2330)*D*Y(13)*Y(161)+RATE(2368)*D*Y(13)*Y(232)&
    &+RATE(2569)*D*Y(28)*Y(161)+RATE(2605)*D*Y(104)*Y(161)+RATE(2665)*D*Y(36)&
    &*Y(161)+RATE(2834)*D*Y(46)*Y(63)+RATE(2848)*D*Y(46)*Y(57)+RATE(2930)*D&
    &*Y(48)*Y(161)+RATE(2941)*D*Y(48)*Y(232)+RATE(2951)*D*Y(48)*Y(56)&
    &+RATE(2954)*D*Y(161)*Y(187)+RATE(2955)*D*Y(161)*Y(319)+RATE(3027)*D*Y(57&
    &)*Y(161)
    YDOT(162) = PROD-LOSS
    LOSS = RATE(188)*Y(163)+RATE(590)*D*Y(163)+RATE(943)*Y(163)&
    &+RATE(1283)*D*Y(16)*Y(163)+RATE(1346)*D*Y(18)*Y(163)+RATE(1364)*D*Y(67)&
    &*Y(163)+RATE(1367)*D*Y(68)*Y(163)+RATE(1372)*D*Y(68)*Y(163)+RATE(1390)*D&
    &*Y(76)*Y(163)+RATE(1495)*D*Y(20)*Y(163)+RATE(1496)*D*Y(20)*Y(163)&
    &+RATE(1503)*D*Y(21)*Y(163)+RATE(1540)*D*Y(21)*Y(163)+RATE(1541)*D*Y(21)&
    &*Y(163)+RATE(1598)*D*Y(24)*Y(163)+RATE(1599)*D*Y(24)*Y(163)+RATE(1616)*D&
    &*Y(25)*Y(163)+RATE(1641)*D*Y(32)*Y(163)+RATE(1663)*D*Y(33)*Y(163)&
    &+RATE(1703)*D*Y(41)*Y(163)+RATE(1738)*D*Y(53)*Y(163)+RATE(1760)*D*Y(82)&
    &*Y(163)+RATE(1768)*D*Y(83)*Y(163)+RATE(1792)*D*Y(100)*Y(163)+RATE(1902)&
    &*D*Y(4)*Y(163)+RATE(1998)*D*Y(6)*Y(163)+RATE(2062)*D*Y(132)*Y(163)&
    &+RATE(2070)*D*Y(132)*Y(163)+RATE(2117)*D*Y(63)*Y(163)+RATE(2133)*D*Y(63)&
    &*Y(163)+RATE(2134)*D*Y(63)*Y(163)+RATE(2197)*D*Y(10)*Y(163)+RATE(2257)*D&
    &*Y(91)*Y(163)+RATE(2266)*D*Y(91)*Y(163)+RATE(2308)*D*Y(117)*Y(163)&
    &+RATE(2458)*D*Y(145)*Y(163)+RATE(2607)*D*Y(104)*Y(163)+RATE(2619)*D&
    &*Y(120)*Y(163)+RATE(2659)*D*Y(35)*Y(163)+RATE(2660)*D*Y(35)*Y(163)&
    &+RATE(2666)*D*Y(36)*Y(163)+RATE(2695)*D*Y(36)*Y(163)+RATE(2696)*D*Y(36)&
    &*Y(163)+RATE(2732)*D*Y(44)*Y(163)+RATE(2753)*D*Y(44)*Y(163)+RATE(2754)*D&
    &*Y(44)*Y(163)+RATE(2820)*D*Y(133)*Y(163)+RATE(2821)*D*Y(133)*Y(163)&
    &+RATE(2962)*D*Y(161)*Y(163)+RATE(2966)*D*Y(162)*Y(163)+RATE(2971)*D&
    &*Y(162)*Y(163)+RATE(3015)*D*Y(56)*Y(163)+RATE(3028)*D*Y(57)*Y(163)&
    &+RATE(3046)*D*Y(57)*Y(163)+RATE(3047)*D*Y(57)*Y(163)+RATE(3051)*D*Y(163)&
    &*Y(18)+RATE(3052)*D*Y(163)*Y(184)+RATE(3053)*D*Y(163)*Y(174)+RATE(3054)&
    &*D*Y(163)*Y(189)+RATE(3055)*D*Y(163)*Y(176)+RATE(3056)*D*Y(163)*Y(236)&
    &+RATE(3057)*D*Y(163)*Y(116)+RATE(3058)*D*Y(163)*Y(116)+RATE(3059)*D&
    &*Y(163)*Y(173)+RATE(3060)*D*Y(163)*Y(320)+RATE(3061)*D*Y(163)*Y(277)&
    &+RATE(3096)*D*Y(136)*Y(163)
    PROD = RATE(142)*Y(233)+RATE(149)*Y(183)+RATE(163)*Y(173)+RATE(164)&
    &*Y(324)+RATE(179)*Y(265)+RATE(186)*Y(300)+RATE(189)*Y(316)+RATE(189)&
    &*Y(316)+RATE(199)*Y(302)+RATE(200)*Y(277)+RATE(272)*Y(164)/safeMantle&
    &+RATE(355)*D*Y(164)/safeMantle*Y(2)+RATE(438)*Y(164)/safeMantle+RATE(878&
    &)*Y(233)+RATE(893)*Y(183)+RATE(911)*Y(173)+RATE(913)*Y(174)+RATE(915)&
    &*Y(324)+RATE(931)*Y(265)+RATE(939)*Y(300)+RATE(945)*Y(316)+RATE(945)&
    &*Y(316)+RATE(964)*Y(302)+RATE(965)*Y(277)+RATE(1120)*Y(164)+RATE(1203)&
    &*Y(171)+RATE(1252)*D*Y(16)*Y(233)+RATE(1256)*D*Y(16)*Y(173)+RATE(1267)*D&
    &*Y(16)*Y(265)+RATE(1273)*D*Y(16)*Y(316)+RATE(1276)*D*Y(16)*Y(277)&
    &+RATE(1337)*D*Y(18)*Y(277)+RATE(1342)*D*Y(18)*Y(302)+RATE(1457)*D*Y(20)&
    &*Y(174)+RATE(1797)*D*Y(234)*Y(333)+RATE(1844)*D*Y(2)*Y(173)+RATE(1853)*D&
    &*Y(2)*Y(265)+RATE(1863)*D*Y(2)*Y(316)+RATE(1865)*D*Y(2)*Y(277)+RATE(2097&
    &)*D*Y(62)*Y(174)+RATE(2142)*D*Y(184)*Y(333)+RATE(2238)*D*Y(189)*Y(333)&
    &+RATE(2251)*D*Y(90)*Y(174)+RATE(2279)*D*Y(116)*Y(165)+RATE(2318)*D*Y(245&
    &)*Y(333)+RATE(2372)*D*Y(13)*Y(233)+RATE(2378)*D*Y(13)*Y(261)+RATE(2416)&
    &*D*Y(13)*Y(265)+RATE(2423)*D*Y(13)*Y(300)+RATE(2425)*D*Y(13)*Y(316)&
    &+RATE(2429)*D*Y(13)*Y(277)+RATE(2444)*D*Y(13)*Y(302)+RATE(2452)*D*Y(92)&
    &*Y(174)+RATE(2463)*D*Y(173)*Y(173)+RATE(2464)*D*Y(174)*Y(333)+RATE(2465)&
    &*D*Y(174)*Y(183)+RATE(2467)*D*Y(325)*Y(333)+RATE(2484)*D*Y(69)*Y(165)&
    &+RATE(2529)*D*Y(27)*Y(233)+RATE(2534)*D*Y(27)*Y(244)+RATE(2537)*D*Y(27)&
    &*Y(173)+RATE(2544)*D*Y(27)*Y(265)+RATE(2549)*D*Y(27)*Y(316)+RATE(2551)*D&
    &*Y(27)*Y(277)+RATE(2745)*D*Y(44)*Y(183)+RATE(2765)*D*Y(54)*Y(165)&
    &+RATE(2785)*D*Y(54)*Y(174)+RATE(2813)*D*Y(133)*Y(165)+RATE(2823)*D*Y(266&
    &)*Y(333)+RATE(2833)*D*Y(46)*Y(234)+RATE(2839)*D*Y(46)*Y(245)+RATE(2846)&
    &*D*Y(46)*Y(266)+RATE(2871)*D*Y(46)*Y(233)+RATE(2887)*D*Y(46)*Y(173)&
    &+RATE(2895)*D*Y(46)*Y(265)+RATE(2900)*D*Y(46)*Y(300)+RATE(2903)*D*Y(46)&
    &*Y(316)+RATE(2905)*D*Y(46)*Y(277)+RATE(2976)*D*Y(301)*Y(333)+RATE(3062)&
    &*D*Y(165)*Y(207)+RATE(3063)*D*Y(165)*Y(302)+RATE(3067)*D*Y(165)*Y(333)&
    &+RATE(3068)*D*Y(318)*Y(333)+RATE(3068)*D*Y(318)*Y(333)+RATE(3075)*D&
    &*Y(105)*Y(165)+RATE(3089)*D*Y(121)*Y(165)+RATE(3106)*D*Y(303)*Y(333)&
    &+RATE(3107)*D*Y(278)*Y(333)+RATE(3115)*D*Y(319)*Y(333)
    YDOT(163) = PROD-LOSS
    LOSS = RATE(272)*Y(164)/safeMantle+RATE(355)*D*Y(164)/safeMantle*Y(2&
    &)+RATE(438)*Y(164)/safeMantle+RATE(666)*Y(164)*Y(3)+RATE(770)*Y(164)*Y(3&
    &)+RATE(1037)*Y(164)*totalSwap/safeMantle+RATE(1120)*Y(164)
    PROD = RATE(70)*Y(171)*bulkLayersReciprocal+RATE(590)*D*Y(163)&
    &+RATE(591)*D*Y(165)
    YDOT(164) = PROD-LOSS
    LOSS = RATE(591)*D*Y(165)+RATE(1282)*D*Y(16)*Y(165)+RATE(1359)*D&
    &*Y(67)*Y(165)+RATE(1421)*D*Y(109)*Y(165)+RATE(1467)*D*Y(20)*Y(165)&
    &+RATE(1570)*D*Y(24)*Y(165)+RATE(1617)*D*Y(32)*Y(165)+RATE(1698)*D*Y(41)&
    &*Y(165)+RATE(1699)*D*Y(41)*Y(165)+RATE(1980)*D*Y(6)*Y(165)+RATE(2004)*D&
    &*Y(6)*Y(165)+RATE(2060)*D*Y(131)*Y(165)+RATE(2061)*D*Y(131)*Y(165)&
    &+RATE(2279)*D*Y(116)*Y(165)+RATE(2286)*D*Y(116)*Y(165)+RATE(2484)*D*Y(69&
    &)*Y(165)+RATE(2641)*D*Y(35)*Y(165)+RATE(2765)*D*Y(54)*Y(165)+RATE(2813)&
    &*D*Y(133)*Y(165)+RATE(2957)*D*Y(161)*Y(165)+RATE(2991)*D*Y(56)*Y(165)&
    &+RATE(3062)*D*Y(165)*Y(207)+RATE(3063)*D*Y(165)*Y(302)+RATE(3064)*D&
    &*Y(165)*Y(183)+RATE(3065)*D*Y(165)*Y(183)+RATE(3066)*D*Y(165)*Y(300)&
    &+RATE(3067)*D*Y(165)*Y(333)+RATE(3075)*D*Y(105)*Y(165)+RATE(3089)*D&
    &*Y(121)*Y(165)+RATE(3090)*D*Y(121)*Y(165)
    PROD = RATE(188)*Y(163)+RATE(879)*Y(234)+RATE(912)*Y(174)+RATE(943)&
    &*Y(163)+RATE(1336)*D*Y(18)*Y(277)+RATE(1367)*D*Y(68)*Y(163)+RATE(1390)*D&
    &*Y(76)*Y(163)+RATE(1503)*D*Y(21)*Y(163)+RATE(1768)*D*Y(83)*Y(163)&
    &+RATE(1792)*D*Y(100)*Y(163)+RATE(1817)*D*Y(2)*Y(174)+RATE(1902)*D*Y(4)&
    &*Y(163)+RATE(1928)*D*Y(4)*Y(183)+RATE(1935)*D*Y(4)*Y(173)+RATE(2043)*D&
    &*Y(8)*Y(183)+RATE(2062)*D*Y(132)*Y(163)+RATE(2117)*D*Y(63)*Y(163)&
    &+RATE(2257)*D*Y(91)*Y(163)+RATE(2371)*D*Y(13)*Y(233)+RATE(2377)*D*Y(13)&
    &*Y(261)+RATE(2384)*D*Y(13)*Y(183)+RATE(2403)*D*Y(13)*Y(324)+RATE(2405)*D&
    &*Y(13)*Y(173)+RATE(2415)*D*Y(13)*Y(265)+RATE(2422)*D*Y(13)*Y(300)&
    &+RATE(2425)*D*Y(13)*Y(316)+RATE(2426)*D*Y(13)*Y(320)+RATE(2428)*D*Y(13)&
    &*Y(277)+RATE(2443)*D*Y(13)*Y(302)+RATE(2585)*D*Y(28)*Y(183)+RATE(2596)*D&
    &*Y(28)*Y(300)+RATE(2607)*D*Y(104)*Y(163)+RATE(2611)*D*Y(104)*Y(183)&
    &+RATE(2613)*D*Y(104)*Y(300)+RATE(2666)*D*Y(36)*Y(163)+RATE(2732)*D*Y(44)&
    &*Y(163)+RATE(2840)*D*Y(46)*Y(174)+RATE(2944)*D*Y(48)*Y(183)+RATE(2950)*D&
    &*Y(48)*Y(300)+RATE(2966)*D*Y(162)*Y(163)+RATE(3028)*D*Y(57)*Y(163)&
    &+RATE(3051)*D*Y(163)*Y(18)+RATE(3052)*D*Y(163)*Y(184)+RATE(3053)*D*Y(163&
    &)*Y(174)
    YDOT(165) = PROD-LOSS
    LOSS = RATE(197)*Y(166)+RATE(608)*D*Y(166)+RATE(958)*Y(166)+RATE(959&
    &)*Y(166)+RATE(960)*Y(166)+RATE(1415)*D*Y(81)*Y(166)+RATE(1416)*D*Y(81)&
    &*Y(166)+RATE(1417)*D*Y(81)*Y(166)+RATE(1418)*D*Y(81)*Y(166)+RATE(1665)*D&
    &*Y(33)*Y(166)+RATE(1739)*D*Y(53)*Y(166)+RATE(1761)*D*Y(82)*Y(166)&
    &+RATE(1911)*D*Y(4)*Y(166)+RATE(1940)*D*Y(4)*Y(166)+RATE(2203)*D*Y(10)&
    &*Y(166)+RATE(2311)*D*Y(117)*Y(166)+RATE(2438)*D*Y(13)*Y(166)+RATE(2439)&
    &*D*Y(13)*Y(166)+RATE(2913)*D*Y(46)*Y(166)
    PROD = RATE(281)*Y(158)/safeMantle+RATE(364)*D*Y(158)/safeMantle*Y(2&
    &)+RATE(447)*Y(158)/safeMantle+RATE(775)*Y(151)*Y(3)+RATE(827)*Y(155)*Y(5&
    &)*bulkLayersReciprocal+RATE(1129)*Y(158)+RATE(1212)*Y(170)+RATE(2106)*D&
    &*Y(62)*Y(177)+RATE(3102)*D*Y(177)*Y(333)
    YDOT(166) = PROD-LOSS
    LOSS = RATE(609)*D*Y(167)+RATE(1781)*D*Y(99)*Y(167)+RATE(1982)*D*Y(6&
    &)*Y(167)+RATE(2105)*D*Y(62)*Y(167)+RATE(3099)*D*Y(167)*Y(333)+RATE(3100)&
    &*D*Y(167)*Y(333)
    PROD = RATE(1911)*D*Y(4)*Y(166)+RATE(2202)*D*Y(10)*Y(146)
    YDOT(167) = PROD-LOSS
    LOSS = RATE(65)*Y(168)*bulkLayersReciprocal+RATE(1198)*Y(168)
    PROD = RATE(715)*Y(51)*Y(51)*bulkLayersReciprocal+RATE(1032)*Y(156)&
    &*totalSwap/safeMantle
    YDOT(168) = PROD-LOSS
    LOSS = RATE(22)*Y(169)*bulkLayersReciprocal+RATE(1155)*Y(169)
    PROD = RATE(696)*Y(5)*Y(154)*bulkLayersReciprocal+RATE(697)*Y(5)&
    &*Y(152)*bulkLayersReciprocal+RATE(701)*Y(152)*Y(228)&
    &*bulkLayersReciprocal+RATE(989)*Y(157)*totalSwap/safeMantle
    YDOT(169) = PROD-LOSS
    LOSS = RATE(79)*Y(170)*bulkLayersReciprocal+RATE(1212)*Y(170)
    PROD = RATE(723)*Y(155)*Y(5)*bulkLayersReciprocal+RATE(1046)*Y(158)&
    &*totalSwap/safeMantle
    YDOT(170) = PROD-LOSS
    LOSS = RATE(70)*Y(171)*bulkLayersReciprocal+RATE(718)*Y(171)*Y(5)&
    &*bulkLayersReciprocal+RATE(822)*Y(171)*Y(5)*bulkLayersReciprocal&
    &+RATE(1203)*Y(171)
    PROD = RATE(1037)*Y(164)*totalSwap/safeMantle
    YDOT(171) = PROD-LOSS
    LOSS = RATE(490)*D*Y(172)+RATE(1673)*D*Y(172)*Y(333)+RATE(1674)*D&
    &*Y(172)*Y(333)+RATE(1675)*D*Y(172)*Y(333)+RATE(1676)*D*Y(172)*Y(333)&
    &+RATE(1677)*D*Y(172)*Y(333)+RATE(1678)*D*Y(172)*Y(54)+RATE(2245)*D*Y(90)&
    &*Y(172)
    PROD = RATE(1509)*D*Y(21)*Y(159)+RATE(1666)*D*Y(33)*Y(62)+RATE(1712)&
    &*D*Y(42)*Y(159)+RATE(2067)*D*Y(132)*Y(159)+RATE(2160)*D*Y(10)*Y(159)&
    &+RATE(2212)*D*Y(142)*Y(159)+RATE(2223)*D*Y(66)*Y(159)+RATE(2297)*D*Y(117&
    &)*Y(159)
    YDOT(172) = PROD-LOSS
    LOSS = RATE(163)*Y(173)+RATE(551)*D*Y(173)+RATE(911)*Y(173)&
    &+RATE(1255)*D*Y(16)*Y(173)+RATE(1256)*D*Y(16)*Y(173)+RATE(1324)*D*Y(18)&
    &*Y(173)+RATE(1656)*D*Y(33)*Y(173)+RATE(1844)*D*Y(2)*Y(173)+RATE(1890)*D&
    &*Y(4)*Y(173)+RATE(1935)*D*Y(4)*Y(173)+RATE(1990)*D*Y(6)*Y(173)+RATE(2182&
    &)*D*Y(10)*Y(173)+RATE(2304)*D*Y(117)*Y(173)+RATE(2405)*D*Y(13)*Y(173)&
    &+RATE(2463)*D*Y(173)*Y(173)+RATE(2463)*D*Y(173)*Y(173)+RATE(2536)*D*Y(27&
    &)*Y(173)+RATE(2537)*D*Y(27)*Y(173)+RATE(2887)*D*Y(46)*Y(173)+RATE(2888)&
    &*D*Y(46)*Y(173)+RATE(3059)*D*Y(163)*Y(173)
    PROD = RATE(150)*Y(330)+RATE(150)*Y(330)+RATE(164)*Y(324)+RATE(252)&
    &*Y(179)/safeMantle+RATE(335)*D*Y(179)/safeMantle*Y(2)+RATE(418)*Y(179&
    &)/safeMantle+RATE(770)*Y(164)*Y(3)+RATE(822)*Y(171)*Y(5)&
    &*bulkLayersReciprocal+RATE(892)*Y(183)+RATE(894)*Y(330)+RATE(894)*Y(330)&
    &+RATE(915)*Y(324)+RATE(1100)*Y(179)+RATE(1183)*Y(181)+RATE(1496)*D*Y(20)&
    &*Y(163)+RATE(1497)*D*Y(20)*Y(277)+RATE(1625)*D*Y(32)*Y(183)+RATE(1703)*D&
    &*Y(41)*Y(163)+RATE(1795)*D*Y(100)*Y(183)+RATE(1821)*D*Y(2)*Y(303)&
    &+RATE(1835)*D*Y(2)*Y(183)+RATE(1852)*D*Y(2)*Y(265)+RATE(1861)*D*Y(2)&
    &*Y(300)+RATE(1863)*D*Y(2)*Y(316)+RATE(1864)*D*Y(2)*Y(277)+RATE(1998)*D&
    &*Y(6)*Y(163)+RATE(2061)*D*Y(131)*Y(165)+RATE(2090)*D*Y(62)*Y(184)&
    &+RATE(2128)*D*Y(63)*Y(183)+RATE(2141)*D*Y(184)*Y(333)+RATE(2143)*D*Y(184&
    &)*Y(183)+RATE(2146)*D*Y(329)*Y(333)+RATE(2146)*D*Y(329)*Y(333)+RATE(2236&
    &)*D*Y(189)*Y(333)+RATE(2237)*D*Y(189)*Y(333)+RATE(2381)*D*Y(13)*Y(330)&
    &+RATE(2403)*D*Y(13)*Y(324)+RATE(2467)*D*Y(325)*Y(333)+RATE(2469)*D*Y(307&
    &)*Y(333)+RATE(2480)*D*Y(69)*Y(174)+RATE(2584)*D*Y(28)*Y(183)+RATE(2659)&
    &*D*Y(35)*Y(163)+RATE(2744)*D*Y(44)*Y(183)+RATE(2761)*D*Y(54)*Y(174)&
    &+RATE(2774)*D*Y(54)*Y(184)+RATE(2802)*D*Y(55)*Y(183)+RATE(2811)*D*Y(133)&
    &*Y(174)+RATE(2876)*D*Y(46)*Y(183)+RATE(2882)*D*Y(46)*Y(244)+RATE(3001)*D&
    &*Y(56)*Y(233)+RATE(3005)*D*Y(56)*Y(183)+RATE(3053)*D*Y(163)*Y(174)&
    &+RATE(3058)*D*Y(163)*Y(116)+RATE(3072)*D*Y(105)*Y(174)
    YDOT(173) = PROD-LOSS
    LOSS = RATE(552)*D*Y(174)+RATE(912)*Y(174)+RATE(913)*Y(174)&
    &+RATE(1234)*D*Y(16)*Y(174)+RATE(1457)*D*Y(20)*Y(174)+RATE(1691)*D*Y(41)&
    &*Y(174)+RATE(1817)*D*Y(2)*Y(174)+RATE(1968)*D*Y(6)*Y(174)+RATE(2003)*D&
    &*Y(6)*Y(174)+RATE(2097)*D*Y(62)*Y(174)+RATE(2251)*D*Y(90)*Y(174)&
    &+RATE(2452)*D*Y(92)*Y(174)+RATE(2464)*D*Y(174)*Y(333)+RATE(2465)*D*Y(174&
    &)*Y(183)+RATE(2466)*D*Y(174)*Y(183)+RATE(2480)*D*Y(69)*Y(174)+RATE(2501)&
    &*D*Y(27)*Y(174)+RATE(2761)*D*Y(54)*Y(174)+RATE(2785)*D*Y(54)*Y(174)&
    &+RATE(2811)*D*Y(133)*Y(174)+RATE(2840)*D*Y(46)*Y(174)+RATE(2841)*D*Y(46)&
    &*Y(174)+RATE(3053)*D*Y(163)*Y(174)+RATE(3072)*D*Y(105)*Y(174)
    PROD = RATE(1541)*D*Y(21)*Y(163)+RATE(1738)*D*Y(53)*Y(163)+RATE(1815&
    &)*D*Y(2)*Y(184)+RATE(1890)*D*Y(4)*Y(173)+RATE(1927)*D*Y(4)*Y(183)&
    &+RATE(1937)*D*Y(4)*Y(300)+RATE(1980)*D*Y(6)*Y(165)+RATE(2042)*D*Y(8)&
    &*Y(183)+RATE(2070)*D*Y(132)*Y(163)+RATE(2133)*D*Y(63)*Y(163)+RATE(2197)&
    &*D*Y(10)*Y(163)+RATE(2266)*D*Y(91)*Y(163)+RATE(2286)*D*Y(116)*Y(165)&
    &+RATE(2308)*D*Y(117)*Y(163)+RATE(2381)*D*Y(13)*Y(330)+RATE(2383)*D*Y(13)&
    &*Y(183)+RATE(2458)*D*Y(145)*Y(163)+RATE(2583)*D*Y(28)*Y(183)+RATE(2610)&
    &*D*Y(104)*Y(183)+RATE(2619)*D*Y(120)*Y(163)+RATE(2695)*D*Y(36)*Y(163)&
    &+RATE(2743)*D*Y(44)*Y(183)+RATE(2754)*D*Y(44)*Y(163)+RATE(2835)*D*Y(46)&
    &*Y(184)+RATE(2943)*D*Y(48)*Y(183)+RATE(3046)*D*Y(57)*Y(163)+RATE(3055)*D&
    &*Y(163)*Y(176)
    YDOT(174) = PROD-LOSS
    LOSS = RATE(183)*Y(175)+RATE(583)*D*Y(175)+RATE(935)*Y(175)+RATE(936&
    &)*Y(175)+RATE(1489)*D*Y(20)*Y(175)+RATE(1490)*D*Y(20)*Y(175)+RATE(1635)&
    &*D*Y(32)*Y(175)+RATE(1786)*D*Y(99)*Y(175)+RATE(1855)*D*Y(2)*Y(175)&
    &+RATE(1856)*D*Y(2)*Y(175)+RATE(1857)*D*Y(2)*Y(175)+RATE(2293)*D*Y(116)&
    &*Y(175)+RATE(2546)*D*Y(27)*Y(175)+RATE(2897)*D*Y(46)*Y(175)+RATE(3013)*D&
    &*Y(56)*Y(175)
    PROD = RATE(268)*Y(178)/safeMantle+RATE(351)*D*Y(178)/safeMantle*Y(2&
    &)+RATE(434)*Y(178)/safeMantle+RATE(1116)*Y(178)+RATE(1199)*Y(180)&
    &+RATE(1420)*D*Y(89)*Y(161)+RATE(1634)*D*Y(32)*Y(161)+RATE(1701)*D*Y(41)&
    &*Y(161)+RATE(1994)*D*Y(6)*Y(161)+RATE(2069)*D*Y(132)*Y(161)+RATE(2292)*D&
    &*Y(116)*Y(161)+RATE(2970)*D*Y(162)*Y(267)
    YDOT(175) = PROD-LOSS
    LOSS = RATE(584)*D*Y(176)+RATE(1238)*D*Y(16)*Y(176)+RATE(1358)*D&
    &*Y(67)*Y(176)+RATE(1382)*D*Y(75)*Y(176)+RATE(1465)*D*Y(20)*Y(176)&
    &+RATE(1568)*D*Y(24)*Y(176)+RATE(1746)*D*Y(82)*Y(176)+RATE(1779)*D*Y(99)&
    &*Y(176)+RATE(1978)*D*Y(6)*Y(176)+RATE(2059)*D*Y(131)*Y(176)+RATE(2102)*D&
    &*Y(62)*Y(176)+RATE(2254)*D*Y(90)*Y(176)+RATE(2285)*D*Y(116)*Y(176)&
    &+RATE(2454)*D*Y(92)*Y(176)+RATE(2600)*D*Y(103)*Y(176)+RATE(2639)*D*Y(35)&
    &*Y(176)+RATE(2720)*D*Y(43)*Y(176)+RATE(2789)*D*Y(54)*Y(176)+RATE(2816)*D&
    &*Y(133)*Y(176)+RATE(2847)*D*Y(46)*Y(176)+RATE(2972)*D*Y(176)*Y(333)&
    &+RATE(2973)*D*Y(176)*Y(232)+RATE(2990)*D*Y(56)*Y(176)+RATE(3055)*D*Y(163&
    &)*Y(176)
    PROD = RATE(2050)*D*Y(8)*Y(161)+RATE(2191)*D*Y(10)*Y(161)+RATE(2284)&
    &*D*Y(116)*Y(162)+RATE(2692)*D*Y(36)*Y(161)
    YDOT(176) = PROD-LOSS
    LOSS = RATE(610)*D*Y(177)+RATE(2106)*D*Y(62)*Y(177)+RATE(3101)*D&
    &*Y(177)*Y(333)+RATE(3102)*D*Y(177)*Y(333)
    PROD = RATE(1982)*D*Y(6)*Y(167)+RATE(2007)*D*Y(6)*Y(147)+RATE(2203)&
    &*D*Y(10)*Y(166)+RATE(2311)*D*Y(117)*Y(166)
    YDOT(177) = PROD-LOSS
    LOSS = RATE(268)*Y(178)/safeMantle+RATE(351)*D*Y(178)/safeMantle*Y(2&
    &)+RATE(434)*Y(178)/safeMantle+RATE(1033)*Y(178)*totalSwap/safeMantle&
    &+RATE(1116)*Y(178)
    PROD = RATE(66)*Y(180)*bulkLayersReciprocal+RATE(583)*D*Y(175)&
    &+RATE(584)*D*Y(176)
    YDOT(178) = PROD-LOSS
    LOSS = RATE(252)*Y(179)/safeMantle+RATE(335)*D*Y(179)/safeMantle*Y(2&
    &)+RATE(418)*Y(179)/safeMantle+RATE(652)*Y(179)*Y(3)+RATE(756)*Y(179)*Y(3&
    &)+RATE(1017)*Y(179)*totalSwap/safeMantle+RATE(1100)*Y(179)
    PROD = RATE(50)*Y(181)*bulkLayersReciprocal+RATE(92)*Y(182)+RATE(551&
    &)*D*Y(173)+RATE(552)*D*Y(174)+RATE(666)*Y(164)*Y(3)
    YDOT(179) = PROD-LOSS
    LOSS = RATE(66)*Y(180)*bulkLayersReciprocal+RATE(1199)*Y(180)
    PROD = RATE(1033)*Y(178)*totalSwap/safeMantle
    YDOT(180) = PROD-LOSS
    LOSS = RATE(50)*Y(181)*bulkLayersReciprocal+RATE(704)*Y(181)*Y(5)&
    &*bulkLayersReciprocal+RATE(808)*Y(181)*Y(5)*bulkLayersReciprocal&
    &+RATE(1183)*Y(181)
    PROD = RATE(718)*Y(171)*Y(5)*bulkLayersReciprocal+RATE(1017)*Y(179)&
    &*totalSwap/safeMantle
    YDOT(181) = PROD-LOSS
    LOSS = RATE(92)*Y(182)+RATE(237)*Y(182)/safeMantle+RATE(320)*D*Y(182&
    &)/safeMantle*Y(2)+RATE(403)*Y(182)/safeMantle+RATE(1002)*Y(182)&
    &*totalSwap/safeMantle+RATE(1085)*Y(182)
    PROD = RATE(35)*Y(185)*bulkLayersReciprocal+RATE(517)*D*Y(183)&
    &+RATE(518)*D*Y(184)+RATE(527)*D*Y(189)+RATE(652)*Y(179)*Y(3)
    YDOT(182) = PROD-LOSS
    LOSS = RATE(148)*Y(183)+RATE(149)*Y(183)+RATE(517)*D*Y(183)+RATE(891&
    &)*Y(183)+RATE(892)*Y(183)+RATE(893)*Y(183)+RATE(1287)*D*Y(18)*Y(183)&
    &+RATE(1318)*D*Y(18)*Y(183)+RATE(1403)*D*Y(81)*Y(183)+RATE(1411)*D*Y(81)&
    &*Y(183)+RATE(1522)*D*Y(21)*Y(183)+RATE(1523)*D*Y(21)*Y(183)+RATE(1607)*D&
    &*Y(25)*Y(183)+RATE(1608)*D*Y(25)*Y(183)+RATE(1609)*D*Y(25)*Y(183)&
    &+RATE(1625)*D*Y(32)*Y(183)+RATE(1654)*D*Y(33)*Y(183)+RATE(1706)*D*Y(42)&
    &*Y(183)+RATE(1718)*D*Y(42)*Y(183)+RATE(1732)*D*Y(53)*Y(183)+RATE(1788)*D&
    &*Y(100)*Y(183)+RATE(1795)*D*Y(100)*Y(183)+RATE(1835)*D*Y(2)*Y(183)&
    &+RATE(1885)*D*Y(4)*Y(183)+RATE(1927)*D*Y(4)*Y(183)+RATE(1928)*D*Y(4)&
    &*Y(183)+RATE(2018)*D*Y(8)*Y(183)+RATE(2042)*D*Y(8)*Y(183)+RATE(2043)*D&
    &*Y(8)*Y(183)+RATE(2111)*D*Y(63)*Y(183)+RATE(2127)*D*Y(63)*Y(183)&
    &+RATE(2128)*D*Y(63)*Y(183)+RATE(2136)*D*Y(183)*Y(301)+RATE(2137)*D*Y(183&
    &)*Y(199)+RATE(2138)*D*Y(183)*Y(325)+RATE(2139)*D*Y(183)*Y(307)+RATE(2140&
    &)*D*Y(183)*Y(278)+RATE(2143)*D*Y(184)*Y(183)+RATE(2172)*D*Y(10)*Y(183)&
    &+RATE(2213)*D*Y(142)*Y(183)+RATE(2226)*D*Y(66)*Y(183)+RATE(2274)*D*Y(102&
    &)*Y(183)+RATE(2275)*D*Y(102)*Y(183)+RATE(2301)*D*Y(117)*Y(183)+RATE(2327&
    &)*D*Y(13)*Y(183)+RATE(2383)*D*Y(13)*Y(183)+RATE(2384)*D*Y(13)*Y(183)&
    &+RATE(2465)*D*Y(174)*Y(183)+RATE(2466)*D*Y(174)*Y(183)+RATE(2561)*D*Y(28&
    &)*Y(183)+RATE(2583)*D*Y(28)*Y(183)+RATE(2584)*D*Y(28)*Y(183)+RATE(2585)&
    &*D*Y(28)*Y(183)+RATE(2602)*D*Y(104)*Y(183)+RATE(2610)*D*Y(104)*Y(183)&
    &+RATE(2611)*D*Y(104)*Y(183)+RATE(2728)*D*Y(44)*Y(183)+RATE(2742)*D*Y(44)&
    &*Y(183)+RATE(2743)*D*Y(44)*Y(183)+RATE(2744)*D*Y(44)*Y(183)+RATE(2745)*D&
    &*Y(44)*Y(183)+RATE(2802)*D*Y(55)*Y(183)+RATE(2876)*D*Y(46)*Y(183)&
    &+RATE(2926)*D*Y(48)*Y(183)+RATE(2943)*D*Y(48)*Y(183)+RATE(2944)*D*Y(48)&
    &*Y(183)+RATE(2965)*D*Y(162)*Y(183)+RATE(3005)*D*Y(56)*Y(183)+RATE(3023)&
    &*D*Y(57)*Y(183)+RATE(3037)*D*Y(57)*Y(183)+RATE(3064)*D*Y(165)*Y(183)&
    &+RATE(3065)*D*Y(165)*Y(183)
    PROD = RATE(237)*Y(182)/safeMantle+RATE(320)*D*Y(182)/safeMantle*Y(2&
    &)+RATE(403)*Y(182)/safeMantle+RATE(756)*Y(179)*Y(3)+RATE(808)*Y(181)*Y(5&
    &)*bulkLayersReciprocal+RATE(1085)*Y(182)+RATE(1168)*Y(185)+RATE(1990)*D&
    &*Y(6)*Y(173)+RATE(2056)*D*Y(131)*Y(189)+RATE(2099)*D*Y(62)*Y(307)&
    &+RATE(2144)*D*Y(184)*Y(333)+RATE(2235)*D*Y(189)*Y(333)+RATE(2248)*D*Y(90&
    &)*Y(189)+RATE(2277)*D*Y(116)*Y(184)+RATE(2449)*D*Y(92)*Y(189)+RATE(2463)&
    &*D*Y(173)*Y(173)+RATE(2478)*D*Y(69)*Y(184)+RATE(2759)*D*Y(54)*Y(184)&
    &+RATE(2777)*D*Y(54)*Y(189)+RATE(2809)*D*Y(133)*Y(184)+RATE(3052)*D*Y(163&
    &)*Y(184)+RATE(3071)*D*Y(105)*Y(184)
    YDOT(183) = PROD-LOSS
    LOSS = RATE(518)*D*Y(184)+RATE(1228)*D*Y(16)*Y(184)+RATE(1815)*D*Y(2&
    &)*Y(184)+RATE(1965)*D*Y(6)*Y(184)+RATE(2090)*D*Y(62)*Y(184)+RATE(2141)*D&
    &*Y(184)*Y(333)+RATE(2142)*D*Y(184)*Y(333)+RATE(2143)*D*Y(184)*Y(183)&
    &+RATE(2144)*D*Y(184)*Y(333)+RATE(2277)*D*Y(116)*Y(184)+RATE(2478)*D*Y(69&
    &)*Y(184)+RATE(2500)*D*Y(27)*Y(184)+RATE(2759)*D*Y(54)*Y(184)+RATE(2774)&
    &*D*Y(54)*Y(184)+RATE(2809)*D*Y(133)*Y(184)+RATE(2835)*D*Y(46)*Y(184)&
    &+RATE(2836)*D*Y(46)*Y(184)+RATE(3052)*D*Y(163)*Y(184)+RATE(3071)*D*Y(105&
    &)*Y(184)
    PROD = RATE(148)*Y(183)+RATE(891)*Y(183)+RATE(1287)*D*Y(18)*Y(183)&
    &+RATE(1403)*D*Y(81)*Y(183)+RATE(1706)*D*Y(42)*Y(183)+RATE(1788)*D*Y(100)&
    &*Y(183)+RATE(1816)*D*Y(2)*Y(189)+RATE(1885)*D*Y(4)*Y(183)+RATE(1968)*D&
    &*Y(6)*Y(174)+RATE(2004)*D*Y(6)*Y(165)+RATE(2018)*D*Y(8)*Y(183)+RATE(2060&
    &)*D*Y(131)*Y(165)+RATE(2111)*D*Y(63)*Y(183)+RATE(2136)*D*Y(183)*Y(301)&
    &+RATE(2182)*D*Y(10)*Y(173)+RATE(2304)*D*Y(117)*Y(173)+RATE(2327)*D*Y(13)&
    &*Y(183)+RATE(2561)*D*Y(28)*Y(183)+RATE(2602)*D*Y(104)*Y(183)+RATE(2728)&
    &*D*Y(44)*Y(183)+RATE(2926)*D*Y(48)*Y(183)+RATE(2965)*D*Y(162)*Y(183)&
    &+RATE(3023)*D*Y(57)*Y(183)
    YDOT(184) = PROD-LOSS
    LOSS = RATE(35)*Y(185)*bulkLayersReciprocal+RATE(1168)*Y(185)
    PROD = RATE(704)*Y(181)*Y(5)*bulkLayersReciprocal+RATE(1002)*Y(182)&
    &*totalSwap/safeMantle
    YDOT(185) = PROD-LOSS
    LOSS = RATE(100)*Y(186)+RATE(137)*Y(186)+RATE(494)*D*Y(186)+RATE(872&
    &)*Y(186)+RATE(1740)*D*Y(186)*Y(4)+RATE(1741)*D*Y(186)*Y(6)+RATE(2168)*D&
    &*Y(10)*Y(186)
    PROD = RATE(153)*Y(193)+RATE(226)*Y(188)/safeMantle+RATE(309)*D&
    &*Y(188)/safeMantle*Y(2)+RATE(392)*Y(188)/safeMantle+RATE(900)*Y(193)&
    &+RATE(1074)*Y(188)+RATE(1157)*Y(190)+RATE(1742)*D*Y(187)*Y(2)+RATE(1743)&
    &*D*Y(187)*Y(333)+RATE(2053)*D*Y(196)*Y(333)+RATE(2239)*D*Y(193)*Y(2)&
    &+RATE(2240)*D*Y(194)*Y(333)+RATE(2954)*D*Y(161)*Y(187)
    YDOT(186) = PROD-LOSS
    LOSS = RATE(495)*D*Y(187)+RATE(1742)*D*Y(187)*Y(2)+RATE(1743)*D&
    &*Y(187)*Y(333)+RATE(1963)*D*Y(6)*Y(187)+RATE(2954)*D*Y(161)*Y(187)
    PROD = RATE(100)*Y(186)+RATE(137)*Y(186)+RATE(872)*Y(186)+RATE(1740)&
    &*D*Y(186)*Y(4)+RATE(2397)*D*Y(13)*Y(193)
    YDOT(187) = PROD-LOSS
    LOSS = RATE(226)*Y(188)/safeMantle+RATE(309)*D*Y(188)/safeMantle*Y(2&
    &)+RATE(392)*Y(188)/safeMantle+RATE(624)*Y(188)*Y(3)+RATE(728)*Y(188)*Y(3&
    &)+RATE(991)*Y(188)*totalSwap/safeMantle+RATE(1074)*Y(188)
    PROD = RATE(24)*Y(190)*bulkLayersReciprocal+RATE(494)*D*Y(186)&
    &+RATE(495)*D*Y(187)
    YDOT(188) = PROD-LOSS
    LOSS = RATE(527)*D*Y(189)+RATE(1816)*D*Y(2)*Y(189)+RATE(2056)*D&
    &*Y(131)*Y(189)+RATE(2235)*D*Y(189)*Y(333)+RATE(2236)*D*Y(189)*Y(333)&
    &+RATE(2237)*D*Y(189)*Y(333)+RATE(2238)*D*Y(189)*Y(333)+RATE(2248)*D*Y(90&
    &)*Y(189)+RATE(2449)*D*Y(92)*Y(189)+RATE(2777)*D*Y(54)*Y(189)+RATE(3054)&
    &*D*Y(163)*Y(189)
    PROD = RATE(1411)*D*Y(81)*Y(183)+RATE(1522)*D*Y(21)*Y(183)+RATE(1608&
    &)*D*Y(25)*Y(183)+RATE(1718)*D*Y(42)*Y(183)+RATE(1732)*D*Y(53)*Y(183)&
    &+RATE(1965)*D*Y(6)*Y(184)+RATE(2003)*D*Y(6)*Y(174)+RATE(2127)*D*Y(63)&
    &*Y(183)+RATE(2138)*D*Y(183)*Y(325)+RATE(2139)*D*Y(183)*Y(307)+RATE(2143)&
    &*D*Y(184)*Y(183)+RATE(2172)*D*Y(10)*Y(183)+RATE(2213)*D*Y(142)*Y(183)&
    &+RATE(2226)*D*Y(66)*Y(183)+RATE(2274)*D*Y(102)*Y(183)+RATE(2275)*D*Y(102&
    &)*Y(183)+RATE(2301)*D*Y(117)*Y(183)+RATE(2465)*D*Y(174)*Y(183)+RATE(2742&
    &)*D*Y(44)*Y(183)+RATE(3037)*D*Y(57)*Y(183)
    YDOT(189) = PROD-LOSS
    LOSS = RATE(24)*Y(190)*bulkLayersReciprocal+RATE(676)*Y(190)*Y(5)&
    &*bulkLayersReciprocal+RATE(780)*Y(190)*Y(5)*bulkLayersReciprocal&
    &+RATE(1157)*Y(190)
    PROD = RATE(991)*Y(188)*totalSwap/safeMantle
    YDOT(190) = PROD-LOSS
    LOSS = RATE(242)*Y(191)/safeMantle+RATE(325)*D*Y(191)/safeMantle*Y(2&
    &)+RATE(408)*Y(191)/safeMantle+RATE(1007)*Y(191)*totalSwap/safeMantle&
    &+RATE(1090)*Y(191)
    PROD = RATE(40)*Y(195)*bulkLayersReciprocal+RATE(508)*D*Y(196)&
    &+RATE(529)*D*Y(193)+RATE(530)*D*Y(194)+RATE(624)*Y(188)*Y(3)
    YDOT(191) = PROD-LOSS
    LOSS = RATE(469)*D*Y(192)+RATE(1427)*D*Y(192)*Y(333)+RATE(2244)*D&
    &*Y(90)*Y(192)
    PROD = RATE(1222)*D*Y(16)*Y(76)+RATE(1303)*D*Y(18)*Y(75)+RATE(1320)&
    &*D*Y(18)*Y(288)+RATE(1369)*D*Y(68)*Y(67)+RATE(1443)*D*Y(20)*Y(68)&
    &+RATE(1506)*D*Y(21)*Y(67)+RATE(1507)*D*Y(21)*Y(75)+RATE(2353)*D*Y(13)&
    &*Y(206)
    YDOT(192) = PROD-LOSS
    LOSS = RATE(153)*Y(193)+RATE(529)*D*Y(193)+RATE(900)*Y(193)+RATE(901&
    &)*Y(193)+RATE(1735)*D*Y(53)*Y(193)+RATE(1888)*D*Y(4)*Y(193)+RATE(2178)*D&
    &*Y(10)*Y(193)+RATE(2239)*D*Y(193)*Y(2)+RATE(2397)*D*Y(13)*Y(193)
    PROD = RATE(242)*Y(191)/safeMantle+RATE(325)*D*Y(191)/safeMantle*Y(2&
    &)+RATE(408)*Y(191)/safeMantle+RATE(728)*Y(188)*Y(3)+RATE(780)*Y(190)*Y(5&
    &)*bulkLayersReciprocal+RATE(1090)*Y(191)+RATE(1173)*Y(195)+RATE(1741)*D&
    &*Y(186)*Y(6)+RATE(1775)*D*Y(99)*Y(196)+RATE(2054)*D*Y(196)*Y(333)&
    &+RATE(2089)*D*Y(62)*Y(196)
    YDOT(193) = PROD-LOSS
    LOSS = RATE(530)*D*Y(194)+RATE(1967)*D*Y(6)*Y(194)+RATE(2240)*D&
    &*Y(194)*Y(333)
    PROD = RATE(901)*Y(193)+RATE(1888)*D*Y(4)*Y(193)+RATE(1963)*D*Y(6)&
    &*Y(187)+RATE(2168)*D*Y(10)*Y(186)
    YDOT(194) = PROD-LOSS
    LOSS = RATE(40)*Y(195)*bulkLayersReciprocal+RATE(1173)*Y(195)
    PROD = RATE(676)*Y(190)*Y(5)*bulkLayersReciprocal+RATE(1007)*Y(191)&
    &*totalSwap/safeMantle
    YDOT(195) = PROD-LOSS
    LOSS = RATE(508)*D*Y(196)+RATE(1775)*D*Y(99)*Y(196)+RATE(2053)*D&
    &*Y(196)*Y(333)+RATE(2054)*D*Y(196)*Y(333)+RATE(2089)*D*Y(62)*Y(196)
    PROD = RATE(1735)*D*Y(53)*Y(193)+RATE(1967)*D*Y(6)*Y(194)+RATE(2178)&
    &*D*Y(10)*Y(193)
    YDOT(196) = PROD-LOSS
    LOSS = RATE(119)*Y(197)+RATE(120)*Y(197)+RATE(465)*D*Y(197)+RATE(842&
    &)*Y(197)+RATE(843)*Y(197)+RATE(1244)*D*Y(16)*Y(197)+RATE(1875)*D*Y(4)&
    &*Y(197)+RATE(2151)*D*Y(10)*Y(197)+RATE(2345)*D*Y(13)*Y(197)+RATE(2516)*D&
    &*Y(27)*Y(197)+RATE(2864)*D*Y(46)*Y(197)
    PROD = RATE(210)*Y(198)/safeMantle+RATE(293)*D*Y(198)/safeMantle*Y(2&
    &)+RATE(376)*Y(198)/safeMantle+RATE(1058)*Y(198)+RATE(1141)*Y(202)&
    &+RATE(1253)*D*Y(16)*Y(101)+RATE(1258)*D*Y(16)*Y(291)+RATE(1424)*D*Y(290)&
    &*Y(333)+RATE(1426)*D*Y(204)*Y(333)+RATE(1429)*D*Y(309)*Y(333)+RATE(2515)&
    &*D*Y(27)*Y(75)+RATE(2518)*D*Y(27)*Y(284)+RATE(2538)*D*Y(27)*Y(291)&
    &+RATE(2865)*D*Y(46)*Y(284)
    YDOT(197) = PROD-LOSS
    LOSS = RATE(210)*Y(198)/safeMantle+RATE(293)*D*Y(198)/safeMantle*Y(2&
    &)+RATE(376)*Y(198)/safeMantle+RATE(975)*Y(198)*totalSwap/safeMantle&
    &+RATE(1058)*Y(198)
    PROD = RATE(8)*Y(202)*bulkLayersReciprocal
    YDOT(198) = PROD-LOSS
    LOSS = RATE(466)*D*Y(199)+RATE(1422)*D*Y(199)*Y(333)+RATE(1423)*D&
    &*Y(199)*Y(333)+RATE(2082)*D*Y(62)*Y(199)+RATE(2083)*D*Y(62)*Y(199)&
    &+RATE(2137)*D*Y(183)*Y(199)+RATE(2769)*D*Y(54)*Y(199)
    PROD = RATE(1290)*D*Y(18)*Y(291)+RATE(1319)*D*Y(18)*Y(288)+RATE(1323&
    &)*D*Y(18)*Y(92)+RATE(1514)*D*Y(21)*Y(82)+RATE(1524)*D*Y(21)*Y(90)&
    &+RATE(1875)*D*Y(4)*Y(197)+RATE(2386)*D*Y(13)*Y(288)+RATE(2491)*D*Y(27)&
    &*Y(76)+RATE(2493)*D*Y(27)*Y(81)+RATE(2587)*D*Y(28)*Y(291)+RATE(2625)*D&
    &*Y(35)*Y(68)+RATE(2669)*D*Y(36)*Y(67)
    YDOT(199) = PROD-LOSS
    LOSS = RATE(470)*D*Y(200)+RATE(1245)*D*Y(16)*Y(200)+RATE(2517)*D&
    &*Y(27)*Y(200)
    PROD = RATE(211)*Y(201)/safeMantle+RATE(294)*D*Y(201)/safeMantle*Y(2&
    &)+RATE(377)*Y(201)/safeMantle+RATE(1059)*Y(201)+RATE(1142)*Y(203)&
    &+RATE(1242)*D*Y(16)*Y(89)+RATE(1471)*D*Y(20)*Y(80)
    YDOT(200) = PROD-LOSS
    LOSS = RATE(211)*Y(201)/safeMantle+RATE(294)*D*Y(201)/safeMantle*Y(2&
    &)+RATE(377)*Y(201)/safeMantle+RATE(976)*Y(201)*totalSwap/safeMantle&
    &+RATE(1059)*Y(201)
    PROD = RATE(9)*Y(203)*bulkLayersReciprocal+RATE(470)*D*Y(200)
    YDOT(201) = PROD-LOSS
    LOSS = RATE(8)*Y(202)*bulkLayersReciprocal+RATE(1141)*Y(202)
    PROD = RATE(975)*Y(198)*totalSwap/safeMantle
    YDOT(202) = PROD-LOSS
    LOSS = RATE(9)*Y(203)*bulkLayersReciprocal+RATE(1142)*Y(203)
    PROD = RATE(976)*Y(201)*totalSwap/safeMantle
    YDOT(203) = PROD-LOSS
    LOSS = RATE(468)*D*Y(204)+RATE(1426)*D*Y(204)*Y(333)
    PROD = RATE(1525)*D*Y(21)*Y(90)+RATE(2082)*D*Y(62)*Y(199)+RATE(2151)&
    &*D*Y(10)*Y(197)+RATE(2494)*D*Y(27)*Y(81)+RATE(2704)*D*Y(43)*Y(68)
    YDOT(204) = PROD-LOSS
    LOSS = RATE(220)*Y(205)/safeMantle+RATE(303)*D*Y(205)/safeMantle*Y(2&
    &)+RATE(386)*Y(205)/safeMantle+RATE(985)*Y(205)*totalSwap/safeMantle&
    &+RATE(1068)*Y(205)
    PROD = RATE(18)*Y(210)*bulkLayersReciprocal+RATE(471)*D*Y(213)&
    &+RATE(484)*D*Y(206)
    YDOT(205) = PROD-LOSS
    LOSS = RATE(484)*D*Y(206)+RATE(1307)*D*Y(18)*Y(206)+RATE(2221)*D&
    &*Y(66)*Y(206)+RATE(2295)*D*Y(117)*Y(206)+RATE(2353)*D*Y(13)*Y(206)
    PROD = RATE(220)*Y(205)/safeMantle+RATE(303)*D*Y(205)/safeMantle*Y(2&
    &)+RATE(386)*Y(205)/safeMantle+RATE(1068)*Y(205)+RATE(1151)*Y(210)&
    &+RATE(1243)*D*Y(16)*Y(123)+RATE(1428)*D*Y(213)*Y(333)+RATE(1472)*D*Y(20)&
    &*Y(109)
    YDOT(206) = PROD-LOSS
    LOSS = RATE(191)*Y(207)+RATE(596)*D*Y(207)+RATE(947)*Y(207)&
    &+RATE(1299)*D*Y(18)*Y(207)+RATE(1338)*D*Y(18)*Y(207)+RATE(1908)*D*Y(4)&
    &*Y(207)+RATE(2432)*D*Y(13)*Y(207)+RATE(2433)*D*Y(13)*Y(207)+RATE(2552)*D&
    &*Y(27)*Y(207)+RATE(2908)*D*Y(46)*Y(207)+RATE(2909)*D*Y(46)*Y(207)&
    &+RATE(3062)*D*Y(165)*Y(207)
    PROD = RATE(192)*Y(292)+RATE(275)*Y(209)/safeMantle+RATE(358)*D&
    &*Y(209)/safeMantle*Y(2)+RATE(441)*Y(209)/safeMantle+RATE(949)*Y(314)&
    &+RATE(1123)*Y(209)+RATE(1206)*Y(211)+RATE(1277)*D*Y(16)*Y(121)+RATE(2906&
    &)*D*Y(46)*Y(292)+RATE(3086)*D*Y(293)*Y(333)+RATE(3088)*D*Y(313)*Y(333)
    YDOT(207) = PROD-LOSS
    LOSS = RATE(597)*D*Y(208)+RATE(2507)*D*Y(27)*Y(208)+RATE(2849)*D&
    &*Y(46)*Y(208)+RATE(3084)*D*Y(208)*Y(333)
    PROD = RATE(1240)*D*Y(16)*Y(122)+RATE(1299)*D*Y(18)*Y(207)+RATE(1339&
    &)*D*Y(18)*Y(135)+RATE(1340)*D*Y(18)*Y(121)+RATE(1342)*D*Y(18)*Y(302)&
    &+RATE(1360)*D*Y(67)*Y(236)+RATE(1468)*D*Y(20)*Y(106)+RATE(1908)*D*Y(4)&
    &*Y(207)+RATE(3062)*D*Y(165)*Y(207)
    YDOT(208) = PROD-LOSS
    LOSS = RATE(275)*Y(209)/safeMantle+RATE(358)*D*Y(209)/safeMantle*Y(2&
    &)+RATE(441)*Y(209)/safeMantle+RATE(1040)*Y(209)*totalSwap/safeMantle&
    &+RATE(1123)*Y(209)
    PROD = RATE(73)*Y(211)*bulkLayersReciprocal+RATE(596)*D*Y(207)&
    &+RATE(597)*D*Y(208)
    YDOT(209) = PROD-LOSS
    LOSS = RATE(18)*Y(210)*bulkLayersReciprocal+RATE(1151)*Y(210)
    PROD = RATE(985)*Y(205)*totalSwap/safeMantle
    YDOT(210) = PROD-LOSS
    LOSS = RATE(73)*Y(211)*bulkLayersReciprocal+RATE(1206)*Y(211)
    PROD = RATE(1040)*Y(209)*totalSwap/safeMantle
    YDOT(211) = PROD-LOSS
    LOSS = RATE(222)*Y(212)/safeMantle+RATE(305)*D*Y(212)/safeMantle*Y(2&
    &)+RATE(388)*Y(212)/safeMantle+RATE(987)*Y(212)*totalSwap/safeMantle&
    &+RATE(1070)*Y(212)
    PROD = RATE(20)*Y(215)*bulkLayersReciprocal+RATE(486)*D*Y(214)
    YDOT(212) = PROD-LOSS
    LOSS = RATE(471)*D*Y(213)+RATE(1428)*D*Y(213)*Y(333)
    PROD = RATE(1409)*D*Y(81)*Y(214)+RATE(1649)*D*Y(33)*Y(109)+RATE(1683&
    &)*D*Y(41)*Y(81)+RATE(2221)*D*Y(66)*Y(206)+RATE(2295)*D*Y(117)*Y(206)
    YDOT(213) = PROD-LOSS
    LOSS = RATE(133)*Y(214)+RATE(486)*D*Y(214)+RATE(862)*Y(214)&
    &+RATE(1409)*D*Y(81)*Y(214)+RATE(1410)*D*Y(81)*Y(214)+RATE(1651)*D*Y(33)&
    &*Y(214)+RATE(1669)*D*Y(214)*Y(243)+RATE(1919)*D*Y(4)*Y(214)+RATE(2158)*D&
    &*Y(10)*Y(214)+RATE(2222)*D*Y(66)*Y(214)+RATE(2270)*D*Y(102)*Y(214)&
    &+RATE(2271)*D*Y(102)*Y(214)+RATE(2296)*D*Y(117)*Y(214)+RATE(2355)*D*Y(13&
    &)*Y(214)+RATE(2356)*D*Y(13)*Y(214)+RATE(2616)*D*Y(120)*Y(214)
    PROD = RATE(222)*Y(212)/safeMantle+RATE(305)*D*Y(212)/safeMantle*Y(2&
    &)+RATE(388)*Y(212)/safeMantle+RATE(1070)*Y(212)+RATE(1153)*Y(215)&
    &+RATE(1642)*D*Y(32)*Y(82)+RATE(1670)*D*Y(220)*Y(333)
    YDOT(214) = PROD-LOSS
    LOSS = RATE(20)*Y(215)*bulkLayersReciprocal+RATE(1153)*Y(215)
    PROD = RATE(987)*Y(212)*totalSwap/safeMantle
    YDOT(215) = PROD-LOSS
    LOSS = RATE(217)*Y(216)/safeMantle+RATE(300)*D*Y(216)/safeMantle*Y(2&
    &)+RATE(383)*Y(216)/safeMantle+RATE(982)*Y(216)*totalSwap/safeMantle&
    &+RATE(1065)*Y(216)
    PROD = RATE(15)*Y(223)*bulkLayersReciprocal+RATE(480)*D*Y(218)
    YDOT(216) = PROD-LOSS
    LOSS = RATE(223)*Y(217)/safeMantle+RATE(306)*D*Y(217)/safeMantle*Y(2&
    &)+RATE(389)*Y(217)/safeMantle+RATE(988)*Y(217)*totalSwap/safeMantle&
    &+RATE(1071)*Y(217)
    PROD = RATE(21)*Y(224)*bulkLayersReciprocal+RATE(465)*D*Y(197)&
    &+RATE(466)*D*Y(199)+RATE(468)*D*Y(204)+RATE(487)*D*Y(219)+RATE(488)*D&
    &*Y(220)
    YDOT(217) = PROD-LOSS
    LOSS = RATE(127)*Y(218)+RATE(480)*D*Y(218)+RATE(854)*Y(218)&
    &+RATE(1826)*D*Y(2)*Y(218)+RATE(2350)*D*Y(13)*Y(218)+RATE(2351)*D*Y(13)&
    &*Y(218)
    PROD = RATE(217)*Y(216)/safeMantle+RATE(300)*D*Y(216)/safeMantle*Y(2&
    &)+RATE(383)*Y(216)/safeMantle+RATE(1065)*Y(216)+RATE(1148)*Y(223)&
    &+RATE(2856)*D*Y(46)*Y(89)+RATE(2858)*D*Y(46)*Y(109)+RATE(2994)*D*Y(56)&
    &*Y(80)
    YDOT(218) = PROD-LOSS
    LOSS = RATE(487)*D*Y(219)
    PROD = RATE(223)*Y(217)/safeMantle+RATE(306)*D*Y(217)/safeMantle*Y(2&
    &)+RATE(389)*Y(217)/safeMantle+RATE(1071)*Y(217)+RATE(1154)*Y(224)
    YDOT(219) = PROD-LOSS
    LOSS = RATE(488)*D*Y(220)+RATE(1670)*D*Y(220)*Y(333)+RATE(1671)*D&
    &*Y(220)*Y(333)
    PROD = RATE(1410)*D*Y(81)*Y(214)+RATE(1667)*D*Y(33)*Y(90)+RATE(1669)&
    &*D*Y(214)*Y(243)+RATE(2158)*D*Y(10)*Y(214)+RATE(2222)*D*Y(66)*Y(214)&
    &+RATE(2245)*D*Y(90)*Y(172)+RATE(2270)*D*Y(102)*Y(214)+RATE(2271)*D*Y(102&
    &)*Y(214)+RATE(2296)*D*Y(117)*Y(214)+RATE(2616)*D*Y(120)*Y(214)
    YDOT(220) = PROD-LOSS
    LOSS = RATE(184)*Y(221)+RATE(585)*D*Y(221)+RATE(937)*Y(221)&
    &+RATE(1269)*D*Y(16)*Y(221)+RATE(1331)*D*Y(18)*Y(221)+RATE(1858)*D*Y(2)&
    &*Y(221)+RATE(1859)*D*Y(2)*Y(221)+RATE(1860)*D*Y(2)*Y(221)+RATE(2418)*D&
    &*Y(13)*Y(221)+RATE(2419)*D*Y(13)*Y(221)+RATE(2819)*D*Y(133)*Y(221)&
    &+RATE(2898)*D*Y(46)*Y(221)+RATE(2899)*D*Y(46)*Y(221)+RATE(2960)*D*Y(161)&
    &*Y(221)+RATE(2961)*D*Y(161)*Y(221)
    PROD = RATE(269)*Y(222)/safeMantle+RATE(352)*D*Y(222)/safeMantle*Y(2&
    &)+RATE(435)*Y(222)/safeMantle+RATE(726)*Y(34)*Y(226)+RATE(753)*Y(143)&
    &*Y(226)+RATE(761)*Y(37)*Y(226)+RATE(762)*Y(45)*Y(226)+RATE(769)*Y(58)&
    &*Y(226)+RATE(778)*Y(38)*Y(228)*bulkLayersReciprocal+RATE(805)*Y(152)&
    &*Y(228)*bulkLayersReciprocal+RATE(813)*Y(39)*Y(228)*bulkLayersReciprocal&
    &+RATE(814)*Y(50)*Y(228)*bulkLayersReciprocal+RATE(821)*Y(60)*Y(228)&
    &*bulkLayersReciprocal+RATE(1117)*Y(222)+RATE(1200)*Y(225)+RATE(1484)*D&
    &*Y(20)*Y(133)+RATE(1755)*D*Y(82)*Y(264)+RATE(1757)*D*Y(82)*Y(133)&
    &+RATE(1759)*D*Y(82)*Y(161)+RATE(2533)*D*Y(27)*Y(116)+RATE(2873)*D*Y(46)&
    &*Y(101)+RATE(2879)*D*Y(46)*Y(90)+RATE(2999)*D*Y(56)*Y(82)+RATE(3010)*D&
    &*Y(56)*Y(291)
    YDOT(221) = PROD-LOSS
    LOSS = RATE(269)*Y(222)/safeMantle+RATE(352)*D*Y(222)/safeMantle*Y(2&
    &)+RATE(435)*Y(222)/safeMantle+RATE(621)*Y(17)*Y(222)+RATE(646)*Y(3)&
    &*Y(222)+RATE(664)*Y(47)*Y(222)+RATE(725)*Y(17)*Y(222)+RATE(750)*Y(3)&
    &*Y(222)+RATE(768)*Y(47)*Y(222)+RATE(1034)*Y(222)*totalSwap/safeMantle&
    &+RATE(1117)*Y(222)
    PROD = RATE(67)*Y(225)*bulkLayersReciprocal+RATE(585)*D*Y(221)&
    &+RATE(622)*Y(34)*Y(226)+RATE(649)*Y(143)*Y(226)+RATE(657)*Y(37)*Y(226)&
    &+RATE(658)*Y(45)*Y(226)+RATE(665)*Y(58)*Y(226)
    YDOT(222) = PROD-LOSS
    LOSS = RATE(15)*Y(223)*bulkLayersReciprocal+RATE(1148)*Y(223)
    PROD = RATE(982)*Y(216)*totalSwap/safeMantle
    YDOT(223) = PROD-LOSS
    LOSS = RATE(21)*Y(224)*bulkLayersReciprocal+RATE(1154)*Y(224)
    PROD = RATE(988)*Y(217)*totalSwap/safeMantle
    YDOT(224) = PROD-LOSS
    LOSS = RATE(67)*Y(225)*bulkLayersReciprocal+RATE(673)*Y(19)*Y(225)&
    &*bulkLayersReciprocal+RATE(698)*Y(5)*Y(225)*bulkLayersReciprocal&
    &+RATE(716)*Y(51)*Y(225)*bulkLayersReciprocal+RATE(777)*Y(19)*Y(225)&
    &*bulkLayersReciprocal+RATE(802)*Y(5)*Y(225)*bulkLayersReciprocal&
    &+RATE(820)*Y(51)*Y(225)*bulkLayersReciprocal+RATE(1200)*Y(225)
    PROD = RATE(674)*Y(38)*Y(228)*bulkLayersReciprocal+RATE(701)*Y(152)&
    &*Y(228)*bulkLayersReciprocal+RATE(709)*Y(39)*Y(228)*bulkLayersReciprocal&
    &+RATE(710)*Y(50)*Y(228)*bulkLayersReciprocal+RATE(717)*Y(60)*Y(228)&
    &*bulkLayersReciprocal+RATE(1034)*Y(222)*totalSwap/safeMantle
    YDOT(225) = PROD-LOSS
    LOSS = RATE(109)*Y(226)+RATE(250)*Y(226)/safeMantle+RATE(333)*D&
    &*Y(226)/safeMantle*Y(2)+RATE(416)*Y(226)/safeMantle+RATE(622)*Y(34)&
    &*Y(226)+RATE(647)*Y(3)*Y(226)+RATE(649)*Y(143)*Y(226)+RATE(657)*Y(37)&
    &*Y(226)+RATE(658)*Y(45)*Y(226)+RATE(665)*Y(58)*Y(226)+RATE(726)*Y(34)&
    &*Y(226)+RATE(751)*Y(3)*Y(226)+RATE(753)*Y(143)*Y(226)+RATE(761)*Y(37)&
    &*Y(226)+RATE(762)*Y(45)*Y(226)+RATE(769)*Y(58)*Y(226)+RATE(829)*Y(226)&
    &+RATE(1015)*Y(226)*totalSwap/safeMantle+RATE(1098)*Y(226)
    PROD = RATE(48)*Y(228)*bulkLayersReciprocal+RATE(545)*D*Y(227)&
    &+RATE(646)*Y(3)*Y(222)+RATE(656)*Y(37)*Y(97)
    YDOT(226) = PROD-LOSS
    LOSS = RATE(161)*Y(227)+RATE(545)*D*Y(227)+RATE(909)*Y(227)&
    &+RATE(1933)*D*Y(4)*Y(227)+RATE(2455)*D*Y(227)*Y(16)
    PROD = RATE(250)*Y(226)/safeMantle+RATE(333)*D*Y(226)/safeMantle*Y(2&
    &)+RATE(416)*Y(226)/safeMantle+RATE(750)*Y(3)*Y(222)+RATE(760)*Y(37)*Y(97&
    &)+RATE(802)*Y(5)*Y(225)*bulkLayersReciprocal+RATE(812)*Y(39)*Y(111)&
    &*bulkLayersReciprocal+RATE(1098)*Y(226)+RATE(1181)*Y(228)+RATE(1585)*D&
    &*Y(24)*Y(133)
    YDOT(227) = PROD-LOSS
    LOSS = RATE(48)*Y(228)*bulkLayersReciprocal+RATE(674)*Y(38)*Y(228)&
    &*bulkLayersReciprocal+RATE(699)*Y(5)*Y(228)*bulkLayersReciprocal&
    &+RATE(701)*Y(152)*Y(228)*bulkLayersReciprocal+RATE(709)*Y(39)*Y(228)&
    &*bulkLayersReciprocal+RATE(710)*Y(50)*Y(228)*bulkLayersReciprocal&
    &+RATE(717)*Y(60)*Y(228)*bulkLayersReciprocal+RATE(778)*Y(38)*Y(228)&
    &*bulkLayersReciprocal+RATE(803)*Y(5)*Y(228)*bulkLayersReciprocal&
    &+RATE(805)*Y(152)*Y(228)*bulkLayersReciprocal+RATE(813)*Y(39)*Y(228)&
    &*bulkLayersReciprocal+RATE(814)*Y(50)*Y(228)*bulkLayersReciprocal&
    &+RATE(821)*Y(60)*Y(228)*bulkLayersReciprocal+RATE(1181)*Y(228)
    PROD = RATE(698)*Y(5)*Y(225)*bulkLayersReciprocal+RATE(708)*Y(39)&
    &*Y(111)*bulkLayersReciprocal+RATE(1015)*Y(226)*totalSwap/safeMantle
    YDOT(228) = PROD-LOSS
    LOSS = RATE(282)*Y(229)/safeMantle+RATE(365)*D*Y(229)/safeMantle*Y(2&
    &)+RATE(448)*Y(229)/safeMantle+RATE(672)*Y(229)*Y(3)+RATE(776)*Y(229)*Y(3&
    &)+RATE(1047)*Y(229)*totalSwap/safeMantle+RATE(1130)*Y(229)
    PROD = RATE(80)*Y(239)*bulkLayersReciprocal+RATE(611)*D*Y(235)&
    &+RATE(612)*D*Y(236)
    YDOT(229) = PROD-LOSS
    LOSS = RATE(230)*Y(230)/safeMantle+RATE(313)*D*Y(230)/safeMantle*Y(2&
    &)+RATE(396)*Y(230)/safeMantle+RATE(627)*Y(230)*Y(3)+RATE(731)*Y(230)*Y(3&
    &)+RATE(995)*Y(230)*totalSwap/safeMantle+RATE(1078)*Y(230)
    PROD = RATE(28)*Y(240)*bulkLayersReciprocal+RATE(95)*Y(246)+RATE(501&
    &)*D*Y(233)+RATE(502)*D*Y(234)
    YDOT(230) = PROD-LOSS
    LOSS = RATE(229)*Y(231)/safeMantle+RATE(312)*D*Y(231)/safeMantle*Y(2&
    &)+RATE(395)*Y(231)/safeMantle+RATE(994)*Y(231)*totalSwap/safeMantle&
    &+RATE(1077)*Y(231)
    PROD = RATE(27)*Y(241)*bulkLayersReciprocal+RATE(500)*D*Y(232)&
    &+RATE(536)*D*Y(243)+RATE(626)*Y(97)*Y(58)+RATE(648)*Y(130)*Y(47)
    YDOT(231) = PROD-LOSS
    LOSS = RATE(140)*Y(232)+RATE(500)*D*Y(232)+RATE(876)*Y(232)&
    &+RATE(1312)*D*Y(18)*Y(232)+RATE(1475)*D*Y(20)*Y(232)+RATE(1515)*D*Y(21)&
    &*Y(232)+RATE(1604)*D*Y(25)*Y(232)+RATE(1714)*D*Y(42)*Y(232)+RATE(1728)*D&
    &*Y(53)*Y(232)+RATE(1830)*D*Y(2)*Y(232)+RATE(1924)*D*Y(4)*Y(232)&
    &+RATE(2037)*D*Y(8)*Y(232)+RATE(2164)*D*Y(10)*Y(232)+RATE(2259)*D*Y(91)&
    &*Y(232)+RATE(2366)*D*Y(13)*Y(232)+RATE(2367)*D*Y(13)*Y(232)+RATE(2368)*D&
    &*Y(13)*Y(232)+RATE(2369)*D*Y(13)*Y(232)+RATE(2457)*D*Y(145)*Y(232)&
    &+RATE(2528)*D*Y(27)*Y(232)+RATE(2579)*D*Y(28)*Y(232)+RATE(2617)*D*Y(120)&
    &*Y(232)+RATE(2673)*D*Y(36)*Y(232)+RATE(2674)*D*Y(36)*Y(232)+RATE(2675)*D&
    &*Y(36)*Y(232)+RATE(2870)*D*Y(46)*Y(232)+RATE(2941)*D*Y(48)*Y(232)&
    &+RATE(2973)*D*Y(176)*Y(232)+RATE(3033)*D*Y(57)*Y(232)+RATE(3077)*D*Y(105&
    &)*Y(232)
    PROD = RATE(229)*Y(231)/safeMantle+RATE(312)*D*Y(231)/safeMantle*Y(2&
    &)+RATE(395)*Y(231)/safeMantle+RATE(730)*Y(97)*Y(58)+RATE(752)*Y(130)&
    &*Y(47)+RATE(782)*Y(111)*Y(60)*bulkLayersReciprocal+RATE(804)*Y(139)*Y(51&
    &)*bulkLayersReciprocal+RATE(1077)*Y(231)+RATE(1160)*Y(241)+RATE(1232)*D&
    &*Y(16)*Y(243)+RATE(1485)*D*Y(20)*Y(161)+RATE(1586)*D*Y(24)*Y(161)&
    &+RATE(1587)*D*Y(24)*Y(161)+RATE(1669)*D*Y(214)*Y(243)+RATE(1689)*D*Y(41)&
    &*Y(243)+RATE(1776)*D*Y(99)*Y(243)+RATE(1780)*D*Y(99)*Y(319)+RATE(1782)*D&
    &*Y(99)*Y(236)+RATE(1783)*D*Y(99)*Y(144)+RATE(1784)*D*Y(99)*Y(264)&
    &+RATE(1785)*D*Y(99)*Y(161)+RATE(1786)*D*Y(99)*Y(175)+RATE(1796)*D*Y(100)&
    &*Y(320)+RATE(2094)*D*Y(62)*Y(243)+RATE(2291)*D*Y(116)*Y(161)+RATE(2315)&
    &*D*Y(243)*Y(333)+RATE(2782)*D*Y(54)*Y(243)+RATE(2819)*D*Y(133)*Y(221)&
    &+RATE(2880)*D*Y(46)*Y(116)+RATE(2900)*D*Y(46)*Y(300)+RATE(2950)*D*Y(48)&
    &*Y(300)+RATE(2960)*D*Y(161)*Y(221)+RATE(3000)*D*Y(56)*Y(99)+RATE(3114)*D&
    &*Y(278)*Y(300)
    YDOT(232) = PROD-LOSS
    LOSS = RATE(141)*Y(233)+RATE(142)*Y(233)+RATE(501)*D*Y(233)+RATE(877&
    &)*Y(233)+RATE(878)*Y(233)+RATE(1252)*D*Y(16)*Y(233)+RATE(1880)*D*Y(4)&
    &*Y(233)+RATE(2167)*D*Y(10)*Y(233)+RATE(2224)*D*Y(66)*Y(233)+RATE(2298)*D&
    &*Y(117)*Y(233)+RATE(2371)*D*Y(13)*Y(233)+RATE(2372)*D*Y(13)*Y(233)&
    &+RATE(2529)*D*Y(27)*Y(233)+RATE(2871)*D*Y(46)*Y(233)+RATE(2872)*D*Y(46)&
    &*Y(233)+RATE(3001)*D*Y(56)*Y(233)+RATE(3002)*D*Y(56)*Y(233)
    PROD = RATE(146)*Y(261)+RATE(230)*Y(230)/safeMantle+RATE(313)*D&
    &*Y(230)/safeMantle*Y(2)+RATE(396)*Y(230)/safeMantle+RATE(887)*Y(261)&
    &+RATE(1078)*Y(230)+RATE(1161)*Y(240)+RATE(1255)*D*Y(16)*Y(173)+RATE(1266&
    &)*D*Y(16)*Y(265)+RATE(1270)*D*Y(16)*Y(300)+RATE(1273)*D*Y(16)*Y(316)&
    &+RATE(1275)*D*Y(16)*Y(277)+RATE(1283)*D*Y(16)*Y(163)+RATE(1315)*D*Y(18)&
    &*Y(261)+RATE(1364)*D*Y(67)*Y(163)+RATE(1493)*D*Y(20)*Y(300)+RATE(1495)*D&
    &*Y(20)*Y(163)+RATE(1598)*D*Y(24)*Y(163)+RATE(1839)*D*Y(2)*Y(244)&
    &+RATE(2072)*D*Y(262)*Y(333)+RATE(2214)*D*Y(274)*Y(333)+RATE(2319)*D&
    &*Y(245)*Y(333)+RATE(2396)*D*Y(13)*Y(244)+RATE(2421)*D*Y(13)*Y(300)&
    &+RATE(2461)*D*Y(306)*Y(333)+RATE(2476)*D*Y(69)*Y(234)+RATE(2783)*D*Y(54)&
    &*Y(245)+RATE(2975)*D*Y(301)*Y(333)+RATE(3069)*D*Y(105)*Y(234)
    YDOT(233) = PROD-LOSS
    LOSS = RATE(502)*D*Y(234)+RATE(879)*Y(234)+RATE(1685)*D*Y(41)*Y(234)&
    &+RATE(1797)*D*Y(234)*Y(333)+RATE(1962)*D*Y(6)*Y(234)+RATE(2476)*D*Y(69)&
    &*Y(234)+RATE(2833)*D*Y(46)*Y(234)+RATE(2956)*D*Y(161)*Y(234)+RATE(3069)&
    &*D*Y(105)*Y(234)
    PROD = RATE(141)*Y(233)+RATE(877)*Y(233)+RATE(1234)*D*Y(16)*Y(174)&
    &+RATE(1282)*D*Y(16)*Y(165)+RATE(1324)*D*Y(18)*Y(173)+RATE(1328)*D*Y(18)&
    &*Y(265)+RATE(1332)*D*Y(18)*Y(300)+RATE(1335)*D*Y(18)*Y(277)+RATE(1346)*D&
    &*Y(18)*Y(163)+RATE(1359)*D*Y(67)*Y(165)+RATE(1372)*D*Y(68)*Y(163)&
    &+RATE(1467)*D*Y(20)*Y(165)+RATE(1540)*D*Y(21)*Y(163)+RATE(1880)*D*Y(4)&
    &*Y(233)+RATE(1932)*D*Y(4)*Y(244)+RATE(2376)*D*Y(13)*Y(261)+RATE(2395)*D&
    &*Y(13)*Y(244)+RATE(2420)*D*Y(13)*Y(300)+RATE(2595)*D*Y(28)*Y(300)
    YDOT(234) = PROD-LOSS
    LOSS = RATE(198)*Y(235)+RATE(611)*D*Y(235)+RATE(961)*Y(235)+RATE(962&
    &)*Y(235)+RATE(1341)*D*Y(18)*Y(235)+RATE(1913)*D*Y(4)*Y(235)+RATE(2205)*D&
    &*Y(10)*Y(235)+RATE(2234)*D*Y(66)*Y(235)+RATE(2313)*D*Y(117)*Y(235)&
    &+RATE(2441)*D*Y(13)*Y(235)+RATE(2442)*D*Y(13)*Y(235)+RATE(3050)*D*Y(57)&
    &*Y(235)
    PROD = RATE(151)*Y(263)+RATE(282)*Y(229)/safeMantle+RATE(365)*D&
    &*Y(229)/safeMantle*Y(2)+RATE(448)*Y(229)/safeMantle+RATE(895)*Y(263)&
    &+RATE(896)*Y(263)+RATE(1130)*Y(229)+RATE(1213)*Y(239)+RATE(2280)*D*Y(116&
    &)*Y(236)+RATE(2487)*D*Y(69)*Y(236)+RATE(2791)*D*Y(54)*Y(247)+RATE(2815)&
    &*D*Y(133)*Y(236)+RATE(2909)*D*Y(46)*Y(207)+RATE(2910)*D*Y(46)*Y(135)&
    &+RATE(2911)*D*Y(46)*Y(135)+RATE(2914)*D*Y(46)*Y(121)+RATE(2918)*D*Y(46)&
    &*Y(105)+RATE(2958)*D*Y(161)*Y(303)+RATE(3017)*D*Y(56)*Y(105)+RATE(3077)&
    &*D*Y(105)*Y(232)+RATE(3078)*D*Y(105)*Y(99)+RATE(3079)*D*Y(105)*Y(133)&
    &+RATE(3080)*D*Y(105)*Y(161)+RATE(3105)*D*Y(247)*Y(333)
    YDOT(235) = PROD-LOSS
    LOSS = RATE(612)*D*Y(236)+RATE(963)*Y(236)+RATE(1241)*D*Y(16)*Y(236)&
    &+RATE(1360)*D*Y(67)*Y(236)+RATE(1470)*D*Y(20)*Y(236)+RATE(1571)*D*Y(24)&
    &*Y(236)+RATE(1782)*D*Y(99)*Y(236)+RATE(1983)*D*Y(6)*Y(236)+RATE(2280)*D&
    &*Y(116)*Y(236)+RATE(2487)*D*Y(69)*Y(236)+RATE(2508)*D*Y(27)*Y(236)&
    &+RATE(2509)*D*Y(27)*Y(236)+RATE(2815)*D*Y(133)*Y(236)+RATE(2853)*D*Y(46)&
    &*Y(236)+RATE(3056)*D*Y(163)*Y(236)+RATE(3103)*D*Y(236)*Y(333)
    PROD = RATE(962)*Y(235)+RATE(1913)*D*Y(4)*Y(235)+RATE(2849)*D*Y(46)&
    &*Y(208)+RATE(2850)*D*Y(46)*Y(122)+RATE(2917)*D*Y(46)*Y(106)+RATE(2959)*D&
    &*Y(161)*Y(303)+RATE(2992)*D*Y(56)*Y(106)
    YDOT(236) = PROD-LOSS
    LOSS = RATE(131)*Y(237)+RATE(132)*Y(237)+RATE(485)*D*Y(237)+RATE(860&
    &)*Y(237)+RATE(861)*Y(237)+RATE(1650)*D*Y(33)*Y(237)+RATE(2155)*D*Y(10)&
    &*Y(237)+RATE(2156)*D*Y(10)*Y(237)+RATE(2157)*D*Y(10)*Y(237)+RATE(2354)*D&
    &*Y(13)*Y(237)
    PROD = RATE(221)*Y(238)/safeMantle+RATE(304)*D*Y(238)/safeMantle*Y(2&
    &)+RATE(387)*Y(238)/safeMantle+RATE(727)*Y(34)*Y(118)+RATE(779)*Y(38)&
    &*Y(126)*bulkLayersReciprocal+RATE(1069)*Y(238)+RATE(1152)*Y(242)&
    &+RATE(1473)*D*Y(20)*Y(159)+RATE(2861)*D*Y(46)*Y(123)+RATE(3018)*D*Y(56)&
    &*Y(89)
    YDOT(237) = PROD-LOSS
    LOSS = RATE(221)*Y(238)/safeMantle+RATE(304)*D*Y(238)/safeMantle*Y(2&
    &)+RATE(387)*Y(238)/safeMantle+RATE(986)*Y(238)*totalSwap/safeMantle&
    &+RATE(1069)*Y(238)
    PROD = RATE(19)*Y(242)*bulkLayersReciprocal+RATE(485)*D*Y(237)&
    &+RATE(623)*Y(34)*Y(118)
    YDOT(238) = PROD-LOSS
    LOSS = RATE(80)*Y(239)*bulkLayersReciprocal+RATE(724)*Y(239)*Y(5)&
    &*bulkLayersReciprocal+RATE(828)*Y(239)*Y(5)*bulkLayersReciprocal&
    &+RATE(1213)*Y(239)
    PROD = RATE(1047)*Y(229)*totalSwap/safeMantle
    YDOT(239) = PROD-LOSS
    LOSS = RATE(28)*Y(240)*bulkLayersReciprocal+RATE(679)*Y(240)*Y(5)&
    &*bulkLayersReciprocal+RATE(783)*Y(240)*Y(5)*bulkLayersReciprocal&
    &+RATE(1161)*Y(240)
    PROD = RATE(995)*Y(230)*totalSwap/safeMantle
    YDOT(240) = PROD-LOSS
    LOSS = RATE(27)*Y(241)*bulkLayersReciprocal+RATE(1160)*Y(241)
    PROD = RATE(678)*Y(111)*Y(60)*bulkLayersReciprocal+RATE(700)*Y(139)&
    &*Y(51)*bulkLayersReciprocal+RATE(994)*Y(231)*totalSwap/safeMantle
    YDOT(241) = PROD-LOSS
    LOSS = RATE(19)*Y(242)*bulkLayersReciprocal+RATE(1152)*Y(242)
    PROD = RATE(675)*Y(38)*Y(126)*bulkLayersReciprocal+RATE(986)*Y(238)&
    &*totalSwap/safeMantle
    YDOT(242) = PROD-LOSS
    LOSS = RATE(536)*D*Y(243)+RATE(1232)*D*Y(16)*Y(243)+RATE(1669)*D&
    &*Y(214)*Y(243)+RATE(1689)*D*Y(41)*Y(243)+RATE(1776)*D*Y(99)*Y(243)&
    &+RATE(2094)*D*Y(62)*Y(243)+RATE(2315)*D*Y(243)*Y(333)+RATE(2316)*D*Y(243&
    &)*Y(333)+RATE(2317)*D*Y(243)*Y(333)+RATE(2782)*D*Y(54)*Y(243)+RATE(2837)&
    &*D*Y(46)*Y(243)
    PROD = RATE(1714)*D*Y(42)*Y(232)+RATE(1728)*D*Y(53)*Y(232)+RATE(2037&
    &)*D*Y(8)*Y(232)+RATE(2164)*D*Y(10)*Y(232)+RATE(2259)*D*Y(91)*Y(232)&
    &+RATE(2457)*D*Y(145)*Y(232)+RATE(2617)*D*Y(120)*Y(232)+RATE(2673)*D*Y(36&
    &)*Y(232)+RATE(2970)*D*Y(162)*Y(267)+RATE(2973)*D*Y(176)*Y(232)+RATE(2986&
    &)*D*Y(56)*Y(117)+RATE(3033)*D*Y(57)*Y(232)
    YDOT(243) = PROD-LOSS
    LOSS = RATE(158)*Y(244)+RATE(539)*D*Y(244)+RATE(907)*Y(244)&
    &+RATE(1839)*D*Y(2)*Y(244)+RATE(1932)*D*Y(4)*Y(244)+RATE(2177)*D*Y(10)&
    &*Y(244)+RATE(2395)*D*Y(13)*Y(244)+RATE(2396)*D*Y(13)*Y(244)+RATE(2534)*D&
    &*Y(27)*Y(244)+RATE(2882)*D*Y(46)*Y(244)+RATE(2883)*D*Y(46)*Y(244)
    PROD = RATE(247)*Y(246)/safeMantle+RATE(330)*D*Y(246)/safeMantle*Y(2&
    &)+RATE(413)*Y(246)/safeMantle+RATE(731)*Y(230)*Y(3)+RATE(783)*Y(240)*Y(5&
    &)*bulkLayersReciprocal+RATE(1095)*Y(246)+RATE(1178)*Y(253)+RATE(1599)*D&
    &*Y(24)*Y(163)+RATE(2073)*D*Y(262)*Y(333)+RATE(3109)*D*Y(278)*Y(80)&
    &+RATE(3112)*D*Y(278)*Y(109)
    YDOT(244) = PROD-LOSS
    LOSS = RATE(540)*D*Y(245)+RATE(2318)*D*Y(245)*Y(333)+RATE(2319)*D&
    &*Y(245)*Y(333)+RATE(2783)*D*Y(54)*Y(245)+RATE(2838)*D*Y(46)*Y(245)&
    &+RATE(2839)*D*Y(46)*Y(245)
    PROD = RATE(158)*Y(244)+RATE(907)*Y(244)+RATE(1228)*D*Y(16)*Y(184)&
    &+RATE(1318)*D*Y(18)*Y(183)+RATE(1421)*D*Y(109)*Y(165)+RATE(1523)*D*Y(21)&
    &*Y(183)+RATE(1537)*D*Y(21)*Y(300)+RATE(1570)*D*Y(24)*Y(165)+RATE(1609)*D&
    &*Y(25)*Y(183)+RATE(1615)*D*Y(25)*Y(300)+RATE(1616)*D*Y(25)*Y(163)&
    &+RATE(1663)*D*Y(33)*Y(163)+RATE(1685)*D*Y(41)*Y(234)+RATE(1699)*D*Y(41)&
    &*Y(165)+RATE(1962)*D*Y(6)*Y(234)+RATE(2137)*D*Y(183)*Y(199)+RATE(2167)*D&
    &*Y(10)*Y(233)+RATE(2224)*D*Y(66)*Y(233)+RATE(2298)*D*Y(117)*Y(233)&
    &+RATE(3110)*D*Y(278)*Y(80)
    YDOT(245) = PROD-LOSS
    LOSS = RATE(95)*Y(246)+RATE(247)*Y(246)/safeMantle+RATE(330)*D*Y(246&
    &)/safeMantle*Y(2)+RATE(413)*Y(246)/safeMantle+RATE(651)*Y(246)*Y(3)&
    &+RATE(755)*Y(246)*Y(3)+RATE(1012)*Y(246)*totalSwap/safeMantle+RATE(1095)&
    &*Y(246)
    PROD = RATE(45)*Y(253)*bulkLayersReciprocal+RATE(539)*D*Y(244)&
    &+RATE(540)*D*Y(245)+RATE(627)*Y(230)*Y(3)
    YDOT(246) = PROD-LOSS
    LOSS = RATE(613)*D*Y(247)+RATE(2791)*D*Y(54)*Y(247)+RATE(3104)*D&
    &*Y(247)*Y(333)+RATE(3105)*D*Y(247)*Y(333)
    PROD = RATE(254)*Y(248)/safeMantle+RATE(337)*D*Y(248)/safeMantle*Y(2&
    &)+RATE(420)*Y(248)/safeMantle+RATE(1102)*Y(248)+RATE(1185)*Y(254)&
    &+RATE(1929)*D*Y(4)*Y(263)+RATE(1983)*D*Y(6)*Y(236)+RATE(2099)*D*Y(62)&
    &*Y(307)+RATE(2103)*D*Y(62)*Y(106)+RATE(2205)*D*Y(10)*Y(235)+RATE(2234)*D&
    &*Y(66)*Y(235)+RATE(2313)*D*Y(117)*Y(235)+RATE(2385)*D*Y(13)*Y(263)&
    &+RATE(2851)*D*Y(46)*Y(136)+RATE(2852)*D*Y(46)*Y(147)+RATE(3050)*D*Y(57)&
    &*Y(235)+RATE(3081)*D*Y(106)*Y(159)+RATE(3095)*D*Y(136)*Y(161)
    YDOT(247) = PROD-LOSS
    LOSS = RATE(254)*Y(248)/safeMantle+RATE(337)*D*Y(248)/safeMantle*Y(2&
    &)+RATE(420)*Y(248)/safeMantle+RATE(654)*Y(248)*Y(3)+RATE(758)*Y(248)*Y(3&
    &)+RATE(1019)*Y(248)*totalSwap/safeMantle+RATE(1102)*Y(248)
    PROD = RATE(52)*Y(254)*bulkLayersReciprocal+RATE(613)*D*Y(247)&
    &+RATE(672)*Y(229)*Y(3)
    YDOT(248) = PROD-LOSS
    LOSS = RATE(570)*D*Y(249)
    PROD = RATE(261)*Y(250)/safeMantle+RATE(344)*D*Y(250)/safeMantle*Y(2&
    &)+RATE(427)*Y(250)/safeMantle+RATE(763)*Y(45)*Y(118)+RATE(765)*Y(45)&
    &*Y(130)+RATE(815)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(817)*Y(50)&
    &*Y(139)*bulkLayersReciprocal+RATE(1109)*Y(250)+RATE(1192)*Y(255)&
    &+RATE(2727)*D*Y(43)*Y(131)
    YDOT(249) = PROD-LOSS
    LOSS = RATE(261)*Y(250)/safeMantle+RATE(344)*D*Y(250)/safeMantle*Y(2&
    &)+RATE(427)*Y(250)/safeMantle+RATE(1026)*Y(250)*totalSwap/safeMantle&
    &+RATE(1109)*Y(250)
    PROD = RATE(59)*Y(255)*bulkLayersReciprocal+RATE(570)*D*Y(249)&
    &+RATE(659)*Y(45)*Y(118)+RATE(661)*Y(45)*Y(130)
    YDOT(250) = PROD-LOSS
    LOSS = RATE(537)*D*Y(251)
    PROD = RATE(245)*Y(252)/safeMantle+RATE(328)*D*Y(252)/safeMantle*Y(2&
    &)+RATE(411)*Y(252)/safeMantle+RATE(1093)*Y(252)+RATE(1176)*Y(256)
    YDOT(251) = PROD-LOSS
    LOSS = RATE(245)*Y(252)/safeMantle+RATE(328)*D*Y(252)/safeMantle*Y(2&
    &)+RATE(411)*Y(252)/safeMantle+RATE(1010)*Y(252)*totalSwap/safeMantle&
    &+RATE(1093)*Y(252)
    PROD = RATE(43)*Y(256)*bulkLayersReciprocal+RATE(94)*Y(268)+RATE(537&
    &)*D*Y(251)
    YDOT(252) = PROD-LOSS
    LOSS = RATE(45)*Y(253)*bulkLayersReciprocal+RATE(703)*Y(253)*Y(5)&
    &*bulkLayersReciprocal+RATE(807)*Y(253)*Y(5)*bulkLayersReciprocal&
    &+RATE(1178)*Y(253)
    PROD = RATE(679)*Y(240)*Y(5)*bulkLayersReciprocal+RATE(1012)*Y(246)&
    &*totalSwap/safeMantle
    YDOT(253) = PROD-LOSS
    LOSS = RATE(52)*Y(254)*bulkLayersReciprocal+RATE(706)*Y(254)*Y(5)&
    &*bulkLayersReciprocal+RATE(810)*Y(254)*Y(5)*bulkLayersReciprocal&
    &+RATE(1185)*Y(254)
    PROD = RATE(724)*Y(239)*Y(5)*bulkLayersReciprocal+RATE(1019)*Y(248)&
    &*totalSwap/safeMantle
    YDOT(254) = PROD-LOSS
    LOSS = RATE(59)*Y(255)*bulkLayersReciprocal+RATE(1192)*Y(255)
    PROD = RATE(711)*Y(50)*Y(126)*bulkLayersReciprocal+RATE(713)*Y(50)&
    &*Y(139)*bulkLayersReciprocal+RATE(1026)*Y(250)*totalSwap/safeMantle
    YDOT(255) = PROD-LOSS
    LOSS = RATE(43)*Y(256)*bulkLayersReciprocal+RATE(1176)*Y(256)
    PROD = RATE(1010)*Y(252)*totalSwap/safeMantle
    YDOT(256) = PROD-LOSS
    LOSS = RATE(235)*Y(257)/safeMantle+RATE(318)*D*Y(257)/safeMantle*Y(2&
    &)+RATE(401)*Y(257)/safeMantle+RATE(1000)*Y(257)*totalSwap/safeMantle&
    &+RATE(1083)*Y(257)
    PROD = RATE(33)*Y(269)*bulkLayersReciprocal+RATE(512)*D*Y(261)&
    &+RATE(513)*D*Y(262)+RATE(525)*D*Y(274)+RATE(651)*Y(246)*Y(3)
    YDOT(257) = PROD-LOSS
    LOSS = RATE(264)*Y(258)/safeMantle+RATE(347)*D*Y(258)/safeMantle*Y(2&
    &)+RATE(430)*Y(258)/safeMantle+RATE(1029)*Y(258)*totalSwap/safeMantle&
    &+RATE(1112)*Y(258)
    PROD = RATE(62)*Y(270)*bulkLayersReciprocal+RATE(576)*D*Y(264)
    YDOT(258) = PROD-LOSS
    LOSS = RATE(265)*Y(259)/safeMantle+RATE(348)*D*Y(259)/safeMantle*Y(2&
    &)+RATE(431)*Y(259)/safeMantle+RATE(1030)*Y(259)*totalSwap/safeMantle&
    &+RATE(1113)*Y(259)
    PROD = RATE(63)*Y(271)*bulkLayersReciprocal+RATE(548)*D*Y(275)&
    &+RATE(577)*D*Y(265)+RATE(578)*D*Y(266)
    YDOT(259) = PROD-LOSS
    LOSS = RATE(239)*Y(260)/safeMantle+RATE(322)*D*Y(260)/safeMantle*Y(2&
    &)+RATE(405)*Y(260)/safeMantle+RATE(1004)*Y(260)*totalSwap/safeMantle&
    &+RATE(1087)*Y(260)
    PROD = RATE(37)*Y(272)*bulkLayersReciprocal+RATE(521)*D*Y(263)&
    &+RATE(654)*Y(248)*Y(3)
    YDOT(260) = PROD-LOSS
    LOSS = RATE(146)*Y(261)+RATE(512)*D*Y(261)+RATE(887)*Y(261)&
    &+RATE(1315)*D*Y(18)*Y(261)+RATE(1882)*D*Y(4)*Y(261)+RATE(2170)*D*Y(10)&
    &*Y(261)+RATE(2300)*D*Y(117)*Y(261)+RATE(2376)*D*Y(13)*Y(261)+RATE(2377)&
    &*D*Y(13)*Y(261)+RATE(2378)*D*Y(13)*Y(261)
    PROD = RATE(235)*Y(257)/safeMantle+RATE(318)*D*Y(257)/safeMantle*Y(2&
    &)+RATE(401)*Y(257)/safeMantle+RATE(755)*Y(246)*Y(3)+RATE(807)*Y(253)*Y(5&
    &)*bulkLayersReciprocal+RATE(1083)*Y(257)+RATE(1166)*Y(269)+RATE(1641)*D&
    &*Y(32)*Y(163)+RATE(2074)*D*Y(262)*Y(333)+RATE(2215)*D*Y(274)*Y(333)
    YDOT(261) = PROD-LOSS
    LOSS = RATE(513)*D*Y(262)+RATE(2072)*D*Y(262)*Y(333)+RATE(2073)*D&
    &*Y(262)*Y(333)+RATE(2074)*D*Y(262)*Y(333)
    PROD = RATE(1614)*D*Y(25)*Y(300)+RATE(1617)*D*Y(32)*Y(165)+RATE(1656&
    &)*D*Y(33)*Y(173)+RATE(1882)*D*Y(4)*Y(261)+RATE(2177)*D*Y(10)*Y(244)&
    &+RATE(3108)*D*Y(278)*Y(80)+RATE(3111)*D*Y(278)*Y(109)
    YDOT(262) = PROD-LOSS
    LOSS = RATE(151)*Y(263)+RATE(521)*D*Y(263)+RATE(895)*Y(263)+RATE(896&
    &)*Y(263)+RATE(1929)*D*Y(4)*Y(263)+RATE(2385)*D*Y(13)*Y(263)
    PROD = RATE(239)*Y(260)/safeMantle+RATE(322)*D*Y(260)/safeMantle*Y(2&
    &)+RATE(405)*Y(260)/safeMantle+RATE(758)*Y(248)*Y(3)+RATE(810)*Y(254)*Y(5&
    &)*bulkLayersReciprocal+RATE(1087)*Y(260)+RATE(1170)*Y(272)+RATE(2912)*D&
    &*Y(46)*Y(146)
    YDOT(263) = PROD-LOSS
    LOSS = RATE(178)*Y(264)+RATE(576)*D*Y(264)+RATE(930)*Y(264)&
    &+RATE(1582)*D*Y(24)*Y(264)+RATE(1630)*D*Y(32)*Y(264)+RATE(1755)*D*Y(82)&
    &*Y(264)+RATE(1784)*D*Y(99)*Y(264)+RATE(1849)*D*Y(2)*Y(264)+RATE(1936)*D&
    &*Y(4)*Y(264)+RATE(2188)*D*Y(10)*Y(264)+RATE(2540)*D*Y(27)*Y(264)&
    &+RATE(2541)*D*Y(27)*Y(264)+RATE(2542)*D*Y(27)*Y(264)+RATE(2649)*D*Y(35)&
    &*Y(264)+RATE(2893)*D*Y(46)*Y(264)+RATE(2949)*D*Y(48)*Y(264)
    PROD = RATE(264)*Y(258)/safeMantle+RATE(347)*D*Y(258)/safeMantle*Y(2&
    &)+RATE(430)*Y(258)/safeMantle+RATE(1112)*Y(258)+RATE(1195)*Y(270)&
    &+RATE(2818)*D*Y(133)*Y(161)+RATE(2884)*D*Y(46)*Y(144)+RATE(2961)*D*Y(161&
    &)*Y(221)+RATE(3012)*D*Y(56)*Y(133)
    YDOT(264) = PROD-LOSS
    LOSS = RATE(179)*Y(265)+RATE(577)*D*Y(265)+RATE(931)*Y(265)&
    &+RATE(1266)*D*Y(16)*Y(265)+RATE(1267)*D*Y(16)*Y(265)+RATE(1293)*D*Y(18)&
    &*Y(265)+RATE(1328)*D*Y(18)*Y(265)+RATE(1852)*D*Y(2)*Y(265)+RATE(1853)*D&
    &*Y(2)*Y(265)+RATE(1896)*D*Y(4)*Y(265)+RATE(2190)*D*Y(10)*Y(265)&
    &+RATE(2305)*D*Y(117)*Y(265)+RATE(2415)*D*Y(13)*Y(265)+RATE(2416)*D*Y(13)&
    &*Y(265)+RATE(2544)*D*Y(27)*Y(265)+RATE(2895)*D*Y(46)*Y(265)+RATE(2896)*D&
    &*Y(46)*Y(265)
    PROD = RATE(265)*Y(259)/safeMantle+RATE(348)*D*Y(259)/safeMantle*Y(2&
    &)+RATE(431)*Y(259)/safeMantle+RATE(1113)*Y(259)+RATE(1196)*Y(271)&
    &+RATE(1760)*D*Y(82)*Y(163)+RATE(2459)*D*Y(275)*Y(333)+RATE(2536)*D*Y(27)&
    &*Y(173)+RATE(2549)*D*Y(27)*Y(316)+RATE(2550)*D*Y(27)*Y(277)+RATE(2660)*D&
    &*Y(35)*Y(163)+RATE(2820)*D*Y(133)*Y(163)
    YDOT(265) = PROD-LOSS
    LOSS = RATE(578)*D*Y(266)+RATE(2823)*D*Y(266)*Y(333)+RATE(2846)*D&
    &*Y(46)*Y(266)
    PROD = RATE(1293)*D*Y(18)*Y(265)+RATE(1896)*D*Y(4)*Y(265)+RATE(2500)&
    &*D*Y(27)*Y(184)+RATE(2501)*D*Y(27)*Y(174)+RATE(2506)*D*Y(27)*Y(278)&
    &+RATE(2641)*D*Y(35)*Y(165)+RATE(2696)*D*Y(36)*Y(163)
    YDOT(266) = PROD-LOSS
    LOSS = RATE(157)*Y(267)+RATE(538)*D*Y(267)+RATE(906)*Y(267)&
    &+RATE(2175)*D*Y(10)*Y(267)+RATE(2176)*D*Y(10)*Y(267)+RATE(2970)*D*Y(162)&
    &*Y(267)
    PROD = RATE(246)*Y(268)/safeMantle+RATE(329)*D*Y(268)/safeMantle*Y(2&
    &)+RATE(412)*Y(268)/safeMantle+RATE(754)*Y(118)*Y(58)+RATE(806)*Y(126)&
    &*Y(60)*bulkLayersReciprocal+RATE(1094)*Y(268)+RATE(1177)*Y(273)&
    &+RATE(3004)*D*Y(56)*Y(131)
    YDOT(267) = PROD-LOSS
    LOSS = RATE(94)*Y(268)+RATE(246)*Y(268)/safeMantle+RATE(329)*D*Y(268&
    &)/safeMantle*Y(2)+RATE(412)*Y(268)/safeMantle+RATE(1011)*Y(268)&
    &*totalSwap/safeMantle+RATE(1094)*Y(268)
    PROD = RATE(44)*Y(273)*bulkLayersReciprocal+RATE(538)*D*Y(267)&
    &+RATE(650)*Y(118)*Y(58)
    YDOT(268) = PROD-LOSS
    LOSS = RATE(33)*Y(269)*bulkLayersReciprocal+RATE(1166)*Y(269)
    PROD = RATE(703)*Y(253)*Y(5)*bulkLayersReciprocal+RATE(1000)*Y(257)&
    &*totalSwap/safeMantle
    YDOT(269) = PROD-LOSS
    LOSS = RATE(62)*Y(270)*bulkLayersReciprocal+RATE(1195)*Y(270)
    PROD = RATE(1029)*Y(258)*totalSwap/safeMantle
    YDOT(270) = PROD-LOSS
    LOSS = RATE(63)*Y(271)*bulkLayersReciprocal+RATE(1196)*Y(271)
    PROD = RATE(1030)*Y(259)*totalSwap/safeMantle
    YDOT(271) = PROD-LOSS
    LOSS = RATE(37)*Y(272)*bulkLayersReciprocal+RATE(1170)*Y(272)
    PROD = RATE(706)*Y(254)*Y(5)*bulkLayersReciprocal+RATE(1004)*Y(260)&
    &*totalSwap/safeMantle
    YDOT(272) = PROD-LOSS
    LOSS = RATE(44)*Y(273)*bulkLayersReciprocal+RATE(1177)*Y(273)
    PROD = RATE(702)*Y(126)*Y(60)*bulkLayersReciprocal+RATE(1011)*Y(268)&
    &*totalSwap/safeMantle
    YDOT(273) = PROD-LOSS
    LOSS = RATE(525)*D*Y(274)+RATE(2214)*D*Y(274)*Y(333)+RATE(2215)*D&
    &*Y(274)*Y(333)
    PROD = RATE(1607)*D*Y(25)*Y(183)+RATE(1654)*D*Y(33)*Y(183)+RATE(1661&
    &)*D*Y(33)*Y(300)+RATE(1691)*D*Y(41)*Y(174)+RATE(1698)*D*Y(41)*Y(165)&
    &+RATE(2170)*D*Y(10)*Y(261)+RATE(2300)*D*Y(117)*Y(261)+RATE(3113)*D*Y(278&
    &)*Y(109)
    YDOT(274) = PROD-LOSS
    LOSS = RATE(548)*D*Y(275)+RATE(2459)*D*Y(275)*Y(333)
    PROD = RATE(2190)*D*Y(10)*Y(265)+RATE(2305)*D*Y(117)*Y(265)&
    &+RATE(2753)*D*Y(44)*Y(163)
    YDOT(275) = PROD-LOSS
    LOSS = RATE(284)*Y(276)/safeMantle+RATE(367)*D*Y(276)/safeMantle*Y(2&
    &)+RATE(450)*Y(276)/safeMantle+RATE(1049)*Y(276)*totalSwap/safeMantle&
    &+RATE(1132)*Y(276)
    PROD = RATE(82)*Y(279)*bulkLayersReciprocal+RATE(556)*D*Y(281)&
    &+RATE(616)*D*Y(277)+RATE(617)*D*Y(278)
    YDOT(276) = PROD-LOSS
    LOSS = RATE(200)*Y(277)+RATE(201)*Y(277)+RATE(616)*D*Y(277)+RATE(965&
    &)*Y(277)+RATE(966)*Y(277)+RATE(1275)*D*Y(16)*Y(277)+RATE(1276)*D*Y(16)&
    &*Y(277)+RATE(1295)*D*Y(18)*Y(277)+RATE(1335)*D*Y(18)*Y(277)+RATE(1336)*D&
    &*Y(18)*Y(277)+RATE(1337)*D*Y(18)*Y(277)+RATE(1497)*D*Y(20)*Y(277)&
    &+RATE(1498)*D*Y(20)*Y(277)+RATE(1664)*D*Y(33)*Y(277)+RATE(1864)*D*Y(2)&
    &*Y(277)+RATE(1865)*D*Y(2)*Y(277)+RATE(1904)*D*Y(4)*Y(277)+RATE(2199)*D&
    &*Y(10)*Y(277)+RATE(2309)*D*Y(117)*Y(277)+RATE(2428)*D*Y(13)*Y(277)&
    &+RATE(2429)*D*Y(13)*Y(277)+RATE(2550)*D*Y(27)*Y(277)+RATE(2551)*D*Y(27)&
    &*Y(277)+RATE(2905)*D*Y(46)*Y(277)+RATE(2916)*D*Y(46)*Y(277)+RATE(2963)*D&
    &*Y(161)*Y(277)+RATE(3016)*D*Y(56)*Y(277)+RATE(3061)*D*Y(163)*Y(277)
    PROD = RATE(202)*Y(320)+RATE(284)*Y(276)/safeMantle+RATE(367)*D&
    &*Y(276)/safeMantle*Y(2)+RATE(450)*Y(276)/safeMantle+RATE(967)*Y(320)&
    &+RATE(1132)*Y(276)+RATE(1215)*Y(279)+RATE(1274)*D*Y(16)*Y(320)+RATE(2471&
    &)*D*Y(281)*Y(333)+RATE(2473)*D*Y(326)*Y(333)+RATE(2474)*D*Y(326)*Y(333)&
    &+RATE(2485)*D*Y(69)*Y(278)+RATE(2766)*D*Y(54)*Y(278)+RATE(2821)*D*Y(133)&
    &*Y(163)+RATE(2872)*D*Y(46)*Y(233)+RATE(2888)*D*Y(46)*Y(173)+RATE(2896)*D&
    &*Y(46)*Y(265)+RATE(2901)*D*Y(46)*Y(300)+RATE(2903)*D*Y(46)*Y(316)&
    &+RATE(2904)*D*Y(46)*Y(320)+RATE(2959)*D*Y(161)*Y(303)+RATE(2962)*D*Y(161&
    &)*Y(163)+RATE(2974)*D*Y(301)*Y(333)+RATE(3015)*D*Y(56)*Y(163)+RATE(3056)&
    &*D*Y(163)*Y(236)+RATE(3060)*D*Y(163)*Y(320)+RATE(3060)*D*Y(163)*Y(320)&
    &+RATE(3116)*D*Y(319)*Y(333)
    YDOT(277) = PROD-LOSS
    LOSS = RATE(617)*D*Y(278)+RATE(2140)*D*Y(183)*Y(278)+RATE(2485)*D&
    &*Y(69)*Y(278)+RATE(2506)*D*Y(27)*Y(278)+RATE(2766)*D*Y(54)*Y(278)&
    &+RATE(3107)*D*Y(278)*Y(333)+RATE(3108)*D*Y(278)*Y(80)+RATE(3109)*D*Y(278&
    &)*Y(80)+RATE(3110)*D*Y(278)*Y(80)+RATE(3111)*D*Y(278)*Y(109)+RATE(3112)&
    &*D*Y(278)*Y(109)+RATE(3113)*D*Y(278)*Y(109)+RATE(3114)*D*Y(278)*Y(300)
    PROD = RATE(201)*Y(277)+RATE(966)*Y(277)+RATE(1295)*D*Y(18)*Y(277)&
    &+RATE(1334)*D*Y(18)*Y(320)+RATE(1780)*D*Y(99)*Y(319)+RATE(1796)*D*Y(100)&
    &*Y(320)+RATE(1819)*D*Y(2)*Y(319)+RATE(1904)*D*Y(4)*Y(277)+RATE(2427)*D&
    &*Y(13)*Y(320)+RATE(2836)*D*Y(46)*Y(184)+RATE(2841)*D*Y(46)*Y(174)&
    &+RATE(2952)*D*Y(48)*Y(320)+RATE(2957)*D*Y(161)*Y(165)+RATE(2958)*D*Y(161&
    &)*Y(303)+RATE(2971)*D*Y(162)*Y(163)+RATE(2991)*D*Y(56)*Y(165)+RATE(3047)&
    &*D*Y(57)*Y(163)
    YDOT(278) = PROD-LOSS
    LOSS = RATE(82)*Y(279)*bulkLayersReciprocal+RATE(1215)*Y(279)
    PROD = RATE(1049)*Y(276)*totalSwap/safeMantle
    YDOT(279) = PROD-LOSS
    LOSS = RATE(122)*Y(280)+RATE(473)*D*Y(280)+RATE(845)*Y(280)&
    &+RATE(2347)*D*Y(13)*Y(280)+RATE(2519)*D*Y(27)*Y(280)
    PROD = RATE(213)*Y(282)/safeMantle+RATE(296)*D*Y(282)/safeMantle*Y(2&
    &)+RATE(379)*Y(282)/safeMantle+RATE(1061)*Y(282)+RATE(1144)*Y(283)&
    &+RATE(1245)*D*Y(16)*Y(200)+RATE(1361)*D*Y(67)*Y(80)
    YDOT(280) = PROD-LOSS
    LOSS = RATE(556)*D*Y(281)+RATE(2471)*D*Y(281)*Y(333)
    PROD = RATE(2134)*D*Y(63)*Y(163)+RATE(2199)*D*Y(10)*Y(277)+RATE(2309&
    &)*D*Y(117)*Y(277)
    YDOT(281) = PROD-LOSS
    LOSS = RATE(213)*Y(282)/safeMantle+RATE(296)*D*Y(282)/safeMantle*Y(2&
    &)+RATE(379)*Y(282)/safeMantle+RATE(978)*Y(282)*totalSwap/safeMantle&
    &+RATE(1061)*Y(282)
    PROD = RATE(11)*Y(283)*bulkLayersReciprocal+RATE(473)*D*Y(280)
    YDOT(282) = PROD-LOSS
    LOSS = RATE(11)*Y(283)*bulkLayersReciprocal+RATE(1144)*Y(283)
    PROD = RATE(978)*Y(282)*totalSwap/safeMantle
    YDOT(283) = PROD-LOSS
    LOSS = RATE(121)*Y(284)+RATE(472)*D*Y(284)+RATE(844)*Y(284)&
    &+RATE(2346)*D*Y(13)*Y(284)+RATE(2518)*D*Y(27)*Y(284)+RATE(2865)*D*Y(46)&
    &*Y(284)
    PROD = RATE(212)*Y(285)/safeMantle+RATE(295)*D*Y(285)/safeMantle*Y(2&
    &)+RATE(378)*Y(285)/safeMantle+RATE(1060)*Y(285)+RATE(1143)*Y(286)&
    &+RATE(1362)*D*Y(67)*Y(90)+RATE(1430)*D*Y(309)*Y(333)+RATE(2520)*D*Y(27)&
    &*Y(308)+RATE(2866)*D*Y(46)*Y(308)
    YDOT(284) = PROD-LOSS
    LOSS = RATE(212)*Y(285)/safeMantle+RATE(295)*D*Y(285)/safeMantle*Y(2&
    &)+RATE(378)*Y(285)/safeMantle+RATE(977)*Y(285)*totalSwap/safeMantle&
    &+RATE(1060)*Y(285)
    PROD = RATE(10)*Y(286)*bulkLayersReciprocal
    YDOT(285) = PROD-LOSS
    LOSS = RATE(10)*Y(286)*bulkLayersReciprocal+RATE(1143)*Y(286)
    PROD = RATE(977)*Y(285)*totalSwap/safeMantle
    YDOT(286) = PROD-LOSS
    LOSS = RATE(241)*Y(287)/safeMantle+RATE(324)*D*Y(287)/safeMantle*Y(2&
    &)+RATE(407)*Y(287)/safeMantle+RATE(1006)*Y(287)*totalSwap/safeMantle&
    &+RATE(1089)*Y(287)
    PROD = RATE(39)*Y(289)*bulkLayersReciprocal+RATE(472)*D*Y(284)&
    &+RATE(528)*D*Y(288)
    YDOT(287) = PROD-LOSS
    LOSS = RATE(152)*Y(288)+RATE(528)*D*Y(288)+RATE(899)*Y(288)&
    &+RATE(1319)*D*Y(18)*Y(288)+RATE(1320)*D*Y(18)*Y(288)+RATE(1321)*D*Y(18)&
    &*Y(288)+RATE(2386)*D*Y(13)*Y(288)+RATE(2387)*D*Y(13)*Y(288)
    PROD = RATE(241)*Y(287)/safeMantle+RATE(324)*D*Y(287)/safeMantle*Y(2&
    &)+RATE(407)*Y(287)/safeMantle+RATE(1089)*Y(287)+RATE(1172)*Y(289)&
    &+RATE(1384)*D*Y(75)*Y(90)+RATE(1385)*D*Y(75)*Y(92)+RATE(1386)*D*Y(75)&
    &*Y(291)+RATE(1388)*D*Y(75)*Y(82)+RATE(1747)*D*Y(82)*Y(80)+RATE(2084)*D&
    &*Y(62)*Y(309)+RATE(2517)*D*Y(27)*Y(200)
    YDOT(288) = PROD-LOSS
    LOSS = RATE(39)*Y(289)*bulkLayersReciprocal+RATE(1172)*Y(289)
    PROD = RATE(1006)*Y(287)*totalSwap/safeMantle
    YDOT(289) = PROD-LOSS
    LOSS = RATE(467)*D*Y(290)+RATE(1399)*D*Y(80)*Y(290)+RATE(1424)*D&
    &*Y(290)*Y(333)+RATE(1425)*D*Y(290)*Y(333)+RATE(1809)*D*Y(2)*Y(290)&
    &+RATE(2243)*D*Y(90)*Y(290)
    PROD = RATE(918)*Y(291)+RATE(1771)*D*Y(83)*Y(90)+RATE(2588)*D*Y(28)&
    &*Y(291)
    YDOT(290) = PROD-LOSS
    LOSS = RATE(168)*Y(291)+RATE(565)*D*Y(291)+RATE(918)*Y(291)+RATE(919&
    &)*Y(291)+RATE(1258)*D*Y(16)*Y(291)+RATE(1290)*D*Y(18)*Y(291)+RATE(1386)&
    &*D*Y(75)*Y(291)+RATE(1845)*D*Y(2)*Y(291)+RATE(2407)*D*Y(13)*Y(291)&
    &+RATE(2538)*D*Y(27)*Y(291)+RATE(2587)*D*Y(28)*Y(291)+RATE(2588)*D*Y(28)&
    &*Y(291)+RATE(3010)*D*Y(56)*Y(291)
    PROD = RATE(258)*Y(294)/safeMantle+RATE(341)*D*Y(294)/safeMantle*Y(2&
    &)+RATE(424)*Y(294)/safeMantle+RATE(1106)*Y(294)+RATE(1189)*Y(296)&
    &+RATE(1399)*D*Y(80)*Y(290)+RATE(1751)*D*Y(82)*Y(90)+RATE(1753)*D*Y(82)&
    &*Y(92)+RATE(2243)*D*Y(90)*Y(290)
    YDOT(291) = PROD-LOSS
    LOSS = RATE(192)*Y(292)+RATE(598)*D*Y(292)+RATE(948)*Y(292)&
    &+RATE(1297)*D*Y(18)*Y(292)+RATE(1906)*D*Y(4)*Y(292)+RATE(2430)*D*Y(13)&
    &*Y(292)+RATE(2906)*D*Y(46)*Y(292)
    PROD = RATE(193)*Y(314)+RATE(276)*Y(295)/safeMantle+RATE(359)*D&
    &*Y(295)/safeMantle*Y(2)+RATE(442)*Y(295)/safeMantle+RATE(950)*Y(314)&
    &+RATE(1124)*Y(295)+RATE(1207)*Y(297)+RATE(1401)*D*Y(80)*Y(105)+RATE(2907&
    &)*D*Y(46)*Y(314)+RATE(3087)*D*Y(313)*Y(333)
    YDOT(292) = PROD-LOSS
    LOSS = RATE(599)*D*Y(293)+RATE(3085)*D*Y(293)*Y(333)+RATE(3086)*D&
    &*Y(293)*Y(333)
    PROD = RATE(1297)*D*Y(18)*Y(292)+RATE(1383)*D*Y(75)*Y(106)+RATE(1414&
    &)*D*Y(81)*Y(105)+RATE(1906)*D*Y(4)*Y(292)+RATE(2431)*D*Y(13)*Y(314)
    YDOT(293) = PROD-LOSS
    LOSS = RATE(258)*Y(294)/safeMantle+RATE(341)*D*Y(294)/safeMantle*Y(2&
    &)+RATE(424)*Y(294)/safeMantle+RATE(1023)*Y(294)*totalSwap/safeMantle&
    &+RATE(1106)*Y(294)
    PROD = RATE(56)*Y(296)*bulkLayersReciprocal+RATE(467)*D*Y(290)&
    &+RATE(565)*D*Y(291)
    YDOT(294) = PROD-LOSS
    LOSS = RATE(276)*Y(295)/safeMantle+RATE(359)*D*Y(295)/safeMantle*Y(2&
    &)+RATE(442)*Y(295)/safeMantle+RATE(1041)*Y(295)*totalSwap/safeMantle&
    &+RATE(1124)*Y(295)
    PROD = RATE(74)*Y(297)*bulkLayersReciprocal+RATE(598)*D*Y(292)&
    &+RATE(599)*D*Y(293)
    YDOT(295) = PROD-LOSS
    LOSS = RATE(56)*Y(296)*bulkLayersReciprocal+RATE(1189)*Y(296)
    PROD = RATE(1023)*Y(294)*totalSwap/safeMantle
    YDOT(296) = PROD-LOSS
    LOSS = RATE(74)*Y(297)*bulkLayersReciprocal+RATE(1207)*Y(297)
    PROD = RATE(1041)*Y(295)*totalSwap/safeMantle
    YDOT(297) = PROD-LOSS
    LOSS = RATE(283)*Y(298)/safeMantle+RATE(366)*D*Y(298)/safeMantle*Y(2&
    &)+RATE(449)*Y(298)/safeMantle+RATE(1048)*Y(298)*totalSwap/safeMantle&
    &+RATE(1131)*Y(298)
    PROD = RATE(81)*Y(304)*bulkLayersReciprocal+RATE(555)*D*Y(307)&
    &+RATE(614)*D*Y(302)+RATE(615)*D*Y(303)
    YDOT(298) = PROD-LOSS
    LOSS = RATE(270)*Y(299)/safeMantle+RATE(353)*D*Y(299)/safeMantle*Y(2&
    &)+RATE(436)*Y(299)/safeMantle+RATE(1035)*Y(299)*totalSwap/safeMantle&
    &+RATE(1118)*Y(299)
    PROD = RATE(68)*Y(305)*bulkLayersReciprocal+RATE(550)*D*Y(306)&
    &+RATE(586)*D*Y(300)+RATE(587)*D*Y(301)
    YDOT(299) = PROD-LOSS
    LOSS = RATE(185)*Y(300)+RATE(186)*Y(300)+RATE(586)*D*Y(300)+RATE(938&
    &)*Y(300)+RATE(939)*Y(300)+RATE(1270)*D*Y(16)*Y(300)+RATE(1294)*D*Y(18)&
    &*Y(300)+RATE(1332)*D*Y(18)*Y(300)+RATE(1493)*D*Y(20)*Y(300)+RATE(1537)*D&
    &*Y(21)*Y(300)+RATE(1538)*D*Y(21)*Y(300)+RATE(1614)*D*Y(25)*Y(300)&
    &+RATE(1615)*D*Y(25)*Y(300)+RATE(1661)*D*Y(33)*Y(300)+RATE(1709)*D*Y(42)&
    &*Y(300)+RATE(1720)*D*Y(42)*Y(300)+RATE(1861)*D*Y(2)*Y(300)+RATE(1899)*D&
    &*Y(4)*Y(300)+RATE(1937)*D*Y(4)*Y(300)+RATE(2116)*D*Y(63)*Y(300)&
    &+RATE(2194)*D*Y(10)*Y(300)+RATE(2306)*D*Y(117)*Y(300)+RATE(2420)*D*Y(13)&
    &*Y(300)+RATE(2421)*D*Y(13)*Y(300)+RATE(2422)*D*Y(13)*Y(300)+RATE(2423)*D&
    &*Y(13)*Y(300)+RATE(2570)*D*Y(28)*Y(300)+RATE(2595)*D*Y(28)*Y(300)&
    &+RATE(2596)*D*Y(28)*Y(300)+RATE(2606)*D*Y(104)*Y(300)+RATE(2613)*D*Y(104&
    &)*Y(300)+RATE(2900)*D*Y(46)*Y(300)+RATE(2901)*D*Y(46)*Y(300)+RATE(2931)&
    &*D*Y(48)*Y(300)+RATE(2950)*D*Y(48)*Y(300)+RATE(3066)*D*Y(165)*Y(300)&
    &+RATE(3082)*D*Y(106)*Y(300)+RATE(3114)*D*Y(278)*Y(300)
    PROD = RATE(270)*Y(299)/safeMantle+RATE(353)*D*Y(299)/safeMantle*Y(2&
    &)+RATE(436)*Y(299)/safeMantle+RATE(1118)*Y(299)+RATE(1201)*Y(305)&
    &+RATE(1498)*D*Y(20)*Y(277)+RATE(2096)*D*Y(62)*Y(306)+RATE(2136)*D*Y(183)&
    &*Y(301)+RATE(2462)*D*Y(306)*Y(333)+RATE(2764)*D*Y(54)*Y(301)+RATE(2883)&
    &*D*Y(46)*Y(244)+RATE(3002)*D*Y(56)*Y(233)+RATE(3057)*D*Y(163)*Y(116)
    YDOT(300) = PROD-LOSS
    LOSS = RATE(587)*D*Y(301)+RATE(2136)*D*Y(183)*Y(301)+RATE(2764)*D&
    &*Y(54)*Y(301)+RATE(2974)*D*Y(301)*Y(333)+RATE(2975)*D*Y(301)*Y(333)&
    &+RATE(2976)*D*Y(301)*Y(333)
    PROD = RATE(185)*Y(300)+RATE(938)*Y(300)+RATE(1294)*D*Y(18)*Y(300)&
    &+RATE(1709)*D*Y(42)*Y(300)+RATE(1899)*D*Y(4)*Y(300)+RATE(2116)*D*Y(63)&
    &*Y(300)+RATE(2570)*D*Y(28)*Y(300)+RATE(2606)*D*Y(104)*Y(300)+RATE(2838)&
    &*D*Y(46)*Y(245)+RATE(2931)*D*Y(48)*Y(300)+RATE(2956)*D*Y(161)*Y(234)
    YDOT(301) = PROD-LOSS
    LOSS = RATE(199)*Y(302)+RATE(614)*D*Y(302)+RATE(964)*Y(302)&
    &+RATE(1302)*D*Y(18)*Y(302)+RATE(1342)*D*Y(18)*Y(302)+RATE(1914)*D*Y(4)&
    &*Y(302)+RATE(2206)*D*Y(10)*Y(302)+RATE(2314)*D*Y(117)*Y(302)+RATE(2443)&
    &*D*Y(13)*Y(302)+RATE(2444)*D*Y(13)*Y(302)+RATE(3063)*D*Y(165)*Y(302)
    PROD = RATE(283)*Y(298)/safeMantle+RATE(366)*D*Y(298)/safeMantle*Y(2&
    &)+RATE(449)*Y(298)/safeMantle+RATE(1131)*Y(298)+RATE(1214)*Y(304)&
    &+RATE(2139)*D*Y(183)*Y(307)+RATE(2252)*D*Y(90)*Y(307)+RATE(2470)*D*Y(307&
    &)*Y(333)+RATE(2787)*D*Y(54)*Y(307)
    YDOT(302) = PROD-LOSS
    LOSS = RATE(615)*D*Y(303)+RATE(1821)*D*Y(2)*Y(303)+RATE(2958)*D&
    &*Y(161)*Y(303)+RATE(2959)*D*Y(161)*Y(303)+RATE(3106)*D*Y(303)*Y(333)
    PROD = RATE(1302)*D*Y(18)*Y(302)+RATE(1914)*D*Y(4)*Y(302)+RATE(3063)&
    &*D*Y(165)*Y(302)+RATE(3082)*D*Y(106)*Y(300)+RATE(3090)*D*Y(121)*Y(165)
    YDOT(303) = PROD-LOSS
    LOSS = RATE(81)*Y(304)*bulkLayersReciprocal+RATE(1214)*Y(304)
    PROD = RATE(1048)*Y(298)*totalSwap/safeMantle
    YDOT(304) = PROD-LOSS
    LOSS = RATE(68)*Y(305)*bulkLayersReciprocal+RATE(1201)*Y(305)
    PROD = RATE(1035)*Y(299)*totalSwap/safeMantle
    YDOT(305) = PROD-LOSS
    LOSS = RATE(550)*D*Y(306)+RATE(2096)*D*Y(62)*Y(306)+RATE(2461)*D&
    &*Y(306)*Y(333)+RATE(2462)*D*Y(306)*Y(333)
    PROD = RATE(1538)*D*Y(21)*Y(300)+RATE(1664)*D*Y(33)*Y(277)+RATE(1720&
    &)*D*Y(42)*Y(300)+RATE(2194)*D*Y(10)*Y(300)+RATE(2306)*D*Y(117)*Y(300)
    YDOT(306) = PROD-LOSS
    LOSS = RATE(555)*D*Y(307)+RATE(2099)*D*Y(62)*Y(307)+RATE(2139)*D&
    &*Y(183)*Y(307)+RATE(2252)*D*Y(90)*Y(307)+RATE(2469)*D*Y(307)*Y(333)&
    &+RATE(2470)*D*Y(307)*Y(333)+RATE(2787)*D*Y(54)*Y(307)
    PROD = RATE(2206)*D*Y(10)*Y(302)+RATE(2314)*D*Y(117)*Y(302)&
    &+RATE(3096)*D*Y(136)*Y(163)
    YDOT(307) = PROD-LOSS
    LOSS = RATE(474)*D*Y(308)+RATE(2520)*D*Y(27)*Y(308)+RATE(2866)*D&
    &*Y(46)*Y(308)
    PROD = RATE(214)*Y(310)/safeMantle+RATE(297)*D*Y(310)/safeMantle*Y(2&
    &)+RATE(380)*Y(310)/safeMantle+RATE(1062)*Y(310)+RATE(1145)*Y(311)&
    &+RATE(2519)*D*Y(27)*Y(280)
    YDOT(308) = PROD-LOSS
    LOSS = RATE(475)*D*Y(309)+RATE(1429)*D*Y(309)*Y(333)+RATE(1430)*D&
    &*Y(309)*Y(333)+RATE(2084)*D*Y(62)*Y(309)
    PROD = RATE(1321)*D*Y(18)*Y(288)+RATE(2244)*D*Y(90)*Y(192)
    YDOT(309) = PROD-LOSS
    LOSS = RATE(214)*Y(310)/safeMantle+RATE(297)*D*Y(310)/safeMantle*Y(2&
    &)+RATE(380)*Y(310)/safeMantle+RATE(979)*Y(310)*totalSwap/safeMantle&
    &+RATE(1062)*Y(310)
    PROD = RATE(12)*Y(311)*bulkLayersReciprocal+RATE(474)*D*Y(308)&
    &+RATE(475)*D*Y(309)
    YDOT(310) = PROD-LOSS
    LOSS = RATE(12)*Y(311)*bulkLayersReciprocal+RATE(1145)*Y(311)
    PROD = RATE(979)*Y(310)*totalSwap/safeMantle
    YDOT(311) = PROD-LOSS
    LOSS = RATE(277)*Y(312)/safeMantle+RATE(360)*D*Y(312)/safeMantle*Y(2&
    &)+RATE(443)*Y(312)/safeMantle+RATE(1042)*Y(312)*totalSwap/safeMantle&
    &+RATE(1125)*Y(312)
    PROD = RATE(75)*Y(321)*bulkLayersReciprocal+RATE(600)*D*Y(314)&
    &+RATE(601)*D*Y(313)
    YDOT(312) = PROD-LOSS
    LOSS = RATE(601)*D*Y(313)+RATE(3087)*D*Y(313)*Y(333)+RATE(3088)*D&
    &*Y(313)*Y(333)
    PROD = RATE(1298)*D*Y(18)*Y(314)+RATE(1907)*D*Y(4)*Y(314)
    YDOT(313) = PROD-LOSS
    LOSS = RATE(193)*Y(314)+RATE(600)*D*Y(314)+RATE(949)*Y(314)+RATE(950&
    &)*Y(314)+RATE(1298)*D*Y(18)*Y(314)+RATE(1907)*D*Y(4)*Y(314)+RATE(2431)*D&
    &*Y(13)*Y(314)+RATE(2907)*D*Y(46)*Y(314)
    PROD = RATE(277)*Y(312)/safeMantle+RATE(360)*D*Y(312)/safeMantle*Y(2&
    &)+RATE(443)*Y(312)/safeMantle+RATE(1125)*Y(312)+RATE(1208)*Y(321)
    YDOT(314) = PROD-LOSS
    LOSS = RATE(285)*Y(315)/safeMantle+RATE(368)*D*Y(315)/safeMantle*Y(2&
    &)+RATE(451)*Y(315)/safeMantle+RATE(1050)*Y(315)*totalSwap/safeMantle&
    &+RATE(1133)*Y(315)
    PROD = RATE(83)*Y(322)*bulkLayersReciprocal+RATE(557)*D*Y(326)&
    &+RATE(618)*D*Y(320)+RATE(619)*D*Y(319)
    YDOT(315) = PROD-LOSS
    LOSS = RATE(189)*Y(316)+RATE(592)*D*Y(316)+RATE(944)*Y(316)+RATE(945&
    &)*Y(316)+RATE(1273)*D*Y(16)*Y(316)+RATE(1863)*D*Y(2)*Y(316)+RATE(1901)*D&
    &*Y(4)*Y(316)+RATE(2196)*D*Y(10)*Y(316)+RATE(2230)*D*Y(66)*Y(316)&
    &+RATE(2307)*D*Y(117)*Y(316)+RATE(2425)*D*Y(13)*Y(316)+RATE(2549)*D*Y(27)&
    &*Y(316)+RATE(2903)*D*Y(46)*Y(316)
    PROD = RATE(273)*Y(317)/safeMantle+RATE(356)*D*Y(317)/safeMantle*Y(2&
    &)+RATE(439)*Y(317)/safeMantle+RATE(1121)*Y(317)+RATE(1204)*Y(323)&
    &+RATE(2138)*D*Y(183)*Y(325)+RATE(2468)*D*Y(325)*Y(333)+RATE(2814)*D&
    &*Y(133)*Y(318)+RATE(3059)*D*Y(163)*Y(173)+RATE(3061)*D*Y(163)*Y(277)
    YDOT(316) = PROD-LOSS
    LOSS = RATE(273)*Y(317)/safeMantle+RATE(356)*D*Y(317)/safeMantle*Y(2&
    &)+RATE(439)*Y(317)/safeMantle+RATE(667)*Y(317)*Y(3)+RATE(771)*Y(317)*Y(3&
    &)+RATE(1038)*Y(317)*totalSwap/safeMantle+RATE(1121)*Y(317)
    PROD = RATE(71)*Y(323)*bulkLayersReciprocal+RATE(592)*D*Y(316)&
    &+RATE(593)*D*Y(318)
    YDOT(317) = PROD-LOSS
    LOSS = RATE(593)*D*Y(318)+RATE(1672)*D*Y(159)*Y(318)+RATE(2814)*D&
    &*Y(133)*Y(318)+RATE(3068)*D*Y(318)*Y(333)
    PROD = RATE(944)*Y(316)+RATE(1901)*D*Y(4)*Y(316)+RATE(2140)*D*Y(183)&
    &*Y(278)+RATE(2404)*D*Y(13)*Y(324)+RATE(3065)*D*Y(165)*Y(183)+RATE(3066)&
    &*D*Y(165)*Y(300)+RATE(3114)*D*Y(278)*Y(300)
    YDOT(318) = PROD-LOSS
    LOSS = RATE(619)*D*Y(319)+RATE(1780)*D*Y(99)*Y(319)+RATE(1819)*D*Y(2&
    &)*Y(319)+RATE(1981)*D*Y(6)*Y(319)+RATE(2955)*D*Y(161)*Y(319)+RATE(3115)&
    &*D*Y(319)*Y(333)+RATE(3116)*D*Y(319)*Y(333)
    PROD = RATE(1903)*D*Y(4)*Y(320)+RATE(2331)*D*Y(13)*Y(320)+RATE(2933)&
    &*D*Y(48)*Y(320)
    YDOT(319) = PROD-LOSS
    LOSS = RATE(202)*Y(320)+RATE(618)*D*Y(320)+RATE(967)*Y(320)&
    &+RATE(1274)*D*Y(16)*Y(320)+RATE(1334)*D*Y(18)*Y(320)+RATE(1796)*D*Y(100)&
    &*Y(320)+RATE(1903)*D*Y(4)*Y(320)+RATE(2135)*D*Y(63)*Y(320)+RATE(2198)*D&
    &*Y(10)*Y(320)+RATE(2331)*D*Y(13)*Y(320)+RATE(2426)*D*Y(13)*Y(320)&
    &+RATE(2427)*D*Y(13)*Y(320)+RATE(2904)*D*Y(46)*Y(320)+RATE(2933)*D*Y(48)&
    &*Y(320)+RATE(2952)*D*Y(48)*Y(320)+RATE(3060)*D*Y(163)*Y(320)
    PROD = RATE(285)*Y(315)/safeMantle+RATE(368)*D*Y(315)/safeMantle*Y(2&
    &)+RATE(451)*Y(315)/safeMantle+RATE(1133)*Y(315)+RATE(1216)*Y(322)&
    &+RATE(2098)*D*Y(62)*Y(326)+RATE(2472)*D*Y(326)*Y(333)+RATE(2786)*D*Y(54)&
    &*Y(326)+RATE(2916)*D*Y(46)*Y(277)+RATE(2955)*D*Y(161)*Y(319)+RATE(2963)&
    &*D*Y(161)*Y(277)+RATE(3016)*D*Y(56)*Y(277)
    YDOT(320) = PROD-LOSS
    LOSS = RATE(75)*Y(321)*bulkLayersReciprocal+RATE(1208)*Y(321)
    PROD = RATE(1042)*Y(312)*totalSwap/safeMantle
    YDOT(321) = PROD-LOSS
    LOSS = RATE(83)*Y(322)*bulkLayersReciprocal+RATE(1216)*Y(322)
    PROD = RATE(1050)*Y(315)*totalSwap/safeMantle
    YDOT(322) = PROD-LOSS
    LOSS = RATE(71)*Y(323)*bulkLayersReciprocal+RATE(719)*Y(323)*Y(5)&
    &*bulkLayersReciprocal+RATE(823)*Y(323)*Y(5)*bulkLayersReciprocal&
    &+RATE(1204)*Y(323)
    PROD = RATE(1038)*Y(317)*totalSwap/safeMantle
    YDOT(323) = PROD-LOSS
    LOSS = RATE(164)*Y(324)+RATE(553)*D*Y(324)+RATE(914)*Y(324)+RATE(915&
    &)*Y(324)+RATE(1889)*D*Y(4)*Y(324)+RATE(2181)*D*Y(10)*Y(324)+RATE(2229)*D&
    &*Y(66)*Y(324)+RATE(2303)*D*Y(117)*Y(324)+RATE(2403)*D*Y(13)*Y(324)&
    &+RATE(2404)*D*Y(13)*Y(324)
    PROD = RATE(253)*Y(327)/safeMantle+RATE(336)*D*Y(327)/safeMantle*Y(2&
    &)+RATE(419)*Y(327)/safeMantle+RATE(771)*Y(317)*Y(3)+RATE(823)*Y(323)*Y(5&
    &)*bulkLayersReciprocal+RATE(1101)*Y(327)+RATE(1184)*Y(328)+RATE(2145)*D&
    &*Y(329)*Y(333)
    YDOT(324) = PROD-LOSS
    LOSS = RATE(554)*D*Y(325)+RATE(2138)*D*Y(183)*Y(325)+RATE(2467)*D&
    &*Y(325)*Y(333)+RATE(2468)*D*Y(325)*Y(333)
    PROD = RATE(914)*Y(324)+RATE(1889)*D*Y(4)*Y(324)+RATE(2196)*D*Y(10)&
    &*Y(316)+RATE(2230)*D*Y(66)*Y(316)+RATE(2307)*D*Y(117)*Y(316)+RATE(2382)&
    &*D*Y(13)*Y(330)+RATE(2466)*D*Y(174)*Y(183)+RATE(3064)*D*Y(165)*Y(183)
    YDOT(325) = PROD-LOSS
    LOSS = RATE(557)*D*Y(326)+RATE(2098)*D*Y(62)*Y(326)+RATE(2472)*D&
    &*Y(326)*Y(333)+RATE(2473)*D*Y(326)*Y(333)+RATE(2474)*D*Y(326)*Y(333)&
    &+RATE(2786)*D*Y(54)*Y(326)
    PROD = RATE(1981)*D*Y(6)*Y(319)+RATE(2135)*D*Y(63)*Y(320)+RATE(2198)&
    &*D*Y(10)*Y(320)
    YDOT(326) = PROD-LOSS
    LOSS = RATE(253)*Y(327)/safeMantle+RATE(336)*D*Y(327)/safeMantle*Y(2&
    &)+RATE(419)*Y(327)/safeMantle+RATE(653)*Y(327)*Y(3)+RATE(757)*Y(327)*Y(3&
    &)+RATE(1018)*Y(327)*totalSwap/safeMantle+RATE(1101)*Y(327)
    PROD = RATE(51)*Y(328)*bulkLayersReciprocal+RATE(553)*D*Y(324)&
    &+RATE(554)*D*Y(325)+RATE(667)*Y(317)*Y(3)
    YDOT(327) = PROD-LOSS
    LOSS = RATE(51)*Y(328)*bulkLayersReciprocal+RATE(705)*Y(328)*Y(5)&
    &*bulkLayersReciprocal+RATE(809)*Y(328)*Y(5)*bulkLayersReciprocal&
    &+RATE(1184)*Y(328)
    PROD = RATE(719)*Y(323)*Y(5)*bulkLayersReciprocal+RATE(1018)*Y(327)&
    &*totalSwap/safeMantle
    YDOT(328) = PROD-LOSS
    LOSS = RATE(520)*D*Y(329)+RATE(2145)*D*Y(329)*Y(333)+RATE(2146)*D&
    &*Y(329)*Y(333)
    PROD = RATE(1672)*D*Y(159)*Y(318)+RATE(1884)*D*Y(4)*Y(330)+RATE(2181&
    &)*D*Y(10)*Y(324)+RATE(2229)*D*Y(66)*Y(324)+RATE(2303)*D*Y(117)*Y(324)&
    &+RATE(3054)*D*Y(163)*Y(189)
    YDOT(329) = PROD-LOSS
    LOSS = RATE(150)*Y(330)+RATE(519)*D*Y(330)+RATE(894)*Y(330)&
    &+RATE(1884)*D*Y(4)*Y(330)+RATE(2381)*D*Y(13)*Y(330)+RATE(2382)*D*Y(13)&
    &*Y(330)
    PROD = RATE(238)*Y(331)/safeMantle+RATE(321)*D*Y(331)/safeMantle*Y(2&
    &)+RATE(404)*Y(331)/safeMantle+RATE(757)*Y(327)*Y(3)+RATE(809)*Y(328)*Y(5&
    &)*bulkLayersReciprocal+RATE(1086)*Y(331)+RATE(1169)*Y(332)
    YDOT(330) = PROD-LOSS
    LOSS = RATE(238)*Y(331)/safeMantle+RATE(321)*D*Y(331)/safeMantle*Y(2&
    &)+RATE(404)*Y(331)/safeMantle+RATE(1003)*Y(331)*totalSwap/safeMantle&
    &+RATE(1086)*Y(331)
    PROD = RATE(36)*Y(332)*bulkLayersReciprocal+RATE(519)*D*Y(330)&
    &+RATE(520)*D*Y(329)+RATE(653)*Y(327)*Y(3)
    YDOT(331) = PROD-LOSS
    LOSS = RATE(36)*Y(332)*bulkLayersReciprocal+RATE(1169)*Y(332)
    PROD = RATE(705)*Y(328)*Y(5)*bulkLayersReciprocal+RATE(1003)*Y(331)&
    &*totalSwap/safeMantle
    YDOT(332) = PROD-LOSS
    LOSS = RATE(503)*D*Y(333)+RATE(1347)*D*Y(18)*Y(333)+RATE(1368)*D&
    &*Y(68)*Y(333)+RATE(1391)*D*Y(76)*Y(333)+RATE(1392)*D*Y(76)*Y(333)&
    &+RATE(1406)*D*Y(81)*Y(333)+RATE(1407)*D*Y(81)*Y(333)+RATE(1408)*D*Y(81)&
    &*Y(333)+RATE(1422)*D*Y(199)*Y(333)+RATE(1423)*D*Y(199)*Y(333)+RATE(1424)&
    &*D*Y(290)*Y(333)+RATE(1425)*D*Y(290)*Y(333)+RATE(1426)*D*Y(204)*Y(333)&
    &+RATE(1427)*D*Y(192)*Y(333)+RATE(1428)*D*Y(213)*Y(333)+RATE(1429)*D&
    &*Y(309)*Y(333)+RATE(1430)*D*Y(309)*Y(333)+RATE(1505)*D*Y(21)*Y(333)&
    &+RATE(1601)*D*Y(25)*Y(333)+RATE(1602)*D*Y(25)*Y(333)+RATE(1603)*D*Y(25)&
    &*Y(333)+RATE(1646)*D*Y(33)*Y(333)+RATE(1647)*D*Y(33)*Y(333)+RATE(1648)*D&
    &*Y(33)*Y(333)+RATE(1668)*D*Y(33)*Y(333)+RATE(1670)*D*Y(220)*Y(333)&
    &+RATE(1671)*D*Y(220)*Y(333)+RATE(1673)*D*Y(172)*Y(333)+RATE(1674)*D&
    &*Y(172)*Y(333)+RATE(1675)*D*Y(172)*Y(333)+RATE(1676)*D*Y(172)*Y(333)&
    &+RATE(1677)*D*Y(172)*Y(333)+RATE(1710)*D*Y(42)*Y(333)+RATE(1711)*D*Y(42)&
    &*Y(333)+RATE(1721)*D*Y(53)*Y(333)+RATE(1722)*D*Y(53)*Y(333)+RATE(1723)*D&
    &*Y(53)*Y(333)+RATE(1724)*D*Y(53)*Y(333)+RATE(1725)*D*Y(53)*Y(333)&
    &+RATE(1743)*D*Y(187)*Y(333)+RATE(1769)*D*Y(83)*Y(333)+RATE(1793)*D*Y(100&
    &)*Y(333)+RATE(1797)*D*Y(234)*Y(333)+RATE(1944)*D*Y(4)*Y(333)+RATE(1951)&
    &*D*Y(6)*Y(333)+RATE(2027)*D*Y(8)*Y(333)+RATE(2053)*D*Y(196)*Y(333)&
    &+RATE(2054)*D*Y(196)*Y(333)+RATE(2063)*D*Y(132)*Y(333)+RATE(2064)*D&
    &*Y(132)*Y(333)+RATE(2065)*D*Y(132)*Y(333)+RATE(2066)*D*Y(132)*Y(333)&
    &+RATE(2071)*D*Y(132)*Y(333)+RATE(2072)*D*Y(262)*Y(333)+RATE(2073)*D&
    &*Y(262)*Y(333)+RATE(2074)*D*Y(262)*Y(333)+RATE(2075)*D*Y(160)*Y(333)&
    &+RATE(2076)*D*Y(160)*Y(333)+RATE(2119)*D*Y(63)*Y(333)+RATE(2120)*D*Y(63)&
    &*Y(333)+RATE(2121)*D*Y(63)*Y(333)+RATE(2141)*D*Y(184)*Y(333)+RATE(2142)&
    &*D*Y(184)*Y(333)+RATE(2144)*D*Y(184)*Y(333)+RATE(2145)*D*Y(329)*Y(333)&
    &+RATE(2146)*D*Y(329)*Y(333)+RATE(2147)*D*Y(10)*Y(333)+RATE(2148)*D*Y(10)&
    &*Y(333)+RATE(2207)*D*Y(142)*Y(333)+RATE(2208)*D*Y(142)*Y(333)+RATE(2209)&
    &*D*Y(142)*Y(333)+RATE(2210)*D*Y(142)*Y(333)+RATE(2211)*D*Y(142)*Y(333)&
    &+RATE(2214)*D*Y(274)*Y(333)+RATE(2215)*D*Y(274)*Y(333)+RATE(2216)*D*Y(66&
    &)*Y(333)+RATE(2217)*D*Y(66)*Y(333)+RATE(2218)*D*Y(66)*Y(333)+RATE(2219)&
    &*D*Y(66)*Y(333)+RATE(2235)*D*Y(189)*Y(333)+RATE(2236)*D*Y(189)*Y(333)&
    &+RATE(2237)*D*Y(189)*Y(333)+RATE(2238)*D*Y(189)*Y(333)+RATE(2240)*D&
    &*Y(194)*Y(333)+RATE(2258)*D*Y(91)*Y(333)+RATE(2267)*D*Y(102)*Y(333)&
    &+RATE(2268)*D*Y(102)*Y(333)+RATE(2269)*D*Y(102)*Y(333)+RATE(2294)*D&
    &*Y(117)*Y(333)+RATE(2315)*D*Y(243)*Y(333)+RATE(2316)*D*Y(243)*Y(333)&
    &+RATE(2317)*D*Y(243)*Y(333)+RATE(2318)*D*Y(245)*Y(333)+RATE(2319)*D&
    &*Y(245)*Y(333)+RATE(2445)*D*Y(13)*Y(333)+RATE(2446)*D*Y(15)*Y(333)&
    &+RATE(2456)*D*Y(145)*Y(333)+RATE(2459)*D*Y(275)*Y(333)+RATE(2460)*D&
    &*Y(119)*Y(333)+RATE(2461)*D*Y(306)*Y(333)+RATE(2462)*D*Y(306)*Y(333)&
    &+RATE(2464)*D*Y(174)*Y(333)+RATE(2467)*D*Y(325)*Y(333)+RATE(2468)*D&
    &*Y(325)*Y(333)+RATE(2469)*D*Y(307)*Y(333)+RATE(2470)*D*Y(307)*Y(333)&
    &+RATE(2471)*D*Y(281)*Y(333)+RATE(2472)*D*Y(326)*Y(333)+RATE(2473)*D&
    &*Y(326)*Y(333)+RATE(2474)*D*Y(326)*Y(333)+RATE(2488)*D*Y(70)*Y(333)&
    &+RATE(2598)*D*Y(28)*Y(333)+RATE(2608)*D*Y(104)*Y(333)+RATE(2614)*D*Y(120&
    &)*Y(333)+RATE(2615)*D*Y(120)*Y(333)+RATE(2667)*D*Y(36)*Y(333)+RATE(2733)&
    &*D*Y(44)*Y(333)+RATE(2734)*D*Y(44)*Y(333)+RATE(2797)*D*Y(55)*Y(333)&
    &+RATE(2798)*D*Y(55)*Y(333)+RATE(2805)*D*Y(64)*Y(333)+RATE(2806)*D*Y(64)&
    &*Y(333)+RATE(2807)*D*Y(64)*Y(333)+RATE(2822)*D*Y(134)*Y(333)+RATE(2823)&
    &*D*Y(266)*Y(333)+RATE(2953)*D*Y(48)*Y(333)+RATE(2967)*D*Y(162)*Y(333)&
    &+RATE(2972)*D*Y(176)*Y(333)+RATE(2974)*D*Y(301)*Y(333)+RATE(2975)*D&
    &*Y(301)*Y(333)+RATE(2976)*D*Y(301)*Y(333)+RATE(3029)*D*Y(57)*Y(333)&
    &+RATE(3067)*D*Y(165)*Y(333)+RATE(3068)*D*Y(318)*Y(333)+RATE(3083)*D&
    &*Y(106)*Y(333)+RATE(3084)*D*Y(208)*Y(333)+RATE(3085)*D*Y(293)*Y(333)&
    &+RATE(3086)*D*Y(293)*Y(333)+RATE(3087)*D*Y(313)*Y(333)+RATE(3088)*D&
    &*Y(313)*Y(333)+RATE(3091)*D*Y(122)*Y(333)+RATE(3092)*D*Y(136)*Y(333)&
    &+RATE(3093)*D*Y(136)*Y(333)+RATE(3094)*D*Y(136)*Y(333)+RATE(3097)*D&
    &*Y(147)*Y(333)+RATE(3098)*D*Y(147)*Y(333)+RATE(3099)*D*Y(167)*Y(333)&
    &+RATE(3100)*D*Y(167)*Y(333)+RATE(3101)*D*Y(177)*Y(333)+RATE(3102)*D&
    &*Y(177)*Y(333)+RATE(3103)*D*Y(236)*Y(333)+RATE(3104)*D*Y(247)*Y(333)&
    &+RATE(3105)*D*Y(247)*Y(333)+RATE(3106)*D*Y(303)*Y(333)+RATE(3107)*D&
    &*Y(278)*Y(333)+RATE(3115)*D*Y(319)*Y(333)+RATE(3116)*D*Y(319)*Y(333)
    PROD = RATE(99)*Y(16)+RATE(100)*Y(186)+RATE(101)*Y(99)+RATE(102)*Y(2&
    &)+RATE(103)*Y(6)+RATE(104)*Y(6)+RATE(106)*Y(11)+RATE(107)*Y(27)+RATE(108&
    &)*Y(46)+RATE(110)*Y(16)+RATE(113)*Y(75)+RATE(114)*Y(80)+RATE(125)*Y(24)&
    &+RATE(129)*Y(32)+RATE(137)*Y(186)+RATE(141)*Y(233)+RATE(143)*Y(2)&
    &+RATE(148)*Y(183)+RATE(156)*Y(116)+RATE(158)*Y(244)+RATE(159)*Y(11)&
    &+RATE(165)*Y(69)+RATE(166)*Y(27)+RATE(170)*Y(35)+RATE(171)*Y(43)&
    &+RATE(174)*Y(54)+RATE(176)*Y(133)+RATE(180)*Y(46)+RATE(181)*Y(161)&
    &+RATE(185)*Y(300)+RATE(188)*Y(163)+RATE(190)*Y(105)+RATE(201)*Y(277)&
    &+RATE(830)*Y(16)+RATE(831)*Y(67)+RATE(835)*Y(75)+RATE(837)*Y(80)&
    &+RATE(847)*Y(20)+RATE(849)*Y(24)+RATE(856)*Y(32)+RATE(864)*Y(159)&
    &+RATE(868)*Y(41)+RATE(872)*Y(186)+RATE(877)*Y(233)+RATE(885)*Y(131)&
    &+RATE(886)*Y(131)+RATE(888)*Y(62)+RATE(891)*Y(183)+RATE(901)*Y(193)&
    &+RATE(904)*Y(116)+RATE(907)*Y(244)+RATE(914)*Y(324)+RATE(916)*Y(69)&
    &+RATE(918)*Y(291)+RATE(921)*Y(35)+RATE(923)*Y(43)+RATE(926)*Y(54)&
    &+RATE(928)*Y(133)+RATE(932)*Y(161)+RATE(938)*Y(300)+RATE(941)*Y(56)&
    &+RATE(943)*Y(163)+RATE(944)*Y(316)+RATE(946)*Y(105)+RATE(953)*Y(135)&
    &+RATE(956)*Y(146)+RATE(962)*Y(235)+RATE(966)*Y(277)+RATE(1431)*D*Y(20)&
    &*Y(46)+RATE(1951)*D*Y(6)*Y(333)
    YDOT(333) = PROD-LOSS
    PROD = YDOT(5)+YDOT(9)+YDOT(14)+YDOT(19)+YDOT(23)+YDOT(30)+YDOT(31)&
    &+YDOT(38)+YDOT(39)+YDOT(49)+YDOT(50)+YDOT(51)+YDOT(59)+YDOT(60)+YDOT(65)&
    &+YDOT(73)+YDOT(74)+YDOT(78)+YDOT(85)+YDOT(86)+YDOT(94)+YDOT(95)+YDOT(96)&
    &+YDOT(111)+YDOT(112)+YDOT(113)+YDOT(114)+YDOT(115)+YDOT(126)+YDOT(127)&
    &+YDOT(128)+YDOT(138)+YDOT(139)+YDOT(140)+YDOT(152)+YDOT(153)+YDOT(154)&
    &+YDOT(155)+YDOT(168)+YDOT(169)+YDOT(170)+YDOT(171)+YDOT(180)+YDOT(181)&
    &+YDOT(185)+YDOT(190)+YDOT(195)+YDOT(202)+YDOT(203)+YDOT(210)+YDOT(211)&
    &+YDOT(215)+YDOT(223)+YDOT(224)+YDOT(225)+YDOT(228)+YDOT(239)+YDOT(240)&
    &+YDOT(241)+YDOT(242)+YDOT(253)+YDOT(254)+YDOT(255)+YDOT(256)+YDOT(269)&
    &+YDOT(270)+YDOT(271)+YDOT(272)+YDOT(273)+YDOT(279)+YDOT(283)+YDOT(286)&
    &+YDOT(289)+YDOT(296)+YDOT(297)+YDOT(304)+YDOT(305)+YDOT(311)+YDOT(321)&
    &+YDOT(322)+YDOT(323)+YDOT(328)+YDOT(332)
    YDOT(334) = PROD
    PROD = YDOT(3)+YDOT(7)+YDOT(12)+YDOT(17)+YDOT(22)+YDOT(26)+YDOT(29)&
    &+YDOT(34)+YDOT(37)+YDOT(40)+YDOT(45)+YDOT(47)+YDOT(52)+YDOT(58)+YDOT(61)&
    &+YDOT(71)+YDOT(72)+YDOT(77)+YDOT(79)+YDOT(84)+YDOT(87)+YDOT(88)+YDOT(93)&
    &+YDOT(97)+YDOT(98)+YDOT(107)+YDOT(108)+YDOT(110)+YDOT(118)+YDOT(124)&
    &+YDOT(125)+YDOT(129)+YDOT(130)+YDOT(137)+YDOT(143)+YDOT(148)+YDOT(150)&
    &+YDOT(151)+YDOT(156)+YDOT(157)+YDOT(158)+YDOT(164)+YDOT(178)+YDOT(179)&
    &+YDOT(182)+YDOT(188)+YDOT(191)+YDOT(198)+YDOT(201)+YDOT(205)+YDOT(209)&
    &+YDOT(212)+YDOT(216)+YDOT(217)+YDOT(222)+YDOT(226)+YDOT(229)+YDOT(230)&
    &+YDOT(231)+YDOT(238)+YDOT(246)+YDOT(248)+YDOT(250)+YDOT(252)+YDOT(257)&
    &+YDOT(258)+YDOT(259)+YDOT(260)+YDOT(268)+YDOT(276)+YDOT(282)+YDOT(285)&
    &+YDOT(287)+YDOT(294)+YDOT(295)+YDOT(298)+YDOT(299)+YDOT(310)+YDOT(312)&
    &+YDOT(315)+YDOT(317)+YDOT(327)+YDOT(331)
    YDOT(335) = PROD
!Update surface species for bulk growth, replace surfaceCoverage with alpha_des
!Since ydot(surface_index) is negative, bulk is lost and surface forms
IF (YDOT(335) .lt. 0) THEN
    surfaceCoverage = MIN(1.0,safeBulk/safeMantle)
    YDOT(3)=YDOT(3)-YDOT(335)*surfaceCoverage*Y(5)/safeBulk
    YDOT(5)=YDOT(5)+YDOT(335)*surfaceCoverage*Y(5)/safeBulk
    YDOT(7)=YDOT(7)-YDOT(335)*surfaceCoverage*Y(9)/safeBulk
    YDOT(9)=YDOT(9)+YDOT(335)*surfaceCoverage*Y(9)/safeBulk
    YDOT(12)=YDOT(12)-YDOT(335)*surfaceCoverage*Y(14)/safeBulk
    YDOT(14)=YDOT(14)+YDOT(335)*surfaceCoverage*Y(14)/safeBulk
    YDOT(17)=YDOT(17)-YDOT(335)*surfaceCoverage*Y(19)/safeBulk
    YDOT(19)=YDOT(19)+YDOT(335)*surfaceCoverage*Y(19)/safeBulk
    YDOT(22)=YDOT(22)-YDOT(335)*surfaceCoverage*Y(23)/safeBulk
    YDOT(23)=YDOT(23)+YDOT(335)*surfaceCoverage*Y(23)/safeBulk
    YDOT(26)=YDOT(26)-YDOT(335)*surfaceCoverage*Y(30)/safeBulk
    YDOT(29)=YDOT(29)-YDOT(335)*surfaceCoverage*Y(31)/safeBulk
    YDOT(30)=YDOT(30)+YDOT(335)*surfaceCoverage*Y(30)/safeBulk
    YDOT(31)=YDOT(31)+YDOT(335)*surfaceCoverage*Y(31)/safeBulk
    YDOT(34)=YDOT(34)-YDOT(335)*surfaceCoverage*Y(38)/safeBulk
    YDOT(37)=YDOT(37)-YDOT(335)*surfaceCoverage*Y(39)/safeBulk
    YDOT(38)=YDOT(38)+YDOT(335)*surfaceCoverage*Y(38)/safeBulk
    YDOT(39)=YDOT(39)+YDOT(335)*surfaceCoverage*Y(39)/safeBulk
    YDOT(40)=YDOT(40)-YDOT(335)*surfaceCoverage*Y(49)/safeBulk
    YDOT(45)=YDOT(45)-YDOT(335)*surfaceCoverage*Y(50)/safeBulk
    YDOT(47)=YDOT(47)-YDOT(335)*surfaceCoverage*Y(51)/safeBulk
    YDOT(49)=YDOT(49)+YDOT(335)*surfaceCoverage*Y(49)/safeBulk
    YDOT(50)=YDOT(50)+YDOT(335)*surfaceCoverage*Y(50)/safeBulk
    YDOT(51)=YDOT(51)+YDOT(335)*surfaceCoverage*Y(51)/safeBulk
    YDOT(52)=YDOT(52)-YDOT(335)*surfaceCoverage*Y(59)/safeBulk
    YDOT(58)=YDOT(58)-YDOT(335)*surfaceCoverage*Y(60)/safeBulk
    YDOT(59)=YDOT(59)+YDOT(335)*surfaceCoverage*Y(59)/safeBulk
    YDOT(60)=YDOT(60)+YDOT(335)*surfaceCoverage*Y(60)/safeBulk
    YDOT(61)=YDOT(61)-YDOT(335)*surfaceCoverage*Y(65)/safeBulk
    YDOT(65)=YDOT(65)+YDOT(335)*surfaceCoverage*Y(65)/safeBulk
    YDOT(71)=YDOT(71)-YDOT(335)*surfaceCoverage*Y(73)/safeBulk
    YDOT(72)=YDOT(72)-YDOT(335)*surfaceCoverage*Y(74)/safeBulk
    YDOT(73)=YDOT(73)+YDOT(335)*surfaceCoverage*Y(73)/safeBulk
    YDOT(74)=YDOT(74)+YDOT(335)*surfaceCoverage*Y(74)/safeBulk
    YDOT(77)=YDOT(77)-YDOT(335)*surfaceCoverage*Y(78)/safeBulk
    YDOT(78)=YDOT(78)+YDOT(335)*surfaceCoverage*Y(78)/safeBulk
    YDOT(79)=YDOT(79)-YDOT(335)*surfaceCoverage*Y(85)/safeBulk
    YDOT(84)=YDOT(84)-YDOT(335)*surfaceCoverage*Y(86)/safeBulk
    YDOT(85)=YDOT(85)+YDOT(335)*surfaceCoverage*Y(85)/safeBulk
    YDOT(86)=YDOT(86)+YDOT(335)*surfaceCoverage*Y(86)/safeBulk
    YDOT(87)=YDOT(87)-YDOT(335)*surfaceCoverage*Y(94)/safeBulk
    YDOT(88)=YDOT(88)-YDOT(335)*surfaceCoverage*Y(95)/safeBulk
    YDOT(93)=YDOT(93)-YDOT(335)*surfaceCoverage*Y(96)/safeBulk
    YDOT(94)=YDOT(94)+YDOT(335)*surfaceCoverage*Y(94)/safeBulk
    YDOT(95)=YDOT(95)+YDOT(335)*surfaceCoverage*Y(95)/safeBulk
    YDOT(96)=YDOT(96)+YDOT(335)*surfaceCoverage*Y(96)/safeBulk
    YDOT(97)=YDOT(97)-YDOT(335)*surfaceCoverage*Y(111)/safeBulk
    YDOT(98)=YDOT(98)-YDOT(335)*surfaceCoverage*Y(112)/safeBulk
    YDOT(107)=YDOT(107)-YDOT(335)*surfaceCoverage*Y(113)/safeBulk
    YDOT(108)=YDOT(108)-YDOT(335)*surfaceCoverage*Y(114)/safeBulk
    YDOT(110)=YDOT(110)-YDOT(335)*surfaceCoverage*Y(115)/safeBulk
    YDOT(111)=YDOT(111)+YDOT(335)*surfaceCoverage*Y(111)/safeBulk
    YDOT(112)=YDOT(112)+YDOT(335)*surfaceCoverage*Y(112)/safeBulk
    YDOT(113)=YDOT(113)+YDOT(335)*surfaceCoverage*Y(113)/safeBulk
    YDOT(114)=YDOT(114)+YDOT(335)*surfaceCoverage*Y(114)/safeBulk
    YDOT(115)=YDOT(115)+YDOT(335)*surfaceCoverage*Y(115)/safeBulk
    YDOT(118)=YDOT(118)-YDOT(335)*surfaceCoverage*Y(126)/safeBulk
    YDOT(124)=YDOT(124)-YDOT(335)*surfaceCoverage*Y(127)/safeBulk
    YDOT(125)=YDOT(125)-YDOT(335)*surfaceCoverage*Y(128)/safeBulk
    YDOT(126)=YDOT(126)+YDOT(335)*surfaceCoverage*Y(126)/safeBulk
    YDOT(127)=YDOT(127)+YDOT(335)*surfaceCoverage*Y(127)/safeBulk
    YDOT(128)=YDOT(128)+YDOT(335)*surfaceCoverage*Y(128)/safeBulk
    YDOT(129)=YDOT(129)-YDOT(335)*surfaceCoverage*Y(138)/safeBulk
    YDOT(130)=YDOT(130)-YDOT(335)*surfaceCoverage*Y(139)/safeBulk
    YDOT(137)=YDOT(137)-YDOT(335)*surfaceCoverage*Y(140)/safeBulk
    YDOT(138)=YDOT(138)+YDOT(335)*surfaceCoverage*Y(138)/safeBulk
    YDOT(139)=YDOT(139)+YDOT(335)*surfaceCoverage*Y(139)/safeBulk
    YDOT(140)=YDOT(140)+YDOT(335)*surfaceCoverage*Y(140)/safeBulk
    YDOT(143)=YDOT(143)-YDOT(335)*surfaceCoverage*Y(152)/safeBulk
    YDOT(148)=YDOT(148)-YDOT(335)*surfaceCoverage*Y(153)/safeBulk
    YDOT(150)=YDOT(150)-YDOT(335)*surfaceCoverage*Y(154)/safeBulk
    YDOT(151)=YDOT(151)-YDOT(335)*surfaceCoverage*Y(155)/safeBulk
    YDOT(152)=YDOT(152)+YDOT(335)*surfaceCoverage*Y(152)/safeBulk
    YDOT(153)=YDOT(153)+YDOT(335)*surfaceCoverage*Y(153)/safeBulk
    YDOT(154)=YDOT(154)+YDOT(335)*surfaceCoverage*Y(154)/safeBulk
    YDOT(155)=YDOT(155)+YDOT(335)*surfaceCoverage*Y(155)/safeBulk
    YDOT(156)=YDOT(156)-YDOT(335)*surfaceCoverage*Y(168)/safeBulk
    YDOT(157)=YDOT(157)-YDOT(335)*surfaceCoverage*Y(169)/safeBulk
    YDOT(158)=YDOT(158)-YDOT(335)*surfaceCoverage*Y(170)/safeBulk
    YDOT(164)=YDOT(164)-YDOT(335)*surfaceCoverage*Y(171)/safeBulk
    YDOT(168)=YDOT(168)+YDOT(335)*surfaceCoverage*Y(168)/safeBulk
    YDOT(169)=YDOT(169)+YDOT(335)*surfaceCoverage*Y(169)/safeBulk
    YDOT(170)=YDOT(170)+YDOT(335)*surfaceCoverage*Y(170)/safeBulk
    YDOT(171)=YDOT(171)+YDOT(335)*surfaceCoverage*Y(171)/safeBulk
    YDOT(178)=YDOT(178)-YDOT(335)*surfaceCoverage*Y(180)/safeBulk
    YDOT(179)=YDOT(179)-YDOT(335)*surfaceCoverage*Y(181)/safeBulk
    YDOT(180)=YDOT(180)+YDOT(335)*surfaceCoverage*Y(180)/safeBulk
    YDOT(181)=YDOT(181)+YDOT(335)*surfaceCoverage*Y(181)/safeBulk
    YDOT(182)=YDOT(182)-YDOT(335)*surfaceCoverage*Y(185)/safeBulk
    YDOT(185)=YDOT(185)+YDOT(335)*surfaceCoverage*Y(185)/safeBulk
    YDOT(188)=YDOT(188)-YDOT(335)*surfaceCoverage*Y(190)/safeBulk
    YDOT(190)=YDOT(190)+YDOT(335)*surfaceCoverage*Y(190)/safeBulk
    YDOT(191)=YDOT(191)-YDOT(335)*surfaceCoverage*Y(195)/safeBulk
    YDOT(195)=YDOT(195)+YDOT(335)*surfaceCoverage*Y(195)/safeBulk
    YDOT(198)=YDOT(198)-YDOT(335)*surfaceCoverage*Y(202)/safeBulk
    YDOT(201)=YDOT(201)-YDOT(335)*surfaceCoverage*Y(203)/safeBulk
    YDOT(202)=YDOT(202)+YDOT(335)*surfaceCoverage*Y(202)/safeBulk
    YDOT(203)=YDOT(203)+YDOT(335)*surfaceCoverage*Y(203)/safeBulk
    YDOT(205)=YDOT(205)-YDOT(335)*surfaceCoverage*Y(210)/safeBulk
    YDOT(209)=YDOT(209)-YDOT(335)*surfaceCoverage*Y(211)/safeBulk
    YDOT(210)=YDOT(210)+YDOT(335)*surfaceCoverage*Y(210)/safeBulk
    YDOT(211)=YDOT(211)+YDOT(335)*surfaceCoverage*Y(211)/safeBulk
    YDOT(212)=YDOT(212)-YDOT(335)*surfaceCoverage*Y(215)/safeBulk
    YDOT(215)=YDOT(215)+YDOT(335)*surfaceCoverage*Y(215)/safeBulk
    YDOT(216)=YDOT(216)-YDOT(335)*surfaceCoverage*Y(223)/safeBulk
    YDOT(217)=YDOT(217)-YDOT(335)*surfaceCoverage*Y(224)/safeBulk
    YDOT(222)=YDOT(222)-YDOT(335)*surfaceCoverage*Y(225)/safeBulk
    YDOT(223)=YDOT(223)+YDOT(335)*surfaceCoverage*Y(223)/safeBulk
    YDOT(224)=YDOT(224)+YDOT(335)*surfaceCoverage*Y(224)/safeBulk
    YDOT(225)=YDOT(225)+YDOT(335)*surfaceCoverage*Y(225)/safeBulk
    YDOT(226)=YDOT(226)-YDOT(335)*surfaceCoverage*Y(228)/safeBulk
    YDOT(228)=YDOT(228)+YDOT(335)*surfaceCoverage*Y(228)/safeBulk
    YDOT(229)=YDOT(229)-YDOT(335)*surfaceCoverage*Y(239)/safeBulk
    YDOT(230)=YDOT(230)-YDOT(335)*surfaceCoverage*Y(240)/safeBulk
    YDOT(231)=YDOT(231)-YDOT(335)*surfaceCoverage*Y(241)/safeBulk
    YDOT(238)=YDOT(238)-YDOT(335)*surfaceCoverage*Y(242)/safeBulk
    YDOT(239)=YDOT(239)+YDOT(335)*surfaceCoverage*Y(239)/safeBulk
    YDOT(240)=YDOT(240)+YDOT(335)*surfaceCoverage*Y(240)/safeBulk
    YDOT(241)=YDOT(241)+YDOT(335)*surfaceCoverage*Y(241)/safeBulk
    YDOT(242)=YDOT(242)+YDOT(335)*surfaceCoverage*Y(242)/safeBulk
    YDOT(246)=YDOT(246)-YDOT(335)*surfaceCoverage*Y(253)/safeBulk
    YDOT(248)=YDOT(248)-YDOT(335)*surfaceCoverage*Y(254)/safeBulk
    YDOT(250)=YDOT(250)-YDOT(335)*surfaceCoverage*Y(255)/safeBulk
    YDOT(252)=YDOT(252)-YDOT(335)*surfaceCoverage*Y(256)/safeBulk
    YDOT(253)=YDOT(253)+YDOT(335)*surfaceCoverage*Y(253)/safeBulk
    YDOT(254)=YDOT(254)+YDOT(335)*surfaceCoverage*Y(254)/safeBulk
    YDOT(255)=YDOT(255)+YDOT(335)*surfaceCoverage*Y(255)/safeBulk
    YDOT(256)=YDOT(256)+YDOT(335)*surfaceCoverage*Y(256)/safeBulk
    YDOT(257)=YDOT(257)-YDOT(335)*surfaceCoverage*Y(269)/safeBulk
    YDOT(258)=YDOT(258)-YDOT(335)*surfaceCoverage*Y(270)/safeBulk
    YDOT(259)=YDOT(259)-YDOT(335)*surfaceCoverage*Y(271)/safeBulk
    YDOT(260)=YDOT(260)-YDOT(335)*surfaceCoverage*Y(272)/safeBulk
    YDOT(268)=YDOT(268)-YDOT(335)*surfaceCoverage*Y(273)/safeBulk
    YDOT(269)=YDOT(269)+YDOT(335)*surfaceCoverage*Y(269)/safeBulk
    YDOT(270)=YDOT(270)+YDOT(335)*surfaceCoverage*Y(270)/safeBulk
    YDOT(271)=YDOT(271)+YDOT(335)*surfaceCoverage*Y(271)/safeBulk
    YDOT(272)=YDOT(272)+YDOT(335)*surfaceCoverage*Y(272)/safeBulk
    YDOT(273)=YDOT(273)+YDOT(335)*surfaceCoverage*Y(273)/safeBulk
    YDOT(276)=YDOT(276)-YDOT(335)*surfaceCoverage*Y(279)/safeBulk
    YDOT(279)=YDOT(279)+YDOT(335)*surfaceCoverage*Y(279)/safeBulk
    YDOT(282)=YDOT(282)-YDOT(335)*surfaceCoverage*Y(283)/safeBulk
    YDOT(283)=YDOT(283)+YDOT(335)*surfaceCoverage*Y(283)/safeBulk
    YDOT(285)=YDOT(285)-YDOT(335)*surfaceCoverage*Y(286)/safeBulk
    YDOT(286)=YDOT(286)+YDOT(335)*surfaceCoverage*Y(286)/safeBulk
    YDOT(287)=YDOT(287)-YDOT(335)*surfaceCoverage*Y(289)/safeBulk
    YDOT(289)=YDOT(289)+YDOT(335)*surfaceCoverage*Y(289)/safeBulk
    YDOT(294)=YDOT(294)-YDOT(335)*surfaceCoverage*Y(296)/safeBulk
    YDOT(295)=YDOT(295)-YDOT(335)*surfaceCoverage*Y(297)/safeBulk
    YDOT(296)=YDOT(296)+YDOT(335)*surfaceCoverage*Y(296)/safeBulk
    YDOT(297)=YDOT(297)+YDOT(335)*surfaceCoverage*Y(297)/safeBulk
    YDOT(298)=YDOT(298)-YDOT(335)*surfaceCoverage*Y(304)/safeBulk
    YDOT(299)=YDOT(299)-YDOT(335)*surfaceCoverage*Y(305)/safeBulk
    YDOT(304)=YDOT(304)+YDOT(335)*surfaceCoverage*Y(304)/safeBulk
    YDOT(305)=YDOT(305)+YDOT(335)*surfaceCoverage*Y(305)/safeBulk
    YDOT(310)=YDOT(310)-YDOT(335)*surfaceCoverage*Y(311)/safeBulk
    YDOT(311)=YDOT(311)+YDOT(335)*surfaceCoverage*Y(311)/safeBulk
    YDOT(312)=YDOT(312)-YDOT(335)*surfaceCoverage*Y(321)/safeBulk
    YDOT(315)=YDOT(315)-YDOT(335)*surfaceCoverage*Y(322)/safeBulk
    YDOT(317)=YDOT(317)-YDOT(335)*surfaceCoverage*Y(323)/safeBulk
    YDOT(321)=YDOT(321)+YDOT(335)*surfaceCoverage*Y(321)/safeBulk
    YDOT(322)=YDOT(322)+YDOT(335)*surfaceCoverage*Y(322)/safeBulk
    YDOT(323)=YDOT(323)+YDOT(335)*surfaceCoverage*Y(323)/safeBulk
    YDOT(327)=YDOT(327)-YDOT(335)*surfaceCoverage*Y(328)/safeBulk
    YDOT(328)=YDOT(328)+YDOT(335)*surfaceCoverage*Y(328)/safeBulk
    YDOT(331)=YDOT(331)-YDOT(335)*surfaceCoverage*Y(332)/safeBulk
    YDOT(332)=YDOT(332)+YDOT(335)*surfaceCoverage*Y(332)/safeBulk
ELSE
    YDOT(3)=YDOT(3)-YDOT(335)*surfaceCoverage*Y(3)
    YDOT(5)=YDOT(5)+YDOT(335)*surfaceCoverage*Y(3)
    YDOT(7)=YDOT(7)-YDOT(335)*surfaceCoverage*Y(7)
    YDOT(9)=YDOT(9)+YDOT(335)*surfaceCoverage*Y(7)
    YDOT(12)=YDOT(12)-YDOT(335)*surfaceCoverage*Y(12)
    YDOT(14)=YDOT(14)+YDOT(335)*surfaceCoverage*Y(12)
    YDOT(17)=YDOT(17)-YDOT(335)*surfaceCoverage*Y(17)
    YDOT(19)=YDOT(19)+YDOT(335)*surfaceCoverage*Y(17)
    YDOT(22)=YDOT(22)-YDOT(335)*surfaceCoverage*Y(22)
    YDOT(23)=YDOT(23)+YDOT(335)*surfaceCoverage*Y(22)
    YDOT(26)=YDOT(26)-YDOT(335)*surfaceCoverage*Y(26)
    YDOT(29)=YDOT(29)-YDOT(335)*surfaceCoverage*Y(29)
    YDOT(30)=YDOT(30)+YDOT(335)*surfaceCoverage*Y(26)
    YDOT(31)=YDOT(31)+YDOT(335)*surfaceCoverage*Y(29)
    YDOT(34)=YDOT(34)-YDOT(335)*surfaceCoverage*Y(34)
    YDOT(37)=YDOT(37)-YDOT(335)*surfaceCoverage*Y(37)
    YDOT(38)=YDOT(38)+YDOT(335)*surfaceCoverage*Y(34)
    YDOT(39)=YDOT(39)+YDOT(335)*surfaceCoverage*Y(37)
    YDOT(40)=YDOT(40)-YDOT(335)*surfaceCoverage*Y(40)
    YDOT(45)=YDOT(45)-YDOT(335)*surfaceCoverage*Y(45)
    YDOT(47)=YDOT(47)-YDOT(335)*surfaceCoverage*Y(47)
    YDOT(49)=YDOT(49)+YDOT(335)*surfaceCoverage*Y(40)
    YDOT(50)=YDOT(50)+YDOT(335)*surfaceCoverage*Y(45)
    YDOT(51)=YDOT(51)+YDOT(335)*surfaceCoverage*Y(47)
    YDOT(52)=YDOT(52)-YDOT(335)*surfaceCoverage*Y(52)
    YDOT(58)=YDOT(58)-YDOT(335)*surfaceCoverage*Y(58)
    YDOT(59)=YDOT(59)+YDOT(335)*surfaceCoverage*Y(52)
    YDOT(60)=YDOT(60)+YDOT(335)*surfaceCoverage*Y(58)
    YDOT(61)=YDOT(61)-YDOT(335)*surfaceCoverage*Y(61)
    YDOT(65)=YDOT(65)+YDOT(335)*surfaceCoverage*Y(61)
    YDOT(71)=YDOT(71)-YDOT(335)*surfaceCoverage*Y(71)
    YDOT(72)=YDOT(72)-YDOT(335)*surfaceCoverage*Y(72)
    YDOT(73)=YDOT(73)+YDOT(335)*surfaceCoverage*Y(71)
    YDOT(74)=YDOT(74)+YDOT(335)*surfaceCoverage*Y(72)
    YDOT(77)=YDOT(77)-YDOT(335)*surfaceCoverage*Y(77)
    YDOT(78)=YDOT(78)+YDOT(335)*surfaceCoverage*Y(77)
    YDOT(79)=YDOT(79)-YDOT(335)*surfaceCoverage*Y(79)
    YDOT(84)=YDOT(84)-YDOT(335)*surfaceCoverage*Y(84)
    YDOT(85)=YDOT(85)+YDOT(335)*surfaceCoverage*Y(79)
    YDOT(86)=YDOT(86)+YDOT(335)*surfaceCoverage*Y(84)
    YDOT(87)=YDOT(87)-YDOT(335)*surfaceCoverage*Y(87)
    YDOT(88)=YDOT(88)-YDOT(335)*surfaceCoverage*Y(88)
    YDOT(93)=YDOT(93)-YDOT(335)*surfaceCoverage*Y(93)
    YDOT(94)=YDOT(94)+YDOT(335)*surfaceCoverage*Y(87)
    YDOT(95)=YDOT(95)+YDOT(335)*surfaceCoverage*Y(88)
    YDOT(96)=YDOT(96)+YDOT(335)*surfaceCoverage*Y(93)
    YDOT(97)=YDOT(97)-YDOT(335)*surfaceCoverage*Y(97)
    YDOT(98)=YDOT(98)-YDOT(335)*surfaceCoverage*Y(98)
    YDOT(107)=YDOT(107)-YDOT(335)*surfaceCoverage*Y(107)
    YDOT(108)=YDOT(108)-YDOT(335)*surfaceCoverage*Y(108)
    YDOT(110)=YDOT(110)-YDOT(335)*surfaceCoverage*Y(110)
    YDOT(111)=YDOT(111)+YDOT(335)*surfaceCoverage*Y(97)
    YDOT(112)=YDOT(112)+YDOT(335)*surfaceCoverage*Y(98)
    YDOT(113)=YDOT(113)+YDOT(335)*surfaceCoverage*Y(107)
    YDOT(114)=YDOT(114)+YDOT(335)*surfaceCoverage*Y(108)
    YDOT(115)=YDOT(115)+YDOT(335)*surfaceCoverage*Y(110)
    YDOT(118)=YDOT(118)-YDOT(335)*surfaceCoverage*Y(118)
    YDOT(124)=YDOT(124)-YDOT(335)*surfaceCoverage*Y(124)
    YDOT(125)=YDOT(125)-YDOT(335)*surfaceCoverage*Y(125)
    YDOT(126)=YDOT(126)+YDOT(335)*surfaceCoverage*Y(118)
    YDOT(127)=YDOT(127)+YDOT(335)*surfaceCoverage*Y(124)
    YDOT(128)=YDOT(128)+YDOT(335)*surfaceCoverage*Y(125)
    YDOT(129)=YDOT(129)-YDOT(335)*surfaceCoverage*Y(129)
    YDOT(130)=YDOT(130)-YDOT(335)*surfaceCoverage*Y(130)
    YDOT(137)=YDOT(137)-YDOT(335)*surfaceCoverage*Y(137)
    YDOT(138)=YDOT(138)+YDOT(335)*surfaceCoverage*Y(129)
    YDOT(139)=YDOT(139)+YDOT(335)*surfaceCoverage*Y(130)
    YDOT(140)=YDOT(140)+YDOT(335)*surfaceCoverage*Y(137)
    YDOT(143)=YDOT(143)-YDOT(335)*surfaceCoverage*Y(143)
    YDOT(148)=YDOT(148)-YDOT(335)*surfaceCoverage*Y(148)
    YDOT(150)=YDOT(150)-YDOT(335)*surfaceCoverage*Y(150)
    YDOT(151)=YDOT(151)-YDOT(335)*surfaceCoverage*Y(151)
    YDOT(152)=YDOT(152)+YDOT(335)*surfaceCoverage*Y(143)
    YDOT(153)=YDOT(153)+YDOT(335)*surfaceCoverage*Y(148)
    YDOT(154)=YDOT(154)+YDOT(335)*surfaceCoverage*Y(150)
    YDOT(155)=YDOT(155)+YDOT(335)*surfaceCoverage*Y(151)
    YDOT(156)=YDOT(156)-YDOT(335)*surfaceCoverage*Y(156)
    YDOT(157)=YDOT(157)-YDOT(335)*surfaceCoverage*Y(157)
    YDOT(158)=YDOT(158)-YDOT(335)*surfaceCoverage*Y(158)
    YDOT(164)=YDOT(164)-YDOT(335)*surfaceCoverage*Y(164)
    YDOT(168)=YDOT(168)+YDOT(335)*surfaceCoverage*Y(156)
    YDOT(169)=YDOT(169)+YDOT(335)*surfaceCoverage*Y(157)
    YDOT(170)=YDOT(170)+YDOT(335)*surfaceCoverage*Y(158)
    YDOT(171)=YDOT(171)+YDOT(335)*surfaceCoverage*Y(164)
    YDOT(178)=YDOT(178)-YDOT(335)*surfaceCoverage*Y(178)
    YDOT(179)=YDOT(179)-YDOT(335)*surfaceCoverage*Y(179)
    YDOT(180)=YDOT(180)+YDOT(335)*surfaceCoverage*Y(178)
    YDOT(181)=YDOT(181)+YDOT(335)*surfaceCoverage*Y(179)
    YDOT(182)=YDOT(182)-YDOT(335)*surfaceCoverage*Y(182)
    YDOT(185)=YDOT(185)+YDOT(335)*surfaceCoverage*Y(182)
    YDOT(188)=YDOT(188)-YDOT(335)*surfaceCoverage*Y(188)
    YDOT(190)=YDOT(190)+YDOT(335)*surfaceCoverage*Y(188)
    YDOT(191)=YDOT(191)-YDOT(335)*surfaceCoverage*Y(191)
    YDOT(195)=YDOT(195)+YDOT(335)*surfaceCoverage*Y(191)
    YDOT(198)=YDOT(198)-YDOT(335)*surfaceCoverage*Y(198)
    YDOT(201)=YDOT(201)-YDOT(335)*surfaceCoverage*Y(201)
    YDOT(202)=YDOT(202)+YDOT(335)*surfaceCoverage*Y(198)
    YDOT(203)=YDOT(203)+YDOT(335)*surfaceCoverage*Y(201)
    YDOT(205)=YDOT(205)-YDOT(335)*surfaceCoverage*Y(205)
    YDOT(209)=YDOT(209)-YDOT(335)*surfaceCoverage*Y(209)
    YDOT(210)=YDOT(210)+YDOT(335)*surfaceCoverage*Y(205)
    YDOT(211)=YDOT(211)+YDOT(335)*surfaceCoverage*Y(209)
    YDOT(212)=YDOT(212)-YDOT(335)*surfaceCoverage*Y(212)
    YDOT(215)=YDOT(215)+YDOT(335)*surfaceCoverage*Y(212)
    YDOT(216)=YDOT(216)-YDOT(335)*surfaceCoverage*Y(216)
    YDOT(217)=YDOT(217)-YDOT(335)*surfaceCoverage*Y(217)
    YDOT(222)=YDOT(222)-YDOT(335)*surfaceCoverage*Y(222)
    YDOT(223)=YDOT(223)+YDOT(335)*surfaceCoverage*Y(216)
    YDOT(224)=YDOT(224)+YDOT(335)*surfaceCoverage*Y(217)
    YDOT(225)=YDOT(225)+YDOT(335)*surfaceCoverage*Y(222)
    YDOT(226)=YDOT(226)-YDOT(335)*surfaceCoverage*Y(226)
    YDOT(228)=YDOT(228)+YDOT(335)*surfaceCoverage*Y(226)
    YDOT(229)=YDOT(229)-YDOT(335)*surfaceCoverage*Y(229)
    YDOT(230)=YDOT(230)-YDOT(335)*surfaceCoverage*Y(230)
    YDOT(231)=YDOT(231)-YDOT(335)*surfaceCoverage*Y(231)
    YDOT(238)=YDOT(238)-YDOT(335)*surfaceCoverage*Y(238)
    YDOT(239)=YDOT(239)+YDOT(335)*surfaceCoverage*Y(229)
    YDOT(240)=YDOT(240)+YDOT(335)*surfaceCoverage*Y(230)
    YDOT(241)=YDOT(241)+YDOT(335)*surfaceCoverage*Y(231)
    YDOT(242)=YDOT(242)+YDOT(335)*surfaceCoverage*Y(238)
    YDOT(246)=YDOT(246)-YDOT(335)*surfaceCoverage*Y(246)
    YDOT(248)=YDOT(248)-YDOT(335)*surfaceCoverage*Y(248)
    YDOT(250)=YDOT(250)-YDOT(335)*surfaceCoverage*Y(250)
    YDOT(252)=YDOT(252)-YDOT(335)*surfaceCoverage*Y(252)
    YDOT(253)=YDOT(253)+YDOT(335)*surfaceCoverage*Y(246)
    YDOT(254)=YDOT(254)+YDOT(335)*surfaceCoverage*Y(248)
    YDOT(255)=YDOT(255)+YDOT(335)*surfaceCoverage*Y(250)
    YDOT(256)=YDOT(256)+YDOT(335)*surfaceCoverage*Y(252)
    YDOT(257)=YDOT(257)-YDOT(335)*surfaceCoverage*Y(257)
    YDOT(258)=YDOT(258)-YDOT(335)*surfaceCoverage*Y(258)
    YDOT(259)=YDOT(259)-YDOT(335)*surfaceCoverage*Y(259)
    YDOT(260)=YDOT(260)-YDOT(335)*surfaceCoverage*Y(260)
    YDOT(268)=YDOT(268)-YDOT(335)*surfaceCoverage*Y(268)
    YDOT(269)=YDOT(269)+YDOT(335)*surfaceCoverage*Y(257)
    YDOT(270)=YDOT(270)+YDOT(335)*surfaceCoverage*Y(258)
    YDOT(271)=YDOT(271)+YDOT(335)*surfaceCoverage*Y(259)
    YDOT(272)=YDOT(272)+YDOT(335)*surfaceCoverage*Y(260)
    YDOT(273)=YDOT(273)+YDOT(335)*surfaceCoverage*Y(268)
    YDOT(276)=YDOT(276)-YDOT(335)*surfaceCoverage*Y(276)
    YDOT(279)=YDOT(279)+YDOT(335)*surfaceCoverage*Y(276)
    YDOT(282)=YDOT(282)-YDOT(335)*surfaceCoverage*Y(282)
    YDOT(283)=YDOT(283)+YDOT(335)*surfaceCoverage*Y(282)
    YDOT(285)=YDOT(285)-YDOT(335)*surfaceCoverage*Y(285)
    YDOT(286)=YDOT(286)+YDOT(335)*surfaceCoverage*Y(285)
    YDOT(287)=YDOT(287)-YDOT(335)*surfaceCoverage*Y(287)
    YDOT(289)=YDOT(289)+YDOT(335)*surfaceCoverage*Y(287)
    YDOT(294)=YDOT(294)-YDOT(335)*surfaceCoverage*Y(294)
    YDOT(295)=YDOT(295)-YDOT(335)*surfaceCoverage*Y(295)
    YDOT(296)=YDOT(296)+YDOT(335)*surfaceCoverage*Y(294)
    YDOT(297)=YDOT(297)+YDOT(335)*surfaceCoverage*Y(295)
    YDOT(298)=YDOT(298)-YDOT(335)*surfaceCoverage*Y(298)
    YDOT(299)=YDOT(299)-YDOT(335)*surfaceCoverage*Y(299)
    YDOT(304)=YDOT(304)+YDOT(335)*surfaceCoverage*Y(298)
    YDOT(305)=YDOT(305)+YDOT(335)*surfaceCoverage*Y(299)
    YDOT(310)=YDOT(310)-YDOT(335)*surfaceCoverage*Y(310)
    YDOT(311)=YDOT(311)+YDOT(335)*surfaceCoverage*Y(310)
    YDOT(312)=YDOT(312)-YDOT(335)*surfaceCoverage*Y(312)
    YDOT(315)=YDOT(315)-YDOT(335)*surfaceCoverage*Y(315)
    YDOT(317)=YDOT(317)-YDOT(335)*surfaceCoverage*Y(317)
    YDOT(321)=YDOT(321)+YDOT(335)*surfaceCoverage*Y(312)
    YDOT(322)=YDOT(322)+YDOT(335)*surfaceCoverage*Y(315)
    YDOT(323)=YDOT(323)+YDOT(335)*surfaceCoverage*Y(317)
    YDOT(327)=YDOT(327)-YDOT(335)*surfaceCoverage*Y(327)
    YDOT(328)=YDOT(328)+YDOT(335)*surfaceCoverage*Y(327)
    YDOT(331)=YDOT(331)-YDOT(335)*surfaceCoverage*Y(331)
    YDOT(332)=YDOT(332)+YDOT(335)*surfaceCoverage*Y(331)
ENDIF
!Update total rate of change of bulk and surface for bulk growth
    PROD = YDOT(5)+YDOT(9)+YDOT(14)+YDOT(19)+YDOT(23)+YDOT(30)+YDOT(31)&
    &+YDOT(38)+YDOT(39)+YDOT(49)+YDOT(50)+YDOT(51)+YDOT(59)+YDOT(60)+YDOT(65)&
    &+YDOT(73)+YDOT(74)+YDOT(78)+YDOT(85)+YDOT(86)+YDOT(94)+YDOT(95)+YDOT(96)&
    &+YDOT(111)+YDOT(112)+YDOT(113)+YDOT(114)+YDOT(115)+YDOT(126)+YDOT(127)&
    &+YDOT(128)+YDOT(138)+YDOT(139)+YDOT(140)+YDOT(152)+YDOT(153)+YDOT(154)&
    &+YDOT(155)+YDOT(168)+YDOT(169)+YDOT(170)+YDOT(171)+YDOT(180)+YDOT(181)&
    &+YDOT(185)+YDOT(190)+YDOT(195)+YDOT(202)+YDOT(203)+YDOT(210)+YDOT(211)&
    &+YDOT(215)+YDOT(223)+YDOT(224)+YDOT(225)+YDOT(228)+YDOT(239)+YDOT(240)&
    &+YDOT(241)+YDOT(242)+YDOT(253)+YDOT(254)+YDOT(255)+YDOT(256)+YDOT(269)&
    &+YDOT(270)+YDOT(271)+YDOT(272)+YDOT(273)+YDOT(279)+YDOT(283)+YDOT(286)&
    &+YDOT(289)+YDOT(296)+YDOT(297)+YDOT(304)+YDOT(305)+YDOT(311)+YDOT(321)&
    &+YDOT(322)+YDOT(323)+YDOT(328)+YDOT(332)
    YDOT(334) = PROD
    PROD = YDOT(3)+YDOT(7)+YDOT(12)+YDOT(17)+YDOT(22)+YDOT(26)+YDOT(29)&
    &+YDOT(34)+YDOT(37)+YDOT(40)+YDOT(45)+YDOT(47)+YDOT(52)+YDOT(58)+YDOT(61)&
    &+YDOT(71)+YDOT(72)+YDOT(77)+YDOT(79)+YDOT(84)+YDOT(87)+YDOT(88)+YDOT(93)&
    &+YDOT(97)+YDOT(98)+YDOT(107)+YDOT(108)+YDOT(110)+YDOT(118)+YDOT(124)&
    &+YDOT(125)+YDOT(129)+YDOT(130)+YDOT(137)+YDOT(143)+YDOT(148)+YDOT(150)&
    &+YDOT(151)+YDOT(156)+YDOT(157)+YDOT(158)+YDOT(164)+YDOT(178)+YDOT(179)&
    &+YDOT(182)+YDOT(188)+YDOT(191)+YDOT(198)+YDOT(201)+YDOT(205)+YDOT(209)&
    &+YDOT(212)+YDOT(216)+YDOT(217)+YDOT(222)+YDOT(226)+YDOT(229)+YDOT(230)&
    &+YDOT(231)+YDOT(238)+YDOT(246)+YDOT(248)+YDOT(250)+YDOT(252)+YDOT(257)&
    &+YDOT(258)+YDOT(259)+YDOT(260)+YDOT(268)+YDOT(276)+YDOT(282)+YDOT(285)&
    &+YDOT(287)+YDOT(294)+YDOT(295)+YDOT(298)+YDOT(299)+YDOT(310)+YDOT(312)&
    &+YDOT(315)+YDOT(317)+YDOT(327)+YDOT(331)
    YDOT(335) = PROD
