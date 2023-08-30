import numpy

from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles, Particle
from amuse.community.kepler.interface import Kepler

from amuse.ext.solarsystem import solar_system_in_time

#mass in units of 10^16kg
_lunar_data = [
    ["Earth", "Moon",         7.34767309e+6, 1.737, 384400., 0.0554, 318.15, 135.27, 5.16, 125.08],
    ["Mars", "Phobos",        1.08, 11.1, 9376., 0.0151, 150.057, 91.059, 1.075, 207.784],
    ["Mars", "Deimos",        0.2, 6.3, 23458., 0.0002, 260.729, 325.329, 1.788, 24.525],
    ["Jupiter", "Io",         8932000., 1821.5, 421800., 0.0041, 84.129, 342.021, 0.036, 43.977],	
    ["Jupiter", "Europa",     4800000., 1560.8,	671100., 0.0094, 88.970, 171.016, 0.466, 219.106],	
    ["Jupiter", "Ganymede",   14819000., 2631.2, 1070400., 0.0013, 192.417, 317.540, 0.177, 63.552],	
    ["Jupiter", "Callisto",   10759000., 2410.3, 1882700., 0.0074, 52.643, 181.408, 0.192, 298.848],	
    ["Jupiter", "Amalthea",       1, 50, 181400., 0.0032, 155.873, 185.194, 0.380, 108.946], 
    ["Jupiter", "Thebe",          1, 50, 221900., 0.0176, 234.269, 135.956, 1.080, 235.694],	 
    ["Jupiter", "Adrastea",       1, 50, 129000., 0.0018, 328.047, 135.673, 0.054, 228.378],	 
    ["Jupiter", "Metis",          1, 50, 128000., 0.0012, 297.177, 276.047, 0.019, 146.912], 
    ["Jupiter", "Himalia",        1, 50, 11460000., 0.1586, 331.234, 66.874, 28.612, 64.798], 
    ["Jupiter", "Elara",          1, 50, 11740000., 0.2108, 142.001, 330.985, 27.945, 115.511], 
    ["Jupiter", "Pasiphae",       1, 50, 23629000., 0.4062, 169.226, 279.769, 151.413, 314.223], 
    ["Jupiter", "Sinope",         1, 50, 23942000., 0.2552, 354.541, 165.352, 158.189, 309.199], 
    ["Jupiter", "Lysithea",       1, 50, 11717000., 0.1161, 49.670, 330.475, 27.663, 5.326], 
    ["Jupiter", "Carme",          1, 50, 23401000., 0.2546, 26.416, 233.375, 164.994, 114.854], 
    ["Jupiter", "Ananke",         1, 50, 21254000., 0.2332, 95.772, 253.384, 148.693, 15.959], 
    ["Jupiter", "Leda",           1, 50, 11164000., 0.1624, 269.393, 230.352, 27.882, 219.181], 
    ["Jupiter", "Callirrhoe",     1, 50, 24099000., 0.2796, 23.909, 107.962, 147.080, 283.104], 
    ["Jupiter", "Themisto",       1, 50, 7504000.,  0.2435, 217.147, 313.051, 42.977, 192.288], 
    ["Jupiter", "Megaclite",      1, 50, 23814000., 0.4156, 288.882, 135.272, 152.781, 280.575], 
    ["Jupiter", "Taygete",        1, 50, 23363000., 0.2523, 231.540, 94.756, 165.253, 305.114], 
    ["Jupiter", "Chaldene",       1, 50, 23181000., 0.2503, 243.878, 267.454, 165.155, 134.240], 
    ["Jupiter", "Harpalyke",      1, 50, 21106000., 0.2296, 134.505, 215.956, 148.759, 29.834], 
    ["Jupiter", "Kalyke",         1, 50, 23565000., 0.2466, 218.934, 255.702, 165.121, 43.864], 
    ["Jupiter", "Iocaste",        1, 50, 21272000., 0.2152, 64.727, 213.675, 149.411, 269.613], 
    ["Jupiter", "Erinome",        1, 50, 23286000., 0.2655, 10.274, 267.136, 164.914, 317.497], 
    ["Jupiter", "Isonoe",         1, 50, 23231000., 0.2471, 116.879, 124.941, 165.250, 130.961], 
    ["Jupiter", "Praxidike",      1, 50, 21148000., 0.2274, 190.862, 117.480, 148.885, 280.956], 
    ["Jupiter", "Autonoe",        1, 50, 24037000., 0.3152, 54.793, 142.035, 152.364, 272.817], 
    ["Jupiter", "Thyone",         1, 50, 21197000., 0.2307, 97.023, 238.786, 148.595, 233.022], 
    ["Jupiter", "Hermippe",       1, 50, 21297000., 0.2095, 300.836, 131.854, 150.740, 330.393], 
    ["Jupiter", "Aitne",          1, 50, 23317000., 0.2627, 99.401, 105.000, 165.048, 8.679], 
    ["Jupiter", "Eurydome",       1, 50, 23146000., 0.2755, 223.631, 287.689, 150.271, 302.470], 
    ["Jupiter", "Euanthe",        1, 50, 21039000., 0.2320, 320.635, 333.101, 148.915, 254.297], 
    ["Jupiter", "Euporie",        1, 50, 19336000., 0.1438, 89.904, 70.243, 145.740, 60.143], 
    ["Jupiter", "Orthosie",       1, 50, 21158000., 0.2807, 216.805, 204.517, 146.004, 221.949], 
    ["Jupiter", "Sponde",         1, 50, 23790000., 0.3112, 61.885, 174.044, 150.997, 116.363], 
    ["Jupiter", "Kale",           1, 50, 23306000., 0.2597, 44.233, 212.853, 164.944, 60.170], 
    ["Jupiter", "Pasithee",       1, 50, 23091000., 0.2682, 231.920, 215.443, 165.117, 327.729], 
    ["Jupiter", "Hegemone",       1, 50, 23575000., 0.3445, 197.144, 236.950, 154.164, 318.902], 
    ["Jupiter", "Mneme",          1, 50, 21033000., 0.2258, 40.542, 256.860, 148.585, 13.467], 
    ["Jupiter", "Aoede",          1, 50, 23974000., 0.4325, 59.739, 197.676, 158.272, 173.392], 
    ["Jupiter", "Thelxinoe",      1, 50, 21160000., 0.2201, 313.183, 268.013, 151.390, 169.962], 
    ["Jupiter", "Arche",          1, 50, 23352000., 0.2495, 171.632, 39.713, 165.015, 339.210], 
    ["Jupiter", "Kallichore",     1, 50, 23276000., 0.2509, 9.836, 55.937, 165.102, 30.339], 
    ["Jupiter", "Helike",         1, 50, 21065000., 0.1498, 299.482, 43.659, 154.842, 89.749], 
    ["Jupiter", "Carpo",          1, 50, 17056000., 0.4317, 90.372, 337.062, 51.624, 50.597], 
    ["Jupiter", "Eukelade",       1, 50, 23323000., 0.2619, 309.685, 204.846, 165.265, 193.558], 
    ["Jupiter", "Cyllene",        1, 50, 23800000., 0.4155, 187.429, 128.345, 150.336, 252.611], 
    ["Jupiter", "Kore",           1, 50, 24482000., 0.3313, 138.071, 33.416, 145.173, 313.355], 
    ["Jupiter", "Herse",          1, 50, 23408000., 0.2541, 330.295, 141.667, 164.964, 295.702], 
    ["Jupiter", "S/2000 J11",     1, 50, 12297000., 0.2320, 173.544, 309.734, 28.631, 294.497], 
    ["Jupiter", "S/2003 J2",      1, 50, 28347000., 0.4100, 165.201, 237.932, 157.291, 344.782], 
    ["Jupiter", "S/2003 J3",      1, 50, 20221000., 0.1969, 66.338, 311.780, 147.547, 231.489], 
    ["Jupiter", "S/2003 J4",      1, 50, 23929000., 0.3624, 197.401, 260.480, 149.589, 179.131], 
    ["Jupiter", "S/2003 J5",      1, 50, 23495000., 0.2476, 90.066, 336.636, 165.248, 176.683], 
    ["Jupiter", "S/2003 J9",      1, 50, 23385000., 0.2632, 292.662, 348.415, 165.047, 44.321], 
    ["Jupiter", "S/2003 J10",     1, 50, 23042000., 0.4299, 170.833, 258.937, 165.073, 151.911], 
    ["Jupiter", "S/2003 J12",     1, 50, 17830000., 0.4904, 13.288, 38.543, 151.003, 65.530], 
    ["Jupiter", "S/2003 J15",     1, 50, 22627000., 0.1899, 18.405, 58.865, 146.492, 236.674], 
    ["Jupiter", "S/2003 J16",     1, 50, 21097000., 0.2281, 57.681, 307.563, 148.683, 16.883], 
    ["Jupiter", "S/2003 J18",     1, 50, 20508000., 0.0895, 130.894, 202.160, 146.077, 158.247], 
    ["Jupiter", "S/2003 J19",     1, 50, 23533000., 0.2552, 176.668, 223.035, 165.116, 27.442], 
    ["Jupiter", "S/2003 J23",     1, 50, 23567000., 0.2746, 255.114, 144.222, 146.424, 41.706], 
    ["Jupiter", "S/2010 J1",      1, 50, 23449000., 0.2491, 189.230, 160.525, 165.100, 282.871], 
    ["Jupiter", "S/2010 J2",      1, 50, 21004000., 0.2267, 18.252, 312.074, 148.673, 5.802], 
    ["Jupiter", "S/2011 J1",      1, 50, 23446000., 0.2534, 31.514, 256.027, 165.318, 250.728], 
    ["Jupiter", "S/2011 J2",      1, 50, 23124000., 0.3493, 270.154, 285.597, 153.597, 24.866],
    ["Saturn", "Mimas", 3790, 208, 185539., 0.0196, 332.499, 14.848, 1.574, 173.027], 
    ["Saturn", "Enceladus", 10800, 257, 238042., 0.0000, 0.076, 199.686, 0.003, 342.507],        
    ["Saturn", "Tethys", 61800,	538, 294672., 0.0001, 45.202, 243.367, 1.091, 259.842], 
    ["Saturn", "Dione", 110000, 563, 377415., 0.0022, 284.315, 322.232, 0.028, 290.415], 
    ["Saturn", "Rhea", 231000, 765, 527068., 0.0002, 241.619, 179.781, 0.333, 351.042], 
    ["Saturn", "Titan", 13455000,  2575, 1221865., 0.0288, 180.532, 163.310, 0.306, 28.060], 
    ["Saturn", "Hyperion", 560, 180, 1500933., 0.0232, 303.178, 86.342, 0.615, 263.847], 
    ["Saturn", "Iapetus", 181000, 746, 3560854., 0.0293, 271.606, 201.789, 8.298, 81.105], 
    ["Saturn", "Phoebe",     1, 50, 12947918., 0.1634, 342.500, 53.038, 175.243, 241.086], 
    ["Saturn", "Janus", 190, 102, 151450., 0.0098, 16.012, 17.342, 0.165, 154.175], 
    ["Saturn", "Epimetheus", 53, 65, 151450., 0.0161, 88.975, 80.377, 0.353, 192.762], 
    ["Saturn", "Helene", 3, 99, 377444., 0.0000, 33.134, 43.186, 0.213, 163.112], 
    ["Saturn", "Telesto",     1, 50, 294720., 0.0002, 119.135, 260.157, 1.180, 229.182], 
    ["Saturn", "Calypso",     1, 50, 294721., 0.0005, 17.470, 156.660, 1.500, 314.226], 
    ["Saturn", "Atlas",     1, 50, 137774., 0.0011, 210.851, 283.282, 0.003, 236.422], 
    ["Saturn", "Prometheus", 16, 50, 139429., 0.0022, 37.514, 96.886, 0.007, 319.176], 
    ["Saturn", "Pandora", 15, 50, 141810., 0.0042, 66.248, 125.112, 0.050, 147.272], 
    ["Saturn", "Pan",     1, 50, 133585., 0.0000, 103.331, 351.187, 0.000, 52.076], 
    ["Saturn", "Methone",     1, 50, 194402., 0.0000, 134.636, 71.189, 0.013, 313.562], 
    ["Saturn", "Pallene",     1, 50, 212282., 0.0040, 16.074, 356.229, 0.001, 123.180], 
    ["Saturn", "Polydeuces",     1, 50, 377222., 0.0191, 311.847, 89.307, 0.175, 67.936],
    ["Saturn", "Daphnis",     1, 50, 136504., 0.0000, 266.931, 113.790, 0.003, 132.867], 
    ["Saturn", "Anthe",     1, 50, 196888., 0.0011, 138.902, 190.473, 0.015, 287.852], 
    ["Saturn", "Aegaeon",     1, 50, 167425., 0.0002, 152.905, 322.771, 0.001, 317.202], 
    ["Saturn", "Ymir",     1, 50, 23128000., 0.3338, 21.352, 228.673, 173.497, 192.937], 
    ["Saturn", "Paaliaq",     1, 50, 15204000., 0.3325, 237.522, 321.654, 46.228, 330.022], 
    ["Saturn", "Tarvos",     1, 50, 18243000., 0.5382, 274.104, 265.783, 33.725, 102.504], 
    ["Saturn", "Ijiraq",     1, 50, 11408000., 0.2717, 92.899, 17.328, 47.485, 130.779], 
    ["Saturn", "Suttungr",     1, 50, 19468000., 0.1139, 34.281, 321.133, 175.815, 227.259], 
    ["Saturn", "Kiviuq",     1, 50, 11384000., 0.3325, 90.205, 172.018, 46.764, 353.584], 
    ["Saturn", "Mundilfari",     1, 50, 18653000., 0.2097, 309.694, 92.821, 167.439, 82.856], 
    ["Saturn", "Albiorix",     1, 50, 16393000., 0.4797, 55.932, 32.828, 34.060, 102.512], 
    ["Saturn", "Skathi",     1, 50, 15635000., 0.2718, 203.517, 114.689, 152.633, 286.599], 
    ["Saturn", "Erriapus",     1, 50, 17602000., 0.4722, 282.522, 294.829, 34.481, 150.985], 
    ["Saturn", "Siarnaq",     1, 50, 18182000., 0.2802, 65.929, 201.288, 45.809, 47.826], 
    ["Saturn", "Thrymr",     1, 50, 20418000., 0.4659, 125.404, 30.075, 177.659, 285.762], 
    ["Saturn", "Narvi",     1, 50, 19349000., 0.4296, 169.959, 114.172, 145.731, 174.435], 
    ["Saturn", "Aegir",     1, 50, 20751000., 0.2524, 242.651, 26.017, 166.668, 179.064], 
    ["Saturn", "Bebhionn",     1, 50, 17116000., 0.4682, 358.141, 168.045, 35.101, 199.128], 
    ["Saturn", "Bergelmir",     1, 50, 19336000., 0.1420, 133.400, 306.494, 158.557, 202.164], 
    ["Saturn", "Bestla",     1, 50, 20209000., 0.5145, 81.185, 239.156, 145.136, 288.308], 
    ["Saturn", "Farbauti",     1, 50, 20390000., 0.2414, 342.995, 282.813, 156.520, 135.109], 
    ["Saturn", "Fenrir",     1, 50, 22454000., 0.1347, 120.982, 131.678, 164.963, 226.595], 
    ["Saturn", "Fornjot",     1, 50, 25146000., 0.2077, 324.787, 214.499, 170.372, 259.946], 
    ["Saturn", "Hati",     1, 50, 19868000., 0.3710, 21.286, 163.640, 165.808, 324.380], 
    ["Saturn", "Hyrrokkin",     1, 50, 18440000., 0.3359, 273.076, 291.841, 151.536, 45.402], 
    ["Saturn", "Kari",     1, 50, 22093000., 0.4756, 163.935, 286.021, 156.067, 281.211], 
    ["Saturn", "Loge",     1, 50, 23059000., 0.1862, 32.821, 337.237, 167.689, 343.811], 
    ["Saturn", "Skoll",     1, 50, 17668000., 0.4636, 193.115, 44.965, 161.010, 296.623], 
    ["Saturn", "Surtur",     1, 50, 22941000., 0.4459, 303.662, 136.191, 169.688, 236.537], 
    ["Saturn", "Jarnsaxa",     1, 50, 19354000., 0.2178, 237.422, 198.750, 163.649, 22.519],
    ["Saturn", "Greip",     1, 50, 18457000., 0.3146, 152.160, 314.541, 174.800, 349.350], 
    ["Saturn", "Tarqeq",     1, 50, 17962000., 0.1675, 34.767, 161.020, 46.291, 83.291], 
    ["Saturn", "S/2004 S7",     1, 50, 21000000., 0.5290, 84.036, 79.762, 165.693, 341.236], 
    ["Saturn", "S/2004 S12",     1, 50, 19886000., 0.3268, 87.128, 1.599, 165.261, 307.942], 
    ["Saturn", "S/2004 S13",     1, 50, 18406000., 0.2591, 346.186, 41.077, 168.798, 205.701], 
    ["Saturn", "S/2004 S17",     1, 50, 19448000., 0.1795, 180.792, 228.545, 168.239, 26.664], 
    ["Saturn", "S/2006 S1",     1, 50, 18780000., 0.1412, 154.950, 96.596, 156.180, 336.641], 
    ["Saturn", "S/2006 S3",     1, 50, 22428000., 0.3792, 188.728, 167.147, 158.631, 206.993], 
    ["Saturn", "S/2007 S2",     1, 50, 16718000., 0.1791, 57.720, 84.066, 174.057, 111.277], 
    ["Saturn", "S/2007 S3",     1, 50, 18938000., 0.1853, 111.854, 292.691, 177.595, 276.824],
    ["Uranus", "Ariel",     129000, 581.1, 190900., 0.0012, 115.349, 39.481, 0.041, 22.394], 
    ["Uranus", "Umbriel",     122000, 584.7, 266000., 0.0039, 84.709, 12.469, 0.128, 33.485], 
    ["Uranus", "Titania",     342000, 788.9, 436300., 0.0011, 284.400, 24.614, 0.079, 99.771], 
    ["Uranus", "Oberon",     288000, 761.4, 583500., 0.0014, 104.400, 283.088, 0.068, 279.771], 
    ["Uranus", "Miranda",     6600, 240, 129900., 0.0013, 68.312, 311.330, 4.338, 326.438], 
    ["Uranus", "Cordelia",     1, 50, 49800., 0.0003, 136.827, 254.805, 0.085, 38.374], 
    ["Uranus", "Ophelia",     1, 50, 53800., 0.0099, 17.761, 116.259, 0.104, 164.048], 
    ["Uranus", "Bianca",     1, 50, 59200., 0.0009, 8.293, 138.486, 0.193, 93.220], 
    ["Uranus", "Cressida",     1, 50, 61800., 0.0004, 44.236, 233.795, 0.006, 99.403], 
    ["Uranus", "Desdemona",     1, 50, 62700., 0.0001, 183.285, 184.627, 0.113, 306.089], 
    ["Uranus", "Juliet",     1, 50, 64400., 0.0007, 223.819, 244.696, 0.065, 200.155], 
    ["Uranus", "Portia",     1, 50, 66100., 0.0001, 222.433, 218.312, 0.059, 260.067], 
    ["Uranus", "Rosalind",     1, 50, 69900., 0.0001, 140.477, 136.181, 0.279, 12.847], 
    ["Uranus", "Belinda",     1, 50, 75300., 0.0001, 42.406, 357.224, 0.031, 279.337], 
    ["Uranus", "Puck",     1, 50, 86000., 0.0001, 177.094, 245.796, 0.319, 268.734], 
    ["Uranus", "Perdita",     1, 50, 76417., 0.0116, 253.925, 192.405, 0.470, 309.376], 
    ["Uranus", "Mab",     1, 50, 97736., 0.0025, 249.565, 273.769, 0.134, 350.737], 
    ["Uranus", "Cupid",     1, 50, 74392., 0.0013, 247.608, 163.830, 0.099, 182.793], 
    ["Uranus", "Caliban",     1, 50, 7231100., 0.1812, 354.339, 7.271, 141.529, 171.189], 
    ["Uranus", "Sycorax",     1, 50, 12179400., 0.5219, 20.103, 266.583, 159.420, 263.034], 
    ["Uranus", "Prospero",     1, 50, 16276800., 0.4445, 174.152, 233.586, 151.830, 319.003], 
    ["Uranus", "Setebos",     1, 50, 17420400., 0.5908, 359.953, 179.449, 158.235, 250.235], 
    ["Uranus", "Stephano",     1, 50, 8007400., 0.2248, 14.956, 270.163, 143.819, 191.411], 
    ["Uranus", "Trinculo",     1, 50, 8505200., 0.2194, 158.688, 180.374, 166.971, 193.755], 
    ["Uranus", "Francisco",     1, 50, 4282900., 0.1324, 140.644, 3.202, 147.250, 100.738], 
    ["Uranus", "Margaret",     1, 50, 14146700., 0.6772, 90.017, 322.187, 57.367, 7.067], 
    ["Uranus", "Ferdinand",     1, 50, 20430000., 0.3993, 156.298, 26.163, 169.793, 217.350],
    ["Neptune", "Triton", 2140000, 1353, 354759., 0.0000, 66.142, 352.257, 156.865, 177.608],
    ["Neptune", "Nereid",    3000, 170, 5513818., 0.7507, 281.117, 216.692, 7.090, 335.570], 
    ["Neptune", "Naiad",       20, 48, 48227., 0.0003, 2.045, 30.035, 4.691, 42.279], 
    ["Neptune", "Thalassa",    40, 54, 50074., 0.0002, 237.065, 262.923, 0.135, 145.980],
    ["Neptune", "Despina",    200, 90, 52526., 0.0002, 176.857, 230.812, 0.068, 77.060],
    ["Neptune", "Galatea",    400, 102, 61953., 0.0001, 343.011, 65.999, 0.034, 37.247], 
    ["Neptune", "Larissa",    500, 108, 73548., 0.0014, 249.891, 166.246, 0.205, 308.127], 
    ["Neptune", "Proteus",   5000, 220, 117646., 0.0005, 67.968, 250.938, 0.075, 315.131], 
    ["Neptune", "S/2004 N1",    1, 10, 105284., 0.0000, 0.000, 302.652, 0.000, 0.000], 
    ["Neptune", "Halimede",    20, 30, 16681000., 0.2909, 162.119, 105.258, 112.898, 217.288], 
    ["Neptune", "Psamathe",     2, 20, 46705000., 0.4617, 144.158, 190.027, 137.679, 298.074], 
    ["Neptune", "Sao",         10, 20, 22619000., 0.2827, 65.047, 168.139, 49.907, 60.354], 
    ["Neptune", "Laomedeia",   10, 20, 23613000., 0.4339, 140.107, 285.863, 34.049, 59.124], 
    ["Neptune", "Neso",        20, 30, 50258000., 0.4243, 86.441, 260.648, 131.265, 49.151]
] 


def new_kepler():
  converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  kepler = Kepler(converter)
  kepler.initialize_code()
  return kepler

def get_position(mass_sun, mass_planet, ecc, semi, mean_anomaly, incl, argument, longitude, delta_t=0.|units.day):
  """
  cartesian position and velocity from orbital elements,
  where the orbit is evolved from given mean_anomaly 
  by time delta_t
  argument -- argument of perihelion
  longitude -- longitude of ascending node
  """
  kepler = new_kepler()
  kepler.initialize_from_elements(mass=(mass_sun+mass_planet),
                                  semi=semi,
                                  ecc=ecc,
                                  mean_anomaly=mean_anomaly)
  kepler.transform_to_time(time=delta_t)
  r = kepler.get_separation_vector()
  v = kepler.get_velocity_vector()
  
  kepler.stop()
  
  a1 = ([numpy.cos(longitude), -numpy.sin(longitude), 0.0], [numpy.sin(longitude), numpy.cos(longitude), 0.0], [0.0, 0.0, 1.0])
  a2 = ([1.0, 0.0, 0.0], [0.0, numpy.cos(incl), -numpy.sin(incl)], [0.0, numpy.sin(incl), numpy.cos(incl)])
  a3 = ([numpy.cos(argument), -numpy.sin(argument), 0.0], [numpy.sin(argument), numpy.cos(argument), 0.0], [0.0, 0.0, 1.0])
  A = numpy.dot(numpy.dot(a1,a2),a3)
  print(A, r)
  
  # old version from P2.7
  #  r_vec = numpy.dot(A,numpy.reshape(r,3,1))
  #  v_vec = numpy.dot(A,numpy.reshape(v,3,1))
  r_vec = numpy.dot(A,numpy.reshape(r,3,'F'))
  v_vec = numpy.dot(A,numpy.reshape(v,3,'F'))
  
  # for relative vectors
  r[0] = r_vec[0]
  r[1] = r_vec[1]
  r[2] = r_vec[2]
  v[0] = v_vec[0]
  v[1] = v_vec[1]
  v[2] = v_vec[2]
  
  return r,v

def get_moons_for_planet(planet, delta_JD=0.|units.day):
  """
  The Earth's moon
  as for JD = 2457099.500000000 = A.D. 2015-Mar-18 00:00:00.0000 (CT)
  https://ssd.jpl.nasa.gov/?sat_elem
  """

  data = numpy.array([tuple(entry) for entry in _lunar_data],
        dtype=[('planet_name','S10'), ('name','S10'), 
        ('mass','<f8'), ('radius','<f8'), ('semimajor_axis','<f8'), 
        ('eccentricity','<f8'), ('argument_of_peri','<f8'),
        ('mean_anomaly','<f8'), ('inclination','<f8'), ('longitude_oan','<f8')])
  moon_data = data[data['planet_name']==planet.name.encode('UTF-8')]
  print("Planet=", planet.name, "moon=", moon_data["name"])
  moons = Particles()
  if len(moon_data["name"]):
      print(len(moon_data["name"]))
      for moon in moon_data:
          #print moon
          r, v = get_position(planet.mass,
                              moon["mass"] * 1.e+16 | units.kg,
                              moon["eccentricity"],
                              moon["semimajor_axis"]|units.km,
                              numpy.deg2rad(moon["mean_anomaly"]),
                              numpy.deg2rad(moon["inclination"]),
                              numpy.deg2rad(moon["longitude_oan"]),
                              numpy.deg2rad(moon["argument_of_peri"]),
                              delta_t=delta_JD)
          single_moon = Particle()
          single_moon.type = "moon"
          single_moon.name = moon["name"]
          single_moon.mass = moon["mass"] * 1.e+16 | units.kg
          single_moon.hostname = moon["planet_name"]
          single_moon.radius = moon["radius"] | units.km
          single_moon.position = r
          single_moon.position += planet.position
          single_moon.velocity = v
          single_moon.velocity += planet.velocity
          moons.add_particle(single_moon)
    
  return moons

def new_lunar_system_in_time(time_JD=2457099.5|units.day):
  """
  Initial conditions of Solar system --
  particle set with the sun + eight moons,
  at the center-of-mass reference frame.

  Defined attributes: 
  name, mass, radius, x, y, z, vx, vy, vz
  """
  time_0 = 2457099.5 | units.day
  delta_JD = time_JD-time_0
  solar_system = solar_system_in_time(time_JD)
  solar_system[0].type = "star"
  solar_system[0].name = "sun"
  solar_system[1:].type = "planet"
  for pi in solar_system:
      moons = get_moons_for_planet(pi, delta_JD=delta_JD)
      solar_system.add_particles(moons)
  solar_system.move_to_center()
  
  ### to compare with JPL, relative positions and velocities need to be corrected for the
  # Sun's vectors with respect to the barycenter
  #r_s = (3.123390770608490E-03, -4.370830943817017E-04, -1.443425433116342E-04) | units.AU
  #v_s = (3.421633816761503E-06,  5.767414405893875E-06, -8.878039607570240E-08) | (units.AU / units.day)
  #print sun
  #print moons.position.in_(units.AU) + r_s
  #print moons.velocity.in_(units.AU/units.day) + v_s
  
  return solar_system

def new_lunar_system(Julian_date=-1|units.day):
    return new_lunar_system_in_time(Julian_date)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-d", dest="Julian_date", unit=units.day,
                      type=float, default = 2438871.5|units.day,
                      help="julian date [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
  o, arguments  = new_option_parser().parse_args()
  lunar_system = new_lunar_system(o.Julian_date)
  print(lunar_system)


