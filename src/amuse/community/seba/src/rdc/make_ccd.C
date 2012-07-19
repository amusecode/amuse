//// mk_ccd:     create CCD image in FITS format from input snapshot
////
//// Options:    -A    number of electrons to ADU [3]
////             -a    arcsecond per pixel scale [0.3]
////             -B    filter of observation [3]
////                   (options are: U=0, B=1, V=2, R=3, I=4) 
////             -b    add background stars to ccd [false]
////             -c    add cosmics rays to ccd [false]
////             -d    distance to cluster in parsec [1000]
////             -F    output filename [ccd.fits]
////             -f    full with half maximum of the PSF [0.7]
////             -L    upper luminosity limit in Lsun [1]
////             -N    add number of bad columns [0]
////             -P    Kolmogorov point spread function beta [5]
////             -p    give projecttion [1]
////                       -3: project on xy plane viewed from -z
////                       -2: project on xz plane viewed from -y
////                       -1: project on yz plane viewed from -x
////                        1: project on yz plane viewed from  x
////                        2: project on xz plane viewed from  y
////                        3: project on xy plane viewed from  z
////             -R    reddening [0]
////             -S    random number seed [clock]
////             -s    add standard stars [false];
////             -x    x-offset from density center [0]
////             -y    y-offset from density center [0]
////             -v    verbose, prints stars and positions to cerr [false]
////
//// Latest version (SPZ:3.0) October 1999.
////
//// Version 3.0 and later requires fits library
//
//
//   version 1.0:  September 1999   Simon F. Portegies Zwart
//                                  spz@komodo.bu.edu
//
//++ version 1.0: basic ccd creation routines.
//++         2.0: output file written in FITS format
//++         3.0: added: overexposure bleading
//++                     dask current
//++                     read out noise
//++                     scyntilation
//++                     thermal luminosity
//

#include "single_star.h"
#include "sstar_to_dyn.h"
#include "dyn.h"

#include "fitsio.h"

enum wavelength {X=-1, U=0, B, F, R, J};

#define XCCD_MAX 1024
#define YCCD_MAX 1024

//#define ARCSEC_PER_PIXEL 0.0455  // HST WFPC
//#define ARCSEC_PER_PIXEL 0.0966  // HST Planetary camera
//#define ARCSEC_PER_PIXEL 0.076
//#define ARCSEC_PER_PIXEL 0.49    // 0.9m Cassegrain telescope
//#define ARCSEC_PER_PIXEL 0.15
//#define ARCSEC_PER_PIXEL 0.30      // [arcsec] Ken Janes
//#define FWHM             0.70      // [arcsec]

#define SATURATION_LIMIT 65536

#ifdef TOOLBOX

local real hst_psf(real r) {

    real a, b, c, d, e;
	
    a = 6.66183E-01;
    b = 3.51739E-01;
    c = -1.17196E-01;
    d = 1.50760E-02;
    e = -6.55790E-04;

    real f=1;
    if(r<=9)
    f = a + r*(b + r*(c + r*(d + r*e))) - 0.3; 

    return f;
} 

local real Kolmogorov_psf(real r, real beta, real fwhm) {

    real hw = 0.5 * fwhm;

    real first_part = 1./(cnsts.mathematics(pi)*pow(hw, 2));
    real teller = (pow(2, 1./beta)-1) * (beta - 1);
    real noemer = pow(1 + (pow(2, 1./beta)-1)*pow(r/hw, 2), beta);

    real I = first_part * teller/noemer;

    return I;
}

local void print_ccd(real **ccd, char filename[], real t_exposure) { 

    fitsfile *fptr;
    int status, ii, jj;
    long  fpixel, nelements, exposure;
    //    unsigned short array[512][512];  
    unsigned short array[XCCD_MAX][YCCD_MAX];  

//    char filename[] = "ccd.fit";
    int bitpix   =  USHORT_IMG; 
    long naxis    =   2;
    //    long naxes[2] = { 512, 512 }; 
    long naxes[2] = { XCCD_MAX, YCCD_MAX }; 

    remove(filename); 
    status = 0;


    if (fits_create_file(&fptr, filename, &status)) 
         printf("%s\n", status);       

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printf("%s\n", status);     

    float val;
    for (jj = 0; jj < naxes[1]; jj++) {   
      for (ii = 0; ii < naxes[0]; ii++) {
//          cin >> val;
            array[jj][ii] = (unsigned short)ccd[jj][ii];
        }
    }

    fpixel = 1;
    nelements = naxes[0] * naxes[1]; 

    if ( fits_write_img(fptr, TUSHORT, fpixel, nelements, array[0], &status) )
         printf("%s\n", status);

    exposure = long(t_exposure);
    if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Exposure time in seconds", &status) )
         printf("%s\n", status);  

    if ( fits_close_file(fptr, &status) )  
         printf("%s\n", status);        
}

real **mkarray(int rows, int cols) {
    real ** ret = new real *[rows];
    for (int i=0; i<cols; i++)
	ret[i] = new real[cols];

    for (int i=0; i<rows; i++)
	for (int j=0; j<cols; j++) {
	ret[i][j] = 0;
    }

    cerr << "Datafile read"<<endl;
    return ret;
}

local bool add_star_to_ccd(real** ccd, real luminosity, 
			   real xpos, real ypos, real beta_PSF,
			   real arcsec_per_pixel, real fwhm) {

  int nx = (int)floor(xpos/arcsec_per_pixel);
  int ny = (int)floor(ypos/arcsec_per_pixel);
  
//  PRC(nx);PRC(ny);PRC(luminosity);
  if (nx<-10 || nx>XCCD_MAX+10 || ny<-10 || ny>YCCD_MAX+10) {

    //cerr << ": not added"<<endl;
    return false;
  }

//  if(luminosity>0) {
    real magn = 4.76 - 2.5*log10(luminosity);
//    PRC(nx);PRC(ny);PRL(magn);
//  }

  real r_locus;
  real norm = Kolmogorov_psf(0, beta_PSF, fwhm);
  int n_psf = int(10*fwhm/arcsec_per_pixel);
  
  for(int ix = Starlab::max(0, nx-n_psf); ix<Starlab::min(XCCD_MAX, nx+n_psf); ix++) 
    for(int iy = Starlab::max(0, ny-n_psf); iy<Starlab::min(YCCD_MAX, ny+n_psf); iy++) {

      r_locus = sqrt(pow(xpos-ix*arcsec_per_pixel, 2) 
		     +      pow(ypos-iy*arcsec_per_pixel, 2)); 
      
      ccd[ix][iy] += luminosity 
	           * Kolmogorov_psf(r_locus, beta_PSF, fwhm)/norm;
    }

    return true;
}

// The exposure time is calubrated to the upper luminosity limit and the
// cluster distance.
// We assume that a 1Lsun star at 10pc just gets over exposed in 1. second.
local real calibrate_exposure_time(real upper_Llimit,
				   real distance_modulus) {

  real Mbol = 4.76;
  real magn = Mbol - 2.5*log10(upper_Llimit);
  real luminosity = pow(10, 0.4*(Mbol-magn));

  real t_exposure = 2./luminosity;

  return t_exposure;              // seconds

}

local void add_standard_stars(real **ccd, wavelength band, real beta_PSF,
			      real arcsec_per_pixel, real fwhm, bool verbose) {

  real Lu, Lb, Lv, Lr, Li;
  Lu= Lb= Lv= Lr= Li = VERY_LARGE_NUMBER;
  real logl = 0;    // 1 Lsun
  real logt = log10(cnsts.parameters(solar_temperature));
  real mass = 1;    // Msun

  ltm_to_ubvri(logl, logt, mass, Lu, Lb, Lv, Lr, Li);

  real luminosity = 0;
  switch(band) {
    case U: luminosity = pow(10, (4.76-Lu)/2.5);
	    break;
    case B: luminosity = pow(10, (4.76-Lb)/2.5);
	    break;
    case F: luminosity = pow(10, (4.76-Lv)/2.5);
	    break;
    case R: luminosity = pow(10, (4.76-Lr)/2.5);
	    break;
    case J: luminosity = pow(10, (4.76-Li)/2.5);
	    break;
    default:  cerr << "No Default value for magnitude"<<endl;
  }

  real xpos = 3*arcsec_per_pixel;
  real ypos = 3*arcsec_per_pixel;

  if (verbose) 
    cerr << "Add starndard star at pos (x, y): " 
         << xpos << " " << ypos << endl;

  add_star_to_ccd(ccd, luminosity, xpos, ypos, beta_PSF, 
		  arcsec_per_pixel, fwhm);

}

local real get_luminosity(dyn* bi, wavelength band,
			  real distance_modulus,
			  real luminosity_norm) {

  real Lu, Lb, Lv, Lr, Li;
  Lu= Lb= Lv= Lr= Li = VERY_LARGE_NUMBER;
  stellar_type stype = NAS;
  get_ubvri_star(bi, stype, Lu, Lb, Lv, Lr, Li);

  //	Lu= Lb= Lv= Lr= Li = VERY_LARGE_NUMBER;
  //	get_Lubvri_star(bi, stype, Lu, Lb, Lv, Lr, Li);
  // PRC(Lu);PRC(Lb);PRC(Lv);PRC(Lr);PRL(Li);
  //cerr << type_string(stype) << endl;

  real magn;
  real Msun;
  real air_mass = 0;
  switch(band) {
    case X: magn = Lu;
	    Msun = 5.61;
	    switch(stype) {
	      //	    case Carbon_Dwarf:
	      //	    case Helium_Dwarf:
	    case Oxygen_Dwarf:
	    case Xray_Pulsar:
	    case Radio_Pulsar:
	    case Neutron_Star:
	    case Black_Hole:
	      magn = Lu - 10;
	      Msun = 5.61;
	      PRC(magn);
		break;
	    default: 
	      magn = VERY_LARGE_NUMBER; // cut out all stars.
	    }

	    break;
    case U: magn = Lu;
	    Msun = 5.61;
	    break;
    case B: magn = Lb + 0.11*(Lb-Lv) - 0.05 + 0 * air_mass;
	    Msun = 5.48;
	    break;
    case F: magn = Lv + 0.05*(Lb-Lv) - 0.02 + 0 * air_mass; 
	    Msun = 4.83;
	    break;
    case R: magn = Lr;
	    break;
    case J: magn = Li;
	    break;
    default:  cerr << "No Default value for magnitude"<<endl;
  }

  real Mbol = 4.76;
  real luminosity = pow(10, 0.4*(Mbol-magn - distance_modulus));
  //cerr << "Magnitude = " << magn << " ";PRL(luminosity);

  real dL = gauss()*sqrt(luminosity*luminosity_norm);
  luminosity += dL/luminosity_norm;
  //  PRC(luminosity);PRL(dL);

  return luminosity;
}

local int add_background_stars_to_ccd(real** ccd, wavelength band,
				      real beta_PSF,
				      real arcsec_per_pixel, real fwhm,
				      bool verbose) {

  real luminosity = 0;
  real xpos, ypos;
  real ccd_size = XCCD_MAX*YCCD_MAX * pow(arcsec_per_pixel/3600., 2);
  PRL(ccd_size);

  real a = -5.86697E+00;
  real b = 9.95179E-01;
  real c = -3.74649E-02;
  real d = 6.47974E-04;
  real e = -4.38781E-06;

  int n_added=0;
  int ns_mag;
  real Mv;
  for(int imag=4; imag<31; imag++) {
    Mv = imag;
    ns_mag = (int)(ccd_size * 2.5
           * pow(10., a + b*Mv + c*pow(Mv, 2) + d*pow(Mv, 3) + e*pow(Mv, 4)));

    PRC(Mv);PRC(ns_mag);
    
    for(int i=0; i<ns_mag; i++) {
      Mv = randinter((real)imag, (real)imag+1);
      if(band!=2) {
	Mv += randinter(-0.2, 1.5);
      }
      luminosity = pow(10, (4.76-Mv)/2.5);
      xpos = randinter(0., XCCD_MAX*arcsec_per_pixel);
      ypos = randinter(0., YCCD_MAX*arcsec_per_pixel);
      if (verbose) {
	cerr << "Background star # " << n_added+1
	     << " pos (x, y): " 
	     << xpos << " " << ypos << " pixels"; 
	cerr << "         magnitudes (M): " << Mv << endl;
      }

      if(add_star_to_ccd(ccd, luminosity, xpos, ypos, beta_PSF, 
			 arcsec_per_pixel, fwhm)) 
	n_added++;
    }
  }
  return n_added;
} 

local void add_bad_column(real **ccd, int n_bc, bool verbose) {

  // reset input seed.
  int input_seed = 11111111;
  int actual_seed = srandinter(input_seed);
  if (verbose) {
    cerr << "Add " << n_bc << " bad columns: ";
    PRC(input_seed);PRL(actual_seed);
  }

  for(int iy = 0; iy<n_bc; iy++) {

    int iyy = int(randinter(0., YCCD_MAX));
    if (verbose)
      cerr << "Bad column at: " << iyy << endl;

    for(int ix = 0; ix<XCCD_MAX; ix++) 
      ccd[ix][iyy] = 0;
  }

}

local int add_cosmic_rays(real **ccd, real t_exposure, bool verbose) {

  if (verbose) 
    cerr << "Cosmic rays: " << endl;

  int n_cosmics = (int)abs(0.05*t_exposure * gauss());

  int x, y, cosm, n_cosm=0;
  int nx, ny;
  for (int i=0; i<n_cosmics; i++) {
      cosm = (int)abs(SATURATION_LIMIT*gauss());
      nx = int(randinter(0, 4));
      ny = int(randinter(0, 4));
      x = int(randinter(nx, XCCD_MAX-nx));
      y = int(randinter(ny, YCCD_MAX-ny));
      
      if (verbose) {
	cerr << "    Cosmic " << cosm << " from (x, y): " << x << " " << y 
             << " to: " << x+nx << " " << y+ny << endl; 
      }

      for(int ix = x-nx; ix<x+nx; ix++) 
	for(int iy = y-ny; iy<y+ny; iy++) {
	  ccd[ix][iy] = cosm;
	  n_cosm++;
	}
    }

  return n_cosm;  
}

local void saturate_pixel(real **ccd, int ix, int iy, int inc, 
			  real maximum_electrons) { 

  if(ccd[ix][iy]>maximum_electrons && ix<XCCD_MAX-1 && ix>0) {

    real excess = ccd[ix][iy]-maximum_electrons;
    ccd[ix][iy] -= excess;
    ccd[ix+inc][iy] += excess;
    saturate_pixel(ccd, ix+inc, iy, inc, maximum_electrons);
  }
  else 
    return;

}

local void create_overexposure_marks(real **ccd, real maximum_electrons) { 

  real excess = 0;
  for(int ix = 0; ix<XCCD_MAX; ix++) 
    for(int iy = 0; iy<YCCD_MAX; iy++)  {
//      excess = ccd[ix][iy]-maximum_electrons;
//      ccd[ix][iy] -= 0.9999*excess;
//      saturate_pixel(ccd, ix, iy, 1, maximum_electrons);
//      ccd[ix][iy] += 0.9999*excess;
      saturate_pixel(ccd, ix, iy, -1, maximum_electrons);
    }
}

local void mk_ccd(dyn* b, vec dc_pos, int project, wavelength band,
		  real distance, real reddening,
		  real upper_Llimit,
		  real xoffset, real yoffset,
		  real electrons_to_ADU, real beta_PSF,
		  real arcsec_per_pixel, real fwhm,
		  char* filename, bool add_background_stars,
		  bool add_standard_star, int input_seed,
		  bool add_cosmics, int nbad_column,
		  bool verbose) {

  PRL(band);

  real maximum_electrons = electrons_to_ADU*SATURATION_LIMIT-1;

  real half_xccd = 0.5*XCCD_MAX*arcsec_per_pixel;
  real half_yccd = 0.5*YCCD_MAX*arcsec_per_pixel;

  real **ccd = mkarray(XCCD_MAX, YCCD_MAX);
  cerr << "CCD allocated"<<endl;
  for(int ix = 0; ix<XCCD_MAX; ix++) 
    for(int iy = 0; iy<YCCD_MAX; iy++)  
      ccd[ix][iy] = 0;

  PRL(distance);

  real distance_modulus = 5 * log10(distance/10.);
  PRC(distance_modulus);
  PRL(distance);

//  magnitude_limit -= distance_modulus;
//  real magnitude_limit = t_obs/10. - distance_modulus;

//  real luminosity_limit = pow(10, (4.76-magnitude_limit)/5.5); // in V
//  PRC(magnitude_limit);
//  PRL(luminosity_limit);
  
  real pos_to_parsec = b->get_starbase()->conv_r_dyn_to_star(1)
                     * cnsts.parameters(Rsun)/cnsts.parameters(parsec);
  real pos_to_arcsec = b->get_starbase()->conv_r_dyn_to_star(1)
                     * 3600 * cnsts.parameters(Rsun)
                     / (distance*cnsts.parameters(parsec));
  

  PRC(pos_to_parsec);PRL(pos_to_arcsec);

  real xpos, ypos, zpos;
  real M_max = VERY_LARGE_NUMBER;
  real local_Mmm=0, magn;
  real luminosity, L_max = 0;
  real xmax=0, ymax=0;
  vec r;
  if (b->get_oldest_daughter()) {

    vec r_com  = b->get_pos() - dc_pos;

    for_all_leaves(dyn, b, bb) {

//      cerr << " new star: ";
//      cerr << bb->get_index()<<endl;

      r = bb->get_pos()-bb->get_parent()->get_pos();
//      PRL(r);

      switch(project) {
        case -3: xpos =  r[0];
		 ypos =  r[1];
		 zpos = -r[2];
		 break;
	case -2: xpos = -r[0];
		 ypos =  r[2];
		 zpos = -r[1];
		 break;
	case -1: xpos = -r[1];
		 ypos =  r[2];
		 zpos = -r[0];
		 break;
	case 1: xpos = r[1];
		ypos = r[2];
		zpos = r[0];
		break;
	case 2: xpos = r[0];
		ypos = r[2];
		zpos = r[1];
		break;
        case 3: xpos = r[0];
		ypos = r[1];
		zpos = r[2];
		break;
      }
      xpos = xpos * pos_to_arcsec + half_xccd - xoffset;
      ypos = ypos * pos_to_arcsec + half_yccd - yoffset;
      zpos *= pos_to_parsec;
//      PRC(zpos);
      if(distance+zpos>0) {
	real local_distance_modulus = 5 * log10((distance+zpos)/10.);
//      PRL(local_distance_modulus);

	luminosity = get_luminosity(bb, band, local_distance_modulus,
				    maximum_electrons/upper_Llimit);

//	PRL(luminosity);
	if(L_max<luminosity) {
	  L_max = luminosity;
	  M_max = 4.76 - 2.5*log10(luminosity) - local_distance_modulus;
	  local_Mmm = local_distance_modulus;
	  xmax = xpos;    //ascsec
	  ymax = ypos;
	}

	if (luminosity<=0) {
	  luminosity = 0;
	  //	  cerr << "Stellar luminosity = " << luminosity << endl;
	  //	  cerr << "        star Not added to ccd"<<endl;
	}
	else {
	  if (verbose) {
	    real nx = xpos/arcsec_per_pixel;
	    real ny = ypos/arcsec_per_pixel;
	    if (nx>=-10 && nx<=XCCD_MAX && 
		ny>=-10 && ny<=YCCD_MAX+10) {
	      cerr << "Cluster star #" << bb->get_index() 
		   << " pos (x, y): " 
	           << nx << " " << ny << " pixels" << endl; 
	      cerr << "             magnitudes (M, m): "  
		   << 4.76 - 2.5*log10(luminosity)  - local_distance_modulus
 	           << " " << 4.76 - 2.5*log10(luminosity) <<endl;
	    }
	  }
	  add_star_to_ccd(ccd, luminosity, xpos, ypos, beta_PSF, 
			  arcsec_per_pixel, fwhm);
	}
      }
    }

    cerr << "standard (brightest) star: \n"; 
    cerr << "L (app) = " << L_max << ", Mag (abs) = " << M_max << endl;
    cerr << "          local (M-m) = " << local_Mmm << endl;
    cerr << "    position in CCD = (" << floor(xmax/arcsec_per_pixel) << ", "
                                      << floor(ymax/arcsec_per_pixel) 
         << ") pixels"<<endl;
    cerr << "    position on sky = (" << xmax << ", "<<ymax 
         << ") arcsec."<<endl; 

    real t_exposure = calibrate_exposure_time(upper_Llimit,
					      distance_modulus);
    PRL(t_exposure);

    if(add_background_stars) {
      cerr << "add background " <<endl;
      int n_backgr = add_background_stars_to_ccd(ccd, band, beta_PSF,
						 arcsec_per_pixel, fwhm,
						 verbose); 
      PRL(n_backgr);
    }
    else {
      cerr <<"no background stars added"<<endl;
    }

    if(add_standard_star) {
      cerr << "add standard star " <<endl;
      add_standard_stars(ccd, band, beta_PSF, arcsec_per_pixel, fwhm, verbose);
    }
    else {
      cerr <<"no standard star added" <<endl;
    }

    // initialize random processes.
    int actual_seed;
    if(input_seed==0) actual_seed = 0;
    actual_seed = srandinter(input_seed);
    PRC(input_seed);PRL(actual_seed);

    real Mmax = 0;
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  
	Mmax = Starlab::max(Mmax, ccd[ix][iy]);

    PRC(Mmax);PRL(upper_Llimit);
//    real upper_limit_counts = calibrate_upper_limit(beta_PSF, upper_Llimit);
//    PRL(upper_limit_counts);
    // Normalize and account for dynamic range of log10(4)
    real lmax = 0;
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  {
	ccd[ix][iy] = ccd[ix][iy]/upper_Llimit;
	lmax = Starlab::max(lmax, ccd[ix][iy]);
      }
    PRL(lmax);

    Mmax = 0;
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  {
	ccd[ix][iy] *= maximum_electrons;
	Mmax = Starlab::max(Mmax, ccd[ix][iy]);
      }
    PRL(Mmax);

    // add Scintilation noise of exposure duration
    real scintilation = 10./sqrt(t_exposure); // electrons
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  
	ccd[ix][iy] += scintilation;

    // add dark current and readout
    real dark_current = 50*t_exposure;  // electrons
    real sky_background = 150;    // Gilliland et al. 1995, ApJ 447, 191
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  {
	ccd[ix][iy] += dark_current;
	ccd[ix][iy] += sky_background;
      }

    // add Electroluminescence only on very cheap ccds
//    real r_origin, electroluminescence = 10*t_exposure;
//  real diagonal = sqrt(pow(XCCD_MAX, 2) + pow(YCCD_MAX, 2));
//    for(int ix = 0; ix<XCCD_MAX; ix++) 
//      for(int iy = 0; iy<YCCD_MAX; iy++)  {
//	r_origin = sqrt(pow(ix, 2) + pow(iy, 2))
//                 / diagonal;
//	ccd[ix][iy] += electroluminescence*pow(1-r_origin, 4);
//      }

    // add Thermal luminosity
//    real readout = 200;
//    for(int ix = 0; ix<XCCD_MAX; ix++) 
//      for(int iy = 0; iy<YCCD_MAX; iy++)  {
//	r_origin = (YCCD_MAX-iy)/YCCD_MAX;
//	ccd[ix][iy] += readout*(1-r_origin);
//      }

    // add Poissonian noise
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  
	ccd[ix][iy] += sqrt(ccd[ix][iy])*gauss();
//	ccd[ix][iy] += abs(sqrt(ccd[ix][iy])*gauss());

    // add Bleeding and Blooming
    create_overexposure_marks(ccd, maximum_electrons);

    // convert electrons to ADUs.
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  
	ccd[ix][iy] /= electrons_to_ADU;

    // cut off too low and too high values.
    for(int ix = 0; ix<XCCD_MAX; ix++) 
      for(int iy = 0; iy<YCCD_MAX; iy++)  
	ccd[ix][iy] = Starlab::min(Starlab::max(0., ccd[ix][iy]), 1.*SATURATION_LIMIT);

    if(add_cosmics) {
      int n_cosm = add_cosmic_rays(ccd, t_exposure, verbose);
      PRL(n_cosm);
    }
    else
      cerr << "No cosmics added"<<endl;

    int n_bc = 0;
    if(nbad_column>0) {
      add_bad_column(ccd, nbad_column, verbose);
    }
    else
      cerr << "No bad columns added"<<endl;

    print_ccd(ccd, filename, t_exposure);
  }
}

// void mk_ccd.C
//
// Options:    -A    number of electrons to ADU [3]
//             -B    filter of observation [3]
//                   (options are: U=0, B=1, V=2, R=3, I=4) 
//             -b    add background stars to ccd [false]
//             -c    add cosmics rays to ccd [false]
//             -d    distance to cluster in parsec [1000]
//             -F    output filename [ccd.fits]
//             -L    upper luminosity limit in Lsun [1]
//             -N    add number of bad columns [0]
//             -P    Kolmogorov point spread function beta [5]
//             -p    give projecttion [1]
//                       -3: project on xy plane viewed from -z
//                       -2: project on xz plane viewed from -y
//                       -1: project on yz plane viewed from -x
//                        1: project on yz plane viewed from  x
//                        2: project on xz plane viewed from  y
//                        3: project on xy plane viewed from  z
//             -R    reddening [0]
//             -S    random number seed [clock]
//             -s    add standard stars [false];
//             -x    x-offset from density center [0]
//             -y    y-offset from density center [0]
//             -v    verbose, prints stars and positions to cerr [false]
//
//
main(int argc, char ** argv) {
    check_help();

    wavelength band = F;

    bool verbose = false;

    int input_seed=0;
    int nbad_column = 0;

    //    int nx_bin = 512; // one pixel = 0.076 arcsec
    //    int ny_bin = 512;
    int nx_bin = XCCD_MAX; // one pixel = 0.076 arcsec
    int ny_bin = YCCD_MAX;

    bool add_background_stars = false;
    bool add_standard_star = false;
    bool add_cosmics = false;

    real arcsec_per_pixel = 0.3;
    real fwhm  = 0.7;

    real distance = 1000;   //pc
    real reddening = 0;
    real electrons_to_ADU = 3;     // standard for modern CCD
    real beta_PSF = 5;             // standard PSF vary between 3 and 10

    real upper_luminosity_limit = 1; //Lsun

    real xoffset = 0;  // offset from cluster density center in arcsec
    real yoffset = 0;  // offset from cluster density center in arcsec

    int project = 1; // project along x(1), y(2), z(3)

    char* filename = "ccd.fits";
    extern char *poptarg;
    int c;
    char* param_string = "A:a:B:bcx:y:d:F:f:p:P:R:sS:L:x:y:N:v";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
        switch(c)
            {
            case 'A': electrons_to_ADU = atof(poptarg);
                      break;
            case 'a': arcsec_per_pixel = atof(poptarg);
                      break;
            case 'f': fwhm = atof(poptarg);
                      break;
            case 'B': band = (wavelength)atoi(poptarg);
                      break;
            case 'b': add_background_stars = true;
                      break;
            case 'c': add_cosmics = true;
                      break;
            case 'p': project = atoi(poptarg);
                      break;
            case 'd': distance = atof(poptarg);
                      break;
            case 'P': beta_PSF = atof(poptarg);
                      break;
            case 'R': reddening = atof(poptarg);
                      break;
            case 'S': input_seed = atoi(poptarg);
                      break;
            case 's': add_standard_star = true;
                      break;
            case 'N': nbad_column = atoi(poptarg);
                      break;
            case 'L': upper_luminosity_limit = atof(poptarg);
                      break;
            case 'F': filename = poptarg;
                      break;
            case 'x': xoffset = atof(poptarg);   
                      break;
            case 'y': yoffset = atof(poptarg);
                      break;
            case 'v': verbose = true;
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
                      get_help();
                      exit(1);
            }

    real distance_modulus = 5 * log10(distance/10.);
    real Mbol = 4.76;
    real magn = Mbol - 2.5*log10(upper_luminosity_limit) + distance_modulus;
    upper_luminosity_limit = pow(10, 0.4*(Mbol-magn));
    PRC(distance_modulus);PRC(magn);PRL(upper_luminosity_limit);


//    ofstream outfile(filename, ios::out);
//    ofstream outfile(filename, ios::binary);
//    if (!outfile) {
//      cerr << "error: couldn't crate file "<< filename<<endl;
//      exit(0);
//    }

    dyn *b = NULL;
    int count = 0;
    bool cod, try_com = false;
    vec dc_pos = 0;
    while (b = get_dyn()) {

	b->flatten_node();

        cout << "\n\nTime = "
             << b->get_starbase()->conv_t_dyn_to_star(b->get_system_time())
             << endl;

        if (find_qmatch(b->get_dyn_story(), "density_center_pos")) {
	  cerr << "found density_center_pos"<<endl;

           if (getrq(b->get_dyn_story(), "density_center_time")
               != b->get_system_time()) {
               warning("mkpovfile: neglecting out-of-date density center");
               try_com = true;
            } else
            cod = true;

            dc_pos = getvq(b->get_dyn_story(), "density_center_pos");

	  PRL(dc_pos);
         }

         if (try_com && find_qmatch(b->get_dyn_story(), "com_pos")) {
	  cerr << "found com_pos"<<endl;

            if (getrq(b->get_dyn_story(), "com_time")
               != b->get_system_time()) {

               warning("lagrad: neglecting out-of-date center of mass");
            } else
                dc_pos = getvq(b->get_dyn_story(), "com_pos");

	  PRL(dc_pos);
         }

	cerr << "Create CCD"<<endl;
        mk_ccd(b, dc_pos, project, band, distance, reddening, upper_luminosity_limit, xoffset, yoffset, electrons_to_ADU, beta_PSF, arcsec_per_pixel, fwhm, filename, add_background_stars, add_standard_star, input_seed, add_cosmics, nbad_column, verbose);

       rmtree(b);
    }
}

#endif
