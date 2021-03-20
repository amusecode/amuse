/*MIT License

Copyright (c) 2021 Keigo Nitadori

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

/*
-added support for QD in addition to DD by Thomas Schano  QD is available at http://crd-legacy.lbl.gov/~dhbailey/mpdist/ (https://github.com/GFTwrt/QD)
-added support for additional step acceleration by Thomas Schano
*/

template <class SYS, typename real_type, typename vect_type>
void calc_force_on_i_p8(SYS &sys, const int nbody, const int i, const real_type eps2) {
	vect_type acc0(real_type(0));
	vect_type acc1(real_type(0));
	vect_type acc2(real_type(0));
	vect_type acc3(real_type(0));
	vect_type acc4(real_type(0));
	vect_type acc5(real_type(0));
	vect_type acc6(real_type(0));
	vect_type acc7(real_type(0));
	vect_type icoord0 = sys.template get_coord<0>(i);
	vect_type icoord1 = sys.template get_coord<1>(i);
	vect_type icoord2 = sys.template get_coord<2>(i);
	vect_type icoord3 = sys.template get_coord<3>(i);
	vect_type icoord4 = sys.template get_coord<4>(i);
	vect_type icoord5 = sys.template get_coord<5>(i);
	vect_type icoord6 = sys.template get_coord<6>(i);
	vect_type icoord7 = sys.template get_coord<7>(i);
	for(int j=0; j<nbody; j++){
		vect_type dr0 = sys.template get_coord<0>(j) - icoord0;
		vect_type dr1 = sys.template get_coord<1>(j) - icoord1;
		vect_type dr2 = sys.template get_coord<2>(j) - icoord2;
		vect_type dr3 = sys.template get_coord<3>(j) - icoord3;
		vect_type dr4 = sys.template get_coord<4>(j) - icoord4;
		vect_type dr5 = sys.template get_coord<5>(j) - icoord5;
		vect_type dr6 = sys.template get_coord<6>(j) - icoord6;
		vect_type dr7 = sys.template get_coord<7>(j) - icoord7;

		real_type s0 = eps2 + (dr0*dr0);
		if(eps2==s0) continue;
		real_type s1 = (dr0*dr1);
		real_type s2 = (dr0*dr2) +(dr1*dr1);
		real_type s3 = (dr0*dr3) +3.0*(dr1*dr2);
		real_type s4 = (dr0*dr4) +4.0*(dr1*dr3) +3.0*(dr2*dr2);
		real_type s5 = (dr0*dr5) +5.0*(dr1*dr4) +10.0*(dr2*dr3);
		real_type s6 = (dr0*dr6) +6.0*(dr1*dr5) +15.0*(dr2*dr4) +10.0*(dr3*dr3);
		real_type s7 = (dr0*dr7) +7.0*(dr1*dr6) +21.0*(dr2*dr5) +35.0*(dr3*dr4);

#if 1
		real_type rinv1 = rsqrt(s0);
		real_type rinv2 = rinv1 * rinv1;
#else
		real_type rinv2 = 1.0 / s0;
		real_type rinv1 = sqrt(rinv2);
#endif
		real_type rinv3 = rinv1 * rinv2;
		rinv3 *= sys.get_mass(j);

		real_type q1 = rinv2 * (-3.0*s1);
		real_type q2 = rinv2 * (-3.0*s2 - (5.0*s1)*q1);
		real_type q3 = rinv2 * (-3.0*s3 - (8.0*s2)*q1 - (7.0*s1)*q2);
		real_type q4 = rinv2 * (-3.0*s4 - (11.0*s3)*q1 - (15.0*s2)*q2 - (9.0*s1)*q3);
		real_type q5 = rinv2 * (-3.0*s5 - (14.0*s4)*q1 - (26.0*s3)*q2 - (24.0*s2)*q3 - (11.0*s1)*q4);
		real_type q6 = rinv2 * (-3.0*s6 - (17.0*s5)*q1 - (40.0*s4)*q2 - (50.0*s3)*q3 - (35.0*s2)*q4 - (13.0*s1)*q5);
		real_type q7 = rinv2 * (-3.0*s7 - (20.0*s6)*q1 - (57.0*s5)*q2 - (90.0*s4)*q3 - (85.0*s3)*q4 - (48.0*s2)*q5 - (15.0*s1)*q6);

		acc0 += rinv3 * (dr0);
		acc1 += rinv3 * (dr1 +(q1)*dr0);
		acc2 += rinv3 * (dr2 +(2.0*q1)*dr1 +(q2)*dr0);
		acc3 += rinv3 * (dr3 +(3.0*q1)*dr2 +(3.0*q2)*dr1 +(q3)*dr0);
		acc4 += rinv3 * (dr4 +(4.0*q1)*dr3 +(6.0*q2)*dr2 +(4.0*q3)*dr1 +(q4)*dr0);
		acc5 += rinv3 * (dr5 +(5.0*q1)*dr4 +(10.0*q2)*dr3 +(10.0*q3)*dr2 +(5.0*q4)*dr1 +(q5)*dr0);
		acc6 += rinv3 * (dr6 +(6.0*q1)*dr5 +(15.0*q2)*dr4 +(20.0*q3)*dr3 +(15.0*q4)*dr2 +(6.0*q5)*dr1 +(q6)*dr0);
		acc7 += rinv3 * (dr7 +(7.0*q1)*dr6 +(21.0*q2)*dr5 +(35.0*q3)*dr4 +(35.0*q4)*dr3 +(21.0*q5)*dr2 +(7.0*q6)*dr1 +(q7)*dr0);
	} // for(j)
	acc0 += sys.get_add_acc(i);
	sys.template set_force<0>(i, acc0);
	sys.template set_force<1>(i, acc1);
	sys.template set_force<2>(i, acc2);
	sys.template set_force<3>(i, acc3);
	sys.template set_force<4>(i, acc4);
	sys.template set_force<5>(i, acc5);
	sys.template set_force<6>(i, acc6);
	sys.template set_force<7>(i, acc7);
} // end of function
