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

#include "vector3.h"
#include "taylor.h"


#ifdef _QD_QD_REAL_H
static inline m_real div(double a, double b){
    return  m_real::accurate_div(m_real(a),m_real(b));
}
#endif // _QD_QD_REAL_H

struct Particle{
	enum
	{
		NFORCE = 8,
		ORDER  = 2*NFORCE,
	};
	unsigned long id;
	m_real dt;
	m_real tlast;
	m_real mass;
	m_real radius;

	qvec3 coord[2 + ORDER]; // pos-vel, forces, interpolates
	qvec3 add_acc;
};

struct Predictor{
	m_real mass;
	qvec3 coord[8];

	void predict(const m_real tnext, const Particle &p){
		const m_real dt = tnext - p.tlast;
		mass = p.mass;
		coord[0] = taylor<15>(dt, &p.coord[0]);
		coord[1] = taylor<14>(dt, &p.coord[1]);
		coord[2] = taylor<13>(dt, &p.coord[2]);
		coord[3] = taylor<12>(dt, &p.coord[3]);
		coord[4] = taylor<11>(dt, &p.coord[4]);
		coord[5] = taylor<10>(dt, &p.coord[5]);
		coord[6] = taylor< 9>(dt, &p.coord[6]);
		coord[7] = taylor< 8>(dt, &p.coord[7]);
	}
	void zero_predict(const Particle &p){
		mass     = p.mass;
		coord[0] = p.coord[0];
		coord[1] = p.coord[1];
		coord[2] = p.coord[2];
		coord[3] = p.coord[3];
		coord[4] = p.coord[4];
		coord[5] = p.coord[5];
		coord[6] = p.coord[6];
		coord[7] = p.coord[7];
	}
};

struct Force{
	qvec3 force[8]; // acc, jerk, snap, crackle, pop, pdot, pdot2, pdot3
	void init_assign(Particle &p){
		for(int m=0; m<Particle::NFORCE; m++){
			p.coord[2+m] = force[m];
		}
	}
};

struct Corrector{
	qvec3 pos;
	qvec3 vel;
	qvec3 force[8];
	qvec3 fintp[8];

	void correct(const Particle &p, const Force &f){
		const m_real h = 0.5 * p.dt;
		const m_real c0 = 0.5;
		const m_real c1 = c0 * h;
		const m_real c2 = c1 * (div(1., 2.) * h);
		const m_real c3 = c2 * (div(1., 3.) * h);
		const m_real c4 = c3 * (div(1., 4.) * h);
		const m_real c5 = c4 * (div(1., 5.) * h);
		const m_real c6 = c5 * (div(1., 6.) * h);
		const m_real c7 = c6 * (div(1., 7.) * h);

		const m_real d0 = c0;
		const m_real d1 = div( 7.,   15.) * c1;
		const m_real d2 = div(12.,   45.) * c2;
		const m_real d3 = div( 2.,   13.) * c3;
		const m_real d4 = div(16.,  195.) * c4;
		const m_real d5 = div(16.,  429.) * c5;
		const m_real d6 = div(64., 5005.) * c6;
		const m_real d7 = div(16., 6435.) * c7;

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// correct velocity
		const qvec3 f0pl = d0 * (fr[0] + fl[0]);
		const qvec3 f1mn = d1 * (fr[1] - fl[1]);
		const qvec3 f2pl = d2 * (fr[2] + fl[2]);
		const qvec3 f3mn = d3 * (fr[3] - fl[3]);
		const qvec3 f4pl = d4 * (fr[4] + fl[4]);
		const qvec3 f5mn = d5 * (fr[5] - fl[5]);
		const qvec3 f6pl = d6 * (fr[6] + fl[6]);
		const qvec3 f7mn = d7 * (fr[7] - fl[7]);

		this->vel = p.coord[1] + p.dt * (f0pl - (f1mn - (f2pl - (f3mn - (f4pl - (f5mn - (f6pl - f7mn)))))));

		// correct position
		const qvec3 v0pl = d0 * (this->vel + p.coord[1]);
		const qvec3 v1mn = d1 * (fr[0] - fl[0]);
		const qvec3 v2pl = d2 * (fr[1] + fl[1]);
		const qvec3 v3mn = d3 * (fr[2] - fl[2]);
		const qvec3 v4pl = d4 * (fr[3] + fl[3]);
		const qvec3 v5mn = d5 * (fr[4] - fl[4]);
		const qvec3 v6pl = d6 * (fr[5] + fl[5]);
		const qvec3 v7mn = d7 * (fr[6] - fl[6]);

		this->pos = p.coord[0] + p.dt * (v0pl - (v1mn - (v2pl - (v3mn - (v4pl - (v5mn - (v6pl - v7mn)))))));

		// copy forces
		for(int i=0; i<8; i++) force[i] = f.force[i];
	}

	void interpolate(const Particle &p, const Force &f){
		qvec3 fpl[8];
		qvec3 fmn[8];
		qvec3 fmid[16];
		// qvec3 fleft [16];
		qvec3 fright[16];

		const m_real h = 0.5 * p.dt;
		const m_real c0 = 0.5;
		const m_real c1 = c0 * h;
		const m_real c2 = c1 * (div(1., 2.) * h);
		const m_real c3 = c2 * (div(1., 3.) * h);
		const m_real c4 = c3 * (div(1., 4.) * h);
		const m_real c5 = c4 * (div(1., 5.) * h);
		const m_real c6 = c5 * (div(1., 6.) * h);
		const m_real c7 = c6 * (div(1., 7.) * h);

		const qvec3 *fr = f.force;
		const qvec3 *fl = p.coord + 2;

		// f+, f-
		{
			fpl[0] = c0 * (fr[0] + fl[0]);
			fmn[0] = c0 * (fr[0] - fl[0]);
			fpl[1] = c1 * (fr[1] + fl[1]);
			fmn[1] = c1 * (fr[1] - fl[1]);
			fpl[2] = c2 * (fr[2] + fl[2]);
			fmn[2] = c2 * (fr[2] - fl[2]);
			fpl[3] = c3 * (fr[3] + fl[3]);
			fmn[3] = c3 * (fr[3] - fl[3]);
			fpl[4] = c4 * (fr[4] + fl[4]);
			fmn[4] = c4 * (fr[4] - fl[4]);
			fpl[5] = c5 * (fr[5] + fl[5]);
			fmn[5] = c5 * (fr[5] - fl[5]);
			fpl[6] = c6 * (fr[6] + fl[6]);
			fmn[6] = c6 * (fr[6] - fl[6]);
			fpl[7] = c7 * (fr[7] + fl[7]);
			fmn[7] = c7 * (fr[7] - fl[7]);
		}
#if 0
		// fmid
		{
			qvec3 evn[8], odd[8];

			evn[1] = fmn[1] - 2.*fpl[2] + 3.*fmn[3] - 4.*fpl[4] +  5.*fmn[5] -  6.*fpl[6] +  7.*fmn[7];
			evn[3] =                         fmn[3] - 4.*fpl[4] + 10.*fmn[5] - 20.*fpl[6] + 35.*fmn[7];
			evn[5] =                                                  fmn[5] -  6.*fpl[6] + 21.*fmn[7];
			evn[7] =                                                                            fmn[7];

			odd[0] = fmn[0] - fpl[1] +    fmn[2] -    fpl[3] +    fmn[4] -     fpl[5] +     fmn[6] -     fpl[7];
			odd[2] =                      fmn[2] - 3.*fpl[3] + 6.*fmn[4] - 10.*fpl[5] + 15.*fmn[6] - 21.*fpl[7];
			odd[4] =                                              fmn[4] -  5.*fpl[5] + 15.*fmn[6] - 35.*fpl[7];
			odd[6] =                                                                        fmn[6] -  7.*fpl[7];

			fmid[ 8] = (1./2048.) * (-2145.*evn[1] + 693.*evn[3] - 525.*evn[5] + 1225.*evn[7]);
			fmid[10] = (1./2048.) * ( 1001.*evn[1] - 297.*evn[3] + 189.*evn[5] -  245.*evn[7]);
			fmid[12] = (1./2048.) * ( -273.*evn[1] +  77.*evn[3] -  45.*evn[5] +   49.*evn[7]);
			fmid[14] = (1./2048.) * (   33.*evn[1] -   9.*evn[3] +   5.*evn[5] -    5.*evn[7]);

			fmid[ 9] = (1./2048.) * ( 25025.*odd[0] - 2145.*odd[2] + 693.*odd[4] - 525.*odd[6]);
			fmid[11] = (1./2048.) * (-12285.*odd[0] + 1001.*odd[2] - 297.*odd[4] + 189.*odd[6]);
			fmid[13] = (1./2048.) * (  3465.*odd[0] -  273.*odd[2] +  77.*odd[4] -  45.*odd[6]);
			fmid[15] = (1./2048.) * (  -429.*odd[0] +   33.*odd[2] -   9.*odd[4] +   5.*odd[6]);
		}
#else
		// even
		{
			qvec3 tmp  = fpl[2] - (1./2.) * fmn[1];

			// lower triangle
			qvec3 evn0 = (1./16.)  * ((5./4.)   * tmp + (-3./2.) * fmn[3] + fpl[4]);
			qvec3 evn1 = (1./32.)  * ((-7./4.)  * tmp + (9./4)   * fmn[3] + (-2.0)  * fpl[4] +            fmn[5]);
			qvec3 evn2 = (1./64.)  * ((21./8.)  * tmp + (-7./2.) * fmn[3] + (7./2.) * fpl[4] + (-5./2.) * fmn[5] +          fpl[6]);
			qvec3 evn3 = (1./128.) * ((-33./8.) * tmp + (45./8.) * fmn[3] + (-6.0)  * fpl[4] + (5.0)    * fmn[5] + (-3.0) * fpl[6] + fmn[7]);

			// upper triangle
			fmid[8]  = evn0 - 5.0 * evn1 + 15.0 * evn2 - 35.0 * evn3;
			fmid[10] =              evn1 -  6.0 * evn2 + 21.0 * evn3;
			fmid[12] =                            evn2 -  7.0 * evn3;
			fmid[14] =                                          evn3;
		}
		// odd
		{
			qvec3 tmp = fpl[1] - fmn[0];

			// lower triangle
			qvec3 odd0 = (1./16.)  * ((-35./8.)   * tmp + (15./4.)  * fmn[2] + (-5./2.)  * fpl[3] +           fmn[4]);
			qvec3 odd1 = (1./32.)  * ((63./8)     * tmp + (-7.0)    * fmn[2] + (21./4)   * fpl[3] + (-3.0)  * fmn[4] +            fpl[5]);
			qvec3 odd2 = (1./64.)  * ((-231./16.) * tmp + (105./8.) * fmn[2] + (-21./2.) * fpl[3] + (7.0)   * fmn[4] + (-7./2.) * fpl[5] +          fmn[6]);
			qvec3 odd3 = (1./128.) * ((429./16.)  * tmp + (-99./4.) * fmn[2] + (165./8.) * fpl[3] + (-15.0) * fmn[4] + (9.0)    * fpl[5] + (-4.0) * fmn[6] + fpl[7]);

			// upper triangle
			fmid[9]  = odd0 - 5.0 * odd1 + 15.0 * odd2 - 35.0 * odd3;
			fmid[11] =              odd1 -  6.0 * odd2 + 21.0 * odd3;
			fmid[13] =                            odd2 -  7.0 * odd3;
			fmid[15] =                                          odd3;
		}
#endif
		fright[ 8] = fmid[8] + 9.*fmid[9] +  45.*fmid[10] + 165.*fmid[11] + 495.*fmid[12] + 1287.*fmid[13] + 3003.*fmid[14] + 6435.*fmid[15];
		fright[ 9] =              fmid[9] +  10.*fmid[10] +  55.*fmid[11] + 220.*fmid[12] +  715.*fmid[13] + 2002.*fmid[14] + 5005.*fmid[15];
		fright[10] =                             fmid[10] +  11.*fmid[11] +  66.*fmid[12] +  286.*fmid[13] + 1001.*fmid[14] + 3003.*fmid[15];
		fright[11] =                                             fmid[11] +  12.*fmid[12] +   78.*fmid[13] +  364.*fmid[14] + 1365.*fmid[15];
		fright[12] =                                                             fmid[12] +   13.*fmid[13] +   91.*fmid[14] +  455.*fmid[15];
		fright[13] =                                                                              fmid[13] +   14.*fmid[14] +  105.*fmid[15];
		fright[14] =                                                                                               fmid[14] +   15.*fmid[15];
		fright[15] =                                                                                                                fmid[15];
		// rescale
		{
			const m_real hi = 1.0 / h;
			const m_real s1 = hi;
			const m_real s2 = s1 * (2.*hi);
			const m_real s3 = s2 * (3.*hi);
			const m_real s4 = s3 * (4.*hi);
			const m_real s5 = s4 * (5.*hi);
			const m_real s6 = s5 * (6.*hi);
			const m_real s7 = s6 * (7.*hi);
			const m_real s8 = s7 * (8.*hi);
			const m_real s9 = s8 * (9.*hi);
			const m_real s10= s9 *(10.*hi);
			const m_real s11= s10*(11.*hi);
			const m_real s12= s11*(12.*hi);
			const m_real s13= s12*(13.*hi);
			const m_real s14= s13*(14.*hi);
			const m_real s15= s14*(15.*hi);

			fintp[0] = s8 * fright[8];
			fintp[1] = s9 * fright[9];
			fintp[2] = s10* fright[10];
			fintp[3] = s11* fright[11];
			fintp[4] = s12* fright[12];
			fintp[5] = s13* fright[13];
			fintp[6] = s14* fright[14];
			fintp[7] = s15* fright[15];
		}
	}
	void commit(Particle &p){
		p.coord[0] = pos;
		p.coord[1] = vel;

		p.coord[2] = force[0];
		p.coord[3] = force[1];
		p.coord[4] = force[2];
		p.coord[5] = force[3];
		p.coord[6] = force[4];
		p.coord[7] = force[5];
		p.coord[8] = force[6];
		p.coord[9] = force[7];

		p.coord[10] = fintp[0];
		p.coord[11] = fintp[1];
		p.coord[12] = fintp[2];
		p.coord[13] = fintp[3];
		p.coord[14] = fintp[4];
		p.coord[15] = fintp[5];
		p.coord[16] = fintp[6];
		p.coord[17] = fintp[7];

		p.tlast += p.dt;
	}
	void feedback(Predictor &pr){
		pr.coord[0] = pos;
		pr.coord[1] = vel;
		pr.coord[2] = force[0]; // acc
		pr.coord[3] = force[1]; // jrk
		pr.coord[4] = force[2]; // snp
		pr.coord[5] = force[3]; // crk
		pr.coord[6] = force[4]; // pop
		pr.coord[7] = force[5]; // pdot
	}
	void taylor_test(const Particle &p){
		dvec3 df0 = p.coord[2] - taylor<15>(-p.dt, &force[0]);
		dvec3 df1 = p.coord[3] - taylor<14>(-p.dt, &force[1]);
		dvec3 df2 = p.coord[4] - taylor<13>(-p.dt, &force[2]);
		dvec3 df3 = p.coord[5] - taylor<12>(-p.dt, &force[3]);
		dvec3 df4 = p.coord[6] - taylor<11>(-p.dt, &force[4]);
		dvec3 df5 = p.coord[7] - taylor<10>(-p.dt, &force[5]);
		dvec3 df6 = p.coord[8] - taylor< 9>(-p.dt, &force[6]);
		dvec3 df7 = p.coord[9] - taylor< 8>(-p.dt, &force[7]);
		printf("[0]: %e %e %e\n", df0.x, df0.y, df0.z);
		printf("[1]: %e %e %e\n", df1.x, df1.y, df1.z);
		printf("[2]: %e %e %e\n", df2.x, df2.y, df2.z);
		printf("[3]: %e %e %e\n", df3.x, df3.y, df3.z);
		printf("[4]: %e %e %e\n", df4.x, df4.y, df4.z);
		printf("[5]: %e %e %e\n", df5.x, df5.y, df5.z);
		printf("[6]: %e %e %e\n", df6.x, df6.y, df6.z);
		printf("[7]: %e %e %e\n", df7.x, df7.y, df7.z);
	}
};

