void initialize(int n);

void set_jp(int add, double mass, double pos[3], double vel[3], double acc[3],
	    double jrk[3], double snp[3], double crk[3], double time, int id,
	    double eps2);

void predict_all(double time, int nj);
//void no_predict_all(in nj); additional

void pick_up_predictor(int add, predictor &pred);

void calc_force_on_predictors(int ni, predictor pred[], Force force[], int nj);
void calc_force_on_predictors_epsj(int ni, predictor pred[], Force force[], int nj);

// the definitions for sse2 are in "yesisu6-sse2.h"

struct predictor{
                double pos [3]; //  6
                double vel [3]; // 12
                double acc [3]; // 18
                double mass;    // 20
                double id  ;    // 22
                //double pad [2];
                double eps2;
                    predictor(){
                            assert(sizeof(*this) == 12*8);
                    }
            };  

struct Force{
    //static const int nword = 9;
    vec acc;
    vec jrk;
    vec snp;
    double phi;
    int nnb_id;
    float  nnb_r2;
    Force() : acc(0.0), jrk(0.0), snp(0.0), phi(0.0), nnb_id(-1), nnb_r2(HUGE){}    void operator += (Force &rhs){
        acc += rhs.acc;
        jrk += rhs.jrk;
        snp += rhs.snp;
        phi += rhs.phi;
    }
    void clear(){
        acc[0] = acc[1] = acc[2] = 0.0;
        jrk[0] = jrk[1] = jrk[2] = 0.0;
        snp[0] = snp[1] = snp[2] = 0.0;
        phi = 0.0;
        nnb_id = -1;
        nnb_r2 = HUGE;

    }
    
    
  //// for BH particles ////
    void add_BH_force(double eps2, double mBH, vec pos_BH, vec vel_BH, vec acc_BH,
                      vec pos_i, vec vel_i, vec acc_i, int iBH){

            vec pos0 = pos_BH - pos_i;
            vec vel0 = vel_BH - vel_i;
            vec acc0 = acc_BH - acc_i;
            double dist2 = pos0*pos0;
            double r2 = dist2 + eps2;
            double rinv2 = 1.0 / r2;
            double rinv = sqrt(rinv2);
            double rinv3 = rinv*rinv2;

            double alpha = rinv2*(vel0*pos0);
            double beta = rinv2*(vel0*vel0 + pos0*acc0) + alpha*alpha;

            vec acc1 = mBH * rinv3 * pos0;
            vec jerk1 = mBH*rinv3*vel0 - 3.0*alpha*acc1;
            double phi1 = -mBH * rinv;

            snp += mBH*rinv3*acc0 - 6.0*alpha*jerk1 - 3.0*beta*acc1;
            acc += acc1;
            jrk += jerk1;
            phi += phi1;

            if(dist2<nnb_r2){
                nnb_id = iBH;
                nnb_r2 = dist2;
                //cout << index << endl;
            }
    }
    
};  