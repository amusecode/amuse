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
-ported Keigo Nitadori "eight-16th" to AMUSE by Thomas Schano.
-added support for QD in addition to DD by Thomas Schano.
        QD is available at http://crd-legacy.lbl.gov/~dhbailey/mpdist/ (https://github.com/GFTwrt/QD)
        eight-16th is available at https://github.com/nitadori/HighOrder
-added support for additional step acceleration by Thomas Schano
*/

#include <vector>
#include <limits>
#include <ostream>


#include <cstdio>
#include "vector3.h"

#include <qd/qd_real.h>
#include <qd/qd_inline.h>

typedef qd_real m_real;

std::string arg0;
std::string arg1;
std::string arg2;
std::string arg3;
std::string arg4;
std::string arg5;
std::string arg6;
std::string arg7;
std::string arg8;
std::string arg9;
std::string result;

#ifdef _QD_QD_REAL_H
inline m_real rsqrt(const m_real &x)
{
    double y_app = 1.0 / sqrt(to_double(x));
    m_real x2 = mul_pwr2(x, 0.5);
    return y_app * (m_real(1.5) - x2 * dd_real::sqr(y_app));
}
#else
inline m_real rsqrt(const m_real &x)
{
    double y_app = 1.0 / sqrt(to_double(x));
    m_real x2 = mul_pwr2(x, 0.5);
    return y_app * (m_real(1.5) - x2 * sqr(y_app));
}
#endif // _QD_QD_REAL_H*/

typedef vector3<m_real> qvec3;

template <> template <>
vector3<m_real>::operator dvec3() const
{
    return dvec3(to_double(x), to_double(y), to_double(z));
}


#include "hermite16.h"
#include "test8.h"

#define calc_force_on_i calc_force_on_i_p8

#include "nbodysystem.h"

NbodySystem sys;
long last_id_of_parts;
bool new_data;
m_real i_dt;
int numBits;


/*
 * Interface code
 */



int initialize_code()
{
    std::cout <<"initialize_code"<<std::endl;
    last_id_of_parts = sys.nbody =0;
    sys.tsys = 0.0;
    i_dt= 1.0/1024.0*8;
    sys.ptcl .clear();
    sys.pred .clear();
    sys.force.clear();
    std::numeric_limits<m_real> a;
    numBits=a.digits;
    return 0;
}

int get_word_length(int *mynumBits)
{
    std::cout <<"get_word_length"<<std::endl;
    *mynumBits = numBits;
    return 0;
}

int set_word_length(int mynumBits)
{
    std::cout <<"set_word_length"<<std::endl;
    return 0;
}

int map_id(long id)
{
    for (size_t i=0; i<sys.ptcl.size(); i++)
        if (sys.ptcl[i].id == id) return i;
    return -1;
}

int cleanup_code()
{
    std::cout <<"cleanup_code"<<std::endl;
    last_id_of_parts = sys.nbody =0;
    sys.ptcl .clear();
    sys.pred .clear();
    sys.force.clear();
    return 0;
}

int commit_parameters()
{
    std::cout <<"commit_parameters"<<std::endl;
    return 0;
}

int recommit_parameters()
{
    commit_parameters();
    std::cout <<"recommit_parameters"<<std::endl;
    return 0;
}

int commit_particles()
{
    std::cout <<"commit_particles"<<std::endl;
    sys.pred .resize(sys.ptcl.size());
    sys.force.resize(sys.ptcl.size());
    sys.nbody=sys.ptcl.size();
    new_data=true;
    sys.init_force();
    return 0;
}

int recommit_particles()
{
    std::cout <<"recommit_particles"<<std::endl;
    return commit_particles();
}

int synchronize_model()
{
    std::cout <<"synchronize_model"<<std::endl;
    return 0;
}

int new_particle_float64(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius)
{
    std::cout <<"new_particle_float64"<<std::endl;
    last_id_of_parts++;
    *id=last_id_of_parts;
    Particle part;
    part.id= *id;
    part.mass=mass;
    part.radius=radius;
    part.coord[0]= qvec3( x, y, z);
    part.coord[1]= qvec3( vx, vy, vz);
    part.add_acc= qvec3( 0, 0, 0);
    sys.ptcl.push_back(part);
    new_data=true;
    return 0;
}

int new_particle_string(int *id, char *mass, char *x, char *y, char *z, char *vx, char *vy, char *vz, char *radius)
{
    std::cout <<"new_particle"<<std::endl;
    last_id_of_parts++;
    *id=last_id_of_parts;
    Particle part;
    part.id= *id;
    part.mass=*mass;
    part.radius=*radius;
    part.coord[0]= qvec3( *x, *y, *z);
    part.coord[1]= qvec3( *vx, *vy, *vz);
    part.add_acc= qvec3( 0, 0, 0);
    sys.ptcl.push_back(part);
    new_data=true;
    return 0;
}

int delete_particle (int id)
{
    std::cout <<"delete_particle"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl.erase(sys.ptcl.begin()+part);
        new_data=true;
        return 0;
    }
    return part;
}

int get_number_of_particles(int* N)
{
    std::cout <<"get_number_of_particles"<<std::endl;
    *N=sys.ptcl.size();
    return 0;
}

int get_index_of_first_particle(int* id)
{
    std::cout <<"get_index_of_first_particle"<<std::endl;
    *id=sys.ptcl[0].id;
    return 0;
}

int get_index_of_next_particle(int ida, int* idn)
{
    std::cout <<"get_index_of_next_particle"<<std::endl;
    return -2;
}

int npec =1;
size_t i=0;
m_real e0;

int evolve_model(double t_end)
{
    std::cout <<"evolve_model"<<std::endl;
    m_real ldt = t_end-sys.tsys;
    if(ldt > i_dt)
    {
        ldt=i_dt;
    }
    sys.set_fixed_dt(ldt);
    if(new_data==true)
    {
//        sys.init_force();
//        i=0;
        new_data=false;
        e0 = sys.calc_energy_from_ptcl();
    }

    sys.init_force();
    i=0;
    m_real e1;
    fprintf(stderr, "e0 : %24.63f\n", to_double(e0));
    bool run=true;
    do
    {
        const int n = i ? npec : npec+4;
        sys.tsys += ldt;
        if (sys.tsys >t_end)
        {
            run=false;
            sys.tsys=t_end;
            sys.set_fixed_dt(sys.tsys-sys.ptcl[0].tlast);
        }
        sys.predict_all();

        // P(EC)^n iteration
        for(int nn=0; nn<n; nn++)
        {
            sys.calc_force_on_first_nact(sys.nbody );
            sys.correct_and_feedback();
        }
        sys.calc_force_on_first_nact(sys.nbody );
        sys.correct_and_commit();
        e1 = sys.calc_energy_from_ptcl();


        i++;
    }
    while(run);
    fprintf(stderr, "e1 : %24.63f\n", to_double(e1));
    fprintf(stderr, "de : %24.63f\n", to_double(e0-e1));
    {
        int N = sys.ptcl.size();
        for(int i=0; i<N; i++)
        {
            sys.ptcl[i].add_acc = qvec3( 0, 0, 0);
        }
    }

    return 0;
}

int get_time(double* time)
{
    std::cout <<"get_time"<<std::endl;
    *time = to_double(sys.tsys);
    return 0;
}

int get_mass(int id, double* mass)
{
    std::cout <<"get_mass"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        *mass=to_double(sys.ptcl[part].mass);
        return 0;
    }
    return part;
}

int set_mass(int id, double mass)
{
    std::cout <<"set_mass"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].mass=mass;
        new_data=true;
        return 0;
    }
    return part;
}

int get_radius(int id, double* radius)
{
    std::cout <<"get_radius"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        *radius=to_double(sys.ptcl[part].radius);
        return 0;
    }
    return part;
}

int set_radius(int id, double radius)
{
    std::cout <<"set_radius"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].radius=radius;
        new_data=true;
        return 0;
    }
    return part;
}

int get_position(int id, double* x, double* y, double* z)
{
    std::cout <<"get_position"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        *x=to_double(sys.ptcl[part].coord[0].x);
        *y=to_double(sys.ptcl[part].coord[0].y);
        *z=to_double(sys.ptcl[part].coord[0].z);
        return 0;
    }
    return part;
}

int set_position(int id, double x, double y, double z)
{
    std::cout <<"set_position"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].coord[0].x=x;
        sys.ptcl[part].coord[0].y=y;
        sys.ptcl[part].coord[0].z=z;
        return 0;
    }
    return part;
}

int get_velocity(int id, double* vx, double* vy, double* vz)
{
    std::cout <<"get_velocity"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        *vx=to_double(sys.ptcl[part].coord[1].x);
        *vy=to_double(sys.ptcl[part].coord[1].y);
        *vz=to_double(sys.ptcl[part].coord[1].z);
        return 0;
    }
    return part;
}

int set_velocity(int id, double vx, double vy, double vz)
{
    std::cout <<"set_velocity"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].coord[1].x=vx;
        sys.ptcl[part].coord[1].y=vy;
        sys.ptcl[part].coord[1].z=vz;
        return 0;
    }
    return part;
}

int get_state(int id, double* m, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* radius)
{
    std::cout <<"get_state"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        *m=to_double(sys.ptcl[part].mass);
        *x=to_double(sys.ptcl[part].coord[0].x);
        *y=to_double(sys.ptcl[part].coord[0].y);
        *z=to_double(sys.ptcl[part].coord[0].z);
        *vx=to_double(sys.ptcl[part].coord[1].x);
        *vy=to_double(sys.ptcl[part].coord[1].y);
        *vz=to_double(sys.ptcl[part].coord[1].z);
        *radius=to_double(sys.ptcl[part].radius);
        return 0;
    }
    return part;
}

int set_state(int id, double m, double x, double y, double z, double vx, double vy, double vz, double radius)
{
    std::cout <<"set_state"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].mass=m;
        sys.ptcl[part].coord[0].x=x;
        sys.ptcl[part].coord[0].y=y;
        sys.ptcl[part].coord[0].z=z;
        sys.ptcl[part].coord[1].x=vx;
        sys.ptcl[part].coord[1].y=vy;
        sys.ptcl[part].coord[1].z=vz;
        sys.ptcl[part].radius=radius;
        return 0;
    }
    return part;
}

int get_eps2(double *eps2)
{
    std::cout <<"get_eps2"<<std::endl;
    return -2;
}

int set_eps2(double eps2)
{
    std::cout <<"set_eps2"<<std::endl;
    return -2;
}

int get_time_step(double* dtt)
{
    std::cout <<"get_time_step"<< i_dt.to_string() << std::endl;
    *dtt = to_double(i_dt);
    return 0;
}

int get_const_time_step(double* dtt)
{
    std::cout <<"get_const_time_step"<< i_dt.to_string() << std::endl;
    *dtt = to_double(i_dt);
    return 0;
}

int set_const_time_step(double dtt)
{
    i_dt=dtt;
    std::cout <<"set_const_time_step"<< i_dt.to_string() << std::endl;
    return 0;
}

int get_begin_time(double * output)
{
    std::cout <<"get_begin_time"<<std::endl;
    *output = to_double(sys.ptcl[0].tlast);
    return 0;
}

int set_begin_time(double input)
{
    std::cout <<"set_begin_time"<< input <<std::endl;
    new_data=true;
    sys.tsys = input;
    return 0;
}

int get_total_mass(double* M)
{
    std::cout <<"get_total_mass" <<std::endl;
    int N = sys.ptcl.size();
    m_real Mtot =0;
    for(int i=0; i<N; i++)
    {
        Mtot += sys.ptcl[i].mass;
    }
    *M = to_double(Mtot);
    return 0;
}

int get_center_of_mass_position(double* x, double* y, double* z)
{
    std::cout <<"get_center_of_mass_position" <<std::endl;
    return -2;
}

int get_center_of_mass_velocity(double* vx, double* vy, double* vz)
{
    std::cout <<"get_center_of_mass_velocity" <<std::endl;
    return -2;
}

int get_total_radius(double* R)
{
    std::cout <<"get_total_radius" <<std::endl;
    return -2;
}

int get_acceleration(int id, double* ax, double* ay, double* az)
{
    std::cout <<"get_acceleration" <<std::endl;
    return -2;
}

int set_acceleration(int id, double ax, double ay, double az)
{
    std::cout <<"set_acceleration" <<std::endl;
    return -2;
}

int get_potential(int id, double* pot)
{
    std::cout <<"get_potential" <<std::endl;
    return -2;
}

m_real get_kinetic_energy_m()
{
    int N = sys.ptcl.size();
    m_real ektot = 0;
    m_real m, vx, vy, vz, v2;
    for(int i=0; i<N; i++)
    {
        m_real m  = sys.ptcl[i].mass;
        m_real vx = sys.ptcl[i].coord[1].x;
        m_real vy = sys.ptcl[i].coord[1].y;
        m_real vz = sys.ptcl[i].coord[1].z;
        m_real v2 = vx*vx + vy*vy + vz*vz;
        ektot +=  m * v2 ;
    }
    ektot /=2;
    std::cout <<"ektot " <<ektot.to_string()<<std::endl;
    return ektot;
}

int get_kinetic_energy_string( char **ep)
{
    arg9=get_kinetic_energy_m().to_string();
    *ep =(char*) arg9.c_str();
    return 0;
}

int get_kinetic_energy(double* ek)
{
    *ek = to_double(get_kinetic_energy_m());
    return 0;
}

m_real get_potential_energy_m()
{
    int N = sys.ptcl.size();
    m_real eptot = 0;
    m_real mi,xi, yi, zi, mj, xj, yj, zj, dx, dy, dz;
    m_real dr2;
    for(int i=0; i<N-1; i++)
    {
        mi = sys.ptcl[i].mass;
        xi = sys.ptcl[i].coord[0].x;
        yi = sys.ptcl[i].coord[0].y;
        zi = sys.ptcl[i].coord[0].z;
        for(int j=i+1; j<N; j++)
        {
            mj = sys.ptcl[j].mass;
            xj = sys.ptcl[j].coord[0].x;
            yj = sys.ptcl[j].coord[0].y;
            zj = sys.ptcl[j].coord[0].z;
            dx = xj - xi;
            dy = yj - yi;
            dz = zj - zi;
            dr2 = dx*dx + dy*dy + dz*dz;
            eptot -= mi*mj/sqrt(dr2);
        }
    }
    std::cout <<"eptot " <<eptot.to_string()<<std::endl;
    return eptot;
}

int get_potential_energy_string( char **ep)
{
    arg9=get_potential_energy_m().to_string();
    *ep =(char*) arg9.c_str();
    return 0;
}

int get_potential_energy(double* ep)
{
    *ep = to_double(get_potential_energy_m());
    return 0;
}

int get_total_energy_string( char **ep)
{
    m_real etot = get_kinetic_energy_m() + get_potential_energy_m();
    arg9=etot.to_string();
    std::cout <<"etot " <<arg9<<std::endl;
    std::cout <<"etot_na " <<sys.calc_energy_from_ptcl().to_string()<<std::endl;
    *ep =(char*) arg9.c_str();
    return 0;
}

int get_mass_string(int id, char **mass)
{
    std::cout <<"get_mass_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        arg0 = sys.ptcl[part].mass.to_string();
        *mass = (char*) arg0.c_str();
        return 0;
    }
    return part;
}

int set_mass_string(int id, char *mass)
{
    std::cout <<"set_mass_string"<<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].mass=*mass;
        return 0;
    }
    return part;
}

int get_position_string(int id, char **x, char **y, char **z)
{
    std::cout <<"get_position_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        arg1 = sys.ptcl[part].coord[0].x.to_string();
        arg2 = sys.ptcl[part].coord[0].y.to_string();
        arg3 = sys.ptcl[part].coord[0].z.to_string();
        *x = (char*) arg1.c_str();
        *y = (char*) arg2.c_str();
        *z = (char*) arg3.c_str();
        return 0;
    }
    return part;
}

int set_position_string(int id, char *x, char *y, char *z)
{
    std::cout <<"set_position_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].coord[0].x=*x;
        sys.ptcl[part].coord[0].y=*y;
        sys.ptcl[part].coord[0].z=*z;
        return 0;
    }
    return part;
}
int get_velocity_string(int id, char **vx, char **vy, char **vz)
{
    std::cout <<"get_velocity_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        arg4 = sys.ptcl[part].coord[1].x.to_string();
        arg5 = sys.ptcl[part].coord[1].y.to_string();
        arg6 = sys.ptcl[part].coord[1].z.to_string();
        *vx = (char*) arg4.c_str();
        *vy = (char*) arg5.c_str();
        *vz = (char*) arg6.c_str();
        return 0;
    }
    return part;
}

int set_velocity_string(int id, char* vx, char* vy, char* vz)
{
    std::cout <<"set_velocity_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].coord[1].x=*vx;
        sys.ptcl[part].coord[1].y=*vy;
        sys.ptcl[part].coord[1].z=*vz;
        return 0;
    }
    return part;
}

int get_state_string(int id, char** m, char** x, char** y, char** z, char** vx, char** vy, char** vz, char** radius)
{
    std::cout <<"get_state_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        arg0 = sys.ptcl[part].mass.to_string();
        arg1 = sys.ptcl[part].coord[0].x.to_string();
        arg2 = sys.ptcl[part].coord[0].y.to_string();
        arg3 = sys.ptcl[part].coord[0].z.to_string();
        arg4 = sys.ptcl[part].coord[1].x.to_string();
        arg5 = sys.ptcl[part].coord[1].y.to_string();
        arg6 = sys.ptcl[part].coord[1].z.to_string();
        arg7 = sys.ptcl[part].radius.to_string();
        *m = (char*) arg0.c_str();
        *x = (char*) arg1.c_str();
        *y = (char*) arg2.c_str();
        *z = (char*) arg3.c_str();
        *vx = (char*) arg4.c_str();
        *vy = (char*) arg5.c_str();
        *vz = (char*) arg6.c_str();
        *radius = (char*) arg7.c_str();
        return 0;
    }
    return part;
}

int set_state_string(int id, char* m, char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius)
{
    std::cout <<"set_state_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].mass=*m;
        sys.ptcl[part].coord[0].x=*x;
        sys.ptcl[part].coord[0].y=*y;
        sys.ptcl[part].coord[0].z=*z;
        sys.ptcl[part].coord[1].x=*vx;
        sys.ptcl[part].coord[1].y=*vy;
        sys.ptcl[part].coord[1].z=*vz;
        sys.ptcl[part].radius=*radius;
        return 0;
    }
    return part;
}

int get_radius_string(int id, char** radius)
{
    std::cout <<"get_radius_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        arg7 = sys.ptcl[part].radius.to_string();
        *radius = (char*) arg7.c_str();
        return 0;
    }
    return part;
}

int set_radius_string(int id, char* radius)
{
    std::cout <<"set_radius_string" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].radius=*radius;
        return 0;
    }
    return part;
}

int get_total_mass_string(char **M)
{
    std::cout <<"get_total_mass_string" <<std::endl;
    int N = sys.ptcl.size();
    m_real Mtot =0;
    for(int i=0; i<N; i++)
    {
        Mtot += sys.ptcl[i].mass;
    }
    arg8=Mtot.to_string();
    *M = (char*) arg8.c_str();
    return 0;
}

int add_step_acceleration_float64(int id, double a_step_x, double a_step_y, double a_step_z)
{
    std::cout <<"add_step_acceleration_float64" <<std::endl;
    int part=map_id(id);
    if (part >=0)
    {
        sys.ptcl[part].add_acc.x=a_step_x;
        sys.ptcl[part].add_acc.y=a_step_y;
        sys.ptcl[part].add_acc.z=a_step_z;
        return 0;
    }
    return part;
}


