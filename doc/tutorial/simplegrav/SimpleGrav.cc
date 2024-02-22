#include <math.h>

// Compute a single Eulerian gravity step
int gravity_step (double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int N, double dt, double eps_sq, double G) {

    double *new_x = new double[N];
    double *new_y = new double[N];
    double *new_z = new double[N];

    double *new_vx = new double[N];
    double *new_vy = new double[N];
    double *new_vz = new double[N];

    double ax, ay, az;
    double dR, dx, dy, dz;

    for (int i = 0; i < N; i++) {
        ax = 0.;
        ay = 0.;
        az = 0.;

        for (int j = 0; j < N; j++) {
            if (i != j) {
                dx = x[j] - x[i];
                dy = y[j] - y[i];
                dz = z[j] - z[i];

                dR = sqrt( dx*dx + dy*dy + dz*dz + eps_sq );

                ax += G*mass[j]*dx/pow(dR, 3);
                ay += G*mass[j]*dy/pow(dR, 3);
                az += G*mass[j]*dz/pow(dR, 3);
            }
        }

        new_x[i] = x[i] + vx[i]*dt;
        new_y[i] = y[i] + vy[i]*dt;
        new_z[i] = z[i] + vz[i]*dt;

        new_vx[i] = vx[i] + ax*dt;
        new_vy[i] = vy[i] + ay*dt;
        new_vz[i] = vz[i] + az*dt;
    }

    for (int i = 0; i < N; i++) {
        x[i] = new_x[i];
        y[i] = new_y[i];
        z[i] = new_z[i];

        vx[i] = new_vx[i];
        vy[i] = new_vy[i];
        vz[i] = new_vz[i];
    }

    delete new_x;
    delete new_y;
    delete new_z;
    delete new_vx;
    delete new_vy;
    delete new_vz;

    return 0;
}
