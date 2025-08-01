void evolve_split_pass(int clevel,struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt,
                         int calc_timestep);
void evolve_split_naive(int clevel,struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt,
                         int calc_timestep);
void evolve_split_bridge(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_split_bridge_dkd(int clevel,struct sys s, DOUBLE stime, DOUBLE etime,
                         DOUBLE dt, int calc_timestep);
void evolve_split_hold(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_split_hold_dkd(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_split_pass_dkd(int clevel,struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt,
                         int calc_timestep);
void evolve_split_ppass_dkd(int clevel,struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, 
                         int calc_timestep);
void evolve_sf_4m4(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_sf_4m5(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);

                         
