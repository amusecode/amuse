void evolve_split_pass(struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt,
                         int calc_timestep);
void evolve_split_naive(struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt,
                         int calc_timestep);
void evolve_split_bridge(struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_split_bridge_dkd(struct sys s, DOUBLE stime, DOUBLE etime,
                         DOUBLE dt, int calc_timestep);
void evolve_split_hold(struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_split_hold_dkd(struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_split_pass_dkd(struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt,
                         int calc_timestep);
void evolve_split_ppass_dkd(struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, 
                         int calc_timestep);
void evolve_sf_4m4(struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);
void evolve_sf_4m5(struct sys s, DOUBLE stime, DOUBLE etime, 
                         DOUBLE dt, int calc_timestep);

                         
