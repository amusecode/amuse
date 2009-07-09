
// Accessors used only for communication between the C++ main driver
// and the local dynamics-.cc module.  Probably all the variables
// listed in parameters should be accessible here.

void set_t(double tt);
void set_dt_param(double dt);
void set_dt_dia(double dt);
void set_eps(double e);
