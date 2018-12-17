#pragma once

//-------------------------------------------------------------------------------
//    Macro defs. for Galaxy simulation test
//-------------------------------------------------------------------------------
#if (INITIAL_CONDITION == 0)
#define ENABLE_GRAVITY_INTERACT
#define ENABLE_HYDRO_INTERACT
#define SET_CIRCULAR_VELOCITY
//------------------------------------------
//    Macro defs. for cold collapse test
//------------------------------------------
#elif (INITIAL_CONDITION == 1)
#define ENABLE_GRAVITY_INTERACT
//-----------------------------------
//    Macro defs. for Evrard test
//-----------------------------------
#elif (INITIAL_CONDITION == 2)
#define ENABLE_GRAVITY_INTERACT
#define ENABLE_HYDRO_INTERACT
#if defined(ISOTHERMAL_EOS)
#undef ISOTHERMAL_EOS
#endif
//-------------------------------------------------------------------------------
//    Macro defs. for the mode to make a glass-like distribution of particles
//-------------------------------------------------------------------------------
#elif (INITIAL_CONDITION == 3)
#define ENABLE_HYDRO_INTERACT
#if !defined(ISOTHERMAL_EOS)
#define ISOTHERMAL_EOS
#endif
#define DUMP_VELOCITY_OF_SPH_PARTICLE
#define CHECK_DENSITY_FLUCTUATION
//----------------------
//    Error handling
//----------------------
#else
#error An invalid value is specified for the macro INITIAL_CONDITION.
#endif
