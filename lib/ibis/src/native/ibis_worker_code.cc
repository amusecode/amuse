#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "bhtree_code.h"
#include "worker_code.h"
#include "stopcond.h"

class message_header {
public:
  int tag;
  int len;
  int number_of_doubles;
  int number_of_ints;
  int number_of_floats;
  int number_of_strings;
  int number_of_bools;
  int number_of_longs;

  message_header(): tag(0), len(1), number_of_doubles(0), number_of_ints(0), number_of_floats(0), number_of_strings(0), number_of_bools(0), number_of_longs(0){}

};




  bool must_run_loop = true;
  char * characters = 0, * output_characters = 0;
  
  int max_len = 10;
  int * ints_in = new int[ max_len * 2];
  int * ints_out = new int[ max_len * 3];
  
  int * strings_in = new int[ max_len * 2];
  
  double * doubles_in = new double[ max_len * 8];
  double * doubles_out = new double[ max_len * 8];

  
  
    
    message_header request_header;
    message_header reply_header;
    
    request_header.recv(parent,rank);
    if (request_header.len > max_len) {
      max_len = request_header.len + 255;
      delete[] ints_in;
      delete[] ints_out;
      delete[] strings_in;
      delete[] doubles_in;
      delete[] doubles_out;
      ints_in = new int[ max_len * 2];
      ints_out = new int[ max_len * 3];
      
      strings_in = new int[ max_len * 2];
      
      doubles_in = new double[ max_len * 8];
      doubles_out = new double[ max_len * 8];
      
      
      
      
    }
    if(request_header.number_of_doubles > 0) {
      parent.Bcast(doubles_in, request_header.number_of_doubles * request_header.len, MPI_DOUBLE, 0);
    }
    if(request_header.number_of_ints > 0) {
      parent.Bcast(ints_in, request_header.number_of_ints * request_header.len, MPI_INT, 0);
    }
    
    if(request_header.number_of_strings > 0) {
      parent.Bcast(strings_in, request_header.number_of_strings * request_header.len, MPI_INTEGER, 0);
      characters = new char[strings_in[request_header.number_of_strings * request_header.len - 1] + 1];
      parent.Bcast(characters,  strings_in[request_header.number_of_strings * request_header.len- 1] + 1, MPI_CHARACTER, 0);
    }
    
    
    
    reply_header.tag = request_header.tag;
    
    reply_header.len = request_header.len;
    
    switch(request_header.tag) {
      case 0:
        must_run_loop = false;
        break;
      case 20284990:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_mass(
            ints_in[i] ,
            &doubles_out[i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 20920053:
        ints_out[0] = commit_particles();
        reply_header.number_of_ints = 1;
        break;
      
      case 37921492:
        ints_out[0] = set_stopping_condition_timeout_parameter(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 44188957:
        ints_out[0] = get_time(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 55746553:
        ints_out[0] = get_theta_for_tree(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 104547857:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = set_mass(
            ints_in[i] ,
            doubles_in[i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 111364973:
        ints_out[0] = evolve(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 128926247:
        ints_out[0] = get_index_of_first_particle(
          &ints_out[1]
        );
        reply_header.number_of_ints = 2;
        break;
      
      case 154188853:
        ints_out[0] = get_dt_dia(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 159095171:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = enable_stopping_condition(
            ints_in[i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 205426934:
        ints_out[0] = get_total_radius(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 210995141:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_potential_at_point(
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i] ,
            doubles_in[( 3 * request_header.len) + i] ,
            &doubles_out[i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 219539042:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_number_of_stopping_conditions_set(
            &ints_out[( 1 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 2;
        break;
      
      case 271440754:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = is_stopping_condition_set(
            ints_in[i] ,
            &ints_out[( 1 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 2;
        break;
      
      case 290264013:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = new_particle(
            &ints_out[( 1 * request_header.len) + i] ,
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i] ,
            doubles_in[( 3 * request_header.len) + i] ,
            doubles_in[( 4 * request_header.len) + i] ,
            doubles_in[( 5 * request_header.len) + i] ,
            doubles_in[( 6 * request_header.len) + i] ,
            doubles_in[( 7 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 2;
        break;
      
      case 384567015:
        ints_out[0] = get_total_mass(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 508605261:
        ints_out[0] = set_stopping_condition_out_of_box_parameter(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 509916277:
        ints_out[0] = reinitialize_particles();
        reply_header.number_of_ints = 1;
        break;
      
      case 542058817:
        ints_out[0] = set_eps2(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 632979349:
        ints_out[0] = set_stopping_condition_number_of_steps_parameter(
          ints_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 639951605:
        ints_out[0] = get_stopping_condition_timeout_parameter(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 643372263:
        ints_out[0] = set_theta_for_tree(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 658631024:
        ints_out[0] = get_eps2(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 662019495:
        ints_out[0] = set_dt_dia(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 678380482:
        ints_out[0] = get_index_of_next_particle(
          ints_in[0] ,
          &ints_out[1]
        );
        reply_header.number_of_ints = 2;
        break;
      
      case 728786188:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = delete_particle(
            ints_in[i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 733749514:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = is_stopping_condition_enabled(
            ints_in[i] ,
            &ints_out[( 1 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 2;
        break;
      
      case 835969050:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_potential(
            ints_in[i] ,
            &doubles_out[i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 887125873:
        ints_out[0] = synchronize_model();
        reply_header.number_of_ints = 1;
        break;
      
      case 929181341:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = set_state(
            ints_in[i] ,
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i] ,
            doubles_in[( 3 * request_header.len) + i] ,
            doubles_in[( 4 * request_header.len) + i] ,
            doubles_in[( 5 * request_header.len) + i] ,
            doubles_in[( 6 * request_header.len) + i] ,
            doubles_in[( 7 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 967950880:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_state(
            ints_in[i] ,
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i] ,
            &doubles_out[( 3 * request_header.len) + i] ,
            &doubles_out[( 4 * request_header.len) + i] ,
            &doubles_out[( 5 * request_header.len) + i] ,
            &doubles_out[( 6 * request_header.len) + i] ,
            &doubles_out[( 7 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 8;
        break;
      
      case 1024680297:
        ints_out[0] = get_time_step(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 1039508780:
        ints_out[0] = set_use_self_gravity(
          ints_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 1050085724:
        ints_out[0] = recommit_particles();
        reply_header.number_of_ints = 1;
        break;
      
      case 1071152125:
        ints_out[0] = get_kinetic_energy(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 1082457792:
        ints_out[0] = get_number_of_particles(
          &ints_out[1]
        );
        reply_header.number_of_ints = 2;
        break;
      
      case 1098838617:
        ints_out[0] = get_stopping_condition_number_of_steps_parameter(
          &ints_out[1]
        );
        reply_header.number_of_ints = 2;
        break;
      
      case 1137702459:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = disable_stopping_condition(
            ints_in[i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 1141573512:
        ints_out[0] = internal__redirect_outputs(
          characters + ( 0- 1 < 0 ? 0 :strings_in[0 - 1] + 1) ,
          characters + ( 1- 1 < 0 ? 0 :strings_in[1 - 1] + 1)
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 1214559355:
        ints_out[0] = get_epsilon_squared(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 1221190952:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = set_acceleration(
            ints_in[i] ,
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 1221524787:
        ints_out[0] = get_indices_of_colliding_particles(
          &ints_out[1] ,
          &ints_out[2]
        );
        reply_header.number_of_ints = 3;
        break;
      
      case 1231060790:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_center_of_mass_position(
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 3;
        break;
      
      case 1236118821:
        ints_out[0] = set_time_step(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 1303477613:
        ints_out[0] = set_epsilon_squared(
          doubles_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 1315680918:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_center_of_mass_velocity(
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 3;
        break;
      
      case 1317242279:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_radius(
            ints_in[i] ,
            &doubles_out[i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 1452242270:
        ints_out[0] = set_ncrit_for_tree(
          ints_in[0]
        );
        reply_header.number_of_ints = 1;
        break;
      
      case 1623630901:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = set_radius(
            ints_in[i] ,
            doubles_in[i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 1625942852:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = has_stopping_condition(
            ints_in[i] ,
            &ints_out[( 1 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 2;
        break;
      
      case 1644113439:
        ints_out[0] = cleanup_code();
        reply_header.number_of_ints = 1;
        break;
      
      case 1744145122:
        ints_out[0] = recommit_parameters();
        reply_header.number_of_ints = 1;
        break;
      
      case 1768994498:
        ints_out[0] = initialize_code();
        reply_header.number_of_ints = 1;
        break;
      
      case 1814108848:
        ints_out[0] = get_use_self_gravity(
          &ints_out[1]
        );
        reply_header.number_of_ints = 2;
        break;
      
      case 1852958273:
        ints_out[0] = get_potential_energy(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 1877555269:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_gravity_at_point(
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i] ,
            doubles_in[( 3 * request_header.len) + i] ,
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 3;
        break;
      
      case 1892689129:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_velocity(
            ints_in[i] ,
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 3;
        break;
      
      case 1937183958:
        ints_out[0] = get_stopping_condition_out_of_box_parameter(
          &doubles_out[0]
        );
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 1;
        break;
      
      case 1938095680:
        ints_out[0] = get_ncrit_for_tree(
          &ints_out[1]
        );
        reply_header.number_of_ints = 2;
        break;
      
      case 2010900811:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_position(
            ints_in[i] ,
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 3;
        break;
      
      case 2026192840:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = set_position(
            ints_in[i] ,
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      case 2046699524:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_stopping_condition_info(
            ints_in[i] ,
            &ints_out[( 1 * request_header.len) + i] ,
            &ints_out[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 3;
        break;
      
      case 2061473599:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_acceleration(
            ints_in[i] ,
            &doubles_out[i] ,
            &doubles_out[( 1 * request_header.len) + i] ,
            &doubles_out[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        reply_header.number_of_doubles = 3;
        break;
      
      case 2069478464:
        ints_out[0] = commit_parameters();
        reply_header.number_of_ints = 1;
        break;
      
      case 2129795713:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = get_stopping_condition_particle_index(
            ints_in[i] ,
            ints_in[( 1 * request_header.len) + i] ,
            &ints_out[( 1 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 2;
        break;
      
      case 2144268908:
        for (int i = 0 ; i < request_header.len; i++){
          ints_out[i] = set_velocity(
            ints_in[i] ,
            doubles_in[i] ,
            doubles_in[( 1 * request_header.len) + i] ,
            doubles_in[( 2 * request_header.len) + i]
          );
        }
        reply_header.number_of_ints = 1;
        break;
      
      default:
        reply_header.tag = -1;
    }
    
    MPI::COMM_WORLD.Barrier();
    
    
    if(rank == 0) {
      
      
      reply_header.send(parent, rank);
      
      
      if(reply_header.number_of_doubles > 0) {
        parent.Send(doubles_out, reply_header.number_of_doubles * request_header.len, MPI_DOUBLE, 0, 999);
      }
      if(reply_header.number_of_ints > 0) {
        parent.Send(ints_out, reply_header.number_of_ints * request_header.len, MPI_INT, 0, 999);
      }
    
    }
    if (characters) { delete[] characters; characters = 0;}
    if (output_characters) { delete[] output_characters; output_characters = 0;}
  }
  delete[] ints_in;
  delete[] ints_out;
  delete[] strings_in;
  delete[] doubles_in;
  delete[] doubles_out;
  
  parent.Disconnect();
}

int main(int argc, char *argv[])
{
  
  MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
  atexit(onexit);
  
  run_loop();
  
  MPI_Finalize();
  return 0;
}
