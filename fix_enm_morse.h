/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(enm_morse,FixENM_Morse)

#else

#ifndef LMP_FIX_ENM_MORSE_H
#define LMP_FIX_ENM_MORSE_H


#include "fix.h"
#include "atom.h"

namespace LAMMPS_NS {

class FixENM_Morse : public Fix {
 public:
  //bond_type
  double *v0;
  double *alpha;
  double *r0;
  
  
  struct Bond_List{
    int tag1;
    int tag2;//tag1<tag2
    int v0_type;
    int alpha_type;
    int r0_type;
    double fragment;// For partial bonds range = (0,1)
    int ca_id;//only for partial bond
    int saved_id;//only for partial bond
  }*bond_list;

  struct Saved_Partial_Bond_Data{
    double fx;
    double fy;
    double fz;
    double bond_energy;
    double viral[6];
    bool calculated; 
  }*saved_partial_bond_data;
  
  int v0_type_num;
  int alpha_type_num;
  int r0_type_num;
  int bond_num;
  int partial_bond_num;
  
  double eenm;
  FixENM_Morse(class LAMMPS *, int, char **);
  ~FixENM_Morse();
  int *local_index;
  int setmask();
  double compute_scalar();
  //virtual void init();
  virtual void post_force(int);
  //virtual void final_integrate();
  //virtual void initial_integrate_respa(int, int, int);
  //virtual void final_integrate_respa(int, int);
  //virtual void reset_dt();


};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
