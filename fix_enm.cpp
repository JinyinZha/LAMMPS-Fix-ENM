/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "fix_enm.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "comm.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include "modify.h"
#include "memory.h"
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixENM::FixENM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int i,j,tmp;
  bool not_find;
  FILE *fq;  
  struct Partial_Bond{
  int ca_index;
  int tag_two;
  float sum;
  }*partial_bond;
  int natoms = atom->natoms;
  int me = comm->me;
  char str[500];
   
  if (me == 0) printf("Now fix ENM!\n");  
  if (me == 0 && logfile) fprintf(logfile,"Now fix ENM!\n");
  //1. Check Args
  if ( narg != 6)
    error->all(FLERR,"Invalid ENM command!\nToo few args!\nUsage: fix fix_id group_id enm enm_bond_type_file enm_bond_file");
  
  //2, Read k_styles
  //2.1 Open the file
  if (me ==0) printf("ENM now reading k style file.\n");
  if (me == 0 && logfile) fprintf(logfile,"ENM now reading k style file.\n");
  fq = fopen(arg[3],"r");
  if (fq == NULL)   error->all(FLERR,"Cannot open enm k style file\n");
  //2.2 Read the number of k type
  fscanf(fq,"%d\n",&k_type_num);
  //2.3 Read k type )
  memory->create(k,k_type_num+1,"fix:enm_k");
  //k = new double[k_type_num+1];//first is null
  for(i=1;i<=k_type_num;i++) fscanf(fq,"%d,%lf\n",&tmp,&k[i]);
  //2.4 Close the file 
  fclose(fq);

   //3, Read r0_styles
   //3.1 Open the file
   if (me ==0) printf("ENM now reading r0 style file.\n");
   if (me == 0 && logfile) fprintf(logfile,"ENM now reading r0 style file.\n");
   fq = fopen(arg[4],"r");
   if (fq == NULL) error->all(FLERR,"Cannot open enm r0 style file\n");
   //3.2 Read the number of r0 type
   fscanf(fq,"%d\n",&r0_type_num);
   //3.3 Read r0 type
   memory->create(r0,r0_type_num+1,"fix:enm_r0");
   //r0 = new double[r0_type_num+1];//first is null
   for(i=1;i<=r0_type_num;i++) fscanf(fq,"%d,%lf\n",&tmp,&r0[i]);
   //3.4 Close the file 
   fclose(fq);


  
  //4,Read Bonds
  //4.1 Open the file
  if (me==0) printf("ENM now reading bond file.\n");
  if (me == 0 && logfile) fprintf(logfile,"ENM now reading bond file.\n");
  fq = fopen(arg[5],"r");
     if (fq == NULL)
       error->all(FLERR,"Cannot open enm bond file\n");
  //4.2 Read the number of bonds and get atom_num
  fscanf(fq,"%d\n",&bond_num);  
  //4.3 Read bond data
  bond_list = new Bond_List[bond_num];
  partial_bond = new Partial_Bond[bond_num];
  partial_bond_num = 0;
  for(i=0;i<bond_num;i++){
    //4.3.1 Read
    fscanf(fq,"%d,%d,%d,%d,%lf,%d\n",&bond_list[i].tag1,&bond_list[i].tag2,&bond_list[i].k_type,&bond_list[i].r0_type,&bond_list[i].fragment,&bond_list[i].ca_id);
    //4.3.2 Checking
    if (bond_list[i].k_type <= 0 || bond_list[i].k_type > k_type_num || bond_list[i].r0_type <= 0 || bond_list[i].r0_type > r0_type_num){
      sprintf(str,"Error found in line %d! Invalid bond type!\nk type and r0 type should belongs to [1,type_num]\n",i+1);
      error->all(FLERR,str);
    }
    if (bond_list[i].fragment <= 0.0 || bond_list[i].fragment > 1.0){
      sprintf(str,"Error found in line %d! Invalid Fragment: %lf\nFragement should belongs to (0,1]\n",i+1,bond_list[i].fragment);
      error->all(FLERR,str);
    }
    if (bond_list[i].fragment < 1.0){
      not_find = true;
      for(j=0;j<partial_bond_num;j++){
        if (partial_bond[j].ca_index == bond_list[i].ca_id && partial_bond[j].tag_two == bond_list[i].tag2){
          not_find = false;
          partial_bond[j].sum += bond_list[i].fragment;
          bond_list[i].saved_id = j;
          break;
        }
      }
      if (not_find){
        bond_list[i].saved_id = partial_bond_num;
        partial_bond[partial_bond_num].ca_index = bond_list[i].ca_id;
        partial_bond[partial_bond_num].sum = bond_list[i].fragment;
        partial_bond[partial_bond_num].tag_two = bond_list[i].tag2;
        partial_bond_num++;
      }
      if (bond_list[i].ca_id <= 0 || bond_list[i].ca_id > natoms){
        sprintf(str,"Error found in line %d! Invalid CA id: %d\n CA id should belongs to [1,natoms]\n",i+1,bond_list[i].ca_id);
        error->all(FLERR,str);
      }
    }
    if (bond_list[i].tag1 >= bond_list[i].tag2){
      sprintf(str,"Error found in line %d!\nFirst atom number should be smaller than the second!\n",i+1);
      error->all(FLERR,str);
    }
    if (bond_list[i].tag1 < 0){
      sprintf(str,"Error found in line %d!\nAtom index should be positive!\n",i+1);
      error->all(FLERR,str);
    }
  }
  //4.5 Close the file 
  fclose(fq);
  //5, Make saved partial bond data
  if (partial_bond_num > 0) saved_partial_bond_data = new Saved_Partial_Bond_Data[partial_bond_num];
  //6, Check the validity of partial bonds
  for(i=0;i<partial_bond_num;i++){
    if (partial_bond[i].sum < 0.900 || partial_bond[i].sum > 1.100) {
      sprintf(str,"Error found in partial bonds!(ca_ids:%d and %d) Sum is %f.\nCheck your partial bonds!\nThe sum of fragment belonging to same CA should be 1!\n",partial_bond[i].ca_index,partial_bond[i].tag_two,partial_bond[i].sum);
      error->all(FLERR,str);
    }
  }
if (me==0) printf("Fix ENM Summary:\nENM k Type Number: %d\nENM r0 Type Number: %d\nENM Bond Number: %d\nPartial Bond Number: %d\n",k_type_num,r0_type_num,bond_num,partial_bond_num);
if (me==0 && logfile) fprintf(logfile,"Fix ENM Summary:\nENM k Type Number: %d\nENM r0 Type Number: %d\nENM Bond Number: %d\nPartial Bond Number: %d\n",k_type_num,r0_type_num,bond_num,partial_bond_num);
delete [] partial_bond;
//7, Initialize local_index
memory->create(local_index,natoms+1,"fix:enm_local_index");
}


/* ---------------------------------------------------------------------- */

int FixENM::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
//  mask |= FINAL_INTEGRATE;
//  mask |= INITIAL_INTEGRATE_RESPA;
//  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixENM::~FixENM(){
memory->destroy(k);
memory->destroy(r0);
delete [] bond_list;
delete [] saved_partial_bond_data;
memory->destroy(local_index);


}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */



double FixENM::compute_scalar(){
   double all;
   MPI_Allreduce(&eenm,&all,1,MPI_DOUBLE,MPI_SUM,world);
   
   int me = comm -> me;
   
   return all;
 }






void FixENM::post_force(int vflag)
{
  int i,n,i1,i2,ica,li_in_local_num;
  double delx,dely,delz,ebond,fbond,f_x,f_y,f_z;
  double rsq,r,dr,rk;
  int nlocal = atom->nlocal;
  int nwithin = nlocal + atom->nghost;
  long natoms = atom->natoms;
  
  bool i1_in_local,i2_in_local;
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag; 
  ebond = 0.0;
  eenm = 0.0;
  int me = comm->me;
  double v[6];
  int li[2];
  char str[500];



  //1, Make local_index
  for(i=0;i<=natoms;i++) local_index[i] = -1;
  for(i=0;i<nwithin;i++) local_index[tag[i]] = i;

  //2, ev_setup
  int eflag = 3;//I guess.
  vflag=6;//I guess.
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag=0; 

  //3, Initialize saved partial bond
    for(i=0;i<partial_bond_num;i++)  saved_partial_bond_data[i].calculated = false;

  //4, Looping bondlist. 
  //4.1, Judge. Calculate when at least one atom is in local.   
  for(n=0;n<bond_num;n++){

    i1 = local_index[bond_list[n].tag1];
    i2 = local_index[bond_list[n].tag2];
    if (i1 == -1 or i2 == -1) continue;//two not in local or ghost
    i1_in_local = i1 < nlocal;
    i2_in_local = i2 < nlocal;
    if (i1_in_local || i2_in_local){
      //4.1.1 If local, Judge whether it is a partial bond
      if (bond_list[n].fragment < 1.0){
        //4.1.1.1 If partial, Judge whether this bond is already calculated.
        if (saved_partial_bond_data[bond_list[n].saved_id].calculated){
          //4.1.1.1.1 If calculated, get force, energy and viral
          f_x =  saved_partial_bond_data[bond_list[n].saved_id].fx*bond_list[n].fragment;
          f_y =  saved_partial_bond_data[bond_list[n].saved_id].fy*bond_list[n].fragment;
          f_z =  saved_partial_bond_data[bond_list[n].saved_id].fz*bond_list[n].fragment;
          ebond = saved_partial_bond_data[bond_list[n].saved_id].bond_energy*bond_list[n].fragment;
          v[0] = saved_partial_bond_data[bond_list[n].saved_id].viral[0]*bond_list[n].fragment;
          v[1] = saved_partial_bond_data[bond_list[n].saved_id].viral[1]*bond_list[n].fragment;
          v[2] = saved_partial_bond_data[bond_list[n].saved_id].viral[2]*bond_list[n].fragment;
          v[3] = saved_partial_bond_data[bond_list[n].saved_id].viral[3]*bond_list[n].fragment;
          v[4] = saved_partial_bond_data[bond_list[n].saved_id].viral[4]*bond_list[n].fragment;
          v[5] = saved_partial_bond_data[bond_list[n].saved_id].viral[5]*bond_list[n].fragment;
        }else{
          //4.1.1.2.1 If Not, Find CA first
          ica = local_index[bond_list[n].ca_id];
          //4.1.1.2.2 Then, calculate basic values regardless of fragment
          delx = x[ica][0] - x[i2][0];
          dely = x[ica][1] - x[i2][1];
          delz = x[ica][2] - x[i2][2];
          rsq = delx*delx + dely*dely + delz*delz;
          r = sqrt(rsq);
          dr = r - r0[bond_list[n].r0_type];
          rk = k[bond_list[n].k_type] * dr;
          if (r > 0.0) fbond = -2.0*rk/r;
          else fbond = 0.0;
          ebond = rk*dr; 
          f_x = delx*fbond;
          f_y = dely*fbond;
          f_z = delz*fbond;
          v[0] = delx*delx*fbond;
          v[1] = dely*dely*fbond;
          v[2] = delz*delz*fbond;
          v[3] = delx*dely*fbond;
          v[4] = delx*delz*fbond;
          v[5] = dely*delz*fbond;
          //4.1.1.1.4 Save data
          saved_partial_bond_data[bond_list[n].saved_id].fx = f_x;
          saved_partial_bond_data[bond_list[n].saved_id].fy = f_y;
          saved_partial_bond_data[bond_list[n].saved_id].fz = f_z;
          saved_partial_bond_data[bond_list[n].saved_id].bond_energy = ebond;
          saved_partial_bond_data[bond_list[n].saved_id].viral[0] = v[0];
          saved_partial_bond_data[bond_list[n].saved_id].viral[1] = v[1];
          saved_partial_bond_data[bond_list[n].saved_id].viral[2] = v[2];
          saved_partial_bond_data[bond_list[n].saved_id].viral[3] = v[3];
          saved_partial_bond_data[bond_list[n].saved_id].viral[4] = v[4];
          saved_partial_bond_data[bond_list[n].saved_id].viral[5] = v[5];
          saved_partial_bond_data[bond_list[n].saved_id].calculated = true;
          //4.1.1.2.3 Consider fragment
          f_x *= bond_list[n].fragment;
          f_y *= bond_list[n].fragment;
          f_z *= bond_list[n].fragment;
          ebond *= bond_list[n].fragment;
          v[0] *= bond_list[n].fragment;
          v[1] *= bond_list[n].fragment;
          v[2] *= bond_list[n].fragment;
          v[3] *= bond_list[n].fragment;
          v[4] *= bond_list[n].fragment;
          v[5] *= bond_list[n].fragment;
        }
      }else{
        //4.1.2 If not partial, calculate basic values
        delx = x[i1][0] - x[i2][0];
        dely = x[i1][1] - x[i2][1];
        delz = x[i1][2] - x[i2][2];
        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);
        dr = r - r0[bond_list[n].r0_type];
        rk = k[bond_list[n].k_type] * dr;
        if (r > 0.0) fbond = -2.0*rk/r;
        else fbond = 0.0;
        ebond = rk*dr; 
        f_x = delx*fbond;
        f_y = dely*fbond;
        f_z = delz*fbond;
        v[0] = delx*delx*fbond;
        v[1] = dely*dely*fbond;
        v[2] = delz*delz*fbond;
        v[3] = delx*dely*fbond;
        v[4] = delx*delz*fbond;
        v[5] = dely*delz*fbond;

      }
      //4.2 Apply force and construct li if local
      if (i1_in_local && i2_in_local){
        f[i1][0] += f_x;
        f[i1][1] += f_y;
        f[i1][2] += f_z;
        f[i2][0] -= f_x;
        f[i2][1] -= f_y;
        f[i2][2] -= f_z;
        li[0] = i1;
        li[1] = i2;
        li_in_local_num = 2;
        eenm += ebond;
      }else if (i1_in_local && (not i2_in_local)){
        f[i1][0] += f_x;
        f[i1][1] += f_y;
        f[i1][2] += f_z;
        li[0] = i1;
        li_in_local_num = 1;
        eenm += 0.5 * ebond;
      }else{
        f[i2][0] -= f_x;
        f[i2][1] -= f_y;
        f[i2][2] -= f_z;
        li[0] = i2;
        li_in_local_num = 1;
        eenm += 0.5 * ebond;
      }


      //4.3 ev_tally
      ev_tally(li_in_local_num,li,2.0,ebond,v);
    }
  }
}
   


/* ---------------------------------------------------------------------- */


