#ifndef pumi_mesh1D_h
#define pumi_mesh1D_h

typedef enum pumi_meshflag{
  unassigned = 0x00,
  uniform = 0x01,
  leftBL = 0x02,
  rightBL = 0x04,
} pumi_meshflag_t;

typedef struct pumi_submesh1D{
  double x_left; //coordinate of left end of submesh
  double x_right; //coordinate of right end of submesh

  int uniform_Nel; // number of elements in the uniform mesh
  double uniform_x_left; // (dependent variable) coordinate of left end of uniform mesh region inside submesh
  double uniform_x_right; // (dependent variable)coordinate of right end of uniform mesh reion inside submesh
  double uniform_t0; // (dependent variable) cell sizes in uniform mesh region inside submesh

  double left_T; // left BL thickness
  double left_r; // growth ratio for left BL mesh
  int left_Nel; // number of elements in the graded mesh from left
  double lBL_x_right; // (dependent variable) coordinate of right end of left graded mesh region inside submesh
  double lBL_t0; // (dependent variable) size of first (leftmost) cell in the left graded mesh region inside submesh

  double right_T; // right BL thickness
  double right_r; //  growth ratio for right BL mesh
  int right_Nel; // number of elements in the graded mesh from right
  double rBL_x_left; // (dependent variable) coordinate of left end of right graded mesh region inside submesh
  double rBL_t0; // (dependent variable) size of first (rightmost) cell in the right graded mesh region inside submesh

  int submesh_total_Nel; // (dependent variable) total number of elements/cells in the submesh

  pumi_meshflag_t pumi_flag; // flag for mesh type (uniform mesh, right BL graded mesh or left BL graded mesh)
} pumi_submesh1D_t;

typedef struct pumi_mesh{
  int nsubmeshes; // number of submeshes
  int ndim; // number of dimensions
  void *pumi_submeshes;
} pumi_mesh_t;

#include "pumi_initiate.h"
#include "pumi_routines.h"

#endif /* pumi_mesh1D_h */
