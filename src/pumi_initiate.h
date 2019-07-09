#ifndef pumi_initiate_h
#define pumi_initiate_h

#include "pumiMBBL.h"

#define SUBMESH_FLAGSTRING_LENGTH 23
#define SUBMESH_MAX_SEGMENTS 3
#define SEGMENT_STRING_LENGTH 8

#define DEFAULT_PARAM_VAL 0.0

typedef enum pumi_initiate_flag{
  initiate_from_terminal, //inputs supplied by user in the terminal. The user will be prompted to input mesh parameters
  initiate_from_commandline_inputs, //mesh parameters will be read from the 'command line inputs'
  // of the hpic code and assigned to struct pumi_initiate_input members
} pumi_initiate_flag_t;


typedef struct pumi_initiate_input{
  int ndim; // numnber of dimensions of physical space
  int nsubmeshes; // number of submeshes required
  char **type_flag; //submesh flag as a string
  double *x_left; //coordinate of left end of submesh
  double *x_right; //coordinate of right end of submesh
  int *uniform_Nel; // number of elements in the uniform mesh
  double *left_T; // left BL thickness
  double *left_r; // growth ratio for left BL mesh
  int *left_Nel; // number of elements in the graded mesh from left
  double *right_T; // right BL thickness
  double *right_r; //  growth ratio for right BL mesh
  int *right_Nel; // number of elements in the graded mesh from right
} pumi_initiate_input_t;

pumi_mesh_t* pumi_initiate(pumi_initiate_flag_t pumi_input_initiate_flag, pumi_initiate_input_t *pumi_inputs);
void pumi_initiate_allocate(pumi_initiate_input_t *pumi_inputs, int nsubmeshes);
void pumi_initiate_deallocate(pumi_initiate_input_t *pumi_inputs, int nsubmeshes);
void pumi_getmeshparameters_from_terminal(int *dimension, int *submesh_num, double ***submesh_params, unsigned int **submesh_flag);
void pumi_setsubmesh(pumi_mesh_t *pumi_mesh, int isubmesh, double xleft, double xright, unsigned int submeshflag,
  int N_uniform, double T_left, double r_left, int N_left, double T_right, double r_right, int N_right);
unsigned int pumi_getsubmeshflag(char flagstring[30]);
void pumi_finalize(pumi_mesh_t* pumi_mesh);
void deallocate_submesh_flag(unsigned int *submesh_flag);
void deallocate_submesh_param(int dim, double **submesh_param);

#endif /* pumi_initiate_h */
