#ifndef pumi_initiate_h
#define pumi_initiate_h

#include "pumiMBBL.h"

#define SUBMESH_FLAGSTRING_LENGTH 23 //!< Maximum length of the mesh flag string input for each submesh block (corresponds to "uniform&leftBL&rightBL")
#define SUBMESH_MAX_SEGMENTS 3 //!< Maximum number of segments allowed in submesh block
#define SEGMENT_STRING_LENGTH 8 //!< Maximum length of string that defines the mesh type in a segment of a submesh block
#define MAX_SUBMESHES 100 //!< Maximum number of submeshes allowed in MSBL
#define DEFAULT_PARAM_VAL 0.0 //!< Default values assigned to submesh parameters

/*!
* \brief Initiate flag enum, possible ways to set the mesh parameters
*/
typedef enum pumi_initiate_flag{
  initiate_from_terminal, //!< Inputs supplied by user in the terminal. The user will be prompted to input mesh parameters
  initiate_from_commandline_inputs, //!< Mesh parameters will be read from the 'command line inputs' of the hpic code and assigned to struct pumi_initiate_input members
} pumi_initiate_flag_t;

/*!
* \brief Contains the parameters inputs to the mesh which will be passed as arguments to pumi_initiate()
*/
typedef struct pumi_initiate_input{
  int ndim; //!< number of physical dimensions of the problem space
  int nsubmeshes; //!< number of submesh blocks in the domain
  int *Nel_max; //!< maximum number of elements in a submesh mesh
  double *alpha; //!< Multiplicative factor to determine Nel_max for a submesh
  int *Nel_max_FLAG; //!< Flag to specify type of input for Nel_max i.e with or without alpha for a submesh
  int *p1_i;//! Number of debye lenghts in a submesh
  int *Nel_i;//!< Number of cells in a submesh block
  double *p2max_i;//!< Number of maximum size cells in a Debye Length
  double *p2min_i;//!< Number of minimum size cells in a Debye Length
  int p1_l;//!< Number of debye lengths in leftBL segment
  int p1_r;//!< Number of debye lengths in rightBL segment
  char **type_flag; //!< pointer to array of mesh flag strings of each submesh block
  double *x_left; //!< pointer to array of left end coordinates of each submesh block
  double *x_right; //!< pointer to array of right end coordinates of each submesh block
  int *uniform_Nel; //!< pointer to array of number of elements in the uniform mesh segment for each submesh block
  double *left_T; //!< pointer to array of left BL segment thickness for each submesh block
  double *left_r; //!< pointer to array of growth ratios in the left BL segment for each block
  int *left_Nel; //!< pointer to array of number of elements in the left BL segment for each block
  double *right_T; //!< pointer to array of right BL segment thickness for each submesh block
  double *right_r; //!<  pointer to array of growth ratios in the right BL segment for each block
  int *right_Nel; //!< pointer to array of number of elements in the right BL segment for each block
} pumi_initiate_input_t;
// remove unwanted variables and do more optimizations
pumi_mesh_t* pumi_initiate(pumi_initiate_flag_t pumi_input_initiate_flag, pumi_initiate_input_t *pumi_inputs);
void pumi_inputs_allocate(pumi_initiate_input_t *pumi_inputs, int nsubmeshes);
void pumi_inputs_deallocate(pumi_initiate_input_t *pumi_inputs, int nsubmeshes);
void pumi_getmeshparameters_from_terminal(int *dimension, int *submesh_num, double ***submesh_params, unsigned int **submesh_flag);
void pumi_freemeshparameters_from_terminal(int nsubmeshes, double **submesh_params, unsigned int *submesh_flag);
void pumi_setsubmesh(pumi_mesh_t *pumi_mesh, int isubmesh, double xleft, double xright, unsigned int submeshflag,
  int N_uniform, double T_left, double r_left, int N_left, double T_right, double r_right, int N_right);
unsigned int pumi_getsubmeshflag(char flagstring[SUBMESH_FLAGSTRING_LENGTH]);
void pumi_finalize(pumi_mesh_t* pumi_mesh);
double pumi_compute_grading_ratio_new(double BL_T, double BL_t0, int BL_Nel);
double pumi_compute_grading_ratio(int p1_lr, int p2, int BL_Nel);
void pumi_verify_params(pumi_mesh_t *pumi_mesh);
void pumi_print_node_coordinates(pumi_mesh_t *pumi_mesh);

#endif /* pumi_initiate_h */
