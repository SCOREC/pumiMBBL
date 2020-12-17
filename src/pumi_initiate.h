#ifndef pumi_initiate_h
#define pumi_initiate_h

#include "pumiMBBL.h"

#define SUBMESH_FLAGSTRING_LENGTH 23 //!< Maximum length of the mesh flag string input for each submesh block (corresponds to "uniform&leftBL&rightBL")
#define SUBMESH_MAX_SEGMENTS 3 //!< Maximum number of segments allowed in submesh block
#define SUBMESH_MAXTYPES 5 //!< Number of types of meshes available
#define SEGMENT_STRING_LENGTH 9 //!< Maximum length of string that defines the mesh type in a segment of a submesh block
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
* \brief options to calculate the BL element sizes
*/
typedef enum pumi_cache_BL_elemsize{
  pumi_cache_BL_elemsize_ON = 0, //!< no caching of BL element sizes (will be calculated on-the-fly)
  pumi_cache_BL_elemsize_OFF = 1, //!< caching of BL element sizes (precomputed while mesh initialization)
} pumi_cache_BL_elemsize_t;

typedef enum pumi_cache_nodeoffset{
  pumi_cache_nodeoffset_OFF = 0, //!< no caching of BL element sizes (will be calculated on-the-fly)
  pumi_cache_nodeoffset_ON = 1, //!< caching of BL element sizes (precomputed while mesh initialization)
} pumi_cache_nodeoffset_t;

typedef enum pumi_use_bspline{
  pumi_bspline_OFF = 0, //!< no bspline based charge distribution
  pumi_bspline_ON = 1, //!< initiate routines to allow bspline based charge distribution
} pumi_use_bspline_t;

typedef enum pumi_periodic_mesh{
  pumi_periodic_mesh_OFF = 0, //!< regular mesh (no periodicity)
  pumi_periodic_mesh_ON = 1, //!< periodicity present in mesh
} pumi_periodic_mesh_t;

typedef struct pumi_initiate_mesh_options{
    pumi_cache_BL_elemsize_t BL_cache_flag;
    pumi_cache_nodeoffset_t nodeoffset_cache_flag;
    pumi_use_bspline_t bspline_flag;
    pumi_periodic_mesh_t periodic_mesh_flag;
} pumi_initiate_mesh_options_t;

/*!
* \brief Contains the parameters inputs to the mesh which will be passed as arguments to pumi_initiate()
*/
typedef struct pumi_initiate_input{
  int ndim; //!< number of physical dimensions of the problem space
  int P_spline; //!< order of b-spline used for charge distribution
  //1D params
  int nsubmeshes; //!< number of submesh blocks in the domain
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
  // 2D params
  int nsubmeshes_x1; //!< number of submesh blocks in the domain
  int nsubmeshes_x2; //!< number of submesh blocks in the domain
  bool isactive[MAX_SUBMESHES][MAX_SUBMESHES];
  int *p1_i_x1;//! Number of debye lenghts in a submesh
  int *p1_i_x2;//! Number of debye lenghts in a submesh
  int *Nel_i_x1;//!< Number of cells in a submesh block
  int *Nel_i_x2;//!< Number of cells in a submesh block
  double *p2max_i_x1;//!< Number of maximum size cells in a Debye Length
  double *p2max_i_x2;//!< Number of maximum size cells in a Debye Length
  double *p2min_i_x1;//!< Number of minimum size cells in a Debye Length
  double *p2min_i_x2;//!< Number of minimum size cells in a Debye Length
  double *y_bottom; //!< pointer to array of left end coordinates of each submesh block
  double *y_top; //!< pointer to array of right end coordinates of each submesh block
  int *uniform_Nel_x1; //!< pointer to array of number of elements in the uniform mesh segment for each submesh block
  int *uniform_Nel_x2; //!< pointer to array of number of elements in the uniform mesh segment for each submesh block
  double *bottom_T; //!< pointer to array of left BL segment thickness for each submesh block
  double *bottom_r; //!< pointer to array of growth ratios in the left BL segment for each block
  int *bottom_Nel; //!< pointer to array of number of elements in the left BL segment for each block
  double *top_T; //!< pointer to array of right BL segment thickness for each submesh block
  double *top_r; //!<  pointer to array of growth ratios in the right BL segment for each block
  int *top_Nel; //!< pointer to array of number of elements in the right BL segment for each block
} pumi_initiate_input_t;

pumi_mesh_t* pumi_initiate(pumi_initiate_flag_t pumi_input_initiate_flag, pumi_initiate_input_t *pumi_inputs, pumi_initiate_mesh_options_t pumi_initiate_options);
pumi_initiate_input_t* pumi_inputs_allocate(int nsubmeshes);
void pumi_inputs_deallocate(pumi_initiate_input_t *pumi_inputs);
void pumi_getmeshparameters_from_terminal(int *dimension, int *submesh_num, double ***submesh_params, unsigned int **submesh_flag);
void pumi_freemeshparameters_from_terminal(int nsubmeshes, double **submesh_params, unsigned int *submesh_flag);
void pumi_setsubmesh_x1(pumi_mesh_t *pumi_mesh, int isubmesh, double xmin, double xmax, unsigned int submeshflag,
    int N_uniform, double T_minBL, double r_minBL, int N_minBL, double T_maxBL, double r_maxBL, int N_maxBL);
void pumi_setsubmesh_x2(pumi_mesh_t *pumi_mesh, int isubmesh, double xmin, double xmax, unsigned int submeshflag,
    int N_uniform, double T_minBL, double r_minBL, int N_minBL, double T_maxBL, double r_maxBL, int N_maxBL);
void pumi_setsubmesh_elemoffsets(pumi_mesh_t *pumi_mesh);
void pumi_setsubmesh_nodeoffsets(pumi_mesh_t *pumi_mesh);
unsigned int pumi_getsubmeshflag(char flagstring[SUBMESH_FLAGSTRING_LENGTH]);
void pumi_finalize(pumi_mesh_t* pumi_mesh);
double pumi_compute_grading_ratio_new(double BL_T, double BL_t0, int BL_Nel);
//double pumi_compute_grading_ratio(int p1_lr, int p2, int BL_Nel);
void pumi_verify_params(pumi_mesh_t *pumi_mesh);
void pumi_verify_params_1D(pumi_mesh_t *pumi_mesh);
void pumi_verify_params_2D(pumi_mesh_t *pumi_mesh);
void pumi_print_node_coordinates(pumi_mesh_t *pumi_mesh);
void pumi_print_node_coordinates_1D(pumi_mesh_t *pumi_mesh);
void pumi_print_node_coordinates_2D(pumi_mesh_t *pumi_mesh);
int nchoosek(int n, int k);
pumi_bezier_extractor_t* pumi_bezier_extraction(pumi_mesh_t *pumi_mesh, int dir);
pumi_bezier_extractor_t* pumi_bezier_extraction_periodic(pumi_mesh_t *pumi_mesh, int dir);
void pumi_initiate_bsplines(pumi_mesh_t *pumi_mesh, int dir);
pumi_bezier_extractor_t* pumi_unique_bezier_extractor_matrices(pumi_mesh_t* pumi_mesh, int dir, pumi_bezier_extractor_t *pumi_bez_ex_full);
void pumi_finalize_bsplines(pumi_mesh_t* pumi_mesh);
#endif /* pumi_initiate_h */
