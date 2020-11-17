/*!
* \author Vignesh Vittal-Srinivasaragavan
* \date 07-11-2019
* \mainpage Multi-block Boundary Layer PUMI mesh
*/

#ifndef pumiMBBL_h
#define pumiMBBL_h

#include <stdbool.h>

/*!
* \brief Mesh flag enum that defines the types of meshing in each segment of the submesh block
*/
typedef enum pumi_meshflag{
  unassigned = 0x00, //!< Default flag for each submesh block
  uniform = 0x01, //!< Inidicates the presence of uniform mesh segment in a submesh block
  leftBL = 0x02, //!< Inidicates the presence of left Boundary Layer (BL) segment in a submesh block
  bottomBL = 0x02, //!< Inidicates the presence of bottom Boundary Layer (BL) segment in a submesh block
  rightBL = 0x04, //!< Inidicates the presence of rigth BL segment in a submesh block
  topBL = 0x04, //!< Inidicates the presence of top BL segment in a submesh block
} pumi_meshflag_t;

typedef enum pumi_2D_blocktype_for_nodeoffset{
    type_O = 0, //type for inactive blocks
    type_A = 1, //nodeoffset_skip is same for all nodes
    type_B = 2, //nodeoffset_skip is different for first 2 nodes
    type_C = 3, //nodeoffset_skip is different for last node
    type_D = 4, //nodeoffset_skip is diffrent for first 2 nodes and last node
}pumi_2D_blocktype_for_nodeoffset_t;

typedef struct pumi_bezier_extractor{
    double **C; // local element bspline extraction operator
}pumi_bezier_extractor_t;

typedef struct pumi_bspline{
    int N_spline;
    pumi_bezier_extractor_t *pumi_bez_ex_x1;
    pumi_bezier_extractor_t *pumi_bez_ex_x2;
    double *bernstein_vector;
    double *cov_coeffs;
    double *Q_coeffs;
    double *E_coeffs;
    int *nCk4spline;
}pumi_bspline_t;

/*!
* \brief Contains the parameters used to define a submesh
*/
typedef struct pumi_submesh{
  double coord_min; //!< coordinate of min end of submesh block
  double coord_max; //!< coordinate of max end of submesh block

  double submesh_T; //!< thickness of the submesh block
  int submesh_Nel; //!< number of elements in the submesh block
  int submesh_Nel_minus_1;
  double t0; //!< (dependent variable) smallest element size inside the submesh block
  double r; //!< growth ratio for elements inside the submesh block
  double log_r; //!< (dependent variable) log(r) -- value used in analytic locate routines in BL
  double r_t0_ratio; //!< (dependent variable) (left_r-1)/lBL_t0 -- -- value used in analytic locate routines in BL
  double *BL_elemsize; //!< pointer to array that stores elem size in BL blocks
  double *BL_coords; //!< pointer to array that stores node coords in BL blocks

  int Nel_cumulative; //!< (dependent variable) total number of elements in the previous submeshes

  pumi_meshflag_t pumi_flag; //!< flag for types of mesh segments(i.e. uniform mesh segment, right BL segment or left BL segment) available in the submesh block
} pumi_submesh_t;

/*!
* \brief Contains parameters that defines the mesh
*/
typedef struct pumi_mesh{
  int nsubmeshes_x1; //!< number of submesh blocks in the domain in x-drection
  int nsubmeshes_x2; //!< number of submesh blocks in the domain in y-drection
  int ndim; //!< number of physical dimensions of the problem space
  void *pumi_submeshes_x1; //!< pointer object to access members of the struct pumi_submesh (in x-direction)
  void *pumi_submeshes_x2; //!< pointer object to access members of the struct pumi_submesh (in y-direction)
  int pumi_Nel_total_x1; //!< total number of elements in the mesh (along x-directiom)
  int pumi_Nnp_total_x1;
  int pumi_Nel_total_x2; //!< total number of elements in the mesh (along y-directiom)
  int pumi_Nnp_total_x2;
  int pumi_Nel_total_2D;
  int pumi_Nnp_total_2D;
  int BL_elem_coords_cache_flag;// !< BL elem size and coords precompute flag -- 0=>BL elemsize and node coords array not precomputed, 1=>BL elemsize and node coords array precomputed
  bool **isactive;
  int nodeoffset_cache_flag;
  int **nodeoffset_start;
  int **nodeoffset_skip_top;
  int **nodeoffset_skip_mid;
  int **nodeoffset_skip_bottom;
  int **elemoffset_start;
  int *elemoffset_skip;
  int **global_nodeoffset;
  pumi_2D_blocktype_for_nodeoffset_t **blocktype;
  int bspline_flag;
  pumi_bspline_t pumi_bspl;
  int P_spline;
} pumi_mesh_t;

#include "pumi_initiate.h"
#include "pumi_routines.h"

#endif /* pumiMBBL_h */
