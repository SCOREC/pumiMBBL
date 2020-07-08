/*!
* \author Vignesh Vittal-Srinivasaragavan
* \date 07-11-2019
* \mainpage Multi-block Boundary Layer PUMI mesh
*/

#ifndef pumiMBBL_h
#define pumiMBBL_h

/*!
* \brief Mesh flag enum that defines the types of meshing in each segment of the submesh block
*/
typedef enum pumi_meshflag{
  unassigned = 0x00, //!< Default flag for each submesh block
  uniform = 0x01, //!< Inidicates the presence of uniform mesh segment in a submesh block
  leftBL = 0x02, //!< Inidicates the presence of left Boundary Layer (BL) segment in a submesh block
  rightBL = 0x04, //!< Inidicates the presence of rigth BL segment in a submesh block
} pumi_meshflag_t;

/*!
* \brief Contains the parameters used to define a submesh
*/
typedef struct pumi_submesh{
  double x_min; //!< coordinate of min end of submesh block
  double x_max; //!< coordinate of max end of submesh block

  double submesh_T; //!< thickness of the submesh block
  int submesh_Nel; //!< number of elements in the submesh block
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
  int nsubmeshes_x; //!< number of submesh blocks in the domain
  int ndim; //!< number of physical dimensions of the problem space
  void *pumi_submeshes_x; //!< pointer object to access members of the structs pumi_submesh1D and pumi_submesh2D
  int pumi_Nel_total_x; //!< total number of elements in the mesh
  int BL_elem_coords_cache_flag;// !< BL elem size and coords precompute flag -- 0=>BL elemsize and node coords array not precomputed, 1=>BL elemsize and node coords array precomputed
} pumi_mesh_t;

#include "pumi_initiate.h"
#include "pumi_routines.h"

#endif /* pumiMBBL_h */
