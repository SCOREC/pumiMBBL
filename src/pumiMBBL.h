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
typedef struct pumi_submesh1D{
  double x_left; //!< coordinate of left end of submesh block
  double x_right; //!< coordinate of right end of submesh block

  int uniform_Nel; //!< number of elements in the uniform mesh segment block
  double uniform_x_left; //!< (dependent variable) coordinate of left end of uniform mesh segment inside the submesh block
  double uniform_x_right; //!< (dependent variable) coordinate of right end of uniform mesh segment inside the submesh block
  double uniform_t0; //!< (dependent variable) element sizes in uniform mesh segment inside the submesh block

  double left_T; //!< left Boundary Layer (BL) thickness
  double left_r; //!< growth ratio for left BL mesh
  int left_Nel; //!< number of elements in the left BL segment
  double lBL_x_right; //!< (dependent variable) coordinate of right end of left BL segment inside the submesh block
  double lBL_t0; //!< (dependent variable) size of first (leftmost) element in the left BL segment inside the submesh block
  double log_left_r;//!< (dependent variable) log(left_r) -- value used in pumi_locatepoint_BL_1D algo
  double left_r_lBL_t0_ratio;//!< (dependent variable) (left_r-1)/lBL_t0 -- value used in pumi_locatepoint_1D algo
  double *leftBL_elemsize;//!< pointer to array that stores elem size in leftBL
  int leftBL_elemsize_calc_flag; // !< leftBL elem size calculation flag variable. 0=>BL elemsize array not computed, 1=>BL elemsize array computed

  double right_T; //!< right BL thickness
  double right_r; //!<  growth ratio for right BL segment
  int right_Nel; //!< number of elements in the right BL segment
  double rBL_x_left; //!< (dependent variable) coordinate of left end of right graded mesh region inside the submesh block
  double rBL_t0; //!< (dependent variable) size of first (rightmost) element in the right BL segment inside the submesh block
  double log_right_r;//!< (dependent variable) log(right_r) -- value used in pumi_locatepoint_BL_1D algo
  double right_r_rBL_t0_ratio;//!< (dependent variable) (right_r-1)/rBL_t0 -- value used in pumi_locatepoint_1D algo
  double *rightBL_elemsize;//!< pointer to array that stores elem size in rightBL
  int rightBL_elemsize_calc_flag; // !< rightBL elem size calculation flag variable. 0=>BL elemsize array not computed, 1=>BL elemsize array computed

  int submesh_total_Nel; //!< (dependent variable) total number of elements in the submesh block
  int Nel_cumulative; //!< (dependent variable) total number of elements in the previous submeshes

  pumi_meshflag_t pumi_flag; //!< flag for types of mesh segments(i.e. uniform mesh segment, right BL segment or left BL segment) available in the submesh block
} pumi_submesh1D_t;

/*!
* \brief Contains parameters that defines the mesh
*/
typedef struct pumi_mesh{
  int nsubmeshes; //!< number of submesh blocks in the domain
  int ndim; //!< number of physical dimensions of the problem space
  void *pumi_submeshes; //!< pointer object to access members of the structs pumi_submesh1D and pumi_submesh2D
  int pumi_Nel_total; //!< total number of elements in the mesh
} pumi_mesh_t;

#include "pumi_initiate.h"
#include "pumi_routines.h"

#endif /* pumiMBBL_h */
