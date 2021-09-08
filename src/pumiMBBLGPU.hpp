/*!
* \author Vignesh Vittal-Srinivasaragavan
* \date 05-19-2021
* \mainpage Multi-block Boundary Layer PUMI mesh (with Kokkos)
*/

#ifndef pumiMBBLGPU_hpp
#define pumiMBBLGPU_hpp

#include <iostream>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <ctime>

#include <Kokkos_Core.hpp>
#include "pumiMBBL_utils.hpp"

#define MAX_DIM 2
#define MAX_SUBMESHES 100

namespace pumi {
///////// Mesh-Initiate-Structs  ///////////////////////////////////////

/*!
* \brief enum type for directions
*/
enum directions{
    x1_dir = 0, //!< x1-direction
    x2_dir = 1, //!< x2-direction
    x3_dir = 2, //!< x3-direction
};


/*!
* \brief Mesh type enum that defines the types of meshing in each submesh block
*/
enum Meshtype{
    unassigned = 0x00, //!< Default type for each block
    uniform    = 0x01, //!< Uniform mesh
    minBL      = 0x02, //!< geometrically graded BL mesh biased towards the min-size (i.e left-side or bottom-side)
    maxBL      = 0x04, //!< geometrically graded BL mesh biased towards the max-size (i.e right-side or top-side)
};

/**
 * @brief Mesh input struct
 *
 * Object of this class will be passed mesh_initialize API to initiate the mesh
 */
struct Mesh_Inputs{
    int ndim; //!< number of physical dimensions of the problem space
    // 2D params (input from commandline)
    int nsubmesh_x1; //!< number of x1-submesh blocks in the domain
    int nsubmesh_x2; //!< number of x2-submesh blocks in the domain
    bool isactive[MAX_SUBMESHES][MAX_SUBMESHES];
    std::vector<double> block_length_x1;//! Number of debye lenghts in a x1-submesh
    std::vector<double> block_length_x2;//! Number of debye lenghts in a x2-submesh
    std::vector<double> max_elem_size_x1;//!< Maximum size cells in Debye Length (along x1-direction)
    std::vector<double> max_elem_size_x2;//!< Maximum size cells in Debye Length (along x2-direction)
    std::vector<double> min_elem_size_x1;//!< Minimum size cells in Debye Length (along x1-direction)
    std::vector<double> min_elem_size_x2;//!< Minimum size cells in Debye Length (along x2-direction)
    std::vector<std::string> meshtype_x1; //!< Type of mesh as string (uniform/minBL/maxBL)
    std::vector<std::string> meshtype_x2; //!< Type of mesh as string (uniform/minBL/maxBL)
};

/*!
* \brief enum type for option to store BL coords
*/
enum store_BL_coords{
    store_BL_coords_OFF = 0, //!< option to disable storing explicit node coords for BL blocks
    store_BL_coords_ON  = 1, //!< option to enable storing explicit node coords for BL blocks
};

/*!
* \brief struct of mesh options
*/
struct Mesh_Options{
    store_BL_coords BL_storage_option;
    /*!
    * \brief Struct default constructor
    */
    Mesh_Options():BL_storage_option(store_BL_coords_ON){};
};


///////// PUMI-MBBL-GPU Data-Structures ////////////////////////////////////////////////////

using DoubleViewPtr = Kokkos::View<double*>;

/**
 * @brief Base class for submesh-blocks which stores parameters defining the mesh in a block
 *
 * Different meshtypes are derived from this class
 */
class Submesh{
public:
    // SUBMESH PARAMS
    double xmin; //!< Min-side coordinates
    double xmax; //!< Max-side coordinates
    int Nel; //!< Number of elements in the block
    double t0; //!< Smallest element size in the block (or element size for uniform blocks)
    Meshtype meshtype; //<! Type of meshing in the block
    // DEPENDENT PARAMETERS
    double r; //!< Grading ratio for element sizes in the block
    double length; //!< Length of the block
    int Nel_cumulative; //!< Number of elements in the preceding blocks

    double r_by_t0; //!< value of (r-1.0)/t0 -- value needed in analytical cell-locate functions
    double log_r; //!< value of log(r) -- value needed in analtyical cell-locate functions
    DoubleViewPtr BL_coords; //<! BL coords to be stored (optional) for BL blocks for faster particle search
    /**
    * @brief Default constructor.
    */
    Submesh(){};
    /**
    * @brief Class constructor.
    *
    * \param[in] submesh min-side coords
    * \param[in] submesh max-side coords
    * \param[in] submesh number of elements
    * \param[in] submesh smallest element size
    * \param[in] submesh grading ratio
    * \param[in] submesh mesh-type
    * \param[in] submesh length
    * \param[in] submesh preceding cumulative elements
    * \param[in] submesh (r-1.0)/t0 value
    * \param[in] submesh log(r) value
    * \param[in] submesh BL coordinates (explicitly stored)
    */
    Submesh(double xmin_,
            double xmax_,
            int Nel_,
            double t0_,
            double r_,
            Meshtype meshtype_,
            double length_,
            int Nel_cumulative_,
            double r_by_t0_,
            double log_r_,
            DoubleViewPtr BL_coords_):
            xmin(xmin_),
            xmax(xmax_),
            Nel(Nel_),
            t0(t0_),
            r(r_),
            meshtype(meshtype_),
            length(length_),
            Nel_cumulative(Nel_cumulative_),
            r_by_t0(r_by_t0_),
            log_r(log_r_),
            BL_coords(BL_coords_){};

    KOKKOS_INLINE_FUNCTION
    virtual int locate_cell(double q) { return -1; }

    KOKKOS_INLINE_FUNCTION
    virtual int update_cell(double q, int icell) { return -1; }

    KOKKOS_INLINE_FUNCTION
    virtual double elem_size(int icell) { return -999.0; }

    KOKKOS_INLINE_FUNCTION
    virtual void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *global_cell = -1;
        *Wgh2 = -999.0;
    }

};

using SubmeshDeviceViewPtr = Kokkos::View<DevicePointer<Submesh>*>; //!< Pointer to array of Submesh objects (in device space) poniting to address of derived class
using SubmeshHostViewPtr = Submesh*; //<! Pointer to array of submesh objects (in CPU memory)
/**
 * @brief Uniform submesh class derived from submesh class
 *
 */
class Uniform_Submesh : public Submesh{
public:
    /**
    * @brief Class constructor.
    *
    * \param[in] submesh min-side coords
    * \param[in] submesh max-side coords
    * \param[in] submesh number of elements
    * \param[in] submesh smallest element size
    * \param[in] submesh grading ratio
    * \param[in] submesh length
    * \param[in] submesh preceding cumulative elements
    * \param[in] submesh (r-1.0)/t0 value
    * \param[in] submesh log(r) value
    * \param[in] submesh BL coordinates (explicitly stored)
    */
    Uniform_Submesh(double xmin_,
                    double xmax_,
                    int Nel_,
                    double t0_,
                    double r_,
                    double length_,
                    int Nel_cumulative_,
                    double r_by_t0_,
                    double log_r_,
                    DoubleViewPtr BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,uniform,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_){};

    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        return (q - this->xmin)/this->t0;
    }

    KOKKOS_INLINE_FUNCTION
    int update_cell(double q, int icell){
        return (q - this->xmin)/this->t0;
    }

    KOKKOS_INLINE_FUNCTION
    double elem_size(int icell){
        return this->t0;
    }

    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = ((q - this->xmin) - local_cell*this->t0)/this->t0;
        *global_cell = local_cell + this->Nel_cumulative;
    }
};
/**
 * @brief MinBL submesh class derived from submesh class
 *
 */
class MinBL_Submesh : public Submesh{
public:
    /**
    * @brief Class constructor.
    *
    * \param[in] submesh min-side coords
    * \param[in] submesh max-side coords
    * \param[in] submesh number of elements
    * \param[in] submesh smallest element size
    * \param[in] submesh grading ratio
    * \param[in] submesh length
    * \param[in] submesh preceding cumulative elements
    * \param[in] submesh (r-1.0)/t0 value
    * \param[in] submesh log(r) value
    * \param[in] submesh BL coordinates (explicitly stored)
    */
    MinBL_Submesh(double xmin_,
                    double xmax_,
                    int Nel_,
                    double t0_,
                    double r_,
                    double length_,
                    int Nel_cumulative_,
                    double r_by_t0_,
                    double log_r_,
                    DoubleViewPtr BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,minBL,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_){};

    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        int cell = log(1.0 + (q - this->xmin)*this->r_by_t0)/this->log_r;
        return cell;
    }

    KOKKOS_INLINE_FUNCTION
    int update_cell(double q, int icell){
        while (q < this->BL_coords(icell)){
            icell--;
        }
        while (q > this->BL_coords(icell+1)){
            icell++;
        }
        return icell;
    }

    KOKKOS_INLINE_FUNCTION
    double elem_size(int icell){
        return (this->BL_coords(icell+1)-this->BL_coords(icell));
        // return this->t0*pow(this->r,icell);
    }

    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->BL_coords(local_cell))/(this->BL_coords(local_cell+1)-this->BL_coords(local_cell));
        // double r_pow_cell = pow(this->r, local_cell);
        // *Wgh2 = (q - (this->xmin + (r_pow_cell-1.0)/this->r_by_t0))/(r_pow_cell*this->t0);
        *global_cell = local_cell + this->Nel_cumulative;
    }

};

/**
 * @brief MaxBL submesh class derived from submesh class
 *
 */
class MaxBL_Submesh : public Submesh{
public:
    /**
    * @brief Class constructor.
    *
    * \param[in] submesh min-side coords
    * \param[in] submesh max-side coords
    * \param[in] submesh number of elements
    * \param[in] submesh smallest element size
    * \param[in] submesh grading ratio
    * \param[in] submesh length
    * \param[in] submesh preceding cumulative elements
    * \param[in] submesh (r-1.0)/t0 value
    * \param[in] submesh log(r) value
    * \param[in] submesh BL coordinates (explicitly stored)
    */
    MaxBL_Submesh(double xmin_,
                    double xmax_,
                    int Nel_,
                    double t0_,
                    double r_,
                    double length_,
                    int Nel_cumulative_,
                    double r_by_t0_,
                    double log_r_,
                    DoubleViewPtr BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,maxBL,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_){};

    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        int cell = log(1.0 + (this->xmax - q)*this->r_by_t0)/this->log_r;
        return this->Nel - cell - 1;
    }

    KOKKOS_INLINE_FUNCTION
    int update_cell(double q, int icell){
        while (q < this->BL_coords(icell)){
            icell--;
        }
        while (q > this->BL_coords(icell+1)){
            icell++;
        }
        return icell;
    }

    KOKKOS_INLINE_FUNCTION
    double elem_size(int icell){
        return (this->BL_coords(icell+1)-this->BL_coords(icell));
        // return this->t0*pow(this->r,this->Nel-icell-1);
    }

    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->BL_coords(local_cell))/(this->BL_coords(local_cell+1)-this->BL_coords(local_cell));
        *global_cell = local_cell + this->Nel_cumulative;
        // local_cell = this->Nel-local_cell-1;
        // double r_pow_cell = pow(this->r, local_cell);
        // *Wgh2 = 1.0 - ((this->xmax - (r_pow_cell-1.0)/this->r_by_t0) - q)/(this->t0*r_pow_cell);
    }

};

/**
 * @brief MaxBL submesh class derived from submesh class
 *
 */
class Unassigned_Submesh : public Submesh{
public:
    /**
    * @brief Class constructor.
    *
    * \param[in] submesh min-side coords
    * \param[in] submesh max-side coords
    * \param[in] submesh number of elements
    * \param[in] submesh smallest element size
    * \param[in] submesh grading ratio
    * \param[in] submesh length
    * \param[in] submesh preceding cumulative elements
    * \param[in] submesh (r-1.0)/t0 value
    * \param[in] submesh log(r) value
    * \param[in] submesh BL coordinates (explicitly stored)
    */
    Unassigned_Submesh(double xmin_,
                    double xmax_,
                    int Nel_,
                    double t0_,
                    double r_,
                    double length_,
                    int Nel_cumulative_,
                    double r_by_t0_,
                    double log_r_,
                    DoubleViewPtr BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,unassigned,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_){};

    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        return -1;
    }

    KOKKOS_INLINE_FUNCTION
    int update_cell(double q, int icell){
        return -1;
    }

    KOKKOS_INLINE_FUNCTION
    double elem_size(int icell){
        return -999.0;
    }

    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = -999.0;
        *global_cell = -1;
    }

};


/**
 * @brief Mesh class
 *
 * Object of this class can access all properties of the mesh and its submesh
 */
class Mesh{
public:
    int ndim; //!< dimensions of the domain
    int nsubmesh_x1; //!< number of blocks in x1-direction
    int nsubmesh_x2; //!< number of blocks in x2-direction
    int nsubmesh_x3; //!< number of blocks in x3-direction

    Kokkos::View<bool**> isactive; //!< 2D bool-array defining the activity of blocks
    bool **host_isactive;

    // Kokkos::View<int**> nodeoffset;
    Kokkos::View<int**> nodeoffset_start; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> nodeoffset_skip_bot; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> nodeoffset_skip_mid; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> nodeoffset_skip_top; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> elemoffset_start; //!< aux data structure to compute elemtent offset
    Kokkos::View<int*> elemoffset_skip; //!< aux data structure to compute element offset

    Kokkos::View<bool*> is_bdry; //!< bool value stores if an edge is on boundary
    Kokkos::View<double*[3]> bdry_normal; //!< boundary normal direction 

    int Nel_tot_x1; //!< Total number of elements in x1-direction
    int Nel_tot_x2; //!< Total number of elements in x2-direction
    int Nel_tot_x3; //!< Total number of elements in x3-direction

    int Nel_total; //!< Total active elements in domain
    int Nnp_total; //!< Total active nodes in domain

    /**
    * @brief Default constructor.
    */
    KOKKOS_INLINE_FUNCTION
    Mesh(){};
    /**
    * @brief Constructor for 1D Mesh
    * \param[in] number of x1-submesh blocks
    * \param[in] pointer object for x1-submesh blocks
    * \param[in] total number of elements along x1-direction
    */
    KOKKOS_INLINE_FUNCTION
    Mesh(int nsubmesh_x1_,
         int Nel_tot_x1_):
         ndim(1),
         nsubmesh_x1(nsubmesh_x1_),
         Nel_tot_x1(Nel_tot_x1_){
             Nel_total = Nel_tot_x1_;
             Nnp_total = Nel_tot_x1_+1;
             nsubmesh_x2 = 0;
             nsubmesh_x3 = 0;
             Nel_tot_x2 = 0;
             Nel_tot_x3 = 0;
         };
     /**
     * @brief Constructor for 2D Mesh
     * \param[in] number of x1-submesh blocks
     * \param[in] pointer object for x1-submesh blocks
     * \param[in] total number of elements along x1-direction
     * \param[in] number of x2-submesh blocks
     * \param[in] pointer object for x2-submesh blocks
     * \param[in] total number of elements along x2-direction
     * \param[in] submesh activity info
     */
     KOKKOS_INLINE_FUNCTION
     Mesh(int nsubmesh_x1_,
         int Nel_tot_x1_,
         int nsubmesh_x2_,
         int Nel_tot_x2_,
         Kokkos::View<bool**> isactive_,
         Kokkos::View<int**> elemoffset_start_,
         Kokkos::View<int*> elemoffset_skip_,
         Kokkos::View<int**> nodeoffset_start_,
         Kokkos::View<int**> nodeoffset_skip_bot_,
         Kokkos::View<int**> nodeoffset_skip_mid_,
         Kokkos::View<int**> nodeoffset_skip_top_,
         Kokkos::View<bool*> is_bdry_,
         Kokkos::View<double*[3]> bdry_normal_,
         int Nel_total_,
         int Nnp_total_,
         bool** host_isactive_):
         ndim(2),
         nsubmesh_x1(nsubmesh_x1_),
         Nel_tot_x1(Nel_tot_x1_),
         nsubmesh_x2(nsubmesh_x2_),
         Nel_tot_x2(Nel_tot_x2_),
         isactive(isactive_),
         elemoffset_start(elemoffset_start_),
         elemoffset_skip(elemoffset_skip_),
         nodeoffset_start(nodeoffset_start_),
         nodeoffset_skip_bot(nodeoffset_skip_bot_),
         nodeoffset_skip_mid(nodeoffset_skip_mid_),
         nodeoffset_skip_top(nodeoffset_skip_top_),
         is_bdry(is_bdry_),
         bdry_normal(bdry_normal_),
         Nel_total(Nel_total_),
         Nnp_total(Nnp_total_),
         host_isactive(host_isactive_)
         {
             nsubmesh_x3 = 0;
             Nel_tot_x3 = 0;
         };
};

using MeshDeviceViewPtr = Kokkos::View<Mesh*>;

/**
 * @brief Wrapper structure containing mesh and submesh objects
 *
 * Object of this class will be used to call all mesh related APIs
 */
struct MBBL{
    MeshDeviceViewPtr mesh; //!< Mesh object allocated in device space
    SubmeshDeviceViewPtr submesh_x1;//!< X1-Submesh object allocated in device space
    SubmeshHostViewPtr host_submesh_x1;//!< COPY of X1-Submesh object allocated in host space
    SubmeshDeviceViewPtr submesh_x2;//!< X1-Submesh object allocated in device space
    SubmeshHostViewPtr host_submesh_x2;//!< COPY of X2-Submesh object allocated in host space

    /**
    * @brief Default constructor.
    */
    MBBL(){};

    /**
    * @brief Constructor for 1D Wrapper structure
    * \param[in] Mesh object in allocated in GPU
    * \param[in] x1-submesh object in GPU
    * \param[in] copy of x1-submesh object in CPU
    */
    MBBL(MeshDeviceViewPtr mesh_,
         SubmeshDeviceViewPtr submesh_x1_,
         SubmeshHostViewPtr host_submesh_x1_):
         mesh(mesh_),
         submesh_x1(submesh_x1_),
         host_submesh_x1(host_submesh_x1_){};

    /**
    * @brief Constructor for 2D Wrapper structure
    * \param[in] Mesh object in allocated in GPU
    * \param[in] x1-submesh object in GPU
    * \param[in] copy of x1-submesh object in CPU
    * \param[in] x2-submesh object in GPU
    * \param[in] copy of x2-submesh object in CPU
    */
    MBBL(MeshDeviceViewPtr mesh_,
         SubmeshDeviceViewPtr submesh_x1_,
         SubmeshHostViewPtr host_submesh_x1_,
         SubmeshDeviceViewPtr submesh_x2_,
         SubmeshHostViewPtr host_submesh_x2_):
         mesh(mesh_),
         submesh_x1(submesh_x1_),
         host_submesh_x1(host_submesh_x1_),
         submesh_x2(submesh_x2_),
         host_submesh_x2(host_submesh_x2_){};
};


/**
* @brief Locate the submesh ID and local cell ID for a given x1-coordinate
* Uses analytical formulae to locate the input coordinate
* \param[in] Object of the wrapper mesh structure
* \param[in] x1-coordinate to be located
* \param[out] located x1-submesh ID
* \param[out] located x1-localcell ID
*/
KOKKOS_INLINE_FUNCTION
void locate_submesh_and_cell_x1(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    // int nsubmesh = pumi_obj.submesh_x1.extent(0);
    int nsubmesh = pumi_obj.mesh(0).nsubmesh_x1;
    for (isubmesh=2; isubmesh<=nsubmesh; isubmesh++){
     if (q < (pumi_obj.submesh_x1(isubmesh)()->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh;
    }
    *cellID  = pumi_obj.submesh_x1(*submeshID)()->locate_cell(q);
}

/**
* @brief Locate the submesh ID and local cell ID for a given x2-coordinate
* Uses analytical formulae to locate the input coordinate
* \param[in] Object of the wrapper mesh structure
* \param[in] x2-coordinate to be located
* \param[out] located x2-submesh ID
* \param[out] located x2-localcell ID
*/
KOKKOS_INLINE_FUNCTION
void locate_submesh_and_cell_x2(MBBL pumi_obj, double q, int* submeshID, int *cellID){
    int isubmesh;
    int submesh_located = 0;
    // int nsubmesh = pumi_obj.submesh_x2.extent(0);
    int nsubmesh = pumi_obj.mesh(0).nsubmesh_x2;
    for (isubmesh=2; isubmesh<=nsubmesh; isubmesh++){
     if (q < (pumi_obj.submesh_x2(isubmesh)()->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh;
    }
    *cellID  = pumi_obj.submesh_x2(*submeshID)()->locate_cell(q);
}

/**
 * @brief Update the submesh ID and local cell ID for a given x1-coordinate
 * based on previous submesh and cell IDs.
 * Uses adjacency search to update the IDs
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x1-coordinate
 * \param[in] old x1-submesh ID
 * \param[in] old x1-localcell ID
 * \param[out] updated x1-submesh ID
 * \param[out] updated x1-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_cell_x1(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x1(*submeshID)()->xmin)){
        *submeshID -= 1;
        prev_cellID = pumi_obj.submesh_x1(*submeshID)()->Nel - 1;
    }

    while(q > (pumi_obj.submesh_x1(*submeshID)()->xmax)){
        *submeshID += 1;
        prev_cellID = 0;
    }
    *cellID = pumi_obj.submesh_x1(*submeshID)()->update_cell(q, prev_cellID);
}


/**
 * @brief Update the submesh ID and local cell ID for a given x2-coordinate
 * based on previous submesh and cell IDs.
 * Uses adjacency search to update the IDs
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x2-coordinate
 * \param[in] old x2-submesh ID
 * \param[in] old x2-localcell ID
 * \param[out] updated x2-submesh ID
 * \param[out] updated x2-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_cell_x2(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x2(*submeshID)()->xmin)){
        *submeshID -= 1;
        prev_cellID = pumi_obj.submesh_x2(*submeshID)()->Nel - 1;
    }

    while(q > (pumi_obj.submesh_x2(*submeshID)()->xmax)){
        *submeshID += 1;
        prev_cellID = 0;
    }
    *cellID = pumi_obj.submesh_x2(*submeshID)()->update_cell(q, prev_cellID);
}

/**
 * @brief Update the submesh ID based on previous submesh and
 * locate the local cell ID
 * Uses adjacency search to update the submesh IDs and analytical
 * forumlae for cell IDs. Use this function when BL coords are
 * not stored
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x1-coordinate
 * \param[in] old x1-submesh ID
 * \param[out] updated x1-submesh ID
 * \param[out] updated x1-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_locate_cell_x1(MBBL pumi_obj, double q, int prev_submeshID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x1(*submeshID)()->xmin)){
        *submeshID -= 1;
    }

    while(q > (pumi_obj.submesh_x1(*submeshID)()->xmax)){
        *submeshID += 1;
    }

    *cellID = pumi_obj.submesh_x1(*submeshID)()->locate_cell(q);
}


/**
 * @brief Update the submesh ID based on previous submesh and
 * locate the local cell ID
 * Uses adjacency search to update the submesh IDs and analytical
 * forumlae for cell IDs. Use this function when BL coords are
 * not stored
 * \param[in] Object of the wrapper mesh structure
 * \param[in] new x2-coordinate
 * \param[in] old x2-submesh ID
 * \param[out] updated x2-submesh ID
 * \param[out] updated x2-localcell ID
 */
KOKKOS_INLINE_FUNCTION
void update_submesh_and_locate_cell_x2(MBBL pumi_obj, double q, int prev_submeshID, int *submeshID, int *cellID){
    *submeshID = prev_submeshID;
    while(q < (pumi_obj.submesh_x2(*submeshID)()->xmin)){
        *submeshID -= 1;
    }

    while(q > (pumi_obj.submesh_x2(*submeshID)()->xmax)){
        *submeshID += 1;
    }

    *cellID = pumi_obj.submesh_x2(*submeshID)()->locate_cell(q);
}


/**
 * @brief Computes the partial weights (correspoding to node on the max-side i.e right side)
 * for a located particle coordinate and the global directional cell ID
 * \param[in] Object of the wrapper mesh structure
 * \param[in] x1-coordinate of the particle
 * \param[in] x1-submesh ID of the particle
 * \param[in] x1-localcell ID of the particle
 * \param[out] global cell ID in x1 direction
 * \param[out] partial weight
 */
KOKKOS_INLINE_FUNCTION
void calc_weights_x1(MBBL pumi_obj, double q, int isubmesh, int icell, int *x1_global_cell, double *Wgh2){
    pumi_obj.submesh_x1(isubmesh)()->calc_weights(q, icell, x1_global_cell, Wgh2);
}

/**
 * @brief Computes the partial weights (correspoding to node on the max-side i.e top side)
 * for a located particle coordinate and the global directional cell ID
 * \param[in] Object of the wrapper mesh structure
 * \param[in] x2-coordinate of the particle
 * \param[in] x2-submesh ID of the particle
 * \param[in] x2-localcell ID of the particle
 * \param[out] global cell ID in x2 direction
 * \param[out] partial weight
 */
 KOKKOS_INLINE_FUNCTION
 void calc_weights_x2(MBBL pumi_obj, double q, int isubmesh, int icell, int *x2_global_cell, double *Wgh2){
     pumi_obj.submesh_x2(isubmesh)()->calc_weights(q, icell, x2_global_cell, Wgh2);
 }

 /**
  * @brief Computes the gloabl cell ID and node ID in 2D for a full Mesh
  * with no-inactive blocks (mesh with inactive blocks will need separate implementations)
  * \param[in] global cell ID in x1-direction
  * \param[in] global cell ID in x2-direction
  * \param[out] global cell ID in 2D
  * \param[out] global node ID of the node in left-bottom corner
  * \param[out] global node ID of the node in left-top coner
  */
  KOKKOS_INLINE_FUNCTION
  void calc_global_cellID_and_nodeID_fullmesh(MBBL pumi_obj, int kcell_x1, int kcell_x2, int *global_cell_2D, int *bottomleft_node, int *topleft_node){
      *global_cell_2D = kcell_x1 + kcell_x2*pumi_obj.mesh(0).Nel_tot_x1;
      *bottomleft_node = *global_cell_2D + kcell_x2;
      *topleft_node = *bottomleft_node + pumi_obj.mesh(0).Nel_tot_x1 + 1;
  }

/**
* @brief Computes the gloabl cell ID and node ID in 2D for a full Mesh
* with no-inactive blocks (mesh with inactive blocks will need separate implementations)
* \param[in] global cell ID in x1-direction
* \param[in] global cell ID in x2-direction
* \param[out] global cell ID in 2D
* \param[out] global node ID of the node in left-bottom corner
* \param[out] global node ID of the node in left-top coner
*/
KOKKOS_INLINE_FUNCTION
void calc_global_cellID_and_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node){
    // int nodeoffset_bottom = pumi_obj.mesh(0).nodeoffset(isubmesh,kcell_x2);
    // int nodeoffset_top = pumi_obj.mesh(0).nodeoffset(isubmesh,kcell_x2+1);
    int icell_x2 = kcell_x2 - pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative;
    int elemoffset = pumi_obj.mesh(0).elemoffset_start(isubmesh,jsubmesh) + icell_x2*pumi_obj.mesh(0).elemoffset_skip(jsubmesh);
    *global_cell_2D = kcell_x1 + kcell_x2*pumi_obj.mesh(0).Nel_tot_x1 - elemoffset;
    int nodeoffset_bottom = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh(0).nodeoffset_skip_bot(isubmesh,jsubmesh)
                    +(icell_x2-1)*pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh);
    int nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh);
    if (icell_x2==0){
        nodeoffset_bottom = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh);
        nodeoffset_top = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh(0).nodeoffset_skip_bot(isubmesh,jsubmesh);
    }
    if (icell_x2==pumi_obj.submesh_x2(jsubmesh)()->Nel-1){
        nodeoffset_top = nodeoffset_bottom + pumi_obj.mesh(0).nodeoffset_skip_top(isubmesh,jsubmesh);
    }
    *bottomleft_node = *global_cell_2D + kcell_x2 - nodeoffset_bottom;
    *topleft_node = *bottomleft_node + pumi_obj.mesh(0).Nel_tot_x1 + 1 - nodeoffset_top;
}


KOKKOS_INLINE_FUNCTION
int calc_global_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int Inp, int Jnp){
    int nodeID = Jnp*(pumi_obj.mesh(0).Nel_tot_x1+1) + Inp;
    int jnp = Jnp - pumi_obj.submesh_x2(jsubmesh)()->Nel_cumulative;
    // int nodeoffset = pumi_obj.mesh(0).nodeoffset(isubmesh,Jnp);
    int nodeoffset;
    nodeoffset = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh) + pumi_obj.mesh(0).nodeoffset_skip_bot(isubmesh,jsubmesh)
                    +(jnp-1)*pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh);
    if (jnp==0){
        nodeoffset = pumi_obj.mesh(0).nodeoffset_start(isubmesh,jsubmesh);
    }
    if (jnp==pumi_obj.submesh_x2(jsubmesh)()->Nel){
        nodeoffset +=  (pumi_obj.mesh(0).nodeoffset_skip_top(isubmesh,jsubmesh)-pumi_obj.mesh(0).nodeoffset_skip_mid(isubmesh,jsubmesh));
    }
    return nodeID-nodeoffset;
}

KOKKOS_INLINE_FUNCTION
void push_particle(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain, int *bdry_hit){

    double q1_new = q1+dq1;
    double q2_new = q2+dq2;
    int Nx = pumi_obj.mesh(0).nsubmesh_x1;
    int Nxx = 2*Nx+1;
    // int Ny = pumi_obj.mesh(0).nsubmesh_x2;
    int case_id = (dq2>=0.0)+2*(dq1>=0.0);
    int isub = *isubmesh;
    int jsub = *jsubmesh;
    *in_domain = true;
    if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
        && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){

        *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
        *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
        *bdry_hit = -1;
        return;
    }
    else{
        double del1, del2;
        bool located = false;
        switch (case_id) {
            case 0:
                del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                    isub--;
                }
                else{
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                    jsub--;
                }


                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                            del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                                isub--;
                            }
                            else{
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                                jsub--;
                            }
                        }
                    }
                }

            case 1:
                del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                    isub--;
                }
                else{
                    *bdry_hit = jsub*(Nxx)+isub-1;
                    jsub++;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (q1-pumi_obj.submesh_x1(isub)()->xmin);
                            del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1+Nx;
                                isub--;
                            }
                            else{
                                *bdry_hit = jsub*(Nxx)+isub-1;
                                jsub++;
                            }
                        }
                    }
                }

            case 2:
                del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                    isub++;
                }
                else{
                    *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                    jsub--;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                            del2 = (q2-pumi_obj.submesh_x2(jsub)()->xmin);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                                isub++;
                            }
                            else{
                                *bdry_hit = (jsub-1)*(Nxx)+isub-1;
                                jsub--;
                            }
                        }
                    }
                }

            case 3:
                del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                if (del2/del1 > fabs(dq2/dq1)){
                    *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                    isub++;
                }
                else{
                    *bdry_hit = jsub*(Nxx)+isub-1;
                    jsub++;
                }

                while (!located && in_domain){
                    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                        *in_domain = false;
                        *isubmesh=-1;
                        *icell=-1;
                        *jsubmesh=-1;
                        *jcell=-1;
                        return;
                    }
                    else{
                        if ( q1_new > pumi_obj.submesh_x1(isub)()->xmin && q1_new < pumi_obj.submesh_x1(isub)()->xmax
                            && q2_new > pumi_obj.submesh_x2(jsub)()->xmin && q2_new < pumi_obj.submesh_x2(jsub)()->xmax){
                            *isubmesh = isub;
                            *jsubmesh = jsub;
                            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
                            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
                            *bdry_hit = -1;
                            located = true;
                            return;
                        }
                        else{
                            del1 = (pumi_obj.submesh_x1(isub)()->xmax-q1);
                            del2 = (pumi_obj.submesh_x2(jsub)()->xmax-q2);

                            if (del2/del1 > fabs(dq2/dq1)){
                                *bdry_hit = (jsub-1)*(Nxx)+isub+Nx;
                                isub++;
                            }
                            else{
                                *bdry_hit = jsub*(Nxx)+isub-1;
                                jsub++;
                            }
                        }
                    }
                }
        }

    }

}

KOKKOS_INLINE_FUNCTION
void push_particle_v2(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain, int *bdry_hit){

    double q1_new = q1+dq1;
    double q2_new = q2+dq2;
    int Nx = pumi_obj.mesh(0).nsubmesh_x1;
    int Nxx = 2*Nx+1;
    // int Ny = pumi_obj.mesh(0).nsubmesh_x2;
    int isub = *isubmesh;
    int jsub = *jsubmesh;

    int num_x1_crossed = 0;
    int num_x2_crossed = 0;
    int x1_sub_move = 0;
    int x2_sub_move = 0;
    while(q1_new < (pumi_obj.submesh_x1(isub)()->xmin)){
        isub--;
        num_x1_crossed++;
        x1_sub_move = 2;
        *icell = pumi_obj.submesh_x1(isub)()->Nel-1;
    }
    while(q1_new > (pumi_obj.submesh_x1(isub)()->xmax)){
        isub++;
        num_x1_crossed++;
        x1_sub_move = 1;
        *icell = 0;
    }
    while(q2_new < (pumi_obj.submesh_x2(jsub)()->xmin)){
        jsub--;
        num_x2_crossed++;
        x2_sub_move = 2;
        *jcell = pumi_obj.submesh_x2(jsub)()->Nel-1;
    }
    while(q2_new > (pumi_obj.submesh_x2(jsub)()->xmax)){
        jsub++;
        num_x2_crossed++;
        x2_sub_move = 1;
        *jcell = 0;
    }


    int case_id = x1_sub_move + 3*x2_sub_move;
    // printf("case=%d    \n",case_id );
    double del1, del2;
    int i;

    int isub_tmp = *isubmesh;
    int jsub_tmp = *jsubmesh;

    *isubmesh = isub;
    *jsubmesh = jsub;

    switch (case_id) {
        case 0:
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *isubmesh = isub;
            *jsubmesh = jsub;
            *in_domain = true;
            *bdry_hit = -1;
            return;

        case 1:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
            i=0;
            while (i<num_x1_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d (%d,%d)\n", *bdry_hit,*isubmesh+1,*jsubmesh);
                    isub_tmp += 1;
                    *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *bdry_hit = -1;
            *in_domain = true;
            return;

        case 2:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
            i=0;
            while (i<num_x1_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh-1,*jsubmesh);
                    isub_tmp -= 1;
                    *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *in_domain = true;
            *bdry_hit = -1;
            return;

        case 3:
            *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
            i=0;
            while (i<num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh+1);
                    jsub_tmp += 1;
                    *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *in_domain = true;
            *bdry_hit = -1;
            return;

        case 4:
            del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
            del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                isub_tmp += 1;
            }
            else{
                *bdry_hit = jsub_tmp*(Nxx)+isub_tmp-1;
                jsub_tmp += 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
                    del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                        isub_tmp += 1;
                    }
                    else{
                        *bdry_hit = jsub_tmp*(Nxx)+isub_tmp-1;
                        jsub_tmp += 1;
                    }
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *in_domain = true;
            *bdry_hit = -1;
            return;

        case 5:
            del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
            del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                isub_tmp -= 1;
            }
            else{
                *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
                jsub_tmp += 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
                    del2 = (pumi_obj.submesh_x2(jsub_tmp)()->xmax-q2);

                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                        isub_tmp -= 1;
                    }
                    else{
                        *bdry_hit = (jsub_tmp)*(Nxx)+isub_tmp-1;
                        jsub_tmp += 1;
                    }
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *in_domain = true;
            *bdry_hit = -1;
            return;

        case 6:
            *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
            i=0;
            while (i<num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh-1);
                    jsub_tmp -= 1;
                    *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *in_domain = true;
            *bdry_hit = -1;
            return;

        case 7:
            del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
            del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);
            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                isub_tmp += 1;
            }
            else{
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                jsub_tmp -= 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (pumi_obj.submesh_x1(isub_tmp)()->xmax-q1);
                    del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);
                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp+Nx;
                        isub_tmp += 1;
                    }
                    else{
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                        jsub_tmp -= 1;
                    }
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *in_domain = true;
            *bdry_hit = -1;
            return;

        case 8:
            del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
            del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);

            if (del2/del1 > fabs(dq2/dq1)){
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                isub_tmp -= 1;
            }
            else{
                *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                jsub_tmp -= 1;
            }

            i=0;
            while (i<num_x1_crossed+num_x2_crossed){
                i++;
                if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
                    // printf("hit_bdry = %d\n",*bdry_hit );
                    *in_domain = false;
                    *isubmesh = -1;
                    *jsubmesh = -1;
                    *icell = -1;
                    *jcell = -1;
                    return;
                }
                else{
                    // printf("crossed_edge = %d   (%d,%d)\n", *bdry_hit,*isubmesh,*jsubmesh);
                    del1 = (q1-pumi_obj.submesh_x1(isub_tmp)()->xmin);
                    del2 = (q2-pumi_obj.submesh_x2(jsub_tmp)()->xmin);

                    if (del2/del1 > fabs(dq2/dq1)){
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1+Nx;
                        isub_tmp -= 1;
                    }
                    else{
                        *bdry_hit = (jsub_tmp-1)*(Nxx)+isub_tmp-1;
                        jsub_tmp -= 1;
                    }
                }
            }
            *icell = pumi_obj.submesh_x1(isub)()->update_cell(q1_new, *icell);
            *jcell = pumi_obj.submesh_x2(jsub)()->update_cell(q2_new, *jcell);
            *in_domain = true;
            *bdry_hit = -1;
            return;

    }

}

void check_is_pumi_working(){
    printf("Yes, pumiMBBL-GPU is working\n\n");
}

} // namespace pumi
#include "pumiMBBL_initiate.hpp"
#include "pumiMBBL_routines.hpp"
#endif
