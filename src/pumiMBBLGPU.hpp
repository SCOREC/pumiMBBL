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

    // KOKKOS_INLINE_FUNCTION
    // virtual int locate_cell(double q) { return -1; }
    //
    // KOKKOS_INLINE_FUNCTION
    // virtual int update_cell(double q, int icell) { return -1; }
    //
    // KOKKOS_INLINE_FUNCTION
    // virtual double elem_size(int icell) { return -999.0; }
    //
    // KOKKOS_INLINE_FUNCTION
    // virtual void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
    //     *global_cell = -1;
    //     *Wgh2 = -999.0;
    // }

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
    int** host_nodeoffset_start;
    Kokkos::View<int**> nodeoffset_skip_bot; //!< aux data structure to compute nodeoffset
    int** host_nodeoffset_skip_bot;
    Kokkos::View<int**> nodeoffset_skip_mid; //!< aux data structure to compute nodeoffset
    int** host_nodeoffset_skip_mid;
    Kokkos::View<int**> nodeoffset_skip_top; //!< aux data structure to compute nodeoffset
    int** host_nodeoffset_skip_top;
    Kokkos::View<int**> elemoffset_start; //!< aux data structure to compute elemtent offset
    int** host_elemoffset_start;
    Kokkos::View<int*> elemoffset_skip; //!< aux data structure to compute element offset
    int* host_elemoffset_skip;

    Kokkos::View<bool*> is_bdry; //!< bool value stores if an edge is on boundary
    bool* host_is_bdry;
    Kokkos::View<double*[3]> bdry_normal; //!< boundary normal direction
    double** host_bdry_normal;

    int Nbdry_faces; //!< int value storing number of boundary element faces
    int *edge_to_face;

    int Nel_tot_x1; //!< Total number of elements in x1-direction
    int Nel_tot_x2; //!< Total number of elements in x2-direction
    int Nel_tot_x3; //!< Total number of elements in x3-direction

    int Nel_total; //!< Total active elements in domain
    int Nnp_total; //!< Total active nodes in domain

    /**
    * @brief Default constructor.
    */
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
             Nbdry_faces = 2;
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
         int** host_elemoffset_start_,
         int* host_elemoffset_skip_,
         int** host_nodeoffset_start_,
         int** host_nodeoffset_skip_bot_,
         int** host_nodeoffset_skip_mid_,
         int** host_nodeoffset_skip_top_,
         Kokkos::View<bool*> is_bdry_,
         Kokkos::View<double*[3]> bdry_normal_,
         bool* host_is_bdry_,
         double** host_bdry_normal_,
         int Nbdry_faces_,
         int* edge_to_face_,
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
         host_elemoffset_start(host_elemoffset_start_),
         host_elemoffset_skip(host_elemoffset_skip_),
         host_nodeoffset_start(host_nodeoffset_start_),
         host_nodeoffset_skip_bot(host_nodeoffset_skip_bot_),
         host_nodeoffset_skip_mid(host_nodeoffset_skip_mid_),
         host_nodeoffset_skip_top(host_nodeoffset_skip_top_),
         is_bdry(is_bdry_),
         bdry_normal(bdry_normal_),
         host_is_bdry(host_is_bdry_),
         host_bdry_normal(host_bdry_normal_),
         Nbdry_faces(Nbdry_faces_),
         edge_to_face(edge_to_face_),
         Nel_total(Nel_total_),
         Nnp_total(Nnp_total_),
         host_isactive(host_isactive_)
         {
             nsubmesh_x3 = 0;
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
         Mesh(int nsubmesh_x1_,
             int Nel_tot_x1_,
             int nsubmesh_x2_,
             int Nel_tot_x2_,
             int Nel_total_,
             int Nnp_total_,
             bool** host_isactive_,
             bool* host_is_bdry_,
             double** host_bdry_normal_,
             int Nbdry_faces_,
             int* edge_to_face_,
             int** host_elemoffset_start_,
             int* host_elemoffset_skip_,
             int** host_nodeoffset_start_,
             int** host_nodeoffset_skip_bot_,
             int** host_nodeoffset_skip_mid_,
             int** host_nodeoffset_skip_top_):
             ndim(2),
             nsubmesh_x1(nsubmesh_x1_),
             Nel_tot_x1(Nel_tot_x1_),
             nsubmesh_x2(nsubmesh_x2_),
             Nel_tot_x2(Nel_tot_x2_),
             Nel_total(Nel_total_),
             Nnp_total(Nnp_total_),
             host_isactive(host_isactive_),
             host_is_bdry(host_is_bdry_),
             host_bdry_normal(host_bdry_normal_),
             Nbdry_faces(Nbdry_faces_),
             edge_to_face(edge_to_face_),
             host_elemoffset_start(host_elemoffset_start_),
             host_elemoffset_skip(host_elemoffset_skip_),
             host_nodeoffset_start(host_nodeoffset_start_),
             host_nodeoffset_skip_bot(host_nodeoffset_skip_bot_),
             host_nodeoffset_skip_mid(host_nodeoffset_skip_mid_),
             host_nodeoffset_skip_top(host_nodeoffset_skip_top_)
             {
                 nsubmesh_x3 = 0;
                 Nel_tot_x3 = 0;
             };
};

using MeshDeviceViewPtr = Kokkos::View<Mesh*>;
using MeshHostViewPtr = Mesh*;
/**
 * @brief Wrapper structure containing mesh and submesh objects
 *
 * Object of this class will be used to call all mesh related APIs
 */
struct MBBL{
    MeshDeviceViewPtr mesh; //!< Mesh object allocated in device space
    MeshHostViewPtr host_mesh; //!< Mesh object allocated in host space
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
         host_submesh_x1(host_submesh_x1_){
             MeshDeviceViewPtr::HostMirror h_mesh_ = Kokkos::create_mirror_view(mesh_);
             Kokkos::deep_copy(h_mesh_, mesh_);
             host_mesh = new Mesh(h_mesh_(0).nsubmesh_x1,h_mesh_(0).Nel_tot_x1);
         };

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
         host_submesh_x2(host_submesh_x2_){
             MeshDeviceViewPtr::HostMirror h_mesh_ = Kokkos::create_mirror_view(mesh_);
             Kokkos::deep_copy(h_mesh_, mesh_);
             host_mesh = new Mesh(h_mesh_(0).nsubmesh_x1,h_mesh_(0).Nel_tot_x1,h_mesh_(0).nsubmesh_x2,h_mesh_(0).Nel_tot_x2,
                                    h_mesh_(0).Nel_total,h_mesh_(0).Nnp_total,h_mesh_(0).host_isactive,h_mesh_(0).host_is_bdry,
                                    h_mesh_(0).host_bdry_normal,h_mesh_(0).Nbdry_faces,h_mesh_(0).edge_to_face, h_mesh_(0).host_elemoffset_start,
                                    h_mesh_(0).host_elemoffset_skip,h_mesh_(0).host_nodeoffset_start,h_mesh_(0).host_nodeoffset_skip_bot,
                                    h_mesh_(0).host_nodeoffset_skip_mid,h_mesh_(0).host_nodeoffset_skip_top);
         };
};

KOKKOS_FUNCTION
int locate_cell(DevicePointer<Submesh> submesh, double q);

KOKKOS_FUNCTION
int update_cell(DevicePointer<Submesh> submesh, double q, int icell);

KOKKOS_FUNCTION
double elem_size(DevicePointer<Submesh> submesh, int icell);

KOKKOS_FUNCTION
void calc_weights(DevicePointer<Submesh> submesh, double q, int local_cell, int *global_cell, double *Wgh2);

KOKKOS_FUNCTION
void locate_submesh_and_cell_x1(MBBL pumi_obj, double q, int* submeshID, int *cellID);

KOKKOS_FUNCTION
void locate_submesh_and_cell_x2(MBBL pumi_obj, double q, int* submeshID, int *cellID);

KOKKOS_FUNCTION
void update_submesh_and_cell_x1(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID);

KOKKOS_FUNCTION
void update_submesh_and_cell_x2(MBBL pumi_obj, double q, int prev_submeshID, int prev_cellID, int *submeshID, int *cellID);

KOKKOS_FUNCTION
void update_submesh_and_locate_cell_x1(MBBL pumi_obj, double q, int prev_submeshID, int *submeshID, int *cellID);

KOKKOS_FUNCTION
void update_submesh_and_locate_cell_x2(MBBL pumi_obj, double q, int prev_submeshID, int *submeshID, int *cellID);

KOKKOS_FUNCTION
void calc_weights_x1(MBBL pumi_obj, double q, int isubmesh, int icell, int *x1_global_cell, double *Wgh2);

KOKKOS_FUNCTION
void calc_weights_x2(MBBL pumi_obj, double q, int isubmesh, int icell, int *x2_global_cell, double *Wgh2);

KOKKOS_FUNCTION
void calc_global_cellID_and_nodeID_fullmesh(MBBL pumi_obj, int kcell_x1, int kcell_x2, int *global_cell_2D, int *bottomleft_node, int *topleft_node);

KOKKOS_FUNCTION
void calc_global_cellID_and_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int kcell_x1, int kcell_x2,
                                    int *global_cell_2D, int *bottomleft_node, int *topleft_node);


KOKKOS_FUNCTION
int calc_global_nodeID(MBBL pumi_obj, int isubmesh, int jsubmesh, int Inp, int Jnp);

KOKKOS_FUNCTION
void push_particle(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain, int *bdry_hit);

KOKKOS_FUNCTION
void push_particle_v2(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int *isubmesh, int *jsubmesh, int *icell, int *jcell, bool *in_domain, int *bdry_hit);

KOKKOS_FUNCTION
void get_directional_submeshID_and_cellID(MBBL pumi_obj, int submeshID, int cellID, int* isub, int *icell, int* jsub, int *jcell);

KOKKOS_FUNCTION
int get_x1_cellID(MBBL pumi_obj, int isub, int icell);

KOKKOS_FUNCTION
int get_x2_cellID(MBBL pumi_obj, int isub, int icell);

KOKKOS_FUNCTION
void flatten_submeshID_and_cellID(MBBL pumi_obj, int isub, int icell, int jsub, int jcell, int* submeshID, int* cellID);

KOKKOS_FUNCTION
double get_x1_elem_size_in_submesh(MBBL pumi_obj, int isub, int icell);

KOKKOS_FUNCTION
double get_x2_elem_size_in_submesh(MBBL pumi_obj, int isub, int icell);

//For 1D and 2D tests ONLY -- NOT TO BE USED IN HPIC2
class ParticleData{
public:
    double x1;
    double x2;
    int cellID;
    int submeshID;

    ParticleData(){};

    ParticleData(double x1_,double x2_):
                x1(x1_),x2(x2_){};

    KOKKOS_INLINE_FUNCTION
    ParticleData(double x1_, double x2_, int submeshID_, int cellID_):
                x1(x1_),x2(x2_),cellID(cellID_),submeshID(submeshID_){};

};

} // namespace pumi
#include "pumiMBBL_initiate.hpp"
#include "pumiMBBL_routines.hpp"
#include "pumiMBBLGPU_impl.hpp"
#endif
