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
#include "pumi_utils.hpp"

#define MAX_DIM 2
#define MAX_SUBMESHES 100

namespace pumi {
///////// Mesh-Initiate-Structs and Function declarations ///////////////////////////////////////

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
    double *p1_i_x1;//! Number of debye lenghts in a x1-submesh
    double *p1_i_x2;//! Number of debye lenghts in a x2-submesh
    double *p2max_i_x1;//!< Maximum size cells in Debye Length (along x1-direction)
    double *p2max_i_x2;//!< Maximum size cells in Debye Length (along x2-direction)
    double *p2min_i_x1;//!< Minimum size cells in Debye Length (along x1-direction)
    double *p2min_i_x2;//!< Minimum size cells in Debye Length (along x2-direction)
    std::vector<std::string> meshtype; //!< Type of mesh as string (uniform/minBL/maxBL)

    double *x1_min; //!< Min-side coordinate (along x1-direction)
    double *x1_max; //!< Max-side coordinate (along x1-direction)
    double *x2_min; //!< Min-side coordinate (along x2-direction)
    double *x2_max; //!< Max-side coordinate (along x2-direction)
    int *Nel_i_x1;//!< Number of cells in a x1-submesh block
    int *Nel_i_x2;//!< Number of cells in a x2-submesh block
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

/**
 * @brief Reads the meshtype as string and converts it to unsigned int meshtype
 *
 * \param[in] mesh type in string format "uniform" or "minBL" or "maxBL"
 * \return unsigned int mesh type -- uniform=0x01 , minBL=0x02 and maxBL=0x04
 */
unsigned int get_submesh_type(std::string meshtype_string){
    std::vector<std::string> defined_types {"uniform","minBL","maxBL"};
    unsigned int type = unassigned;

    if ( meshtype_string.compare(defined_types[0]) ==0 ){
        type = uniform;
    }
    else if ( meshtype_string.compare(defined_types[1]) ==0 ){
        type = minBL;
    }
    else if ( meshtype_string.compare(defined_types[2]) ==0 ){
        type = maxBL;
    }
    else{
        std::cout << meshtype_string << " is an invalid meshtype\n";
        std::cout << "Valid types are \"uniform\", \"minBL\", \"maxBL\"\n";
        std::cout << "Terminating\n";
        exit(0);
    }
    return type;
}

/**
 * @brief computes grading ratio using newton raphson
 *
 * \param[in] block Length
 * \param[in] smallest element size in the block
 * \param[in] number of elements in the block
 * \return computed grading ratio
 */
double compute_grading_ratio(double BL_T, double BL_t0, int BL_Nel){
    double tol = 1e-5;
    double r = 1.5;
    double del_r = 1.0;
    int iter = 0;
    int max_iter = 10000;
    double f, f_r;
    while (fabs(del_r) > tol){
        f = BL_t0*pow(r,BL_Nel) - BL_T*r + (BL_T-BL_t0);
        f_r = BL_Nel*BL_t0*pow(r,BL_Nel-1) - BL_T;
        del_r = f/f_r;
        r -= del_r;
        iter++;
        if (iter > max_iter){
          printf("Cannot compute grading ratio. Solver not converging\n" );
          exit(0);
        }
    }
    return r;
}

/**
 * @brief Allocates memory for mesh inputs struct members
 *
 * \param[in] sum of number of submesh-blocks in each direction
 */
Mesh_Inputs* inputs_allocate(int nsubmeshes){
    Mesh_Inputs* pumi_inputs = new Mesh_Inputs;
    pumi_inputs->p1_i_x1 = new double[nsubmeshes];
    pumi_inputs->p1_i_x2 = new double[nsubmeshes];
    pumi_inputs->p2max_i_x1 = new double[nsubmeshes];
    pumi_inputs->p2max_i_x2 = new double[nsubmeshes];
    pumi_inputs->p2min_i_x1 = new double[nsubmeshes];
    pumi_inputs->p2min_i_x2 = new double[nsubmeshes];

    return pumi_inputs;
}

/**
 * @brief Frees up the alocated memory of mesh inputs struct members
 *
 * \param[in] mesh inputs struct pointer
 */
void inputs_deallocate(Mesh_Inputs* pumi_inputs){
    delete[] pumi_inputs->p1_i_x1;
    delete[] pumi_inputs->p1_i_x2;
    delete[] pumi_inputs->p2max_i_x1;
    delete[] pumi_inputs->p2max_i_x2;
    delete[] pumi_inputs->p2min_i_x1;
    delete[] pumi_inputs->p2min_i_x2;
    delete pumi_inputs;
}
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
    Submesh(double submesh_xmin,
            double submesh_xmax,
            int submesh_Nel,
            double submesh_t0,
            double submesh_r,
            Meshtype submesh_type,
            double submesh_length,
            int submesh_Nel_cumulative,
            double r_t0_ratio,
            double logr,
            DoubleViewPtr BLcoords):
            xmin(submesh_xmin),
            xmax(submesh_xmax),
            Nel(submesh_Nel),
            t0(submesh_t0),
            r(submesh_r),
            meshtype(submesh_type),
            length(submesh_length),
            Nel_cumulative(submesh_Nel_cumulative),
            r_by_t0(r_t0_ratio),
            log_r(logr),
            BL_coords(BLcoords){};

    KOKKOS_INLINE_FUNCTION
    virtual int locate_cell(double q) { return 0; }

    KOKKOS_INLINE_FUNCTION
    virtual int update_cell(double q, int icell) { return 0; }

    KOKKOS_INLINE_FUNCTION
    virtual double elem_size(int icell) { return 0.0; }

    KOKKOS_INLINE_FUNCTION
    virtual void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *global_cell = 0;
        *Wgh2 = 0.0;
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
    Uniform_Submesh(double submesh_xmin,
                    double submesh_xmax,
                    int submesh_Nel,
                    double submesh_t0,
                    double submesh_r,
                    double submesh_length,
                    int submesh_Nel_cumulative,
                    double r_t0_ratio,
                    double logr,
                    DoubleViewPtr BLcoords):
                    Submesh(submesh_xmin,submesh_xmax,submesh_Nel,submesh_t0,submesh_r,uniform,submesh_length,submesh_Nel_cumulative,r_t0_ratio,logr,BLcoords){};

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
    MinBL_Submesh(double submesh_xmin,
                  double submesh_xmax,
                  int submesh_Nel,
                  double submesh_t0,
                  double submesh_r,
                  double submesh_length,
                  int submesh_Nel_cumulative,
                  double r_t0_ratio,
                  double logr,
                  DoubleViewPtr BLcoords):
                  Submesh(submesh_xmin,submesh_xmax,submesh_Nel,submesh_t0,submesh_r,minBL,submesh_length,submesh_Nel_cumulative,r_t0_ratio,logr,BLcoords){};

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
    MaxBL_Submesh(double submesh_xmin,
                  double submesh_xmax,
                  int submesh_Nel,
                  double submesh_t0,
                  double submesh_r,
                  double submesh_length,
                  int submesh_Nel_cumulative,
                  double r_t0_ratio,
                  double logr,
                  DoubleViewPtr BLcoords):
                  Submesh(submesh_xmin,submesh_xmax,submesh_Nel,submesh_t0,submesh_r,maxBL,submesh_length,submesh_Nel_cumulative,r_t0_ratio,logr,BLcoords){};

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
 * @brief Mesh class
 *
 * Object of this class can access all properties of the mesh and its submesh
 */
class Mesh{
public:
    int ndim; //!< dimensions of the domain
    int nsubmesh_x1; //!< number of blocks in x1-direction
    int nsubmesh_x2; //!< number of blocks in x2-direction
    int nsubmesh_x3; //1< number of blocks in x3-direction

    Kokkos::View<bool**> isactive; //!< 2D bool-array defining the activity of blocks
    bool **host_isactive;

    Kokkos::View<int**> nodeoffset;
    Kokkos::View<int**> nodeoffset_start;
    Kokkos::View<int**> nodeoffset_skip_bot;
    Kokkos::View<int**> nodeoffset_skip_mid;
    Kokkos::View<int**> nodeoffset_skip_top;
    Kokkos::View<int**> elemoffset_start;
    Kokkos::View<int*> elemoffset_skip;

    Kokkos::View<bool*> is_bdry;

    int Nel_tot_x1; //!< Total number of elements in x1-direction
    int Nel_tot_x2; //!< Total number of elements in x2-direction
    int Nel_tot_x3; //!< Total number of elements in x3-direction

    int Nel_total;
    int Nnp_total;

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
    Mesh(int numsubmesh_x1,
         int Nel_total_x1):
         ndim(1),
         nsubmesh_x1(numsubmesh_x1),
         Nel_tot_x1(Nel_total_x1){
             Nel_total = Nel_total_x1;
             Nnp_total = Nel_total_x1+1;
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
    Mesh(int numsubmesh_x1,
         int Nel_total_x1,
         int numsubmesh_x2,
         int Nel_total_x2,
         Kokkos::View<bool**> submesh_activity,
         Kokkos::View<int**> elem_offset_start,
         Kokkos::View<int*> elem_offset_skip,
         Kokkos::View<int**> node_offset,
         Kokkos::View<int**> node_offset_start,
         Kokkos::View<int**> node_offset_skip_bot,
         Kokkos::View<int**> node_offset_skip_mid,
         Kokkos::View<int**> node_offset_skip_top,
         Kokkos::View<bool*> isbdry,
         int Nel_tot,
         int Nnp_tot,
         bool** h_isactive):
         ndim(2),
         nsubmesh_x1(numsubmesh_x1),
         Nel_tot_x1(Nel_total_x1),
         nsubmesh_x2(numsubmesh_x2),
         Nel_tot_x2(Nel_total_x2),
         isactive(submesh_activity),
         elemoffset_start(elem_offset_start),
         elemoffset_skip(elem_offset_skip),
         nodeoffset(node_offset),
         nodeoffset_start(node_offset_start),
         nodeoffset_skip_bot(node_offset_skip_bot),
         nodeoffset_skip_mid(node_offset_skip_mid),
         nodeoffset_skip_top(node_offset_skip_top),
         is_bdry(isbdry),
         Nel_total(Nel_tot),
         Nnp_total(Nnp_tot),
         host_isactive(h_isactive)
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
    * @brief Constructor for 2D Wrapper structure
    * \param[in] Mesh object in allocated in GPU
    * \param[in] x1-submesh object in GPU
    * \param[in] copy of x1-submesh object in CPU
    * \param[in] x2-submesh object in GPU
    * \param[in] copy of x2-submesh object in CPU
    */
    MBBL(MeshDeviceViewPtr pumi_mesh,
         SubmeshDeviceViewPtr pumi_submesh_x1,
         SubmeshHostViewPtr pumi_host_submesh_x1,
         SubmeshDeviceViewPtr pumi_submesh_x2,
         SubmeshHostViewPtr pumi_host_submesh_x2):
         mesh(pumi_mesh),
         submesh_x1(pumi_submesh_x1),
         host_submesh_x1(pumi_host_submesh_x1),
         submesh_x2(pumi_submesh_x2),
         host_submesh_x2(pumi_host_submesh_x2){};
};

/**
 * @brief Prints the all relevant 1D-mesh details
 *
 * \param[in] mesh object pointer
 * \param[in] x1-submesh object pointer
 */
KOKKOS_INLINE_FUNCTION
void print_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshDeviceViewPtr submesh_x1){
    printf("\n\nPUMI mesh parameter info [X1-Direction] :\n\n");
    printf("\tTotal elements along X1-direction = %d\n\n", pumi_mesh(0).Nel_tot_x1);
    for (int i=0; i<pumi_mesh(0).nsubmesh_x1; i++){
        printf("\tSUBMESH %d  parameters:\n", i+1);
        printf("\n\t submesh-type   = ");
        if (submesh_x1(i)()->meshtype & minBL){
            printf("leftBL\n");
            printf("\t left_t0        = %2.4e \t [m] Cell size of first/leftmost cell in left BL block\n", submesh_x1(i)()->t0);
            printf("\t left_tN        = %2.4e \t [m] Cell size of last/rightmost cell in left BL block\n",
                (submesh_x1(i)()->t0)*pow(submesh_x1(i)()->r,submesh_x1(i)()->Nel-1));
            printf("\t left_T         = %2.4e \t [m] Left boundary layer (left BL) thickness\n", submesh_x1(i)()->length);
            printf("\t left_r         = %2.4e \t Grading ratio in left BL block\n", submesh_x1(i)()->r);
            printf("\t left_Nel       = %d    \t Number of Cells in left BL block\n\n", submesh_x1(i)()->Nel);
            if (submesh_x1(i)()->BL_coords.extent(0)){
                printf("\t %d leftBL node coords stored in a array\n\n", submesh_x1(i)()->BL_coords.extent(0));
            }
        }
        if (submesh_x1(i)()->meshtype & maxBL){
            printf("rightBL\n");
            printf("\t right_t0       = %2.4e \t [m] Cell size of last/rightmost cell in right BL block\n", submesh_x1(i)()->t0);
            printf("\t right_tN       = %2.4e \t [m] Cell size of first/leftmost cell in right BL block\n",
            (submesh_x1(i)()->t0)*pow(submesh_x1(i)()->r,submesh_x1(i)()->Nel-1));
            printf("\t right_T        = %2.4e \t [m] Left boundary layer (right BL) thickness\n", submesh_x1(i)()->length);
            printf("\t right_r        = %2.4e \t Grading ratio in right BL mesh\n", submesh_x1(i)()->r);
            printf("\t right_Nel      = %d    \t Number of Cells in left BL mesh region\n\n", submesh_x1(i)()->Nel);
            if (submesh_x1(i)()->BL_coords.extent(0)){
                printf("\t %d rightBL node coords stored in a array\n\n", submesh_x1(i)()->BL_coords.extent(0));
            }
        }
        if (submesh_x1(i)()->meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_dx1    = %2.4e \t [m] Cell size in uniform block\n", submesh_x1(i)()->t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", submesh_x1(i)()->Nel);
        }
    }
}

/**
 * @brief Prints the all relevant 2D-mesh details
 *
 * \param[in] mesh object pointer
 * \param[in] x1-submesh object pointer
 * \param[in] x2-submesh object pointer
 */
KOKKOS_INLINE_FUNCTION
void print_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshDeviceViewPtr submesh_x1, SubmeshDeviceViewPtr submesh_x2){
    printf("\n\nPUMI mesh parameter info [X1-Direction] :\n\n");
    printf("\tTotal elements along X1-direction = %d\n\n", pumi_mesh(0).Nel_tot_x1);
    for (int i=0; i<pumi_mesh(0).nsubmesh_x1; i++){
        printf("\tSUBMESH %d  parameters:\n", i+1);
        printf("\n\t submesh-type   = ");
        if (submesh_x1(i)()->meshtype & minBL){
            printf("leftBL\n");
            printf("\t left_t0        = %2.4e \t [m] Cell size of first/leftmost cell in left BL block\n", submesh_x1(i)()->t0);
            printf("\t left_tN        = %2.4e \t [m] Cell size of last/rightmost cell in left BL block\n",
                (submesh_x1(i)()->t0)*pow(submesh_x1(i)()->r,submesh_x1(i)()->Nel-1));
            printf("\t left_T         = %2.4e \t [m] Left boundary layer (left BL) thickness\n", submesh_x1(i)()->length);
            printf("\t left_r         = %2.4e \t Grading ratio in left BL block\n", submesh_x1(i)()->r);
            printf("\t left_Nel       = %d    \t Number of Cells in left BL block\n\n", submesh_x1(i)()->Nel);
            if (submesh_x1(i)()->BL_coords.extent(0)){
                printf("\t %d leftBL node coords stored in a array\n\n", submesh_x1(i)()->BL_coords.extent(0));
            }
        }
        if (submesh_x1(i)()->meshtype & maxBL){
            printf("rightBL\n");
            printf("\t right_t0       = %2.4e \t [m] Cell size of last/rightmost cell in right BL block\n", submesh_x1(i)()->t0);
            printf("\t right_tN       = %2.4e \t [m] Cell size of first/leftmost cell in right BL block\n",
            (submesh_x1(i)()->t0)*pow(submesh_x1(i)()->r,submesh_x1(i)()->Nel-1));
            printf("\t right_T        = %2.4e \t [m] Left boundary layer (right BL) thickness\n", submesh_x1(i)()->length);
            printf("\t right_r        = %2.4e \t Grading ratio in right BL mesh\n", submesh_x1(i)()->r);
            printf("\t right_Nel      = %d    \t Number of Cells in left BL mesh region\n\n", submesh_x1(i)()->Nel);
            if (submesh_x1(i)()->BL_coords.extent(0)){
                printf("\t %d rightBL node coords stored in a array\n\n", submesh_x1(i)()->BL_coords.extent(0));
            }
        }
        if (submesh_x1(i)()->meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_dx1    = %2.4e \t [m] Cell size in uniform block\n", submesh_x1(i)()->t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", submesh_x1(i)()->Nel);
        }
    }

    printf("PUMI mesh parameter info [X2-Direction] :\n\n");
    printf("\tTotal elements along X2-direction = %d\n\n", pumi_mesh(0).Nel_tot_x2);
    for (int i=0; i<pumi_mesh(0).nsubmesh_x2; i++){
        printf("\tSUBMESH %d  parameters:\n", i+1);
        printf("\n\t submesh-type   = ");
        if (submesh_x2(i)()->meshtype & minBL){
            printf("bottomBL\n");
            printf("\t bottom_t0      = %2.4e \t [m] Cell size of first/bottom-most cell in left BL block\n", submesh_x2(i)()->t0);
            printf("\t bottom_tN      = %2.4e \t [m] Cell size of last /   top-most cell in left BL block\n",
                (submesh_x2(i)()->t0)*pow(submesh_x2(i)()->r,submesh_x2(i)()->Nel-1));
            printf("\t bottom_T       = %2.4e \t [m] Left boundary layer (left BL) thickness\n", submesh_x2(i)()->length);
            printf("\t bottom_r       = %2.4e \t Grading ratio in left BL block\n", submesh_x2(i)()->r);
            printf("\t bottom_Nel     = %d    \t Number of Cells in left BL block\n\n", submesh_x2(i)()->Nel);
            if (submesh_x2(i)()->BL_coords.extent(0)){
                printf("\t %d bottomBL node coords stored in a array\n\n", submesh_x2(i)()->BL_coords.extent(0));
            }
        }
        if (submesh_x2(i)()->meshtype & maxBL){
            printf("topBL\n");
            printf("\t top_t0         = %2.4e \t [m] Cell size of last /   top-most cell in right BL block\n", submesh_x2(i)()->t0);
            printf("\t top_tN         = %2.4e \t [m] Cell size of first/bottom-most cell in right BL block\n",
            (submesh_x2(i)()->t0)*pow(submesh_x2(i)()->r,submesh_x2(i)()->Nel-1));

            printf("\t top_T          = %2.4e \t [m] Left boundary layer (right BL) thickness\n", submesh_x2(i)()->length);
            printf("\t top_r          = %2.4e \t Grading ratio in right BL mesh\n", submesh_x2(i)()->r);
            printf("\t top_Nel        = %d    \t Number of Cells in left BL mesh region\n\n", submesh_x2(i)()->Nel);
            if (submesh_x2(i)()->BL_coords.extent(0)){
                printf("\t %d topBL node coords stored in a array\n\n", submesh_x2(i)()->BL_coords.extent(0));
            }
        }
        if (submesh_x2(i)()->meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_dx2    = %2.4e \t [m] Cell size in uniform block\n", submesh_x2(i)()->t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", submesh_x2(i)()->Nel);
        }
    }

    printf("PUMI submesh activity info :\n\n");
    for (int jsubmesh=pumi_mesh(0).nsubmesh_x2-1; jsubmesh>=0; jsubmesh--){
        if (jsubmesh != pumi_mesh(0).nsubmesh_x2-1){
            for (int isubmesh=0; isubmesh<pumi_mesh(0).nsubmesh_x1-1; isubmesh++ ){
                printf("_____________");
            }
            printf("__________________\n\n");
        }
        for (int isubmesh=0; isubmesh<pumi_mesh(0).nsubmesh_x1; isubmesh++ ){
            if (isubmesh){
                printf("|");
            }
            if(pumi_mesh(0).isactive(isubmesh,jsubmesh)){
                printf("    ACTIVE    ");
            }
            else{
                printf("    XXXXXX    ");
            }
        }
        printf("\n");
    }
    printf("\n\n");

    printf("Total active elements in 2D Mesh = %d\n",pumi_mesh(0).Nel_total);
    printf("Total active nodes in 2D Mesh = %d\n",pumi_mesh(0).Nnp_total);

}

/**
 * @brief Prints node coordinates of 1D mesh to the terminal and also
 * writes the individual submesh coords to a file as well as
 * full mesh coordinates for all directions
 * \param[in] mesh object pointer
 * \param[in] CPU copy of x1-submesh object pointer
 */
void print_mesh_nodes(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);
    printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
    FILE *mesh_coords_file;
    char mesh_coords_filename[30];
    sprintf(mesh_coords_filename,"X1_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    int inode=0;
    for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
        printf("X1-SUBMESH %d:\n", isubmesh+1 );
        FILE *submesh_coords_file;
        char submesh_coords_filename[30];
        sprintf(submesh_coords_filename,"X1_submesh_%d_coords.dat",isubmesh+1);
        submesh_coords_file = fopen(submesh_coords_filename,"w");
        double icoord = h_submesh_x1[isubmesh].xmin;
        double cell_size = h_submesh_x1[isubmesh].t0;
        double r = h_submesh_x1[isubmesh].r;
        if (h_submesh_x1[isubmesh].meshtype & maxBL){
            cell_size = (h_submesh_x1[isubmesh].t0)*
                        pow(h_submesh_x1[isubmesh].r,h_submesh_x1[isubmesh].Nel-1);
            r = 1.0/r;
        }

        if (isubmesh==0){
            fprintf(mesh_coords_file, "%.16e\n", icoord);
        }
        fprintf(submesh_coords_file, "%.16e\n", icoord);
        printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x1[isubmesh].Nel; icell++){
            inode++;
            icoord += cell_size;
            cell_size *= r;
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        }
        printf("\n\t Coordinates written to file %s\n\n", submesh_coords_filename);
        fclose(submesh_coords_file);
    }
    fclose(mesh_coords_file);
}

/**
 * @brief Prints node coordinates of 2D mesh to the terminal and also
 * writes the individual submesh coords to a file as well as
 * full mesh coordinates for all directions
 * \param[in] mesh object pointer
 * \param[in] CPU copy of x1-submesh object pointer
 * \param[in] CPU copy of x2-submesh object pointer
 */
void print_mesh_nodes(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);

    printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
    FILE *mesh_coords_file;
    char mesh_coords_filename[30];
    sprintf(mesh_coords_filename,"X1_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    int inode=0;
    for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
        printf("X1-SUBMESH %d:\n", isubmesh+1 );
        FILE *submesh_coords_file;
        char submesh_coords_filename[30];
        sprintf(submesh_coords_filename,"X1_submesh_%d_coords.dat",isubmesh+1);
        submesh_coords_file = fopen(submesh_coords_filename,"w");
        double icoord = h_submesh_x1[isubmesh].xmin;
        double cell_size = h_submesh_x1[isubmesh].t0;
        double r = h_submesh_x1[isubmesh].r;
        if (h_submesh_x1[isubmesh].meshtype & maxBL){
            cell_size = (h_submesh_x1[isubmesh].t0)*
                        pow(h_submesh_x1[isubmesh].r,h_submesh_x1[isubmesh].Nel-1);
            r = 1.0/r;
        }

        if (isubmesh==0){
            fprintf(mesh_coords_file, "%.16e\n", icoord);
        }
        fprintf(submesh_coords_file, "%.16e\n", icoord);
        printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x1[isubmesh].Nel; icell++){
            inode++;
            icoord += cell_size;
            cell_size *= r;
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        }
        printf("\n\t Coordinates written to file %s\n\n", submesh_coords_filename);
        fclose(submesh_coords_file);
    }
    fclose(mesh_coords_file);

    sprintf(mesh_coords_filename,"X2_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    inode = 0;
    for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x2; isubmesh++){
        printf("X2-SUBMESH %d:\n", isubmesh+1 );
        FILE *submesh_coords_file;
        char submesh_coords_filename[30];
        sprintf(submesh_coords_filename,"X2_submesh_%d_coords.dat",isubmesh+1);
        submesh_coords_file = fopen(submesh_coords_filename,"w");
        double icoord = h_submesh_x2[isubmesh].xmin;
        double cell_size = h_submesh_x2[isubmesh].t0;
        double r = h_submesh_x2[isubmesh].r;
        if (h_submesh_x2[isubmesh].meshtype & maxBL){
            cell_size = (h_submesh_x2[isubmesh].t0)*
                        pow(h_submesh_x2[isubmesh].r,h_submesh_x2[isubmesh].Nel-1);
            r = 1.0/r;
        }

        if (isubmesh==0){
            fprintf(mesh_coords_file, "%.16e\n", icoord);
        }
        fprintf(submesh_coords_file, "%.16e\n", icoord);
        printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x2[isubmesh].Nel; icell++){
            inode++;
            icoord += cell_size;
            cell_size *= r;
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        }
        printf("\n\t Coordinates written to file %s\n\n", submesh_coords_filename);
        fclose(submesh_coords_file);
    }
    fclose(mesh_coords_file);
}

/**
 * @brief Verifies the computed mesh and submesh parameters
 *
 * \param[in] mesh object pointer
 * \param[in] x1-submesh object pointer
 * \param[out] Boolean value indicating mesh validity (True -- valid mesh, False -- invalid mesh)
 */
KOKKOS_INLINE_FUNCTION
void verify_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshDeviceViewPtr submesh_x1, Kokkos::View<bool*> mesh_verified){
    printf("\n\nNow verifying valdity of pumi mesh parameters for\n");
    int flag = 0;
    for (int isubmesh=0; isubmesh<pumi_mesh(0).nsubmesh_x1; isubmesh++){
        printf("\tX1-SUBMESH %d:\n", isubmesh+1 );
        if (submesh_x1(isubmesh)()->meshtype & minBL){
            if (!(submesh_x1(isubmesh)()->Nel > 0)){
                printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x1(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t left_Nel    -- verified...\n");
            }
            if (!(submesh_x1(isubmesh)()->r >= 1.0)){
                printf("\t\t left_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", submesh_x1(isubmesh)()->r);
                flag++;
            }
            else{
                printf("\t\t left_r      -- verified...\n");
            }
            double min_computed_T = submesh_x1(isubmesh)()->t0*submesh_x1(isubmesh)()->Nel;
            if (!(submesh_x1(isubmesh)()->length > 0.0)){
                printf("\t\t left_T   = %2.4f is not a valid input. It has to be postive\n", submesh_x1(isubmesh)()->length);
                flag++;
            }
            else if (!(submesh_x1(isubmesh)()->length > min_computed_T)){
                printf("\t\t left_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                submesh_x1(isubmesh)()->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t left_T      -- verified...\n");
            }
        }
        if (submesh_x1(isubmesh)()->meshtype & maxBL){
            if (!(submesh_x1(isubmesh)()->Nel > 0)){
                printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x1(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t right_Nel   -- verified...\n");
            }
            if (!(submesh_x1(isubmesh)()->r >= 1.0)){
                printf("\t\t right_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", submesh_x1(isubmesh)()->r);
                flag++;
            }
            else{
                printf("\t\t right_r     -- verified...\n");
            }
            double min_computed_T = submesh_x1(isubmesh)()->t0*submesh_x1(isubmesh)()->Nel;
            if (!(submesh_x1(isubmesh)()->length > 0.0)){
                printf("\t\t right_T   = %2.4f is not a valid input. It has to be postive\n", submesh_x1(isubmesh)()->length);
                flag++;
            }
            else if (!(submesh_x1(isubmesh)()->length > min_computed_T)){
                printf("\t\t right_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                submesh_x1(isubmesh)()->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t right_T     -- verified...\n");
            }
        }
        if (submesh_x1(isubmesh)()->meshtype & uniform){
            if (!(submesh_x1(isubmesh)()->Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x1(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
    }

    if (flag == 0){
        printf("\n\tThe input mesh parameters and the calculated mesh parameters are all valid and verified\n\n");
        mesh_verified(0) = true;
    }
    else{
        printf("\t\nERROR: One or more input/calculated mesh paramater is not valid. Abort\n");
        mesh_verified(0) = false;
    }
}

/**
 * @brief Verifies the computed mesh and submesh parameters
 *
 * \param[in] mesh object pointer
 * \param[in] x1-submesh object pointer
 * \param[in] x2-submesh object pointer
 * \param[out] Boolean value indicating mesh validity (True -- valid mesh, False -- invalid mesh )
 */
KOKKOS_INLINE_FUNCTION
void verify_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshDeviceViewPtr submesh_x1, SubmeshDeviceViewPtr submesh_x2, Kokkos::View<bool*> mesh_verified){
    printf("\n\nNow verifying valdity of pumi mesh parameters for\n");
    int flag = 0;
    for (int isubmesh=0; isubmesh<pumi_mesh(0).nsubmesh_x1; isubmesh++){
        printf("\tX1-SUBMESH %d:\n", isubmesh+1 );
        if (submesh_x1(isubmesh)()->meshtype & minBL){
            if (!(submesh_x1(isubmesh)()->Nel > 0)){
                printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x1(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t left_Nel    -- verified...\n");
            }
            if (!(submesh_x1(isubmesh)()->r >= 1.0)){
                printf("\t\t left_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", submesh_x1(isubmesh)()->r);
                flag++;
            }
            else{
                printf("\t\t left_r      -- verified...\n");
            }
            double min_computed_T = submesh_x1(isubmesh)()->t0*submesh_x1(isubmesh)()->Nel;
            if (!(submesh_x1(isubmesh)()->length > 0.0)){
                printf("\t\t left_T   = %2.4f is not a valid input. It has to be postive\n", submesh_x1(isubmesh)()->length);
                flag++;
            }
            else if (!(submesh_x1(isubmesh)()->length > min_computed_T)){
                printf("\t\t left_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                submesh_x1(isubmesh)()->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t left_T      -- verified...\n");
            }
        }
        if (submesh_x1(isubmesh)()->meshtype & maxBL){
            if (!(submesh_x1(isubmesh)()->Nel > 0)){
                printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x1(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t right_Nel   -- verified...\n");
            }
            if (!(submesh_x1(isubmesh)()->r >= 1.0)){
                printf("\t\t right_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", submesh_x1(isubmesh)()->r);
                flag++;
            }
            else{
                printf("\t\t right_r     -- verified...\n");
            }
            double min_computed_T = submesh_x1(isubmesh)()->t0*submesh_x1(isubmesh)()->Nel;
            if (!(submesh_x1(isubmesh)()->length > 0.0)){
                printf("\t\t right_T   = %2.4f is not a valid input. It has to be postive\n", submesh_x1(isubmesh)()->length);
                flag++;
            }
            else if (!(submesh_x1(isubmesh)()->length > min_computed_T)){
                printf("\t\t right_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                submesh_x1(isubmesh)()->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t right_T     -- verified...\n");
            }
        }
        if (submesh_x1(isubmesh)()->meshtype & uniform){
            if (!(submesh_x1(isubmesh)()->Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x1(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
    }


    for (int isubmesh=0; isubmesh<pumi_mesh(0).nsubmesh_x2; isubmesh++){
        printf("\tX2-SUBMESH %d:\n", isubmesh+1 );
        if (submesh_x2(isubmesh)()->meshtype & minBL){
            if (!(submesh_x2(isubmesh)()->Nel > 0)){
                printf("\t\t bottom_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x2(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t bottom_Nel  -- verified...\n");
            }
            if (!(submesh_x2(isubmesh)()->r >= 1.0)){
                printf("\t\t bottom_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", submesh_x2(isubmesh)()->r);
                flag++;
            }
            else{
                printf("\t\t bottom_r    -- verified...\n");
            }
            double min_computed_T = submesh_x2(isubmesh)()->t0*submesh_x2(isubmesh)()->Nel;
            if (!(submesh_x2(isubmesh)()->length > 0.0)){
                printf("\t\t bottom_T   = %2.4f is not a valid input. It has to be postive\n", submesh_x2(isubmesh)()->length);
                flag++;
            }
            else if (!(submesh_x2(isubmesh)()->length > min_computed_T)){
                printf("\t\t bottom_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                submesh_x2(isubmesh)()->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t bottom_T    -- verified...\n");
            }
        }
        if (submesh_x2(isubmesh)()->meshtype & maxBL){
            if (!(submesh_x2(isubmesh)()->Nel > 0)){
                printf("\t\t top_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x2(isubmesh)()->Nel);
                flag++;
            }
            else{
                printf("\t\t top_Nel     -- verified...\n");
            }
            if (!(submesh_x2(isubmesh)()->r >= 1.0)){
                printf("\t\t top_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", submesh_x2(isubmesh)()->r);
                flag++;
            }
            else{
                printf("\t\t top_r       -- verified...\n");
            }
            double min_computed_T = submesh_x2(isubmesh)()->t0*submesh_x2(isubmesh)()->Nel;
            if (!(submesh_x2(isubmesh)()->length > 0.0)){
                printf("\t\t top_T   = %2.4f is not a valid input. It has to be postive\n", submesh_x2(isubmesh)()->length);
                flag++;
            }
            else if (!(submesh_x2(isubmesh)()->length > min_computed_T)){
                printf("\t\t top_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                submesh_x2(isubmesh)()->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t top_T       -- verified...\n");
            }
        }
        if (submesh_x2(isubmesh)()->meshtype & uniform){
            if (!(submesh_x2(isubmesh)()->Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", submesh_x2(isubmesh)()->Nel);
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
    }


    if (flag == 0){
        printf("\n\tThe input mesh parameters and the calculated mesh parameters are all valid and verified\n\n");
        mesh_verified(0) = true;
    }
    else{
        printf("\t\nERROR: One or more input/calculated mesh paramater is not valid. Abort\n");
        mesh_verified(0) = false;
    }
}

/**
 * @brief Performs necessary submesh parameter calculations
 * and initiates the submesh object and returns the submesh object
 *
 * \param[in] pointer object to pumi inputs structure which containts
 * \param[in] object to structure containing user options for the mesh
 * \param[in] direction along which submesh object is being initialized
 * \param[out] Host-copy of array of submeshes
 * \return Final mesh object pointer
 */
SubmeshDeviceViewPtr submesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, int dir, SubmeshHostViewPtr* hc_submesh){
    int nsubmesh;
    SubmeshDeviceViewPtr submesh;
    SubmeshDeviceViewPtr::HostMirror h_submesh;

    if (dir == x1_dir){
        nsubmesh = pumi_inputs->nsubmesh_x1;
        submesh = SubmeshDeviceViewPtr("x1-submesh-obj",nsubmesh);
        h_submesh = Kokkos::create_mirror_view(submesh);
    }
    else if (dir == x2_dir){
        nsubmesh = pumi_inputs->nsubmesh_x2;
        submesh = SubmeshDeviceViewPtr("x2-submesh-obj",nsubmesh);
        h_submesh = Kokkos::create_mirror_view(submesh);
    }

    SubmeshHostViewPtr submesh_host_copy = new Submesh[nsubmesh];

    if (pumi_options.BL_storage_option == store_BL_coords_OFF){
        pumi_options.BL_storage_option = store_BL_coords_ON;
    }

    double xmin, xmax, xlength, t0, tN, r, r_t0_ratio, logr;
    int Nel;
    unsigned int type;
    int Nel_cumulative = 0;

    for (int isubmesh=0; isubmesh<nsubmesh; isubmesh++){
        DoubleViewPtr BLcoords;
        std::string BLcoordsname;

        if (dir == x1_dir){
            type = get_submesh_type(pumi_inputs->meshtype[isubmesh]);

            xlength = *(pumi_inputs->p1_i_x1 + isubmesh);
            if (isubmesh==0){
                xmin = 0.0;
            }
            else{
                xmin = xmax;
            }
            xmax = xmin + xlength;
            t0 = *(pumi_inputs->p2min_i_x1 + isubmesh);
            tN = *(pumi_inputs->p2max_i_x1 + isubmesh);

            BLcoordsname = "x1-BLcoords-submesh-";
            BLcoordsname += std::to_string(isubmesh);
        }
        else if (dir == x2_dir){
            type = get_submesh_type(pumi_inputs->meshtype[isubmesh+pumi_inputs->nsubmesh_x1]);

            xlength = *(pumi_inputs->p1_i_x2 + isubmesh);
            if (isubmesh==0){
                xmin = 0.0;
            }
            else{
                xmin = xmax;
            }
            xmax = xmin + xlength;
            t0 = *(pumi_inputs->p2min_i_x2 + isubmesh);
            tN = *(pumi_inputs->p2max_i_x2 + isubmesh);

            BLcoordsname = "x2-BLcoords-submesh-";
            BLcoordsname += std::to_string(isubmesh);
        }

        if (isubmesh > 0){
            Nel_cumulative += Nel;
        }

        if (type & uniform){
            r = 1.0;
            Nel = floor(xlength/t0);
            t0 = xlength/Nel;
            r_t0_ratio = (r-1.0)/t0;
            logr = log(r);

            Uniform_Submesh tmp_obj(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords);
            submesh_host_copy[isubmesh] = tmp_obj;
            h_submesh(isubmesh) = copyForDevice<Submesh, Uniform_Submesh> (tmp_obj);
        }
        if (type & minBL){
            r = (xlength-t0)/(xlength-tN);
            if (fabs(r-1.0) < 1e-5){
                Nel = floor(xlength/t0);
                r = 1.0;
            }
            else{
                Nel = 1 + ceil(log(tN/t0)/log(r));
                r = compute_grading_ratio(xlength, t0, Nel);
            }
            r_t0_ratio = (r-1.0)/t0;
            logr = log(r);
            if (pumi_options.BL_storage_option){
                BLcoords = DoubleViewPtr (BLcoordsname, Nel+1);
                DoubleViewPtr::HostMirror h_BLcoords = Kokkos::create_mirror_view(BLcoords);
                h_BLcoords(0) = xmin;
                double cell_size = t0;
                for (int icell=1; icell<Nel; icell++){
                    h_BLcoords(icell) = h_BLcoords(icell-1) + cell_size;
                    cell_size *= r;
                }
                h_BLcoords(Nel) = xmax;

                Kokkos::deep_copy(BLcoords, h_BLcoords);
            }
            MinBL_Submesh tmp_obj(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords);
            submesh_host_copy[isubmesh] = tmp_obj;
            h_submesh(isubmesh) = copyForDevice<Submesh, MinBL_Submesh> (tmp_obj);
        }
        if (type & maxBL){
            r = (xlength-t0)/(xlength-tN);
            if (fabs(r-1.0) < 1e-5){
                Nel = floor(xlength/t0);
                r = 1.0;
            }
            else{
                Nel = 1 + ceil(log(tN/t0)/log(r));
                r = compute_grading_ratio(xlength, t0, Nel);
            }
            r_t0_ratio = (r-1.0)/t0;
            logr = log(r);
            if (pumi_options.BL_storage_option){
                BLcoords = DoubleViewPtr (BLcoordsname, Nel+1);
                DoubleViewPtr::HostMirror h_BLcoords = Kokkos::create_mirror_view(BLcoords);
                h_BLcoords(0) = xmin;
                double cell_size = t0*pow(r,Nel-1);
                for (int icell=1; icell<Nel; icell++){
                    h_BLcoords(icell) = h_BLcoords(icell-1) + cell_size;
                    cell_size /= r;
                }
                h_BLcoords(Nel) = xmax;

                Kokkos::deep_copy(BLcoords, h_BLcoords);
            }
            MaxBL_Submesh tmp_obj(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords);
            submesh_host_copy[isubmesh] = tmp_obj;
            h_submesh(isubmesh) = copyForDevice<Submesh, MaxBL_Submesh> (tmp_obj);
        }

    }

    Kokkos::deep_copy(submesh, h_submesh);
    *hc_submesh = submesh_host_copy;

    return submesh;
}

/**
 * @brief Initiates the mesh and returns the final 1D-mesh object
 *
 * \param[in] pointer object to pumi inputs structure which containts
 * \param[in] x1-submesh object pointer
 * \param[in] Copy of x1-submesh object pointer on CPU
 * \return Final mesh object pointer
 */
MeshDeviceViewPtr mesh_initialize(Mesh_Inputs *pumi_inputs, SubmeshDeviceViewPtr submesh_x1, SubmeshHostViewPtr hc_submesh_x1){
    MeshDeviceViewPtr pumi_mesh("1D-meshobj",1);

    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int Nel_tot_x1;

    Kokkos::View<int*> Nel_total_x1("total-x1-elems",1);
    Kokkos::View<int*>::HostMirror h_Nel_total_x1 = Kokkos::create_mirror_view(Nel_total_x1);

    Kokkos::parallel_for("1D-total-nel-init", 1, KOKKOS_LAMBDA (const int) {
        Nel_total_x1(0) = submesh_x1(nsubmesh_x1-1)()->Nel + submesh_x1(nsubmesh_x1-1)()->Nel_cumulative;
    });
    Kokkos::deep_copy(h_Nel_total_x1,Nel_total_x1);
    Nel_tot_x1 = h_Nel_total_x1(0);

    Kokkos::parallel_for("1D-meshobj-init", 1, KOKKOS_LAMBDA (const int) {
        pumi_mesh(0) = Mesh(nsubmesh_x1, Nel_tot_x1);
        print_mesh_params(pumi_mesh, submesh_x1);
    });

    Kokkos::View<bool*> mesh_verified("mesh-verified",1);
    Kokkos::View<bool*>::HostMirror h_mesh_verified = Kokkos::create_mirror_view(mesh_verified);
    Kokkos::parallel_for("1D-mesh-verifications", 1, KOKKOS_LAMBDA (const int) {
        verify_mesh_params(pumi_mesh, submesh_x1, mesh_verified);
    });
    Kokkos::deep_copy(h_mesh_verified, mesh_verified);

    if (!(h_mesh_verified(0))){
        Kokkos::finalize();
        exit(0);
    }
    else{
        print_mesh_nodes(pumi_mesh, hc_submesh_x1);
    }

    return pumi_mesh;
}

/**
 * @brief Initiates the mesh and returns the final 2D mesh object
 *
 * \param[in] pointer object to pumi inputs structure which containts
 * \param[in] x1-submesh object pointer
 * \param[in] Copy of x1-submesh object pointer on CPU
 * \param[in] x2-submesh object pointer
 * \param[in] Copy of x2-submesh object pointer on CPU
 * \return Final mesh object pointer
 */
MeshDeviceViewPtr mesh_initialize(Mesh_Inputs *pumi_inputs, SubmeshDeviceViewPtr submesh_x1, SubmeshHostViewPtr hc_submesh_x1,
                            SubmeshDeviceViewPtr submesh_x2, SubmeshHostViewPtr hc_submesh_x2){
    MeshDeviceViewPtr pumi_mesh("2D-meshobj",1);

    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int nsubmesh_x2 = pumi_inputs->nsubmesh_x2;
    int Nel_tot_x1, Nel_tot_x2, Nel_total_2D, Nnp_total_2D;

    Kokkos::View<int*> Nel_total_x1("total-x1-elems",1);
    Kokkos::View<int*>::HostMirror h_Nel_total_x1 = Kokkos::create_mirror_view(Nel_total_x1);
    Kokkos::View<int*> Nel_total_x2("total-x2-elems",1);
    Kokkos::View<int*>::HostMirror h_Nel_total_x2 = Kokkos::create_mirror_view(Nel_total_x2);
    Kokkos::parallel_for("2D-total-nel-init", 1, KOKKOS_LAMBDA (const int) {
        Nel_total_x1(0) = submesh_x1(nsubmesh_x1-1)()->Nel + submesh_x1(nsubmesh_x1-1)()->Nel_cumulative;
        Nel_total_x2(0) = submesh_x2(nsubmesh_x2-1)()->Nel + submesh_x2(nsubmesh_x2-1)()->Nel_cumulative;
    });
    Kokkos::deep_copy(h_Nel_total_x1,Nel_total_x1);
    Kokkos::deep_copy(h_Nel_total_x2,Nel_total_x2);
    Nel_tot_x1 = h_Nel_total_x1(0);
    Nel_tot_x2 = h_Nel_total_x2(0);
    Nel_total_2D = Nel_tot_x1*Nel_tot_x2;
    Nnp_total_2D = (Nel_tot_x1+1)*(Nel_tot_x2+1);

    Kokkos::View<bool**> submesh_activity("submesh-isactive",nsubmesh_x1,nsubmesh_x2);
    bool** host_isactive = new bool*[nsubmesh_x1];
    for (int i=0; i<nsubmesh_x1; i++){
        host_isactive[i] = new bool[nsubmesh_x2];
    }
    for (int i=0; i<nsubmesh_x1; i++){
        for (int j=0; j<nsubmesh_x2; j++){
            host_isactive[i][j] = pumi_inputs->isactive[i][j];
        }
    }

    Kokkos::View<bool**>::HostMirror h_submesh_activity = Kokkos::create_mirror_view(submesh_activity);
    for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++ ){
        for (int jsubmesh=0; jsubmesh<nsubmesh_x2; jsubmesh++){
            h_submesh_activity(isubmesh,jsubmesh) = pumi_inputs->isactive[isubmesh][jsubmesh];
        }
    }
    Kokkos::deep_copy(submesh_activity, h_submesh_activity);

    bool is_fullmesh = true;
    for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++ ){
        for (int jsubmesh=0; jsubmesh<nsubmesh_x2; jsubmesh++){
            if (!(h_submesh_activity(isubmesh, jsubmesh))){
                is_fullmesh = false;
            }
        }
    }

    Kokkos::View<int**> elemoffset_start("elemoffset_start", nsubmesh_x1, nsubmesh_x2);
    Kokkos::View<int*> elemoffset_skip("elemoffset_skip", nsubmesh_x2);
    Kokkos::View<int**> nodeoffset("nodeoffset", nsubmesh_x1, Nel_tot_x2+1);
    Kokkos::View<int**> nodeoffset_start("nodeoffset_start", nsubmesh_x1, nsubmesh_x2);
    Kokkos::View<int**> nodeoffset_skip_bot("nodeoffset_skip_bot", nsubmesh_x1, nsubmesh_x2);
    Kokkos::View<int**> nodeoffset_skip_mid("nodeoffset_skip_mid", nsubmesh_x1, nsubmesh_x2);
    Kokkos::View<int**> nodeoffset_skip_top("nodeoffset_skip_top", nsubmesh_x1, nsubmesh_x2);


    if (!is_fullmesh){
        Kokkos::View<int**>::HostMirror h_elemoffset_start = Kokkos::create_mirror_view(elemoffset_start);
        Kokkos::View<int*>::HostMirror h_elemoffset_skip = Kokkos::create_mirror_view(elemoffset_skip);
        Kokkos::View<int**>::HostMirror h_nodeoffset = Kokkos::create_mirror_view(nodeoffset);
        Kokkos::View<int**>::HostMirror h_nodeoffset_start = Kokkos::create_mirror_view(nodeoffset_start);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_bot = Kokkos::create_mirror_view(nodeoffset_skip_bot);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_mid = Kokkos::create_mirror_view(nodeoffset_skip_mid);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_top = Kokkos::create_mirror_view(nodeoffset_skip_top);

        // elemoffsets
        int elemstart_init, elemskip;
        int elemstart = 0;

        for (int jsubmesh=0; jsubmesh<nsubmesh_x2; jsubmesh++){
            elemstart_init = elemstart;

            for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
                if(h_submesh_activity(isubmesh, jsubmesh)){
                    h_elemoffset_start(isubmesh,jsubmesh) = elemstart;
                }
                else{
                    elemstart +=  hc_submesh_x1[isubmesh].Nel;
                    h_elemoffset_start(isubmesh,jsubmesh) = -1;
                }
            }
            elemskip = elemstart-elemstart_init;
            h_elemoffset_skip(jsubmesh) = elemskip;
            elemstart = elemstart_init + elemskip*hc_submesh_x2[jsubmesh].Nel;
        }
        Nel_total_2D -= elemstart;
        // printf("Full Mesh = %d\n", Nel_tot_x1*Nel_tot_x2);
        // printf("Actual Mesh = %d\n", Nel_total_2D);

        // nodeoffsets
        int jnp;
        int jsubmesh = 0;
        int nodestart = 0;
        for (jnp=0; jnp<hc_submesh_x2[jsubmesh].Nel; jnp++){
            int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;
            for (int isubmesh=0; isubmesh<nsubmesh_x2; isubmesh++){
                if(h_submesh_activity(isubmesh, jsubmesh)){
                    h_nodeoffset(isubmesh,Jnp) = nodestart;
                }
                else{
                    h_nodeoffset(isubmesh,Jnp) = -1;
                    if (isubmesh==0){
                        if (h_submesh_activity(isubmesh+1,jsubmesh)){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                        }
                    }
                    else if (isubmesh == nsubmesh_x1-1){
                        nodestart += hc_submesh_x1[isubmesh].Nel;
                    }
                    else{
                        if (h_submesh_activity(isubmesh+1,jsubmesh)){
                            nodestart += (hc_submesh_x1[isubmesh].Nel-1);
                        }
                        else{
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                    }
                }
            }
        }

        jnp = hc_submesh_x2[jsubmesh].Nel;
        for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
            if (jsubmesh == nsubmesh_x2-1){
                int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;
                if (h_submesh_activity(isubmesh, jsubmesh)){
                    h_nodeoffset(isubmesh,Jnp) = nodestart;
                }
                else{
                    h_nodeoffset(isubmesh,Jnp) = -1;
                    if (isubmesh==0){
                        if (h_submesh_activity(isubmesh+1,jsubmesh)){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                        }
                    }
                    else if (isubmesh==nsubmesh_x1-1){
                        nodestart += hc_submesh_x1[isubmesh].Nel;
                    }
                    else{
                        if (h_submesh_activity(isubmesh+1,jsubmesh)){
                            nodestart += (hc_submesh_x1[isubmesh].Nel-1);
                        }
                        else{
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                    }
                }
            }
        }

        for (jsubmesh=1; jsubmesh<nsubmesh_x2; jsubmesh++){
            jnp = 0;
            int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;

            for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
                if (h_submesh_activity(isubmesh,jsubmesh) || h_submesh_activity(isubmesh,jsubmesh-1)){
                    nodestart += 0;
                    h_nodeoffset(isubmesh,Jnp) = nodestart;
                }
                else{
                    h_nodeoffset(isubmesh,Jnp) = -1;
                    if (isubmesh==0){
                        if (h_submesh_activity(isubmesh+1,jsubmesh) || h_submesh_activity(isubmesh+1,jsubmesh-1)){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                        }
                    }
                    else if (isubmesh==nsubmesh_x1-1){
                        nodestart += (hc_submesh_x1[isubmesh].Nel);
                    }
                    else{
                        if (h_submesh_activity(isubmesh+1,jsubmesh)){
                            nodestart += (hc_submesh_x1[isubmesh].Nel-1);
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh].Nel);
                        }
                    }
                }
            }

            for (jnp=1; jnp<hc_submesh_x2[jsubmesh].Nel; jnp++){
                int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;

                for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        h_nodeoffset(isubmesh,Jnp) = nodestart;
                    }
                    else{
                        h_nodeoffset(isubmesh,Jnp) = -1;
                        if (isubmesh==0){
                            if (h_submesh_activity(isubmesh+1,jsubmesh)){
                                nodestart += hc_submesh_x1[isubmesh].Nel;
                            }
                            else{
                                nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                            }
                        }
                        else if (isubmesh == nsubmesh_x1-1){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            if (h_submesh_activity(isubmesh+1,jsubmesh)){
                                nodestart += (hc_submesh_x1[isubmesh].Nel-1);
                            }
                            else{
                                nodestart += hc_submesh_x1[isubmesh].Nel;
                            }
                        }
                    }
                }
            }

            if (jsubmesh==nsubmesh_x2-1){
                jnp = hc_submesh_x2[jsubmesh].Nel;
                int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;

                for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        h_nodeoffset(isubmesh,Jnp) = nodestart;
                    }
                    else{
                        h_nodeoffset(isubmesh,Jnp) = -1;
                        if (isubmesh==0){
                            if (h_submesh_activity(isubmesh+1,jsubmesh)){
                                nodestart += hc_submesh_x1[isubmesh].Nel;
                            }
                            else{
                                nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                            }
                        }
                        else if (isubmesh==nsubmesh_x1-1){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            if (h_submesh_activity(isubmesh+1,jsubmesh)){
                                nodestart += (hc_submesh_x1[isubmesh].Nel-1);
                            }
                            else{
                                nodestart += hc_submesh_x1[isubmesh].Nel;
                            }
                        }
                    }
                }
            }

        }

        Nnp_total_2D -= nodestart;

        for (int jsubmesh=0; jsubmesh<nsubmesh_x2; jsubmesh++){
            int Jnp = hc_submesh_x2[jsubmesh].Nel_cumulative;
            for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
                if (h_submesh_activity(isubmesh,jsubmesh)){
                    h_nodeoffset_start(isubmesh,jsubmesh) = h_nodeoffset(isubmesh,Jnp);
                    h_nodeoffset_skip_bot(isubmesh,jsubmesh) = h_nodeoffset(isubmesh,Jnp+1)-h_nodeoffset(isubmesh,Jnp);
                    h_nodeoffset_skip_mid(isubmesh,jsubmesh) = h_nodeoffset(isubmesh,Jnp+2)-h_nodeoffset(isubmesh,Jnp+1);
                    h_nodeoffset_skip_top(isubmesh,jsubmesh) = h_nodeoffset(isubmesh,Jnp+hc_submesh_x2[jsubmesh].Nel) -
                                                                h_nodeoffset(isubmesh,Jnp+hc_submesh_x2[jsubmesh].Nel-1);
                }
                else{
                    h_nodeoffset_start(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_bot(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_mid(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_top(isubmesh,jsubmesh) = -1;
                }
            }
        }
        // printf("Full Mesh = %d\n", (Nel_tot_x1+1)*(Nel_tot_x2+1));
        // printf("Actual Mesh = %d\n", Nnp_total_2D);
        Kokkos::deep_copy(elemoffset_skip, h_elemoffset_skip);
        Kokkos::deep_copy(elemoffset_start, h_elemoffset_start);
        Kokkos::deep_copy(nodeoffset, h_nodeoffset);
        Kokkos::deep_copy(nodeoffset_start, h_nodeoffset_start);
        Kokkos::deep_copy(nodeoffset_skip_bot, h_nodeoffset_skip_bot);
        Kokkos::deep_copy(nodeoffset_skip_mid, h_nodeoffset_skip_mid);
        Kokkos::deep_copy(nodeoffset_skip_top, h_nodeoffset_skip_top);

    }
    else{
        Kokkos::View<int**>::HostMirror h_elemoffset_start = Kokkos::create_mirror_view(elemoffset_start);
        Kokkos::View<int*>::HostMirror h_elemoffset_skip = Kokkos::create_mirror_view(elemoffset_skip);
        Kokkos::View<int**>::HostMirror h_nodeoffset = Kokkos::create_mirror_view(nodeoffset);
        Kokkos::View<int**>::HostMirror h_nodeoffset_start = Kokkos::create_mirror_view(nodeoffset_start);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_bot = Kokkos::create_mirror_view(nodeoffset_skip_bot);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_mid = Kokkos::create_mirror_view(nodeoffset_skip_mid);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_top = Kokkos::create_mirror_view(nodeoffset_skip_top);

        for (int jsubmesh=0; jsubmesh<nsubmesh_x2; jsubmesh++){
            h_elemoffset_skip(jsubmesh) = 0;
            for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
                h_elemoffset_start(isubmesh,jsubmesh) = 0;
                h_nodeoffset_start(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_bot(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_mid(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_top(isubmesh,jsubmesh) = 0;
            }
        }

        Kokkos::deep_copy(elemoffset_skip, h_elemoffset_skip);
        Kokkos::deep_copy(elemoffset_start, h_elemoffset_start);
        Kokkos::deep_copy(nodeoffset, h_nodeoffset);
        Kokkos::deep_copy(nodeoffset_start, h_nodeoffset_start);
        Kokkos::deep_copy(nodeoffset_skip_bot, h_nodeoffset_skip_bot);
        Kokkos::deep_copy(nodeoffset_skip_mid, h_nodeoffset_skip_mid);
        Kokkos::deep_copy(nodeoffset_skip_top, h_nodeoffset_skip_top);
    }

    Kokkos::View<bool*> is_bdry("is_bdry", 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2);
    Kokkos::View<bool*>::HostMirror h_is_bdry = Kokkos::create_mirror_view(is_bdry);

    for (int jsubmesh=0; jsubmesh<nsubmesh_x2; jsubmesh++){
        for (int isubmesh=0; isubmesh<nsubmesh_x1; isubmesh++){
            if (jsubmesh==0){
                if (h_submesh_activity(isubmesh,jsubmesh)){
                    h_is_bdry(isubmesh) = true;
                    if (isubmesh==0){
                        h_is_bdry(nsubmesh_x1) = true;
                    }
                    if (isubmesh==nsubmesh_x1-1){
                        h_is_bdry(isubmesh+nsubmesh_x1+1) = true;
                    }
                }
                else{
                    h_is_bdry(isubmesh) = false;
                    if (isubmesh==0){
                        h_is_bdry(nsubmesh_x1) = false;
                    }
                    if (isubmesh==nsubmesh_x1-1){
                        h_is_bdry(isubmesh+nsubmesh_x1+1) = false;
                    }
                }

                if (isubmesh>0){
                    int sum = h_submesh_activity(isubmesh-1,jsubmesh)+h_submesh_activity(isubmesh,jsubmesh);
                    if (sum == 1){
                        h_is_bdry(isubmesh+nsubmesh_x1) = true;
                    }
                    else{
                        h_is_bdry(isubmesh+nsubmesh_x1) = false;
                    }
                }
            }
            else{

                int sum = h_submesh_activity(isubmesh,jsubmesh-1)+h_submesh_activity(isubmesh,jsubmesh);
                if (sum == 1){
                    h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+isubmesh) = true;
                }
                else{
                    h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+isubmesh) = false;
                }

                if (isubmesh==0){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+nsubmesh_x1) = true;
                    }
                    else{
                        h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+nsubmesh_x1) = false;
                    }
                }
                else{
                    sum = h_submesh_activity(isubmesh-1,jsubmesh)+h_submesh_activity(isubmesh,jsubmesh);
                    if (sum == 1){
                        h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh) = true;
                    }
                    else{
                        h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh) = false;
                    }
                }

                if (isubmesh==nsubmesh_x1-1){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+isubmesh+nsubmesh_x1+1) = true;
                    }
                    else{
                        h_is_bdry(jsubmesh*(2*nsubmesh_x1+1)+isubmesh+nsubmesh_x1+1) = false;
                    }
                }
            }

            if (jsubmesh == nsubmesh_x2-1){
                if (h_submesh_activity(isubmesh,jsubmesh)){
                    h_is_bdry((jsubmesh+1)*(2*nsubmesh_x1+1)+isubmesh) = true;
                }
                else{
                    h_is_bdry((jsubmesh+1)*(2*nsubmesh_x1+1)+isubmesh) = false;
                }
            }
        }
    }

    Kokkos::deep_copy(is_bdry, h_is_bdry);

    // for (int i=0; i<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2; i++){
    //     if (h_is_bdry(i)){
    //         printf("Edge - %3d is boundary edge\n", i);
    //     }
    // }

    Kokkos::parallel_for("2D-meshobj-init", 1, KOKKOS_LAMBDA (const int) {
        pumi_mesh(0) = Mesh(nsubmesh_x1, Nel_tot_x1, nsubmesh_x2, Nel_tot_x2,
                        submesh_activity, elemoffset_start, elemoffset_skip, nodeoffset,
                        nodeoffset_start, nodeoffset_skip_bot, nodeoffset_skip_mid, nodeoffset_skip_top,
                        is_bdry, Nel_total_2D, Nnp_total_2D, host_isactive);
        print_mesh_params(pumi_mesh, submesh_x1, submesh_x2);
    });

    Kokkos::View<bool*> mesh_verified("mesh-verified",1);
    Kokkos::View<bool*>::HostMirror h_mesh_verified = Kokkos::create_mirror_view(mesh_verified);
    Kokkos::parallel_for("2D-mesh-verifications", 1, KOKKOS_LAMBDA (const int) {
        verify_mesh_params(pumi_mesh, submesh_x1, submesh_x2, mesh_verified);
    });
    Kokkos::deep_copy(h_mesh_verified, mesh_verified);

    if (!(h_mesh_verified(0))){
        Kokkos::finalize();
        exit(0);
    }
    else{
        print_mesh_nodes(pumi_mesh, hc_submesh_x1, hc_submesh_x2);
    }

    return pumi_mesh;

}


// KOKKOS_INLINE_FUNCTION
// void print_BL_coords(MeshDeviceViewPtr pumi_mesh, SubmeshDeviceViewPtr submesh_x1, SubmeshDeviceViewPtr submesh_x2){
//     for (int isubmesh=0; isubmesh < pumi_mesh(0).nsubmesh_x1; isubmesh++){
//         if (submesh_x1(isubmesh)()->meshtype & minBL){
//             printf("\n");
//             for (int i=0; i<=submesh_x1(isubmesh)()->Nel; i++){
//                 printf("%2.4f\n",submesh_x1(isubmesh)()->BL_coords(i));
//             }
//         }
//         if (submesh_x1(isubmesh)()->meshtype & maxBL){
//             printf("\n");
//             for (int i=0; i<=submesh_x1(isubmesh)()->Nel; i++){
//                 printf("%2.4f\n",submesh_x1(isubmesh)()->BL_coords(i));
//             }
//         }
//     }
//
//     for (int isubmesh=0; isubmesh < pumi_mesh(0).nsubmesh_x2; isubmesh++){
//         if (submesh_x2(isubmesh)()->meshtype & minBL){
//             printf("\n");
//             for (int i=0; i<=submesh_x2(isubmesh)()->Nel; i++){
//                 printf("%2.4f\n",submesh_x2(isubmesh)()->BL_coords(i));
//             }
//         }
//         if (submesh_x2(isubmesh)()->meshtype & maxBL){
//             printf("\n");
//             for (int i=0; i<=submesh_x2(isubmesh)()->Nel; i++){
//                 printf("%2.4f\n",submesh_x2(isubmesh)()->BL_coords(i));
//             }
//         }
//     }
// }

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
    for (isubmesh=1; isubmesh<nsubmesh; isubmesh++){
     if (q < (pumi_obj.submesh_x1(isubmesh)()->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh-1;
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
    for (isubmesh=1; isubmesh<nsubmesh; isubmesh++){
     if (q < (pumi_obj.submesh_x2(isubmesh)()->xmin)){
         *submeshID = isubmesh-1;
         submesh_located++;
         break;
     }
    }
    if (!(submesh_located)){
     *submeshID = nsubmesh-1;
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
        nodeoffset +=  pumi_obj.mesh(0).nodeoffset_skip_top(isubmesh,jsubmesh);
    }
    return nodeID-nodeoffset;
}

KOKKOS_INLINE_FUNCTION
void push_particle(MBBL pumi_obj, double q1, double q2, double dq1, double dq2,
                    int isubmesh, int jsubmesh, int icell, int jcell,
                    int *new_isubmesh, int *new_jsubmesh, int *new_icell, int *new_jcell,
                    bool *in_domain, int *bdry_hit){

    bool crossed_right = false;
    bool crossed_left = false;
    bool crossed_up = false;
    bool crossed_down = false;

    int num_x1_blocks_crossed=0;
    int num_x2_blocks_crossed=0;

    double q1_new = q1 + dq1;
    double q2_new = q2 + dq2;
    int Nx = pumi_obj.mesh(0).nsubmesh_x1;
    int Ny = pumi_obj.mesh(0).nsubmesh_x2;

    *new_isubmesh = isubmesh;
    *new_icell = icell;
    while(q1_new < (pumi_obj.submesh_x1(*new_isubmesh)()->xmin)){
        *new_isubmesh -= 1;
        crossed_left = true;
        num_x1_blocks_crossed++;
        if (*new_isubmesh == -1){
            *new_icell = -1;
            *in_domain = false;
            break;
        }
        *new_icell = pumi_obj.submesh_x1(*new_isubmesh)()->Nel-1;
    }

    while(q1_new > (pumi_obj.submesh_x1(*new_isubmesh)()->xmax)){
        *new_isubmesh += 1;
        crossed_right = true;
        num_x1_blocks_crossed++;
        if (*new_isubmesh == Nx){
            *new_icell = -1;
            *in_domain = false;
            break;
        }
        *new_icell = 0;
    }

    *new_jsubmesh = jsubmesh;
    *new_jcell = jcell;
    while(q2_new < (pumi_obj.submesh_x2(*new_jsubmesh)()->xmin)){
        *new_jsubmesh -= 1;
        crossed_down = true;
        num_x2_blocks_crossed++;
        if (*new_jsubmesh == -1){
            *new_jcell = -1;
            *in_domain = false;
            break;
        }
        *new_jcell = pumi_obj.submesh_x2(*new_jsubmesh)()->Nel-1;
    }

    while(q2_new > (pumi_obj.submesh_x2(*new_jsubmesh)()->xmax)){
        *new_jsubmesh += 1;
        crossed_up = true;
        num_x2_blocks_crossed++;
        if (*new_jsubmesh == Ny){
            *new_jcell = -1;
            *in_domain = false;
            break;
        }
        *new_jcell = 0;
    }

    if (!crossed_right & !crossed_left & !crossed_down & !crossed_up){
        *in_domain = true;
        *new_icell = pumi_obj.submesh_x1(*new_isubmesh)()->update_cell(q1_new, *new_icell);
        *new_jcell = pumi_obj.submesh_x2(*new_jsubmesh)()->update_cell(q2_new, *new_jcell);
        *bdry_hit = -1;
        return;
    }

    if (crossed_right & !crossed_up & !crossed_up){
        *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh+Nx+1;
    }

    if (crossed_left & !crossed_up & !crossed_up){
        *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh+Nx;
    }

    if (crossed_up & !crossed_right & !crossed_left){
        *bdry_hit = (jsubmesh+1)*(2*Nx+1)+isubmesh;
    }

    if (crossed_down & !crossed_right & !crossed_left){
        *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh;
    }

    double del_x1, del_x2;
    if (crossed_right & crossed_down){
        del_x1 = (pumi_obj.submesh_x1(isubmesh)()->xmax - q1);
        del_x2 = (q2 - pumi_obj.submesh_x2(jsubmesh)()->xmin);

        if (del_x2/del_x1 > fabs(dq2/dq1)){
            *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh+Nx+1;
        }
        else{
            *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh;
        }

    }

    if (crossed_right & crossed_up){
        del_x1 = (pumi_obj.submesh_x1(isubmesh)()->xmax - q1);
        del_x2 = (pumi_obj.submesh_x2(jsubmesh)()->xmax - q2);

        if (del_x2/del_x1 > fabs(dq2/dq1)){
            *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh+Nx+1;
        }
        else{
            *bdry_hit = (jsubmesh+1)*(2*Nx+1)+isubmesh;
        }
    }

    if (crossed_left & crossed_down){
        del_x1 = (q1 - pumi_obj.submesh_x1(isubmesh)()->xmin);
        del_x2 = (q2 - pumi_obj.submesh_x2(jsubmesh)()->xmin);

        if (del_x2/del_x1 > fabs(dq2/dq1)){
            *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh+Nx;
        }
        else{
            *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh;
        }

    }

    if (crossed_left & crossed_up){
        del_x1 = (q1 - pumi_obj.submesh_x1(isubmesh)()->xmin);
        del_x2 = (pumi_obj.submesh_x2(jsubmesh)()->xmax - q2);

        if (del_x2/del_x1 > fabs(dq2/dq1)){
            *bdry_hit = jsubmesh*(2*Nx+1)+isubmesh+Nx;
        }
        else{
            *bdry_hit = (jsubmesh+1)*(2*Nx+1)+isubmesh;
        }
    }

    if (pumi_obj.mesh(0).is_bdry(*bdry_hit)){
        *in_domain = false;
        return;
    }
    else {
        *in_domain = true;
        *bdry_hit = -1;
        *new_icell = pumi_obj.submesh_x1(*new_isubmesh)()->update_cell(q1_new, *new_icell);
        *new_jcell = pumi_obj.submesh_x2(*new_jsubmesh)()->update_cell(q2_new, *new_jcell);
        return;
    }



}

///////Field-related data structures and routines //////////////////////////////////////////////

/*!
* \brief enum type for element index offsets
*
* To be used while querying element size using
* either element ID or associated node ID
*/
enum elemsize_index_offset{
    elem_input_offset       =  0, //!< zero offset for direct element ID input
    elem_on_max_side_offset =  0, //!< offset for node ID input and querying element to the max side
    elem_on_min_side_offset = -1, //!< offset for node ID input and querying element to the min side
};

/**
 * @brief Returns the grading ratio along a queried direction about a queried node
 *
 * \param[in] Object of the wrapper mesh structure
 * \param[in] direction along which grading ratio is needed
 * \param[in] global directional node ID
 * \return grading ratio about the node in the requested direction
 */
double return_gradingratio(MBBL pumi_obj, int dir, int node){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);
    SubmeshHostViewPtr h_submesh;
    int nsubmesh, Nel_total;
    if (dir == x1_dir){
        nsubmesh = h_pumi_mesh(0).nsubmesh_x1;
        Nel_total = h_pumi_mesh(0).Nel_tot_x1;
        h_submesh = pumi_obj.host_submesh_x1;
    }
    else{
        nsubmesh = h_pumi_mesh(0).nsubmesh_x2;
        Nel_total = h_pumi_mesh(0).Nel_tot_x2;
        h_submesh = pumi_obj.host_submesh_x2;
    }
    if (node == 0 || node == Nel_total){
        printf("Grading ratio not defined for the first and last node of the domain -- Terminating \n");
        exit(0);
    }
    else{
        for (int isubmesh=0; isubmesh<nsubmesh; isubmesh++){
            int submesh_min_node =  h_submesh[isubmesh].Nel_cumulative;
            int submesh_max_node = submesh_min_node + h_submesh[isubmesh].Nel;
            if (node > submesh_min_node && node < submesh_max_node){
                if (h_submesh[isubmesh].meshtype & maxBL){
                    return 1.0/h_submesh[isubmesh].r;
                }
                else{
                    return h_submesh[isubmesh].r;
                }
            }
            else if (node == submesh_min_node){
                double max_elem;
                double min_elem;
                if (h_submesh[isubmesh-1].meshtype & minBL){
                    min_elem = h_submesh[isubmesh-1].t0 * pow(h_submesh[isubmesh-1].r, h_submesh[isubmesh-1].Nel-1);
                }
                else{
                    min_elem = h_submesh[isubmesh-1].t0;
                }

                if (h_submesh[isubmesh].meshtype & maxBL){
                    max_elem = h_submesh[isubmesh].t0 * pow(h_submesh[isubmesh].r, h_submesh[isubmesh].Nel-1);
                }
                else{
                    max_elem = h_submesh[isubmesh].t0;
                }
                return max_elem/min_elem;
            }
        }
    }
    return -999.0;
}

/**
 * @brief Returns the element size along a queried direction for a queried node/element
 *
 * \param[in] Object of the wrapper mesh structure
 * \param[in] direction along which element-size is needed
 * \param[in] global directional node ID or element ID index
 * \param[in] relevant offset enum based on provied index
 * \return element size for the queried element in the queried direction
 */
double return_elemsize(MBBL pumi_obj, int dir, int index, int offset){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);
    SubmeshHostViewPtr h_submesh;
    int nsubmesh, Nel_total, elem;
    if (dir == x1_dir){
        nsubmesh = h_pumi_mesh(0).nsubmesh_x1;
        Nel_total = h_pumi_mesh(0).Nel_tot_x1;
        h_submesh = pumi_obj.host_submesh_x1;
    }
    else{
        nsubmesh = h_pumi_mesh(0).nsubmesh_x2;
        Nel_total = h_pumi_mesh(0).Nel_tot_x2;
        h_submesh = pumi_obj.host_submesh_x2;
    }

    elem = index + offset;
    if (elem >= Nel_total){
        elem = Nel_total;
    }
    else if (elem < 0){
        elem = 0;
    }

    for (int isubmesh = 0; isubmesh < nsubmesh; isubmesh++){
        int submesh_min_elem = h_submesh[isubmesh].Nel_cumulative;
        int submesh_max_elem = submesh_min_elem + h_submesh[isubmesh].Nel-1;

        if (elem >= submesh_min_elem && elem <= submesh_max_elem){
            if (h_submesh[isubmesh].meshtype & uniform){
                return h_submesh[isubmesh].t0;
            }
            if (h_submesh[isubmesh].meshtype & minBL){
                int local_cell = elem - submesh_min_elem;
                return (h_submesh[isubmesh].t0 * pow(h_submesh[isubmesh].r , local_cell));
            }
            if (h_submesh[isubmesh].meshtype & maxBL){
                int local_cell = h_submesh[isubmesh].Nel - (elem - submesh_min_elem) - 1;
                return (h_submesh[isubmesh].t0 * pow(h_submesh[isubmesh].r , local_cell));
            }
        }
    }
    return -999.0;
}

/**
 * @brief Returns the element size along a queried direction for a queried node/element
 *
 * \param[in] Object of the wrapper mesh structure
 * \param[in] global node IDs along x1-direction
 * \return covolume for the requested node
 */
double return_covolume(MBBL pumi_obj, int inode_x1){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);

    double covolume, dx1_min, dx1_max;

    if (inode_x1 == 0){
        dx1_min = 0.0;
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
    }
    else if (inode_x1 == h_pumi_mesh(0).Nel_tot_x1){
        dx1_max = 0.0;
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }
    else{
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }

    covolume = (dx1_min+dx1_max)/2.0;
    return covolume;
}


/**
 * @brief Returns the element size along a queried direction for a queried node/element
 *
 * \param[in] Object of the wrapper mesh structure
 * \param[in] global node IDs along x1-direction
 * \param[in] global node IDs along x2-direction
 * \return covolume for the requested node
 */
double return_covolume(MBBL pumi_obj, int inode_x1, int inode_x2){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);

    double covolume, dx1_min, dx1_max, dx2_min, dx2_max;

    if (inode_x1 == 0){
        dx1_min = 0.0;
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
    }
    else if (inode_x1 == h_pumi_mesh(0).Nel_tot_x1){
        dx1_max = 0.0;
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }
    else{
        dx1_max = return_elemsize(pumi_obj, x1_dir, elem_on_max_side_offset, inode_x1);
        dx1_min = return_elemsize(pumi_obj, x1_dir, elem_on_min_side_offset, inode_x1);
    }

    if (inode_x2 == 0){
        dx2_min = 0.0;
        dx2_max = return_elemsize(pumi_obj, x2_dir, elem_on_max_side_offset, inode_x2);
    }
    else if (inode_x2 == h_pumi_mesh(0).Nel_tot_x2){
        dx2_max = 0.0;
        dx2_min = return_elemsize(pumi_obj, x2_dir, elem_on_min_side_offset, inode_x2);
    }
    else{
        dx2_max = return_elemsize(pumi_obj, x2_dir, elem_on_max_side_offset, inode_x2);
        dx2_min = return_elemsize(pumi_obj, x2_dir, elem_on_min_side_offset, inode_x2);
    }

    covolume = (dx1_min*dx2_min + dx1_max*dx2_min + dx1_min*dx2_max + dx1_max*dx2_max)/4.0;
    return covolume;
}

void where_is_node(MBBL pumi_obj, int knode_x1, int knode_x2, bool* on_bdry, bool* in_domain, int* bdry_tag, int* bdry_dim){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);

    int isubmesh, jsubmesh, inp, jnp;
    bool left_edge, right_edge, bottom_edge, top_edge;

    for (isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
        int submesh_min_node = pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
        int submesh_max_node = pumi_obj.host_submesh_x1[isubmesh].Nel + submesh_min_node;
        left_edge =  false;
        right_edge = false;
        if (knode_x1 >= submesh_min_node && knode_x1 <= submesh_max_node){
            inp = knode_x1 - submesh_min_node;
            if (inp == 0){
                left_edge = true;
            }
            if (inp == pumi_obj.host_submesh_x1[isubmesh].Nel){
                right_edge = true;
            }
            break;
        }
    }

    for (jsubmesh=0; jsubmesh<h_pumi_mesh(0).nsubmesh_x2; jsubmesh++){
        int submesh_min_node = pumi_obj.host_submesh_x2[jsubmesh].Nel_cumulative;
        int submesh_max_node = pumi_obj.host_submesh_x2[jsubmesh].Nel + submesh_min_node;
        bottom_edge =  false;
        top_edge = false;
        if (knode_x2 >= submesh_min_node && knode_x2 <= submesh_max_node){
            jnp = knode_x2 - submesh_min_node;
            if (jnp == 0){
                bottom_edge = true;
            }
            if (jnp == pumi_obj.host_submesh_x2[jsubmesh].Nel){
                top_edge = true;
            }
            break;
        }
    }

    if (!left_edge && !right_edge && !bottom_edge && !top_edge){
        *on_bdry = false;

        if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
            *in_domain = true;
            *bdry_tag = -1;
            *bdry_dim = -1;
            return;
        }
        else{
            *in_domain = false;
            *bdry_tag = -999;
            *bdry_dim = -1;
            return;
        }
    }

    if (left_edge & !top_edge & !bottom_edge){

        if (isubmesh==0){
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = jsubmesh*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + h_pumi_mesh(0).nsubmesh_x1;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh]){
                *in_domain = true;
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] + h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh] == 1){
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + h_pumi_mesh(0).nsubmesh_x1;
                    *bdry_dim = 1;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (left_edge & top_edge){
        if (jsubmesh==h_pumi_mesh(0).nsubmesh_x2-1){
            if (isubmesh==0){
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==0){
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh+1] |
                    h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh] + h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh+1] +
                        h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1] + h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }

    if (top_edge & !left_edge & !right_edge){

        if (jsubmesh==h_pumi_mesh(0).nsubmesh_x2-1){
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = (jsubmesh+1)*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1]){
                *in_domain = true;
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] + h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1] == 1){
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh+1)*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 1;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (top_edge & right_edge){
        if (jsubmesh==h_pumi_mesh(0).nsubmesh_x2-1){
            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1] | h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh+1] |
                    h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh+1] + h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh+1] +
                            h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh] + h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = (jsubmesh+1)*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }

    if (right_edge & !top_edge & !bottom_edge){

        if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = jsubmesh*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + h_pumi_mesh(0).nsubmesh_x1 + 1;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh]){
                *in_domain = true;
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] + h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh] == 1){
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + h_pumi_mesh(0).nsubmesh_x1 + 1;
                    *bdry_dim = 1;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (right_edge & bottom_edge){
        if (jsubmesh==0){
            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1] | h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh] |
                    h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh-1] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1] + h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh] +
                            h_pumi_mesh(0).host_isactive[isubmesh+1][jsubmesh-1] + h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh + 1;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }

    if (bottom_edge & !left_edge & !right_edge){

        if (jsubmesh==0){
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                *in_domain = true;
                *on_bdry = true;
                *bdry_tag = jsubmesh*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                *bdry_dim = 1;
                return;
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }
        else{
            if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1]){
                *in_domain = true;
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] + h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1] == 1){
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(2*h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *on_bdry = false;
                    *bdry_tag = -1;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                *in_domain = false;
                *on_bdry = false;
                *bdry_tag = -999;
                *bdry_dim = -1;
                return;
            }
        }

    }

    if (bottom_edge & left_edge){
        if (jsubmesh==0){
            if (isubmesh==0){
                if (h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
        else{
            if (isubmesh==0){
                if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1]){
                    *in_domain = true;
                    *on_bdry = true;
                    *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                    *bdry_dim = 0;
                    return;
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
            else{
                if (h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh-1] | h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh] |
                    h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1] | h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                    *in_domain = true;
                    int sum = h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh-1] + h_pumi_mesh(0).host_isactive[isubmesh-1][jsubmesh] +
                            h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh-1] + h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh];
                    if (sum < 4){
                        *on_bdry = true;
                        *bdry_tag = jsubmesh*(h_pumi_mesh(0).nsubmesh_x1 + 1) + isubmesh;
                        *bdry_dim = 0;
                        return;
                    }
                }
                else{
                    *in_domain = false;
                    *on_bdry = false;
                    *bdry_tag = -999;
                    *bdry_dim = -1;
                    return;
                }
            }
        }
    }
}

double get_global_left_coord(MBBL pumi_obj){
    return pumi_obj.host_submesh_x1[0].xmin;
}

double get_global_right_coord(MBBL pumi_obj){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);
    int nsubmesh = h_pumi_mesh(0).nsubmesh_x1;
    return pumi_obj.host_submesh_x1[nsubmesh-1].xmax;
}

double get_global_bottom_coord(MBBL pumi_obj){
    return pumi_obj.host_submesh_x2[0].xmin;
}

double get_global_top_coord(MBBL pumi_obj){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);
    int nsubmesh = h_pumi_mesh(0).nsubmesh_x2;
    return pumi_obj.host_submesh_x2[nsubmesh-1].xmax;
}

void print_mesh_skeleton(MBBL pumi_obj){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_obj.mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_obj.mesh);
    bool on_bdry, in_domain;
    int bdry_tag, bdry_dim;
    printf("\n\nPrinting the skeleton of the mesh\n");
    printf("E --> Boundary block-edge\n");
    printf("V --> Boundary block-vertex\n\n");
    for (int jsubmesh=h_pumi_mesh(0).nsubmesh_x2-1; jsubmesh>=0; jsubmesh--){
        int Jnp = pumi_obj.host_submesh_x2[jsubmesh].Nel + pumi_obj.host_submesh_x2[jsubmesh].Nel_cumulative;
        for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                if (bdry_tag+1){
                    printf("%3dV",bdry_tag );
                }
                else{
                    printf("----",bdry_tag );
                }
            }
            else{
                printf("    ");
            }
            Inp++;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                if (bdry_tag+1){
                    printf("----%3dE----",bdry_tag );
                }
                else{
                    printf("------------",bdry_tag );
                }
            }
            else{
                printf("            ");
            }
            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                Inp = pumi_obj.host_submesh_x1[isubmesh].Nel + pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
                pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    printf("%3dV",bdry_tag );
                }
                else{
                    printf("    ");
                }
            }
        }
        printf("\n");
        Jnp--;
        for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                printf("   |");
            }
            else{
                printf("    ");
            }
            Inp++;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                printf("            ");
            }
            else{
                printf("            ");
            }

            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                Inp = pumi_obj.host_submesh_x1[isubmesh].Nel + pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
                pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    printf("   |");
                }
                else{
                    printf("    ");
                }
            }
        }
        printf("\n");
        for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                if (bdry_tag+1){
                    printf("%3dE",bdry_tag);
                }
                else{
                    printf("   |",bdry_tag);
                }

            }
            else{
                printf("    ");
            }
            Inp++;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                printf("            ");
            }
            else{
                printf("            ");
            }

            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                Inp = pumi_obj.host_submesh_x1[isubmesh].Nel + pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
                pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    printf("%3dE",bdry_tag);
                }
                else{
                    printf("    ");
                }
            }
        }
        printf("\n");
        for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                printf("   |");
            }
            else{
                printf("    ");
            }
            Inp++;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                printf("            ");
            }
            else{
                printf("            ");
            }

            if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                Inp = pumi_obj.host_submesh_x1[isubmesh].Nel + pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
                pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    printf("   |");
                }
                else{
                    printf("    ");
                }
            }
        }

        if (jsubmesh==0){
            printf("\n");
            Jnp = 0;
            for (int isubmesh=0; isubmesh<h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
                int Inp = pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
                pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    printf("%3dV",bdry_tag );
                }
                else{
                    printf("    ");
                }
                Inp++;
                pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    if (bdry_tag+1){
                        printf("----%3dE----",bdry_tag );
                    }
                    else{
                        printf("------------",bdry_tag );
                    }
                }
                else{
                    printf("            ");
                }
                if (isubmesh==h_pumi_mesh(0).nsubmesh_x1-1){
                    Inp = pumi_obj.host_submesh_x1[isubmesh].Nel + pumi_obj.host_submesh_x1[isubmesh].Nel_cumulative;
                    pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                    if (in_domain){
                        printf("%3dV",bdry_tag );
                    }
                    else{
                        printf("    ");
                    }
                }
            }
        }
        printf("\n");
    }
}

} // namespace pumi

#endif
