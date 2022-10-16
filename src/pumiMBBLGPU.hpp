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
    arbitrary  = 0x08, //!< mesh with arbitrary element sizing
};

///////// PUMI-MBBL-GPU Data-Structures ////////////////////////////////////////////////////

using DoubleView = Kokkos::View<double*>;

/**
 * @brief Base class for submesh-blocks which stores parameters defining the mesh in a block
 *
 * Different meshtypes are derived from this class
 */
class Submesh{
public:

    double xmin; //!< Min-side coordinates
    double xmax; //!< Max-side coordinates
    int Nel; //!< Number of elements in the block
    double t0; //!< Smallest element size in the block (or element size for uniform blocks)
    double r; //!< Grading ratio for element sizes in the block
    Meshtype meshtype; //<! Type of meshing in the block
    double length; //!< Length of the block
    int Nel_cumulative; //!< Number of elements in the preceding blocks

    double r_by_t0; //!< value of (r-1.0)/t0 -- value needed in analytical cell-locate functions
    double log_r; //!< value of log(r) -- value needed in analtyical cell-locate functions
    DoubleView BL_coords; //<! node coords to be stored for non-uniform blocks for faster particle search
    double *host_BL_coords; //<! host copy of node coords to be stored for non-uniform blocks for faster particle search
    /**
    * @brief Default constructor.
    */
    Submesh(){};
    /**
    * @brief Class constructor.
    *
    * @param[in] submesh min-side coords
    * @param[in] submesh max-side coords
    * @param[in] submesh number of elements
    * @param[in] submesh smallest element size
    * @param[in] submesh grading ratio
    * @param[in] submesh mesh-type
    * @param[in] submesh length
    * @param[in] submesh preceding cumulative elements
    * @param[in] submesh (r-1.0)/t0 value
    * @param[in] submesh log(r) value
    * @param[in] submesh BL coordinates (explicitly stored)
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
            DoubleView BL_coords_,
            double* host_BL_coords_):
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
            BL_coords(BL_coords_),
            host_BL_coords(host_BL_coords_){};

    virtual int locate_cell_host(double ) = 0;
    virtual int update_cell_host(double , int) = 0;
    virtual double elem_size_host(int ) = 0;
    virtual void calc_weights_host(double , int , int *global_cell, double *Wgh2) = 0;
    virtual double node_coords_host(int) = 0;
    virtual double grading_ratio_host(int) = 0;

    KOKKOS_INLINE_FUNCTION
    virtual ~Submesh(){};

};

using SubmeshDeviceViewPtr = Kokkos::View<DevicePointer<Submesh>*>; //!< Pointer to array of Submesh objects (in device space) poniting to address of derived class
using SubmeshHostViewPtr = Submesh**; //<! Pointer to array of submesh objects (in CPU memory)
/**
 * @brief Uniform submesh class derived from submesh class
 *
 */
class Uniform_Submesh : public Submesh{
public:
    /**
    * @brief Class constructor.
    *
    * @param[in] submesh min-side coords
    * @param[in] submesh max-side coords
    * @param[in] submesh number of elements
    * @param[in] submesh smallest element size
    * @param[in] submesh grading ratio
    * @param[in] submesh length
    * @param[in] submesh preceding cumulative elements
    * @param[in] submesh (r-1.0)/t0 value
    * @param[in] submesh log(r) value
    * @param[in] submesh node coordinates (explicitly stored)
    * @param[in] submesh node coordinates (host copy)
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
                    DoubleView BL_coords_,
                    double *host_BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,uniform,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_,host_BL_coords_){};

    /**
     * @brief Locates local cell ID of given point inside a block
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        if (q==this->xmax)
            return this->Nel-1;
        else
            return (q - this->xmin)/this->t0;
    }

    /**
     * @brief updates local cell ID of given particle coordinate
     * @note for uniform mesh update is same as locating particle
     * @param[in] point coordinate to be located
     * @param[in] previous local cell ID of particle
     * @return upddated local cell ID
     */
    KOKKOS_INLINE_FUNCTION
    int update_cell(double q, int){
        if (q==this->xmax)
            return this->Nel-1;
        else
            return (q - this->xmin)/this->t0;
    }

    /**
     * @brief Computes element size inside a block
     * @param[in] local cell ID for which elemetn size is needed
     * @return Value of element size
     */
    KOKKOS_INLINE_FUNCTION
    double elem_size(int){
        return this->t0;
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = ((q - this->xmin) - local_cell*this->t0)/this->t0;
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Computes node coordinate inside submesh block
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    KOKKOS_INLINE_FUNCTION
    double node_coords(int inode) {
        return this->xmin + inode*this->t0;
    }

    /**
     * @brief Fetches grading ration inside submesh block
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    KOKKOS_INLINE_FUNCTION
    double grading_ratio(int){
        return 1.0;
    }

    /**
     * @brief Locates local cell ID of given point inside a block
     * @note executable on host-only
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    int locate_cell_host(double q) {
        if (q==this->xmax)
            return this->Nel-1;
        else
            return (q - this->xmin)/this->t0;
    }

    /**
     * @brief updates local cell ID of given particle coordinate
     * @note executable on host-only
     * @param[in] point coordinate to be located
     * @param[in] previous local cell ID of particle
     * @return upddated local cell ID
     */
    int update_cell_host(double q, int ) {
        if (q==this->xmax)
            return this->Nel-1;
        else
            return (q - this->xmin)/this->t0;
    }

    /**
     * @brief Computes element size inside a block
     * @note executable on host-only
     * @param[in] local cell ID for which elemetn size is needed
     * @return Value of element size
     */
    double elem_size_host(int ) {
        return this->t0;
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @note executable on host-only
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    void calc_weights_host(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = ((q - this->xmin) - local_cell*this->t0)/this->t0;
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Computes node coordinate inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    double node_coords_host(int inode) {
        return this->xmin + inode*this->t0;
    }

    /**
     * @brief Fetches grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    double grading_ratio_host(int){
        return 1.0;
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
    * @param[in] submesh min-side coords
    * @param[in] submesh max-side coords
    * @param[in] submesh number of elements
    * @param[in] submesh smallest element size
    * @param[in] submesh grading ratio
    * @param[in] submesh length
    * @param[in] submesh preceding cumulative elements
    * @param[in] submesh (r-1.0)/t0 value
    * @param[in] submesh log(r) value
    * @param[in] submesh BL coordinates (explicitly stored)
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
                    DoubleView BL_coords_,
                    double *host_BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,minBL,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_,host_BL_coords_){};

    /**
     * @brief Locates local cell ID of given point inside a block
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        if (r != 1.0){
            if (q==this->xmax){
                return this->Nel-1;
            }
            else{
                int cell = log(1.0 + (q - this->xmin)*this->r_by_t0)/this->log_r;
                return cell;
            }
        }
        else{
            if (q==this->xmax)
                return this->Nel-1;
            else
                return (q - this->xmin)/this->t0;
        }
    }

    /**
     * @brief updates local cell ID of given particle coordinate
     * @note for non-uniform mesh update is done with adjacency search on stored node coords
     * @param[in] point coordinate to be located
     * @param[in] previous local cell ID of particle
     * @return upddated local cell ID
     */
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

    /**
     * @brief Computes element size inside a block
     * @param[in] local cell ID for which elemetn size is needed
     * @return Value of element size
     */
    KOKKOS_INLINE_FUNCTION
    double elem_size(int icell){
        return (this->BL_coords(icell+1)-this->BL_coords(icell));
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->BL_coords(local_cell))/(this->BL_coords(local_cell+1)-this->BL_coords(local_cell));
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Fetches node coordinate inside submesh block
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    KOKKOS_INLINE_FUNCTION
    double node_coords(int inode) {
        return this->BL_coords(inode);
    }

    /**
     * @brief Fetches grading ration inside submesh block
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    KOKKOS_INLINE_FUNCTION
    double grading_ratio(int){
        return this->r;
    }

    /**
     * @brief Locates local cell ID of given point inside a block
     * @note executable on host-only
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    int locate_cell_host(double q) {
        if (r != 1.0){
            if (q==this->xmax){
                return this->Nel-1;
            }
            else{
                int cell = log(1.0 + (q - this->xmin)*this->r_by_t0)/this->log_r;
                return cell;
            }
        }
        else{
            if (q==this->xmax)
                return this->Nel-1;
            else
                return (q - this->xmin)/this->t0;
        }
    }

    /**
     * @brief updates local cell ID of given particle coordinate
     * @note executable on host-only
     * @param[in] point coordinate to be located
     * @param[in] previous local cell ID of particle
     * @return upddated local cell ID
     */
    int update_cell_host(double q, int icell) {
        while (q < this->host_BL_coords[icell]){
            icell--;
        }
        while (q > this->host_BL_coords[icell+1]){
            icell++;
        }
        return icell;
    }

    /**
     * @brief Computes element size inside a block
     * @note executable on host-only
     * @param[in] local cell ID for which elemetn size is needed
     * @return Value of element size
     */
    double elem_size_host(int icell) {
        return (this->host_BL_coords[icell+1]-this->host_BL_coords[icell]);
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @note executable on host-only
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    void calc_weights_host(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->host_BL_coords[local_cell])/(this->host_BL_coords[local_cell+1]-this->host_BL_coords[local_cell]);
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Fetches node coordinate inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    double node_coords_host(int inode) {
        return this->host_BL_coords[inode];
    }

    /**
     * @brief Fetches grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    double grading_ratio_host(int){
        return this->r;
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
    * @param[in] submesh min-side coords
    * @param[in] submesh max-side coords
    * @param[in] submesh number of elements
    * @param[in] submesh smallest element size
    * @param[in] submesh grading ratio
    * @param[in] submesh length
    * @param[in] submesh preceding cumulative elements
    * @param[in] submesh (r-1.0)/t0 value
    * @param[in] submesh log(r) value
    * @param[in] submesh BL coordinates (explicitly stored)
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
                    DoubleView BL_coords_,
                    double *host_BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,maxBL,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_,host_BL_coords_){};

    /**
     * @brief Locates local cell ID of given point inside a block
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        if (r != 1.0){
            if (q==this->xmin){
                return 0;
            }
            else{
                int cell = log(1.0 + (this->xmax - q)*this->r_by_t0)/this->log_r;
                return this->Nel - cell - 1;
            }
        }
        else{
            if (q==this->xmax)
                return this->Nel-1;
            else
                return (q - this->xmin)/this->t0;
        }
    }

    /**
     * @brief updates local cell ID of given particle coordinate
     * @note for non-uniform mesh update is done with adjacency search on stored node coords
     * @param[in] point coordinate to be located
     * @param[in] previous local cell ID of particle
     * @return upddated local cell ID
     */
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

    /**
     * @brief Computes element size inside a block
     * @param[in] local cell ID for which elemetn size is needed
     * @return Value of element size
     */
    KOKKOS_INLINE_FUNCTION
    double elem_size(int icell){
        return (this->BL_coords(icell+1)-this->BL_coords(icell));
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->BL_coords(local_cell))/(this->BL_coords(local_cell+1)-this->BL_coords(local_cell));
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Fetches node coordinate inside submesh block
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    KOKKOS_INLINE_FUNCTION
    double node_coords(int inode) {
        return this->BL_coords(inode);
    }

    /**
     * @brief Fetches grading ration inside submesh block
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    KOKKOS_INLINE_FUNCTION
    double grading_ratio(int){
        return 1.0/this->r;
    }

    /**
     * @brief Locates local cell ID of given point inside a block
     * @note executable on host-only
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    int locate_cell_host(double q) {
        if (r != 1.0){
            if (q==this->xmin){
                return 0;
            }
            else{
                int cell = log(1.0 + (this->xmax - q)*this->r_by_t0)/this->log_r;
                return this->Nel - cell - 1;
            }
        }
        else{
            if (q==this->xmax)
                return this->Nel-1;
            else
                return (q - this->xmin)/this->t0;
        }
    }

    /**
     * @brief Computes grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    int update_cell_host(double q, int icell) {
        while (q < this->host_BL_coords[icell]){
            icell--;
        }
        while (q > this->host_BL_coords[icell+1]){
            icell++;
        }
        return icell;
    }

    /**
     * @brief Computes grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    double elem_size_host(int icell) {
        return (this->host_BL_coords[icell+1]-this->host_BL_coords[icell]);
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @note executable on host-only
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    void calc_weights_host(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->host_BL_coords[local_cell])/(this->host_BL_coords[local_cell+1]-this->host_BL_coords[local_cell]);
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Fetches node coordinate inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    double node_coords_host(int inode) {
        return this->host_BL_coords[inode];
    }

    /**
     * @brief Fetches grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    double grading_ratio_host(int){
        return 1.0/this->r;
    }

};


/**
 * @brief MaxBL submesh class derived from submesh class
 *
 */
class Arbitrary_Submesh : public Submesh{
public:
    /**
    * @brief Class constructor.
    *
    * @param[in] submesh min-side coords
    * @param[in] submesh max-side coords
    * @param[in] submesh number of elements
    * @param[in] submesh smallest element size
    * @param[in] submesh grading ratio
    * @param[in] submesh length
    * @param[in] submesh preceding cumulative elements
    * @param[in] submesh (r-1.0)/t0 value
    * @param[in] submesh log(r) value
    * @param[in] submesh BL coordinates (explicitly stored)
    */
    Arbitrary_Submesh(double xmin_,
                    double xmax_,
                    int Nel_,
                    double t0_,
                    double r_,
                    double length_,
                    int Nel_cumulative_,
                    double r_by_t0_,
                    double log_r_,
                    DoubleView BL_coords_,
                    double *host_BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,arbitrary,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_,host_BL_coords_){};

    /**
     * @brief binary search to locate local cell ID of given point
     * @param[in] array of node coordinates
     * @param[in] first index of window in which bst is performed
     * @param[in] last index of window in which bst is performed
     * @param[in] point coordinate to be located
     */
    KOKKOS_INLINE_FUNCTION
    int bst_search(DoubleView arr, int first, int last, double val){
        int mid = (first+last)/2;

        if (last == first+1){
            if ((arr(first) <= val) && (arr(last) >= val)){
                return first;
            }
        }
        else{
            if (arr(mid) >= val){
                last = mid;
                return bst_search(arr, first, last, val);
            }
            else if (arr(mid) < val){
                first = mid;
                return bst_search(arr, first, last, val);
            }
        }
        return -1;
    }

    /**
     * @brief Locates local cell ID of given point inside a block
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    KOKKOS_INLINE_FUNCTION
    int locate_cell(double q){
        int iel = bst_search(this->BL_coords,0,this->Nel,q);
        return iel;
    }

    /**
     * @brief updates local cell ID of given particle coordinate
     * @note for non-uniform mesh update is done with adjacency search on stored node coords
     * @param[in] point coordinate to be located
     * @param[in] previous local cell ID of particle
     * @return upddated local cell ID
     */
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

    /**
     * @brief Computes element size inside a block
     * @param[in] local cell ID for which elemetn size is needed
     * @return Value of element size
     */
    KOKKOS_INLINE_FUNCTION
    double elem_size(int icell){
        return (this->BL_coords(icell+1)-this->BL_coords(icell));
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    KOKKOS_INLINE_FUNCTION
    void calc_weights(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->BL_coords(local_cell))/(this->BL_coords(local_cell+1)-this->BL_coords(local_cell));
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Fetches node coordinate inside submesh block
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    KOKKOS_INLINE_FUNCTION
    double node_coords(int inode) {
        return this->BL_coords(inode);
    }

    /**
     * @brief Computes grading ration inside submesh block
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    KOKKOS_INLINE_FUNCTION
    double grading_ratio(int inode){
        double max = this->BL_coords(inode+1)-this->BL_coords(inode);
        double min = this->BL_coords(inode)-this->BL_coords(inode-1);
        return max/min;
    }

    /**
     * @brief binary search to locate local cell ID of given point
     * @note executable on host-only
     * @param[in] array of node coordinates (host copy)
     * @param[in] first index of window in which bst is performed
     * @param[in] last index of window in which bst is performed
     * @param[in] point coordinate to be located
     */
    int bst_search_host(double* arr, int first, int last, double val){
        int mid = (first+last)/2;

        if (last == first+1){
            if ((arr[first] <= val) && (arr[last] >= val)){
                return first;
            }
        }
        else{
            if (arr[mid] >= val){
                last = mid;
                return bst_search_host(arr, first, last, val);
            }
            else if (arr[mid] < val){
                first = mid;
                return bst_search_host(arr, first, last, val);
            }
        }
        return -1;
    }

    /**
     * @brief Locates local cell ID of given point inside a block thru binsary search
     * @note executable on host-only
     * @param[in] point coordinate to be located
     * @return local cell ID
     */
    int locate_cell_host(double q) {
        int iel = bst_search_host(this->host_BL_coords,0,this->Nel,q);
        return iel;
    }

    /**
     * @brief Computes grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    int update_cell_host(double q, int icell) {
        while (q < this->host_BL_coords[icell]){
            icell--;
        }
        while (q > this->host_BL_coords[icell+1]){
            icell++;
        }
        return icell;
    }

    /**
     * @brief Computes grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    double elem_size_host(int icell) {
        return (this->host_BL_coords[icell+1]-this->host_BL_coords[icell]);
    }

    /**
     * @brief Computes linear 1D weights for gather/scatter operations and
     *        directional global cell ID
     * @note executable on host-only
     * @param[in] point coordinate inside submesh block
     * @param[in] local cell ID of point inside submesh block
     * @param[out] directional global cell ID
     * @param[out] linear 1D weight (correspoding to max-side node)
     */
    void calc_weights_host(double q, int local_cell, int *global_cell, double *Wgh2){
        *Wgh2 = (q - this->host_BL_coords[local_cell])/(this->host_BL_coords[local_cell+1]-this->host_BL_coords[local_cell]);
        *global_cell = local_cell + this->Nel_cumulative;
    }

    /**
     * @brief Fetches node coordinate inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of node coordinate
     */
    double node_coords_host(int inode) {
        return this->host_BL_coords[inode];
    }

    /**
     * @brief Fetches grading ration inside submesh block
     * @note executable on host-only
     * @param[in] local node ID inside block
     * @return Value of grading ratio at given node
     */
    double grading_ratio_host(int inode){
        double max = this->host_BL_coords[inode+1]-this->host_BL_coords[inode];
        double min = this->host_BL_coords[inode]-this->host_BL_coords[inode-1];
        return max/min;
    }

};

/**
 * @brief Unassigned submesh class derived from submesh class for padded blocks
 *
 */
class Unassigned_Submesh : public Submesh{
public:
    /**
    * @brief Class constructor.
    *
    * @param[in] submesh min-side coords
    * @param[in] submesh max-side coords
    * @param[in] submesh number of elements
    * @param[in] submesh smallest element size
    * @param[in] submesh grading ratio
    * @param[in] submesh length
    * @param[in] submesh preceding cumulative elements
    * @param[in] submesh (r-1.0)/t0 value
    * @param[in] submesh log(r) value
    * @param[in] submesh BL coordinates (explicitly stored)
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
                    DoubleView BL_coords_,
                    double *host_BL_coords_):
                    Submesh(xmin_,xmax_,Nel_,t0_,r_,unassigned,length_,Nel_cumulative_,r_by_t0_,log_r_,BL_coords_,host_BL_coords_){};

    KOKKOS_INLINE_FUNCTION
    int locate_cell(double){
        return -1;
    }

    KOKKOS_INLINE_FUNCTION
    int update_cell(double, int ){
        return -1;
    }

    KOKKOS_INLINE_FUNCTION
    double elem_size(int ){
        return -999.0;
    }

    KOKKOS_INLINE_FUNCTION
    void calc_weights(double , int , int *global_cell, double *Wgh2){
        *Wgh2 = -999.0;
        *global_cell = -1;
    }

    KOKKOS_INLINE_FUNCTION
    double node_coords(int) {
        return -999.0;
    }

    int locate_cell_host(double ) { return -1; }
    int update_cell_host(double , int) { return -1; }
    double elem_size_host(int ) { return -999.0; }
    void calc_weights_host(double , int , int *global_cell, double *Wgh2){
        *global_cell = -1;
        *Wgh2 = -999.0;
    }
    double node_coords_host(int ){
        return -999.0;
    }
    double grading_ratio_host(int){
        return -999.0;
    }

};

/**
 * @brief Mesh offsets class storing auxiliary data to compute node/element ID
 *        for a 2D mesh with inactive submesh blocks
 *
 */
class MeshOffsets{
public:
    bool is_fullmesh;

    Kokkos::View<int**> nodeoffset_start; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> nodeoffset_skip_bot; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> nodeoffset_skip_mid; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> nodeoffset_skip_top; //!< aux data structure to compute nodeoffset
    Kokkos::View<int**> elemoffset_start; //!< aux data structure to compute element offset
    Kokkos::View<int*> elemoffset_skip; //!< aux data structure to compute element offset

    int** host_nodeoffset_start; //!< aux data structure to compute nodeoffset (host copy)
    int** host_nodeoffset_skip_bot; //!< aux data structure to compute nodeoffset (host copy)
    int** host_nodeoffset_skip_mid; //!< aux data structure to compute nodeoffset (host copy)
    int** host_nodeoffset_skip_top; //!< aux data structure to compute nodeoffset (host copy)
    int** host_elemoffset_start; //!< aux data structure to compute element offset (host copy)
    int* host_elemoffset_skip; //!< aux data structure to compute element offset (host copy)

    int Nel_total; //!< Total active elements in domain
    int Nnp_total; //!< Total active nodes in domain

    /**
    * @brief Default Class constructor.
    */
    MeshOffsets(){
        host_nodeoffset_start = nullptr;
        host_nodeoffset_skip_bot = nullptr;
        host_nodeoffset_skip_mid = nullptr;
        host_nodeoffset_skip_top = nullptr;
        host_elemoffset_start = nullptr;
        host_elemoffset_skip = nullptr;
    };

    /**
    * @brief Class constructor.
    *
    * @param[in] Array of x1-submesh object pointers on host
    * @param[in] Number of x1-blocks in mesh
    * @param[in] Array of x2-submesh object pointers on host
    * @param[in] Number of x1-blocks in mesh
    */
    MeshOffsets(SubmeshHostViewPtr , int, SubmeshHostViewPtr, int, bool**);
};
/**
 * @brief Class to store classification information of mesh entities
 *
 */
class MeshBdry{
public:

    Kokkos::View<bool*> is_bdry_edge; //!< bool value stores if an edge is on boundary
    Kokkos::View<Vector3*> bdry_edge_normal; //!< boundary normal direction for each edge
    Kokkos::View<bool*> is_bdry_vert; //!< bool value stores if an vertex is on boundary
    Kokkos::View<Vector3*> bdry_vert_normal; //!< boundary normal direction for each vertex
    Kokkos::View<int*> edge_to_face; //!< bdry face ID of first element on edge

    bool* host_is_bdry_edge; //!< bool value stores if an edge is on boundary (host copy)
    Vector3* host_bdry_edge_normal; //!< boundary normal direction for each edge (host copy)
    bool* host_is_bdry_vert; //!< bool value stores if an vertex is on boundary (host copy)
    Vector3* host_bdry_vert_normal; //!< boundary normal direction for each vertex (host copy)
    int *host_edge_to_face; //!< bdry face ID of first element on edge  (host copy)

    int Nbdry_faces; //!< int value storing number of boundary element faces

    /**
    * @brief Default Class constructor.
    */
    MeshBdry(){};

    /**
    * @brief Class constructor.
    *
    * @param[in] Number of bdry faces
    */
    MeshBdry(int Nbdry_faces_):Nbdry_faces(Nbdry_faces_){
        host_is_bdry_edge = nullptr;
        host_bdry_edge_normal = nullptr;
        host_is_bdry_vert = nullptr;
        host_bdry_vert_normal = nullptr;
        host_edge_to_face = nullptr;
    };

    /**
    * @brief Class constructor.
    *
    * @param[in] Array of x1-submesh object pointers on host
    * @param[in] Number of x1-blocks in mesh
    * @param[in] Array of x2-submesh object pointers on host
    * @param[in] Number of x1-blocks in mesh
    */
    MeshBdry(SubmeshHostViewPtr , int, SubmeshHostViewPtr, int, bool**);
};

/**
 * @brief Class to store informations regarding block interface entities
 *
 */
class BlockInterface{
public:
    // 1D
    Kokkos::View<double*> if_x1_r; //!< x1-block-interface grading ratio
    Kokkos::View<int*> if_x1_node; //!< x1-block-interface global x1-node ID
    Kokkos::View<double*> if_x2_r; //!< x2-block-interface grading ratio
    Kokkos::View<int*> if_x2_node; //!< x2-block-interface global x2-node ID

    double *host_if_x1_r; //!< x1-block-interface grading ratio (host copy)
    int *host_if_x1_node; //!< x1-block-interface global x1-node ID (host copy)
    double *host_if_x2_r; //!< x2-block-interface grading ratio (host copy)
    int *host_if_x2_node; //!< x2-block-interface global x2-node ID (host copy)

    // 2D
    Kokkos::View<int*> vert_nodeID; //!< global node ID of vertex
    Kokkos::View<int*> vert_subID; //!< flattedned active submesh ID to which vertex belongs to
    Kokkos::View<int*> edge_first_nodeID; //!< global node ID of first node on edge
    Kokkos::View<int*> edge_subID; //!< flattedned active submesh ID to which edge belongs to

    int *host_vert_nodeID; //!< global node ID of vertex (host copy)
    int *host_vert_subID; //!< flattedned active submesh ID to which vertex belongs to (host copy)
    int *host_edge_first_nodeID; //!< global node ID of first node on edge (host copy)
    int *host_edge_subID; //!< flattedned active submesh ID to which edge belongs to (host copy)

    /**
    * @brief Default Class constructor.
    */
    BlockInterface(){};

    /**
    * @brief Class constructor.
    *
    * @param[in] Array of x1-submesh object pointers on host
    * @param[in] Number of x1-blocks in mesh
    */
    BlockInterface(SubmeshHostViewPtr , int);

    /**
    * @brief Class constructor.
    *
    * @param[in] Array of x1-submesh object pointers on host
    * @param[in] Number of x1-blocks in mesh
    * @param[in] Array of x2-submesh object pointers on host
    * @param[in] Number of x1-blocks in mesh
    */
    BlockInterface(SubmeshHostViewPtr , int, SubmeshHostViewPtr, int, bool**);

    /**
     * @brief sets global node IDs for vertex nodes and first nodes on edge
     * @param[in] array of vertex node IDs
     * @param[in] array of edge first node IDs
     * @param[in] number of x1-blocks in mesh
     * @param[in] number of x2-blocks in mesh
     */
    void set_interface_nodeIDs(int *vert_nodeIDs, int *edge_first_nodeIDs, int Nx, int Ny);
};

/**
 * @brief Class to store data structures needed to perform faster
 *        binary search to locate active blocks to which different
 *        classified nodes belong to
 *
 */
class MeshBST{
public:
    int total_active_blocks; //!< total active blocks in domain
    int total_block_nodes; //!< total active block-interior nodes
    int total_active_edges; //!< total active edges in domain
    int total_edge_nodes; //!< total active edge nodes

    Kokkos::View<int*> active_blockID; //!< mapping b/w active block to flatened submesh ID
    Kokkos::View<int*> block_nodes_cumulative; //!< total active block-interior nodes in all preceding active blocks
    Kokkos::View<int*> block_elems_cumulative; //!< total active block elements in all preceding active blocks
    Kokkos::View<int*> active_edgeID; //!< mapping b/w active edges to actual edge ID
    Kokkos::View<int*> edge_nodes_cumulative; //!< total active edge nodes in all preceding active blocks

    int* host_active_blockID; //!< mapping b/w active block to flatened submesh ID
    int* host_block_nodes_cumulative; //!< total active block-interior nodes in all preceding active blocks (host copy)
    int* host_block_elems_cumulative; //!< total active block elements in all preceding active blocks (host copy)
    int* host_active_edgeID; //!< mapping b/w active edges to actual edge ID (host copy)
    int* host_edge_nodes_cumulative; //!< total active edge nodes in all preceding active blocks (host copy)

    /**
    * @brief Default Class constructor.
    */
    MeshBST(){
        host_active_blockID = nullptr;
        host_block_nodes_cumulative = nullptr;
        host_block_elems_cumulative = nullptr;
        host_active_edgeID = nullptr;
        host_edge_nodes_cumulative = nullptr;
    };

    /**
     * @brief intialize data stuctures in class for BST searches
     * @param[in] object to BlockInterface class
     * @param[in] Array of x1-submesh object pointers on host
     * @param[in] Number of x1-blocks in mesh
     * @param[in] Array of x2-submesh object pointers on host
     * @param[in] Number of x1-blocks in mesh
     * @param[in] 2D array of block activity
     */
    void initialize_MeshBST(BlockInterface blkif, SubmeshHostViewPtr, int, SubmeshHostViewPtr, int, bool**);
};

/**
 * @brief Mesh class
 *
 * Object of this class can access all mesh-level properties and other
 * stored auxiliary data structures
 */
class Mesh{
public:
    int ndim; //!< dimensions of the domain
    int nsubmesh_x1; //!< number of blocks in x1-direction
    int Nel_tot_x1; //!< Total number of elements in x1-direction
    int nsubmesh_x2; //!< number of blocks in x2-direction
    int Nel_tot_x2; //!< Total number of elements in x2-direction

    Kokkos::View<bool**> isactive; //!< 2D bool-array defining the activity of blocks
    bool **host_isactive; //!< 2D bool-array defining the activity of blocks (host copy)

    MeshOffsets offsets; //!< Object to class MeshOffsets class
    MeshBdry bdry; //!< Object to class MeshBdry class
    BlockInterface blkif; //!< Object to class BlockInterface class
    MeshBST bst; //!< Object to class MeshBST class


    int Nel_total; //!< Total active elements in domain
    int Nnp_total; //!< Total active nodes in domain

    int nsubmesh_x3; //!< number of blocks in x3-direction
    int Nel_tot_x3; //!< Total number of elements in x3-direction

    /**
    * @brief Default constructor.
    */
    Mesh(){};

    /**
    * @brief Constructor for 1D Mesh
    * @param[in] number of x1-submesh blocks
    * @param[in] total number of elements along x1-direction
    * @param[in] Object to MeshBdry class
    * @param[in] Object to BlockInterface class
    */
    Mesh(int nsubmesh_x1_,
         int Nel_tot_x1_,
         MeshBdry bdry_,
         BlockInterface blkif_):
         ndim(1),
         nsubmesh_x1(nsubmesh_x1_),
         Nel_tot_x1(Nel_tot_x1_),
         bdry(bdry_),
         blkif(blkif_){
             host_isactive = nullptr;
             Nel_total = Nel_tot_x1_;
             Nnp_total = Nel_tot_x1_+1;
             nsubmesh_x2 = 0;
             nsubmesh_x3 = 0;
             Nel_tot_x2 = 0;
             Nel_tot_x3 = 0;
         };
     /**
     * @brief Constructor for 2D Mesh
     * @param[in] number of x1-submesh blocks
     * @param[in] total number of elements along x1-direction
     * @param[in] number of x2-submesh blocks
     * @param[in] total number of elements along x2-direction
     * @param[in] submesh activity info
     * @param[in] submesh activity info (host copy)
     * @param[in] Object to MeshOffsets class
     * @param[in] Object to MeshBdry class
     * @param[in] Object to BlockInterface class
     * @param[in] Object to MeshBST class
     * @param[in] total elements in mesh
     * @param[in] total nodes in mesh
     */
     Mesh(int nsubmesh_x1_,
         int Nel_tot_x1_,
         int nsubmesh_x2_,
         int Nel_tot_x2_,
         Kokkos::View<bool**> isactive_,
         bool** host_isactive_,
         MeshOffsets offsets_,
         MeshBdry bdry_,
         BlockInterface blkif_,
         MeshBST bst_,
         int Nel_total_,
         int Nnp_total_):
         ndim(2),
         nsubmesh_x1(nsubmesh_x1_),
         Nel_tot_x1(Nel_tot_x1_),
         nsubmesh_x2(nsubmesh_x2_),
         Nel_tot_x2(Nel_tot_x2_),
         isactive(isactive_),
         host_isactive(host_isactive_),
         offsets(offsets_),
         bdry(bdry_),
         blkif(blkif_),
         bst(bst_),
         Nel_total(Nel_total_),
         Nnp_total(Nnp_total_)
         {
             nsubmesh_x3 = 0;
             Nel_tot_x3 = 0;
         };

};

/**
 * @brief Wrapper structure containing mesh and submesh objects
 *
 * Object of this class will be used to call all mesh related APIs
 */
struct MBBL{
    Mesh mesh;
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
    * @param[in] Mesh object in allocated in GPU
    * @param[in] x1-submesh object in GPU
    * @param[in] copy of x1-submesh object in CPU
    */
    MBBL(Mesh mesh_,
         SubmeshDeviceViewPtr submesh_x1_,
         SubmeshHostViewPtr host_submesh_x1_):
         mesh(mesh_),
         submesh_x1(submesh_x1_),
         host_submesh_x1(host_submesh_x1_){
             host_submesh_x2 = nullptr;
         };

    /**
    * @brief Constructor for 2D Wrapper structure
    * @param[in] Mesh object in allocated in GPU
    * @param[in] x1-submesh object in GPU
    * @param[in] copy of x1-submesh object in CPU
    * @param[in] x2-submesh object in GPU
    * @param[in] copy of x2-submesh object in CPU
    */
    MBBL(Mesh mesh_,
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

//For 1D and 2D tests ONLY -- NOT TO BE USED IN HPIC2
class ParticleData{
public:
    double x1;
    double x2;
    int cellID;
    int submeshID;
    bool part_active;
    int exit_faceID;
    double w1;
    double w2;
    double w3;
    double w4;

    ParticleData(){};

    ParticleData(double x1_,double x2_):
                x1(x1_),x2(x2_){};

    KOKKOS_INLINE_FUNCTION
    ParticleData(double x1_, double x2_, int submeshID_, int cellID_, bool part_active_, int exit_faceID_):
                x1(x1_),x2(x2_),cellID(cellID_),submeshID(submeshID_),part_active(part_active_),exit_faceID(exit_faceID_){};

    KOKKOS_INLINE_FUNCTION
    ParticleData(double x1_,
                double x2_,
                int submeshID_,
                int cellID_,
                bool part_active_,
                int exit_faceID_,
                double w1_,
                double w2_,
                double w3_,
                double w4_):
                x1(x1_),
                x2(x2_),
                cellID(cellID_),
                submeshID(submeshID_),
                part_active(part_active_),
                exit_faceID(exit_faceID_),
                w1(w1_),
                w2(w2_),
                w3(w3_),
                w4(w4_){};

};

class ParticleDataCPU{
public:
    double x1;
    double x2;
    int cellID;
    int submeshID;
    bool part_active;
    int exit_faceID;
    double w1;
    double w2;
    double w3;
    double w4;

    ParticleDataCPU(){};

    ParticleDataCPU(double x1_,double x2_):
                x1(x1_),x2(x2_){};

    ParticleDataCPU(double x1_, double x2_, int submeshID_, int cellID_, bool part_active_, int exit_faceID_):
                x1(x1_),x2(x2_),cellID(cellID_),submeshID(submeshID_),part_active(part_active_),exit_faceID(exit_faceID_){};

    ParticleDataCPU(double x1_,
                    double x2_,
                    int submeshID_,
                    int cellID_,
                    bool part_active_,
                    int exit_faceID_,
                    double w1_,
                    double w2_,
                    double w3_,
                    double w4_):
                    x1(x1_),
                    x2(x2_),
                    cellID(cellID_),
                    submeshID(submeshID_),
                    part_active(part_active_),
                    exit_faceID(exit_faceID_),
                    w1(w1_),
                    w2(w2_),
                    w3(w3_),
                    w4(w4_){};

};

} // namespace pumi
#include "pumiMBBL_initiate.hpp"
#include "pumiMBBL_meshinfo.hpp"
#include "pumiMBBL_meshops.hpp"
#include "pumiMBBL_meshutils.hpp"
#include "pumiMBBLGPU_impl.hpp"
#include "pumiMBBL_finalize.hpp"
#endif
