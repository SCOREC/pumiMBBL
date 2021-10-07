#ifndef pumiMBBLGPU_initiate_hpp
#define pumiMBBLGPU_initiate_hpp

#include "pumiMBBLGPU.hpp"

namespace pumi {
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
    double domain_x1_min;
    double domain_x2_min;
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
* \brief enum type for printing node coordinates to terminal
*/
enum print_node_coords{
    print_node_coords_OFF = 0, //!< option to disable printing node coords
    print_node_coords_ON  = 1, //!< option to enable printing node coords
};

enum print_mesh_connectivity{
    print_mesh_connectivity_OFF = 0,
    print_mesh_connectivity_ON  = 1,
};

/*!
* \brief struct of mesh options
*/
struct Mesh_Options{
    store_BL_coords BL_storage_option;
    print_node_coords print_node_option;
    print_mesh_connectivity print_mesh_connectivity_option;
    /*!
    * \brief Struct default constructor
    */
    Mesh_Options():BL_storage_option(store_BL_coords_ON),
                   print_node_option(print_node_coords_OFF),
                   print_mesh_connectivity_option(print_mesh_connectivity_ON){};
};

class SubmeshInit{
public:
    SubmeshDeviceViewPtr submesh;
    SubmeshHostViewPtr host_submesh;

    SubmeshInit(){};

    SubmeshInit(SubmeshDeviceViewPtr submesh_, SubmeshHostViewPtr host_submesh_):
                submesh(submesh_), host_submesh(host_submesh_){};
};

unsigned int get_submesh_type(std::string meshtype_string);
double compute_grading_ratio(double BL_T, double BL_t0, int BL_Nel);
Mesh_Inputs* inputs_allocate();
void inputs_deallocate(Mesh_Inputs* pumi_inputs);
void print_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1);
void print_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2);
void print_mesh_nodes(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, Mesh_Options pumi_options);
void print_mesh_nodes(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2, Mesh_Options pumi_options);
bool verify_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1);
bool verify_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2);
SubmeshInit submesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, int dir);
Mesh mesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, SubmeshDeviceViewPtr submesh_x1, SubmeshHostViewPtr hc_submesh_x1);
Mesh mesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, SubmeshDeviceViewPtr submesh_x1, SubmeshHostViewPtr hc_submesh_x1,
                            SubmeshDeviceViewPtr submesh_x2, SubmeshHostViewPtr hc_submesh_x2);

MBBL initialize_MBBL_mesh(Mesh_Inputs* pumi_inputs, Mesh_Options pumi_options);

} // namespace pumi
#endif
