#ifndef pumiMBBL_initiate_hpp
#define pumiMBBL_initiate_hpp

#include "pumiMBBLGPU.hpp"

namespace pumi {
///////// Mesh-Initiate Function declarations ///////////////////////////////////////

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
Mesh_Inputs* inputs_allocate(){
    Mesh_Inputs* pumi_inputs = new Mesh_Inputs;
    // pumi_inputs->block_length_x1 = new double[nsubmeshes];
    // pumi_inputs->block_length_x2 = new double[nsubmeshes];
    // pumi_inputs->max_elem_size_x1 = new double[nsubmeshes];
    // pumi_inputs->max_elem_size_x2 = new double[nsubmeshes];
    // pumi_inputs->min_elem_size_x1 = new double[nsubmeshes];
    // pumi_inputs->min_elem_size_x2 = new double[nsubmeshes];

    return pumi_inputs;
}

/**
 * @brief Frees up the alocated memory of mesh inputs struct members
 *
 * \param[in] mesh inputs struct pointer
 */
void inputs_deallocate(Mesh_Inputs* pumi_inputs){
    // delete[] pumi_inputs->block_length_x1;
    // delete[] pumi_inputs->block_length_x2;
    // delete[] pumi_inputs->max_elem_size_x1;
    // delete[] pumi_inputs->max_elem_size_x2;
    // delete[] pumi_inputs->min_elem_size_x1;
    // delete[] pumi_inputs->min_elem_size_x2;
    delete pumi_inputs;
}

/**
 * @brief Prints the all relevant 1D-mesh details
 *
 * \param[in] mesh object pointer
 * \param[in] host copy of x1-submesh object pointer
 */
void print_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);

    printf("\n\nPUMI mesh parameter info [X1-Direction] :\n\n");
    printf("\tTotal elements along X1-direction = %d\n\n", h_pumi_mesh(0).Nel_tot_x1);
    for (int i=1; i<=h_pumi_mesh(0).nsubmesh_x1; i++){
        printf("\tSUBMESH %d  parameters:\n", i);
        printf("\n\t submesh-type   = ");
        if (h_submesh_x1[i].meshtype & minBL){
            printf("leftBL\n");
            printf("\t left_t0        = %2.4e \t [m] Cell size of first/leftmost cell in left BL block\n", h_submesh_x1[i].t0);
            printf("\t left_tN        = %2.4e \t [m] Cell size of last/rightmost cell in left BL block\n",
                (h_submesh_x1[i].t0)*pow(h_submesh_x1[i].r,h_submesh_x1[i].Nel-1));
            printf("\t left_T         = %2.4e \t [m] Left boundary layer (left BL) thickness\n", h_submesh_x1[i].length);
            printf("\t left_r         = %2.4e \t Grading ratio in left BL block\n", h_submesh_x1[i].r);
            printf("\t left_Nel       = %d    \t Number of Cells in left BL block\n\n", h_submesh_x1[i].Nel);
            // if (h_submesh_x1[i].BL_coords.extent(0)){
                printf("\t %d leftBL node coords stored in a array\n\n", h_submesh_x1[i].Nel+1);
            // }
        }
        if (h_submesh_x1[i].meshtype & maxBL){
            printf("rightBL\n");
            printf("\t right_t0       = %2.4e \t [m] Cell size of last/rightmost cell in right BL block\n", h_submesh_x1[i].t0);
            printf("\t right_tN       = %2.4e \t [m] Cell size of first/leftmost cell in right BL block\n",
            (h_submesh_x1[i].t0)*pow(h_submesh_x1[i].r,h_submesh_x1[i].Nel-1));
            printf("\t right_T        = %2.4e \t [m] Left boundary layer (right BL) thickness\n", h_submesh_x1[i].length);
            printf("\t right_r        = %2.4e \t Grading ratio in right BL mesh\n", h_submesh_x1[i].r);
            printf("\t right_Nel      = %d    \t Number of Cells in left BL mesh region\n\n", h_submesh_x1[i].Nel);
            // if (h_submesh_x1[i].BL_coords.extent(0)){
                printf("\t %d rightBL node coords stored in a array\n\n", h_submesh_x1[i].Nel+1);
            // }
        }
        if (h_submesh_x1[i].meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_dx1    = %2.4e \t [m] Cell size in uniform block\n", h_submesh_x1[i].t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", h_submesh_x1[i].Nel);
        }
    }
}

/**
 * @brief Prints the all relevant 2D-mesh details
 *
 * \param[in] mesh object pointer
 * \param[in] host copy of x1-submesh object pointer
 * \param[in] host copy of x2-submesh object pointer
 */
void print_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);

    printf("\n\nPUMI mesh parameter info [X1-Direction] :\n\n");
    printf("\tTotal elements along X1-direction = %d\n\n", h_pumi_mesh(0).Nel_tot_x1);
    for (int i=1; i<=h_pumi_mesh(0).nsubmesh_x1; i++){
        printf("\tSUBMESH %d  parameters:\n", i);
        printf("\n\t submesh-type   = ");
        if (h_submesh_x1[i].meshtype & minBL){
            printf("leftBL\n");
            printf("\t left_t0        = %2.4e \t [m] Cell size of first/leftmost cell in left BL block\n", h_submesh_x1[i].t0);
            printf("\t left_tN        = %2.4e \t [m] Cell size of last/rightmost cell in left BL block\n",
                (h_submesh_x1[i].t0)*pow(h_submesh_x1[i].r,h_submesh_x1[i].Nel-1));
            printf("\t left_T         = %2.4e \t [m] Left boundary layer (left BL) thickness\n", h_submesh_x1[i].length);
            printf("\t left_r         = %2.4e \t Grading ratio in left BL block\n", h_submesh_x1[i].r);
            printf("\t left_Nel       = %d    \t Number of Cells in left BL block\n\n", h_submesh_x1[i].Nel);
            if (h_submesh_x1[i].BL_coords.extent(0)){
                printf("\t %d leftBL node coords stored in a array\n\n", h_submesh_x1[i].Nel+1);
            }
        }
        if (h_submesh_x1[i].meshtype & maxBL){
            printf("rightBL\n");
            printf("\t right_t0       = %2.4e \t [m] Cell size of last/rightmost cell in right BL block\n", h_submesh_x1[i].t0);
            printf("\t right_tN       = %2.4e \t [m] Cell size of first/leftmost cell in right BL block\n",
            (h_submesh_x1[i].t0)*pow(h_submesh_x1[i].r,h_submesh_x1[i].Nel-1));
            printf("\t right_T        = %2.4e \t [m] Left boundary layer (right BL) thickness\n", h_submesh_x1[i].length);
            printf("\t right_r        = %2.4e \t Grading ratio in right BL mesh\n", h_submesh_x1[i].r);
            printf("\t right_Nel      = %d    \t Number of Cells in left BL mesh region\n\n", h_submesh_x1[i].Nel);
            if (h_submesh_x1[i].BL_coords.extent(0)){
                printf("\t %d rightBL node coords stored in a array\n\n", h_submesh_x1[i].Nel+1);
            }
        }
        if (h_submesh_x1[i].meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_dx1    = %2.4e \t [m] Cell size in uniform block\n", h_submesh_x1[i].t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", h_submesh_x1[i].Nel);
        }
    }

    printf("PUMI mesh parameter info [X2-Direction] :\n\n");
    printf("\tTotal elements along X2-direction = %d\n\n", h_pumi_mesh(0).Nel_tot_x2);
    for (int i=1; i<=h_pumi_mesh(0).nsubmesh_x2; i++){
        printf("\tSUBMESH %d  parameters:\n", i);
        printf("\n\t submesh-type   = ");
        if (h_submesh_x2[i].meshtype & minBL){
            printf("bottomBL\n");
            printf("\t bottom_t0      = %2.4e \t [m] Cell size of first/bottom-most cell in left BL block\n", h_submesh_x2[i].t0);
            printf("\t bottom_tN      = %2.4e \t [m] Cell size of last /   top-most cell in left BL block\n",
                (h_submesh_x2[i].t0)*pow(h_submesh_x2[i].r,h_submesh_x2[i].Nel-1));
            printf("\t bottom_T       = %2.4e \t [m] Left boundary layer (left BL) thickness\n", h_submesh_x2[i].length);
            printf("\t bottom_r       = %2.4e \t Grading ratio in left BL block\n", h_submesh_x2[i].r);
            printf("\t bottom_Nel     = %d    \t Number of Cells in left BL block\n\n", h_submesh_x2[i].Nel);
            if (h_submesh_x2[i].BL_coords.extent(0)){
                printf("\t %d bottomBL node coords stored in a array\n\n", h_submesh_x2[i].Nel+1);
            }
        }
        if (h_submesh_x2[i].meshtype & maxBL){
            printf("topBL\n");
            printf("\t top_t0         = %2.4e \t [m] Cell size of last /   top-most cell in right BL block\n", h_submesh_x2[i].t0);
            printf("\t top_tN         = %2.4e \t [m] Cell size of first/bottom-most cell in right BL block\n",
            (h_submesh_x2[i].t0)*pow(h_submesh_x2[i].r,h_submesh_x2[i].Nel-1));

            printf("\t top_T          = %2.4e \t [m] Left boundary layer (right BL) thickness\n", h_submesh_x2[i].length);
            printf("\t top_r          = %2.4e \t Grading ratio in right BL mesh\n", h_submesh_x2[i].r);
            printf("\t top_Nel        = %d    \t Number of Cells in left BL mesh region\n\n", h_submesh_x2[i].Nel);
            if (h_submesh_x2[i].BL_coords.extent(0)){
                printf("\t %d topBL node coords stored in a array\n\n", h_submesh_x2[i].Nel+1);
            }
        }
        if (h_submesh_x2[i].meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_dx2    = %2.4e \t [m] Cell size in uniform block\n", h_submesh_x2[i].t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", h_submesh_x2[i].Nel);
        }
    }

    printf("PUMI submesh activity info :\n\n");
    for (int jsubmesh=h_pumi_mesh(0).nsubmesh_x2; jsubmesh>=1; jsubmesh--){
        if (jsubmesh != h_pumi_mesh(0).nsubmesh_x2){
            for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x1-1; isubmesh++ ){
                printf("_____________");
            }
            printf("__________________\n\n");
        }
        for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x1; isubmesh++ ){
            if (isubmesh-1){
                printf("|");
            }
            if(h_pumi_mesh(0).host_isactive[isubmesh][jsubmesh]){
                printf("    ACTIVE    ");
            }
            else{
                printf("    XXXXXX    ");
            }
        }
        printf("\n");
    }
    printf("\n\n");

    printf("Total active elements in 2D Mesh = %d\n",h_pumi_mesh(0).Nel_total);
    printf("Total active nodes in 2D Mesh = %d\n",h_pumi_mesh(0).Nnp_total);

}

/**
 * @brief Prints node coordinates of 1D mesh to the terminal and also
 * writes the individual submesh coords to a file as well as
 * full mesh coordinates for all directions
 * \param[in] mesh object pointer
 * \param[in] CPU copy of x1-submesh object pointer
 */
void print_mesh_nodes(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1, Mesh_Options pumi_options){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);
    if (pumi_options.print_node_options)
        printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
    FILE *mesh_coords_file;
    char mesh_coords_filename[30];
    sprintf(mesh_coords_filename,"X1_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    int inode=0;
    for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
        if (pumi_options.print_node_options)
            printf("X1-SUBMESH %d:\n", isubmesh );
        FILE *submesh_coords_file;
        char submesh_coords_filename[33];
        sprintf(submesh_coords_filename,"X1_submesh_%d_coords.dat",isubmesh);
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
        if (pumi_options.print_node_options)
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x1[isubmesh].Nel; icell++){
            inode++;
            icoord += cell_size;
            cell_size *= r;
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            if (pumi_options.print_node_options)
                printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        }
        printf("\n\t X1-SUBMESH %d -- Coordinates written to file %s\n\n",isubmesh, submesh_coords_filename);
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
void print_mesh_nodes(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2, Mesh_Options pumi_options){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);
    if (pumi_options.print_node_options)
        printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
    FILE *mesh_coords_file;
    char mesh_coords_filename[30];
    sprintf(mesh_coords_filename,"X1_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    int inode=0;
    for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
        if (pumi_options.print_node_options)
            printf("X1-SUBMESH %d:\n", isubmesh );
        FILE *submesh_coords_file;
        char submesh_coords_filename[33];
        sprintf(submesh_coords_filename,"X1_submesh_%d_coords.dat",isubmesh);
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
        if (pumi_options.print_node_options)
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x1[isubmesh].Nel; icell++){
            inode++;
            icoord += cell_size;
            cell_size *= r;
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            if (pumi_options.print_node_options)
                printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        }
        printf("\n\t X1-SUBMESH %d -- Coordinates written to file %s\n\n", isubmesh, submesh_coords_filename);
        fclose(submesh_coords_file);
    }
    fclose(mesh_coords_file);

    sprintf(mesh_coords_filename,"X2_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    inode = 0;
    for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x2; isubmesh++){
        if (pumi_options.print_node_options)
            printf("X2-SUBMESH %d:\n", isubmesh );
        FILE *submesh_coords_file;
        char submesh_coords_filename[33];
        sprintf(submesh_coords_filename,"X2_submesh_%d_coords.dat",isubmesh);
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
        if (pumi_options.print_node_options)
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x2[isubmesh].Nel; icell++){
            inode++;
            icoord += cell_size;
            cell_size *= r;
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            if (pumi_options.print_node_options)
                printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        }
        printf("\n\t X2-SUBMESH %d -- Coordinates written to file %s\n\n",isubmesh, submesh_coords_filename);
        fclose(submesh_coords_file);
    }
    fclose(mesh_coords_file);
}

/**
 * @brief Verifies the computed mesh and submesh parameters and returns true if verified
 *
 * \param[in] mesh object pointer
 * \param[in] host copy of x1-submesh object pointer
 */
bool verify_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);

    printf("\n\nNow verifying valdity of pumi mesh parameters for\n");
    int flag = 0;
    for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
        printf("\tX1-SUBMESH %d:\n", isubmesh );
        if (h_submesh_x1[isubmesh].meshtype & minBL){
            if (!(h_submesh_x1[isubmesh].Nel > 0)){
                printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t left_Nel    -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh].r >= 1.0)){
                printf("\t\t left_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh].r);
                flag++;
            }
            else{
                printf("\t\t left_r      -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh].t0*h_submesh_x1[isubmesh].Nel;
            if (!(h_submesh_x1[isubmesh].length > 0.0)){
                printf("\t\t left_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh].length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh].length > min_computed_T)){
                printf("\t\t left_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh].length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t left_T      -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh].meshtype & maxBL){
            if (!(h_submesh_x1[isubmesh].Nel > 0)){
                printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t right_Nel   -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh].r >= 1.0)){
                printf("\t\t right_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh].r);
                flag++;
            }
            else{
                printf("\t\t right_r     -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh].t0*h_submesh_x1[isubmesh].Nel;
            if (!(h_submesh_x1[isubmesh].length > 0.0)){
                printf("\t\t right_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh].length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh].length > min_computed_T)){
                printf("\t\t right_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh].length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t right_T     -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh].meshtype & uniform){
            if (!(h_submesh_x1[isubmesh].Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
    }

    if (flag == 0){
        printf("\n\tThe input mesh parameters and the calculated mesh parameters are all valid and verified\n\n");
        return true;
    }
    else{
        printf("\t\nERROR: One or more input/calculated mesh paramater is not valid. Abort\n");
        return false;
    }
}

/**
 * @brief Verifies the computed mesh and submesh parameters and returns true if verified
 *
 * \param[in] mesh object pointer
 * \param[in] host copy of x1-submesh object pointer
 * \param[in] host copy of x2-submesh object pointer
 */
bool verify_mesh_params(MeshDeviceViewPtr pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2){
    MeshDeviceViewPtr::HostMirror h_pumi_mesh = Kokkos::create_mirror_view(pumi_mesh);
    Kokkos::deep_copy(h_pumi_mesh, pumi_mesh);

    printf("\n\nNow verifying valdity of pumi mesh parameters for\n");
    int flag = 0;
    for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x1; isubmesh++){
        printf("\tX1-SUBMESH %d:\n", isubmesh );
        if (h_submesh_x1[isubmesh].meshtype & minBL){
            if (!(h_submesh_x1[isubmesh].Nel > 0)){
                printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t left_Nel    -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh].r >= 1.0)){
                printf("\t\t left_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh].r);
                flag++;
            }
            else{
                printf("\t\t left_r      -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh].t0*h_submesh_x1[isubmesh].Nel;
            if (!(h_submesh_x1[isubmesh].length > 0.0)){
                printf("\t\t left_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh].length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh].length > min_computed_T)){
                printf("\t\t left_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh].length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t left_T      -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh].meshtype & maxBL){
            if (!(h_submesh_x1[isubmesh].Nel > 0)){
                printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t right_Nel   -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh].r >= 1.0)){
                printf("\t\t right_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh].r);
                flag++;
            }
            else{
                printf("\t\t right_r     -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh].t0*h_submesh_x1[isubmesh].Nel;
            if (!(h_submesh_x1[isubmesh].length > 0.0)){
                printf("\t\t right_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh].length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh].length > min_computed_T)){
                printf("\t\t right_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh].length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t right_T     -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh].meshtype & uniform){
            if (!(h_submesh_x1[isubmesh].Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
    }


    for (int isubmesh=1; isubmesh<=h_pumi_mesh(0).nsubmesh_x2; isubmesh++){
        printf("\tX2-SUBMESH %d:\n", isubmesh );
        if (h_submesh_x2[isubmesh].meshtype & minBL){
            if (!(h_submesh_x2[isubmesh].Nel > 0)){
                printf("\t\t bottom_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x2[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t bottom_Nel  -- verified...\n");
            }
            if (!(h_submesh_x2[isubmesh].r >= 1.0)){
                printf("\t\t bottom_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x2[isubmesh].r);
                flag++;
            }
            else{
                printf("\t\t bottom_r    -- verified...\n");
            }
            double min_computed_T = h_submesh_x2[isubmesh].t0*h_submesh_x2[isubmesh].Nel;
            if (!(h_submesh_x2[isubmesh].length > 0.0)){
                printf("\t\t bottom_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x2[isubmesh].length);
                flag++;
            }
            else if (!(h_submesh_x2[isubmesh].length > min_computed_T)){
                printf("\t\t bottom_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x2[isubmesh].length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t bottom_T    -- verified...\n");
            }
        }
        if (h_submesh_x2[isubmesh].meshtype & maxBL){
            if (!(h_submesh_x2[isubmesh].Nel > 0)){
                printf("\t\t top_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x2[isubmesh].Nel);
                flag++;
            }
            else{
                printf("\t\t top_Nel     -- verified...\n");
            }
            if (!(h_submesh_x2[isubmesh].r >= 1.0)){
                printf("\t\t top_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x2[isubmesh].r);
                flag++;
            }
            else{
                printf("\t\t top_r       -- verified...\n");
            }
            double min_computed_T = h_submesh_x2[isubmesh].t0*h_submesh_x2[isubmesh].Nel;
            if (!(h_submesh_x2[isubmesh].length > 0.0)){
                printf("\t\t top_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x2[isubmesh].length);
                flag++;
            }
            else if (!(h_submesh_x2[isubmesh].length > min_computed_T)){
                printf("\t\t top_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x2[isubmesh].length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t top_T       -- verified...\n");
            }
        }
        if (h_submesh_x2[isubmesh].meshtype & uniform){
            if (!(h_submesh_x2[isubmesh].Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x2[isubmesh].Nel);
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
    }


    if (flag == 0){
        printf("\n\tThe input mesh parameters and the calculated mesh parameters are all valid and verified\n\n");
        return true;
    }
    else{
        printf("\t\nERROR: One or more input/calculated mesh paramater is not valid. Abort\n");
        return false;
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
    int nsubmesh = 0;
    SubmeshDeviceViewPtr submesh;
    SubmeshDeviceViewPtr::HostMirror h_submesh;

    if (dir == x1_dir){
        nsubmesh = pumi_inputs->nsubmesh_x1;
        submesh = SubmeshDeviceViewPtr("x1-submesh-obj",nsubmesh+2);
        h_submesh = Kokkos::create_mirror_view(submesh);
    }
    else if (dir == x2_dir){
        nsubmesh = pumi_inputs->nsubmesh_x2;
        submesh = SubmeshDeviceViewPtr("x2-submesh-obj",nsubmesh+2);
        h_submesh = Kokkos::create_mirror_view(submesh);
    }

    SubmeshHostViewPtr submesh_host_copy = new Submesh[nsubmesh+2];

    if (pumi_options.BL_storage_option == store_BL_coords_OFF){
        pumi_options.BL_storage_option = store_BL_coords_ON;
    }

    double xmin, xmax, xlength, t0, tN, r, r_t0_ratio, logr;
    int Nel = 0;
    unsigned int type = 0x00;
    int Nel_cumulative = 0;
    xmin = 0.0;
    xmax = 0.0;
    xlength = 0.0;
    t0 = 0.0;
    tN = 0.0;
    r = 0.0;
    r_t0_ratio = 0.0;
    logr = 0.0;

    if (dir == x1_dir){
        double total_length = 0.0;
        for (int isubmesh=0; isubmesh<nsubmesh; isubmesh++){
            total_length += pumi_inputs->block_length_x1[isubmesh];
        }

        type = unassigned;
        xlength = 1000.0*total_length;
        DoubleViewPtr BLcoords;
        // padding submesh to min-side
        xmax = 0.0;
        xmin = xmax - xlength;
        Unassigned_Submesh tmp_obj_min(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords);
        submesh_host_copy[0] = tmp_obj_min;
        h_submesh(0) = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_min);
        // padding submesh to max-side
        xmin = total_length;
        xmax = xmin + xlength;
        Unassigned_Submesh tmp_obj_max(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords);
        submesh_host_copy[nsubmesh+1] = tmp_obj_max;
        h_submesh(nsubmesh+1) = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_max);
    }
    else if (dir == x2_dir){
        double total_length = 0.0;
        for (int isubmesh=0; isubmesh<nsubmesh; isubmesh++){
            total_length += pumi_inputs->block_length_x2[isubmesh];
        }
        type = unassigned;
        xlength = 1000.0*total_length;
        DoubleViewPtr BLcoords;
        // padding submesh to min-side
        xmax = 0.0;
        xmin = xmax - xlength;
        Unassigned_Submesh tmp_obj_min(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords);
        submesh_host_copy[0] = tmp_obj_min;
        h_submesh(0) = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_min);
        // padding submesh to max-side
        xmin = total_length;
        xmax = xmin + xlength;
        Unassigned_Submesh tmp_obj_max(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords);
        submesh_host_copy[nsubmesh+1] = tmp_obj_max;
        h_submesh(nsubmesh+1) = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_max);
    }


    for (int isubmesh=1; isubmesh<=nsubmesh; isubmesh++){
        DoubleViewPtr BLcoords;
        std::string BLcoordsname;

        if (dir == x1_dir){
            type = get_submesh_type(pumi_inputs->meshtype_x1[isubmesh-1]);

            xlength = pumi_inputs->block_length_x1[isubmesh-1];
            if (isubmesh==1){
                xmin = 0.0;
            }
            else{
                xmin = xmax;
            }
            xmax = xmin + xlength;
            t0 = pumi_inputs->min_elem_size_x1[isubmesh-1];
            tN = pumi_inputs->max_elem_size_x1[isubmesh-1];

            BLcoordsname = "x1-BLcoords-submesh-";
            BLcoordsname += std::to_string(isubmesh);
        }
        else if (dir == x2_dir){
            type = get_submesh_type(pumi_inputs->meshtype_x2[isubmesh-1]);

            xlength = pumi_inputs->block_length_x2[isubmesh-1];
            if (isubmesh==1){
                xmin = 0.0;
            }
            else{
                xmin = xmax;
            }
            xmax = xmin + xlength;
            t0 = pumi_inputs->min_elem_size_x2[isubmesh-1];
            tN = pumi_inputs->max_elem_size_x2[isubmesh-1];

            BLcoordsname = "x2-BLcoords-submesh-";
            BLcoordsname += std::to_string(isubmesh);
        }

        if (isubmesh > 1){
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
MeshDeviceViewPtr mesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, SubmeshDeviceViewPtr submesh_x1, SubmeshHostViewPtr hc_submesh_x1){
    MeshDeviceViewPtr pumi_mesh("1D-meshobj",1);

    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int Nel_tot_x1;

    // Kokkos::View<int*> Nel_total_x1("total-x1-elems",1);
    // Kokkos::View<int*>::HostMirror h_Nel_total_x1 = Kokkos::create_mirror_view(Nel_total_x1);
    //
    // Kokkos::parallel_for("1D-total-nel-init", 1, KOKKOS_LAMBDA (const int) {
    //     Nel_total_x1(0) = submesh_x1(nsubmesh_x1)()->Nel + submesh_x1(nsubmesh_x1)()->Nel_cumulative;
    // });
    // Kokkos::deep_copy(h_Nel_total_x1,Nel_total_x1);
    // Nel_tot_x1 = h_Nel_total_x1(0);
    Nel_tot_x1 = hc_submesh_x1[nsubmesh_x1].Nel + hc_submesh_x1[nsubmesh_x1].Nel_cumulative;

    Kokkos::parallel_for("1D-meshobj-init", 1, KOKKOS_LAMBDA (const int) {
        pumi_mesh(0) = Mesh(nsubmesh_x1, Nel_tot_x1);
    });
    print_mesh_params(pumi_mesh, hc_submesh_x1);

    bool mesh_verified = verify_mesh_params(pumi_mesh, hc_submesh_x1);
    if (!mesh_verified){
        Kokkos::finalize();
        exit(0);
    }
    else{
        print_mesh_nodes(pumi_mesh, hc_submesh_x1, pumi_options);
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
MeshDeviceViewPtr mesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, SubmeshDeviceViewPtr submesh_x1, SubmeshHostViewPtr hc_submesh_x1,
                            SubmeshDeviceViewPtr submesh_x2, SubmeshHostViewPtr hc_submesh_x2){
    MeshDeviceViewPtr pumi_mesh("2D-meshobj",1);

    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int nsubmesh_x2 = pumi_inputs->nsubmesh_x2;
    int Nel_tot_x1, Nel_tot_x2, Nel_total_2D, Nnp_total_2D;

    // Kokkos::View<int*> Nel_total_x1("total-x1-elems",1);
    // Kokkos::View<int*>::HostMirror h_Nel_total_x1 = Kokkos::create_mirror_view(Nel_total_x1);
    // Kokkos::View<int*> Nel_total_x2("total-x2-elems",1);
    // Kokkos::View<int*>::HostMirror h_Nel_total_x2 = Kokkos::create_mirror_view(Nel_total_x2);
    // Kokkos::parallel_for("2D-total-nel-init", 1, KOKKOS_LAMBDA (const int) {
    //     Nel_total_x1(0) = submesh_x1(nsubmesh_x1)()->Nel + submesh_x1(nsubmesh_x1)()->Nel_cumulative;
    //     Nel_total_x2(0) = submesh_x2(nsubmesh_x2)()->Nel + submesh_x2(nsubmesh_x2)()->Nel_cumulative;
    // });
    // Kokkos::deep_copy(h_Nel_total_x1,Nel_total_x1);
    // Kokkos::deep_copy(h_Nel_total_x2,Nel_total_x2);
    // Nel_tot_x1 = h_Nel_total_x1(0);
    // Nel_tot_x2 = h_Nel_total_x2(0);
    Nel_tot_x1 = hc_submesh_x1[nsubmesh_x1].Nel + hc_submesh_x1[nsubmesh_x1].Nel_cumulative;
    Nel_tot_x2 = hc_submesh_x2[nsubmesh_x2].Nel + hc_submesh_x2[nsubmesh_x2].Nel_cumulative;
    Nel_total_2D = Nel_tot_x1*Nel_tot_x2;
    Nnp_total_2D = (Nel_tot_x1+1)*(Nel_tot_x2+1);

    Kokkos::View<bool**> submesh_activity("submesh-isactive",nsubmesh_x1+2,nsubmesh_x2+2);
    bool** host_isactive = new bool*[nsubmesh_x1+2];
    for (int i=0; i<nsubmesh_x1+2; i++){
        host_isactive[i] = new bool[nsubmesh_x2+2];
    }
    for (int i=0; i<nsubmesh_x1+2; i++){
        for (int j=0; j<nsubmesh_x2+2; j++){
            if (i==0 || i==nsubmesh_x1+1 || j==0 || j==nsubmesh_x2+1){
                host_isactive[i][j] = false;
            }
            else{
                host_isactive[i][j] = pumi_inputs->isactive[i-1][j-1];
            }
        }
    }

    Kokkos::View<bool**>::HostMirror h_submesh_activity = Kokkos::create_mirror_view(submesh_activity);
    for (int isubmesh=0; isubmesh<nsubmesh_x1+2; isubmesh++ ){
        for (int jsubmesh=0; jsubmesh<nsubmesh_x2+2; jsubmesh++){
            if (isubmesh==0 || isubmesh==nsubmesh_x1+1 || jsubmesh==0 || jsubmesh==nsubmesh_x2+1){
                h_submesh_activity(isubmesh,jsubmesh) = false;
            }
            else{
                h_submesh_activity(isubmesh,jsubmesh) = pumi_inputs->isactive[isubmesh-1][jsubmesh-1];
            }
        }
    }
    Kokkos::deep_copy(submesh_activity, h_submesh_activity);

    bool is_fullmesh = true;
    for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++ ){
        for (int jsubmesh=1; jsubmesh<=nsubmesh_x2; jsubmesh++){
            if (!(h_submesh_activity(isubmesh, jsubmesh))){
                is_fullmesh = false;
            }
        }
    }

    Kokkos::View<int**> elemoffset_start("elemoffset_start", nsubmesh_x1+2, nsubmesh_x2+2);
    Kokkos::View<int*> elemoffset_skip("elemoffset_skip", nsubmesh_x2+2);
    // Kokkos::View<int**> nodeoffset("nodeoffset", nsubmesh_x1+2, Nel_tot_x2+1);
    Kokkos::View<int**> nodeoffset_start("nodeoffset_start", nsubmesh_x1+2, nsubmesh_x2+2);
    Kokkos::View<int**> nodeoffset_skip_bot("nodeoffset_skip_bot", nsubmesh_x1+2, nsubmesh_x2+2);
    Kokkos::View<int**> nodeoffset_skip_mid("nodeoffset_skip_mid", nsubmesh_x1+2, nsubmesh_x2+2);
    Kokkos::View<int**> nodeoffset_skip_top("nodeoffset_skip_top", nsubmesh_x1+2, nsubmesh_x2+2);

    int **nodeoffset = new int*[nsubmesh_x1+2];
    for (int i=0; i<nsubmesh_x1+2; i++){
        nodeoffset[i] = new int[Nel_tot_x2+1];
    }

    if (!is_fullmesh){
        Kokkos::View<int**>::HostMirror h_elemoffset_start = Kokkos::create_mirror_view(elemoffset_start);
        Kokkos::View<int*>::HostMirror h_elemoffset_skip = Kokkos::create_mirror_view(elemoffset_skip);
        // Kokkos::View<int**>::HostMirror h_nodeoffset = Kokkos::create_mirror_view(nodeoffset);
        Kokkos::View<int**>::HostMirror h_nodeoffset_start = Kokkos::create_mirror_view(nodeoffset_start);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_bot = Kokkos::create_mirror_view(nodeoffset_skip_bot);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_mid = Kokkos::create_mirror_view(nodeoffset_skip_mid);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_top = Kokkos::create_mirror_view(nodeoffset_skip_top);

        // elemoffsets
        int elemstart_init, elemskip;
        int elemstart = 0;

        for (int jsubmesh=1; jsubmesh<=nsubmesh_x2; jsubmesh++){
            elemstart_init = elemstart;

            for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
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
        int jsubmesh = 1;
        int nodestart = 0;
        for (jnp=0; jnp<hc_submesh_x2[jsubmesh].Nel; jnp++){
            int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;
            for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
                if(h_submesh_activity(isubmesh, jsubmesh)){
                    nodeoffset[isubmesh][Jnp] = nodestart;
                }
                else{
                    nodeoffset[isubmesh][Jnp] = -1;
                    if (isubmesh==1){
                        if (h_submesh_activity(isubmesh+1,jsubmesh)){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                        }
                    }
                    else if (isubmesh == nsubmesh_x1){
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
        for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
            if (jsubmesh == nsubmesh_x2){
                int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;
                if (h_submesh_activity(isubmesh, jsubmesh)){
                    nodeoffset[isubmesh][Jnp] = nodestart;
                }
                else{
                    nodeoffset[isubmesh][Jnp] = -1;
                    if (isubmesh==1){
                        if (h_submesh_activity(isubmesh+1,jsubmesh)){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                        }
                    }
                    else if (isubmesh==nsubmesh_x1){
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

        for (jsubmesh=2; jsubmesh<=nsubmesh_x2; jsubmesh++){
            jnp = 0;
            int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;

            for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
                if (h_submesh_activity(isubmesh,jsubmesh) || h_submesh_activity(isubmesh,jsubmesh-1)){
                    nodestart += 0;
                    nodeoffset[isubmesh][Jnp] = nodestart;
                }
                else{
                    nodeoffset[isubmesh][Jnp] = -1;
                    if (isubmesh==1){
                        if (h_submesh_activity(isubmesh+1,jsubmesh) || h_submesh_activity(isubmesh+1,jsubmesh-1)){
                            nodestart += hc_submesh_x1[isubmesh].Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                        }
                    }
                    else if (isubmesh==nsubmesh_x1){
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

                for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        nodeoffset[isubmesh][Jnp] = nodestart;
                    }
                    else{
                        nodeoffset[isubmesh][Jnp] = -1;
                        if (isubmesh==1){
                            if (h_submesh_activity(isubmesh+1,jsubmesh)){
                                nodestart += hc_submesh_x1[isubmesh].Nel;
                            }
                            else{
                                nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                            }
                        }
                        else if (isubmesh == nsubmesh_x1){
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

            if (jsubmesh==nsubmesh_x2){
                jnp = hc_submesh_x2[jsubmesh].Nel;
                int Jnp = jnp + hc_submesh_x2[jsubmesh].Nel_cumulative;

                for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        nodeoffset[isubmesh][Jnp] = nodestart;
                    }
                    else{
                        nodeoffset[isubmesh][Jnp] = -1;
                        if (isubmesh==1){
                            if (h_submesh_activity(isubmesh+1,jsubmesh)){
                                nodestart += hc_submesh_x1[isubmesh].Nel;
                            }
                            else{
                                nodestart += (hc_submesh_x1[isubmesh].Nel+1);
                            }
                        }
                        else if (isubmesh==nsubmesh_x1){
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

        for (int jsubmesh=1; jsubmesh<=nsubmesh_x2; jsubmesh++){
            int Jnp = hc_submesh_x2[jsubmesh].Nel_cumulative;
            for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
                if (h_submesh_activity(isubmesh,jsubmesh)){
                    h_nodeoffset_start(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp];
                    h_nodeoffset_skip_bot(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp+1]-nodeoffset[isubmesh][Jnp];
                    h_nodeoffset_skip_mid(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp+2]-nodeoffset[isubmesh][Jnp+1];
                    h_nodeoffset_skip_top(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp+hc_submesh_x2[jsubmesh].Nel] -
                                                                nodeoffset[isubmesh][Jnp+hc_submesh_x2[jsubmesh].Nel-1];
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
        // Kokkos::deep_copy(nodeoffset, h_nodeoffset);
        Kokkos::deep_copy(nodeoffset_start, h_nodeoffset_start);
        Kokkos::deep_copy(nodeoffset_skip_bot, h_nodeoffset_skip_bot);
        Kokkos::deep_copy(nodeoffset_skip_mid, h_nodeoffset_skip_mid);
        Kokkos::deep_copy(nodeoffset_skip_top, h_nodeoffset_skip_top);

    }
    else{
        Kokkos::View<int**>::HostMirror h_elemoffset_start = Kokkos::create_mirror_view(elemoffset_start);
        Kokkos::View<int*>::HostMirror h_elemoffset_skip = Kokkos::create_mirror_view(elemoffset_skip);
        // Kokkos::View<int**>::HostMirror h_nodeoffset = Kokkos::create_mirror_view(nodeoffset);
        Kokkos::View<int**>::HostMirror h_nodeoffset_start = Kokkos::create_mirror_view(nodeoffset_start);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_bot = Kokkos::create_mirror_view(nodeoffset_skip_bot);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_mid = Kokkos::create_mirror_view(nodeoffset_skip_mid);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_top = Kokkos::create_mirror_view(nodeoffset_skip_top);

        for (int jsubmesh=0; jsubmesh<nsubmesh_x2+2; jsubmesh++){
            h_elemoffset_skip(jsubmesh) = 0;
            for (int isubmesh=0; isubmesh<nsubmesh_x1+2; isubmesh++){
                h_elemoffset_start(isubmesh,jsubmesh) = 0;
                h_nodeoffset_start(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_bot(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_mid(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_top(isubmesh,jsubmesh) = 0;
            }
        }

        Kokkos::deep_copy(elemoffset_skip, h_elemoffset_skip);
        Kokkos::deep_copy(elemoffset_start, h_elemoffset_start);
        // Kokkos::deep_copy(nodeoffset, h_nodeoffset);
        Kokkos::deep_copy(nodeoffset_start, h_nodeoffset_start);
        Kokkos::deep_copy(nodeoffset_skip_bot, h_nodeoffset_skip_bot);
        Kokkos::deep_copy(nodeoffset_skip_mid, h_nodeoffset_skip_mid);
        Kokkos::deep_copy(nodeoffset_skip_top, h_nodeoffset_skip_top);
    }

    Kokkos::View<bool*> is_bdry("is_bdry", 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2);
    Kokkos::View<double*[3]> bdry_normal("bdry_normal", 2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2);
    Kokkos::View<bool*>::HostMirror h_is_bdry = Kokkos::create_mirror_view(is_bdry);
    Kokkos::View<double*[3]>::HostMirror h_bdry_normal = Kokkos::create_mirror_view(bdry_normal);

    for (int jsubmesh=1; jsubmesh<=nsubmesh_x2; jsubmesh++){
        for (int isubmesh=1; isubmesh<=nsubmesh_x1; isubmesh++){
            if (jsubmesh==1){
                if (h_submesh_activity(isubmesh,jsubmesh)){
                    h_is_bdry(isubmesh-1) = true;
                    h_bdry_normal(isubmesh-1,0) = 0.0;
                    h_bdry_normal(isubmesh-1,1) = -1.0;
                    h_bdry_normal(isubmesh-1,2) = 0.0;
                    if (isubmesh==1){
                        h_is_bdry(nsubmesh_x1) = true;
                        h_bdry_normal(nsubmesh_x1,0) = -1.0;
                        h_bdry_normal(nsubmesh_x1,1) = 0.0;
                        h_bdry_normal(nsubmesh_x1,2) = 0.0;
                    }
                    if (isubmesh==nsubmesh_x1){
                        h_is_bdry(isubmesh-1+nsubmesh_x1+1) = true;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1+1,0) = 1.0;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1+1,1) = 0.0;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1+1,2) = 0.0;
                    }
                }
                else{
                    h_is_bdry(isubmesh-1) = false;
                    h_bdry_normal(isubmesh-1,0) = 0.0;
                    h_bdry_normal(isubmesh-1,1) = 0.0;
                    h_bdry_normal(isubmesh-1,2) = 0.0;
                    if (isubmesh==1){
                        h_is_bdry(nsubmesh_x1) = false;
                        h_bdry_normal(nsubmesh_x1,0) = 0.0;
                        h_bdry_normal(nsubmesh_x1,1) = 0.0;
                        h_bdry_normal(nsubmesh_x1,2) = 0.0;
                    }
                    if (isubmesh==nsubmesh_x1){
                        h_is_bdry(isubmesh-1+nsubmesh_x1+1) = false;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1+1,0) = 0.0;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1+1,1) = 0.0;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1+1,2) = 0.0;
                    }
                }

                if (isubmesh>1){
                    int sum = h_submesh_activity(isubmesh-1,jsubmesh)+h_submesh_activity(isubmesh,jsubmesh);
                    if (sum == 1){
                        h_is_bdry(isubmesh-1+nsubmesh_x1) = true;
                        if (h_submesh_activity(isubmesh,jsubmesh)){
                            h_bdry_normal(isubmesh-1+nsubmesh_x1,0) = -1.0;
                            h_bdry_normal(isubmesh-1+nsubmesh_x1,1) = 0.0;
                            h_bdry_normal(isubmesh-1+nsubmesh_x1,2) = 0.0;
                        }
                        else{
                            h_bdry_normal(isubmesh-1+nsubmesh_x1,0) = 1.0;
                            h_bdry_normal(isubmesh-1+nsubmesh_x1,1) = 0.0;
                            h_bdry_normal(isubmesh-1+nsubmesh_x1,2) = 0.0;
                        }
                    }
                    else{
                        h_is_bdry(isubmesh-1+nsubmesh_x1) = false;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1,0) = 0.0;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1,1) = 0.0;
                        h_bdry_normal(isubmesh-1+nsubmesh_x1,2) = 0.0;
                    }
                }
            }
            else{

                int sum = h_submesh_activity(isubmesh,jsubmesh-1)+h_submesh_activity(isubmesh,jsubmesh);
                if (sum == 1){
                    h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1) = true;
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,0) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,1) = -1.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,2) = 0.0;
                    }
                    else{
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,0) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,1) = 1.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,2) = 0.0;
                    }
                }
                else{
                    h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1) = false;
                    h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,0) = 0.0;
                    h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,1) = 0.0;
                    h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1,2) = 0.0;
                }

                if (isubmesh==1){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1) = true;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1,0) = -1.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1,1) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1,2) = 0.0;
                    }
                    else{
                        h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1) = false;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1,0) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1,1) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1,2) = 0.0;
                    }
                }
                else{
                    sum = h_submesh_activity(isubmesh-1,jsubmesh)+h_submesh_activity(isubmesh,jsubmesh);
                    if (sum == 1){
                        h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1) = true;
                        if (h_submesh_activity(isubmesh,jsubmesh)){
                            h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,0) = -1.0;
                            h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,1) = 0.0;
                            h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,2) = 0.0;
                        }
                        else{
                            h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,0) = 1.0;
                            h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,1) = 0.0;
                            h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,2) = 0.0;
                        }
                    }
                    else{
                        h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1) = false;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,0) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,1) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+nsubmesh_x1+isubmesh-1,2) = 0.0;
                    }
                }

                if (isubmesh==nsubmesh_x1){
                    if (h_submesh_activity(isubmesh,jsubmesh)){
                        h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1) = true;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1,0) = 1.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1,1) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1,2) = 0.0;
                    }
                    else{
                        h_is_bdry((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1) = false;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1,0) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1,1) = 0.0;
                        h_bdry_normal((jsubmesh-1)*(2*nsubmesh_x1+1)+isubmesh-1+nsubmesh_x1+1,2) = 0.0;
                    }
                }
            }

            if (jsubmesh == nsubmesh_x2){
                if (h_submesh_activity(isubmesh,jsubmesh)){
                    h_is_bdry((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1) = true;
                    h_bdry_normal((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1,0) = 0.0;
                    h_bdry_normal((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1,1) = 1.0;
                    h_bdry_normal((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1,2) = 0.0;
                }
                else{
                    h_is_bdry((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1) = false;
                    h_bdry_normal((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1,0) = 0.0;
                    h_bdry_normal((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1,1) = 0.0;
                    h_bdry_normal((jsubmesh)*(2*nsubmesh_x1+1)+isubmesh-1,2) = 0.0;
                }
            }
        }
    }

    Kokkos::deep_copy(is_bdry, h_is_bdry);
    Kokkos::deep_copy(bdry_normal, h_bdry_normal);

    // for (int i=0; i<2*nsubmesh_x1*nsubmesh_x2+nsubmesh_x1+nsubmesh_x2; i++){
    //     if (h_is_bdry(i)){
    //         printf("Edge - %3d is boundary edge\n", i);
    //     }
    // }

    Kokkos::parallel_for("2D-meshobj-init", 1, KOKKOS_LAMBDA (const int) {
        pumi_mesh(0) = Mesh(nsubmesh_x1, Nel_tot_x1, nsubmesh_x2, Nel_tot_x2,
                        submesh_activity, elemoffset_start, elemoffset_skip,
                        nodeoffset_start, nodeoffset_skip_bot, nodeoffset_skip_mid, nodeoffset_skip_top,
                        is_bdry, bdry_normal, Nel_total_2D, Nnp_total_2D, host_isactive);
    });

    print_mesh_params(pumi_mesh, hc_submesh_x1, hc_submesh_x2);

    bool mesh_verified = verify_mesh_params(pumi_mesh, hc_submesh_x1, hc_submesh_x2);
    if (!mesh_verified){
        Kokkos::finalize();
        exit(0);
    }
    else{
        print_mesh_nodes(pumi_mesh, hc_submesh_x1, hc_submesh_x2, pumi_options);
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


} // namespace pumi

#endif
