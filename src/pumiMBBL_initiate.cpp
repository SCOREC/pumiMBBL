#include "pumiMBBLGPU.hpp"

namespace pumi {
///////// Mesh-Initiate Function declarations ///////////////////////////////////////

/**
 * @brief Reads the meshtype as string and converts it to unsigned int meshtype
 *
 * @param[in] mesh type in string format "uniform" or "minBL" or "maxBL"
 * \return unsigned int mesh type -- uniform=0x01 , minBL=0x02 and maxBL=0x04
 */
unsigned int get_submesh_type(std::string meshtype_string){
    std::vector<std::string> defined_types {"uniform","minBL","maxBL","arbitrary"};
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
    else if ( meshtype_string.compare(defined_types[3]) ==0 ){
        type = arbitrary;
    }
    else{
        std::cout << meshtype_string << " is an invalid meshtype\n";
        std::cout << "Valid types are \"uniform\", \"minBL\", \"maxBL\", \"arbitrary\"\n";
        std::cout << "Terminating\n";
        exit(0);
    }
    return type;
}

/**
 * @brief computes grading ratio using newton raphson
 *
 * @param[in] block Length
 * @param[in] smallest element size in the block
 * @param[in] number of elements in the block
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
 * @param[in] sum of number of submesh-blocks in each direction
 */
Mesh_Inputs* inputs_allocate(){
    Mesh_Inputs* pumi_inputs = new Mesh_Inputs;
    return pumi_inputs;
}

/**
 * @brief Frees up the alocated memory of mesh inputs struct members
 *
 * @param[in] mesh inputs struct pointer
 */
void inputs_deallocate(Mesh_Inputs* pumi_inputs){
    delete pumi_inputs;
}

std::vector<double> read_elemsize(double xmin, double xmax, std::string elemsize_file){
    std::vector<double> dx;
    double xlength = xmax-xmin;
    double dx_val;
    std::ifstream fn(elemsize_file);

    if (fn.is_open()){
        while (fn >> dx_val){
            dx.push_back(dx_val);
        }
    }
    else{
        std::cout << "ERROR: Element-size input file " << elemsize_file << " NOT FOUND -- TERMINATING\n";
        exit(0);
    }

    int nel = dx.size();
    double xlength_computed = 0.0;
    for (int i=0; i<nel; i++){
        xlength_computed += dx[i];
    }

    double diff = fabs(xlength_computed-xlength);
    if (diff/xlength > 1e-5){
        std::cout << "ERROR: Element-sizes from file " << elemsize_file << " does not match with input block-size\n";
        printf("Required block-size : %2.5e\n",xlength);
        printf("Block-size from file: %2.5e\n",xlength_computed);
        exit(0);
    }
    return dx;
}

/**
 * @brief Prints the all relevant 1D-mesh details
 *
 * @param[in] mesh object pointer
 * @param[in] host copy of x1-submesh object pointer
 */
void print_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1){

    printf("\n\nPUMI mesh parameter info [X1-Direction] :\n\n");
    printf("\tTotal elements along X1-direction = %d\n\n", pumi_mesh.Nel_tot_x1);
    for (int i=1; i<=pumi_mesh.nsubmesh_x1; i++){
        printf("\tSUBMESH %d  parameters:\n", i);
        printf("\n\t submesh-type   = ");
        if (h_submesh_x1[i]->meshtype & minBL){
            printf("leftBL\n");
            printf("\t left_t0        = %2.4e \t [m] Cell size of first/leftmost cell in left BL block\n", h_submesh_x1[i]->t0);
            printf("\t left_tN        = %2.4e \t [m] Cell size of last/rightmost cell in left BL block\n",
                (h_submesh_x1[i]->t0)*pow(h_submesh_x1[i]->r,h_submesh_x1[i]->Nel-1));
            printf("\t left_T         = %2.4e \t [m] Left boundary layer (left BL) thickness\n", h_submesh_x1[i]->length);
            printf("\t left_r         = %2.4e \t Grading ratio in left BL block\n", h_submesh_x1[i]->r);
            printf("\t left_Nel       = %d    \t Number of Cells in left BL block\n\n", h_submesh_x1[i]->Nel);
            // if (h_submesh_x1[i]->BL_coords.extent(0)){
                printf("\t %d leftBL node coords stored in a array\n\n", h_submesh_x1[i]->Nel+1);
            // }
        }
        if (h_submesh_x1[i]->meshtype & maxBL){
            printf("rightBL\n");
            printf("\t right_t0       = %2.4e \t [m] Cell size of last/rightmost cell in right BL block\n", h_submesh_x1[i]->t0);
            printf("\t right_tN       = %2.4e \t [m] Cell size of first/leftmost cell in right BL block\n",
            (h_submesh_x1[i]->t0)*pow(h_submesh_x1[i]->r,h_submesh_x1[i]->Nel-1));
            printf("\t right_T        = %2.4e \t [m] Left boundary layer (right BL) thickness\n", h_submesh_x1[i]->length);
            printf("\t right_r        = %2.4e \t Grading ratio in right BL mesh\n", h_submesh_x1[i]->r);
            printf("\t right_Nel      = %d    \t Number of Cells in left BL mesh region\n\n", h_submesh_x1[i]->Nel);
            // if (h_submesh_x1[i]->BL_coords.extent(0)){
                printf("\t %d rightBL node coords stored in a array\n\n", h_submesh_x1[i]->Nel+1);
            // }
        }
        if (h_submesh_x1[i]->meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_T      = %2.4e \t [m] Uniform block thickness\n", h_submesh_x1[i]->length);
            printf("\t uniform_dx1    = %2.4e \t [m] Cell size in uniform block\n", h_submesh_x1[i]->t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", h_submesh_x1[i]->Nel);
        }
        if (h_submesh_x1[i]->meshtype & arbitrary){
            printf("arbitrary\n");
            printf("\t arbitrary_T    = %2.4e \t [m] Arbitrary block thickness\n", h_submesh_x1[i]->length);
            printf("\t arbitrary_Nel  = %d    \t Number of Cells in arbitrary block\n\n", h_submesh_x1[i]->Nel);
        }
    }
}

/**
 * @brief Prints the all relevant 2D-mesh details
 *
 * @param[in] mesh object pointer
 * @param[in] host copy of x1-submesh object pointer
 * @param[in] host copy of x2-submesh object pointer
 */
void print_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2){

    printf("\n\nPUMI mesh parameter info [X1-Direction] :\n\n");
    printf("\tTotal elements along X1-direction = %d\n\n", pumi_mesh.Nel_tot_x1);
    for (int i=1; i<=pumi_mesh.nsubmesh_x1; i++){
        printf("\tSUBMESH %d  parameters:\n", i);
        printf("\n\t submesh-type   = ");
        if (h_submesh_x1[i]->meshtype & minBL){
            printf("leftBL\n");
            printf("\t left_t0        = %2.4e \t [m] Cell size of first/leftmost cell in left BL block\n", h_submesh_x1[i]->t0);
            printf("\t left_tN        = %2.4e \t [m] Cell size of last/rightmost cell in left BL block\n",
                (h_submesh_x1[i]->t0)*pow(h_submesh_x1[i]->r,h_submesh_x1[i]->Nel-1));
            printf("\t left_T         = %2.4e \t [m] Left boundary layer (left BL) thickness\n", h_submesh_x1[i]->length);
            printf("\t left_r         = %2.4e \t Grading ratio in left BL block\n", h_submesh_x1[i]->r);
            printf("\t left_Nel       = %d    \t Number of Cells in left BL block\n\n", h_submesh_x1[i]->Nel);
            if (h_submesh_x1[i]->BL_coords.extent(0)){
                printf("\t %d leftBL node coords stored in a array\n\n", h_submesh_x1[i]->Nel+1);
            }
        }
        if (h_submesh_x1[i]->meshtype & maxBL){
            printf("rightBL\n");
            printf("\t right_t0       = %2.4e \t [m] Cell size of last/rightmost cell in right BL block\n", h_submesh_x1[i]->t0);
            printf("\t right_tN       = %2.4e \t [m] Cell size of first/leftmost cell in right BL block\n",
            (h_submesh_x1[i]->t0)*pow(h_submesh_x1[i]->r,h_submesh_x1[i]->Nel-1));
            printf("\t right_T        = %2.4e \t [m] Left boundary layer (right BL) thickness\n", h_submesh_x1[i]->length);
            printf("\t right_r        = %2.4e \t Grading ratio in right BL mesh\n", h_submesh_x1[i]->r);
            printf("\t right_Nel      = %d    \t Number of Cells in left BL mesh region\n\n", h_submesh_x1[i]->Nel);
            if (h_submesh_x1[i]->BL_coords.extent(0)){
                printf("\t %d rightBL node coords stored in a array\n\n", h_submesh_x1[i]->Nel+1);
            }
        }
        if (h_submesh_x1[i]->meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_T      = %2.4e \t [m] Uniform block thickness\n", h_submesh_x1[i]->length);
            printf("\t uniform_dx1    = %2.4e \t [m] Cell size in uniform block\n", h_submesh_x1[i]->t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", h_submesh_x1[i]->Nel);
        }
        if (h_submesh_x1[i]->meshtype & arbitrary){
            printf("arbitrary\n");
            printf("\t arbitrary_T    = %2.4e \t [m] Arbitrary block thickness\n", h_submesh_x1[i]->length);
            printf("\t arbitrary_Nel  = %d    \t Number of Cells in arbitrary block\n\n", h_submesh_x1[i]->Nel);
        }
    }

    printf("PUMI mesh parameter info [X2-Direction] :\n\n");
    printf("\tTotal elements along X2-direction = %d\n\n", pumi_mesh.Nel_tot_x2);
    for (int i=1; i<=pumi_mesh.nsubmesh_x2; i++){
        printf("\tSUBMESH %d  parameters:\n", i);
        printf("\n\t submesh-type   = ");
        if (h_submesh_x2[i]->meshtype & minBL){
            printf("bottomBL\n");
            printf("\t bottom_t0      = %2.4e \t [m] Cell size of first/bottom-most cell in left BL block\n", h_submesh_x2[i]->t0);
            printf("\t bottom_tN      = %2.4e \t [m] Cell size of last /   top-most cell in left BL block\n",
                (h_submesh_x2[i]->t0)*pow(h_submesh_x2[i]->r,h_submesh_x2[i]->Nel-1));
            printf("\t bottom_T       = %2.4e \t [m] Left boundary layer (left BL) thickness\n", h_submesh_x2[i]->length);
            printf("\t bottom_r       = %2.4e \t Grading ratio in left BL block\n", h_submesh_x2[i]->r);
            printf("\t bottom_Nel     = %d    \t Number of Cells in left BL block\n\n", h_submesh_x2[i]->Nel);
            if (h_submesh_x2[i]->BL_coords.extent(0)){
                printf("\t %d bottomBL node coords stored in a array\n\n", h_submesh_x2[i]->Nel+1);
            }
        }
        if (h_submesh_x2[i]->meshtype & maxBL){
            printf("topBL\n");
            printf("\t top_t0         = %2.4e \t [m] Cell size of last /   top-most cell in right BL block\n", h_submesh_x2[i]->t0);
            printf("\t top_tN         = %2.4e \t [m] Cell size of first/bottom-most cell in right BL block\n",
            (h_submesh_x2[i]->t0)*pow(h_submesh_x2[i]->r,h_submesh_x2[i]->Nel-1));

            printf("\t top_T          = %2.4e \t [m] Left boundary layer (right BL) thickness\n", h_submesh_x2[i]->length);
            printf("\t top_r          = %2.4e \t Grading ratio in right BL mesh\n", h_submesh_x2[i]->r);
            printf("\t top_Nel        = %d    \t Number of Cells in left BL mesh region\n\n", h_submesh_x2[i]->Nel);
            if (h_submesh_x2[i]->BL_coords.extent(0)){
                printf("\t %d topBL node coords stored in a array\n\n", h_submesh_x2[i]->Nel+1);
            }
        }
        if (h_submesh_x2[i]->meshtype & uniform){
            printf("uniform\n");
            printf("\t uniform_T      = %2.4e \t [m] Uniform block thickness\n", h_submesh_x2[i]->length);
            printf("\t uniform_dx2    = %2.4e \t [m] Cell size in uniform block\n", h_submesh_x2[i]->t0);
            printf("\t uniform_Nel    = %d    \t Number of Cells in uniform block\n\n", h_submesh_x2[i]->Nel);
        }
        if (h_submesh_x2[i]->meshtype & arbitrary){
            printf("arbitrary\n");
            printf("\t arbitrary_T    = %2.4e \t [m] Arbitrary block thickness\n", h_submesh_x2[i]->length);
            printf("\t arbitrary_Nel  = %d    \t Number of Cells in arbitrary block\n\n", h_submesh_x2[i]->Nel);
        }
    }

    printf("PUMI submesh activity info :\n\n");
    for (int jsubmesh=pumi_mesh.nsubmesh_x2; jsubmesh>=1; jsubmesh--){
        if (jsubmesh != pumi_mesh.nsubmesh_x2){
            for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x1-1; isubmesh++ ){
                printf("_____________");
            }
            printf("__________________\n\n");
        }
        for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x1; isubmesh++ ){
            if (isubmesh-1){
                printf("|");
            }
            if(pumi_mesh.host_isactive[isubmesh][jsubmesh]){
                printf("    ACTIVE    ");
            }
            else{
                printf("    XXXXXX    ");
            }
        }
        printf("\n");
    }
    printf("\n\n");

    printf("Total active elements in 2D Mesh = %d\n",pumi_mesh.Nel_total);
    printf("Total active nodes in 2D Mesh    = %d\n",pumi_mesh.Nnp_total);
    printf("Total boundary faces in 2D Mesh  = %d\n",pumi_mesh.bdry.Nbdry_faces );

}

/**
 * @brief Prints node coordinates of 1D mesh to the terminal and also
 * writes the individual submesh coords to a file as well as
 * full mesh coordinates for all directions
 * @param[in] mesh object pointer
 * @param[in] CPU copy of x1-submesh object pointer
 */
void print_mesh_nodes(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, Mesh_Options pumi_options){

    if (pumi_options.print_node_option)
        printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
    FILE *mesh_coords_file;
    char mesh_coords_filename[30];
    sprintf(mesh_coords_filename,"X1_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    int inode=0;
    for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x1; isubmesh++){
        if (pumi_options.print_node_option)
            printf("X1-SUBMESH %d:\n", isubmesh );
        FILE *submesh_coords_file;
        char submesh_coords_filename[33];
        sprintf(submesh_coords_filename,"X1_submesh_%d_coords.dat",isubmesh);
        submesh_coords_file = fopen(submesh_coords_filename,"w");
        double icoord = h_submesh_x1[isubmesh]->xmin;
        double cell_size = h_submesh_x1[isubmesh]->t0;
        double r = h_submesh_x1[isubmesh]->r;
        if (h_submesh_x1[isubmesh]->meshtype & maxBL){
            cell_size = (h_submesh_x1[isubmesh]->t0)*
                        pow(h_submesh_x1[isubmesh]->r,h_submesh_x1[isubmesh]->Nel-1);
            r = 1.0/r;
        }
        if (h_submesh_x1[isubmesh]->meshtype & arbitrary){
            cell_size = h_submesh_x1[isubmesh]->host_BL_coords[1]-h_submesh_x1[isubmesh]->host_BL_coords[0];
        }

        if (isubmesh==1){
            fprintf(mesh_coords_file, "%.16e\n", icoord);
        }
        fprintf(submesh_coords_file, "%.16e\n", icoord);
        if (pumi_options.print_node_option)
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x1[isubmesh]->Nel; icell++){
            inode++;
            icoord += cell_size;
            if (h_submesh_x1[isubmesh]->meshtype &arbitrary){
                if (icell < h_submesh_x1[isubmesh]->Nel-1)
                    cell_size = h_submesh_x1[isubmesh]->host_BL_coords[icell+2]-h_submesh_x1[isubmesh]->host_BL_coords[icell+1];
            }
            else{
                cell_size *= r;
            }
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            if (pumi_options.print_node_option)
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
 * @param[in] mesh object pointer
 * @param[in] CPU copy of x1-submesh object pointer
 * @param[in] CPU copy of x2-submesh object pointer
 */
void print_mesh_nodes(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2, Mesh_Options pumi_options){

    if (pumi_options.print_node_option)
        printf("\nPrinting the coordinates of the nodes in the pumi mesh...\n\n");
    FILE *mesh_coords_file;
    char mesh_coords_filename[30];
    sprintf(mesh_coords_filename,"X1_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    int inode=0;
    for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x1; isubmesh++){
        if (pumi_options.print_node_option)
            printf("X1-SUBMESH %d:\n", isubmesh );
        FILE *submesh_coords_file;
        char submesh_coords_filename[33];
        sprintf(submesh_coords_filename,"X1_submesh_%d_coords.dat",isubmesh);
        submesh_coords_file = fopen(submesh_coords_filename,"w");
        double icoord = h_submesh_x1[isubmesh]->xmin;
        double cell_size = h_submesh_x1[isubmesh]->t0;
        double r = h_submesh_x1[isubmesh]->r;
        if (h_submesh_x1[isubmesh]->meshtype & maxBL){
            cell_size = (h_submesh_x1[isubmesh]->t0)*
                        pow(h_submesh_x1[isubmesh]->r,h_submesh_x1[isubmesh]->Nel-1);
            r = 1.0/r;
        }
        if (h_submesh_x1[isubmesh]->meshtype & arbitrary){
            cell_size = h_submesh_x1[isubmesh]->host_BL_coords[1]-h_submesh_x1[isubmesh]->host_BL_coords[0];
        }

        if (isubmesh==1){
            fprintf(mesh_coords_file, "%.16e\n", icoord);
        }
        fprintf(submesh_coords_file, "%.16e\n", icoord);
        if (pumi_options.print_node_option)
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x1[isubmesh]->Nel; icell++){
            inode++;
            icoord += cell_size;
            if (h_submesh_x1[isubmesh]->meshtype &arbitrary){
                if (icell < h_submesh_x1[isubmesh]->Nel-1)
                    cell_size = h_submesh_x1[isubmesh]->host_BL_coords[icell+2]-h_submesh_x1[isubmesh]->host_BL_coords[icell+1];
            }
            else{
                cell_size *= r;
            }
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            if (pumi_options.print_node_option)
                printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        }
        printf("\n\t X1-SUBMESH %d -- Coordinates written to file %s\n\n", isubmesh, submesh_coords_filename);
        fclose(submesh_coords_file);
    }
    fclose(mesh_coords_file);

    sprintf(mesh_coords_filename,"X2_fullmesh_coords.dat");
    mesh_coords_file = fopen(mesh_coords_filename,"w");
    inode = 0;
    for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x2; isubmesh++){
        if (pumi_options.print_node_option)
            printf("X2-SUBMESH %d:\n", isubmesh );
        FILE *submesh_coords_file;
        char submesh_coords_filename[33];
        sprintf(submesh_coords_filename,"X2_submesh_%d_coords.dat",isubmesh);
        submesh_coords_file = fopen(submesh_coords_filename,"w");
        double icoord = h_submesh_x2[isubmesh]->xmin;
        double cell_size = h_submesh_x2[isubmesh]->t0;
        double r = h_submesh_x2[isubmesh]->r;
        if (h_submesh_x2[isubmesh]->meshtype & maxBL){
            cell_size = (h_submesh_x2[isubmesh]->t0)*
                        pow(h_submesh_x2[isubmesh]->r,h_submesh_x2[isubmesh]->Nel-1);
            r = 1.0/r;
        }
        if (h_submesh_x2[isubmesh]->meshtype & arbitrary){
            cell_size = h_submesh_x2[isubmesh]->host_BL_coords[1]-h_submesh_x2[isubmesh]->host_BL_coords[0];
        }

        if (isubmesh==1){
            fprintf(mesh_coords_file, "%.16e\n", icoord);
        }
        fprintf(submesh_coords_file, "%.16e\n", icoord);
        if (pumi_options.print_node_option)
            printf("\t\tNode %6d: %2.8e\n", inode, icoord );
        for (int icell=0; icell<h_submesh_x2[isubmesh]->Nel; icell++){
            inode++;
            icoord += cell_size;
            if (h_submesh_x2[isubmesh]->meshtype &arbitrary){
                if (icell < h_submesh_x2[isubmesh]->Nel-1)
                    cell_size = h_submesh_x2[isubmesh]->host_BL_coords[icell+2]-h_submesh_x2[isubmesh]->host_BL_coords[icell+1];
            }
            else{
                cell_size *= r;
            }
            fprintf(mesh_coords_file, "%.16e\n", icoord);
            fprintf(submesh_coords_file, "%.16e\n", icoord);
            if (pumi_options.print_node_option)
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
 * @param[in] mesh object pointer
 * @param[in] host copy of x1-submesh object pointer
 */
bool verify_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1){

    printf("\n\nNow verifying valdity of pumi mesh parameters for\n");
    int flag = 0;
    for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x1; isubmesh++){
        printf("\tX1-SUBMESH %d:\n", isubmesh );
        if (h_submesh_x1[isubmesh]->meshtype & minBL){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t left_Nel    -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh]->r >= 1.0)){
                printf("\t\t left_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh]->r);
                flag++;
            }
            else{
                printf("\t\t left_r      -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh]->t0*h_submesh_x1[isubmesh]->Nel;
            if (!(h_submesh_x1[isubmesh]->length > 0.0)){
                printf("\t\t left_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh]->length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh]->length >= min_computed_T)){
                printf("\t\t left_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh]->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t left_T      -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh]->meshtype & maxBL){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t right_Nel   -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh]->r >= 1.0)){
                printf("\t\t right_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh]->r);
                flag++;
            }
            else{
                printf("\t\t right_r     -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh]->t0*h_submesh_x1[isubmesh]->Nel;
            if (!(h_submesh_x1[isubmesh]->length > 0.0)){
                printf("\t\t right_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh]->length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh]->length >= min_computed_T)){
                printf("\t\t right_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh]->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t right_T     -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh]->meshtype & uniform){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh]->meshtype & arbitrary){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t arbitrary_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t arbitrary_Nel -- verified...\n");
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
 * @param[in] mesh object pointer
 * @param[in] host copy of x1-submesh object pointer
 * @param[in] host copy of x2-submesh object pointer
 */
bool verify_mesh_params(Mesh pumi_mesh, SubmeshHostViewPtr h_submesh_x1, SubmeshHostViewPtr h_submesh_x2){

    printf("\n\nNow verifying valdity of pumi mesh parameters for\n");
    int flag = 0;
    for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x1; isubmesh++){
        printf("\tX1-SUBMESH %d:\n", isubmesh );
        if (h_submesh_x1[isubmesh]->meshtype & minBL){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t left_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t left_Nel    -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh]->r >= 1.0)){
                printf("\t\t left_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh]->r);
                flag++;
            }
            else{
                printf("\t\t left_r      -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh]->t0*h_submesh_x1[isubmesh]->Nel;
            if (!(h_submesh_x1[isubmesh]->length > 0.0)){
                printf("\t\t left_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh]->length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh]->length >= min_computed_T)){
                printf("\t\t left_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh]->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t left_T      -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh]->meshtype & maxBL){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t right_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t right_Nel   -- verified...\n");
            }
            if (!(h_submesh_x1[isubmesh]->r >= 1.0)){
                printf("\t\t right_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x1[isubmesh]->r);
                flag++;
            }
            else{
                printf("\t\t right_r     -- verified...\n");
            }
            double min_computed_T = h_submesh_x1[isubmesh]->t0*h_submesh_x1[isubmesh]->Nel;
            if (!(h_submesh_x1[isubmesh]->length > 0.0)){
                printf("\t\t right_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x1[isubmesh]->length);
                flag++;
            }
            else if (!(h_submesh_x1[isubmesh]->length >= min_computed_T)){
                printf("\t\t right_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x1[isubmesh]->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t right_T     -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh]->meshtype & uniform){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
        if (h_submesh_x1[isubmesh]->meshtype & arbitrary){
            if (!(h_submesh_x1[isubmesh]->Nel > 0)){
                printf("\t\t arbitrary_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x1[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t arbitrary_Nel -- verified...\n");
            }
        }
    }


    for (int isubmesh=1; isubmesh<=pumi_mesh.nsubmesh_x2; isubmesh++){
        printf("\tX2-SUBMESH %d:\n", isubmesh );
        if (h_submesh_x2[isubmesh]->meshtype & minBL){
            if (!(h_submesh_x2[isubmesh]->Nel > 0)){
                printf("\t\t bottom_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x2[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t bottom_Nel  -- verified...\n");
            }
            if (!(h_submesh_x2[isubmesh]->r >= 1.0)){
                printf("\t\t bottom_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x2[isubmesh]->r);
                flag++;
            }
            else{
                printf("\t\t bottom_r    -- verified...\n");
            }
            double min_computed_T = h_submesh_x2[isubmesh]->t0*h_submesh_x2[isubmesh]->Nel;
            if (!(h_submesh_x2[isubmesh]->length > 0.0)){
                printf("\t\t bottom_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x2[isubmesh]->length);
                flag++;
            }
            else if (!(h_submesh_x2[isubmesh]->length >= min_computed_T)){
                printf("\t\t bottom_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x2[isubmesh]->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t bottom_T    -- verified...\n");
            }
        }
        if (h_submesh_x2[isubmesh]->meshtype & maxBL){
            if (!(h_submesh_x2[isubmesh]->Nel > 0)){
                printf("\t\t top_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x2[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t top_Nel     -- verified...\n");
            }
            if (!(h_submesh_x2[isubmesh]->r >= 1.0)){
                printf("\t\t top_r   = %2.4f is not a valid input. It has to be a greater than 1.0\n", h_submesh_x2[isubmesh]->r);
                flag++;
            }
            else{
                printf("\t\t top_r       -- verified...\n");
            }
            double min_computed_T = h_submesh_x2[isubmesh]->t0*h_submesh_x2[isubmesh]->Nel;
            if (!(h_submesh_x2[isubmesh]->length > 0.0)){
                printf("\t\t top_T   = %2.4f is not a valid input. It has to be postive\n", h_submesh_x2[isubmesh]->length);
                flag++;
            }
            else if (!(h_submesh_x2[isubmesh]->length >= min_computed_T)){
                printf("\t\t top_T   = %2.4f is not a valid input. It is smaller than the minimum computed BL length %2.4f\n",
                h_submesh_x2[isubmesh]->length, min_computed_T);
                flag++;
            }
            else{
                printf("\t\t top_T       -- verified...\n");
            }
        }
        if (h_submesh_x2[isubmesh]->meshtype & uniform){
            if (!(h_submesh_x2[isubmesh]->Nel > 0)){
                printf("\t\t uniform_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x2[isubmesh]->Nel);
            }
            else{
                printf("\t\t uniform_Nel -- verified...\n");
            }
        }
        if (h_submesh_x2[isubmesh]->meshtype & arbitrary){
            if (!(h_submesh_x2[isubmesh]->Nel > 0)){
                printf("\t\t arbitrary_Nel = %d is not a valid input. It has to be a positive integer.\n", h_submesh_x2[isubmesh]->Nel);
                flag++;
            }
            else{
                printf("\t\t arbitrary_Nel -- verified...\n");
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
 * @param[in] pointer object to pumi inputs structure which containts
 * @param[in] object to structure containing user options for the mesh
 * @param[in] direction along which submesh object is being initialized
 * \return object to SubmeshInit class
 */
SubmeshInit submesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, int dir){
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

    SubmeshHostViewPtr submesh_host_copy = new Submesh*[nsubmesh+2];

    if (pumi_options.BL_storage_option == store_BL_coords_OFF){
        pumi_options.BL_storage_option = store_BL_coords_ON;
    }

    double xmin, xmax, xlength, t0, tN, r, r_t0_ratio, logr;
    std::string elemsize_file;
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
        DoubleView BLcoords;
        double *host_BLcoords = nullptr;
        // padding submesh to min-side
        xmax = pumi_inputs->domain_x1_min;
        xmin = xmax - xlength;
        Unassigned_Submesh tmp_obj_min(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        Unassigned_Submesh *tmp_obj_min_host = new Unassigned_Submesh(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        submesh_host_copy[0] = tmp_obj_min_host;
        auto dvc_ptr_1 = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_min);
        h_submesh(0) = DevicePointer<Submesh> (dvc_ptr_1);
        // padding submesh to max-side
        xmin = pumi_inputs->domain_x1_min+total_length;
        xmax = xmin + xlength;
        Unassigned_Submesh tmp_obj_max(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        Unassigned_Submesh *tmp_obj_max_host = new Unassigned_Submesh(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        submesh_host_copy[nsubmesh+1] = tmp_obj_max_host;
        auto dvc_ptr_2 = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_max);
        h_submesh(nsubmesh+1) = DevicePointer<Submesh>(dvc_ptr_2);
    }
    else if (dir == x2_dir){
        double total_length = 0.0;
        for (int isubmesh=0; isubmesh<nsubmesh; isubmesh++){
            total_length += pumi_inputs->block_length_x2[isubmesh];
        }
        type = unassigned;
        xlength = 1000.0*total_length;
        DoubleView BLcoords;
        double *host_BLcoords = nullptr;
        // padding submesh to min-side
        xmax = pumi_inputs->domain_x2_min;
        xmin = xmax - xlength;
        Unassigned_Submesh tmp_obj_min(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        Unassigned_Submesh *tmp_obj_min_host = new Unassigned_Submesh(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        submesh_host_copy[0] = tmp_obj_min_host;
        auto dvc_ptr_1 = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_min);
        h_submesh(0) = DevicePointer<Submesh>(dvc_ptr_1);
        // padding submesh to max-side
        xmin = pumi_inputs->domain_x2_min+total_length;
        xmax = xmin + xlength;
        Unassigned_Submesh tmp_obj_max(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        Unassigned_Submesh *tmp_obj_max_host = new Unassigned_Submesh(xmin,xmax,0,0.0,0.0,xlength,0,0.0,0.0,BLcoords,host_BLcoords);
        submesh_host_copy[nsubmesh+1] = tmp_obj_max_host;
        auto dvc_ptr_2 = copyForDevice<Submesh, Unassigned_Submesh> (tmp_obj_max);
        h_submesh(nsubmesh+1) = DevicePointer<Submesh>(dvc_ptr_2);
    }


    for (int isubmesh=1; isubmesh<=nsubmesh; isubmesh++){
        DoubleView BLcoords;
        double *host_BLcoords;
        std::string BLcoordsname;

        if (dir == x1_dir){
            type = get_submesh_type(pumi_inputs->meshtype_x1[isubmesh-1]);

            xlength = pumi_inputs->block_length_x1[isubmesh-1];
            if (isubmesh==1){
                xmin = pumi_inputs->domain_x1_min;
            }
            else{
                xmin = xmax;
            }
            xmax = xmin + xlength;
            t0 = pumi_inputs->min_elem_size_x1[isubmesh-1];
            tN = pumi_inputs->max_elem_size_x1[isubmesh-1];

            BLcoordsname = "x1-BLcoords-submesh-";
            BLcoordsname += std::to_string(isubmesh);

            if (type & arbitrary)
                elemsize_file = pumi_inputs->arbitrary_x1_elemsize_file[isubmesh-1];
        }
        else if (dir == x2_dir){
            type = get_submesh_type(pumi_inputs->meshtype_x2[isubmesh-1]);

            xlength = pumi_inputs->block_length_x2[isubmesh-1];
            if (isubmesh==1){
                xmin = pumi_inputs->domain_x2_min;
            }
            else{
                xmin = xmax;
            }
            xmax = xmin + xlength;
            t0 = pumi_inputs->min_elem_size_x2[isubmesh-1];
            tN = pumi_inputs->max_elem_size_x2[isubmesh-1];

            BLcoordsname = "x2-BLcoords-submesh-";
            BLcoordsname += std::to_string(isubmesh);

            if (type & arbitrary)
                elemsize_file = pumi_inputs->arbitrary_x2_elemsize_file[isubmesh-1];
        }

        if (isubmesh > 1){
            Nel_cumulative += Nel;
        }

        if (type & uniform){
            r = 1.0;
            double Nel_dbl = xlength/t0;
            int Nel_int = (int) Nel_dbl;
            double diff = Nel_dbl-Nel_int;
            if (1.0-diff < 1e-3)
                Nel_int++;
            Nel = Nel_int;
            t0 = xlength/Nel;
            r_t0_ratio = (r-1.0)/t0;
            logr = log(r);
            host_BLcoords = nullptr;

            Uniform_Submesh tmp_obj(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            Uniform_Submesh *tmp_obj_host = new Uniform_Submesh(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            submesh_host_copy[isubmesh] = tmp_obj_host;
            auto dvc_ptr = copyForDevice<Submesh, Uniform_Submesh> (tmp_obj);
            h_submesh(isubmesh) = DevicePointer<Submesh>(dvc_ptr);
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
                BLcoords = DoubleView (BLcoordsname, Nel+1);
                host_BLcoords = new double[Nel+1];
                DoubleView::HostMirror h_BLcoords = Kokkos::create_mirror_view(BLcoords);
                h_BLcoords(0) = xmin;
                host_BLcoords[0] = xmin;
                double cell_size = t0;
                for (int icell=1; icell<Nel; icell++){
                    h_BLcoords(icell) = h_BLcoords(icell-1) + cell_size;
                    host_BLcoords[icell] = host_BLcoords[icell-1] + cell_size;
                    cell_size *= r;
                }
                h_BLcoords(Nel) = xmax;
                host_BLcoords[Nel] = xmax;

                Kokkos::deep_copy(BLcoords, h_BLcoords);
            }
            MinBL_Submesh tmp_obj(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            MinBL_Submesh *tmp_obj_host = new MinBL_Submesh(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            submesh_host_copy[isubmesh] = tmp_obj_host;
            auto dvc_ptr = copyForDevice<Submesh, MinBL_Submesh> (tmp_obj);
            h_submesh(isubmesh) = DevicePointer<Submesh>(dvc_ptr);
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
                BLcoords = DoubleView (BLcoordsname, Nel+1);
                host_BLcoords = new double[Nel+1];
                DoubleView::HostMirror h_BLcoords = Kokkos::create_mirror_view(BLcoords);
                h_BLcoords(0) = xmin;
                host_BLcoords[0] = xmin;
                double cell_size = t0*pow(r,Nel-1);
                for (int icell=1; icell<Nel; icell++){
                    h_BLcoords(icell) = h_BLcoords(icell-1) + cell_size;
                    host_BLcoords[icell] = host_BLcoords[icell-1] + cell_size;
                    cell_size /= r;
                }
                h_BLcoords(Nel) = xmax;
                host_BLcoords[Nel] = xmax;

                Kokkos::deep_copy(BLcoords, h_BLcoords);
            }
            MaxBL_Submesh tmp_obj(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            MaxBL_Submesh *tmp_obj_host = new MaxBL_Submesh(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            submesh_host_copy[isubmesh] = tmp_obj_host;
            auto dvc_ptr = copyForDevice<Submesh, MaxBL_Submesh> (tmp_obj);
            h_submesh(isubmesh) = DevicePointer<Submesh>(dvc_ptr);
        }
        if (type & arbitrary){

            std::vector<double> BL_elemsize = read_elemsize(xmin, xmax, elemsize_file);
            Nel = BL_elemsize.size();
            r = -1.0;
            r_t0_ratio = -1.0;
            logr = -1.0;

            BLcoords = DoubleView (BLcoordsname, Nel+1);
            host_BLcoords = new double[Nel+1];
            DoubleView::HostMirror h_BLcoords = Kokkos::create_mirror_view(BLcoords);
            h_BLcoords(0) = xmin;
            host_BLcoords[0] = xmin;
            for (int icell=1; icell<Nel; icell++){
                h_BLcoords(icell) = h_BLcoords(icell-1) + BL_elemsize[icell-1];
                host_BLcoords[icell] = host_BLcoords[icell-1] + BL_elemsize[icell-1];
            }
            h_BLcoords(Nel) = xmax;
            host_BLcoords[Nel] = xmax;

            Kokkos::deep_copy(BLcoords, h_BLcoords);

            Arbitrary_Submesh tmp_obj(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            Arbitrary_Submesh *tmp_obj_host = new Arbitrary_Submesh(xmin,xmax,Nel,t0,r,xlength,Nel_cumulative,r_t0_ratio,logr,BLcoords,host_BLcoords);
            submesh_host_copy[isubmesh] = tmp_obj_host;
            auto dvc_ptr = copyForDevice<Submesh, Arbitrary_Submesh> (tmp_obj);
            h_submesh(isubmesh) = DevicePointer<Submesh>(dvc_ptr);
        }

    }

    Kokkos::deep_copy(submesh, h_submesh);
    // *hc_submesh = submesh_host_copy;

    SubmeshInit init(submesh, submesh_host_copy);

    return init;
}

/**
 * @brief Initiates the mesh and returns the final 1D-mesh object
 *
 * @param[in] pointer object to pumi inputs structure which containts
 * @param[in] object to structure containing user options for the mesh
 * @param[in] x1-submesh object pointer
 * @param[in] Copy of x1-submesh object pointer on CPU
 * \return Final mesh object
 */
Mesh mesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, SubmeshDeviceViewPtr , SubmeshHostViewPtr hc_submesh_x1){

    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int Nel_tot_x1;

    Nel_tot_x1 = hc_submesh_x1[nsubmesh_x1]->Nel + hc_submesh_x1[nsubmesh_x1]->Nel_cumulative;
    int num_bdry_faces = 2;
    MeshBdry bdry = MeshBdry(num_bdry_faces);

    BlockInterface blk_if = BlockInterface(hc_submesh_x1, nsubmesh_x1);

    // Kokkos::parallel_for("1D-meshobj-init", 1, KOKKOS_LAMBDA (const int) {
    //     pumi_mesh(0) = Mesh(nsubmesh_x1, Nel_tot_x1, bdry);
    // });
    Mesh pumi_mesh = Mesh(nsubmesh_x1, Nel_tot_x1, bdry, blk_if);
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
 * @param[in] pointer object to pumi inputs structure which containts
 * @param[in] object to structure containing user options for the mesh
 * @param[in] x1-submesh object pointer
 * @param[in] Copy of x1-submesh object pointer on CPU
 * @param[in] x2-submesh object pointer
 * @param[in] Copy of x2-submesh object pointer on CPU
 * \return Final mesh object
 */
Mesh mesh_initialize(Mesh_Inputs *pumi_inputs, Mesh_Options pumi_options, SubmeshDeviceViewPtr , SubmeshHostViewPtr hc_submesh_x1,
                            SubmeshDeviceViewPtr , SubmeshHostViewPtr hc_submesh_x2){

    int nsubmesh_x1 = pumi_inputs->nsubmesh_x1;
    int nsubmesh_x2 = pumi_inputs->nsubmesh_x2;
    int Nel_tot_x1, Nel_tot_x2;

    Nel_tot_x1 = hc_submesh_x1[nsubmesh_x1]->Nel + hc_submesh_x1[nsubmesh_x1]->Nel_cumulative;
    Nel_tot_x2 = hc_submesh_x2[nsubmesh_x2]->Nel + hc_submesh_x2[nsubmesh_x2]->Nel_cumulative;

    Kokkos::View<bool**> isactive("isactive",nsubmesh_x1+2,nsubmesh_x2+2);
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

    Kokkos::View<bool**>::HostMirror h_isactive = Kokkos::create_mirror_view(isactive);
    for (int isubmesh=0; isubmesh<nsubmesh_x1+2; isubmesh++ ){
        for (int jsubmesh=0; jsubmesh<nsubmesh_x2+2; jsubmesh++){
            if (isubmesh==0 || isubmesh==nsubmesh_x1+1 || jsubmesh==0 || jsubmesh==nsubmesh_x2+1){
                h_isactive(isubmesh,jsubmesh) = false;
            }
            else{
                h_isactive(isubmesh,jsubmesh) = pumi_inputs->isactive[isubmesh-1][jsubmesh-1];
            }
        }
    }
    Kokkos::deep_copy(isactive, h_isactive);

    MeshOffsets offsets = MeshOffsets(hc_submesh_x1, nsubmesh_x1, hc_submesh_x2, nsubmesh_x2, host_isactive);
    MeshBdry bdry = MeshBdry(hc_submesh_x1, nsubmesh_x1, hc_submesh_x2, nsubmesh_x2, host_isactive);

    BlockInterface blk_if = BlockInterface(hc_submesh_x1, nsubmesh_x1, hc_submesh_x2, nsubmesh_x2, host_isactive);
    MeshBST bst;
    bst.initialize_MeshBST(blk_if,hc_submesh_x1,nsubmesh_x1,hc_submesh_x2,nsubmesh_x2,host_isactive);
    // Kokkos::parallel_for("2D-meshobj-init", 1, KOKKOS_LAMBDA (const int) {
    //     pumi_mesh(0) = Mesh(nsubmesh_x1, Nel_tot_x1, nsubmesh_x2, Nel_tot_x2,
    //                         isactive, host_isactive, offsets, bdry, offsets.Nel_total, offsets.Nnp_total);
    // });
    Mesh pumi_mesh = Mesh(nsubmesh_x1, Nel_tot_x1, nsubmesh_x2, Nel_tot_x2,
                        isactive, host_isactive, offsets, bdry, blk_if, bst, offsets.Nel_total, offsets.Nnp_total);

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

/**
* @brief Class constructor.
*
* @param[in] Array of x1-submesh object pointers on host
* @param[in] Number of x1-blocks in mesh
* @param[in] Array of x2-submesh object pointers on host
* @param[in] Number of x1-blocks in mesh
*/
MeshOffsets::MeshOffsets(SubmeshHostViewPtr hc_submesh_x1,
                         int Nx,
                         SubmeshHostViewPtr hc_submesh_x2,
                         int Ny,
                         bool** host_isactive){
    is_fullmesh = true;
    for (int isubmesh=1; isubmesh<=Nx; isubmesh++ ){
        for (int jsubmesh=1; jsubmesh<=Ny; jsubmesh++){
            if (!host_isactive[isubmesh][jsubmesh]){
                is_fullmesh = false;
            }
        }
    }

    int Nel_tot_x1 = hc_submesh_x1[Nx]->Nel + hc_submesh_x1[Nx]->Nel_cumulative;
    int Nel_tot_x2 = hc_submesh_x2[Ny]->Nel + hc_submesh_x2[Ny]->Nel_cumulative;
    Nel_total = Nel_tot_x1*Nel_tot_x2;
    Nnp_total = (Nel_tot_x1+1)*(Nel_tot_x2+1);

    elemoffset_start = Kokkos::View<int**>("elemoffset_start", Nx+2, Ny+2);
    elemoffset_skip = Kokkos::View<int*>("elemoffset_skip", Ny+2);
    nodeoffset_start = Kokkos::View<int**>("nodeoffset_start", Nx+2, Ny+2);
    nodeoffset_skip_bot = Kokkos::View<int**>("nodeoffset_skip_bot", Nx+2, Ny+2);
    nodeoffset_skip_mid = Kokkos::View<int**>("nodeoffset_skip_mid", Nx+2, Ny+2);
    nodeoffset_skip_top = Kokkos::View<int**>("nodeoffset_skip_top", Nx+2, Ny+2);

    int **nodeoffset = new int*[Nx+2];
    host_elemoffset_start = new int*[Nx+2];
    host_elemoffset_skip = new int[Ny+2];
    host_nodeoffset_start = new int*[Nx+2];
    host_nodeoffset_skip_bot = new int*[Nx+2];
    host_nodeoffset_skip_mid = new int*[Nx+2];
    host_nodeoffset_skip_top = new int*[Nx+2];

    for (int i=0; i<Nx+2; i++){
        nodeoffset[i] = new int[Nel_tot_x2+1];
        host_elemoffset_start[i] = new int[Ny+2];
        host_nodeoffset_start[i] = new int[Ny+2];
        host_nodeoffset_skip_bot[i] = new int[Ny+2];
        host_nodeoffset_skip_mid[i] = new int[Ny+2];
        host_nodeoffset_skip_top[i] = new int[Ny+2];
    }

    if (!is_fullmesh){
        Kokkos::View<int**>::HostMirror h_elemoffset_start = Kokkos::create_mirror_view(elemoffset_start);
        Kokkos::View<int*>::HostMirror h_elemoffset_skip = Kokkos::create_mirror_view(elemoffset_skip);
        Kokkos::View<int**>::HostMirror h_nodeoffset_start = Kokkos::create_mirror_view(nodeoffset_start);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_bot = Kokkos::create_mirror_view(nodeoffset_skip_bot);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_mid = Kokkos::create_mirror_view(nodeoffset_skip_mid);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_top = Kokkos::create_mirror_view(nodeoffset_skip_top);

        // elemoffsets
        int elemstart_init, elemskip;
        int elemstart = 0;

        for (int jsubmesh=1; jsubmesh<=Ny; jsubmesh++){
            elemstart_init = elemstart;

            for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
                if(host_isactive[isubmesh][jsubmesh]){
                    h_elemoffset_start(isubmesh,jsubmesh) = elemstart;
                }
                else{
                    elemstart +=  hc_submesh_x1[isubmesh]->Nel;
                    h_elemoffset_start(isubmesh,jsubmesh) = -1;
                }
            }
            elemskip = elemstart-elemstart_init;
            h_elemoffset_skip(jsubmesh) = elemskip;
            elemstart = elemstart_init + elemskip*hc_submesh_x2[jsubmesh]->Nel;
        }
        Nel_total -= elemstart;

        int jnp;
        int jsubmesh = 1;
        int nodestart = 0;
        for (jnp=0; jnp<hc_submesh_x2[jsubmesh]->Nel; jnp++){
            int Jnp = jnp + hc_submesh_x2[jsubmesh]->Nel_cumulative;
            for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
                if(host_isactive[isubmesh][jsubmesh]){
                    nodeoffset[isubmesh][Jnp] = nodestart;
                }
                else{
                    nodeoffset[isubmesh][Jnp] = -1;
                    if (isubmesh==1){
                        if (host_isactive[isubmesh+1][jsubmesh]){
                            nodestart += hc_submesh_x1[isubmesh]->Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh]->Nel+1);
                        }
                    }
                    else if (isubmesh == Nx){
                        nodestart += hc_submesh_x1[isubmesh]->Nel;
                    }
                    else{
                        if (host_isactive[isubmesh+1][jsubmesh]){
                            nodestart += (hc_submesh_x1[isubmesh]->Nel-1);
                        }
                        else{
                            nodestart += hc_submesh_x1[isubmesh]->Nel;
                        }
                    }
                }
            }
        }

        jnp = hc_submesh_x2[jsubmesh]->Nel;
        for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
            if (jsubmesh == Ny){
                int Jnp = jnp + hc_submesh_x2[jsubmesh]->Nel_cumulative;
                if (host_isactive[isubmesh][jsubmesh]){
                    nodeoffset[isubmesh][Jnp] = nodestart;
                }
                else{
                    nodeoffset[isubmesh][Jnp] = -1;
                    if (isubmesh==1){
                        if (host_isactive[isubmesh+1][jsubmesh]){
                            nodestart += hc_submesh_x1[isubmesh]->Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh]->Nel+1);
                        }
                    }
                    else if (isubmesh==Nx){
                        nodestart += hc_submesh_x1[isubmesh]->Nel;
                    }
                    else{
                        if (host_isactive[isubmesh+1][jsubmesh]){
                            nodestart += (hc_submesh_x1[isubmesh]->Nel-1);
                        }
                        else{
                            nodestart += hc_submesh_x1[isubmesh]->Nel;
                        }
                    }
                }
            }
        }

        for (jsubmesh=2; jsubmesh<=Ny; jsubmesh++){
            jnp = 0;
            int Jnp = jnp + hc_submesh_x2[jsubmesh]->Nel_cumulative;

            for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
                if (host_isactive[isubmesh][jsubmesh] || host_isactive[isubmesh][jsubmesh-1]){
                    nodeoffset[isubmesh][Jnp] = nodestart;
                }
                else{
                    nodeoffset[isubmesh][Jnp] = -1;
                    if (isubmesh==1){
                        if (host_isactive[isubmesh+1][jsubmesh] || host_isactive[isubmesh+1][jsubmesh-1]){
                            nodestart += hc_submesh_x1[isubmesh]->Nel;
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh]->Nel+1);
                        }
                    }
                    else if (isubmesh==Nx){
                        nodestart += (hc_submesh_x1[isubmesh]->Nel);
                    }
                    else{
                        if (host_isactive[isubmesh+1][jsubmesh] || host_isactive[isubmesh+1][jsubmesh-1]){
                            nodestart += (hc_submesh_x1[isubmesh]->Nel-1);
                        }
                        else{
                            nodestart += (hc_submesh_x1[isubmesh]->Nel);
                        }
                    }
                }
            }

            for (jnp=1; jnp<hc_submesh_x2[jsubmesh]->Nel; jnp++){
                int Jnp = jnp + hc_submesh_x2[jsubmesh]->Nel_cumulative;

                for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
                    if (host_isactive[isubmesh][jsubmesh]){
                        nodeoffset[isubmesh][Jnp] = nodestart;
                    }
                    else{
                        nodeoffset[isubmesh][Jnp] = -1;
                        if (isubmesh==1){
                            if (host_isactive[isubmesh+1][jsubmesh]){
                                nodestart += hc_submesh_x1[isubmesh]->Nel;
                            }
                            else{
                                nodestart += (hc_submesh_x1[isubmesh]->Nel+1);
                            }
                        }
                        else if (isubmesh == Nx){
                            nodestart += hc_submesh_x1[isubmesh]->Nel;
                        }
                        else{
                            if (host_isactive[isubmesh+1][jsubmesh]){
                                nodestart += (hc_submesh_x1[isubmesh]->Nel-1);
                            }
                            else{
                                nodestart += hc_submesh_x1[isubmesh]->Nel;
                            }
                        }
                    }
                }
            }

            if (jsubmesh==Ny){
                jnp = hc_submesh_x2[jsubmesh]->Nel;
                int Jnp = jnp + hc_submesh_x2[jsubmesh]->Nel_cumulative;

                for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
                    if (host_isactive[isubmesh][jsubmesh]){
                        nodeoffset[isubmesh][Jnp] = nodestart;
                    }
                    else{
                        nodeoffset[isubmesh][Jnp] = -1;
                        if (isubmesh==1){
                            if (host_isactive[isubmesh+1][jsubmesh]){
                                nodestart += hc_submesh_x1[isubmesh]->Nel;
                            }
                            else{
                                nodestart += (hc_submesh_x1[isubmesh]->Nel+1);
                            }
                        }
                        else if (isubmesh==Nx){
                            nodestart += hc_submesh_x1[isubmesh]->Nel;
                        }
                        else{
                            if (host_isactive[isubmesh+1][jsubmesh]){
                                nodestart += (hc_submesh_x1[isubmesh]->Nel-1);
                            }
                            else{
                                nodestart += hc_submesh_x1[isubmesh]->Nel;
                            }
                        }
                    }
                }
            }

        }

        Nnp_total -= nodestart;

        for (int jsubmesh=1; jsubmesh<=Ny; jsubmesh++){
            int Jnp = hc_submesh_x2[jsubmesh]->Nel_cumulative;
            for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
                if (host_isactive[isubmesh][jsubmesh]){
                    h_nodeoffset_start(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp];
                    h_nodeoffset_skip_bot(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp+1]-nodeoffset[isubmesh][Jnp];
                    h_nodeoffset_skip_mid(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp+2]-nodeoffset[isubmesh][Jnp+1];
                    h_nodeoffset_skip_top(isubmesh,jsubmesh) = nodeoffset[isubmesh][Jnp+hc_submesh_x2[jsubmesh]->Nel] -
                                                                nodeoffset[isubmesh][Jnp+hc_submesh_x2[jsubmesh]->Nel-1];
                }
                else{
                    h_nodeoffset_start(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_bot(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_mid(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_top(isubmesh,jsubmesh) = -1;
                }
            }
        }

        for (int jsubmesh=0; jsubmesh<Ny+2; jsubmesh++){
            for (int isubmesh=0; isubmesh<Nx+2; isubmesh++){
                if (isubmesh==0 || jsubmesh==0 || isubmesh==Nx+1 || jsubmesh==Ny+1){
                    h_nodeoffset_start(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_bot(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_mid(isubmesh,jsubmesh) = -1;
                    h_nodeoffset_skip_top(isubmesh,jsubmesh) = -1;
                    h_elemoffset_start(isubmesh,jsubmesh) = -1;
                }
            }
            if (jsubmesh==0 || jsubmesh== Ny+1){
                h_elemoffset_skip(jsubmesh) = -1;
            }
        }

        for (int jsubmesh=0; jsubmesh<Ny+2; jsubmesh++){
            host_elemoffset_skip[jsubmesh] = h_elemoffset_skip(jsubmesh);
            for (int isubmesh=0; isubmesh<Nx+2; isubmesh++){
                host_elemoffset_start[isubmesh][jsubmesh] = h_elemoffset_start(isubmesh,jsubmesh);
                host_nodeoffset_start[isubmesh][jsubmesh] = h_nodeoffset_start(isubmesh,jsubmesh);
                host_nodeoffset_skip_bot[isubmesh][jsubmesh] = h_nodeoffset_skip_bot(isubmesh,jsubmesh);
                host_nodeoffset_skip_mid[isubmesh][jsubmesh] = h_nodeoffset_skip_mid(isubmesh,jsubmesh);
                host_nodeoffset_skip_top[isubmesh][jsubmesh] = h_nodeoffset_skip_top(isubmesh,jsubmesh);
            }
        }

        Kokkos::deep_copy(elemoffset_skip, h_elemoffset_skip);
        Kokkos::deep_copy(elemoffset_start, h_elemoffset_start);
        Kokkos::deep_copy(nodeoffset_start, h_nodeoffset_start);
        Kokkos::deep_copy(nodeoffset_skip_bot, h_nodeoffset_skip_bot);
        Kokkos::deep_copy(nodeoffset_skip_mid, h_nodeoffset_skip_mid);
        Kokkos::deep_copy(nodeoffset_skip_top, h_nodeoffset_skip_top);

    }
    else{
        Kokkos::View<int**>::HostMirror h_elemoffset_start = Kokkos::create_mirror_view(elemoffset_start);
        Kokkos::View<int*>::HostMirror h_elemoffset_skip = Kokkos::create_mirror_view(elemoffset_skip);
        Kokkos::View<int**>::HostMirror h_nodeoffset_start = Kokkos::create_mirror_view(nodeoffset_start);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_bot = Kokkos::create_mirror_view(nodeoffset_skip_bot);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_mid = Kokkos::create_mirror_view(nodeoffset_skip_mid);
        Kokkos::View<int**>::HostMirror h_nodeoffset_skip_top = Kokkos::create_mirror_view(nodeoffset_skip_top);

        for (int jsubmesh=0; jsubmesh<Ny+2; jsubmesh++){
            h_elemoffset_skip(jsubmesh) = 0;
            for (int isubmesh=0; isubmesh<Nx+2; isubmesh++){
                h_elemoffset_start(isubmesh,jsubmesh) = 0;
                h_nodeoffset_start(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_bot(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_mid(isubmesh,jsubmesh) = 0;
                h_nodeoffset_skip_top(isubmesh,jsubmesh) = 0;
            }
        }

        for (int jsubmesh=0; jsubmesh<Ny+2; jsubmesh++){
            host_elemoffset_skip[jsubmesh] = h_elemoffset_skip(jsubmesh);
            for (int isubmesh=0; isubmesh<Nx+2; isubmesh++){
                host_elemoffset_start[isubmesh][jsubmesh] = h_elemoffset_start(isubmesh,jsubmesh);
                host_nodeoffset_start[isubmesh][jsubmesh] = h_nodeoffset_start(isubmesh,jsubmesh);
                host_nodeoffset_skip_bot[isubmesh][jsubmesh] = h_nodeoffset_skip_bot(isubmesh,jsubmesh);
                host_nodeoffset_skip_mid[isubmesh][jsubmesh] = h_nodeoffset_skip_mid(isubmesh,jsubmesh);
                host_nodeoffset_skip_top[isubmesh][jsubmesh] = h_nodeoffset_skip_top(isubmesh,jsubmesh);
            }
        }

        Kokkos::deep_copy(elemoffset_skip, h_elemoffset_skip);
        Kokkos::deep_copy(elemoffset_start, h_elemoffset_start);
        Kokkos::deep_copy(nodeoffset_start, h_nodeoffset_start);
        Kokkos::deep_copy(nodeoffset_skip_bot, h_nodeoffset_skip_bot);
        Kokkos::deep_copy(nodeoffset_skip_mid, h_nodeoffset_skip_mid);
        Kokkos::deep_copy(nodeoffset_skip_top, h_nodeoffset_skip_top);
    }

    for (int i=0; i<Nx+2; i++){
        delete[] nodeoffset[i];
    }
    delete[] nodeoffset;
}

/**
* @brief Class constructor.
*
* @param[in] Array of x1-submesh object pointers on host
* @param[in] Number of x1-blocks in mesh
* @param[in] Array of x2-submesh object pointers on host
* @param[in] Number of x1-blocks in mesh
*/
MeshBdry::MeshBdry(SubmeshHostViewPtr hc_submesh_x1,
                   int Nx,
                   SubmeshHostViewPtr hc_submesh_x2,
                   int Ny,
                   bool** host_isactive){

    is_bdry_edge = Kokkos::View<bool*> ("is_bdry_edge", 2*Nx*Ny+Nx+Ny);
    bdry_edge_normal = Kokkos::View<Vector3*> ("bdry_edge_normal", 2*Nx*Ny+Nx+Ny);
    is_bdry_vert = Kokkos::View<bool*> ("is_bdry_vertex", Nx*Ny+Nx+Ny+1);
    bdry_vert_normal = Kokkos::View<Vector3*> ("bdry_vert_normal", Nx*Ny+Nx+Ny+1);
    edge_to_face = Kokkos::View<int*> ("Edge2Face", 2*Nx*Ny+Nx+Ny);

    Kokkos::View<bool*>::HostMirror h_is_bdry_edge = Kokkos::create_mirror_view(is_bdry_edge);
    Kokkos::View<Vector3*>::HostMirror h_bdry_edge_normal = Kokkos::create_mirror_view(bdry_edge_normal);
    Kokkos::View<bool*>::HostMirror h_is_bdry_vert = Kokkos::create_mirror_view(is_bdry_vert);
    Kokkos::View<Vector3*>::HostMirror h_bdry_vert_normal = Kokkos::create_mirror_view(bdry_vert_normal);
    Kokkos::View<int*>::HostMirror h_edge_to_face = Kokkos::create_mirror_view(edge_to_face);

    host_is_bdry_edge = new bool[2*Nx*Ny+Nx+Ny];
    host_bdry_edge_normal = new Vector3[2*Nx*Ny+Nx+Ny];
    host_is_bdry_vert = new bool[Nx*Ny+Nx+Ny+1];
    host_bdry_vert_normal = new Vector3[Nx*Ny+Nx+Ny+1];
    host_edge_to_face = new int[2*Nx*Ny+Nx+Ny];

    for (int jsubmesh=1; jsubmesh<=Ny; jsubmesh++){
        for (int isubmesh=1; isubmesh<=Nx; isubmesh++){
            if (jsubmesh==1){
                if (host_isactive[isubmesh][jsubmesh]){
                    h_is_bdry_edge(isubmesh-1) = true;
                    h_bdry_edge_normal(isubmesh-1) = Vector3(0.0,-1.0,0.0);
                    if (isubmesh==1){
                        h_is_bdry_edge(Nx) = true;
                        h_bdry_edge_normal(Nx) = Vector3(-1.0,0.0,0.0);
                    }
                    if (isubmesh==Nx){
                        h_is_bdry_edge(isubmesh+Nx) = true;
                        h_bdry_edge_normal(isubmesh+Nx) = Vector3(1.0,0.0,0.0);
                    }
                }
                else{
                    h_is_bdry_edge(isubmesh-1) = false;
                    h_bdry_edge_normal(isubmesh-1) = Vector3(0.0,0.0,0.0);
                    if (isubmesh==1){
                        h_is_bdry_edge(Nx) = false;
                        h_bdry_edge_normal(Nx) = Vector3(0.0,0.0,0.0);
                    }
                    if (isubmesh==Nx){
                        h_is_bdry_edge(isubmesh+Nx) = false;
                        h_bdry_edge_normal(isubmesh+Nx) = Vector3(0.0,0.0,0.0);
                    }
                }

                if (isubmesh>1){
                    int sum = host_isactive[isubmesh-1][jsubmesh]+host_isactive[isubmesh][jsubmesh];
                    if (sum == 1){
                        h_is_bdry_edge(isubmesh-1+Nx) = true;
                        if (host_isactive[isubmesh][jsubmesh]){
                            h_bdry_edge_normal(isubmesh-1+Nx) = Vector3(-1.0,0.0,0.0);

                        }
                        else{
                            h_bdry_edge_normal(isubmesh-1+Nx) = Vector3(1.0,0.0,0.0);
                        }
                    }
                    else{
                        h_is_bdry_edge(isubmesh-1+Nx) = false;
                        h_bdry_edge_normal(isubmesh-1+Nx) = Vector3(0.0,0.0,0.0);
                    }
                }
            }
            else{

                int sum = host_isactive[isubmesh][jsubmesh-1]+host_isactive[isubmesh][jsubmesh];
                if (sum == 1){
                    h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+isubmesh-1) = true;
                    if (host_isactive[isubmesh][jsubmesh]){
                        h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+isubmesh-1) = Vector3(0.0,-1.0,0.0);
                    }
                    else{
                        h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+isubmesh-1) = Vector3(0.0,1.0,0.0);
                    }
                }
                else{
                    h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+isubmesh-1) = false;
                    h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+isubmesh-1) = Vector3(0.0,0.0,0.0);
                }

                if (isubmesh==1){
                    if (host_isactive[isubmesh][jsubmesh]){
                        h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+Nx) = true;
                        h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+Nx) = Vector3(-1.0,0.0,0.0);

                    }
                    else{
                        h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+Nx) = false;
                        h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+Nx) = Vector3(0.0,0.0,0.0);
                    }
                }
                else{
                    sum = host_isactive[isubmesh-1][jsubmesh]+host_isactive[isubmesh][jsubmesh];
                    if (sum == 1){
                        h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+Nx+isubmesh-1) = true;
                        if (host_isactive[isubmesh][jsubmesh]){
                            h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+Nx+isubmesh-1) = Vector3(-1.0,0.0,0.0);
                        }
                        else{
                            h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+Nx+isubmesh-1) = Vector3(1.0,0.0,0.0);
                        }
                    }
                    else{
                        h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+Nx+isubmesh-1) = false;
                        h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+Nx+isubmesh-1) = Vector3(0.0,0.0,0.0);
                    }
                }

                if (isubmesh==Nx){
                    if (host_isactive[isubmesh][jsubmesh]){
                        h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+isubmesh+Nx) = true;
                        h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+isubmesh+Nx) = Vector3(1.0,0.0,0.0);
                    }
                    else{
                        h_is_bdry_edge((jsubmesh-1)*(2*Nx+1)+isubmesh+Nx) = false;
                        h_bdry_edge_normal((jsubmesh-1)*(2*Nx+1)+isubmesh+Nx) = Vector3(0.0,0.0,0.0);
                    }
                }
            }

            if (jsubmesh == Ny){
                if (host_isactive[isubmesh][jsubmesh]){
                    h_is_bdry_edge((jsubmesh)*(2*Nx+1)+isubmesh-1) = true;
                    h_bdry_edge_normal((jsubmesh)*(2*Nx+1)+isubmesh-1) = Vector3(0.0,1.0,0.0);
                }
                else{
                    h_is_bdry_edge((jsubmesh)*(2*Nx+1)+isubmesh-1) = false;
                    h_bdry_edge_normal((jsubmesh)*(2*Nx+1)+isubmesh-1) = Vector3(0.0,0.0,0.0);
                }
            }
        }
    }

    for (int jvert=0; jvert<Ny+1; jvert++){
        for (int ivert=0; ivert<Nx+1; ivert++){
            int VertID = jvert*(Nx+1)+ivert;
            int jsub = VertID/(Nx+1) + 1;
            int isub = VertID - (jsub-1)*(Nx+1) + 1;
            if (ivert==0){
                if (jvert==0){
                    if (host_isactive[isub][jsub]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(-1.0,-1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
                else if (jvert==Ny){
                    if (host_isactive[isub][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(-1.0,1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
                else{
                    if (host_isactive[isub][jsub] && host_isactive[isub][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(-1.0,0.0,0.0);
                    }
                    else if (host_isactive[isub][jsub] && !host_isactive[isub][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(-1.0,-1.0,0.0);
                    }
                    else if (!host_isactive[isub][jsub] && host_isactive[isub][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(-1.0,1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
            }
            else if (ivert==Nx){
                if (jvert==0){
                    if (host_isactive[isub-1][jsub]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(1.0,-1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
                else if (jvert==Ny){
                    if (host_isactive[isub-1][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(1.0,1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
                else{
                    if (host_isactive[isub-1][jsub] && host_isactive[isub-1][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(1.0,0.0,0.0);
                    }
                    else if (host_isactive[isub-1][jsub] && !host_isactive[isub-1][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(1.0,-1.0,0.0);
                    }
                    else if (!host_isactive[isub-1][jsub] && host_isactive[isub-1][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(1.0,1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
            }
            else{
                if (jvert==0){
                    if (host_isactive[isub][jsub] && host_isactive[isub-1][jsub]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,-1.0,0.0);
                    }
                    else if (host_isactive[isub][jsub] && !host_isactive[isub-1][jsub]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(-1.0,-1.0,0.0);
                    }
                    else if (!host_isactive[isub][jsub] && host_isactive[isub-1][jsub]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(1.0,-1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
                else if (jvert==Ny){
                    if (host_isactive[isub][jsub-1] && host_isactive[isub-1][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,1.0,0.0);
                    }
                    else if (host_isactive[isub][jsub-1] && !host_isactive[isub-1][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(-1.0,1.0,0.0);
                    }
                    else if (!host_isactive[isub][jsub-1] && host_isactive[isub-1][jsub-1]){
                        h_is_bdry_vert(VertID) = true;
                        h_bdry_vert_normal(VertID) = Vector3(1.0,1.0,0.0);
                    }
                    else{
                        h_is_bdry_vert(VertID) = false;
                        h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);
                    }
                }
                else{
                    bool tr = host_isactive[isub][jsub];
                    bool tl = host_isactive[isub-1][jsub];
                    bool bl = host_isactive[isub-1][jsub-1];
                    bool br = host_isactive[isub][jsub-1];
                    int sum = bl+br+tl+tr;
                    if (sum==0 || sum==4) {
                        h_is_bdry_vert(VertID) = false;
                    }
                    else{
                        h_is_bdry_vert(VertID) = true;
                    }

                    if      ( tr &&  tl &&  bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//1
                    else if (!tr &&  tl &&  bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//2
                    else if ( tr && !tl &&  bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//3
                    else if ( tr &&  tl && !bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//4
                    else if ( tr &&  tl &&  bl && !br ) h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//5
                    else if ( tr && !tl &&  bl && !br ) h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//6
                    else if (!tr &&  tl && !bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//7
                    else if ( tr &&  tl && !bl && !br ) h_bdry_vert_normal(VertID) = Vector3(0.0,-1.0,0.0);//8
                    else if (!tr && !tl &&  bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(0.0,1.0,0.0);//9
                    else if (!tr &&  tl &&  bl && !br ) h_bdry_vert_normal(VertID) = Vector3(1.0,0.0,0.0);//10
                    else if ( tr && !tl && !bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(-1.0,0.0,0.0);//11
                    else if ( tr && !tl && !bl && !br ) h_bdry_vert_normal(VertID) = Vector3(-1.0,-1.0,0.0);//12
                    else if (!tr &&  tl && !bl && !br ) h_bdry_vert_normal(VertID) = Vector3(1.0,-1.0,0.0);//13
                    else if (!tr && !tl &&  bl && !br ) h_bdry_vert_normal(VertID) = Vector3(1.0,1.0,0.0);//14
                    else if (!tr && !tl && !bl &&  br ) h_bdry_vert_normal(VertID) = Vector3(-1.0,1.0,0.0);//15
                    else h_bdry_vert_normal(VertID) = Vector3(0.0,0.0,0.0);//16
                }
            }
        }
    }


    Kokkos::deep_copy(is_bdry_edge, h_is_bdry_edge);
    Kokkos::deep_copy(bdry_edge_normal, h_bdry_edge_normal);
    Kokkos::deep_copy(is_bdry_vert, h_is_bdry_vert);
    Kokkos::deep_copy(bdry_vert_normal, h_bdry_vert_normal);

    Nbdry_faces = 0;
    for (int iedge=0; iedge<2*Nx*Ny+Nx+Ny; iedge++){
        int Nx2p1 = 2*Nx+1;
        if (h_is_bdry_edge(iedge)){
            h_edge_to_face(iedge) = Nbdry_faces;
            int num = iedge/Nx2p1;
            int rem = iedge - num*Nx2p1;
            if (rem < Nx){
                Nbdry_faces += hc_submesh_x1[rem+1]->Nel;
            }
            else{
                Nbdry_faces += hc_submesh_x2[num+1]->Nel;
            }
        }
        else{
            h_edge_to_face(iedge) = -1;
        }
        host_edge_to_face[iedge] = h_edge_to_face(iedge);
        host_is_bdry_edge[iedge] = h_is_bdry_edge(iedge);
        host_bdry_edge_normal[iedge] = h_bdry_edge_normal(iedge);
    }
    Kokkos::deep_copy(edge_to_face, h_edge_to_face);

    for (int ivert=0; ivert<Nx*Ny+Nx+Ny+1; ivert++){
        host_is_bdry_vert[ivert] = h_is_bdry_vert(ivert);
        host_bdry_vert_normal[ivert] = Vector3(h_bdry_vert_normal(ivert)[0],
                                               h_bdry_vert_normal(ivert)[1],
                                               h_bdry_vert_normal(ivert)[2]);
    }
}

/**
* @brief Class constructor.
*
* @param[in] Array of x1-submesh object pointers on host
* @param[in] Number of x1-blocks in mesh
*/
BlockInterface::BlockInterface(SubmeshHostViewPtr hc_submesh_x1,
                               int Nx){

    host_if_x1_r = new double[Nx-1];
    host_if_x1_node = new int[Nx+1];
    host_if_x2_r = nullptr;
    host_if_x2_node = nullptr;
    host_vert_nodeID = nullptr;
    host_vert_subID = nullptr;
    host_edge_first_nodeID = nullptr;
    host_edge_subID = nullptr;

    if_x1_r = Kokkos::View<double*> ("interface_x1_grading", Nx-1);
    if_x1_node = Kokkos::View<int*> ("interface_x1_node", Nx+1);

    host_if_x1_node[0] = 0;
    for (int i=1; i<Nx; i++){
        double min_elem = 999.0;
        double max_elem = 999.0;
        if (hc_submesh_x1[i]->meshtype & minBL)
            min_elem = (hc_submesh_x1[i]->t0)*pow(hc_submesh_x1[i]->r,hc_submesh_x1[i]->Nel-1);
        if (hc_submesh_x1[i]->meshtype & maxBL)
            min_elem = hc_submesh_x1[i]->t0;
        if (hc_submesh_x1[i]->meshtype & uniform)
            min_elem = hc_submesh_x1[i]->t0;
        if (hc_submesh_x1[i]->meshtype & arbitrary){
            int Nel = hc_submesh_x1[i]->Nel;
            min_elem = hc_submesh_x1[i]->host_BL_coords[Nel]-hc_submesh_x1[i]->host_BL_coords[Nel-1];
        }

        if (hc_submesh_x1[i+1]->meshtype & minBL)
            max_elem = hc_submesh_x1[i+1]->t0;
        if (hc_submesh_x1[i+1]->meshtype & maxBL)
            max_elem = (hc_submesh_x1[i+1]->t0)*pow(hc_submesh_x1[i+1]->r,hc_submesh_x1[i+1]->Nel-1);
        if (hc_submesh_x1[i+1]->meshtype & uniform)
            max_elem = hc_submesh_x1[i+1]->t0;
        if (hc_submesh_x1[i+1]->meshtype & arbitrary)
            max_elem = hc_submesh_x1[i+1]->host_BL_coords[1]-hc_submesh_x1[i+1]->host_BL_coords[0];

        host_if_x1_r[i-1] = max_elem/min_elem;
        host_if_x1_node[i] = hc_submesh_x1[i+1]->Nel_cumulative;
    }
    host_if_x1_node[Nx] = host_if_x1_node[Nx-1] + hc_submesh_x1[Nx]->Nel;

    Kokkos::View<double*>::HostMirror h_if_x1_r = Kokkos::create_mirror_view(if_x1_r);
    Kokkos::View<int*>::HostMirror h_if_x1_node = Kokkos::create_mirror_view(if_x1_node);
    h_if_x1_node(0) = host_if_x1_node[0];
    for (int i=1; i<Nx; i++){
        h_if_x1_node(i) = host_if_x1_node[i];
        h_if_x1_r(i-1) = host_if_x1_r[i-1];
    }
    h_if_x1_node(Nx) = host_if_x1_node[Nx];

    Kokkos::deep_copy(if_x1_r,h_if_x1_r);
    Kokkos::deep_copy(if_x1_node,h_if_x1_node);
}

/**
* @brief Class constructor.
*
* @param[in] Array of x1-submesh object pointers on host
* @param[in] Number of x1-blocks in mesh
* @param[in] Array of x2-submesh object pointers on host
* @param[in] Number of x1-blocks in mesh
*/
BlockInterface::BlockInterface(SubmeshHostViewPtr hc_submesh_x1,
                               int Nx,
                               SubmeshHostViewPtr hc_submesh_x2,
                               int Ny,
                               bool** host_isactive){
    host_if_x1_r = new double[Nx-1];
    host_if_x1_node = new int[Nx+1];
    host_if_x2_r = new double[Ny-1];
    host_if_x2_node = new int[Ny+1];

    if_x1_r = Kokkos::View<double*> ("interface_x1_grading", Nx-1);
    if_x1_node = Kokkos::View<int*> ("interface_x1_node", Nx+1);
    if_x2_r = Kokkos::View<double*> ("interface_x2_grading", Ny-1);
    if_x2_node = Kokkos::View<int*> ("interface_x2_node", Ny+1);

    host_if_x1_node[0] = 0;
    for (int i=1; i<Nx; i++){
       double min_elem = 999.0;
       double max_elem = 999.0;
       if (hc_submesh_x1[i]->meshtype & minBL)
           min_elem = (hc_submesh_x1[i]->t0)*pow(hc_submesh_x1[i]->r,hc_submesh_x1[i]->Nel-1);
       if (hc_submesh_x1[i]->meshtype & maxBL)
           min_elem = hc_submesh_x1[i]->t0;
       if (hc_submesh_x1[i]->meshtype & uniform)
           min_elem = hc_submesh_x1[i]->t0;
       if (hc_submesh_x1[i]->meshtype & arbitrary){
           int Nel = hc_submesh_x1[i]->Nel;
           min_elem = hc_submesh_x1[i]->host_BL_coords[Nel]-hc_submesh_x1[i]->host_BL_coords[Nel-1];
       }

       if (hc_submesh_x1[i+1]->meshtype & minBL)
           max_elem = hc_submesh_x1[i+1]->t0;
       if (hc_submesh_x1[i+1]->meshtype & maxBL)
           max_elem = (hc_submesh_x1[i+1]->t0)*pow(hc_submesh_x1[i+1]->r,hc_submesh_x1[i+1]->Nel-1);
       if (hc_submesh_x1[i+1]->meshtype & uniform)
           max_elem = hc_submesh_x1[i+1]->t0;
       if (hc_submesh_x1[i+1]->meshtype & arbitrary)
           max_elem = hc_submesh_x1[i+1]->host_BL_coords[1]-hc_submesh_x1[i+1]->host_BL_coords[0];

       host_if_x1_r[i-1] = max_elem/min_elem;
       host_if_x1_node[i] = hc_submesh_x1[i+1]->Nel_cumulative;
    }
    host_if_x1_node[Nx] = host_if_x1_node[Nx-1] + hc_submesh_x1[Nx]->Nel;

    host_if_x2_node[0] = 0;
    for (int i=1; i<Ny; i++){
       double min_elem = 999.0;
       double max_elem = 999.0;
       if (hc_submesh_x2[i]->meshtype & minBL)
           min_elem = (hc_submesh_x2[i]->t0)*pow(hc_submesh_x2[i]->r,hc_submesh_x2[i]->Nel-1);
       if (hc_submesh_x2[i]->meshtype & maxBL)
           min_elem = hc_submesh_x2[i]->t0;
       if (hc_submesh_x2[i]->meshtype & uniform)
           min_elem = hc_submesh_x2[i]->t0;
       if (hc_submesh_x2[i]->meshtype & arbitrary){
           int Nel = hc_submesh_x2[i]->Nel;
           min_elem = hc_submesh_x2[i]->host_BL_coords[Nel]-hc_submesh_x2[i]->host_BL_coords[Nel-1];
       }

       if (hc_submesh_x2[i+1]->meshtype & minBL)
           max_elem = hc_submesh_x2[i+1]->t0;
       if (hc_submesh_x2[i+1]->meshtype & maxBL)
           max_elem = (hc_submesh_x2[i+1]->t0)*pow(hc_submesh_x2[i+1]->r,hc_submesh_x2[i+1]->Nel-1);
       if (hc_submesh_x2[i+1]->meshtype & uniform)
           max_elem = hc_submesh_x2[i+1]->t0;
       if (hc_submesh_x2[i+1]->meshtype & arbitrary)
           max_elem = hc_submesh_x2[i+1]->host_BL_coords[1]-hc_submesh_x2[i+1]->host_BL_coords[0];

       host_if_x2_r[i-1] = max_elem/min_elem;
       host_if_x2_node[i] = hc_submesh_x2[i+1]->Nel_cumulative;
    }
    host_if_x2_node[Ny] = host_if_x2_node[Ny-1] + hc_submesh_x2[Ny]->Nel;

    Kokkos::View<double*>::HostMirror h_if_x1_r = Kokkos::create_mirror_view(if_x1_r);
    Kokkos::View<int*>::HostMirror h_if_x1_node = Kokkos::create_mirror_view(if_x1_node);
    h_if_x1_node(0) = host_if_x1_node[0];
    for (int i=1; i<Nx; i++){
       h_if_x1_node(i) = host_if_x1_node[i];
       h_if_x1_r(i-1) = host_if_x1_r[i-1];
    }
    h_if_x1_node(Nx) = host_if_x1_node[Nx];

    Kokkos::deep_copy(if_x1_r,h_if_x1_r);
    Kokkos::deep_copy(if_x1_node,h_if_x1_node);

    Kokkos::View<double*>::HostMirror h_if_x2_r = Kokkos::create_mirror_view(if_x2_r);
    Kokkos::View<int*>::HostMirror h_if_x2_node = Kokkos::create_mirror_view(if_x2_node);
    h_if_x2_node(0) = host_if_x2_node[0];
    for (int i=1; i<Nx; i++){
       h_if_x2_node(i) = host_if_x2_node[i];
       h_if_x2_r(i-1) = host_if_x2_r[i-1];
    }
    h_if_x2_node(Nx) = host_if_x2_node[Nx];

    Kokkos::deep_copy(if_x2_r,h_if_x2_r);
    Kokkos::deep_copy(if_x2_node,h_if_x2_node);

    // block-vertices
    host_vert_subID = new int[Nx*Ny+Nx+Ny+1];
    vert_subID = Kokkos::View<int*>("vert_subID",Nx*Ny+Nx+Ny+1);
    Kokkos::View<int*>::HostMirror h_vert_subID = Kokkos::create_mirror_view(vert_subID);
    // host_vert_nodeID = new int[Nx*Ny+Nx+Ny+1];

    for (int vertID=0; vertID<Nx*Ny+Nx+Ny+1; vertID++){
        int jsub = 1 + vertID/(Nx+1);
        int isub = vertID - (Nx+1)*(jsub-1) + 1;

        if (host_isactive[isub-1][jsub-1]){ // bottom-left
            host_vert_subID[vertID] = (isub-2) + (jsub-2)*Nx;
        }
        else if (host_isactive[isub][jsub-1]){// bottom-right
            host_vert_subID[vertID] = (isub-1) + (jsub-2)*Nx;
        }
        else if (host_isactive[isub-1][jsub]){// top-left
            host_vert_subID[vertID] = (isub-2) + (jsub-1)*Nx;
        }
        else if (host_isactive[isub][jsub]){
            host_vert_subID[vertID] = (isub-1) + (jsub-1)*Nx;
        }
        else {
            host_vert_subID[vertID] = -1;
        }

        h_vert_subID(vertID) = host_vert_subID[vertID];
    }
    Kokkos::deep_copy(vert_subID,h_vert_subID);

    // block-edges
    host_edge_subID = new int[2*Nx*Ny+Nx+Ny];
    edge_subID = Kokkos::View<int*>("edge_subID",2*Nx*Ny+Nx+Ny);
    Kokkos::View<int*>::HostMirror h_edge_subID = Kokkos::create_mirror_view(edge_subID);
    // host_edge_nodeID = new int[2*Nx*Ny+Nx+Ny];

    for (int edgeID=0; edgeID<2*Nx*Ny+Nx+Ny; edgeID++){
        int num = edgeID/(2*Nx+1);
        int rem = edgeID - num*(2*Nx+1);

        if (rem < Nx){
            int jsub = 1 + edgeID/(2*Nx+1);
            int isub = edgeID - (2*Nx+1)*(jsub-1) + 1;

            if (host_isactive[isub][jsub-1]){ // bottom
                host_edge_subID[edgeID] = (isub-1) + (jsub-2)*Nx;
            }
            else if (host_isactive[isub][jsub]){// top
                host_edge_subID[edgeID] = (isub-1) + (jsub-1)*Nx;
            }
            else {
                host_edge_subID[edgeID] = -1;
            }
        }
        else{
            int jsub = 1 + edgeID/(2*Nx+1);
            int isub = edgeID - (2*Nx+1)*(jsub-1) - Nx + 1;

            if (host_isactive[isub-1][jsub]){ // left
                host_edge_subID[edgeID] = (isub-2) + (jsub-1)*Nx;
            }
            else if (host_isactive[isub][jsub]){// right
                host_edge_subID[edgeID] = (isub-1) + (jsub-1)*Nx;
            }
            else {
                host_edge_subID[edgeID] = -1;
            }
        }

        h_edge_subID(edgeID) = host_edge_subID[edgeID];
    }
    Kokkos::deep_copy(edge_subID,h_edge_subID);

}

/**
 * @brief sets global node IDs for vertex nodes and first nodes on edge
 * @param[in] array of vertex node IDs
 * @param[in] array of edge first node IDs
 * @param[in] number of x1-blocks in mesh
 * @param[in] number of x2-blocks in mesh
 */
void BlockInterface::set_interface_nodeIDs(int *vert_nodeIDs,
                                           int *edge_first_nodeIDs,
                                           int Nx,
                                           int Ny){
    host_vert_nodeID = new int[Nx*Ny+Nx+Ny+1];
    host_edge_first_nodeID = new int[2*Nx*Ny+Nx+Ny];
    vert_nodeID = Kokkos::View<int*>("vert_nodeID",Nx*Ny+Nx+Ny+1);
    Kokkos::View<int*>::HostMirror h_vert_nodeID = Kokkos::create_mirror_view(vert_nodeID);
    edge_first_nodeID = Kokkos::View<int*>("edge_first_nodeID",2*Nx*Ny+Nx+Ny);
    Kokkos::View<int*>::HostMirror h_edge_first_nodeID = Kokkos::create_mirror_view(edge_first_nodeID);

    for (int vertID=0; vertID < Nx*Ny+Nx+Ny+1; vertID++){
        host_vert_nodeID[vertID] = vert_nodeIDs[vertID];
        h_vert_nodeID(vertID) = vert_nodeIDs[vertID];
    }
    Kokkos::deep_copy(vert_nodeID,h_vert_nodeID);
    for (int edgeID=0; edgeID < 2*Nx*Ny+Nx+Ny; edgeID++){
        host_edge_first_nodeID[edgeID] = edge_first_nodeIDs[edgeID];
        h_edge_first_nodeID(edgeID) = edge_first_nodeIDs[edgeID];
    }
    Kokkos::deep_copy(edge_first_nodeID,h_edge_first_nodeID);

}

/**
 * @brief intialize all interface node IDs (for edge/vertex nodes)
 * @param[in] Object of the wrapper mesh structure
 * @return  Updated object of the wrapper mesh structure
 */
MBBL initialize_interface_nodeIDs(MBBL pumi_obj){
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int Ny = pumi_obj.mesh.nsubmesh_x2;

    int *vert_nodeIDs = new int[Nx*Ny+Nx+Ny+1];
    int *edge_nodeIDs = new int[2*Nx*Ny+Nx+Ny];

    for (int jvert=0; jvert<Ny+1; jvert++){
        for (int ivert=0; ivert<Nx+1; ivert++){
            int vertID = ivert + jvert*(Nx+1);

            int inp = get_x1_nodeID_at_interface_host(pumi_obj,ivert);
            int jnp = get_x2_nodeID_at_interface_host(pumi_obj,jvert);

            if (pumi_obj.mesh.blkif.host_vert_subID[vertID]+1){
                vert_nodeIDs[vertID] = get_global_nodeID_2D(pumi_obj,inp,jnp);
            }
            else{
                vert_nodeIDs[vertID] = -1;
            }
        }
    }

    //horizontal block-edges
    for (int jvert=0; jvert<Ny+1; jvert++){
        for (int iedge=0; iedge<Nx; iedge++){
            int edgeID = iedge + jvert*(2*Nx+1);
            int inp = get_x1_nodeID_at_interface_host(pumi_obj,iedge);
            int jnp = get_x2_nodeID_at_interface_host(pumi_obj,jvert);
            if (pumi_obj.mesh.blkif.host_edge_subID[edgeID]+1){
                edge_nodeIDs[edgeID] = get_global_nodeID_2D(pumi_obj,inp+1,jnp);
            }
            else{
                edge_nodeIDs[edgeID] = -1;
            }
        }
    }

    //vertical block-edges
    for (int jedge=0; jedge<Ny; jedge++){
        for (int ivert=0; ivert<Nx+1; ivert++){
            int edgeID = Nx + jedge*(2*Nx+1) + ivert;
            int inp = get_x1_nodeID_at_interface_host(pumi_obj,ivert);
            int jnp = get_x2_nodeID_at_interface_host(pumi_obj,jedge);
            if (pumi_obj.mesh.blkif.host_edge_subID[edgeID]+1){
                edge_nodeIDs[edgeID] = get_global_nodeID_2D(pumi_obj,inp,jnp+1);
            }
            else{
                edge_nodeIDs[edgeID] = -1;
            }
        }
    }

    pumi_obj.mesh.blkif.set_interface_nodeIDs(vert_nodeIDs,edge_nodeIDs,Nx,Ny);

    delete vert_nodeIDs;
    delete edge_nodeIDs;

    return pumi_obj;
}

/**
 * @brief intialize data stuctures in class for BST searches
 * @param[in] object to BlockInterface class
 * @param[in] Array of x1-submesh object pointers on host
 * @param[in] Number of x1-blocks in mesh
 * @param[in] Array of x2-submesh object pointers on host
 * @param[in] Number of x1-blocks in mesh
 * @param[in] 2D array of block activity
 */
void MeshBST::initialize_MeshBST(BlockInterface blkif,
                                    SubmeshHostViewPtr hc_submesh_x1,
                                    int Nx,
                                    SubmeshHostViewPtr hc_submesh_x2,
                                    int Ny,
                                    bool **host_isactive){

    total_active_blocks = 0;
    for (int jsub=1; jsub<=Ny; jsub++){
        for (int isub=1; isub<=Nx; isub++){
            if (host_isactive[isub][jsub]){
                total_active_blocks++;
            }
        }
    }
    host_active_blockID = new int[total_active_blocks];
    host_block_nodes_cumulative = new int[total_active_blocks];
    host_block_elems_cumulative = new int[total_active_blocks];
    active_blockID = Kokkos::View<int*>("active-blockIDs",total_active_blocks);
    block_nodes_cumulative = Kokkos::View<int*>("block-nodes-cumulative",total_active_blocks);
    block_elems_cumulative = Kokkos::View<int*>("block-elems-cumulative",total_active_blocks);
    Kokkos::View<int*>::HostMirror h_active_blockID = Kokkos::create_mirror_view(active_blockID);
    Kokkos::View<int*>::HostMirror h_block_nodes_cumulative = Kokkos::create_mirror_view(block_nodes_cumulative);
    Kokkos::View<int*>::HostMirror h_block_elems_cumulative = Kokkos::create_mirror_view(block_elems_cumulative);

    int blkID = 0;
    for (int jsub=1; jsub<=Ny; jsub++){
        for (int isub=1; isub<=Nx; isub++){
            if (host_isactive[isub][jsub]){
                host_active_blockID[blkID] = (isub-1) + (jsub-1)*Nx;
                if (blkID==0){
                    host_block_nodes_cumulative[blkID] = (hc_submesh_x1[isub]->Nel-1)*(hc_submesh_x2[jsub]->Nel-1);
                    host_block_elems_cumulative[blkID] = (hc_submesh_x1[isub]->Nel)*(hc_submesh_x2[jsub]->Nel);
                }
                else{
                    host_block_nodes_cumulative[blkID] = host_block_nodes_cumulative[blkID-1] +
                                                        (hc_submesh_x1[isub]->Nel-1)*(hc_submesh_x2[jsub]->Nel-1);
                    host_block_elems_cumulative[blkID] = host_block_elems_cumulative[blkID-1] +
                                                        (hc_submesh_x1[isub]->Nel)*(hc_submesh_x2[jsub]->Nel);
                }
                h_block_nodes_cumulative(blkID) = host_block_nodes_cumulative[blkID];
                h_block_elems_cumulative(blkID) = host_block_elems_cumulative[blkID];
                h_active_blockID(blkID) = (isub-1) + (jsub-1)*Nx;
                blkID++;
            }
        }
    }
    Kokkos::deep_copy(active_blockID, h_active_blockID);
    Kokkos::deep_copy(block_nodes_cumulative, h_block_nodes_cumulative);
    Kokkos::deep_copy(block_elems_cumulative, h_block_elems_cumulative);
    total_block_nodes = h_block_nodes_cumulative(total_active_blocks-1);
    total_active_edges = 0;
    for(int iEdge=0; iEdge<2*Nx*Ny+Nx+Ny; iEdge++){
        if (blkif.host_edge_subID[iEdge]+1){
            total_active_edges++;
        }
    }

    host_active_edgeID = new int[total_active_edges];
    host_edge_nodes_cumulative = new int[total_active_edges];
    active_edgeID = Kokkos::View<int*>("active-edgeIDs",total_active_edges);
    edge_nodes_cumulative = Kokkos::View<int*>("edge-nodes-cumulative",total_active_edges);
    Kokkos::View<int*>::HostMirror h_active_edgeID = Kokkos::create_mirror_view(active_edgeID);
    Kokkos::View<int*>::HostMirror h_edge_nodes_cumulative = Kokkos::create_mirror_view(edge_nodes_cumulative);


    int edgeID = 0;
    for(int iEdge=0; iEdge<2*Nx*Ny+Nx+Ny; iEdge++){
        if (blkif.host_edge_subID[iEdge]+1){
            host_active_edgeID[edgeID] = iEdge;

            int num = iEdge/(2*Nx+1);
            int rem = iEdge - num*(2*Nx+1);
            int nel_edge;
            if (rem < Nx){
                nel_edge = hc_submesh_x1[rem+1]->Nel-1;
            }
            else{
                nel_edge = hc_submesh_x2[num+1]->Nel-1;
            }
            if (edgeID==0){
                host_edge_nodes_cumulative[edgeID] = nel_edge;
            }
            else{
                host_edge_nodes_cumulative[edgeID] = host_edge_nodes_cumulative[edgeID-1] +
                                                    nel_edge;
            }
            h_edge_nodes_cumulative(edgeID) = host_edge_nodes_cumulative[edgeID];
            h_active_edgeID(edgeID) = iEdge;
            edgeID++;
        }
    }
    Kokkos::deep_copy(active_edgeID,h_active_edgeID);
    Kokkos::deep_copy(edge_nodes_cumulative,h_edge_nodes_cumulative);
    total_edge_nodes = h_edge_nodes_cumulative(total_active_edges-1);
}

/**
 * @brief intializes MBBL mesh from user inputs and options
 * @param[in] Object of the wrapper mesh structure
 * @param[in] object to structure containing user options for the mesh
 * @return  Object of the wrapper mesh structure
 */
MBBL initialize_MBBL_mesh(Mesh_Inputs* pumi_inputs, Mesh_Options pumi_options){
    pumi::Mesh mesh;
    pumi::SubmeshInit x1_sub_obj;
    pumi::SubmeshInit x2_sub_obj;

    MBBL pumi_obj;

    if (pumi_inputs->ndim == 1){
        x1_sub_obj = pumi::submesh_initialize(pumi_inputs, pumi_options, pumi::x1_dir);
        mesh = pumi::mesh_initialize(pumi_inputs, pumi_options, x1_sub_obj.submesh, x1_sub_obj.host_submesh);
        pumi_obj = MBBL(mesh, x1_sub_obj.submesh, x1_sub_obj.host_submesh);
    }
    else if (pumi_inputs->ndim == 2){
        x1_sub_obj = pumi::submesh_initialize(pumi_inputs, pumi_options, pumi::x1_dir);
        x2_sub_obj = pumi::submesh_initialize(pumi_inputs, pumi_options, pumi::x2_dir);
        mesh = pumi::mesh_initialize(pumi_inputs, pumi_options, x1_sub_obj.submesh, x1_sub_obj.host_submesh, x2_sub_obj.submesh, x2_sub_obj.host_submesh);

        pumi_obj = MBBL(mesh, x1_sub_obj.submesh, x1_sub_obj.host_submesh, x2_sub_obj.submesh, x2_sub_obj.host_submesh);

        if (pumi_options.print_mesh_connectivity_option){
            pumi::print_2D_node_coordinates(pumi_obj);
            pumi::print_2D_node_elem_connectivity(pumi_obj);
        }
        pumi_obj = initialize_interface_nodeIDs(pumi_obj);
    }
    return pumi_obj;
}


} // namespace pumi
