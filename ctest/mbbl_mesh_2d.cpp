#include "pumiMBBLGPU.hpp"

int main( int argc, char* argv[] )
{
    bool test_passed = true;
    Kokkos::initialize( argc, argv );
    {
        pumi::Mesh_Inputs *pumi_inputs = pumi::inputs_allocate();
        pumi_inputs->ndim = 2;
        pumi_inputs->nsubmesh_x1 = 3;
        pumi_inputs->domain_x1_min = 1.5;

        std::string dummyfile;
        // submesh-1
        (pumi_inputs->meshtype_x1).push_back("maxBL");
        (pumi_inputs->block_length_x1).push_back(10.0);
        (pumi_inputs->max_elem_size_x1).push_back(3.0);
        (pumi_inputs->min_elem_size_x1).push_back(1.0);
        (pumi_inputs->arbitrary_x1_elemsize_file).push_back(dummyfile);
        // submesh-2
        (pumi_inputs->meshtype_x1).push_back("minBL");
        (pumi_inputs->block_length_x1).push_back(7.0);
        (pumi_inputs->max_elem_size_x1).push_back(2.0);
        (pumi_inputs->min_elem_size_x1).push_back(0.5);
        (pumi_inputs->arbitrary_x1_elemsize_file).push_back(dummyfile);
        // submesh-3
        (pumi_inputs->meshtype_x1).push_back("uniform");
        (pumi_inputs->block_length_x1).push_back(15.0);
        (pumi_inputs->max_elem_size_x1).push_back(3.0);
        (pumi_inputs->min_elem_size_x1).push_back(3.0);
        (pumi_inputs->arbitrary_x1_elemsize_file).push_back(dummyfile);

        pumi_inputs->nsubmesh_x2 = 4;
        pumi_inputs->domain_x2_min = -10.5;
        // std::vector<double> es = {-10.5, -8.75,-6.15,-5.09,-4.71,-2.66,-1.0,0.01};
        std::vector<double> es = {1.75, 2.6,1.06,0.38,2.05,1.66,1.01};

        FILE *coord_file;
        char coord_filename[30];
        sprintf(coord_filename,"arb-x2-elemsize.dat");
        coord_file = fopen(coord_filename,"w");
        for (int i=0; i< (int) es.size(); i++){
            fprintf(coord_file, "%.16e \n", es[i]);
        }
        fclose(coord_file);
        // submesh-1
        (pumi_inputs->meshtype_x2).push_back("arbitrary");
        (pumi_inputs->block_length_x2).push_back(10.51);
        (pumi_inputs->max_elem_size_x2).push_back(3.0);
        (pumi_inputs->min_elem_size_x2).push_back(1.0);
        (pumi_inputs->arbitrary_x2_elemsize_file).push_back("arb-x2-elemsize.dat");
        // submesh-2
        (pumi_inputs->meshtype_x2).push_back("minBL");
        (pumi_inputs->block_length_x2).push_back(6.99);
        (pumi_inputs->max_elem_size_x2).push_back(2.0);
        (pumi_inputs->min_elem_size_x2).push_back(0.5);
        (pumi_inputs->arbitrary_x2_elemsize_file).push_back(dummyfile);
        // submesh-3
        (pumi_inputs->meshtype_x2).push_back("uniform");
        (pumi_inputs->block_length_x2).push_back(5.5);
        (pumi_inputs->max_elem_size_x2).push_back(0.5);
        (pumi_inputs->min_elem_size_x2).push_back(0.5);
        (pumi_inputs->arbitrary_x2_elemsize_file).push_back(dummyfile);
        // submesh-4
        (pumi_inputs->meshtype_x2).push_back("uniform");
        (pumi_inputs->block_length_x2).push_back(5.0);
        (pumi_inputs->max_elem_size_x2).push_back(0.5);
        (pumi_inputs->min_elem_size_x2).push_back(0.5);
        (pumi_inputs->arbitrary_x2_elemsize_file).push_back(dummyfile);

        pumi_inputs->isactive[0][0] = true;
        pumi_inputs->isactive[1][0] = false;
        pumi_inputs->isactive[2][0] = false;
        pumi_inputs->isactive[0][1] = true;
        pumi_inputs->isactive[1][1] = true;
        pumi_inputs->isactive[2][1] = true;
        pumi_inputs->isactive[0][2] = false;
        pumi_inputs->isactive[1][2] = false;
        pumi_inputs->isactive[2][2] = true;
        pumi_inputs->isactive[0][3] = false;
        pumi_inputs->isactive[1][3] = true;
        pumi_inputs->isactive[2][3] = true;


        pumi::Mesh_Options pumi_options;

        pumi::MBBL pumi_obj = pumi::initialize_MBBL_mesh(pumi_inputs, pumi_options);
        pumi::inputs_deallocate(pumi_inputs);

        int knode=0;
        int error=0;
        for (int jnp=0; jnp<=pumi_obj.mesh.Nel_tot_x2; jnp++){
            for (int inp=0; inp<=pumi_obj.mesh.Nel_tot_x1; inp++){
                bool on_bdry, in_domain;
                int bdry_tag, bdry_dim;
                pumi::where_is_node(pumi_obj, inp, jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    int nodeID = get_global_nodeID_2D(pumi_obj, inp, jnp);
                    error += (nodeID-knode);
                    knode++;
                }
            }
        }

        if (error != 0) test_passed = false;
        pumi::free_mbbl(pumi_obj);
    }

    Kokkos::finalize();
    if (test_passed) return 0;
    else return 1;
}
