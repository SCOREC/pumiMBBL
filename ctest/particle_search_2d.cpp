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

        pumi_inputs->nsubmesh_x2 = 2;
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

        int ksubmesh=0;
        for(int jsubmesh=0; jsubmesh<pumi_inputs->nsubmesh_x2; jsubmesh++){
            for(int isubmesh=0; isubmesh<pumi_inputs->nsubmesh_x1; isubmesh++){
                ksubmesh++;
                pumi_inputs->isactive[isubmesh][jsubmesh] = true;
            }
        }

        pumi::Mesh_Options pumi_options;

        pumi::MBBL pumi_obj = pumi::initialize_MBBL_mesh(pumi_inputs, pumi_options);
        pumi::inputs_deallocate(pumi_inputs);

        int N_part = 5;
        Kokkos::View<pumi::ParticleData*> Partdata("particle-data",N_part);
        Kokkos::View<pumi::ParticleData*>::HostMirror h_Partdata = Kokkos::create_mirror_view(Partdata);

        h_Partdata(0) = pumi::ParticleData(1.677,0.52);
        h_Partdata(1) = pumi::ParticleData(10.0,-10.0);
        h_Partdata(2) = pumi::ParticleData(17.0,-5.3);
        h_Partdata(3) = pumi::ParticleData(12.0,6.99);
        h_Partdata(4) = pumi::ParticleData(30.0,0.0);
        Kokkos::deep_copy(Partdata, h_Partdata);
        Kokkos::parallel_for("particle-locate", N_part, KOKKOS_LAMBDA (int ipart) {
            int isub, jsub, icell, jcell, submeshID, cellID;
            double q1 = Partdata(ipart).x1;
            double q2 = Partdata(ipart).x2;
            pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
            pumi::locate_submesh_and_cell_x2(pumi_obj, q2, &jsub, &jcell);
            pumi::flatten_submeshID_and_cellID(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
            Partdata(ipart) = pumi::ParticleData(q1,q2,submeshID,cellID,true,-1);
        });
        Kokkos::deep_copy(h_Partdata,Partdata);

        int isub, icell;
        // particle-1
        isub = h_Partdata(0).submeshID;
        icell = h_Partdata(0).cellID;
        if (isub != 3) test_passed = false;
        if (icell != 6) test_passed = false;
        // particle-2
        isub = h_Partdata(1).submeshID;
        icell = h_Partdata(1).cellID;
        if (isub != 0) test_passed = false;
        if (icell != 4) test_passed = false;
        // particle-3
        isub = h_Partdata(2).submeshID;
        icell = h_Partdata(2).cellID;
        if (isub != 1) test_passed = false;
        if (icell != 20) test_passed = false;
        // particle-4
        isub = h_Partdata(3).submeshID;
        icell = h_Partdata(3).cellID;
        if (isub != 4) test_passed = false;
        if (icell != 43) test_passed = false;
        // particle-5
        isub = h_Partdata(4).submeshID;
        icell = h_Partdata(4).cellID;
        if (isub != 2) test_passed = false;
        if (icell != 33) test_passed = false;

        pumi::free_mbbl(pumi_obj);

    }
    Kokkos::finalize();
    if (test_passed) return 0;
    else return 1;
}
