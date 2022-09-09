#include "pumiMBBLGPU.hpp"

int main( int argc, char* argv[] )
{
    bool test_passed = true;
    Kokkos::initialize( argc, argv );
    {
        pumi::Mesh_Inputs *pumi_inputs = pumi::inputs_allocate();
        pumi_inputs->ndim = 1;
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

        pumi::Mesh_Options pumi_options;

        pumi::MBBL pumi_obj = pumi::initialize_MBBL_mesh(pumi_inputs, pumi_options);
        pumi::inputs_deallocate(pumi_inputs);

        int N_part = 5;
        Kokkos::View<pumi::ParticleData*> Partdata("particle-data",N_part);
        Kokkos::View<pumi::ParticleData*>::HostMirror h_Partdata = Kokkos::create_mirror_view(Partdata);

        h_Partdata(0) = pumi::ParticleData(1.677,0.0);
        h_Partdata(1) = pumi::ParticleData(10.0,0.0);
        h_Partdata(2) = pumi::ParticleData(17.0,0.0);
        h_Partdata(3) = pumi::ParticleData(12.0,0.0);
        h_Partdata(4) = pumi::ParticleData(30.0,0.0);
        Kokkos::deep_copy(Partdata, h_Partdata);
        Kokkos::parallel_for("particle-locate", N_part, KOKKOS_LAMBDA (int ipart) {
            int isub, icell;
            double q1 = Partdata(ipart).x1;
            pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
            printf("q=%2.2f isub=%d icell=%d\n",q1,isub,icell );
            Partdata(ipart) = pumi::ParticleData(q1,0.0,isub,icell,true,-1);
        });
        Kokkos::deep_copy(h_Partdata,Partdata);

        int isub, icell;
        // particle-1
        isub = h_Partdata(0).submeshID;
        icell = h_Partdata(0).cellID;
        if (isub != 1) test_passed = false;
        if (icell != 0) test_passed = false;
        // particle-2
        isub = h_Partdata(1).submeshID;
        icell = h_Partdata(1).cellID;
        if (isub != 1) test_passed = false;
        if (icell != 4) test_passed = false;
        // particle-3
        isub = h_Partdata(2).submeshID;
        icell = h_Partdata(2).cellID;
        if (isub != 2) test_passed = false;
        if (icell != 6) test_passed = false;
        // particle-4
        isub = h_Partdata(3).submeshID;
        icell = h_Partdata(3).cellID;
        if (isub != 2) test_passed = false;
        if (icell != 1) test_passed = false;
        // particle-5
        isub = h_Partdata(4).submeshID;
        icell = h_Partdata(4).cellID;
        if (isub != 3) test_passed = false;
        if (icell != 3) test_passed = false;

    }
    Kokkos::finalize();

    if (test_passed) return 0;
    else return 1;
}
