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
        double x1_min = pumi::get_global_x1_min_coord(pumi_obj);
        double x1_max = pumi::get_global_x1_max_coord(pumi_obj);

        int N_part = 5;
        Kokkos::View<pumi::ParticleData*> Partdata("particle-data",N_part);
        Kokkos::View<pumi::ParticleData*>::HostMirror h_Partdata = Kokkos::create_mirror_view(Partdata);

        h_Partdata(0) = pumi::ParticleData(1.677,28.7);
        h_Partdata(1) = pumi::ParticleData(10.0,-8.49);
        h_Partdata(2) = pumi::ParticleData(17.0,2.0);
        h_Partdata(3) = pumi::ParticleData(12.0,-0.501);
        h_Partdata(4) = pumi::ParticleData(30.0,3.5001);
        Kokkos::deep_copy(Partdata, h_Partdata);
        Kokkos::parallel_for("particle-locate", N_part, KOKKOS_LAMBDA (int ipart) {
            int isub, icell;
            double q1 = Partdata(ipart).x1;
            double dq1 = Partdata(ipart).x2;
            pumi::locate_submesh_and_cell_x1(pumi_obj, q1, &isub, &icell);
            Partdata(ipart) = pumi::ParticleData(q1,dq1,isub,icell,true,-1);
        });
        Kokkos::parallel_for("particle-push-test-0", N_part, KOKKOS_LAMBDA (const int ipart) {
            int isub, icell, kcell_x1;
            double q1 = Partdata(ipart).x1;
            double dq1 = Partdata(ipart).x2;
            isub = Partdata(ipart).submeshID;
            icell = Partdata(ipart).cellID;

            double Wgh2_x1, Wgh1_x1;
            pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
            Wgh1_x1 = 1.0-Wgh2_x1;

            q1 += dq1;

            if (q1 < x1_min || q1 > x1_max){
                int exit_faceID = (q1 < x1_min) + 2*(q1>x1_max) - 1;
                isub = -1;
                icell = -1;
                Partdata(ipart) = pumi::ParticleData(q1,0.0,isub,icell,false,exit_faceID);
                // printf("q1=%2.3f isub=%d, icell=%d, exitface=%d status=%d\n",q1,isub,icell,exit_faceID,false );
            }
            else {
                pumi::update_submesh_and_cell_x1(pumi_obj, q1, isub, icell, &isub, &icell);
                pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                Wgh1_x1 = 1.0-Wgh2_x1;
                Partdata(ipart) = pumi::ParticleData(q1,0.0,isub,icell,true,-1,Wgh1_x1,Wgh2_x1,0.0,0.0);
                // printf("q1=%2.3f isub=%d, icell=%d, exitface=%d status=%d\n",q1,isub,icell,-1,true );
            }

        });
        Kokkos::deep_copy(h_Partdata,Partdata);

        int isub, icell, exitface;
        bool status;
        // particle-1
        isub = h_Partdata(0).submeshID;
        icell = h_Partdata(0).cellID;
        status = h_Partdata(0).part_active;
        exitface = h_Partdata(0).exit_faceID;
        if (isub != 3) test_passed = false;
        if (icell != 3) test_passed = false;
        if (!status) test_passed = false;
        // particle-2
        isub = h_Partdata(1).submeshID;
        icell = h_Partdata(1).cellID;
        status = h_Partdata(1).part_active;
        exitface = h_Partdata(1).exit_faceID;
        if (isub != 1) test_passed = false;
        if (icell != 0) test_passed = false;
        if (!status) test_passed = false;
        // particle-3
        isub = h_Partdata(2).submeshID;
        icell = h_Partdata(2).cellID;
        status = h_Partdata(2).part_active;
        exitface = h_Partdata(2).exit_faceID;
        if (isub != 3) test_passed = false;
        if (icell != 0) test_passed = false;
        if (!status) test_passed = false;
        // particle-4
        isub = h_Partdata(3).submeshID;
        icell = h_Partdata(3).cellID;
        status = h_Partdata(3).part_active;
        exitface = h_Partdata(3).exit_faceID;
        if (isub != 1) test_passed = false;
        if (icell != 5) test_passed = false;
        if (!status) test_passed = false;
        // particle-5
        isub = h_Partdata(4).submeshID;
        icell = h_Partdata(4).cellID;
        status = h_Partdata(4).part_active;
        exitface = h_Partdata(4).exit_faceID;
        if (isub != -1) test_passed = false;
        if (icell != -1) test_passed = false;
        if (exitface != 1) test_passed = false;
        if (status) test_passed = false;

    }
    Kokkos::finalize();

    if (test_passed) return 0;
    else return 1;
}
