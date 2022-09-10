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
        for (int i=0; i<es.size(); i++){
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
                if (isubmesh==2 && jsubmesh==1)
                    pumi_inputs->isactive[isubmesh][jsubmesh] = false;
                else
                    pumi_inputs->isactive[isubmesh][jsubmesh] = true;
            }
        }

        pumi::Mesh_Options pumi_options;

        pumi::MBBL pumi_obj = pumi::initialize_MBBL_mesh(pumi_inputs, pumi_options);
        pumi::inputs_deallocate(pumi_inputs);

        int N_part = 5;
        Kokkos::View<pumi::ParticleData*> Partdata("particle-data",N_part);
        Kokkos::View<double*[2]> dx("particle-disp-data",N_part);
        Kokkos::View<pumi::ParticleData*>::HostMirror h_Partdata = Kokkos::create_mirror_view(Partdata);
        Kokkos::View<double*[2]>::HostMirror h_dx = Kokkos::create_mirror_view(dx);

        h_Partdata(0) = pumi::ParticleData(1.677,0.52);
        h_dx(0,0) = -0.18 ; h_dx(0,1) = -20.0;
        h_Partdata(1) = pumi::ParticleData(10.0,-10.0);
        h_dx(1,0) = 5.18 ; h_dx(1,1) = 13.87;
        h_Partdata(2) = pumi::ParticleData(17.0,-5.3);
        h_dx(2,0) = 14.22 ; h_dx(2,1) = 12.0;
        h_Partdata(3) = pumi::ParticleData(12.0,6.99);
        h_dx(3,0) = 0.01 ; h_dx(3,1) = -15.0;
        h_Partdata(4) = pumi::ParticleData(30.0,0.0);
        h_dx(4,0) = -14.77 ; h_dx(4,1) = 3.5;
        Kokkos::deep_copy(Partdata, h_Partdata);
        Kokkos::deep_copy(dx, h_dx);
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

        Kokkos::parallel_for("particle-push", N_part, KOKKOS_LAMBDA (const int ipart) {
            int isub, jsub, icell, jcell, kcell_x1, kcell_x2, bdry_hit, submeshID, cellID, bdry_faceID;
            bool in_domain;
            double fraction_done;
            double q1 = Partdata(ipart).x1;
            double q2 = Partdata(ipart).x2;
            double dq1 = dx(ipart,0);
            double dq2 = dx(ipart,1);
            submeshID = Partdata(ipart).submeshID;
            cellID = Partdata(ipart).cellID;

            pumi::get_directional_submeshID_and_cellID(pumi_obj,submeshID,cellID,&isub,&icell,&jsub,&jcell);

            double Wgh2_x1, Wgh2_x2, Wgh1_x1, Wgh1_x2;
            pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
            pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
            Wgh1_x1 = 1.0-Wgh2_x1;
            Wgh1_x2 = 1.0-Wgh2_x2;

            pumi::Vector3 qnew = pumi::push_particle(pumi_obj, pumi::Vector3(q1,q2,0.0), pumi::Vector3(dq1,dq2,0.0), &isub, &jsub, &icell, &jcell,
                                &in_domain, &bdry_hit, &fraction_done, &bdry_faceID);
            if (!in_domain){
                q1 = qnew[0];
                q2 = qnew[1];
                pumi::flatten_submeshID_and_cellID(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
                Partdata(ipart) = pumi::ParticleData(q1,q2,submeshID,cellID,false,bdry_faceID);
                printf("%d q=(%2.2f,%2.2f) isub=%d icell=%d bdry=%d face=%d\n",ipart,q1,q2,submeshID,cellID,bdry_hit,bdry_faceID);
            }
            else{
                q1 = qnew[0];
                q2 = qnew[1];
                pumi::calc_weights_x1(pumi_obj, q1, isub, icell, &kcell_x1, &Wgh2_x1);
                pumi::calc_weights_x2(pumi_obj, q2, jsub, jcell, &kcell_x2, &Wgh2_x2);
                Wgh1_x1 = 1.0-Wgh2_x1;
                Wgh1_x2 = 1.0-Wgh2_x2;
                pumi::flatten_submeshID_and_cellID(pumi_obj,isub,icell,jsub,jcell,&submeshID,&cellID);
                Partdata(ipart) = pumi::ParticleData(q1,q2,submeshID,cellID,true,-1,Wgh1_x1,Wgh2_x1,Wgh1_x2,Wgh2_x2);
                printf("%d q=(%2.2f,%2.2f) isub=%d icell=%d\n",ipart,q1,q2,submeshID,cellID);
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
        if (isub != 0) test_passed = false;
        if (icell != 0) test_passed = false;
        if (exitface != 0) test_passed = false;
        if (status) test_passed = false;
        // particle-2
        isub = h_Partdata(1).submeshID;
        icell = h_Partdata(1).cellID;
        status = h_Partdata(1).part_active;
        exitface = h_Partdata(1).exit_faceID;
        if (isub != 4) test_passed = false;
        if (icell != 32) test_passed = false;
        if (!status) test_passed = false;
        // particle-3
        isub = h_Partdata(2).submeshID;
        icell = h_Partdata(2).cellID;
        status = h_Partdata(2).part_active;
        exitface = h_Partdata(2).exit_faceID;
        if (isub != 2) test_passed = false;
        if (icell != 31) test_passed = false;
        if (status) test_passed = false;
        if (exitface != 33) test_passed = false;
        // particle-4
        isub = h_Partdata(3).submeshID;
        icell = h_Partdata(3).cellID;
        status = h_Partdata(3).part_active;
        exitface = h_Partdata(3).exit_faceID;
        if (isub != 1) test_passed = false;
        if (icell != 8) test_passed = false;
        if (!status) test_passed = false;
        // particle-5
        isub = h_Partdata(4).submeshID;
        icell = h_Partdata(4).cellID;
        status = h_Partdata(4).part_active;
        exitface = h_Partdata(4).exit_faceID;
        if (isub != 2) test_passed = false;
        if (icell != 33) test_passed = false;
        if (exitface != 35) test_passed = false;
        if (status) test_passed = false;

    }
    Kokkos::finalize();
    if (test_passed) return 0;
    else return 1;
}
