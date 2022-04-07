#include "pumiMBBL_meshutils.hpp"

namespace pumi {
void check_is_pumi_working(){
    printf("Yes, pumiMBBL-GPU is working\n\n");
}

/**
 * @brief Prints the block skeleton (along with block edge tags and block vertex tags )
 *
 * @param[in] Object of the wrapper mesh structure
 */
void print_mesh_skeleton(MBBL pumi_obj){

    bool on_bdry, in_domain;
    int bdry_tag, bdry_dim;
    printf("\n\nPrinting the skeleton of the mesh\n");
    printf("E --> Boundary block-edge\n");
    printf("V --> Boundary block-vertex\n\n");
    for (int jsubmesh=pumi_obj.mesh.nsubmesh_x2; jsubmesh>=1; jsubmesh--){
        int Jnp = pumi_obj.host_submesh_x2[jsubmesh]->Nel + pumi_obj.host_submesh_x2[jsubmesh]->Nel_cumulative;
        for (int isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                if (bdry_tag+1){
                    printf("%3dV",bdry_tag );
                }
                else{
                    printf("----");
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
                    printf("------------");
                }
            }
            else{
                printf("            ");
            }
            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel + pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
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
        for (int isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
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

            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel + pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
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
        for (int isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
            pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
            if (in_domain){
                if (bdry_tag+1){
                    printf("%3dE",bdry_tag);
                }
                else{
                    printf("   |");
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

            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel + pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
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
        for (int isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
            int Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
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

            if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel + pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
                pumi::where_is_node(pumi_obj, Inp, Jnp, &on_bdry, &in_domain, &bdry_tag, &bdry_dim);
                if (in_domain){
                    printf("   |");
                }
                else{
                    printf("    ");
                }
            }
        }

        if (jsubmesh==1){
            printf("\n");
            Jnp = 0;
            for (int isubmesh=1; isubmesh<=pumi_obj.mesh.nsubmesh_x1; isubmesh++){
                int Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
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
                        printf("------------");
                    }
                }
                else{
                    printf("            ");
                }
                if (isubmesh==pumi_obj.mesh.nsubmesh_x1){
                    Inp = pumi_obj.host_submesh_x1[isubmesh]->Nel + pumi_obj.host_submesh_x1[isubmesh]->Nel_cumulative;
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

void print_blockwise_nodeIDs(MBBL pumi_obj){
    printf("\nPrinting NodeIDs in each active blocks\n" );
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int Ny = pumi_obj.mesh.nsubmesh_x2;
    for (int jsub=Ny; jsub>0; jsub--){
        int Nel_y = pumi_obj.host_submesh_x2[jsub]->Nel+1;
        for (int jnp=Nel_y-1; jnp>=0; jnp--){
            int Jnp = jnp+pumi_obj.host_submesh_x2[jsub]->Nel_cumulative;
            for (int isub=1; isub<=Nx; isub++){
                int Nel_x = pumi_obj.host_submesh_x1[isub]->Nel+1;
                for (int inp=0; inp<Nel_x; inp++){
                    int Inp = inp+pumi_obj.host_submesh_x1[isub]->Nel_cumulative;
                    if (pumi_obj.mesh.host_isactive[isub][jsub]){
                        // int Knp = Jnp*(pumi_obj.mesh.Nel_tot_x1+1)+Inp;
                        // int Ksub = (jsub-1)*Nx+isub-1;
                        // int nodeID = get_global_nodeID(pumi_obj,Ksub,Knp);
                        int nodeID = get_global_nodeID_2D(pumi_obj, Inp, Jnp);
                        printf("%5d",nodeID);
                    }
                    else{
                        printf("     ");
                    }
                }
                if (isub<Nx)
                    printf(" || " );
            }
            printf("\n");
        }
        printf("\n\n");
    }
}

void print_node_submeshID(MBBL pumi_obj){
    printf("\nPrinting NodeIDs in each active blocks\n" );
    int Nx = pumi_obj.mesh.nsubmesh_x1;
    int Ny = pumi_obj.mesh.nsubmesh_x2;
    for (int jsub=Ny; jsub>0; jsub--){
        int Nel_y = pumi_obj.host_submesh_x2[jsub]->Nel+1;
        for (int jnp=Nel_y-1; jnp>=0; jnp--){
            int Jnp = jnp+pumi_obj.host_submesh_x2[jsub]->Nel_cumulative;
            for (int isub=1; isub<=Nx; isub++){
                int Nel_x = pumi_obj.host_submesh_x1[isub]->Nel+1;
                for (int inp=0; inp<Nel_x; inp++){
                    int Inp = inp+pumi_obj.host_submesh_x1[isub]->Nel_cumulative;
                    if (pumi_obj.mesh.host_isactive[isub][jsub]){
                        int subID = get_node_submeshID(pumi_obj,Inp,Jnp);
                        printf("%5d",subID);
                    }
                    else{
                        printf("     ");
                    }
                }
                if (isub<Nx)
                    printf(" || " );
            }
            printf("\n");
        }
        printf("\n\n");
    }
}

void print_fullmesh_nodeIDs(MBBL pumi_obj){
    int Nel_y = pumi_obj.mesh.Nel_tot_x2;
    int Nel_x = pumi_obj.mesh.Nel_tot_x1;
    for (int jnode=Nel_y; jnode>=0; jnode--){
        for (int inode=0; inode<=Nel_x; inode++){
            int nodeID = get_global_nodeID_2D(pumi_obj, inode, jnode);
            if (nodeID+1){
                printf("%5d", nodeID);
            }
            else{
                printf("     ");
            }
        }
        printf("\n");
    }
}

void print_2D_node_coordinates(MBBL pumi_obj){
    int Nel_y = pumi_obj.mesh.Nel_tot_x2;
    int Nel_x = pumi_obj.mesh.Nel_tot_x1;

    FILE *node_coords_file;
    char node_coords_filename[30];
    sprintf(node_coords_filename,"node_coordinates.dat");
    node_coords_file = fopen(node_coords_filename,"w");

    for (int jnode=0; jnode<=Nel_y; jnode++){
        for (int inode=0; inode<=Nel_x; inode++){
            int subID = get_node_submeshID(pumi_obj,inode,jnode);
            if (subID+1){
                int jsub = subID/pumi_obj.mesh.nsubmesh_x1 + 1;
                int isub = subID - (jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
                int inp = inode - pumi_obj.host_submesh_x1[isub]->Nel_cumulative;
                int jnp = jnode - pumi_obj.host_submesh_x2[jsub]->Nel_cumulative;
                double icoord = pumi_obj.host_submesh_x1[isub]->node_coords_host(inp);
                double jcoord = pumi_obj.host_submesh_x2[jsub]->node_coords_host(jnp);
                fprintf(node_coords_file, "%.16e %.16e\n",icoord, jcoord );
            }
        }
    }
    fclose(node_coords_file);
}

void print_2D_node_elem_connectivity(MBBL pumi_obj){
    int Nel_y = pumi_obj.mesh.Nel_tot_x2;
    int Nel_x = pumi_obj.mesh.Nel_tot_x1;

    FILE *elem_conn_file;
    char elem_conn_filename[30];
    sprintf(elem_conn_filename,"elem_node_connectivity.dat");
    elem_conn_file = fopen(elem_conn_filename,"w");
    for (int jcell=0; jcell<Nel_y; jcell++){
        for (int icell=0; icell<Nel_x; icell++){
            int subID = get_elem_submeshID(pumi_obj, icell, jcell);
            if (subID+1){
                int jsub = subID/pumi_obj.mesh.nsubmesh_x1 + 1;
                int isub = subID - (jsub-1)*pumi_obj.mesh.nsubmesh_x1 + 1;
                int inp1, inp2, inp3, inp4, global_cell;
                calc_global_cellID_and_nodeID_host(pumi_obj, isub, jsub, icell, jcell,
                                                    &global_cell, &inp1, &inp4);
                inp2 = inp1+1;
                inp3 = inp4+1;
                fprintf(elem_conn_file, "%4d %4d %4d %4d\n",inp1, inp2, inp3, inp4 );
            }
        }
    }
    fclose(elem_conn_file);

}

Vector3View compute_2D_field_gradient(MBBL pumi_obj, DoubleView phi){
    int nnp_total = phi.extent(0);
    Vector3View phi_grad = Vector3View("phi_grad",nnp_total);

    // Loop #1 -- Loop over all block interior nodes
    int nsubmesh_tot = get_total_submesh_blocks(pumi_obj);
    for(int submeshID=0; submeshID<nsubmesh_tot; submeshID++){
        int isub, jsub;
        get_directional_submeshID_host(pumi_obj, submeshID, &isub, &jsub);
        if (is_block_active_host(pumi_obj,isub,jsub)){
            int blk_interior_nnp = get_num_interior_nodes_on_block(pumi_obj, isub, jsub);

            Kokkos::parallel_for("grad_blk_interior", blk_interior_nnp, KOKKOS_LAMBDA (const int inode){
                int inp, jnp;
                get_directional_interior_nodeIDs(pumi_obj, isub, jsub, inode, &inp, &jnp);
                int inode_curr, inode_north, inode_south, inode_east, inode_west;
                double dx1_max, dx1_min, dx2_max, dx2_min;
                inode_curr = calc_global_nodeID(pumi_obj, isub, jsub, inp, jnp);
                inode_east = inode_curr+1;
                dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inp);
                inode_west = inode_curr-1;
                dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inp-1);
                inode_north = calc_global_nodeID(pumi_obj, isub, jsub, inp, jnp+1);
                dx2_max = get_x2_elem_size_in_submesh(pumi_obj, jsub, jnp);
                inode_south = calc_global_nodeID(pumi_obj, isub, jsub, inp, jnp-1);
                dx2_min = get_x2_elem_size_in_submesh(pumi_obj, jsub, jnp-1);

                phi_grad(inode_curr)[0] = - (phi(inode_east) - phi(inode_west))/(dx1_min+dx1_max);
                phi_grad(inode_curr)[1] = - (phi(inode_north) - phi(inode_south))/(dx2_min+dx2_max);
            });
        }
    }

    // Loop #2 -- Loop over all block-edge-interior nodes
    int nedges_tot = get_total_mesh_block_edges(pumi_obj);
    for (int iEdge=0; iEdge<nedges_tot; iEdge++){
        int submeshID = get_block_edge_submeshID_host(pumi_obj, iEdge);

        if (submeshID+1){
            int isub, jsub;
            get_directional_submeshID_host(pumi_obj, submeshID, &isub, &jsub);
            int edge_interior_nnp = get_num_interior_nodes_on_edge(pumi_obj, iEdge);
            Vector3 edge_nrml = get_edge_normal_host(pumi_obj,iEdge);
            if (is_horizontal_edge(pumi_obj,iEdge)){
                if (edge_nrml[1]>0.0){
                    Kokkos::parallel_for("grad_horizontal_edge_interior_north_normal",edge_interior_nnp,KOKKOS_LAMBDA (const int inode){
                        int inode_curr = calc_global_nodeID_on_horizontal_edge(pumi_obj,iEdge,inode);
                        int inode_east = inode_curr+1;
                        int inode_west = inode_curr-1;
                        double dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inode);
                        double dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inode+1);
                        phi_grad(inode_curr)[0] = - (phi(inode_east) - phi(inode_west))/(dx1_min+dx1_max);

                        int inode_south_first = calc_first_south_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                        int inode_south_second = calc_second_south_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                        double r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                        r = 1.0/r;
                        int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub);
                        double dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,nel_blk-1);
                        phi_grad(inode_curr)[1] = -( -(1.0+r)*(1.0+r)*phi(inode_south_first) + phi(inode_south_second) + r*(r+2.0)*phi(inode_curr) )/
                                                    (r*(r+1.0)*dx2);
                    });
                }
                else if (edge_nrml[1]<0.0){
                    Kokkos::parallel_for("grad_horizontal_edge_interior_south_normal",edge_interior_nnp,KOKKOS_LAMBDA (const int inode){
                        int inode_curr = calc_global_nodeID_on_horizontal_edge(pumi_obj,iEdge,inode);
                        int inode_east = inode_curr+1;
                        int inode_west = inode_curr-1;
                        double dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inode);
                        double dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inode+1);
                        phi_grad(inode_curr)[0] = - (phi(inode_east) - phi(inode_west))/(dx1_min+dx1_max);

                        int inode_north_first = calc_first_north_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                        int inode_north_second = calc_second_north_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                        double r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                        double dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,0);
                        phi_grad(inode_curr)[1] = -( (1.0+r)*(1.0+r)*phi(inode_north_first) - phi(inode_north_second) - r*(r+2.0)*phi(inode_curr) ) /
                                                     (r*(r+1.0)*dx2);
                    });
                }
                else{
                    Kokkos::parallel_for("grad_horizontal_edge_interior_nonbdry",edge_interior_nnp,KOKKOS_LAMBDA (const int inode){
                        int inode_curr = calc_global_nodeID_on_horizontal_edge(pumi_obj,iEdge,inode);
                        int inode_east = inode_curr+1;
                        int inode_west = inode_curr-1;
                        double dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inode);
                        double dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inode+1);
                        phi_grad(inode_curr)[0] = - (phi(inode_east) - phi(inode_west))/(dx1_min+dx1_max);

                        int inode_north_first = calc_first_north_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                        int inode_south_first = calc_first_south_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                        int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub);
                        double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub,nel_blk-1);
                        double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub+1,0);
                        phi_grad(inode_curr)[1] = -(phi(inode_north_first) - phi(inode_south_first))/(dx2_min+dx2_max);
                    });
                }
            }
            else{
                if (edge_nrml[0] > 0.0){
                    Kokkos::parallel_for("grad_vertical_edge_interior_east_normal",edge_interior_nnp,KOKKOS_LAMBDA (const int inode){
                        int inode_curr = calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,inode);
                        int inode_north = calc_first_north_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
                        int inode_south = calc_first_south_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
                        double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode);
                        double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode+1);
                        phi_grad(inode_curr)[1] = -(phi(inode_north) - phi(inode_south))/(dx2_min+dx2_max);

                        int inode_west_first = inode_curr-1;
                        int inode_west_second = inode_curr-2;
                        double r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                        r = 1.0/r;
                        int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub);
                        double dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,nel_blk-1);
                        phi_grad(inode_curr)[0] = -( -(1.0+r)*(1.0+r)*phi(inode_west_first) + phi(inode_west_second) + r*(r+2.0)*phi(inode_curr) ) /
                                                     (r*(r+1.0)*dx1);
                    });
                }
                else if (edge_nrml[0] < 0.0){
                    Kokkos::parallel_for("grad_vertical_edge_interior_west_normal",edge_interior_nnp,KOKKOS_LAMBDA (const int inode){
                        int inode_curr = calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,inode);
                        int inode_north = calc_first_north_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
                        int inode_south = calc_first_south_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
                        double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode);
                        double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode+1);
                        phi_grad(inode_curr)[1] = -(phi(inode_north) - phi(inode_south))/(dx2_min+dx2_max);

                        int inode_east_first = inode_curr+1;
                        int inode_east_second = inode_curr+2;
                        double r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                        double dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,0);
                        phi_grad(inode_curr)[0] = -( (1.0+r)*(1.0+r)*phi(inode_east_first) - phi(inode_east_second)- r*(r+2.0)*phi(inode_curr) ) /
                                                     (r*(r+1.0)*dx1);
                    });
                }
                else{
                    Kokkos::parallel_for("grad_vertical_edge_interior_nonbdry",edge_interior_nnp,KOKKOS_LAMBDA (const int inode){
                        int inode_curr = calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,inode);
                        int inode_north = calc_first_north_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
                        int inode_south = calc_first_south_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
                        double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode);
                        double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode+1);
                        phi_grad(inode_curr)[1] = -(phi(inode_north) - phi(inode_south))/(dx2_min+dx2_max);

                        int inode_east_first = inode_curr+1;
                        int inode_west_first = inode_curr-1;
                        int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub);
                        double dx1_min = get_x1_elem_size_in_submesh(pumi_obj,isub,nel_blk-1);
                        double dx1_max = get_x1_elem_size_in_submesh(pumi_obj,isub+1,0);
                        phi_grad(inode_curr)[0] = - (phi(inode_east_first)-phi(inode_west_first))/(dx1_max+dx1_min);
                    });
                }
            }
        }
    }

    // Loop #3 --  Loop over all block-vertex nodes
    int nverts_tot = get_total_mesh_block_verts(pumi_obj);
    Kokkos::parallel_for("grad_block_verts",nverts_tot,KOKKOS_LAMBDA (const int inode){
        int inode_curr = get_block_vert_nodeID(pumi_obj,inode);
        Vector3 vert_nrml = get_vert_normal(pumi_obj,inode);
        int submeshID = get_block_vert_submeshID(pumi_obj, inode);

        if (submeshID+1){
            int isub, jsub;
            get_directional_submeshID(pumi_obj, submeshID, &isub, &jsub);
            if (vert_nrml[0]==1.0){
                int inode_west_first = calc_first_west_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_west_second = calc_second_west_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx1;
                r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub);
                dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,nel_blk-1);
                r = 1.0/r;
                phi_grad(inode_curr)[0] = -( -(1.0+r)*(1.0+r)*phi(inode_west_first) + phi(inode_west_second) + r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx1);
            }
            else if (vert_nrml[0]==-1.0){
                int inode_east_first = calc_first_east_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_east_second = calc_second_east_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx1;
                r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,0);
                phi_grad(inode_curr)[0] = -( (1.0+r)*(1.0+r)*phi(inode_east_first) - phi(inode_east_second)- r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx1);
            }
            else{
                int inode_east_first = calc_first_east_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_west_first = calc_first_west_global_nodeID_to_vertex(pumi_obj,inode);
                int isub_west = get_x1_submeshID_west_to_vertex(pumi_obj,inode);
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub_west);
                double dx1_min = get_x1_elem_size_in_submesh(pumi_obj,isub_west,nel_blk-1);
                double dx1_max = get_x1_elem_size_in_submesh(pumi_obj,isub_west+1,0);
                phi_grad(inode_curr)[0] = - (phi(inode_east_first)-phi(inode_west_first))/(dx1_max+dx1_min);
            }

            if (vert_nrml[1]==1.0){
                int inode_south_first = calc_first_south_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_south_second = calc_second_south_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx2;
                r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub);
                dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,nel_blk-1);
                r = 1.0/r;
                phi_grad(inode_curr)[1] = -( -(1.0+r)*(1.0+r)*phi(inode_south_first) + phi(inode_south_second) + r*(r+2.0)*phi(inode_curr) )/
                                            (r*(r+1.0)*dx2);
            }
            else if (vert_nrml[1]==-1.0){
                int inode_north_first = calc_first_north_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_north_second = calc_second_north_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx2;
                r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,0);
                phi_grad(inode_curr)[1] = -( (1.0+r)*(1.0+r)*phi(inode_north_first) - phi(inode_north_second) - r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx2);
            }
            else{
                int inode_north_first = calc_first_north_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_south_first = calc_first_south_global_nodeID_to_vertex(pumi_obj,inode);
                int jsub_south = get_x2_submeshID_south_to_vertex(pumi_obj,inode);
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub_south);
                double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub_south,nel_blk-1);
                double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub_south+1,0);
                phi_grad(inode_curr)[1] = -(phi(inode_north_first) - phi(inode_south_first))/(dx2_min+dx2_max);
            }
        }
    });

    return phi_grad;
}

Vector3View compute_2D_field_gradient_v2(MBBL pumi_obj, DoubleView phi){
    int nnp_total = phi.extent(0);
    Vector3View phi_grad = Vector3View("phi_grad",nnp_total);

    // Loop #1 -- Loop over all block interior nodes
    int tot_blk_nodes = get_num_block_interior_nodes(pumi_obj);
    Kokkos::parallel_for("grad_blk_interior",tot_blk_nodes, KOKKOS_LAMBDA (const int ignode){
        int isub, jsub, inp, jnp;
        get_submeshIDs_and_localnodeIDs_of_block_interior_nodes(pumi_obj, ignode, &isub, &jsub, &inp, &jnp);
        int inode_curr, inode_north, inode_south, inode_east, inode_west;
        double dx1_max, dx1_min, dx2_max, dx2_min;
        inode_curr = calc_global_nodeID(pumi_obj, isub, jsub, inp, jnp);
        inode_east = inode_curr+1;
        dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inp);
        inode_west = inode_curr-1;
        dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inp-1);
        inode_north = calc_global_nodeID(pumi_obj, isub, jsub, inp, jnp+1);
        dx2_max = get_x2_elem_size_in_submesh(pumi_obj, jsub, jnp);
        inode_south = calc_global_nodeID(pumi_obj, isub, jsub, inp, jnp-1);
        dx2_min = get_x2_elem_size_in_submesh(pumi_obj, jsub, jnp-1);

        phi_grad(inode_curr)[0] = - (phi(inode_east) - phi(inode_west))/(dx1_min+dx1_max);
        phi_grad(inode_curr)[1] = - (phi(inode_north) - phi(inode_south))/(dx2_min+dx2_max);
    });

    // Loop #2 -- Loop over all block-edge-interior nodes
    int tot_edg_nodes = get_num_block_edge_interior_nodes(pumi_obj);
    Kokkos::parallel_for("grad_edge_interior",tot_edg_nodes,KOKKOS_LAMBDA (const int ignode){
        int inode, isub, jsub, iEdge;
        get_edgeIDs_submeshIDs_and_localnodeIDs_of_block_edge_interior_nodes
                                    (pumi_obj, ignode, &iEdge, &isub, &jsub, &inode);
        if (is_horizontal_edge(pumi_obj,iEdge)){
            Vector3 edge_nrml = get_edge_normal(pumi_obj, iEdge);
            int inode_curr = calc_global_nodeID_on_horizontal_edge(pumi_obj,iEdge,inode);
            int inode_east = inode_curr+1;
            int inode_west = inode_curr-1;
            double dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inode);
            double dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inode+1);
            phi_grad(inode_curr)[0] = - (phi(inode_east) - phi(inode_west))/(dx1_min+dx1_max);
            if (edge_nrml[1]==1.0){
                int inode_south_first = calc_first_south_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                int inode_south_second = calc_second_south_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                // double r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                // r = 1.0/r;
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub);
                double r = get_x2_gradingratio_in_submesh(pumi_obj,jsub,nel_blk-1);
                r = 1.0/r;
                double dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,nel_blk-1);
                phi_grad(inode_curr)[1] = -( -(1.0+r)*(1.0+r)*phi(inode_south_first) + phi(inode_south_second) + r*(r+2.0)*phi(inode_curr) )/
                                            (r*(r+1.0)*dx2);
            }
            else if (edge_nrml[1]==-1.0){
                int inode_north_first = calc_first_north_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                int inode_north_second = calc_second_north_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                // double r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                double r = get_x2_gradingratio_in_submesh(pumi_obj,jsub,1);
                double dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,0);
                phi_grad(inode_curr)[1] = -( (1.0+r)*(1.0+r)*phi(inode_north_first) - phi(inode_north_second) - r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx2);
            }
            else{
                int inode_north_first = calc_first_north_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                int inode_south_first = calc_first_south_global_nodeID_to_horizontal_edge(pumi_obj,iEdge,inode);
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub);
                double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub,nel_blk-1);
                double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub+1,0);
                phi_grad(inode_curr)[1] = -(phi(inode_north_first) - phi(inode_south_first))/(dx2_min+dx2_max);
            }
        }
        else{
            Vector3 edge_nrml = get_edge_normal(pumi_obj, iEdge);
            int inode_curr = calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,inode);
            int inode_north = calc_first_north_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
            int inode_south = calc_first_south_global_nodeID_to_vertical_edge(pumi_obj,iEdge,inode);
            double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode);
            double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode+1);
            phi_grad(inode_curr)[1] = -(phi(inode_north) - phi(inode_south))/(dx2_min+dx2_max);

            if (edge_nrml[0]==1.0){
                int inode_west_first = inode_curr-1;
                int inode_west_second = inode_curr-2;
                // double r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                // r = 1.0/r;
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub);
                double r = get_x1_gradingratio_in_submesh(pumi_obj,isub,nel_blk-1);
                r = 1.0/r;
                double dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,nel_blk-1);
                phi_grad(inode_curr)[0] = -( -(1.0+r)*(1.0+r)*phi(inode_west_first) + phi(inode_west_second) + r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx1);
            }
            else if (edge_nrml[0]==-1.0){
                int inode_east_first = inode_curr+1;
                int inode_east_second = inode_curr+2;
                // double r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                double r = get_x1_gradingratio_in_submesh(pumi_obj,isub,1);
                double dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,0);
                phi_grad(inode_curr)[0] = -( (1.0+r)*(1.0+r)*phi(inode_east_first) - phi(inode_east_second)- r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx1);
            }
            else{
                int inode_east_first = inode_curr+1;
                int inode_west_first = inode_curr-1;
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub);
                double dx1_min = get_x1_elem_size_in_submesh(pumi_obj,isub,nel_blk-1);
                double dx1_max = get_x1_elem_size_in_submesh(pumi_obj,isub+1,0);
                phi_grad(inode_curr)[0] = - (phi(inode_east_first)-phi(inode_west_first))/(dx1_max+dx1_min);
            }
        }
    });

    // Loop #3 --  Loop over all block-vertex nodes
    int nverts_tot = get_total_mesh_block_verts(pumi_obj);
    Kokkos::parallel_for("grad_block_verts",nverts_tot,KOKKOS_LAMBDA (const int inode){
        int inode_curr = get_block_vert_nodeID(pumi_obj,inode);
        Vector3 vert_nrml = get_vert_normal(pumi_obj,inode);
        int submeshID = get_block_vert_submeshID(pumi_obj, inode);

        if (submeshID+1){
            int isub, jsub;
            get_directional_submeshID(pumi_obj, submeshID, &isub, &jsub);
            if (vert_nrml[0]==1.0){
                int inode_west_first = calc_first_west_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_west_second = calc_second_west_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx1;
                // r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub);
                r = get_x1_gradingratio_in_submesh(pumi_obj,isub,nel_blk-1);
                dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,nel_blk-1);
                r = 1.0/r;
                phi_grad(inode_curr)[0] = -( -(1.0+r)*(1.0+r)*phi(inode_west_first) + phi(inode_west_second) + r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx1);
            }
            else if (vert_nrml[0]==-1.0){
                int inode_east_first = calc_first_east_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_east_second = calc_second_east_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx1;
                // r = get_x1_gradingratio_in_submesh(pumi_obj,isub);
                r = get_x1_gradingratio_in_submesh(pumi_obj,isub,1);
                dx1 = get_x1_elem_size_in_submesh(pumi_obj,isub,0);
                phi_grad(inode_curr)[0] = -( (1.0+r)*(1.0+r)*phi(inode_east_first) - phi(inode_east_second)- r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx1);
            }
            else{
                int inode_east_first = calc_first_east_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_west_first = calc_first_west_global_nodeID_to_vertex(pumi_obj,inode);
                int isub_west = get_x1_submeshID_west_to_vertex(pumi_obj,inode);
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub_west);
                double dx1_min = get_x1_elem_size_in_submesh(pumi_obj,isub_west,nel_blk-1);
                double dx1_max = get_x1_elem_size_in_submesh(pumi_obj,isub_west+1,0);
                phi_grad(inode_curr)[0] = - (phi(inode_east_first)-phi(inode_west_first))/(dx1_max+dx1_min);
            }

            if (vert_nrml[1]==1.0){
                int inode_south_first = calc_first_south_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_south_second = calc_second_south_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx2;
                // r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub);
                r = get_x2_gradingratio_in_submesh(pumi_obj,jsub,nel_blk-1);
                dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,nel_blk-1);
                r = 1.0/r;
                phi_grad(inode_curr)[1] = -( -(1.0+r)*(1.0+r)*phi(inode_south_first) + phi(inode_south_second) + r*(r+2.0)*phi(inode_curr) )/
                                            (r*(r+1.0)*dx2);
            }
            else if (vert_nrml[1]==-1.0){
                int inode_north_first = calc_first_north_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_north_second = calc_second_north_global_nodeID_to_vertex(pumi_obj,inode);
                double r, dx2;
                // r = get_x2_gradingratio_in_submesh(pumi_obj,jsub);
                r = get_x2_gradingratio_in_submesh(pumi_obj,jsub,1);
                dx2 = get_x2_elem_size_in_submesh(pumi_obj,jsub,0);
                phi_grad(inode_curr)[1] = -( (1.0+r)*(1.0+r)*phi(inode_north_first) - phi(inode_north_second) - r*(r+2.0)*phi(inode_curr) ) /
                                             (r*(r+1.0)*dx2);
            }
            else{
                int inode_north_first = calc_first_north_global_nodeID_to_vertex(pumi_obj,inode);
                int inode_south_first = calc_first_south_global_nodeID_to_vertex(pumi_obj,inode);
                int jsub_south = get_x2_submeshID_south_to_vertex(pumi_obj,inode);
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub_south);
                double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub_south,nel_blk-1);
                double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub_south+1,0);
                phi_grad(inode_curr)[1] = -(phi(inode_north_first) - phi(inode_south_first))/(dx2_min+dx2_max);
            }
        }
    });
    Kokkos::fence();

    return phi_grad;
}

Vector3View compute_2D_field_gradient_fulluniform(MBBL pumi_obj, DoubleView phi){
    int num_interior_nodes = (pumi_obj.mesh.Nel_tot_x1-1)*(pumi_obj.mesh.Nel_tot_x2-1);
    int num_boundary_nodes = (pumi_obj.mesh.Nel_tot_x1+1)*(pumi_obj.mesh.Nel_tot_x2+1)-num_interior_nodes;
    int num_x1_nodes = pumi_obj.mesh.Nel_tot_x1+1;
    int num_x2_nodes = pumi_obj.mesh.Nel_tot_x2+1;
    double dx1 = pumi_obj.host_submesh_x1[1]->t0;
    double dx2 = pumi_obj.host_submesh_x2[1]->t0;

    int nnp_total = num_x1_nodes*num_x2_nodes;
    Vector3View phi_grad = Vector3View("phi_grad",nnp_total);

    Kokkos::parallel_for(
        "UniformMesh2D::gradientNodalScalarField::interior_nodes",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {num_x1_nodes-1, num_x2_nodes-1}),
        KOKKOS_LAMBDA (const int inode, const int jnode)
    {
        // Need to find the global index corresponding to this interior node
        int gnode = inode + jnode * num_x1_nodes;
        phi_grad(gnode)[0] = -( phi(gnode+1) - phi(gnode-1) ) / dx1 / 2.0;
        // Under our dictionary ordering,
        // adding num_x1_nodes gives the index of the neighboring node
        // in the positive x2 direction.
        phi_grad(gnode)[1] = -( phi(gnode+num_x1_nodes) - phi(gnode-num_x1_nodes) ) / dx2 / 2.0;
    });

    Kokkos::parallel_for(
        "UniformMesh2D::gradientNodalScalarField::boundary_nodes",
        num_boundary_nodes,
        KOKKOS_LAMBDA (const int boundary_node)
    {
        short north_flag = boundary_node / (num_x1_nodes + 2*num_x2_nodes - 4);
        short south_flag = (boundary_node / num_x1_nodes) == 0;

        int inode = 0;
        int jnode = 0;
        if (south_flag) {
            inode = boundary_node;
        } else if (north_flag) {
            inode = boundary_node + num_interior_nodes - num_x1_nodes*(num_x2_nodes - 1);
            jnode = num_x2_nodes - 1;
        } else {
            inode = ((boundary_node - num_x1_nodes) % 2) * (num_x1_nodes - 1);
            jnode =  (boundary_node - num_x1_nodes) / 2 + 1;
        }

        // Need to find the global index corresponding to this boundary node
        int gnode = inode + jnode * num_x1_nodes;

        if (inode == 0) {
            phi_grad(gnode)[0] = -( -0.5*phi(gnode+2) + 2.0*phi(gnode+1) - 1.5*phi(gnode  ) ) / dx1;
        }
        else if (inode == num_x1_nodes-1) {
            phi_grad(gnode)[0] = -(  1.5*phi(gnode)   - 2.0*phi(gnode-1) + 0.5*phi(gnode-2) ) / dx1;
        }
        else {
            phi_grad(gnode)[0] = -( phi(gnode+1) - phi(gnode-1) ) / dx1 / 2.0;
        }

        if (south_flag) {
            phi_grad(gnode)[1] = -( -0.5*phi(gnode+2*num_x1_nodes) + 2.0*phi(gnode+num_x1_nodes) - 1.5*phi(gnode) ) / dx2;
        }
        else if (north_flag) {
            phi_grad(gnode)[1] = -(  1.5*phi(gnode) - 2.0*phi(gnode-num_x1_nodes) + 0.5*phi(gnode-2*num_x1_nodes) ) / dx2;
        }
        else {
            phi_grad(gnode)[1] = -( phi(gnode+num_x1_nodes) - phi(gnode-num_x1_nodes) ) / dx2 / 2.0;
        }
    });
    Kokkos::fence();

    return phi_grad;
}

DoubleView compute_2D_field_density(MBBL pumi_obj, DoubleView Q){
    int nnp_total = Q.extent(0);
    DoubleView rho = DoubleView("Q-density",nnp_total);
    // Loop #1 -- Loop over all block interior nodes
    int tot_blk_nodes = get_num_block_interior_nodes(pumi_obj);
    Kokkos::parallel_for("grad_blk_interior",tot_blk_nodes, KOKKOS_LAMBDA (const int ignode){
        int isub, jsub, inp, jnp;
        get_submeshIDs_and_localnodeIDs_of_block_interior_nodes(pumi_obj, ignode, &isub, &jsub, &inp, &jnp);
        int inode_curr;
        double dx1_max, dx1_min, dx2_max, dx2_min, cov;
        inode_curr = calc_global_nodeID(pumi_obj, isub, jsub, inp, jnp);
        dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inp);
        dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inp-1);
        dx2_max = get_x2_elem_size_in_submesh(pumi_obj, jsub, jnp);
        dx2_min = get_x2_elem_size_in_submesh(pumi_obj, jsub, jnp-1);
        cov = 0.25*(dx1_min+dx1_max)*(dx2_min+dx2_max);
        rho(inode_curr) = Q(inode_curr)/cov;
    });

    // Loop #2 -- Loop over all block-edge-interior nodes
    int tot_edg_nodes = get_num_block_edge_interior_nodes(pumi_obj);
    Kokkos::parallel_for("grad_edge_interior",tot_edg_nodes,KOKKOS_LAMBDA (const int ignode){
        int inode, isub, jsub, iEdge;
        get_edgeIDs_submeshIDs_and_localnodeIDs_of_block_edge_interior_nodes
                                    (pumi_obj, ignode, &iEdge, &isub, &jsub, &inode);
        if (is_horizontal_edge(pumi_obj,iEdge)){
            int inode_curr = calc_global_nodeID_on_horizontal_edge(pumi_obj,iEdge,inode);
            int jsub_north = get_x2_submeshID_north_to_horizontal_edge(pumi_obj,iEdge);
            int jsub_south = jsub_north-1;
            double cov = 0.0;
            double dx2_max, dx2_min;
            double dx1_min = get_x1_elem_size_in_submesh(pumi_obj, isub, inode);
            double dx1_max = get_x1_elem_size_in_submesh(pumi_obj, isub, inode+1);

            if (pumi_obj.mesh.isactive(isub,jsub_north)){
                dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub_north,0);
                cov += 0.25*(dx2_max*(dx1_min+dx1_max));
            }
            if (pumi_obj.mesh.isactive(isub,jsub_south)){
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub_south);
                dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub_south,nel_blk-1);
                cov += 0.25*(dx2_min*(dx1_min+dx1_max));
            }

            rho(inode_curr) = Q(inode_curr)/cov;
        }
        else{
            int inode_curr = calc_global_nodeID_on_vertical_edge(pumi_obj,iEdge,inode);
            int isub_east = get_x1_submeshID_east_to_vertical_edge(pumi_obj,iEdge);
            int isub_west = isub_east-1;
            double cov = 0.0;
            double dx1_max, dx1_min;
            double dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode);
            double dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub,inode+1);

            if (pumi_obj.mesh.isactive(isub_east,jsub)){
                dx1_max = get_x1_elem_size_in_submesh(pumi_obj,isub_east,0);
                cov += 0.25*(dx1_max*(dx2_min+dx2_max));
            }
            if (pumi_obj.mesh.isactive(isub_west,jsub)){
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub_west);
                dx1_min = get_x1_elem_size_in_submesh(pumi_obj,isub_west,nel_blk-1);
                cov += 0.25*(dx1_min*(dx2_min+dx2_max));
            }

            rho(inode_curr) = Q(inode_curr)/cov;
        }
    });

    // Loop #3 --  Loop over all block-vertex nodes
    int nverts_tot = get_total_mesh_block_verts(pumi_obj);
    Kokkos::parallel_for("grad_block_verts",nverts_tot,KOKKOS_LAMBDA (const int inode){
        int inode_curr = get_block_vert_nodeID(pumi_obj,inode);
        int submeshID = get_block_vert_submeshID(pumi_obj, inode);
        if (submeshID+1){
            int isub_east = get_x1_submeshID_east_to_vertex(pumi_obj, inode);
            int isub_west = get_x1_submeshID_west_to_vertex(pumi_obj, inode);
            int jsub_north = get_x2_submeshID_north_to_vertex(pumi_obj, inode);
            int jsub_south = get_x2_submeshID_south_to_vertex(pumi_obj, inode);
            double dx1_min=0.0, dx1_max=0.0, dx2_min=0.0, dx2_max=0.0, cov=0.0;

            if (pumi_obj.mesh.isactive(isub_west,jsub_south)){
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub_west);
                dx1_min = get_x1_elem_size_in_submesh(pumi_obj,isub_west,nel_blk-1);
                nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub_south);
                dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub_south,nel_blk-1);
                cov += 0.25*dx1_min*dx2_min;
            }

            if (pumi_obj.mesh.isactive(isub_west,jsub_north)){
                int nel_blk = get_num_x1_elems_in_submesh(pumi_obj,isub_west);
                dx1_min = get_x1_elem_size_in_submesh(pumi_obj,isub_west,nel_blk-1);
                dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub_north,0);
                cov += 0.25*dx1_min*dx2_max;
            }

            if (pumi_obj.mesh.isactive(isub_east,jsub_south)){
                dx1_max = get_x1_elem_size_in_submesh(pumi_obj,isub_east,0);
                int nel_blk = get_num_x2_elems_in_submesh(pumi_obj,jsub_south);
                dx2_min = get_x2_elem_size_in_submesh(pumi_obj,jsub_south,nel_blk-1);
                cov += 0.25*dx1_max*dx2_min;
            }

            if (pumi_obj.mesh.isactive(isub_east,jsub_north)){
                dx1_max = get_x1_elem_size_in_submesh(pumi_obj,isub_east,0);
                dx2_max = get_x2_elem_size_in_submesh(pumi_obj,jsub_north,0);
                cov += 0.25*dx1_max*dx2_max;
            }

            rho(inode_curr) = Q(inode_curr)/cov;
        }
    });

    return rho;
}
} // namespace pumi
