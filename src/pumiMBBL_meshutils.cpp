#include "pumiMBBL_meshutils.hpp"

namespace pumi {
void check_is_pumi_working(){
    printf("Yes, pumiMBBL-GPU is working\n\n");
}

/**
 * @brief Prints the block skeleton (along with block edge tags and block vertex tags )
 *
 * \param[in] Object of the wrapper mesh structure
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
        int nnp_y = pumi_obj.host_submesh_x2[jsub]->Nel+1;
        for (int jnp=nnp_y-1; jnp>=0; jnp--){
            int Jnp = jnp+pumi_obj.host_submesh_x2[jsub]->Nel_cumulative;
            for (int isub=1; isub<=Nx; isub++){
                int nnp_x = pumi_obj.host_submesh_x1[isub]->Nel+1;
                for (int inp=0; inp<nnp_x; inp++){
                    int Inp = inp+pumi_obj.host_submesh_x1[isub]->Nel_cumulative;
                    if (pumi_obj.mesh.host_isactive[isub][jsub]){
                        int Knp = Jnp*(pumi_obj.mesh.Nel_tot_x1+1)+Inp;
                        int Ksub = (jsub-1)*Nx+isub-1;
                        int nodeID = get_global_nodeID(pumi_obj,Ksub,Knp);
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
        int nnp_y = pumi_obj.host_submesh_x2[jsub]->Nel+1;
        for (int jnp=nnp_y-1; jnp>=0; jnp--){
            int Jnp = jnp+pumi_obj.host_submesh_x2[jsub]->Nel_cumulative;
            for (int isub=1; isub<=Nx; isub++){
                int nnp_x = pumi_obj.host_submesh_x1[isub]->Nel+1;
                for (int inp=0; inp<nnp_x; inp++){
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
    int Nnp_y = pumi_obj.mesh.Nel_tot_x2;
    int Nnp_x = pumi_obj.mesh.Nel_tot_x1;
    for (int jnode=Nnp_y; jnode>=0; jnode--){
        for (int inode=0; inode<=Nnp_x; inode++){
            int subID = get_node_submeshID(pumi_obj,inode,jnode);
            if (subID+1){
                int Knp = jnode*(Nnp_x+1)+inode;
                int nodeID = get_global_nodeID(pumi_obj,subID,Knp);
                printf("%5d", nodeID);
            }
            else{
                printf("     ");
            }
        }
        printf("\n");
    }
}

} // namespace pumi
