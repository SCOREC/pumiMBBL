#include "pumiMBBLGPU.hpp"

namespace pumi {

    void free_mbbl(MBBL pumi_obj){

        Mesh mesh = pumi_obj.mesh;
        SubmeshHostViewPtr host_submesh_x1 = pumi_obj.host_submesh_x1;
        SubmeshHostViewPtr host_submesh_x2 = pumi_obj.host_submesh_x2;

        if (host_submesh_x1){
            for (int i=0; i<mesh.nsubmesh_x1+2; i++){
                if (host_submesh_x1[i]->host_BL_coords) {
                    delete[] host_submesh_x1[i]->host_BL_coords;
                }

                delete host_submesh_x1[i];

            }
            delete[] host_submesh_x1;
        }
        if (host_submesh_x2){
            for (int i=0; i<mesh.nsubmesh_x2+2; i++){
                if (host_submesh_x2[i]->host_BL_coords){
                    delete[] host_submesh_x2[i]->host_BL_coords;
                }
                delete host_submesh_x2[i];
            }
            delete[] host_submesh_x2;
        }

        if (mesh.host_isactive){
            for (int i=0; i<mesh.nsubmesh_x1+2; i++){
                delete[] mesh.host_isactive[i];
            }
            delete[] mesh.host_isactive;
        }

        free_mesh(mesh);

    }

    void free_mesh(Mesh mesh){

        free_offsets(mesh.offsets, mesh.nsubmesh_x1);
        free_bdry(mesh.bdry);
        free_blkif(mesh.blkif);
        free_bst(mesh.bst);

    }

    void free_offsets(MeshOffsets offsets, int Nx){
        if (offsets.host_nodeoffset_start){
            for (int i=0; i<Nx+2; i++){
                delete[] offsets.host_nodeoffset_start[i];
            }
            delete[] offsets.host_nodeoffset_start;
        }

        if (offsets.host_nodeoffset_skip_bot){
            for (int i=0; i<Nx+2; i++){
                delete[] offsets.host_nodeoffset_skip_bot[i];
            }
            delete[] offsets.host_nodeoffset_skip_bot;
        }

        if (offsets.host_nodeoffset_skip_mid){
            for (int i=0; i<Nx+2; i++){
                delete[] offsets.host_nodeoffset_skip_mid[i];
            }
            delete[] offsets.host_nodeoffset_skip_mid;
        }

        if (offsets.host_nodeoffset_skip_top){
            for (int i=0; i<Nx+2; i++){
                delete[] offsets.host_nodeoffset_skip_top[i];
            }
            delete[] offsets.host_nodeoffset_skip_top;
        }

        if (offsets.host_elemoffset_start){
            for (int i=0; i<Nx+2; i++){
                delete[] offsets.host_elemoffset_start[i];
            }
            delete[] offsets.host_elemoffset_start;
        }

        if (offsets.host_elemoffset_skip){
            delete[] offsets.host_elemoffset_skip;
        }
    }

    void free_bdry(MeshBdry bdry){
        if (bdry.host_is_bdry_edge){
            delete[] bdry.host_is_bdry_edge;
        }

        if (bdry.host_bdry_edge_normal){
            delete[] bdry.host_bdry_edge_normal;
        }

        if (bdry.host_is_bdry_vert){
            delete[] bdry.host_is_bdry_vert;
        }

        if (bdry.host_bdry_vert_normal){
            delete[] bdry.host_bdry_vert_normal;
        }

        if (bdry.host_edge_to_face){
            delete[] bdry.host_edge_to_face;
        }
    }

    void free_blkif(BlockInterface blkif){
        if (blkif.host_if_x1_r){
            delete[] blkif.host_if_x1_r;
        }

        if (blkif.host_if_x1_node){
            delete[] blkif.host_if_x1_node;
        }

        if (blkif.host_if_x2_r){
            delete[] blkif.host_if_x2_r;
        }

        if (blkif.host_if_x2_node){
            delete[] blkif.host_if_x2_node;
        }

        if (blkif.host_vert_nodeID){
            delete[] blkif.host_vert_nodeID;
        }

        if (blkif.host_vert_subID){
            delete[] blkif.host_vert_subID;
        }

        if (blkif.host_edge_first_nodeID){
            delete[] blkif.host_edge_first_nodeID;
        }

        if (blkif.host_edge_subID){
            delete[] blkif.host_edge_subID;
        }
    }

    void free_bst(MeshBST bst){
        if (bst.host_active_blockID){
            delete[] bst.host_active_blockID;
        }

        if (bst.host_block_nodes_cumulative){
            delete[] bst.host_block_nodes_cumulative;
        }

        if (bst.host_block_elems_cumulative){
            delete[] bst.host_block_elems_cumulative;
        }

        if (bst.host_active_edgeID){
            delete[] bst.host_active_edgeID;
        }

        if (bst.host_edge_nodes_cumulative){
            delete[] bst.host_edge_nodes_cumulative;
        }
    }
}
