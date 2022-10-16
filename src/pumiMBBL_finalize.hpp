#ifndef pumiMBBLGPU_finalize_hpp
#define pumiMBBLGPU_finalize_hpp

#include "pumiMBBLGPU.hpp"

namespace pumi{

    void free_mbbl(MBBL pumi_obj);
    void free_mesh(Mesh mesh);
    void free_offsets(MeshOffsets offsets, int Nx);
    void free_bdry(MeshBdry bdry);
    void free_blkif(BlockInterface blkif);
    void free_bst(MeshBST bst);

}

#endif
