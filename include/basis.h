/**
 * basis.h — SAFB basis construction
 *
 * Translated from: Sg_init/space_group_plane_family.py
 *
 * Core function: build_basis() — assembles SAFBBasis from symmetry ops,
 * lattice parameters, and desired mode count.
 */

#ifndef BASIS_H
#define BASIS_H

#include "symmetry_ops.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * Build SAFB basis from symmetry operations, lattice, and mode count
 *
 * Parameters:
 *   ops         — symmetry operation group
 *   sg_symbol   — space group symbol (e.g. "Ia-3d")
 *   lattice     — LatticeInfo
 *   N           — number of modes to generate
 *   result      — output SAFBBasis (caller must ensure safb is NULL or
 *                 basis_build frees any previous safb)
 *
 * Returns: 0 on success, -1 on error.
 * ======================================================================== */

int basis_build(const SymmGroup *ops,
                const char *sg_symbol,
                const LatticeInfo *lattice,
                int N,
                SAFBBasis *result);

#ifdef __cplusplus
}
#endif

#endif /* BASIS_H */
