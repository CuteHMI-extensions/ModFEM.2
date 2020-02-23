#ifndef DDH_INTF_H_
#define DDH_INTF_H_

#ifdef __cplusplus
extern "C"{
#endif

/** @defgroup DDD Domain decomposition
 *
 *  @{
 */

enum ddt_mesh_handling {
    DDC_LOCAL_MESH = 0,
    DDC_GLOBAL_MESH = 1,
    DDC_DEFAULT = 2
};

/**
 * @brief ddr_create_subdomains to create partitioning scheme with (or without) overlap
 * Depending on the param Mesh_handling partitioning is done locally or globally
 * @param IN Mesh_id
 * @param IN Mesh_handling
 *      DDC_LOCAL_MESH - creates subdomains only from local mesh
 *      DDC_GLOBAL_MESH - creates subdomains from all procs' meshes
 * @param IN N_subdomains
 * @param OUT N_subdomains_elems - sizes (in core elements) for each subdomain
 *  preferred size: N_subdomains
 * @param OUT Subdomains_elems - core elements stored in continogus manner for each subdomain
 *  preferred size: number of all elements (LOCAL or GLOBAL)
 * @param OUT Overlap_sizes - number of overlap element for each subdomain
 * preferred size: N_subdomains
 * @param OUT Subdomains_elems_overlap - elements creating overlap for each subdomain
 *  preferred size: same as Subdomains_elems
 * @return number of created subdomains
 */
int ddr_create_subdomains_scheme(const int Mesh_id,
              ddt_mesh_handling Mesh_handling,
              int N_subdomains,
              int * N_subdomains_elems,
              int *Subdomains_elems,
              int * Overlap_sizes,
              int * Subdomains_elems_overlap,
              int ** Part_ptr);



/**
 * @brief ddr_balance_load to improve existing partitioning created by ddr_create_subdomains
 * according to changes in mesh structures.
 *
 * @param IN Mesh_id
 * @param OUT N_subdomains
 * @param OUT N_subdomains_elems
 * @param OUT Subdomains_elems
 * @param OUT Subdomains_elems_overlap
 * @return
 */
int ddr_balance_subdomains(const int Mesh_id,
                     int * N_subdomains_elems,
                     int *Subdomains_elems,
                     int ** Subdomains_elems_overlap
                     );

/**
 *
 *  @}
 */

#ifdef __cplusplus
}
#endif

#endif //DDH_INTF_H_
