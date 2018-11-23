from __future__ import division
import re
import numpy as np
import itertools
import cython
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
cimport numpy as np
from phyilpd.stats_utils import qn

cdef int K
cdef object leaf_node_ids
cdef object leaf_to_ancestors
cdef np.ndarray leafpair_to_distance
cdef np.ndarray leaf_dist_to_node

def get_closely_related_wcl(K, leaf_node_ids, leafpair_to_distance, leaf_to_ancestors, leaf_dist_to_node):

    cdef unsigned int _
    cdef unsigned int leaf 
    cdef unsigned int closest_leaf
    cdef unsigned int mrca_node
    cdef unsigned int decimals_to_round
    cdef float mrca_leaf_dist
    cdef float mrca_neighbour_dist
    cdef float median_closest_distance_diff
    cdef np.ndarray j_array
    cdef np.ndarray closest_distance_diff_distribution = np.zeros(len(leaf_node_ids), dtype='f4')
    cdef object k_range

    # determine distance distribution of closely-related tips
    if K:
        k_range = [K]
        print ('WARNING: k fixed at {}.'.format(K))
    else:
        k_range = range(2, 6)
        from phyilpx import p_hypotest
        print ('Auto-scaling k...')

    for _, leaf in enumerate(leaf_node_ids):
        # get closest neighbouring leaf of leaf
        j_array = leafpair_to_distance[leafpair_to_distance['leaf_i'] == leaf][['leaf_j', 'dist']]
        closest_leaf = np.sort(j_array, order='dist')['leaf_j'][0]

        # get diff of distance of mrca to leaf and neighbour leaf
        mrca_node = np.max(list(set(leaf_to_ancestors[leaf])&set(leaf_to_ancestors[closest_leaf])))
        mrca_leaf_dist = leaf_dist_to_node[(leaf_dist_to_node["leaf"] == leaf) & (leaf_dist_to_node["node"] == mrca_node)]["dist"].sum()
        mrca_neighbour_dist = leaf_dist_to_node[(leaf_dist_to_node["leaf"] == closest_leaf) & (leaf_dist_to_node["node"] == mrca_node)]["dist"].sum()

        closest_distance_diff_distribution[_] = np.abs(mrca_leaf_dist - mrca_neighbour_dist)

    # get median difference of closest pairwise distances
    median_closest_distance_diff = np.median(closest_distance_diff_distribution)

    # if median closest pair distance < 1, then number of demical place to the 2 significant digit will be the deimcal place to round
    if median_closest_distance_diff < 1:
        decimals_to_round = np.int(np.abs(np.log10(median_closest_distance_diff))) + 1
    else:
        decimals_to_round = 1

    cdef int k_strains
    cdef int i
    cdef int found_leaf
    cdef int neighbor_leaf
    cdef float dist
    cdef float neighbor_dist
    cdef float med_x
    cdef float mad_x
    cdef float p_val
    cdef float wcl
    cdef np.ndarray core_member_pairwise_leafdist
    cdef np.ndarray prev_core_member_pairwise_distance
    cdef np.ndarray j_array_dist
    cdef np.ndarray unique_rounded_k_distances
    cdef np.ndarray neighbor_leaves_with_dist
    cdef np.ndarray neighbor_j_array
    cdef np.ndarray neighbor_j_array_dist
    cdef np.ndarray neighbor_unique_rounded_k_distances
    cdef np.ndarray additional_closely_related_distances

    for k_strains in k_range:
        core_member_pairwise_leafdist = np.zeros(len(leaf_node_ids)*k_strains, dtype='f4')
        i = 0

        for leaf in leaf_node_ids:
            j_array = leafpair_to_distance[leafpair_to_distance['leaf_i'] == leaf]['leaf_j']
            j_array_dist = np.around(leafpair_to_distance[leafpair_to_distance['leaf_i'] == leaf]['dist'], decimals=decimals_to_round)

            unique_rounded_k_distances = np.sort(np.unique(j_array_dist))[:k_strains]
            found_leaf = 0

            for _, dist in enumerate(unique_rounded_k_distances):
                neighbor_leaves_with_dist = np.extract(j_array_dist == dist, j_array)

                for neighbor_leaf in neighbor_leaves_with_dist:
                    neighbor_j_array = leafpair_to_distance[leafpair_to_distance['leaf_i'] == neighbor_leaf]['leaf_j']
                    neighbor_j_array_dist = np.around(leafpair_to_distance[leafpair_to_distance['leaf_i'] == neighbor_leaf]['dist'], decimals=decimals_to_round)
                    neighbor_unique_rounded_k_distances = np.sort(np.unique(neighbor_j_array_dist))[:k_strains]

                    for neighbor_dist in neighbor_unique_rounded_k_distances:
                        if leaf in np.extract(neighbor_j_array_dist == neighbor_dist, neighbor_j_array):
                            found_leaf = 1
                            break

                    if found_leaf == 1:
                        break

                if found_leaf == 1:
                    core_member_pairwise_leafdist[i:i+(_+1)] = unique_rounded_k_distances[:_+1]
                    i += (_ + 1)
                    break

        core_member_pairwise_leafdist = np.sort(core_member_pairwise_leafdist[:i])

        med_x = np.median(core_member_pairwise_leafdist)
        mad_x = qn(core_member_pairwise_leafdist)
        core_member_pairwise_leafdist = core_member_pairwise_leafdist[np.where(core_member_pairwise_leafdist <= med_x + mad_x)]

        if len(k_range) > 1: # auto-scaling of k
            if k_strains > 2:
                additional_closely_related_distances = np.setdiff1d(core_member_pairwise_leafdist, prev_core_member_pairwise_distance)
                if additional_closely_related_distances.size == 0:
                    if k_strains == 5:
                        wcl = np.amax(core_member_pairwise_leafdist)
                    continue
                else:
                    p_val = p_hypotest(prev_core_member_pairwise_distance, additional_closely_related_distances, 1)
                    if p_val < 0.05:
                        wcl = np.amax(prev_core_member_pairwise_distance)
                        k_strains -= 1
                        break
                    elif k_strains == 5:
                        wcl = np.amax(core_member_pairwise_leafdist)
                        break
            prev_core_member_pairwise_distance = core_member_pairwise_leafdist[:]
        else:
            wcl = np.amax(core_member_pairwise_leafdist)

    return wcl, k_strains