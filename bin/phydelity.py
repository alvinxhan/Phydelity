#!/usr/bin/env python

# Phydelity
# Author: Alvin X. Han & Edyth Parker

from __future__ import division
import re
import argparse
import numpy as np
import itertools
import time

import pyximport; pyximport.install()
from phyilpd import node_leaves_reassociation, clean_up_modules, phydelity_output
from phyilpd.tree_utils import parse_newick_tree
from phyilpd.stats_utils import qn, get_cluster_size_distribution

if __name__ == '__main__':
    # parse parameters
    version = 1.0
    parser = argparse.ArgumentParser(description='Phydelity v{}'.format(version))

    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-t', '--tree', type=str, help='Input phylogenetic tree in NEWICK format.')

    analyses_options = parser.add_argument_group('Analysis options')
    analyses_options.add_argument('--k', type=int, help='Custom k neighbours (optional).')
    analyses_options.add_argument('--outgroup', type=str, default=False, help='Taxon (name as appeared in tree) to be set as outgroup OR type \'midpoint\' for midpoint-rooting.')
    analyses_options.add_argument('--collapse_zero_branch_length', action='store_true', help='Collapse internal nodes with zero branch length of tree before running Phydelity.')
    analyses_options.add_argument('--equivalent_zero_length', default=0.000001, type=float, help='Maximum branch length to be rounded to zero if the --collapse_zero_branch_length flag is passed (default = %(default)s).')

    solver_options = parser.add_argument_group('Solver options')
    """
    solver_options.add_argument('--solver', default='gurobi', choices=['glpk', 'gurobi'], type=str, help='Preferred ILP solver IF more than one solvers are available (default: %(default)s).')
    solver_options.add_argument('--solver_verbose', default=0, choices=[0, 1], type=int, help='ILP solver verbose (default: %(default)s)')
    solver_options.add_argument('--solver_check', action='store_true', help='Check available ILP solver(s) installed.')
    """
    solver_options.add_argument('--solver_verbose', default=0, choices=[0, 1], type=int, help='Gurobi solver verbose (default: %(default)s)')
    solver_options.add_argument('--solver_check', action='store_true', help='Check if Gurobi is installed.')

    output_options = parser.add_argument_group('Output options')
    output_options.add_argument('--pdf_tree', action='store_true', help='PDF tree output annotated with cluster results (X server required).')

    # solver_options.add_argument('--threads', type=int, help='Number of threads to use (default = all).')

    params = parser.parse_args()
    # need to change this if we ever update with glpk support
    params.solver = 'gurobi'

    print ('{}\n\n{:^72}\n{:^72}\n\n{}'.format(''.join(['-']*72), 'Phydelity', 'v{}'.format(version), ''.join(['-']*72)))

    # check solver availability
    available_solvers = {}
    """
    try:
        # check for glpsol
        cmd = ['glpsol', '--version']
        solver_version = 'glpk_{}'.format(re.search('v\d+\.\d+', subprocess.check_output(cmd)).group())
        available_solvers['glpk'] = solver_version
    except:
        pass
    """

    try:
        from gurobipy import gurobi
        solver_version = 'gurobi_v{}'.format('.'.join(map(str, gurobi.version())))
        available_solvers['gurobi'] = solver_version
    except:
        pass

    if len(available_solvers) > 0:
        # exit if only to check available ILP solver(s) installed
        if params.solver_check:
            print ('\nAvailable solvers...{}\n'.format(', '.join(available_solvers.values())))
            exit(0)
    else:
        raise Exception('\nNo supported solvers installed. See manual for details on how to download and install supported ILP solvers.\n')

    # check if tree input is given minimally
    if not params.tree:
        raise Exception('\nTree input missing. Check --help for options.\n')

    # preferred solver
    if params.solver not in available_solvers:
        print ('\nWARNING: {} is not installed.'.format(params.solver))
        params.solver = available_solvers.keys()[0]

    print('\nILP solver...{}'.format(available_solvers[params.solver]))

    """
    # limit number of threads for parallelisation
    try:
        ncpu = int(params.threads)
    except:
        from pathos.helpers import mp
        ncpu = mp.cpu_count()
    print ('Threads...{}'.format(ncpu))
    """

    print ('\nParsing tree...')
    # filenames
    treefname = re.sub('([^/]+/|\.[^.]+$)', '', params.tree)

    # parse newick tree file
    try:
        newick_tree_string = parse_newick_tree(params.tree)
    except:
        raise Exception('\nInvalid tree file.\n')

    # parse for tree distance/toplogical information
    from phyilpx import phyilpx_treeinfo
    phyilpx_obj = phyilpx_treeinfo(newick_tree_string, treefname, params.outgroup, params.collapse_zero_branch_length, params.equivalent_zero_length)

    global_leaf_node_id_to_leafname, original_tree_string = phyilpx_obj.properties()

    # get pairwise distance array (nodepair = all nodes including leaves, leafpair = structured, just leaves)
    global_nodepair_to_dist, global_leafpair_to_distance = phyilpx_obj.get_nodepair_distance()

    # get structured array of leaf dist to node / array of node to leaves (reverse sorted by distance to node)
    global_leaf_dist_to_node, global_node_to_leaves = phyilpx_obj.get_leaf_dist_to_node()

    # structured array of node to parent node
    global_node_to_parent_node = phyilpx_obj.get_node_to_parent_node()

    # ancestry relations/global_node_to_mean_pwdist are python dictionaries
    # global_node_to_mean_child_dist2anc = np.array(N,N)
    # global_node_to_pwdist = dictionary of array
    global_node_to_ancestral_nodes, global_node_to_descendant_nodes, global_leaf_to_ancestors, global_node_to_mean_child_dist2anc, global_node_to_pwdist, global_node_to_mean_pwdist = phyilpx_obj.get_ancestral_relations()

    # !-- phydelity specific --#
    print ('\nCalculating distance distribution of closely-related tips...')

    # determine distance distribution of closely-related tips
    if params.k:
        k_range = [params.k]
        print ('WARNING: k fixed at {}.'.format(params.k))
    else:
        k_range = range(2, 6)
        from phyilpx import p_hypotest
        print ('Auto-scaling k...')

    closest_distance_diff_distribution = np.zeros(len(global_leaf_node_id_to_leafname), dtype='f4')

    for _, leaf in enumerate(global_leaf_node_id_to_leafname.keys()):
        # get closest neighbouring leaf of leaf
        j_array = global_leafpair_to_distance[global_leafpair_to_distance['leaf_i'] == leaf][['leaf_j', 'dist']]
        closest_leaf = np.sort(j_array, order='dist')['leaf_j'][0]

        # get diff of distance of mrca to leaf and neighbour leaf
        mrca_node = np.max(list(set(global_leaf_to_ancestors[leaf])&set(global_leaf_to_ancestors[closest_leaf])))
        mrca_leaf_dist = global_leaf_dist_to_node[(leaf, mrca_node)]
        mrca_neighbour_dist = global_leaf_dist_to_node[(closest_leaf, mrca_node)]

        closest_distance_diff_distribution[_] = np.abs(mrca_leaf_dist - mrca_neighbour_dist)

    # get median difference of closest pairwise distances
    median_closest_distance_diff = np.median(closest_distance_diff_distribution)

    # if median closest pair distance < 1, then number of demical place to the 2 significant digit will be the deimcal place to round
    if median_closest_distance_diff < 1:
        decimals_to_round = np.int(np.abs(np.log10(median_closest_distance_diff))) + 2
    else:
        decimals_to_round = 2

    leaf_to_dist_to_csleaf = {}
    for leaf_i, leaf_j in itertools.combinations(global_leaf_node_id_to_leafname.keys(), 2):
        dist = np.round(global_nodepair_to_dist[(leaf_i, leaf_j)], decimals=decimals_to_round)

        try:
            leaf_to_dist_to_csleaf[leaf_i][dist].append(leaf_j)
        except:
            try:
                leaf_to_dist_to_csleaf[leaf_i][dist] = [leaf_j]
            except:
                leaf_to_dist_to_csleaf[leaf_i] = {dist:[leaf_j]}

        try:
            leaf_to_dist_to_csleaf[leaf_j][dist].append(leaf_i)
        except:
            try:
                leaf_to_dist_to_csleaf[leaf_j][dist] = [leaf_i]
            except:
                leaf_to_dist_to_csleaf[leaf_j] = {dist:[leaf_i]}

    for k_strains in k_range:

        leaf_to_kth_sorted_cs_leaves = {}
        for leaf, dist_to_csleaf in leaf_to_dist_to_csleaf.items():
            sorted_pw_distances = sorted(dist_to_csleaf.keys())[:k_strains]
            leaf_to_kth_sorted_cs_leaves[leaf] = [(distance, tuple(dist_to_csleaf[distance])) for distance in sorted_pw_distances]

        core_member_pairwise_leafdist = []
        for leaf, sorted_kth_dist_cleaves in leaf_to_kth_sorted_cs_leaves.items():

            dist_to_add = []

            add_to_core_dist_binary = 0
            for distance, cleaves_tuple in sorted_kth_dist_cleaves:
                dist_to_add.append(distance)

                found_leaf_binary = 0
                for cleaf in cleaves_tuple:
                    for (_distance, _cleaves_tuple) in leaf_to_kth_sorted_cs_leaves[cleaf]:
                        if leaf in _cleaves_tuple:
                            found_leaf_binary = 1
                            break

                    if found_leaf_binary == 1:
                        add_to_core_dist_binary = 1
                        break

                if add_to_core_dist_binary == 1:
                    core_member_pairwise_leafdist += dist_to_add

        med_x = np.median(core_member_pairwise_leafdist)
        mad_x = qn(core_member_pairwise_leafdist)

        core_member_pairwise_leafdist = [_ for _ in core_member_pairwise_leafdist if _ <= med_x + mad_x]

        if len(k_range) > 1: # auto-scaling of k
            if k_strains > 2:

                p_val = p_hypotest(np.array(sorted(set(core_member_pairwise_leafdist)), dtype='f4'),
                                   np.array(sorted(set(prev_core_member_pairwise_distance)), dtype='f4'), 1)

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

    # distal dissociation
    # level-order sorted list of nodes with leaves >= always 2 for transmission clusters (only reassociate such nodes)
    print ('\nDistal dissociation...')
    cs = 2 # min cluster size = pair
    curr_list_of_ancestral_node = np.sort(np.array([node for node, leaves in global_node_to_leaves.items() if len(leaves) >= cs], dtype='i4'))
    
    nla_object = node_leaves_reassociation(cs, wcl, curr_list_of_ancestral_node, global_node_to_leaves, global_node_to_ancestral_nodes, global_node_to_descendant_nodes, global_node_to_mean_pwdist, global_node_to_mean_child_dist2anc, global_node_to_parent_node, global_nodepair_to_dist, global_leaf_dist_to_node, global_leaf_to_ancestors)

    curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist = nla_object.nla_main()

    # remove any nodes with len(leaves) < cs after distal dissociation
    for node, leaves in curr_node_to_leaves.items():
        if len(leaves) < cs:
            del curr_node_to_leaves[node]
            continue

    # update pairwise distance distributions and ancestry relations
    print ('\nUpdating tree info to ILP model...')
    curr_leaves = [x for y in curr_node_to_leaves.values() for x in y]
    curr_list_of_ancestral_node = curr_node_to_leaves.keys()[:]
    curr_node_to_descendant_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}

    # Build ILP model and solve
    from phyilpd.gurobi_solver import gurobi_solver
    all_solutions = gurobi_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_node_to_mean_pwdist, wcl, cs, params.solver_verbose)

    if all_solutions == 'na':
        # continue to next parameter set if no solution
        print ('\nProgram EXIT: No clustering solution found.\n')
        exit(1)

    elif len(all_solutions) > 1:
        print ('\nMultiple ({}) solutions found..'.format(len(all_solutions)))

    # analyse solution and print outputs
    for sol_index, curr_taxon_to_clusterid in all_solutions.items():
        # failed integrality
        if curr_taxon_to_clusterid == False:
            if len(all_solutions) == 1:
                print ('\nProgram EXIT: No optimal solution found (Solver failed integrality).\n')
                exit(1)
            else:
                print ('\nWARNING: No optimal solution found for solution nr. {} (Solver failed integrality).'.format(sol_index))
                continue

        curr_outfname = 'phydelity_k{}_sol{}_{}'.format(str(k_strains), sol_index, treefname)

        curr_clusterid_to_taxa = {}
        for taxon, clusterid in curr_taxon_to_clusterid.items():
            try:
                curr_clusterid_to_taxa[clusterid].append(taxon)
            except:
                curr_clusterid_to_taxa[clusterid] = [taxon]

        print ('\nCleaning up clusters...')
        cleanup_object = clean_up_modules(curr_node_to_descendant_nodes, global_node_to_leaves, curr_node_to_leaves, wcl, cs, global_leaf_dist_to_node, global_leaf_to_ancestors, global_node_to_parent_node, global_nodepair_to_dist)

        # ensure that the most descendant-possible node-id is subtending each cluster
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.transmission_cleanup(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # leave-one-out clean up for clusters violating wcl
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.loo_wcl_violation(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # determine distribution of clusters
        curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa)

        # remove cluster-size sensitivity-induced subclusters
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, 25)

        # de-cluster outliers of descendent cluster that got clustered in ancestor cluster
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.decluster_outlying_taxa_clustered_to_anc_clusters(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # delete any clusters < cs
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.remove_clusters_below_cs(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # print outputs
        print ('\nWriting outputs...')

        output_obj = phydelity_output(original_tree_string, global_leaf_node_id_to_leafname, curr_taxon_to_clusterid, curr_clusterid_to_taxa, curr_outfname)

        # cluster file
        curr_modified_tree_string = output_obj.cluster_output()
        
        # figtree annotated tree file
        output_obj.figtree_output(curr_modified_tree_string)

        # output pdf tree
        if params.pdf_tree:
            output_obj.ete3_pdf_tree_output()

    print ('\n...done.\n')