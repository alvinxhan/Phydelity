from __future__ import division
import re
import itertools
import numpy as np
from copy import deepcopy as dc
from phyilpd.stats_utils import qn
import time

def parse_newick_tree(tree_file):

    fhandle = filter(None, [_.strip() for _ in open(tree_file, 'rU')])
    fhandle = fhandle.pop(0)
    new_fhandle = []
    prev_end = 0
    for expr in re.finditer('(\(|,)([^(,:]+):', fhandle):
        new_fhandle.append(fhandle[prev_end:expr.start()+1])
        new_fhandle.append("'{}'".format(re.sub("(^'|'$)", "", expr.group(2))))
        prev_end = expr.end()-1
    new_fhandle.append(fhandle[prev_end:])

    return ''.join(new_fhandle)

def collapse_zero_branch_length(tree_object, retain_length_bound, treefname):
    '''
    Collapse nodes with zero branch lengths and reorder tree
    '''
    no_node_collapse = 1
    diff = max([retain_length_bound, 1e-6])
    for node in tree_object.traverse():
        if node.is_leaf():
            continue
        else:
            if np.absolute(np.float64(node.dist) - 2*retain_length_bound) <= 2*diff:
                no_node_collapse = 0
                node.delete()

    if no_node_collapse == 1:
        return False # no branches collapsed
    else:
        tree_object.ladderize()

        # write zero branch length collapsed tree as a newick file
        treefname = 'zero-branch-length-collapsed_{}'.format(treefname)
        tree_object.write(format=5, outfile='{}.nwk'.format(treefname))
        print ('Collapsed newick tree file...{}.nwk'.format(treefname))

        return tree_object

# class of modules to reassociate node and leaves by current min cluster size and within-cluster limit
class node_leaves_reassociation():

    def __init__(self, min_cluster_size, within_cluster_limit=None, nodes_list=None, node_to_leaves=None, node_to_ancestral_nodes=None, node_to_descendant_nodes=None, node_to_mean_pwdist=None, node_to_mean_child_dist2anc=None, node_to_parent_node=None, nodepair_to_dist=None, leaf_dist_to_node=None, leaf_to_ancestors=None, equivalent_zero_length=None, leaf_node_id_to_name=None, node_to_child_leaves=None, leaf_node_id_to_leafname=None):

        self.min_cluster_size = min_cluster_size
        self.within_cluster_limit = within_cluster_limit
        self.nodes_list = nodes_list
        self.node_to_leaves = node_to_leaves
        self.node_to_ancestral_nodes = node_to_ancestral_nodes
        self.node_to_descendant_nodes = node_to_descendant_nodes
        self.node_to_mean_pwdist = node_to_mean_pwdist
        self.node_to_mean_child_dist2anc = node_to_mean_child_dist2anc
        self.node_to_parent_node = node_to_parent_node
        self.nodepair_to_dist = nodepair_to_dist
        self.leaf_dist_to_node = leaf_dist_to_node
        self.leaf_to_ancestors = leaf_to_ancestors
        self.equivalent_zero_length = equivalent_zero_length
        self.leaf_node_id_to_name = leaf_node_id_to_name
        self.node_to_child_leaves = node_to_child_leaves
        self.leaf_node_id_to_leafname = leaf_node_id_to_leafname

    def get_pwdist_from_leaf_distances_to_node(self, leaves_dist_to_node, desc_node_to_leaves):
        n_i = len(leaves_dist_to_node)

        term_a = (n_i-1)*sum(leaves_dist_to_node)
        term_b = 0

        for desc_node, desc_node_leaves in desc_node_to_leaves.items():
            n_desc = len(desc_node_leaves)
            parent_node = self.node_to_parent_node[desc_node]
            term_b += (n_desc)*(n_desc-1)*self.nodepair_to_dist[(parent_node, desc_node)]

        return np.float64(2*(term_a-term_b)/(n_i*(n_i-1)))

    def leave_one_out_leaf_reduction(self, sorted_leaves, main_node):
        # note that sorted_leaves is already reverse sorted by distance to main_node

        if len(sorted_leaves) < self.min_cluster_size:
            # immediate return if len(sorted_leaves) < self.min_cluster_size
            return False

        for l_index in range(-1, len(sorted_leaves), 1):
            # return if len(sorted_leaves) after reduction == self.min_cluster_size
            if l_index + 1 == len(sorted_leaves) - self.min_cluster_size:
                return False
                break

            remaining_leaves_to_node_dist = {}
            remaining_descendant_nodes_to_leaves = {}

            # start from not having any leaves removed
            for leaf in sorted_leaves[l_index+1:]:
                remaining_leaves_to_node_dist[leaf] = self.leaf_dist_to_node[(leaf, main_node)]

                # get all descendant nodes of main node subtending leaf
                try:
                    desc_nodes_subtending_leaf = list(set(self.node_to_descendant_nodes[main_node])&set(self.leaf_to_ancestors[leaf]))
                except:
                    desc_nodes_subtending_leaf = []

                # save as desc node to leaves
                for rdn in desc_nodes_subtending_leaf:
                    try:
                        remaining_descendant_nodes_to_leaves[rdn].append(leaf)
                    except:
                        remaining_descendant_nodes_to_leaves[rdn] = [leaf]

            reduced_mean_pwdist = self.get_pwdist_from_leaf_distances_to_node(remaining_leaves_to_node_dist.values(), remaining_descendant_nodes_to_leaves)

            # break if <= self.within_cluster_limit and dissociate
            if reduced_mean_pwdist <= self.within_cluster_limit:
                # find nodes of dissociated leaves
                leaves_dissociated = sorted_leaves[:l_index+1]

                if len(leaves_dissociated) == 0:
                    old_descendant_nodes_to_dissociate = []
                else:
                    try:
                        nodes_of_dissociated_leaves = list(set([x for y in [list(set(self.node_to_descendant_nodes[main_node])&set(self.leaf_to_ancestors[leaf])) for leaf in leaves_dissociated] for x in y]))
                        old_descendant_nodes_to_dissociate = [desc_node for desc_node in list(set(nodes_of_dissociated_leaves) - set(remaining_descendant_nodes_to_leaves.keys())) if set(self.node_to_leaves[desc_node]) <= set(leaves_dissociated)]
                    except:
                        old_descendant_nodes_to_dissociate = []

                # return leaves to keep for node AND descendant nodes (which all of its subtended leaves are to be removed) to dissociate
                return sorted_leaves[l_index+1:], old_descendant_nodes_to_dissociate, reduced_mean_pwdist

        return False

    def remove_distantly_related_nodes(self, list_to_analyse, main_node, desc_node_remaining=None, pwdist=None):

        N = len(list_to_analyse) # node list length
        if pwdist != None and N == 2: # leaves analysis and and we have a pair
            if pwdist <= self.within_cluster_limit:
                return True, list_to_analyse, desc_node_remaining, pwdist
            else:
                return False, False

        """
        if pwdist != None:
            for _ in list_to_analyse:
                if re.search('IMH0028', self.leaf_node_id_to_leafname[_]):
                    print main_node
                    print map(lambda __: self.leaf_node_id_to_leafname[__], list_to_analyse)
        """

        # get pairwise distance array of node list
        pairwise_d_array = np.full((N, N), np.nan, dtype=np.float64)

        for _i, _j in itertools.combinations(range(N), 2):
            leaf_i = list_to_analyse[_i]
            leaf_j = list_to_analyse[_j]
            dist = self.nodepair_to_dist[(leaf_i, leaf_j)]
            pairwise_d_array[(_i, _j)] = pairwise_d_array[(_j, _i)] = dist

        # get i-th wise mean pairwise distance
        i_wise_mean_pd_array = np.zeros(N, dtype=np.float64)
        for _i in xrange(N):
            i_wise_pd = pairwise_d_array[(_i,)]
            i_wise_mean_pd_array[_i] = np.mean(i_wise_pd[~np.isnan(i_wise_pd)])

        # calculate median and mad of pairwise nodes distance distribution
        tril_d_array = np.tril(pairwise_d_array)
        tril_d_array = tril_d_array[(~np.isnan(tril_d_array))&(tril_d_array > 0.)]
        try:
            mad_x = qn(tril_d_array)
        except:
            mad_x = 0.
        med_x = np.median(tril_d_array)
        dispersal_threshold = med_x + 3 * mad_x

        if np.absolute(dispersal_threshold - 2 * self.equivalent_zero_length) <= self.equivalent_zero_length:
            dispersal_nodes_to_stay = [_i for _i in range(N) if i_wise_mean_pd_array[_i] <= self.within_cluster_limit]
        else:
            dispersal_nodes_to_stay = [_i for _i in range(N) if i_wise_mean_pd_array[_i] <= min([dispersal_threshold, self.within_cluster_limit])]

        if len(dispersal_nodes_to_stay) == 0:
            #print main_node
            return False, False

        set_distances = list(set([pairwise_d_array[(_i, _j)] for _i, _j in itertools.combinations(dispersal_nodes_to_stay, 2)]))
        if all(np.absolute(dist - 2 * self.equivalent_zero_length) <= self.equivalent_zero_length for dist in set_distances):
            dispersal_nodes_to_stay = range(N)

        # nodes = internal nodes
        if pwdist == None:

            # mean distance between the i-th internal node and all others must be <= min(med_x+3*mad_x, wcl)
            """
            if np.isclose(dispersal_threshold, 2 * self.equivalent_zero_length):
                nodes_to_remove = np.array([list_to_analyse[_i] for _i in range(N) if i_wise_mean_pd_array[_i] > self.within_cluster_limit])
            else:
                nodes_to_remove = np.array([list_to_analyse[_i] for _i in range(N) if i_wise_mean_pd_array[_i] > min([dispersal_threshold, self.within_cluster_limit])])
            """
            nodes_to_remove = np.array([list_to_analyse[_i] for _i in list(set(range(N)) - set(dispersal_nodes_to_stay))])

            # delete dispersal nodes to stay if the main_node is the odd node
            if main_node in nodes_to_remove:
                return True, np.array([list_to_analyse[_i] for _i in dispersal_nodes_to_stay])

            # no nodes to remove
            if len(nodes_to_remove) == 0:
                return True, False
            else:
                # return internal nodes to remove
                return True, nodes_to_remove

        # nodes = leaves
        else:
            """
            if N != len(dispersal_nodes_to_stay): # if there are leaves removed by dispersal
                # get children tips of main node with the shortest distance to said node.
                children_leaves_of_main_node_to_dist = {_i:self.nodepair_to_dist[(list_to_analyse[_i], main_node)] for _i in range(N) if self.node_to_parent_node[list_to_analyse[_i]] == main_node}
                children_leaves_of_main_node = children_leaves_of_main_node_to_dist.keys()
                children_leaves_of_main_node_w_min_dist = [_i for _i, dist in children_leaves_of_main_node_to_dist.items() if dist == min(children_leaves_of_main_node_to_dist.values())]

                dispersal_nodes_removed = list(set(range(N)) - set(dispersal_nodes_to_stay))

                if set(children_leaves_of_main_node) <= set(dispersal_nodes_removed) : # all children of main_node removed by dispersal
                    return False, False

                _i_to_j = {}
                for _i, _j in itertools.product(dispersal_nodes_removed, dispersal_nodes_to_stay):
                    if pairwise_d_array[(_i, _j)] > np.min([dispersal_threshold, self.within_cluster_limit]):
                        try:
                            _i_to_j[_i].append(_j)
                        except:
                            _i_to_j[_i] = [_j]

                for _i, _j_list in _i_to_j.items():
                    # if _i is a child of the main_node, remove it
                    if _i in children_leaves_of_main_node:
                        continue

                    # if any (_i, _j) pair violates the distance threshold and that _j is a child of the main_node, _i should be removed
                    set_of_main_node_child_j_with_longest_dist_to_i = [self.leaf_node_id_to_name[list_to_analyse[_j]] for _j in set(_j_list)&set(children_leaves_of_main_node_w_min_dist)]

                    #print self.leaf_node_id_to_name[list_to_analyse[_i]], set_of_main_node_child_j_with_longest_dist_to_i
                    if len(set_of_main_node_child_j_with_longest_dist_to_i) == 0:
                        dispersal_nodes_to_stay.append(_i)
                        continue

            parent_node_of_remaining_leaves = [self.node_to_parent_node[list_to_analyse[_i]] for _i in dispersal_nodes_to_stay]

            if main_node not in parent_node_of_remaining_leaves: # all direct children are removed
                return False, False

            
            # get direct children node of main node
            direct_children_leaves_indices = [_i for _, _i in enumerate(dispersal_nodes_to_stay) if parent_node_of_remaining_leaves[_] == main_node]

            # cut away any remaining leaves which pwdist is > wcl from the direct children leaves of the main_node
            non_child_leaves_indices = list(set(dispersal_nodes_to_stay) - set(direct_children_leaves_indices))
            for _i, _j in itertools.product(direct_children_leaves_indices, non_child_leaves_indices):
                if pairwise_d_array[(_i, _j)] > self.within_cluster_limit:
                    remove_from_dispersal.append(_j)
            """
            #final_nodes_to_stay = list(set(dispersal_nodes_to_stay) - set(remove_from_dispersal))
            #nodes_to_remove = np.array([list_to_analyse[_i] for _i in list(set(range(N)) - set(final_nodes_to_stay))])

            nodes_to_remove = np.array([list_to_analyse[_i] for _i in list(set(range(N)) - set(dispersal_nodes_to_stay))])

            """
            if main_node == 275:
                print map(lambda _: self.leaf_node_id_to_leafname[_], nodes_to_remove), main_node
            """

            # if there are leaves removed
            if len(nodes_to_remove) > 0:

                # remove node if all direct child leaves of node were removed...
                if set(self.node_to_child_leaves[main_node]) <= set(nodes_to_remove):
                    #print main_node
                    return False, False

                # ...or if the number of remaining leaves < min cluster size
                list_to_analyse = np.array(list(set(list_to_analyse) - set(nodes_to_remove)))
                if len(list_to_analyse) < self.min_cluster_size:
                    #print main_node
                    return False, False

                # remove any descendant nodes that wholly subtend the leaves to be removed.
                if len(desc_node_remaining) > 0:
                    desc_node_remaining = list(set(desc_node_remaining) - set([dnode for dnode in desc_node_remaining if set(self.node_to_leaves[dnode])&set(list_to_analyse) <= set(nodes_to_remove)]))

                # re-calculate mean pairwise distance with remaining leaves
                if len(list_to_analyse) == 2:
                    pwdist = self.nodepair_to_dist[(list_to_analyse[0], list_to_analyse[1])]
                else:
                    pwdist = np.mean([self.nodepair_to_dist[(leaf_i, leaf_j)] for leaf_i, leaf_j in itertools.combinations(list_to_analyse, 2)])

            return True, list_to_analyse, desc_node_remaining, pwdist

    def network_conservation(self, tcct_leaves, main_node, desc_nodes):

        # check that mean distances between internal nodes within a cluster would not exceed wcl
        #nodes_to_visit = np.array([main_node] + desc_nodes, dtype=np.int64)
        nodes_to_visit = np.array([main_node], dtype=np.int64)

        # check for if a center node exist that is equidistant to all leaves (non-cherries)
        node_to_distance_dist = np.zeros((len(nodes_to_visit), len(tcct_leaves)), dtype=np.float64)

        for _n, node in enumerate(nodes_to_visit):
            for _l, leaf in enumerate(tcct_leaves):

                parent_node_of_leaf = self.node_to_parent_node[leaf]
                distance_to_node = self.leaf_dist_to_node[(leaf, parent_node_of_leaf)]

                if node != parent_node_of_leaf:
                    try:
                        mrca_node = max(set(self.node_to_ancestral_nodes[parent_node_of_leaf])&set(self.node_to_ancestral_nodes[node]))
                    except:
                        if parent_node_of_leaf == 0 or node == 0:
                            mrca_node = 0 # root
                        else:
                            raise Exception('No common ancestor found between subtrees. Contact Alvin at hanxc@bii.a-star.edu.sg with your tree/paramters.')

                    # CA of parent node and node == main_node
                    if mrca_node == main_node:
                        distance_to_node += self.nodepair_to_dist[(parent_node_of_leaf, node)] + self.nodepair_to_dist[(main_node, node)]
                    else:
                        distance_to_node += abs(self.nodepair_to_dist[(main_node, parent_node_of_leaf)] - self.nodepair_to_dist[(main_node, node)])

                node_to_distance_dist[(_n, _l)] = distance_to_node

        node_to_median_distance_dist = np.apply_along_axis(np.median, 1, node_to_distance_dist)
        min_distance_dist = np.min(node_to_median_distance_dist)
        nodes_with_min_dist = nodes_to_visit[np.where(node_to_median_distance_dist <= min_distance_dist)]

        for center_node in nodes_with_min_dist:
            #if main_node != center_node:
            center_node_distance_dist = node_to_distance_dist[np.where(nodes_to_visit == center_node)][0]
            #med_x = np.median(center_node_distance_dist)
            #mad_x = qn(center_node_distance_dist)

            if len(np.where(center_node_distance_dist > self.within_cluster_limit)[0]) > 0:
            #if len(np.where(center_node_distance_dist > max([med_x + 3*mad_x, self.within_cluster_limit]))[0]) > 0:
                return False
        return True

    def nla_main(self):

        # create deep copies to be edited
        sd_node_to_leaves = dc(self.node_to_leaves)
        sd_node_to_descendant_nodes = dc(self.node_to_descendant_nodes)
        sd_node_to_mean_pwdist = {} # mean pairwise distance dictionary for current set of parameters

        # reassociate subtrees - starting from the root by level order
        for node in self.nodes_list:  # nodes_list already sorted

            # leaves to node
            leaves_to_node = sd_node_to_leaves[node]

            """
            if any(re.search("NUH0018", self.leaf_node_id_to_leafname[leaf]) for leaf in leaves_to_node):
                print node, "*"
            """

            try:
                # get descendant nodes if any
                desc_nodes_list = sd_node_to_descendant_nodes[node]
            except:
                desc_nodes_list = []

            # delete_node_binary
            delete_node_binary = 0

            # current node <= self.within_cluster_limit
            if self.node_to_mean_pwdist[node] <= self.within_cluster_limit:
                mean_pwdist = self.node_to_mean_pwdist[node]

            # current node > self.within_cluster_limit
            else:
                if len(desc_nodes_list) > 0: # nodes with descendant nodes

                    # perform leave one out for remaining leaves
                    loo_output = self.leave_one_out_leaf_reduction(leaves_to_node, node)

                    # leave-one-out-wcl approach removes leaves
                    if loo_output != False:

                        leaves_to_node, loo_descendant_nodes_to_dissociate, mean_pwdist = loo_output
                        # dissociate descendant nodes to be removed from node
                        desc_nodes_list = list(set(desc_nodes_list) - set(loo_descendant_nodes_to_dissociate))

                    # leave-one-out-wcl approach did not remove any of the remaining leaves
                    else:
                        # perform leave-one-out-wcl approach again based on entire leaves_to_node
                        loo_output = self.leave_one_out_leaf_reduction(leaves_to_node, node)

                        # leave-one-out-wcl approach still did not remove any leaves
                        if loo_output == False:
                            mean_pwdist = self.node_to_mean_pwdist[node]

                        # leave-one-out-wcl approach removes leaves
                        else:
                            leaves_to_node, loo_descendant_nodes_to_dissociate, mean_pwdist = loo_output
                            # dissociate descendant nodes to be removed from node
                            desc_nodes_list = list(set(desc_nodes_list) - set(loo_descendant_nodes_to_dissociate))

                else: # dead-end node with no descendants but > self.within_cluster_limit
                    # see if we can reduce mean mean_pwdist by leave one out approach
                    loo_output = self.leave_one_out_leaf_reduction(leaves_to_node, node)

                    # no leaves to remove by leave-one-out-wcl approach
                    if loo_output == False:
                        mean_pwdist = self.node_to_mean_pwdist[node]

                    # leave-one-out-wcl approach removes some leaves from dead-end node
                    else:
                        leaves_to_node, mean_pwdist = loo_output[0], loo_output[-1]

            # remove any distantly related leaves that violates wcl
            rdrn_leaves_output = self.remove_distantly_related_nodes(leaves_to_node, node, desc_nodes_list, mean_pwdist)

            '''
            if node == 275:
                print rdrn_leaves_output
            '''

            if rdrn_leaves_output[0] == False:
                # delete node if either all child leaves of node have been removed or node membership < min cluster size
                delete_node_binary = 1
            else:
                # leaves remaining
                remaining_leaves, desc_nodes_list, mean_pwdist = rdrn_leaves_output[1:]

                nodes_to_visit = [node] + desc_nodes_list
                if len(nodes_to_visit) > 2:
                    # check for distantly related descendant nodes which mean inter-nodal distance > min([med_x + 3*mad_x, self.within_cluster_limit])
                    rdrn_nodes_output = self.remove_distantly_related_nodes(nodes_to_visit, node)

                    if rdrn_nodes_output[0] == False: # main_node is a node to be removed
                        delete_node_binary = 1
                    else:
                        rdrn_descendant_nodes_to_dissociate = rdrn_nodes_output[1]
                        if isinstance(rdrn_descendant_nodes_to_dissociate, bool) and rdrn_descendant_nodes_to_dissociate == False:
                            pass
                        else:
                            # remove any descendant nodes flagged
                            desc_nodes_list = list(set(desc_nodes_list) - set(rdrn_descendant_nodes_to_dissociate))

                            # get remaining leaves subtended by node
                            leaves_to_remove = list(set([x for y in [self.node_to_leaves[desc_node] for desc_node in rdrn_descendant_nodes_to_dissociate] for x in y]))
                            remaining_leaves = np.array(list(set(remaining_leaves) - set(leaves_to_remove)))

                            if len(remaining_leaves) < self.min_cluster_size: # del node if number of remaining leaves < min cluster size
                                delete_node_binary = 1
                            else:
                                if len(remaining_leaves) == 2:
                                    mean_pwdist = self.nodepair_to_dist[(remaining_leaves[0], remaining_leaves[1])]
                                else:
                                    mean_pwdist = np.mean([self.nodepair_to_dist[(leaf_i, leaf_j)] for leaf_i, leaf_j in itertools.combinations(remaining_leaves, 2)])

            if delete_node_binary == 0 and len(remaining_leaves) > 2:

                if self.network_conservation(remaining_leaves, node, desc_nodes_list) == False:
                    delete_node_binary = 1

            # delete node if binary == 1 or else save/update node information
            if delete_node_binary == 1:
                del sd_node_to_leaves[node]
                try:
                    del sd_node_to_descendant_nodes[node]
                except:
                    pass
            else:
                sd_node_to_mean_pwdist[node] = mean_pwdist
                sd_node_to_leaves[node] = remaining_leaves
                sd_node_to_descendant_nodes[node] = desc_nodes_list

        return sd_node_to_leaves, sd_node_to_descendant_nodes, sd_node_to_mean_pwdist

# clean-up modules
class clean_up_modules():

    def __init__(self, current_node_to_descendant_nodes=None, node_to_leaves=None, current_node_to_leaves=None, within_cluster_limit=None, min_cluster_size=None, leaf_dist_to_node=None, leaf_to_ancestors=None, node_to_parent_node=None, nodepair_to_dist=None, node_to_child_leaves=None, leaf_node_id_to_leafname=None):

        self.current_node_to_descendant_nodes = current_node_to_descendant_nodes
        self.node_to_leaves = node_to_leaves
        self.current_node_to_leaves = current_node_to_leaves
        self.within_cluster_limit = within_cluster_limit
        self.min_cluster_size = min_cluster_size
        self.leaf_dist_to_node = leaf_dist_to_node
        self.leaf_to_ancestors = leaf_to_ancestors
        self.node_to_parent_node = node_to_parent_node
        self.nodepair_to_dist = nodepair_to_dist
        self.node_to_child_leaves = node_to_child_leaves
        self.leaf_node_id_to_leafname = leaf_node_id_to_leafname

    def transmission_cleanup(self, clusterid_to_taxa, taxon_to_clusterid):
        for cluster in sorted(clusterid_to_taxa.keys()): # clusters are ordered by ancestry
            # get descendant internal nodes of cluster which are clusters as well
            try:
                desc_nodes_clustered = list(set(self.current_node_to_descendant_nodes[cluster])&set(clusterid_to_taxa.keys()))
            except:
                continue

            # if there are desc clusters subtending from cluster node
            if len(desc_nodes_clustered) > 0:
                for desc_cluster in sorted(desc_nodes_clustered): # desc clusters are ordered by ancestry
                    # check to see if leaves of cluster overlap with leaves subtended by desc cluster
                    taxa_clustered_in_anc_but_subtended_by_desc = list(set(clusterid_to_taxa[cluster])&set(self.current_node_to_leaves[desc_cluster]))
                    if len(taxa_clustered_in_anc_but_subtended_by_desc) > 0:
                        # let desc cluster be the node clustering anc cluster leaves overlapped
                        clusterid_to_taxa[cluster] = list(set(clusterid_to_taxa[cluster]) - set(taxa_clustered_in_anc_but_subtended_by_desc))

                        clusterid_to_taxa[desc_cluster] = clusterid_to_taxa[desc_cluster] + taxa_clustered_in_anc_but_subtended_by_desc

                        for taxon in taxa_clustered_in_anc_but_subtended_by_desc:
                            taxon_to_clusterid[taxon] = desc_cluster

        # remove if cluster falls below min cluster size
        for cluster, taxa in clusterid_to_taxa.items():
            if len(taxa) < self.min_cluster_size:
                del clusterid_to_taxa[cluster]
                for taxon in taxa:
                    del taxon_to_clusterid[taxon]

        return clusterid_to_taxa, taxon_to_clusterid

    def remove_odd_leaf(self, clusterid_to_taxa, taxon_to_clusterid):
        from phyilpd.stats_utils import qn

        for clusterid, taxa in clusterid_to_taxa.items():

            if len(taxa) <= 3: # skip anything less than/equals to a trio
                continue

            # check mean pairwise distance
            sorted_taxa = sorted(taxa)
            N = len(sorted_taxa)
            new_pwdist = np.full((N, N), np.nan, dtype=np.float64)
            for _i, _j in itertools.combinations(range(N), 2):
                new_pwdist[(_i, _j)] = new_pwdist[(_j, _i)] = self.nodepair_to_dist[(sorted_taxa[_i], sorted_taxa[_j])]

            leaves_flagged = {}
            # get the median distance of pairwise leaf distance for each leaf to all of its cluster mates
            for i in xrange(N):
                pwdist_of_leaf = new_pwdist[(i,)]
                pwdist_of_leaf = pwdist_of_leaf[~np.isnan(pwdist_of_leaf)]
                med_x = np.median(pwdist_of_leaf)
                mad_x = qn([_ for _ in pwdist_of_leaf if _ >= med_x])
                wc_tol = med_x + 5*mad_x # set to 5x to be CLEAR outliers

                for j in xrange(N):
                    if j == i:
                        continue
                    if new_pwdist[(i,j)] > wc_tol:
                        try:
                            leaves_flagged[sorted_taxa[j]] += 1
                        except:
                            leaves_flagged[sorted_taxa[j]] = 1

            leaves_flagged = [k for k,v in leaves_flagged.items() if v/(N-1) > 0.9]

            if len(leaves_flagged) > 0:
                clusterid_to_taxa[clusterid] = list(set(taxa) - set(leaves_flagged))
                for leaf in leaves_flagged:
                    del taxon_to_clusterid[leaf]

        clusterid_to_taxa, taxon_to_clusterid = self.remove_clusters_below_cs(clusterid_to_taxa, taxon_to_clusterid)

        return clusterid_to_taxa, taxon_to_clusterid

    def get_pwdist_from_leaf_distances_to_node_cleanup(self, leaves_dist_to_node, desc_node_to_leaves):
        n_i = len(leaves_dist_to_node)

        term_a = (n_i-1)*sum(leaves_dist_to_node)
        term_b = 0
        for desc_node, desc_node_leaves in desc_node_to_leaves.items():
            n_desc = len(desc_node_leaves)
            parent_node = self.node_to_parent_node[desc_node]
            term_b += (n_desc)*(n_desc-1)*self.nodepair_to_dist[(parent_node, desc_node)]

        return np.float64(2*(term_a-term_b)/(n_i*(n_i-1)))

    def leave_one_out_leaf_reduction_cleanup(self, sorted_leaves, main_node):
        # note that sorted_leaves is already reverse sorted by distance to main_node
        # immediate return if len(sorted_leaves) < self.min_cluster_size
        if len(sorted_leaves) < self.min_cluster_size:
            return False, None

        for l_index in range(-1, len(sorted_leaves), 1):
            if len(sorted_leaves[l_index+1:]) < self.min_cluster_size:
                return False, None

            remaining_leaves_to_node_dist = {}
            remaining_descendant_nodes_to_leaves = {}

            # start from not having any leaves removed
            for leaf in sorted_leaves[l_index + 1:]:
                remaining_leaves_to_node_dist[leaf] = self.leaf_dist_to_node[(leaf, main_node)]
                try:
                    desc_nodes_subtending_leaf = list(set(self.current_node_to_descendant_nodes[main_node])&set(self.leaf_to_ancestors[leaf]))
                except:
                    desc_nodes_subtending_leaf = []

                for rdn in desc_nodes_subtending_leaf:
                    try:
                        remaining_descendant_nodes_to_leaves[rdn].append(leaf)
                    except:
                        remaining_descendant_nodes_to_leaves[rdn] = [leaf]

            reduced_mean_pwdist = self.get_pwdist_from_leaf_distances_to_node_cleanup(remaining_leaves_to_node_dist.values(), remaining_descendant_nodes_to_leaves)

            # break if <= self.within_cluster_limit and dissociate
            if reduced_mean_pwdist <= self.within_cluster_limit:
                # return leaves to keep for node
                return True, sorted_leaves[l_index+1:]

        return False, None

    def loo_wcl_violation(self, clusterid_to_taxa, taxon_to_clusterid):

        for clusterid, taxa in clusterid_to_taxa.items():
            # check mean pairwise distance
            leafpairs_list = list(itertools.combinations(taxa, 2))
            new_pwdist = np.zeros(len(leafpairs_list), dtype=np.float64)
            for _, (i,j) in enumerate(leafpairs_list):
                new_pwdist[_] = self.nodepair_to_dist[(i, j)]

            if np.mean(new_pwdist) > self.within_cluster_limit:
                # reverse sort clustered taxa by distance to node
                leaf_dist_of_cluster = {}
                for leaf in taxa:
                    leaf_dist_of_cluster[leaf] = self.leaf_dist_to_node[(leaf, clusterid)]
                rsorted_taxa = sorted(leaf_dist_of_cluster.keys(), key=leaf_dist_of_cluster.get, reverse=True)

                loo_output_binary, loo_output = self.leave_one_out_leaf_reduction_cleanup(rsorted_taxa, clusterid)

                if loo_output_binary == False:
                    # remvoe entire cluster since it entirely violates the within cluster limit
                    del clusterid_to_taxa[clusterid]
                    for taxon in taxa:
                        del taxon_to_clusterid[taxon]
                else:
                    # if we could still have a cluster after removing "outlying" taxa
                    for taxon in list(set(taxa) - set(loo_output)):
                        del taxon_to_clusterid[taxon]
                        clusterid_to_taxa[clusterid].remove(taxon)

        return clusterid_to_taxa, taxon_to_clusterid

    def subsume_subclusters_under_x_percentile(self, clusterid_to_taxa, taxon_to_clusterid, clusterlen_distribution, percentile):

        #subsumed_taxa_to_clusterid = {}

        # regard any clusters of length <= 25 (default) percentile lacking evidence to be a potential trajectory
        clusterlen_cutoff = np.percentile(clusterlen_distribution, percentile)
        if percentile < 100:
            print ('Subsuming cluster-size sensitivity-induced subclusters <= {} taxa (if any)...'.format(int(clusterlen_cutoff)))

        # determine ancestry of clusters
        cluster_to_desc_clusters = {}
        for cluster in sorted(clusterid_to_taxa.keys()):
            try:
                desc_clusters = list(set(self.current_node_to_descendant_nodes[cluster])&set(clusterid_to_taxa.keys()))
                if len(desc_clusters) > 0:
                    cluster_to_desc_clusters[cluster] = desc_clusters
            except:
                continue

        cluster_to_anc_clusters  = {}
        for cluster, desc_clusters in cluster_to_desc_clusters.items():
            for desc in desc_clusters:
                try:
                    cluster_to_anc_clusters[desc].append(cluster)
                except:
                    cluster_to_anc_clusters[desc] = [cluster]

        # check nodes which are <= x-percentile and do not have any descending clusters
        for cluster in sorted(clusterid_to_taxa.keys()):
            taxa = clusterid_to_taxa[cluster]
            if len(taxa) <= clusterlen_cutoff and cluster not in cluster_to_desc_clusters:
                try:
                    parent_cluster = sorted(cluster_to_anc_clusters[cluster])[-1]
                    parent_taxa = clusterid_to_taxa[parent_cluster][:]
                except:
                    continue

                # check that taxa set <= self.current_node_to_leaves[parent_cluster] (in other words, the parent cluster, a statistically significant cluster, is only so with the inclusion of the taxa set as well)
                if set(taxa) < set(self.current_node_to_leaves[parent_cluster]):
                    # check that if we subsume cluster into immediate parent cluster, it would still fulfill self.within_cluster_limit
                    #combined_mean_pwd = np.mean([self.nodepair_to_dist[(x, y)] for x, y in itertools.combinations(parent_taxa + taxa, 2)])
                    #if combined_mean_pwd <= self.within_cluster_limit:
                    if all(self.nodepair_to_dist[(x, y)] <= self.within_cluster_limit for x, y in itertools.product(taxa, parent_taxa)):
                        for taxon in taxa:
                            taxon_to_clusterid[taxon] = parent_cluster
                            #subsumed_taxa_to_clusterid[taxon] = cluster

                        clusterid_to_taxa[parent_cluster] = list(set(parent_taxa)|set(taxa))
                        del clusterid_to_taxa[cluster]

        return clusterid_to_taxa, taxon_to_clusterid #, subsumed_taxa_to_clusterid

    def remove_clusters_below_cs(self, clusterid_to_taxa, taxon_to_clusterid):
        # remove any clusters < min cluster size
        for clusterid, taxa in clusterid_to_taxa.items():
            if len(taxa) < self.min_cluster_size:
                for taxon in taxa:
                    del taxon_to_clusterid[taxon]
                del clusterid_to_taxa[clusterid]
        return clusterid_to_taxa, taxon_to_clusterid

    def decluster_outlying_taxa_clustered_to_anc_clusters(self, clusterid_to_taxa, taxon_to_clusterid):
        clusterid_to_skip = []
        # start from the most ancestral to the most descendant clusters
        for clusterid in sorted(clusterid_to_taxa.keys()):
            taxa = clusterid_to_taxa[clusterid]
            # skip if all leaves subtended by cluster node is of the same set as those parsed in the tree
            if clusterid in clusterid_to_skip or set(taxa) == set(self.node_to_leaves[clusterid]):
                continue

            # get descendant cluster nodes of current cluster
            try:
                descendant_clusterids = list(set(self.current_node_to_descendant_nodes[clusterid])&set(clusterid_to_taxa.keys()))
            except:
                continue

            if len(descendant_clusterids) > 0:
                for desc_clusterid in descendant_clusterids:
                    # skip if all leaves subtended by cluster node is of the same set as those parsed in the tree
                    if set(self.node_to_leaves[desc_clusterid]) == set(clusterid_to_taxa[desc_clusterid]):
                        clusterid_to_skip.append(desc_clusterid)
                        continue

                    # leaf outliers of nodes appears in global node_to_leaves dict but not in current
                    outliers_of_desc_clusterid = list(set(self.node_to_leaves[desc_clusterid])-set(self.current_node_to_leaves[desc_clusterid]))
                    # decluster any outlying leaves of descending nodes that are clustered in current node
                    outliers_of_desc_clusterid_clustered_in_clusterid = list(set(taxa)&set(outliers_of_desc_clusterid))

                    if len(outliers_of_desc_clusterid_clustered_in_clusterid) > 0:
                        for outlier in outliers_of_desc_clusterid_clustered_in_clusterid:
                            del taxon_to_clusterid[outlier]
                            clusterid_to_taxa[clusterid].remove(outlier)
                            # delete any empty clusters
                            if len(clusterid_to_taxa[clusterid]) == 0:
                                del clusterid_to_taxa[clusterid]

        return clusterid_to_taxa, taxon_to_clusterid