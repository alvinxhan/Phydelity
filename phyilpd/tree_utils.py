from __future__ import division
import re
import itertools
import numpy as np
from copy import deepcopy as dc
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
    for node in tree_object.traverse():
        if node.is_leaf():
            continue
        elif np.float32(node.dist) <= retain_length_bound:
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

    def __init__(self, min_cluster_size, within_cluster_limit=None, nodes_list=None, node_to_leaves=None, node_to_ancestral_nodes=None, node_to_descendant_nodes=None, node_to_mean_pwdist=None, node_to_mean_child_dist2anc=None, node_to_parent_node=None, nodepair_to_dist=None, leaf_dist_to_node=None, leaf_to_ancestors=None, leafpair_to_distance=None):

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
        self.leafpair_to_distance = leafpair_to_distance

    def get_pwdist_from_leaf_distances_to_node(self, leaves_dist_to_node, desc_node_to_leaves):
        n_i = len(leaves_dist_to_node)

        term_a = (n_i-1)*sum(leaves_dist_to_node)
        term_b = 0

        for desc_node, desc_node_leaves in desc_node_to_leaves.items():
            n_desc = len(desc_node_leaves)
            parent_node = self.node_to_parent_node[self.node_to_parent_node["node"] == desc_node]["parent"].sum()
            term_b += (n_desc)*(n_desc-1)*self.nodepair_to_dist[(parent_node, desc_node)]

        return np.float32(2*(term_a-term_b)/(n_i*(n_i-1)))

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
                remaining_leaves_to_node_dist[leaf] = self.leaf_dist_to_node[(self.leaf_dist_to_node['leaf'] == leaf) & (self.leaf_dist_to_node['node'] == main_node)]['dist'].sum()

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

    def network_conservation(self, tcct_leaves, main_node, desc_nodes):
        from phyilpd.stats_utils import qn

        nodes_to_visit = np.array([main_node] + desc_nodes, dtype='i4')
        node_to_distance_dist = np.zeros((len(nodes_to_visit), len(tcct_leaves)), dtype='f4')

        for _n, node in enumerate(nodes_to_visit):
            for _l, leaf in enumerate(tcct_leaves):

                parent_node_of_leaf = self.node_to_parent_node[self.node_to_parent_node["node"] == leaf]["parent"].sum()
                distance_to_node = self.leaf_dist_to_node[(self.leaf_dist_to_node["leaf"] == leaf) & (self.leaf_dist_to_node["node"] == parent_node_of_leaf)]["dist"].sum()

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
            if main_node != center_node:
                center_node_distance_dist = node_to_distance_dist[np.where(nodes_to_visit == center_node)][0]
                med_x = np.median(center_node_distance_dist)
                mad_x = qn(center_node_distance_dist)

                if len(np.where(center_node_distance_dist > med_x+3*mad_x)[0]) > 0:
                    return False
                
        return True

    def nla_main(self):

        # create deep copies to be edited
        sd_node_to_leaves = dc(self.node_to_leaves)
        sd_node_to_descendant_nodes = dc(self.node_to_descendant_nodes)
        sd_node_to_mean_pwdist = {} # mean pairwise distance dictionary for current set of parameters

        # reassociate subtrees - starting from the root by level order
        for node in self.nodes_list:  # nodes_list already sorted
            # current node <= self.within_cluster_limit
            if self.node_to_mean_pwdist[node] <= self.within_cluster_limit:
                try:
                    # descendant nodes if any
                    tmp_desc_nodes = self.node_to_descendant_nodes[node]
                except:
                    tmp_desc_nodes = []

                # apply network conservation test
                if self.network_conservation(self.node_to_leaves[node], node, tmp_desc_nodes) == True:
                    sd_node_to_mean_pwdist[node] = self.node_to_mean_pwdist[node]
                else:
                    # del node if failed
                    del sd_node_to_leaves[node]
                    try:
                        del sd_node_to_descendant_nodes[node]
                    except:
                        pass
                continue

            # current node > self.within_cluster_limit
            descendant_nodes_to_dissociate = []  # list to save nodes for dissociation
            try:
                # descendants nodes are already sorted by mean child-nodes' distance to node
                for desc_node in self.node_to_descendant_nodes[node]:
                    if desc_node in descendant_nodes_to_dissociate:
                        continue

                    elif self.node_to_mean_child_dist2anc[(desc_node, node)] > self.within_cluster_limit: # node_to_mean_child_dist2anc indexed by desc to anc

                        # append not only desc_node but all descendant nodes of desc_node itself
                        descendant_nodes_to_dissociate.append(desc_node)
                        try:
                            descendant_nodes_to_dissociate = list(set(descendant_nodes_to_dissociate)|set(self.node_to_descendant_nodes[desc_node]))
                        except:
                            pass
                # !--- code con't below

            except:
                # dead-end node with no descendants but > self.within_cluster_limit
                loo_output = self.leave_one_out_leaf_reduction(self.node_to_leaves[node], node)

                # no leaves to remove by leave-one-out-wcl approach
                if loo_output == False:
                    try:
                        # dissociate descendant nodes from node if any
                        tmp_desc_nodes = self.node_to_descendant_nodes[node]
                    except:
                        tmp_desc_nodes = []
                    if self.network_conservation(self.node_to_leaves[node], node, tmp_desc_nodes) == True:
                        sd_node_to_mean_pwdist[node] = self.node_to_mean_pwdist[node]
                    else:
                        del sd_node_to_leaves[node]
                        try:
                            del sd_node_to_descendant_nodes[node]
                        except:
                            pass

                # leave-one-out-wcl approach removes some leaves from dead-end node
                else:
                    leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = loo_output
                    # update node_to_leaves as per output of leave-one-out-wcl approach first
                    sd_node_to_leaves[node] = np.array(leaves_to_keep[:], dtype = 'i4')

                    try:
                        # dissociate descendant nodes from node if any
                        tmp_desc_nodes = sd_node_to_descendant_nodes[node]
                    except:
                        tmp_desc_nodes = []

                    if self.network_conservation(sd_node_to_leaves[node], node, tmp_desc_nodes) == True:
                        # update mean pwd dist
                        sd_node_to_mean_pwdist[node] = mean_pwdist
                    else:
                        del sd_node_to_leaves[node]
                        try:
                            del sd_node_to_descendant_nodes[node]
                        except:
                            pass

                continue

            # code resumed from above ---!#
            leaves_to_remove = list(set([x for y in [self.node_to_leaves[desc_node] for desc_node in descendant_nodes_to_dissociate] for x in y]))
            # remove all leaves from nodes that could be potentially dissociated (leaves in self.node_to_leaves[node] already reverse-sorted by distance to node)
            remaining_leaves = [leaf for leaf in self.node_to_leaves[node] if leaf not in leaves_to_remove]

            loo_output = self.leave_one_out_leaf_reduction(remaining_leaves, node)

            # leave-one-out-wcl approach removes leaves
            if loo_output != False:

                leaves_to_keep, loo_descendant_nodes_to_dissociate, mean_pwdist = loo_output
                # update as per output of leave-one-out-wcl approach
                # update node to remaining leaves
                sd_node_to_leaves[node] = np.array(leaves_to_keep[:], dtype = 'i4')
                # dissociate descendant nodes from node
                descendant_nodes_to_dissociate = list(set(descendant_nodes_to_dissociate)|set(loo_descendant_nodes_to_dissociate))

                sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))
                try:
                    # dissociate descendant nodes from node if any
                    tmp_desc_nodes = sd_node_to_descendant_nodes[node]
                except:
                    tmp_desc_nodes = []
                if self.network_conservation(sd_node_to_leaves[node], node, tmp_desc_nodes) == True:
                    # update mean pwd dist
                    sd_node_to_mean_pwdist[node] = mean_pwdist
                else:
                    del sd_node_to_leaves[node]
                    try:
                        del sd_node_to_descendant_nodes[node]
                    except:
                        pass

            # leave-one-out-wcl approach did not remove any leaves
            else:
                # perform leave-one-out-wcl approach again, now based on self.node_to_leaves[node]
                loo_output = self.leave_one_out_leaf_reduction(self.node_to_leaves[node], node)

                # leave-one-out-wcl approach still did not remove any leaves
                if loo_output == False:
                    try:
                        # dissociate descendant nodes from node if any
                        tmp_desc_nodes = self.node_to_descendant_nodes[node]
                    except:
                        tmp_desc_nodes = []

                    if self.network_conservation(self.node_to_leaves[node], node, tmp_desc_nodes) == True:
                        # update mean pwd dist
                        sd_node_to_mean_pwdist[node] = self.node_to_mean_pwdist[node]
                    else:
                        del sd_node_to_leaves[node]
                        try:
                            del sd_node_to_descendant_nodes[node]
                        except:
                            pass
                # leave-one-out-wcl approach removes leaves
                else:
                    leaves_to_keep, descendant_nodes_to_dissociate, mean_pwdist = loo_output
                    # update as per output of leave-one-out-wcl approach
                    # update node to remaining leaves
                    sd_node_to_leaves[node] = np.array(leaves_to_keep[:], dtype='i4')
                    # dissociate descendant nodes from node
                    sd_node_to_descendant_nodes[node] = list(set(sd_node_to_descendant_nodes[node])-set(descendant_nodes_to_dissociate))

                    try:
                        # dissociate descendant nodes from node if any
                        tmp_desc_nodes = sd_node_to_descendant_nodes[node]
                    except:
                        tmp_desc_nodes = []
                    if self.network_conservation(sd_node_to_leaves[node], node, tmp_desc_nodes) == True:
                        # update mean pwd dist
                        sd_node_to_mean_pwdist[node] = mean_pwdist
                    else:
                        del sd_node_to_leaves[node]
                        try:
                            del sd_node_to_descendant_nodes[node]
                        except:
                            pass

        return sd_node_to_leaves, sd_node_to_descendant_nodes, sd_node_to_mean_pwdist

# clean-up modules
class clean_up_modules():

    def __init__(self, current_node_to_descendant_nodes=None, node_to_leaves=None, leafpair_to_distance=None, current_node_to_leaves=None, within_cluster_limit=None, min_cluster_size=None, leaf_dist_to_node=None, leaf_to_ancestors=None, node_to_parent_node=None, nodepair_to_dist=None):

        self.current_node_to_descendant_nodes = current_node_to_descendant_nodes
        self.node_to_leaves = node_to_leaves
        self.leafpair_to_distance = leafpair_to_distance
        self.current_node_to_leaves = current_node_to_leaves
        self.within_cluster_limit = within_cluster_limit
        self.min_cluster_size = min_cluster_size
        self.leaf_dist_to_node = leaf_dist_to_node
        self.leaf_to_ancestors = leaf_to_ancestors
        self.node_to_parent_node = node_to_parent_node
        self.nodepair_to_dist = nodepair_to_dist

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

    def get_pwdist_from_leaf_distances_to_node_cleanup(self, leaves_dist_to_node, desc_node_to_leaves):
        n_i = len(leaves_dist_to_node)

        term_a = (n_i-1)*sum(leaves_dist_to_node)
        term_b = 0
        for desc_node, desc_node_leaves in desc_node_to_leaves.items():
            n_desc = len(desc_node_leaves)
            parent_node = self.node_to_parent_node[self.node_to_parent_node["node"] == desc_node]["parent"].sum()
            term_b += (n_desc)*(n_desc-1)*self.nodepair_to_dist[(parent_node, desc_node)]

        return np.float32(2*(term_a-term_b)/(n_i*(n_i-1)))

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
                remaining_leaves_to_node_dist[leaf] = self.leaf_dist_to_node[(self.leaf_dist_to_node["leaf"] == leaf) & (self.leaf_dist_to_node["node"] == main_node)]["dist"].sum()
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
            new_pwdist = np.zeros(len(leafpairs_list), dtype='f4')
            for _, (i,j) in enumerate(leafpairs_list):
                new_pwdist[_] = self.nodepair_to_dist[(i, j)]

            if np.mean(new_pwdist) > self.within_cluster_limit:
                # reverse sort clustered taxa by distance to node
                leaf_dist_of_cluster = self.leaf_dist_to_node[self.leaf_dist_to_node['node'] == clusterid][["leaf", "dist"]]
                rsorted_taxa = np.sort(leaf_dist_of_cluster, order="dist")["leaf"][::-1]

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
                    combined_mean_pwd = np.mean([self.nodepair_to_dist[(x, y)] for x, y in itertools.combinations(parent_taxa + taxa, 2)])
                    if combined_mean_pwd <= self.within_cluster_limit:
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