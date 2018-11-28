from __future__ import division
import re
import numpy as np
import itertools
import cython
from scipy.misc import comb as nCr
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
cimport numpy as np

from ete3 import Tree
import time

from phyilpd.tree_utils import collapse_zero_branch_length

cdef extern from "math.h":
    float exp "exp" (float)
    float sqrt "sqrt" (float)
    float fabs "fabs" (float)

def p_hypotest(data1, data2, method):
    if method == 1:
        return kuiper(data1, data2)

@cython.boundscheck(False)
cdef float hypotest(float[:] data1, float[:] data2, int method):
    if method == 1:
        return kuiper(data1, data2)

@cython.boundscheck(False)
cdef float _kuiper_dist(float q):

    cdef float EPS1 = 0.001
    cdef float EPS2 = 0.00000001

    cdef float a2
    cdef float term
    cdef float sum = 0.
    cdef float termbf = 0.
    cdef float dist = 1.

    cdef unsigned int j

    a2 = -2*q*q

    for j in range(1, 10001):
        term = 2*((4*j*j*q*q)-1)*exp(a2*j*j)
        sum += term

        if fabs(term) <= EPS1*termbf or fabs(term) <= EPS2*sum:
            dist = sum
            break

        termbf = fabs(term)
    return dist

@cython.boundscheck(False)
cdef float kuiper(float [:] data1, float [:] data2):
    cdef int j1 = 0
    cdef int j2 = 0
    cdef int n1
    cdef int n2
    cdef float eff_n
    cdef float fn1 = 0.
    cdef float fn2 = 0.
    cdef float d1
    cdef float d2
    cdef float dt
    cdef float lam
    cdef float d = 0.

    n1 = len(data1)  # number of elements in data1
    n2 = len(data2)  # number of elements in data2
    eff_n = sqrt((n1*n2)/(n1+n2))

    while j1 < n1 and j2 < n2:
        d1 = data1[j1]
        d2 = data2[j2]

        if d1 <= d2:
            j1 += 1
            fn1 = j1/n1

        if d2 <= d1:
            j2 += 1
            fn2 = j2/n2

        dt = fabs(fn2-fn1)
        if dt > d:
            d = dt

    lam = (eff_n + 0.155 + 0.24/eff_n)*d

    return _kuiper_dist(lam)

    """
    def _ks_2samp(self, data1, data2):
        from scipy.stats import kstwobign

        # change data list to array
        data1 = np.asarray(data1, dtype=np.float32)
        data2 = np.asarray(data2, dtype=np.float32)

        n1 = len(data1) # number of elements in data1
        n2 = len(data2) # number of elements in data2

        data_all = np.concatenate([data1,data2]) # concatenate both data arrays and sort - this is basically the array of all possible values

        cdf1 = np.searchsorted(data1,data_all,side='right')/n1
        cdf2 = np.searchsorted(data2,data_all,side='right')/n2

        d = np.max(np.absolute(cdf1-cdf2), dtype=np.float32)
        # Note: d absolute not signed distance
        en = np.sqrt((n1*n2)/(n1+n2), dtype=np.float32)

        try:
            prob = kstwobign.sf((en + 0.12 + 0.11 / en) * d)
        except:
            prob = 1.0
        return prob
    """

cdef struct Node:
    int parent
    int* children
    int children_length
    float edge_distance

@cython.no_gc_clear
cdef class phyilpx_treeinfo:

    cdef Node* data
    cdef unsigned int depth
    cdef unsigned int total_nr_nodes
    cdef object leafname_to_leaf_nodeid
    cdef object leaf_nodeid_to_leafname
    cdef object internalnodes
    cdef object np_buffer
    cdef object node_to_leaves
    cdef object node_to_pwdist
    cdef str original_tree_string
    #cdef np.ndarray node_to_pwdist
    cdef np.ndarray leaf_dist_to_node
    cdef np.ndarray nodepair_to_distance

    #cdef object ete_nodeid_to_node

    def __init__(self, newick_tree, treefname, outgroup, collapse_zero_branch_length_binary, eq_zero_branch_length):

        cdef int node_id

        self.np_buffer = None # numpy memory buffer

        # read newick tree as ete tree
        tree = Tree(newick_tree)

        # changing tree root
        if outgroup != False:

            if outgroup == 'midpoint': # midpoint rooting
                root_node = tree.get_midpoint_outgroup()
            else:
                # search for outgroup by name
                try:
                    root_node = tree&"'{}'".format(outgroup)
                except:
                    raise Exception('\nGiven outgroup {} did not match any taxon node in tree.\n'.format(outgroup))

            tree.set_outgroup(root_node)
            print ('\nWARNING: Tree outgroup changed to {}.'.format(outgroup))

        tree.ladderize()  # ladderize tree

        # collapse zero branch length
        if collapse_zero_branch_length_binary:
            tree = collapse_zero_branch_length(tree, np.float32(eq_zero_branch_length), treefname)

            if tree == False:
                raise Exception('No branches were collapsed. Check the upper limit of zero-length branches on your tree and adjust accordingly using --equivalent_zero_length')

        # get original tree string
        self.original_tree_string = tree.write(format=5)

        treesize = len(tree.get_descendants()) + 1 # size of tree (all nodes + leaves)
        self.data = <Node*> PyMem_Malloc(treesize*sizeof(Node)) # create data array and allocate memory
        if not self.data:
            raise MemoryError()

        self.leaf_nodeid_to_leafname = {}
        self.leafname_to_leaf_nodeid = {}
        self.internalnodes = []

        #self.ete_nodeid_to_node = {}

        for node_id, node in enumerate(tree.traverse(strategy='levelorder')):

            #self.ete_nodeid_to_node[node_id] = node
            node.add_feature('node_id', node_id)

            if node.is_leaf() :
                self.leaf_nodeid_to_leafname[node_id] = node.name
                self.leafname_to_leaf_nodeid[node.name] = node_id
            else:
                self.internalnodes.append(node_id)

        self.total_nr_nodes = node_id + 1 # total number of nodes

        for node_id, node in enumerate(tree.traverse(strategy='levelorder')):
            try:
                self.data[node_id].parent = node.up.node_id
            except:
                self.data[node_id].parent = -1 # root

            self.data[node_id].edge_distance = node.dist

            if node.is_leaf():
                children_list = []
            else:
                children_list = [child_node.node_id for child_node in node.get_children()]

            self.data[node_id].children_length = len(children_list)

            self.data[node_id].children = <int*> PyMem_Malloc(len(children_list) * sizeof(int))
            if not self.data[node_id].children:
                raise MemoryError()

            for c, child in enumerate(children_list):
                self.data[node_id].children[c] = child

        # determine the maximum depth of the tree
        for node_id in self.leaf_nodeid_to_leafname.keys():
            n = 1
            while True :
                if self.data[node_id].parent == -1: # break when parent is root
                    break
                node_id = self.data[node_id].parent
                n += 1
            if n > self.depth :
                self.depth = n

    def get_children(self, node_id):
        # return node ids of child nodes for a given node
        return [self.data[node_id].children[c] for c in xrange(self.data[node_id].children_length)]

    def get_leaves(self, node_id):
        # return array of leaves subtended by a given node
        cdef unsigned int i
        cdef unsigned int n = 0

        if self.np_buffer is None:
            self.np_buffer = np.ndarray(len(self.leaf_nodeid_to_leafname), dtype=int)

        to_visit = [node_id]
        for i in to_visit :
            children_list = self.get_children(i)
            if len(children_list) == 0: # is leaf
                self.np_buffer[n] = i
                n += 1
            else:
                to_visit += children_list
        return np.array(self.np_buffer[:n])

    def mrca(self, a, b) :
        # return mrca of 2 nodes
        visited = np.zeros(self.depth, dtype=int)

        return self._mrca(visited, a, b)

    @cython.boundscheck(False)
    cdef int _mrca(self, long[:] visited, int a, int b) nogil :
        cdef int n
        cdef int i
        cdef int mrca = -1
        cdef int a_depth

        n = a
        i = 0
        while True :
            visited[i] = n
            n = self.data[n].parent
            i += 1
            if n == -1 : break
        a_depth = i

        n = b
        while True :
            i = 0
            while True :
                if i >= a_depth : break
                if visited[i] == n :
                    mrca = visited[i]
                    break
                i += 1
            if mrca != -1 : break
            n = self.data[n].parent
            if n == -1 :
                mrca = n
                break

        for i in xrange(self.depth) :
            visited[i] = -1

        return mrca

    @cython.boundscheck(False)
    cdef float get_distance(self, int a, int b):
        cdef int mrca
        cdef float d = 0
        cdef int n

        mrca = self.mrca(a, b)

        n = a
        while n != mrca :
            d += self.data[n].edge_distance
            n =  self.data[n].parent
        n = b
        while n != mrca :
            d += self.data[n].edge_distance
            n =  self.data[n].parent
        return d

    def is_ancestor(self, a, b):
        # returns 1 if a is ancestor of b, -1 if b is ancestor of a, 0 if otherwise 
        return self._is_ancestor( a, b )

    @cython.boundscheck(False)
    cdef int _is_ancestor(self, int a, int b) nogil :
        cdef int i
        cdef int n

        # is a an ancestor of b?
        i = b
        while True :
            n = self.data[i].parent
            if n == -1 :
                break
            if n == a :
                return 1
            i = n

        # is b an ancestor of a?
        i = a
        while True :
            n = self.data[i].parent
            if n == -1 :
                break
            if n == b :
                return -1
            i = n

        # or neither?
        return 0

    def properties(self):
        return self.leaf_nodeid_to_leafname, self.original_tree_string

    def get_nodepair_distance(self):
        cdef unsigned int i
        cdef unsigned int j
        cdef unsigned int k
        cdef unsigned int N = self.total_nr_nodes # length of all nodes
        cdef float dist

        self.nodepair_to_distance = np.zeros((N,N), dtype='f4')

        cdef unsigned int N_leafpairs = 2*nCr(len(self.leaf_nodeid_to_leafname), 2)
        cdef np.ndarray leafpair_to_distance = np.zeros(N_leafpairs,
                                                        dtype = {'names':('leaf_i', 'leaf_j', 'dist'), 'formats':('i4', 'i4', 'f4')}) # structured array

        k = 0
        # pairwise patristic distance array of all nodes present
        for (i, j) in itertools.combinations(range(N), 2):
            dist = self.get_distance(i, j)

            if self.data[i].children_length == 0 and self.data[j].children_length == 0: # leaf pairs
                leafpair_to_distance[k] = (i,j,dist)
                leafpair_to_distance[k+1] = (j,i,dist)
                k += 2

            self.nodepair_to_distance[(i,j)] = dist
            self.nodepair_to_distance[(j,i)] = dist

        return self.nodepair_to_distance, leafpair_to_distance

    def get_leaf_dist_to_node(self):
        cdef unsigned int i, node_id
        cdef unsigned int k = 0
        #cdef unsigned int N = 0
        cdef unsigned int N = self.total_nr_nodes
        cdef float dist

        cdef object leaf_to_ancestors

        # get leaf to internal node distances
        leaf_to_ancestors = {}
        for i, node_id in itertools.product(self.leaf_nodeid_to_leafname.keys(), self.internalnodes):
            if self.is_ancestor(node_id, i) == 1:
                try:
                    leaf_to_ancestors[i].append(node_id)
                except:
                    leaf_to_ancestors[i] = [node_id]
                #N += 1

        # structured array of node-leaf distance
        #self.leaf_dist_to_node = np.zeros(N, dtype={'names':('leaf', 'node', 'dist'), 'formats':('i4', 'i4', 'f4')})
        self.leaf_dist_to_node = np.full((N, N), np.nan, dtype = 'f4') # leaf v node
        for i in leaf_to_ancestors.keys():
            for node_id in leaf_to_ancestors[i]:
                dist = self.nodepair_to_distance[(i, node_id)]
                #self.leaf_dist_to_node[k] = (i, node_id, dist)
                self.leaf_dist_to_node[(i, node_id)] = dist
                k += 1

        # dictionary of leaf arrays reverse sorted by distance to node
        self.node_to_leaves = {}
        cdef np.ndarray curr_struc_array, leafdist_to_node_arr, leaves_of_node
        for node_id in self.internalnodes:
            #curr_struc_array = self.leaf_dist_to_node[self.leaf_dist_to_node['node'] == node_id][['leaf', 'dist']]
            #self.node_to_leaves[node_id] = np.sort(curr_struc_array, order='dist')['leaf'][::-1]

            curr_struc_array = self.leaf_dist_to_node[:,node_id]
            leafdist_to_node_arr = curr_struc_array[np.isnan(curr_struc_array) == False]
            leaves_of_node = np.where(np.isnan(curr_struc_array) == False)[0]
            self.node_to_leaves[node_id] = leaves_of_node[leafdist_to_node_arr.argsort()][::-1]

        return self.leaf_dist_to_node, self.node_to_leaves

    def get_node_to_parent_node(self):
        cdef unsigned int node_id, parent, i
        cdef unsigned int N = self.total_nr_nodes - 1
        #cdef np.ndarray node_to_parent_node = np.zeros(N, dtype={'names':('node', 'parent'), 'formats':('i4', 'i4')})
        cdef object node_to_parent_node = {}

        for i, node_id in enumerate(range(1, self.total_nr_nodes)):
            parent = self.data[node_id].parent
            #node_to_parent_node[i] = (node_id, parent)
            node_to_parent_node[node_id] = parent

        return node_to_parent_node

    def get_node_to_mean_child_dist2root(self):

        cdef unsigned int node_id
        cdef unsigned int child_node_id
        cdef unsigned int i
        cdef unsigned int k
        cdef float dist

        cdef object children_of_node

        cdef unsigned int N = len(self.internalnodes)
        cdef np.ndarray child_dist_to_node
        cdef np.ndarray node_to_mean_child_dist2root = np.zeros(N, dtype={'names':('node', 'dist'), 'formats':('i4', 'f4')})

        for k, node_id in enumerate(self.internalnodes):
            children_of_node = self.get_children(node_id)

            child_dist_to_node = np.zeros(len(children_of_node), dtype='f4')
            for i, child_node_id in enumerate(children_of_node):

                dist = self.nodepair_to_distance[(child_node_id, 0)]
                child_dist_to_node[i] = dist

            node_to_mean_child_dist2root[k] = (node_id, np.mean(child_dist_to_node))

        return node_to_mean_child_dist2root

    def get_ancestral_relations(self):

        cdef object entry
        cdef unsigned int i
        cdef unsigned int j
        cdef unsigned int k
        cdef unsigned int node_id
        cdef unsigned int N = len(self.internalnodes) + len(self.leaf_nodeid_to_leafname)
        cdef float mean_dist_to_root

        cdef object leaf_to_ancestors
        cdef object node_to_ancestral_nodes
        cdef object node_to_descendant_nodes
        cdef object ancestors_to_n
        cdef object leafpairs_to_node
        cdef object node_to_mean_pwdist

        cdef np.ndarray leaves_to_node
        cdef np.ndarray node_to_mean_child_dist2anc = np.zeros((N,N), dtype='f4')

        node_to_mean_child_dist2root = self.get_node_to_mean_child_dist2root() # get mean child node distances to root for every internal node

        leaf_to_ancestors = {}
        self.node_to_pwdist = {}
        node_to_mean_pwdist = {}
        node_to_ancestral_nodes = {}
        node_to_descendant_nodes = {}

        for entry in np.sort(node_to_mean_child_dist2root, order='dist'): # starting node with the lowest mean children distance to root
            node_id, mean_dist_to_root = entry

            leaves_to_node = self.node_to_leaves[node_id]
            for i in leaves_to_node:
                try:
                    leaf_to_ancestors[i].append(node_id)
                except:
                    leaf_to_ancestors[i] = [node_id]

            ancestors_to_n = [i for i in self.internalnodes[::-1] if self.is_ancestor(i, node_id) == 1]

            node_to_ancestral_nodes[node_id] = ancestors_to_n[:]

            for i in ancestors_to_n:
                try:
                    node_to_descendant_nodes[i].append(node_id)
                except:
                    node_to_descendant_nodes[i] = [node_id]

                node_to_mean_child_dist2anc[(node_id, i)] = mean_dist_to_root - self.nodepair_to_distance[i][0] # indexed by node to anc

            leafpairs_to_node = list(itertools.combinations(leaves_to_node, 2))
            pwdist_to_node = np.zeros(len(leafpairs_to_node), dtype='f4')
            for k, (i,j) in enumerate(leafpairs_to_node):
                pwdist_to_node[k] = self.nodepair_to_distance[(i,j)]

            self.node_to_pwdist[node_id] = np.sort(pwdist_to_node)
            node_to_mean_pwdist[node_id] = np.mean(pwdist_to_node)

        return node_to_ancestral_nodes, node_to_descendant_nodes, leaf_to_ancestors, node_to_mean_child_dist2anc, self.node_to_pwdist, node_to_mean_pwdist

    def get_global_pval(self, hytest_method):
        cdef unsigned int i
        cdef unsigned int j
        cdef unsigned int x
        cdef unsigned int y
        cdef unsigned int k
        cdef unsigned int N = len(self.internalnodes) + len(self.leaf_nodeid_to_leafname)
        cdef float pval
        cdef float pval0
        cdef float pval1
        cdef np.ndarray ij_product_pwdist
        cdef np.ndarray ij_pwdist
        cdef np.ndarray nodepair_to_pval = np.zeros((N, N), dtype='f4')
        cdef object leaves_product

        for i, j in itertools.combinations(self.internalnodes, 2):

            if self.is_ancestor(i, j) == 0: # neither i or j are ancestors of each other's

                leaves_product = list(itertools.product(self.node_to_leaves[i], self.node_to_leaves[j]))

                ij_product_pwdist = np.zeros(len(leaves_product), dtype='f4')
                for k, (x,y) in enumerate(leaves_product):
                    ij_product_pwdist[k] = self.nodepair_to_distance[(x,y)]
                ij_pwdist = np.concatenate((self.node_to_pwdist[i], self.node_to_pwdist[j], ij_product_pwdist))
                ij_pwdist = np.sort(ij_pwdist)

                # take the conservative (max) p-value comparing node i/j individually to i+j
                pval0 = hypotest(self.node_to_pwdist[i], ij_pwdist, hytest_method)
                pval1 = hypotest(self.node_to_pwdist[j], ij_pwdist, hytest_method)
                if pval0 > pval1:
                    pval = pval0
                else:
                    pval = pval1
            else:
                pval = hypotest(self.node_to_pwdist[i], self.node_to_pwdist[j], hytest_method)

            nodepair_to_pval[(i,j)] = nodepair_to_pval[(j,i)] = pval

    def __dealloc__(self):
        PyMem_Free(self.data)  # no-op if self.data is NULL