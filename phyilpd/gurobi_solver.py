from __future__ import division
from gurobipy import *
import itertools
import numpy as np
import time

def gurobi_solver(node_to_leaves, all_leaves, list_of_ancestral_node, node_to_mean_pwdist, within_cluster_limit, min_cluster_size, verbose):

    print ('Solving with gurobi...')

    # set up indices
    if verbose == 1:
        print ('Setting up indices...')

    leaf_binary_indices = []
    leaf_to_ancestral_nodes = {}

    for n, leaves in node_to_leaves.items():
        for leaf in leaves:
            leaf_binary_indices.append((leaf, n))

            try:
                leaf_to_ancestral_nodes[leaf].append(n)
            except:
                leaf_to_ancestral_nodes[leaf] = [n]

    # model
    model = Model('NewModel')

    # set method
    model.Params.Seed = 666 # Fixed seed to maintain same search path
    model.Params.Method = 1 # always solve by dual simplex (avoid numerical issues)
    model.Params.NumericFocus = 3 # greatest care on numerics (suffer on speed)
    # model.Params.Threads = ncpu_to_use # set number of threads to use

    # verbose
    model.Params.LogToConsole = verbose

    # variables
    if verbose == 1:
        print ('Setting up variables...')
    node_decision = model.addVars(list_of_ancestral_node, vtype=GRB.BINARY)
    leaf_decision = model.addVars(leaf_binary_indices, vtype=GRB.BINARY)

    # update model
    model.update()

    # constraints
    if verbose == 1:
        print ('Setting up constraints...')
    # leaf binary should be related to node binary
    model.addConstrs(leaf_decision[(leaf,n)] <= node_decision[n] for (leaf, n) in leaf_binary_indices)

    # each leaf can only select one node to cluster to
    model.addConstrs(quicksum(leaf_decision[(leaf, n)] for n in leaf_to_ancestral_nodes[leaf]) <= 1 for leaf in all_leaves)

    # leaf always chooses the oldest cluster-able node (transmission clustering)
    for leaf, ancestral_nodes in leaf_to_ancestral_nodes.items():
        for n, m in itertools.combinations(ancestral_nodes, 2):
            if n < m:
                model.addConstr(leaf_decision[(leaf, m)] <= 2 - node_decision[n] - node_decision[m])
            else:
                model.addConstr(leaf_decision[(leaf, n)] <= 2 - node_decision[m] - node_decision[n])

    bigM = 999
    # check that bigM is > min_cluster_size
    if bigM <= min_cluster_size:
        bigM = min_cluster_size + 1

    for n in list_of_ancestral_node:
        # cluster size constraint
        model.addConstr(bigM*(node_decision[n]-1) + min_cluster_size <= quicksum(leaf_decision[(leaf, n)] for leaf in node_to_leaves[n]))
        # within-cluster constraint
        model.addConstr(bigM*(node_decision[n]-1) <= within_cluster_limit - node_to_mean_pwdist[n])

    # update model
    model.update()

    # objective function - maximize number of strains clustered
    if verbose == 1:
        print ('Solving...')

    model.ModelSense = GRB.MAXIMIZE

    model.setObjective(0.5*quicksum(leaf_decision[(leaf, n)] for (leaf, n) in leaf_binary_indices) + 0.5*quicksum(-node_decision[n] for n in list_of_ancestral_node))

    # update model
    model.update()

    # optimize
    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        optimal_solutions = {}

        intfeastol = model.Params.IntFeasTol

        # query optimal solution objective
        optsol_obj = model.getAttr('PoolObjVal')

        # get solution pool size
        solution_pool_size = model.getAttr('SolCount')

        for sol_index in range(solution_pool_size):

            # get solution
            model.params.solutionNumber = sol_index
            curr_sol_obj = model.getAttr('PoolObjVal')

            if curr_sol_obj == optsol_obj:
                fail_integral = 0  # binary to denote if all solution is considered integral

                primaryobj_solution = model.getAttr('Xn', leaf_decision)

                taxon_to_clusterid = {}

                for (leaf, n) in leaf_binary_indices:
                    sol_bin = primaryobj_solution[(leaf, n)]

                    # check for integrality
                    if abs(sol_bin - np.floor(sol_bin+0.5)) > intfeastol:
                        fail_integral = 1
                        # if any solution fails integrality, we can't trust solution (numerical issues)
                        break

                    # round solution to integer
                    if int(round(sol_bin)) == 1:
                        taxon_to_clusterid[leaf] = n

                if fail_integral == 0:
                    optimal_solutions[sol_index] = taxon_to_clusterid
                else:
                    # failed integrality
                    optimal_solutions[sol_index] = False

        # save to solution only if len(optimal_solutions) > 0
        if len(optimal_solutions) > 0:
            return optimal_solutions
    else:
        return 'na'