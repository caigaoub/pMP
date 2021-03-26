import pulp
import numpy as np
import sys
import DataHandler as dh


def _get_solver(param_gap):
    cplex = pulp.CPLEX_CMD(msg=True,
                           path='/home/cplex',
                           options=['set mip strategy probe 3'])
    if cplex.available():
        return cplex
    return pulp.PULP_CBC_CMD(msg=True, fracGap=param_gap)


def _solve_pMp(distances, n_dcs, decimals):
    distances = np.around(distances, decimals)
    prob = pulp.LpProblem("MIPOptimization", pulp.LpMinimize)

    customers = list(range(distances.shape[0]))
    sites = list(range(distances.shape[1]))

    # The assignment to each DC
    x = pulp.LpVariable.dicts('X',
                              ((c, s) for s in sites for c in customers),
                              lowBound=0,
                              cat=pulp.LpContinuous)

    # If the site is used
    y = pulp.LpVariable.dicts('Y',
                              (s for s in sites),
                              0, 1, cat=pulp.LpBinary)

    # Objective function to minimize - the flow to any DC 
    prob += (pulp.lpSum(x[(c, s)] * distances[c, s] for s in sites for c in customers),
             'Objective')

    # Each customer's flows should add to 1
    for c in customers:
        prob += (pulp.lpSum(x[(c, s)] for s in sites) == 1)

    # A site should only have flow if used
#     for s in sites:
#         prob += (pulp.lpSum(x[(c, s)] for c in customers) <= 999999999999 * y[s])
    for c in customers:
        for s in sites:
            prob += (x[c,s] <= y[s])
    # number of sites constraint
    prob += (pulp.lpSum((y[s] for s in sites)) == n_dcs)
    solver = pulp.PULP_CBC_CMD()
    solver.tmpDir = '/opt/ibm/ILOG/CPLEX_Studio201/cplex/'
    prob.solve(solver)

    # decision = [bool(y[s].varValue) for s in sites]
    decision = [y[s].varValue for s in sites]
# np.array(decision)
    return decision

if __name__ == "__main__":
  # distances = [[1, 3, 4,],[5, 3, 1]]
  # param_gap = 0.001
  # y = _solve_pMp(distances, 2, decimals)
  # print(y)
  cus_file = sys.argv[1]
  czLats, czLongts, czDemands = dh._read_cz_info(cus_file)
  fac_file = sys.argv[2]
  facLats, facLongts = dh._read_fac_info(fac_file)
  dist_matx = dh._calc_dist_matx(czLats, czLongts, czDemands, facLats, facLongts)
  # print(dist_matx)
  decimals = 2 
  nb_dcs = int(sys.argv[3])
  y = _solve_pMp(dist_matx, nb_dcs, decimals)
  for idx, e in enumerate(y):
    if e == 1:
      print(idx, end= ', ')