{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "round-watts",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pulp\n",
    "import numpy as np"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "vanilla-syndicate",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Collecting pulp\n",
      "  Downloading PuLP-2.4-py3-none-any.whl (40.6 MB)\n",
      "Collecting amply>=0.1.2\n",
      "  Downloading amply-0.1.4-py3-none-any.whl (16 kB)\n",
      "Collecting docutils>=0.3\n",
      "  Downloading docutils-0.16-py2.py3-none-any.whl (548 kB)\n",
      "Requirement already satisfied: pyparsing in c:\\python38\\scg_64\\lib\\site-packages (from amply>=0.1.2->pulp) (2.4.7)\n",
      "Installing collected packages: docutils, amply, pulp\n",
      "Successfully installed amply-0.1.4 docutils-0.16 pulp-2.4\n"
     ]
    }
   ],
   "source": [
    "# !pip3 install pulp"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "altered-father",
   "metadata": {},
   "outputs": [],
   "source": [
    "def _solve_msd_only_problem(distances, n_dcs, solver, decimals):\n",
    "    distances = np.around(distances, decimals)\n",
    "    prob = pulp.LpProblem(\"MSDOptimization\", pulp.LpMinimize)\n",
    "\n",
    "    customers = list(range(distances.shape[0]))\n",
    "    sites = list(range(distances.shape[1]))\n",
    "\n",
    "    # The assignment to each DC\n",
    "    x = pulp.LpVariable.dicts('X',\n",
    "                              ((c, s) for s in sites for c in customers),\n",
    "                              lowBound=0,\n",
    "                              cat=pulp.LpContinuous)\n",
    "\n",
    "    # If the site is used\n",
    "    y = pulp.LpVariable.dicts('Y',\n",
    "                              (s for s in sites),\n",
    "                              0, 1, cat=pulp.LpBinary)\n",
    "\n",
    "    # probt = sum(x[(c, s)] * violation[c, s])\n",
    "    # min # of violations or min weighted vol of violations\n",
    "                                                                                                                          41,5          17%     # probt = sum(x[(c, s)] * violation[c, s]\n",
    "    # Objective function to minimize - the flow to any DC \n",
    "    prob += (pulp.lpSum(x[(c, s)] * distances[c, s] for s in sites for c in customers),\n",
    "             'Objective')\n",
    "\n",
    "    # Each customer's flows should add to 1\n",
    "    for c in customers:\n",
    "        prob += (pulp.lpSum(x[(c, s)] for s in sites) == 1)\n",
    "\n",
    "    # A site should only have flow if used\n",
    "#     for s in sites:\n",
    "#         prob += (pulp.lpSum(x[(c, s)] for c in customers) <= 999999999999 * y[s])\n",
    "    for c in customers:\n",
    "        for s in sites:\n",
    "            prob += (x[c,s] <= y[s])\n",
    "    # number of sites constraint\n",
    "    prob += (pulp.lpSum((y[s] for s in sites)) == n_dcs)\n",
    "\n",
    "    prob.solve(solver)\n",
    "\n",
    "    decision = [bool(y[s].varValue) for s in sites]\n",
    "    return np.array(decision)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
