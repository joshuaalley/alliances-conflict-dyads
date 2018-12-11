# alliances-conflict-dyads

This project replicates three prior studies of how alliance participation affects the likelihood of militarized conflict. I assess the robustness of each paper's results to corrections for dyadic clustering of observations, using [cluster-robust variance estimation](https://arxiv.org/abs/1312.3398), and varying intercept models.

The three studies are:
1. Benson 2011: *Unpacking Alliances: Deterrent and Compellent Alliances and Their Relationship with Conflict, 1816--2000*
2. Johnson and Leeds 2011: *Defense Pacts: A Prescription for Peace?*
3. Kenwick and Vasquez 2017: *Defense Pacts and Deterrence: Caveat Emptor* 

Each paper has a corresponding directory with the data and R scripts. `VI STAN Model.stan` provides the varying intercepts model for all three replications- all three re-analysis scripts will compile this model. 
