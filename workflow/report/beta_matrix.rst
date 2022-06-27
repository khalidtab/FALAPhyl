The dissimilarity matrix which was used for the analysis.

Some notes about some of the dissimilarity matrix options:

- PhILR: Compostionally-aware distance matrix. Uses ILR and a Ward-hierarchical tree that is automatically generated.
- bray: Bray-Curtis. Good for detecting unlying ecological gradients
- manhattan: Manhattan distance is a distance metric between two points in a N dimensional vector space measured as the sum of the lengths of the projections of the line segment between the points onto the coordinate axes.
	- Euclidean and Manhattan dissimilarities are not good in gradient separation without proper standardization but are still included for comparison and special needs.
- canberra
- clark
- kulczynski: good in detecting underlying ecological gradients
- Jaccard: good in detecting underlying ecological gradients.
- gower: Good for detecting unlying ecological gradients
- altGower
- morisita: Can only use count data (integers)
- horn: Moritisa-horn. Can handle any abundance data
- raup: Raup-Crick index. It is a probabilistic index based on presence/absence data. It is defined as $1 - prob(j)$, or based on the probability of observing at least $j$ species in shared in compared communities. The current function uses analytic result from hypergeometric distribution (phyper) to find the probabilities. This probability (and the index) is dependent on the number of species missing in both sites, and adding all-zero species to the data or removing missing species from the data will influence the index. The probability (and the index) may be almost zero or almost one for a wide range of parameter values. The index is nonmetric: two communities with no shared species may have a dissimilarity slightly below one, and two identical communities may have dissimilarity slightly above zero. The index uses equal occurrence probabilities for all species, but Raup and Crick originally suggested that sampling probabilities should be proportional to species frequencies (Chase et al. 2011). A simulation approach with unequal species sampling probabilities is implemented in raupcrick function following Chase et al. (2011). The index can be also used for transposed data to give a probabilistic dissimilarity index of species co-occurrence (identical to Veech 2013).
- binomial: Binomial index is derived from Binomial deviance under null hypothesis that the two compared communities are equal. It should be able to handle variable sample sizes. The index does not have a fixed upper limit, but can vary among sites with no shared species.
- chao: Chao index tries to take into account the number of unseen species pairs
- cao: Cao index or CYd index (Cao et al. 1997) was suggested as a minimally biased index for high beta diversity and variable sampling intensity. Cao index does not have a fixed upper limit, but can vary among sites with no shared species. The index is intended for count (integer) data, and it is undefined for zero abundances; these are replaced with arbitrary value $0.1$ following Cao et al. (1997). Cao et al. (1997) used $log10$, but the current function uses natural logarithms so that the values are approximately $2.30$ times higher than with 10-based logarithms.
- mahalanobis: a Euclidean matrix where columns are centred, have unit variance, and are uncorrelated. The index is not commonly used for community data, but it is sometimes used for environmental variables. The calculation is based on transforming data matrix and then using Euclidean distances following Mardia et al. (1979).

