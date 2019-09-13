does not answer at all!
implemented!

eps=1e-15


A 4 by (n−1) matrix Z is returned. At the i-th iteration, clusters with indices Z[i, 0] and Z[i, 1]
are combined to form cluster n+i. A cluster with an index less than n corresponds to one of the
n original observations. The distance between clusters Z[i, 0] and Z[i, 1] is given by Z[i, 2].
 The fourth value Z[i, 3] represents the number of original observations in the newly formed cluster

cut heights should be specified in cut_heights.txt file

look at https://docs.scipy.org/doc/scipy-0.9.0/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
inconsistency is a concept about cutting dendograms!

cut heights be better to be converted to percent