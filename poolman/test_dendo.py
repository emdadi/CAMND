import matplotlib
# matplotlib.use('PDF')
import scipy as sp
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from numpy import array

if __name__ == '__main__':
    mat = array([[1, 0.5, 0.9],
                 [0.5, 1, -0.5],
                 [0.9, -0.5, 1]
    ])
    plt.subplot(1,2,1)
    plt.title("mat")
    dist_mat = mat
    linkage_matrix = sp.cluster.hierarchy.linkage(dist_mat,
                                                  "single")
    print "linkage2:"
    print linkage(1-dist_mat, "single")
    dendrogram(linkage_matrix,
               color_threshold=1,
               labels=["a", "b", "c"],
               show_leaf_counts=True)
    plt.subplot(1,2,2)
    plt.title("1 - mat")
    dist_mat = 1 - mat
    linkage_matrix = linkage(dist_mat,
                             "single")
    dendrogram(linkage_matrix,
               color_threshold=1,
               labels=["a", "b", "c"],
               show_leaf_counts=True)

    plt.show()