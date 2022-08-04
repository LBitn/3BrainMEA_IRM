module Pystuff

using PyCall

export clus_1D
export TwoGauss

function __init__()
    py"""

    import numpy as np
    from sklearn.cluster import MeanShift, estimate_bandwidth

    def clus(x,quant):
        X = np.array(list(zip(x,np.zeros(len(x)))))
        bandwidth = estimate_bandwidth(X, quantile=quant)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(X)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_

        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
        clusters = []
        for k in range(n_clusters_):
            my_members = labels == k
            print( "cluster {0}: {1}".format(k, X[my_members, 0]))
        return labels
    """
end


function __init__()
    py"""

    import numpy as np
    from sklearn.mixture import GaussianMixture

    def separar2Gauss( samples ):

        ok = [ ]
        mixture = GaussianMixture(n_components=2).fit(samples.reshape(-1, 1))
        means_hat = mixture.means_.flatten()
        weights_hat = mixture.weights_.flatten()
        sds_hat = np.sqrt(mixture.covariances_).flatten()

        ok.append( mixture.converged_ )
        ok.append( means_hat )
        ok.append( sds_hat )
        ok.append( weights_hat )

        print(mixture.converged_)
        print(means_hat)
        print(sds_hat)
        print(weights_hat)
        return ok
    """
end

TwoGauss( samples ) = py"separar2Gauss"( samples );
clus_1D( x, y ) = py"clus"( x, y );


end


