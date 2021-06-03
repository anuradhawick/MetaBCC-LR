import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from kneed import KneeLocator
from sklearn.neighbors import NearestNeighbors
from random import sample
from sklearn.metrics.cluster import adjusted_rand_score
import matplotlib
import logging
from tabulate import tabulate
import umap
from song.song import SONG

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger('MetaBCC-LR')

class Read:
    def __init__(self, i, p3, p15):
        self.i = i
        self.p3 = p3
        self.p15 = p15
        self.embededValue = None

class Cluster:
    def __init__(self, name, embedding):
        self.name = name
        self.reads = []
        self.sampledReads = None
        self.embedding = embedding

    def sample(self, count = 1000):
        self.sampledReads = sample(self.reads, count)

    def addRead(self, read):
        self.reads.append(read)

    def getDataP3(self):
        arr = []
        for read in self.reads:
            arr.append(read.p3)

        return np.array(arr)

    def getDataP15(self):
        arr = []
        for read in self.reads:
            arr.append(read.p15)

        return np.array(arr)

    def getDataP3sampled(self):
        arr = []
        for read in self.sampledReads:
            arr.append(read.p3)

        return np.array(arr)

    def getDataP15sampled(self):
        arr = []
        for read in self.sampledReads:
            arr.append(read.p15)

        return np.array(arr)

    def embedP15(self):
        if self.embedding == 'song':
            return SONG().fit_transform(self.getDataP15())
        if self.embedding == 'umap':
            return umap.UMAP().fit_transform(self.getDataP15())
        return TSNE(n_components=2,
                    init='pca').fit_transform(self.getDataP15())

    def embedP3(self):
        if self.embedding == 'song':
            return SONG().fit_transform(self.getDataP3())
        if self.embedding == 'umap':
            return umap.UMAP().fit_transform(self.getDataP3())
        return TSNE(n_components=2, init='pca').fit_transform(self.getDataP3())

    def embedP15sampled(self):
        if self.embedding == 'song':
            return SONG().fit_transform(self.getDataP15sampled())
        if self.embedding == 'umap':
            return umap.UMAP().fit_transform(self.getDataP15sampled())
        return TSNE(n_components=2,
                    init='pca').fit_transform(self.getDataP15sampled())

    def embedP3sampled(self):
        if self.embedding == 'song':
            return SONG().fit_transform(self.getDataP3sampled())
        if self.embedding == 'umap':
            return umap.UMAP().fit_transform(self.getDataP3sampled())
        return TSNE(n_components=2,
                    init='pca').fit_transform(self.getDataP3sampled())

    def getMeanP15(self):
        return np.mean(self.getDataP15(), axis=0)

    def getMeanP3(self):
        return np.mean(self.getDataP3(), axis=0)

    def getStdP15(self):
        return np.std(self.getDataP15(), axis=0)

    def getStdP3(self):
        return np.std(self.getDataP3(), axis=0)

    def getLabels(self):
        return list(map(lambda r: r.i, self.reads))

    def getProbabilityP3(self, read):
        x = np.array(read.p3)
        mu = self.getMeanP3()
        s = self.getStdP3()

        p = -0.5 * (((x - mu)/s) ** 2.0) - np.log(((2.0 * np.pi) ** 0.5) * s)

        return np.sum(p)

    def getProbabilityP15(self, read):
        x = np.array(read.p15)
        mu = self.getMeanP15()
        s = self.getStdP15()

        p = -0.5 * ((x - mu)/s ** 2.0) - np.log(((2.0 * np.pi) ** 0.5) * s)

        return np.sum(p)

def plot_cluster(X, Y, title, labels, output):
    fig = plt.figure(figsize=(10, 10))
    fig.suptitle(title, fontsize=20)

    if len(labels) == len(X):
        sns.scatterplot(x=X, y=Y, hue=labels).plot()
    else:
        sns.scatterplot(x=X, y=Y).plot()

    plt.savefig(f"{output}/images/{title}.png", dpi=200, format="png")
    plt.close()

def estimate_epsilon(X_embedded, sensitivity, plot_name, output):
    neigh = NearestNeighbors(n_neighbors=2)
    nbrs = neigh.fit(X_embedded)
    distances, indices = nbrs.kneighbors(X_embedded)
    distances = np.sort(distances, axis=0)
    distances = distances[int(distances.shape[0]/2):, 1]
    i = np.arange(distances.shape[0])

    # get the elbow
    kneedle = KneeLocator(i,
                        distances,
                        S=1.0,
                        curve='convex',
                        direction='increasing',
                        interp_method='polynomial')

    if plot_name:
        plt.figure()
        plt.plot(distances)
        kneedle.plot_knee_normalized()
        plt.savefig(f"{output}/images/{plot_name}.png", dpi=200, format="png")
        plt.close()

    return sensitivity * distances[kneedle.knee]

def evaluate_clusters(clusters, species):
    all_clusters = list(clusters.keys())
    all_species = species

    s_map = {}
    c_map = {}

    labels = []
    truth = []

    for s in all_species:
        s_map[s] = len(s_map)

    for c in all_clusters:
        c_map[c] = len(c_map)

    mat = [[0 for c in range(len(c_map))] for s in range(len(s_map))]

    for c in all_clusters:
        cc = clusters[c]
        for s in cc.getLabels():
            mat[s_map[s]][c_map[c]] += 1
            labels.append(c_map[c])
            truth.append(s_map[s])

    rMax = 0
    cMax = 0
    tot = 0

    for x in mat:
        tot += sum(x)
        rMax += max(x)

    matT = [[mat[j][i] for j in range(len(mat))] for i in range(len(mat[0]))]

    for x in matT:
        cMax += max(x)

    for n, s in enumerate(all_species):
        mat[n] = [s] + mat[n]

    logger.info("\n" + tabulate(mat, headers=[""] + [c for c in all_clusters], tablefmt="fancy_grid"))

    logger.info(f"Total reads  {tot:15}")
    logger.info(f"Precision    {rMax / tot:15.3f}")
    logger.info(f"Recall       {cMax / tot:15.3f}")
    logger.info(f"ARI          {adjusted_rand_score(truth, labels):15.3f}")

def cluster_reads(cluster, t, sensitivity, threads, output, embedding, ground_truth, sample=False, plot=False):
    if len(cluster.reads) < 200:
        return {cluster.name + "-N": cluster}
    elif len(cluster.reads) > 2000 and sample:
        cluster.sample()
        if t == "composition":
            embeddedRoot = cluster.embedP3sampled()
        else:
            embeddedRoot = cluster.embedP15sampled()
    else:
        if t == "composition":
            embeddedRoot = cluster.embedP3()
        else:
            embeddedRoot = cluster.embedP15()
    
    try:
        eps = estimate_epsilon(embeddedRoot, sensitivity, "Kneedle for " + t  + " " + cluster.name, output)
    except:
        logger.error("Unable to locate an epsilon for clustering. Using 0.5")
        eps = 0.5

    # eps cannot be zero
    if eps==0:
        eps = 1

    clustering = DBSCAN(eps=eps, n_jobs=threads).fit(embeddedRoot)
    
    if t == "composition":
        labels = list(map(lambda x: f"{cluster.name}-c-{x+1}", list(clustering.labels_)))
    else:
        labels = list(map(lambda x: f"{cluster.name}-x-{x+1}", list(clustering.labels_)))

    if plot:
        plot_cluster(embeddedRoot[:, 0],
                    embeddedRoot[:, 1],
                    "Separation by " + t  + " " + cluster.name,
                    labels, output)

    if ground_truth is not None:
        if sample and len(cluster.reads) > 2000:
            reads = cluster.sampledReads
        else:
            reads = cluster.reads
        labels1 = list(map(lambda r: r.i, reads))
        if plot:
            plot_cluster(embeddedRoot[:, 0],
                        embeddedRoot[:, 1],
                        "Separation by " + t + " Truth " + cluster.name,
                        labels1, output)

    if sample and len(cluster.reads) > 2000:
        reads = cluster.sampledReads
    else:
        reads = cluster.reads
    new_clusters = {}
    
    for cname, read in zip(labels, reads):
        # discard noise
        if "--1" in cname:
            continue
        if cname not in new_clusters:
            new_clusters[cname] = Cluster(cname, embedding)
        new_clusters[cname].addRead(read)

    # reassign all reads back
    if sample and len(cluster.reads) > 2000:
        new_clustersAllReads = {}
        reads = cluster.reads
        for r in reads:
            maxP = float("-inf")
            bestC = None
            for cn, c in new_clusters.items():
                if t == "composition":
                    p = c.getProbabilityP3(r)
                else:
                    p = c.getProbabilityP15(r)
                if p > maxP:
                    maxP = p
                    bestC = cn
            
            if bestC not in new_clustersAllReads:
                new_clustersAllReads[bestC] = Cluster(bestC, embedding)

            new_clustersAllReads[bestC].addRead(r)
        return new_clustersAllReads
    return new_clusters

def cluster_composition(cluster, sensitivity, threads, output, embedding, ground_truth):
    r1 = cluster_reads(cluster, "composition", sensitivity, threads, output, embedding, ground_truth, True, True)
    result = {}

    if len(r1) > 100:
        for _, c2 in r1.items():
            r2 = cluster_composition(c2, sensitivity, threads, output, embedding, ground_truth)
            result.update(r2)
    else:
        result.update(r1)
    
    return result

def run_binner(output, ground_truth, threads, sensitivity, embedding):
    sensitivity = 11 - sensitivity
    p3 = np.load(f"{output}/profiles/3mers_sampled.npy")
    p15 = np.load(f"{output}/profiles/15mers_sampled.npy")
    output_binning = f"{output}/misc/"
    
    if ground_truth is not None:
        ground_truth = np.load(f"{output}/misc/filtered_truth_sampled.npy")

    cluster_init = Cluster("Root", embedding)
    all_species = []

    if ground_truth is not None:
        for p3_row, p15_row, truth_read in zip(p3, p15, ground_truth):
            all_species.append(truth_read)
            cluster_init.addRead(Read(truth_read, p3_row, p15_row))
        all_species = list(set(all_species))
    else:
        for p3_row, p15_row in zip(p3, p15):
            cluster_init.addRead(Read(0, p3_row, p15_row))

    stats = open(f"{output_binning}/cluster-stats.txt", "w+")

    logger.debug("Clustering using coverage")
    coverage_based_clusters = cluster_reads(cluster_init, "coverage", sensitivity, threads, output, embedding, ground_truth, False, True)
    logger.debug(f"Identified number of coverage clusters - {len(coverage_based_clusters)}")

    # for each cov cluster cluster by composition
    logger.debug("Clustering using composition")
    final_clusters = {}
    bin_name = 1
    final_binning_result = {}

    if sensitivity < 8:
        logger.debug(f"Discarding small clusters (< 500 reads in the sampled set)")

    for _, coverage_based_cluster in coverage_based_clusters.items():
        composition_based_clusters = cluster_composition(coverage_based_cluster, sensitivity, threads, output, embedding, ground_truth)
            
        final_clusters.update(composition_based_clusters)

        for _, composition_based_cluster in composition_based_clusters.items():
            if len(composition_based_cluster.reads) > 500 and sensitivity < 8:
                stats.write(f"Bin-{bin_name}")
                final_binning_result[f"Bin-{bin_name}"] = composition_based_cluster
                stats.write("\n")
                stats.write(" ".join(list(map(str, coverage_based_cluster.getMeanP15()))))
                stats.write("\n")
                stats.write(" ".join(list(map(str, composition_based_cluster.getMeanP3()))))
                stats.write("\n")
                stats.write(" ".join(list(map(str, coverage_based_cluster.getStdP15()))))
                stats.write("\n")
                stats.write(" ".join(list(map(str, composition_based_cluster.getStdP3()))))
                stats.write("\n")
                bin_name += 1
    
    logger.debug(f"Identified number of coverage and composition clusters - {bin_name}")

    if sensitivity < 8:
        logger.debug(f"Discarding small clusters")
        # discarding poor bins
        for k in list(final_clusters.keys()):
            if len(final_clusters[k].reads) < 100:
                del final_clusters[k]

    stats.close()

    if ground_truth is not None:
        evaluate_clusters(final_binning_result, all_species)
