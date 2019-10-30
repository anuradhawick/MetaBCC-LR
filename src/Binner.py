import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from collections import defaultdict
import sys
from kneed import KneeLocator
from sklearn.neighbors import NearestNeighbors
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re
from random import sample
from sklearn.metrics.cluster import adjusted_rand_score

parser = argparse.ArgumentParser(description='Identify Initial Bins.')

parser.add_argument('-p3',
                    help="Path to trimer profiles",
                    type=str,
                    required=True)
parser.add_argument('-p15',
                    help="Path to 15-mer coverage histogram profiles",
                    type=str,
                    required=True)
parser.add_argument('-ids',
                    help="Read ids of reads (For dry runs with ground truth)",
                    type=str,
                    required=False,
                    default=None)
parser.add_argument('-o', help="Output directory", type=str, required=True)

args = parser.parse_args()

p3 = args.p3
p15 = args.p15
o = args.o
ids = args.ids

if not os.path.exists(o):
    os.makedirs(o)


class Read:
    def __init__(self, i, p3, p15):
        self.i = i
        self.p3 = p3
        self.p15 = p15
        self.embededValue = None


class Cluster:
    def __init__(self, name):
        self.name = name
        self.reads = []
        self.sampledReads = None

    def sample(self):
        self.sampledReads = sample(self.reads, 1000)

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
        return TSNE(n_components=2,
                    init='pca').fit_transform(self.getDataP15())

    def embedP3(self):
        return TSNE(n_components=2, init='pca').fit_transform(self.getDataP3())

    def embedP15sampled(self):
        return TSNE(n_components=2,
                    init='pca').fit_transform(self.getDataP15sampled())

    def embedP3sampled(self):
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


cluster_0 = Cluster("R")
all_species = []
# using coverage
with open(p15) as fp15, open(p3) as fp3:
    if ids:
        with open(ids) as fids:
            for p15line, p3line, idline in zip(fp15, fp3, fids):
                all_species.append(idline.strip())
                cluster_0.addRead(
                    Read(idline.strip(),
                         list(map(float,
                                  p3line.strip().split())),
                         list(map(float,
                                  p15line.strip().split()))))
        all_species = list(set(all_species))
    else:
        for p15line, p3line in zip(fp15, fp3):
            cluster_0.addRead(
                Read("read", list(map(float,
                                      p3line.strip().split())),
                     list(map(float,
                              p15line.strip().split()))))


def scatterPlot(X, Y, title="Untitled", labels=[]):
    fig = plt.figure()
    fig.suptitle(title, fontsize=20)
    if len(labels) == len(X):
        sns.scatterplot(X, Y, hue=labels, alpha=0.6).plot()
    else:
        sns.scatterplot(X, Y, alpha=0.6).plot()

    plt.savefig(o + "/" + title + ".eps", dpi=1200, format="eps")


def estimateEpsilon(X_embedded, plot=False):
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
                
    if plot:
        plt.figure()
        plt.plot(distances)
        kneedle.plot_knee_normalized()
    return 10 * distances[kneedle.knee]


def evaluateClusters(clusters, species):
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

    print("\n")
    print(all_clusters)
    print("\n")
    for x in range(len(mat)):
        print(all_species[x], " " * (30-len(all_species[x])), " ".join(list(map(lambda x: " " * (20- len(str(x))) + str(x),mat[x]))))

    print("Total reads = ", tot)
    print("\nPrecision ", rMax / tot)
    print("Recall ", cMax / tot)
    print("ARI ", adjusted_rand_score(truth, labels))


def clusterReads(cluster, t, sample=False, plot=False):
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
    
    clustering = DBSCAN(eps=estimateEpsilon(embeddedRoot)).fit(embeddedRoot)
    if t == "composition":
        labels = list(map(lambda x: cluster.name + "-c-" + str(x), list(clustering.labels_)))
    else:
        labels = list(map(lambda x: cluster.name + "-x-" + str(x), list(clustering.labels_)))

    if plot:
        scatterPlot(embeddedRoot[:, 0],
                    embeddedRoot[:, 1],
                    title="Separation by " + t  + " " + cluster.name,
                    labels=labels)
    label_set = list(set(labels))

    if ids:
        if sample and len(cluster.reads) > 2000:
            reads = cluster.sampledReads
        else:
            reads = cluster.reads
        labels1 = list(map(lambda r: r.i, reads))
        if plot:
            scatterPlot(embeddedRoot[:, 0],
                        embeddedRoot[:, 1],
                        title="Separation by " + t + " Truth " + cluster.name,
                        labels=labels1)
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
            new_clusters[cname] = Cluster(cname)
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
                new_clustersAllReads[bestC] = Cluster(bestC)

            new_clustersAllReads[bestC].addRead(r)
        return new_clustersAllReads
    return new_clusters

stats = open(o + "/cluster-stats.txt", "w+")

clx = clusterReads(cluster_0, "coverage", False, True)

# for each cov cluster cluster by composition
print("\n")
cly = {}

def clusterByComposition(cluster):
    r1 = clusterReads(cluster, "composition", True, True)
    result = {}

    if len(r1) > 1:
        for cn, c2 in r1.items():
            r2 = clusterByComposition(c2)
            result.update(r2)
    else:
        result.update(r1)
    
    return result

for cn, c in clx.items():
    t = clusterByComposition(c)
        
    cly.update(t)

    for cn2, cc in t.items():
        stats.write(cn2)
        stats.write("\n")
        stats.write(" ".join(list(map(str, c.getMeanP15()))))
        stats.write("\n")
        stats.write(" ".join(list(map(str, cc.getMeanP3()))))
        stats.write("\n")
        stats.write(" ".join(list(map(str, c.getStdP15()))))
        stats.write("\n")
        stats.write(" ".join(list(map(str, cc.getStdP3()))))
        stats.write("\n")

# discarding poor bins
for k in list(cly.keys()):
    if len(cly[k].reads) < 100:
        del cly[k]

stats.close()

if ids:
    evaluateClusters(cly, all_species)


