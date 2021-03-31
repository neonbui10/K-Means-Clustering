from Bio import SeqIO
import numpy
numpy.set_printoptions(threshold=numpy.inf)
from sklearn.manifold import MDS
from sklearn.cluster import KMeans as kmean
from matplotlib import pyplot as plt

#method for hamming distance
def hammingDistance(str1, str2):
    i=0
    count=0
    while(i < len(str1)):
        if(str1[i] != str2[i]):
            count+=1
        i+=1
    return count

#reading file and storing into array
array = []
for seq_record in SeqIO.parse("seq.fas", "fasta"):
    array.append(seq_record.seq)

size = len(array)-1

#finding hamming distance and store into list
hamList = list()
for a in range(0, size, 1):
   for b in range(0, size, 1):
        hamList.append(hammingDistance(array[a], array[b]))

#creating 2D matrix
arr = numpy.array(hamList).reshape(size, size)

#MDS
multiDimen = MDS(dissimilarity="precomputed")
data = multiDimen.fit_transform(arr)

#plotting graph
points = numpy.array(data)
a,b = points.T
plt.scatter(a,b, color='forestgreen')

#clustering data using k-means algorithm
kcenter = kmean(n_clusters=3, init='k-means++', max_iter=300, n_init=10, random_state=0)
pred_y = kcenter.fit_predict(points)
plt.scatter(points[:,0], points[:,1], c=kcenter.labels_.astype(float))
plt.scatter(kcenter.cluster_centers_[:,0], kcenter.cluster_centers_[:,1], s=50, color='red')
plt.show()
