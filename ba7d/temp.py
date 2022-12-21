import re
with open("ba7d.in","r") as f:
    context=f.read()
lines=context.split("\n")
n=int(lines[0])
dist_mat=[[0]*n for _ in range(n)]
i=0
for i,line in enumerate(lines[1:]):
    for j,num in enumerate(line.split()):
        dist_mat[i][j]=int(num)

class Cluster:
    def __init__(self,name):
        self.name=name
        self.age=0
        self.nodes=[name]
def calculate_clusters_distance(dist_mat,c1,c2):
    dis=0
    for i in c1.nodes:
        for j in c2.nodes:
            dis+=dist_mat[i][j]
    return dis/float(len(c1.nodes)*len(c2.nodes))
def find_closest(dist_mat,current_clusters,clusters):
    clus_dist_list=[]
    p=0
    tmp_dic={}
    for i in current_clusters:
        for j in current_clusters:
            if i==j: continue
            clus_dist=calculate_clusters_distance(dist_mat, clusters[i],clusters[j])
            clus_dist_list.append(clus_dist)
            tmp_dic[p]=(i,j)
            p+=1
    min_value=min(clus_dist_list)
    argmin=clus_dist_list.index(min_value)
    return tmp_dic[argmin]

def merge_clusters(dist_mat,clusters,i,j,name):
    C_new=Cluster(name)
    C_new.nodes=clusters[i].nodes+clusters[j].nodes
    distance=calculate_clusters_distance(dist_mat, clusters[i], clusters[j])
    age= distance/2
    C_new.age=age
    return C_new

def add_edge(edge_list,C_new,c1,c2):
    distance1=C_new.age-c1.age
    distance2=C_new.age-c2.age

    if C_new.name not in edge_list.keys():
        edge_list[C_new.name]=[]
    if c1.name not in edge_list.keys():
        edge_list[c1.name]=[]
    if c2.name not in edge_list.keys():
        edge_list[c2.name]=[]       

    edge_list[C_new.name].append((c1.name,distance1))
    edge_list[C_new.name].append((c2.name,distance2))
    # print(f"{C_new.name}->{c1.name}:{distance1}")
    # print(f"{C_new.name}->{c2.name}:{distance2}")
    edge_list[c1.name].append((C_new.name,distance1))
    edge_list[c2.name].append((C_new.name,distance2))
def update_dis_mat(C_new,clusters,dist_mat):
    distance_list=[]
    for cluster in clusters:
        distance=calculate_clusters_distance(dist_mat, C_new, cluster)
        distance_list.append(distance)
    print(distance_list)
    for i in range(len(dist_mat)):
        dist_mat[i].append(distance_list[i])
    distance_list.append(0)

    dist_mat.append(distance_list)


def UPGMA(dist_mat,n):
    clusters=[Cluster(i) for i in range(n)]
    current_clusters=set([i for i in range(n)])
    edge_list={}
    name=n
    while len(current_clusters)>1:
        i,j=find_closest(dist_mat,current_clusters,clusters)
        print(i, j)
        print(current_clusters)
        C_new=merge_clusters(dist_mat,clusters,i,j,name)
        name+=1
        add_edge(edge_list,C_new,clusters[i],clusters[j])

        current_clusters.remove(i)
        current_clusters.remove(j)
        current_clusters.add(C_new.name)

        clusters.append(C_new)

        update_dis_mat(C_new,clusters,dist_mat)
    return edge_list
edge_list=UPGMA(dist_mat, n)
for i in range(len(edge_list)):
    for j in range(len(edge_list[i])):
        print(f"{i}->{edge_list[i][j][0]}:{edge_list[i][j][1]:.3f}")