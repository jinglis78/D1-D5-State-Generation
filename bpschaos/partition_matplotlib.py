import numpy as np
import matplotlib.pyplot as plt
import math
from array import array

'''
# lexicographic ordering
def partitions(N):
    result = []

    def helper(remaining, max_k, vec):
        if remaining == 0:
            result.append(vec[:])
            return
        for k in range(1, max_k+1):
            if k <= remaining:
                vec[k-1] += 1
                helper(remaining - k, k, vec)
                vec[k-1] -= 1

    helper(N, N, [0]*N)

    # sort partitions: descending lexicographic order
    result.sort(reverse=True)
    return result
'''

# get list of partitions for integer N
def partitions(N):
    result = []

    def helper(remaining, max_k, vec):
        if remaining == 0:
            result.append(vec[:])
            return
        for k in range(1, max_k+1):
            if k <= remaining:
                vec[k-1] += 1
                helper(remaining - k, k, vec)
                vec[k-1] -= 1

    helper(N, N, [0]*N)

    # sort: first by cycle count (descending), then lex order (descending)
    result.sort(key=lambda v: (sum(v)))
    return result


N = 30
parts = partitions(N)
print ("There are",len(parts),"partitions for N = ",N)
    
# emulate the splitting of a state
def split(V):
    split_parts = []
    for k, N_k in enumerate(V):   
        length = k+1 
        #print ("length: ", length)
        # length-1 strands cannot be split
        if length == 1:
            continue
        if N_k == 0:
            continue
        # there are N_k * floor(k/2) possible splittings.
        for j in range(1,math.floor(length/2)+1): 
            #print (j)
            # produce integer pairs (j, length-j) = (1, length-1), (2, length-2), (3, length-3), ...
            pair = (j, length-j)
            # splitting creates a new partition vector V_split with N_k' = N_k - 1, N_j' = N_k + 1,  N_(k-j)' = N_(k-j) + 1 
            V_split = V[:]  
            #print (V_split)
            # length-1 is in positon 0, so subtract by -1. 
            V_split[length-1] -= 1
            V_split[j-1] += 1
            V_split[length-j-1] += 1
            # add this vector to the split_parts
            split_parts.append(V_split)
            #print (V_split)
    return split_parts

# emulate the joining of two strands
def join(V,N):
    joint_parts = []
    for k, N_k in enumerate(V):   
        # length = k+1, N_k = number of strands of length = k+1
        length = k+1
        if length == N:
            continue
        if N_k == 0:
            continue
        N_k_new = N_k - 1
        # loop over all integers j where length+j <= N
        for j in range(k+1, N+1-k):   
            V_join = V[:]  
            if length == j:
                N_j = V[j-1] - 1
            else:
                N_j = V[j-1]
        
            # joining is possible 
            if  N_j != 0:
                V_join[length-1] -= 1
                V_join[j-1] -= 1
                V_join[length+j-1] += 1
                # add this vector to the split_parts
                joint_parts.append(V_join)
                #print (V_join)
    return joint_parts     

# build a mapping from partition tuple -> position
part_to_index = {tuple(p): i for i, p in enumerate(parts)}

print ("Converted partitions to tuples")

pN = len(parts)




"""# Test 
V = [0, 0, 3, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
split_parts = split(V)
joint_parts = join(V, N)
    
print (V, " has the following splittings: ")
for p in split_parts:
    print(p)

print (V, " has the following joinings: ")
for p in joint_parts:
    print(p)"""


rows = []
cols = []

print("Collecting nonzero entries...")

for row_idx, V in enumerate(parts):
    comb_list = split(V) + join(V, N)
    
    for P in comb_list:
        col_idx = part_to_index[tuple(P)]
        
        rows.append(len(parts) - (row_idx+1))
        cols.append(col_idx) 

print("Plotting...")

plt.figure(figsize=(8, 8))

#plt.scatter(cols, rows, s=0.2, color='blue')
plt.scatter(cols, rows, s=0.2, color='blue', rasterized=False) # change rasterized=True if hard to open, but lose detail in pdf

plt.xticks([])
plt.yticks([])

plt.title(f"Partition Matrix N={N}")

plt.savefig(f"partition_matrix_N{N}.pdf", bbox_inches='tight', dpi=200) 