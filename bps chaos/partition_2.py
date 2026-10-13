import math
from array import array
import ROOT

# define a simple palette: one color (blue) for all values
ncontours = 1
stops = array('d', [0.0, 1.0])          # normalized positions in palette
red   = array('d', [0.0, 0.0])
green = array('d', [0.0, 0.0])
blue  = array('d', [1.0, 1.0])          # full blue
ROOT.TColor.CreateGradientColorTable(2, stops, red, green, blue, ncontours)
ROOT.gStyle.SetNumberContours(ncontours)


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
    result.sort(key=lambda v: (sum(v), v), reverse=True)
    return result
'''

# Example usage
N = 30
parts = partitions(N)
print ("There are",len(parts),"partitions for N = ",N)
for p in parts:
    print(p)
    
def split(V):
    split_parts = []
    for k, N_k in enumerate(V):   
        length = k+1 
        # print ("length: ", length)
        # length-1 strands cannot be split
        if length == 1:
            continue
        if N_k == 0:
            continue
        # there are N_k * floor(k/2) possible splittings.
        for j in range(1,math.floor(length/2)+1): 
            # print (j)
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

def join(V):
    joint_parts = []
    for k, N_k in enumerate(V):   
        length = k+1 
        #print ("length: ", length)
        if length == 10:
            continue
        if N_k == 0:
            continue
        N_k_new = N_k - 1
        # loop over all integers j where length+j <= 10
        for j in range(k+1, 11-k):   
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
                


#Test 
"""
V = [0,2,0,0,0,0,0,1,0,0]
split_parts = split(V)
joint_parts = join(V)
    
print (V, " has the following splittings: ")
for p in split_parts:
    print(p)

print (V, " has the following joinings: ")
for p in joint_parts:
    print(p)
"""

# build a mapping from partition tuple -> position
part_to_index = {tuple(p): i for i, p in enumerate(parts)}
pN = len(parts)

M = [[0 for _ in range(pN)] for _ in range(pN)]

# loop over partitions
for row_idx, V in enumerate(parts):
    split_parts = split(V)
    joint_parts = join(V)
    comb_list = split_parts + joint_parts

    # create indicator vector
    for P in comb_list:
        col_idx = part_to_index[tuple(P)]
        M[row_idx] [col_idx] = 1    # or +=1 if counting multiplicities
    
    

hist = ROOT.TH2I("M", "", pN, 0.5, pN+0.5, pN, 0.5, pN+0.5)
hist.SetStats(0)
hist.GetZaxis().SetLabelSize(0)    # remove labels
hist.GetZaxis().SetTitle("")       # remove title
hist.GetZaxis().SetTickLength(0)   # remove ticks
hist.SetFillColor(ROOT.kBlue)

for row_idx, V in enumerate(parts):
    split_parts = split(V)
    joint_parts = join(V)
    comb_list = split_parts + joint_parts

    for P in comb_list:
        col_idx = part_to_index[tuple(P)]
        hist.SetBinContent(row_idx + 1, col_idx + 1, 1)
# draw
c = ROOT.TCanvas("c", "Partition Matrix", 800, 800)
hist.SetMinimum(0)
hist.SetMaximum(1)
hist.SetFillColor(ROOT.kBlue)
hist.Draw("COL")
c.SaveAs("partition_matrix_lexocographic_N30.pdf")
            
            