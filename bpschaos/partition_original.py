import math
from array import array
import ROOT

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

print ("Initalising matrix")

M = [[0 for _ in range(pN)] for _ in range(pN)]

print ("Initiailised matrix M")

# loop over partitions and either split or join
for row_idx, V in enumerate(parts):
    #print (V) 
    split_parts = split(V)
    #print ("\tsplit: ", split_parts)
    joint_parts = join(V, N)
    #print ("\tjoin: ", joint_parts)
    comb_list = split_parts + joint_parts

    # create indicator vector
    for P in comb_list:
        col_idx = part_to_index[tuple(P)]
        M[row_idx] [col_idx] = 1    # or +=1 if counting multiplicities

print ("Looped over partitions and filled matrix")


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



# ___ Produce the plot ____
        
# set up colors and canvas
ncontours = 1
stops = array('d', [0.0, 1.0])          # normalized positions in palette
red   = array('d', [0.0, 0.0])
green = array('d', [0.0, 0.0])
blue  = array('d', [1.0, 1.0])          # full blue
ROOT.TColor.CreateGradientColorTable(2, stops, red, green, blue, ncontours)
ROOT.gStyle.SetNumberContours(ncontours)

# make histogram
hist = ROOT.TH2I("M", "", pN, 0.5, pN+0.5, pN, 0.5, pN+0.5)
hist.SetStats(0)
hist.GetZaxis().SetLabelSize(0)    # remove labels
hist.GetZaxis().SetTitle("")       # remove title
hist.GetZaxis().SetTickLength(0)   # remove ticks
hist.SetFillColor(ROOT.kBlue)

print ("About to fill histogram")

count = 0
for row_idx, V in enumerate(parts):
    sum = 0
    N_1 = 0
    for k,N_k in enumerate(V):
        if k == 0:
            N_1 = N_k
        sum += N_k


    """if sum == 8:
        count += 1
        if count > 600:"""
    #print (count, V)
    split_parts = split(V)
    joint_parts = join(V,N)
    comb_list = split_parts + joint_parts
    
    for P in comb_list:
        col_idx = part_to_index[tuple(P)]
        hist.SetBinContent(row_idx + 1, len(parts) - (col_idx + 1), 1)

print (count)

print ("Filled histogram")

# draw
c = ROOT.TCanvas("c", "Partition Matrix", 800, 800)
hist.SetMinimum(0)
hist.SetMaximum(1)
hist.SetFillColor(ROOT.kBlue)
hist.Draw("COL")
c.SaveAs(f"partition_matrix_N{N}.pdf")
            
            