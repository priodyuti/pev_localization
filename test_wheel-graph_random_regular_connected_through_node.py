'''
It simulates the optimized structure: two graph components (wheel and random regular) connected via a node.
'''

import networkx as nx                                              
import matplotlib.pyplot as plt 
import numpy as np
import scipy as sp
import codecs
import sys
from operator import itemgetter
from numpy import linalg as LA
from functions_base import eigval_bound, store_undirected_network


def get_ipr(A):
  e_val, e_vec = LA.eigh(A)
  #print e_val
  d = np.shape(e_vec)   
  n = d[0]   
  lam1 = e_val[n-1]  
  lam2 = e_val[n-2] 
  lam3 = e_val[n-3]   
  Ev1 = e_vec[:,n-1]
  Ev2 = e_vec[:,n-2]
            
  IPR1 = 0.0;
  IPR2 = 0.0;
  for i in range(n):
    IPR1 = IPR1 + pow(Ev1[i],4);
    IPR2 = IPR2 + pow(Ev2[i],4); 
   
  del e_val
  del e_vec
  del Ev1
  del Ev2
  return IPR1,IPR2,lam1,lam2,lam3

path = 'wheel_RR_graph_1.txt'
#  N = 520;
#  E = 2630; 37 n2: 986
r = 45 
s1 = 1937
s2 = 470
print 's1: ',s1
p = r/float(s1)
Flag = False
while not Flag:
  G2 = nx.random_regular_graph(r,s2)
  Flag = nx.is_connected(G2)
  print Flag

G1 = nx.wheel_graph(s1)
#G1 = nx.star_graph(s1-1)
#print 'Number of nodes in Wheel graph: ', len(nx.nodes(G1))
deg_seq = nx.degree(G1)
max_deg = 0
node_index = 0
for i in list(nx.degree(G1).keys()):
  if deg_seq[i] > max_deg:
     max_deg = deg_seq[i];
     node_index = i

print 'max deg: ',max_deg
print 'Node label with max deg: ',node_index


N = s1 + s2 + 1
M = np.zeros((N,N), dtype=int)
M1 = np.zeros((s1,s1), dtype=int)
M2 = np.zeros((s2,s2), dtype=int)

Edge_list1 = nx.edges(G1)
for edge in Edge_list1:
  x = edge[0]
  y = edge[1]
  M1[x][y] = 1
  M1[y][x] = 1

Edge_list2 = nx.edges(G2)
for edge in Edge_list2:
  x = edge[0]
  y = edge[1]
  M2[x][y] = 1
  M2[y][x] = 1


#print 'Number of nodes in random regular graph: %d' %s1
Edges1 = nx.number_of_edges(G1)
avg_deg1 = Edges1*2/float(s1)
print 'Number of edges in wheel graph: ', Edges1
print 'Average degree of wheel graph : %f' %avg_deg1
ipr1,ipr2,lam1,lam2,lam3 = get_ipr(M1)                
print 's1: %d ipr1: %0.4f ipr2: %0.4f lam1: %0.16f lam2: %0.4f lam3: %0.4f\n'%(s1,ipr1,ipr2,lam1,lam2,lam3)                                                                                                        


#print 'Number of nodes in wheel graph: %d' %s2
Edges2 = nx.number_of_edges(G2)
avg_deg2 = Edges2*2/float(s2)
print 'Number of edges in random regular graph: ', Edges2
print 'Average degree of random regular graph : %f' %avg_deg2
ipr1,ipr2,lam1,lam2,lam3 = get_ipr(M2)                
print 's2: %d ipr1: %0.4f ipr2: %0.4f lam1: %0.16f lam2: %0.4f lam3: %0.4f\n'%(s2,ipr1,ipr2,lam1,lam2,lam3)                                                                                                        

for i in range(s1):
 for j in range(s1):
    M[i][j] = M1[i][j]
#print(len(M2))
#print(len(M))

M[s1-1][s1] = 1
M[s1][s1-1] = 1 
M[s1][s1+1] = 1
M[s1+1][s1] = 1

i = s1 + 1
j = s1 + 1
m = 0
n = 0

while i<N and m<s2:
   j = s1+1
   n = 0
   while j<N and n<s2:
     M[i][j] = M2[m][n]
     #print i,j,m,n
     j = j + 1
     n = n + 1
   i = i + 1
   m = m + 1

M[0][4] = 0
M[4][0] = 0

M[2000][2010] = 1
M[2010][2000] = 1

gp = nx.from_numpy_matrix(M)
Edges = nx.number_of_edges(gp)
avg_deg = Edges*2/float(N)
print 'Number of nodes in wheel-random_regular_graph_combined: %d' %N
print 'Number of edges in wheel-random_regular_graph_combined: %d' %Edges
print 'Average degree of wheel-random_regular_graph_combined : %f' %avg_deg

ipr1,ipr2,lam1,lam2,lam3 = get_ipr(M)                
print 'N: %d \n ipr1: %0.4f \n ipr2: %0.4f \n lam1: %0.16f \n lam2: %0.16f \n lam3: %0.16f\n'%(N,ipr1,ipr2,lam1,lam2,lam3)                                                                                                       
print(nx.is_connected(gp))
#eigval_bound(M)
store_undirected_network(M, N, path)


'''def determine_roots(N,E):

  if E < (N-1)+ceil(((N+2)*sqrt(3*(N+2)))/9):
      print'Correct!!!'

  eps = 0.0 
  sigma = 1.0 - eps
  p = -2*sigma - 4
  q = 8*sigma + sigma^2 + 1 - N
  r = 2*E - 4*sigma^2

  a = (1/3) * (3*q - p*p)
  b = (1/27) * (2*p*p*p - 9*p*q + 27*r)

  T1 = b*b/4+a*a*a/27

  if T1==0.0 || T1>0.0
     A = nthroot(((-b/2) + sqrt(T1)),3)
     B = nthroot(((-b/2) - sqrt(T1)),3)
  else
     A = ((-b/2) + sqrt(T1))^(1/3)
     B = ((-b/2) - sqrt(T1))^(1/3)
  end

  y1 = A + B;
  y2 = -(1/2)*(A+B) - (i*sqrt(3)/2)*(A-B)
  #y3 = -(1/2)*(A+B) + (i*sqrt(3)/2)*(A-B);

  ##y2 = -(1/2)*(A+B) + ((sqrt(-1)*sqrt(3))/2)*(A-B)
  ##y3 = -(1/2)*(A+B) - ((sqrt(-1)*sqrt(3))/2)*(A-B)


  x1 = y1 - p/3;
  x2 = y2 - p/3;
  %x3 = y3 - p/3;
 
  print'x1: %0.2f x2: %0.2f\n',x1,x2

  n1 = (ceil(x1)-1)^2
  n2 = N-n1-1
  print'n1: %d n2: %d\n',n1,n2

  n1 = (ceil(x2)-1)^2
  n2 = N-n1-1
  print'n1: %d n2: %d\n',n1,n2
'''
