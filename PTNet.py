import numpy as np
import pandas as pd
import sys
from scipy import stats
import matplotlib.pyplot as plt


######adjacency matrix and input modifications 

def adj_matrix2(mRNA_path,miRNA_path,network):
  miRNA = pd.read_csv(miRNA_path,delimiter=',',low_memory=False,index_col=0)
  mRNA = pd.read_csv(mRNA_path,delimiter=',',low_memory=False,index_col=0)
  adj = pd.read_csv(network,index_col=0)
  adj=adj.astype(np.int8).T
  adj[adj == 1] = -1
  adj[adj==0]=1
  print('data loaded')

  xy, x_ind, y_ind = np.intersect1d(miRNA.columns,mRNA.columns,return_indices=True)
  _, x_ind1, y_ind1 = np.intersect1d(miRNA.index,adj.index,return_indices=True)
  xy1, x_ind2, y_ind2 = np.intersect1d(mRNA.index.astype('<U15'),adj.columns.astype('<U15'),return_indices=True)
  mRNA=mRNA.iloc[:,y_ind].values.astype(float)
  miRNA=miRNA.iloc[:,x_ind].values.astype(float)
  mRNA=mRNA[x_ind2,:]
  adj=adj.iloc[:,y_ind2].values
  miRNA=miRNA[x_ind1,:]
  adj=adj[y_ind1,:]
  mRNA=np.exp2(mRNA)-.001
  miRNA=np.nan_to_num(miRNA)
  
  return adj,mRNA,miRNA,xy,xy1


######network propagation


def NP(mRNA,miRNA,adj,alpha,uni_mRNA_names,sample_name,maxiter):
    
  print('Predicting protein expression for alpha=',alpha)
  miRNA1=miRNA
  mRNA1=mRNA
  for j in range(maxiter):
    old_miRNA = miRNA1
    old_mRNA = mRNA1

    miRNA1 = alpha*(np.matmul(adj,old_mRNA)) + (1-alpha)*miRNA
    mRNA1 = alpha*(np.matmul(adj.transpose(),old_miRNA)) + (1-alpha)*mRNA

    if np.logical_and((np.amax(np.absolute(miRNA1 - old_miRNA)))<10**-11, (np.amax(np.absolute(mRNA1-old_mRNA)))<10**-11):
      break

  print('optimized') 
  print(mRNA1.shape)
  print(uni_mRNA_names.shape)
  print(new.shape)
  mRNA1 = np.concatenate((uni_mRNA_names,mRNA1),axis=1)
  mRNA1 = np.concatenate((new,mRNA1),axis=0)
  np.savetxt('predicted_protein.csv',mRNA1, delimiter=',',fmt='%s')
  
  return 
 

  

if len(sys.argv)!=5:
  print('wrong number of inputs. please specify 5 parameters: code, mRNA, miRNA, network, alpha')
  sys.exit()

mRNA_path = sys.argv[1]
miRNA_path = sys.argv[2]
network = sys.argv[3]
alpha = float(sys.argv[4])

##### processing mRNA, miRNA and interaction netowrk to find compatible matrices 
adj,mRNA,miRNA,sample_name,uni_mRNA_names=adj_matrix2(mRNA_path,miRNA_path,network)
C=np.sqrt(np.outer(np.sum(np.absolute(adj),0),np.sum(np.absolute(adj),1)))
adj=np.divide(adj,C.transpose())
print('adj normalized')

##### calculating protein expression

uni_mRNA_names = np.expand_dims(uni_mRNA_names,1)
new = []
new.append('sample_name')
new.extend(sample_name)
new = np.expand_dims(new,0)
NP(mRNA,miRNA,adj,alpha,uni_mRNA_names,new,1000)






