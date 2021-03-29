import numpy as np
import pandas as pd
import sys
from scipy import stats
import matplotlib.pyplot as plt

#alpha=np.array([.6])
ROC=np.array([50,100])

######adjacency matrix and input modifications 

def adj_matrix2(mRNA_path,miRNA_path):
  miRNA = pd.read_csv(miRNA_path,delimiter=',',header=None)
  miRNA=np.array(miRNA)
  data = pd.read_csv(mRNA_path,delimiter=',',header=None)
  data=np.array(data)
  mRNA=data[1:,1:]
  feature = data[1:,0]
  sample = data[0,1:]

  data = pd.read_excel('bipartite_targetscan.xlsx',header=None)
  data=np.transpose(np.array(data))
  final_miRNA=data[1:,0]
  final_mRNA=data[0,1:]
  adj=data[1:,1:].astype(np.int8)

  adj[adj == 1] = -1
  adj[adj==0]=1
  print('data loaded')

  xy, x_ind, y_ind = np.intersect1d(miRNA[0,:],sample,return_indices=True)
  _, x_ind1, y_ind1 = np.intersect1d(miRNA[:,0],final_miRNA,return_indices=True)
  xy1, x_ind2, y_ind2 = np.intersect1d(feature.astype('<U15'),final_mRNA.astype('<U15'),return_indices=True)
  mRNA=mRNA[:,y_ind].astype(float)
  miRNA=miRNA[:,x_ind]
  mRNA=mRNA[x_ind2,:]
  adj=adj[:,y_ind2]
  miRNA=miRNA[x_ind1,:].astype(float)
  adj=adj[y_ind1,:]
  mRNA=np.exp2(mRNA)-.001
  miRNA=np.nan_to_num(miRNA)

  return adj,mRNA,miRNA,xy,xy1


######network propagation


def NP(mRNA,miRNA,adj,alpha,maxiter):
    
  print('Predicting protein expression for alpha=',alpha)
  miRNA1=miRNA
  mRNA1=mRNA
  for j in range(maxiter):
    #print(j)
    old_miRNA = miRNA1
    old_mRNA = mRNA1

    miRNA1 = alpha*(np.matmul(adj,old_mRNA)) + (1-alpha)*miRNA
    mRNA1 = alpha*(np.matmul(adj.transpose(),old_miRNA)) + (1-alpha)*mRNA

    if np.logical_and((np.amax(np.absolute(miRNA1 - old_miRNA)))<10**-11, (np.amax(np.absolute(mRNA1-old_mRNA)))<10**-11):
      break

  print('optimized') 
  np.savetxt('predicted_protein.csv',mRNA1, delimiter=',')
  return mRNA1
 


#correlation


def corr(protein,sample_name,mRNA_name,ROC,ground_truth):
  print('calculating correlations')
  gencode = pd.read_csv('gencode23.csv',delimiter=',',usecols=[0,10])
  gencode=np.array(gencode)
  spectral_count = pd.read_csv('spectral_count.tsv',delimiter='\t',header=None)
  spectral_count=np.array(spectral_count)

  gene=[]
  c=[]
  j=0
  gencode=gencode[gencode[:,0].argsort()].astype('<U15')

  for i,x in enumerate(mRNA_name): 
    a=np.searchsorted(gencode[:,0],x.astype('<U15'))
    if a<=np.size(gencode,0)-1 and gencode[a,0]==x:
      gene.append(gencode[a,1])
      j+=1
    else:
      c.append(i)

  protein=np.delete(protein,c,0)  


  gene=np.array(gene)
  uni_gene=np.unique(gene)

  gene_exp=[]

  for i,x in enumerate(uni_gene):
    b= np.where(gene == x)[0]
    gene_exp.append(np.sum(protein[b], axis=0))

  j=0     
  k=0
  d=[]
  gene_matched_spectral_count=[]
  for i,x in enumerate(uni_gene):   
    c= np.where(spectral_count[:,0] == x)[0]
    if c.size == 0:
      d.append(i)
      j+=1
      continue
    gene_matched_spectral_count.append(spectral_count[c,:][0])
    k+=1
  gene_exp=np.delete(gene_exp,d,0)        
    
  print(j,' genes deleted from gene_exp')
  print(k,' genes added to gene_matched_spectral_count')
  print('dimesnion of gene_matched_spectral_count is: ', np.shape(gene_matched_spectral_count))
  print('dimesnion of gene_exp is: ', np.shape(gene_exp))
        

  final_protein=[]
  d=[]
  for i,x in enumerate(spectral_count[0,:]):

    if 'Spectral Counts' in x:
      names=x.split(':')
      names = [x[0:10] for x in names]
      names= ['TCGA-'+str(x) for x in names]
      _,x_ind,_ = np.intersect1d(sample_name,names,return_indices=True)
      if x_ind.size==0:
        d.append(i)
        continue
      final_protein.append(np.sum(gene_exp[:,x_ind], axis=1)/np.size(x_ind))  
    else:
      d.append(i)
  gene_matched_spectral_count=np.delete(gene_matched_spectral_count,d,1)     
  final_protein=np.transpose(final_protein)

  final_protein=np.array(final_protein,dtype=float)
  gene_matched_spectral_count=np.array(gene_matched_spectral_count,dtype=float)
  final_protein=np.log2(final_protein+1)
  gene_matched_spectral_count=np.log2(gene_matched_spectral_count+1)
  #final_protein=np.nan_to_num(final_protein)

  pear_coef = [stats.pearsonr(final_protein[:, i], gene_matched_spectral_count[:, i])[0] for i in range(np.size(final_protein,1))]

  spearman_coef =[stats.spearmanr(final_protein[:, i], gene_matched_spectral_count[:, i])[0] for i in range(np.size(final_protein,1))]
  sorted_prediction=np.argsort(final_protein,0)
  sorted_prediction=sorted_prediction[::-1]
  sorted_ground=np.argsort(gene_matched_spectral_count,0)
  sorted_ground=sorted_ground[::-1]

  res=[]
  for x in ROC:
    for i in range(np.size(final_protein,1)):
      xy= np.intersect1d(sorted_prediction[0:x,i],sorted_ground[0:x,i])  
      res.append(np.size(xy))
  roc11=np.reshape(res,(-1, np.size(final_protein,1)))
  
  print('correlations calculated')

  return pear_coef,spearman_coef,roc11
  

if len(sys.argv)!=5:
  print('wrong number of inputs. please specify 5 parameters: code, mRNA, miRNA, ground truth, alpha')
  sys.exit()

mRNA_path = sys.argv[1]
miRNA_path = sys.argv[2]
ground_truth = sys.argv[3]
alpha = float(sys.argv[4])

##### processing mRNA, miRNA and interaction netowrk to find compatible matrices 
adj,mRNA,miRNA,sample_name,uni_mRNA_names=adj_matrix2(mRNA_path,miRNA_path)
C=np.sqrt(np.outer(np.sum(np.absolute(adj),0),np.sum(np.absolute(adj),1)))
adj=np.divide(adj,C.transpose())
print('adj normalized')

pearson=[]
spearman=[]
roc12=[]

##### calculating protein expression and correlation for different alphas
# correlation is for evaluation purpose only. The variable "protein" below can be 
# directly saved as predicted protein expression without running the code further. 

protein=NP(mRNA,miRNA,adj,alpha,1000)
pear_coef,spearman_coef,roc=corr(protein,sample_name,uni_mRNA_names,ROC,ground_truth)
pearson.append(pear_coef)
spearman.append(spearman_coef)
roc12.append(roc)

##### calculating correlation for baseline (mRNA and ground truth)
pear_coef,spearman_coef,roc =corr(mRNA,sample_name,uni_mRNA_names,ROC)
pearson.append(pear_coef)
spearman.append(spearman_coef)
roc12.append(roc)
roc12=[e for sl in roc12 for e in sl]


#np.savetxt('pearson_coef.csv',pearson, delimiter=',')
#np.savetxt('spearman_coef.csv',spearman, delimiter=',')
#np.savetxt('ROC.csv',roc12, delimiter=',',fmt='%1.1d')





