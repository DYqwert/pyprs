import pandas as pd
import sys

phe=pd.read_csv(sys.argv[1],sep='\t')
prs=pd.read_csv(sys.argv[2],sep='\t')
line=pd.read_csv(sys.argv[3],sep=',')
F1=pd.read_csv(sys.argv[4],sep=',')
#print(phe)
res=pd.merge(prs,phe)
res=pd.merge(res,line)
res['gamete']=res['prs']/2
F1['PRS_F1']=pd.Series() 
for i in range(F1.shape[0]):
    Parent1=F1['Parent1'][i]
    Parent2=F1['Parent2'][i]
    if Parent1 in list(res['Parent']) and Parent2 in list(res['Parent']):
        p2=float(res[res['Parent']==Parent2].loc[:,'gamete'])
        p1=float(res[res['Parent']==Parent1].loc[:,'gamete'])
        PRS=p2+p1
        print(PRS)
        F1['PRS_F1'][i]=str(PRS)
    else:
        F1['PRS_F1'][i]='nan'

F1.to_csv(sys.argv[5])
'''
classinformation=res['Lineage'].unique()
print(res.corr('spearman'))
print(res.corr('pearson'))
print(res.corr('kendall'))
for item in classinformation:
    item_data=res[res['Lineage'].isin([item])]
    print(item,item_data)
    print(item_data.corr('spearman'))
    print(item_data.corr('pearson'))
    print(item_data.corr('kendall'))
'''
