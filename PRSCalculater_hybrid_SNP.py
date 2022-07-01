import os,re
import argparse
from collections import Counter
def indexgene(df):
    # genome annotation 
    res={}
    for line in open(df,'r'):
        if line[0]=='#':continue
        l=line.split()
        if l[2] != 'gene':continue
        x=l[-1].split(';')
        for i in x:
            if 'Name=' in i:
                gene=re.findall("Name=(.*)",i)[0]
            elif 'locus_tag=' in i:
                index=re.findall("locus_tag=(.*)",i)[0]
        res[gene]=index
    return(res)
def CNV266(df,GCFgene):
    # genome CNV message
    a=0
    res={}
    out=''
    probij={}
    for line in open(df,'r'):
        l=line.strip().split(',')
        if a==0:
            a+=1
            cc=l[2:]
        else:
            res[l[1]]=[]
            for i in l[2:]:
                if float(i)<0.73:
                    res[l[1]].append('0')
                elif 0.73<=float(i)<=1.74:
                    res[l[1]].append('1')
                elif float(i)>1.74:
                    res[l[1]].append('2')
               # if 0.34<float(i)<=0.73:
               #     res[l[1]].append('0.5')
               # elif 0.73<float(i)<=1.2:
               #     res[l[1]].append('1')
               # elif float(i)<=0.34:
               #     res[l[1]].append('0')
               # elif 1.2<float(i)<=1.74:
               #     res[l[1]].append('2')
               # elif float(i)>1.74:
#                    res[l[1]].append('3')
    xx=[]
    for i in cc:
        if i in GCFgene:
            xx.append(GCFgene[i])
        else:
            xx.append(i)
    for n in res.values():
        for i,j in enumerate(xx):
            if j not in probij:probij[j]={}
            if n[i] not in probij[j]: 
                probij[j][n[i]]=1
            else:
                probij[j][n[i]]+=1
    return(res,xx,probij)

def risklocis(df):
    #risk loci for GWAS result
    a=0
    risk={}
    for line in open(df,'r'):
        l=line.split(',')
        if a==0:a+=1;continue
        if l[0] in risk:
            risk[l[0]][l[2]]=[l[1],l[7]]
        else:
            risk[l[0]]={l[2]:[l[1],l[7]]}
    return(risk)
def snpgeno(df,probwij):
    #genome SNP
    snp={}
    for line in open(df,'r'):
        if line[0] =='#':continue
        l=line.strip().split()
        j=0
        if float(l[-2])==0:#genetic heterogeneity
            j=1
        elif float(l[-2])==1:
            j=2
        chrom=str(int(l[0][-3:-1])-32)
        if chrom not in snp:
            snp[chrom]={l[1]:j}
        else:
            snp[chrom][l[1]]=j
        if chrom in probwij:
            if l[1] in probwij[chrom]: 
                probwij[chrom][l[1]].append(j)
            else:
                probwij[chrom][l[1]]=[j]
        else:
            probwij[chrom]={l[1]:[j]}
    return(snp,probwij)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A parser for calculating the individual polygenic risk score that accepts input from the user")
    parser.add_argument('--anno', type=str, required=True,help="GCF file for finding gene name")
 #   parser.add_argument('--cnv', type=str, required=True, help="CNV file of individuals")
    parser.add_argument('--riskloci', type=str, required=False, help="GWAS following sig loci")
    parser.add_argument('--snpdir', type=str,required=False,help="a directory of snpfiles ")
    parser.add_argument('--output', type=str, required=True, help="filename for output")
    args = parser.parse_args()
    trancnv,cnvloc,probcnv={},{},{}
    Gene=indexgene(args.anno)    
 #   trancnv,cnvloc,probcnv=CNV266(args.cnv,Gene)
    risklocii=risklocis(args.riskloci)
    snploc={}
    probwij={}
    for snpfile in os.listdir(args.snpdir):
        if snpfile.endswith('.snp'):
            print(snpfile)
            snp,probwij=snpgeno(args.snpdir+snpfile,probwij)
            snploc[snpfile[:-4]]=snp
    prob={}
    for i,j in probwij.items():
        prob[i]={}
        for m,n in j.items():
            prob[i][m]=dict(Counter(n))            
    a=0
    out='prs\tParent\n'
    FJ11=''
    for stra,snpll in snploc.items():
        a+=1
       
        prs=0
        if stra=='FJ11':print(stra)
        if stra in trancnv:
            for i,j in risklocii.items():

                for m,n in j.items():
                    if i != '17':
                        if m not in snpll[i]:
                            wij=1
                            geno=0
                        else:
                            geno=1#snpll[i][m]
                            if geno in prob[i][m]:wij=float(prob[i][m][geno])
                            else:wij=0
                        FJ11+=str(wij)+'\t'+str(float(n[-1]))+'\t'+str(float(geno))+'\t'+n[0]+'\n'
                        prs+=wij*float(n[-1])*float(geno)
                    else:
                        if n[0][2:] not in cnvloc:
                            geno=0
                            wij=1
                        else:
                            #geno=trancnv[stra][cnvloc.index(n[0][2:])]
                            gene=float(trancnv[stra][cnvloc.index(n[0][2:])])/2
                            if geno not in probcnv[n[0][2:]]:wij=0
                            else:wij=float(probcnv[n[0][2:]][geno])
                        FJ11+=str(wij)+'\t'+str(float(n[-1]))+'\t'+str(float(geno))+'\t'+n[0]+'\n'
                        prs+=wij*float(n[-1])*float(geno)
        else:
            for i,j in risklocii.items():
                for m,n in j.items():
                    if i != '17':
                        if m not in snpll[i]:
                            wij=1
                            geno=0
                        else:
                            geno=1#snpll[i][m]
                            if geno in prob[i][m]:wij=float(prob[i][m][geno])
                            else:wij=0
                        prs+=wij*float(n[-1])*float(geno)
        out+=str(prs)+'\t'+stra+'\n'


    open(args.output,'w').write(out)
    print(FJ11)
