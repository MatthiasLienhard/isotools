from umap import UMAP # pylint: disable-msg=E0611
from sklearn.decomposition import PCA
from scipy.stats import binom,norm, chi2, betabinom, fisher_exact,beta # pylint: disable-msg=E0611
from scipy.special import gammaln, polygamma,gamma# pylint: disable-msg=E0611
from scipy.optimize import minimize
import statsmodels.stats.multitest as multi
import logging
import numpy as np
import pandas as pd
import itertools
from tqdm import tqdm
import warnings
from ._utils import overlap
logger=logging.getLogger('isotools')

# differential splicing
def proportion_test(x,n):
    # Normal approximation
    #x,n should be lenght 2(the two groups)
    #tests H0: proportions are equal vs H1: proportions are different (two sided)
    x=[xi.sum() for xi in x]
    n=[ni.sum() for ni in n]
    p1=[x[i]/n[i] for i in range(2)]
    p0=(x[0]+x[1])/(n[0]+n[1])
    z=abs(p1[0]-p1[1])/np.sqrt(p0*(1-p0)*(1/n[0]+1/n[1]))
    return(2*norm.sf(z)), (p1[0],0,p1[1],0,p0,0)#two sided alternative
    
def binom_lr_test(x,n):
    # likelihood ratio test
    # x,n should be length 2 (the two groups)
    # principle: log likelihood ratio of M0/M1 is chi2 distributed
    x=[xi.sum() for xi in x]
    n=[ni.sum() for ni in n]
    p1=[x[i]/n[i] for i in range(2)]
    p0=(x[0]+x[1])/(n[0]+n[1])
    # calculate the log likelihoods
    l0 = binom.logpmf(x, n, p0).sum() 
    l1 = binom.logpmf(x,n, p1).sum()
    # calculate the pvalue (sf=1-csf(), 1df)
    return chi2.sf(2*(l1-l0),1),(p1[0],0,p1[1],0,p0,0)

def loglike_betabinom(params, k,n):
    '''returns  log likelihood of betabinomial and its partial derivatives'''
    a, b = params
    logpdf = gammaln(n+1) + gammaln(k+a) + gammaln(n-k+b) + gammaln(a+b) - \
     (gammaln(k+1) + gammaln(n-k+1) + gammaln(a) + gammaln(b) + gammaln(n+a+b))
    e=polygamma(0,a+b) -  polygamma(0,n+a+b)
    da= e+polygamma(0,k+a)    - polygamma(0,a) 
    db= e+polygamma(0,n-k+b)  - polygamma(0,b) 
    return -np.sum(logpdf), np.array((-np.sum(da),-np.sum(db)))

def betabinom_lr_test(x,n):
    ''' likelihood ratio test with random-effects betabinomial model
     x,n should be 2 sets (the two groups)
     x is betabinomial(n,a,b), eg a binomial distribution, where p follows beta ditribution with parameters a,b>0
     mean m=a/(a+b) overdispersion d=ab/((a+b+1)(a+b)^2) --> a=-m(m^2-m+d)/d b=(m-1)(m^2-m+d)/d
     principle: log likelihood ratio of M0/M1 is chi2 distributed
     '''    
    params=list()
    success=True
    if any(ni.sum()==0 for ni in n):
        return (np.nan,[None, None]) #one group is not covered at all - no test possible. Checking this to avoid RuntimeWarnings (Mean of empty slice)
    x_all, n_all=(np.concatenate(x),np.concatenate(n))
    for xi,ni in itertools.chain(zip(x,n),((x_all,n_all),)):
        xi, ni=xi[ni>0], ni[ni>0] #avoid div by 0
        #x and n must be np arrays
        
        prob=xi/ni    
        m=prob.mean() #estimate initial parameters
        d=prob.var()
        if d==0: #just one sample? or all exactly the same proportion
            params.append((m,None)) # in this case the betabinomial reduces to the binomial
        else:
            d=max(d,1e-6)# to avoid division by 0
            e=(m**2-m+d) #helper          
            #find ml estimates for a and b
            mle = minimize(loglike_betabinom, x0=[-m*e/d,((m-1)*e)/d],bounds=((1e-6,None),(1e-6,None)),  args=(xi,ni),options={'maxiter': 250}, method='L-BFGS-B', jac=True)
            params.append(mle.x)                
            #mle = minimize(loglike_betabinom2, x0=[-d/(m*e),d/((m-1)*e)],bounds=((1e-9,None),(1e-9,None)),  args=(xi,ni),options={'maxiter': 250}, method='L-BFGS-B', tol=1e-6)
            #params.append([1/p for p in mle.x])                
            if not mle.success:
                logger.debug(f'no convergence in betabinomial fit: k={xi}\nn={ni}\nparams={params}\nmessage={mle.message}') #should not happen to often, mainly with mu close to boundaries
                success=False #prevent calculation of p-values based on non optimal parameters
    # calculate the log likelihoods
    params_alt=[(a/(a+b),  a*b/((a+b)**2*(a+b+1))) if b is not None else (a,None) for a,b in params ] #get alternative parametrization (mu and disp)
    if not success:
        return np.nan, [v for p in params_alt for v in p]
    try:
        l0 = betabinom_ll(x_all,n_all, *params[2]).sum() 
        l1 = betabinom_ll(x[0],n[0], *params[0]).sum()+betabinom_ll(x[1],n[1], *params[1]).sum()
    except (ValueError, TypeError):
        logger.critical(f'betabinom error: x={x}\nn={n}\nparams={params}')#should not happen
        raise
    return chi2.sf(2*(l1-l0),2), [v for p in params_alt for v in p] #note that we need two degrees of freedom here as h0 hsa two parameters, h1 has 4

def betabinom_ll(x,n,a,b):
    if b is None:
        return binom.logpmf(x, n, a).sum() 
    else:
        return betabinom.logpmf(x,n,a,b).sum() 


TESTS={ 'betabinom_lr':betabinom_lr_test,
        'binom_lr':binom_lr_test,
        'proportions': proportion_test}


def altsplice_test(self,groups, min_total=100,min_alt_fraction=.1, min_n=10, min_sa=.51, test='auto',padj_method='fdr_bh'):
    #min_total: minimum total coverage over splice bubble for both groups combined
    #min_alt_fraction: each alternative must contribute at least this fraction in both groups combined
    #min_n: for each group min_sa % of the samples must have coverage> min_n over splice_bubble
    
    assert len(groups) == 2 , "length of groups should be 2, but found %i" % len(groups)
    #find groups and sample indices
    if isinstance(groups, dict):
        groupnames=list(groups)
        groups=list(groups.values())
    elif all (isinstance(gn,str) and gn in self.groups() for gn in groups):
        groupnames=list(groups)
        groups=[self.groups()[gn] for gn in groupnames]
    elif all( isinstance(grp,list) for grp in groups):
        groupnames=['group1','group2']
    else:
        raise ValueError('groups not found in dataset')
    notfound=[sa for grp in groups for sa in grp if sa not in self.samples]
    if notfound:
        raise ValueError(f"Cannot find the following samples: {notfound}")    

    if isinstance(test,str):
        if test=='auto':
            test='betabinom_lr' if min(len(g) for g in groups)>1 else 'proportions'
        test_name=test
        try:
            test=TESTS[test]
        except KeyError:
            raise ValueError('test must be one of %s', str(list(TESTS)))
    else:
        test_name='custom'
    
    logger.info('testing differential splicing for %s using %s test',' vs '.join(f'{groupnames[i]} ({len(groups[i])})' for i in range(2)) ,test_name)
    sa_idx={sa:idx[0] for sa,idx in self._get_sample_idx().items()}
    grp_idx=[[sa_idx[sa] for sa in grp] for grp in groups]
    sidx=grp_idx[0]+grp_idx[1]
    if min_sa<1:
        min_sa*=sum(len(gr) for gr in groups)
    res=[]
    for g in tqdm(self):
        if g.coverage[sidx,:].sum()<min_total:
            continue
        known={} #check for known events
        if g.is_annotated and g.n_transcripts: 
            sg=g.ref_splice_graph
            for _,_,nX,nY, splice_type in sg.find_splice_bubbles():#find annotated alternatives (known)
                if splice_type in ("TSS","PAS"):
                    if (splice_type=="TSS")==(g.strand=="+"):
                        known.setdefault(splice_type,set()).add((sg[nY].end))
                    else:
                        known.setdefault(splice_type,set()).add((sg[nX].start))
                else:
                    known.setdefault(splice_type,set()).add((sg[nX].end,sg[nY].start))
        sg=g.splice_graph
        for setA,setB,nX, nY, splice_type in sg.find_splice_bubbles():
            
            junction_cov=g.coverage[:,setB].sum(1)
            total_cov=g.coverage[:,setA].sum(1)+junction_cov
            if total_cov[sidx].sum() < min_total:
                continue
            alt_fraction=junction_cov[sidx].sum()/total_cov[sidx].sum() 
            if alt_fraction<min_alt_fraction or alt_fraction > 1-min_alt_fraction:
                continue
            x=[junction_cov[grp] for grp in grp_idx]
            n=[total_cov[grp] for grp in grp_idx]
            if sum((ni>=min_n).sum() for ni in n)<min_sa:
                continue
            pval, params=test(x,n)
            if splice_type in ['TSS', 'PAS']:
                start, end=sg[nX].start, sg[nY].end
                if (splice_type=="TSS")==(g.strand=="+"):
                    novel=end not in known.get(splice_type,set())
                else:
                    novel=start not in known.get(splice_type,set())
            else:
                start, end=sg[nX].end, sg[nY].start
                novel=(start, end) not in known.get(splice_type,set())
            res.append(tuple(itertools.chain((g.name,g.id,g.chrom,g.strand, start, end,splice_type,novel,pval),params ,
                (val for lists in zip(x,n) for pair in zip(*lists) for val in pair ))))
        
    df=pd.DataFrame(res, columns= (['gene','gene_id','chrom','strand', 'start', 'end','splice_type','novel','pvalue']+ 
            [gn+part for gn in groupnames+['total'] for part in ['_PSI', '_disp'] ]+  
            [f'{sa}_{gn}_{w}' for gn,grp in zip(groupnames, groups) for sa in grp for w in ['in_cov', 'total_cov'] ]))
    try:
        mask = np.isfinite(df['pvalue'])
        padj = np.empty(mask.shape)
        padj.fill(np.nan) 
        padj[mask] = multi.multipletests(df.loc[mask,'pvalue'],method=padj_method)[1]
        df.insert(8,'padj',padj)
    except TypeError as e: #apparently this happens if df is empty...
        logger.error('unexpected error during calculation of adjusted p-values: {e}' ,e=e)
    return df

def splice_dependence_test(self,samples=None, min_cov=20,padj_method='fdr_bh',region=None):
    logger.warning('splice_dependence_test is depreciated - revise based on splice bubbles')
    if isinstance(samples, str):
        samples=[samples]
    elif samples is None:
        samples=self.samples
    sa_dict=self._get_sample_idx()
    sa_idx=[sa_dict[sa][0] for sa in samples]
    res=[]
    for g in tqdm(self.iter_genes(region=region)):
        for ma,j1,j2 in g.splice_graph.splice_dependence(sa_idx, min_cov):
            oddsratio, pval = fisher_exact(ma)    
            res.append((g.name,g.id,g.chrom, j1, j2,pval,oddsratio,ma))
    df=pd.DataFrame(res, columns= (['gene','gene_id','chrom', 'junction1', 'junction2','pvalue', 'oddsratio','counts']))
    try:
        mask = np.isfinite(df['pvalue'])
        padj = np.empty(mask.shape)
        padj.fill(np.nan) 
        padj[mask] = multi.multipletests(df.loc[mask,'pvalue'],method=padj_method)[1]
        df.insert(5,'padj',padj)
    except TypeError as e: #apparently this happens if df is empty...
        logger.error('unexpected error during calculation of adjusted p-values: {e}' ,e=e)
    return df            
        
def find_splice_bubbles(self, min_total=100, min_alt_fraction=.1,samples=None, region=None,include=None, remove=None):
    bubbles=[]
    if samples is None:
        samples=self.samples
    assert all(s in self.samples for s in samples), 'not all specified samples found'
    sa_dict={sa:i for i,sa in enumerate(self.samples)}
    sidx=np.array([sa_dict[sa] for sa in samples])

    assert 0<min_alt_fraction<.5, 'min_alt_fraction must be > 0 and < 0.5'
    for g in self.iter_genes(region, include, remove):
        if g.coverage[sidx,:].sum()<min_total:
            continue
        known={} #check for known events
        if g.is_annotated and g.n_transcripts: 
            sg=g.ref_splice_graph
            for _,_,nX,nY, splice_type in sg.find_splice_bubbles():#find annotated alternatives (known)
                if splice_type in ("TSS","PAS"):
                    if (splice_type=="TSS")==(g.strand=="+"):
                        known.setdefault(splice_type,set()).add((sg[nY].end))
                    else:
                        known.setdefault(splice_type,set()).add((sg[nX].start))
                else:
                    known.setdefault(splice_type,set()).add((sg[nX].end,sg[nY].start))
        sg=g.splice_graph
        for setA,setB,nX, nY, splice_type in sg.find_splice_bubbles():
            junction_cov=g.coverage[np.ix_(sidx,setA)].sum(1)
            total_cov=g.coverage[np.ix_(sidx,setB)].sum(1)+junction_cov
            if total_cov.sum()>=min_total and min_alt_fraction < junction_cov.sum()/total_cov.sum() < 1-min_alt_fraction:
                if splice_type in ['TSS', 'PAS']:
                    start, end=sg[nX].start, sg[nY].end
                    if (splice_type=="TSS")==(g.strand=="+"):
                        novel=end not in known.get(splice_type,set())
                    else:
                        novel=start not in known.get(splice_type,set())
                else:
                    start, end=sg[nX].end, sg[nY].start
                    novel=(start, end) not in known.get(splice_type,set())
                bubbles.append([g.id, g.chrom, start, end, splice_type, novel]+list(junction_cov)+list(total_cov))
    return pd.DataFrame(bubbles,columns=['gene', 'chr', 'start', 'end', 'splice_type', 'novel']+[f'{sa}_{what}' for what in ['in_cov','total_cov'] for sa in samples ] )

# summary tables (can be used as input to plot_bar / plot_dist)
def altsplice_stats(self, groups=None , weight_by_coverage=True, min_coverage=2, tr_filter={}):    
    weights=dict()
    #if groups is not None:
    #    gi={r:i for i,r in enumerate(runs)}
    #    groups={gn:[gi[r] for r in gr] for gn,gr in groups.items()}
    current=None
    if groups is not None:
        sidx={sa:i for i,sa in enumerate(self.samples)} #idx
        groups={gn:[sidx[sa] for sa in gr] for gn,gr in groups.items()}

    for g,trid,tr in self.iter_transcripts(**tr_filter):
        if g!=current:
            current=g
            w=g.coverage.copy() if groups is None else np.array([g.coverage[grp,:].sum(0) for grp in groups.values()])
            w[w<min_coverage]=0
            if not weight_by_coverage:
                w[w>0]=1
        if 'annotation' not in tr or tr['annotation'] is None:
            weights['unknown']=weights.get('unknown',np.zeros(w.shape[0]))+w[:,trid]
        else:
            for stype in tr['annotation'][1]:
                weights[stype]=weights.get(stype,np.zeros(w.shape[0]))+w[:,trid]
        weights['total']=weights.get('total',np.zeros(w.shape[0]))+w[:,trid]

    df=pd.DataFrame(weights, index=self.samples if groups is None else groups.keys()).T
    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0) #sort by row mean
    if weight_by_coverage:
        title='Expressed Transcripts'
        ylab='fraction of reads'
    else:
        title='Different Transcripts'
        ylab='fraction of  different transcripts'
        if min_coverage>1:
            title+=f' > {min_coverage} reads'

    return df, {'ylabel':ylab,'title':title}
    #
    
def filter_stats(self, groups=None, weight_by_coverage=True, min_coverage=2,consider=None, tr_filter={}):   

    weights=dict()
    if groups is not None:
        sidx={sa:i for i,sa in enumerate(self.samples)} #idx
        groups={gn:[sidx[sa] for sa in gr] for gn,gr in groups.items()}
    current=None
    for g,trid,tr in self.iter_transcripts(**tr_filter):
        if g!=current:
            current=g
            w=g.coverage.copy() if groups is None else np.array([g.coverage[grp,:].sum(0) for grp in groups.values()])
            w[w<min_coverage]=0
            if not weight_by_coverage:
                w[w>0]=1
        relevant_filter=[f for f in tr['filter'] if  consider is None or f in consider]
        for f in relevant_filter:
            weights[f]=weights.get(f,np.zeros(w.shape[0]))+w[:,trid]
        if not relevant_filter:
            weights['PASS']=weights.get('PASS',np.zeros(w.shape[0]))+w[:,trid]
        weights['total']=weights.get('total',np.zeros(w.shape[0]))+w[:,trid]
            
    df=pd.DataFrame(weights, index=self.samples if groups is None else groups.keys()).T

    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)
    ylab='fraction of reads' if weight_by_coverage else 'fraction of different transcripts'
    if weight_by_coverage:
        title='Expressed Transcripts'
    else:
        title='Different Transcripts'
    if min_coverage>1:
        title+=f' > {min_coverage} reads'
    return df, {'ylabel':ylab,'title':title}

def transcript_length_hist(self=None,  groups=None, add_reference=False,bins=50,x_range=(0,10000),weight_by_coverage=True,min_coverage=2,use_alignment=True, tr_filter={}, ref_filter={}):
    trlen=[]
    cov=[]
    current=None
    for g,trid,tr in self.iter_transcripts(**tr_filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
        trlen.append(sum(e[1]-e[0] for e in tr['exons']) if use_alignment else tr['source_len']) #source_len is not set in the current version
    cov=pd.DataFrame(cov,columns=self.samples)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0]-.5,x_range[1]-.5,bins+1)
    cov[cov<min_coverage]=0
    if not weight_by_coverage:    
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(trlen, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if add_reference:
        ref_len=[sum(e[1]-e[0] for e in tr['exons']) for _,_,tr in self.iter_ref_transcripts(**ref_filter)]
        counts['reference']=np.histogram(ref_len, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='linear', title='transcript length',xlabel='transcript length [bp]', density=True)
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def transcript_coverage_hist(self,  groups=None,bins=50,x_range=(1,1001), tr_filter={}):
    # get the transcript coverage in bins for groups
    # return count dataframe and suggested default parameters for plot_distr
    cov=[]
    current=None
    for g,trid,_ in self.iter_transcripts(**tr_filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
    cov=pd.DataFrame(cov, columns=self.samples)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0]-.5,x_range[1]-.5,bins+1)
    counts=pd.DataFrame({gn:np.histogram(g_cov, bins=bins)[0] for gn,g_cov in cov.items()})
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='log', title='transcript coverage',xlabel='reads per transcript')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    #plot histogram
    # cov.mask(cov.lt(x_range[0]) | cov.gt(x_range[1])).plot.hist(ax=ax, alpha=0.5, bins=n_bins)
    # ax=counts.plot.bar() 
    # ax.plot(x, counts)
   
def transcripts_per_gene_hist(self,   groups=None, add_reference=False, bins=49,x_range=(1,50), min_coverage=2, tr_filter={}, ref_filter={}):
    ntr=[]
    current=None
    if groups is None:
        group_names=self.samples        
    else:
        group_names=groups.keys()
        sidx={sa:i for i,sa in enumerate(self.samples)} #idx
        groups={gn:[sidx[sa] for sa in gr] for gn,gr in groups.items()}
    n_sa=len(group_names)
    for g,trid,_ in self.iter_transcripts(**tr_filter):
        if g!=current:
            current=g
            current_cov=g.coverage if groups is None else np.array([g.coverage[grp,:].sum(0) for grp in groups.values()])
            ntr.append(np.zeros(n_sa))
        ntr[-1]+= current_cov[:,trid]>= min_coverage
        
    ntr=pd.DataFrame((n for n in ntr if n.sum()>0), columns=group_names)    
    if isinstance(bins,int):
        bins=np.linspace(x_range[0]-.5,x_range[1]-.5,bins+1)
    counts=pd.DataFrame({gn:np.histogram(n, bins=bins)[0] for gn,n in ntr.items()})   
    if add_reference:
        if ref_filter:
            logger.warning('reference filter not implemented')
        ref_ntr=[g.n_ref_transcripts for g in self] #todo: add reference filter
        counts['reference']=np.histogram(ref_ntr, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    sub=f'counting transcripts covered by >= {min_coverage} reads'
    if 'include' in tr_filter:
        sub+=f', only including {", ".join(tr_filter["include"])}'
    if 'remove' in tr_filter:
        sub+=f', excluding {", ".join(tr_filter["remove"])}'
    params=dict(yscale='log',title='transcript per gene\n'+sub,xlabel='transcript per gene')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    
def exons_per_transcript_hist(self,  groups=None, add_reference=False, bins=34,x_range=(1,69),weight_by_coverage=True,  min_coverage=2, tr_filter={}, ref_filter={}):

    n_exons=[]
    cov=[]
    current=None
    for g,trid,tr in self.iter_transcripts(**tr_filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
        n_exons.append(len(tr['exons']))
    cov=pd.DataFrame(cov, columns=self.samples) 
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0]-.5,x_range[1]-.5,bins+1)
    cov[cov<min_coverage]=0
    if not weight_by_coverage:        
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(n_exons, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if add_reference:
        ref_n_exons=[len(tr['exons']) for _,_,tr in self.iter_ref_transcripts(**ref_filter)]
        counts['reference']=np.histogram(ref_n_exons, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    sub=f'counting transcripts covered by >= {min_coverage} reads'
    if 'include' in tr_filter:
        sub+=f', only including {", ".join(tr_filter["include"])}'
    if 'remove' in tr_filter:
        sub+=f', excluding {", ".join(tr_filter["remove"])}'
    params=dict(yscale='log', title='exons per transcript\n'+sub,xlabel='number of exons per transcript')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def downstream_a_hist(self, groups=None,add_reference=False,bins=30,x_range=(0,1),weight_by_coverage=True,  min_coverage=2, tr_filter={}, ref_filter={}):
    acontent=[]
    cov=[]
    current=None
    for g,trid,tr in self.iter_transcripts(**tr_filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
        try:
            acontent.append(tr['downstream_A_content'])
        except KeyError:
            acontent.append(-1)
    cov=pd.DataFrame(cov, columns=self.samples)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins+1)
    cov[cov<min_coverage]=0
    if not weight_by_coverage:        
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(acontent, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if add_reference:
        ref_acontent=[tr['downstream_A_content'] for _,_,tr in self.iter_ref_transcripts(**ref_filter) if 'downstream_A_content' in tr]
        counts['reference']=np.histogram(ref_acontent, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict( title='downstream genomic A content',xlabel='fraction of A downstream the transcript')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def direct_repeat_hist(self, groups=None, bins=10, x_range=(0,10), weight_by_coverage=True, min_coverage=2, tr_filter={}):
    # find the direct repeat length distribution in FSM transcripts and putative RTTS
    # putative RTTS are identified by introns where both splice sites are novel but within annotated exons
    novel_rl=[]
    known_rl=[]
    for g,trid,tr in self.iter_transcripts():
        if 'annotation' in tr and 'novel exonic splice donor' in tr['annotation'][1] and 'novel exonic splice acceptor' in tr['annotation'][1]:
            novel1, novel2=(tr['annotation'][1][k] for k in ('novel exonic splice donor','novel exonic splice acceptor'))
            if g.strand=='-':
                novel1, novel2=novel2,novel1 # sort according to genomic position
            #find the exon numbers
            e1=[next(i for i,e in enumerate(tr['exons']) if e[1]==alt[0]) for alt in novel1]    #alt[0] is the genomic position
            e2=[next(i for i,e in enumerate(tr['exons']) if e[0]==alt[0]) for alt in novel2]    
            candidates=[e for e in e1 if e+1 in e2]
            nc={v[0] for v in tr['noncanonical_splicing']} if 'noncanonical_splicing' in tr else {}
            #splicesite=dict(tr.get('noncanonical_splicing',[]))
            #novel_rl.extend((tr['direct_repeat_len'][c],g.coverage[:,trid],splicesite.get(c,'GTAG')) for c in candidates)
            novel_rl.extend((tr['direct_repeat_len'][c],g.coverage[:,trid]) for c in candidates if c in nc)
        if 'annotation' in tr and tr['annotation'][0]==0: #e.g. FSM
            known_rl.extend((l,g.coverage[:,trid]) for l in tr['direct_repeat_len'])

    known_rl_cov=pd.DataFrame((v[1] for v in known_rl), columns=self.samples)
    novel_rl_cov=pd.DataFrame((v[1] for v in novel_rl), columns=self.samples)
    if groups is not None:
        known_rl_cov=pd.DataFrame({grn:known_rl_cov[grp].sum(1) for grn, grp in groups.items()})
        novel_rl_cov=pd.DataFrame({grn:novel_rl_cov[grp].sum(1) for grn, grp in groups.items()})
    known_rl_cov[known_rl_cov<min_coverage]=0
    novel_rl_cov[novel_rl_cov<min_coverage]=0
    if not weight_by_coverage:        
        known_rl_cov[known_rl_cov>0]=1
        novel_rl_cov[novel_rl_cov>0]=1
    rl={'novel noncanonical':novel_rl, 'known':known_rl}
    cov={'novel noncanonical':novel_rl_cov, 'known':known_rl_cov}
    if isinstance(bins,int):
        bins=np.linspace(x_range[0]-.5,x_range[1]-.5,bins+1)
    counts=pd.DataFrame({f'{sa} {k}':np.histogram([val[0] for val in v],weights=cov[k][sa], bins=bins)[0]  for k,v in rl.items() for sai,sa in enumerate(groups)}   )

    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict( title='direct repeat length',xlabel='length of direct repeats at splice junctons', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
