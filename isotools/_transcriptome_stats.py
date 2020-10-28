from umap import UMAP
from sklearn.decomposition import PCA
from scipy.stats import binom,norm, chi2, betabinom, beta, fisher_exact
from scipy.special import gammaln, gamma, polygamma
from scipy.optimize import minimize
import statsmodels.stats.multitest as multi
import logging
import numpy as np
import pandas as pd
import itertools
from tqdm import tqdm

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
    return(2*norm.sf(z)), (*p1,p0)#two sided alternative
    
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
    return chi2.sf(2*(l1-l0),1),(*p1,p0)

# Differential splicing
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
    for xi,ni in itertools.chain(zip(x,n),((np.concatenate(x),np.concatenate(n)),)):
        xi, ni=xi[ni>0], ni[ni>0] #avoid div by 0
        #x and n must be np arrays
        prob=xi/ni    
        m=prob.mean() #estimate initial parameters
        d=max(prob.var(),1e-6) #max to avoid division by 0
        e=(m**2-m+d) #helper          
        #find ml estimates for a and b
        mle = minimize(loglike_betabinom, x0=[-m*e/d,((m-1)*e)/d],bounds=((1e-6,None),(1e-6,None)),  args=(xi,ni),options={'maxiter': 250}, method='L-BFGS-B', jac=True)
        params.append(mle.x)                
        #mle = minimize(loglike_betabinom2, x0=[-d/(m*e),d/((m-1)*e)],bounds=((1e-9,None),(1e-9,None)),  args=(xi,ni),options={'maxiter': 250}, method='L-BFGS-B', tol=1e-6)
        #params.append([1/p for p in mle.x])                
        if not mle.success:
            logging.debug(f'no convergence in betabinomial fit: k={xi}\nn={ni}\nparams={params}\nmessage={mle.message}') #should not happen to often, mainly with mu close to boundaries
            success=False #prevent calculation of p-values based on non optimal parameters
    # calculate the log likelihoods
    params_alt=[(a/(a+b),  a*b/((a+b)**2*(a+b+1))) for a,b in params] #get alternative parametrization (mu and disp)
    if not success:
        return np.nan, params_alt
    try:
        l0 = betabinom.logpmf(np.concatenate(x), np.concatenate(n), *params[2]).sum() 
        l1 = betabinom.logpmf(x[0],n[0], *params[0]).sum()+betabinom.logpmf(x[1],n[1], *params[1]).sum()
    except ValueError:
        logging.critical(f'betabinom error: x={x}\nn={n}\nparams={params}')#should not happen
        raise
    return chi2.sf(2*(l1-l0),2), params_alt #note that we need two degrees of freedom here as h0 hsa two parameters, h1 has 4

def loglike_betabinom_alt(params, k,n):
    #alternative parametrization with mean, disp - less efficient and seems to have numerical issues
    m,v=params #0<m<1 and 0<v<m-m^2
    e=(m**2-m+v) #helper
    a=-m*e/v
    b=((m-1)*e)/v
    logpdf = gammaln(n+1) + gammaln(k+a) + gammaln(n-k+b) + gammaln(a+b) - \
     (gammaln(k+1) + gammaln(n-k+1) + gammaln(a) + gammaln(b) + gammaln(n+a+b))
    return - np.sum(logpdf) 

def betabinom_lr_test_alt(x,n, disp=None): #alternative formulation
    ''' alternative parametrization of betabinomial log likelihood function L(mu, disp), mainly for testing purposes
     '''
    
    params=list()
    d_range=(1e-6,None) if disp is None else (disp,disp)
    for xi,ni in itertools.chain(zip(x,n),((np.concatenate(x),np.concatenate(n)),)):
        xi, ni=xi[ni>0], ni[ni>0] #avoid div by 0
        #x and n must be np arrays
        #find good initialization parameters for a and b
        prob=xi/ni
        m=prob.mean()
        d=max(prob.var(),1e-6) #to avoid division by 0
        mle = minimize(loglike_betabinom_alt, x0=(m,d),bounds=((1e-6,1-1e-6),d_range),  args=(xi,ni),options={'maxiter': 250})
        if not mle.success:
            logging.debug(f'{mle.message}')
        params.append(mle.x)                
    try:
        params_alt=[((-m*(m**2-m+v))/v,((m-1)*(m**2-m+v))/v)  for m,v in params]
        l0 = betabinom.logpmf(np.concatenate(x), np.concatenate(n), *params_alt[2]).sum() 
        l1 = betabinom.logpmf(x[0],n[0], *params_alt[0]).sum()+betabinom.logpmf(x[1],n[1], *params_alt[1]).sum()
        
    except ValueError:
        logging.error(f'betabinom error: x={x}\nn={n}\nparams={params}')
        raise
    return chi2.sf(2*(l1-l0),2), params #note that we need two degrees of freedom here as h0 hsa two parameters, h1 has 4
    
def altsplice_test(self,groups, min_cov=20, min_n=10, min_sa=.5, test=betabinom_lr_test,padj_method='fdr_bh'):
    #multitest_default={}
    assert len(groups) == 2 , "length of groups should be 2, but found %i" % len(groups)
    if isinstance(groups, dict):
        groupnames=list(groups)
        groups=list(groups.values())
    elif all (isinstance(gn,str) and gn in self.groups for gn in groups):
        groupnames=list(groups)
        groups=[self.groups[gn] for gn in groupnames]
    elif all( isinstance(grp,list) for grp in groups):
        groupnames=['group1','group2']
    else:
        raise ValueError('groups not found in dataset')
    notfound=[sa for grp in groups for sa in grp if sa not in self.samples]
    if notfound:
        raise ValueError(f"Cannot find the following samples: {notfound}")    
    logging.info('testing differential splicing for {g}',g="vs".join(f'{groupnames[i]} ({len(groups[i])})' for i in range(2)) )
    sa_idx={sa:idx[0] for sa,idx in self._get_sample_idx().items()}
    grp_idx=[[sa_idx[sa] for sa in grp] for grp in groups]
    if min_sa<1:
        min_sa*=sum(len(gr) for gr in groups)
    res=[]
    for g in tqdm(self):
        splice_type=['ES','3AS','5AS','IR','ME'] if g.strand=='+' else ['ES','5AS','3AS','IR','ME']
        for junction_cov,total_cov,start,end,sti in g.splice_graph.get_splice_coverage():
            try:
                x=[junction_cov[grp] for grp in grp_idx]
            except IndexError:
                logging.error(f'error at {g.name}:{start}-{end}')
                raise
            n=[total_cov[grp] for grp in grp_idx]
            if sum((ni>=min_n).sum() for ni in n)<min_sa:
                continue
            x_sum=sum(xi.sum() for xi in x)     
            n_sum=sum(ni.sum() for ni in n)       
            if x_sum < min_cov or n_sum-x_sum < min_cov:
                continue
            pval, params=test(x,n)
            res.append(tuple(itertools.chain((g.name,g.id,g.chrom, start, end,splice_type[sti],pval),params ,
                (val for lists in zip(x,n) for pair in zip(*lists) for val in pair ))))
    df=pd.DataFrame(res, columns= (['gene','gene_id','chrom', 'start', 'end','splice_type','pvalue']+ 
            ['prop_'+gn for gn in groupnames+['total']]+  
            [f'{w}_{sa}_{gn}' for gn,grp in zip(groupnames, groups) for sa in grp for w in ['cov', 'span_cov'] ]))
    try:
        mask = np.isfinite(df['pvalue'])
        padj = np.empty(mask.shape)
        padj.fill(np.nan) 
        padj[mask] = multi.multipletests(df.loc[mask,'pvalue'],method=padj_method)[1]
        df.insert(5,'padj',padj)
    except TypeError as e: #apparently this happens if df is empty...
        logging.error('unexpected error during calculation of adjusted p-values: {e}' ,e=e)
    return df

def splice_dependence_test(self,samples=None, min_cov=20,padj_method='fdr_bh',region=None):
    if samples is None:
        samples=self.samples
    sa_dict=self._get_sample_idx()
    sa_idx=[sa_dict[sa][0] for sa in samples]
    res=[]
    for g in tqdm(self.iter_genes(region=region)):
        for ma,j1,j2 in g.splice_graph.splice_dependence(sa_idx, min_cov):
            oddsratio, pval = fisher_exact(ma)    
            res.append((g.name,g.id,g.chrom, j1, j2,pval,oddsratio,ma))
    df=pd.DataFrame(res, columns= (['gene','gene_id','chrom', 'junction1', 'junction1','pvalue', 'oddsratio','counts']))
    try:
        mask = np.isfinite(df['pvalue'])
        padj = np.empty(mask.shape)
        padj.fill(np.nan) 
        padj[mask] = multi.multipletests(df.loc[mask,'pvalue'],method=padj_method)[1]
        df.insert(5,'padj',padj)
    except TypeError as e: #apparently this happens if df is empty...
        logging.error('unexpected error during calculation of adjusted p-values: {e}' ,e=e)
    return df            
        

# PCA and UMAP
def embedding(self, genes=None, reducer='PCA',n_top_var=500, n_subsample=None,samples=None, min_cov=20):
    ''' Compute embedding of most variable junction for each gene if genes are provided as a list (e.g. genes=[gene_name1, gene_name2 ...]). 
        Alternatively, if genes is a dict, junctions can be provided as a list (e.g. genes[gene_name]=[(jstart,jend)]).
        If genes is None, all genes with coverage > min_cov are considered. 
    '''
    joi=dict()
    if samples is not None:
        s_filter=[True if sn in samples else False for sn in self.samples]
    else:
        s_filter=[True]*len(self.samples)
    if genes is None:        
        genes=[g for i,g in enumerate(self) if all(g.splice_graph.weights[s_filter].sum(1)>=min_cov)] #all covered genes
        logging.info(f'found {len(genes)} genes > {min_cov} reads in all {sum(s_filter)} selected samples.')
    else:
        n_in=len(genes)
        if isinstance(genes, dict):
            joi=genes
        genes=[self[g] for i,g in enumerate(genes) if g in self and all(self[g].splice_graph.weights[s_filter].sum(1)>=min_cov)] #all covered genes
        logging.info(f'{len(genes)} of {n_in} genes have > {min_cov} reads in all {sum(s_filter)} selected samples.')
    if n_subsample is not None and n_subsample<len(genes):
        genes=[genes[i] for i in np.random.choice(range(len(genes)), n_subsample)]
        logging.info(f'sampled {n_subsample} random genes')        
    data={}
    # since junctions of a single gene are often correlated, I select one top variable junction per gene
    for g in tqdm(genes):   
        sg=g.splice_graph
        gene_coverage=sg.weights[s_filter].sum(1)
        if any(cov<min_cov for cov in gene_coverage):
            continue
        if g.name in joi and joi[g.name]:#
            for j in joi[g.name]:
                start=next(n for n in sg if n.end==j[0])            
                jcov=sum(sg.weights[s_filter,trid] for trid,suc in start.suc.items() if sg[suc].start==j[1])/gene_coverage
                data[f'{g.name}_{j[0]}_{j[1]}']=jcov
        else:
            jcov={}                 
            for n in sg:
                #print(n)
                for trid,suc in n[3].items():
                    #if all(sg.weights[s_filter,trid] < min_cov):
                    #    continue
                    j_name=f'{g.name}_{n[1]}_{sg[suc][0]}'
                    jcov.setdefault(j_name,np.zeros(len(gene_coverage)))
                    jcov[j_name]+=(sg.weights[s_filter,trid]) /gene_coverage
            if jcov:
                top_idx=pd.DataFrame(jcov).var(0).idxmax()
                data[top_idx]=jcov[top_idx]
    data=pd.DataFrame(data,index=self.infos['sample_table'].loc[s_filter,'name'])
    #filter for highest variance
    logging.info(f'selecting the most variable junction from {min(data.shape[1],n_top_var)} genes.')
    topvar=data[data.var(0).nlargest(n_top_var).index]

    #compute embedding
    if reducer=='PCA':
        # Linear dimensionality reduction using Singular Value Decomposition of the data to project it to a lower dimensional space. 
        # The input data is centered but not scaled for each feature before applying the SVD.
        reducer=PCA()
    elif reducer=='UMAP':
        reducer==UMAP()
    elif reducer is None:
        return topvar    
    return pd.DataFrame(reducer.fit_transform(topvar), index=topvar.index)

# summary tables (can be used as input to plot_bar / plot_dist)
def altsplice_stats(self, coverage=True,groups=None,  min_coverage=1, filter={}):    #todo: filter like in make table
    runs=self.samples

    weights=dict()
    #if groups is not None:
    #    gi={r:i for i,r in enumerate(runs)}
    #    groups={gn:[gi[r] for r in gr] for gn,gr in groups.items()}
    for _,_,tr in self.iter_transcripts(**filter):
        w=tr['coverage'] 
        if groups is not None:
            w=[sum(w[gi] for gi in g) for g in groups.values()]
        if not coverage:
            w=[1 if wi>=min_coverage else 0 for wi in w]
        if tr['annotation'] is None:
            weights['novel/unknown']=weights.get('novel/unknown',np.zeros(len(w)))+w
        else:
            for stype in tr['annotation']['as']:
                weights[stype]=weights.get(stype,np.zeros(len(w)))+w
        weights['total']=weights.get('total',np.zeros(len(w)))+w

    df=pd.DataFrame(weights, index=runs if groups is None else groups.keys()).T
    #if groups is not None:
    #    df=pd.DataFrame({grn:df[grp].sum(1) for grn, grp in groups.items()})
        # sum per group
    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)
    if coverage:
        title='Expressed Transcripts'
        ylab='fraction of reads'
    else:
        title='Different Transcripts'
        ylab='fraction of  different transcripts'
        if min_coverage>1:
            title+=f' > {min_coverage} reads'

    return df, {'ylabel':ylab,'title':title}
    #
    
def filter_stats(self, groups=None, coverage=True, min_coverage=1,consider=None, filter={}):    #todo: filter like in make table
    try:
        runs=self.samples
    except KeyError:
        runs=['self']
    weights=dict()
    #if groups is not None:
    #    gi={r:i for i,r in enumerate(runs)}
    #    groups={gn:[gi[r] for r in gr] for gn,gr in groups.items()}
    current=None
    for g,trid,tr in self.iter_transcripts(**filter):
        if g!=current:
            current=g
            current_cov=g.coverage       
        w= current_cov[:,trid]
        if groups is not None:
            w=[sum(w[gi] for gi in g) for g in groups.values()]
        if not coverage:
            w=[1 if wi>=min_coverage else 0 for wi in w]        
        relevant_filter=[f for f in tr['filter']+g.data['filter'] if f in consider]
        for f in relevant_filter:
            weights[f]=weights.get(f,np.zeros(len(w)))+w
        if not relevant_filter:
            weights['PASS']=weights.get('PASS',np.zeros(len(w)))+w
        weights['total']=weights.get('total',np.zeros(len(w)))+w
            
    df=pd.DataFrame(weights, index=runs if groups is None else groups.keys()).T

    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)
    ylab='fraction of reads' if coverage else 'fraction of different transcripts'
    if coverage:
        title='Expressed Transcripts'
    else:
        title='Different Transcripts'
        if min_coverage>1:
            title+=f' > {min_coverage} reads'
    return df, {'ylabel':ylab,'title':title}

def transcript_length_hist(self=None,  groups=None, add_reference=False,bins=50,x_range=(100,10000),coverage=True,min_coverage=1,use_alignment=False, filter={}, ref_filter={}):
    trlen=[]
    cov=[]
    current=None
    for g,trid,tr in self.iter_transcripts(**filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
        trlen.append(sum(e[1]-e[0] for e in tr['exons']) if use_alignment else tr['source_len']) #source_len is not set in the current version
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(trlen, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if add_reference is not None:
        ref_len=[sum(e[1]-e[0] for e in tr['exons']) for _,_,tr in self.iter_ref_transcripts(**ref_filter)]
        counts['reference']=np.histogram(ref_len, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='linear', title='transcript length',xlabel='transcript length [bp]', ylabel='density', density=True)
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def transcript_coverage_hist(self,  groups=None,bins=50,x_range=(0.5,1000), filter={}):
    # get the transcript coverage in bins for groups
    # return count dataframe and suggested default parameters for plot_distr
    cov=[]
    current=None
    for g,trid,_ in self.iter_transcripts(**filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    counts=pd.DataFrame({gn:np.histogram(g_cov, bins=bins)[0] for gn,g_cov in cov.items()})
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='log', title='transcript coverage',xlabel='reads per transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    #plot histogram
    # cov.mask(cov.lt(x_range[0]) | cov.gt(x_range[1])).plot.hist(ax=ax, alpha=0.5, bins=n_bins)
    # ax=counts.plot.bar() 
    # ax.plot(x, counts)
   
def transcripts_per_gene_hist(self,   groups=None, add_reference=False, bins=50,x_range=(.5,49.5), min_coverage=1, filter={}, ref_filter={}):
    ntr=[]
    current=None
    for g,trid,_ in self.iter_transcripts(**filter):
        if g!=current:
            current=g
            current_cov=g.coverage  
            ntr.append(0)      
        if current_cov[:,trid]>= min_coverage:
            ntr[-1]+=1

    ntr=pd.DataFrame([n for n in ntr if n>0])
    if groups is not None:
        ntr=pd.DataFrame({grn:ntr[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    counts=pd.DataFrame({gn:np.histogram(trl, bins=bins)[0] for gn,trl in ntr.items()})   
    if add_reference:
        ref_ntr=[g.n_ref_transcripts for g in self] #todo: add reference filter
        counts['reference']=np.histogram(ref_ntr, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    sub=f'counting transcripts covered by >= {min_coverage} reads'
    if 'include' in filter:
        sub+=f'; only the following categories:{filter["include"]}'
    if 'remove' in filter:
        sub+=f'; without the following categories:{filter["remove"]}'
    params=dict(yscale='log',title='transcript per gene\n'+sub,xlabel='transcript per gene', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    
def exons_per_transcript_hist(self,  groups=None, add_reference=False, bins=35,x_range=(0.5,69.5),coverage=True,  min_coverage=2, filter={}, ref_filter={}):
    n_exons=[]
    cov=[]
    current=None
    for g,trid,tr in self.iter_transcripts(**filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
        n_exons.append(len(tr['exons']))
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(n_exons, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if add_reference:
        ref_n_exons=[len(tr['exons']) for _,_,tr in self.iter_ref_transcripts(**ref_filter)]
        counts['reference']=np.histogram(ref_n_exons, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='log', title='exons per transcript',xlabel='number of exons per transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def downstream_a_hist(self, groups=None,add_reference=False,bins=30,x_range=(0,1),coverage=True,  min_coverage=2, filter={}, ref_filter={}):
    acontent=[]
    cov=[]
    current=None
    for g,trid,tr in self.iter_transcripts(**filter):
        if g!=current:
            current=g
            current_cov=g.coverage        
        cov.append(current_cov[:,trid])
        acontent.append(tr['downstream_A_content'])
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(acontent, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if add_reference is not None:
        try:
            ref_acontent=[tr['downstream_A_content'] for _,_,tr in self.iter_ref_transcripts(**ref_filter)]
        except KeyError:
            logging.error('No A content for reference found. Run add_biases for reference first')
        else:
            counts['reference']=np.histogram(ref_acontent, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict( title='downstream genomic A content',xlabel='fraction of A downstream the transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
