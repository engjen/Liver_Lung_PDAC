import os
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import math
import re
import warnings
from scipy.stats import pearsonr
import networkx as nx
from bokeh.palettes import Colorblind
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm, gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve, RocCurveDisplay
mpl.rc('figure', max_open_warning = 0)

import sklearn
from sklearn.preprocessing import minmax_scale, scale, FunctionTransformer
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import lifelines
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines import exceptions
warnings.filterwarnings("ignore",category = exceptions.ApproximationWarning)

import scipy
from scipy import stats
from scipy.stats import entropy, norm
from scipy.spatial import cKDTree

import skimage
from skimage.filters import unsharp_mask
from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma)
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage import color, morphology
from skimage.transform import rescale
import tifffile
from scipy.ndimage import median_filter
from skimage.util import img_as_ubyte,  img_as_float
from math import sqrt
import statsmodels
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.api as sm
import statannotations
from statannotations.Annotator import Annotator
from itertools import combinations

import anndata
from anndata import AnnData
import plotly.express as px

codedir = os.getcwd()
def plot_violins(plotting,ls_groups,b_correct=False,figsize=(2.4,3)):
    '''
    '''
    if figsize == 0:
        fig, ax = plt.subplots(figsize=(1.3*len(ls_groups),3.2),dpi=300)
    else:
        fig, ax = plt.subplots(figsize=figsize,dpi=300)
    pairs = [item for item in itertools.combinations(ls_groups,2)]
    sns.violinplot(**plotting,ax=ax,alpha=0.2,linewidth=1,cut=0,inner=None,hue_order=ls_groups,color='white')
    sns.boxplot(**plotting,ax=ax,showmeans=True,medianprops={'visible': False},whiskerprops={'visible': False},meanline=True,showcaps=False,meanprops={'color': 'k', 'ls': '-', 'lw': 2},showfliers=False,showbox=False)#
    sns.stripplot(**plotting,s=4,dodge=True,ax=ax,jitter=0.2,alpha=0.8,hue_order=ls_groups) #
    #dummy
    pvalues, pairs2, test_results = extract_pvals(plotting,pairs,test='t-test_ind')
    #correct
    reject, corrected, __, __ = statsmodels.stats.multitest.multipletests(pvalues,method='fdr_bh')
    formatted_pvalues = [f'P={pvalue:.2}' for pvalue in list(pvalues)]
    if b_correct:
        formatted_pvalues = [f'FDR={pvalue:.2}' for pvalue in list(corrected)]
    # create real annotator
    ax = real_annotator(plotting,pairs2,formatted_pvalues,ax)
    ax = labels_with_n(plotting,ls_groups,ax)
    return(fig,ax,test_results)

def extract_pvals(plotting,pairs,test='t-test_ind'):
    #dummy
    f,a = plt.subplots()
    annotd = Annotator(a,pairs,**plotting,verbose=False)
    annotd.configure(test=test,text_format="star")
    annotd.apply_test()
    a, test_results = annotd.apply_test().annotate()
    plt.close(f)
    #extract annotations
    pvalues = []
    pairs2 = []
    for res in test_results:
        pairs2.append((res.data.group1,res.data.group2))
        pvalues.append(res.data.pvalue)
    return(pvalues, pairs2, test_results)
    
def real_annotator(plotting,pairs2,formatted_pvalues,ax):
    annot = Annotator(ax,pairs2,**plotting,verbose=False)
    annot.configure(fontsize='small')
    annot.set_custom_annotations(formatted_pvalues)
    annot.annotate()
    return(ax)

def labels_with_n(plotting,ls_groups,ax):
    # labels with N
    linebreak = '\n'
    df_mean = plotting['data']
    s_group = plotting['x']
    if len(ls_groups)==2:
        xlabels = [f'{ls_groups[0]}{linebreak}(N={(df_mean.loc[:,s_group].value_counts()[ls_groups[0]])})',
               f'{ls_groups[1]}{linebreak}(N={(df_mean.loc[:,s_group].value_counts()[ls_groups[1]])})']
    elif len(ls_groups)==3:
        xlabels = [f'{ls_groups[0]}{linebreak}(N={(df_mean.loc[:,s_group].value_counts()[ls_groups[0]])})',
               f'{ls_groups[1]}{linebreak}(N={(df_mean.loc[:,s_group].value_counts()[ls_groups[1]])})',
                  f'{ls_groups[2]}{linebreak}(N={(df_mean.loc[:,s_group].value_counts()[ls_groups[2]])})']
    elif len(ls_groups)==4:
        xlabels = [f'{ls_groups[0]}{linebreak}N={(df_mean.loc[:,s_group].value_counts()[ls_groups[0]])}',
               f'{ls_groups[1]}{linebreak}N={(df_mean.loc[:,s_group].value_counts()[ls_groups[1]])}',
                  f'{ls_groups[2]}{linebreak}N={(df_mean.loc[:,s_group].value_counts()[ls_groups[2]])}',
               f'{ls_groups[3]}{linebreak}N={(df_mean.loc[:,s_group].value_counts()[ls_groups[3]])}']
    else:
        return(ax)
    ax.set_xticklabels(xlabels,fontsize='small')
    ax.set_xlabel(ax.get_xlabel().replace("_"," "))
    ax.set_ylabel(ax.get_ylabel().replace("_"," "))
    return(ax)
        
def annotated_stripplot(plotting,ls_groups,b_correct=False,figsize=(2.4,3),s=4):
    '''
    no hue, show pvalues
    default no correction
    '''
    if figsize == 0:
        fig, ax = plt.subplots(figsize=(1.3*len(ls_groups),3.2),dpi=300)
    else:
        fig, ax = plt.subplots(figsize=figsize,dpi=300)
    sns.stripplot(**plotting,ax=ax,alpha=0.8,label='__',s=s)
    sns.boxplot(**plotting,ax=ax,showmeans=True,medianprops={'visible': False},
                                   whiskerprops={'visible': False},meanline=True,showcaps=False,
                           meanprops={'color': 'k', 'ls': '-', 'lw': 2},showfliers=False,showbox=False)
    pairs = [item for item in itertools.combinations(ls_groups,2)]
    pvalues, pairs2, test_results = extract_pvals(plotting,pairs)
    #correct
    reject, corrected, __, __ = statsmodels.stats.multitest.multipletests(pvalues,method='fdr_bh')
    formatted_pvalues = [f'P={pvalue:.2}' for pvalue in list(pvalues)]
    if b_correct:
        formatted_pvalues = [f'FDR={pvalue:.2}' for pvalue in list(corrected)]
    # create real annotator
    ax = real_annotator(plotting,pairs2,formatted_pvalues,ax)
    ax = labels_with_n(plotting,ls_groups,ax)
    plt.tight_layout()
    return(fig,ax,test_results)
    
def annotated_stripplot_hue(df,x,y,hue,order,hue_order,figsize,b_correct=True,loc='upper right',s=2,b_box=False):
    '''
    with hue, show pvalues, default FDR corrected 
    '''

    plotting = {"data":df,"x":x,"y":y,"order":order,
               "hue":hue,"hue_order":hue_order}
    pairs = [((item,hue_order[0]),(item,hue_order[1])) for item in order]

    fig,ax = plt.subplots(dpi=200,figsize=figsize)
    f,a = plt.subplots()
    sns.stripplot(**plotting,dodge=True,ax=ax,s=s,alpha=0.7)
    if b_box:
        sns.boxplot(**plotting,ax=ax,showfliers=False,saturation=0.6,linewidth=0.5)
    else:
        sns.boxplot(**plotting,ax=ax,showmeans=True,medianprops={'visible': False},
                                whiskerprops={'visible': False},meanline=True,showcaps=False,
                           meanprops={'color': 'k', 'ls': '-', 'lw': 2},showfliers=False,showbox=False) 
    annot = Annotator(a,pairs,**plotting,verbose=False)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[0:2],l[0:2],loc=loc,title=hue)
    annot.configure(test='t-test_ind',text_format="star")
    annot.apply_test()
    a, test_results = annot.apply_test().annotate()
    plt.close(f)
    annot = Annotator(ax,pairs,**plotting,verbose=False)
    d_pval = {}
    for res in test_results:
        d_pval.update({res.data.group1[0]:res.data.pvalue})
    pvalues = [d_pval[item] for item in order]
    reject, corrected, __, __ = statsmodels.stats.multitest.multipletests(pvalues,method='fdr_bh')
    formatted_pvalues = [f'p={pvalue:.2}' for pvalue in list(pvalues)]
    if b_correct:
        formatted_pvalues = [f'FDR={pvalue:.2}' for pvalue in list(corrected)]
    annot.set_custom_annotations(formatted_pvalues)
    annot.annotate()
    return(fig, ax)

def annotated_stripplot_pairs(plotting,pairs,order,test='Mann-Whitney'):
    '''
    no hue, mann whitney, no corrections
    '''
    fig,ax=plt.subplots(dpi=300,figsize=(4+ len(pairs)*.7,3.5))
    #create dummy annotator
    f, a = plt.subplots()
    annotator = Annotator(ax=a, pairs=pairs,**plotting)
    a, test_results  = annotator.configure(test=test, text_format='full',verbose=0).apply_test().annotate()
    plt.close(f)
    # extract pvalues
    pairs2 = []
    #order = []
    for res in test_results:
        pairs2.append((res.data.group1,res.data.group2))
        # order.append(res.data.group1)
        # order.append(res.data.group2)
    print(pairs2)
    print(order)
    sns.stripplot(**plotting,ax=ax,alpha=0.8,order=order)
    sns.boxplot(**plotting,ax=ax,order=order,showmeans=True,medianprops={'visible': False}, whiskerprops={'visible': False},meanline=True,showcaps=False,meanprops={'color': 'k', 'ls': '-', 'lw': 2},showfliers=False,showbox=False)
    # create real annotator
    annot = Annotator(ax=ax, pairs=pairs2,**plotting,order=order)
    formatted_pvalues = [f"P={res.data.pvalue:.2}" for res in test_results]
    annot.configure(fontsize='large',verbose=0)
    annot.set_custom_annotations(formatted_pvalues)
    annot.annotate()
    plt.tight_layout()
    return(fig,ax, test_results,order,pairs2)
    
def annotated_stripplot3(plotting,ls_groups,label,hue):
    '''
    with hue, default ttest FDR corrected, stars
    '''
    fig, ax = plt.subplots(figsize=(1.3*len(ls_groups),3.2),dpi=300)
    sns.stripplot(**plotting,hue=hue,ax=ax,alpha=0.8,label=label,dodge=False)
    sns.boxplot(**plotting,ax=ax,
                       whiskerprops={'visible': False},showcaps=False,
                       #meanline=True,showmeans=True,medianprops={'visible': False},
                       medianprops={'color': 'k', 'ls': '-', 'lw': 1},#meanprops
                showfliers=False,showbox=False)
    pairs = [item for item in itertools.combinations(ls_groups,2)]
    annot = Annotator(ax,pairs,**plotting)
    annot.configure(test='t-test_ind',comparisons_correction="fdr_bh",text_format='star', verbose=False)
    ax, test_results = annot.apply_test().annotate()
    return(fig,ax,test_results)

#excel func
def open_write_excel(d_result_fig,filename='Source_Data_Figure3.xlsx'):
    if not os.path.exists(filename):
        with pd.ExcelWriter(filename,mode='w',engine="openpyxl") as writer:
            for sheet_name, df in d_result_fig.items():
                print(f'saving {sheet_name[0:30]}')
                df.to_excel(writer, sheet_name=sheet_name[0:30])
    else:
        with pd.ExcelWriter(filename,mode='a',if_sheet_exists='replace',engine="openpyxl") as writer:
            for sheet_name, df in d_result_fig.items():
                print(f'saving {sheet_name[0:30]}')
                df.to_excel(writer, sheet_name=sheet_name[0:30])

def excel_dict_from_csvs(s_find,s_prefix='ED_Fig2',s_dir='.'):
    df = pd.DataFrame(sorted(os.listdir(s_dir)),columns=['Test'])
    df_test = df[df.Test.str.contains(s_find)].copy()
    df_test['Panel'] = [item.split(s_find)[1].split('.csv')[0] for item in df_test.Test]
    df_test['Sheet'] = [f'{s_prefix}{item}' for item in df_test.Panel]
    d_fig = {}
    for ids, row in df_test.iterrows():
        print(row.Test)
        df_sheet = pd.read_csv(row.Test,index_col=0)
        d_fig.update({row.Sheet:df_sheet})
    return(d_fig)
    
def get_blobs2(image_gray,min_sigma,max_sigma,threshold,exclude_border):

    blobs_dog = blob_dog(image_gray,  min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold,exclude_border=exclude_border)
    blobs_dog[:, 2] = blobs_dog[:, 2] * sqrt(2)

    blobs_list = [image_gray,  blobs_dog] #blobs_doh ,
    colors = ['red','red', ]
    titles = ['Original','Difference of Gaussian',#
              ]
    sequence = zip(blobs_list, colors, titles)

    fig, axes = plt.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)
    ax = axes.ravel()

    for idx, (blobs, color, title) in enumerate(sequence):
        if idx == 1:
            ax[idx].set_title(f'{title}\nmin={min_sigma} max={max_sigma} thresh={threshold}')
        else:
            ax[idx].set_title(f'{title}')
        ax[idx].imshow(image)
        if not title == 'Original':
            for blob in blobs:
                y, x, r = blob
                c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
                ax[idx].add_patch(c)
        #ax[idx].set_axis_off()

    plt.tight_layout()
    plt.close(fig)
    return(blobs_dog,fig)


def correlation_regplot(df_merge,s_tum,s_marker,alpha=0.05,b_legend=False):    
    '''
    df_merge: dataframe with two columns to correlate
    s_tum: x value
    s_marker: y value
    '''
    try:
        r, pvalue = stats.spearmanr(a=df_merge.loc[:,s_tum], b=df_merge.loc[:,s_marker])
        #print(pvalue)
    except:
        print(f'skipping {s_marker}')
        pvalue=1
    if pvalue < alpha:
        fig, ax = plt.subplots(figsize=(3,3),dpi=300)
        sns.regplot(x=s_tum, y=s_marker, data=df_merge,ax=ax)  
        ax.set_title(f'{s_marker}\n vs {s_tum}\nSpearman r={r:.2f} p={pvalue:.2f} n={len(df_merge)}') 
        if b_legend:
            ax.legend(bbox_to_anchor=(1.01,.3))
        plt.tight_layout()
    else:
        fig=None
        ax=None
    return(fig,ax)

def km_plot(df,s_col,s_time,s_censor,loc='upper center',fontsize='small',figsize=(3,3),b_pval=True):
    df = df.loc[:,[s_col,s_time,s_censor]].dropna().copy()
    results = multivariate_logrank_test(event_durations=df.loc[:,s_time],
                                    groups=df.loc[:,s_col], event_observed=df.loc[:,s_censor])        
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=figsize,dpi=400)
    ls_order = sorted(df.loc[:,s_col].dropna().unique())
    for s_group in ls_order:
        #print(s_group)
        df_abun = df[df.loc[:,s_col]==s_group]
        durations = df_abun.loc[:,s_time]
        event_observed = df_abun.loc[:,s_censor]
        kmf.fit(durations,event_observed,label=s_group)
        kmf.plot(ax=ax,ci_show=True,show_censors=True,ci_alpha=0.08,censor_styles={"marker": "|"})
        print(f'{s_col} {s_group} Median = {kmf.median_survival_time_}')
    if b_pval:#len(ls_order)==2:
        pvalue = f"{results.summary.p[0]:.2}"
        demo_con_style(ax,0.74,0.9,0.73,0.83, f"P={pvalue}")
        ax.set_title(f'{s_col.replace("_"," ")}')
    else:
        ax.set_title(f'{s_col.replace("_"," ")}\nP={results.summary.p[0]:.2}')
    pvalue_int = results.summary.p[0]
    ls_n = [df.loc[:,s_col].value_counts()[item] for item in ls_order]
    print(ls_n)
    h,l = ax.get_legend_handles_labels()
    print(l)
    labels = [f'{item}\n(N={ls_n[idx]})' for idx, item in enumerate(l)]
    ax.legend(labels=labels,handles=h,fancybox=False,frameon=False,loc=loc,fontsize=fontsize)
    ax.set_ylim(-0.05,1.05)
    # Hide the right and top spines
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_ylabel(f'Probability of {s_censor}')
    ax.set_xlabel(s_time)
    plt.tight_layout()
    return(fig,ax,pvalue_int)
    
def demo_con_style(ax,x,top, bottom,ytext,pvalue):
    connectionstyle = "bar,fraction=0.2"
    x1 = x #x1, y1 = 0.8, 0.9
    x2 = x  #x2, y2 = 0.8, 0.8
    y1 = top
    y2 = bottom
    #ax.plot([x1, x2], [y1, y2], ".",color='white')
    ax.annotate("",#pvalue,#
                xy=(x1, y1), xycoords=ax.transAxes,#'data',
                xytext=(x2, y2), textcoords=ax.transAxes,#'data',
                arrowprops=dict(arrowstyle="-", color="k",
                                #shrinkA=5, shrinkB=5,
                                patchA=None, patchB=None,
                                connectionstyle=connectionstyle,
                                ),
                )
    ax.text(x+0.06, ytext, pvalue,transform=ax.transAxes,
            ha="left", va="top",fontsize='medium')

def km_plot_overlay(df_km,ls_col,s_time,s_censor,loc='upper center',fontsize='small',figsize=(4,4)):
    fig, ax = plt.subplots(figsize=figsize,dpi=300)
    ls_title = []
    nl = '\n'
    for s_col in ls_col:
        df = df_km.loc[:,[s_col,s_time,s_censor]].dropna().copy()
        results = multivariate_logrank_test(event_durations=df.loc[:,s_time],
                                        groups=df.loc[:,s_col], event_observed=df.loc[:,s_censor])        
        kmf = KaplanMeierFitter()
        
        ls_order = sorted(df.loc[:,s_col].dropna().unique())
        for s_group in ls_order:
            #print(s_group)
            df_abun = df[df.loc[:,s_col]==s_group]
            durations = df_abun.loc[:,s_time]
            event_observed = df_abun.loc[:,s_censor]
            kmf.fit(durations,event_observed,label=s_group)
            kmf.plot(ax=ax,ci_show=True,show_censors=True,ci_alpha=0.15,censor_styles={"marker": "|"})
            print(f'{s_col} {s_group} Median = {kmf.median_survival_time_}')
        s_title = f'{s_col.replace("_"," ")}: P={results.summary.p[0]:.2}' #N={[df.loc[:,s_col].value_counts()[item] for item in ls_order]}'
        ls_n = [df.loc[:,s_col].value_counts()[item] for item in ls_order]
        ls_title.append(ls_n)
    ls_title = flatten(ls_title)
    print(ls_title)
    h,l = ax.get_legend_handles_labels()
    print(l)
    labels = [f'{item}\n(N={ls_title[idx]})' for idx, item in enumerate(l)]
    ax.legend(labels=labels,handles=h,fancybox=False,frameon=False,loc=loc,fontsize=fontsize)
    ax.set_ylim(-0.05,1.05)
    # Hide the right and top spines
    ax.spines[['right', 'top']].set_visible(False)
    return(fig,ax,ls_order)

def flatten(xss):
    return [x for xs in xss for x in xs]
    
def cph_plot(df,s_multi,s_time,s_censor,figsize=(3,3)):
    cph = CoxPHFitter()  #penalizer=0.1
    if df.columns.isin(['Stage']).any():
        df.Stage = df.Stage.replace({'I':1,'II':2,'III':3,'IV':4})
    cph.fit(df.dropna(), s_time, event_col=s_censor) 
    fig, ax = plt.subplots(figsize=figsize,dpi=300)
    cph.plot(ax=ax)
    pvalue = cph.summary.p[s_multi]
    hr = cph.summary.loc[:,'exp(coef)'][s_multi]
    ax.set_title(f'{s_multi}\nHR={hr:.3},P={pvalue:.2} N={(len(df.dropna()))}')
    plt.tight_layout()
    return(fig,cph)
        
def plot_pearson(df_pri,s_porg,s_foci,s_stats,ls_plots=['Primaries','Mets','Both']):
    pvalues = []
    fig,ax=plt.subplots(1,len(ls_plots),figsize=(len(ls_plots)*3.3,3.5), sharex='col',sharey='row',squeeze=False,dpi=200)
    ax = ax.ravel()
    for idx,s_met_pri in enumerate(ls_plots):
        #print(s_met_pri)
        if s_met_pri == 'Mets':
            b_met = df_pri.Tumor_Type=='Met'#.str.contains('-M',na=False)
        elif s_met_pri == 'Primaries':
            b_met = df_pri.Tumor_Type=='Primary'#~(df_pri.Patient_Specimen_ID.str.contains('-M',na=True))
        else:
            b_met = df_pri.Public_Patient_ID.str.contains('ST-',na=False)
        sns.regplot(data=df_pri.loc[b_met,[s_porg,s_foci]],x=s_foci,y=s_porg,label=s_met_pri,ax=ax[idx],color="tab:blue",scatter_kws={'s':3})
        y = df_pri.loc[b_met,[s_foci,s_porg]].dropna().loc[:,s_foci]
        x = df_pri.loc[y.index,s_porg]
        if s_stats == 'non-parametric':
            statistic, pvalue = stats.spearmanr(x, y)
        else:
            statistic, pvalue = stats.pearsonr(x, y)
        pvalues.append(pvalue)
        if idx > 0:
            ax[idx].set_ylabel('')
        else:
            ax[idx].set_ylabel(s_porg.replace('trim_padj_0.2_','').replace("_"," "))
        ax[idx].set_title(f'{s_met_pri} r={statistic:.3} P={pvalue:.3} N={len(x)}')
    plt.tight_layout()
    return(fig, pvalues)

# QQ-plot
import statsmodels.api as sm
import matplotlib.pyplot as plt
from bioinfokit.analys import stat
def qq_plot_hist(df_pri,s_cat,s_foci):
    fig, ax = plt.subplots(2,1)
    res = stat()
    df_melt = df_pri.loc[~df_pri.loc[:,s_cat].isna(),[s_cat,s_foci]]
    df_melt.columns = ['treatments', 'value']
    res.anova_stat(df=df_melt, res_var='value', anova_model='value ~ C(treatments)')
    # res.anova_std_residuals are standardized residuals obtained from ANOVA (check above)
    sm.qqplot(res.anova_std_residuals, line='45',ax=ax[0])
    ax[0].set_title(f'{s_foci} ({s_cat})')
    ax[0].set_xlabel("Theoretical Quantiles")
    ax[0].set_ylabel("Standardized Residuals")

    # histogram
    ax[1].hist(res.anova_model_out.resid, bins='auto', histtype='bar', ec='k') 
    ax[1].set_title(f'')
    ax[1].set_xlabel("Residuals")
    ax[1].set_ylabel('Frequency')
    plt.tight_layout()
    
def add_quartiles(df_merge,s_porg):
    x = df_merge.loc[:,s_porg].dropna()
    b_cut = df_merge.loc[:,s_porg].dropna().index
    #print(len(x))
    d_cut = {'quartiles':(4,['low','med-low','med-high','high']),
             'tertiles' : (3,['low','med','high']),
            'medians' : (2,['low','high'])}
    for s_col, tu_cut in d_cut.items():
        i_cut = tu_cut[0]
        labels = tu_cut[1]
        q = pd.qcut(x, q=i_cut,labels=labels) 
        if s_col == 'quartiles':
            df_merge[s_col] = pd.NA #np.NaN
            df_merge.loc[b_cut,s_col] = q.replace({'med-low':pd.NA,'med-high':pd.NA})
        elif s_col == 'tertiles':
            df_merge[s_col] = pd.NA
            df_merge.loc[b_cut,s_col] = q.replace({'med':pd.NA})#'high'
        else:
            df_merge[s_col] = pd.NA
            df_merge.loc[b_cut,s_col] = q
        #print(df_merge[s_col].value_counts())
    return(df_merge)


# Liver vermillion
# Lung blue
# High pORG orange
# Low pORG sky blue
# Basal black
# Classical reddish purple
# High pSUB bluish green
# Low pSUB yellow

def violin_stats2(df_pri,d_order,s_foci,s_stats):
    order = []
    d_pval = {}
    df_both = pd.DataFrame()
    for idx, s_cat in enumerate(d_order.keys()):
        ls_order = d_order[s_cat]
        s_bad = ls_order[0]
        s_good = ls_order[1]
        d_replace = {s_bad:'bad',s_good:'good'}
        a = df_pri.loc[df_pri.loc[:,s_cat]==ls_order[0],s_foci].dropna()
        b = df_pri.loc[df_pri.loc[:,s_cat]==ls_order[1],s_foci].dropna()
        if s_stats == 'mean':
            statistic, pvalue = stats.f_oneway(b,a)
            #statistic, pvalue = stats.ttest_ind(b,a)
        elif s_stats == 'non-parametric':
            statistic, pvalue = stats.kruskal(b,a)
        df_pri['hue'] = df_pri.loc[:,s_cat].replace(d_replace)
        df_pri['x'] = s_cat
        df_both=pd.concat([df_both,df_pri.loc[df_pri.hue.isin(['bad','good']),['x','hue',s_foci,s_cat]].rename({s_cat:'color'},axis=1)])
        for s_test in ls_order:
            order.append(s_test)
        d_pval.update({s_cat:pvalue})
    return(df_both,d_pval,order)

def plot_violins2(df_both,d_pval,d_order,s_stats,s_foci,order,d_colorblind,s_porg,b_correct=False,figsize=(3,3),s=4,b_title=True,alpha=0.8):
    fig,ax=plt.subplots(dpi=300,figsize=figsize)
    hue_order = df_both.sort_values(by=['x','hue']).color.unique()
    if s_stats == 'non-parametric':
        sns.violinplot(data=df_both,y=s_foci,x='color',ax=ax,alpha=1,linewidth=1,cut=0,inner='quartile',order=hue_order,color='white')#
    elif s_stats == 'mean':
        sns.violinplot(data=df_both,y=s_foci,x='color',ax=ax,alpha=1,linewidth=1,cut=0,inner=None,
                       order=hue_order,color='white')
        sns.boxplot(data=df_both,y=s_foci,x='color',ax=ax,showmeans=True,medianprops={'visible': False},
                       whiskerprops={'visible': False},meanline=True,showcaps=False,order=hue_order,
                       meanprops={'color': 'k', 'ls': '-', 'lw': 2},showfliers=False,showbox=False)#
    sns.stripplot(data=df_both,y=s_foci,x='color',s=s,dodge=True,ax=ax,palette=d_colorblind,jitter=0.2,alpha=alpha,order=hue_order) #
    #annotate
    if len(order) == 6:
        pairs = [(order[0],order[1]),(order[2],order[3]),(order[4],order[5])]
        pvalues = [d_pval[list(d_order.keys())[0]],d_pval[list(d_order.keys())[1]],d_pval[list(d_order.keys())[2]]]
    elif len(order) == 4:
        pairs = [(order[0],order[1]),(order[2],order[3])]
        pvalues = [d_pval[list(d_order.keys())[0]],d_pval[list(d_order.keys())[1]]]
    elif len(order) == 2:
        pairs = [(order[0],order[1])]
        pvalues = [d_pval[list(d_order.keys())[0]]]
    elif len(order) == 8:
        pairs = [(order[0],order[1]),(order[2],order[3]),(order[4],order[5]),(order[6],order[7])]
        pvalues = [d_pval[list(d_order.keys())[0]],d_pval[list(d_order.keys())[1]],d_pval[list(d_order.keys())[2]],d_pval[list(d_order.keys())[3]]]
    reject, corrected, __, __ = statsmodels.stats.multitest.multipletests(pvalues,method='fdr_bh')
    formatted_pvalues = [f'P={pvalue:.2}' for pvalue in list(pvalues)]
    if b_correct:
        formatted_pvalues = [f'P={pvalue:.2}' for pvalue in list(corrected)]
    annotator = Annotator(ax, pairs=pairs, data=df_both,y=s_foci,x='color',verbose=False)
    annotator.set_custom_annotations(formatted_pvalues)
    annotator.annotate()
    #ax.legend().remove()
    df_label = df_both.dropna().groupby('x').count().hue.loc[d_order.keys()]
    s_replace = s_porg.split('_')[0]
    ls_labs = [f'{item.replace("quartiles",f"{s_replace} quartiles")} N={df_label[item]}' for item in df_label.index]
    ax.set_xlabel(" | ".join(ls_labs),fontsize=8)
    ax.set_ylabel(ax.get_ylabel(),fontsize=8)
    if b_title:
        ax.set_title(f"{s_foci.replace('_',' ')} ({s_porg.split('_')[-1][0:3]})", fontsize='large',loc='right',pad=10) #
    plt.tight_layout()
    return(fig,pvalues,corrected)

def plot_violins3(df_both,s_stats,s_foci,s_comp,s_porg,hue_order,hue='TCR_Met_Site',figsize=(3,3)):
    fig,ax=plt.subplots(dpi=300,figsize=figsize)
    if s_stats == 'non-parametric':
        sns.violinplot(data=df_both,y=s_foci,x=s_comp,ax=ax,alpha=0.2,linewidth=1,cut=0,
                       inner='quartile',color='white')#
    elif s_stats == 'mean':
        sns.violinplot(data=df_both,y=s_foci,x=s_comp,ax=ax,alpha=0.2,linewidth=1,cut=0,inner=None,
                       color='white',scale='width')
        sns.boxplot(data=df_both,y=s_foci,x=s_comp,ax=ax,showmeans=True,medianprops={'visible': False},
                       whiskerprops={'visible': False},meanline=True,showcaps=False,
                       meanprops={'color': 'k', 'ls': '-', 'lw': 2},showfliers=False,showbox=False)#
    sns.stripplot(data=df_both,y=s_foci,x=s_comp,hue=hue,s=8,dodge=False,ax=ax,jitter=0.2,alpha=0.8,hue_order=hue_order,
                 palette='tab20') #
    ax.legend(bbox_to_anchor=(1,1),fontsize='small')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=45,fontsize=8)
    #ax.set_title(f"{s_foci.replace('_',' ')} ({s_porg.split('_')[-1][0:3]})", loc='left',fontsize='large',pad=10) #
    plt.tight_layout()
    return(fig, ax)


def violin_stats(df_pri,d_order,s_foci,s_stats):
    order = []
    ls_ticks = []
    d_pval = {}
    df_both = pd.DataFrame()
    for idx, s_cat in enumerate(d_order.keys()):
        ls_order = d_order[s_cat]
        s_bad = ls_order[0]
        s_good = ls_order[1]
        d_replace = {s_bad:'bad',s_good:'good'}
        a = df_pri.loc[df_pri.loc[:,s_cat]==ls_order[0],s_foci].dropna()
        b = df_pri.loc[df_pri.loc[:,s_cat]==ls_order[1],s_foci].dropna()
        if s_stats == 'mean':
            statistic, pvalue = stats.f_oneway(b,a)
        elif s_stats == 'non-parametric':
            statistic, pvalue = stats.kruskal(b,a)
        df_pri['hue'] = df_pri.loc[:,s_cat].replace(d_replace)
        df_pri['x'] = s_cat
        df_both=pd.concat([df_both,df_pri.loc[df_pri.hue.isin(['bad','good']),['x','hue',s_foci,s_cat]].rename({s_cat:'color'},axis=1)])
        for s_test in ls_order:
            order.append((s_cat,d_replace[s_test]))
            ls_ticks.append(s_test)
        d_pval.update({s_cat:pvalue})
    return(df_both,d_pval,order,ls_ticks)
    
def plot_violins_OLD(df_both,d_pval,d_order,s_stats,s_foci,order,ls_ticks,b_correct=False):
    figsize=(3,3)
    fig,ax=plt.subplots(dpi=300,figsize=figsize)
    if s_stats == 'non-parametric':
        sns.violinplot(data=df_both,y=s_foci,x='x',hue='hue',ax=ax,alpha=0.2,linewidth=1,cut=0,inner='quartile',hue_order=['bad','good'],color='white')#
    elif s_stats == 'mean':
        sns.violinplot(data=df_both,y=s_foci,x='x',hue='hue',ax=ax,alpha=0.2,linewidth=1,cut=0,inner=None,
                       hue_order=['bad','good'],color='white')
        sns.boxplot(data=df_both,y=s_foci,x='x',hue='hue',ax=ax,showmeans=True,medianprops={'visible': False},
                       whiskerprops={'visible': False},meanline=True,showcaps=False,
                       meanprops={'color': 'k', 'ls': '-', 'lw': 2},showfliers=False,showbox=False)#
    sns.stripplot(data=df_both,y=s_foci,x='x',hue='hue',s=4,dodge=True,ax=ax,palette=Colorblind[8],jitter=0.2,alpha=0.8,hue_order=['bad','good']) #
    #annotate
    if len(order) == 6:
        pairs = [(order[0],order[1]),(order[2],order[3]),(order[4],order[5])]
        pvalues = [d_pval[list(d_order.keys())[0]],d_pval[list(d_order.keys())[1]],d_pval[list(d_order.keys())[2]]]
    elif len(order) == 4:
        pairs = [(order[0],order[1]),(order[2],order[3])]
        pvalues = [d_pval[list(d_order.keys())[0]],d_pval[list(d_order.keys())[1]]]
    else:
        pairs = [(order[0],order[1]),(order[2],order[3]),(order[4],order[5]),(order[6],order[7])]
        pvalues = [d_pval[list(d_order.keys())[0]],d_pval[list(d_order.keys())[1]],d_pval[list(d_order.keys())[2]],d_pval[list(d_order.keys())[3]]]
    reject, corrected, __, __ = statsmodels.stats.multitest.multipletests(pvalues,method='fdr_bh')
    formatted_pvalues = [f'p={pvalue:.2}' for pvalue in list(pvalues)]
    if b_correct:
        formatted_pvalues = [f'p={pvalue:.2}' for pvalue in list(corrected)]
    annotator = Annotator(ax, pairs=pairs, data=df_both,y=s_foci,x='x',hue='hue',verbose=False)
    annotator.set_custom_annotations(formatted_pvalues)
    annotator.annotate()
    ax.legend().remove()
    if len(order) == 6:
        ax.set_xticks([-0.2,0.2, 0.8,1.2,1.8,2.2])
    elif len(order) == 4:
        ax.set_xticks([-0.2,0.2, 0.8,1.2])
    else:
        ax.set_xticks([-0.2,0.2, 0.8,1.2,1.8,2.2,2.8,3.2])
    ax.set_xticklabels(ls_ticks,rotation=45)
    df_label = df_both.dropna().groupby('x').count().hue.loc[d_order.keys()]
    ls_labs = [f'{item.replace("Subtype","")} n={df_label[item]}' for item in df_label.index]
    ax.set_xlabel(" | ".join(ls_labs),fontsize=8)
    ax.set_title(f"{s_foci}", fontsize='large') #
    plt.tight_layout()
    return(fig,pvalues,corrected)
    
############################################# for survival/spatial notebook ########################
def more_plots(adata,df_p,s_subtype,s_type,s_partition,s_cell,n_neighbors,resolution,z_score,linkage,
               s_color_p='Platform',d_color_p = {'cycIF':'gold','IMC':'darkblue'},
               savedir=f'{codedir}/20220222/Survival_Plots_Both',figsize=(7,6)):
    #more plots
    #color by platform/leiden
    from matplotlib.pyplot import gcf
    d_color = dict(zip(sorted(adata.obs.leiden.unique()),sns.color_palette()[0:len(adata.obs.leiden.unique())]))
    
    network_colors = df_p.leiden.astype('str').map(d_color)#
    network_colors.name = 'cluster'
    node_colors  = df_p.loc[:,s_color_p].astype('str').map(d_color_p)
    network_node_colors = pd.DataFrame(node_colors).join(pd.DataFrame(network_colors))

    g = sns.clustermap(df_p.loc[:,ls_col].dropna(),figsize=figsize,cmap='viridis',
            row_colors=network_node_colors,method=linkage,dendrogram_ratio=0.16)
    for label,color in d_color_p.items():
        g.ax_col_dendrogram.bar(0, 0, color=color,label=label, linewidth=0)
    l1 = g.ax_col_dendrogram.legend(loc="right", ncol=1,bbox_to_anchor=(-0.1, 0.72),bbox_transform=gcf().transFigure)
    for label,color in d_color.items():
        g.ax_row_dendrogram.bar(0, 0, color=color,label=label, linewidth=0)
    l2 = g.ax_row_dendrogram.legend(loc="right", ncol=1,bbox_to_anchor=(-0.1, 0.5),bbox_transform=gcf().transFigure)
    g.savefig(f'{savedir}/clustermap_PlatformandSubtype_{s_sample}_{s_type}_{s_partition}_{s_cell}_{s_type}_{n_neighbors}_{resolution}.png',dpi=200)

    #subtypes' mean
    d_replace = {}
    df_plot = df_p.loc[:,ls_col.tolist()+['leiden']].dropna().groupby('leiden').mean()
    df_plot.index.name = f'leiden {resolution}'
    g = sns.clustermap(df_plot.dropna().T,z_score=z_score,figsize=(4,len(ls_col)*.25+1),cmap='viridis',vmin=-2,vmax=2,method='ward')
    g.fig.suptitle(f'leiden {resolution}',x=.9) 
    g.savefig(f'{savedir}/clustermap_subtypes_{s_sample}_{s_type}_{s_partition}_{s_cell}_{s_type}_{n_neighbors}_{resolution}.png',dpi=200)
    marker_genes = df_plot.dropna().T.iloc[:,g.dendrogram_col.reordered_ind].columns.tolist()
    categories_order = df_plot.dropna().T.iloc[g.dendrogram_row.reordered_ind,:].index.tolist()
    #barplot
    fig,ax=plt.subplots(figsize=(2.5,2.5),dpi=200)
    df_p.groupby(['leiden','Platform','Subtype']).count().iloc[:,0].unstack().loc[marker_genes].plot(kind='barh',title='Patient Count',ax=ax)
    plt.tight_layout()
    fig.savefig(f'{savedir}/barplot_subtyping_{s_sample}_{s_type}_{s_partition}_{s_cell}_{s_type}_{n_neighbors}_{resolution}.png')

def group_mean_diff(df_marker,s_group,s_marker):
    lls_result = []
    for s_test in df_marker.loc[:,s_group].dropna().unique():
        ls_result = df_marker.loc[df_marker.loc[:,s_group] == s_test,s_marker].values
        lls_result.append(ls_result)
    if len(lls_result)==2:
        try:
            statistic, pvalue = stats.ttest_ind(lls_result[0],lls_result[1])
        except:
            print('error in ttest_inde')
            pvalue = 1.0
            statistic = None
    elif len(lls_result) > 2:
        try:
            statistic, pvalue = stats.f_oneway(*lls_result,nan_policy='omit')
        except:
            print('error in f_oneway')
            pvalue = 1.0
            statistic = None
    else:
        pvalue = 1.0
        statistic = None
    return(statistic, pvalue)

def quartile_km(df,s_col,s_title_str='',savedir='',alpha=0.05,i_cut=4,labels=['low','med-low','med-high','high'],s_time='Survival_time',s_censor='Survival'):
    '''
    make sure labels has a high and low, other names will be ignored
    s_title_str additional label for title
    '''
    df = df.loc[:,[s_col,s_time,s_censor]].dropna()
    if len(df) > 1:
        #KM
        q = pd.qcut(df.loc[:,s_col], q=i_cut,labels=labels,duplicates='drop').str.replace('X','') 
        s_title1 = f'{s_col}'
        s_title2 = f'quantiles={i_cut}'
        df.loc[q=='low','abundance'] = 'low'
        df.loc[q=='high','abundance'] = 'high'
        #log rank
        results = multivariate_logrank_test(event_durations=df.loc[:,s_time],
                                            groups=df.abundance, event_observed=df.loc[:,s_censor])
        pvalue = results.summary.p[0]
        if np.isnan(pvalue):
            print(f'{s_col}: pvalue is na')
            pvalue = 1
        if pvalue < alpha:
            kmf = KaplanMeierFitter()
            fig, ax = plt.subplots(figsize=(3,3),dpi=300)
            for s_group in ['high','low']:
                df_abun = df[df.abundance==s_group]
                durations = df_abun.loc[:,s_time]
                event_observed = df_abun.loc[:,s_censor]
                try:
                    kmf.fit(durations, event_observed,label=s_group)
                    kmf.plot(ax=ax,ci_show=False,show_censors=True)
                except:
                    results.summary.p[0] = 1
            s_pval = f'{results.summary.p[0]:.2}'
            ax.set_title(f'{s_title1}\n{s_title_str} {s_title2}\np={s_pval} n={len(df)} [{sum(df.abundance=="high")}, {sum(df.abundance=="low")}]',fontsize=10)
            ax.set_xlabel(s_time)
            ax.legend(loc='upper right')
            plt.tight_layout()
            fig.savefig(f"{savedir}/KM_{s_title1.replace(' ','_')}_{s_title_str}_{s_title2.replace(' ','_')}_{i_cut}_{s_censor}_{s_pval}.png",dpi=300)
            #plt.close(fig)
    else:
        print(f'{s_col}: too many nas')
        df = pd.DataFrame(columns=[s_time,s_censor,'abundance'],index=[0],data=np.nan)
        pvalue = 1
    return(df, pvalue)


## find best cutpoint 
def single_km(df_all,s_cell,s_subtype,s_plat,s_col,savedir,alpha=0.05,cutp=0.5,s_time='Survival_time',s_censor='Survival',s_propo='in'):
    df_all.index = df_all.index.astype('str')
    if s_subtype != '':
        df = df_all[(df_all.subtype==s_subtype)].copy() #(df_all.Platform==s_plat) &
    else:
        df = df_all.copy()
    df = df.loc[:,[s_col,s_time,s_censor]].dropna()
    if len(df) > 1:
        #KM
        i_cut = np.quantile(df.loc[:,s_col],cutp)
        b_low = df.loc[:,s_col] <= i_cut
        s_title1 = f'{s_col} {s_plat}'
        s_title2 = f'{s_propo} {s_cell}'
        if i_cut == 0:
            b_low = df.loc[:,s_col] <= 0
        df.loc[b_low,'abundance'] = 'low'
        df.loc[~b_low,'abundance'] = 'high'
        #log rank
        
        results = multivariate_logrank_test(event_durations=df.loc[:,s_time],
                                            groups=df.abundance, event_observed=df.loc[:,s_censor])
        pvalue = results.summary.p[0]
        if np.isnan(pvalue):
            print(f'{s_col}: pvalue is na')
            pvalue = 1
        if pvalue < alpha:
            kmf = KaplanMeierFitter()
            fig, ax = plt.subplots(figsize=(3,3),dpi=300)
            for s_group in ['high','low']:
                df_abun = df[df.abundance==s_group]
                durations = df_abun.loc[:,s_time]
                event_observed = df_abun.loc[:,s_censor]
                try:
                    kmf.fit(durations, event_observed,label=s_group)
                    kmf.plot(ax=ax,ci_show=False,show_censors=True)
                except:
                    results.summary.p[0] = 1
            s_pval = f'{results.summary.p[0]:.2}'
            ax.set_title(f'{s_title1}\n{s_title2} {s_subtype}\np={s_pval} n={len(df)} [{len(df.loc[~b_low,:])}, {len(df.loc[b_low,:])}]',fontsize=10)
            ax.set_xlabel(s_time)
            ax.legend(loc='upper right',title=f'<={i_cut:.2}')#,title=f'{len(df.loc[~b_low,:])}, {len(df.loc[b_low,:])}'
            plt.tight_layout()
            fig.savefig(f"{savedir}/Survival_Plots/KM_{s_title1.replace(' ','_')}_{s_title2.replace(' ','_')}_{s_subtype}_{cutp}_{s_censor}_{s_pval}.pdf",dpi=300)
            #plt.close(fig)
    else:
        print(f'{s_col}: too many nas')
        df = pd.DataFrame(columns=[s_time,s_censor,'abundance'],index=[0],data=np.nan)
        pvalue = 1
    return(df, pvalue)

warnings.filterwarnings("default",category = exceptions.ApproximationWarning)

# def cph_plot(df,s_multi,s_time,s_censor,figsize=(3,3)):
#     cph = CoxPHFitter()  #penalizer=0.1
#     if df.columns.isin(['Stage']).any():
#         df.Stage = df.Stage.replace({'I':1,'II':2,'III':3,'IV':4})
#     cph.fit(df.dropna(), s_time, event_col=s_censor) 
#     fig, ax = plt.subplots(figsize=figsize,dpi=300)
#     cph.plot(ax=ax)
#     pvalue = cph.summary.p[s_multi]
#     ax.set_title(f'{s_multi}\np={pvalue:.2} n={(len(df.dropna()))}')
#     plt.tight_layout()
#     return(fig,cph)

# def km_plot(df,s_col,s_time,s_censor):
#     results = multivariate_logrank_test(event_durations=df.loc[:,s_time],
#                                     groups=df.loc[:,s_col], event_observed=df.loc[:,s_censor])        
#     kmf = KaplanMeierFitter()
#     fig, ax = plt.subplots(figsize=(4,4),dpi=300)
#     ls_order = sorted(df.loc[:,s_col].dropna().unique())
#     for s_group in ls_order:
#         print(s_group)
#         df_abun = df[df.loc[:,s_col]==s_group]
#         durations = df_abun.loc[:,s_time]
#         event_observed = df_abun.loc[:,s_censor]
#         kmf.fit(durations,event_observed,label=s_group)
#         kmf.plot(ax=ax,ci_show=True,show_censors=True)
#     ax.set_title(f'{s_col}\np={results.summary.p[0]:.2} n={[df.loc[:,s_col].value_counts()[item] for item in ls_order]}')
#     ax.set_ylim(-0.05,1.05)
#     return(fig,ls_order)

def patient_heatmap(df_p,ls_col,ls_annot,figsize=(7,6),method="complete",metric="euclidean", standard_scale=0,
                    ls_color=[mpl.cm.tab10.colors,mpl.cm.Set1.colors,mpl.cm.Set2.colors,mpl.cm.Set3.colors,mpl.cm.Paired.colors,mpl.cm.Pastel1.colors],
                    ):
    #more plots
    #color by platform/leiden
    from matplotlib.pyplot import gcf

    df_annot = pd.DataFrame()
    dd_color = {}
    for idx, s_annot in enumerate(ls_annot):
        color_palette = ls_color[idx]
        d_color = dict(zip(sorted(df_p.loc[:,s_annot].dropna().unique()),color_palette[0:len(df_p.loc[:,s_annot].dropna().unique())]))
        network_colors = df_p.loc[:,s_annot].map(d_color) 
        df_annot[s_annot] = pd.DataFrame(network_colors)
        dd_color.update({s_annot:d_color})
    #try:
    g = sns.clustermap(df_p.loc[:,ls_col],figsize=figsize,cmap='viridis',method=method,metric=metric,standard_scale=standard_scale,#z_score=z_score,
        row_colors=df_annot,dendrogram_ratio=0.16,xticklabels=1,yticklabels=1,
        cbar_pos= (-.05, .8, 0.05, 0.1),
        cbar_kws= {'orientation':'vertical','label':'z-score'}#'fraction':0.01, 'pad':0.04,'anchor':(-.01,-.01),
                   #'aspect':20,'fraction':.03,'shrink':0.5}
                   )
    for idx, (s_annot, d_color) in enumerate(dd_color.items()):
        g.ax_col_dendrogram.bar(0, 0, color='w',label=' ', linewidth=0)
        for label,color in d_color.items():
            g.ax_col_dendrogram.bar(0, 0, color=color,label=label, linewidth=0)
    
    l1 = g.ax_col_dendrogram.legend(loc="right", ncol=1,bbox_to_anchor=(0, 0.6),bbox_transform=gcf().transFigure)
    # except:
    #     print('clustermap error')
    #     g= df_p.loc[:,ls_col].dropna(how='any')
    return(g,df_annot)
    #g.savefig(f'{savedir}/clustermap_PlatformandSubtype_{s_sample}_{s_type}_{s_partition}_{s_cell}_{s_type}_{n_neighbors}_{resolution}.png',dpi=200)

    
#functions

import matplotlib
def df_from_mcomp(m_comp):
    df_test = pd.DataFrame.from_records(m_comp.summary().data,coerce_float=True)
    df_test.columns=df_test.loc[0].astype('str')
    df_test.drop(0,inplace=True)
    df_test =df_test.apply(pd.to_numeric, errors='ignore')
    ls_order = (pd.concat([df_test.group1,df_test.group2])).unique()
    return(df_test, ls_order)
def plt_sig(df_test,ax,ax_factor=5):
    ls_order = pd.concat([df_test.group1,df_test.group2]).unique()
    props = {'connectionstyle':matplotlib.patches.ConnectionStyle.Bar(armA=0.0, armB=0.0, fraction=0.0, angle=None),
             'arrowstyle':'-','linewidth':.5}
    #draw on axes
    y_lim = ax.get_ylim()[1]
    y_lim_min = ax.get_ylim()[0]
    y_diff = y_lim-y_lim_min
    for count, s_index in enumerate(df_test[df_test.reject].index):
        text =f"p = {df_test.loc[s_index,'p-adj']:.1}"
        #text = "*"
        one = df_test.loc[s_index,'group1']
        two = df_test.loc[s_index,'group2']
        x_one = np.argwhere(ls_order == one)[0][0]
        x_two = np.argwhere(ls_order == two)[0][0]
        ax.annotate(text, xy=(np.mean([x_one,x_two]),y_lim - (y_diff+count)/ax_factor),fontsize=6)
        ax.annotate('', xy=(x_one,y_lim - (y_diff+count)/ax_factor), xytext=(x_two,y_lim - (y_diff+count)/ax_factor), arrowprops=props)
        #break
    return(ax)

def plot_chi(df_all,s_clin,s_x_var,figsize=(4,4),alpha=0.05,annot=None):
    df_patient = df_all.loc[:,[s_clin,s_x_var]].dropna()
    crosstab = pd.crosstab(df_patient.loc[:,s_clin],df_patient.loc[:,s_x_var])
    df = crosstab.T.unstack().reset_index().rename({0:'Total Pts'},axis=1)
    statistic,pvalue, dof, expected_freq = stats.chi2_contingency(crosstab)
    if pvalue < alpha:
        fig, ax = plt.subplots(figsize=figsize,dpi=300)
        if not annot is None:
            sns.heatmap(crosstab/crosstab.sum(),ax=ax,annot=annot,cbar_kws={'label':f'Percent Pts. in {s_x_var}'},fmt='')
        else:
            sns.heatmap(crosstab/crosstab.sum(),ax=ax,annot=True,cbar_kws={'label':f'Percent Pts. in {s_x_var}'})
        ax.set_title(f'{s_clin} vs. {s_x_var}\np={pvalue:.3} n={len(df_patient)}')
    else:
        fig = None
        ax = None
    return(fig,ax,crosstab)
    
def post_hoc(confusion_matrix):
    chi2, pvalue, dof, expected  = stats.chi2_contingency(confusion_matrix)
    observed_vals = confusion_matrix
    expected_vals = pd.DataFrame(expected,index=confusion_matrix.index,columns=confusion_matrix.columns)
    result_val = pd.DataFrame(data='',index=confusion_matrix.index,columns=confusion_matrix.columns)
    col_sum = observed_vals.sum(axis=1)
    row_sum = observed_vals.sum(axis=0)

    for indx in confusion_matrix.index:
        for cols in confusion_matrix.columns:
            observed = float(observed_vals.loc[indx,cols])
            expected = float(expected_vals.loc[indx,cols])
            col_total = float(col_sum[indx])
            row_total = float(row_sum[cols])
            expected_row_prop = expected/row_total
            expected_col_prop = expected/col_total
            std_resid = (observed - expected) / (math.sqrt(expected * (1-expected_row_prop) * (1-expected_col_prop)))
            p_val = norm.sf(abs(std_resid))
            if p_val < 0.05/(len(confusion_matrix.index)*len(confusion_matrix.columns)):
                print(indx,cols, "***", p_val)
                result_val.loc[indx,cols] = '***'
            elif p_val < 0.05:
                print (indx,cols, '*', p_val)
                result_val.loc[indx,cols] = '*'
            else:
                print (indx,cols, 'not sig', p_val)
    print('cutoff')
    print(0.05/(len(confusion_matrix.index)*len(confusion_matrix.columns)))
    return(result_val)

def single_var_km_cph(df_all,df_surv,s_subtype,s_platform,s_cell,alpha=0.05,min_cutoff=0.003,savedir=f"/home/groups/graylab_share/OMERO.rdsStore/engje/Data/20200000/20200406_JP-TMAs/20220408/Survival_Plots"):
    df_all.index = df_all.index.astype('str')
    df_surv.index = df_surv.index.astype('str')
    df_all = df_all.merge(df_surv.loc[:,['Survival','Survival_time','subtype','Platform']],left_index=True,right_index=True)
    if s_platform == 'IMC':
        df = df_all[(df_all.Platform==s_platform) & (~df_all.index.str.contains('Z')) & (df_all.subtype==s_subtype)].copy()
    elif s_platform == 'cycIF':
        df = df_all[(df_all.Platform==s_platform) & (~df_all.index.str.contains('JP-TMA2')) & (df_all.subtype==s_subtype)].copy()
    else:
        df = df_all[(df_all.Platform==s_platform) & (df_all.subtype==s_subtype)].copy()
    df = df.dropna() #df.dropna(axis=1).dropna()
    #KM
    for s_col in df.columns.drop(['Survival','Survival_time','subtype','Platform']):
        b_low = df.loc[:,s_col] <= df.loc[:,s_col].median()
        s_title1 = f'{s_subtype} {s_platform}'
        s_title2 = f'{s_cell} {s_col.replace(".","")}'
        if df.loc[:,s_col].median() < min_cutoff:
            continue
        elif len(df) < 1:
            continue
        df.loc[b_low,'abundance'] = 'low'
        df.loc[~b_low,'abundance'] = 'high'
        #log rank
        results = multivariate_logrank_test(event_durations=df.Survival_time,
                                            groups=df.abundance, event_observed=df.Survival)
        if results.summary.p[0] < alpha:
            print(s_col)
            #kaplan meier plotting
            kmf = KaplanMeierFitter()
            fig, ax = plt.subplots(figsize=(3,3),dpi=300)
            for s_group in ['high','low']:
                df_abun = df[df.abundance==s_group]
                durations = df_abun.Survival_time
                event_observed = df_abun.Survival
                try:
                    kmf.fit(durations, event_observed,label=s_group)
                    kmf.plot(ax=ax,ci_show=False,show_censors=True)
                except:
                    print('.')
            ax.set_title(f'{s_title1}\n{s_title2}\np={results.summary.p[0]:.2} (n={len(df)})',fontsize=10)
            ax.legend(loc='upper right',title=f'{df.loc[:,s_col].median():.2}')
            plt.tight_layout()
            fig.savefig(f"{savedir}/KM_{s_title1.replace(' ','_')}_{s_title2.replace(' ','_')}.png",dpi=300)
        #CPH
        cph2 = CoxPHFitter(penalizer=0.1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                cph2.fit(df.loc[:,[s_col,'Survival_time','Survival']], duration_col='Survival_time', event_col='Survival')
                if cph2.summary.p[0] < alpha:
                    print(s_col)
                    fig, ax = plt.subplots(figsize=(2.5,2),dpi=300)
                    cph2.plot(ax=ax)
                    ax.set_title(f'{s_title1} (n={len(df)})\n{s_title2}\np={cph2.summary.p[0]:.2} ({df.loc[:,s_col].median():.2})',fontsize=10)
                    ax.set_ylabel(f'{s_col}')
                    ax.set_yticklabels([])
                    plt.tight_layout()
                    fig.savefig(f"{savedir}/CPH_{s_title1.replace(' ','_')}_{s_title2.replace(' ','_')}.png",dpi=300)
            except:
                print(f'skipped {s_col}')   
    return(df)


                
def make_adata(df,ls_col,n_neighbors):
    print('making adata')
    adata = sc.AnnData(df.loc[:,ls_col].fillna(0))
    adata.raw = adata
    #reduce dimensionality
    sc.tl.pca(adata, svd_solver='auto')
    print('scaling')
    sc.pp.scale(adata, zero_center=False, max_value=20)
    print('calc umap')
    # calculate neighbors 
    sc.pp.neighbors(adata, n_neighbors=n_neighbors) 
    sc.tl.umap(adata)
    return(adata)

def km_cph(df_st):
    print(len(df_st)) 
    T = df_st['Survival_time']     ## time to event
    E = df_st['Survival']      ## event occurred or censored
    groups = df_st.loc[:,'leiden'] 
    kmf1 = KaplanMeierFitter() ## instantiate the class to create an object
    fig1, ax1 = plt.subplots(figsize=(3,3),dpi=200)
    for idx, s_group in enumerate(sorted(df_st.leiden.unique())):
        i1 = (groups == s_group)
        if sum(i1) > 0:
            kmf1.fit(T[i1], E[i1], label=s_group)    ## fit thedata
            kmf1.plot(ax=ax1,ci_show=False,color=mpl.cm.Set1.colors[idx],show_censors=True)
            print(f'{s_group}: {kmf1.median_survival_time_}, ({i1.sum()})') #{kmf1.percentile(.75)} 
    results = multivariate_logrank_test(event_durations=T, groups=groups, event_observed=E)
    p1=results.summary.p[0]
    ax1.set_title(f'p={p1:.2}')
    #CPH
    df_dummy = pd.get_dummies(df_st.loc[:,['Survival_time','Survival','leiden']])
    df_dummy = df_dummy.loc[:,df_dummy.sum() != 0]
    cph = CoxPHFitter(penalizer=0.1)  ## Instantiate the class to create a cph object
    cph.fit(df_dummy, 'Survival_time', event_col='Survival')
    fig2, ax2 = plt.subplots(figsize=(2.5,3),dpi=200)
    cph.plot(ax=ax2)
    p2 = cph.summary.loc[:,'p'].min()
    ax2.set_title(f'p={p2:.2}')
    return(fig1, p1, fig2, p2)

def make_adata_old(df, ls_col,df_surv, n_neighbors, s_subtype, s_type, s_partition, s_cell,ncols=4):
    print('making adata')
    adata = sc.AnnData(df.loc[:,ls_col].fillna(0))
    adata.raw = adata
    #reduce dimensionality
    sc.tl.pca(adata, svd_solver='auto')
    print('scaling')
    sc.pp.scale(adata, zero_center=False, max_value=20)
    print('calc umap')
    # calculate neighbors 
    sc.pp.neighbors(adata, n_neighbors=n_neighbors) 
    sc.tl.umap(adata)
    #color by markers   
    figname = f"Umapboth_markers_{s_subtype}_{s_type}_{s_partition}_{s_cell}_{n_neighbors}neigh.png"
    title=figname.split('.png')[0].replace('_',' ')
    sc.pl.umap(adata, color=ls_col,vmin='p1.5',vmax='p99.5',ncols=ncols,save=figname,size=250)
    #platform
    adata.obs['Platform'] = adata.obs.index.astype('str').map(dict(zip(df_surv.index.astype('str'),df_surv.Platform)))
    figname = f"Umapboth_Platform_{s_subtype}_{s_type}_{s_partition}_{s_cell}_{n_neighbors}neigh.png"
    title=figname.split('.png')[0].replace('_',' ')
    sc.pl.umap(adata, color='Platform',save=figname,size=250)
    #subtype
    adata.obs['subtype'] = adata.obs.index.astype('str').map(dict(zip(df_surv.index.astype('str'),df_surv.subtype)))
    #CAREFUL
    adata.obs['subtype'] = adata.obs['subtype'].fillna('TNBC')
    figname = f"Umapboth_subtype_{s_subtype}_{s_type}_{s_partition}_{s_cell}_{n_neighbors}neigh.png"
    title=figname.split('.png')[0].replace('_',' ')
    sc.pl.umap(adata, color='subtype',save=figname,size=250)
    return(adata)
def cluster_leiden_old(adata, resolution,n_neighbors, s_subtype, s_type, s_partition, s_cell):
    sc.tl.leiden(adata,resolution=resolution)
    fig,ax = plt.subplots(figsize=(2.5,2),dpi=200)
    figname=f'both_{s_subtype}_{s_partition}_{s_cell}_{n_neighbors}_{resolution}.png'
    sc.pl.umap(adata, color='leiden',ax=ax,title=figname.split('.png')[0].replace('_',' '),wspace=.25,save=figname,size=40)
    return(adata)
def km_cph_old(adata,df_surv,s_subtype,s_plat,s_type,s_partition,s_cell,savedir=f'{codedir}/20220222/Survival_Plots_Both'):
    if type(adata) == anndata._core.anndata.AnnData:
        df_p = pd.DataFrame(data=adata.raw.X, index=adata.obs.index, columns=adata.var.index) #adata.to_df()
        df_p['Subtype'] = adata.obs.subtype
        df_p['leiden'] = adata.obs.leiden
        df_p['Platform'] = adata.obs.Platform
    else:
        df_p = adata
    df_p.index = df_p.index.astype('str')
    df_p['Survival'] = df_p.index.map(dict(zip(df_surv.index,df_surv.Survival)))
    df_p['Survival_time'] = df_p.index.map(dict(zip(df_surv.index,df_surv.Survival_time)))
    df_st = df_p[(df_p.Subtype==s_subtype)].dropna()
    if s_plat != 'Both':
        df_st = df_p[(df_p.Platform==s_plat) & (df_p.Subtype==s_subtype)].dropna()
    if not len(df_st) < 1:
        
        print(len(df_st)) 
        T = df_st['Survival_time']     ## time to event
        E = df_st['Survival']      ## event occurred or censored
        groups = df_st.loc[:,'leiden'] 
        kmf1 = KaplanMeierFitter() ## instantiate the class to create an object
        fig, ax = plt.subplots(figsize=(3,3),dpi=200)
        for idx, s_group in enumerate(sorted(df_p.leiden.unique())):
            i1 = (groups == s_group)
            if sum(i1) > 0:
                kmf1.fit(T[i1], E[i1], label=s_group)    ## fit thedata
                kmf1.plot(ax=ax,ci_show=False,color=f'C{idx}',show_censors=True)
                print(f'{s_group}: {kmf1.median_survival_time_}, {kmf1.percentile(.75)} ({i1.sum()})')
        results = multivariate_logrank_test(event_durations=T, groups=groups, event_observed=E)
        ax.set_title(f'{s_subtype} {s_plat} {s_cell} \nk={resolution}  p={results.summary.p[0]:.1} n={len(df_st)}')
        ax.legend(loc='upper right')
        ax.set_ylim(-0.05,1.05)
        plt.tight_layout()
        fig.savefig(f'{savedir}/KM_{s_subtype}_{s_plat}_{s_type}_{s_partition}_{s_cell}_{n_neighbors}_{resolution}.png',dpi=300)
        #CPH
        df_dummy = pd.get_dummies(df_st.loc[:,['Survival_time','Survival','leiden']])
        df_dummy = df_dummy.loc[:,df_dummy.sum() != 0]
        cph = CoxPHFitter(penalizer=0.1)  ## Instantiate the class to create a cph object
        cph.fit(df_dummy, 'Survival_time', event_col='Survival')
        fig, ax = plt.subplots(figsize=(2.5,3),dpi=200)
        cph.plot(ax=ax)
        pvalue = cph.summary.loc[:,'p'].min()
        ax.set_title(f'CPH: {s_subtype} {s_plat} {s_cell}\np={pvalue:.2}')
        plt.tight_layout()
        fig.savefig(f'{savedir}/CoxPH_{s_subtype}_{s_plat}_{s_type}_{s_partition}_{s_cell}_{n_neighbors}_{resolution}.png',dpi=300)
    else:
        cph = 'zero'
    return(df_p, cph)

def km_cph_entropy(df_p,df,ls_col,s_subtype,s_plat,s_cell,savedir=f'{codedir}/20220222/Survival_Plots_Both'):
    df_p['entropy'] = entropy(df_p.loc[:,df_p.columns[df_p.dtypes=='float32']].fillna(0),axis=1,base=2)
    df_st = df_p[(df_p.Subtype==s_subtype)].dropna()
    if s_plat != 'Both':
        df_st = df_p[(df_p.Platform==s_plat) & (df_p.Subtype==s_subtype)].dropna()
    #######3 Entropy
    s_col = 'entropy'
    # no df and ls_col variable
    df_st = df.loc[:,ls_col].merge(df_st.loc[:,['Subtype','Platform','Survival','Survival_time','entropy']],left_index=True,right_index=True)
    if not len(df_st) < 1:
        b_low = df_st.loc[:,s_col] <= df_st.loc[:,s_col].median()
        if df_st.loc[:,s_col].median() == 0:
            b_low = df.loc[:,s_col] <= 0
        df_st.loc[b_low,'abundance'] = 'low'
        df_st.loc[~b_low,'abundance'] = 'high'
        kmf = KaplanMeierFitter()
        results = multivariate_logrank_test(event_durations=df_st.Survival_time, groups=df_st.abundance, event_observed=df_st.Survival)
        print(f'entropy {results.summary.p[0]}')
        if results.summary.p[0] < 0.2:
            fig, ax = plt.subplots(figsize=(3,3),dpi=200)
            for s_group in ['high','low']:
                    df_abun = df_st[df_st.abundance==s_group]
                    durations = df_abun.Survival_time
                    event_observed = df_abun.Survival
                    kmf.fit(durations, event_observed,label=s_group)
                    kmf.plot(ax=ax,ci_show=False,show_censors=True)
            s_title1 = f'{s_subtype} {s_plat}'
            s_title2 = f'{s_cell} {s_col}'
            ax.set_title(f'{s_title1}\n{s_title2}\np={results.summary.p[0]:.2}',fontsize=10)
            ax.legend(loc='upper right')
            plt.tight_layout()
            fig.savefig(f"{savedir}/KM_{s_title1.replace(' ','_')}_{s_title2.replace(' ','_')}.png",dpi=300)
            cph = CoxPHFitter(penalizer=0.1)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    cph.fit(df_st.loc[:,[s_col,'Survival','Survival_time']], duration_col='Survival_time', event_col='Survival')
                    if cph.summary.p[0] < 0.1:
                        print(s_col)
                        fig, ax = plt.subplots(figsize=(2.5,2),dpi=200)
                        cph.plot(ax=ax)
                        s_title1 = f'{s_subtype} {s_plat}'
                        s_title2 = f'{s_cell} {s_col}'
                        ax.set_title(f'{s_title1}\n{s_title2}\np={cph.summary.p[0]:.2}',fontsize=10)
                        plt.tight_layout()
                        fig.savefig(f"{savedir}/CPH_{s_title1.replace(' ','_')}_{s_title2.replace(' ','_')}.png",dpi=300)
                except:
                    print(f'skipped {s_col}')   


def group_median_diff(df_marker,s_group,s_marker):
    lls_result = []
    for s_test in df_marker.loc[:,s_group].unique():
        ls_result = df_marker.loc[df_marker.loc[:,s_group] == s_test,s_marker].values
        lls_result.append(ls_result)
    if len(lls_result)==2:
        try:
            statistic, pvalue = stats.mannwhitneyu(lls_result[0],lls_result[1])
        except:
            pvalue = 1
            statistic = None
    elif len(lls_result) > 2:
        try:
            statistic, pvalue = stats.kruskal(*lls_result,nan_policy='omit')
        except:
            pvalue = 1
            statistic = None
    else:
        pvalue = None
        statistic = None
    #print(pvalue)
    return(statistic,pvalue)


#functions
def silheatmap(adata,clust,marker_list,sil_key):
    cluster_list = [str(item) for item in adata.uns[f'dendrogram_{clust}']['categories_ordered']]
    #dataframe
    df = adata.to_df()
    df[clust] = adata.obs[clust]
    #sort by sil
    df[sil_key] = adata.obs[sil_key]
    df = df.sort_values(by=sil_key)
    #sort by cluster, markers
    df['old_index'] = df.index
    obs_tidy = df.set_index(clust)
    obs_tidy.index = obs_tidy.index.astype('str')
    obs_tidy = obs_tidy.loc[cluster_list,:]
    df = df.loc[obs_tidy.old_index]
    obs_tidy = obs_tidy.loc[:,marker_list]
    #scale
    obs_tidy = pd.DataFrame(data=minmax_scale(obs_tidy),index=obs_tidy.index,columns=obs_tidy.columns)
    # define a layout of 3 rows x 3 columns
    # The first row is for the dendrogram (if not dendrogram height is zero)
    # second row is for main content. This col is divided into three axes:
    #   first ax is for the heatmap
    #   second ax is for 'brackets' if any (othwerise width is zero)
    #   third ax is for colorbar
    colorbar_width = 0.2
    var_names = marker_list
    width = 10
    dendro_height = 0.8 #if dendrogram else 0
    groupby_height = 0.13 #if categorical else 0
    heatmap_height = len(var_names) * 0.18 + 1.5
    height = heatmap_height + dendro_height + groupby_height + groupby_height
    height_ratios = [dendro_height, heatmap_height, groupby_height,groupby_height]
    width_ratios = [width, 0, colorbar_width, colorbar_width]
    fig = plt.figure(figsize=(width, height),dpi=200)
    axs = gridspec.GridSpec(
        nrows=4,
        ncols=4,
        wspace=1 / width,
        hspace=0.3 / height,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )
    norm = mpl.colors.Normalize(vmin=0, vmax=1, clip=False)
    norm2 = mpl.colors.Normalize(vmin=-1, vmax=1, clip=False)

    # plot heatmap
    heatmap_ax = fig.add_subplot(axs[1, 0])
    im = heatmap_ax.imshow(obs_tidy.T.values, aspect='auto',norm=norm,interpolation='nearest') # ,interpolation='nearest'
    heatmap_ax.set_xlim(0 - 0.5, obs_tidy.shape[0] - 0.5)
    heatmap_ax.set_ylim(obs_tidy.shape[1] - 0.5, -0.5)
    heatmap_ax.tick_params(axis='x', bottom=False, labelbottom=False)
    heatmap_ax.set_xlabel('')
    heatmap_ax.grid(False)
    heatmap_ax.tick_params(axis='y', labelsize='small', length=1)
    heatmap_ax.set_yticks(np.arange(len(var_names)))
    heatmap_ax.set_yticklabels(var_names, rotation=0)

    #colors
    value_sum = 0
    ticks = []  # list of centered position of the labels
    labels = []
    label2code = {}  # dictionary of numerical values asigned to each label
    for code, (label, value) in enumerate(
            obs_tidy.index.value_counts().loc[cluster_list].iteritems()
        ):
            ticks.append(value_sum + (value / 2))
            labels.append(label)
            value_sum += value
            label2code[label] = code

    groupby_cmap = mpl.colors.ListedColormap(adata.uns[f'{clust}_colors'])
    groupby_ax = fig.add_subplot(axs[3, 0])
    groupby_ax.imshow(
                np.array([[label2code[lab] for lab in obs_tidy.index]]),
                aspect='auto',
                cmap=groupby_cmap,
            )
    groupby_ax.grid(False)
    groupby_ax.yaxis.set_ticks([])
    groupby_ax.set_xticks(ticks,labels,fontsize='xx-small',rotation=90)
    groupby_ax.set_ylabel('Cluster',fontsize='x-small',rotation=0,ha='right',va='center')


    #sil
    sil_ax = fig.add_subplot(axs[2, 0])
    #max_index = df[sil_key].idxmax()    #df.loc[max_index,sil_key] = 1    #min_index = df[sil_key].idxmin()    #df.loc[min_index,sil_key] = -1 #not needed
    a=np.array([df[sil_key]]) #f'{clust}_silhuette'
    a_tile = np.tile(a,(int(len(df)/80),1))
    sil_ax.imshow(a_tile,cmap='bwr',norm=norm2)
    sil_ax.xaxis.set_ticks([])
    sil_ax.yaxis.set_ticks([])
    sil_ax.set_ylabel('Silhouette',fontsize='x-small',rotation=0,ha='right',va='center')
    sil_ax.grid(False)

    #dendrogram
    dendro_ax = fig.add_subplot(axs[0, 0], sharex=heatmap_ax)
    #_plot_dendrogram(dendro_ax, adata, groupby, dendrogram_key=dendrogram,ticks=ticks, orientation='top', )
    dendro_info = adata.uns[f'dendrogram_{clust}']['dendrogram_info']
    leaves = dendro_info["ivl"]
    icoord = np.array(dendro_info['icoord'])
    dcoord = np.array(dendro_info['dcoord'])
    orig_ticks = np.arange(5, len(leaves) * 10 + 5, 10).astype(float)
    for xs, ys in zip(icoord, dcoord):
        if ticks is not None:
            xs = translate_pos(xs, ticks, orig_ticks)
        dendro_ax.plot(xs, ys, color='#555555')
    dendro_ax.tick_params(bottom=False, top=False, left=False, right=False)
    ticks = ticks if ticks is not None else orig_ticks
    dendro_ax.set_xticks(ticks)
    #dendro_ax.set_xticklabels(leaves, fontsize='small', rotation=90)
    dendro_ax.set_xticklabels([])
    dendro_ax.tick_params(labelleft=False, labelright=False)
    dendro_ax.grid(False)
    dendro_ax.spines['right'].set_visible(False)
    dendro_ax.spines['top'].set_visible(False)
    dendro_ax.spines['left'].set_visible(False)
    dendro_ax.spines['bottom'].set_visible(False)

    # plot colorbar
    cbar_ax = fig.add_subplot(axs[1, 2])
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap='viridis')
    cbar = plt.colorbar(mappable=mappable, cax=cbar_ax)
    cbar_ax.tick_params(axis='both', which='major', labelsize='xx-small',rotation=90,length=.1)
    cbar_ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(locs=[0,1]))
    cbar.set_label('Expression', fontsize='xx-small',labelpad=-5)

    # plot colorbar2
    cbar_ax = fig.add_subplot(axs[1, 3])
    mappable = mpl.cm.ScalarMappable(norm=norm2, cmap='bwr')
    cbar = plt.colorbar(mappable=mappable, cax=cbar_ax)
    cbar_ax.tick_params(axis='both', which='major', labelsize='xx-small',rotation=90,length=.1)
    cbar_ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(locs=[-1,0,1]))
    cbar.set_label('Silhouette Score', fontsize='xx-small',labelpad=0)

    #return dict
    return_ax_dict = {'heatmap_ax': heatmap_ax}
    return_ax_dict['groupby_ax'] = groupby_ax
    return_ax_dict['dendrogram_ax'] = dendro_ax
    return(fig)

def translate_pos(pos_list, new_ticks, old_ticks):
    """
    transforms the dendrogram coordinates to a given new position.
    """
    # of given coordinates.

    if not isinstance(old_ticks, list):
        # assume that the list is a numpy array
        old_ticks = old_ticks.tolist()
    new_xs = []
    for x_val in pos_list:
        if x_val in old_ticks:
            new_x_val = new_ticks[old_ticks.index(x_val)]
        else:
            # find smaller and bigger indices
            idx_next = np.searchsorted(old_ticks, x_val, side="left")
            idx_prev = idx_next - 1
            old_min = old_ticks[idx_prev]
            old_max = old_ticks[idx_next]
            new_min = new_ticks[idx_prev]
            new_max = new_ticks[idx_next]
            new_x_val = ((x_val - old_min) / (old_max - old_min)) * (
                new_max - new_min
            ) + new_min
        new_xs.append(new_x_val)
    return new_xs

#functions

# count the neighbors.
class NeighborsCounter:

    def __init__(self, rad, xy=['CentroidX', 'CentroidY']):
        self.rad = rad
        self.xy = xy

    def query_balltree_vanilla(self, coords_np):
        """
        input coords_np:
            these are coordinates. possible shape: (N,2)

        output neighbor_indices:
            this is a list of lists.
            there is one list per row in coords_np (i.e. there are N)
            the i'th list contains the indices of the neighbors of i,
            not including itself.
        """
        n_points = coords_np.shape[0]
        print(f'Counting neighbors for {n_points} points.')

        tree = cKDTree(coords_np)
        neighbor_indices = tree.query_ball_tree(tree, self.rad)
        for i in range(n_points):
            neighbor_indices[i].remove(i)
        return neighbor_indices

    def run(self, dataframe):
        """
        Splits the input dataframe into cell types and coordinates
        Runs query_balltree_vanilla on the coordinates
        Uses the neighbor indices to get cell type neighbor counts.

        Input:
            a dataframe with boolean cell type columns and coordinate columns
            the coordinate columns by default are named ['CentroidX', 'CentroidY']
            (coordinate column names are stored in attribute self.xy)

        Output:
            a dataframe with the same shape and index as the input dataframe.
        """

        types = [c for c in dataframe.columns if c not in self.xy]
        #why do we have to do this?
        types.remove('slide')
        g = self.query_balltree_vanilla(dataframe[self.xy].to_numpy())
        counts = np.zeros((len(g), len(types)))
        df_arra = dataframe[types].to_numpy()

        #return(counts)
        
        for n in range(dataframe.shape[0]):
            idx = np.array(g[n])
            if idx.size:
                counts[n, :] = df_arra[idx, :].sum(axis=0)

        return pd.DataFrame(counts, index=dataframe.index, columns=types)
        

def plot_sil(d_sil,s_name='Tumor'):
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(dpi=200)
    pd.Series(d_sil).plot(ax=ax)
    ax.set_title(f'{s_name}: Mean Silhoutte Scores')
    ax.set_xlabel('k')
    plt.tight_layout()
    fig.savefig(f'{s_name}_Silhouette.png')
    
def single_km_cat(df_all,s_col,savedir,alpha=0.05,s_time='Survival_time',s_censor='Survival'):
    df_all.index = df_all.index.astype('str')
    df = df_all.copy()
    df = df.loc[:,[s_col,s_time,s_censor]].dropna()
    if len(df) > 0:
        #log rank
        results = multivariate_logrank_test(event_durations=df.loc[:,s_time],
                                            groups=df.loc[:,s_col], event_observed=df.loc[:,s_censor])
        #kaplan meier plotting
        if results.summary.p[0] < alpha:
            kmf = KaplanMeierFitter()
            fig, ax = plt.subplots(figsize=(3,3),dpi=300)
            for s_group in sorted(df.loc[:,s_col].unique()):
                df_abun = df[df.loc[:,s_col]==s_group]
                durations = df_abun.loc[:,s_time]
                event_observed = df_abun.loc[:,s_censor]
                try:
                    kmf.fit(durations, event_observed,label=s_group)
                    kmf.plot(ax=ax,ci_show=False,show_censors=True)
                except:
                    results.summary.p[0] = 1
            s_title1 = f'{s_col}'
            ax.set_title(f'{s_title1}\np={results.summary.p[0]:.2} (n={len(df)})',fontsize=10)
            ax.set_xlabel(s_time)
            ax.legend(loc='upper right')
            plt.tight_layout()
            #fig.savefig(f"{savedir}/Survival_Plots/KM_{s_title1.replace(' ','_')}_{s_censor}.png",dpi=300)
        else:
            fig = None
            ax = None
    else:
        fig = None
        ax = None
    return(fig, ax)

def add_patient_results(df_file, s_file,s_sample):
    if s_file.find(f'results_{s_sample}_GatedCellTypes') > -1:
        s_type = 'GatedCellTypes'
        s_subtype = s_file.split('.csv')[0].split('_')[-1]
        s_partition = 'gating'
        s_cell = s_file.split('.csv')[0].split('_')[-2].split('by')[1]
    elif s_file.find(f'results_{s_sample}_LeidenClustering_') > -1:
        s_type = 'LeidenClustering'
        s_subtype = s_file.split('.csv')[0].split('_')[-1]
        s_partition = s_file.split('.csv')[0].split('_')[-3].split('by')[1]
        s_cell = s_file.split('.csv')[0].split('_')[-2].split('in')[1]   
    elif s_file.find(f'results_{s_sample}_Abundance_') > -1:
        s_type = 'Abundance'
        s_subtype = s_file.split('.csv')[0].split('_')[-1]
        s_partition = s_file.split('.csv')[0].split('_')[-3].split('by')[1]
        s_cell = s_file.split('.csv')[0].split('_')[-2].split('in')[1]   
    elif s_file.find(f'results_{s_sample}_BGSubtractedMeanIntensity') > -1:
        s_type = 'MeanIntensity'
        s_subtype = s_file.split('.csv')[0].split('_')[-1]
        s_partition = s_file.split('.csv')[0].split('_')[-3].split('by')[1]
        s_cell = s_file.split('.csv')[0].split('_')[-2].split('in')[1]   
    elif s_file.find(f'results_{s_sample}_Density_') > -1:
        s_type = 'Density'
        s_subtype = s_file.split('.csv')[0].split('_')[-1]
        s_partition = s_file.split('.csv')[0].split('_')[-2].split('by')[1]
        s_cell = 'density'
    elif s_file.find(f'results_{s_sample}_FractionPositive') > -1:
        s_type = 'FractionPositive'
        s_subtype = s_file.split('.csv')[0].split('_')[-1]
        s_partition = s_file.split('.csv')[0].split('_')[-3].split('by')[1]
        s_cell = s_file.split('.csv')[0].split('_')[-2].split('in')[1]   
    elif s_file.find(f'Neighbors') > -1 and s_file.find(s_sample.split('_')[0] + '_' + s_sample.split('_')[1]) > -1:
            print('neighbors')
            s_type = s_file.split('_')[3]
            s_subtype = s_file.split('_')[3] #s_file.split('.csv')[0].split('_')[-1]
            s_partition = s_file.split('.csv')[0].split('_')[-2].split('by')[1]
            s_cell = s_file.split('.csv')[0].split('_')[-1].split('in')[1]   
    elif s_file.find(f'results_{s_sample}_') > -1:
        s_type = s_file.split('.csv')[0].split('_')[3]
        s_subtype = s_file.split('.csv')[0].split('_')[3]
        # results/results_20230307_U54-TMA-9_LeidenPatientSubtypes_Ripleys_Gcross.csv
        s_partition = s_file.split('.csv')[0].split('_')[4]#'epi'
        s_cell = s_file.split('.csv')[0].split('_')[-1]
        #try:
        #    s_cell = s_file.split('.csv')[0].split('_')[-1].split('in')[1]
        #except:
        #    s_cell = 'epi'
           
    else:
        return(df_file)
    df_file.loc[s_file,'subtype'] = s_subtype
    df_file.loc[s_file,'type'] = s_type
    df_file.loc[s_file,'partition'] = s_partition
    df_file.loc[s_file,'cell'] = s_cell
    return(df_file)

def load_patient_results(df_file, df_surv, s_index,  s_time, s_censor):
    s_type_one = df_file.loc[s_index,'type']
    #ls_type = re.findall('[A-Z][a-z]*', s_type_one)
    s_type = uppercase_to_space(s_type_one)
    s_cell = df_file.loc[s_index,'cell']
    s_measure = s_index.split('.csv')[0].split('_PDAC')[0].split('_')[1]
    if s_measure == 'Fraction':
        s_measure = s_cell
        s_cell = 'Fraction Positive in'
    if not s_measure.isnumeric():
        s_cell = f'{s_cell} {s_measure}'
    df_all=pd.read_csv(f'results/{s_index}',index_col=0)
    ls_marker = (df_all.columns[(df_all.dtypes=='float64') & (~df_all.columns.isin([s_time,s_censor]))]).tolist() #
    df_all['Old_Pt_ID'] = df_all.index
    df_all = df_all.merge(df_surv,on='Old_Pt_ID',how='left',suffixes=('_1',''))
    df_all = df_all[~df_all.Old_Pt_ID.duplicated(keep='first')]
    d_replace={'Ovary_Met':'Mets', 'Liver_Met':'Mets', 'Lung_Met':'Mets', 'Lymph node':'Mets',
       'Peritoneum_Met':'Mets'}
    df_all['subtype'] = df_all.subtype.replace(d_replace)
    return(ls_marker, df_all, s_type, s_cell)

def correlation_heatmap(df_all, s_title, dim = (8,7)):
    #just the dendrogram
    mask = (np.ones_like(df_all.corr().fillna(0))) 
    g = sns.clustermap(df_all.corr().fillna(0),yticklabels=1,figsize=dim,
                      cbar=False,dendrogram_ratio=0.1,cbar_pos=None,tree_kws={'linewidths':3},
                   mask=mask)
    plt.close()
    categories_order = df_all.corr().iloc[g.dendrogram_col.reordered_ind,:].index.tolist()
    df_all = df_all.loc[:,categories_order]
    rho = df_all.corr()
    pval = df_all.corr(method=lambda x, y: pearsonr(x, y)[1]) - np.eye(*rho.shape)
    p_vals = pval.applymap(lambda x: ''.join(['*' for t in [0.001,0.005,0.05] if x<=t]))
    fig, ax = plt.subplots(figsize=dim,dpi=300)
    sns.heatmap(rho, vmin=-1, vmax=1, annot=p_vals, fmt = '', cmap='RdBu_r',ax=ax,
               yticklabels=1,xticklabels=1,cbar_kws={'label':'Pearson Correlation','aspect':30,'shrink':0.8})
    ax.set_title(f'{s_title} Correlation n={len(df_all)}', fontdict={'fontsize':16}, pad=12);
    return(fig, g, rho, pval)

def caps_to_space(s_type_one):
    ls_type = re.findall('[A-Z][a-z]*', s_type_one)
    s_type = " ".join(ls_type)
    return(s_type)

def uppercase_to_space(inp):
    out = re.sub(r'(?<![A-Z\W])(?=[A-Z])', ' ', inp)
    if out[0] == ' ':
        out = out[1::]
    return out

def merge_patient_df(df_tissue_type,df_surv,s_hue,s_patient='Old_Pt_ID'):
    df_tissue_type[s_patient] = df_tissue_type.index
    df_merge = df_tissue_type.merge(df_surv.loc[:,[s_patient,s_hue]],on=s_patient)
    df_merge.set_index(s_patient,inplace=True)
    df_merge = df_merge[~df_merge.index.duplicated()]
    return(df_merge)

def correlation_scatterplot(df_merge,s_hue,ls_scores,s_type_plot,alpha = 0.05,b_legend=False,b_in=False):
    ls_colors = mpl.cm.tab10.colors
    ls_color = [ls_colors[7],ls_colors[0],ls_colors[1],ls_colors[2],ls_colors[3],ls_colors[4],ls_colors[5],ls_colors[6],ls_colors[8],ls_colors[9]] #
    d_fig = {}
    df_merge[s_hue] = df_merge.loc[:,s_hue].fillna('__')
    d_color = dict(zip(sorted(df_merge.loc[:,s_hue].unique()),ls_color))
    for s_marker in ls_scores:
        for s_tum in df_merge.drop(s_hue,axis=1).columns:
            ls_index = sorted(set(df_merge.loc[:,s_marker].dropna().index).intersection(set(df_merge.loc[:,s_tum].dropna().index)))
            if not s_tum == s_marker:
                try:
                    r, pvalue = stats.pearsonr(y=df_merge.loc[ls_index,s_tum], x=df_merge.loc[ls_index,s_marker])
                except:
                    print(f'skipping {s_tum}')
                    pvalue=1
                if pvalue < alpha:
                    fig, ax = plt.subplots(figsize=(2.5,2.6),dpi=300)
                    sns.scatterplot(x=s_tum, y=s_marker, data=df_merge.loc[ls_index],
                                    ax=ax,s=15,legend=b_legend,hue=s_hue,palette=d_color)  
                    ax.set_ylabel(f'{s_marker}') 
                    if b_in:
                        ax.set_xlabel(f'{s_type_plot} {s_tum.replace("_"," in ")}') 
                    else:
                        ax.set_xlabel(f'{s_tum}') 
                    ax.set_title(f'{s_type_plot} {s_tum.split("_")[0]}\n vs {s_marker}\n r = {r:.2f} \np = {pvalue:.3f} n={len(df_merge.loc[ls_index])}',fontsize=12) 
                    if b_legend:
                        ax.legend(bbox_to_anchor=(1.01,.3))
                    plt.tight_layout()
                    d_fig.update({f'scatterplot_{s_type_plot}_{s_tum}_vs_{s_marker}':fig})
    return(d_fig)

def multi_plots(df_merge,df_surv,s_find='all',n_neighbors=5,resolution=0.5,ls_hue=['Cohort'],s_title='Abundance',figsize=(14,10),figsize2=(12,4)):
    ls_col = df_merge.columns[(df_merge.columns.str.contains(s_find))] #
    adata = make_adata(df_merge,ls_col,n_neighbors=6)
    #umap
    for s_hue in ls_hue:
        adata.obs[s_hue] = adata.obs.index.map(dict(zip(df_surv.Old_Pt_ID,df_surv.loc[:,s_hue])))
        sc.pl.umap(adata, color=s_hue,title=s_title)
    # clusters
    sc.tl.leiden(adata,resolution=resolution)
    sc.pl.umap(adata, color='leiden',title=s_title,palette=mpl.cm.Set1.colors)
    # patient heatmap
    ls_annot = adata.obs.columns.tolist()
    df_p = df_merge.loc[:,ls_col.tolist()].fillna(0).merge(adata.obs.astype('object'),left_index=True,right_index=True,suffixes=("","_y")) #[ls_hue[0]]
    gg, df_annot = patient_heatmap(df_p,ls_col,ls_annot,figsize=figsize,linkage='average')
    #survival
    df_st = df_p.merge(df_surv.loc[:,['Survival','Survival_time','Old_Pt_ID']].set_index('Old_Pt_ID'),left_index=True,right_index=True).dropna()
    df_st = df_st[~df_st.index.duplicated()]
    fig1, p1, fig2, p2 = km_cph(df_st)
    print(p1)
    g = sns.clustermap(df_p.drop(s_hue,axis=1).groupby('leiden').mean(),xticklabels=1,figsize=figsize2,cmap='viridis',
                  dendrogram_ratio=0.1, cbar_pos=(.04, 0.92, 0.03, 0.10))
    g.ax_heatmap.set_xticklabels([item.get_text().replace(s_find,'').replace('_',' ').replace('in','') for item in g.ax_heatmap.get_xticklabels()])
    return(adata)


# networkx
def make_links(rho,pval,s_split='_in',s_split2='  '):
    pvals = pval.stack().reset_index()
    pvals.columns = ['var1', 'var2', 'pvalue']
    links = rho.stack().reset_index()
    links.columns = ['var1', 'var2', 'value']
    links['var1c'] = [f'{item.split(s_split)[0]}' for item in links.var1]
    links['var2c'] = [f'{item.split(s_split)[0]}' for item in links.var2]
    links['var1c'] = [f'{item.split(s_split2)[-1]}' for item in links.var1c]
    links['var2c'] = [f'{item.split(s_split2)[-1]}' for item in links.var2c]
    links['pval'] = pvals.pvalue
    return(links)

def drop_list(links,ls_drop):
    es_drop_out = set()
    for s_drop in ls_drop:
        es_drop_out.update(set(links[links.var1c.str.contains(s_drop)].var1c.unique()))
        es_drop_out.update(set(links[links.var2c.str.contains(s_drop)].var2c.unique()))
    ls_drop_out = sorted(es_drop_out)
    return(ls_drop_out)

def filter_links(links, ls_drop, thresh=0.5,b_greater=True):
    if b_greater:
        links_filtered = links.loc[(links['pval'] < 0.05) & (links['value'] > thresh) 
                                   & (links['var1'] != links['var2']) & (links['var1c'] != links['var2c'])
                                   & (~links['var2c'].isin(ls_drop)) & (~links['var1c'].isin(ls_drop))]
    else:
        links_filtered = links.loc[(links['pval'] < 0.05) & (links['value'] < thresh) 
                                   & (links['var1'] != links['var2']) & (links['var1c'] != links['var2c'])
                                   & (~links['var2c'].isin(ls_drop)) & (~links['var1c'].isin(ls_drop))]
    links_filtered_copy = links_filtered.copy()
    return(links_filtered_copy)

def make_graph_layout(links_filtered,source,target):
    G=nx.from_pandas_edgelist(df=links_filtered, source=source, target=target,edge_attr=['weight','color','style'])
    df = pd.DataFrame(index=G.nodes(), columns=G.nodes())
    for row, data in nx.shortest_path_length(G):
        for col, dist in data.items():
            df.loc[row,col] = dist
    df = df.fillna(df.max().max())
    layout = nx.kamada_kawai_layout(G, dist=df.to_dict())
    return(G, layout)

def concat_links(links_filtered,links_filtered_neg):
    links_filtered_neg['color'] = 'blue'
    links_filtered['color'] = 'red'
    links_filtered_neg['style'] = 'dashed'
    links_filtered['style'] = 'solid'
    links_filtered2 = pd.concat([links_filtered_neg,links_filtered]) 
    links_filtered2['weight'] = abs(links_filtered2.value)
    links_filtered2['names'] = links_filtered2.var1c+links_filtered2.var2c
    return(links_filtered2)

def get_G_attributes(G):
    ls_color = [n3['color'] for n1, n2, n3 in list(G.edges(data=True))]
    ls_style = [n3['style'] for n1, n2, n3 in list(G.edges(data=True))]
    ls_weight = [n3['weight']*2 for n1, n2, n3 in list(G.edges(data=True))]
    return(ls_color,ls_style,ls_weight)

def label_hubs(G, hubs):
    labels = {}    
    for node in G.nodes():
        if node in hubs:
            #set the node name as the key and the label as its value 
            labels[node] = node
    return(labels)

def plt_sig3(df_test,ls_order,ax):
    props = {'connectionstyle':matplotlib.patches.ConnectionStyle.Bar(armA=0.0, armB=0.0, fraction=0.0, angle=None),
             'arrowstyle':'-','linewidth':1}
    #draw on axes
    y_lim = ax.get_ylim()[1]
    y_lim_min = ax.get_ylim()[0]
    y_diff = (y_lim-y_lim_min)/10
    for count, s_index in enumerate(df_test[df_test.reject].index):
        y_test = (y_diff+count*y_diff)
        text =f"p-adj={df_test.loc[s_index,'p-adj']:.1}"
        one = df_test.loc[s_index,'group1']
        two = df_test.loc[s_index,'group2']
        x_one = np.argwhere(np.array(ls_order) == one)[0][0]
        x_two = np.argwhere(np.array(ls_order) == two)[0][0]
        ax.annotate(text, xy=(np.mean([x_one,x_two]),y_lim - y_test),fontsize=8)
        ax.annotate('', xy=(x_one,y_lim - y_test), xytext=(x_two,y_lim - y_test), arrowprops=props)
        #break
    return(ax)

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])          

def categorical_correlation_boxplot(df_all,s_group,s_marker,s_type,s_cell,alpha=0.05,s_propo='in',b_ttest=False,figsize=(3,3)):
        df_group = df_all.loc[:,[s_group,s_marker]]
        df_group = df_group.dropna()
        df_group = df_group.loc[~df_group.index.duplicated(),~df_group.columns.duplicated()]
        #print(len(df_group))
        if s_group==s_marker:
            return(None,1,None,None)
        if df_group.loc[:,s_group].nunique() > 2:
            statistic, pvalue = group_median_diff(df_group,s_group=s_group,s_marker=s_marker)
        elif df_group.loc[:,s_group].nunique() == 2:
            s_high = df_group.loc[:,s_group].unique()[0]
            s_low = df_group.loc[:,s_group].unique()[1]
            n_high = sum(df_group.loc[:,s_group]==s_high)
            n_low = sum(df_group.loc[:,s_group]==s_low)

            try:
                statistic,pvalue = stats.mannwhitneyu(df_group.loc[df_group.loc[:,s_group]==s_high,s_marker],
                                               df_group.loc[df_group.loc[:,s_group]==s_low,s_marker])
            except:
                print(f'mann whitney u error: {s_group} vs {s_marker}')
                pvalue = 1
            if b_ttest:
                statistic,pvalue = stats.ttest_ind(df_group.loc[df_group.loc[:,s_group]==s_high,s_marker],
                                               df_group.loc[df_group.loc[:,s_group]==s_low,s_marker])
        else:
            fig = None
            pvalue = 1
        # perform multiple pairwise comparison (Tukey HSD)
        if np.isnan(pvalue):
            pvalue = 1
        if pvalue <= alpha:
            m_comp = pairwise_tukeyhsd(endog=df_group.loc[:,s_marker].fillna(0), groups=df_group.loc[:,s_group], alpha=0.1)
            df_test, ls_order = df_from_mcomp(m_comp)
            if pd.Series(['Templates_binary','Rearrangements_binary','TE_Clones_binary','Recurrence']).isin([s_group]).any(): 
                ls_order = np.flip(ls_order) 
            # if df_group.loc[:,s_group].nunique() > 2:
            #     figsize=(6,3)
            # else:
            #     figsize=(3,3)
            fig, ax = plt.subplots(figsize=figsize,dpi=300)
            sns.boxplot(data=df_group,x=s_group,y=s_marker,showfliers=False,ax=ax,order=[(item) for item in ls_order],palette=[lighten_color(color,0.7) for color in sns.color_palette()])
            sns.stripplot(data=df_group,x=s_group,y=s_marker,ax=ax,order=[(item) for item in ls_order],alpha=1,linewidth=1)#,palette='dark'
            if df_group.loc[:,s_group].nunique() > 2:
                plt_sig3(df_test,ls_order,ax)
                ax.set_title(f'{s_group} versus\n {s_marker} {s_propo} {s_cell}\n p={pvalue:.2f}')
            else:
                ax.set_title(f'{s_group} versus\n {s_marker} {s_propo} {s_cell}\n p={pvalue:.2f} (n={n_low}, {n_high})')
            #ax.set_ylim(ax.get_ylim()[0],ax.get_ylim()[1])
            ax.set_ylabel(f'{s_marker} {s_type}')
            plt.tight_layout()
        else:
            fig = None
            df_test = None
        return(fig, pvalue,df_test,df_group) 
    
def add_old_pt_ID(df_all,df_surv,ls_group,ls_Annot,s_primary_met='PDAC'):
    ls_need = ['Tissue','Old_Pt_ID']
    df_all['Old_Pt_ID'] = df_all.index
    df_all.index.name = None
    df_all = df_all.merge(df_surv.loc[:,(ls_group + ls_Annot + ls_need)],on='Old_Pt_ID',suffixes=('','_also'))
    df_all = df_all[~df_all.Old_Pt_ID.duplicated()]

    if not s_primary_met=='Mets':
        df_all = df_all[df_all.Tissue=='PDAC']

    if s_primary_met=='Mets':
        s_hue = 'Cohort'
        print('updating Mets Cohort') #necessary?
        d_mets = {'Ovary_Met':np.nan, 'Peritoneum_Met':np.nan, 'Lung_Met':'lung_cohort', 'Liver_Met':'liver_cohort',  'Lymph node':np.nan}
        df_all[s_hue] = df_all.Old_Pt_ID.map(dict(zip(df_all.Old_Pt_ID,df_all.Tissue))).map(d_mets)
    return(df_all)

def plot_youden_tcr(df_patient,s_tcr,pos_label,figsize=(5,6)):
    d_result = {}
    ls_col_thresh = df_patient.columns[(df_patient.columns.str.contains('_day_survival')) & (df_patient.columns.str.contains(s_tcr))]
    if len(ls_col_thresh)<3:
        fig,ax = plt.subplots(2,1,dpi=200,figsize=figsize,sharey=True)
    else:
        fig,ax = plt.subplots(3,2,dpi=200,figsize=figsize,sharex=True,sharey=True)
    ax=ax.ravel()
    for ax_idx, s_col in enumerate(ls_col_thresh):
        #print(s_col)
        #s_tcr = '_'.join(s_col.split('_day_survival')[0].split('_')[0:-1])
        df_patient.loc[:,[s_col,s_tcr]]
        df_test = df_patient.loc[:,[s_col,s_tcr]].dropna().copy()
        fpr, tpr, thresholds = roc_curve(y_true=df_test.loc[:,s_col].values,
                                         y_score=df_test.loc[:,s_tcr].values,
                                         pos_label=pos_label)
        idx = np.argmax(tpr - fpr)
        youden = np.max(tpr - fpr)
        cutoff = thresholds[idx]
        d_result.update({s_col:cutoff})
        print(f'{s_col} {youden:.2}')
        display = RocCurveDisplay(fpr=fpr, tpr=tpr)
        ax[ax_idx].plot([0, 1], [0, 1], transform=ax[ax_idx].transAxes,ls='--',color='lightgray')
        display.plot(ax=ax[ax_idx])
        ax[ax_idx].text(fpr[idx]+.01,tpr[idx]-.075,f'Youden={youden:.2}\ncutoff={cutoff:.2}',fontsize='small')
        ax[ax_idx].scatter(fpr[idx],tpr[idx])
        ax[ax_idx].set_title(s_col.replace(f'{s_tcr}_',''),fontsize='medium')
        ax[ax_idx].set_xlabel('')
        ax[ax_idx].set_ylabel('')#ax[ax_idx].get_ylabel(),fontsize='small'#
    if len(ls_col_thresh) % 2:
        ax[-1].axis('off')
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(s_tcr,pad=20)
    plt.tight_layout()
    return(fig,d_result)
    
#GSEA

#plotting functions

def compare_dataframe(ls_plot_items,d_en,sorter_combined,ls_columns):
    df_all = pd.DataFrame()
    es_marker = set()
    for s_comp in ls_plot_items: 
        df_plot_long = pd.DataFrame()
        for s_direction in ['UP','DN']:
            s_compare =f'{s_comp}_{s_direction}'
            #print(s_compare)
            #manual sorting
            df_plot = d_en[s_compare].loc[d_en[s_compare].NAME.isin(sorter_combined),ls_columns]
            df_plot['direction'] = s_direction
            df_plot_long = pd.concat([df_plot_long,df_plot])
            es_marker = es_marker.union(set(df_plot.NAME))#.union(es_add)
            #sorter_all.update({s_compare:sorter})
        df_plot_long = df_plot_long.sort_values(by='NES')
        df_plot_long['comparison'] = s_comp
        df_all = pd.concat([df_all,df_plot_long])
    df_all.NAME = df_all.NAME.astype('category')
    df_all.NAME = df_all.NAME.cat.set_categories(sorter_combined) #order manually
    df_plot_bar = df_all[df_all.NAME.isin(sorter_combined)].sort_values(by=['comparison','NAME'])
    return(df_plot_bar,es_marker)

def plot_double_bars(df_plot_bar,ls_plot_items,d_colorblind,sorter_combined,hatchbar='',figsize=(5,4)):
    # plot figure
    fig, ax = plt.subplots(dpi=200,figsize=figsize)
    height=0.4
    for idx, s_comp in enumerate(ls_plot_items):
        df_comp = df_plot_bar[df_plot_bar.comparison==s_comp]
        df_comp.set_index('NAME',inplace=True)
        if idx == 0:
            indices = np.arange(len(df_comp.index))
            ax.barh(y=indices+height/2, width=df_comp.loc[sorter_combined,'NES'],height=height,
                    edgecolor='k', linewidth=1,
                color=[d_colorblind[item] for item in df_comp.color],label=s_comp,alpha=0.8,
               )
        else: 
            ax.barh(y=indices-height/2, width=df_comp.loc[sorter_combined,'NES'],height=height,
                    label=s_comp,alpha=0.8,edgecolor='k',hatch=hatchbar,
                color=[d_colorblind[item] for item in df_comp.color],
               )
        ax.set_yticks(range(len(df_comp.index)))
        ax.set_yticklabels(df_comp.index)
    ax.set_title('',loc='left')
    ax.set_xlabel('NES')
    return(fig, ax)

def twin_pvalue_axis(df_plot_bar,ls_plot_items,sorter_combined,fig,ax,height=0.4,hatch=''):
    #pvalue axis
    ax2 = ax.twiny()
    for idx, s_comp in enumerate(ls_plot_items):
        print(s_comp)
        df_comp = df_plot_bar[df_plot_bar.comparison==s_comp]
        df_comp.set_index('NAME',inplace=True)
        if idx == 0:
            indices = np.arange(len(df_comp.index))
            ax2.scatter(y=indices+height/2, x=df_comp.loc[sorter_combined,'FDR.q.val'],
                        facecolor='gray',linewidth=0.5,edgecolor='k',alpha=0.8)
        else: 
            ax2.scatter(y=indices-height/2, x=df_comp.loc[sorter_combined,'FDR.q.val'],
                       facecolor='gray',linewidth=0.5,edgecolor='k',alpha=0.8,hatch=3*hatch,)#hatch=3*'//',
        #break
    ax2.set_xlabel('FDR.Q', color='gray') 
    ax2.tick_params(axis='x', labelcolor='gray')
    return(fig,ax2)

from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
class HandlerEllipse(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Ellipse(xy=center, width=height + xdescent,
                             height=height + ydescent)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]
    
def plot_double_bars_heat(ls_plot_items,df_plot_bar,d_labels,sorter_combined,cmap,norm,hatch='//',height=0.4,figsize=(4,4),anchor=(1.25,1),x_label='NES',ncol=2):#mappable
    # Choose colormap
    N = len(norm.boundaries)
    cmap2 = cm.get_cmap('RdYlBu_r', N-1)
    colors = [mpl.colors.to_hex(cmap2(i)) for i in range(N-1)]
    fig, ax = plt.subplots(dpi=300,figsize=figsize)
    for idx, s_comp in enumerate(ls_plot_items):
        df_comp = df_plot_bar[df_plot_bar.comparison==s_comp]
        df_comp.set_index('NAME',inplace=True)
        if idx == 0:
            indices = np.arange(len(df_comp.index))
            ax.barh(y=indices+height/2, width=df_comp.loc[sorter_combined,'NES'],height=height,
               color=[colors[item[0]] for item in df_comp.FDR_color],
                    label=s_comp,linewidth=1,
                    edgecolor='k',
               )
        else: 
            ax.barh(y=indices-height/2, width=df_comp.loc[sorter_combined,'NES'],height=height,
                    hatch=hatch,edgecolor='k',
                color=[colors[item[0]] for item in df_comp.FDR_color],label=s_comp,linewidth=1,
               )
        ax.set_yticks(range(len(df_comp.index)))
        ax.set_yticklabels(df_comp.index)
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),#mappable=mappable,
                 ax=ax,label='FDR.q.val')
    if len(ls_plot_items) == 2:
        handles = [mpl.patches.Patch(facecolor='white', edgecolor='black',
                             label=d_labels[ls_plot_items[0]]),
                      mpl.patches.Patch(facecolor='white', edgecolor='black',
                             label=d_labels[ls_plot_items[1]],
                                        hatch=hatch)]
    else:
        handles = [mpl.patches.Patch(facecolor='white', edgecolor='black',
                             label=d_labels[ls_plot_items[0]])]
    ax.legend(handles=handles,bbox_to_anchor=anchor,markerscale=1,title='Comparison',ncol=ncol,loc='upper center')
    ax.set_xlabel(x_label)
    return(fig, ax)

#add FDRQ
def add_fdrq_legend(texts,hatch,anchor=(1.01,0.6)):
    if hatch == '//':
        c1 = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='gray', edgecolor="k")
        c2 = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='gray', edgecolor="k",hatch=3*hatch)
        
        plt.legend([c1,c2],texts,loc='upper left', bbox_to_anchor=anchor, ncol=1,title='FDR.Q',
                   handler_map={mpatches.Circle: HandlerEllipse()}).get_frame()
    else:
        c1 = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='gray', edgecolor="k")
        plt.legend([c1],['FDR.Q'],loc='upper left', bbox_to_anchor=anchor, ncol=1,#title='FDR.Q',
                   handler_map={mpatches.Circle: HandlerEllipse()}).get_frame()
        
def plot_double_bars_grid(df_plot_bar,ls_plot_items,d_colorblind,sorter_combined,hatchbar='',figsize=(5,4),xlabel='NES'):
    # plot figure
    fig = plt.figure(figsize=figsize,dpi=200)#, layout="constrained"
    spec = fig.add_gridspec(1,2,wspace=0,)
    ax1 = fig.add_subplot(spec[:, 0])
    height=0.4
    for idx, s_comp in enumerate(ls_plot_items):
        df_comp = df_plot_bar[df_plot_bar.comparison==s_comp]
        df_comp.set_index('NAME',inplace=True)
        if idx == 0:
            indices = np.arange(len(df_comp.index))
            ax1.barh(y=indices+height/2, width=df_comp.loc[sorter_combined,'NES'],height=height,
                    edgecolor='k', linewidth=1,
                color=[d_colorblind[item] for item in df_comp.color],label=s_comp,alpha=0.8,
               )
        else: 
            ax1.barh(y=indices-height/2, width=df_comp.loc[sorter_combined,'NES'],height=height,
                    label=s_comp,alpha=0.8,edgecolor='k',hatch=hatchbar,
                color=[d_colorblind[item] for item in df_comp.color],
               )
        ax1.set_yticks(range(len(df_comp.index)))
        ax1.set_yticklabels(df_comp.index)
    ax1.set_title('',loc='left')
    ax1.set_xlabel(xlabel)
    return(fig, ax1)

def add_pvalue_axis(df_plot_bar,ls_plot_items,sorter_combined,fig,ax,height=0.4,hatch=''):
    spec = fig.add_gridspec(1,2)
    ax2 = fig.add_subplot(spec[:, 1])
    for idx, s_comp in enumerate(ls_plot_items):
        print(s_comp)
        df_comp = df_plot_bar[df_plot_bar.comparison==s_comp]
        df_comp.set_index('NAME',inplace=True)
        if idx == 0:
            indices = np.arange(len(df_comp.index))
            ax2.scatter(y=indices+height/2, x=df_comp.loc[sorter_combined,'FDR.q.val'],
                        facecolor='gray',linewidth=0.5,edgecolor='k',alpha=0.8)
        else: 
            ax2.scatter(y=indices-height/2, x=df_comp.loc[sorter_combined,'FDR.q.val'],
                       facecolor='gray',linewidth=0.5,edgecolor='k',alpha=0.8,hatch=3*hatch,)#hatch=3*'//',
        #break
    ax2.set_xlabel('FDR.Q', color='gray') 
    ax2.tick_params(axis='x', labelcolor='gray')
    ax2.set_yticks([])
    return(fig,ax2)
def add_fdrq_legend_ax2(texts,hatch,ax2,anchor=(1.01,0.6)):
    if hatch == '//':
        c1 = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='gray', edgecolor="k")
        c2 = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='gray', edgecolor="k",hatch=3*hatch)
    
        legend2 = ax2.legend([c1,c2],texts,loc='upper left', bbox_to_anchor=anchor, ncol=1,title='FDR.Q',frameon=True,
                   handler_map={mpatches.Circle: HandlerEllipse()}).get_frame()
    else:
        c1 = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor='gray', edgecolor="k")
        legend2 = ax2.legend([c1],['FDR.Q'],loc='upper left', bbox_to_anchor=anchor, ncol=1,frameon=True,
                   handler_map={mpatches.Circle: HandlerEllipse()}).get_frame()
    return(legend2)

### youden ###

#df_counts = pd.DataFrame(columns=['Long','Short','Short_Censored'])
#     df_counts['total'] = df_counts.sum(axis=1)
#     df_counts['included_ROC'] = df_counts.loc[:,['Short','Long']].sum(axis=1)
    #break
#df_counts.to_csv('patient_counts_TCR_ROC.csv')
            #Num_censored = (df_patient2.loc[:,s_time] <= SurvivalThreshold) & (df_patient2.loc[:,s_tcr].notna()) & (df_patient2.loc[:,s_censor] == 0) & (df_patient2.Alive_30_days_post_surgery)
            #df_counts.loc[f'{s_tcr}_{SurvivalThreshold}'] = [GoodIdx.sum(),BadIdx.sum(),Num_censored.sum()]
        
def youden_high_good(df_patient,b_primary,s_time,s_censor,s_tcr,figsize=(4,4)):
    ls_foci = pd.Series(['Shannon_Entropy_Tumor','Templates_per_ng', 
               'Productive_Rearrangements','Simpsons_Diversity_Tumor',#'Fraction Shared Clones 10',
               'Clonality_Tumor', 'Clonality_Blood','Productive_Rearrangements_Blood',
            'Fraction Tumor Distinct TCRs','Percent Tumor Distinct Clones',
            'Hill_1', 'Hill_2', 'Hill_3','Fraction Shared Clones 10',
            #'Fraction Shared Clones 5','Fraction Shared Clones 2',   
       'Hill_4', 'Hill_5', 'Hill_6',
    'propshared_proportion_shared_clones_blood','propshared_proportion_shared_templates_blood',
           'Hill_1_blood','d50_Clones', 'd50_Clones_blood',
           'chao_Estimator','true_Value','Fraction Shared Clones 0to2',
      'Fraction Shared Clones 2to5','Fraction Shared Clones 5to10','Fraction Shared Clones >=10',
     'Fraction Tumor Distinct Clones 0to2','jaccard_Tumor','public_Tumor','jaccard_Tumor_blood','public_Tumor_blood',#'morisita_Tumor',
     'Fraction Tumor Distinct Clones 2to5','Fraction Tumor Distinct Clones 5to10',
     'Fraction Tumor Distinct Clones >=10','jaccard_Met','jaccard_Primary',
              'jaccard_Blood','public_Blood','public_Met','public_Primary',])
    pal_porg_r = ('#E69F00','#56B4E9')
    sns.set_palette(pal_porg_r)
    pos_label='long'
    dd_result = {}
    d_fig={}
    #for s_tcr in ls_foci:
        #print(s_tcr)
    if ls_foci.isin([s_tcr]).any():
        print(f'{pos_label} = good')
        fig, d_result = plot_youden_tcr(df_patient,s_tcr,pos_label)
        fig.savefig(f'figures/youden_{s_tcr}_{pos_label}_{b_primary}.png')
        plt.close(fig)    
        dd_result.update(d_result)
    else:
        return(d_fig)
    for s_surv in ['545_day_survival']:#'180_day_survival',
        keys = df_patient.columns[(df_patient.columns.str.contains(s_surv))]
        for key in keys:
            s_tcr = '_'.join(key.split('_day_survival')[0].split('_')[0:-1])
            try:
                thresh = dd_result[key]
            except:
                continue
            df_km = df_patient.loc[df_patient.Alive_30_days_post_surgery,[s_tcr,s_time,s_censor]].dropna()
            df_km[key] = (df_km.loc[:,s_tcr] > thresh).replace({True:'High',False:'Low'})
            fig,ax,ls_order = km_plot(df_km,key,s_time,s_censor,fontsize='small',figsize=figsize)
            d_fig.update({f'{b_primary}_{key}':fig})
    return(d_fig)
#low = good
def youden_low_good(df_patient,b_primary,s_time,s_censor,s_tcr,b_include=False,figsize=(4,4)):
    ls_foci = pd.Series(['Shannon_Entropy_Blood','Simpsons_Diversity_Blood',
                         'Fraction Shared Clones 10',
           # 'Fraction Shared Clones 5','Fraction Shared Clones 2',   
     "Simpson's Evenness tumor", "Simpson's Evenness blood",
               'Hill_2_blood', 'Hill_3_blood','chao_Estimator_blood','true_Value_blood',
        'Hill_4_blood', 'Hill_5_blood','Hill_6_blood',
              'Simpsons_Evenness_Blood','Simpsons_Evenness_Tumor',
                         'ginisimp_Value','ginisimp_Value_blood',
                  'propshared_proportion_shared_clones','morisita_Tumor_blood',
                   'propshared_proportion_shared_clones_blood',
       'propshared_proportion_shared_templates_blood',
      'propshared_proportion_shared_templates',])

    pal_porg_r = ('#E69F00','#56B4E9')
    sns.set_palette(pal_porg_r)
    pal_porg = ('#56B4E9','#E69F00')
    sns.set_palette(pal_porg)
    pos_label='short'
    dd_result = {}
    d_fig = {}
    #for s_tcr in ls_foci:
    if ls_foci.isin([s_tcr]).any():
            fig, d_result = plot_youden_tcr(df_patient,s_tcr,pos_label)
            print(f'{pos_label} = good')
            fig.savefig(f'figures/youden_{s_tcr}_{pos_label}_{b_primary}.png')
            plt.close(fig)
            dd_result.update(d_result)
    elif b_include:
        fig, d_result = plot_youden_tcr(df_patient,s_tcr,pos_label)
        plt.close(fig)
        dd_result.update(d_result)
    else:
        return(d_fig)
    sns.set_palette(pal_porg_r)
    for s_surv in ['545_day_survival']:#'180_day_survival',
        keys = df_patient.columns[(df_patient.columns.str.contains(s_surv))]
        for key in keys:
            s_tcr = '_'.join(key.split('_day_survival')[0].split('_')[0:-1])
            try:
                thresh = dd_result[key]
            except:
                continue
            df_km = df_patient.loc[df_patient.Alive_30_days_post_surgery,['Public_Patient_ID',s_tcr,s_time,s_censor]].dropna()
            df_km[key] = (df_km.loc[:,s_tcr] > thresh).replace({True:'High',False:'Low'})
            fig,ax,ls_order = km_plot(df_km,key,s_time,s_censor,fontsize='small',figsize=figsize)
            d_fig.update({f'{b_primary}_{key}':fig})
    return(d_fig)

def process_overlap(df,df_meta,ls_site=['Blood_Type','Site','Sample_Type']):
    df = df.copy()
    df_meta['Cohort_Blood'] = df_meta.Cohort + '_' + df_meta.Blood_Type
    df_meta['Cohort_Tumor'] = df_meta.Cohort + '_'+ df_meta.Site.replace('Blood',np.nan)
    #ls_site = [#'Cohort_Blood','Cohort_Tumor',
         #'pORG_0.2_Primary_quartiles_Site','pORG_0.2_Met_quartiles_Site' ]   
    for s_site in ls_site:
        #print(f'\n{s_site}:')
        for s_group in df_meta.loc[:,s_site].dropna().unique():
            df_test = df_meta[df_meta.loc[:,s_site]==s_group]
            if s_site.find('pORG') > -1:
                s_rep = s_site.split('_')[2]
                #print(s_rep)
                s_group = s_group.replace('high_Blood',f'high_{s_rep}_Bld').replace('low_Blood',f'low_{s_rep}_Bld')
            for s_test in ['Blood','Primary','Met']:
                if (len(df_test[df_test.Site==s_test])) > 0:
                    #print(len(df))
                    ls_patient = df_test.loc[df_test.Site==s_test,'Sample']
                    se_sum = df.loc[:,ls_patient].sum(axis=1)
                    #print(len(se_sum))
                    df.loc[se_sum.index,s_group] = se_sum
                    #break
            #break
        #break
    df_out = df.loc[:,~df.columns.str.contains('ST-')]
    #df_out = df_out.reset_index().rename({'index':'Sample'},axis=1)
    return(df_out)

def process_overlap2(df,df_meta,ls_site=['Blood_Type','Site','Sample_Type']):
    df = df.copy()
    df_meta['Cohort_Blood'] = df_meta.Cohort + '_' + df_meta.Blood_Type
    df_meta['Cohort_Tumor'] = df_meta.Cohort + '_'+ df_meta.Site.replace('Blood',np.nan)
    #ls_site = [#'Cohort_Blood','Cohort_Tumor',
         #'pORG_0.2_Primary_quartiles_Site','pORG_0.2_Met_quartiles_Site' ]   
    for s_site in ls_site:
        #print(f'\n{s_site}:')
        for s_group in df_meta.loc[:,s_site].dropna().unique():
            df_test = df_meta[df_meta.loc[:,s_site]==s_group]
            #print(len(df_test))
            if s_site.find('pORG') > -1:
                s_rep = s_site.split('_')[2]
                #print(s_rep)
                s_group = s_group.replace('high_Blood',f'high_{s_rep}_Bld').replace('low_Blood',f'low_{s_rep}_Bld')
            if len(df_test) > 0:
                ls_patient = df_test.loc[:,'Sample']
                se_sum = df.loc[:,ls_patient].sum(axis=1)
                #print(len(se_sum))
                df.loc[se_sum.index,s_group] = se_sum
                #break
            #break
        #break
    df_out = df.loc[:,~df.columns.str.contains('ST-')]
    #df_out = df_out.reset_index().rename({'index':'Sample'},axis=1)
    return(df_out)

# function to calculate Cohen's d for independent samples
def cohend(df_sp,x,y,ls_groups):
    d1 = df_sp.loc[df_sp.loc[:,x] == ls_groups[0],y]
    d2 = df_sp.loc[df_sp.loc[:,x] == ls_groups[1],y]
    # calculate the size of samples
    n1, n2 = len(d1), len(d2)
    # calculate the variance of the samples
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    # calculate the pooled standard deviation
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    # calculate the means of the samples
    u1, u2 = np.mean(d1), np.mean(d2)
    # calculate the effect size
    return (u1 - u2) / s

#https://github.com/sflury/kendall/blob/main/kendall.py
import numpy as np
from scipy.special import erf,ndtr,ndtri

'''
Name
    kendall

Purpose
    Compute Kendall's tau correlation coefficient following Akritas & Siebert
    (1996) MNRAS 278, 919-924, accounting for left-censored data (i.e., data
    with upper limits).

Arguments
    :x (*np.ndarray*): 1xN array of independent variable, containing either
                            measured values or detection limits
    :y (*np.ndarray*): 1xN array of dependent variable, containing either
                            measured values or detection limits

Keyword Arguments:
    :censors (*np.ndarray*): 2xN array containing censors of `x` and `y` with 1
                            representing an uncensored datum and 0 representing
                            a left-censored (upper-limit) datum. Default is
                            uncensored (2xN array of 1s).

Returns
    :tau (*float*): Kendall's tau correlation coefficient, accounting for
                            upper limits on data
    :p (*float*): probability that the null hypothesis (no correlation) is true
'''
def kendall(x,y,censors=None):
    # check that x and y are identical lengths
    if len(x) != len(y):
        print('X and Y different lengths!')
        return None
    else:
        # length of data
        n = len(x)
    # check if censors exist, and if not, assume no censoring
    if np.any(censors==None):
        censors = np.ones((2,n))
    # otherwise, check that censor has same shape as T
    elif len(censors) != 2 or len(censors[0]) != n:
        print('Censor length must match data length!')
        return
    # array of rank differences (J in A&S Eqns 1, 3, 5, and 6)
    J = np.zeros((2,n,n))
    # array of data (T in A&S Eqn 1)
    # switching from left- to right-censored
    T = np.array([-x,-y])
    tau = 0.
    # compute Kendall's tau
    for j in range(n):
        for i in range(j):
            for k in range(2):
                I_kij = 0.
                if T[k,i] < T[k,j]: I_kij = 1.
                I_kji = 0.
                if T[k,j] < T[k,i]: I_kji = 1.
                J[k,i,j] = censors[k,i]*I_kij - censors[k,j]*I_kji
    # A&S Eqn 3, 7
    tau = np.sum(np.prod(J,axis=0))
    # normalize by binomial coefficient (n choose 2)
    tau *= 2/(float(n)*(float(n)-1))
    # simple p-value
    var = 2*(2*n+5)/(9*n*(n-1))
    p = (1-erf(abs(tau)/np.sqrt(2*var)))
    return tau,p

'''
Name
    tau_conf

Purpose
    Determine confidence intervals on Kendall's tau with left-censored data by
    bootstrapping values of x and y. Kendall's tau calculated following
    Akritas & Siebert (1996) MNRAS 278, 919-924.

Arguments
    :x (*np.ndarray*): 1xN array of independent variable, containing either
                            measured values or detection limits
    :y (*np.ndarray*): 1xN array of dependent variable, containing either
                            measured values or detection limits

Keyword Arguments:
    :censors (*np.ndarray*): 2xN array containing censors of `x` and `y` with 1
                            representing an uncensored datum and 0 representing
                            a left-censored (upper-limit) datum. Default is
                            uncensored (2xN array of 1s).
    :p_conf (*float*): two-sided  probability interval of the desired confidence
                            interval. Default is 0.6826 (1-sigma).
    :n_samp (*int*): number of samples to draw. Default is 10^4.
    :method (*str*): 'montecarlo' for Monte Carlo sampling of uncertainties or
                            'bootstrap' for bootstrapping of measurements.
                            Default is 'montecarlo'.

'''
def tau_conf(x,y,x_err=None,y_err=None,censors=None,p_conf=0.6826,n_samp=int(1e4),method='montecarlo'):
    # check if censors exist, and if not, assume no censoring
    if np.any(censors==None):
        censors = np.ones((2,len(x)))
    # bootstrap sampling
    if method == 'bootstrap':
        tau_boot = []
        for i in range(n_samp):
            inds = np.unique(np.random.randint(0,high=len(x)-1,size=len(x)))
            tau_i,p_i = kendall(x[inds],y[inds],censors[:,inds])
            tau_boot += [tau_i]
        quants = np.quantile(tau_boot,[0.5-0.5*p_conf,0.5,0.5+0.5*p_conf])
        tau_lower,tau_upper = np.diff(quants)
    # MC sampling
    else:
        tau_mc = []
        x_mc = np.random.normal(size=(n_samp,len(x)),loc=x,scale=x_err)
        y_mc = np.random.normal(size=(n_samp,len(y)),loc=y,scale=y_err)
        for i in range(n_samp):
            tau_i,p_i = kendall(x_mc[i,:],y_mc[i,:],censors)
            tau_mc += [tau_i]
        quants = np.quantile(tau_mc,[0.5-0.5*p_conf,0.5,0.5+0.5*p_conf])
        tau_lower,tau_upper = np.diff(quants)
    return tau_lower,tau_upper

def patient_kendall(df_patient,sx,sy,s_censors):
    x = df_patient.loc[:,sx]
    y = df_patient.loc[:,sy]
    # gen random "upper limits"
    censors = np.ones((2,len(x)))
    uplims = df_patient.loc[:,s_censors]==0#np.unique(np.random.randint(0,high=len(x)-1,size=60))
    censors[1,uplims]=0
    #uncens = [i for i in range(len(x)) if i not in uplims]
    # faux uncertainties, accounting for upper limits
    x_err = 0.1*np.random.rand(len(x))*censors[0]
    y_err = 0.1*np.random.rand(len(y))*censors[1]
    # Kendall tau
    tau,p = kendall(x,y,censors=censors)
    #tau_lo,tau_up = util.tau_conf(x,y,censors,method='bootstrap')
    print(f'censored Kendall tau: {tau:.3f}, p: {p:.3}')#-{tau_lo:.3f}+{tau_up:.3f}
    return(tau,p)

def patient_pearson(df_patient,sx,sy):
    df = df_patient.loc[:,[sx,sy]].dropna()
    x = df.loc[:,sx]
    y = df.loc[:,sy]
    # gen random "upper limits"
    statistic,p=stats.pearsonr(x, y)
    print(f'pearson r: {statistic:.3f}, p: {p:.3}')#-{tau_lo:.3f}+{tau_up:.3f}
    return(statistic,p)
# # organotropism and pORG (side by side axis)
# sorter_combined =  ['HALLMARK_MYOGENESIS' ,'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
#   'HALLMARK_ANDROGEN_RESPONSE','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
#   'HALLMARK_DNA_REPAIR','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_E2F_TARGETS',
#   'HALLMARK_MYC_TARGETS_V1','HALLMARK_G2M_CHECKPOINT','HALLMARK_MITOTIC_SPINDLE',
#     'HALLMARK_GLYCOLYSIS','HALLMARK_MTORC1_SIGNALING','HALLMARK_PROTEIN_SECRETION',
#                   ]
# ls_plot_items = ['Top4th_pORG.20_vs_Bottom4th_pORG.20',
#     'LiverCohort_vs_LungNotLiverCohort']
# #generate dataframe with comparisons
# df_plot_bar,es_marker = util.compare_dataframe(ls_plot_items,d_en,sorter_combined,ls_columns)
# #add colors
# b_purist = df_plot_bar.comparison=='LiverCohort_vs_LungNotLiverCohort'
# b_psub = df_plot_bar.comparison=='Top4th_pORG.20_vs_Bottom4th_pORG.20'
# df_plot_bar.loc[b_purist & (df_plot_bar.direction=='UP'),'color'] = 'Liver'
# df_plot_bar.loc[b_purist & (df_plot_bar.direction=='DN'),'color'] = 'Lung'
# df_plot_bar.loc[b_psub & (df_plot_bar.direction=='UP'),'color'] = 'high pORG'
# df_plot_bar.loc[b_psub & (df_plot_bar.direction=='DN'),'color'] = 'low pORG'
# # plot figure
# fig, ax1 = plot_double_bars_grid(df_plot_bar,ls_plot_items,d_colorblind,sorter_combined,hatchbar,
#                                 figsize=(5.5,4.5)) #util.plot_double_bars
# #custom legend
# handles = []
# for s_color in ['high pORG', 'low pORG']:
#     facecolor = d_colorblind[s_color]
#     handles.append(mpl.patches.Patch(facecolor=facecolor, edgecolor='black',label=s_color,alpha=0.8,))
# #hatch
# for s_color in ['Liver', 'Lung']:
#     facecolor = d_colorblind[s_color]
#     handles.append(mpl.patches.Patch(facecolor=facecolor, edgecolor='black',label=s_color,
#                                      alpha=0.8,hatch=hatchbar))

# ax1.set_yticklabels([item.replace('HALLMARK_','').replace('EPITHELIAL_MESENCHYMAL_TRANSITION','EMT') for item in sorter_combined])
# ax1.set_xlim(-2.2,2.2)
# fig, ax2 = add_pvalue_axis(df_plot_bar,ls_plot_items,sorter_combined,fig,ax1,
#                                  height=height,hatch=hatch) #util.twin_pvalue_axis
# legend1 = ax2.legend(handles=handles,bbox_to_anchor=(1.1,1.01),markerscale=1,title='NES',frameon=False)
# ax2.set_xlim((-0.1, 0.5))
# ax2.axvline(0.2,linestyle='--',color='gray')
# plt.tight_layout()
# legend2 = add_fdrq_legend_ax2(texts= ['pORG', 'Cohort'],hatch=hatch,ax2=ax2)
# ax2.add_artist(legend1)
# ax2.add_artist(legend2)

# #pSUB and purist  (shared axis)
# hatch = '//'#''#
# hatchbar = '//'#''
# height=0.4
# sorter_combined = ['HALLMARK_XENOBIOTIC_METABOLISM','HALLMARK_PEROXISOME', 'HALLMARK_FATTY_ACID_METABOLISM',
#   'HALLMARK_BILE_ACID_METABOLISM', 'HALLMARK_PANCREAS_BETA_CELLS','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
#   'HALLMARK_APICAL_JUNCTION', 'HALLMARK_HYPOXIA', 'HALLMARK_GLYCOLYSIS']
# ls_plot_items = ['Top4th_pSUB.1eNeg4_vs_Bottom4th_pSUB.1eNeg4', 'Top4th_PurIST.Score_vs_Bottom4th_PurIST.Score']
# #generate dataframe with comparisons
# df_plot_bar,es_marker = util.compare_dataframe(ls_plot_items,d_en,sorter_combined,ls_columns)
# #add colors
# b_purist = df_plot_bar.comparison=='Top4th_PurIST.Score_vs_Bottom4th_PurIST.Score'
# b_psub = df_plot_bar.comparison=='Top4th_pSUB.1eNeg4_vs_Bottom4th_pSUB.1eNeg4'
# df_plot_bar.loc[b_purist & (df_plot_bar.direction=='UP'),'color'] = 'high PurIST'
# df_plot_bar.loc[b_purist & (df_plot_bar.direction=='DN'),'color'] = 'low PurIST'
# df_plot_bar.loc[b_psub & (df_plot_bar.direction=='UP'),'color'] = 'high pSUB'
# df_plot_bar.loc[b_psub & (df_plot_bar.direction=='DN'),'color'] = 'low pSUB'
# # plot figure
# fig, ax = util.plot_double_bars(df_plot_bar,ls_plot_items,d_colorblind,sorter_combined,hatchbar)
# #custom legend
# handles = []
# for s_color in ['high pSUB', 'low pSUB']:
#     facecolor = d_colorblind[s_color]
#     handles.append(mpl.patches.Patch(facecolor=facecolor, edgecolor='black',label=s_color,alpha=0.8,))
# #hatch
# for s_color in [ 'high PurIST', 'low PurIST']:
#     facecolor = d_colorblind[s_color]
#     handles.append(mpl.patches.Patch(facecolor=facecolor, edgecolor='black',label=s_color,
#                                      alpha=0.8,hatch=hatchbar))
# ax.legend(handles=handles,bbox_to_anchor = (1.01,1),markerscale=1,title='NES')
# ax.set_yticklabels([item.replace('HALLMARK_','').replace('EPITHELIAL_MESENCHYMAL_TRANSITION','EMT') for item in sorter_combined])
# #add pvalue axis
# fig, ax2 = util.twin_pvalue_axis(df_plot_bar,ls_plot_items,sorter_combined,fig,ax,
#                                  height=height,hatch=hatch)
# ax2.set_xlim((-0.0092, 0.21))
# plt.tight_layout()
# util.add_fdrq_legend(texts= ['pSUB', 'PurIST'],hatch=hatch)


#erase low FDR.q
#df_plot_bar.loc[(abs(df_plot_bar.loc[:,'NES']) < 1.5) | (df_plot_bar.loc[:,'FDR.q.val'] > 0.15),'NES'] = 0

# # check the quartiles code
# #s_porg ='pORG_0.2_Primary'#'pORG_78_Primary'#'pORG_0.2_Met'
# x = df_pri.loc[:,s_porg].dropna()
# b_cut = df_pri.loc[:,[s_porg,'Public_Patient_ID']].dropna().Public_Patient_ID
# b_cut = df_pri.loc[:,s_porg].dropna().index
# print(len(x))
# d_cut = {'quartiles':(4,['low','med-low','med-high','high']),
#          'tertiles' : (3,['low','med','high']),
#         'medians' : (2,['low','high'])}
# for s_col, tu_cut in d_cut.items():
#     i_cut = tu_cut[0]
#     labels = tu_cut[1]
#     q = pd.qcut(x, q=i_cut,labels=labels) 
#     break
# df_pri['test'] = q
# fig,ax=plt.subplots()
# sns.stripplot(data=df_pri,x='test',y=s_porg,ax=ax)
# ax.axhline(df_pri.loc[df_pri.test=='low',s_porg].max())
# ax.axhline(df_pri.loc[df_pri.test=='high',s_porg].min())