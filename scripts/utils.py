import numpy as np
import pandas as pd
import sys

def neut_format_parser(x, upper=None, lower=0.0005):
    try:
        _ = float(x)
        if _ < lower:
            return lower
        if upper is not None and _ > upper:
            return upper
        return _
    except:
        if x[0] == '>' or x[-1] == "*":
            return np.inf if not upper else upper
        if x[0] == '<':
            return lower
        if x == '--':
            return pd.NA
        sys.stderr.write(f'{x}\n')
        raise ValueError


# DMS calculation
import logomaker
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
rcParams['pdf.fonttype'] = 42

def plot_res_logo(res, prefix, shownames={}, rownames=None, site_thres=0.0, force_plot_sites = None, force_ylim = None, width=None):
    flat_res = res.pivot(index=['antibody', 'site'], columns='mutation', values='mut_escape').fillna(0)
    sites_mean_score = flat_res.mean(axis=1)
    sites_total_score = flat_res.sum(axis=1)
    _ = sites_total_score[sites_total_score>site_thres].index
    strong_sites = np.unique(np.array(sorted([i[1] for i in _])))
    print(strong_sites)

    plot_sites = strong_sites
    plot_sites = plot_sites[plot_sites < 520].astype(int)
    print(plot_sites)

    if force_plot_sites is not None:
        plot_sites = force_plot_sites

    flat_res = flat_res[flat_res.index.isin(plot_sites, level=1)]

    _ = pd.DataFrame(sites_total_score)
    _.columns = ['value']
    _['site'] = [i[1] for i in _.index]
    _['antibody'] = [i[0] for i in _.index]

    if rownames is not None:
        Abs = rownames
    else:
        Abs = np.unique([i[0] for i in flat_res.index])
    print(Abs)
    Npages = len(Abs)//10 + 1
    if width is None:
        width=30
    with PdfPages(prefix+'_aa_logo.pdf') as pdf:
        for p in range(Npages):
            Abs_p = Abs[p*10:min(len(Abs),(p+1)*10)]
            fig = plt.figure(figsize=(width,len(Abs_p)*4.6)).subplots_adjust(wspace=0.2,hspace=0.5)
            site2pos = {}
            for i in range(len(plot_sites)):
                site2pos[plot_sites[i]] = i

            for i in range(len(Abs_p)):
                ab = Abs_p[i]
                _ = flat_res.query('antibody == @ab').droplevel(0)
                add_sites = np.setdiff1d(plot_sites, _.index)
                for _site in add_sites:
                    _.loc[_site,:] = 0.0
                _ = _.sort_index()
                _.index = range(len(_))
                ax = plt.subplot(len(Abs_p), 1, i+1)
                logo = logomaker.Logo(_,
                               ax=ax, 
                               color_scheme='dmslogo_funcgroup', 
                               vpad=.1, 
                               width=.8)
                logo.style_xticks(anchor=1, spacing=1, rotation=90, fontsize=16)
                _max = np.sum(_.to_numpy(), axis=1).max()
                # ax.set_xticklabels(plot_sites[1::2])
                ax.set_xticklabels(plot_sites)
                
                # ax.set_yticks([])
                ax.tick_params(axis='both', which='both', length=0)
                if force_ylim is not None:
                    ax.set_ylim(0.0,force_ylim)
                ax.yaxis.set_tick_params(labelsize=20)
                if ab in shownames:
                    ax.set_title(shownames[ab], fontsize=8, fontweight="bold")
                else:
                    ax.set_title(ab, fontsize=8, fontweight="bold")
            pdf.savefig()
            plt.close()

def do_calc(scores_r, ag, targets, neut2se, abinfo, A_adv = True, A_codon = True, A_neut = True,
            E=1.0, B=1.0, use_log=False, use_max=False, use_norm=False,
            logo=False, return_df=False, use_codon_weights=False,
            title=None,blacklist=None
           ):
    # scores_r = targets[ag]['dms']
    neut_data = abinfo[targets[ag]['neut']].to_dict()
    src_dict = abinfo['source'].to_dict()
    
    single_mut_effects = pd.read_csv(neut2se[ag]).assign(
        coef=lambda x: [y for y in (np.tanh(x['expr_avg']*(x['expr_avg']<0)*E)+1)*(np.tanh(x['bind_avg']*B)+1)]
    ) 
    single_mut_effects.index = single_mut_effects['site'].astype('str') + single_mut_effects['mutation']
    single_mut_effects = single_mut_effects['coef'].to_dict()

    if use_codon_weights is False:
        use_codon = pd.read_csv(targets[ag]['codon_mut1'])
        _umuts = {}
        for i in range(len(use_codon)):
            _ms = use_codon['mut1'][i]
            for x in _ms:
                _umuts[str(use_codon['pos'][i])+x] = 1.0
    else:
        _umuts = pd.read_csv(targets[ag]['codon_weights']).assign(
            mutation = lambda x: x['pos'].astype(str)+x['mut']).set_index('mutation')['weight'].to_dict()
    
    
    scores = scores_r.assign(site_mut = lambda x: x['site'].astype(str)+x['mutation']).assign(
        adv_weight = (lambda x: [single_mut_effects[y] if (y in single_mut_effects and not np.isnan(single_mut_effects[y])) else 0.0 for y in x['site_mut']]) if A_adv else 1.0,
        codon_weight = (lambda x: [(_umuts[y] if y in _umuts else 0.0) for y in x['site_mut'].to_list()]) if A_codon else 1.0
    )
    if use_norm:
        scores = scores.assign(escape_max = lambda x: x.groupby('antibody')['mut_escape'].transform('max')).assign(
            mut_escape = lambda x: x['mut_escape']/x['escape_max']).drop(columns=['escape_max'])
    if use_log:
        scores = scores.assign(neut_weight = lambda x: [(0.0 if np.isnan(neut_data[y]) else max(0.0,np.log10(1/min(1,neut_data[y])))) if A_neut else 1.0 for y in x['antibody']])
    else:
        scores = scores.assign(neut_weight = lambda x: [(0.0 if np.isnan(neut_data[y]) else 1.0/neut_data[y]) if A_neut else 1.0 for y in x['antibody']])
    
    scores = scores.assign(
        mut_escape_adj = lambda x: x['mut_escape'] * x['neut_weight'] * x['adv_weight'] * x['codon_weight']
    )
    _title = (
              'weight: '+ag+' expr_bind:'+str(A_adv)+
              ' codon:'+str(A_codon)+' log:'+str(use_log)+
              ' norm:'+str(use_norm)+' max:'+str(use_max)+
              ' Expr:'+str(E)+' Bind:'+str(B)) if title is None else title
    
    if blacklist is not None:
        scores = scores.assign(site_mut = lambda x: x['site'].astype(str)+x['mutation']).query('site_mut not in @blacklist')
    
    if logo:
        scores = scores.groupby(['site','mutation']).sum()['mut_escape_adj'].reset_index().assign(antibody=_title)
        scores['mut_escape_adj'] = scores['mut_escape_adj']/scores['mut_escape_adj'].max()
        return scores
    
    if use_max:
        site_avg = scores.groupby(['site', 'antibody']).max()['mut_escape_adj'].reset_index().groupby('site').sum().reset_index()
    else:
        site_avg = scores.groupby(['site', 'mutation']).sum()['mut_escape_adj'].reset_index().groupby('site').sum().reset_index()
    site_avg['mut_escape_adj'] = site_avg['mut_escape_adj']/site_avg['mut_escape_adj'].max()
    
    if return_df:
        return site_avg.assign(
            absrc = ag, 
            # absrc = '+'.join(use_ab_src) if use_ab_src is not None else 'ALL', 
            weight = ag, is_expr_bind = A_adv, is_codon = A_codon, 
            is_neut_log = use_log, is_norm = use_norm, is_max = use_max, expr_coef = E, bind_coef = B
        )
    else:
        raise ValueError("Not implemented")
    
def get_enrichment_stat(y_true, y_pred, true_thres=0.05):
    res = {
        'thres':[],
        'TP':[],
        'FP':[],
        'TN':[],
        'FN':[],
        'precision':[],
        'recall':[],
        'FPR': [],
    }
    N = len(y_true)
    y_true_bin = (y_true < true_thres)
    for thres in sorted(pd.unique(y_pred)):
        y_pred_bin = (y_pred < thres)
        
        res['thres'].append(thres)
        res['TP'].append(sum([y_true_bin[i] and y_pred_bin[i] for i in range(N)]))
        res['FP'].append(sum([not y_true_bin[i] and y_pred_bin[i] for i in range(N)]))
        res['TN'].append(sum([not y_true_bin[i] and not y_pred_bin[i] for i in range(N)]))
        res['FN'].append(sum([y_true_bin[i] and not y_pred_bin[i] for i in range(N)]))
        res['precision'].append(res['TP'][-1]/(res['TP'][-1]+res['FP'][-1]) if (res['TP'][-1]+res['FP'][-1] > 0) else pd.NA)
        res['recall'].append(res['TP'][-1]/(res['TP'][-1]+res['FN'][-1]) if (res['TP'][-1]+res['FN'][-1] > 0) else pd.NA)
        res['FPR'].append(res['FP'][-1]/(res['FP'][-1]+res['TN'][-1]) if (res['FP'][-1]+res['TN'][-1] > 0) else pd.NA)
    y_pred_bin = [True]*N
    res['thres'].append(thres)
    res['TP'].append(sum([y_true_bin[i] and y_pred_bin[i] for i in range(N)]))
    res['FP'].append(sum([not y_true_bin[i] and y_pred_bin[i] for i in range(N)]))
    res['TN'].append(sum([not y_true_bin[i] and not y_pred_bin[i] for i in range(N)]))
    res['FN'].append(sum([y_true_bin[i] and not y_pred_bin[i] for i in range(N)]))
    res['precision'].append(res['TP'][-1]/(res['TP'][-1]+res['FP'][-1]) if (res['TP'][-1]+res['FP'][-1] > 0) else pd.NA)
    res['recall'].append(res['TP'][-1]/(res['TP'][-1]+res['FN'][-1]) if (res['TP'][-1]+res['FN'][-1] > 0) else pd.NA)
    res['FPR'].append(res['FP'][-1]/(res['FP'][-1]+res['TN'][-1]) if (res['FP'][-1]+res['TN'][-1] > 0) else pd.NA)
    res = pd.DataFrame(res).dropna().reset_index(drop=True)
    auc = 0
    _rec = 0
    _y = 0
    for i in res.sort_values('FPR').index:
        auc += (res.loc[i, 'FPR']-_rec)*res.loc[i, 'recall']-(res.loc[i, 'FPR']-_rec)*(res.loc[i, 'recall']-_y)
        _rec = res.loc[i, 'FPR']
        _y = res.loc[i, 'recall']
    auc += (1-_rec)*res.loc[res.index[-1],'recall']
    return auc, res