import pandas as pd
from utils import neut_format_parser
from pathlib import Path

srcdir = Path('../source_data')
tgdir = Path('../processed_source_data')

def prepare_fig1():
    tg = tgdir/'Fig1'
    tg.mkdir(exist_ok=True, parents=True)
    data = pd.read_csv(srcdir/'abinfo_table.csv', index_col=0)
    for s in data.columns:
        if '_IC50' in s:
            data[s] = [neut_format_parser(x) for x in data[s]]

    # rename source cohort
    src_map = {
        'WT convalescents':'WT', 'WT vaccinees':'WT', 'SARS exposure':'SARS+WT', 'BA.1 BTI':'BA.1 BTI',
        'BA.2 BTI':'BA.2 BTI', 'BA.5 BTI':'BA.5/BF.7 BTI', 'BF.7 BTI':'BA.5/BF.7 BTI',
        'BA.1 BTI + BA.5/BF.7 infection':'BA.1 BTI + BA.5/BF.7', 'BA.2 BTI + BA.5/BF.7 infection':'BA.2 BTI + BA.5/BF.7',
        'BA.1 + BA.5/BF.7 without vaccine':'BA.1 + BA.5/BF.7', 'BA.2 + BA.5/BF.7 infection':'BA.2 + BA.5/BF.7',
            
        'BA.5 BTI + XBB infection': 'BA.5 BTI + XBB', 
        'BA.5 + XBB infection': 'BA.5 + XBB',
    }
    use_src = list(src_map.keys())

    data_subs = data.query("source in @use_src").copy()
    data_subs['source'] = [src_map[i] for i in data_subs['source']]

    auto_IC50 = {
        'WT': 'D614G_IC50', 
        'SARS+WT': 'D614G_IC50', 
        'BA.1 BTI': 'BA1_IC50', 
        'BA.2 BTI': 'BA2_IC50', 
        'BA.5/BF.7 BTI': 'BA5_IC50', 
        
        'BA.1 BTI + BA.5/BF.7': 'BA5_IC50', 
        'BA.2 BTI + BA.5/BF.7': 'BA5_IC50',
        'BA.1 + BA.5/BF.7': 'BA5_IC50', 
        'BA.2 + BA.5/BF.7': 'BA5_IC50',
        
        'BA.5 BTI + XBB': 'XBB1_5_IC50', 
        'BA.5 + XBB': 'XBB1_5_IC50',
        'XBB BTI': 'XBB1_5_IC50', 
        'XBB (no vaccine)': 'XBB1_5_IC50',
    }


    data_subs = data_subs.assign(
        # best_IC50=lambda df: [df.loc[i,terms].min() for i in df.index],
        auto_IC50=lambda df: [df.loc[i,auto_IC50[df.loc[i,'source']]] for i in df.index],
        # num_effective=lambda df: [(df.loc[i, terms] < 0.05).sum() for i in df.index],
        # tier=lambda df: [get_tier(df.loc[i, terms].to_dict()) for i in df.index]
    ).query('not auto_IC50.isna()')

    terms = ['D614G_IC50', 'BA1_IC50', 'BA2_IC50', 'BA5_IC50', 'BQ1_1_IC50', 'XBB1_5_IC50', 'HK3_1_IC50', 'JN1_IC50', 'KP3_IC50', 'SARS_IC50']
    data_subs[['source','auto_IC50']+terms].to_csv(tg/"all_mAbs_IC50.csv")
    data_subs = data_subs.query(
        'auto_IC50 < 0.05')
    
    x = data_subs[['source']+terms].reset_index()
    x.columns = ['id', 'source', 'D614G', 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5', 'HK.3.1', 'JN.1', 'KP.3','SARS']
    # groups = pd.read_csv("/gshare/xielab/jianfc/COVID/temp_tasks/202311_XBB_response_MS/DMS/dms_groups.csv", index_col=0)[['WT_group','BA.5_group','XBB.1.5_group']]
    # x.merge(groups, left_on='id', right_index=True, how='left').to_csv("use_neutralization_table.csv", index=None)

    x = x.melt(id_vars=['id','source'])
    x.columns = ['id', 'source', 'variant', 'IC50']
    x.to_csv(tg/"selected_Ab_clean_IC50.csv", index=None)

import numpy as np
def prepare_fig2():
    tg = tgdir/'Fig2'
    plot_dir = tg/'..'/'..'/'plot_scripts'/'plots'/'Fig2'
    plot_dir.mkdir(exist_ok=True, parents=True)
    tg.mkdir(exist_ok=True, parents=True)

    WT_DMS_CSV = pd.read_csv(srcdir/"dms_scores_clean.csv")
    ABINFO = pd.read_csv(srcdir/'abinfo_table.csv', index_col=0)
    targets = {
        "WT": {
            "dms": WT_DMS_CSV,
            "neut": "D614G_IC50",
            "src": ["WT convalescents", "WT vaccinees"],
            "codon_weights": srcdir/"mutation_weights_SARSCoV2_WuhanHu1_Spike.csv",
            "codon_mut1": srcdir/"mut1_SARSCoV2_WuhanHu1_Spike.csv"
        },
    }
    neut2se = {
        "WT": srcdir/"JBloom_new_WT_single_mut_effects.csv",
    }

    strains = list(np.unique([targets[x]['neut'] for x in targets]))

    data = ABINFO[[
        "source", *strains
    ]].assign(antibody=lambda x: x.index)
    for s in strains:
        data[s] = [neut_format_parser(x, 10) for x in data[s]]

    eff_coef = {
        'expr': {"WT":1.0},
        'bind': {"WT":1.0},
    }

    exclude_mutations = {
        'WT': set(),
    }

    # add mutations causing glycan change to exclusion
    for ag in targets:
        wt_res = pd.read_csv(targets['WT']['codon_mut1'], index_col=1)['wt'].to_dict()
        ag_res = pd.read_csv(targets[ag]['codon_mut1'], index_col=1)['wt'].to_dict()
        for i in ag_res:
            if wt_res[i] != ag_res[i]:
                exclude_mutations[ag].add(str(i)+wt_res[i])
        if not (ag in exclude_mutations):
            continue
        wt = pd.read_csv(targets[ag]['codon_mut1'])['wt'].to_list()
        for i in range(len(wt)):
            if i+2 < len(wt) and wt[i] == 'N' and (wt[i+2] != 'S' and wt[i+2] != 'T'):
                exclude_mutations[ag].add(str(331+i+2)+'S')
                exclude_mutations[ag].add(str(331+i+2)+'T')
            if i+2 < len(wt) and (wt[i+2] == 'S' or wt[i+2] == 'T') and wt[i] != 'N':
                exclude_mutations[ag].add(str(331+i)+'N')

    from utils import do_calc, plot_res_logo
    ag = 'WT'
    blacklist = exclude_mutations[ag]
    use_ab_src = targets[ag]['src']
    _uabs = set(data.query('source in @use_ab_src').index.to_list())
    scores_orig = targets[ag]['dms'].query('antibody in @_uabs')

    plot_df = []

    df = do_calc(scores_orig, ag, targets, neut2se, data,
                 A_adv = True, A_codon = True, A_neut=True, use_log=True,
                    E=eff_coef['expr'][ag], B=eff_coef['bind'][ag], 
            use_norm=True, use_max=False, use_codon_weights=True, logo=True, blacklist=blacklist).rename(
        columns={'mut_escape_adj':'mut_escape'})
    plot_res_logo(df, f"{plot_dir}/Fig2a-calc-logo-WT", site_thres=0.15, width=10)
    df.to_csv(tg/"SourceData-calc-logo-WT.csv")

    ts = do_calc(scores_orig, ag, targets, neut2se, data,
                 A_adv = True, A_codon = True, A_neut=True, use_log=True,
                        E=eff_coef['expr'][ag], B=eff_coef['bind'][ag], 
                use_norm=True, use_max=False, use_codon_weights=True, logo=False, return_df=True, blacklist=blacklist).assign(
        target=ag,
    )
    ts.to_csv(tg/"site-calc-sum-WT.csv", index=None)

    data_subs = pd.read_csv(srcdir/"Fig2-WT-NAbs-neut.csv")
    from utils import get_enrichment_stat
    from statsmodels.stats.proportion import proportion_confint
    preds = ['D614G_IC50', 'filter_worst']
    targs = ['BA5_IC50', 'XBB1_5_IC50', "JN1_IC50"]
    thres = 0.05
    res = []
    for targ in targs:
        for pred in preds:
            N_effect = 0
            N_total = 0
            ratio = []
            Eff = []
            _tmp_data = data_subs.sort_values(pred).reset_index(drop=True)
            for x in _tmp_data.index:
                N_total += 1
                N_effect += (_tmp_data.loc[x, targ] < thres)
                ratio.append(N_effect / N_total)
                Eff.append(N_effect)
            res.append(
                pd.DataFrame({'id':range(len(ratio)), 'N_effect':Eff, 'N_total':range(1,N_total+1), 'ratio':ratio}).assign(metric = pred, target=targ)
            )
    res = pd.concat(res)
    CI = proportion_confint(res['N_effect'], res['N_total'], alpha=0.05, method='normal')
    res = res.assign(low=CI[0], high=CI[1])
    res.to_csv(tg/"enrichment_curve_data_wt_sars.csv", index=None)
    
    targs = ['BA5_IC50', 'XBB1_5_IC50', "JN1_IC50"]
    res = []
    for targ in targs:
        for pred in preds:
            N_effect = 0
            N_total = 0
            ratio = []
            Eff = []
            _tmp_data = data_subs.query('source != "SARS"').sort_values(pred).reset_index(drop=True)
            for x in _tmp_data.index:
                N_total += 1
                N_effect += (_tmp_data.loc[x, targ] < thres)
                ratio.append(N_effect / N_total)
                Eff.append(N_effect)
            res.append(
                pd.DataFrame({'id':range(len(ratio)), 'N_effect':Eff, 'N_total':range(1,N_total+1), 'ratio':ratio}).assign(metric = pred, target=targ)
            )
    res = pd.concat(res)
    CI = proportion_confint(res['N_effect'], res['N_total'], alpha=0.05, method='normal')
    res = res.assign(low=CI[0], high=CI[1])
    res.to_csv(tg/"enrichment_curve_data_wt.csv", index=None)

    for WT_only in [False, True]:
        x_name = []
        y_name = []
        val = []
        for metric in [["D614G_IC50"], ["D614G+E484K"], 
                    ["D614G+E484K", "D614G-S1"], 
                    ["D614G+E484K", "D614G-S1", "D614G-S2"], 
                    ["D614G+E484K", "D614G-S1", "D614G-S2", "D614G-S3"],
                    ["D614G+E484K", "D614G-S1", "D614G-S2", "D614G-S3", "D614G-S4", "D614G-S5"],
                    ]:
            if len(metric) > 0:
                _tmp_data = data_subs.assign(metric_best = lambda x: [np.max(x.loc[y,metric]) for y in x.index]).query('metric_best < 0.05')
            else:
                _tmp_data = data_subs.copy()

            if WT_only:
                _tmp_data = _tmp_data.query("source == 'WT'")
            for target in [[], ["BA1_IC50", "BA2_IC50"], ["BA1_IC50", "BA2_IC50", "BA5_IC50"], ["BA1_IC50", "BA2_IC50", "BA5_IC50", "BQ1_1_IC50", "XBB1_5_IC50"]]:
                if len(target) > 0:
                    _tg_data = _tmp_data.assign(target_best = lambda x: [np.max(x.loc[y, target]) for y in x.index]).query("target_best < 0.05")
                else:
                    _tg_data = _tmp_data
                x_name.append('_'.join(metric))
                y_name.append('_'.join(target) if len(target) > 0 else "all")
                val.append(len(_tg_data))

        pd.DataFrame({'indicator':x_name, 'target':y_name, 'value': val}).to_csv(tg/f"fig2d_WT_only_{WT_only}.csv")

    from scipy.stats import fisher_exact

    for WT_only in [False, True]:
        x_name = []
        y_name = []
        pvalue = []
        for metric in [["D614G+E484K"], ["D614G-S1"], ["D614G-S2"], ["D614G-S3"], ["D614G-S4"], ["D614G-S5"], ["D614G-S1", "D614G-S2", "D614G-S3", "D614G-S4", "D614G-S5"]]:
            for target in ['BA1_IC50', 'BA2_IC50', 'BA5_IC50', 'BQ1_1_IC50', 'XBB1_5_IC50']:
                if WT_only:
                    use_pred = data_subs.query("source == 'WT'")
                else:
                    use_pred = data_subs
                use_pred = use_pred.assign(
                    is_pred_pass=lambda x: [np.max(x.loc[y,metric]) < 0.05 for y in x.index],
                    is_broad_pass=lambda df: df[target] < 0.05,
                ).groupby(['is_pred_pass', 'is_broad_pass'])['source'].count().reset_index().pivot(columns="is_broad_pass", index="is_pred_pass", values="source")

                odd_r, pval = fisher_exact(use_pred.fillna(0))

                x_name.append('_'.join(metric))
                y_name.append(target)
                pvalue.append(pval)

        pd.DataFrame({'filter':x_name, 'target':y_name, 'pval': pvalue}).to_csv(tg/f"fig2e_WT_only_{WT_only}.csv")

def prepare_fig3():
    plot_dir = srcdir/'..'/'plot_scripts'/'plots'/'Fig3'
    plot_dir.mkdir(exist_ok=True, parents=True)

    dms_profiles = {
        'BA.5' : pd.read_csv(srcdir/"Ext-plot_dms_scores_BA.5.csv"),
        'XBB.1.5' : pd.read_csv(srcdir/"Fig3-plot_dms_scores_XBB.1.5.csv"),
    }

    plot_sites = {
        'BA.5': None,
        'XBB.1.5': [403, 404, 420, 439, 447, 453, 455, 456, 472, 473,474,475,476, 485, 487, 489, 491, 498, 499, 500, 501, 502, 503, 504, 505],
    }
    from utils import plot_res_logo
    for ag in dms_profiles:
        plot_res_logo(dms_profiles[ag], prefix=f"{plot_dir}/plot_logo_mAbs_{ag}", site_thres=1, width=6, force_plot_sites = plot_sites[ag])

def prepare_fig5():
    tg = tgdir/'Fig5'
    tg.mkdir(exist_ok=True, parents=True)
    plot_dir = tg/'..'/'..'/'plot_scripts'/'plots'/'Fig5'
    plot_dir.mkdir(exist_ok=True, parents=True)

    from scipy.optimize import curve_fit
    # fit ELISA kinetics to bi-exponential model
    def func(t, a, k3, k):
        return a*np.exp(-k3*t)*(1.0-np.exp(-k*t))

    def func_log(t, a, k3, k):
        return -k3*t+np.log(1.0-np.exp(-k*t))+a

    data = pd.read_csv(srcdir/"Fig5-mice-data-elisa.csv").dropna()
    data['group'] = data['type']+'-'+data['sex']
    animals = pd.unique(data['mouse'])

    animal_groups = data[['mouse', 'group']].drop_duplicates()
    x = np.linspace(min(data['hours'])-1, max(data['hours'])+1, 1000)
    fit_data = []
    fit_params = []

    for animal in animals:
        datam = data.query('mouse == @animal')
        popt, pcov = curve_fit(func_log, datam['hours'], np.log(datam['conc']), p0=[10, 0.001, 0.001], maxfev=20000)
        fit_data.append(pd.DataFrame({'hours': x, 'conc': [np.exp(func_log(i, *popt)) for i in x], 'mouse': animal}))
        # popt, pcov = curve_fit(func, datam['hours'], datam['conc'], p0=[2000, 0.001, 0.001], maxfev=20000)
        # fit_data.append(pd.DataFrame({'hours': x, 'conc': [func(i, *popt) for i in x], 'mouse': animal}))
        fit_params.append(popt)

    fit_data = pd.concat(fit_data)
    fit_params = pd.DataFrame(fit_params, index=animals, columns=['fit_a', 'fit_k3', 'fit_k']).merge(animal_groups, left_index=True, right_on='mouse', how='outer')

    fit_data.to_csv(tg/"ELISA_fit_data.csv", index=None)
    fit_params.to_csv(tg/"ELISA_fit_params.csv", index=None)
    def gmean(x):
        return np.exp(np.mean(np.log(x)))

    def gstd(x):
        return np.exp(np.std(np.log(x)))

    data_agg = data.groupby(["hours", "group"])['conc'].agg([np.mean, np.std, gmean, gstd]).reset_index()

    fit_data = []
    fit_params = []
    groups = pd.unique(data_agg['group'])

    for _g in groups:
        datam = data_agg.query("group == @_g")
        popt, pcov = curve_fit(func_log, datam['hours'], np.log(datam['gmean']), p0=[10, 0.001, 0.001], maxfev=20000)
        fit_data.append(pd.DataFrame({'hours': x, 'conc': [np.exp(func_log(i, *popt)) for i in x], 'group': _g}))
        # popt, pcov = curve_fit(func, datam['hours'], datam['gmean'], p0=[2000, 0.001, 0.001], maxfev=20000)
        # fit_data.append(pd.DataFrame({'hours': x, 'conc': [func(i, *popt) for i in x], 'group': _g}))
        fit_params.append(popt)

    fit_data = pd.concat(fit_data)
    fit_params = pd.DataFrame(fit_params, index=groups, columns=['fit_a', 'fit_k3', 'fit_k'])
    fit_params.index.name = 'group'

    data_agg.to_csv(tg/"ELISA_data_agg.csv", index=None)
    fit_data.to_csv(tg/"ELISA_fit_data_agg.csv", index=None)
    fit_params.to_csv(tg/"ELISA_fit_params_agg.csv")

def prepare_extfigs():
    # ExtFig2
    tg = tgdir/'ExtFig2'
    tg.mkdir(exist_ok=True, parents=True)
    
    neut = pd.read_csv(srcdir/'Fig2-WT-NAbs-neut.csv')
    single = pd.read_csv(srcdir/'ExtFig2-single_muts_neut.csv')
    filters = single.columns[single.columns.str.startswith('D614G+')].tolist() + ['D614G+E484K']
    neut = neut.merge(
        single, on='id', how='left'
    ).assign(
        filter_worst_single = lambda x: x[filters].max(axis=1)
    )

    neut = neut.query('source == "WT"').sort_values('filter_worst_single').rename(
        columns={k: k.replace('D614G+', 'B.1+') for k in filters}
    )

    neut.to_csv(tg/'neutralization_single.csv', index=False)

    # enrich
    x_name = []
    y_name = []
    val = []
    for metric in [["D614G_IC50"], ["B.1+E484K"], 
                    ["filter_worst_single"],
                    ]:
        if len(metric) > 0:
            _tmp_data = neut.assign(metric_best = lambda x: [np.max(x.loc[y,metric]) for y in x.index]).query('metric_best < 0.05')
        else:
            _tmp_data = neut.copy()

        for target in [[], ["BA1_IC50", "BA2_IC50"], ["BA1_IC50", "BA2_IC50", "BA5_IC50"], ["BA1_IC50", "BA2_IC50", "BA5_IC50", "BQ1_1_IC50", "XBB1_5_IC50"]]:
            if len(target) > 0:
                _tg_data = _tmp_data.assign(target_best = lambda x: [np.max(x.loc[y, target]) for y in x.index]).query("target_best < 0.05")
            else:
                _tg_data = _tmp_data
            x_name.append('_'.join(metric))
            y_name.append('_'.join(target) if len(target) > 0 else "all")
            val.append(len(_tg_data))

    pd.DataFrame({'indicator':x_name, 'target':y_name, 'value': val}).to_csv(tg/"enrich_bars_data.csv")
    # ExtFig5
    tg = tgdir/'ExtFig5'
    plot_dir = tg/'..'/'..'/'plot_scripts'/'plots'/'ExtFig5'
    plot_dir.mkdir(exist_ok=True, parents=True)
    tg.mkdir(exist_ok=True, parents=True)
    from scipy.optimize import curve_fit

    def hill_func(x, EC50, n):
        return 1/(1+(EC50/x)**n)
        
    data = pd.read_csv(srcdir/"ExtFig5-ACE2_neutralization.csv")
    x = np.exp(np.linspace(np.log(min(data['conc'])/0.9), np.log(max(data['conc'])*1.1), 500))
    variants = pd.unique(data['variant'])
    fit_params = []
    fit_data = []

    for _var in variants:
        datam = data.groupby(["variant","conc"])['inhibition'].mean().reset_index().query("variant == @_var")
        popt, _ = curve_fit(hill_func, datam['conc'], datam['inhibition'], p0=[0.1, 1], maxfev=20000)
        fit_params.append([_var, *popt])
        fit_data.append(pd.DataFrame({'variant': _var, 'conc':x, 'inhibition': [hill_func(i, *popt) for i in x]}))
        
    fit_params = pd.DataFrame(fit_params, columns=['variant','EC50', 'n'])
    fit_data = pd.concat(fit_data)

    fit_params.to_csv(tg/"ACE2_neut_fit_params.csv", index=None)
    fit_data.to_csv(tg/"ACE2_neut_fit_data.csv", index=None)
    # fit for each rep
    x = np.exp(np.linspace(np.log(min(data['conc'])/0.9), np.log(max(data['conc'])*1.1), 500))

    variants = pd.unique(data['variant'])
    fit_params = []
    fit_data = []

    for _var in variants:
        datam = data.query("variant == @_var")
        reps = pd.unique(datam['rep'])
        for _rep in reps:
            _data = datam.query("rep == @_rep")
            popt, _ = curve_fit(hill_func, _data['conc'], _data['inhibition'], p0=[0.1, 1], maxfev=20000)
            fit_params.append([_var, _rep, *popt])
            fit_data.append(pd.DataFrame({'variant': _var, 'rep': _rep, 'conc':x, 'inhibition': [hill_func(i, *popt) for i in x]}))
        
    fit_params = pd.DataFrame(fit_params, columns=['variant','rep','EC50', 'n'])
    fit_data = pd.concat(fit_data)

    fit_params.to_csv(tg/"ACE2_neut_fit_params_ind.csv", index=None)
    fit_data.to_csv(tg/"ACE2_neut_fit_data_ind.csv", index=None)
    
    tg = tgdir/'ExtFig6'
    tg.mkdir(exist_ok=True, parents=True)
    def hill_func(x, EC50, n):
        return 1/(1+(EC50/x)**n)

    data = pd.read_csv(srcdir/"ExtFig6d-S490Y_neut_raw.csv")
    VC = data.groupby('variant')['VC'].mean().reset_index()
    data = data.merge(VC, left_on='variant', right_on='variant', how='left')

    data['rep1'] = 1.0 - data['rep1'] / data['VC_y']
    data['rep2'] = 1.0 - data['rep2'] / data['VC_y']


    data = data[['variant', 'conc', 'rep1', 'rep2']].melt(id_vars=['variant', 'conc'])

    data.columns = ['variant', 'conc', 'rep', 'inhibition']
    data.to_csv(tg/"490_neut_data_clean.csv", index=None)
    x = np.exp(np.linspace(np.log(min(data['conc'])/0.9), np.log(max(data['conc'])*1.1), 500))
    variants = pd.unique(data['variant'])
    fit_params = []
    fit_data = []

    for _var in variants:
        datam = data.groupby(["variant","conc"])['inhibition'].mean().reset_index().query("variant == @_var")
        popt, _ = curve_fit(hill_func, datam['conc'], datam['inhibition'], p0=[0.1, 1], maxfev=20000)
        fit_params.append([_var, *popt])
        fit_data.append(pd.DataFrame({'variant': _var, 'conc':x, 'inhibition': [hill_func(i, *popt) for i in x]}))
        
    fit_params = pd.DataFrame(fit_params, columns=['variant','EC50', 'n'])
    fit_data = pd.concat(fit_data)

    fit_params.to_csv(tg/"490_neut_fit_params.csv", index=None)
    fit_data.to_csv(tg/"490_neut_fit_data.csv", index=None)

    x = np.exp(np.linspace(np.log(min(data['conc'])/0.9), np.log(max(data['conc'])*1.1), 500))

    variants = pd.unique(data['variant'])
    fit_params = []
    fit_data = []

    for _var in variants:
        datam = data.query("variant == @_var")
        reps = pd.unique(datam['rep'])
        for _rep in reps:
            _data = datam.query("rep == @_rep")
            popt, _ = curve_fit(hill_func, _data['conc'], _data['inhibition'], p0=[0.1, 1], maxfev=20000)
            fit_params.append([_var, _rep, *popt])
            fit_data.append(pd.DataFrame({'variant': _var, 'rep': _rep, 'conc':x, 'inhibition': [hill_func(i, *popt) for i in x]}))
        
    fit_params = pd.DataFrame(fit_params, columns=['variant','rep','EC50', 'n'])
    fit_data = pd.concat(fit_data)

    fit_params.to_csv(tg/"490_neut_fit_params_ind.csv", index=None)
    fit_data.to_csv(tg/"490_neut_fit_data_ind.csv", index=None)

def main():
    prepare_fig1()
    prepare_fig2()
    prepare_fig3()
    prepare_fig5()
    prepare_extfigs()

if __name__ == "__main__":
    main()