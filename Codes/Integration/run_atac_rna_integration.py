#!/usr/bin/env python3
"""Integrate ATAC/RNA changes using representative peaks and direction-specific gene groups."""
import argparse
from pathlib import Path
import warnings
import matplotlib
matplotlib.use("Agg")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_atac_csv", type=Path)
    parser.add_argument("input_rna_csv", type=Path)
    parser.add_argument("input_mapping_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    for path in (args.input_atac_csv, args.input_rna_csv, args.input_mapping_tsv):
        if not path.is_file():
            parser.error(f"Input missing: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    correlation_rows, group_rows, count_rows, test_rows = [], [], [], []
    import pandas as pd, numpy as np

    atac = pd.read_csv(args.input_atac_csv)
    rna  = pd.read_csv(args.input_rna_csv)
    mapping = pd.read_csv(args.input_mapping_tsv, sep='\t')

    # RNA: FBgn -> symbol
    rna = rna.merge(mapping, left_on='FBgn', right_on='FLYBASE', how='left')
    rna = rna.dropna(subset=['SYMBOL']).drop_duplicates('SYMBOL')
    rna = rna.rename(columns={'log2FoldChange':'rna_l2fc','padj':'rna_padj'})[['SYMBOL','rna_l2fc','rna_padj']]

    # TSS peak-gene (matched_gene)
    tss = atac[(atac['TSS']=='O') & (atac['matched_gene']!='-') & (atac['matched_gene'].notna())].copy()
    tss = tss.assign(gene=tss['matched_gene'].str.split(';')).explode('gene')
    tss['gene'] = tss['gene'].str.strip(); tss = tss[tss['gene']!='']
    tss = tss[['PeakID','gene','log2FoldChange','padj','ClusterLabel']].rename(
        columns={'log2FoldChange':'atac_l2fc','padj':'atac_padj'})

    # Enhancer peak-gene (enhancer column)
    enh = atac[(atac['TSS']=='X') & (atac['enhancer']!='-') & (atac['enhancer'].notna())].copy()
    enh = enh.assign(gene=enh['enhancer'].str.split(';')).explode('gene')
    enh['gene'] = enh['gene'].str.strip(); enh = enh[enh['gene']!='']
    enh = enh[['PeakID','gene','log2FoldChange','padj','ClusterLabel']].rename(
        columns={'log2FoldChange':'atac_l2fc','padj':'atac_padj'})

    # Select the peak with the smallest adjusted p-value for each gene.
    tss_gene = tss.sort_values('atac_padj').drop_duplicates('gene').reset_index(drop=True)
    enh_gene = enh.sort_values('atac_padj').drop_duplicates('gene').reset_index(drop=True)

    tss_gene.to_csv(args.output_dir / 'output_tss_representative_peaks.csv', index=False)
    enh_gene.to_csv(args.output_dir / 'output_enhancer_representative_peaks.csv', index=False)
    rna.to_csv(args.output_dir / 'output_rna_gene_symbols.csv', index=False)
    ATAC_PADJ=0.05; ATAC_L2FC=np.log2(1.5)
    RNA_PADJ=0.05;  RNA_L2FC=np.log2(1.5)


    import matplotlib.pyplot as plt
    from scipy import stats
    from adjustText import adjust_text

    plt.rcParams['pdf.fonttype']=42; plt.rcParams['font.family']='DejaVu Sans'
    COL_UP='#B2182B'; COL_DOWN='#2166AC'; COL_NS='#B0B0B0'

    def make_scatter(atac_gene, title, out_pdf):
        d = atac_gene.merge(rna, left_on='gene', right_on='SYMBOL', how='inner').dropna(subset=['rna_l2fc'])
        atac_sig = (d['atac_padj']<ATAC_PADJ)&(d['atac_l2fc'].abs()>ATAC_L2FC)
        rna_sig  = d['rna_padj'].notna()&(d['rna_padj']<RNA_PADJ)&(d['rna_l2fc'].abs()>RNA_L2FC)
        both = atac_sig & rna_sig
        d['cat']='ns'
        d.loc[both&(d['atac_l2fc']>0),'cat']='up'
        d.loc[both&(d['atac_l2fc']<0),'cat']='down'

        # Pearson correlation and regression use genes significant in both assays.
        sig_pts = d[both]
        if len(sig_pts) < 2 or sig_pts['atac_l2fc'].nunique() < 2 or sig_pts['rna_l2fc'].nunique() < 2:
            raise ValueError(f'{title}: correlation needs at least two nonconstant significant pairs')
        r_val, p_val = stats.pearsonr(sig_pts['atac_l2fc'], sig_pts['rna_l2fc'])
        d.to_csv(args.output_dir / f'output_{title.lower()}_scatter_data.csv', index=False)
        correlation_rows.append({'region': title, 'n_joined': len(d), 'n_both_significant': len(sig_pts),
                                 'pearson_r': r_val, 'pvalue': p_val})

        fig, ax = plt.subplots(figsize=(5,4.5))
        ns=d[d['cat']=='ns']; up=d[d['cat']=='up']; dn=d[d['cat']=='down']
        ax.scatter(ns['atac_l2fc'],ns['rna_l2fc'],c=COL_NS,s=6,alpha=0.4,edgecolor='none',rasterized=True)
        ax.scatter(dn['atac_l2fc'],dn['rna_l2fc'],c=COL_DOWN,s=12,alpha=0.85,edgecolor='none',label='ATAC(dec) & RNA (sig.)')
        ax.scatter(up['atac_l2fc'],up['rna_l2fc'],c=COL_UP,s=12,alpha=0.85,edgecolor='none',label='ATAC(inc) & RNA (sig.)')
        # regression = 유의 점만
        m,b_=np.polyfit(sig_pts['atac_l2fc'],sig_pts['rna_l2fc'],1)
        xr=np.array([d['atac_l2fc'].min(),d['atac_l2fc'].max()])
        ax.plot(xr,m*xr+b_,'k-',linewidth=0.8)
        ax.axhline(0,ls='--',color='grey',linewidth=0.5); ax.axvline(0,ls='--',color='grey',linewidth=0.5)
        # Q1(++)/Q3(--) 라벨 각 5개
        d['score']=d['atac_l2fc'].abs()+d['rna_l2fc'].abs()
        q1=d[(d['atac_l2fc']>0)&(d['rna_l2fc']>0)&(d['cat']!='ns')].nlargest(5,'score')
        q3=d[(d['atac_l2fc']<0)&(d['rna_l2fc']<0)&(d['cat']!='ns')].nlargest(5,'score')
        texts=[ax.text(r['atac_l2fc'],r['rna_l2fc'],r['gene'],fontsize=7.5) for _,r in pd.concat([q1,q3]).iterrows()]
        try: adjust_text(texts,arrowprops=dict(arrowstyle='-',color='grey',lw=0.4),ax=ax)
        except Exception as exc:
            warnings.warn(f"Label placement failed: {exc}")
        ax.text(0.02,0.98,f'r = {r_val:.3f}\np = {p_val:.2e}',transform=ax.transAxes,ha='left',va='top',fontsize=9)
        ax.set_title(title,fontsize=11,fontweight='bold')
        ax.set_xlabel('ATAC log2 FC (0h vs 3h)',fontsize=10); ax.set_ylabel('RNA log2 FC',fontsize=10)
        ax.legend(loc='lower right',fontsize=6,frameon=False)
        for s in ['top','right']: ax.spines[s].set_visible(False)
        plt.tight_layout(); plt.savefig(out_pdf,dpi=300); plt.close()
        print(f'{title}: r={r_val:.3f}, p={p_val:.2e}, n_sig={len(sig_pts)}')

    make_scatter(tss_gene,'TSS',args.output_dir / 'output_tss_scatter.pdf')
    make_scatter(enh_gene,'Enhancer',args.output_dir / 'output_enhancer_scatter.pdf')


    import matplotlib.pyplot as plt
    from scipy import stats
    plt.rcParams['pdf.fonttype']=42; plt.rcParams['font.family']='DejaVu Sans'

    tss_sig = tss_gene[(tss_gene['atac_padj']<ATAC_PADJ)&(tss_gene['atac_l2fc'].abs()>ATAC_L2FC)]
    enh_sig = enh_gene[(enh_gene['atac_padj']<ATAC_PADJ)&(enh_gene['atac_l2fc'].abs()>ATAC_L2FC)]

    def sig_label(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'n.s.'

    def compare_direction_groups(direction, out_pdf):
        if direction=='up':
            tss_g=set(tss_sig[tss_sig['atac_l2fc']>0]['gene']); enh_g=set(enh_sig[enh_sig['atac_l2fc']>0]['gene'])
            rna_dir=rna[(rna['rna_padj'].notna())&(rna['rna_padj']<RNA_PADJ)&(rna['rna_l2fc']>RNA_L2FC)]
            fill='#B2182B'
        else:
            tss_g=set(tss_sig[tss_sig['atac_l2fc']<0]['gene']); enh_g=set(enh_sig[enh_sig['atac_l2fc']<0]['gene'])
            rna_dir=rna[(rna['rna_padj'].notna())&(rna['rna_padj']<RNA_PADJ)&(rna['rna_l2fc']<-RNA_L2FC)]
            fill='#2166AC'
        rset=set(rna_dir['SYMBOL'])
        only_tss=(tss_g-enh_g)&rset; only_enh=(enh_g-tss_g)&rset; both=(tss_g&enh_g)&rset
        rl=rna.set_index('SYMBOL')['rna_l2fc'].to_dict()
        cats=['TSS only','Enh only','TSS+Enh']
        data=[[rl[g] for g in sorted(s) if g in rl] for s in [only_tss,only_enh,both]]
        counts=[len(x) for x in data]
        for category, genes in zip(cats, [only_tss, only_enh, both]):
            group_rows.extend({'direction': direction, 'group': category, 'gene': g, 'rna_l2fc': rl[g]}
                              for g in sorted(genes) if g in rl)
            count_rows.append({'direction': direction, 'group': category, 'n': len(genes)})
        if any(len(arr) < 2 or np.ptp(arr) == 0 for arr in data):
            raise ValueError(f'{direction}: violin plots need at least two nonconstant values per group')

        fig,ax=plt.subplots(figsize=(4.5,4.2)); pos=np.arange(3)
        parts=ax.violinplot(data,positions=pos,widths=0.7,showmeans=False,showmedians=False,showextrema=False)
        for pc in parts['bodies']: pc.set_facecolor(fill); pc.set_alpha(0.8); pc.set_edgecolor('#4a4a4a'); pc.set_linewidth(0.8)
        for i,arr in enumerate(data):
            if arr:
                j=np.random.default_rng(42+i).uniform(-0.12,0.12,len(arr))
                ax.scatter(np.full(len(arr),i)+j,arr,s=4,c='#333',alpha=0.3,edgecolor='none',rasterized=True)
                ax.plot([i-0.2,i+0.2],[np.median(arr)]*2,'w-',linewidth=2)
        # pairwise Mann-Whitney
        allv=np.concatenate([a for a in data if a]); ymax=allv.max(); ymin=allv.min(); yr=ymax-ymin
        hs=[ymax+yr*0.05,ymax+yr*0.13,ymax+yr*0.21] if direction=='up' else [ymin-yr*0.09,ymin-yr*0.17,ymin-yr*0.25]
        for (i,j),h in zip([(0,1),(1,2),(0,2)],hs):
            if len(data[i])>=2 and len(data[j])>=2:
                statistic,p=stats.mannwhitneyu(data[i],data[j],alternative='two-sided')
                test_rows.append({'direction': direction, 'group_1': cats[i], 'group_2': cats[j],
                                  'n_1': len(data[i]), 'n_2': len(data[j]), 'U': statistic, 'pvalue': p})
                tick=yr*0.02*(1 if direction=='up' else -1)
                ax.plot([i,i,j,j],[h,h+tick,h+tick,h],'k-',linewidth=0.6)
                ax.text((i+j)/2,h+tick,sig_label(p),ha='center',va='bottom' if direction=='up' else 'top',fontsize=9)
        ax.set_xticks(pos); ax.set_xticklabels([f'{c}\n(n={n})' for c,n in zip(cats,counts)],fontsize=8)
        ax.set_ylabel('Log2 FC (RNA)',fontsize=10)
        for s in ['top','right']: ax.spines[s].set_visible(False)
        plt.tight_layout(); plt.savefig(out_pdf,dpi=300); plt.close()
        print(f'{direction}: n={counts}')

    compare_direction_groups('up',args.output_dir / 'output_groups_up.pdf')
    compare_direction_groups('down',args.output_dir / 'output_groups_down.pdf')


    pd.DataFrame(correlation_rows).to_csv(args.output_dir / 'output_correlations.csv', index=False)
    pd.DataFrame(group_rows).to_csv(args.output_dir / 'output_group_data.csv', index=False)
    pd.DataFrame(count_rows).to_csv(args.output_dir / 'output_group_counts.csv', index=False)
    pd.DataFrame(test_rows).to_csv(args.output_dir / 'output_group_tests.csv', index=False)


if __name__ == "__main__":
    main()
