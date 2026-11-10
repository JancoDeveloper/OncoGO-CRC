import pandas as pd,numpy as np,json,matplotlib
matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
O='/workspace/tesis-oncogo-entrega/';T='results/tesis_cierre/'
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
fi=pd.read_csv(T+'fichas_proteinas.tsv',sep='\t');u2g=dict(zip(fi.uniprot,fi.gene))
m=pd.read_csv('plasma_raw.tsv',sep='\t').drop_duplicates('protein_id').set_index('protein_id')
m=np.log2(m.where(m>0))
for _u in ['O94812','P02751','Q9BY43']:
    _h=[i for i in m.index if i.split(';')[0].split('-')[0]==_u]
    if _h: m.loc[_u]=m.loc[_h[0]]
grp={'CP':'Cáncer','CN':'Colonoscopia negativa','Normal':'Sanos'}
cols={k:[c for c in m.columns if c.startswith(k+'_')] for k in grp}
LFC=0.585
def volcano(ax,df,title,lab):
    df=df.dropna(subset=['log2fc','padj']).copy();y=-np.log10(df.padj.clip(lower=1e-12))
    s=df.significant.astype(str)=='True'
    ax.scatter(df.log2fc[~s],y[~s],s=6,c='#b0b0b0',lw=0);ax.scatter(df.log2fc[s],y[s],s=14,c='#c0392b',lw=0)
    ax.axhline(-np.log10(.05),ls='--',c='k',lw=.6);[ax.axvline(v,ls=':',c='k',lw=.6) for v in(-LFC,LFC)]
    for _,r in df[s].iterrows():
        if lab: ax.annotate(lab(r.protein_id),(r.log2fc,-np.log10(max(r.padj,1e-12))),fontsize=7,xytext=(3,2),textcoords='offset points')
    ax.set_title(f'{title}\n{s.sum()} de {len(df)} con FDR<0,05 y |log2FC|>0,585',fontsize=9);ax.set_xlabel('log2FC');ax.set_ylabel('−log10(FDR)')
# F1 serum volcanoes
f,ax=plt.subplots(1,3,figsize=(13,4.2),sharey=True)
for a,(fn,t) in zip(ax,[('crc_vs_enfermedad','Cáncer vs. colonoscopia negativa'),('crc_vs_sano','Cáncer vs. sanos'),('enfermedad_vs_sano','Colonoscopia negativa vs. sanos')]):
    volcano(a,pd.read_csv(T+fn+'.tsv',sep='\t'),t,lambda p:u2g.get(p.split(';')[0].split('-')[0],p.split(';')[0]))
f.tight_layout();f.savefig(O+'Fig3-1_volcanes_suero_PXD055218.png');plt.close()
# F2 boxplots
genes=[('SAA1','P0DJI8'),('LBP','P18428'),('CHMP4A','Q9BY43'),('APOA4','P06727')]
f,ax=plt.subplots(1,4,figsize=(13,3.8));rows=[]
for a,(g,u) in zip(ax,genes):
    v=[m.loc[u,cols[k]].dropna().values if u in m.index else [] for k in grp]
    a.boxplot(v,widths=.5,showfliers=False)
    for i,x in enumerate(v): a.scatter(np.random.default_rng(1).normal(i+1,.06,len(x)),x,s=8,c='#34495e',alpha=.7)
    a.set_xticks([1,2,3],[f'{n}\n(n={len(x)})' for n,x in zip(['Cáncer','Col. neg.','Sanos'],v)],fontsize=8);a.set_title(f'{g} ({u})',fontsize=9);a.set_ylabel('log2 intensidad')
    rows.append([g,u]+[len(x) for x in v]+[round(float(np.median(x)),2) if len(x) else '' for x in v])
f.tight_layout();f.savefig(O+'Fig3-2_cajas_SAA1_LBP_CHMP4A_APOA4.png');plt.close()
pd.DataFrame(rows,columns=['gen','uniprot','n_cancer','n_colneg','n_sanos','mediana_cancer','mediana_colneg','mediana_sanos']).to_csv(O+'Tabla3-2_n_y_medianas_cajas.csv',index=False)
# F3 ROC 8 candidates CRC vs col neg
es=pd.read_csv(T+'especificidad_candidatas.tsv',sep='\t');e8=es[es.veredicto=='especifica'].sort_values('auc_crc_vs_enfermedad',ascending=False)
f,ax=plt.subplots(figsize=(5.6,5.2));chk=[]
for _,r in e8.iterrows():
    u=fi.set_index('gene').uniprot.get(r.gene)
    if u not in m.index: chk.append([r.gene,u,'sin_matriz','',r.auc_crc_vs_enfermedad]);continue
    a=m.loc[u,cols['CP']].dropna().values;b=m.loc[u,cols['CN']].dropna().values
    auc=mannwhitneyu(a,b).statistic/(len(a)*len(b));sgn=1 if auc>=.5 else -1
    sc=np.r_[a,b]*sgn;y=np.r_[np.ones(len(a)),np.zeros(len(b))];th=np.unique(sc)[::-1]
    tpr=[0]+[((sc>=t)&(y==1)).sum()/len(a) for t in th];fpr=[0]+[((sc>=t)&(y==0)).sum()/len(b) for t in th]
    ls='--' if str(r.lista_ruido)=='fase_aguda' else '-'
    ax.plot(fpr,tpr,ls,lw=1.3,label=f'{r.gene} AUC {max(auc,1-auc):.2f} ({len(a)}/{len(b)})'.replace('.',','))
    chk.append([r.gene,u,len(a),len(b),round(max(auc,1-auc),4),round(r.auc_crc_vs_enfermedad,4)])
ax.plot([0,1],[0,1],c='grey',lw=.6);ax.set_xlabel('1 − especificidad (colonoscopia negativa)');ax.set_ylabel('Sensibilidad (cáncer)')
ax.legend(fontsize=7,frameon=False,loc='lower right');ax.set_title('Curvas ROC por proteína, cáncer vs. colonoscopia negativa\n(línea discontinua: marca de fase aguda)',fontsize=9)
f.tight_layout();f.savefig(O+'Fig3-3_ROC_ocho_candidatas.png');plt.close()
pd.DataFrame(chk,columns=['gen','uniprot','n_cancer','n_colneg','auc_recalculado','auc_tesis_cierre']).to_csv(O+'control_auc_recalculado_vs_tesis_cierre.csv',index=False)
# F4 forest AUC
ev=es.dropna(subset=['auc_crc_vs_enfermedad']).sort_values('auc_crc_vs_enfermedad')
f,ax=plt.subplots(figsize=(6.5,6.5));yy=np.arange(len(ev))
c=['#c0392b' if v=='especifica' else '#7f8c8d' for v in ev.veredicto]
ax.errorbar(ev.auc_crc_vs_enfermedad,yy,xerr=[(ev.auc_crc_vs_enfermedad-ev.auc_crc_vs_enfermedad_ci_bajo).clip(lower=0),(ev.auc_crc_vs_enfermedad_ci_alto-ev.auc_crc_vs_enfermedad).clip(lower=0)],fmt='none',ecolor='#bbb')
ax.scatter(ev.auc_crc_vs_enfermedad,yy,c=c,s=25,zorder=3,label='Cáncer vs. col. negativa');ax.scatter(ev.auc_enfermedad_vs_sano,yy,marker='x',c='k',s=20,label='Col. negativa vs. sanos')
ax.axvline(.7,ls='--',lw=.6,c='k');ax.set_yticks(yy,[f'{g}{" *" if isinstance(n,str) else ""}' for g,n in zip(ev.gene,ev.lista_ruido)],fontsize=8)
ax.set_xlabel('AUC (IC 95 % bootstrap)');ax.legend(fontsize=7,frameon=False,loc='lower right');ax.set_title('Regla de especificidad de las candidatas de literatura\n(rojo: cumple; * marca previa de ruido)',fontsize=9)
f.tight_layout();f.savefig(O+'Fig3-4_AUC_candidatas_literatura.png');plt.close()
# T ocho candidatas
t=e8.merge(fi[['gene','uniprot','padj_crc_vs_enfermedad','log2fc_crc_vs_enfermedad']],on='gene',how='left')
t=t[['gene','uniprot','lista_ruido','auc_crc_vs_enfermedad','auc_crc_vs_enfermedad_ci_bajo','auc_crc_vs_enfermedad_ci_alto','auc_enfermedad_vs_sano','log2fc_crc_vs_enfermedad','padj_crc_vs_enfermedad','lectura']].round(3)
t.to_csv(O+'Tabla3-3_ocho_candidatas_AUC.csv',index=False)
cl=pd.read_csv(T+'candidatas_literatura.tsv',sep='\t');cl.to_csv(O+'TablaA_catalogo_37_candidatas_literatura.csv',index=False)
# auditorias heatmap
V={'replicada':3,'no_replicada':1,'no_evaluable':0,'significativa_en_sentido_opuesto':2,'direccion_no_unica_en_la_fuente':0.5};aud={}
for k,l in [('auditoria_cptac','CPTAC tejido'),('auditoria_crc_vs_sano','Suero cáncer vs sanos'),('auditoria_crc_vs_enfermedad','Suero cáncer vs col. neg.'),('auditoria_enfermedad_vs_sano','Suero col. neg. vs sanos')]:
    a=pd.read_csv(T+k+'.tsv',sep='\t');aud[l]=a.set_index('gene').veredicto
A=pd.DataFrame(aud);A.to_csv(O+'Tabla3-5_auditoria_catalogo_por_cohorte.csv')
f,ax=plt.subplots(figsize=(5,9));from matplotlib.colors import ListedColormap
ax.imshow(A.apply(lambda c:c.map(V)).values.astype(float),cmap=ListedColormap(['#eeeeee','#d6dbdf','#f5b7b1','#2e86c1','#c0392b']),vmin=0,vmax=3,aspect='auto')
ax.set_yticks(range(len(A)),A.index,fontsize=7);ax.set_xticks(range(4),A.columns,rotation=35,ha='right',fontsize=8)
ax.set_title('Auditoría del catálogo de literatura\ngris: no evaluable / dirección no única · rosa: no replica\nazul: significativa en sentido opuesto · rojo: replica',fontsize=9);f.tight_layout();f.savefig(O+'Fig3-7_auditoria_catalogo_heatmap.png');plt.close()
# Tanaka
J=json.load(open(T+'cierre_audit.json'));tk={}
f,ax=plt.subplots(1,3,figsize=(13,4.2),sharey=True)
for a,(fn,tt) in zip(ax,[('primario_vs_mucosa','Primario (P) vs. mucosa normal (N)'),('metastasis_vs_mucosa','Metástasis hepática (M) vs. mucosa (N)'),('metastasis_vs_primario','Metástasis (M) vs. primario (P)')]):
    d=pd.read_csv(T+'tanaka_'+fn+'.tsv',sep='\t');tk[fn]=d.set_index('protein_id').log2fc;volcano(a,d,tt,None)
    d[d.significant.astype(str)=='True'].sort_values('padj')[['protein_id','log2fc','padj','n_obs']].round(4).to_csv(O+f'TablaA_tanaka_{fn}_significativas.csv',index=False)
f.tight_layout();f.savefig(O+'Fig3-5_volcanes_tejido_Tanaka_PXD034575.png');plt.close()
x=pd.concat([tk['primario_vs_mucosa'],tk['metastasis_vs_mucosa']],axis=1,keys=['p','m']).dropna()
f,ax=plt.subplots(figsize=(5,5));ax.scatter(x.m,x.p,s=5,c='#34495e',alpha=.5,lw=0);lim=[x.min().min(),x.max().max()];ax.plot(lim,lim,c='grey',lw=.6)
ax.axhline(0,c='k',lw=.4);ax.axvline(0,c='k',lw=.4);ax.set_xlabel('log2FC metástasis vs. mucosa');ax.set_ylabel('log2FC primario vs. mucosa')
ax.set_title(f'Mismo paciente, dos tejidos tumorales\n{len(x)} proteínas con ambos efectos',fontsize=9);f.tight_layout();f.savefig(O+'Fig3-6_tanaka_primario_vs_metastasis.png');plt.close()
# Tables of counts
C=J['contrastes'];pd.DataFrame([[k,v['n_a'],v['n_b'],v['n_probadas'],v['n_significativas']] for k,v in C.items()],columns=['contraste','n_a','n_b','probadas','significativas']).to_csv(O+'Tabla3-1_conteos_suero.csv',index=False)
pd.read_csv(T+'peptidos_neoantigenos.tsv',sep='\t').to_csv(O+'Tabla3-6_peptidos_neoantigenos.csv',index=False)
fi.iloc[:,:16].round(4).to_csv(O+'TablaA_fichas_21_proteinas.csv',index=False)
print(json.dumps({k:v for k,v in J.items() if 'via' in k or 'gsea' in k.lower()},ensure_ascii=False)[:1500]);print(list(J.keys()))
