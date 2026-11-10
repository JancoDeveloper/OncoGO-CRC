import json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,pandas as pd
O='/workspace/tesis-oncogo-entrega/';plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
J=json.load(open('results/tesis_cierre/cierre_audit.json'))['vias_suero']['contrastes']
lab={'crc_vs_sano':'Cáncer vs. sanos','crc_vs_enfermedad':'Cáncer vs. col. negativa','enfermedad_vs_sano':'Col. negativa vs. sanos'}
df=pd.DataFrame([[lab[k],v['n_terminos'],v['pearson_nes'],v['p_nulo'],v['n_permutaciones']] for k,v in J.items()],columns=['contraste_suero','terminos_GO_evaluables','pearson_NES_tejido_suero','p_nulo','permutaciones']);df.round(3).to_csv(O+'Tabla3-4_vias_tejido_en_suero.csv',index=False)
f,ax=plt.subplots(figsize=(6,3.6));b=ax.bar(df.contraste_suero,df.pearson_NES_tejido_suero,color='#7f8c8d')
for r,(n,p) in zip(b,zip(df.terminos_GO_evaluables,df.p_nulo)): ax.text(r.get_x()+r.get_width()/2,r.get_height()-0.03,f'{n} términos\np nulo={p:.2f}'.replace('.',','),ha='center',va='top',fontsize=8)
ax.axhline(0,c='k',lw=.6);ax.set_ylim(-.5,.2);ax.set_ylabel('Pearson de NES (CPTAC vs. suero)')
ax.set_title('Términos GO alterados en tejido CPTAC (FDR≤0,25)\ny su puntuación en el suero PXD055218',fontsize=9);f.tight_layout();f.savefig(O+'Fig3-8_vias_tejido_en_suero.png');plt.close()
# pipeline
f,ax=plt.subplots(figsize=(11,4.4));ax.axis('off')
def box(x,y,t,w=2.3,h=.8,c='#eaf2f8'):
    ax.add_patch(plt.Rectangle((x,y),w,h,fc=c,ec='#2c3e50',lw=.8));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=7.5)
def arr(x1,y1,x2,y2): ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',lw=.8))
rows=[(3.2,'PDC000116 CPTAC\nTMT, tejido pareado'),(1.9,'PXD034575 Tanaka\nLFQ, P / M / N (16 pac.)'),(0.6,'PXD055218 Rai\nDIA suero, 20/20/20')]
for y,t in rows:
    box(0,y,t);box(2.8,y,'log2 de valor observado\ncero = ausencia\nsin imputación');arr(2.3,y+.4,2.8,y+.4)
    box(5.6,y,'t moderada + treat\n|log2FC|>0,585\nBH FDR<0,05');arr(5.1,y+.4,5.6,y+.4)
box(8.4,2.5,'Comparación de\nsigno y vías\n(cohortes NO fusionadas)',h=1.2,c='#fdedec');arr(7.9,3.6,8.4,3.3);arr(7.9,2.3,8.4,2.8)
box(8.4,0.3,'AUC + bootstrap 2000\nregla de especificidad\nmarcas de ruido previas',h=1.2,c='#fdedec');arr(7.9,1.0,8.4,0.9)
ax.set_xlim(-.1,10.9);ax.set_ylim(0,4.2);f.tight_layout();f.savefig(O+'Fig2-1_flujo_de_analisis.png');plt.close()
