"""Draw Figure 1 (primary age tests: donor-level scores and adjusted effect estimates) from the two delivered source files.
Inputs (same folder as this script, read only, never rewritten):
  figure1_donor_level_source_values.tsv  - one row per modelled donor: panel (A,B,C), test_id, donor, age_years, score_fixed_reference_unadjusted
  figure1_summary_source_values.tsv       - saved primary covariate-adjusted estimates (P01,P02,P03)
Outputs: Figure1_age_donor_level.png / .pdf (written to the folder given as first argument, default: this folder) and Figure1_plot_QA.json.
Presentation only: no model is fitted, no trend line, correlation, rescaling, resampling or new P value is computed.
Usage: python draw_figure1_donor_level.py [output_dir]
Requires: matplotlib (see requirements.lock.txt). The upstream analysis pipeline is NOT included in this package."""
from pathlib import Path
import csv,json,hashlib,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
F=Path(__file__).resolve().parent
OUTD=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else F
OUTD.mkdir(parents=True,exist_ok=True)
rd=lambda p:list(csv.DictReader(open(p),delimiter='\t'))
out=rd(F/'figure1_donor_level_source_values.tsv')
summ={r['test_id']:r for r in rd(F/'figure1_summary_source_values.tsv')}
prim=summ
assert set(prim)=={'P01','P02','P03'}
assert [sum(1 for r in out if r['panel']==p) for p in 'ABC']==[43,29,31]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
blue='#28627d';dark='#243747'
fig=plt.figure(figsize=(11.0,9.0))
gs=fig.add_gridspec(1,3,left=.075,right=.975,top=.80,bottom=.465,wspace=.30)
gs2=fig.add_gridspec(1,1,left=.27,right=.975,top=.345,bottom=.085)
titles={'A':'Normal bulk hippocampus (GSE48350)\ncore score, n=43, ages 20–99',
        'B':'SEA-AD CA3 neurons\ncKO-defined program, n=29, ages 69–99',
        'C':'SEA-AD subicular neurons\ncKO-defined program, n=31, ages 69–100'}
xl={'A':(15,104),'B':(64,104),'C':(64,104)}
ylim=[]
for i,pan in enumerate('ABC'):
    ax=fig.add_subplot(gs[0,i])
    rows=[r for r in out if r['panel']==pan]
    x=[float(r['age_years']) for r in rows];y=[float(r['score_fixed_reference_unadjusted']) for r in rows]
    ax.scatter(x,y,s=34,color=blue,alpha=.72,edgecolor='white',linewidth=.6)
    ax.axhline(0,color='#b4bcc2',lw=.8);ax.set_xlim(*xl[pan]);ax.set_title(titles[pan],fontsize=9.6,loc='left',pad=6)
    ax.set_xlabel('Recorded age (years)');ax.grid(alpha=.15)
    if i==0:ax.set_ylabel('Fixed-reference score (SD units)')
    ax.text(-.02,1.27,pan,transform=ax.transAxes,fontsize=14,fontweight='bold',va='bottom',ha='right')
    ylim+= [min(y),max(y)]
lo,hi=min(ylim),max(ylim);pad=.08*(hi-lo)
for ax in fig.axes:ax.set_ylim(lo-pad,hi+pad)
ax=fig.add_subplot(gs2[0,0])
labels=['Normal bulk hippocampus\ncore score (n=43)','CA3 neurons\ncKO-defined program (n=29)','Subicular neurons\ncKO-defined program (n=31)']
for k,t in enumerate(['P01','P02','P03']):
    r=prim[t];y=2-k;b,l,h=float(r['beta']),float(r['CI95_low']),float(r['CI95_high'])
    ax.plot([l,h],[y,y],color=blue,lw=2.2,solid_capstyle='butt');ax.plot([l,l],[y-.1,y+.1],color=blue,lw=1.6);ax.plot([h,h],[y-.1,y+.1],color=blue,lw=1.6)
    ax.scatter([b],[y],s=60,color=blue,zorder=3)
    ax.text(.27,y+.17,f"genes {int(float(r['A_n']))}/{int(float(r['G0_n']))};  unadjusted P={float(r['raw_P']):.4f}",fontsize=9.2,color=dark,va='center')
    ax.text(.27,y-.17,f"Holm-adjusted P={float(r['adjusted_P']):.4f}",fontsize=9.2,color=dark,va='center')
ax.axvline(0,color='#89939a',lw=1);ax.set_yticks([2,1,0],labels);ax.tick_params(axis='y',length=0)
ax.set_xlim(-1.2,1.05);ax.set_ylim(-.55,2.55);ax.set_xticks([-1,-.5,0,.5]);ax.grid(axis='x',alpha=.15)
ax.set_xlabel('Age coefficient from the primary covariate-adjusted model (fixed-reference score SD per decade)')
ax.text(-.30,1.16,'D',transform=ax.transAxes,fontsize=14,fontweight='bold',va='bottom',ha='right')
ax.set_title('Primary age family (Holm m=3): pointwise 95% HC3-t intervals',fontsize=10,loc='left',pad=6,color=dark)
fig.text(.075,.965,'Primary age tests: donor-level scores and adjusted effect estimates',fontsize=13.5,color='black')
fig.text(.075,.937,'A–C: unadjusted saved scores, one point per modelled donor, no fitted line.  D: saved adjusted estimates.  CA3 and subicular models share 24 donors.',fontsize=9,color=dark)
fig.savefig(OUTD/'Figure1_age_donor_level.png',dpi=300,facecolor='white');fig.savefig(OUTD/'Figure1_age_donor_level.pdf',facecolor='white');plt.close(fig)
n={p:sum(1 for r in out if r['panel']==p) for p in 'ABC'}
ages={p:(min(float(r['age_years']) for r in out if r['panel']==p),max(float(r['age_years']) for r in out if r['panel']==p)) for p in 'ABC'}
shared=len({r['donor'] for r in out if r['panel']=='B'}&{r['donor'] for r in out if r['panel']=='C'})
uniq={p:len({r['donor'] for r in out if r['panel']==p}) for p in 'ABC'}
assert n=={'A':43,'B':29,'C':31}==uniq and shared==24
qa={'status':'PASS','n_points':n,'unique_donors':uniq,'age_range':ages,'shared_B_C_donors':shared,
 'all_scores_ages_finite':all(abs(float(r['age_years']))<1e3 and abs(float(r['score_fixed_reference_unadjusted']))<1e3 for r in out),
 'inputs_read_only':['figure1_donor_level_source_values.tsv','figure1_summary_source_values.tsv'],
 'input_sha256':{n_:hashlib.sha256((F/n_).read_bytes()).hexdigest() for n_ in ['figure1_donor_level_source_values.tsv','figure1_summary_source_values.tsv']},
 'no_trend_line_correlation_refit_resampling_or_adjustment':True,'point_weighting':'none','outliers_removed':0,'ylim_shared_across_A_C':[lo-pad,hi+pad]}
(OUTD/'Figure1_plot_QA.json').write_text(json.dumps(qa,indent=1));print(json.dumps(qa)[:300])
