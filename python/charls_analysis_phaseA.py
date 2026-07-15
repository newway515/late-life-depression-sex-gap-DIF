# -*- coding: utf-8 -*-
"""Exploratory CHARLS descriptives and nonparametric DIF screen.
- pick main wave (age>=60, complete 10 items)
- descriptives + sex x education collinearity
- non-parametric DIF screen: polytomous stratified SMD + dichotomized Mantel-Haenszel (ETS delta)
  for SEX (female vs male) and EDUCATION (low vs high)
- raw gap (Cohen d) and 'DIF-purified' gap (drop flagged items) as first-pass Delta_latent direction
Outputs: CSVs + figures + printed summary.
"""
import argparse, os, sqlite3, math, numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--db", default=os.getenv("CESD_DB"), help="SQLite input database (or set CESD_DB)")
parser.add_argument("--out", default=os.getenv("CESD_OUT", "output"), help="Output directory")
args=parser.parse_args()
if not args.db: parser.error("Provide --db or set CESD_DB")
DB=Path(args.db); OUT=Path(args.out)
if not DB.is_file(): parser.error(f"Database not found: {DB}")
OUT.mkdir(parents=True, exist_ok=True)
ITEMS=["depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear"]
COMMON6=["depres","effort","sleepr","whappy","flone","going"]
D=[i+"_d" for i in ITEMS]

con=sqlite3.connect(str(DB))
df=pd.read_sql("SELECT * FROM charls_cesd_items_long", con); con.close()

# ---- pick main wave: age>=60 & all 10 items & gender & education ----
def analytic(dw):
    m=(dw.agey>=60)&dw[D].notna().all(axis=1)&dw.ragender.isin([1,2])&dw.raeducl.notna()
    return dw[m].copy()
cnt={w:len(analytic(df[df.wave==w])) for w in sorted(df.wave.unique())}
main=max(cnt,key=cnt.get)
print("age>=60 complete-case N by wave:",cnt,"-> main wave =",main)

a=analytic(df[df.wave==main]).reset_index(drop=True)
a["total"]=a[D].sum(axis=1)                       # 0-30 depressive direction
a["female"]=(a.ragender==2).astype(int)
# education dichotomy: low = raeducl==1 (lowest) vs high = >=2
a["edu_low"]=(a.raeducl<=1).astype(int)
N=len(a); print(f"\nMain wave {main}: analytic N={N}")

# ---- descriptives ----
def cohend(x,y):
    nx,ny=len(x),len(y); sp=math.sqrt(((nx-1)*x.var(ddof=1)+(ny-1)*y.var(ddof=1))/(nx+ny-2))
    return (x.mean()-y.mean())/sp
f_tot=a.loc[a.female==1,"total"]; m_tot=a.loc[a.female==0,"total"]
draw=cohend(f_tot,m_tot)
print(f"sex: female n={ (a.female==1).sum() }, male n={ (a.female==0).sum() }")
print(f"mean total F={f_tot.mean():.2f} M={m_tot.mean():.2f} | Delta_raw Cohen d(F-M)={draw:.3f}")
print("\nsex x education (row%):")
ct=pd.crosstab(a.female.map({1:'female',0:'male'}), a.raeducl, normalize='index').round(3)
print(ct)
print("edu_low share: female %.3f male %.3f"%(a.loc[a.female==1,'edu_low'].mean(),a.loc[a.female==0,'edu_low'].mean()))

# ---- helper: strata by total score deciles ----
def make_strata(score, k=10):
    # rank-based deciles, collapse to unique edges
    q=np.unique(np.quantile(score,np.linspace(0,1,k+1)))
    lab=np.clip(np.digitize(score,q[1:-1]),0,len(q)-2)
    return lab

# ---- polytomous stratified SMD (Zwick/Dorans STD P-DIF), focal vs ref ----
def poly_smd(item, focal, restscore):
    strata=make_strata(restscore)
    num=0.0; wsum=0.0
    for s in np.unique(strata):
        idx=strata==s
        fo=item[idx & (focal==1)]; re=item[idx & (focal==0)]
        if len(fo)>=1 and len(re)>=1:
            w=len(fo)                     # weight by focal N (standardization on focal)
            num+=w*(fo.mean()-re.mean()); wsum+=w
    smd=num/wsum if wsum>0 else np.nan
    return smd/item.std(ddof=1)           # standardized by item SD

# ---- dichotomized Mantel-Haenszel + ETS delta ----
def mh_dif(item_bin, focal, score, k=10):
    strata=make_strata(score,k)
    numOR=0.0; denOR=0.0; A=0.0; EA=0.0; VA=0.0
    for s in np.unique(strata):
        idx=strata==s
        g=focal[idx]; y=item_bin[idx]
        Nk=idx.sum()
        if Nk<2: continue
        nf=(g==1).sum(); nr=(g==0).sum()
        m1=(y==1).sum(); m0=(y==0).sum()
        if nf==0 or nr==0 or m1==0 or m0==0: continue
        a_=((g==1)&(y==1)).sum(); b_=((g==1)&(y==0)).sum()
        c_=((g==0)&(y==1)).sum(); d_=((g==0)&(y==0)).sum()
        numOR+=a_*d_/Nk; denOR+=b_*c_/Nk
        A+=a_; EA+=nf*m1/Nk; VA+=nf*nr*m1*m0/(Nk*Nk*(Nk-1))
    if denOR==0 or VA==0: return np.nan,np.nan,np.nan
    OR=numOR/denOR
    delta=-2.35*math.log(OR)              # ETS delta scale (ref=male,focal=female)
    chi2=(abs(A-EA)-0.5)**2/VA
    p=math.erfc(math.sqrt(chi2/2))        # chi2 df=1 -> p
    return OR,delta,p

def ets_class(delta):
    ad=abs(delta)
    return "A" if ad<1 else ("B" if ad<1.5 else "C")

# ---- run DIF for SEX and EDUCATION ----
rows=[]
tot=a["total"].values
for it in ITEMS:
    x=a[it+"_d"].astype(float).values
    rest=tot-x
    xb=(x>=1).astype(int)                 # endorsed (any symptom)
    # SEX: focal=female
    smd_s=poly_smd(pd.Series(x),a.female.values,pd.Series(rest))
    OR_s,delt_s,p_s=mh_dif(xb,a.female.values,tot)
    # EDU: focal=low education
    smd_e=poly_smd(pd.Series(x),a.edu_low.values,pd.Series(rest))
    OR_e,delt_e,p_e=mh_dif(xb,a.edu_low.values,tot)
    rows.append(dict(item=it,common6=it in COMMON6,
        sex_SMD=round(smd_s,3), sex_OR=round(OR_s,3), sex_ETSdelta=round(delt_s,3),
        sex_class=ets_class(delt_s), sex_p=round(p_s,4),
        edu_SMD=round(smd_e,3), edu_OR=round(OR_e,3), edu_ETSdelta=round(delt_e,3),
        edu_class=ets_class(delt_e), edu_p=round(p_e,4)))
dif=pd.DataFrame(rows)
pd.set_option("display.width",200,"display.max_columns",30)
print("\n=== DIF screen (ETS delta: >0 => item endorsed MORE by focal at same trait; focal=female / low-edu) ===")
print(dif.to_string(index=False))

# flag substantial sex DIF: ETS B/C
flagged=dif.loc[dif.sex_class.isin(["B","C"]),"item"].tolist()
print("\nItems with substantial SEX DIF (ETS B/C):", flagged)

# ---- Delta_raw vs DIF-purified gap ----
keep=[i+"_d" for i in ITEMS if i not in flagged]
a["total_pure"]=a[keep].sum(axis=1)
dpure=cohend(a.loc[a.female==1,"total_pure"],a.loc[a.female==0,"total_pure"])
print(f"\nDelta_raw (all 10 items) Cohen d = {draw:.3f}")
print(f"Delta_purified (drop {len(flagged)} DIF item(s)) Cohen d = {dpure:.3f}")
print(f"=> change {draw-dpure:+.3f} ({100*(draw-dpure)/draw:.1f}% of raw gap)")

# ---- save outputs ----
dif.to_csv(f"{OUT}/charls_dif_screen_w{main}.csv",index=False)
summ=pd.DataFrame({"metric":["main_wave","N","female_n","male_n","mean_total_F","mean_total_M",
                             "Delta_raw_d","Delta_purified_d","n_sexDIF_items"],
                   "value":[main,N,int((a.female==1).sum()),int((a.female==0).sum()),
                            round(f_tot.mean(),2),round(m_tot.mean(),2),round(draw,3),round(dpure,3),len(flagged)]})
summ.to_csv(f"{OUT}/charls_summary_w{main}.csv",index=False)

# ---- figures ----
# item endorsement (mean 0-3) by sex
im=a.groupby("female")[D].mean().T; im.index=ITEMS; im.columns=["male","female"]
ax=im.plot(kind="barh",figsize=(7,5),color=["#4C72B0","#C44E52"])
ax.set_xlabel("item mean (0-3, depressive direction)"); ax.set_title(f"CHARLS w{main} age>=60: item means by sex")
plt.tight_layout(); plt.savefig(f"{OUT}/charls_item_means_by_sex_w{main}.png",dpi=130); plt.close()
# sex DIF effect sizes
fig,ax=plt.subplots(figsize=(7,5))
colors=["#C44E52" if c in ("B","C") else "#AAAAAA" for c in dif.sex_class]
ax.barh(dif.item,dif.sex_ETSdelta,color=colors)
for x in (-1.5,-1,1,1.5): ax.axvline(x,ls="--",lw=.7,color="k",alpha=.4)
ax.set_xlabel("ETS delta (sex DIF; |>1| = B, |>1.5| = C; red=substantial)")
ax.set_title(f"CHARLS w{main}: sex DIF (MH ETS-delta)")
plt.tight_layout(); plt.savefig(f"{OUT}/charls_sex_dif_w{main}.png",dpi=130); plt.close()

print("\nSaved: charls_dif_screen / charls_summary CSVs + 2 figures to project folder.")
