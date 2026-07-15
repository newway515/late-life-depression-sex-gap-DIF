# -*- coding: utf-8 -*-
"""Exploratory longitudinal CHARLS nonparametric sex-DIF and gap analysis.
Checks stability of (a) Delta_raw, (b) which items show substantial sex DIF (ETS B/C),
(c) DIF-purified gap. Data: cesd_analysis.db / charls_cesd_items_long.
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
D=[i+"_d" for i in ITEMS]

con=sqlite3.connect(str(DB)); df=pd.read_sql("SELECT * FROM charls_cesd_items_long",con); con.close()
df["ragender"]=pd.to_numeric(df.ragender,errors="coerce")
df["raeducl"]=pd.to_numeric(df.raeducl,errors="coerce")
for c in D+["cesd","agey"]: df[c]=pd.to_numeric(df[c],errors="coerce")

def cohend(f,m):
    sp=math.sqrt(((len(f)-1)*f.var(ddof=1)+(len(m)-1)*m.var(ddof=1))/(len(f)+len(m)-2))
    return (f.mean()-m.mean())/sp

def strata(score,k=10):
    q=np.unique(np.quantile(score,np.linspace(0,1,k+1)))
    return np.clip(np.digitize(score,q[1:-1]),0,len(q)-2)

def mh_delta(item_bin, focal, score, k=10):
    s=strata(score,k); numOR=denOR=A=EA=VA=0.0
    for st in np.unique(s):
        idx=s==st; g=focal[idx]; y=item_bin[idx]; Nk=idx.sum()
        if Nk<2: continue
        nf=(g==1).sum(); nr=(g==0).sum(); m1=(y==1).sum(); m0=(y==0).sum()
        if nf==0 or nr==0 or m1==0 or m0==0: continue
        a_=((g==1)&(y==1)).sum(); b_=((g==1)&(y==0)).sum()
        c_=((g==0)&(y==1)).sum(); d_=((g==0)&(y==0)).sum()
        numOR+=a_*d_/Nk; denOR+=b_*c_/Nk
        A+=a_; EA+=nf*m1/Nk; VA+=nf*nr*m1*m0/(Nk*Nk*(Nk-1))
    if denOR==0 or VA==0: return np.nan,np.nan
    OR=numOR/denOR; delta=-2.35*math.log(OR)
    chi2=(abs(A-EA)-0.5)**2/VA; p=math.erfc(math.sqrt(chi2/2))
    return delta,p
def ets(d):
    ad=abs(d); return "A" if ad<1 else ("B" if ad<1.5 else "C")

waves=[1,2,3,4]
summ=[]; delmat={}
for w in waves:
    dw=df[df.wave==w]
    a=dw[(dw.agey>=60)&dw[D].notna().all(axis=1)&dw.ragender.isin([1,2])&dw.raeducl.notna()].copy()
    a["female"]=(a.ragender==2).astype(int); tot=a["cesd"].values
    draw=cohend(a[a.female==1].cesd,a[a.female==0].cesd)
    deltas={}; flagged=[]
    for it in ITEMS:
        xb=(a[it+"_d"].values>=1).astype(int)
        d_,p=mh_delta(xb,a.female.values,tot)
        deltas[it]=d_
        if ets(d_) in ("B","C"): flagged.append(it)
    delmat[w]=deltas
    keep=[i+"_d" for i in ITEMS if i not in flagged]
    a["tp"]=a[keep].sum(axis=1)
    dpure=cohend(a[a.female==1].tp,a[a.female==0].tp)
    summ.append(dict(wave=w,year={1:2011,2:2013,3:2015,4:2018}[w],N=len(a),
                     femalePct=round(a.female.mean(),3),Draw=round(draw,3),
                     Dpure=round(dpure,3),pct_measurement=round(100*(draw-dpure)/draw,1),
                     flagged=";".join(flagged) or "none",sleepr_delta=round(deltas["sleepr"],3)))
S=pd.DataFrame(summ)
pd.set_option("display.width",200)
print("=== CHARLS longitudinal analysis (age>=60) ===")
print(S.to_string(index=False))
print("\n=== sex ETS delta per item across waves (|>1|=B substantial) ===")
DM=pd.DataFrame(delmat).round(2); DM.index.name="item"; DM.columns=[f"w{w}" for w in waves]
print(DM.to_string())

S.to_csv(f"{OUT}/charls_longitudinal_summary.csv",index=False)
DM.to_csv(f"{OUT}/charls_longitudinal_sexDIF_delta.csv")

# figure: ETS delta trajectories per item
fig,ax=plt.subplots(figsize=(8,5))
for it in ITEMS:
    ys=[delmat[w][it] for w in waves]
    lw=2.5 if it=="sleepr" else 1.0
    ax.plot([f"w{w}" for w in waves],ys,marker="o",lw=lw,label=it)
for y in (-1,1): ax.axhline(y,ls="--",lw=.7,color="k",alpha=.4)
ax.axhspan(-1,1,color="grey",alpha=.08)
ax.set_ylabel("sex ETS delta (|>1| = substantial B)"); ax.set_title("CHARLS: sex DIF per item across waves (bold=sleepr)")
ax.legend(fontsize=7,ncol=2,loc="upper right")
plt.tight_layout(); plt.savefig(f"{OUT}/charls_longitudinal_sexDIF.png",dpi=130); plt.close()
print("\nSaved: charls_longitudinal_summary.csv, _sexDIF_delta.csv, _sexDIF.png")
