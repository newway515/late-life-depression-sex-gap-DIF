# -*- coding: utf-8 -*-
"""三库(CHARLS/ELSA/HRS)描述统计 Table 1 + Δ_raw(bootstrap CI) + 非参数性别 DIF(MH ETS Δ)。
主分析波: 各库 >=60 岁完整个案最多之波。沙箱 numpy/pandas 可发表层结果。
"""
import sqlite3, math, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rng=np.random.default_rng(20260709)
DB="/sessions/keen-magical-pasteur/mnt/SQLitedatabase/cesd_analysis.db"; OUT="/sessions/keen-magical-pasteur/mnt/ccProj01"

SPECS={
 "CHARLS":("charls_cesd_items_long","ID",
    ["depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear"],"raeducl"),
 "ELSA":("elsa_cesd_items_long","idauniq",
    ["depres","effort","sleepr","whappy","flone","fsad","going","enlife"],"raeducl"),
 "HRS":("hrs_cesd_items_long","hhidpn",
    ["depres","effort","sleepr","whappy","flone","fsad","going","enlife"],"raeduc"),
}
COMMON6=["depres","effort","sleepr","whappy","flone","going"]

def cohend(f,m):
    s=math.sqrt(((len(f)-1)*np.var(f,ddof=1)+(len(m)-1)*np.var(m,ddof=1))/(len(f)+len(m)-2))
    return (np.mean(f)-np.mean(m))/s
def strata(x,k=10):
    q=np.unique(np.quantile(x,np.linspace(0,1,k+1))); return np.clip(np.digitize(x,q[1:-1]),0,len(q)-2)
def mh_delta(item_bin,focal,score,k=10):
    s=strata(score,k); numOR=denOR=A=EA=VA=0.0
    for st in np.unique(s):
        idx=s==st; g=focal[idx]; y=item_bin[idx]; Nk=idx.sum()
        if Nk<2: continue
        nf=(g==1).sum(); nr=(g==0).sum(); m1=(y==1).sum(); m0=(y==0).sum()
        if nf==0 or nr==0 or m1==0 or m0==0: continue
        a=((g==1)&(y==1)).sum(); b=((g==1)&(y==0)).sum(); c=((g==0)&(y==1)).sum(); d=((g==0)&(y==0)).sum()
        numOR+=a*d/Nk; denOR+=b*c/Nk; A+=a; EA+=nf*m1/Nk; VA+=nf*nr*m1*m0/(Nk*Nk*(Nk-1))
    if denOR==0 or VA==0: return np.nan,np.nan
    OR=numOR/denOR; delta=-2.35*math.log(OR); chi2=(abs(A-EA)-0.5)**2/VA
    return delta, math.erfc(math.sqrt(chi2/2))
def ets(d): ad=abs(d); return "A" if ad<1 else ("B" if ad<1.5 else "C")

con=sqlite3.connect(DB)
table1=[]; difrows=[]; common_delta={}
for coh,(tab,idc,items,educ) in SPECS.items():
    d=pd.read_sql(f"SELECT * FROM {tab}",con)
    D=[i+"_d" for i in items]
    for c in D+["cesd","agey","ragender",educ]: d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d[(d.agey>=60)&d[D].notna().all(axis=1)&d.ragender.isin([1,2])&d[educ].notna()]
    w=d.groupby("wave").size().idxmax(); a=d[d.wave==w].copy()
    a["female"]=(a.ragender==2).astype(int); tot=a["cesd"].values
    f=a.cesd[a.female==1].values; m=a.cesd[a.female==0].values
    draw=cohend(f,m)
    bs=[cohend(rng.choice(f,len(f)),rng.choice(m,len(m))) for _ in range(2000)]
    lo,hi=np.percentile(bs,[2.5,97.5])
    table1.append(dict(cohort=coh,wave=int(w),N=len(a),female_pct=round(a.female.mean(),3),
        age_mean=round(a.agey.mean(),1), cesd_F=round(f.mean(),2), cesd_M=round(m.mean(),2),
        Draw=round(draw,3), Draw_CI=f"[{lo:.3f}, {hi:.3f}]"))
    # DIF per item
    dl={}
    for it in items:
        xb=(a[it+"_d"].values>=1).astype(int)
        de,p=mh_delta(xb,a.female.values,tot)
        difrows.append(dict(cohort=coh,item=it,common6=it in COMMON6,
            sex_ETSdelta=round(de,3),ets=ets(de),p=round(p,4)))
        if it in COMMON6: dl[it]=de
    common_delta[coh]=dl
con.close()

T1=pd.DataFrame(table1); DIF=pd.DataFrame(difrows)
pd.set_option("display.width",200,"display.max_columns",30)
print("===== Table 1 (主分析波, 年龄>=60, 完整个案) =====")
print(T1.to_string(index=False))
print("\n===== 非参数性别 DIF (MH ETS Δ; 负=女性同等抑郁下更易报告; |>1|=B) =====")
print(DIF.to_string(index=False))
print("\n实质性(B/C)条目:")
print(DIF[DIF.ets.isin(["B","C"])][["cohort","item","sex_ETSdelta","ets"]].to_string(index=False))
T1.to_csv(f"{OUT}/table1_threecohort.csv",index=False)
DIF.to_csv(f"{OUT}/dif_threecohort.csv",index=False)

# figure: common-6 sex ETS delta by cohort
cd=pd.DataFrame(common_delta).reindex(COMMON6)
fig,ax=plt.subplots(figsize=(8,5)); x=np.arange(len(COMMON6)); w=0.26
for i,coh in enumerate(["CHARLS","ELSA","HRS"]):
    ax.bar(x+(i-1)*w, cd[coh].values, w, label=coh)
for y in (-1,1): ax.axhline(y,ls="--",lw=.7,color="k",alpha=.4)
ax.axhspan(-1,1,color="grey",alpha=.08)
ax.set_xticks(x); ax.set_xticklabels(COMMON6,rotation=30); ax.set_ylabel("Sex ETS delta")
ax.set_title("Common 6-item sex DIF across cohorts (|delta|>1 = substantial)"); ax.legend()
plt.tight_layout(); plt.savefig(f"{OUT}/dif_common6_bycohort.png",dpi=130); plt.close()
print("\nSaved: table1_threecohort.csv, dif_threecohort.csv, dif_common6_bycohort.png")
