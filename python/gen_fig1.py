import argparse
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams["font.family"]="DejaVu Sans"
fig,ax=plt.subplots(figsize=(8.2,10.6)); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")

def box(cx,cy,w,h,title,body,fc="#EAF1F8",ec="#2B5A86",tsz=10,bsz=8.2):
    ax.add_patch(FancyBboxPatch((cx-w/2,cy-h/2),w,h,boxstyle="round,pad=0.4,rounding_size=1.6",fc=fc,ec=ec,lw=1.3))
    if title and body:
        ax.text(cx,cy+h*0.20,title,ha="center",va="center",fontsize=tsz,fontweight="bold",color="#173754")
        ax.text(cx,cy-h*0.20,body,ha="center",va="center",fontsize=bsz,color="#222",linespacing=1.3)
    else:
        ax.text(cx,cy,title,ha="center",va="center",fontsize=tsz,fontweight="bold",color="#173754")

def arrow(x1,y1,x2,y2,color="#2B5A86"):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=15,lw=1.5,color=color,shrinkA=1,shrinkB=1))

# Row A: three cohorts (title+body inside)
box(20,94,30,9,"CHARLS · China","wave 4 (2018)\nN = 8,143",fc="#FBEFE6",ec="#B5651D",tsz=10,bsz=8)
box(50,94,30,9,"ELSA · England","wave 6\nN = 6,522",fc="#FBEFE6",ec="#B5651D",tsz=10,bsz=8)
box(80,94,30,9,"HRS · United States","wave 4\nN = 12,848",fc="#FBEFE6",ec="#B5651D",tsz=10,bsz=8)

box(50,81,86,8.5,"Harmonized CES-D · adults aged ≥ 60 · complete cases",
    "Pooled analytic sample N = 27,513  (CHARLS-10 graded; ELSA/HRS-8 binary; 6 common items)",fc="#EFEFEF",ec="#555",tsz=10,bsz=8)
arrow(20,89.5,45,85.4); arrow(50,89.5,50,85.4); arrow(80,89.5,55,85.4)

box(50,70,80,9,"Measurement model (GRM)",
    "Single-factor scoring; bifactor confirms a dominant general factor (ECV 0.84, ωH 0.89)",tsz=10,bsz=8)
arrow(50,76.7,50,74.6)
box(50,57.5,80,9.5,"Sex DIF detection",
    "IRT χ² (BH-adjusted) gated by effect size ESSD ≥ 0.20\nSubstantive DIF: sleepr, fear, flone",tsz=10,bsz=8)
arrow(50,65.2,50,62.4)
box(50,44,80,10,"Sex gap under three estimands",
    "Δ_raw 0.346  →  Δ_latent 0.319 [0.266, 0.372]  →  Δ_adj 0.295\nH2 supported: gap shrinks but stays positive (one-sided p < 0.0001)",fc="#E7F0E7",ec="#2E6B34",tsz=10,bsz=8)
arrow(50,51.5,50,49.1)

# validation branches
arrow(50,38.9,50,35.6)
ax.add_patch(FancyArrowPatch((14,34),(86,34),arrowstyle="-",lw=1.3,color="#2B5A86"))
for x in (14,38,62,86): arrow(x,34,x,31.4)
box(14,25.5,20,11,"Longitudinal","w1–w4 stable;\nsleepr biased 4/4",tsz=9.2,bsz=7.6)
box(38,25.5,20,11,"Cross-cohort","alignment R² > 0.90;\nChina largest",tsz=9.2,bsz=7.6)
box(62,25.5,20,11,"External validity","CIDI robust;\ngrip → more aligned",tsz=9.2,bsz=7.6)
box(86,25.5,20,11,"Sensitivity","9 settings all +;\nrange 0.081 ≤ 0.10",tsz=9.2,bsz=7.6)
for x in (14,38,62,86): arrow(x,19.9,50,14.6)

box(50,8.5,86,10,"Conclusion",
    "The female–male late-life depression gap is real and robust to measurement non-invariance,\nwith limited measurement amplification in somatic, fear, and loneliness items.",
    fc="#F3ECF7",ec="#5B3A82",tsz=10.5,bsz=8.2)

parser=argparse.ArgumentParser(description="Generate the study-flow figure")
parser.add_argument("--output",type=Path,default=Path("figures/Figure1_study_flow.png"))
args=parser.parse_args(); args.output.parent.mkdir(parents=True,exist_ok=True)
plt.savefig(args.output,dpi=170,bbox_inches="tight"); print(f"saved {args.output}")
