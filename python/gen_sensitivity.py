import argparse
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
settings=["Full-invariance baseline","Wave 3 (2015)","Main (data-driven anchors)",
          "Leave-one-out (drop sleepr)","w4 recompute (age>=60)","Wave 1 (2011)",
          "Affective-only anchors","Wave 2 (2013)","Common 6 items (binary)"]
vals=[0.367,0.366,0.319,0.319,0.319,0.317,0.306,0.302,0.286]
order=np.argsort(vals); settings=[settings[i] for i in order]; vals=[vals[i] for i in order]
fig,ax=plt.subplots(figsize=(7.2,4.6))
cmap=plt.cm.Reds; norm=plt.Normalize(0.0,0.40)
for i,v in enumerate(vals):
    ax.barh(i,v,color=cmap(norm(v)),edgecolor="k",lw=.5)
    ax.text(v+0.004,i,f"{v:.3f}",va="center",fontsize=9)
ax.axvline(0.319,ls="--",lw=.8,color="#444"); ax.text(0.319,len(vals)-0.4,"  main = 0.319",fontsize=8,color="#444")
ax.set_yticks(range(len(settings))); ax.set_yticklabels(settings,fontsize=9)
ax.set_xlim(0,0.42); ax.set_xlabel("Latent-mean sex gap, female − male (d units)")
ax.set_title("Robustness of the female–male depression gap across analytic settings\n(CHARLS; all directions positive, cross-setting range = 0.081 ≤ 0.10)",fontsize=10)
parser=argparse.ArgumentParser(description="Generate the sensitivity-matrix figure")
parser.add_argument("--output",type=Path,default=Path("figures/Figure5_sensitivity_matrix.png"))
args=parser.parse_args(); args.output.parent.mkdir(parents=True,exist_ok=True)
plt.tight_layout(); plt.savefig(args.output,dpi=150); print(f"saved {args.output}")
