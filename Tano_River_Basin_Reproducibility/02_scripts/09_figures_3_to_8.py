from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"01_data/Tano_Final_Annual_Hydroclimatic_Results.csv"; R=ROOT/"03_results"; F=ROOT/"04_figures"; F.mkdir(exist_ok=True)
df=pd.read_csv(D)
# Figure 3: six annual series. Statistical annotations must come from corrected Table 3.
t=pd.read_csv(R/"Table_3_Modified_MK.csv").set_index("variable")
fig,axs=plt.subplots(3,2,figsize=(12,11),constrained_layout=True)
vars=["rainfall","rx1day","rx5day","cdd","ndvi","lst"]
labels=["Annual rainfall (mm)","Rx1day (mm)","Rx5day (mm)","CDD (days)","NDVI","LST (°C)"]
for ax,v,lab in zip(axs.flat,vars,labels):
 ax.plot(df.year,df[v],marker="o",linewidth=1)
 ax.set_ylabel(lab); ax.set_xlabel("Year"); ax.grid(alpha=.2)
 q=t.loc[v]; ax.text(.02,.96,f"Sen slope={q.sen_slope:.4g}; modified MK p={q.modified_p:.3g}",transform=ax.transAxes,va="top")
fig.savefig(F/"Figure_3.png",dpi=300,bbox_inches="tight"); plt.close(fig)
# Figure 4: hydro-climate temporal variability.
fig,ax=plt.subplots(figsize=(11,5))
for v in ["rainfall","rx1day","rx5day"]: ax.plot(df.year,df[v],marker="o",label=v)
ax.set_xlabel("Year"); ax.set_ylabel("Precipitation metric (mm)"); ax.legend(); ax.grid(alpha=.2)
fig.savefig(F/"Figure_4a.png",dpi=300,bbox_inches="tight"); plt.close(fig)
fig,ax=plt.subplots(figsize=(11,4)); ax.plot(df.year,df.cdd,marker="o"); ax.set(xlabel="Year",ylabel="CDD (days)"); ax.grid(alpha=.2)
fig.savefig(F/"Figure_4b.png",dpi=300,bbox_inches="tight"); plt.close(fig)
# Figure 5: NDVI-LST interaction.
fig,ax=plt.subplots(figsize=(11,5)); ax2=ax.twinx()
ax.plot(df.year,df.ndvi,marker="o",label="NDVI"); ax2.plot(df.year,df.lst,marker="s",linestyle="--",label="LST")
ax.axvline(2013,linestyle=":",linewidth=1); ax.set(xlabel="Year",ylabel="NDVI"); ax2.set_ylabel("LST (°C)")
fig.savefig(F/"Figure_5.png",dpi=300,bbox_inches="tight"); plt.close(fig)
print("Figures 3-5 generated. Figures 6-8 should be generated only after their source result CSVs are validated.")
