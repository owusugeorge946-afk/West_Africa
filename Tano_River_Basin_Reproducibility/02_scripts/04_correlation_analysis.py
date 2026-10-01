from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
df=pd.read_csv(ROOT/"01_data/Tano_Final_Annual_Hydroclimatic_Results.csv")
OUT=ROOT/"03_results"; OUT.mkdir(exist_ok=True)
cols=["rainfall","rx1day","rx5day","cdd","ndvi","lst"]
rho=df[cols].corr(method="spearman"); pval=pd.DataFrame(index=cols,columns=cols,dtype=float)
for a in cols:
 for b in cols: pval.loc[a,b]=spearmanr(df[a],df[b],nan_policy="omit").pvalue
rho.to_csv(OUT/"Spearman_Correlation.csv"); pval.to_csv(OUT/"Spearman_Pvalues.csv"); print(rho)
