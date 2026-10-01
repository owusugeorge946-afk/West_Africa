from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
x=pd.read_csv(ROOT/"03_results/DCI_ERPI_Annual.csv"); OUT=ROOT/"03_results"
w={"rainfall_deficit":0.261396,"cdd":0.323791,"lst":0.159832}; sw=sum(w.values()); w={k:v/sw for k,v in w.items()}
x["dci_no_ndvi"]=sum(w[k]*x[k] for k in w)
x["rank_primary"]=x.dci_entropy.rank(ascending=False,method="min"); x["rank_no_ndvi"]=x.dci_no_ndvi.rank(ascending=False,method="min")
x["rank_shift"]=(x.rank_primary-x.rank_no_ndvi).abs(); rho,p=spearmanr(x.dci_entropy,x.dci_no_ndvi)
x.to_csv(OUT/"DCI_NDVI_Sensitivity.csv",index=False)
pd.DataFrame([{"spearman_rho":rho,"p_value":p,"max_absolute_rank_shift":x.rank_shift.max()}]).to_csv(OUT/"DCI_NDVI_Sensitivity_Summary.csv",index=False)
print("rho=",rho,"p=",p,"max rank shift=",x.rank_shift.max())
