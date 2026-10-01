from pathlib import Path
import pandas as pd
from sklearn.metrics import r2_score,mean_squared_error,mean_absolute_error
from scipy.stats import pearsonr
ROOT=Path(__file__).resolve().parents[1]
d=pd.read_csv(ROOT/"01_data/validation/chirps_imerg_monthly.csv").dropna(subset=["chirps","imerg"]); OUT=ROOT/"03_results"; OUT.mkdir(exist_ok=True)
o=d.chirps.to_numpy(); p=d.imerg.to_numpy(); r=pearsonr(o,p)
res={"n":len(d),"R2":r2_score(o,p),"pearson_r":r.statistic,"pearson_p":r.pvalue,"RMSE":mean_squared_error(o,p)**0.5,"MAE":mean_absolute_error(o,p),"bias_imerg_minus_chirps":(p-o).mean()}
pd.DataFrame([res]).to_csv(OUT/"CHIRPS_IMERG_Validation.csv",index=False); print(res)
