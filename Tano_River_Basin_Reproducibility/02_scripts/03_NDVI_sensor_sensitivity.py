from pathlib import Path
import pandas as pd
import statsmodels.api as sm
import importlib.util
ROOT=Path(__file__).resolve().parents[1]
df=pd.read_csv(ROOT/"01_data/Tano_Final_Annual_Hydroclimatic_Results.csv")
OUT=ROOT/"03_results"; OUT.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location("mkmod",ROOT/"02_scripts/02_modified_MK_SenSlope.py")
mkmod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mkmod)
rows=[]
for label,a,b in [("pre_OLI_2000_2012",2000,2012),("OLI_era_2013_2024",2013,2024),("pooled_2000_2024",2000,2024)]:
 d=df[df.year.between(a,b)]; q=mkmod.hr_test(d.year,d.ndvi); rows.append({"period":label,"n":len(d),**q})
pd.DataFrame(rows).to_csv(OUT/"NDVI_Sensor_Era_Sensitivity.csv",index=False)
t=df.year-df.year.min(); post=(df.year>=2013).astype(int); post_t=(df.year-2013).clip(lower=0)
X=sm.add_constant(pd.DataFrame({"time":t,"post2013":post,"post2013_time":post_t}))
fit=sm.OLS(df.ndvi,X).fit()
pd.DataFrame({"term":fit.params.index,"estimate":fit.params.values,"se":fit.bse.values,"p":fit.pvalues.values}).to_csv(OUT/"NDVI_2013_Segmented_Regression.csv",index=False)
print(pd.DataFrame(rows)[["period","sen_slope","modified_Z","modified_p"]])
