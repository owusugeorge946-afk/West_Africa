from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
src=ROOT/"01_data/Tano_Final_Annual_Hydroclimatic_Results.csv"
out=ROOT/"03_results"; out.mkdir(exist_ok=True)
df=pd.read_csv(src)
cols=["rainfall","rx1day","rx5day","cdd","ndvi","lst"]
missing=[c for c in ["year",*cols] if c not in df]
if missing: raise ValueError(f"Missing columns: {missing}")
if df.year.duplicated().any(): raise ValueError("Duplicate years")
s=df[cols].describe(percentiles=[.25,.5,.75]).T
s=s.rename(columns={"25%":"Q1","50%":"median","75%":"Q3"})
s[["count","mean","std","min","Q1","median","Q3","max"]].to_csv(out/"Table_2_Descriptive_Statistics.csv")
print(s[["count","mean","std","min","Q1","median","Q3","max"]])
