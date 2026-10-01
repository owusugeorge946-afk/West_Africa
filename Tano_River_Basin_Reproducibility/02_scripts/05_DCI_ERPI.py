from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
df=pd.read_csv(ROOT/"01_data/Tano_Final_Annual_Hydroclimatic_Results.csv")
OUT=ROOT/"03_results"; OUT.mkdir(exist_ok=True)
def mm(s):
 d=s.max()-s.min()
 if d==0: raise ValueError("Zero range in "+s.name)
 return (s-s.min())/d
wd={"rainfall_deficit":0.261396,"cdd":0.323791,"lst":0.159832,"ndvi_deficit":0.254981}
we={"rainfall":0.302832,"rx1day":0.398009,"rx5day":0.299158}
z=pd.DataFrame({"year":df.year})
z["rainfall_deficit"]=1-mm(df.rainfall); z["cdd"]=mm(df.cdd); z["lst"]=mm(df.lst); z["ndvi_deficit"]=1-mm(df.ndvi)
z["dci_entropy"]=sum(wd[k]*z[k] for k in wd)
z["rainfall"]=mm(df.rainfall); z["rx1day"]=mm(df.rx1day); z["rx5day"]=mm(df.rx5day)
z["erpi_entropy"]=sum(we[k]*z[k] for k in we)
z.to_csv(OUT/"DCI_ERPI_Annual.csv",index=False); pd.DataFrame([wd,we],index=["DCI","ERPI"]).to_csv(OUT/"DCI_ERPI_Weights.csv")
print(z[["year","dci_entropy","erpi_entropy"]].to_string(index=False))
