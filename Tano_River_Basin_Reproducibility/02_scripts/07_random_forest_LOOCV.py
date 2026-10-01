from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score,mean_squared_error,mean_absolute_error
ROOT=Path(__file__).resolve().parents[1]
d=pd.read_csv(ROOT/"03_results/DCI_ERPI_Annual.csv"); OUT=ROOT/"03_results"
specs={"DCI":(["rainfall_deficit","cdd","ndvi_deficit","lst"],"dci_entropy"),"ERPI":(["rainfall","rx1day","rx5day"],"erpi_entropy")}
pred_rows=[]; imp_rows=[]; summary=[]
for name,(features,target) in specs.items():
 X=d[features].to_numpy(); y=d[target].to_numpy(); pred=np.zeros(len(y))
 for tr,te in LeaveOneOut().split(X):
  m=RandomForestRegressor(n_estimators=300,random_state=42,max_features=1.0,min_samples_leaf=1,bootstrap=True,n_jobs=1); m.fit(X[tr],y[tr]); pred[te]=m.predict(X[te])
 full=RandomForestRegressor(n_estimators=300,random_state=42,max_features=1.0,min_samples_leaf=1,bootstrap=True,n_jobs=1).fit(X,y)
 for f,v in zip(features,full.feature_importances_): imp_rows.append({"index":name,"feature":f,"importance":v})
 for yr,o,pr in zip(d.year,y,pred): pred_rows.append({"index":name,"year":yr,"observed":o,"predicted":pr})
 summary.append({"index":name,"r2":r2_score(y,pred),"rmse":mean_squared_error(y,pred)**0.5,"mae":mean_absolute_error(y,pred)})
pd.DataFrame(imp_rows).to_csv(OUT/"RF_Feature_Importance.csv",index=False); pd.DataFrame(pred_rows).to_csv(OUT/"RF_LOOCV_Predictions.csv",index=False); pd.DataFrame(summary).to_csv(OUT/"RF_LOOCV_Summary.csv",index=False); print(pd.DataFrame(summary))
