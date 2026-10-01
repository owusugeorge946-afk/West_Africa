"""Hamed-Rao variance-corrected Mann-Kendall test + Sen slope.
Autocorrelation correction is applied to ranks after detrending by Sen slope.
Only statistically significant rank ACF lags (two-sided 95% normal threshold)
enter the variance inflation factor. Classical S/tau are retained for transparency.
"""
from pathlib import Path
import math, numpy as np, pandas as pd
from scipy.stats import norm, theilslopes, rankdata

ROOT=Path(__file__).resolve().parents[1]
df=pd.read_csv(ROOT/"01_data/Tano_Final_Annual_Hydroclimatic_Results.csv")
OUT=ROOT/"03_results"; OUT.mkdir(exist_ok=True)

def mk_s(x):
    x=np.asarray(x,float); n=len(x)
    return int(sum(np.sign(x[j]-x[i]) for i in range(n-1) for j in range(i+1,n)))

def tie_var(x):
    n=len(x); _,cnt=np.unique(x,return_counts=True)
    return (n*(n-1)*(2*n+5)-np.sum(cnt*(cnt-1)*(2*cnt+5)))/18.0

def hr_test(year,x):
    year=np.asarray(year,float); x=np.asarray(x,float); n=len(x)
    S=mk_s(x); var0=tie_var(x)
    tau=S/(0.5*n*(n-1))
    slope,intercept,_,_=theilslopes(x,year,0.95)
    detr=x-(intercept+slope*year)
    r=rankdata(detr)
    rc=r-r.mean(); den=np.dot(rc,rc)
    crit=norm.ppf(.975)/math.sqrt(n)
    acc=0.0; used=[]
    for k in range(1,n):
        rho=float(np.dot(rc[:-k],rc[k:])/den) if den else 0.0
        if abs(rho)>crit:
            acc+=(n-k)*(n-k-1)*(n-k-2)*rho
            used.append((k,rho))
    nns=1+(2/(n*(n-1)*(n-2)))*acc if n>2 else 1
    # Guard against pathological negative finite-sample factors.
    nns=max(nns,1e-12)
    var=var0*nns
    z=(S-1)/math.sqrt(var) if S>0 else ((S+1)/math.sqrt(var) if S<0 else 0.0)
    p=2*(1-norm.cdf(abs(z)))
    return dict(S=S,tau=tau,sen_slope=slope,classical_variance=var0,
                variance_correction=nns,modified_variance=var,
                modified_Z=z,modified_p=p,
                significant_acf_lags=";".join(f"{k}:{v:.6f}" for k,v in used))

rows=[]
for v in ["rainfall","rx1day","rx5day","cdd","ndvi","lst"]:
    q=hr_test(df.year,df[v])
    q["variable"]=v; rows.append(q)
res=pd.DataFrame(rows)
res.to_csv(OUT/"Table_3_Modified_MK.csv",index=False)
print(res.to_string(index=False))
