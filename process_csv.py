import pandas as pd

df = pd.read_csv("max_LD_results.csv")

df = df[df["status"]=="ok"]
df = df[df["CL_max"]>1.5]
df = df[df["CL"]<0.4]
df = df[df["t_over_c"]>=0.1]
df_sorted = df.sort_values(by="L_over_D", ascending=False)

print(df_sorted.head(20))
