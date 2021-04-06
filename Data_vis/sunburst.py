# Sunburst plot of PhD data
import plotly.express as px
import pandas as pd

df = pd.read_csv("~/data_for_phd.csv")
df['value'] = 1

fig = px.sunburst(df, path=['Data', 'Diagnosis', 'CN', 'Exome', "RNA", "WGS"], values='value', color='Diagnosis')
fig.show()

fig.write_image("thesis_data.svg")
