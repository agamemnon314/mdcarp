# %%
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
# %%
df = pd.read_csv('solution_kpi.csv', index_col=False)
# %%
fig = go.Figure()
fig.add_trace(go.Scatter(y=df['approx'],
                         mode='markers',
                         name='approx upper bound (from paper)'))
fig.add_trace(go.Scatter(y=df['approx_improved_1'],
                         mode='markers',
                         name='improved approx upper bound 1'))
fig.add_trace(go.Scatter(y=df['approx_improved_2'],
                         mode='markers',
                         name='improved approx upper bound 2'))
fig.add_trace(go.Scatter(y=df['ub_approx'],
                         mode='lines',
                         name='theoretical approx upper bound 2'))
fig.update_layout(title='All Sample',
yaxis_title='Approximation rate')
fig.show()
fig.write_image('all_data.svg')
# %%
df_mean = df[['n_depots', 'n_required_edges', 'approx', 'approx_improved_1', 'approx_improved_2',
              'ub_approx']].groupby(by=['n_depots']).mean()
# %%
fig = go.Figure()
fig.add_trace(go.Bar(x=df_mean.index.astype('str'), y=df_mean['approx'],
                     name='approx upper bound (from paper)'))
fig.add_trace(go.Bar(x=df_mean.index.astype('str'), y=df_mean['approx_improved_1'],
                     name='improved approx upper bound 1'))
fig.add_trace(go.Bar(x=df_mean.index.astype('str'), y=df_mean['approx_improved_2'],
                     name='improved approx upper bound 2'))
fig.add_trace(go.Bar(x=df_mean.index.astype('str'), y=df_mean['ub_approx'],
                     name='theoretical approx upper bound 2'))
fig.show()
# %%
ddf = df[['n_depots', 'n_required_edges', 'approx', 'approx_improved_1', 'approx_improved_2',
          'ub_approx']].groupby(by=['n_depots', 'n_required_edges']).mean()
# %%
for i in df_mean.index:
    df_i = df[df.n_depots==i]
    fig = go.Figure()
    fig.add_trace(go.Scatter(y=df_i['approx'],
                            mode='markers',
                            name='approx upper bound (from paper)'))
    fig.add_trace(go.Scatter(y=df_i['approx_improved_1'],
                            mode='markers',
                            name='improved approx upper bound 1'))
    fig.add_trace(go.Scatter(y=df_i['approx_improved_2'],
                            mode='markers',
                            name='improved approx upper bound 2'))
    fig.add_trace(go.Scatter(y=df_i['ub_approx'],
                            mode='lines',
                            name='theoretical approx upper bound 2'))
    fig.update_layout(title='number of depot = '+str(i),
    yaxis_title='Approximation rate')
    fig.show()
    fig.write_image('sub_data_'+str(i)+'.svg')
# %%

# %%
