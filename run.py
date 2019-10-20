# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objects as go
import numpy
from dash.dependencies import Input, Output

from timecourse import RenderableTC

#Opens and reads file from url -- could be useful on heroku to avoid huge uploads through GIT
#df = pd.read_csv('https://raw.githubusercontent.com/RasmusBurge/DBT-T-Cell/master/rawdata.csv', sep=";")


##################################################################################################################
############# Prepare data #######################################################################################
##################################################################################################################
print("====== Reading data =============\n")

### Read gene symbol mapping tables. Create a lower case version of symbol for fast search later
map_ensembl_genesym_mouse = pd.read_csv('processed/ensembl_mouse.csv', sep=",")

### Read gene symbol mapping tables. Create a lower case version of symbol for fast search later
map_ensembl_genesym_human = pd.read_csv('processed/ensembl_human.csv', sep=",")



heatmap_data = pd.read_csv('processed/ThExpress/meanNormalisedCounts.txt', sep=",")


### Prepare suggestions for gene names
gene_suggestions = [] 
#map_ensembl_genesym['Associated Gene Name']
#!! TODO problem. if giving all the gene symbols the browser crashes. should search and only provide an approximate list. for now, ignore problem







### Mouse time-course RNA-seq
tc_mouse = RenderableTC(
    pd.read_csv('processed/th2crispr/th2crispr_mouse_tcavg_data.csv', sep=",", index_col="id"),
    pd.read_csv('processed/th2crispr/th2crispr_mouse_tcavg_samplemeta.csv', sep=","),
    map_ensembl_genesym_mouse,
    "Mouse gene expression (RNA)"
)

### Human time-course RNA-seq
tc_human = RenderableTC(
    pd.read_csv('processed/th2crispr/th2crispr_human_tcavg_data.csv', sep=",", index_col="id"),
    pd.read_csv('processed/th2crispr/th2crispr_human_tcavg_samplemeta.csv', sep=","),
    map_ensembl_genesym_human,
    "Human gene expression (RNA)"
)


##################################################################################################################
############# Main page ##########################################################################################
##################################################################################################################
print("====== starting server =========\n");

app = dash.Dash()
app.title = "T cell data visualizer"
server = app.server

#Decides the layout of web app with titles and labels
app.layout = html.Div(children=[
    html.H1(children='T cell data visualizer'),
    html.Div(children='''
        Search for example the gene Gata3
    '''),
    html.Datalist(id='list-suggested-genes', children=[html.Option(value=word) for word in gene_suggestions]),
    dcc.Input(id='input-gene-name',
        type='text',
        list='list-suggested-genes',
        value='',
        debounce=True
    ),
    html.Div(
      children=[
          html.Div(id="tc-human-out", style={'width': '49%', 'display': 'inline-block'}),
          html.Div(id="tc-mouse-out", style={'width': '49%', 'display': 'inline-block'}),

      ]
    ),

    dcc.Tabs(id="tabs", value='tab-1', children=[
        dcc.Tab(label='Tab one', value='tab-1'),
        dcc.Tab(label='Tab two', value='tab-2'),
    ]),

    ])




##################################################################################################################
############# Time course visualizer #############################################################################
##################################################################################################################
@app.callback(
    Output("tc-mouse-out", "children"),
    [Input("input-gene-name", "value")]
)
def number_render(name_of_gene):
    return tc_mouse.render_tc(name_of_gene)

@app.callback(
    Output("tc-human-out", "children"),
    [Input("input-gene-name", "value")]
)
def number_render(name_of_gene):
    return tc_human.render_tc(name_of_gene)




################ To do later, good for book viewer
@app.callback(
    Output('tabs', 'children'),  #tabs-content
    [Input('tabs', 'value')])
def render_content(tab):
    if tab == 'tab-1':
        return html.Div([
            html.H3('Tab content 1')
        ])
    elif tab == 'tab-2':
        return html.Div([
            html.H3('Tab content 2')
        ])



if __name__ == '__main__':
    app.run_server(debug=True)


#        raise PreventUpdate


