# -*- coding: utf-8 -*-
import flask
import math
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objects as go
import numpy
from dash.dependencies import Input, Output

from timecourse import RenderableTC
from heatmap_gene import RenderableHeatmap
from mara import RenderableMARA

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

#foo = pd.read_csv('unprocessed/mTORCpaper/processed_raw.csv', sep="\t", index_col="Ensembl Gene ID")
#print(foo.head)
#exit(0)

### Read heatmap, ThExpress RNAseq
heatmap_thexpress = RenderableHeatmap(
    pd.read_csv('processed/ThExpress/meanNormalisedCounts.txt', sep="\t", index_col="Ensembl Gene ID"),
    map_ensembl_genesym_mouse,
    "Stubbington et al 2015 (RNA-seq)",
    ["naive", "Th1", "Th2", "Th17", "iTreg", "nTreg"]
)

heatmap_mtor = RenderableHeatmap(
    pd.read_csv('processed/mTORCpaper/processed_raw.csv', sep="\t", index_col="Ensembl Gene ID"),
    map_ensembl_genesym_mouse,
    "Howden et al 2019 (proteomics MS)",
    ["CD4 naive","CD4 tcr","CD4 Th1","CD8 tcr","CD8 naive"]
)





### Prepare suggestions for gene names
gene_suggestions = [] 
#map_ensembl_genesym['Associated Gene Name']
#!! TODO problem. if giving all the gene symbols the browser crashes. should search and only provide an approximate list. for now, ignore problem



### Prepare suggestions for motif names
motif_suggestions = pd.read_feather("processed/mara/henrikssonHumanTF/motif_ranking.feather")["motif"].to_list()


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

### MARA over tc
tc_mara_mouse = RenderableMARA(
    pd.read_feather("processed/mara/henrikssonMouseTF/activity.feather"),
    pd.read_feather("processed/mara/henrikssonMouseTF/cellcond.feather"),
    "MARA mouse"
)

tc_mara_human = RenderableMARA(
    pd.read_feather("processed/mara/henrikssonHumanTF/activity.feather"),
    pd.read_feather("processed/mara/henrikssonHumanTF/cellcond.feather"),
    "MARA human"
)


##################################################################################################################
############# Main page ##########################################################################################
##################################################################################################################
print("====== starting server =========\n");



prefix="/"    #on heroku
#prefix="/tcell/"
server = flask.Flask(__name__)
app = dash.Dash(__name__, server=server, routes_pathname_prefix=prefix)
app.title = "T cell data visualizer"

#Decides the layout of web app with titles and labels
app.layout = html.Div(children=[
    html.H1(children='T cell data visualizer'),
    html.Div(children='''
        Per-gene data: Search for example the gene Gata3
    '''),
    html.Datalist(id='list-suggested-genes',  children=[html.Option(value=word) for word in gene_suggestions]),
    html.Datalist(id='list-suggested-motifs', children=[html.Option(value=word) for word in motif_suggestions]),
    #############################
    dcc.Input(id='input-gene-name',
        type='text',
        list='list-suggested-genes',
        value='',
        debounce=True
    ),
    #############################
    html.Div(
      children=[
          html.Div(id="tc-human-out", style={'width': '49%', 'display': 'inline-block'}),
          html.Div(id="tc-mouse-out", style={'width': '49%', 'display': 'inline-block'}),
      ]
    ),

    #############################
    html.Div(
        children=[
            dcc.Graph(id="heatmap_the",  style={"width": "49%", 'display': 'inline-block'}),
            dcc.Graph(id="heatmap_mtor", style={"width": "49%", 'display': 'inline-block'}),
        ]
    ),


    ############################# MARA
    html.Div(children='''
        Per-motif data: Search for example the motif RBM6.1
    '''),
    dcc.Input(id='input-motif-name',
        type='text',
        list='list-suggested-motifs',
        value='',
        debounce=True
    ),
    html.Div(
        children=[
            html.Div(id="tc-mara-mouse-out",  style={"width": "49%", 'display': 'inline-block'}),
            html.Div(id="tc-mara-human-out",  style={"width": "49%", 'display': 'inline-block'}),
        ]
    ),


    #dcc.Tabs(id="tabs", value='tab-1', children=[
    #    dcc.Tab(label='Tab one', value='tab-1'),
    #    dcc.Tab(label='Tab two', value='tab-2'),
    #]),

    ])







##################################################################################################################
############# Time course RNA-seq visualizer #####################################################################
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


##################################################################################################################
############# Heatmaps ###########################################################################################
##################################################################################################################
@app.callback(
    Output("heatmap_the", "figure"),
    [Input("input-gene-name", "value")]
)
def update_heatmap_the(name_of_gene):
    return heatmap_thexpress.render_heatmap(name_of_gene)

@app.callback(
    Output("heatmap_mtor", "figure"),
    [Input("input-gene-name", "value")]
)
def update_heatmap_mtor(name_of_gene):
    return heatmap_mtor.render_heatmap(name_of_gene)



##################################################################################################################
############# Time course MARA visualizer ########################################################################
##################################################################################################################
@app.callback(
    Output("tc-mara-mouse-out", "children"),
    [Input("input-motif-name", "value")]
)
def render_mara(name_of_gene):
    return tc_mara_mouse.render_tc(name_of_gene)

@app.callback(
    Output("tc-mara-human-out", "children"),
    [Input("input-motif-name", "value")]
)
def render_mara(name_of_gene):
    return tc_mara_human.render_tc(name_of_gene)



##################################################################################################################
############# Tabs content #######################################################################################
##################################################################################################################
#@app.callback(
#    Output('tabs', 'children'),
#    [Input('tabs', 'value')])
def render_content(tab):
    if tab == 'tab-1':
        return html.Div([
            html.H3('Tab content 1')
        ])
    elif tab == 'tab-2':
        return html.Div([
            html.H3('Tab content 2')
        ])

# To do later, good for book viewer

##################################################################################################################
##################################################################################################################
##################################################################################################################
if __name__ == '__main__':
    app.run_server(debug=True)




@app.route('/wtf')
def index():
    return "wtf"
