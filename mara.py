import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas
import plotly.graph_objects as go
import numpy
from dash.dependencies import Input, Output


class RenderableMARA(object):

    def __init__(self, activity, cellcond, tcname):
        self.activity = activity
        self.cellcond = cellcond
        self.activity = self.activity.set_index('index')
        self.tcname = tcname

    def render_tc(self, name_of_gene):

        #Test a precise match first
        #found_genes = self.map_ensembl_genesym[self.map_ensembl_genesym['Associated Gene Name lc']==name_of_gene.lower()]

        #Otherwise find all genes starting with this name
        #if found_genes.empty or name_of_gene=="":
        #    found_genes = self.map_ensembl_genesym[self.map_ensembl_genesym['Associated Gene Name lc'].str.startswith(name_of_gene.lower())]

        #if found_genes.empty or name_of_gene=="":
        #    return "Please select a gene to display"
        if name_of_gene=="":
            return "Please select a motif to display"

        #selected_gene = found_genes.iloc[0]
        #selected_gene_id = selected_gene["Ensembl Gene ID"]
        #selected_gene_sym = selected_gene['Associated Gene Name']
        selected_gene = name_of_gene

        ################## todo , check in a better way. !!!!!!
        #if not selected_gene_id in self.df_data.index:
        #    return "Gene not in dataset"

        #return "Gene not in dataset"

        activity=self.activity
        cellcond=self.cellcond
        themotif=name_of_gene

        selected_gene_sym = themotif


        if not(themotif in activity.transpose()):
            return "Motif not in dataset"

        ## Prepare activity levels
        df = pandas.DataFrame(
            {'activity': activity.transpose()[themotif].to_numpy().tolist(),
             'time': cellcond["time"].array.to_numpy().tolist(),
             'type': cellcond["type"]
            }, 
            columns = ['activity', 'time','type'])
        df = df.groupby(['time','type']).mean().reset_index()
        df.sort_values(by=['time'])

        df0 = df.loc[df['type'].isin(['Th0','Naive'])]
        df2 = df.loc[df['type'].isin(['Th2','Naive'])]

        return dcc.Graph(
            id='mara',
            figure={
                'data': [
                    go.Scatter(
                        x=df2['time'],
                        y=df2['activity'],
                        name='Th2', mode='lines+markers'),
                    go.Scatter(
                        x=df0['time'],
                        y=df0['activity'],
                        name='Th0', mode='lines+markers')
                ],
                'layout': go.Layout(
                    title=self.tcname+' -- '+selected_gene_sym,
                    xaxis={'title' : 'Time (h)'},
                    yaxis={'title': 'Activity level'},
                )
            })
