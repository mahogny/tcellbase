# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objects as go


app = dash.Dash()
server = app.server
#Opens and reads file from url
df = pd.read_csv('https://raw.githubusercontent.com/RasmusBurge/DBT-T-Cell/master/rawdata.csv', sep=";")

#Decides the layout of web app with titles and labels
app.layout = html.Div(children=[
    html.H1(children='Försök till visualisering'),
    html.Div(children='''
        Grupp nr 6
    '''),
#Grap line code with graph id
    dcc.Graph(
        id='Försök till datavisualisering',
        figure={
            #Data points are chosen by defining which columns are y and x and then plotted from the csv-file through a scatter plot with lines.
            # Names are put to each data series .
            'data': [
                go.Scatter(x=df.Tidpunkter, y=df.ENSG00000283297, name='ENSG00000283297', mode='lines+markers'),
                go.Scatter(x=df.Tidpunkter, y=df.ENSG00000283573, name='ENSG00000283573 ', mode='lines+markers'),
                go.Scatter(x=df.Tidpunkter, y=df.ENSG00000282961, name='ENSG00000282961 ', mode='lines+markers')

            ],
            #Title label for graph, and for y and x axis.
            'layout': go.Layout(
                title='TH0 Genexpression över tid',
                xaxis={'title' : 'Tid(h)'},
                yaxis={'title': 'Expression .../...'},
            )
        })])


