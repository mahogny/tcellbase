# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objects as go
import numpy
import math


#### See for reference
# https://dash-gallery.plotly.host/dash-clinical-analytics/
# https://github.com/plotly/dash-sample-apps/blob/master/apps/dash-clinical-analytics/app.py





class RenderableHeatmap(object):

    def __init__(self, heatmap_data, map_ensembl_genesym, plotname, conditions):
        self.plotname     = plotname
        self.heatmap_data = heatmap_data

        ### Create a lower case version of symbol for fast search later
        map_ensembl_genesym["Associated Gene Name lc"] = map_ensembl_genesym['Associated Gene Name'].str.lower()
        map_ensembl_genesym["Associated Gene Name lc"] = map_ensembl_genesym["Associated Gene Name lc"].fillna(value="")

        self.map_ensembl_genesym = map_ensembl_genesym

        self.conditions=conditions


    def render_heatmap(self, name_of_gene):

        # find all genes starting with this name
        found_genes = self.map_ensembl_genesym[
            self.map_ensembl_genesym['Associated Gene Name lc'].str.startswith(name_of_gene.lower())]
        if found_genes.shape[0] > 300:  ### quick hack, do not give back too long gene lists
            found_genes = self.map_ensembl_genesym[
                self.map_ensembl_genesym['Associated Gene Name lc'].str.startswith("6666666666666666666666666")]

        ## Extract subset of gene expression
        hdata = self.heatmap_data[self.heatmap_data.index.isin(found_genes["Ensembl Gene ID"].tolist())]

        ## Sort by gene symbol
        hdata = hdata.sort_values(by="Associated Gene Name", ascending=False)

        ## Decide axis names and IDs
        x_axis = self.conditions
        y_axis = hdata.index.tolist()
        y_axis_sym = hdata["Associated Gene Name"].tolist()

        ## Build up the heatmap content
        z = numpy.zeros((len(y_axis), len(x_axis)))
        annotations = []
        for ind_y, geneid in enumerate(y_axis):
            filtered_df = hdata[hdata.index == geneid]
            for ind_x, x_val in enumerate(x_axis):
                gene_exp = filtered_df[x_val][0].item()
                gene_exp = math.log10(gene_exp + 1)
                gene_exp_string = "{0:.1f}".format(gene_exp)

                # Currently stored in the file; but this should be changed later. Better look up on the fly
                gene_sym = filtered_df["Associated Gene Name"][0]

                z[ind_y][ind_x] = gene_exp
                annotation_dict = dict(
                    showarrow=False,
                    text="<b> " + gene_exp_string + " <b>",
                    xref="x",
                    yref="y",
                    x=x_val,
                    y=gene_sym,
                    font=dict(family="sans-serif"),
                )

                annotations.append(annotation_dict)

        # Heatmap
        # hovertemplate = "<b> %{y}  %{x} <br><br> %{z} Patient Records"

        data = [
            dict(
                x=x_axis,
                y=y_axis_sym,
                z=z,
                type="heatmap",
                name="",
                # hovertemplate=hovertemplate,
                showscale=False,
                colorscale=[[0, "#caf3ff"], [1, "#2c82ff"]],
            )
        ]

        layout = dict(
            title=self.plotname,
            margin=dict(l=70, b=50, t=50, r=50),
            modebar={"orientation": "v"},
            font=dict(family="Open Sans"),
            annotations=annotations,
            # shapes=shapes,
            xaxis=dict(
                side="top",
                ticks="",
                ticklen=2,
                tickfont=dict(family="sans-serif"),
                tickcolor="#ffffff",
            ),
            yaxis=dict(
                side="left", ticks="", tickfont=dict(family="sans-serif"), ticksuffix=" "
            ),
            hovermode="closest",
            showlegend=False,
        )

        return {"data": data, "layout": layout}


