class RenderableTC(object):

    def __init__(self, df_data, df_meta, map_ensembl_genesym, tcname):
        self.tcname   = tcname
        self.df_data  = df_data
        self.df_meta  = df_meta
        self.df_times = df_meta['hours'][df_meta['Cell Type'] == "Th2"].tolist()

        ### Create a lower case version of symbol for fast search later
        map_ensembl_genesym["Associated Gene Name lc"] = map_ensembl_genesym['Associated Gene Name'].str.lower()
        map_ensembl_genesym["Associated Gene Name lc"] = map_ensembl_genesym["Associated Gene Name lc"].fillna(value="")

        self.map_ensembl_genesym = map_ensembl_genesym

    def render_tc(self, name_of_gene):
        print(self.tcname)
        print(self.map_ensembl_genesym.head())
        print(self.map_ensembl_genesym['Associated Gene Name lc'])
        print(name_of_gene.lower())

        #Test a precise match first
        found_genes = self.map_ensembl_genesym[self.map_ensembl_genesym['Associated Gene Name lc']==name_of_gene.lower()]

        #Otherwise find all genes starting with this name
        if found_genes.empty or name_of_gene=="":
            found_genes = self.map_ensembl_genesym[self.map_ensembl_genesym['Associated Gene Name lc'].str.startswith(name_of_gene.lower())]

        print(found_genes)
        if found_genes.empty or name_of_gene=="":
            return "Please select a gene to display"

        selected_gene = found_genes.iloc[0]
        selected_gene_id = selected_gene["Ensembl Gene ID"]
        selected_gene_sym = selected_gene['Associated Gene Name']

        if not selected_gene_id in self.df_data.index:
            return "Gene not in dataset"

        return dcc.Graph(
            id='tcmouse',
            figure={
                'data': [
                    go.Scatter(
                        x=self.df_times,
                        y=self.df_data.loc[selected_gene_id,numpy.array(self.df_meta['Cell Type']=="Th0")],
                        name='Th0', mode='lines+markers'),
                    go.Scatter(
                        x=self.df_times,
                        y=self.df_data.loc[selected_gene_id,numpy.array(self.df_meta['Cell Type']=="Th2")],
                        name='Th2', mode='lines+markers')
                ],
                'layout': go.Layout(
                    title=self.tcname+' -- '+selected_gene_sym,
                    xaxis={'title' : 'Time (h)'},
                    yaxis={'title': 'Expression level'},
                )
            })
