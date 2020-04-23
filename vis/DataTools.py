from julia import Main

class DataTools():
    def __init__(self, filename):
        Main.using('Distributed')
        Main.eval('@everywhere push!(LOAD_PATH, "/home/wayne/Documents/school/thesis/multicellular-grn")')
        Main.using('DataMod')
        Main.using('ProteinStoreMod')
        Main.using('ProteinPropsMod')
        Main.eval('data = Data("{}")'.format(filename))

    def close(self):
        Main.eval('DataMod.close(data)')

    def get_tree(self, index):
        Main.eval('tree = DataMod.get_tree(data, {}, {}, {})'.format(*index))
        return Main.tree

    def get_indiv(self, index):
        #note: get_indiv does not require the reg_step arg to be passed
        Main.eval('indiv = DataMod.get_indiv(data, {}, {})'.format(*index[:2]))
        return Main.indiv

    def get_run(self):
        Main.eval('run = data.run')
        return Main.run

    def get_protein(self, cell, props):
        get_fcn = Main.eval('ProteinStoreMod.get')
        return get_fcn(cell.proteins, props)

    def get_props_str(self, props):
        to_str_fcn = Main.eval('ProteinPropsMod.to_str')
        return to_str_fcn(props)

    def get_protein_info_for_tree(self, index):
        self.get_tree(index)
        Main.eval('info = DataMod.get_protein_info_for_tree(tree)')
        
        return Main.info

    def get_cell_children(self, cell):
        Main.cell = cell
        Main.eval('children = cell.children')

        return Main.children

    def get_root_cell(self, tree):
        Main.tree = tree
        Main.eval('root = tree.root')
        
        return Main.root

    def get_probs_info_for_cell(self, cell):
        Main.cell = cell
        Main.eval('info = DataMod.get_probs_info_for_cell(cell)')

        return Main.info

    def get_num_genes(self, cell):
        Main.cell = cell
        Main.eval('num_genes = length(cell.gene_states)')

        return Main.num_genes

    def get_sensor_concs(self, cell):
        concs = {}
        if cell is not None:
            Main.cell = cell
            Main.eval('sensor_concs = DataMod.get_sensor_concs(cell)')
            concs = Main.sensor_concs

        return concs

    def get_interaction_graph(self, ea_step, pop_index, cell):
        if cell is None:
            pixmap_data = []
        else:
            Main.cell = cell
            Main.ea_step = ea_step
            Main.pop_index = pop_index
            Main.eval('pixmap_data = DataMod.build_graph_for_cell(data, ea_step, pop_index, cell)')
            pixmap_data = Main.pixmap_data

        return pixmap_data

    def get_gs_table_data(self, cell, index):
        if cell is None:
            table_data = []
        else:
            Main.cell = cell
            Main.ea_step = index[0]
            Main.indiv_index = index[1]
            Main.reg_step = index[2]
            Main.eval('table_data = DataMod.get_gs_table_data(data, cell, ea_step, indiv_index, reg_step)')
            table_data = Main.table_data

        return table_data

    def get_best_fitnesses(self):
        Main.eval('bests, gen_bests = DataMod.get_best_fitnesses(data)')
        return (Main.bests, Main.gen_bests)
