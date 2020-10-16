import socket
if socket.gethostname() == 'ibis':
    from julia.api import Julia
    jl = Julia(compiled_modules=False)
    module_path = '/home/umfranzw/multicellular-grn'
else:
    module_path = '/home/wayne/Documents/school/thesis/multicellular-grn'
    
from julia import Main
from Cache import Cache

class DataTools():
    CACHE_SIZE = 20

    def __init__(self, filename):
        Main.using('Distributed')
        Main.eval('@everywhere push!(LOAD_PATH, "{}")'.format(module_path))
        Main.using('DataMod')
        Main.using('ProteinStoreMod')
        Main.using('ProteinPropsMod')
        Main.using('CellTreeMod')
        Main.using('TrackerMod')
        Main.eval('data = Data("{}")'.format(filename))

        Main.eval('state_time_instances = instances(TrackerMod.StateTime)')
        Main.eval('state_time_names = map(string, state_time_instances)')
        Main.eval('state_time_dict = Dict{String, TrackerMod.StateTime}(zip(state_time_names, state_time_instances))')
        self.state_time_dict = Main.state_time_dict
        self.cache = Cache(DataTools.CACHE_SIZE)

    def close(self):
        Main.eval('DataMod.close(data)')

    def get_tree(self, index, state_time):
        indiv = self.get_indiv(index, state_time)
        Main.indiv = indiv
        Main.eval('tree = indiv.cell_tree')
        
        return Main.tree

    def get_indiv(self, index, state_time):
        key = (index, state_time)
        if key not in self.cache:
            #note: get_indiv does not require the reg_step arg to be passed
            Main.ea_step = index[0]
            Main.pop_index = index[1]
            Main.reg_step = index[2]
            Main.state_time = state_time
            Main.eval('indiv = DataMod.get_indiv(data, ea_step, pop_index, reg_step, state_time)')
            self.cache[key] = Main.indiv
        #else:
            #still need to set this, since other methods rely on it...
            #Main.indiv = self.cache[key]
        
        return self.cache[key]

    def get_run(self):
        Main.eval('run = data.run')
        return Main.run

    def get_protein(self, cell, props):
        #get_fcn = Main.eval('ProteinStoreMod.get')
        #return get_fcn(cell.proteins, props)
        Main.proteins = cell.proteins
        Main.props = props
        Main.eval('protein = ProteinStoreMod.get(proteins, props)')
        
        return Main.protein

    def get_props_str(self, props, is_initial):
        Main.props = props
        Main.is_initial = is_initial
        Main.eval('props_str = ProteinPropsMod.to_str(props, is_initial)')
        
        return Main.props_str

    #checkpoint is a value from self.state_time_dict
    def get_protein_info_for_indiv(self, index, state_time):
        indiv = self.get_indiv(index, state_time)
        Main.indiv = indiv
        Main.eval('info = DataMod.get_protein_info_for_indiv(indiv)')
        
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

    def get_interaction_graph(self, index, cell):
        if cell is None:
            pixmap_data = []
        else:
            Main.cell = cell
            Main.ea_step = index[0]
            Main.pop_index = index[1]
            Main.reg_step = index[2]
            Main.eval('pixmap_data = DataMod.build_graph_for_cell(data, ea_step, pop_index, reg_step, cell)')
            pixmap_data = Main.pixmap_data

        return pixmap_data

    def get_neighbour_graph(self, index):
        Main.ea_step = index[0]
        Main.pop_index = index[1]
        Main.reg_step = index[2]
        Main.eval('pixmap_data = DataMod.build_neighbour_comm_graph(data, ea_step, pop_index, reg_step)')
        pixmap_data = Main.pixmap_data

        return pixmap_data
        
    def save_all_interaction_graphs(self, index, cell, path):
        Main.ea_step = index[0]
        Main.pop_index = index[1]
        Main.cell = cell
        Main.path = path
        Main.eval('DataMod.save_all_graphs_for_cell(data, ea_step, pop_index, cell, path)')

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

    def save_all_gs_table_data(self, cell, index, filename):
        if cell is not None:
            Main.cell = cell
            Main.ea_step = index[0]
            Main.indiv_index = index[1]
            Main.filename = filename
            Main.eval('DataMod.save_all_gs_table_data(data, cell, ea_step, indiv_index, filename)')

    def get_best_fitnesses(self):
        Main.eval('bests, gen_bests, gen_avgs = DataMod.get_best_fitnesses(data)')
        return (Main.bests, Main.gen_bests, Main.gen_avgs)

    def get_indiv_fitness(self, index):
        indiv = self.get_indiv(index, self.state_time_dict['AfterBind'])
        Main.indiv = indiv
        Main.eval('fitness = indiv.fitness')
        
        return Main.fitness

    def get_indiv_code(self, index):
        indiv = self.get_indiv(index, self.state_time_dict['AfterBind'])
        Main.indiv = indiv
        Main.eval('code = CellTreeMod.to_expr_str(indiv.cell_tree)')

        return Main.code

    def get_gene_scores(self, index):
        indiv = self.get_indiv(index, self.state_time_dict['AfterBind'])
        Main.indiv = indiv
        Main.eval('scores = indiv.gene_scores')
        
        return Main.scores

    def get_run_best_index(self):
        Main.eval('best_index = data.run_best.index')
        
        return Main.best_index

    def get_gene_descs(self, pop_index):
        Main.pop_index = pop_index
        Main.eval('table, max_cols = DataMod.get_gene_descs(data, pop_index)')

        return Main.table, Main.max_cols
    
    def export_gene_descs(self, pop_index, filename):
        Main.filename = filename
        Main.pop_index = pop_index
        Main.eval('DataMod.export_gene_descs(data, pop_index, filename)')

    def get_base_seed(self):
        Main.eval('base_seed = DataMod.get_base_seed(data)')

        return Main.base_seed
