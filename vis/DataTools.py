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
        return Main.data.run

    def get_protein(self, cell, props):
        get_fcn = Main.eval('ProteinStoreMod.get')
        return get_fcn(cell.proteins, props)

    def get_props_str(self, props):
        to_str_fcn = Main.eval('ProteinPropsMod.to_str')
        return to_str_fcn(props)

    def get_protein_info_for_indiv(self, index):
        self.get_indiv(index)
        Main.eval('info = DataMod.get_protein_info(indiv)')
        props_list = list(map(lambda row: row, Main.info))
        
        return props_list
    
