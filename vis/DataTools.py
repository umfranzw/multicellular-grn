from julia import Main

class DataTools():
    def __init__(self, filename):
        Main.using('Distributed')
        Main.eval('@everywhere push!(LOAD_PATH, "/home/wayne/Documents/school/thesis/multicellular-grn")')
        Main.using('DataMod')
        Main.eval('data = Data("{}")'.format(filename))

    def close(self):
        Main.eval('DataMod.close(data)')

    def get_tree(self, index):
        Main.eval('tree = DataMod.get_tree(data, {}, {}, {})'.format(*index))
        return Main.tree

    def get_indiv(self, index):
        Main.eval('indiv = DataMod.get_indiv(data, {}, {}, {})'.format(*index))
        return Main.indiv

    def get_run(self):
        return Main.data.run
    
