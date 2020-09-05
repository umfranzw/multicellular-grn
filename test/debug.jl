using DataMod

data = Data("data")
i1 = DataMod.get_indiv(data, 0, 1, 1)
root1 = i1.cell_tree.root

i2 = DataMod.get_indiv(data, 0, 1, 2)
root2 = i2.cell_tree.root
