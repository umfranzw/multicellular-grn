using DataMod
using CellTreeMod

data = Data("data");

# tree = DataMod.get_tree(data, 10, 1, 1)
# cell = tree.root
# png_data = DataMod.build_graph_for_cell(data, 10, 1, cell)
# f = open("test.png", "w")
# write(f, png_data)
# close(f)

#table = DataMod.get_gs_table_data(data, cell, 0, 1, 1)

indiv = DataMod.get_indiv(data, 19, 29, 5)
size = CellTreeMod.size(indiv.cell_tree)
