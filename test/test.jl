using DataMod
using ChainGraphMod

data = Data("data");

tree = DataMod.get_tree(data, 0, 1, 0)
cell = tree.root
# png_data = DataMod.build_graph_for_cell(data, 1, 0, cell)
# f = open("test.png", "w")
# write(f, png_data)
# close(f)

table = DataMod.get_gs_table_data(data, cell, 0, 1, 0)
