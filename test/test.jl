using DataMod
using ChainGraphMod

data = Data("data");

tree = DataMod.get_tree(data, 1, 1, 1)
cell = tree.root
graph = DataMod.build_graph_for_cell(data, 1, 1, cell)
png_data = ChainGraphMod.plot(graph)
f = open("test.png", "w")
write(f, png_data)
close(f)
