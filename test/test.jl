using DataMod
using CellTreeMod

run, ea_pops, reg_trees = DataMod.read_data("data");

ea_step = 1
reg_step = 1
indiv_index = 1

#indiv = ea_pops["before_reg_sim"][ea_step][indiv_index]
tree = reg_trees["after_reg_step"][ea_step][indiv_index][1]
cell1 = tree.root
#cell5 = CellTreeMod.get_bf_node(tree, 5)

tree = reg_trees["after_reg_step"][ea_step][indiv_index][2]
cell2 = tree.root

tree = reg_trees["after_reg_step"][ea_step][indiv_index][3]
cell3 = tree.root

