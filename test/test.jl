using DataMod

data = Data("data");

ea_step = 0
indiv_index = 1

indiv = DataMod.get_indiv(data, ea_step, indiv_index)
info = DataMod.get_protein_info(indiv)
