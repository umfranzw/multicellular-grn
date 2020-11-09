using CellMod
using GeneMod
using ProteinMod
using Serialization

f = open("/home/umfranzw/multicellular-grn/data/test6", "r")

size = read(f, Int64)
data = read(f, size)
protein = Serialization.deserialize(IOBuffer(data))

size = read(f, Int64)
data = read(f, size)
gene = Serialization.deserialize(IOBuffer(data))

size = read(f, Int64)
data = read(f, size)
cell = Serialization.deserialize(IOBuffer(data))

close(f)
