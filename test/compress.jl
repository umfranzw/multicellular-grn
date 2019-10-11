import RunMod
using IndividualMod
using SymMod
using ProteinMod
using ProteinPropsMod
using CellMod
using Serialization
using CodecZlib

pop_size = 1000
num_cells = 10
num_proteins = 50
    
run = RunMod.get_first_run()
pop = Array{Individual, 1}()
for i in 1:pop_size
    indiv = IndividualMod.rand_init(run)
    
    indiv.root_cell.sym = Sym(:+, SymMod.FcnCall, -1)
    for j in 1:num_cells
        Cell(run, indiv.genes, indiv.root_cell, Sym(i, SymMod.IntConst, 0))
    end

    for j in 1:num_proteins
        protein = Protein(run, ProteinProps(ProteinPropsMod.Reg, ProteinPropsMod.Intra, ProteinPropsMod.Activate, ProteinPropsMod.B), true)
        push!(indiv.initial_cell_proteins, protein)
    end

    push!(pop, indiv)
end

buf = IOBuffer()
serialized_time = @elapsed Serialization.serialize(buf, pop)
serialized_size = buf.size / 2^20
println("Serialized size: $(serialized_size) MB")
println("Serialization time: $(serialized_time)")

compression_time = @elapsed compressed = transcode(GzipCompressor, buf.data)
comp_size = length(compressed) / 2^20
println("Compressed size: $(comp_size) MB")
println("Compression time: $(compression_time)")
println()

comp_factor = serialized_size / comp_size
println("Compression factor: $(comp_factor)")
println()

decompression_time = @elapsed decompressed = transcode(GzipDecompressor, compressed)
println("Decompression time: $(decompression_time)")

new_buf = IOBuffer(decompressed; read=true, write=true)
deserialization_time = @elapsed new_pop = Serialization.deserialize(new_buf)
println("Deserialization time: $(deserialization_time)")
    
