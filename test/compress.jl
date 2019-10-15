import RunMod
using IndividualMod
using SymMod
using ProteinMod
using ProteinPropsMod
using CellMod
using Serialization
using CodecZlib
using Printf

pop_size = 1000
num_cells = 10
num_proteins = 50

println("pop_size: $(pop_size)")
println("num_cells: $(num_cells)")
println("num_proteins: $(num_proteins)")
println()

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
@printf("Serialized size: %0.2f MB\n", serialized_size)
@printf("Serialization time: %0.2f sec\n", serialized_time)

compression_time = @elapsed compressed = transcode(GzipCompressor, buf.data)
comp_size = length(compressed) / 2^20
@printf("Compressed size: %0.2f MB\n", comp_size)
@printf("Compression time: %0.2f sec\n", compression_time)
comp_factor = serialized_size / comp_size
@printf("Compression factor: %0.2f\n", comp_factor)
println()

decompression_time = @elapsed decompressed = transcode(GzipDecompressor, compressed)
@printf("Decompression time: %0.2f\n", decompression_time)

new_buf = IOBuffer(decompressed; read=true, write=true)
deserialization_time = @elapsed new_pop = Serialization.deserialize(new_buf)
@printf("Deserialization time: %0.2f\n", deserialization_time)
    
