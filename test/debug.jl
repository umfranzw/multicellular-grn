# using DataMod
# using DeepDiffs
# import TrackerMod
import Serialization
# using CodecXz

# function print_array(array)
#     for item in array
#         println(item)
#     end
# end

# function update_info(first::Int64, removed::Array{UnitRange{Int64}, 1}, added::Array{Tuple{Int64, Array{UInt8, 1}}, 1}, result::VectorDiff{Array{UInt8,1},Array{UInt8,1}}, chunk2::Array{UInt8, 1})
#     offset = first - 1
#     start = nothing
#     stop = nothing
#     for index in result.removed
#         if stop == nothing
#             start = index
#             stop = index
#         elseif index == stop + 1
#             stop = index
#         else #index > stop, start new range
#             push!(removed, start + offset : stop + offset)
#             start = index
#             stop = index
#         end
#     end
#     #push last one, if necessary
#     if stop != nothing
#         push!(removed, start + offset : stop + offset)
#     end

#     start = nothing
#     stop = nothing
#     for index in result.added
#         if stop == nothing
#             start = index
#             stop = index
#         elseif index == stop + 1
#             stop = index
#         else #index > stop, start new range
#             push!(added, (start + offset, chunk2[start:stop]))
#             start = index
#             stop = index
#         end
#     end
#     #push last one, if necessary
#     if stop != nothing
#         push!(added, (start + offset, chunk2[start:stop]))
#     end
# end

# function test(d1::Array{UInt8, 1}, d2::Array{UInt8, 1}, chunk_size::Int64)
#     first = 1
#     last = chunk_size
#     removed = Array{UnitRange{Int64}, 1}()
#     added = Array{Tuple{Int64, Array{UInt8, 1}}, 1}()
#     # i = 1
    
#     while last <= length(d1) && last <= length(d2)
#         chunk1 = d1[first:last]
#         chunk2 = d2[first:last]
#         result = deepdiff(chunk1, chunk2)
#         # println("chunk $(i)")
#         # println("first, last: $(first), $(last)")
#         # println("removed: $(result.removed)")
#         # println("added: $(result.added)")
#         # println()
        
#         update_info(first, removed, added, result, chunk2)

#         first = last + 1
#         last = first + chunk_size - 1
#         # i += 1
#     end

#     #d2 is smaller
#     if length(d2) < length(d1)
#         # println("length(d2) < length(d1)")

#         #finish last (non-full) chunk of d2, if any
#         last = length(d2)
#         if last >= first
#             chunk1 = d1[first:last]
#             chunk2 = d2[first:last]
#             result = deepdiff(chunk1, chunk2)

#             # println("first, last: $(first), $(last)")
#             # println("removed: $(result.removed)")
#             # println("added: $(result.added)")
            
#             update_info(first, removed, added, result, chunk2)
#             first = last + 1
#         end
        
#         #finish remainder of d1, if any
#         if first <= length(d1)
#             push!(removed, first:length(d1))
#         end
#         # println()

#         #d1 is smaller
#     else
#         # println("length(d1) < length(d2)")
        
#         #finish last (non-full) chunk of d1, if any
#         last = length(d1)
#         if last >= first
#             chunk1 = d1[first:last]
#             chunk2 = d2[first:last]
#             result = deepdiff(chunk1, chunk2)

#             # println("first, last: $(first), $(last)")
#             # println("removed: $(result.removed)")
#             # println("added: $(result.added)")
            
#             update_info(first, removed, added, result, chunk2)
#             first = last + 1
#         end

#         #finish remainder of d2, if any
#         if first <= length(d2)
#             push!(added, (first, d2[first:end]))
#         end
#         # println()
#     end

#     removed, added
# end

# function untest(d1::Array{UInt8, 1}, removed::Array{UnitRange{Int64}, 1}, added::Array{Tuple{Int64, Array{UInt8, 1}}, 1})
#     # println()
#     result = deepcopy(d1)

#     #do removals (proceed in reverse order, so the current iteration does not affect parts of the list the next one will operate on (since removed is in sorted order by start index)
#     for range in Iterators.reverse(removed)
#         # println("Deleted: $(range)")

#         deleteat!(result, range)
        
#         # println("d2:\t$(d2)")
#         # println("result:\t$(result)")
#         # println()
#     end

#     #do additions - no need to reverse here
#     for (start, data) in added
#         # println("Added: $((start, data))")

#         i = 0
#         while i < length(data)
#             insert!(result, start + i, data[i + 1])
#             i += 1
#         end

#         # println("d2:\t$(d2)")
#         # println("result:\t$(result)")
#         # println()
#     end
    
#     result
# end

# function orig_compress(d1::Array{UInt8, 1}, d2::Array{UInt8, 1})
#     d1c = CodecXz.transcode(CodecXz.XzCompressor, d1)
#     d2c = CodecXz.transcode(CodecXz.XzCompressor, d2)

#     (d1c, d2c)
# end

# chunk_size = 512

# data = Data("data")
# i1 = DataMod.get_indiv(data, 1, 6, 1, TrackerMod.IndivStateAfterBind)
# i2 = DataMod.get_indiv(data, 1, 6, 51, TrackerMod.IndivStateAfterBind)
# buf1 = IOBuffer()
# buf2 = IOBuffer()
# Serialization.serialize(buf1, i1)
# Serialization.serialize(buf2, i2)
# d1 = take!(buf1)
# d2 = take!(buf2)

# osize = Base.summarysize((d1, d2)) / 2^20
# println("Original size: $(osize)")

# println("Xz Compressing")
# @time d3 = orig_compress(d1, d2)
# ocsize = Base.summarysize(d3) / 2^20
# println("Xz compressed size: $(ocsize)")
# println("ratio: $(ocsize / osize)")

# println("Checkpoint Compressing")
# @time removed, added = test(d1, d2, chunk_size)
# csize = Base.summarysize((d1, removed, added)) / 2^20
# println("size: $(csize)")
# println("ratio: $(csize / osize)")

# println("Xz Compressing the Checkpoint data")
# temp = IOBuffer()
# Serialization.serialize(temp, (d1, removed, added))

# d4 = CodecXz.transcode(CodecXz.XzCompressor, temp.data)
# best_size = Base.summarysize(d4) / 2^20
# println("size: $(best_size)")
# println("ratio: $(best_size / osize)")

# # println("removed:")
# # print_array(removed)
# # println("added:")
# # print_array(added)


# println("decompressing")
# @time reconstructed = untest(d1, removed, added)
# # println("d2:\t\t$(d2)")
# # println("reconstructed:\t$(reconstructed)")
# println("correctness check: $(d2 == reconstructed)")

# buf = IOBuffer(reconstructed)
# indiv = Serialization.deserialize(buf)

mutable struct MutableTest
    x::Array{Int64, 1}
    y::Int64
    z::Dict{Char, Int64}
end

struct ConstTest
    x::Array{Int64, 1}
    y::Int64
    z::Dict{Char, Int64}
end

mt = MutableTest(collect(1:1000), 2, Dict{Char, Int64}('a' => 1, 'b' => 2))
ct = ConstTest(collect(1:1000), 2, Dict{Char, Int64}('a' => 1, 'b' => 2))

mt_buf = IOBuffer()
Serialization.serialize(mt_buf, mt)
mt_bytes = take!(mt_buf)

ct_buf = IOBuffer()
Serialization.serialize(ct_buf, ct)
ct_bytes = take!(ct_buf)

println("mt: $(length(mt_bytes))")
println("ct: $(length(ct_bytes))")
