import RunMod
using IndividualMod
using Serialization
using CodecZlib

function compress()
    run = RunMod.get_first_run()
    pop1 = map(i -> IndividualMod.rand_init(run), 1:run.pop_size)
    #pop2 = map(i -> IndividualMod.rand_init(run), 1:run.pop_size)

    buf = IOBuffer()
    Serialization.serialize(buf, pop1)
    compressed = CodecZlib.transcode(GzipCompressor, buf.data)

    out_stream = open("test.bin", "w")
    write(out_stream, Int64(length(compressed)))
    write(out_stream, compressed)
    close(out_stream)
end

function decompress()
    in_stream = open("test.bin", "r")
    len = read(in_stream, Int64)
    compressed = Array{UInt8, 1}()
    for i in 1:len
        push!(compressed, read(in_stream, UInt8))
    end
    close(in_stream)

    decompressed = CodecZlib.transcode(GzipDecompressor, compressed)
    buf = IOBuffer(decompressed; read=true)
    pop = Serialization.deserialize(buf)
end
