module CompressionMod

import Blosc
import CodecXz
import CodecZlib
import CodecBzip2
import CodecZstd

export compress, decompress

@enum CompressionAlg::UInt8 blosc codecxz codeczlib codecbzip2 codeczstd

# function check_alg(run)
#     if run.compression_alg ∉ CompressionMod.routing_dict
#         throw(ErrorException("Unrecognized compression_alg argument: $(run.compression_alg)"))
#     end
# end

function compress(compression_alg::CompressionAlg, data::Array{UInt8, 1})
    CompressionMod.routing_dict[compression_alg][1](data)
end

function decompress(compression_alg::CompressionAlg, data::Array{UInt8, 1})
    CompressionMod.routing_dict[compression_alg][2](data)
end

function compress_blosc(data::Array{UInt8, 1})
    Blosc.compress(data)
end

function decompress_blosc(data::Array{UInt8, 1})
    Blosc.decompress(UInt8, data)
end

function compress_codecxz(data::Array{UInt8, 1})
    CodecXz.transcode(CodecXz.XzCompressor, data)
end

function decompress_codecxz(data::Array{UInt8, 1})
    CodecXz.transcode(CodecXz.XzDecompressor, data)
end

function compress_codeczlib(data::Array{UInt8, 1})
    CodecXz.transcode(CodecZlib.GzipCompressor, data)
end

function decompress_codeczlib(data::Array{UInt8, 1})
    CodecXz.transcode(CodecZlib.GzipDecompressor, data)
end

function compress_codecbzip2(data::Array{UInt8, 1})
    CodecXz.transcode(CodecBzip2.Bzip2Compressor, data)
end

function decompress_codecbzip2(data::Array{UInt8, 1})
    CodecXz.transcode(CodecBzip2.Bzip2Decompressor, data)
end

function compress_codeczstd(data::Array{UInt8, 1})
    CodecXz.transcode(CodecZstd.ZstdCompressor, data)
end

function decompress_codeczstd(data::Array{UInt8, 1})
    CodecXz.transcode(CodecZstd.ZstdDecompressor, data)
end

#note: this must be defined last
routing_dict = Dict{CompressionAlg, Tuple{Function, Function}}(
    blosc => (CompressionMod.compress_blosc, CompressionMod.decompress_blosc),
    codecxz => (CompressionMod.compress_codecxz, CompressionMod.decompress_codecxz),
    codeczlib => (CompressionMod.compress_codeczlib, CompressionMod.decompress_codeczlib),
    codecbzip2 => (CompressionMod.compress_codecbzip2, CompressionMod.decompress_codecbzip2),
    codeczstd => (CompressionMod.compress_codeczstd, CompressionMod.decompress_codeczstd)
)

end
