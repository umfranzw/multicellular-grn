using DataMod
using IndividualMod
import TrackerMod
import CheckpointMod
import Serialization
import CompressionMod

function ser_comp(indiv::Individual)
    buf = IOBuffer()
    Serialization.serialize(buf, indiv)

    CompressionMod.compress(indiv.config.run.compression_alg, take!(buf))
end

function decomp_deser(data::Array{UInt8, 1}, alg::CompressionMod.CompressionAlg)
    data = CompressionMod.decompress(alg, data)
    
    Serialization.deserialize(IOBuffer(data))
end

data = Data("data")
#problem: 1, 6, 2
i1 = DataMod.get_indiv(data, 1, 6, 1, TrackerMod.AfterBind)
i2 = DataMod.get_indiv(data, 1, 6, 2, TrackerMod.AfterBind)
# compression_alg = i1.config.run.compression_alg
# i2 = DataMod.get_indiv(data, DataMod.IndexKey((1, 1, 1), TrackerMod.AfterBind))
# data1 = ser_comp(i1)
# data2 = ser_comp(i2)
# info = CheckpointMod.compress(data1, data2)

# data4 = CheckpointMod.decompress(data1, info)
# i3 = decomp_deser(data4, compression_alg)
