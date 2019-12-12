module CustomEnumMod

using RunMod

import Base.getproperty
export CustomEnum

mutable struct CustomEnum
    items::Dict{Symbol, UInt8}
    start::UInt8
    
    function CustomEnum(symbols::Array{Symbol}; start::UInt8=UInt8(1))
        items = Dict{Symbol, UInt8}()
        
        for i in 1:length(symbols)
            items[symbols[i]] = i + start - 1
        end

        new(items, start)
    end
end

function getproperty(enum::CustomEnum, name::Symbol)
    #note: using getfield() allows us to avoid the dot syntax (enum.items), which would cause infinite recursion
    items = getfield(enum, :items)
    
    if name in keys(items)
        return items[name]
    else
        throw(UndefVarError(name))
    end
end

function rand_val(config::Config, enum::CustomEnum)
    items = getfield(enum, :items)
    
    Random.rand(config.rng, keys(items))
end
    
end
