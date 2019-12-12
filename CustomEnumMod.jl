module CustomEnumMod

using RunMod
using DataStructures

import Base.getproperty
import Base.iterate
import Random.rand
export CustomEnum

mutable struct CustomEnum
    items::OrderedDict{Symbol, UInt8}
    start::UInt8
    
    function CustomEnum(symbols::Array{Symbol}; start::UInt8=UInt8(1))
        items = OrderedDict{Symbol, UInt8}()
        
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

function get_val(enum::CustomEnum, sym::Symbol)
    getfield(enum, :items)[sym]
end

function rand(config::Config, enum::CustomEnum)
    items = getfield(enum, :items)
    
    rand(config.rng, items) #returns a Pair{Symbol, Int64}
end

iterate(enum::CustomEnum) = Base.iterate(getfield(enum, :items))
iterate(enum::CustomEnum, state) = Base.iterate(getfield(enum, :items), state)
    
end
