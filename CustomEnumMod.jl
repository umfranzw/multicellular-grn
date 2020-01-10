module CustomEnumMod

using RunMod
using DataStructures
using Mustache

import Base.getproperty
import Base.iterate
import Random.rand
import Base.length
import Base.show
import Base.Int64

export CustomEnum, CustomVal, define_enum

abstract type CustomEnum end
abstract type CustomVal <: Integer end

function build_items(val_type::Any, items:: Array{Symbol, 1})
    OrderedDict{Symbol, val_type}([(sym, val_type(sym, i)) for (i, sym) in enumerate(items)])
end

function show(io::IO, cv::CustomVal)
    print(io, "$(typeof(cv))::$(cv.name) = $(cv.val)")
end

function Base.Int64(val::CustomVal)
    val.val
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

function rand(config::Config, enum::CustomEnum)
    items = getfield(enum, :items)
    
    rand(config.rng, items) #returns a Pair{Symbol, Int64}
end

length(enum::CustomEnum) = length(getfield(enum, :items))

#iterate(enum::CustomEnum) = Base.iterate(getfield(enum, :items))
#iterate(enum::CustomEnum, state) = Base.iterate(getfield(enum, :items), state)
function iterate(enum::CustomEnum)
    vals = collect(values(getfield(enum, :items)))
    Base.iterate(vals)
end

function iterate(enum::CustomEnum, state)
    vals = collect(values(getfield(enum, :items)))

    Base.iterate(vals, state)
end

struct_templates = [
    mt"""
struct {{name}}Val <: CustomEnumMod.CustomVal
  name::Symbol
  val::Int64
end
""",
    mt"""
struct {{name}}s <: CustomEnumMod.CustomEnum
  items::OrderedDict{Symbol, {{name}}Val}

  function {{name}}s()
    new(CustomEnumMod.build_items({{name}}Val, {{values}}))
  end
end
"""
]

#note: regardless of where this function is called from, the enum will be defined within the CustomEnumMod module.
#The calling module can then create a local instance of it.
function define_enum(name::Symbol, values::Array{Symbol, 1})
    dict = Dict(
        "name" => name,
        "values" => values
    )

    code_segs = map(t -> Mustache.render(t, dict), struct_templates)
    map(eval ∘ Meta.parse, code_segs)
end

#----------------
# struct RegSiteVal <: CustomVal
#     name::Symbol
#     val::Int64
# end

# #note: these enum values will be used as indices into the corresponding arrays, so they should start at 1
# #accepts, regulates
# struct RegSites <: CustomEnum
#     items::OrderedDict{Symbol, RegSiteVal}
    
#     function RegSites()
#         new(build_items(RegSiteVal, [:IntraIntra, :IntraInter, :InterIntra, :InterInter]))
#     end
# end

# struct ProdSiteVal <: CustomVal
#     name::Symbol
#     val::Int64
# end

# #produces
# struct ProdSites <: CustomEnum
#     items::OrderedDict{Symbol, ProdSiteVal}

#     function ProdSites()
#         new(build_items(ProdSiteVal, [:Intra, :Inter]))
#     end
# end

# struct ProteinTypeVal <: CustomVal
#     name::Symbol
#     val::Int64
# end

# struct ProteinTypes <: CustomEnum
#     items::OrderedDict{Symbol, ProteinTypeVal}

#     function ProdSites()
#         new(build_items(ProteinTypeVal, [:Activate, :Inhibit]))
#     end
# end

# struct ProteinTargetVal <: CustomVal
#     name::Symbol
#     val::Int64
# end

# struct ProteinTargets <: CustomEnum
#     items::OrderedDict{Symbol, ProteinTargetVal}

#     function ProdSites()
#         new(build_items(ProteinTargetVal, [:Intra, :Inter]))
#     end
# end

# struct ProteinRegActionVal <: CustomVal
#     name::Symbol
#     val::Int64
# end

# struct ProteinRegActions <: CustomEnum
#     items::OrderedDict{Symbol, ProteinRegActionVal}

#     function ProdSites()
#         new(build_items(ProteinRegActionVal, [:Intra, :Inter]))
#     end
# end

# struct ProteinAppActionVal <: CustomVal
#     name::Symbol
#     val::Int64
# end

# struct ProteinAppActions <: CustomEnum
#     items::OrderedDict{Symbol, ProteinAppActionVal}

#     function ProdSites()
#         new(build_items(ProteinAppActionVal, [:Intra, :Inter]))
#     end
# end
    
end
