module MiscUtilsMod

export iprint, iprintln

const indent_spaces = 2

function num_enum_vals(enum::Any)
    length(instances(enum))
end

function enum_val_to_str(val::Any)
    buf = IOBuffer()
    show(buf, val)
    seek(buf, 0)

    read(buf, String)
end

function digits_needed(n::Int64)
    needed = 1
    if n > 0
        needed = max(Int64(ceil(log10(n))), 1)
    end

    needed
end

#indented print
function iprint(io::IO, obj::Any, indent_level::Int64=0)
    print(io, repeat(' ', indent_level * indent_spaces))
    print(io, obj)
end

#indented println
function iprintln(io::IO, obj::Any, indent_level::Int64=0)
    iprint(io, obj, indent_level)
    println(io, "")
end

end
