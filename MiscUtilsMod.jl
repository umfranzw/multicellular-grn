module MiscUtilsMod

using Printf

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

function get_time_str(elapsed_sec::Float64)
    #credit to https://stackoverflow.com/a/52360730
    minutes, seconds = fldmod(elapsed_sec, 60)
    hours, minutes = fldmod(minutes, 60)

    if hours > 0
        result = @sprintf("%02d:%02d:%05.2f", Int64(hours), Int64(minutes), seconds)
    elseif minutes > 0
        result = @sprintf("%02d:%05.2f", Int64(minutes), seconds)
    else
        result = @sprintf("%0.2f sec", seconds)
    end

    result
end

end
