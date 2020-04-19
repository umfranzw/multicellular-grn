module GraphVizMod

function plot(dot_code::String)
    stdout_buf = IOBuffer()
    stdin_buf = IOBuffer()
    write(stdin_buf, dot_code)
    seek(stdin_buf, 0)
    cmd = `dot -Tpng`
    run(pipeline(ignorestatus(cmd), stdin=stdin_buf, stdout=stdout_buf))

    stdout_buf.data
end

end
