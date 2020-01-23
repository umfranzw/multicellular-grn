module GraphVizMod

function gen_graph(dot_code::String, filename::String)
    stdout_buf = IOBuffer()
    stdin_buf = IOBuffer()
    write(stdin_buf, dot_code)
    seek(stdin_buf, 0)
    cmd = `dot -Tpng`
    run(pipeline(ignorestatus(cmd), stdin=stdin_buf, stdout=stdout_buf))
    seek(stdout_buf, 0)
    
    file = open(filename, "w")
    write(file, stdout_buf)
    close(file)
end

end
