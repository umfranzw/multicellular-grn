module ListUtilsMod

export lisp_eval

function lisp_eval(code::String)
    stdout_buf = IOBuffer()
    stderr_buf = IOBuffer()
    cmd = `clisp -q -norc -x $(code)`
    run(pipeline(ignorestatus(cmd), stdout=stdout_buf, stderr=stderr_buf))

    stdout_str = rstrip(String(take!(stdout_buf)))
    stderr_str = rstrip(String(take!(stderr_buf)))

    (stdout=stdout_str, stderr=stderr_str)
end

end
