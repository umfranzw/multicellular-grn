using Printf

function print_matrix(M)
    rows, cols = size(M)
    for i in 1:rows
        for j in 1:cols
            @printf("%5.2f", M[i, j])
        end
        @printf("\n")
    end
    @printf("\n")
end

function get(mat, row, col)
    height, width = size(mat)
    if row < 1 || col < 1 || row > height || col > width
        return 0.0
    else
        return mat[row, col]
    end
end

const alpha = 0.1
const h = 0.75
const dt = 1
const rows = 9
const cols = 9
const timesteps = 15

Tk = zeros((rows, cols))
Tk[5, 5] = 1.0
Tk1 = copy(Tk)

print_matrix(Tk)

for k in 1:timesteps
    for i in 1:rows
        for j in 1:cols
            global Tk1[i, j] = (1 - 4 * dt * alpha / h^2) * Tk[i, j] + dt * alpha * ((get(Tk, i, j-1) + get(Tk, i-1, j) + get(Tk, i+1, j) + get(Tk, i, j+1)) / h^2)
        end
    end
    global Tk = copy(Tk1)
    print_matrix(Tk)
end




        
