import Base.~, Base.⊻, Base.|, Base.&
import Base.Int64

hamming_dist(x::BitArray, y::BitArray) = sum(map(xor, x, y))
count_common_bits(x::BitArray, y::BitArray) = sum(map(==, x, y))
~(x::BitArray) = map(~, x)
⊻(x::BitArray, y::BitArray) = map(xor, x, y)
|(x::BitArray, y::BitArray) = map(|, x, y)
(&)(x::BitArray, y::BitArray) = map((&), x, y)

function Int64(bits::BitArray)
    result = 0
    place_val = 1
    for i in length(bits):-1:1
        result += bits[i] * place_val
        place_val <<= 1
    end

    result
end

function BitArray(n::Int64)
    result = BitArray{1}()
    while n != 0
        push!(result, n % 2)
        n >>= 1 #integer division
    end

    reverse!(result)
end

