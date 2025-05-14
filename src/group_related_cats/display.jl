


function Base.show(io::IO, G::Type{<:ProductGroup})
    T = G.parameters[1]
    groups = T.parameters
    for (i, G) in enumerate(groups)
        if i != 1
            print(io, "×")
        end
        print(io, G)
    end
end

# function Base.show(io::IO, G::Type{<:CohomologyGroup})
#     T = G.parameters
#     print(io, "H", superscript(T[1]), "(", T[2], " , ", T[3], ")")
# end

function subscript(n::Integer)
    subs = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
                '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉')
    return join([subs[c] for c in string(n)])
end

function to_subscript(string::String)
    subscript_map = Dict(
        'a' => 'ₐ',
        'e' => 'ₑ',
        'h' => 'ₕ',
        'i' => 'ᵢ',
        'j' => 'ⱼ',
        'k' => 'ₖ',
        'l' => 'ₗ',
        'm' => 'ₘ',
        'n' => 'ₙ',
        'o' => 'ₒ',
        'p' => 'ₚ',
        'r' => 'ᵣ',
        's' => 'ₛ',
        't' => 'ₜ',
        'u' => 'ᵤ',
        'v' => 'ᵥ',
        'x' => 'ₓ'
    )
    return join([subscript_map[c] for c in string])
end

function superscript(n::Integer)
    subs = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹', '-' => '⁻')
    return join([subs[c] for c in string(n)])
end


function Base.show(io::IO, x::ℤ{N}) where {N}
    if x.a == 0
        print(io, "e")
    elseif x.a == 1
        print(io, "a")
    else
        print(io, "a", superscript(x.value))
    end
end

function Base.show(io::IO, x::D{N}) where {N}
    s, r = x.s, x.r
    if s == 0 && r == 0
        print(io, "e")
    elseif s == 0 && r == 1
        print(io, "r")
    elseif s == 0 && r > 0
        print(io, "r", superscript(r))
    elseif s == 1 && r == 0
        print(io, "s")
    elseif s == 1 && r == 1
        print(io, "sr")
    else
        print(io, "sr", superscript(r))
    end
end
function Base.show(io::IO, x::ProductGroup{Gs}) where {Gs<:GroupTuple}
    print(io, "(", join(x.components, ", "), ")")
end

function Base.show(io::IO, x::VecGωIrr{G,ω}) where {G<:Group,ω}
    print(io, x.g)  
end

function Base.show(io::IO, x::TYIrr{A,χ,ϵ}) where {A<:Group, χ,ϵ}
    print(io, x.value)  
end