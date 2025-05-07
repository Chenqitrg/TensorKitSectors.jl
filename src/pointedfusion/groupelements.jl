# ===========================================
# Abstract Group and Group Element Definitions
# ===========================================

# include("../groups.jl")
# include("display.jl")



struct GroupElement{G<:Group}
    value::Any

    function GroupElement{G}(value) where {G<:ℤ}
        N = G.parameters[1]
        if N == Inf
            new{G}(Int(value))
        else
            new{G}(Int(mod(value, N)))
        end
    end
    function GroupElement{G}(values...) where {G<:D}
        s, r = values
        N = G.parameters[1]
        new((mod(s, 2), mod(r, N)))
    end
    function GroupElement{ProductGroup{Gs}}(value...) where {Gs<:GroupTuple}
        new{ProductGroup{Gs}}(value)
    end
    function GroupElement{G}(value) where {G<:CohomologyGroup}
        new{G}(value)
    end
end

function GroupElement{ProductGroup{Gs}}(value) where {Gs<:GroupTuple}
    return GroupElement{ProductGroup{Gs}}(value...)
end


   
#    @show  GroupElement{ℤ{3}}(4)
#    @show GroupElement{D{3}}(2,1)
# ===========================================
# Group Operations
# =========================================== 

function elements(::Type{ℤ{N}}) where {N}
    return ntuple(i -> GroupElement{ℤ{N}}(i - 1), N)
end

# elements(ℤ{3})
function elements(::Type{D{N}}) where {N}
    rotations = ntuple(i -> GroupElement{D{N}}(0, i - 1), N)  # (e, r, r², ...)
    reflections = ntuple(i -> GroupElement{D{N}}(1, i - 1), N)  # (s, sr, sr², ...)
    return (rotations..., reflections...)  # Joining the two tuples
end

# elements(D{3})
function elements(::Type{ProductGroup{Gs}}) where {Gs<:GroupTuple}
    group_elements = map(elements, Gs.parameters)
    cartesian_product = collect(Iterators.product(group_elements...))
    return map(GroupElement{ProductGroup{Gs}},cartesian_product)
end


# @show GroupElement{ℤ{3}×ℤ{3}}(GroupElement{ℤ₃}(0), GroupElement{ℤ₃}(1))
# display(elements(ℤ{3}×D{3}))

function identity_element(::Type{ℤ{N}}) where {N}
    return GroupElement{ℤ{N}}(0)
end
function identity_element(::Type{D{N}}) where {N}
    return GroupElement{D{N}}(0,0)
end
function identity_element(::Type{ProductGroup{Gs}}) where {Gs<:GroupTuple}
    groups = Gs.parameters
    return GroupElement{ProductGroup{Gs}}((identity_element(g) for g in groups), Gs)
end


function inverse(x::GroupElement{ℤ{N}}) where {N}
    return GroupElement{ℤ{N}}(-x.value)
end
function inverse(x::GroupElement{D{N}}) where {N}
    s, r = x.value
    return GroupElement{D{N}}(-s, (-1)^(s + 1) * r)
end
function inverse(x::GroupElement{ProductGroup{Gs}}) where {Gs<:GroupTuple}
    groups = Gs.parameters
    inverse_elements = map(x -> inverse(x), x.value)
    return GroupElement{ProductGroup{Gs}}(inverse_elements)
end

# inverse(GroupElement{ℤ₃}(1))
# inverse(GroupElement{D₃}(1,2))
# inverse(GroupElement{ℤ₃×D₃}(GroupElement{ℤ₃}(1), GroupElement{D₃}(0,1)))

function Base.:*(x::GroupElement{ℤ{N}}, y::GroupElement{ℤ{N}}) where {N}
    return GroupElement{ℤ{N}}(mod(x.value + y.value, N))
end
function Base.:*(x::GroupElement{D{N}}, y::GroupElement{D{N}}) where {N}
    s1, r1 = x.value
    s2, r2 = y.value
    return GroupElement{D{N}}(mod(s1 + s2, 2), mod((-1)^s2 * r1 + r2, N))
end
function Base.:*(x::GroupElement{ProductGroup{Gs}}, y::GroupElement{ProductGroup{Gs}}) where {Gs<:GroupTuple}
    newelement = ()
    for (x1, y1) in zip(x.value, y.value)
        newelement = (newelement..., x1 * y1)
    end
    return GroupElement{ProductGroup{Gs}}(newelement...)
end

# @show GroupElement{ℤ₃}(2) * GroupElement{ℤ₃}(2) * GroupElement{ℤ₃}(1)
# GroupElement{D₃}(1,2) * GroupElement{D₃}(1,1) * GroupElement{D₃}(0,1)
# GroupElement{ℤ₃×D₃}(GroupElement{ℤ₃}(1), GroupElement{D₃}(1,2)) * GroupElement{ℤ₃×D₃}(GroupElement{ℤ₃}(1), GroupElement{D₃}(0,1)) * GroupElement{ℤ₃×D₃}(GroupElement{ℤ₃}(1), GroupElement{D₃}(0,1))

Base.getindex(::Type{ℤ{N}}, i::Int) where {N} = GroupElement{ℤ{N}}(i)
Base.getindex(::Type{D{N}}, i::Int) where {N} = GroupElement{D{N}}(i÷N, i%N)
function Base.getindex(::Type{ProductGroup{Gs}}, i::Int) where {Gs<:GroupTuple}
    groups = Gs.parameters
    elementtuple = ()
    for n in length(groups):-1:1
        group_tail = groups[n]
        order_tail = order(group_tail)
        if n>1 && order_tail == Inf
            throw(ArgumentError("The infinite group $group_tail must be the first one."))
        end
        ind_tail =  Int(i%order_tail)
        i = i÷order_tail
        elementtuple = (group_tail[ind_tail], elementtuple...)
    end
    return GroupElement{ProductGroup{Gs}}(elementtuple...)
end

findindex(g::GroupElement{ℤ{N}}) where {N} = g.value
findindex(g::GroupElement{D{N}}) where {N} = g.value[1] * N + g.value[2]