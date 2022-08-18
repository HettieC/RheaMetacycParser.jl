"""
$(TYPEDEF)

A structure for representing a single reaction taken from the rhea-reactions.txt database `id`
"""
struct RheaReaction
    id::String
    def::String
    stoich::Dict{String,Float64}
    EC::Vector{String}
end

"""
$(TYPEDEF)

A structure for representing a single reaction taken from the BiGG database, used in order to 
give the Rhea reactions a human-readable name.
"""
struct BiGGReaction
    id::String
    name::String
    databaselinks::Dict{String,Vector{String}}
end

"""
    struct MetacycPathway
        id::String
        name::String
        reactions::Vector{String}
    end

A structure for representing a single pathway from metacyc pathways.dat taken from the downloaded 
flat files
"""
struct MetacycPathway
    id::String
    name::String
    reactions::Vector{Pair{String,String}}
end

"""
    struct MetacycReaction
        id::String
        databaselinks::Dict{String,Vector{String}}
        stoich::Dict{String,Float64}
    end
Structure for representing a single reactions from reactions.dat metacyc file.
"""
struct MetacycReaction
    id::String
    databaselinks::Dict{String,Vector{String}}
    stoich::Dict{String,Float64}
end

"""
    struct SuperPathway
        id::String
        name::String
        pathways::Vector{String}
    end
Structure for representing a single super pathway from pathways.dat
"""
struct SuperPathway
    id::String
    name::String
    pathways::Vector{String}
end
