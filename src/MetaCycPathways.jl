"""
$(FUNCTIONNAME)(path_to_metacyc_pathways::String)
Make a dict of pathways from metacyc pathways.tsv, including whether or not they are a super pathway.
"""
function parse_metacyc_pathways(path_to_metacyc_pathways::String)
    pathways_dict = Dict{String,Pair{Vector{String},Vector{String}}}()
    open(path_to_metacyc_pathways) do io
        id = "" 
        type = "normal"
        name = "" 
        reactions_or_pways = String[]
        for ln in eachline(io)
            data = split(ln, " - ")
            if data[1] == "#"
                continue
            elseif ln == "//"
                pathways_dict[id] = Pair([name,type],reactions_or_pways)
                id = "" 
                type = "normal"
                name = "" 
                reactions_or_pways = String[]
            elseif data[1] == "UNIQUE-ID"
                id = data[2]
            elseif data[1] == "TYPES" && data[2] == "Super-Pathways" 
                type = "superpathway"
            elseif data[1] == "COMMON-NAME"
                name = data[2]
            elseif data[1] == "REACTION-LIST"
                append!(reactions_or_pways,[data[2]])
            end
        end
    end
    return pathways_dict
end


"""
$(FUNCTIONNAME)(r_or_p::String)
If r_or_p is itself a reaction, return a list containing this reaction.
Otherwise, if r_or_p is a pathway, collect all reactions involved this pathway. Use recursion to check
and collect all reactions or pathways that are part of r_or_p.
"""
function collect_pathway_reactions(r_or_p::String)
    reactions = Vector{Pair{String,String}}()
    if haskey(metacyc_rhea_dict,r_or_p)
        append!(reactions,[Pair(r_or_p,metacyc_rhea_dict[r_or_p])])
    elseif haskey(pathways_dict,r_or_p)
        for rs_or_ps in pathways_dict[r_or_p].second 
            append!(reactions,collect_pathway_reactions(rs_or_ps))
        end
    else
        append!(reactions,[Pair(r_or_p,"")])
    end
    return reactions
end

"""
$(FUNCTIONNAME)(path_to_metacyc_pathways::String,metacyc_rhea_dict::Dict{String,String})
Make a list of MetacycPathway objects from the pathways_dict, with pathway.reactions being the list of reactions involved
calculated using collect_pathway_reactions(r_or_p::String).
"""
function make_metacyc_pathways(path_to_metacyc_pathways::String,metacyc_rhea_dict::Dict{String,String})
    pathways_dict = parse_metacyc_pathways(path_to_metacyc_pathways)
    metacyc_pathways = Dict{String,MetacycPathway}()
    for (k,v) in pathways_dict
        reactions = Vector{Pair{String,String}}()
        
        # first deal with the "normal" pathways
        if v.first[2] == "normal"
            metacyc_reactions = v.second
            for r in metacyc_reactions
                if haskey(metacyc_rhea_dict,r)
                    append!(reactions,[Pair(r,metacyc_rhea_dict[r])])
                else
                    append!(reactions,[Pair(r,"")])
                end
            end
            metacyc_pathways[k] = MetacycPathway(k,v.first[1],reactions)

            # then deal with "superpathway"
        elseif v.first[2] == "superpathway"
            for r_or_p in v.second
                append!(reactions,collect_pathway_reactions(r_or_p))
            end
            metacyc_pathways[k] = MetacycPathway(k,v.first[1],reactions)
        end
    end
    return metacyc_pathways
end


"""
$(FUNCTIONNAME)(metacyc_pathways)
Make a dict of metacyc pathways for which all reactions are in rhea.
"""
function select_all_rhea_pathways(metacyc_pathways)
    non_rhea_pathways = Vector{MetacycPathway}()
    for (id,pathway) in metacyc_pathways
        rhea_reactions = String[]
        for r in pathway.reactions
            append!(rhea_reactions, [r.second])
        end
        unique!(rhea_reactions)
        if any(rhea_reactions.=="")
            push!(non_rhea_pathways,pathway)
        end
    end
    rhea_pathway_ids = [p for p in keys(metacyc_pathways) if p ∉ [k.id for k in non_rhea_pathways]]
    return Dict(p => metacyc_pathways[p] for p in rhea_pathway_ids)
end 

"""
$(FUNCTIONNAME)
Make the dict of rhea pathways in one function.
    NOTE: update MetaCycRheaMapping.jl so that I can make metacyc_rhea_dict in one function
"""
function make_rhea_pathways(path_to_metacyc_pathways::String,metacyc_rhea_dict::Dict{String,String})
    metacyc_pathways = make_metacyc_pathways(path_to_metacyc_pathways,metacyc_rhea_dict)
    return select_all_rhea_pathways(metacyc_pathways)
end