"""
$(FUNCTIONNAME)(path_to_rhea_directions::String)
Make list of lists of the reaction directions, using the rhea-directions.tsv file
"""
function parse_rhea_directions(path_to_rhea_directions::String)
    directions = Vector{Vector{String}}()
    open(path_to_rhea_directions) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline=false;continue) 
            dirs = split(ln, "\t")
            push!(directions,[String(x) for x in dirs])
        end
    end 
    return directions
end

"""
$(FUNCTIONNAME)(path_to_metacyc_reactions::String)
Make a dict of metacyc reaction against rhea database links, using the metacyc reactions.dat file
Also make list of MetacycReaction objects so that I can compare them to RheaReaction objects and 
match their stoichiometries.
"""
function parse_metacyc_reactions(path_to_metacyc_reactions::String)
    metacyc_rhea_dict = Dict{String,Vector{String}}()
    metacyc_reactions = MetacycReaction[]
    open(path_to_metacyc_reactions) do io 
        id = "" 
        rhea = String[]
        stoich = Dict{String,Float64}()
        databaselinks = Dict{String,Vector{String}}()
        compounds = String[]
        for ln in eachline(io)
            data = split(ln," - ")
            if data[1] == "#"
                continue
            elseif ln == "//"
                metacyc_rhea_dict[id] = rhea
                rhea = String[]
                push!(metacyc_reactions,MetacycReaction(id,databaselinks,stoich))
                id = "" 
                stoich = Dict{String,Float64}()
                databaselinks = Dict{String,Vector{String}}()
                compounds = String[]
            elseif data[1] == "UNIQUE-ID"
                id = data[2]
            elseif data[1] == "DBLINKS" && startswith(data[2],"(RHEA ")
                append!(rhea,[string(chop(split(data[2]," ")[2];head=1,tail=1))])
                databaselinks["Rhea"] = rhea
            elseif data[1] == "LEFT"
                stoich[data[2]] = -1.0
                append!(compounds,[data[2]])
            elseif data[1] == "^COEFFICIENT"
                last_cmpd = last(compounds)
                right_or_left = stoich[last_cmpd]
                if !isnothing(tryparse(Float64,data[2]))
                    stoich[last_cmpd] = tryparse(Float64,data[2])*right_or_left
                else
                    id = "" 
                    stoich = Dict{String,Float64}()
                    databaselinks = Dict{String,Vector{String}}()
                end
            elseif data[1] == "RIGHT"
                stoich[data[2]] = 1.0
                append!(compounds,[data[2]])
            end
        end
    end
    return metacyc_rhea_dict, metacyc_reactions
end

"""
$(FUNCTIONNAME)(path_to_metacyc_compounds::String)
Make a dict of metacyc compnd id againts chebi id, using the metacyc compounds.dat file
"""
function map_metacyc_to_chebi(path_to_metacyc_compounds::String)
    metacyc_chebi_dict = Dict{String,String}()
    open(path_to_metacyc_compounds) do io
        meta = "" 
        chebi = "" 
        for ln in eachline(io)
            data = split(ln, " - ")
            if data[1] == "#"
                continue
            elseif data[1] == "//"
                metacyc_chebi_dict[meta] = chebi 
                meta = "" 
                chebi = "" 
            elseif data[1] == "UNIQUE-ID"
                meta = data[2]
            elseif data[1] == "DBLINKS" && startswith(data[2],"(CHEBI ")
                chebi = string(chop(split(data[2]," ")[2];head=1,tail=1))
            end
        end
    end
    return metacyc_chebi_dict
end

"""
$(FUNCTIONNAME)(metacyc_reactions::Vector{MetacycReaction},metacyc_chebi_dict::Dict{String,String},
metacyc_rhea_dict::Dict{String,Vector{String}},rhea_reactions_list::Vector{RheaReaction})

Match metacyc reactions to rhea reactions based on their stoichiometries. We only use the metacyc reactions
whose compounds have chebi ids.
"""
function match_metacyc_rhea_stoich(metacyc_reactions::Vector{MetacycReaction},metacyc_chebi_dict::Dict{String,String},
        metacyc_rhea_dict::Dict{String,Vector{String}},rhea_reactions_list::Vector{RheaReaction})
    metacyc_reactions_with_chebi_ids = Dict{String,MetacycReaction}()
    for reaction in metacyc_reactions
        compounds = keys(reaction.stoich)
        chebi_ids = [haskey(metacyc_chebi_dict,c) ? metacyc_chebi_dict[c] : "" for c in compounds]
        if any(chebi_ids.=="")
            continue
        else
            new_stoich = Dict{String,Float64}()
            for (k,v) in reaction.stoich
                new_stoich[string("CHEBI:",metacyc_chebi_dict[k])] = v 
            end
            metacyc_reactions_with_chebi_ids[reaction.id] = MetacycReaction(reaction.id,reaction.databaselinks,new_stoich)
            new_stoich = Dict{String,Float64}()
        end
    end
    #Now match stoichiometries and add matches to the metacyc_rhea_dict
    for (m_id,m_rxn) in metacyc_reactions_with_chebi_ids
        for r_rxn in rhea_reactions_list
            if m_rxn.stoich == r_rxn.stoich
                if haskey(metacyc_rhea_dict,m_id)
                    append!(metacyc_rhea_dict[m_id],[r_rxn.id])
                    unique!(metacyc_rhea_dict[m_id])
                else
                    metacyc_rhea_dict[m_id] = [r_rxn.id]
                end
            end
        end
    end
    return Dict(k=>v for (k,v) in metacyc_rhea_dict if v!=[])
end

"""
$(FUNCTIONNAME)(path_to_rhea2metacyc::String,metacyc_rhea_dict::Dict{String,Vector{String}})
Incorporate rhea2metacyc.tsv into the metacyc_rhea_dict
"""
function parse_rhea2metacyc(path_to_rhea2metacyc::String,metacyc_rhea_dict::Dict{String,Vector{String}})
    open(path_to_rhea2metacyc) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline=false;continue)
            data = split(ln,"\t")
            metacyc = data[4]
            rhea = [data[1],data[3]]
            if haskey(metacyc_rhea_dict,metacyc)
                append!(metacyc_rhea_dict[metacyc],rhea)
                unique!(metacyc_rhea_dict[metacyc])
            else
                metacyc_rhea_dict[metacyc] = unique!(rhea)
            end
            metacyc = "" 
            rhea = String[]
        end 
    end
    return metacyc_rhea_dict
end

"""
$(FUNCTIONNAME)(path_to_rhea_directions::String,path_to_metacyc_reactions::String,
path_to_metacyc_compounds::String,path_to_rhea2metacyc::String,path_to_rhea_reactions::String)
Map metacyc reactions to rhea reactions in one step.
Only output the non-directional master rhea reaction
"""
function map_metacyc_reactions_to_rhea(rhea_reactions_list::Vector{RheaReaction},path_to_metacyc_reactions::String,
        path_to_metacyc_compounds::String,path_to_rhea2metacyc::String,path_to_rhea_directions::String)
    directions = parse_rhea_directions(path_to_rhea_directions)
    first_metacyc_rhea_dict, metacyc_reactions = parse_metacyc_reactions(path_to_metacyc_reactions)
    metacyc_chebi_dict = map_metacyc_to_chebi(path_to_metacyc_compounds)
    metacyc_rhea_dict = match_metacyc_rhea_stoich(metacyc_reactions, metacyc_chebi_dict, first_metacyc_rhea_dict, rhea_reactions_list)
    metacyc_rhea_dict = parse_rhea2metacyc(path_to_rhea2metacyc, metacyc_rhea_dict)
    new_metacyc_rhea_dict = Dict{String,String}()
    for (k,v) in metacyc_rhea_dict
        for r in v
            for dir in directions 
                if r ∈ dir
                    new_metacyc_rhea_dict[k] = dir[1]
                    break
                end
            end
        end
    end
    return new_metacyc_rhea_dict
end