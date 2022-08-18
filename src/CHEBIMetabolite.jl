""" 
    make_metabolites(cobrexa_reactions_dict::Dict{String,Reaction},path_to_compounds::String)
Make metabolite objects for the ChEBI metabolites that are involved in reactions, using the ChEBI 
    'compounds.tsv' file. This requires that the cobrexa_reactions_dict has already been made.
"""
function make_metabolites(cobrexa_reactions_dict::Dict{String,Reaction},path_to_compounds::String)
    metabolite_ids = Vector{String}()
    reactions_list = [r for r in values(cobrexa_reactions_dict)]
    for r in reactions_list
        append!(metabolite_ids,String.(keys(r.metabolites)))
    end
    unique!(metabolite_ids)
    id_name_dict = Dict{String,String}()
    metabolites_dict = Dict{String,Metabolite}()
    open(path_to_compounds) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline=false;continue)
            data = split(ln,"\t";limit=7)
            id = string("CHEBI:",data[1])
            id ∈ metabolite_ids
            name = (data[6]=="null" ? nothing : data[6])
            id_name_dict[data[1]]= (data[6]=="null" ? "" : data[6])
            metabolites_dict[id] = Metabolite(string(id);name=name)            
        end
    end
    return Dict(x => metabolites_dict[x] for x in metabolite_ids)
end

"""
    add_chemical_data(metabolites_dict::Dict{String,Metabolite},path_to_chemical_data::String)
Add the chemical data from `chemical_data.tsv` ChEBI file to the metabolites created in `make_metabolites(..)`
"""
function add_chemical_data(metabolites_dict::Dict{String,Metabolite},path_to_chemical_data::String)
    open(path_to_chemical_data) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline=false;continue)
            data = string.(split(ln,"\t"))
            id = string("CHEBI:",data[2])
            if haskey(metabolites_dict,id)
                dtype = data[4]
                v = metabolites_dict[id]
                if dtype == "FORMULA"
                    v.formula = data[5]
                elseif dtype == "CHARGE"
                    v.charge = parse(Int64, data[5])
                end
            end
        end
    end
    return metabolites_dict
end

