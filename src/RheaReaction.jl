""" 
    rhea_reactions(path_to_rhea_reactions::String) 
Create RheaReaction objects out of each entry in the "rhea-reactions.txt" file. This file uses spaces instead of tabular delimiters.
"""
function make_rhea_reactions(path_to_rhea_reactions::String)
    rhea_reactions_list = RheaReaction[]
    open(path_to_rhea_reactions) do io
        id = "" 
        rs = Dict{String,Float64}()
        def = "" 
        EC = String[] 
        # regex expressions to match for the metabolites and the enzymes
        rcoeff = zero_or_more(DIGIT) * maybe(".") * zero_or_more(DIGIT) * exactly(1, " ")
        rmet = exactly(1, "CHEBI:") * one_or_more(DIGIT)
        rmiddle = exactly(1," ") * maybe("<") * exactly(1,"=") * maybe(">") * exactly(1," ")
        renzyme = one_or_more(DIGIT) * exactly(1,".") * one_or_more(DIGIT) * exactly(1,".") * one_or_more(DIGIT) * exactly(1,".") * one_or_more(DIGIT)
        for line in eachline(io)
            if line == "///"
                push!(rhea_reactions_list, RheaReaction(id,def,rs,EC))
                id = ""   
                rs = Dict{String,Float64}()
                def = ""
                EC = String[]
            else
                splitline = split(line," ";limit=2)
                k = first(splitline)
                if first(splitline) == "ENTRY"
                    id = String(strip(last(split(last(splitline), ":"))))
    
                elseif first(splitline) == "DEFINITION"
                    def = String(last(splitline))
    
                elseif first(splitline) == "EQUATION"
                    rxn = last(splitline)
                    # split the equation up using readable regex to get the stoichiometry
                    delim = match(rmiddle, rxn)
                    subs_prods = split(rxn, delim.match)
                    subs = split(first(subs_prods), " + ")
                    prods = split(last(subs_prods), " + ")
    
                    rs = Dict{String, Float64}()
    
                    for sub in subs
                        m = match(rcoeff, sub)
                        c = parse(Float64, first(split(isnothing(m) || m.match == " " ? "1" : m.match)))                    
                        metmatch = match(rmet, sub)
    
                        if isnothing(metmatch)==false
                            met = metmatch.match
                            rs[met] = -c                      
                        end
                    end
    
                    for prod in prods
                        m = match(rcoeff, prod)
                        c = parse(Float64, first(split(isnothing(m) || m.match == " " ? "1" : m.match)))                    
                        metmatch = match(rmet, prod)
    
                        if isnothing(metmatch) == false
                            met = metmatch.match
                            rs[met] = c                        
                        end
                    end
    
                elseif first(splitline) == "ENZYME" 
                    enz_split = split(last(splitline)," ")
                    for x in enz_split[2:end]
                        enz = match(renzyme,x)
                        if isnothing(enz) == false
                        append!(EC,[String(enz.match)])
                        end
                    end
                elseif first(splitline) == "" # if there are too many enzymes associated with the reaction then they are 
                    # displayed over multiple lines in the "rhea-reactions.txt" file
                    enz_split = split(strip(last(splitline))," ")
                    for x in enz_split
                        enz = match(renzyme,x)
                        if isnothing(enz) == false
                            append!(EC,[String(enz.match)])
                        end
                    end
                else
                    @warn("Unhandled key: ", line)
                end   
            end
        end
    end
    return out = rhea_reactions_list
end

"""
    bigg_reactions(path_to_bigg_reactions::String)
Creat BiGGReaction objects of every reaction in the bigg_models_reactions.txt file taken from http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt
which is a tab delimitted file with `bigg_id\tname\treaction_string\tmodel_list\tdatabase_links\told_bigg_ids` as the header. 
"""
function make_bigg_reactions(path_to_bigg_reactions::String)
    bigg_reactions_list = BiGGReaction[]
    open(path_to_bigg_reactions) do io
        firstline = true
        for line in eachline(io)
            # the first line gives the headers for the columns, so do not parse this
            firstline && (firstline=false;continue) 
            databaselinks = Dict{String,Vector{String}}()
            id = ""
            name = ""
            # split each line on "\t" to get the column entries
            splitline = split(line,"\t")
            id = splitline[1]
            name = splitline[2]
            dbvec = split(splitline[5],";") # split the database_links line on ";" to account
            # for multiple links per id
            for entry in dbvec
                if entry == "" # there is no database_link entry for this id
                    databaselinks = Dict()
                else
                    if entry[1] == " " # some links have a " " space in front of the database name
                        dbentry = split(strip(entry),": ";limit=2) # each entry is in the form 
                        # "db_name: http://db/entry" so split on only the first colon to then add the 
                        # db_name as the key to a dictionary of "databaselinks" and the link as the 
                        # dict value
                        db = strip(dbentry[1])
                        db_link = String(strip(dbentry[2]))
                        if db ∉ keys(databaselinks)
                            databaselinks[db] = [db_link]
                        else
                            append!(databaselinks[db],db_link)
                        end
                    else 
                        dbentry = split(entry,": ";limit=2)
                        db = strip(dbentry[1])
                        db_link = String(strip(dbentry[2]))
                        if db ∉ keys(databaselinks)
                            databaselinks[db] = [db_link]
                        else
                            append!(databaselinks[db],[db_link])
                        end
                    end
                end
            end
            push!(bigg_reactions_list, BiGGReaction(id,name,databaselinks))
        end
    end
    return out = bigg_reactions_list
end


""" 
    rhea_bigg_mapping(bigg_reactions_list::Vector{BiGGReaction})
Make a dictionary of `rhea_id => bigg_id`.
"""
function rhea_bigg_mapping(bigg_reactions_list::Vector{BiGGReaction})
    bigg_rhea_dict = Dict{String, Vector{String}}()
    rhea_bigg_dict = Dict{String, Vector{String}}()
    for r in bigg_reactions_list
        if "RHEA" ∈ keys(r.databaselinks)
            links = last.(rsplit.(r.databaselinks["RHEA"],"/"))
            for v in links
                rhea_bigg_dict[v] = [r.id,r.name]
            end
            bigg_rhea_dict[r.id] = links
        end
    end
    return out = rhea_bigg_dict
end

"""
    make_cobrexa_reactions(rhea_reactions_list::Vector{RheaReaction},rhea_bigg_dict::Dict{String,String})
Create COBREXA reactions from the list of rhea reactions. Rhea IDs are used as the `id` and BiGG IDs are used as the `name`.
"""
function make_cobrexa_reactions(rhea_reactions_list::Vector{RheaReaction},rhea_bigg_dict::Dict{String,Vector{String}})
    cobrexa_reactions_list = Reaction[]
    for r in rhea_reactions_list
        if r.id in keys(rhea_bigg_dict)
            append!(cobrexa_reactions_list, [Reaction(r.id; name=rhea_bigg_dict[r.id][1], metabolites = r.stoich, notes = Dict("BiGG name"=>[rhea_bigg_dict[r.id][2]]))])
        else 
            append!(cobrexa_reactions_list, [Reaction(r.id; metabolites = r.stoich)])
        end
    end
    return out = Dict(x.id => x for x in cobrexa_reactions_list)
end

"""
    make_reactions(rhea_path::String,bigg_path::String,compounds_path::String,chemical_data_path::String)
Parse the Rhea and BiGG files to create COBREXA reactions from Rhea reactions, with BiGG ids as the reaction.name.
"""
function make_reactions(rhea_path::String,bigg_path::String)
    rhea_reactions_list = make_rhea_reactions(rhea_path)
    bigg_reactions_list = make_bigg_reactions(bigg_path)
    rhea_bigg_dict = rhea_bigg_mapping(bigg_reactions_list)
    cobrexa_reactions_dict = make_cobrexa_reactions(rhea_reactions_list,rhea_bigg_dict)
    return cobrexa_reactions_dict
end