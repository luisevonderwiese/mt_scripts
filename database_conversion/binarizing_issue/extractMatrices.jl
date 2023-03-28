using Pkg
Pkg.activate(".")
Pkg.instantiate()

##

using DataFrames
using CSV
using Pipe

##

pth = "raw/lexibank"

datasets = readdir(pth)

filter!(x -> x != "README.md", datasets)



##

cog_dbs = []
for db in datasets
    if isfile(joinpath(pth, db, "cldf", "cognates.csv"))
        cognates = CSV.File(joinpath(pth, db, "cldf", "cognates.csv")) |>
            DataFrame
        if "expert" ∈ cognates.Cognate_Detection_Method
            push!(cog_dbs, db)
        end
    end
end

##

function db2charMtx(db)
    cognates = @pipe CSV.File(joinpath(pth, db, "cldf", "cognates.csv")) |>
        DataFrame |>
        filter(x -> x.Cognate_Detection_Method == "expert", _)
    forms = CSV.File(joinpath(pth, db, "cldf", "forms.csv")) |> DataFrame
    d = outerjoin(
        forms,
        cognates,
        on=:ID => :Form_ID,
        makeunique=true,
    )
    df = @pipe d |>
               select(_, [:Language_ID, :Parameter_ID, :Cognateset_ID]) |>
               dropmissing |>
               unique |>
               rename!(_, [:language, :concept, :cc])

    concepts = df.concept |> unique
    dfs_ = []
    for c in concepts
        dfc = @pipe df |>
                    filter(x -> x.concept == c, _) |>
                    unique |>
                    unstack(_, :cc, :language, :concept) |>
                    _[:, 2:end] |>
                    1 .- ismissing.(_)
        push!(dfs_, dfc)
    end
    vcat(dfs_..., cols=:union)
end
##

charMatrices = db2charMtx.(cog_dbs)

##

function writeCharMtx(cm, fn)
    function printChar(ch)
        ismissing(ch) ? "-" : string(ch)
    end

    pad = maximum(length.(names(cm))) + 5

    nex = """#Nexus
    BEGIN DATA;
    DIMENSIONS ntax=$(size(cm, 2)) nchar = $(size(cm, 1));
    FORMAT DATATYPE=Restriction GAP=? MISSING=- interleave=no;
    MATRIX

    """
    for l in names(cm)
        nex *= rpad(l, pad) * join(printChar.(cm[:,l]))*"\n"
    end
    nex *= ";\nEND"
    open(fn, "w") do file
        write(file, nex)
    end
end

##

pth = "julia"

try
    mkdir(pth)
catch e
end

for (db, cm) in zip(cog_dbs, charMatrices)
    writeCharMtx(cm, joinpath(pth, db*".nex"))
end
