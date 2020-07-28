function convertMathMode!(name::String)
    # expects markdown files as input
    @assert getFileEnding(name)=="md" "no markdown file provided"
    lines = Vector{String}
    open(name) do file
        lines = readlines(file)
    end
    rm(name)
    fout = open(name,"w")
    first = true
    for (i,line) in enumerate(lines)
        if line=="\$\$"
            if first==true
                lines[i] = "```math"
                first=false
            else
                lines[i] = "```"
                first=true
            end
        end
        print(fout,lines[i]*"\n")
    end
    close(fout)
    return name
end

function convertJuliaMode!(name::String,envname::String)
    # expects markdown files as input
    @assert getFileEnding(name)=="md" "no markdown file provided"
    lines = Vector{String}
    open(name) do file
        lines = readlines(file)
    end
    rm(name)
    fout = open(name,"w")
    first = true
    for (i,line) in enumerate(lines)
        if line=="```julia"
            lines[i] = "```@example $envname"
        end
        print(fout,lines[i]*"\n")
    end
    close(fout)
    return name
end

function removeComments!(name::String)
    # expects julia file as input, i.e. name = "xxx.jl"
    @assert getFileEnding(name)=="jl" "no julia file provided"
    lines = Vector{String}
    open(name) do file
        lines = readlines(file)
    end
    rm(name)
    fout = open(name,"w")
    remove = zeros(Int64,0)
    keep = zeros(Int64,0)
    for (i,line) in enumerate(lines)
        # display(line)
        if length(line)>0 && Int(line[1])!=35
            # Int()==35 --> character is a hashtag
            print(fout,lines[i]*"\n")
        end
    end
    close(fout)
end

function addSetupAndMakeMarkdown(name::String,envname::String="mysetup")
    # expects julia file as input, i.e. name = "xxx.jl"
    @assert getFileEnding(name)=="jl" "no julia file provided"
    newname = removeFileEnding(name)
    newname = "$(newname)_setup.md"
    # s = split(name,"")
    # nameind = findfirst(x->x==".",s)
    lines = Vector{String}
    open(name) do file
        lines = readlines(file)
    end
    rm(name)
    fout = open(newname,"w")
    print(fout,"```@setup "*envname*"\n")
    [ print(fout,lines[i]*"\n") for i=1:length(lines) ]
    print(fout,"```\n")
    close(fout)
    return newname
end

function createPreamble(name::String)
    fname = "preamble.md"
    f = open(fname,"w")
    print(f,"----------------------------\n")
    print(f,"----------------------------\n")
    print(f,"----------------------------\n")
    print(f, "This file was automatically created from $name.\n")
    print(f,"----------------------------\n")
    print(f,"----------------------------\n")
    print(f,"----------------------------\n")
    close(f)
    return fname
end

"""
    notebook2markdown(fname::String)
Workflow to get markdown file from Julia notebook
"""
function notebook2markdown(fname::String)
    @assert getFileEnding(fname)=="ipynb" "no ipynb file provided"
    run(`jupyter-nbconvert --ClearOutputPreprocessor.enabled=True --inplace $fname --to markdown`) # ipynb to markdown
    # remove file ending (if any)
    name = removeFileEnding(fname)
    file_text = convertMathMode!("$name.md")
end



"""
    notebook2markdown_repl(fname::String;envname::String="mysetup")
Workflow to get REPL-based markdown file from Julia notebook
"""
function notebook2markdown_repl(fname::String;envname::String="mysetup")
    @assert getFileEnding(fname)=="ipynb" "no ipynb file provided"
    # remove file ending (if any)
    name = removeFileEnding(fname)
    # create ```setup mysetup ``` entry
    run(`jupyter-nbconvert $fname --ClearOutputPreprocessor.enabled=True --inplace`)
    run(`jupyter-nbconvert $fname --to python`) # ipynb to python
    run(`mv $name.py $name.jl`) # python to julia
    removeComments!("$name.jl")
    # file_preamble = createPreamble("$name.jl")
    file_setup = addSetupAndMakeMarkdown("$name.jl",envname)

    # convert $$...$$ to ```math...```
    run(`jupyter-nbconvert $fname --to markdown`) # ipynb to markdown
    file_text = convertMathMode!("$name.md")
    file_text = convertJuliaMode!(file_text,envname)
    # append the two files
    text = [ readlines(file_setup); readlines(file_text) ]
    rm(file_text)
    rm(file_setup)
    fout = open(name*".md","w")
    [ print(fout,line*"\n") for line in text]
    close(fout)
end

function removeFileEnding(name::String)
    s = split(name,"")
    dotindex = findfirst(x->x==".",s)
    dotindex==nothing ? (return name) : (return name[1:dotindex-1])
end

function getFileEnding(name::String)
    s = split(name,"")
    dotindex = findfirst(x->x==".",s)
    dotindex==nothing ? (return name) : (return lowercase(name[dotindex+1:end]))
end

# Workflow to get REPL-executed markdown file from Julia notebook
    # 1. Create setup
    #     jupyter-nbconvert file --to python
    #     mv file.py file.jl
    #     removeComments!()
    #     add ```@setup```
    # 2. Add math-converted markdown
    #     jupyter-nbconvert file --to markdown
    #     convert2math
    # 3. concatenate results from 1. and 2.
