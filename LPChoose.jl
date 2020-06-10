#using DataFrames,LinearAlgebra,CSV,DelimitedFiles,SparseArrays,GLPK,JuMP,Statistics
#include("LPseq.jl")
#input="smalldata.txt"
#sol,hap =LPChoose(input,100,0.00,nsteps=5);

#same results as fixed_budget_LP with MAF=0 if hap frequency is calculated correctly in fixed_budget_LP
function LPChoose(hapblock,budget=100,MAF=0.0;nsteps=1) #budget is #of selected animals
    #Get incidence matrix
    A01,A012,freq  = convert_input_to_A(hapblock,MAF)
    animals  = collect(1:size(A01,2))
    nind     = length(animals)
    if budget == "unlimited"
        c=ones(Int64,N)
        model = Model(GLPK.Optimizer)
        @variable(model, select_this_animal[1:nind],Bin)
        @objective(model, Min, sum([select_this_animal[i] for i= 1:nind]))
        @constraint(model,con,A01 * select_this_animal .>= 1)
        @time optimize!(model)
        termination_status(model)#check if an optimal solution is found
        primal_status(model)
        println("The minimum nunber of selected animals is ", objective_value(model))
    else
        A01_temp   =copy(A01)
        #get importance of each individual
        budget_each_step = Int(budget/nsteps)
        importance = Matrix(freq'A01_temp)) #importance = vec(Matrix(freq'A01))
        for stepi = 1:nsteps
            model = Model(GLPK.Optimizer)
            @variable(model, select_this_animal[1:nind],Bin)
            @objective(model, Max, sum([importance[i]*select_this_animal[i] for i= 1:nind]))
            @constraint(model,A01 * select_this_animal .<= 2)
            @constraint(model,sum(select_this_animal) <= budget_each_step) #e.g.,selecte 2 animals at each step
            @time optimize!(model)
            termination_status(model)     #check if the solver found an optimal solution
            primal_status(model)          #the primal only inform that the primal solutions are feasible
            objective_value(model)        #final objective solution

            # solution
            sol     = JuMP.value.(select_this_animal)
            #which haplotypes are identified by selected animals
            haps_identified = (vec(sum(A01_temp[:,sol .== 1.0],dims=2)) .!= 0)
            A01_temp        = A01_temp[haps_identified,sol]
            freq            = freq[sol .== 1.0]#vec(mean(A01_temp,dims=2)/2)
            importance      = Matrix(freq'A01_temp))
            println("Number of animals selected at step", stepi, " is: ",sum(sol))
        end
        #need a line for which animal is selected
        println(sum(A01_temp)/sum(A01), " of the genome for the population is covered.")
    end
end

#the input file is haplotypes blocks
#hapblock is a file where each row representing one individual and
#each column representing one haplotype block, e.g.
#
#column 1: Individual Name
#column 2-3: Maternal and Paternal Haplotypes block 1
#column 4-5: Maternal and Paternal Haplotypes block 2
#
#ind1 hap1_1 hap1_1 hap2_1 hap2_1
#ind2 hap1_2 hap1_1 hap2_1 hap2_2
#ind3 hap1_1 hap1_3 hap2_2 hap2_3
#i.e.,
#1,1,1,1*,1*
#2,2,1,1*,2*
#3,1,3,2*,3*
#
#
#
#convert input file to "incidence matrix" of zeros/ones where each row represents
# one haplotype in the haplotype block,e.g, A'=
#ID,  hap1_1,hap1_2,hap1_3,hap2_1,hap2_2,hap2_3
#ind1,1,0,0,1,0,0
#ind2,1,1,0,1,1,0
#ind3,1,0,1,0,1,1
#
#Test:
#input=[1 1 1 1 1
#       2 2 1 1 2
#       3 1 3 2 3]
#A = convert_input_to_A(input)
#Matrix(A)'

function convert_input_to_A(hapblock,MAF=0.0)
    if typeof(hapblock) == String
        df = readdlm(hapblock,Int64)
    else
        df = Int.(hapblock)
    end
    nind   = size(df,1)
    nblock = Int((size(df,2)-1)/2)

    num_haps = zeros(Int64,nblock)
    breaks   = zeros(Int64,nblock)
    # number of haplotypes for each block
    for i in 1:nblock
        num_haps[i]=length(union(df[:,2*i],df[:,2*i+1]))
    end
    #
    breaks[1]=0
    for i in 2:nblock
        breaks[i]=num_haps[i-1]+breaks[i-1]
    end

    rowindex=[]
    colindex=[]
    value=[]
    for k in 1:2
        push!(rowindex,Int64[])
        push!(colindex,Int64[])
        push!(value,Int64[])
    end

    for j in 1:nblock
        haps=union(df[:,2*j],df[:,2*j+1])
        for i in 1:nind
            val1=df[i,2*j]
            val2=df[i,2*j+1]
            hap1=findall(x -> x==val1,haps)
            hap2=findall(x -> x==val2,haps)
            push!(rowindex[1],breaks[j]+ hap1[1])
            push!(rowindex[2],breaks[j]+ hap2[1])
            push!(colindex[1],i)
            push!(colindex[2],i)
            push!(value[1],1)
            push!(value[2],1)
        end
    end

    rowindex=[rowindex[1];rowindex[2]]
    colindex=[colindex[1];colindex[2]]
    value=[value[1];value[2]]
    XMat=sparse(rowindex,colindex,value)

    #remove haplotypes with low frequency
    freq = vec(mean(XMat,dims=2)/2)
    freq = freq[freq .> MAF]
    A012 = XMat[freq .> MAF,:]
    #To make sparse matrix only has value "1",one unique haplotypes in one animal in one block
    #was calculated twice.
    A01 = replace!(A012, 2=>1)
    return A01,A012,freq
end
