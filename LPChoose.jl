#include("LPChoose.jl"); LPChoose("smalldata.txt",100,0.00);
using DataFrames,LinearAlgebra,CSV,DelimitedFiles,SparseArrays,GLPK,JuMP,Statistics, StatsBase
"""
    LPChoose(hapblock,budget=100,MAF=0.0;nsteps= (budget=="unlimited" ? 1 : Int(ceil(budget/2)))

* Choose animals for sequencing given haplotype information **hapblock** filterd by minor haplotype frequency **MAF** for two applications:
    * identify minimum number of animals containing all unique haplotypes in the population if `budget = "unlimited"`;
    * identify a fixed number of animals whose haplotypes include as large a proportion as possible of the haplotypes
      present in the population given a limited **budget**, defaulting to `100` (100 animals).
* A fast approximation may be used to speed up computation in practice to select a fixed number of animals. This approximation
  is performed by selecting **budget** animals in **nsteps**, defaulting to selecting 2 animals at each step. For example,
  we can select 2 animals in each step to select 100 animals with 100/2=50 steps.
* If a text file is provided for **hapblock**, the file format should be:
    * ```
      1,1,1,1,4       #ind1, hap1_1, hap1_1, hap2_1, hap2_4
      2,2,1,1,2       #ind2, hap1_2, hap1_1, hap2_1, hap2_2
      3,1,3,2,3       #ind3, hap1_1, hap1_3, hap2_2 hap2_3
      ```
    where individual IDs are in 1st column, maternal and paternal haplotypes for haplotype block 1 are in column 2-3,
    maternal and paternal haplotypes for haplotype block 2 are in column 4-5.

"""
function LPChoose(hapblock,budget="unlimited",MAF=0.0;
                  nsteps= (budget=="unlimited" ? 1 : Int(ceil(budget/2)))) #budget is #of selected animals
    #Get incidence matrix
    A01,freq = convert_input_to_A(hapblock,MAF)
    animals  = collect(1:size(A01,2))
    nind     = length(animals)
    if budget == "unlimited" #application 1
        println("---------------1ST APPLICATION--------------------")
        println("-------identify minimum number of animals---------")
        println("-------containing all unique haplotypes----------")
        println("-----------RUN LINEAR PROGRAMMING-----------------\n")
        importance = ones(nind)
        model = Model(GLPK.Optimizer)
        @variable(model, select_this_animal[1:nind],Bin)
        @objective(model, Min, sum([importance[i]*select_this_animal[i] for i= 1:nind]))
        @constraint(model, A01 * select_this_animal .>= 1)
        print("It took")
        @time optimize!(model)
        if JuMP.has_values(model)
            print("\nThe minimum number of selected animals is: ")
            printstyled(Int(objective_value(model)),"\n\n",bold=true,color=:red)
            select_this_animal = (JuMP.value.(select_this_animal) .== 1.0)
            writedlm("identified_animals.txt",animals[select_this_animal])
            println("IDs for identified animals were saved in identified_animals.txt.\n")
        else
            error("No Solutions")
        end
    else #application 2
        println("---------------2ND APPLICATION--------------------")
        println("------------identify best $budget animals--------------")
        println("--representing maximum proportions of haplotypes--")
        println("-----------RUN LINEAR PROGRAMMING-----------------\n")
        animals_now      = copy(animals)
        A01now           = copy(A01)
        budget_each_step = Int(budget/nsteps)
        importance       = Matrix(freq'A01now)#get importance of each individual
        for stepi = 1:nsteps
            nind  = length(animals_now)
            model = Model(GLPK.Optimizer)
            @variable(model, select_this_animal[1:nind],Bin)
            @objective(model, Max, sum([importance[i]*select_this_animal[i] for i= 1:nind]))
            @constraint(model, A01now * select_this_animal .<= 2)
            @constraint(model, sum(select_this_animal) <= budget_each_step) #e.g.,selecte 2 animals at each step
            print("Step $stepi took")
            @time optimize!(model)
            if JuMP.has_values(model) == false
                error("No Solutions")
            end
            select_this_animal = (JuMP.value.(select_this_animal) .== 1.0) #boolean
            haps_identified    = (vec(sum(A01now[:,select_this_animal],dims=2)) .!= 0) #boolean
            A01now             = A01now[.!haps_identified,.!select_this_animal]
            freq               = freq[.!haps_identified]
            importance         = Matrix(freq'A01now)
            animals_now        = animals_now[.!select_this_animal]
            #println("number of animals selected at step", stepi, " is: ",sum(select_this_animal))
        end
        println("\n",1 - size(A01now,1)/size(A01,1), " of the unique haplotypes in the population is covered.")
        println(1 - sum(A01now)/sum(A01), " of the genome in the population is covered.\n")
        writedlm("identified_animals.txt",setdiff(animals,animals_now))
        println("IDs for identified animals were saved in identified_animals.txt.\n")
    end
    println("---------------------DONE-------------------------")
end

#the input file is haplotypes blocks
#hapblock is a file where each row representing one individual and
#each column representing one haplotype block, e.g.
#
#column 1: Individual Name
#column 2-3: Maternal and Paternal Haplotypes block 1
#column 4-5: Maternal and Paternal Haplotypes block 2
#
#ind1 hap1_1 hap1_1 hap2_1 hap2_4
#ind2 hap1_2 hap1_1 hap2_1 hap2_2
#ind3 hap1_1 hap1_3 hap2_2 hap2_3
#i.e.,
#1,1,1,1*,4*
#2,2,1,1*,2*
#3,1,3,2*,3*
#
#
#
#convert input file to "incidence matrix" of zeros/ones where each row represents
# one haplotype in the haplotype block,e.g, A'=
#ID,  hap1_1,hap1_2,hap1_3,hap2_1,hap2_2,hap2_3,hap2_4
#ind1,1,0,0,1,0,0,1
#ind2,1,1,0,1,1,0,0
#ind3,1,0,1,0,1,1,0
#
#Test:
#input=[1 1 1 1 4
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
    A012=sparse(rowindex,colindex,value)
    #To make sparse matrix only has value "1",one unique haplotypes in one animal in one block
    #was calculated twice.
    A01 = replace!(A012, 2=>1)
    freq = vec(mean(A01,dims=2))
    println("--------------INPUT----------------------------")
    println("#Animal:",size(A01,2))
    println("#Unique Haplotypes:",size(A01,1))
    println("Haplotype Frequency ",summarystats(freq))

    #remove haplotypes with low frequency
    freq = freq[freq .> MAF]
    A01  = A01[freq .> MAF,:]
    println("--------------QUALITY CONTROL-------------------")
    println("----------minor haplotype frequency: ",MAF,"--------")
    println("#Animal:",size(A01,2))
    println("#Unique Haplotypes:",size(A01,1))
    println("Haplotype Frequency ",summarystats(freq))
    return A01,freq
end
