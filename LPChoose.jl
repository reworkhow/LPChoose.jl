#include("LPChoose.jl"); LPChoose("smalldata.txt",100,0.00);
using DataFrames,LinearAlgebra,CSV,DelimitedFiles,SparseArrays,GLPK,JuMP,Statistics, StatsBase
using ProgressMeter
"""
    LPChoose(hapblock,budget=100,MAF=0.0;
             nsteps= (budget=="unlimited" ? 1 : Int(ceil(budget/2)),
             preselected_animals    = false,
             weights_for_haplotypes = "haplotype frequency",
             sequencing_homozygous_haplotypes_only = false)

* Choose animals for sequencing given haplotype information **hapblock** filterd by minor haplotype frequency **MAF** for two applications:
    * identify minimum number of animals containing all unique haplotypes in the population if `budget = "unlimited"`;
    * identify a fixed number of animals whose haplotypes include as large a proportion as possible of the haplotypes
      present in the population given a limited **budget**, defaulting to `100` (100 animals).
* If a text file is provided for **hapblock**, the file format should be:
    * ```
      1,1,1,1,4       #ind1, hap1_1, hap1_1, hap2_1, hap2_4
      2,2,1,1,2       #ind2, hap1_2, hap1_1, hap2_1, hap2_2
      3,1,3,2,3       #ind3, hap1_1, hap1_3, hap2_2 hap2_3
      ```
    where individual IDs (they are required to be intergeres) are in 1st column, maternal and paternal haplotypes
    for haplotype block 1 are in column 2-3, maternal and paternal haplotypes for haplotype block 2 are in column 4-5.
* MISC
    * A fast approximation may be used to speed up computation in practice to select a fixed number of animals.
      This approximation is performed by selecting **budget** animals in **nsteps**, defaulting to selecting 2 animals
      at each step. For example, we can select 2 animals in each step to select 100 animals with 100/2=50 steps.
    * A list of preselected animals can be provided as an array of animal IDs for **preselected_animals**.
    * To identify a fixed number of animals, multiple options for `weights_for_haplotypes` are available, including
      "haplotype frequency" (default), "rare haplotype preferred", and "equal".
    * If `sequencing_homozygous_haplotypes_only`=`true`, LPChoose will only focus on sequencing homozygous haplotype
      segments to achieve a reduction in cost with an added benefit of phasing variant calls efficiently (Bickhart et al. 2015).


"""
function LPChoose(hapblock,budget="unlimited",MAF=0.0;
                  nsteps= (budget=="unlimited" ? 1 : Int(ceil(budget/2))),  #budget is #of selected animals
                  preselected_animals = false,
                  weights_for_haplotypes = "haplotype frequency", #"equal", "haplotype frequency", "rare haplotype preferred",
                  sequencing_homozygous_haplotypes_only = false)
    #Get incidence matrix
    A01_all,freq_all, animals_all = convert_input_to_A(hapblock,MAF,sequencing_homozygous_haplotypes_only)
    A01,freq, animals = select_these_animals(A01_all,freq_all,animals_all,preselected_animals)
    nind     = length(animals)
    #User-defined weights for haplotypes used in application 2
    if weights_for_haplotypes == "haplotype frequency"
        weights_for_haplotypes = freq
    elseif weights_for_haplotypes == "equal"
        weights_for_haplotypes = ones(freq)
    elseif weights_for_haplotypes == "rare haplotype preferred"
        weights_for_haplotypes = (freq .- 1).^2 #IWS weights
    elseif length(weights_for_haplotypes) != length(freq)
        error("the length of weights_for_haplotypes has to be equal to the number of haplotypes.")
    end

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
        if budget > length(animals)
            error("try to identify best $budget animals from a population of $(length(animals)) animals.")
        end
        animals_now      = copy(animals)
        A01now           = copy(A01)
        budget_each_step = Int(budget/nsteps)
        importance       = Matrix(weights_for_haplotypes'A01now)#get importance of each individual

        #Create text files to save output at each step
        file_genome_coverage    =open("genome_coverage.txt","w")
        file_hap_coverage       =open("haplotype_coverage.txt","w")
        file_identified_animals =open("identified_animals.txt","w")
        writedlm(file_genome_coverage,["step" "genome_covereage"],',')
        writedlm(file_hap_coverage,["step" "haplotype_coverage"],',')
        writedlm(file_identified_animals, ["step" "ID"],',')
        writedlm(file_genome_coverage,["0" 1-sum(A01now)/sum(A01_all)],      ',')
        writedlm(file_hap_coverage,   ["0" 1-size(A01now,1)/size(A01_all,1)],',')
        if preselected_animals != false
            writedlm(file_identified_animals,[fill("0",length(preselected_animals)) preselected_animals],',')
        else
            writedlm(file_identified_animals,["0" "NA"],',')
        end

        @showprogress "identifying most representative animals ..." for stepi = 1:nsteps
            nind  = length(animals_now)
            model = Model(GLPK.Optimizer)
            @variable(model, select_this_animal[1:nind],Bin)
            @objective(model, Max, sum([importance[i]*select_this_animal[i] for i= 1:nind]))
            @constraint(model, A01now * select_this_animal .<= 2)
            @constraint(model, sum(select_this_animal) <= budget_each_step) #e.g.,selecte 2 animals at each step
            #print("Step $stepi took") @time optimize!(model)
            optimize!(model)
            if JuMP.has_values(model) == false
                error("No Solutions")
            end
            #get identified animals at this step
            select_this_animal = (JuMP.value.(select_this_animal) .== 1.0) #boolean
            writedlm(file_identified_animals,[fill(string(stepi),sum(select_this_animal)) animals_now[select_this_animal]],',')
            #get identified haplotypes at this step
            haps_identified    = (vec(sum(A01now[:,select_this_animal],dims=2)) .!= 0) #boolean
            #update to current remaining animals for next round of selection
            A01now             = A01now[.!haps_identified,.!select_this_animal]

            weights_for_haplotypes  = weights_for_haplotypes[.!haps_identified]
            importance              = Matrix(weights_for_haplotypes'A01now)#get importance of each individual
            animals_now             = animals_now[.!select_this_animal]
            #save output to text files at this step
            writedlm(file_genome_coverage,[string(stepi) 1-sum(A01now)/sum(A01_all)],',')
            writedlm(file_hap_coverage,[string(stepi) 1-size(A01now,1)/size(A01_all,1)],',')
        end
        close(file_genome_coverage)
        close(file_hap_coverage)
        close(file_identified_animals)
        println("\n",1 - size(A01now,1)/size(A01_all,1), " of the unique haplotypes in the population is covered.")
        println(1 - sum(A01now)/sum(A01_all), " of the genome in the population is covered.\n")
        if preselected_animals != false
            println("\n",1 - size(A01now,1)/size(A01,1), " of the unique haplotypes in the population (exclude preselected animals) is covered.")
            println(1 - sum(A01now)/sum(A01), " of the genome in the population (exclude preselected animals) is covered.\n")
        end
        #writedlm("identified_animals.txt",setdiff(animals,animals_now))
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

function convert_input_to_A(hapblock,MAF=0.0,sequencing_homozygous_haplotypes_only=false)
    if typeof(hapblock) == String
        df = readdlm(hapblock,Int64)
    else
        df = Int.(hapblock)
    end
    animalIDs = df[:,1]
    nind      = size(df,1)
    nblock    = Int((size(df,2)-1)/2)

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

    if sequencing_homozygous_haplotypes_only == true
        A012 = replace!(A012, 1=>0) #MAF threshold will remove rows(haps) of all zeros
    end
    #To make sparse matrix only has value "1",one unique haplotypes in one animal
    #in one block was calculated twice.
    A01 = replace!(A012, 2=>1)
    freq = vec(mean(A01,dims=2))
    println("--------------INPUT----------------------------")
    println("#Animal:",size(A01,2))
    println("#Unique Haplotypes:",size(A01,1))
    println("Haplotype Frequency ",summarystats(freq))

    #remove haplotypes with low frequency
    A01  = A01[freq .> MAF,:]
    freq = freq[freq .> MAF]
    println("--------------QUALITY CONTROL-------------------")
    println("----------minor haplotype frequency: ",MAF,"--------")
    println("#Animal:",size(A01,2))
    println("#Unique Haplotypes:",size(A01,1))
    println("Haplotype Frequency ",summarystats(freq))
    return A01,freq,animalIDs
end


"""
    check_genome_coverage(hapblock,preselected_animals,MAF=0.0)

* Check haplotype and genome coverage for a given set of animails in a population
  defined by given haplotype information **hapblock** filterd by minor haplotype frequency **MAF**.
* If a text file is provided for **hapblock**, the file format should be:
    * ```
      1,1,1,1,4       #ind1, hap1_1, hap1_1, hap2_1, hap2_4
      2,2,1,1,2       #ind2, hap1_2, hap1_1, hap2_1, hap2_2
      3,1,3,2,3       #ind3, hap1_1, hap1_3, hap2_2 hap2_3
      ```
    where individual IDs (they are required to be intergeres) are in 1st column, maternal and paternal haplotypes
    for haplotype block 1 are in column 2-3, maternal and paternal haplotypes for haplotype block 2 are in column 4-5.
* If a text file is provided for **preselected_animals**, the file format should be
    * ```
    1
    3
    10
    ```
    where individual IDs (they are required to be intergeres) are in 1st column.
"""
function check_genome_coverage(hapblock,preselected_animals,MAF=0.0)

    A01,freq, animals = convert_input_to_A(hapblock,MAF)

    println("\n\n\n\n\n")
    println("--------------CHECK GENOME COVERAGE---------------")
    println("-------------by preselected animals---------------")

    A01_remaining, freq_remaining, animals_remaining = select_these_animals(A01,freq, animals,preselected_animals)
    println("\n",1 - size(A01_remaining,1)/size(A01,1), " of the unique haplotypes in the population is covered.")
    println(1 - sum(A01_remaining)/sum(A01), " of the genome in the population is covered.\n")

    println("------------Remaining Animals-----------------------")
    println("#Animal:",length(animals))
    println("Haplotype Frequency (remaining animals) ",summarystats(freq_remaining))
end

function select_these_animals(A01,freq,animals,preselected_animals)
    if preselected_animals == false
        return A01,freq,animals
    elseif typeof(preselected_animals) == String
        preselected_animals = readdlm(preselected_animals,Int64)
    else
        preselected_animals = Int.(preselected_animals)
    end
    if !issubset(preselected_animals,animals)
        error("Selected animals IDs are not found!")
    end
    select_this_animal = [i in preselected_animals for i in animals] #boolean
    haps_identified    = (vec(sum(A01[:,select_this_animal],dims=2)) .!= 0) #boolean
    A01_remaining      = A01[.!haps_identified,.!select_this_animal]
    animals_remaining  = animals[.!select_this_animal]
    freq_remaining     = freq[.!haps_identified]
    return A01_remaining, freq_remaining, animals_remaining
end
