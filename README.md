# LPChoose.jl
Optimizing Sequencing Resources in Genotyped Livestock Populations by Linear Programming

> Hao Cheng, Keyu Xu, Kuruvilla Joseph Abraham, Optimizing Sequencing Resources in Genotyped Livestock Populations Using Linear Programming
bioRxiv 2020.06.29.179093; doi: https://doi.org/10.1101/2020.06.29.179093


### LPChoose function

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

### 1st Application
```julia
julia> LPChoose("smalldata.txt","unlimited")

--------------INPUT----------------------------
#Animal:6000
#Unique Haplotypes:27473
Haplotype Frequency Summary Stats:
Length:         27473
Missing Count:  0
Mean:           0.007139
Minimum:        0.000167
1st Quartile:   0.000167
Median:         0.000500
3rd Quartile:   0.002000
Maximum:        0.353000

--------------QUALITY CONTROL-------------------
----------minor haplotype frequency: 0.0--------
#Animal:6000
#Unique Haplotypes:27473
Haplotype Frequency Summary Stats:
Length:         27473
Missing Count:  0
Mean:           0.007139
Minimum:        0.000167
1st Quartile:   0.000167
Median:         0.000500
3rd Quartile:   0.002000
Maximum:        0.353000

---------------1ST APPLICATION--------------------
-------identify minimum number of animals---------
-------containing all unique haplotypes----------
-----------RUN LINEAR PROGRAMMING-----------------

It took 10.922356 seconds (942.22 k allocations: 121.610 MiB)

The minimum number of selected animals is: 4135
IDs for identified animals were saved in identified_animals.txt.

---------------------DONE-------------------------

```

### 2nd application
```julia
julia> LPChoose("smalldata.txt",10);

--------------INPUT----------------------------
#Animal:6000
#Unique Haplotypes:27473
Haplotype Frequency Summary Stats:
Length:         27473
Missing Count:  0
Mean:           0.007139
Minimum:        0.000167
1st Quartile:   0.000167
Median:         0.000500
3rd Quartile:   0.002000
Maximum:        0.353000

--------------QUALITY CONTROL-------------------
----------minor haplotype frequency: 0.0--------
#Animal:6000
#Unique Haplotypes:27473
Haplotype Frequency Summary Stats:
Length:         27473
Missing Count:  0
Mean:           0.007139
Minimum:        0.000167
1st Quartile:   0.000167
Median:         0.000500
3rd Quartile:   0.002000
Maximum:        0.353000

---------------2ND APPLICATION--------------------
------------identify best 10 animals--------------
--representing maximum proportions of haplotypes--
-----------RUN LINEAR PROGRAMMING-----------------

Step 1 took  0.415148 seconds (942.24 k allocations: 121.932 MiB)
Step 2 took  0.434525 seconds (930.95 k allocations: 110.687 MiB, 15.05% gc time)
Step 3 took  0.348062 seconds (921.41 k allocations: 102.636 MiB)
Step 4 took  0.283420 seconds (910.40 k allocations: 95.106 MiB)
Step 5 took  0.289956 seconds (901.46 k allocations: 89.485 MiB, 11.84% gc time)

0.05023113602446039 of the unique haplotypes in the population is covered.
0.5006245666272382 of the genome in the population is covered.

IDs for identified animals were saved in identified_animals.txt.

---------------------DONE-------------------------

```

# allow preselected animals

```julia
julia> animals_selected=[1,10,12,13,20,25]
julia> LPChoose("smalldata.txt",10,preselected_animals=animals_selected);

--------------INPUT----------------------------
#Animal:6000
#Unique Haplotypes:27473
Haplotype Frequency Summary Stats:
Length:         27473
Missing Count:  0
Mean:           0.007139
Minimum:        0.000167
1st Quartile:   0.000167
Median:         0.000500
3rd Quartile:   0.002000
Maximum:        0.353000

--------------QUALITY CONTROL-------------------
----------minor haplotype frequency: 0.0--------
#Animal:6000
#Unique Haplotypes:27473
Haplotype Frequency Summary Stats:
Length:         27473
Missing Count:  0
Mean:           0.007139
Minimum:        0.000167
1st Quartile:   0.000167
Median:         0.000500
3rd Quartile:   0.002000
Maximum:        0.353000

---------------2ND APPLICATION--------------------
------------identify best 10 animals--------------
--representing maximum proportions of haplotypes--
-----------RUN LINEAR PROGRAMMING-----------------

Step 1 took  0.374401 seconds (908.63 k allocations: 101.468 MiB, 5.89% gc time)
Step 2 took  0.283823 seconds (899.96 k allocations: 94.103 MiB)
Step 3 took  0.289010 seconds (891.79 k allocations: 87.604 MiB, 9.55% gc time)
Step 4 took  0.276065 seconds (885.11 k allocations: 83.597 MiB, 12.94% gc time)
Step 5 took  0.274014 seconds (876.84 k allocations: 79.381 MiB, 19.45% gc time)

0.07884104393404434 of the unique haplotypes in the population is covered.
0.619383149107422 of the genome in the population is covered.


0.0440448759113059 of the unique haplotypes in the population (exclude preselected animals) is covered.
0.47576847092310137 of the genome in the population (exclude preselected animals) is covered.

IDs for identified animals were saved in identified_animals.txt.

---------------------DONE-------------------------

```
