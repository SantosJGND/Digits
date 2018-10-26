## KDE estimation of hand-written digits and population genetics.

The creation of this repository was stimulated by observations made in the course of the simulations of 
haplotype samples ([Genetic Data Analysis repository](https://github.com/SantosJGND/Genetic-data-analysis)). 
In summary, as part of a larger project on the description of genetic variation in large genomics data sets, 
we simulated haplotype populations as a means to test the use of different descriptors of genetic variation. 

For this purpose haplotypes were generated from static populations simulated as allele frequency vectors 
drawn from the Beta distribution. Shape parameters of the Beta distribution were made to vary in order to 
approximate arbitrary selection and dominance scenarios. Because the descriptors used did not model demographic 
or historic parameters we did not complicate simulations further.

One of the benefits of the simulation protocol developped is the ability to include a function of genetic differentiation 
by proxy genomic position (Notebook 10, [Genetic Data Analysis repository](https://github.com/SantosJGND/Genetic-data-analysis)).
While the accompanying simulations rest on a number of simplifying assumptions, this development opened an interesting 
avenue of research. For example, we combined it with the inclusion individual admixture proportions to study our 
ability to identify local introgressions, foreign and reference, in admixed genomes as a function of genetic structure. 

This path led to the question: **can we include more complicated models of population admixture in our simulations?**
We wish to answer this question without delving into complicated modelling. That is, by what mechanism would we 
generate admixed genomes if we had access to historical and pedigree data. **In what format must that model come?**

In order to explore this question without confusing our purpose with that of modelling population genetics parameters,
we will extract our prior of individual local genomic admixture proportions from the most obviously synthetic model
possible. 

**What if the historical dynamics of natural populations had produced patterns of introgression which, when organised 
in a certain fasion, would produce recognisable patterns?**

We resorted to the MINST data set of hand-written digits for our sythetic priors.