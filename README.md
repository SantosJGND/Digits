## KDE estimation of hand-written digits and population genetics.

This repository is composed of two main Jupyter notebooks:
- [Digits](https://nbviewer.jupyter.org/github/SantosJGND/Digits/blob/master/Digits.ipynb)

- [Digital_population_genetics](https://nbviewer.jupyter.org/github/SantosJGND/Digits/blob/master/Digital_population_genetics.ipynb)

*Note*: copy-paste the urls of jupyter notebooks hosted on GitHub onto [NBviewer](https://nbviewer.jupyter.org/) to view the analyses within.

## Context

The creation of this repository was stimulated by observations made in the course of the simulations of 
haplotype samples ([Stats Lab](https://github.com/SantosJGND/Stats_Lab)). 
In summary, as part of a larger project on the description of genetic variation in large genomics data sets, 
we simulated haplotype populations as a means to test the use of different descriptors of genetic variation. 

For this purpose haplotypes were generated from static populations simulated as allele frequency vectors 
drawn from the Beta distribution. Shape parameters of the Beta distribution were made to vary in order to 
approximate arbitrary selection and dominance scenarios. Because the descriptors used did not model demographic 
or historic parameters we did not complicate simulations further.

One of the benefits of the simulation protocol developped is the ability to include a function of genetic differentiation 
by proxy genomic position ([Fig. 2](https://github.com/SantosJGND/Stats_Lab/blob/master/Complementary_data/Ideo_comparison.png), 
[Stats Lab](https://github.com/SantosJGND/Stats_Lab), Notebook).
While the accompanying simulations rest on a number of simplifying assumptions, this development opened an interesting 
avenue of research. For example, we combined it with the inclusion individual admixture proportions to study our 
ability to identify local introgressions, foreign and reference, in admixed genomes as a function of genetic structure. 

In these simulations, individual admixture was modelled as a uniform probability of transition between origins associated
to a population ID. The hypothesis of a uniform distribution draws from our single estimate of admixture proportion. In reality,
patterns of admixture are likely to present non-random distributions at the individual and population level. 

In this context, we explored a method of managing local patterns of assignement at a population level that could be 
of use to model and simulate genomic data.

In order to explore this question without confusing our purpose with that of modelling population genetics parameters,
we will extract our prior of individual local genomic admixture proportions from the most obviously synthetic model
possible. 

**What if the historical dynamics of natural populations had produced patterns of introgression which, when organised 
in a certain fasion, would produce recognisable patterns?**

We resorted to the MINST data set of hand-written digits for our sythetic priors.

- Download the MNIST train and test data sets from the [MNIST Database](http://yann.lecun.com/exdb/mnist/)
