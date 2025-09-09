# A Generalized Evolutionary Classifier for Evolutionary Guided Precision Medicine



![publication_abstract](abs-fig.jpg)


**Purpose**

Current precision medicine (CPM) matches patients to therapies using traditional biomarkers, but inevitably resistance develops. Dynamic precision medicine (DPM) is a new evolutionary guided precision medicine (EGPM) approach undergoing translational development. It tracks intratumoral genetic heterogeneity and evolutionary dynamics, adapts as frequently as every 6 weeks, plans proactively for future resistance development, and incorporates multiple therapeutic agents. Simulations indicated DPM can significantly improve long-term survival and cure rates in a cohort of 3 million virtual patients representing a variety of clinical scenarios. Given the cost and invasiveness of monitoring subclones frequently, we sought to determine the value of a short DPM window of only two 6-week adaptations (moves).

**Methods**

In a new simulation, nearly 3 million virtual patients, differing in DPM input parameters of initial subclone compositions, drug sensitivities, and growth and mutational kinetics, were simulated as previously described. Each virtual patient was treated with CPM, DPM, and DPM for two moves followed by CPM.

**Results**

The first two DPM moves provide similar average benefit to a 5-year, 40-move sequence in the full virtual population. If the first two moves are identical for DPM and CPM, patients will not benefit from DPM (65% negative predictive value). A patient subset (20%) in which 2-move DPM and 40-move DPM provide closely similar outcomes has extraordinary predicted benefit (hazard ratio of DPM/CPM 0.03).

**Conclusion**

The first two DPM moves provide most of the clinical benefit of DPM, reducing the duration required for subclone monitoring. This also leads to an evolutionary classifier selecting patients who will benefit: those in whom DPM and CPM recommendations differ early. These advances bring DPM (and potentially other EGPM approaches) closer to potential clinical testing.

[link to publication](https://ascopubs.org/doi/10.1200/PO.23.00714)

## Running Simulation Code

### input parameters

### executing simulations
The code can be run in parallel using the bash scripts found in simulation/code/

## processing/agregating results

### output files
Zenodo archive


### code

### files for downstream analysis

## analysis

### population subset analysis

### indiviual patient analysis

Analysis scripts can be found in the R markdown file simulation/dpmAnalysis_trialSim.Rmd 

The analysis directory contiains code used for ongoing downstream ananlysis of the simulation results. 
