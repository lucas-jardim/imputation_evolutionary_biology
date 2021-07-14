# A Cautionary Note on Phylogenetic Signal Estimation from Imputed Databases

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains codes of the [paper](https://doi.org/10.1007/s11692-021-09534-0) published in Evolutionary Biology. 


## How to cite

> Jardim, L., Bini, L.M., Diniz-Filho, J.A.F. et al. 2021. A Cautionary Note on Phylogenetic Signal Estimation from Imputed Databases. Evolutionary Biology, 48, 246–258. [https://doi.org/10.1007/s11692-021-09534-0](https://doi.org/10.1007/s11692-021-09534-0).

## Data  

This paper evaluated the effects of imputed data on phylogenetic signal estimation using synthetic data and computational simulations. The simulations also applied on primate's brain and body size data extracted from [DeCasien et al. (2017)](https://doi.org/10.1038/s41559-017-0112). 

### Repository structure

```bash
imputation_evolutionary_biology/
├── graphics
│   ├── beta_RT.tif
│   ├── beta.trait1.tif
│   ├── imputation_error.tif
│   ├── imputation_graph.tif
│   ├── k_RT.tif
│   ├── ktrait1.tif
│   ├── mean_RT.tif
│   ├── mean.trait1.tif
│   ├── mechanism_scheme.png
│   ├── moran1.tif
│   ├── moran_RT.tif
│   ├── var_RT.tif
│   └── var.trait1.tif
├── imputation.Rproj
├── LICENSE
├── R
│   ├── graphics.R
│   ├── running_analyses.R
│   ├── run.simulation.R
│   ├── simulacao.2.1.R
│   └── Simulation_script.R
├── README.html
├── README.md
├── real_data
│   ├── data
│   │   ├── analysis_data.csv
│   │   ├── analysis_tree.tree
│   │   ├── DeCasien_data.xls
│   │   ├── primate_data.csv
│   │   └── Upham_phylogeny.nexus
│   ├── graphics
│   │   ├── beta_RT.tif
│   │   ├── beta.trait1.tif
│   │   ├── imputation_error.tif
│   │   ├── k_RT.tif
│   │   ├── ktrait1.tif
│   │   ├── mean_RT.tif
│   │   ├── mean.trait1.tif
│   │   ├── moran1.tif
│   │   ├── moran_RT.tif
│   │   ├── var_RT.tif
│   │   └── var.trait1.tif
│   ├── results
│   │   └── 0.6
│   │       ├── phylo
│   │       │   └── result.RData
│   │       ├── rand
│   │       │   └── result.RData
│   │       └── trait
│   │           └── result.RData
│   └── scripts
│       ├── Downloading_data.R
│       ├── Filtering_data.R
│       ├── graphics.R
│       ├── running_analyses.R
│       ├── run.simulation.R
│       ├── simulacao.2.1.R
│       └── Simulation_script.R
└── Results
    └── 0.6
        ├── phylo
        │   └── result.RData
        ├── rand
        │   └── result.RData
        └── trait
            └── result.RData

```

### Glossary

* **ktrait** = Blomberg's K of complete data      
* **ktrait1** = Blomberg's K of imputed data     
* **moran** = Moran's index of complete data          
* **moran1** = Moran's index of complete data       
* **mean.trait** = Trait's average of complete data
* **mean.trait1** = Trait's average of imputed data
* **var.trait** = Trait's variance of complete data 
* **var.trait1** = Trait's variance of imputed data 
* **beta.trait** = Regression coefficient of complete data 
* **beta.trait1** = Regression coefficient of imputed data
* **NRMSE** = Imputation error     
* **methods** = Imputation method  
* **mechanism** = Missing data mechanism 
* **alpha** = Simulated Ornstein-Uhlenbeck's alpha      
* **percentage** = Percentage of missing data

## Contributing

Contributions are welcome. Open an issue to discuss the changes. 

## Funding

* Coordenação de Aperfeiçoamento de Pessoal de Nível Superior

* Ministério da Ciência, Tecnologia e Inovações 

* Conselho Nacional de Desenvolvimento Científico e Tecnológico

* Fundação de Amparo à Pesquisa do Estado de Goiás

* CONACYT Ciencia Básica (A1-S-34563)