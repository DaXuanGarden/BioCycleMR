# BioCycleMR
![BioCycleMR_logo](https://github.com/DaXuanGarden/BioCycleMR/assets/140375963/f421447f-ccf5-4b3f-a384-b5ea677083a3)


## Introduction

`BioCycleMR` is an R package crafted with the objective of enhancing Mendelian Randomization (MR) analysis in the field of biomedical research. Designed to integrate diverse exposure data types, the tool is an amalgamation of traditional medicinal insights and the dynamism of modern genomics.


## Key Features

### Comprehensive Exposure Data Integration

- **eQTL**: Seamlessly integrate expression quantitative trait loci data to unveil the genetic intricacies behind gene expression variations.
- **Immune Profiling**: Dive deep into the realm of immune cell data, highlighting potential associations and patterns crucial to disease research.
- **pQTL**: Explore protein quantitative trait loci data, shedding light on protein-level genetic nuances.
- **Inflammatory Markers & Gut Microbiota**: A dedicated module focusing on inflammatory factors and gut microbiota, revealing their complex interplay and implications on health.

### Efficient GWAS Data Integration

- **Automated Data Fetching**: Leverage automated mechanisms to gather pertinent GWAS data IDs efficiently, minimizing manual intervention and maximizing accuracy.
- **Iterative Loop Analyses**: A unique feature facilitating repeated cycles of MR analysis, ensuring comprehensive insights.

### Connectivity with Bioinformatics and Single-Cell Data

- **Bioinformatics Bridge**: A built-in framework designed to synchronize MR analyses with relevant bioinformatics datasets and tools.
- **Single-Cell Target Exploration**: Stay ahead by integrating critical targets deduced from single-cell analyses, thus providing a multi-dimensional research perspective.

## Installation and Usage

1. **Installing the Package**:
   ```r
   install.packages("devtools")
   devtools::install_github("DaXuanGarden/BioCycleMR")
   ```

2. **Load and Commence**:
   ```r
   library(BioCycleMR)
   results <- runBioCycleMR(your_exposure_data, your_outcome_data)
   print(results)
   ```

3. **Immune cells**
```r
data(immune)
outcome_ids = c("finn-b-N14_ENDOMETRIOSIS","finn-b-N14_FEMALEINFERT")
outcome_ids=as.vector(outcome_ids)
for (outcome_id in outcome_ids){
  write.table(analyze_outcome_immune(outcome_id=outcome_id,max_retries = 50),file = paste0(outcome_id,".txt"),sep = "\t",quote = F,row.names = F)
}
```

For an in-depth understanding and advanced functionalities, always refer to the package's extensive documentation, which elucidates various modules, functions, and their respective use cases.

## Future Developments

`BioCycleMR` is designed with adaptability and scalability in mind. While it stands as a reflection of the current knowledge and skills of its creators, they envisage it to evolve, incorporating advancements in biomedical research and feedback from the scientific community.

## Community Engagement and Feedback

Open communication channels and collaborations are the lifeblood of `BioCycleMR`. The creators earnestly invite the community to pitch in, share insights, suggest enhancements, or even critique – every interaction is a step towards refinement.

## Creators

- **Xuanyu Wang**: An undergraduate student majoring in Traditional Chinese Medicine from the 2021 cohort at Tianjin University of Traditional Chinese Medicine. With a keen interest in the convergence of age-old medical wisdom and contemporary genetic research, Xuanyu Wang plays an integral role in shaping the essence of `BioCycleMR`.

- **Yangyang Zhang**: An undergraduate student of the 2021 batch majoring in Clinical Medicine at Fudan University. Zhang's expertise and enthusiasm for eQTL and immune cell data have added depth to the capabilities of the package.

## Contact

For discussions, feedback, or potential collaborations:

📧 [Xuanyu Wang](mailto:daxuan111000@163.com)
📧 [Yangyang Zhang (pigudogzyy)](pigudogzyy@gmail.com)

