# Methods And Results Outline

## Methods

### GWAS input preparation

Describe:

- source GWAS traits
- summary statistic fields
- munging and alignment steps
- reference panel and ancestry assumptions

### LDSC analysis

Describe:

- univariate LDSC QC
- multivariable LDSC setup
- thresholds used for trait inclusion

### EFA and CFA

Describe:

- factor-number decision criteria
- EFA rotation and loading threshold
- CFA fitting strategy
- residual constraints used in the refined model

### Factor GWAS

Describe:

- native WSL GenomicSEM environment
- chunking strategy
- `userGWAS` settings
- final standardized exports

### LDSC validation of latent factors

Describe:

- univariate LDSC for `F1` and `F2`
- interpretation of heritability, intercept, and ratio

## Results

### LDSC QC and genetic correlation structure

Report:

- trait-level `h2` quality
- pairwise genetic correlation pattern

### EFA-supported latent structure

Report:

- retained factor number
- main loading pattern

### Refined CFA model

Report:

- model fit
- standardized loadings
- factor correlation

### Factor GWAS outputs

Report:

- total SNP count
- merged result completeness
- lead signals and top loci

### LDSC validation for `F1` and `F2`

Report:

- factor heritability
- intercept and ratio
- interpretation of reliability
