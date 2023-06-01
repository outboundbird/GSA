
# require data from NCBI Gene Expression Omnibus (GEO) [database](https://www.ncbi.nlm.nih.gov/gds)
# https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
# datasets:
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5260 Psoriasis
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4500 AML train set
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4501 AML test set
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4844 CF
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4345 GI cancer muscle tissue
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5022 AMI
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5074 AMI 48 hrs
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5403 RA and OA
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5186 RIA
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3920 MS
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3628 RA
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3713 smoking
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4266 renal transplant
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4270 UC
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4181 AML
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3886 MS
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3898 social isolation
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4193 SLE
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1436 smoking
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1269 smoking
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1412 hormone therapy
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3005 IL1 and IL6
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3115 CHF
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3383 chronic stress
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3318 scikle cell disease
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3503 exercise
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4719 SLE
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS838 CML
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS534 smoking
library(GEOquery)

GSE18723 <- getGEO("GDS3713")
head(Meta(GSE18723))
