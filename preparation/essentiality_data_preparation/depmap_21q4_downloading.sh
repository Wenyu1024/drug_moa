# This DepMap release contains data from CRISPR knockout screens from project Achilles, as well as genomic characterization data from the CCLE project.
cd /home/cloud-user/cluster_scratch/depmap_21q4/

# 1crispr (crispr_gene_effect.csv)
# File Description
# _Post-Chronos_ Combined Achilles and Sanger SCORE Chronos data using Harmonia (the batch correction pipeline described here: https://www.biorxiv.org/content/10.1101/2020.05.22.110247v3) - Columns: genes in the format "HUGO (Entrez)" - Rows: cell lines (Broad IDs)
# pwd
wget "https://ndownloader.figshare.com/files/31315996" -P ./
mv 31315996 crispr_gene_effect.csv


#3 exp_seq
# RNAseq TPM gene expression data for just protein coding genes using RSEM. Log2 transformed, using a pseudo-count of 1. - Rows: cell lines (Broad IDs) - Columns: genes (HGNC symbol and Entrez ID)
# wget "https://ndownloader.figshare.com/files/26261476" -P /home/cloud-user/cluster_scratch/depmap_21q1/
# mv 26261476 exp_seq

# 4cn
# number Gene level copy number data, log2 transformed with a pseudo count of 1. This is generated by mapping genes onto the segment level calls. - Rows: cell lines (Broad IDs) - Columns: genes (HGNC symbol and Entrez ID)
# wget "https://ndownloader.figshare.com/files/26261524" -P /home/cloud-user/cluster_scratch/depmap_21q1/
# mv 26261524 cn

# 5mut
# MAF of gene mutations. For all columns with AC, the allelic ratio is presented as [ALTERNATE:REFERENCE]. - CGA_WES_AC: the allelic ratio for this variant in all our WES/WGS(exon only) using a cell line adapted version of the 2019 CGA pipeline that includes germline filtering. - SangerWES_AC: in Sanger WES (called by sanger) (legacy) - SangerRecalibWES_AC: in Sanger WES after realignment at Broad (legacy) - RNAseq_AC: in Broad RNAseq data from the CCLE2 project (legacy) - HC_AC: in Broad Hybrid capture data from the CCLE2 project (legacy) - RD_AC: in Broad Raindance data from the CCLE2 project (legacy) - legacy_wgs_exon_only: in Broad WGS data from the CCLE2 project (legacy) Additional columns: - isTCGAhotspot: is this mutation commonly found in TCGA - TCGAhsCnt: number of times this mutation is observed in TCGA - isCOSMIChotspot: is this mutation commonly found in COSMIC - COSMIChsCnt: number of samples in COSMIC with this mutation - ExAC_AF: the allelic frequency in the Exome Aggregation Consortium (ExAC) Descriptions of the remaining columns in the MAF can be found here: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
# wget "https://ndownloader.figshare.com/files/26261527" -P /home/cloud-user/cluster_scratch/depmap_21q1/
# mv 26261527 mut

# 6cell info
# Cell line information definitions - DepMap_ID: Static primary key assigned by DepMap to each cell line - stripped_cell_line_name: Cell line name with alphanumeric characters only - CCLE_Name: Previous naming system that used the stripped cell line name followed by the lineage; no longer assigned to new cell lines - alias: Additional cell line identifiers (not a comprehensive list) - COSMIC_ID: Cell line ID used in Cosmic cancer database - lineage, lineage_subtype, lineage_sub_subtype, lineage_molecular_subtype: Cancer type classifications in a standardized form - sex: Sex of tissue donor if known - source: Source of cell line vial used by DepMap - Achilles_n_replicates: Number of replicates used in Achilles CRISPR screen passing QC - cell_line_NNMD: Difference in the means of positive and negative controls normalized by the standard deviation of the negative control distribution - culture_type: Growth pattern of cell line (Adherent, Suspension, Mixed adherent and suspension, 3D, or Adherent (requires laminin coating)) - culture_medium: Medium used to grow cell line - cas9_activity: Percentage of cells remaining GFP positive on days 12-14 of cas9 activity assay as measured by FACs - RRID: Cellosaurus research resource identifier - sample_collection_site: Tissue collection site - primary_or_metastasis: Indicates whether tissue sample is from primary or metastatic site - disease: General cancer lineage category - disease_subtype: Subtype of disease; specific disease name - age: If known, age of tissue donor at time of sample collection - Sanger_model_ID: Sanger Institute Cell Model Passport ID - additional_info: Further information about cell line modifications and drug resistance
# wget "https://ndownloader.figshare.com/files/26261569" -P /home/cloud-user/cluster_scratch/depmap_21q1/
# mv 26261569 sampleinfo



# 2demeter
# DEMETER2 Data v6 (combined)
# Cancer cell line genetic dependencies estimated using the DEMETER2 model. DEMETER2 is applied to three large-scale RNAi screening datasets: the Broad Institute Project Achilles, Novartis Project DRIVE, and the Marcotte et al. breast cell line dataset. The model is also applied to generate a combined dataset of gene dependencies covering a total of 712 unique cancer cell lines. For version history, see the description in the figshare. For more information visit https://depmap.org/R2-D2/.
# wget "https://ndownloader.figshare.com/files/13515395" -P /home/cloud-user/cluster_scratch/depmap_21q1/
# mv 13515395 demeter2