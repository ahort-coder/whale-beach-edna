# whale-beach-edna

Code for **"Environmental DNA distinguishes nor'easter from hurricane overwash deposits"**
(Horton et al., in review). This code takes the raw 16S rRNA amplicon reads and produces
ASV tables, community statistics and the figures, using four salt-marsh sediment cores
from Whale Beach, New Jersey.


## Repository layout

| Folder / script | Purpose | Manuscript item |
|---|---|---|
| `00_install_packages.R` | Installs the R packages used | — |
| `01_sequence_processing/01_trim_polyCG_cutadapt.sh` | cutadapt trimming of reads (run in a terminal) | Methods: Sequence analysis |
| `01_sequence_processing/02_dada2_asv_inference.R` | DADA2 filtering, error learning, ASV inference, chimera removal, SILVA v138.1 taxonomy | Methods: Sequence analysis |
| `01_sequence_processing/03_build_phyloseq.R` | Combines per-sample ASV tables into one phyloseq object per core | Methods: Statistical analysis |
| `02_statistics/01_within_core_permanova_dbrda.R` | PERMANOVA, pairwise PERMANOVA, PERMDISP and dbRDA within each core | Supplementary Tables 4–8 |
| `02_statistics/02_across_core_permanova_permdisp.R` | Differences among cores (phylum, class and ASV level; rarefied sensitivity test) | Supplementary Tables 2–3 |
| `02_statistics/03_stratified_permanova_hurricane_vs_noreaster.R` | Hurricane vs. nor'easter test combined across cores | Results; Supplementary Table 6 |
| `02_statistics/04_alpha_diversity_WB2501.R` | Shannon and Simpson diversity | Results |
| `02_statistics/05_deseq2_WB2501.R`, `05b_…`, `05c_…` | DESeq2 differential abundance (WB2501) and taxonomy annotation | Results, Discussion |
| `02_statistics/06_indicator_species_WB2501.R` | IndVal indicator-species analysis (indicspecies multipatt) | Fig. 5 |
| `03_figures/fig1_field_site.py` | Study-site map | Fig. 1 |
| `03_figures/fig2_core_sedimentology.py` | Lithology, grain size and LOI profiles | Fig. 2 |
| `03_figures/fig3_pcoa.R` | PCoA ordinations | Fig. 3 |
| `03_figures/fig4_dbrda.R` | dbRDA ordinations | Fig. 4 |
| `03_figures/fig5_indicator_taxa.R` | Indicator-taxa dot plot | Fig. 5 |
| `03_figures/ed_fig2_class_composition.py`, `ed_fig3_depth_profiles.py`, `load.py` | Class-level composition figures (read the phyloseq `.rds` files in Python) | Extended Data Figs. 2–3 |
| `04_supplement/build_supplementary_tables.py` | Builds the Supplementary Information tables | Supplementary Tables |

## Running the code

1. Download the raw reads from SRA, or the processed data from Zenodo.
2. **File paths are hard-coded** to the author's computer (`~/Desktop/...`). Edit the path variables at the top of each script to point to your copies.
3. R scripts were run with R v4.6.1 (phyloseq, vegan, DESeq2, indicspecies, dada2, ggplot2, patchwork). Python scripts use Python 3 with numpy, pandas, matplotlib and rdata.

## Citation
If you use this code, please cite the article and the archived release: [Zenodo DOI].

## License
MIT (see `LICENSE`).
