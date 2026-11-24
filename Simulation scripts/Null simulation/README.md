# Null simulation

Scripts for the null simulation that evaluates the type I error of JAGWAS and minP
(Methods, Null simulation; Supplementary Tables 3–5 and Supplementary Figs. 8 and 9).

Phenotypes are simulated independently of genotype for UK Biobank discovery-cohort
participants and tested against the imputed genotypes of the 22,250 unrelated
participants analyzed, so every variant is null. There are 128 phenotypes, 257
replicates and four scenarios:

| Scenario (folder) | Phenotypes |
|---|---|
| `S1_normal` | Multivariate normal with correlation matrix R |
| `S2_lognormal` | Log-normal (σ = 0.5 on the log scale; skewness 1.75), Pearson correlation matched to R |
| `S3_cauchy` | Standard Cauchy through a Gaussian copula with latent correlation R (same rank correlations as S1) |
| `S4_normal_miss1pct` | S1 with each value missing independently with probability 0.01 |

Within a replicate, the four scenarios are built from the same standard normal draws,
so differences between scenarios reflect the marginal distributions and the
missingness rather than Monte Carlo error.

## Not included

- The UK Biobank genotypes and phenotypes (available to approved UK Biobank projects)
  and the UDIP correlation matrix used as R.
- The single-phenotype GWAS of the simulated phenotypes, which were run with
  TorchGWAS 1.0 (Zhao et al., arXiv 2026): ordinary linear regression on allele
  dosage, with the t statistics used as z-scores. The scripts read its output: for
  each replicate, a float32 file of t statistics in variant-major order
  (8,931,083 variants × 128 phenotypes), with the variant list and an allele
  frequency table.

## Requirements

R (we used version 4.0.2) with the package data.table; Matrix and CompQuadForm are
used when installed. The check in `06` also needs the JAGWAS executable from this
repository.

Every input and output path is a command-line option (`--name=value`). The defaults
are the paths on our servers and must be replaced. The header of each script lists
its options.

## Scripts, in the order they were run

| Step | Script | What it does |
|---|---|---|
| Simulate | `00_check_correlation_matrix.R` | Checks R before simulating: eigenvalues, condition number, Li–Ji effective number of tests, and whether the log-normal correlation can be matched. |
| | `01_simulate_null_phenotypes.R` | Simulates the phenotypes; replicate r uses seed 20260902 + r. All four scenarios of a replicate are written together; `--rep_start` and `--rep_end` split the replicates across jobs. |
| | `02_validate_simulation.R` | Checks the simulated correlations, skewness, kurtosis and missingness. |
| | `03_check_scenario_misspecification.R` | Checks whether the realized phenotype correlation alone would miscalibrate the test. |
| GWAS | TorchGWAS (not included) | Single-phenotype association tests for each replicate. |
| JAGWAS | `04_estimate_scenario_R.R` | Estimates R for each replicate as the sample correlation matrix of its simulated phenotypes (pairwise-complete in S4); `--keep` limits it to listed individuals. `--overlap_adjust=fill` writes the overlap-adjusted matrix G/N used for S4 in the paper (`drop` writes the version for GWAS software that drops missing values). |
| | `08_check_layout.R` | Confirms the layout of the TorchGWAS output file. |
| | `05_prepare_jagwas_input.R` | Writes the JAGWAS input files and computes T = z′R⁻¹z directly, for one scenario and a range of replicates (`--rep_start`, `--rep_end`). `--write_sumstats=false` skips the input files when only the direct statistic is needed. |
| | `06_compare_jagwas_direct.R` | Checks that the JAGWAS executable gives the same p-values as the direct computation (done for replicate 1 of every scenario). Run `05` with `--save_direct_p=true`, run JAGWAS on its input files with `--beta_se 1 --logP 1` (the files give the t statistic as BETA with SE = 1) and save the output as `<outbase>/<scenario>/rep###.jagwas.txt.gz`. |
| | `10_maf_class_pass.R` | Splits the statistics into common (MAF > 0.01) and rare (MAF ≤ 0.01) variants. |
| | `12_trait_subsets.R` | The same test for nested random subsets of 64 and 16 phenotypes; replicate r uses seed 20261008 + r. |
| | `13_minp_pass.R` | Bonferroni minP, min(1, k·p_min), and effective-tests minP, min(1, M_eff·p_min), on the same z-scores and subsets. |
| Results | `07_summarise_type1.R` | Empirical type I error rates with 95% confidence intervals from a cluster bootstrap with replicates as clusters (10,000 resamples), and λ_GC (Supplementary Tables 3 and 5). The output also lists exact Poisson intervals, which ignore LD, for reference. |
| | `11_plot_qq.R` | QQ plots (Supplementary Figs. 8 and 9). |
| | `09_maf_and_univariate.R`, `14_summarise_maf.R` | Rates by MAF bin for JAGWAS and for the single-phenotype statistics (Supplementary Table 4). `09` is run in chunks of replicates (`--out_tag`), and `14` pools them. |
| | `15_check_bootstrap.R` | Compares the bootstrap intervals with the spread of the per-replicate counts. |
| | `common.R` | Functions shared by the scripts; not run on its own. |

For S4 with the overlap-adjusted matrix, run `04` with `--overlap_adjust=fill`, then
`05`, `10`, `12` and `13` with `--rhat_prefix=Rhat_adjfill --out_scenario=S4_adjfill`.

## Example

```bash
Rscript 00_check_correlation_matrix.R --corr=R.txt --outdir=sim
Rscript 01_simulate_null_phenotypes.R --corr=R.txt --samples=ids.txt --outdir=sim
Rscript 02_validate_simulation.R --corr=R.txt --outdir=sim --reps=1,2,3
# single-phenotype GWAS of every replicate and scenario with TorchGWAS
Rscript 04_estimate_scenario_R.R --pheno_dir=sim/pheno --outdir=input
Rscript 05_prepare_jagwas_input.R --tstat_root=gwas --outdir=input --scenario=S1_normal \
    --rep_start=1 --rep_end=257 --variants=variants.bim
Rscript 07_summarise_type1.R --inbase=input --variants=all
Rscript 11_plot_qq.R --inbase=input --variants=all
```
