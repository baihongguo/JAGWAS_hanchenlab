# Power simulation

Power of JAGWAS and Bonferroni minP for 128 phenotypes, computed from simulated
z-scores (Methods, Power simulation; Supplementary Table 6 and Supplementary Fig. 10).

## Files

| File | What it does |
|---|---|
| `pz01_zscore_power.R` | The simulation. |
| `pz02_plot_power.R` | Redraws the power figure from the results of `pz01`, with the labels used in the paper. |
| `common.R` | Functions shared with the null simulation (the same file as in `null_simulation/`). `pz01` needs it in this folder. |

## Design

- JAGWAS and minP use only the k single-phenotype z-scores of a variant and their
  correlation matrix R, so the z-scores are drawn directly as z ~ N(δ, R), where δ
  is the vector of expected z-scores. No genotypes or GWAS are simulated.
- The total signal is fixed. δ is scaled so that the noncentrality λ = δ′R⁺δ
  (R⁺ is the Moore–Penrose inverse, equal to R⁻¹ for an invertible R) is 127.7 or
  105.1, which gives JAGWAS with 128 degrees of freedom 80% or 50% power at
  α = 5 × 10⁻⁸.
- Three tests: JAGWAS with the standard inverse on 128 degrees of freedom; JAGWAS
  with the Moore–Penrose generalized inverse (eigenvalues below 10⁻⁵ set to zero,
  degrees of freedom equal to the rank), which is identical to the standard
  inverse unless an eigenvalue is below 10⁻⁵; and Bonferroni minP (max |z_j| > 6.26).
  JAGWAS power is computed from the noncentral chi-square distribution. The power
  of minP is estimated from 1,000 random sets of associated phenotypes with 1,000
  draws each.
- Number of associated phenotypes: m = 1, 4, 8, 16, 32 and 64.
- Effect patterns (`--patterns`):
  - equal effects on the associated phenotypes and none on the others (`direct` in
    the code and the results file);
  - weighted effects, δ = Rw with w equal on the associated phenotypes and zero
    elsewhere (`correlated`).
- Correlation structures (`--scenarios`):
  - `empirical`: the UDIP correlation matrix given with `--corr`;
  - `nearsingular`: 64 random UDIPs and a near-duplicate of each (noise variance
    10⁻⁶), which has 64 eigenvalues of about 5 × 10⁻⁷ and rank 64. A UDIP and its
    near-duplicate have the same effect, so m is even and m = 1 is skipped;
  - `signflip`: mixed effect directions, with each associated phenotype given a
    negative effect with probability 1/2 (the code reverses the sign of each
    phenotype with probability 1/2, R → DRD, and gives the effects one sign in
    the new coding).
- Without signal (δ = 0), 10⁶ draws per matrix check the type I error of each test.
- Seed: 20261028.

## Requirements

R (we used version 4.0.2) with the package data.table.

## Run

```bash
Rscript pz01_zscore_power.R --corr=T1_correlation_matrix.txt --out=zscore_T1
Rscript pz01_zscore_power.R --corr=T2_correlation_matrix.txt --out=zscore_T2
Rscript pz02_plot_power.R --results=zscore_T1/zpower_results.tsv --out=zscore_T1 --targets=0.8
```

`--corr` is a 128 × 128 correlation matrix; the paper used the T1 and T2 UDIP
matrices pooled over the discovery and replication cohorts (not included). The
defaults of `--corr` and `--out` are the paths on our servers, so pass both. A run
takes about 50 minutes with R's reference BLAS and about 10 minutes with an
optimized BLAS; `--n_sets=50 --n_draws=500 --n_null=100000` gives a test run of
about 2 minutes. The header of each script lists all options.

## Output (in `--out`)

| File | Contents |
|---|---|
| `zpower_matrix_summary.tsv` | Eigenvalues, condition number, rank, tr(R⁻¹)/k and the Li–Ji effective number of tests of each matrix. |
| `zpower_results.tsv` | Power by correlation structure, effect pattern, m and target power (Supplementary Table 6a). |
| `zpower_null_check.tsv` | Type I error of each test without signal, with exact binomial 95% confidence intervals (Supplementary Table 6b). |
| `zpower_power_80.png`, `zpower_power_50.png` | Power against m (PDF with `--plot_format=pdf` or when R has no PNG device). Supplementary Fig. 10 is the T1 matrix at 80%. |
| `pz01_zscore_power.log` | The log. |
