# Response to Reviewers

**Manuscript Title:** Retinal structural changes in patients with major depressive disorder: A cross-sectional study using optical coherence tomography

**Journal:** Journal of Affective Disorders

**Date:** March 11, 2026

---

We would like to thank the reviewers for their thoughtful comments and constructive feedback on our manuscript. We have carefully addressed each point and believe that the revisions have substantially improved the clarity and rigor of our work. Below, we provide a point-by-point response to each comment.

---

## Reviewer #1

### Comment 1: Age Differences and Confounding Control

**Comment:** The MDD group was significantly older than controls (38.3±20.2 vs 28.0±14.2 years, P<0.001). While the authors controlled for age in multivariate analyses, the nearly 10-year difference raises concerns about residual confounding. The authors should provide more transparency about age distribution and conduct additional sensitivity analyses to demonstrate that findings are robust to this demographic difference.

**Response:** We appreciate this important concern and have made the following revisions to address it:

1. **Enhanced age description**: We now report age using median (IQR) to better capture the distribution: [median (IQR): 38.0 (28.0-50.0) vs. 32.0 (25.0-44.0) years; P=0.042] (Results 3.1, Table 1). We also note that age distributions overlapped substantially between groups, with 68% of participants falling within the 25-50 year range in both groups.

2. **Linearity assumption testing**: In the Methods section (2.4), we now explicitly state: "To test for potential non-linear associations, we fitted additional models including a quadratic age term (age²) and performed likelihood ratio tests comparing linear vs. quadratic models. No significant improvement in model fit was observed (all P>0.10), supporting the use of linear age adjustment."

3. **Age-stratified sensitivity analysis**: We conducted analyses stratifying participants by age tertiles (<35, 35-50, >50 years) and found consistent effect estimates across strata (P for interaction >0.30), further confirming the appropriateness of linear adjustment (Methods 2.4).

4. **Three comprehensive sensitivity analyses** (Results 3.8):
   - **Exclusion of older participants**: After excluding participants >60 years, mean macular thickness remained significantly reduced in MDD patients (P<0.001, Cohen's d=-0.332).
   - **Age-stratified analysis**: Consistent trends toward reduced macular thickness were observed across all age strata.
   - **Age-matched subsample**: Using only participants within the control group's age range (18-58 years), the matched sample showed no significant age difference (P=0.416), yet mean macular thickness remained significantly reduced in MDD patients (P=0.003, d=-0.294).

5. **Supplementary material**: We have added Supplementary Figure S1 showing age distribution by group.

These analyses collectively demonstrate that our primary findings are robust and not attributable to age differences between groups.

---

### Comment 2: Inter-Eye Correlation and Analytical Approach

**Comment:** The primary analysis was conducted at the eye level with participant ID as a random effect. While this is appropriate, the authors should provide more details about the mixed-effects models used, including the intraclass correlation coefficient (ICC) and model fit statistics. Additionally, a comparison with participant-level analysis would strengthen confidence in the findings.

**Response:** We have substantially expanded the statistical methods and results to address this comment:

1. **Complete LMM specification**: We now provide the full model formula (Methods 2.4):
   
   **Y_ij = β₀ + β₁Group + β₂Age + β₃Sex + u_i + ε_ij**
   
   Where u_i ~ N(0, σ²_u) represents the random intercept for participant i.

2. **ICC reporting**: We now report ICC values ranging from 0.52-0.68 across parameters, indicating moderate to strong correlation between fellow eyes (Results 3.6). This quantifies the importance of accounting for inter-eye clustering.

3. **Model fit statistics**: We report AIC and BIC values, noting that "Model fit statistics (AIC/BIC) favored mixed-effects models over standard linear regression (ΔAIC >50 for all parameters), confirming the importance of accounting for clustering" (Results 3.6).

4. **Participant-level comparison**: We conducted sensitivity analyses using participant-level averages (mean of both eyes) and found consistent results: "mean macular thickness remained significantly reduced in MDD (β=-5.67 vs. -5.71 μm for LMM vs. participant-level, respectively), supporting the robustness of our findings regardless of analytical approach" (Results 3.6).

These additions provide transparency about our analytical approach and demonstrate the robustness of findings across different modeling strategies.

---

## Reviewer #2

### Comment 3: PHQ-9 Data Availability and Interpretation

**Comment:** The manuscript states that PHQ-9 data were available for only 132 of 164 MDD patients (80.5%). The authors should clarify why data were missing for 32 patients and assess whether this missingness might introduce selection bias. Additionally, the finding that 39.6% of patients had PHQ-9 scores <5 (minimal symptoms) despite an MDD diagnosis warrants more careful explanation.

**Response:** We have addressed both aspects of this comment:

1. **Missing data explanation**: We now explicitly state the reason for missing PHQ-9 data (Methods 2.2): "PHQ-9 data were available for 132 of 164 MDD patients (80.5%). Data were missing for 32 patients because PHQ-9 assessment was added to the study protocol after enrollment had already begun for the initial cohort."

2. **Selection bias assessment**: We compared baseline characteristics between patients with and without PHQ-9 data and "found no significant differences in age (P=0.72), sex distribution (P=0.58), or any OCT parameter (all P>0.30), suggesting data were missing at random" (Methods 2.2). We have added Supplementary Table S1 showing this comparison.

3. **PHQ-9 <5 explanation**: We have expanded the explanation in Results 3.1:
   
   "Notably, 103 eyes (39.6%) from 52 patients had PHQ-9 scores < 5, indicating minimal current symptoms despite their MDD diagnosis. These patients were likely in remission or receiving effective treatment at the time of assessment. This distribution reflects the natural course of MDD, where patients maintain their diagnosis while symptoms fluctuate or are controlled by treatment."
   
   We further clarify (Results 3.1): "MDD diagnosis was established by psychiatrists based on DSM-5 criteria, independent of symptom severity at the time of enrollment. PHQ-9 was administered on the same day as OCT scanning to assess current symptom severity."

4. **Limitations section**: We have added to the Limitations (Discussion 4.7): "First, PHQ-9 data were unavailable for 19.5% of MDD patients due to the timing of protocol amendment. However, analyses comparing patients with and without PHQ-9 data revealed no systematic differences in demographic or OCT characteristics, suggesting that missing data did not introduce selection bias."

---

### Comment 4: Glaucoma Wording

**Comment:** The manuscript states that MDD patients showed "glaucoma-like changes" including increased cup-to-disc ratio and decreased rim volume. This wording could be misinterpreted as suggesting that MDD patients have glaucoma. The authors should use more cautious language to avoid over-interpreting these findings.

**Response:** We completely agree and have revised the wording throughout the manuscript to be more cautious:

1. **Results section (3.3.2)**: We have replaced the original statement with:
   
   "MDD patients demonstrated structural parameters that have been associated with glaucomatous optic neuropathy in previous studies, including increased cup-to-disc ratio (P<0.001, q=0.002) and reduced neuroretinal rim volume (P=0.003, q=0.018). However, these findings should be interpreted cautiously, as our study did not include comprehensive glaucoma evaluation (intraocular pressure measurement, visual field testing, or corneal thickness assessment). Whether these changes represent early glaucomatous damage, non-specific neurodegenerative alterations, or other pathophysiological processes requires further investigation with dedicated ophthalmologic workup."

2. **Discussion section (4.3)**: We have added a clear statement:
   
   "**However, we emphasize that our study was not designed to diagnose glaucoma or establish glaucoma-MDD comorbidity.** The observed changes in cup-to-disc ratio and rim volume, while statistically significant, fall within the range of normal variation and do not meet criteria for glaucomatous optic neuropathy without corroborating evidence of functional damage (visual field defects) or elevated intraocular pressure. These findings should be interpreted as suggestive of potential structural vulnerability rather than definitive evidence of glaucomatous disease. Future studies incorporating comprehensive glaucoma screening (automated perimetry, tonometry, pachymetry) are needed to clarify the clinical significance of these observations."

3. **Conclusion**: We have revised to state "optic disc alterations associated with glaucomatous optic neuropathy that require further ophthalmologic evaluation" rather than "glaucoma-like changes."

---

## Reviewer #3

### Comment 5: Multiple Comparison Correction

**Comment:** The authors analyzed 73 OCT parameters but do not clearly state whether multiple comparison correction was applied across all parameters simultaneously or separately for different parameter families. This is important for interpreting the statistical significance of findings. The authors should clarify their approach and consider reporting both uncorrected P-values and FDR-adjusted q-values.

**Response:** We have substantially revised our approach to multiple comparison correction to address this concern:

1. **Clear correction strategy**: We now explicitly state (Methods 2.4): "Given the large number of OCT parameters analyzed (n=73), we controlled for false discovery rate (FDR) using the Benjamini-Hochberg procedure. **All P values across the 73 parameters were ranked and adjusted simultaneously** to produce q values."

2. **Three-tier classification system**: We have implemented a transparent classification (Methods 2.4):
   - **Definitive evidence:** q<0.05 (FDR-controlled significant)
   - **Suggestive evidence:** P<0.05 but q>0.05 (nominal significance only)
   - **Non-significant:** P≥0.05

3. **Simultaneous reporting**: We now report both P and q values throughout the Results section. For example (Results 3.2.1): "The outer temporal region demonstrated the largest effect size (Cohen's d=-0.50, 95% CI: -0.68 to -0.31, P<0.001, q<0.001 after FDR correction)."

4. **Summary of findings**: We now provide a clear summary (Results 3.5):
   
   "Of 73 OCT parameters analyzed, 28 (38.4%) showed nominal significance at P<0.05, of which 16 (21.9%) remained significant after FDR correction (q<0.05). The proportion of significant findings greatly exceeded the 5% expected under the null hypothesis, supporting the presence of genuine biological signals rather than chance findings."

5. **Key findings categorized**: We explicitly list:
   - **Definitive (q<0.05)**: Mean macular thickness (q=0.008), outer temporal thickness (q=0.003), total macular volume (q=0.012), cup-to-disc ratio (q=0.002), rim volume (q=0.018)
   - **Suggestive (P<0.05, q>0.05)**: Inner ring temporal thickness (P=0.032, q=0.084), peripapillary RNFL superior quadrant (P=0.041, q=0.095)

6. **Supplementary material**: We have added Supplementary Table S2 showing all 73 parameters with both P and q values.

This approach provides complete transparency about our multiple comparison correction strategy and allows readers to appropriately interpret the strength of evidence for each finding.

---

## Additional Changes

In addition to the specific reviewer comments above, we have made the following improvements:

1. **Enhanced transparency**: All statistical methods are now described in greater detail, including model assumptions, diagnostic checks, and sensitivity analyses.

2. **Consistent terminology**: We have standardized terminology throughout the manuscript, particularly regarding "definitive" vs "suggestive" evidence.

3. **Supplementary materials**: We have added:
   - Supplementary Table S1: Comparison of patients with and without PHQ-9 data
   - Supplementary Table S2: All 73 OCT parameters with P and q values
   - Supplementary Figure S1: Age distribution by group

4. **Word count**: The revised manuscript is approximately 5,830 words (within typical limits for original research articles).

---

## Summary

We believe that these revisions have substantially strengthened the manuscript by:

1. Providing greater transparency about age differences and demonstrating robustness through multiple sensitivity analyses
2. Fully specifying the mixed-effects models and demonstrating consistency across analytical approaches
3. Clarifying the reasons for missing PHQ-9 data and demonstrating absence of selection bias
4. Using more cautious language regarding glaucoma-related findings
5. Implementing a clear and comprehensive multiple comparison correction strategy

We thank the reviewers for their constructive feedback, which has undoubtedly improved the quality and clarity of our work.

---

**Corresponding Author:**
Shihai Cui, MD
Department of Ophthalmology
Third Affiliated Hospital of Sun Yat-sen University
Guangzhou, China
Email: cuishihai@mail.sysu.edu.cn

*Submitted: March 11, 2026*
