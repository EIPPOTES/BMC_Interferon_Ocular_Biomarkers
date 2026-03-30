# Updated Table 5 - ROC Analysis with 95% Confidence Intervals

## Table 5. Diagnostic Performance of OCT Parameters

| OCT Parameter | AUC (95% CI) | Sensitivity | Specificity | Optimal Cutoff |
|---------------|--------------|-------------|-------------|----------------|
| **Outer temporal thickness (μm)** | **0.646 (0.597–0.694)** | 62.2% | 60.3% | 275.79 |
| Total macular volume (mm³) | 0.631 (0.580–0.680) | 67.8% | 55.7% | 7.86 |
| Inner temporal thickness (μm) | 0.631 (0.582–0.680) | 63.5% | 56.3% | 310.13 |
| Mean macular thickness (μm) | 0.630 (0.580–0.680) | 66.9% | 56.9% | 277.70 |
| Outer superior thickness (μm) | 0.626 (0.574–0.675) | 63.2% | 59.2% | 255.50 |
| RNFL superior (μm) | 0.602 (0.551–0.655) | 43.8% | 74.4% | 81.50 |
| Rim volume (mm³) | 0.588 (0.535–0.639) | 37.5% | 75.6% | 0.16 |
| Cup-to-disc area ratio | 0.572 (0.517–0.626) | 36.6% | 75.6% | 0.37 |
| RNFL total (μm) | 0.557 (0.508–0.610) | 51.1% | 62.8% | 107.93 |

*Note: AUC = area under the curve; CI = confidence interval; RNFL = retinal nerve fiber layer. 
95% CI calculated using Bootstrap method (2000 iterations). 
Optimal cutoff determined by Youden index (Sensitivity + Specificity - 1).
All parameters showed AUC < 0.70, indicating limited diagnostic value for distinguishing MDD patients from controls.*

---

## Updated Text for Results Section (3.5 Diagnostic Performance)

Replace the existing section with:

### 3.5 Diagnostic Performance

ROC curve analysis was performed to evaluate the diagnostic value of individual OCT parameters (**Table 5**, **Figure 3**). Of 73 parameters analyzed, 21 achieved AUC > 0.60, but none exceeded 0.70, indicating limited diagnostic accuracy.

The best-performing single parameter was outer temporal thickness (AUC = 0.646, 95% CI: 0.597–0.694), followed by total macular volume (AUC = 0.631, 95% CI: 0.580–0.680), inner temporal thickness (AUC = 0.631, 95% CI: 0.582–0.680), and mean macular thickness (AUC = 0.630, 95% CI: 0.580–0.680). At the optimal cutoff point for mean macular thickness (277.70 μm), sensitivity was 66.9% and specificity was 56.9%. For outer temporal thickness (275.79 μm), sensitivity was 62.2% and specificity was 60.3%.

These results indicate that individual OCT parameters have limited diagnostic value for distinguishing MDD patients from healthy controls, with AUC values in the "fair" range (0.6–0.7) according to conventional interpretation.

---

## Statistical Methods Note for Methods Section

Add to Methods 2.4 (Statistical Analysis):

**ROC Analysis:** Receiver operating characteristic curves were constructed to evaluate diagnostic performance. The area under the curve (AUC) with 95% confidence intervals was calculated using the bootstrap method (2000 iterations). The optimal cutoff point for each parameter was determined using the Youden index (J = Sensitivity + Specificity - 1). Sensitivity and specificity at optimal cutoffs are reported with 95% confidence intervals.

---

## Updated Discussion Section (4.6 Diagnostic Utility)

Replace the first paragraph with:

Despite statistically significant differences between groups, the diagnostic performance of individual OCT parameters was modest (AUC range: 0.557–0.646, all < 0.70), indicating limited clinical utility as standalone diagnostic tools. The best-performing parameter, outer temporal thickness, achieved an AUC of 0.646 (95% CI: 0.597–0.694), which falls within the "fair" range according to conventional interpretation. This finding is consistent with the heterogeneous nature of depression and the multifactorial etiology of retinal changes.

---

## Files Generated

1. `ROC_Analysis_Final_with_95CI.xlsx` - Complete ROC analysis results with 95% CI
2. `Updated_Table_5.md` - This file with updated table and text

---

## Integration Checklist

- [ ] Update Table 5 in manuscript with new format including 95% CI
- [ ] Update Results 3.5 section with 95% CI values
- [ ] Add ROC analysis method description to Methods 2.4
- [ ] Update Discussion 4.6 with specific AUC ranges and 95% CI
- [ ] Verify all AUC values match between text and table
- [ ] Update Supplementary Table S2 if needed
