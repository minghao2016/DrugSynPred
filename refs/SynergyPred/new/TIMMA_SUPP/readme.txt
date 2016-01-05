
Dataset_S1.xlsx
<div><p><b>Binding affinity data for the CanOS1224 cell line case study.</b></p>
</div>
10.1371/journal.pcbi.1003226.s001

Dataset_S2.xlsx
<div><p><b>TIMMA prediction result for the CanOS1224 data as compared to PKIM prediction.</b></p>
</div>
10.1371/journal.pcbi.1003226.s002

Dataset_S3.xlsx
<div><p><b>Anticancer efficacy measured in activity area for the 15 kinase inhibitors selected from the CCLE cancer cell line drug collection.</b></p>
</div>
10.1371/journal.pcbi.1003226.s003

Dataset_S4.xlsx
<div><p><b>Binding affinity data for the MCF-7, BxPC-3 and MDA-MB-231 studies.</b></p>
</div>
10.1371/journal.pcbi.1003226.s004

Dataset_S5.xlsx
<div><p><b>TIMMA prediction result for MCF-7 cancer cell.</b></p>
</div>
10.1371/journal.pcbi.1003226.s005

Dataset_S6.xlsx
<div><p><b>Synthetic lethality score for the targets paired with AURKB for BxPC-3 cancer cell.</b></p>
</div>
10.1371/journal.pcbi.1003226.s006

Dataset_S7.xlsx
<div><p><b>TIMMA prediction result for MDA-MB-231 cancer cell.</b></p>
</div>
10.1371/journal.pcbi.1003226.s007

Dataset_S8.xlsx
<div><p><b>Single and double siRNA screen data for MDA-MB-231 cancer cell.</b></p>
</div>
10.1371/journal.pcbi.1003226.s008

Figure_S1.tif
<div><p><b>Flowchart for TIMMA model construction.</b></p>
</div>
10.1371/journal.pcbi.1003226.s009

Figure_S2.tif
<div><p><b>Flowchart of the SFFS algorithm for TIMMA model selection.</b></p>
</div>
10.1371/journal.pcbi.1003226.s010

Figure_S3.tif
<div><p><b>Concordance between TIMMA model assumptions and the actual data.</b> The data contains those 12 kinase inhibitors for which the quantitative binding affinity <img src="http://www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1003226.e168&amp;representation=PNG"> profiles across 384 kinase targets can be obtained from <a href="http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1003226#pcbi.1003226-Davis1" target="_blank">[41]</a>. The drug treatment efficacy data was obtained for the CCLE collection of 504 cancer cell lines measured by Activity area <a href="http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1003226#pcbi.1003226-Pal1" target="_blank">[30]</a>. The concordance index was calculated between the actual relationship between two drug efficacies (i.e. greater than or smaller than) and the model prediction using the basic subset and superset rules. The 95% confidence interval at each threshold was derived by summarizing the concordance indices for 504 cell lines. (<b>A</b>) The binary drug-target profiles were determined using a drug-specific <img src="http://www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1003226.e169&amp;representation=PNG"> threshold defined as the n-fold of the minimal <img src="http://www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1003226.e170&amp;representation=PNG"> value for each drug. The feasibility of the TIMMA model assumptions is manifested by significant enhancement of the concordance index compared to random predictions (paired t-test; p-value<10<sup>&#8722;15</sup>). (<b>B</b>) The binary drug target profiles were determined using a fixed level of <img src="http://www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1003226.e171&amp;representation=PNG"> cut-off threshold. When the threshold is lower than 600 nM the model assumptions cannot be tested as none of the binarized drug-target inhibition profiles is totally inclusive of each other. At higher cut-off thresholds the model prediction performs no better than random prediction.</p>
</div>
10.1371/journal.pcbi.1003226.s011

Figure_S4.tif
<div><p><b>Synergy scores do not necessarily correlate with the single drug treatment efficacies.</b> (<b>A</b>) Scatter plot between synergy scores and average single drug sensitivity scores for the selected 68 drug pairs in the MDA-MB-231 study, fitted by a loess smoothing function. (<b>B</b>) Synergy scores for the selected drugs when paired with dasatinib (blue line, right axis) and the corresponding individual drug sensitivity scores (grey bars, left axis).</p>
</div>
10.1371/journal.pcbi.1003226.s012

Table_S1.docx
<div><p><b>The top most synergistic drug pairs based on the TIMMA model predictions in MCF-7 breast cancer cells.</b></p>
</div>
10.1371/journal.pcbi.1003226.s013

Table_S2.docx
<div><p><b>Fold-changes of the average inhibition percentages between the kinase groups in the MDA-MB-231 siRNA screen.</b></p>
</div>
10.1371/journal.pcbi.1003226.s014
