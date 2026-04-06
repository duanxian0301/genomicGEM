# cTWAS 与 SMR 分析模块说明

本模块用于在 ALPS latent factors (`F1`, `F2`) 与神经退行性疾病 (`AD`, `PD`, `LBD`) 之间开展基因层面的 follow-up 分析，重点包括 trait-wise multigroup cTWAS、bulk/cell-type SMR，以及与 `conjFDR`、`coloc`、`PWCoCo` 的综合候选基因优先级整合。

## 方法学命名约定

当前实现的 cTWAS 部分应表述为：

- **trait-wise multigroup cTWAS**
- **cTWAS-based cross-trait gene prioritization/comparison**

而不应表述为：

- full pairwise M-cTWAS joint analysis

原因是当前实现方式是在单一 trait 内联合多个脑组织 context 建模，而不是把 `F1` 与 `AD` 这样的两个 GWAS 在同一 joint model 中一起建模。

## SMR 的阈值与解释

为了兼顾论文主结果和补充材料，SMR 结果建议分三层报告：

1. `Bonferroni-significant + HEIDI-supported`
2. `FDR < 0.05 + HEIDI-supported`
3. `p_SMR < 1e-4 + HEIDI-supported` 作为探索性层

其中：

- `AD/PD/LBD` 可采用标准频率检查的 SMR
- `F1/F2` 因标准化导出中 `FREQ = MAF`，故在 SMR 中采用 `--disable-freq-ck`
- 因此 factor-side SMR 结果应作为较审慎的支持性证据，而不宜与 disease-side SMR 等量齐观

## 建议在论文中保留的核心输出

- cTWAS primary genes
- cTWAS secondary genes
- SMR bulk Bonferroni/FDR/core exploratory hits
- SMR cell-type Bonferroni/FDR/core exploratory hits
- 综合 `conjFDR + coloc/PWCoCo + cTWAS + SMR` 的 candidate gene master/shortlist

## 审稿与复现友好性

为了满足审稿人和后续读者复现需求，代码仓库中应同时提供：

- 入口脚本
- 运行顺序说明
- 输入格式约定
- 关键方法学判断说明
- 不同结果层级的导出策略

而不应只保留零散脚本或仅保留最终 Excel 表格。
