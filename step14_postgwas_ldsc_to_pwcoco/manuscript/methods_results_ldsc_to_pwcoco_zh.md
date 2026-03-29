# 方法

为系统评估脑类淋巴相关潜变量与神经退行性疾病之间的共享遗传结构，我们围绕两个因子 GWAS 表型（F1 和 F2）以及三种神经退行性疾病 GWAS（阿尔茨海默病，AD；帕金森病，PD；路易体痴呆，LBD）开展了多层次跨性状分析。总体分析框架依次包括全局遗传相关分析、双变量混合模型分析、跨性状条件富集与 pleiotropic locus 鉴定，以及区域层面的共定位与条件共定位验证。为评估潜变量相较于原始 ALPS 指标是否具有方法学优势，我们还将 5 个原始 ALPS 表型及其汇总指标（Mean_ALPS、tALPS 和 mALPS）纳入扩展比较分析。

全局遗传相关分析采用 linkage disequilibrium score regression（LDSC）完成。我们首先对 6 个目标 pair（F1/F2 × AD/PD/LBD）估计 SNP 遗传率和遗传相关系数，并将同一流程扩展到 F1/F2、5 个原始 ALPS 指标及 Mean_ALPS/tALPS/mALPS 与三种神经退行性疾病之间的全部配对，以比较潜变量与原始表型在疾病共享遗传结构上的差异。多重检验校正采用 Benjamini-Hochberg 方法。

随后，我们采用 MiXeR 对 F1/F2 与 AD/PD/LBD 的双变量遗传结构进行建模。MiXeR 用于估计各性状的 polygenicity、可检出位点数量及两性状之间的 polygenic overlap。对于病例-对照性状，样本量参数采用有效样本量而非总样本量，以符合 MiXeR 官方推荐。分析按官方 real-data workflow 进行，包括单变量拟合和双变量拟合，并提取遗传相关、共享多基因成分比例及共享可检出位点数量等指标，用于补充 LDSC 对全局相关结构的刻画。

为从位点层面识别共享关联区域，我们进一步应用 pleioFDR 框架对所有 pair 进行 conditional QQ 和 conjunctional false discovery rate（conjFDR）分析。主分析针对 6 个目标 pair（F1/F2 × AD/PD/LBD）开展；扩展分析则将 5 个原始 ALPS 指标及 Mean_ALPS/tALPS/mALPS 与三种疾病的所有 pair 一并纳入，以比较因子表型与原始 ALPS 表型在共享 SNP 和独立 loci 发现上的差异。显著阈值设为 conjFDR < 0.05。对于 conjFDR 结果，我们对 lead SNP 和独立 locus 进行整理，并基于去重后的 sentinel SNP 构建候选共定位区域。

区域层面的共定位分析首先采用 classical coloc 完成。我们以 conjFDR 鉴定到的 factor-level 共享 loci 为起点，将相邻 sentinel SNP 按 ±500 kb 合并为区域级分析窗口，并在每个区域内提取原始 summary statistics 的完整 SNP 集进行共定位分析。coloc 分析输出后验概率 PP.H0–PP.H4，其中 PP.H4 表示两种表型更可能由同一潜在致病变异驱动。由于 classical coloc 基于单因果变异假设，对于可能存在多个独立信号的区域，我们进一步采用 PWCoCo 进行条件共定位分析。PWCoCo 结合 stepwise 条件分析与 coloc 框架，可在复杂 LD 区域中区分“共享关联”与“共享因果变异”。参考 LD 面板采用 1000 Genomes Project Phase 3 欧洲人群 PLINK 数据。考虑到部分 factor GWAS 信号弱于传统 genome-wide significance 阈值，我们对若干高优先级区域实施了预设的 sensitivity analysis，通过放宽 factor 侧和/或疾病侧 stepwise selection 的 P 值阈值，以评估高 H4 信号是否仅由默认选择阈值所掩盖。

# 结果

## F2 与 AD、F1 与 PD 表现出最明确的全局共享遗传结构

LDSC 分析显示，在 6 个主要 pair 中，F2 与 AD 的遗传相关最强，表现为正向相关趋势；相比之下，F1 与 PD 呈负向相关趋势，且方向上较 F2 与 PD 更为一致。LBD 相关结果整体不稳定，其主要原因在于 LBD 自身 SNP 遗传率估计精度较低，因此后续所有涉及 LBD 的结果均仅作为探索性证据解读。将分析扩展至 5 个原始 ALPS 指标及 Mean_ALPS/tALPS/mALPS 后，整体模式进一步得到支持：AD 相关信号主要集中于 F2、tALPS 和 mALPS 所代表的维度，而 PD 则更偏向 F1 及 mean-like/侧化相关指标。这表明潜变量并非简单重复单一 ALPS 指标，而是在疾病相关性上对原始表型进行了更高层次的组织与提炼。

## MiXeR 进一步支持 AD 更偏向 F2、PD 更偏向 F1 的多基因重叠模式

在 F1/F2 与三种神经退行性疾病的双变量 MiXeR 分析中，F2–AD 为最强组合，其双变量遗传相关高于 F1–AD，同时共享多基因成分的比例也更高，提示 AD 的 polygenic architecture 中有较大一部分与 F2 重叠。相比之下，F1–AD 的相关和 overlap 均明显较弱。对于 PD，两组结果均表现为负向，但 F1–PD 的负向相关强于 F2–PD，支持 PD 更偏向 F1 支路。LBD 与 F1 或 F2 的双变量结果均较弱，且共享可检出位点数量有限，与 LDSC 提示的低稳定性一致。整体而言，MiXeR 在“多基因重叠”这一层面强化了 LDSC 的结论，即 AD 更接近 F2 所代表的遗传维度，而 PD 更接近 F1 所代表的遗传维度。

## conjFDR 揭示了显著的共享关联区域，但其优势在不同疾病和不同 ALPS 表型之间具有明显异质性

在 factor-level conjFDR 分析中，F2–AD 检测到最多的共享 SNP 之一，并形成少数但集中的独立 loci；F2–PD 也鉴定出大量共享 SNP，而 F1–PD 则提供了数量较少但可重复的共享 loci。相较之下，F1–AD 的 shared signal 较少，而 F1/F2 与 LBD 在该阈值下几乎未能形成稳定的共享 loci。值得注意的是，将原始 ALPS 表型与汇总 ALPS 指标纳入同一套 conjFDR 比较框架后，我们发现潜变量的优势并不简单表现为“发现更多 SNP”。在 AD 这条线上，若干原始 ALPS 或汇总 ALPS 指标在 SNP 数量上可超过 F2，但 F2 的信号通常更集中于较少的独立 loci，提示其可能更适合提炼结构化、疾病相关的共享遗传模式。相反，在 PD 这条线上，一些原始 ALPS 表型能够检测到比 F1 更多的 conjFDR hits，说明潜变量的主要价值可能在于对共享架构进行压缩和概括，而不一定在每一种 cross-trait 方法下都产生最多的显著位点。

## classical coloc 在多个候选区域提示潜在共定位，但这些信号在条件共定位框架下大多未能保持

基于 factor-level conjFDR 结果，我们共整理出 21 个去重后的 sentinel loci，并进一步合并为 18 个区域级共定位窗口，共形成 22 个实际的 region-pair coloc 任务。classical coloc 显示，多个区域具有较高的 PP.H4，其中最突出的包括 F1–PD 在 chr12:49.5–50.5 Mb（R12）和 chr2:136.2–137.2 Mb（R03）区域，以及 F1–AD 在 chr1:200.5–201.5 Mb（R01）区域。这些结果提示若干 loci 可能存在 shared causal signal。然而，coloc 采用单因果变异假设，因此我们进一步选取最关键及最复杂的区域开展 PWCoCo 条件共定位验证。

## PWCoCo 表明已验证的重点区域大多更符合复杂多信号区或独立信号区，而非稳定的共享因果位点

在 chr17 复杂区域（R17）中，未条件化 coloc 对 F1–AD 给出了很高的 H4，但在 PWCoCo 中，一旦纳入疾病侧或双侧条件分析，H4 显著下降。对 F1–AD 进行 factor 侧阈值放宽的 sensitivity analysis 后，F1 lead SNP 成功进入 stepwise selection，但双侧条件后 H4 仍接近于零，表明该结果并非由默认阈值过严所致，而更可能反映该区域存在多个独立信号。类似地，F2–AD、F1–PD 和 F2–PD 在 R17 中也均未获得稳定的条件后高 H4 证据。另一个优先验证的 chr18 区域（R18）在未条件化时即主要表现为高 H3，条件分析后仍不支持共享因果变异，提示该区域更符合“同一共享关联区域内存在彼此独立的疾病和因子信号”这一模型。

我们随后将 PWCoCo 扩展到若干在 classical coloc 中 H4 较高且结构相对更简单的区域。R12 和 R03 在未条件化时分别表现出极高或较高的 H4，但在对 factor 侧实施 sensitivity analysis 使其进入 stepwise selection 后，条件后 H4 均明显下降，未能支持稳定的 shared causal variant。R01 呈现出相似模式，即便在同时放宽 F1 和 AD 侧的 selection threshold 后，双侧条件结果仍不支持共同因果驱动。R06 亦表现为未条件化 H4 较高，但在让 PD 侧进入条件模型后，双侧条件 H4 降至接近零。R15 是本轮验证中最接近保留共定位信号的区域：在仅对 PD 条件于一个 lead SNP 时 H4 仍保持较高，但一旦进一步允许 F1 侧进入条件模型，该高 H4 信号同样消失。综合这些结果，我们未能在已重点验证的 loci 中获得稳定支持 shared causal variant 的证据。

## 综合结论

综合来看，F1/F2 与 AD/PD 之间存在明确的共享遗传结构和共享关联区域，但这些共享区域多数不对应于简单的一对一共享因果变异。LDSC 和 MiXeR 从全局水平一致支持 AD 与 F2、PD 与 F1 之间存在疾病特异性的共享遗传结构；conjFDR 进一步在位点层面识别出多个共享关联区域；而 classical coloc 与 PWCoCo 的结合则提示，这些区域中的相当一部分更可能代表复杂 LD 背景下的多信号结构，而非由单一 shared causal variant 驱动。因此，本研究更稳妥的结论是：脑类淋巴相关遗传维度与神经退行性疾病之间存在显著但复杂的 shared genetic architecture，且这种共享具有明显的疾病异质性，其中 AD 更偏向 F2，PD 更偏向 F1，而 LBD 的证据目前仍有限。
