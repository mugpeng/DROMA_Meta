  1. Meta_project/workflow/00-Workflow_Common.R
     定义公共配置和读写逻辑。会定位仓库根目录，加载 Meta_project/DROMA_R2/R 里的 helper，并设定本次分析参数：
     drug_names = "Paclitaxel"、tumor_types = "breast cancer"、feature_type = "mRNA"、es_t = 0.1、fdr_t = 0.1。
     同时把模型分成两组：
     cellline = CellLine + PDC
     pdcpdx = PDO + PDX
  2. Meta_project/workflow/01-Setup_And_Project_Grouping.R
     连 droma.sqlite，读所有 project 注释，然后按 dataset_type 分成 cellline 和 pdcpdx 两大组。
     还会为每组构建一个 MultiDromaSet，供后面统一取表达和药敏数据。

  3. Meta_project/workflow/02-Stage1_Coverage_Filter.R
     做 coverage / overlap 初筛。
     对每个 project 检查：某 drug、某 tumor type 下，同时有表达数据和药敏数据的样本数够不够。
     阈值在 helper 里定义：
     cellline 每个 study 至少 overlap 20，且至少 3 个 study 过线；
     pdcpdx 每个 study 至少 overlap 10，且至少 2 个 study 过线。
     这个study 定义为两两组合，类似batchFindSignificantFeatures 的思路（比如ccle exp 对ctrp2 的drug）
     这一步的作用是先确认“这个药物-癌种组合是否有足够研究可做 meta”。
  4. Meta_project/workflow/03-Stage2_WithinStudy_Screen.R
  先把cellline 和 pdcpdx 的交集基因提取出来。
    然后看每个基因在各个pair study 内absR 是不是大于0.2 且 p < .05，只要满足一个pair study 显著，则通过。
  5. Meta_project/workflow/04-Stage3_Group_Meta.R
     上一步通过的基因，全部按照batchFindSignificantFeatures 分别跑cellline 和pdcpdx，通过上述参数阈值后，再取交集。
  6. Meta_project/workflow/05-Stage4_TCGA_Translation_Filter.R
     把前临床候选拿去跟 TCGA 表达分布做转译一致性检查。
     逻辑是：
     把模型系统里同癌种该基因的表达值汇总起来，
     再读取对应 TCGA 队列表达，
     对同一基因，在tcga 和对应cellline 和 pdcpdx 分别跑 Anderson-Darling，只要在其中一个fdr < .01 即通过。
  7. Meta_project/workflow/06-Stage5_Preclinical_Merge.R
     会标记：
     invitro_supported
     pdcpdx_supported
     tcga_ad_cellline_fdr 
     tcga_ad_pdcpdx_fdr
     direction_concordant
     并生成综合 effect_size
     这一步的目标是形成统一的前临床候选列表。
     选择 direction_concordant 的结果
  8. Meta_project/workflow/07-Stage6_CDRTB_Clinical_Validation.R
     用 ctrdb.sqlite 做临床验证。
     这里优先找对应tumor type的，如果没有对应的，再用全部数据，不过最后需要标记。
     对前临床候选基因调用 batchFindClinicalSigResponse()，找临床 drug response 关联，
     再用 getClinicalCandidateFeatures() 按阈值筛临床支持信号。
     最后打上：
     clinical_supported
     direction_concordant
     retained
     其中 retained = clinical_supported 且 临床方向与前临床方向一致。
