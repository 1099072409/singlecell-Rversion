# 项目结构描述（Rversion）

23 例样本（ypN0 / ypN+）10x 单细胞 RNA-seq 的 R/Seurat 分析工程。

## 目录树

```
Rversion/
├── AGENT.md                      # 项目协作规范（AI 必读）
├── PROJECT_STRUCTURE.md          # 本文件：项目结构描述
├── Git使用规范.md                # Git 日常操作手册
├── .gitignore                    # 大文件排除规则
├── .gitattributes                # 换行符与二进制控制
│
├── 01_initiation/                # 步骤1：数据初始化
│   ├── 01_initiation.R           #   脚本
│   └── output/                   #   中间产物（.rds，不入库）
├── 02_quality_control/           # 步骤2：质量控制
│   ├── 02_quality_control.R
│   └── output/
├── 03_extract_cd45/              # 步骤3：提取 CD45+ 细胞
│   ├── 03_extract_cd45.R
│   ├── exclude/                  #   备选方案脚本（.R，入库）
│   ├── logs/
│   └── output/
├── 04_integration_clustering/    # 步骤4：整合与聚类
│   ├── 04_integration_clustering.R
│   ├── exclude/
│   └── output/
├── 05_celltype_annotation/       # 步骤5：细胞类型注释
│   ├── 05_celltype_annotation.R
│   ├── manual_markers.R          #   人工注释（入库）
│   ├── celltype_map.csv          #   注释映射表（入库）
│   ├── 细胞注释.md               #   注释说明（入库）
│   ├── exclude/ logs/ output/
├── 06_differential_expression/   # 步骤6：差异表达
│   └── 06_differential_expression.R
│
├── 00test/                       # 合并测试脚本
│   ├── merged_pipeline_01_05.R   #   01~05 合并版（入库）
│   ├── README_merged.md
│   ├── manual_markers.R / celltype_map.csv
│   └── output/ sweep_output/     #   产物（不入库）
│
├── agent_doc/                    # 各环节注意事项（分子目录）
│   ├── 01_数据与样本对齐/
│   ├── 02_分析流程/
│   ├── 03_git与多agent协作/
│   └── 04_环境与依赖/
│
├── lessons/                      # 以往踩坑沉淀
│
└── Lesson-*.md                   # 学习笔记
    yijian.R                      # 一键流程脚本
```

## 入库 vs 不入库

| 类别 | 扩展名/位置 | 入库 |
|---|---|---|
| R 代码 | `*.R` | ✅ |
| 文档 | `*.md` | ✅ |
| 结果表 | `*.csv`（含 `celltype_map.csv`、`output/` 下小表） | ✅ |
| 配置文件 | `.gitignore` `.gitattributes` | ✅ |
| 中间对象 | `*.rds`（Seurat 对象，单文件最大 2.4GB） | ❌ |
| 图 | `*.pdf` `*.png` | ❌（脚本可重绘） |
| 日志/会话 | `*.log` `.Rhistory` `.RData` | ❌ |

> 完整规则见 `.gitignore`。

## 共享数据（不进本仓库）
- 原始 10x 矩阵：上级目录 `../data/`（约 1.6GB）
- Python 版分析：上级目录 `../pythonver/`（约 22GB）
- 多工作树共享数据方式见 `agent_doc/03_git与多agent协作/`。
