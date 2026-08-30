# AGENT.md — 单细胞分析项目（R/Seurat）协作规范

> 23 例样本（ypN0 vs ypN+）10x Genomics 单细胞 RNA-seq 分析，R/Seurat 工具链。
> 项目根 = 本目录（`Rversion`）。
> **维护者**：AI 助手（WorkBuddy）担任指定维护者，git 操作全部由其执行（见 `MAINTENANCE.md`）。

## 必读顺序
1. **项目结构** → `PROJECT_STRUCTURE.md`
2. **各环节注意事项** → `agent_doc/`（按子目录分类，动手前查对应情况）
3. **以往踩坑** → `lessons/`（动手前先扫一眼，别重蹈覆辙）

## 铁律（违反会出事）
1. **大文件永不入库**：`*.rds` / `*.pdf` / `*.log` 由 `.gitignore` 排除；只提交 `.R` 代码、`.csv` 结果表、`.md` 文档。
2. **数据与代码分离**：原始矩阵在上级 `../data/`，中间 `.rds` 在 `output/`，都不进 git。
3. **提交信息**：`<type>: <描述>`，type ∈ feat / fix / docs / refactor / chore。中文描述务必 UTF-8（Windows 下优先在 Git Bash 提交）。
4. **一个阶段一提交**：做完一步就 commit，出错可随时回退（见 `agent_doc/03_git与多agent协作/`）。
5. **并行任务用 worktree**：多任务/多 agent 并行时用 `git worktree`，别在 main 上直接改。
6. **新踩坑即沉淀**：任何新错误/新坑，记入 `lessons/`。

## 目录速览
| 目录/文件 | 作用 |
|---|---|
| `00test/` ~ `06_*/` | 分析步骤（脚本 + `output/` + `exclude/`） |
| `agent_doc/` | 不同情况的注意事项（分子目录） |
| `lessons/` | 以往犯过的错误沉淀 |
| `PROJECT_STRUCTURE.md` | 完整项目结构描述 |
| `MAINTENANCE.md` | 维护职责说明书（AI 为指定维护者） |
| `Git使用规范.md` | Git 日常操作手册 |
