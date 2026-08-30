# MAINTENANCE.md — 仓库维护职责说明书

> 本仓库由 **AI 助手（WorkBuddy）担任指定维护者**。用户只负责写 R/Seurat 分析代码，git 相关操作（提交、推送、回退、worktree、大文件护栏）交给维护者执行。
> 项目根 = `A:\Workbuddy\singlecell\Rversion`，远程 = `github.com:1099072409/singlecell-Rversion`（`main` 分支，上游已绑定）。

---

## 一、职责边界

| 谁 | 负责 |
|---|---|
| **用户** | 写/改分析代码、运行分析、产出结果表；告诉我"要提交/要回退/要开 worktree" |
| **维护者（AI）** | 执行全部 git 操作、保持 `.gitignore` 护栏、把踩坑记进 `lessons/`、按需更新 `agent_doc/` |

维护者在**每次会话内**响应维护请求；不主动后台扫描仓库。用户把变化讲清楚，维护者按本文档落地。

---

## 二、日常维护流程（维护者执行）

> PowerShell 里 git 默认不在 PATH，每次新窗口先加：
> ```powershell
> $env:PATH = "C:\Users\10990\.workbuddy\binaries\PortableGit\versions\1.2.0\cmd;" + $env:PATH
> cd A:\Workbuddy\singlecell\Rversion
> ```

1. **开工前同步**（多机/多会话时）：`git pull`
2. **暂存**：`git add -A`（或只加改动的文件 `git add 01_initiation/xxx.R`）
3. **提交**：中文信息用 **Git Bash** 提交（`cd /a/Workbuddy/singlecell/Rversion` 后 `git commit -m "..."`），避免 PowerShell 下 GBK 乱码
4. **推送**：`git push`（上游已绑定 `main → origin/main`）

提交信息格式：`<type>: <描述>`，type ∈ `feat` / `fix` / `docs` / `refactor` / `chore`。

---

## 三、回退操作（做错了随时可用）

| 场景 | 命令 |
|---|---|
| 未推送，撤销最近提交**保留改动** | `git reset --soft HEAD~1` |
| 已推送，安全回退（反向提交，历史不丢） | `git revert <哈希>` |
| 单文件改乱，恢复到上次提交 | `git restore <文件>` |
| 查看历史/分支关系 | `git log --oneline` / `git log --graph --oneline -10` |
| 整体回退到某点（**危险，丢后续改动**） | `git reset --hard <哈希>` |

原则：未推送可随意 `reset`；已推送用 `revert`，**禁止** `reset --hard` 后强推（破坏远程历史）。

---

## 四、大文件护栏（最高优先级）

- `.gitignore` 已排除 `*.rds` / `*.pdf` / `*.log` / `data/`，**切勿** `git add` 大文件。
- 新增一类大产物（如 `.h5ad` / `.h5` / 新 `.rds`）→ **立即**往 `.gitignore` 加一行。
- 误 add 大文件补救：
  - 未 commit → `git reset` 取消暂存，补 `.gitignore`
  - 已 commit → `git rm --cached <文件>` + 补 `.gitignore` + 再 commit
- 每次提交前用 `git status --short` 核对无 `.rds`/`.pdf`/`.log` 混入。

---

## 五、并行 / 实验：worktree

完整流程见 `agent_doc/03_git与多agent协作/worktree并行工作流.md`。要点：

```bash
git worktree add ../.worktrees/<分支名> -b <分支名>   # 开独立工作树
# 在 worktree 里改完提交
git merge <分支名>        # 回到主仓库合并回 main
git worktree remove --force ../.worktrees/<分支名>
git branch -d <分支名>
```

合并前用目录联接（junction）共享 `../data`，**不复制**原始数据；清理 worktree 时先解联接再删目录，防误删原始数据。

---

## 六、文档同步

- 每踩一个坑 → `lessons/` 下按日期建 `.md` 或追加，AI 下次开工先读。
- 分析流程有变动 → 更新 `agent_doc/02_分析流程/` 对应文档。
- 本说明书、AGENT.md、PROJECT_STRUCTURE.md 三者保持一致。

---

## 七、新会话启动检查清单

1. 读 `AGENT.md` → `PROJECT_STRUCTURE.md` → `lessons/` 最新一条。
2. `git status --short` 确认工作区干净/有无未跟踪大文件。
3. `git log --oneline -3` 确认本地与远程同步（必要时 `git pull`）。
4. 按用户请求执行维护，结束前再次 `git status` 确认无大文件混入。
