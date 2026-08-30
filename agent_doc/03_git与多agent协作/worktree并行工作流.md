# Git 版本管理 + worktree 并行工作流

> 目标：阶段性成果可保存、可回退；多 agent 并行互不干扰，完成后合并。

## 一、版本管理：把「阶段」存下来

### 什么时候提交
- 每完成一个可独立运行的阶段（一个脚本、一个结果表、一次 bug 修复）就 commit 一次。
- 提交信息：`<type>: <描述>`，type ∈ feat / fix / docs / refactor / chore。
  例：`feat: 完成 CD45+ 细胞提取脚本`、`fix: 修正质控线粒体阈值`。

### 回退（AI 做错了怎么办）
```bash
git log --oneline -10          # 看最近提交
git checkout -- <文件>          # 丢弃工作区某个文件未提交的改动
git reset --soft HEAD~1        # 撤销最近一次提交，改动保留在工作区
git reset --hard <commit-id>   # 彻底回到某提交（丢弃之后所有改动，慎用）
```

## 二、worktree：多 agent 并行

### 概念
- 一个仓库可同时有多个「工作树」（worktree），每个挂在不同分支上。
- 好处：agent A 在分支 `feat/A`，agent B 在 `feat/B`，**互不踩对方的文件**，共用同一 `.git` 对象库（不重复占空间）。

### 创建 worktree
```bash
cd /a/Workbuddy/singlecell/Rversion      # 主仓库（Git Bash）
git worktree add ../.worktrees/agent-A -b feat/agent-A
git worktree add ../.worktrees/agent-B -b feat/agent-B
```
- `../.worktrees/agent-A` 是工作树路径（在仓库外，不会被主仓库跟踪）。
- `-b feat/agent-A` 为它新建分支。

### 共享数据：目录联接（junction）代替复制
- 原始数据（`../data`，约 1.6GB）与中间 `.rds` 不进 git，但 agent 需要访问。
- 在每个 worktree 里建一个**目录联接**指向共享数据，省去复制：
```powershell
New-Item -ItemType Junction -Path "A:\Workbuddy\singlecell\.worktrees\agent-A\data" -Target "A:\Workbuddy\singlecell\data"
```
- 效果：worktree 里的 `data/` 与主数据是**同一份**，改一处处处同步，不额外占空间。

### 并行 → 合并
```bash
# 各 agent 在自己的 worktree 里正常提交
cd /a/Workbuddy/singlecell/.worktrees/agent-A
git add -A && git commit -m "feat: A 的产出"

# 回到主仓库，逐个合并
cd /a/Workbuddy/singlecell/Rversion
git merge feat/agent-A
git merge feat/agent-B
```

### 收尾
```bash
git worktree remove ../.worktrees/agent-A
git branch -d feat/agent-A
```

### 冲突处理
- 两个 agent 改了同一个文件 → merge 报冲突。`git status` 定位，手动改好后 `git add` + `git commit`。
- 原则：**让 agent 各自改不同的文件/目录**，从源头避免冲突。

## 三、查看所有 worktree / 分支
```bash
git worktree list              # 列出所有工作树
git branch -a                  # 列出所有分支
git log --graph --oneline --all # 分支合并图
```
