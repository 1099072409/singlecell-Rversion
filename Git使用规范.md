# R 项目 Git 版本管理使用规范

> 适用对象：`A:\Workbuddy\singlecell\Rversion`（R / Seurat 单细胞分析流程）
> 编写日期：2026-08-30
> 环境：Windows + Git Bash，Git 2.55.0.windows.3
>
> **阅读建议**：第 1 节必读（有会踩的实际坑），第 2～4 节按顺序执行一次即可，第 5 节起是日常查阅。

---

## 1. 前置检查：先看清楚你要管的是什么

### 1.1 为什么仓库根选在 `Rversion/`，而不是 `singlecell/`

整个 `A:\Workbuddy\singlecell` 目录有 **36 GB**，构成如下：

```
A:\Workbuddy\singlecell\                  ← 36 GB，不建议整体建库
├── singlecell.Rproj                       RStudio 项目文件（在仓库外）
├── .Rhistory / .Renviron                  R 会话与环境文件（在仓库外）
├── data\                    1.6 GB        140 个 10x 原始矩阵（在仓库外）
├── pythonver\              22   GB        Python/scanpy 版本（在仓库外）
├── bioinfo_skills\          38   MB
└── Rversion\               14   GB   ★ ← 本仓库的根目录
    ├── 01_initiation\
    ├── 02_quality_control\
    ├── 03_extract_cd45\
    ├── 04_integration_clustering\
    ├── 05_celltype_annotation\
    ├── 06_differential_expression\
    ├── yijian.R
    └── Lesson-*.md
```

把仓库根设在 `Rversion/`，上层的 `data/`、`pythonver/` 就天然落在版本控制之外，**不需要**在 `.gitignore` 里为它们写任何规则（见 3.3 节的常见误解）。

### 1.2 体积分布：你的代码只有 0.4 MB

| 类型 | 数量 | 体积 | 是否入库 |
|---|---:|---:|---|
| `.R` 代码 | 14 | **0.4 MB** | ✅ 入库 |
| `.md` 文档 | 8 | ~0.2 MB | ✅ 入库（含本规范文档） |
| `.csv` 结果表 | 26 | ~0.1 MB | ✅ 入库 |
| **`.rds` Seurat 对象** | **9** | **13 GB**（单文件 248 MB – 2.4 GB） | ❌ 排除 |
| `.pdf` 结果图 | 56 | 206 MB（最大 44.1 MB） | ❌ 排除 |
| `.log` 日志 | 8 | ~0 MB | ❌ 排除 |
| `.Rhistory` | 5 | ~0 MB | ❌ 排除 |

**入仓后仓库 0.71 MB**（实测：52 个文件入库，13.06 GB 被排除）。9 个 `.rds` 全部超过 GitHub 的 100 MB 硬限制，其中 4 个在 2.3 GB 以上——这不是"要不要优化"的问题，是"不排除就根本推不上去"的问题。

### 1.3 你的 Git 环境现状（三项需要修正）

用下面命令自查：

```bash
git config --global --list
```

当前实测结果：

| 配置项 | 现状 | 影响 |
|---|---|---|
| `user.name` / `user.email` | **未配置** | 首次 commit 会报错 `Please tell me who you are` |
| `init.defaultBranch` | **未配置** | `git init` 会建出 `master`，与远程 `main` 不一致 |
| `credential.helper` | 全局未设；`helperselector.selected=<no helper>` | **凭据管理器未启用 → 每次 push 都要手输密码，且极易表现为"认证失败"** |
| `core.autocrlf` | system 级 = `true`（PortableGit 共享，勿改） | 见第 4 节 |
| `core.quotepath` | 未设置（默认 true） | 中文文件名（`细胞注释.md`）在 `git status` 显示为 `\347\273\206...` |

后三项会在第 2 节一次性修好。

### 1.4 ⚠️ PowerShell / Windows Terminal 用户必读（两种报错来源）

本指南命令默认在 **Git Bash** 中执行（WorkBuddy 自带的 Git 只挂在该 shell 的 PATH 里）。若你在 **PowerShell / Windows Terminal** 中执行，会遇到两种错误：

1. **路径写法不同**：Git Bash 用 `/a/Workbuddy/...`（盘符 + 正斜杠），PowerShell 用 `A:\Workbuddy\singlecell\Rversion`（盘符 + 反斜杠）。直接复制 `/a/...` 会被 PowerShell 当成相对路径、解析成 `C:\a\...`，于是报"找不到路径"。
2. **git 不在 PATH**：WorkBuddy 的 git（`PortableGit`）默认没进 PowerShell 的 PATH，报错"git 不是内部或外部命令"或"无法将 git 识别为 cmdlet"。

**解决办法**：第一次执行前，在 PowerShell 里先运行下面一行（仅为当前窗口生效，关闭即失效）：

```powershell
$env:PATH = "C:\Users\10990\.workbuddy\binaries\PortableGit\versions\1.2.0\cmd;" + $env:PATH
```

之后本指南所有 `git` 命令在 PowerShell 中即可直接使用。要**永久生效**，在 PowerShell 中执行一次：

```powershell
[Environment]::SetEnvironmentVariable("PATH", $env:PATH + ";C:\Users\10990\.workbuddy\binaries\PortableGit\versions\1.2.0\cmd", "User")
```

> 若你不确定当前用的是哪个 shell：看提示符。`PS C:\...>` 是 PowerShell；`MINGW64 / ~ $` 之类是 Git Bash。本指南每处涉及 `cd` 的地方都会同时给出两种写法。

---

## 2. 第一步：初始化本地仓库与身份配置

### 2.1 一次性全局配置（复制整段执行）

把尖括号内容替换成你自己的信息：

```bash
git config --global user.name  "<你的姓名>"
git config --global user.email "<你的邮箱@example.com>"

# 今后 git init 默认建 main 分支，避免 master/main 不一致
git config --global init.defaultBranch main

# 让中文文件名正常显示，而不是 \347\273\206 这类八进制转义
git config --global core.quotepath false

# 启用凭据管理器（本机当前未启用，不设这项每次 push 都要输密码）
git config --global credential.helper manager

# 推送大文件时降低超时概率（500 MB 缓冲区）
git config --global http.postBuffer 524288000
```

验证：

```bash
git config --global --get-regexp 'user\.|init\.|core\.quotepath|credential\.'
```

> **PowerShell 用户**：若还没把 git 加进 PATH（见 1.4 节），上面每条都需在前面补上 git 的完整路径，或先执行 `$env:PATH = "C:\Users\10990\.workbuddy\binaries\PortableGit\versions\1.2.0\cmd;" + $env:PATH`。若想永久生效，执行一次 `[Environment]::SetEnvironmentVariable("PATH", $env:PATH + ";C:\Users\10990\.workbuddy\binaries\PortableGit\versions\1.2.0\cmd", "User")`（注意：此命令会修改你的用户环境变量，只增不删，可随时在"系统→关于→高级系统设置→环境变量"里撤销）。

### 2.2 初始化仓库

**Git Bash：**

```bash
cd /a/Workbuddy/singlecell/Rversion
git init -b main
```

**PowerShell：**

```powershell
cd A:\Workbuddy\singlecell\Rversion
git init -b main
```

`-b main` 需要 Git ≥ 2.28，你的是 2.55，可用。

> **如果已经 `git init` 过并生成了 `master` 分支**，用这条改名：
> ```bash
> git branch -M main
> ```
> 注意是**大写 `-M`**。小写 `-m` 在目标分支名已存在时会报错 `fatal: a branch named 'main' already exists`，这是高频踩坑点。

确认状态：

```bash
git status
git branch --show-current      # 应输出 main
```

---

## 3. 第二步：`.gitignore` 设计（本项目的核心）

文件已生成好，位于 `Rversion/.gitignore`。这里说明**为什么这么写**。

### 3.1 完整规则

```gitignore
# ---------- R / RStudio 会话与用户状态文件 ----------
.Rhistory
.RData
.RDataTmp
.Rapp.history
.Rproj.user/
.Ruserdata
.Renviron          # 常含 API key、本地绝对路径，绝不入库

# ---------- R 序列化对象（本项目主要体积来源，共 13 GB）----------
*.rds
*.RDS
*.rda
*.rdata

# ---------- 图形产物（56 个 / 206 MB，均可由脚本重绘）----------
*.pdf
*.png
*.jpg
*.jpeg
*.tiff
*.tif
*.svg

# ---------- 运行日志 ----------
*.log

# ---------- 系统杂项 ----------
.DS_Store
Thumbs.db
desktop.ini

# ---------- 显式白名单（兜底，防止上面规则误伤）----------
!*.R
!*.r
!*.md
!*.csv
!*.tsv
!*.txt
!*.Rproj
!.gitignore
!.gitattributes
```

### 3.2 ⚠️ 为什么**不能**整目录忽略 `output/` 和 `exclude/`

这是本项目最容易犯的错误。看似"输出目录当然要忽略"，实际会丢东西：

**`exclude/` 名字带"排除"，但它不是被排除的文件——它是备选方案脚本目录**，里面有 6 个必须版本控制的 `.R`：

```
03_extract_cd45/exclude/03_extract_cd45_exclude.R          23.3 KB
03_extract_cd45/exclude/04_integration_clustering_exclude.R 43.1 KB
04_integration_clustering/exclude/04_integration_clustering_exclude.R 47.6 KB
04_integration_clustering/exclude/logs/verify_04_exclude.R   1.6 KB
05_celltype_annotation/exclude/05_celltype_annotation_exclude.R 44.3 KB
05_celltype_annotation/exclude/manual_markers.R              1.1 KB
```

**`output/` 里有 25 个 csv 结果表**，全部小于 16 KB，但它们是流程的**数值证据**——QC 汇总、PC 方差表、聚类汇总、细胞类型评分、样本组成。这些内容 PDF 图上看不精确，反而比图更值得入库：

```
04_integration_clustering/output/02_cluster_by_sample.csv        15.7 KB
03_extract_cd45/exclude/02_cluster_by_sample.csv                 13.5 KB
05_celltype_annotation/output/01_celltype_score_by_cluster.csv    4.6 KB
04_integration_clustering/output/01_PC_variance_table.csv         2.8 KB
...
```

另外还有一份**输入级**配置 `05_celltype_annotation/celltype_map.csv`（人工注释映射表），它虽然躺在步骤目录里，但本质是手写的输入，改错了分析就全错，必须入库。

> **结论：按扩展名忽略，不按目录忽略。** 体积大的类型（`*.rds`、`*.pdf`）整类排除，体积小的类型（`.R`、`.csv`、`.md`）无论在哪都保留。

### 3.3 ⚠️ 常见误解：`.gitignore` 管不到上层目录

`.gitignore` **只对它所在目录及其子目录生效**。上层的 `A:\Workbuddy\singlecell\data\`、`pythonver\`、`singlecell.Rproj` 在仓库根之上，`git add -A` 永远扫不到它们，无需也不应该在本文件里处理。

并且 `.gitignore` **不支持 `..` 语法**——写 `../data/` 是无效的，还容易让人误以为已经排除干净了。

### 3.4 大体积数据与输出文件的处理建议

| 内容 | 建议 | 理由 |
|---|---|---|
| 原始数据 `data/`（1.6 GB） | 留在仓库外；仓库内只保留下载脚本或 `data/README.md` 记录 GSE 编号、下载链接、校验值 | 换机器时可重新拉取；原始数据的"版本"应由 GEO 编号而非 Git 管理 |
| 中间产物 `*.rds`（13 GB） | **不入库**；由脚本重新生成 | 二进制、不可 diff、不可增量压缩；一旦入库，即便后续删除，历史中也永久留存 |
| 结果图 `*.pdf`（206 MB） | **不入库**；由脚本重绘 | 同理，且每张图都能从对应的 `.R` 精确复现 |
| 结果表 `*.csv`（<16 KB） | ✅ 入库 | 体积小、可 diff、是结论的数值依据 |
| 真正需要长期留存的个别产物 | 放对象存储 / 网盘 / 实验室 NAS，在 `output/README.md` 里记录路径与生成命令 | 兼顾可追溯与仓库轻量 |

**若确实需要把个别大对象纳入版本控制**，再考虑 Git LFS（需远程平台支持且有存储配额）：

```bash
git lfs install
git lfs track "*.rds"
git add .gitattributes
```
但对本项目的 13 GB Seurat 对象，**不推荐**——重新跑一遍脚本通常比维护 LFS 更省事。

---

## 4. 第二步附：`.gitattributes` 与换行符

文件已生成，位于 `Rversion/.gitattributes`。

**要点**：

- **R 脚本对换行符不敏感**：R 解析器同时接受 LF 与 CRLF，`source()` 无需关心。`.md`、`.csv` 同理。
- **真正敏感的是二进制**：`.rds` / `.pdf` 若被误判为文本并做换行转换，文件会直接损坏、无法读取。
- 本机 `core.autocrlf=true` 是 **system 级**（PortableGit 是共享安装，改动会影响其他项目），因此**不改全局配置**，改用仓库内的 `.gitattributes` 精确控制：

```gitattributes
* text=auto

*.R    text eol=lf
*.md   text eol=lf
*.csv  text eol=lf
...

*.rds  binary
*.rda  binary
*.pdf  binary
*.png  binary
...
```

好处：`git diff` 不会因行尾差异让整个文件飘红；二进制文件被明确标记为不可转换。

---

## 5. 第三步：首次提交

### 5.1 命令序列

```bash
cd /a/Workbuddy/singlecell/Rversion

# 1) 先加配置文件本身
git add .gitignore .gitattributes

# 2) 加入其余文件
git add -A

# 3) 提交前务必复核（关键！）
git status --short
```

### 5.2 ⚠️ 提交前的复核清单

`git status --short` 的输出里，**不应出现**任何 `.rds`、`.pdf`、`.log`、`.Rhistory`。如果出现了，说明 `.gitignore` 没生效——先回到第 3 节排查，不要硬提交。

再看一眼总体积，确保没有几十 MB 以上的文件混入：

**Git Bash：**

```bash
git ls-files | xargs -I{} du -b "{}" 2>/dev/null | awk '{s+=$1} END {printf "待入库总大小: %.2f MB\n", s/1048576}'
```

**PowerShell：**

```powershell
$total = 0; git ls-files | ForEach-Object { $total += (Get-Item $_).Length }; "待入库总大小: {0:N2} MB" -f ($total/1MB)
```

正常情况下应输出 **1 MB 左右**。若显示几百 MB 甚至上 GB，立即停止，用下面命令找出大文件：

**Git Bash：**

```bash
git ls-files | xargs -I{} du -m "{}" 2>/dev/null | sort -rn | head -10
```

**PowerShell：**

```powershell
git ls-files | ForEach-Object { [PSCustomObject]@{MB=([Math]::Round((Get-Item $_).Length/1MB,2)); File=$_} } | Sort-Object MB -Descending | Select-Object -First 10
```

### 5.3 提交

```bash
git commit -m "feat: 初始化单细胞 Seurat 分析流程（01-06 步骤脚本与结果表）"
```

确认结果：

```bash
git log --oneline --stat
```

> **⚠️ PowerShell 中中文提交信息的编码坑（重要）**
> 在 PowerShell 里直接 `git commit -m "中文"` 很容易把中文按 GBK 存进提交对象，`git log` 里会变成 `?????` 或乱码。本项目已踩过这个坑并修好了。两种可靠做法：
> 1. **推荐：在 Git Bash 里提交**（Git Bash 默认 UTF-8 locale，存储最干净）。
> 2. **坚持用 PowerShell**：先在真实环境执行一次（仅需一次，永久生效）
>    ```powershell
>    git config --global i18n.commitEncoding utf-8
>    git config --global i18n.logOutputEncoding utf-8
>    ```
>    再用**文件**传消息、并显式声明编码，避免被系统 ANSI 代码页替换成 `?`：
>    ```powershell
>    $enc = New-Object System.Text.UTF8Encoding($false)
>    [IO.File]::WriteAllText("msg.txt", "feat: 你的中文提交信息", $enc)
>    git -c i18n.commitEncoding=utf-8 commit --amend -F msg.txt
>    Remove-Item msg.txt
>    ```
>    另外，即便存储是 UTF-8，PowerShell 控制台代码页默认是 GBK，`git log` 仍可能显示乱码——这是**显示问题不是存储问题**。想正常显示，先 `chcp 65001` 切到 UTF-8 代码页，或直接在 Git Bash 看。

### 5.4 提交信息规范（Conventional Commits）

格式：

```
<type>: <简短描述，中文即可，不超过 50 字>

[可选正文：说明动机、与上一版的差异]
[可选脚注：Closes #12]
```

常用 `type`：

| type | 用途 |
|---|---|
| `feat` | 新增分析步骤、新脚本、新功能 |
| `fix` | 修正脚本 bug、修正参数 |
| `docs` | 只改文档（`.md`、注释） |
| `refactor` | 重构代码，不改分析逻辑与结果 |
| `perf` | 性能优化（如加速某步骤） |
| `chore` | 配置、`.gitignore`、依赖等杂项 |

本项目实际场景的示例：

```bash
git commit -m "feat: 新增 07 细胞通讯分析步骤"
git commit -m "fix: 修正 04 整合时 harmony 参数导致过度校正的问题"
git commit -m "docs: 补充细胞注释 marker 基因的文献依据"
git commit -m "refactor: 抽出 QC 阈值到独立配置段"
git commit -m "chore: 应用 .gitignore 规则，移除已追踪的 rds 文件"
git commit -m "feat: 更新 05 细胞类型注释，重跑注释后更新结果表

调整了 T 细胞亚群的 marker 权重，NK 细胞比例由 8.2% 变为 6.7%。"
```

**中英混写建议**：`type` 用英文（便于工具解析与筛选），描述用中文（便于团队阅读）。

---

## 6. 第四步：连接远程仓库

以下命令**平台无关**，把 `<REMOTE_URL>` 换成你的仓库地址即可。

### 6.1 场景 A：远程是全新空仓库（无 README、无 LICENSE）

这是最干净的情况：

```bash
git remote add origin <REMOTE_URL>
git push -u origin main
```

### 6.2 场景 B：远程已有内容（带 README / LICENSE / .gitignore 初始提交）

先拉取再变基，避免产生一次无意义的合并提交：

```bash
git remote add origin <REMOTE_URL>
git fetch origin
git rebase origin/main
git push -u origin main
```

如果 `git rebase` 报 `fatal: refusing to merge unrelated histories`（两边各有独立的初始提交，这是**预期行为**，不是错误）：

```bash
git pull origin main --allow-unrelated-histories
# 解决完冲突后：
git add -A
git commit -m "chore: 合并远程初始提交"
git push -u origin main
```

### 6.3 场景 C：本地已经建出 master（本仓库最可能遇到的情况）

如果你忘了加 `-b main`，本地就是 `master`：

```bash
git branch -M main                 # 大写 -M，强制重命名
git push -u origin main
```

若远程**已经存在** `main` 分支且有内容，则走场景 B 的流程：

```bash
git branch -M main
git fetch origin
git rebase origin/main
git push -u origin main
```

### 6.4 关于 `-u`（上游分支）

`git push -u origin main` 中的 `-u` 是 `--set-upstream` 的简写，作用是**把本地 `main` 绑定到远程 `origin/main`**。绑定后，后续只需：

```bash
git push      # 等价 git push origin main
git pull      # 等价 git pull origin main
```

三种等价写法：

```bash
git push -u origin main
git push --set-upstream origin main
git push origin main && git branch --set-upstream-to=origin/main main
```

**已经推送过但忘记绑定上游**时补绑：

```bash
git branch --set-upstream-to=origin/main main
```

验证：

```bash
git remote -v        # 查看远程地址
git branch -vv       # 查看分支与上游的绑定关系，应显示 [origin/main]
```

### 6.5 修改或更换远程地址

```bash
git remote set-url origin <NEW_REMOTE_URL>     # 改地址
git remote remove origin                        # 删除后重新添加
```

---

## 7. 第五步：认证配置（平台差异集中在这一节）

两条路任选其一。**推荐 HTTPS + 凭据管理器**，配置最简单、穿透性最好。

### 7.1 HTTPS + Personal Access Token

1. 在远程平台生成 token（各平台叫法见 7.3）
2. 执行 `git push`，弹窗或命令行提示时：
   - **用户名** = 平台用户名（**不是邮箱**，Gitee 尤其注意）
   - **密码** = 刚才生成的 token（不是账号密码）
3. 由 Git Credential Manager（GCM）存入 Windows 凭据管理器，之后免密

如果没弹出凭据窗口、或每次都索要密码，检查第 2.1 节那条配置是否已执行：

```bash
git config --global credential.helper manager
```

**手动清理旧凭据**（token 更换后必须做，否则一直用旧的失败凭据）：

```bash
git credential-manager erase
```
按提示输入 `protocol=https`、`host=github.com` 后连按两次回车；
或直接在 **控制面板 → 凭据管理器 → Windows 凭据** 里删除 `git:https://xxx` 条目。

### 7.2 SSH 密钥

```bash
# 生成密钥（邮箱换成你的）
ssh-keygen -t ed25519 -C "<你的邮箱@example.com>"

# 复制公钥内容，粘贴到平台的 SSH Keys 设置页
cat ~/.ssh/id_ed25519.pub

# 测试连通性
ssh -T git@github.com
```

使用 SSH 地址：

```bash
git remote set-url origin git@github.com:<USERNAME>/<REPO>.git
```

### 7.3 各平台差异对照

| 平台 | HTTPS 凭据名称 | 生成入口 | 最小必要权限 | 单文件上限 |
|---|---|---|---|---|
| **GitHub** | fine-grained PAT（推荐）/ classic PAT | Settings → Developer settings → Personal access tokens | fine-grained：`Contents: Read and write`（可限定到单个仓库，最长有效期 366 天）；classic：勾选 `repo` | **100 MB 硬拒**（>50 MB 警告）；仓库软建议 <1 GB |
| **GitLab.com** | Project access token（推荐）/ Personal access token（前缀 `glpat-`） | 项目 → Settings → Access tokens | `write_repository` | 仓库含 LFS 10 GB；单次 push 上限 5 GiB |
| **Gitee** | 私人令牌 | 设置 → 安全设置 → 私人令牌 | 勾选 `projects` 相关写权限 | 随套餐变动，平台侧可配 `max_file_size`，**需自行核实** |
| **Bitbucket** | App password | Atlassian account → Security → App passwords | `repository: write` | **需自行核实** |

**平台无关的通用注意点**：

- 所有平台的 token **只在生成时完整显示一次**，务必当场保存。
- **不要把 token 拼进 remote URL**（`https://<TOKEN>@github.com/...`），它会明文落进 `.git/config`，一旦提交或分享即泄露。用 GCM 存凭据。
- 100 MB 单文件限制是行业通行的硬门槛，本项目的 9 个 `.rds` 全部超标。

---

## 8. 第六步：大文件已经误入库了怎么办

如果 `.gitignore` 配置晚了一步，`*.rds` 已经进了提交历史，从工作区删除**没有用**——历史里依然留着，clone 时照样下载。

需要先从历史中彻底剥离（推荐 `git filter-repo`，比 BFG 与 `filter-branch` 更快更安全）：

```bash
# 安装（需要 Python）
pip install git-filter-repo

# 剥离所有 rds 文件（会重写历史）
git filter-repo --path-glob '*.rds' --invert-paths

# 强制推送（仅限个人仓库或已通知协作者）
git push --force-with-lease origin main
```

⚠️ **重写历史是破坏性操作**：所有协作者必须重新 clone。能用 `--force-with-lease` 就不要用 `--force`（前者会在远程有你未知的新提交时拒绝执行，多一道保险）。

**最佳策略仍是预防**：第一次提交前就用第 5.2 节的方法复核体积。

---

## 9. 第七步：日常使用约定

### 9.1 提交粒度

- **一个分析步骤一次提交**：脚本与其产生的结果 csv **同批次提交**，保证"代码版本"与"结果版本"严格对应。
- **不要把一周的改动攒成一次大提交**。攒得越久，出问题时越难定位到具体改动。
- **一个提交只做一件事**。改了 QC 阈值又顺手调了配色，应拆成两次提交。
- **不要把未完成、跑不通的脚本提交到 `main`**——用特性分支（9.3）。

### 9.2 推送前的拉取与变基（每日开工必做）

```bash
# 开工时：先同步
git pull --rebase origin main

# 做完改动
git add -A
git commit -m "feat: ..."

# 推送前：再次同步，避免冲突
git pull --rebase origin main

# 推送
git push
```

用 `--rebase` 而不是默认的 merge，可以让提交历史保持一条直线，不产生 `Merge branch 'main' into main` 这类噪音提交。

设置为默认行为，省得每次都敲：

```bash
git config --global pull.rebase true
```

遇到冲突时：

```bash
# 编辑冲突文件（搜索 <<<<<<< 标记），解决后：
git add <冲突文件>
git rebase --continue

# 想放弃这次变基：
git rebase --abort
```

### 9.3 分支策略（特性分支）

单人或小团队，用最简模型即可：

```
main          ← 受保护，始终可运行、可复现
 └── feature/<步骤>-<内容>
 └── fix/<步骤>-<问题>
```

创建与合并：

```bash
# 从 main 切出特性分支
git switch -c feature/05-marker-refine

# 开发、提交若干次...
git add -A
git commit -m "feat: 调整 T 细胞亚群 marker 权重"

# 合并回 main（先确保 main 是最新的）
git switch main
git pull --rebase origin main
git merge --no-ff feature/05-marker-refine -m "merge: 合并 marker 权重调整"

# 推送并清理
git push
git branch -d feature/05-marker-refine
```

命名建议（`05` 对应步骤目录编号，便于一眼看出改的是哪一步）：

- `feature/05-marker-refine`
- `feature/07-cellchat`
- `fix/04-harmony-params`

### 9.4 几条硬性纪律

1. **绝不 `git push --force` 到 `main`**（真需要时用 `--force-with-lease`）
2. **提交前必看 `git status --short`**，确认没有 `.rds` / `.pdf` 混入
3. **不提交 `.Renviron`、token、密码、绝对路径**
4. **每次跑完分析，csv 结果表随脚本一起提交**——它们是你的实验记录

---

## 附录 A：常见错误排查表

| 症状 / 报错 | 原因 | 解决办法 |
|---|---|---|
| `Please tell me who you are` | 未配置 user.name / user.email | 执行第 2.1 节的 `git config --global user.name/email` |
| `fatal: Authentication failed` / `403` | ① token 过期 ② token 权限（scope）不足 ③ 凭据管理器未启用，把账号密码当 token 提交 | 重新生成 token 并勾选写权限；`git config --global credential.helper manager`；`git credential-manager erase` 清旧凭据后重推 |
| `401 Unauthorized` | 用户名填错（用了邮箱）或 token 无效 | 用户名必须是**平台用户名**；重新生成 token |
| `fatal: refusing to merge unrelated histories` | 本地与远程各有独立的初始提交 | `git pull origin main --allow-unrelated-histories` |
| `error: src refspec main does not match any` | 本地分支还叫 `master`，或一次提交都还没有 | `git branch -M main`；或先 `git commit` |
| `! [rejected] (non-fast-forward)` | 远程有新提交，本地未同步 | `git pull --rebase origin main` 后再 push **（不要用 `--force`）** |
| `fatal: a branch named 'main' already exists` | 用了小写 `-m` 重命名 | 改用大写 `-M`：`git branch -M main` |
| `.gitignore` 加了，文件仍在被追踪 | 规则只对**未追踪**文件生效，已入库的不受影响 | `git rm --cached <file>`；整类清理：`git rm -r --cached . && git add -A && git commit -m "chore: 应用 .gitignore 规则"` |
| `remote: error: File xxx is 2.30 GB; this exceeds GitHub's file size limit of 100 MB` | `.rds` 超过平台硬限制 | ① `git filter-repo --path-glob '*.rds' --invert-paths` 剥离历史；② 或改用 LFS；③ 本项目推荐**两者都不用**——仓库只留代码与结果表 |
| push 大文件超时 / 连接中断 | 单次推送体积或网络超时 | `git config --global http.postBuffer 524288000`；仍失败则改用 SSH |
| `Filename too long` | Windows 路径长度限制 | `git config --global core.longpaths true` |
| 中文文件名显示为 `\344\275\240\345\245\275` | `core.quotepath` 默认为 true（只是显示问题，文件本身正常） | `git config --global core.quotepath false` |
| 提交后 `git status` 仍显示文件被修改（行尾） | CRLF/LF 转换 | 确认 `.gitattributes` 已提交；`git add --renormalize .` |
| `CONFLICT (content)` during rebase | 与远程改了同一处 | 编辑冲突文件 → `git add` → `git rebase --continue`；放弃用 `git rebase --abort` |

---

## 附录 B：命令速查卡

```bash
# ===== 每日流程 =====
git pull --rebase origin main          # 开工同步
git status --short                     # 看改了什么
git diff                               # 看具体改动
git add -A                             # 暂存全部
git commit -m "feat: 描述"              # 提交
git pull --rebase origin main          # 推送前再同步
git push                               # 推送

# ===== 分支 =====
git branch                             # 列出本地分支
git switch -c feature/xxx              # 新建并切换
git switch main                        # 切回主干
git branch -d feature/xxx              # 删除已合并分支

# ===== 撤销（未推送）=====
git restore <file>                     # 丢弃工作区改动
git restore --staged <file>            # 取消暂存
git commit --amend -m "新信息"           # 修改上一次提交信息

# ===== 撤销（已推送，慎用）=====
git revert <commit-id>                 # 生成一个反向提交，安全
git reset --hard <commit-id>           # 回退并丢弃，危险

# ===== 查看 =====
git log --oneline -20                  # 最近 20 条提交
git log --oneline --stat               # 带改动文件统计
git show <commit-id>                   # 看某次提交的具体内容
git remote -v                          # 远程地址
git branch -vv                         # 分支与上游绑定关系
```

---

## 附录 C：本仓库资产清单

**入库（实测 52 个文件 / 0.71 MB，首次提交 `eaf78dc`→经 UTF-8 修正后为最终提交）**

| 类型 | 数量 | 说明 |
|---|---:|---|
| `.R` | 14 | 含 `exclude/` 下 6 个备选方案脚本；`yijian.R` 在根目录；`00test/manual_markers.R` |
| `.md` | 8 | 6 个 `Lesson-*.md` + `05_celltype_annotation/细胞注释.md` + 本规范文档 |
| `.csv` | 26 | 25 个在 `output/`、`exclude/` 下的结果表 + `celltype_map.csv` 人工注释映射表 + `00test/celltype_map.csv` |
| 配置文件 | 2 | `.gitignore`、`.gitattributes` |

> 注：`00test/` 是项目里的一个小型验证目录（1 个 `.R` + 1 个 `.csv`），已一并入库；如不需要可随时 `git rm -r --cached 00test`。另外本仓库已设仓库级 `core.autocrlf=false`，换行符完全交给 `.gitattributes` 的 `eol=lf` 控制，避免 `system` 级 `autocrlf=true` 在 Windows 上带来的 CRLF 噪音。

> 以上已用脚本按 gitignore 语义模拟核对：`.R` / `.md` / `.csv` 误伤数为 0，入库最大单文件仅 0.12 MB（`yijian.R`）。

**排除**

| 类型 | 数量 | 体积 | 排除理由 |
|---|---:|---:|---|
| `.rds` | 9 | 13 GB | 二进制、不可 diff、全部超平台 100 MB 限制 |
| `.pdf` | 56 | 206 MB | 可由脚本重绘 |
| `.log` | 8 | ~0 MB | 运行时日志 |
| `.Rhistory` | 5 | ~0 MB | R 会话历史 |

**在仓库外（位于 `Rversion/` 之上，不受 Git 管辖）**

- `A:\Workbuddy\singlecell\data\` — 1.6 GB 原始 10x 矩阵
- `A:\Workbuddy\singlecell\pythonver\` — 22 GB Python/scanpy 版本
- `A:\Workbuddy\singlecell\singlecell.Rproj`、`.Renviron`
