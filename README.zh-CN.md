# CellRoB

CellRoB 用于单细胞染色质接触数据中的边界检测与跨细胞对齐。安装和输入格式见 [英文 README](README.md)，算法步骤见 [中文方法说明](docs/METHOD.md)。

## 主要文件

- `src/cellrob/core.py`：完整核心算法，原类名 RobHiC 改为 CellRoB。
- `src/cellrob/__main__.py`：通用命令行入口，替代依赖服务器绝对路径的脚本。
- [方法说明](docs/METHOD.md)：补全、绝缘分数、候选边界评分和对齐流程。
- `tests/test_impute.py`：邻居选择、权重和重复调用测试。
- [数据准备与示例](docs/DATA_PREPARATION.md)：接触表、细胞顺序和相似度矩阵准备。

`_source_snapshot/` 是本地原始候选文件快照，包含 SHA-256 校验记录，并被 `.gitignore` 排除；不要把这些历史脚本与正式发布代码混合上传。NAS 原文件没有修改。

## 本地运行

```powershell
.\.venv\Scripts\python.exe -m cellrob --help
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

完整输入格式与参数示例见英文 README。本地已经建立虚拟环境；其他机器请先按 README 安装。

补全时选择最多 impute_n 个不同邻居，排除自身，只纳入相似度严格大于 0.05 的细胞，使用各自相似度作为接触权重，并保持原相似度矩阵不变。

运行独立合成示例：`python examples/run_example.py --output results/example`。真实数据的相似度准备沿用 Higashi 细胞嵌入 → Scanpy 邻居图 → connectivities 权重矩阵的流程，脚本为 `scripts/prepare_similarity.py`。
