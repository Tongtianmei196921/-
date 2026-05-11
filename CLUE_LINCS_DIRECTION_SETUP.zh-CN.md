# DrugReflector 真实方向性证据接入说明

## 原则

DrugReflector 的模型排序不直接等于 Reverse 或 Mimic。Reverse / Mimic 只能来自真实 signed connectivity 证据，例如 CLUE / CMap / LINCS L1000 Touchstone 查询结果。

如果没有真实 signed connectivity 表，系统会继续显示 `No objective evidence`，不会猜测方向。

## 支持的真实证据表格式

后端支持通过环境变量读取 CSV / TSV：

```bash
DRUGREFLECTOR_CONNECTIVITY_TABLE=/path/to/clue_signed_connectivity.csv
```

也支持备用变量：

```bash
CLUE_CONNECTIVITY_TABLE=/path/to/clue_signed_connectivity.csv
LINCS_CONNECTIVITY_TABLE=/path/to/lincs_signed_connectivity.csv
```

表中至少需要两列：

| 必需内容 | 可接受列名示例 |
|---|---|
| 化合物编号 | `compound`, `pert_id`, `brd_id`, `broad_id`, `perturbagen` |
| signed connectivity 分数 | `connectivity_score`, `signed_connectivity_score`, `score`, `cs`, `tau` |

可选列：

| 可选内容 | 可接受列名示例 |
|---|---|
| signature/query 名称 | `signature`, `sample`, `query`, `query_name`, `input_signature` |
| 证据来源 | `source`, `dataset`, `evidence` |

示例：

```csv
signature,compound,connectivity_score,source
FXS_vs_Unaffected_all,BRD-K06854232,-87.3,CLUE Touchstone L1000
FXS_vs_Unaffected_all,BRD-K49685476,76.1,CLUE Touchstone L1000
```

## 判定规则

默认假设 CLUE/LINCS signed connectivity 分数范围为 -100 到 100：

| 分数 | 显示 |
|---|---|
| `score <= -20` | `Reverse` |
| `score >= 20` | `Mimic` |
| `-20 < score < 20` | `No objective evidence` |

如果表内分数范围像相关系数，即绝对值不超过 1.5，则默认阈值为 `0.1`。

可以手动设置阈值：

```bash
DRUGREFLECTOR_CONNECTIVITY_SCORE_THRESHOLD=30
```

## 当前 GSE198138 的 CLUE 查询文件

可以用脚本从 DrugReflector signature 生成 CLUE 查询所需的 up/down gene sets：

```bash
python scripts/export_clue_query_sets.py
```

输出目录：

```text
outputs/GSE198138_processed/clue_query/
```

其中：

| 文件 | 用途 |
|---|---|
| `clue_uptag_entrez.gmt` | 上传给 CLUE 的 up gene set，Entrez ID 格式 |
| `clue_dntag_entrez.gmt` | 上传给 CLUE 的 down gene set，Entrez ID 格式 |
| `clue_uptag_symbols.gmt` | 备查用 HGNC symbol 版本 |
| `clue_dntag_symbols.gmt` | 备查用 HGNC symbol 版本 |
| `clue_query_gene_set_mapping_stats.csv` | symbol 到 Entrez ID 的映射统计 |

CLUE 官方 Query API 通常需要 user key，并使用 up/down gene sets 查询 L1000 Touchstone。拿到真实结果后，把化合物编号和 signed connectivity 分数整理成上面的 CSV/TSV 格式，部署时设置 `DRUGREFLECTOR_CONNECTIVITY_TABLE` 即可。

