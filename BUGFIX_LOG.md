# Bug Fix Log

## 2026-05-09

### Bug 1 — `core.py:49` — assert 消息缺少 `f` 前缀

**影响**：触发断言时打印字面量 `{save_path}`，不会插值实际路径，难以定位问题目录。

```python
# 修复前
assert os.path.isdir(save_path), "{save_path} doesn't exist or not a dir"

# 修复后
assert os.path.isdir(save_path), f"{save_path} doesn't exist or not a dir"
```

---

### Bug 2 — `dist_complex.py:14` — `NameError: mean_pmf` 未定义

**影响**：`combine_complex_distribution_df` 中，当复合体的所有蛋白均不在分布 DataFrame 的列中时，走到空分支会抛出 `NameError: name 'mean_pmf' is not defined`，导致程序崩溃。

```python
# 修复前
return pd.Series(index=mean_pmf.index)

# 修复后
return pd.Series(index=set_distribution_df.index)
```

---

### Bug 3 — `infer_query.py:130` — 不可达的死代码

**影响**：无运行时错误，但代码逻辑错误，`return "The variable is a string."` 在所有 try/except 分支均已 return 的情况下永远不会执行。

```python
# 修复前
except json.JSONDecodeError:
    ...
    return None
return "The variable is a string."   # 永远不可达

# 修复后（删除死代码）
except json.JSONDecodeError:
    ...
    return None
```

---

### Bug 4 — `infer_query.py:402, 404, 511, 513` — `np.max(x, 0)` 不能截断下界为 0

**影响**：`np.max(x, 0)` 的第二个参数被当作 `axis=0` 而非比较值，截断逻辑静默失效，`ligand_low` / `receptor_low` 可能为负数。同时 `receptor_low` 的乘法顺序有误（应先乘再截断）。

```python
# 修复前
ligand_low  = np.max(ref_p1.loc[index, col] * (1-1/k*1.96), 0)
receptor_low = np.max(ref_p2.loc[index, col], 0) * (1-1/k*1.96)

# 修复后
ligand_low  = np.maximum(ref_p1.loc[index, col] * (1-1/k*1.96), 0)
receptor_low = np.maximum(ref_p2.loc[index, col] * (1-1/k*1.96), 0)
```

---

### Bug 5 — `distrib.py:384` — 支撑域边界检查差一个等号

**影响**：`log1p` 值恰好等于 14.0（高表达基因可能触及）时，断言直接抛出 `AssertionError`，但 14.0 实为合法值（直方图 bin edge 包含该值）。

```python
# 修复前
assert np.max(samples) < max_value_4_log1p, "Support domain is not valid."

# 修复后
assert np.max(samples) <= max_value_4_log1p, "Support domain is not valid."
```

---

### Bug 6 — `distrib.py:433, 448` — 调试用 `print()` 残留

**影响**：`get_pvalue_from_complex_pmf` 中两行裸 `print()` 会在正常运行时向 stdout 打印调试信息。

```python
# 删除以下两行
print(pmf.min_cdf_non_zero, pmf.min_cdf_one)   # 第 433 行
print(precise_pmf_array[49:54])                   # 第 448 行
```

---

### Bug 7 — `distrib.py:537` — `get_quantile_pmf_for_n_iid_distribution` 中 numpy 负索引陷阱（核心 bug）

**影响**：PMF 和显著偏离 1（实测约为 2.0），导致 p 值计算完全错误。仅影响 Quantile 类方法（`Q3`、`Quantile_0.9` 等），Mean 方法不受影响。

**触发条件**：某个基因在当前 cell type 所有细胞中均有表达（即 `pmf[0] = 0`，零表达概率为 0）时：

1. `split_integer_by_probability` 对第一个 bin 分配整数 0
2. `section[0] = cumsum[0] = 0`
3. `section[0] - 1 = -1`
4. numpy 负索引：`fxi_cdf[-1] = 1.0`（最后一个元素，而非期望的 `0.0`）
5. `cdf[0]` 错误地等于 `1.0`，后续 diff 产生大量负值（如 `-0.999`）
6. `np.clip` 截掉负值后，sum 约为 2.0

**修复**：在 `fxi_cdf` 前插入 `0.0`，使 `section=0` 时正确返回 `0.0`：

```python
# 修复前
cdf = fxi_cdf[section-1]

# 修复后
fxi_cdf_padded = np.concatenate([[0.0], fxi_cdf])
cdf = fxi_cdf_padded[section]
```

索引边界验证：
- `section[k] = 0` → `fxi_cdf_padded[0] = 0.0` ✓
- `section[k] = j > 0` → `fxi_cdf_padded[j] = fxi_cdf[j-1]`（与原逻辑一致）✓
- `section[-1] = 100000` → `fxi_cdf_padded[100000] = fxi_cdf[99999] = 1.0` ✓

---

## 2026-05-21

### Bug 8 — `distrib.py:165` — `__add__` 将 `is_complex_analytic` 设在操作数而非结果上（核心 bug）

**影响**：`get_pvalue_from_complex_pmf` 高精度解析路径永远不可达，静默退化为低精度的离散 PMF 路径（`1 - np.sum(y[...])`）。仅影响 L-R 其中一方为 min-of-normals 类型的 complex protein 的情形。

**根本原因**：`__add__` 中，当两个操作数均为 analytic 分布时，`is_complex_analytic = True` 被设在了 `self`（左操作数）上，而 `__truediv__` 后续检查的是返回值 `new_pmf.is_complex_analytic`，该值始终为 `False`，导致传播逻辑从不触发：

```python
# 修复前
if self.is_analytic and pmf2.is_analytic:
    self.is_complex_analytic = True   # 设在输入对象上，结果对象不受影响
    new_pmf.ligand = self
    new_pmf.receptor = pmf2

# 修复后
if self.is_analytic and pmf2.is_analytic:
    new_pmf.is_complex_analytic = True  # 设在返回的结果对象上
    new_pmf.ligand = self
    new_pmf.receptor = pmf2
```

---

### Bug 9 — `distrib.py:401` — `get_pvalue_from_pmf` 在 analytic 路径下多余地调用 `get_pmf_array()`

**影响**：`is_analytic` 或 `is_complex_analytic` 为 `True` 时，`y = pmf.get_pmf_array()` 会触发不必要的 PMF 离散化计算（对正态分布为 CDF 数值积分），但结果 `y` 在这两个分支中完全不使用。在大规模 L-R 对批量计算时造成可感知的性能损失。

```python
# 修复前
def get_pvalue_from_pmf(value, pmf):
    y = pmf.get_pmf_array()   # 无条件调用，analytic 分支不需要
    if pmf.is_analytic:
        pvalue = 1 - pmf.cdf_analytic_func(value)
    elif pmf.is_complex_analytic:
        pvalue = get_pvalue_from_complex_pmf(value, pmf)
    else:
        pvalue = 1 - np.sum(y[:int(np.ceil(value / precision))])

# 修复后
def get_pvalue_from_pmf(value, pmf):
    if pmf.is_analytic:
        pvalue = 1 - pmf.cdf_analytic_func(value)
    elif pmf.is_complex_analytic:
        pvalue = get_pvalue_from_complex_pmf(value, pmf)
    else:
        y = pmf.get_pmf_array()   # 仅在需要时调用
        pvalue = 1 - np.sum(y[:int(np.ceil(value / precision))])
```

---

### Bug 10 — `dist_lr.py` — `timeit` 调试残留

**影响**：`start`/`stop` 计时变量被计算但从不使用，`import timeit` 也成为无用导入。

```python
# 删除以下内容
import timeit
...
start = timeit.default_timer()
...
stop = timeit.default_timer()
```
