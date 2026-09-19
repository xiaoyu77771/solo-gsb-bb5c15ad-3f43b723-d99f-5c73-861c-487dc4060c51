"""氧化还原反应自动补全：仅引入 H2O / H+ / OH- 补齐守恒项。

已知物种所在侧固定；按介质把 H2O、H+（酸性）或 H2O、OH-（碱性）
作为可出现在任一侧的辅助物种加入守恒矩阵（含电荷行），用有理数
零空间求出唯一（不计倍数）的最小正整数净反应：

  - 零空间维数 > 1：存在多组最简解，无法唯一判定（422）；
  - 零空间维数 = 0：守恒方程组无非零解，无法补全（422）；
  - 维数 = 1 但已知物种系数异号或为零：无法在不改变已知物种
    所在侧的前提下补全（422）。

碱性介质直接用 OH- 辅助列求解，与教科书“先按酸性配平、再把 H+
按 H+ + OH- -> H2O 折算为 H2O / OH-、约消两侧同种辅助物”的流程
等价（列向量 OH- = H2O - H+，张成同一零空间，净反应相同）。
辅助物种以带符号系数表示（正为反应物侧、负为生成物侧），同种
辅助物天然只出现在一侧，无需额外约消。与已知物种同式的辅助
物种不再重复成列，避免零空间维数虚增。
"""
from __future__ import annotations

from collections import Counter

from .balancer import (
    BalanceError,
    account,
    build_matrix,
    format_equation,
    null_space,
    to_integer_coeffs,
)
from .parser import parse_formula

# 各介质下允许引入的辅助物种（即响应中 added_species 的列举顺序）
MEDIUM_AUX = {
    "acidic": ("H2O", "H+"),
    "basic": ("H2O", "OH-"),
}


def _signature(counts: Counter, charge: int) -> tuple:
    """物种签名（元素组成 + 电荷）：同式不同写法的物种签名一致。"""
    return tuple(sorted(counts.items())), charge


def balance_redox(reactants: list[str], products: list[str], medium: str) -> dict:
    """主流程：辅助列建矩阵 → 零空间 → 定向整数化 → 组装净反应。

    FormulaError / BalanceError 继续向上抛，由路由层映射为 4xx。
    """
    n_r = len(reactants)
    known = list(reactants) + list(products)
    known_parsed = [parse_formula(f) for f in known]
    known_signatures = {_signature(c, q) for c, q in known_parsed}

    # 辅助物种列：跳过与已知物种同式者，避免重复列虚增零空间维数
    aux: list[tuple[str, Counter, int]] = []
    for formula in MEDIUM_AUX[medium]:
        counts, charge = parse_formula(formula)
        if _signature(counts, charge) not in known_signatures:
            aux.append((formula, counts, charge))

    compositions = [c for c, _ in known_parsed] + [a[1] for a in aux]
    charges = [q for _, q in known_parsed] + [a[2] for a in aux]

    # 全部列按“反应物为正”约定建矩阵：解向量中已知反应物须为正、
    # 已知生成物须为负；辅助列正则位于反应物侧、负则位于生成物侧
    _, matrix = build_matrix(compositions, charges, len(compositions))
    basis = null_space(matrix)
    if len(basis) != 1:
        if not basis:
            raise BalanceError(
                "守恒方程组无非零解：仅引入 H2O / H+ / OH- 无法补全该反应",
                nullity=0,
            )
        raise BalanceError(
            f"零空间维数为 {len(basis)}（必须恰好为 1 才能唯一判定），"
            "存在多组最简解，无法唯一补全",
            nullity=len(basis),
        )
    vec = basis[0]

    # 定向：折算到各自所在侧后，已知物种系数必须全部为正
    sided = [vec[j] if j < n_r else -vec[j] for j in range(len(known))]
    if all(s < 0 for s in sided):
        vec = [-x for x in vec]
        sided = [-s for s in sided]
    if not all(s > 0 for s in sided):
        zeros = [f for f, s in zip(known, sided) if s == 0]
        reason = (
            f"物种系数为零（不参与净反应）: {', '.join(zeros)}"
            if zeros
            else "已知物种系数异号（部分物种需位于对侧才能守恒）"
        )
        raise BalanceError(
            f"无法在不改变已知物种所在侧的前提下补全：{reason}",
            nullity=1,
        )

    ints = to_integer_coeffs(vec)
    n_known = len(known)
    reactant_coeffs = ints[:n_r]
    product_coeffs = [-x for x in ints[n_r:n_known]]
    aux_amounts = ints[n_known:]  # 正：反应物侧；负：生成物侧

    final_reactants = list(reactants)
    final_products = list(products)
    added = []
    for (formula, _, _), amount in zip(aux, aux_amounts):
        if amount == 0:
            continue
        if amount > 0:
            final_reactants.append(formula)
            reactant_coeffs.append(amount)
        else:
            final_products.append(formula)
            product_coeffs.append(-amount)
        added.append(
            {
                "formula": formula,
                "side": "reactants" if amount > 0 else "products",
                "coefficient": abs(amount),
            }
        )

    final_formulas = final_reactants + final_products
    final_coeffs = reactant_coeffs + product_coeffs
    parsed = [parse_formula(f) for f in final_formulas]
    final_comps = [c for c, _ in parsed]
    final_charges = [q for _, q in parsed]
    elements, _ = build_matrix(final_comps, final_charges, len(final_reactants))
    per_element, charge = account(
        final_comps, final_charges, final_coeffs, len(final_reactants), elements
    )
    return {
        "medium": medium,
        "equation": format_equation(final_reactants, final_products, final_coeffs),
        "reactants": final_reactants,
        "products": final_products,
        "coefficients": final_coeffs,
        "reactant_coefficients": reactant_coeffs,
        "product_coefficients": product_coeffs,
        "added_species": added,
        "elements": per_element,
        "charge": charge,
        "balanced": charge["balanced"] and all(e["balanced"] for e in per_element),
    }
