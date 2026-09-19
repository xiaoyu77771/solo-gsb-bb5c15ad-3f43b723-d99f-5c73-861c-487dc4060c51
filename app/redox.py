"""氧化还原反应自动补全配平。

在 /api/v1/balance 的精确配平基础上放宽输入：已知物种两侧的元素与
电荷不必守恒，服务仅引入 H2O、H+（酸性）或 H2O、OH-（碱性）作为
辅助物种补齐守恒项，且已知物种所在侧保持不变。

实现要点：
  - 辅助物种统一按生成物列（负号）追加到守恒矩阵后求零空间；
    辅助列系数可正可负：正号表示落在生成物侧，负号表示应移到
    反应物侧（取绝对值），为零表示无需引入。
  - 零空间维数为 0 → 无法补全；大于 1 → 存在多组最简解，无法唯一
    判定；恰为 1 时要求已知物种系数同号（取全正方向），否则同样
    无法在不改变已知物种所在侧的前提下补全。三种情况均抛
    BalanceError（对应 HTTP 422，与 /api/v1/balance 一致）。
  - 碱性模式：先求酸性形式（H2O / H+）的解，再按 H+ + OH- -> H2O
    把 H+ 转换为 H2O / OH-——即辅助系数作代换 w' = w + h、oh = -h
    （w、h 分别为 H2O、H+ 的带符号系数），两侧同种辅助物随之约消。
    该代换与直接以 H2O / OH- 为辅助列求解完全等价（OH- 列恰为
    H2O 列减 H+ 列，两者张成同一 (H, O, 电荷) 调整空间），且
    gcd(knowns, w, h) = gcd(knowns, w + h, h)，代换后系数仍为
    最大公约数约简后的最小正整数。
  - 已知物种不得与 H2O / H+ / OH- 解析等价（无论介质）：否则守恒
    矩阵出现重复列，解被人为放大为不唯一。此类输入抛
    RedoxInputError（对应 HTTP 400），提示改用 /api/v1/balance。
"""
from __future__ import annotations

from .balancer import (
    BalanceError,
    account,
    build_matrix,
    format_equation,
    null_space,
    to_integer_coeffs,
)
from .parser import parse_formula

# 酸性核心统一以 H2O / H+ 为辅助列求解；碱性模式再把 H+ 转换为 H2O / OH-
_ACIDIC_AUX = ("H2O", "H+")
# 已知物种禁止与任一辅助物种等价（否则出现重复列，解不唯一）
_FORBIDDEN_KNOWN = ("H2O", "H+", "OH-")


class RedoxInputError(ValueError):
    """请求本身不合法（已知物种与自动补全物种冲突），对应 HTTP 400。"""


def _species_key(formula: str) -> tuple:
    """化学式的规范键（排序后的元素组成 + 电荷），用于识别辅助物种。"""
    composition, charge = parse_formula(formula)
    return tuple(sorted(composition.items())), charge


def balance_redox(
    reactants: list[str], products: list[str], medium: str
) -> dict:
    """自动补全配平主流程。

    FormulaError / BalanceError / RedoxInputError 继续向上抛，
    由路由层映射为 4xx。
    """
    n_r = len(reactants)
    known_formulas = list(reactants) + list(products)
    parsed = [parse_formula(f) for f in known_formulas]

    # ---- 已知物种不得与辅助物种（H2O / H+ / OH-）等价 ----
    forbidden = {_species_key(f): f for f in _FORBIDDEN_KNOWN}
    for formula, (composition, charge) in zip(known_formulas, parsed):
        key = (tuple(sorted(composition.items())), charge)
        if key in forbidden:
            raise RedoxInputError(
                f"已知物种 {formula!r} 与自动补全物种 {forbidden[key]!r} 重复："
                "H2O、H+、OH- 由服务按需自动引入，请从已知物种中移除；"
                "如需固定其所在侧，请改用 /api/v1/balance"
            )

    # ---- 守恒矩阵：已知物种 + 辅助列（形式生成物列，可正可负） ----
    aux_parsed = [parse_formula(f) for f in _ACIDIC_AUX]
    all_compositions = [c for c, _ in parsed] + [c for c, _ in aux_parsed]
    all_charges = [q for _, q in parsed] + [q for _, q in aux_parsed]
    _, matrix = build_matrix(all_compositions, all_charges, n_r)
    basis = null_space(matrix)

    n_known = len(known_formulas)
    if not basis:
        raise BalanceError(
            "无法补全：仅引入 H2O、H+、OH- 后守恒方程组仍无非零解，"
            "请检查已知物种是否完整（如缺少电子转移的另一对氧化还原电对）",
            nullity=0,
        )
    if len(basis) > 1:
        raise BalanceError(
            f"存在多组最简解，无法唯一判定（零空间维数为 {len(basis)}）："
            "请补充已知物种以约束反应",
            nullity=len(basis),
        )

    # ---- 定向：已知物种系数必须能化为全正（所在侧保持不变） ----
    vec = basis[0]
    if all(x < 0 for x in vec[:n_known]):
        vec = [-x for x in vec]
    if not all(x > 0 for x in vec[:n_known]):
        raise BalanceError(
            "无法补全：零空间基向量中已知物种系数无法化为全正"
            "（存在异号或零分量），在不改变已知物种所在侧的前提下，"
            "仅引入 H2O、H+、OH- 不能配平该反应",
            nullity=1,
        )

    coeffs = to_integer_coeffs(vec)
    known_coeffs = coeffs[:n_known]
    water, proton = coeffs[n_known], coeffs[n_known + 1]

    # ---- 介质：碱性模式把 H+ 转换为 H2O / OH- 并约消同种辅助物 ----
    if medium == "basic":
        aux = (("H2O", water + proton), ("OH-", -proton))
    else:
        aux = (("H2O", water), ("H+", proton))

    # 辅助列按符号归侧：正号留在生成物侧，负号移到反应物侧，零不引入
    left_aux = [(f, -c) for f, c in aux if c < 0]
    right_aux = [(f, c) for f, c in aux if c > 0]

    final_reactants = list(reactants) + [f for f, _ in left_aux]
    final_products = list(products) + [f for f, _ in right_aux]
    final_coeffs = (
        known_coeffs[:n_r]
        + [c for _, c in left_aux]
        + known_coeffs[n_r:]
        + [c for _, c in right_aux]
    )

    # ---- 逐元素与电荷核算（复用 balancer.account） ----
    final_parsed = [parse_formula(f) for f in final_reactants + final_products]
    final_compositions = [c for c, _ in final_parsed]
    final_charges = [q for _, q in final_parsed]
    elements = sorted({el for comp in final_compositions for el in comp})
    per_element, charge = account(
        final_compositions,
        final_charges,
        final_coeffs,
        len(final_reactants),
        elements,
    )

    added_species = [
        {"formula": f, "side": "reactant", "coefficient": c} for f, c in left_aux
    ] + [
        {"formula": f, "side": "product", "coefficient": c} for f, c in right_aux
    ]

    return {
        "equation": format_equation(final_reactants, final_products, final_coeffs),
        "medium": medium,
        "reactants": final_reactants,
        "products": final_products,
        "coefficients": final_coeffs,
        "reactant_coefficients": final_coeffs[: len(final_reactants)],
        "product_coefficients": final_coeffs[len(final_reactants):],
        "added_species": added_species,
        "elements": per_element,
        "charge": charge,
        "balanced": True,
    }
