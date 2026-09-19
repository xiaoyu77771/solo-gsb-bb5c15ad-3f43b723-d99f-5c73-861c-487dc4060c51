"""氧化还原自动补全接口测试。"""
import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def redox(reactants, products, medium="acidic", expect=200):
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": reactants, "products": products, "medium": medium},
    )
    assert resp.status_code == expect, resp.text
    return resp.json()


# ---------- 酸性介质 ----------

def test_redox_acidic_permanganate_iron():
    data = redox(["MnO4-", "Fe^2+"], ["Mn^2+", "Fe^3+"], "acidic")
    assert data["equation"] == "MnO4- + 5Fe^2+ + 8H+ -> Mn^2+ + 5Fe^3+ + 4H2O"
    # 已知物种保持所在侧且顺序不变，辅助物种追加在对应侧末尾
    assert data["reactants"] == ["MnO4-", "Fe^2+", "H+"]
    assert data["products"] == ["Mn^2+", "Fe^3+", "H2O"]
    assert data["coefficients"] == [1, 5, 8, 1, 5, 4]
    assert data["reactant_coefficients"] == [1, 5, 8]
    assert data["product_coefficients"] == [1, 5, 4]
    assert data["added_species"] == [
        {"formula": "H+", "side": "reactant", "coefficient": 8},
        {"formula": "H2O", "side": "product", "coefficient": 4},
    ]
    assert data["medium"] == "acidic"
    assert data["balanced"] is True
    assert all(e["balanced"] for e in data["elements"])
    assert data["charge"] == {"reactants": 17, "products": 17, "balanced": True}


def test_redox_acidic_water_on_reactant_side():
    # Cl2 歧化：H2O 落在反应物侧、H+ 落在生成物侧
    data = redox(["Cl2"], ["Cl-", "ClO-"], "acidic")
    assert data["equation"] == "Cl2 + H2O -> Cl- + ClO- + 2H+"
    assert data["coefficients"] == [1, 1, 1, 1, 2]
    assert data["added_species"] == [
        {"formula": "H2O", "side": "reactant", "coefficient": 1},
        {"formula": "H+", "side": "product", "coefficient": 2},
    ]
    assert data["charge"]["balanced"] is True


def test_redox_acidic_half_reaction_with_electron():
    # 电子作为已知物种保留在原侧，仅补 H+ 与 H2O
    data = redox(["MnO4-", "e-"], ["MnO2"], "acidic")
    assert data["equation"] == "MnO4- + 3e- + 4H+ -> MnO2 + 2H2O"
    assert data["coefficients"] == [1, 3, 4, 1, 2]
    assert data["added_species"] == [
        {"formula": "H+", "side": "reactant", "coefficient": 4},
        {"formula": "H2O", "side": "product", "coefficient": 2},
    ]


def test_redox_default_medium_is_acidic():
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["MnO4-", "Fe^2+"], "products": ["Mn^2+", "Fe^3+"]},
    )
    assert resp.status_code == 200
    assert resp.json()["medium"] == "acidic"
    assert "H+" in resp.json()["reactants"]


# ---------- 碱性介质 ----------

def test_redox_basic_permanganate_sulfite():
    data = redox(["MnO4-", "SO3^2-"], ["MnO2", "SO4^2-"], "basic")
    assert data["equation"] == "2MnO4- + 3SO3^2- + H2O -> 2MnO2 + 3SO4^2- + 2OH-"
    assert data["coefficients"] == [2, 3, 1, 2, 3, 2]
    assert data["added_species"] == [
        {"formula": "H2O", "side": "reactant", "coefficient": 1},
        {"formula": "OH-", "side": "product", "coefficient": 2},
    ]
    assert data["medium"] == "basic"
    assert data["charge"] == {"reactants": -8, "products": -8, "balanced": True}
    assert all(e["balanced"] for e in data["elements"])


def test_redox_basic_disproportionation():
    data = redox(["Cl2"], ["Cl-", "ClO3-"], "basic")
    assert data["equation"] == "3Cl2 + 6OH- -> 5Cl- + ClO3- + 3H2O"
    assert data["coefficients"] == [3, 6, 5, 1, 3]
    assert data["added_species"] == [
        {"formula": "OH-", "side": "reactant", "coefficient": 6},
        {"formula": "H2O", "side": "product", "coefficient": 3},
    ]


def test_redox_basic_half_reaction_with_electron():
    data = redox(["MnO4-", "e-"], ["MnO2"], "basic")
    assert data["equation"] == "MnO4- + 3e- + 2H2O -> MnO2 + 4OH-"
    assert data["coefficients"] == [1, 3, 2, 1, 4]


def test_redox_basic_never_introduces_proton():
    # 碱性模式下新增辅助物种只能是 H2O / OH-
    data = redox(["MnO4-", "Fe^2+"], ["Mn^2+", "Fe^3+"], "basic")
    formulas = [s["formula"] for s in data["added_species"]]
    assert "H+" not in formulas
    assert set(formulas) <= {"H2O", "OH-"}
    # 同一辅助物种不会同时出现在两侧（净反应）
    assert not (set(data["reactants"]) & set(data["products"]))


def test_redox_acidic_and_basic_share_known_coefficients():
    # 同一组已知物种：两种介质下已知物种系数一致，仅辅助物种不同
    acidic = redox(["MnO4-", "SO3^2-"], ["MnO2", "SO4^2-"], "acidic")
    basic = redox(["MnO4-", "SO3^2-"], ["MnO2", "SO4^2-"], "basic")
    assert acidic["coefficients"][:2] == basic["coefficients"][:2] == [2, 3]
    assert acidic["equation"] == "2MnO4- + 3SO3^2- + 2H+ -> 2MnO2 + 3SO4^2- + H2O"
    assert basic["equation"] == "2MnO4- + 3SO3^2- + H2O -> 2MnO2 + 3SO4^2- + 2OH-"


# ---------- 无需补全 / 边界 ----------

def test_redox_no_auxiliary_needed():
    data = redox(["H2", "Cl2"], ["HCl"], "acidic")
    assert data["equation"] == "H2 + Cl2 -> 2HCl"
    assert data["coefficients"] == [1, 1, 2]
    assert data["added_species"] == []


def test_redox_impossible_returns_422():
    # Fe^2+ -> Fe^3+：电荷无法仅靠 H2O/H+/OH- 补平（零空间维数 0）
    data = redox(["Fe^2+"], ["Fe^3+"], "acidic", expect=422)
    assert data["detail"]["null_space_dim"] == 0
    assert "无法补全" in data["detail"]["message"]


def test_redox_known_species_wrong_side_returns_422():
    # H2 -> H2O2：零空间一维但 H2 系数为负，无法改变已知物种所在侧
    data = redox(["H2"], ["H2O2"], "acidic", expect=422)
    assert data["detail"]["null_space_dim"] == 1
    assert "全正" in data["detail"]["message"]


def test_redox_multiple_solutions_returns_422():
    # NO 与 NO2 同为生成物：存在多组最简解，无法唯一判定
    data = redox(["Cu", "HNO3"], ["Cu(NO3)2", "NO", "NO2"], "acidic", expect=422)
    assert data["detail"]["null_space_dim"] == 2
    assert "多组最简解" in data["detail"]["message"]


@pytest.mark.parametrize("aux", ["H2O", "H+", "H^+", "OH-", "HO-"])
def test_redox_known_species_conflicting_aux_returns_400(aux):
    # 已知物种与自动补全物种（H2O/H+/OH-，含等价写法）重复 → 400
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": [aux, "O2"], "products": ["H2O2"], "medium": "acidic"},
    )
    assert resp.status_code == 400
    assert "自动补全物种" in resp.json()["detail"]


def test_redox_bad_formula_returns_400():
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["H2$"], "products": ["H2O2"], "medium": "acidic"},
    )
    assert resp.status_code == 400


def test_redox_invalid_medium_returns_422():
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["H2"], "products": ["H2O2"], "medium": "neutral"},
    )
    assert resp.status_code == 422
