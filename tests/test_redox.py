"""氧化还原自动补全接口测试。"""
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def redox(reactants, products, medium):
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": reactants, "products": products, "medium": medium},
    )
    assert resp.status_code == 200, resp.text
    return resp.json()


# ---------- 酸性介质 ----------

def test_acidic_permanganate_iron():
    data = redox(["MnO4-", "Fe^2+"], ["Mn^2+", "Fe^3+"], "acidic")
    assert data["medium"] == "acidic"
    assert data["equation"] == "MnO4- + 5Fe^2+ + 8H+ -> Mn^2+ + 5Fe^3+ + 4H2O"
    assert data["reactants"] == ["MnO4-", "Fe^2+", "H+"]
    assert data["products"] == ["Mn^2+", "Fe^3+", "H2O"]
    assert data["coefficients"] == [1, 5, 8, 1, 5, 4]
    assert data["reactant_coefficients"] == [1, 5, 8]
    assert data["product_coefficients"] == [1, 5, 4]
    assert data["added_species"] == [
        {"formula": "H2O", "side": "products", "coefficient": 4},
        {"formula": "H+", "side": "reactants", "coefficient": 8},
    ]
    assert data["balanced"] is True
    assert all(e["balanced"] for e in data["elements"])
    assert data["charge"] == {"reactants": 17, "products": 17, "balanced": True}


def test_acidic_water_on_reactant_side():
    data = redox(["MnO2", "Cl-"], ["Mn^2+", "Cl2"], "acidic")
    assert data["equation"] == "MnO2 + 2Cl- + 4H+ -> Mn^2+ + Cl2 + 2H2O"
    assert data["coefficients"] == [1, 2, 4, 1, 1, 2]
    assert data["added_species"] == [
        {"formula": "H2O", "side": "products", "coefficient": 2},
        {"formula": "H+", "side": "reactants", "coefficient": 4},
    ]


def test_acidic_no_auxiliary_needed():
    data = redox(["H2", "O2"], ["H2O"], "acidic")
    assert data["equation"] == "2H2 + O2 -> 2H2O"
    assert data["coefficients"] == [2, 1, 2]
    assert data["added_species"] == []


def test_acidic_known_water_not_duplicated():
    # H2O 已是已知物种：不再重复成列，只补 H+
    data = redox(["MnO4-", "Fe^2+"], ["Mn^2+", "Fe^3+", "H2O"], "acidic")
    assert data["equation"] == "MnO4- + 5Fe^2+ + 8H+ -> Mn^2+ + 5Fe^3+ + 4H2O"
    assert data["coefficients"] == [1, 5, 8, 1, 5, 4]
    assert data["added_species"] == [
        {"formula": "H+", "side": "reactants", "coefficient": 8}
    ]


# ---------- 碱性介质 ----------

def test_basic_permanganate_sulfite():
    data = redox(["MnO4-", "SO3^2-"], ["MnO2", "SO4^2-"], "basic")
    assert data["medium"] == "basic"
    assert data["equation"] == "2MnO4- + 3SO3^2- + H2O -> 2MnO2 + 3SO4^2- + 2OH-"
    assert data["coefficients"] == [2, 3, 1, 2, 3, 2]
    assert data["added_species"] == [
        {"formula": "H2O", "side": "reactants", "coefficient": 1},
        {"formula": "OH-", "side": "products", "coefficient": 2},
    ]
    assert data["balanced"] is True
    assert data["charge"]["balanced"] is True
    # 净反应：同种辅助物种不得同时出现在两侧
    assert not (set(data["reactants"]) & set(data["products"]))


def test_basic_disproportionation_with_known_aux_species():
    # OH- 与 H2O 均已在已知物种中：辅助列全部去重，等价于直接配平
    data = redox(["Cl2", "OH-"], ["Cl-", "ClO3-", "H2O"], "basic")
    assert data["equation"] == "3Cl2 + 6OH- -> 5Cl- + ClO3- + 3H2O"
    assert data["coefficients"] == [3, 6, 5, 1, 3]
    assert data["added_species"] == []


def test_basic_known_water_on_reactant_side():
    # 碱性净反应中 H2O 本就在反应物侧：已知物种所在侧保持不变
    data = redox(["MnO4-", "Fe^2+", "H2O"], ["Mn^2+", "Fe^3+"], "basic")
    assert data["equation"] == "MnO4- + 5Fe^2+ + 4H2O -> Mn^2+ + 5Fe^3+ + 8OH-"
    assert data["coefficients"] == [1, 5, 4, 1, 5, 8]
    assert data["added_species"] == [
        {"formula": "OH-", "side": "products", "coefficient": 8}
    ]


def test_basic_half_reaction_with_electron():
    data = redox(["MnO4-", "e-"], ["MnO2"], "basic")
    assert data["equation"] == "MnO4- + 3e- + 2H2O -> MnO2 + 4OH-"
    assert data["coefficients"] == [1, 3, 2, 1, 4]
    assert data["added_species"] == [
        {"formula": "H2O", "side": "reactants", "coefficient": 2},
        {"formula": "OH-", "side": "products", "coefficient": 4},
    ]


# ---------- 错误映射 ----------

def test_multiple_solutions_returns_422():
    # H2/O2 -> H2O2 在自由引入 H2O/H+ 后存在多组最简解
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O2"], "medium": "acidic"},
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["null_space_dim"] == 2
    assert "多组最简解" in detail["message"]


def test_no_solution_returns_422():
    # KMnO4 -> KCl + MnCl2 + Cl2：氧原子无任何去处，无法补全
    resp = client.post(
        "/api/v1/redox/balance",
        json={
            "reactants": ["KMnO4"],
            "products": ["KCl", "MnCl2", "Cl2"],
            "medium": "acidic",
        },
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["null_space_dim"] == 0
    assert "无法补全" in detail["message"]


def test_side_conflict_returns_422():
    # O2 必须位于反应物侧才能守恒，但其所在侧被固定为生成物
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["H2"], "products": ["O2"], "medium": "acidic"},
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["null_space_dim"] == 1
    assert "所在侧" in detail["message"]


def test_basic_known_water_wrong_side_returns_422():
    # 碱性净反应要求 H2O 在反应物侧，但已知物种把它固定在生成物侧
    resp = client.post(
        "/api/v1/redox/balance",
        json={
            "reactants": ["MnO4-", "Fe^2+"],
            "products": ["Mn^2+", "Fe^3+", "H2O"],
            "medium": "basic",
        },
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 1


def test_bad_formula_returns_400():
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["H2O$"], "products": ["H2"], "medium": "acidic"},
    )
    assert resp.status_code == 400


def test_unknown_element_returns_400():
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["Xx2"], "products": ["Xx"], "medium": "basic"},
    )
    assert resp.status_code == 400


def test_invalid_medium_returns_422():
    resp = client.post(
        "/api/v1/redox/balance",
        json={"reactants": ["H2"], "products": ["H2O"], "medium": "neutral"},
    )
    assert resp.status_code == 422
