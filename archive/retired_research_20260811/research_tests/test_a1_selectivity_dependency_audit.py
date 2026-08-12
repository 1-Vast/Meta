from research.meta_fewshot.a1_selectivity_dependency_audit import UnionFind, union_shared


def test_dependency_closure_is_transitive():
    components = UnionFind(["a", "b", "c"])
    union_shared(components, {"x": {"a", "b"}, "y": {"b", "c"}})
    assert len({components.find(key) for key in ("a", "b", "c")}) == 1


def test_dependency_closure_keeps_unrelated_groups_apart():
    components = UnionFind(["a", "b", "c"])
    union_shared(components, {"x": {"a", "b"}})
    assert components.find("a") == components.find("b")
    assert components.find("a") != components.find("c")
