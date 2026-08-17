from pathlib import Path

from dataset import episode


def test_default_artifact_paths_prefer_processed_layout(monkeypatch, tmp_path):
    cache = tmp_path / "dataset" / "processed" / "cache"
    sealed = tmp_path / "dataset" / "processed" / "sealed" / "DAVIS"
    cache.mkdir(parents=True)
    sealed.mkdir(parents=True)
    expected_cache = cache / "protein_DAVIS.pt"
    expected_cache.touch()
    monkeypatch.setattr(episode, "ROOT", str(tmp_path))

    assert episode.default_sealed_dir("DAVIS") == str(sealed)
    assert episode.default_protein_cache_path("DAVIS") == str(expected_cache)


def test_default_artifact_paths_fall_back_to_legacy_layout(monkeypatch, tmp_path):
    cache = tmp_path / "dataset" / "cache"
    sealed = tmp_path / "dataset" / "sealed" / "DAVIS"
    cache.mkdir(parents=True)
    sealed.mkdir(parents=True)
    expected_cache = cache / "protein_DAVIS.pt"
    expected_cache.touch()
    monkeypatch.setattr(episode, "ROOT", str(tmp_path))

    assert episode.default_sealed_dir("DAVIS") == str(sealed)
    assert episode.default_protein_cache_path("DAVIS") == str(expected_cache)
