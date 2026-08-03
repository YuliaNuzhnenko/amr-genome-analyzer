import pytest
import pandas as pd
import os
from project import (
    find_resistance_genes,
    visualize_results,
    generate_demo_data
)


# ── тести generate_demo_data ──────────────────────

def test_generate_demo_data_returns_dataframe(tmp_path):
    """Перевіряє що генератор повертає DataFrame."""
    df = generate_demo_data(output_dir=str(tmp_path))
    assert isinstance(df, pd.DataFrame)


def test_generate_demo_data_has_correct_columns(tmp_path):
    """Перевіряє що DataFrame має потрібні колонки."""
    df = generate_demo_data(output_dir=str(tmp_path))
    assert "subject" in df.columns
    assert "identity" in df.columns


def test_generate_demo_data_creates_tsv_file(tmp_path):
    """Перевіряє що TSV файл створюється в папці results."""
    generate_demo_data(output_dir=str(tmp_path))
    blast_file = tmp_path / "blast_hits.tsv"
    assert blast_file.exists()


# ── тести find_resistance_genes ───────────────────

def test_find_resistance_genes_filters_correctly(tmp_path):
    """Перевіряє що фільтрація по identity працює правильно."""
    test_file = tmp_path / "test_blast.tsv"
    test_file.write_text(
        "NC_001\tARO:mecA\t95.0\t300\t5\t0\t1\t300\t1\t300\t0.0\t500\n"
        "NC_001\tARO:tetM\t60.0\t200\t10\t0\t1\t200\t1\t200\t1e-20\t300\n"
        "NC_001\tARO:blaZ\t80.0\t250\t7\t0\t1\t250\t1\t250\t1e-30\t400\n"
    )
    hits = find_resistance_genes(str(test_file), min_identity=70.0)
    assert len(hits) == 2


def test_find_resistance_genes_returns_dataframe(tmp_path):
    """Перевіряє що функція повертає DataFrame."""
    test_file = tmp_path / "test_blast.tsv"
    test_file.write_text(
        "NC_001\tARO:mecA\t95.0\t300\t5\t0\t1\t300\t1\t300\t0.0\t500\n"
    )
    hits = find_resistance_genes(str(test_file))
    assert isinstance(hits, pd.DataFrame)


def test_find_resistance_genes_sorted_by_identity(tmp_path):
    """Перевіряє що результати відсортовані від найкращого."""
    test_file = tmp_path / "test_blast.tsv"
    test_file.write_text(
        "NC_001\tARO:tetM\t75.0\t300\t5\t0\t1\t300\t1\t300\t1e-50\t500\n"
        "NC_001\tARO:mecA\t95.0\t200\t10\t0\t1\t200\t1\t200\t1e-20\t300\n"
    )
    hits = find_resistance_genes(str(test_file), min_identity=70.0)
    assert hits.iloc[0]["identity"] == 95.0


def test_find_resistance_genes_empty_result(tmp_path):
    """Перевіряє що повертається порожній DataFrame якщо нічого не знайдено."""
    test_file = tmp_path / "test_blast.tsv"
    test_file.write_text(
        "NC_001\tARO:mecA\t50.0\t300\t5\t0\t1\t300\t1\t300\t1e-50\t500\n"
    )
    hits = find_resistance_genes(str(test_file), min_identity=70.0)
    assert hits.empty


def test_find_resistance_genes_raises_on_missing_file():
    """Перевіряє що функція кидає FileNotFoundError якщо файл не існує."""
    import pytest
    with pytest.raises(FileNotFoundError):
        find_resistance_genes("nonexistent_file.tsv")


def test_find_resistance_genes_raises_on_empty_file(tmp_path):
    """Перевіряє що порожній файл викидає ValueError."""
    import pytest
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("")
    with pytest.raises(ValueError):
        find_resistance_genes(str(empty_file))


# ── тести visualize_results ───────────────────────

def test_visualize_results_creates_file(tmp_path):
    """Перевіряє що функція створює PNG файл."""
    test_data = pd.DataFrame({
        "subject": ["mecA|beta-lactam", "tetM|tetracycline", "blaZ|penicillin"],
        "identity": [99.2, 94.1, 87.8]
    })
    output_path = visualize_results(
        test_data,
        output_dir=str(tmp_path),
        top_n=3
    )
    assert output_path is not None
    assert os.path.exists(output_path)


def test_visualize_results_empty_dataframe(tmp_path):
    """Перевіряє що порожній DataFrame обробляється без помилок."""
    empty_df = pd.DataFrame()
    result = visualize_results(empty_df, output_dir=str(tmp_path))
    assert result is None
