import pytest

import pystemfinder as psf


def test_mouse_cell_cycle_phases():
    s = psf.cell_cycle_genes("mouse", phase="S")
    g2m = psf.cell_cycle_genes("mouse", phase="G2M")

    assert (len(s), len(g2m)) == (40, 51)
    assert s[:4] == ["Rrm2", "Mcm4", "Msh2", "Exo1"]
    assert g2m[:4] == ["Cbx5", "Gtse1", "Dlgap5", "Smc4"]


@pytest.mark.parametrize("species, gene", [("mouse", "E2f8"), ("human", "E2F8")])
def test_e2f8_is_an_s_phase_gene_only(species, gene):
    # the R package also lists it under G2M; Seurat's cc.genes lists it under S only
    assert gene in psf.cell_cycle_genes(species, phase="S")
    assert gene not in psf.cell_cycle_genes(species, phase="G2M")


@pytest.mark.parametrize("species, n_genes", [("mouse", 91), ("human", 97), ("celegans", 130)])
def test_both_phases_list_each_gene_once(species, n_genes):
    # C. elegans: 11 genes (cyclins and others) are on both of the R package's lists
    s = psf.cell_cycle_genes(species, phase="S")
    both = psf.cell_cycle_genes(species)

    assert len(both) == len(set(both)) == n_genes
    assert both[: len(s)] == s


@pytest.mark.parametrize(
    "species, n_s, n_g2m, first_s",
    [("human", 43, 54, "MCM5"), ("celegans", 41, 100, "orc-1")],
)
def test_other_species(species, n_s, n_g2m, first_s):
    s = psf.cell_cycle_genes(species, phase="S")

    assert (len(s), len(psf.cell_cycle_genes(species, phase="G2M"))) == (n_s, n_g2m)
    assert s[0] == first_s


@pytest.mark.parametrize("species, n, first", [("human", 1652, "ZFY"), ("mouse", 1838, "Adnp"), ("celegans", 305, "ppw-1")])
def test_transcription_factors_are_unique(species, n, first):
    tfs = psf.transcription_factors(species)

    assert len(tfs) == len(set(tfs)) == n
    assert tfs[0] == first


def test_unknown_species_raises():
    with pytest.raises(ValueError, match="species"):
        psf.cell_cycle_genes("zebrafish")
    with pytest.raises(ValueError, match="species"):
        psf.transcription_factors("zebrafish")


def test_unknown_phase_raises():
    with pytest.raises(ValueError, match="phase"):
        psf.cell_cycle_genes("mouse", phase="M")
