"""Gene lists bundled with the R stemFinder package (https://github.com/CahanLab/stemfinder, data/).

The files in pystemfinder/data are copies of the R package's lists, except that E2f8 (mouse) and E2F8 (human)
are removed from the G2M lists. The R package lists them under both S and G2M; they are S phase genes, as in
Seurat's cc.genes.
"""
from importlib.resources import files

# species -> (suffix of the R cell cycle lists, name of the R transcription factor list)
_SPECIES = {'human': ('human', 'hsTFs'), 'mouse': ('mouse', 'mmTFs'), 'celegans': ('celeg', 'ceTFs')}


def _read(name):
    return (files('pystemfinder') / 'data' / f'{name}.txt').read_text().split()


def _species(species):
    if species not in _SPECIES:
        raise ValueError(f"Unknown species '{species}'. Expected one of {', '.join(_SPECIES)}.")
    return _SPECIES[species]


def cell_cycle_genes(species='mouse', phase='both'):
    """S and G2M phase cell cycle genes, the standard markers for :func:`stemfinder`.

    With ``phase='both'``, returns the S genes followed by the G2M genes that are not also S genes, so each gene is
    listed once. (The R vignette's ``c(s_genes_mouse, g2m_genes_mouse)`` lists E2f8 twice; the R package's
    C. elegans lists share 11 genes.)

    Args:
        species (str, optional): ``'human'``, ``'mouse'``, or ``'celegans'``. Defaults to ``'mouse'``.
        phase (str, optional): ``'S'``, ``'G2M'``, or ``'both'``. Defaults to ``'both'``.

    Returns:
        list of str: Gene symbols.
    """
    suffix, _ = _species(species)
    if phase not in ['S', 'G2M', 'both']:
        raise ValueError("Invalid value for 'phase'. Expected one of 'S', 'G2M', 'both'.")
    s = _read(f's_genes_{suffix}') if phase in ['S', 'both'] else []
    g2m = _read(f'g2m_genes_{suffix}') if phase in ['G2M', 'both'] else []
    return list(dict.fromkeys(s + g2m))


def transcription_factors(species='mouse'):
    """Reference transcription factor list (``hsTFs``, ``mmTFs``, or ``ceTFs`` in the R package), without repeats.

    Args:
        species (str, optional): ``'human'``, ``'mouse'``, or ``'celegans'``. Defaults to ``'mouse'``.

    Returns:
        list of str: Gene symbols.
    """
    _, name = _species(species)
    return list(dict.fromkeys(_read(name)))
