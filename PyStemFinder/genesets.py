"""Gene lists bundled with the R stemFinder package (https://github.com/CahanLab/stemfinder, data/)."""
from importlib.resources import files

__all__ = ['cell_cycle_genes', 'transcription_factors']

# species -> (suffix of the R cell cycle lists, name of the R transcription factor list)
_SPECIES = {'human': ('human', 'hsTFs'), 'mouse': ('mouse', 'mmTFs'), 'celegans': ('celeg', 'ceTFs')}


def _read(name):
    return (files('PyStemFinder') / 'data' / f'{name}.txt').read_text().split()


def _species(species):
    if species not in _SPECIES:
        raise ValueError(f"Unknown species '{species}'. Expected one of {', '.join(_SPECIES)}.")
    return _SPECIES[species]


def cell_cycle_genes(species='mouse', phase='both'):
    """
    S and G2M phase cell cycle genes, the standard stemFinder markers.
    
    With phase='both', returns the S genes followed by the G2M genes, like c(s_genes_mouse, g2m_genes_mouse) in the 
    R vignette. Genes on both lists (E2f8 in mouse, E2F8 in human, 11 genes in C. elegans) therefore appear twice 
    and count twice in run_stemFinder, as in R. Use list(dict.fromkeys(genes)) to drop the repeats.
    
    Args:
        species (str, optional): 'human', 'mouse', or 'celegans'. Defaults to 'mouse'.
        phase (str, optional): 'S', 'G2M', or 'both'. Defaults to 'both'.
    
    Returns:
        list of str: Gene symbols.
    """
    suffix, _ = _species(species)
    if phase not in ['S', 'G2M', 'both']:
        raise ValueError("Invalid value for 'phase'. Expected one of 'S', 'G2M', 'both'.")
    s = _read(f's_genes_{suffix}') if phase in ['S', 'both'] else []
    g2m = _read(f'g2m_genes_{suffix}') if phase in ['G2M', 'both'] else []
    return s + g2m


def transcription_factors(species='mouse'):
    """
    Reference transcription factor list (hsTFs, mmTFs, or ceTFs in the R package), without repeats.
    
    Args:
        species (str, optional): 'human', 'mouse', or 'celegans'. Defaults to 'mouse'.
    
    Returns:
        list of str: Gene symbols.
    """
    _, name = _species(species)
    return list(dict.fromkeys(_read(name)))
