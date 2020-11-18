import os


#: Absolute path to output directory.
OUT_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', 'out')
)
#: AcrAB-TolC variable name and key in wcEcoli.
ACRAB_TOLC_KEY = 'TRANS-CPLX-201[s]'
#: AmpC variable name and key in wcEcoli.
BETA_LACTAMASE_KEY = 'EG10040-MONOMER[p]'
#: Variable name for antibiotic.
ANTIBIOTIC_KEY = 'nitrocefin'
#: Path to agents store from hierarchy root.
AGENTS_PATH = ('agents',)
#: Path to volume variable from each agent's root store.
VOLUME_PATH = ('boundary', 'volume')
#: Path to mass variable from each agent's root store.
MASS_PATH = ('boundary', 'mass')
#: Path to store of environmental fields.
FIELDS_PATH = ('fields',)
#: Path to environment bounds variable
BOUNDS_PATH = ('dimensions', 'bounds')
