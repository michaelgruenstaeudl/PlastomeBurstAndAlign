from datetime import date
from importlib.metadata import version

__name__ = "plastburstalign"
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "m_gruenstaeudl@fhsu.edu"

try:
    __version__ = version(__name__)
except ModuleNotFoundError:
    __version__ = date.today()

from .user_parameters import UserParameters
from .seqfeature_ops import PlastidData
from .extraction_ops import ExtractAndCollect, DataCleaning
from .alignment_ops import AlignmentCoordination, MAFFT
from .plastome_burst_and_align import PlastomeRegionBurstAndAlign

__all__ = ['PlastomeRegionBurstAndAlign', 'UserParameters', 'PlastidData',
           'ExtractAndCollect', 'DataCleaning', 'AlignmentCoordination', 'MAFFT']
