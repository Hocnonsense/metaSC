from _typeshed import Incomplete

class DefaultValues:
    MARKERS_TO_EXCLUDE: Incomplete
    E_VAL: float
    LENGTH: float
    PSEUDOGENE_LENGTH: float
    SEQ_CONCAT_CHAR: str
    MARKER_FILE: str
    HMM_MODEL_INFO_FILE: str
    BIN_STATS_FILE: str
    HMM_TABLE_FILE: str

class BinMarkerSets:
    markerSets: Incomplete
    binId: Incomplete
    selectedLinageSpecificMarkerSet: Incomplete
    def __init__(self, binId) -> None: ...
    def read(self, line) -> None: ...

class MarkerSet:
    UID: Incomplete
    markerSet: Incomplete
    def __init__(self, UID, markerSet) -> None: ...
    def getMarkerGenes(self) -> None: ...
    def numMarkers(self) -> None: ...
    def genomeCheck(self, hits, bIndividualMarkers): ...

class HmmerHitDOM:
    target_name: Incomplete
    target_accession: Incomplete
    target_length: Incomplete
    query_name: Incomplete
    query_accession: Incomplete
    query_length: Incomplete
    full_e_value: Incomplete
    full_score: Incomplete
    full_bias: Incomplete
    dom: Incomplete
    ndom: Incomplete
    c_evalue: Incomplete
    i_evalue: Incomplete
    dom_score: Incomplete
    dom_bias: Incomplete
    hmm_from: Incomplete
    hmm_to: Incomplete
    ali_from: Incomplete
    ali_to: Incomplete
    env_from: Incomplete
    env_to: Incomplete
    acc: Incomplete
    target_description: Incomplete
    def __init__(self, values) -> None: ...

class ResultsManager:
    binId: Incomplete
    markerHits: Incomplete
    bIgnoreThresholds: Incomplete
    evalueThreshold: Incomplete
    lengthThreshold: Incomplete
    bSkipPseudoGeneCorrection: Incomplete
    models: Incomplete
    binStats: Incomplete
    def __init__(self, binId, models, bIgnoreThresholds: bool = ..., evalueThreshold=..., lengthThreshold=..., bSkipPseudoGeneCorrection: bool = ..., binStats: Incomplete | None = ...) -> None: ...
    def vetHit(self, hit): ...
    def addHit(self, hit) -> None: ...
    def identifyAdjacentMarkerGenes(self) -> None: ...
