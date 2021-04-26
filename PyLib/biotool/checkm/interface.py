# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-02 22:32:18
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-29 12:41:53
 * @FilePath: /HScripts/Python/mylib/biotool/checkm/_checkm.py
 * @Description:
        Just like checkm.
"""

from ast import literal_eval


class DefaultValues():
    """Default values for filenames and common constants."""

    # set of markers recognized to be unreliable. These are often
    # ubiquitous, single-copy genes, but ones which are challenging
    # to correctly annotate with the PFAM and TIGRFAM models.
    MARKERS_TO_EXCLUDE = {'TIGR00398', 'TIGR00399'}

    E_VAL = 1e-10
    LENGTH = 0.7
    PSEUDOGENE_LENGTH = 0.3

    SEQ_CONCAT_CHAR = '&&'

    MARKER_FILE = r"lineage.ms"
    HMM_MODEL_INFO_FILE = r"checkm_hmm_info.pkl.gz"
    BIN_STATS_FILE = r"bin_stats.analyze.tsv"
    HMM_TABLE_FILE = r"hmmer.analyze.txt"


class BinMarkerSets():
    """A collection of one or more marker sets associated with a bin."""

    def __init__(self, binId):
        self.markerSets = []
        self.binId = binId
        self.selectedLinageSpecificMarkerSet = None

    def read(self, line):
        """Construct bin marker set data from line."""
        lineSplit = line.split('\t')
        numMarkerSets = int(lineSplit[1])
        for i in range(0, numMarkerSets):
            uid = lineSplit[i * 4 + 2]
            # lineageStr = lineSplit[i * 4 + 3]
            # numGenomes = int(lineSplit[i * 4 + 4])
            markerSet = literal_eval(lineSplit[i * 4 + 5])
            self.markerSets.append(MarkerSet(uid, markerSet))


class MarkerSet():
    """A collection of marker genes organized into co-located sets."""

    def __init__(self, UID, markerSet):
        self.UID = UID  # unique ID of marker set
        self.markerSet = markerSet  # marker genes organized into co-located sets

    def genomeCheck(self, hits, bIndividualMarkers):
        """Calculate genome completeness and contamination."""
        # self.markerSet: [{'PF04760.10', }, ]
        # hits: {'PF11874.3': [<__main__.HmmerHitDOM object at 0x000001E9CA1332C8>, ], }
        if bIndividualMarkers:
            present = 0
            multiCopyCount = 0
            for marker in self.getMarkerGenes():
                if marker in hits:
                    present += 1
                    multiCopyCount += (len(hits[marker]) - 1)

            percComp = 100 * float(present) / self.numMarkers()
            percCont = 100 * float(multiCopyCount) / self.numMarkers()
        else:
            comp = 0.0
            cont = 0.0
            # ms: {'PF04760.10', }
            for ms in self.markerSet:
                present = 0
                multiCopy = 0
                # marker: 'PF04760.10'
                for marker in ms:
                    # len([<__main__.HmmerHitDOM object at 0x000001E9CA1332C8>, ])
                    count = len(hits.get(marker, []))
                    if count == 1:
                        present += 1
                    elif count > 1:
                        present += 1
                        multiCopy += (count - 1)

                comp += float(present) / len(ms)
                cont += float(multiCopy) / len(ms)

            percComp = 100 * comp / len(self.markerSet)
            percCont = 100 * cont / len(self.markerSet)

        return percComp, percCont


class HmmerHitDOM():
    """Encapsulate a HMMER hit given in domtblout format."""

    def __init__(self, values):
        if len(values) == 23:
            self.target_name = values[0]
            self.target_accession = values[1]
            self.target_length = int(values[2])
            self.query_name = values[3]

            self.query_accession = values[4]
            if self.query_accession == '-':
                self.query_accession = self.query_name

            self.query_length = int(values[5])
            self.full_e_value = float(values[6])
            self.full_score = float(values[7])
            self.full_bias = float(values[8])
            self.dom = int(values[9])
            self.ndom = int(values[10])
            self.c_evalue = float(values[11])
            self.i_evalue = float(values[12])
            self.dom_score = float(values[13])
            self.dom_bias = float(values[14])
            self.hmm_from = int(values[15])
            self.hmm_to = int(values[16])
            self.ali_from = int(values[17])
            self.ali_to = int(values[18])
            self.env_from = int(values[19])
            self.env_to = int(values[20])
            self.acc = float(values[21])
            self.target_description = values[22]


class ResultsManager():
    """Store all the results for a single bin"""

    def __init__(self, binId, models,
                 bIgnoreThresholds=False,
                 evalueThreshold=DefaultValues.E_VAL,
                 lengthThreshold=DefaultValues.LENGTH,
                 bSkipPseudoGeneCorrection=False,
                 binStats=None):
        self.binId = binId
        self.markerHits = {}
        self.bIgnoreThresholds = bIgnoreThresholds
        self.evalueThreshold = evalueThreshold
        self.lengthThreshold = lengthThreshold
        self.bSkipPseudoGeneCorrection = bSkipPseudoGeneCorrection
        self.models = models
        self.binStats = binStats

    def vetHit(self, hit):
        """Check if hit meets required thresholds."""
        model = self.models[hit.query_accession]

        # preferentially use model specific bit score thresholds, before
        # using the user specified e-value and length criteria

        # Give preference to the gathering threshold unless the model
        # is marked as TIGR (i.e., TIGRFAM model)

        if not self.bSkipPseudoGeneCorrection:
            alignment_length = float(hit.ali_to - hit.ali_from)
            length_perc = alignment_length / float(hit.query_length)
            if length_perc < DefaultValues.PSEUDOGENE_LENGTH:
                return False

        if model.nc is not None and not self.bIgnoreThresholds and 'TIGR' in model.acc:
            if model.nc[0] <= hit.full_score and model.nc[1] <= hit.dom_score:
                return True
        elif model.ga is not None and not self.bIgnoreThresholds:
            if model.ga[0] <= hit.full_score and model.ga[1] <= hit.dom_score:
                return True
        elif model.tc is not None and not self.bIgnoreThresholds:
            if model.tc[0] <= hit.full_score and model.tc[1] <= hit.dom_score:
                return True
        elif model.nc is not None and not self.bIgnoreThresholds:
            if model.nc[0] <= hit.full_score and model.nc[1] <= hit.dom_score:
                return True
        else:
            if hit.full_e_value > self.evalueThreshold:
                return False

            alignment_length = float(hit.ali_to - hit.ali_from)
            length_perc = alignment_length / float(hit.query_length)
            if length_perc >= self.lengthThreshold:
                return True

        return False

    def addHit(self, hit):
        """Process hit and add it to the set of markers if it passes filtering criteria."""
        if self.vetHit(hit):
            if hit.query_accession in self.markerHits:
                # retain only the best domain hit for a given marker to a specific ORF
                previousHitToORF = None
                for h in self.markerHits[hit.query_accession]:
                    if h.target_name == hit.target_name:
                        previousHitToORF = h
                        break

                if not previousHitToORF:
                    self.markerHits[hit.query_accession].append(hit)
                else:
                    if previousHitToORF.dom_score < hit.dom_score:
                        self.markerHits[hit.query_accession].append(hit)
                        self.markerHits[hit.query_accession].remove(
                            previousHitToORF)

            else:
                self.markerHits[hit.query_accession] = [hit]

    def identifyAdjacentMarkerGenes(self):
        """Identify adjacent marker genes and exclude these from the contamination estimate."""

        # check for adjacent ORFs with hits to the same marker gene
        for markerId, hits in self.markerHits.items():

            bCombined = True
            while bCombined:
                for i in range(0, len(hits)):
                    orfI = hits[i].target_name
                    scaffoldIdI = orfI[0:orfI.rfind('_')]

                    bCombined = False
                    for j in range(i + 1, len(hits)):
                        orfJ = hits[j].target_name
                        scaffoldIdJ = orfJ[0:orfJ.rfind('_')]

                        # check if hits are on adjacent ORFs
                        if scaffoldIdI == scaffoldIdJ:
                            try:
                                orfNumI = int(orfI[orfI.rfind('_') + 1:])
                                orfNumJ = int(orfJ[orfJ.rfind('_') + 1:])
                            except Exception:
                                # it appears called genes are not labeled
                                # according to the prodigal format, so
                                # it is not possible to perform this correction
                                break

                            if abs(orfNumI - orfNumJ) == 1:
                                # check if hits are to different parts of the HMM
                                sI = hits[i].hmm_from
                                eI = hits[i].hmm_to

                                sJ = hits[j].hmm_from
                                eJ = hits[j].hmm_to

                                if (sI <= sJ and eI > sJ) or (sJ <= sI and eJ > sI):
                                    # models overlap so this could represent contamination,
                                    # but it seems more likely that adjacent genes hitting
                                    # the same marker represent legitimate gene duplication,
                                    # a gene calling error, or an assembly error and thus
                                    # should not be treated as contamination
                                    bCombined = True
                                    break
                                else:
                                    # combine the two hits
                                    bCombined = True
                                    break

                    if bCombined:
                        newHit = hits[i]

                        # produce concatenated label indicating the two genes being combined
                        orfA, orfB = sorted([orfI, orfJ])
                        newHit.target_name = DefaultValues.SEQ_CONCAT_CHAR.join([
                                                                                orfA, orfB])

                        newHit.target_length = hits[i].target_length + \
                            hits[j].target_length
                        newHit.hmm_from = min(
                            hits[i].hmm_from, hits[j].hmm_from)
                        newHit.hmm_to = min(hits[i].hmm_to, hits[j].hmm_to)

                        newHit.ali_from = min(
                            hits[i].ali_from, hits[j].ali_from)
                        newHit.ali_to = min(hits[i].ali_to, hits[j].ali_to)

                        newHit.env_from = min(
                            hits[i].env_from, hits[j].env_from)
                        newHit.env_to = min(hits[i].env_to, hits[j].env_to)

                        hits.remove(hits[j])
                        hits.remove(hits[i])

                        hits.append(newHit)

                        break

            self.markerHits[markerId] = hits


# # ResultsParser.printSummary
# # ResultsManager.printSummary
# # ResultsManager.geneCountsForSelectedMarkerSet
# # : markerSet, markerHits, False
# # : binMarkerSets.selectedMarkerSet(), ResultsManager.markerHits
# # markerSet: binIdToBinMarkerSets[binId].markerSets[?]
# # markerHits: resultsManager.markerHits
# # bIndividualMarkers = False
# for binId, resultsManager in resultsManagers.items():
#     # for markerSet in binIdToBinMarkerSets[binId].markerSets:
#     #     for markerset in markerSet.markerSet:
#     #         print("markerSet.markerSet:", markerset)
#     #     break
#     print("binId:", binId)
#     for markerSet in binIdToBinMarkerSets[binId].markerSets:
#     #     print("UID:", markerSet.UID, end = "\t")
#     #     print(markerSet.genomeCheck(resultsManager.markerHits, bIndividualMarkers))
#     #     print()
#         if markerSet.UID != "UID4444":
#             continue
#         print("markerSet.markerSet:", markerSet.markerSet)
#         print("resultsManager.markerHits['PF11874.3'][0].__dict__:", resultsManager.markerHits['PF11874.3'][0].__dict__)
#     break
