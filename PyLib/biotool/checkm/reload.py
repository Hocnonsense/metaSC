# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-24 12:55:39
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-12-11 20:25:49
 * @FilePath: /HScripts/Python/mylib/biotool/checkm/reload.py
 * @Description:
    Reload from checkm output.
"""


import gzip
import os
import pickle
from ast import literal_eval
from collections import OrderedDict
from re import split as re_split
from sys import stderr

from mylib.biotool.checkm.interface import (
    BinMarkerSets, DefaultValues, HmmerHitDOM, ResultsManager
)


def popBinIdFromSets(pop_dict: dict, pop_list: list) -> None:
    """Pop binId from pop_dict by pop_list"""
    for binId in pop_list:
        pop_dict.pop(binId)


def reload_binIdToBinMarkerSets(checkM_OUTPUT_DIR):
    """ load binIdToBinMarkerSets. """
    # read file like these:
    """
        # [Lineage Marker File]
        25_4	2	UID1	root	5656	[{'PF00411.14', 'PF01000.21'}, {'PF00831.18'}]	UID1	root	5656	[{'PF00411.14'}]
        25_3	1	UID1	root	5656	[{'PF00411.14'}]
        """
    binIdToBinMarkerSets = {}
    with open(os.path.join(checkM_OUTPUT_DIR, DefaultValues.MARKER_FILE)) as f:
        f.readline()  # skip header

        for line in f:
            lineSplit = line.split('\t')
            binId = lineSplit[0]

            binMarkerSets = BinMarkerSets(binId)
            binMarkerSets.read(line)

            binIdToBinMarkerSets[binId] = binMarkerSets

        # check point #
        # print("binIdToBinMarkerSets[25_1].markerSets[0].UID:", binIdToBinMarkerSets[r"25_1"].markerSets[0].UID)
        # print("binIdToBinMarkerSets[25_1].markerSets[0].markerSet[0]:", binIdToBinMarkerSets[r"25_1"].markerSets[0].markerSet[0])
    return binIdToBinMarkerSets


def reload_resultsManagers(checkM_OUTPUT_DIR, binIdToBinMarkerSets):
    # loadBinModels
    binIdToModels = {}
    with gzip.open(os.path.join(checkM_OUTPUT_DIR, 'storage', DefaultValues.HMM_MODEL_INFO_FILE)) as f:
        binIdToModels = pickle.load(f)
        # print("binIdToModels.keys():", binIdToModels.keys())
        # sth in binIdToBinMarkerSets but not in binIdToModels
        popBinIds = []
        for binId in binIdToBinMarkerSets:
            if binId not in binIdToModels:
                popBinIds.append(binId)
        popBinIdFromSets(binIdToBinMarkerSets, popBinIds)

    # parseBinStats
    binStats = {}
    with open(os.path.join(checkM_OUTPUT_DIR, 'storage', DefaultValues.BIN_STATS_FILE)) as f:
        for line in f:
            lineSplit = line.split('\t')
            binStats[lineSplit[0]] = literal_eval(lineSplit[1])
        popBinIds = []
        for binId in binIdToBinMarkerSets:
            if binId not in binIdToModels:
                popBinIds.append(binId)
        popBinIdFromSets(binIdToBinMarkerSets, popBinIds)

    # read file like these:
    """
        #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
        # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
        #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
        scaffold_4712_4      -            141 DUF3897              PF13036.1    180   1.2e-44  150.3   5.8   1   1   1.3e-45   9.4e-43  144.2   5.8    62   180     2   129     1   129 0.89 # 2503 # 2925 # -1 # ID=364_4;partial=01;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.584
        """
    resultsManagers = {}
    popBinIds = []
    for binId in binIdToBinMarkerSets:
        resultsManager = ResultsManager(
            binId, binIdToModels[binId], binStats=binStats[binId])

        hmmTableFile = os.path.join(
            checkM_OUTPUT_DIR, 'bins', binId, DefaultValues.HMM_TABLE_FILE)
        # parseHmmerResults
        # print("hmmTableFile:", hmmTableFile)
        try:
            f = open(hmmTableFile)
            line_count = 0
            while True:
                # HMMERParser.next()
                # readHitsDOM
                while True:
                    line = f.readline().rstrip()
                    try:
                        if line[0] != '#' and len(line) != 0:
                            dMatch = re_split(r'\s+', line.rstrip())
                            if len(dMatch) < 23:
                                raise Exception(
                                    "Error processing line:\n%s" % (line))
                            refined_match = dMatch[0:22] + \
                                [" ".join([str(i) for i in dMatch[22:]])]
                            hit = HmmerHitDOM(refined_match)
                            line_count += 1
                            break
                    except IndexError:
                        hit = None
                        break
                if hit is None:
                    break
                print("read {0} hits: {1}, {2} in bin: {3}".format(
                    line_count, hit.target_name, hit.query_name, binId), end="\r", file=stderr)
                resultsManager.addHit(hit)
        except FileNotFoundError:
            popBinIds.append(binId)
            continue

        resultsManager.identifyAdjacentMarkerGenes()
        resultsManagers[binId] = resultsManager
        print("", file=stderr)
        # print("resultsManager.__dict__:", resultsManager.__dict__.keys())
    popBinIdFromSets(binIdToBinMarkerSets, popBinIds)
    return resultsManagers


def reload_checkMOutput(checkM_OUTPUT_STDOUT):
    """
    * @description: read checkm ouput and get message there
    * @param {str} checkM_OUTPUT_STDOUT
    * @return {OrderedDict} ckmap
        ckmap -> OrderedDict{Bin_Id:[
            Marker_lineage, UID,
            genomes, markers, marker sets, x0, x1, x2, x3, x4, x5plus,
            Completeness, Contamination, Strain_heterogeneity
        ]}
    """
    ckmap = OrderedDict()
    with open(checkM_OUTPUT_STDOUT) as fin:
        for line in fin:
            #print(line)
            if line.strip().split()[0] == "Bin":
                break
        else:
            print("""
                we read the file like:
                    '''
                    [2020-09-19 16:44:36] INFO: CheckM v1.1.2
                    [2020-09-19 16:44:36] INFO: checkm lineage ...
                    ...                                        ...
                    [2020-09-19 17:33:36] INFO: Parsing HMM hi ...
                    ------------------------------------------ ...
                 ->   Bin Id                            Marker ...
                    ------------------------------------------ ...
                      metabat2_90_60.124          o__Cytophaga ...
                    '''
                and we just reach the line of '->' table head

                So, what's your problem?
                """, file=stderr)
            raise Exception("Havn't found expected list")
        fin.readline()  # read the secend line between table body and head
        for line in fin:
            if line[0] == "-":
                break
            # Bin Id | Marker lineage | (UID) | genomes | markers | marker sets | 0-5+ | Completeness | Contamination | Strain heterogeneity
            values = line.strip().split()
            # Bin Id: [Marker lineage (UID), Completeness, Contamination]
            ckmap[values[0]] = [
                values[1], values[2][1:-1],
                *[int(values[i]) for i in range(3, 11)],
                *[float(values[i]) for i in range(12, 15)],
            ]
        else:
            print("""
                the table will end with a line of "-"

                So, what's your problem?
                """, file=stderr)
    return ckmap
