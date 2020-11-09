# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-06 21:57:58
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-11-09 22:28:03
 * @FilePath: /HScripts/Python/mylib/biotool/read_outputs.py
 * @Description:
"""

from io import StringIO
from sys import stderr
from typing import Any, Callable
from numpy import nan
from mylib.biotool.checkm.reload import reload_checkMOutput


checkM = reload_checkMOutput


def gtdbtk(text: StringIO) -> list:
    """ read --out_dir/gtdbtk.bac120.summary.tsv
                     ./gtdbtk.ar122.summary.tsv
       @return:
            user_genome, classification, fastani_reference, fastani_reference_radius, fastani_taxonomy, fastani_ani, fastani_af, closest_placement_reference, closest_placement_radius, closest_placement_taxonomy, closest_placement_ani, closest_placement_af, pplacer_taxonomy, classification_method, note, other_related_references, aa_percent, translation_table, red_value, warnings
       @warnings: some line are need to deal with.
    """
    gtdbtklist = []
    text.readline()
    print("""
    raise FutureWarning("Some lines are still raw. please check carefully")
    #""", file=stderr)
    for line in text:
        user_genome, classification, fastani_reference, fastani_reference_radius, fastani_taxonomy, fastani_ani, fastani_af, closest_placement_reference, closest_placement_radius, closest_placement_taxonomy, closest_placement_ani, closest_placement_af, pplacer_taxonomy, classification_method, note, other_related_references, aa_percent, translation_table, red_value, warnings = line[:-1].split("\t")
        gtdbtklist.append([
            user_genome, classification, fastani_reference, fastani_reference_radius, fastani_taxonomy, fastani_ani, fastani_af, closest_placement_reference, closest_placement_radius, closest_placement_taxonomy, closest_placement_ani, closest_placement_af, pplacer_taxonomy, classification_method, note, other_related_references, aa_percent, translation_table, red_value, warnings
        ])
    return gtdbtklist


def read_nan(dtype: type, numstr: str, *nanstr) -> Any:
    """ read args, if arg is nan, return nan
     * @param dtype     type of return
     * @param numstr    number to transfer
                        if numstr in nan: num = str
     * @return type, nan
    """
    if numstr in nanstr:
        return nan
    return dtype(numstr)


def iRep(text: StringIO) -> list:
    """ read -o.tsv
       @return:
            binId | index of replication (iRep) | un-filtered index of replication (iRep) | raw index of replication (no GC bias correction) | r^2 | coverage | % windows passing filter | fragments/Mbp | GC bias | GC r^2
            n/a will be recorded to np.nan
    """
    def read_values(text: StringIO, head: str, value_type: type,
                    recode_func: Callable[[str, Any], None]) -> str:
        '''
            ../03_modify/7_final//metabat2_90_90.5.fa	n/a
            #
        --> ## r^2
            # genome	E001.sam
            ../03_modify/7_final//17.fa	n/a
            ../03_modify/7_final//4.fa	0.907207353884764
            ../03_modify/7_final//55_sub22.fa	0.9151732818669122
        '''
        line = text.readline()
        assert line[:len(head) + 3] + "## " + head
        line = text.readline()  # genome name
        for line in text:
            if line[0] == "#":
                return line
            genome_dir, raw_value = line.strip().split()
            # ../03_modify/7_final//17.fa	n/a
            binId = genome_dir[genome_dir.rfind("/")+1:genome_dir.rfind(".")]
            value = read_nan(value_type, raw_value, "n/a", "False")
            recode_func(binId, value)
    # collect data to this map
    irmap = {}
    irvalues = {
        'index of replication'            : float,
        'un-filtered index of replication': float,
        'raw index of replication'        : float,
        'r^2'                             : float,
        'coverage'                        : float,
        '% windows passing filter'        : float,
        'fragments/Mbp'                   : int  ,
        'GC bias'                         : float,
        'GC r^2'                          : float,
    }
    heads_num = len(irvalues)
    def index_lambda(binId: str, value: Any, index: int) -> None:
        if binId not in irmap:
            irmap[binId] = [nan for _ in range(heads_num)]
        irmap[binId][index] = value
    for index, (name, dtype) in enumerate(irvalues.items()):
        read_values(text, name, dtype, lambda p0, p1: index_lambda(p0, p1, index))
    return [[binId] + values for binId, values in irmap.items()]
