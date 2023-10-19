import xml.dom.minidom
from _typeshed import Incomplete as Incomplete
from io import FileIO
from typing import Generator, List, Tuple

BLAST_TSV_specifiers: Incomplete
BLAST_TSV_DEFAULT_keywords: Incomplete
BLAST_XML_Hit_TITLE: Incomplete
BLAST_XML_Hsp_TITLE: Incomplete

def getEle(xmlEle: xml.dom.minidom.Element, TagName): ...
def blast_xml_iter(blast_file: FileIO): ...
def blast_unmatch_iter(values: List[str], wind: int = ...) -> Generator[Tuple[List[str], List[int], List[str]], None, None]: ...
