# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-02 22:11:16
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-06-14 19:53:45
 * @FilePath: /metaSC/PyLib/biotool/kegg/amino_acid_metabolism.py
 * @Description:
"""


from PyLib.biotool.kegg.kegg import _load_module, KModule


amino_acid_metabolism = """
    Serine and threonine metabolism
    Entry	M00020
    Name	Serine biosynthesis, glycerate-3P => serine
    Definition	K00058 K00831 (K01079,K02203,K22305)
    Entry	M00018
    Name	Threonine biosynthesis, aspartate => homoserine => threonine
    Definition	(K00928,K12524,K12525,K12526) K00133 (K00003,K12524,K12525) (K00872,K02204,K02203) K01733
    Entry	M00555
    Name	Betaine biosynthesis, choline => betaine
    Definition	(K17755,((K00108,K11440,K00499) (K00130,K14085)))
    Entry	M00033
    Name	Ectoine biosynthesis, aspartate => ectoine
    Definition	K00928 K00133 K00836 K06718 K06720

    Cysteine and methionine metabolism
    Entry	M00021
    Name	Cysteine biosynthesis, serine => cysteine
    Definition	(K00640,K23304) (K01738,K13034,K17069)
    Entry	M00338
    Name	Cysteine biosynthesis, homocysteine + serine => cysteine
    Definition	(K01697,K10150) K01758
    Entry	M00609
    Name	Cysteine biosynthesis, methionine => cysteine
    Definition	K00789 K17462 K01243 K07173 K17216 K17217
    Entry	M00017
    Name	Methionine biosynthesis, apartate => homoserine => methionine
    Definition	(K00928,K12524,K12525) K00133 (K00003,K12524,K12525) (K00651,K00641) K01739 (K01760,K14155) (K00548,K24042,K00549)
    Entry	M00034
    Name	Methionine salvage pathway
    Definition	K00789 K01611 K00797 ((K01243,K01244) K00899,K00772) K08963 (K16054,K08964 (K09880,K08965 K08966)) K08967 (K00815,K08969,K23977,K00832,K00838)
    Entry	M00035
    Name	Methionine degradation
    Definition	K00789 (K00558,K17398,K17399) K01251 (K01697,K10150)
    Entry	M00368
    Name	Ethylene biosynthesis, methionine => ethylene
    Definition	K00789 (K01762,K20772) K05933

    Branched-chain amino acid metabolism
    Entry	M00019
    Name	Valine/isoleucine biosynthesis, pyruvate => valine / 2-oxobutanoate => isoleucine
    Definition	K01652+(K01653,K11258) K00053 K01687 K00826
    Entry	M00535
    Name	Isoleucine biosynthesis, pyruvate => 2-oxobutanoate
    Definition	K09011 K01703+K01704 K00052
    Entry	M00570
    Name	Isoleucine biosynthesis, threonine => 2-oxobutanoate => isoleucine
    Definition	(K17989,K01754) K01652+(K01653,K11258) K00053 K01687 K00826
    Entry	M00432
    Name	Leucine biosynthesis, 2-oxoisovalerate => 2-oxoisocaproate
    Definition	K01649 (K01702,K01703+K01704) K00052
    Entry	M00036
    Name	Leucine degradation, leucine => acetoacetate + acetyl-CoA
    Definition	K00826 ((K00166+K00167,K11381)+K09699+K00382) (K00253,K00249) (K01968+K01969) (K05607,K13766) K01640

    Lysine metabolism
    Entry	M00016
    Name	Lysine biosynthesis, succinyl-DAP pathway, aspartate => lysine
    Definition	(K00928,K12524,K12525,K12526) K00133 K01714 K00215 K00674 (K00821,K14267) K01439 K01778 (K01586,K12526)
    Entry	M00525
    Name	Lysine biosynthesis, acetyl-DAP pathway, aspartate => lysine
    Definition	K00928 K00133 K01714 K00215 K05822 K00841 K05823 K01778 K01586
    Entry	M00526
    Name	Lysine biosynthesis, DAP dehydrogenase pathway, aspartate => lysine
    Definition	(K00928,K12524,K12525,K12526) K00133 K01714 K00215 K03340 (K01586,K12526)
    Entry	M00527
    Name	Lysine biosynthesis, DAP aminotransferase pathway, aspartate => lysine
    Definition	(K00928,K12524,K12525,K12526) K00133 K01714 K00215 K10206 K01778 (K01586,K12526)
    Entry	M00030
    Name	Lysine biosynthesis, AAA pathway, 2-oxoglutarate => 2-aminoadipate => lysine
    Definition	K01655 K17450 K01705 K05824 K00838 K00143 (K00293,K24034) K00290
    Entry	M00433
    Name	Lysine biosynthesis, 2-oxoglutarate => 2-oxoadipate
    Definition	K01655 (K17450 K01705,K16792+K16793) K05824
    Entry	M00031
    Name	Lysine biosynthesis, mediated by LysW, 2-aminoadipate => lysine
    Definition	K05826,(K05827 K05828 K05829 K05830 K05831)
    Entry	M00032
    Name	Lysine degradation, lysine => saccharopine => acetoacetyl-CoA
    Definition	K14157 K14085 K00825 (K15791+K00658+K00382) K00252 (K07514,(K07515,K07511) K00022)

    Arginine and proline metabolism
    Entry	M00028
    Name	Ornithine biosynthesis, glutamate => ornithine
    Definition	(K00618,K00619,K14681,K14682,K00620,K22477,K22478) ((K00930,K22478) K00145,K12659) (K00818,K00821) (K01438,K14677,K00620)
    Entry	M00763
    Name	Ornithine biosynthesis, mediated by LysW, glutamate => ornithine
    Definition	K05826,(K19412 K05828 K05829 K05830 K05831)
    Entry	M00844
    Name	Arginine biosynthesis, ornithine => arginine
    Definition	K00611 K01940 (K01755,K14681)
    Entry	M00845
    Name	Arginine biosynthesis, glutamate => acetylcitrulline => arginine
    Definition	K22478 K00145 K00821 K09065 K01438 K01940 K01755
    Entry	M00029
    Name	Urea cycle
    Definition	K01948 K00611 K01940 (K01755,K14681) K01476
    Entry	M00015
    Name	Proline biosynthesis, glutamate => proline
    Definition	((K00931 K00147),K12657) K00286
    Entry	M00047
    Name	Creatine pathway
    Definition	K00613 K00542 K00933
    Entry	M00879
    Name	Arginine succinyltransferase pathway, arginine => glutamate
    Definition	K00673 K01484 K00840 K06447 K05526

    Polyamine biosynthesis
    Entry	M00133
    Name	Polyamine biosynthesis, arginine => agmatine => putrescine => spermidine
    Definition	((K01583,K01584,K01585,K02626) K01480),K01611,K00797
    Entry	M00134
    Name	Polyamine biosynthesis, arginine => ornithine => putrescine
    Definition	K01476 K01581
    Entry	M00135
    Name	GABA biosynthesis, eukaryotes, putrescine => GABA
    Definition	K00657 K00274 (K00128,K14085,K00149) --
    Entry	M00136
    Name	GABA biosynthesis, prokaryotes, putrescine => GABA
    Definition	K09470 K09471 K09472 K09473

    Histidine metabolism
    Entry	M00026
    Name	Histidine biosynthesis, PRPP => histidine
    Definition	(K00765-K02502) (K01523 K01496,K11755,K14152) (K01814,K24017) (K02501+K02500,K01663) ((K01693 K00817 (K04486,K05602,K18649)),(K01089 K00817)) (K00013,K14152)
    Entry	M00045
    Name	Histidine degradation, histidine => N-formiminoglutamate => glutamate
    Definition	K01745 K01712 K01468 (K01479,K00603,K13990,(K05603 K01458))

    Aromatic amino acid metabolism
    Entry	M00022
    Name	Shikimate pathway, phosphoenolpyruvate + erythrose-4P => chorismate
    Definition	(K01626,K03856,K13853) (((K01735,K13829) ((K03785,K03786) K00014,K13832)),K13830) ((K00891,K13829) (K00800,K24018),K13830) K01736
    Entry	M00023
    Name	Tryptophan biosynthesis, chorismate => tryptophan
    Definition	(((K01657+K01658,K13503,K13501,K01656) K00766),K13497) (((K01817,K24017) (K01656,K01609)),K13498,K13501) (K01695+(K01696,K06001),K01694)
    Entry	M00024
    Name	Phenylalanine biosynthesis, chorismate => phenylpyruvate => phenylalanine
    Definition	((K01850,K04092,K14187,K04093,K04516,K06208,K06209)
                (K01713,K04518,K05359),K14170) (K00832,K00838)
    Entry	M00910
    Name	Phenylalanine biosynthesis, chorismate => arogenate => phenylalanine
    Definition	K01850 K15849 K05359
    Entry	M00025
    Name	Tyrosine biosynthesis, chorismate => HPP => tyrosine
    Definition	(((K01850,K04092,K14170,K04093,K04516,K06208,K06209)
                (K04517,K00211)),K14187) (K00832,K00838)
    Entry	M00040
    Name	Tyrosine biosynthesis, chorismate => arogenate => tyrosine
    Definition	(K01850,K04092,K14170) (K00832,K15849) (K00220,K24018,K15227)
    Entry	M00042
    Name	Catecholamine biosynthesis, tyrosine => dopamine => noradrenaline => adrenaline
    Definition	(K00505,K00501) K01593 K00503 K00553
    Entry	M00043
    Name	Thyroid hormone biosynthesis, tyrosine => triiodothyronine/thyroxine
    Definition	K00431
    Entry	M00044
    Name	Tyrosine degradation, tyrosine => homogentisate
    Definition	(K00815,K00838,K00832,K03334) K00457 K00451 K01800 (K01555,K16171)
    Entry	M00533
    Name	Homoprotocatechuate degradation, homoprotocatechuate => 2-oxohept-3-enedioate
    Definition	K00455 K00151 K01826 K05921
    Entry	M00545
    Name	Trans-cinnamate degradation, trans-cinnamate => acetyl-CoA
    Definition	(((K05708+K05709+K05710+K00529) K05711),K05712) K05713 K05714 K02554 K01666 K04073
    Entry	M00037
    Name	Melatonin biosynthesis, tryptophan => serotonin => melatonin
    Definition	K00502 K01593 K00669 K00543
    Entry	M00038
    Name	Tryptophan metabolism, tryptophan => kynurenine => 2-aminomuconate
    Definition	(K00453,K00463) (K01432,K14263,K07130) K00486 K01556 K00452 K03392 (K10217,K23234)

    Other amino acid metabolism
    Entry	M00027
    Name	GABA (gamma-Aminobutyrate) shunt
    Definition	K01580 (K13524,K07250,K00823,K16871) (K00135,K00139,K17761)
    Entry	M00369
    Name	Cyanogenic glycoside biosynthesis, tyrosine => dhurrin
    Definition	K13027 K13029 K13030
    Entry	M00118
    Name	Glutathione biosynthesis, glutamate => glutathione
    Definition	(K11204+K11205,K01919) (K21456,K01920)
"""


AAm = _load_module(amino_acid_metabolism)
del amino_acid_metabolism


def Amino_acid_metabolism(ko_match):
    """A interface. """
    out_data = {}
    for pathway in AAm:
        metabolism_data = {}
        for metabolism, element in AAm[pathway].items():
            value = element.completeness(ko_match)
            if value:
                metabolism_data[metabolism] = value
        out_data.update(metabolism_data)
    return out_data


def test():
    def echo(express):
        """repeat the express"""
        print(express)
        km = KModule(express)
        print(km)
        print(len(km))
    print("test begin:")
    echo("K00058 K00831 (K01079,K02203,K22305)")
    echo("(K00928,K12524,K12525,K12526) K00133 (K00003,K12524,K12525) (K00872,K02204,K02203) K01733")
    echo("(K17755,((K00108,K11440,K00499) (K00130,K14085)))")
    echo("(K00640,K23304) (K01738,K13034,K17069)")
    echo("K00826 ((K00166+K00167,K11381)+K09699+K00382) (K00253,K00249) (K01968+K01969) (K05607,K13766) K01640")
    echo("K09011 K01703+K01704 K00052")
    echo("(K00765-K02502) (K01523 K01496,K11755,K14152) (K01814,K24017) (K02501+K02500,K01663) ((K01693 K00817 (K04486,K05602,K18649)),(K01089 K00817)) (K00013,K14152)")
    echo("K00455 K00151 K01826 K05921")
    print(KModule("K00058 K00831 (K01079,K02203,K22305)")["K00831"])
    print(KModule("K00058 K00831 (K01079,K02203,K22305)")["K02079"])
    print(KModule(
        "K00826 ((K00166+K00167,K11381)+K09699+K00382) (K00253,K00249) (K01968+K01969) (K05607,K13766) K01640"
    )["K00382"])
    print("test ends. please check the echo")

    print(Amino_acid_metabolism({"K00831": 1, "K00133": 1, "K12524": 1}))


if __name__ == "__main__":
    test()
