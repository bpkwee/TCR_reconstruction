import numpy as np

# Nomenclature and Symbolism for Amino Acids and Peptides. IUPAC-IUB Joint Commission on Biochemical Nomenclature. 1983.
AA_DICT = ["A", "C", "D", "E",
           "F", "G", "H", "I",
           "K", "L", "M", "N",
           "P", "Q", "R", "S",
           "T", "V", "W", "Y",
           "B", "X", "Z", "J",
           "U", "O", "*", "-"]

TRANSLATION_TABLE = {
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'TGG': 'W',
    'TGT': 'C', 'TGC': 'C',
    'TTT': 'F', 'TTC': 'F',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y',
    'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'ATG': 'M',
    # Stop codons
    'TAG': '*', 'TAA': '*', 'TGA': '*'}
VDJ_NAME = ['TRAV', 'TRAJ', 'TRBV', 'TRBD', 'TRBJ']

# cellranger
START_CODON = 'ATG'
STOP_CODONS = ['TAG', 'TAA', 'TGA']
AMBIGUOUS_AA_CODE = 'X'
VDJ_CDR3_COMMON_END_MOTIFS = ['FGXG', 'WGXG']
VDJ_CDR3_RARE_END_MOTIFS = ['XGXG', 'FXXG']
VDJ_CDR3_ALL_END_MOTIFS = VDJ_CDR3_COMMON_END_MOTIFS + VDJ_CDR3_RARE_END_MOTIFS


# Max possible CDR3 length (in nts)
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4998859/
VDJ_MAX_CDR3_LEN = 80
VDJ_MIN_CDR3_LEN = 26


# TODO: ask mireille for the publication explaining why subtype is important
# TODO: check all annotations and re-check all comments made on the entries of the HLA dict
# TODO: make an option to save and load the dicts instead of defining here
# translation (from either 'HLA-allele' or 'expert assigned type' to [expert assigned type, 4 digit, MHCI or MHCII])
# via https://www.ebi.ac.uk/ipd/imgt/hla/dictionary.html
HLA_DICT = {'HLA-A*02': ['A2', np.nan, 'MHCI'],  # A*02, HLA-A*02 or just A2
            'HLA-A*01': ['A1', np.nan, 'MHCI'],
            'HLA-DQw1': [np.nan, np.nan, 'MHCII'],
            'HLA-B*14': ['B14', np.nan, 'MHCI'],  # B64, B65 or  B14
            'nan': [np.nan, np.nan, np.nan],
            'HLA-A2': ['A2', np.nan, 'MHCI'],
            'HLA-Cw* 16:01': [np.nan, np.nan, 'MHCI'],  # not assigned/ undefined
            'HLA-A*011': ['A11', np.nan, 'MHCI'],  # A*11? A*011 not in database
            'DR3*02:02': ['DR52', 'DRB3*02:02', 'MHCII'],
            # DR3*02:02. Deduced it must have been DRB3 and DR52, but not sure
            'HLA-A*2:01': ['A2', 'A*02:01', 'MHCI'],
            'HLA-B*8': ['B8', np.nan, 'MHCI'],  # B81 or B82 or B8
            'HLA-A*02:01': ['A2', 'A*02:01', 'MHCI'],
            'DRB1*04:01': ['DR4', 'DRB1*04:01', 'MHCII'],
            'HLA-B*44:05': ['B44', 'B*44:05', 'MHCI'],
            'HLA-B*57:01': ['B57', 'B*57:01', 'MHCI'],
            'HLA-B*57:03': ['B57', 'B*57:03', 'MHCI'],
            'HLA-B*07:02': ['B7', 'B*07:02', 'MHCI'],
            'HLA-A*01:01': ['A1', 'A*01:01', 'MHCI'],
            'HLA-A*4:04,*13:01': ['B13', np.nan, 'MHCI'],  # A4 does not exist
            'HLA-A*2:01,*6:801': ['A2', np.nan, 'MHCI'],  # A2 or A6: are almost same
            'Complete HLA typing is provided in the paper for each patient': [np.nan, np.nan, np.nan],
            'HLA- DQ5': ['DQ5', np.nan, 'MHCII'],  # not in database
            'HLA-B08': ['B8', np.nan, 'MHCI'],  # can only find with B*08
            'non HLA-A*02, non HLA-B08': [np.nan, np.nan, np.nan],
            'HLA-DPB1*02:01': [None, None, 'MHCII'],  # not in database
            'HLA-B*35:02': ['B35', 'B*35:02', 'MHCI'],
            'HLA-DRB1*04:01': ['DR4', 'DRB1*04:01', 'MHCII'],
            'HLA-DQ2': ['DQ2', np.nan, 'MHCII'],
            'H-2Kb': [None, np.nan, 'MHCII'],  # not in database
            'HLA-A1': ['A1', np.nan, 'MHCI'],  # cannot find A1, A*1 gives A11, A*01 gives A1
            'HLA-B7': ['B7', np.nan, 'MHCI'],  # cannot find B7, B*7 gives B73, only B*07 gives B7
            'DRB1*15:03': ['DR15', 'DRB1*15:03', 'MHCII'],
            'HLA-DQ8': ['DQ8', np.nan, 'MHCII'],
            'DQ8-trans': [None, np.nan, 'MHCII'],  # not in database, Do I use DQ8 ?
            'HLA-A*02:02': ['A2', 'A*02:02', 'MHCI'],
            'HLA-A*02:03': ['A2', 'A*02:03', 'MHCI'],
            'HLA-A*02:04': ['A2', 'A*02:04', 'MHCI'],
            'HLA-A*02:05': ['A2', 'A*02:05', 'MHCI'],
            'HLA-A*02:06': ['A2', 'A*02:06', 'MHCI'],
            'HLA-A*02:07': ['A2', 'A*02:07', 'MHCI'],
            'HLA-A*02:08': ['A2', 'A*02:08', 'MHCI'],
            'HLA-A*02:09': ['A2', 'A*02:09', 'MHCI'],
            'HLA-A*02:10': ['A2', 'A*02:10', 'MHCI'],
            'HLA-A*02:11': ['A2', 'A*02:11', 'MHCI'],
            'HLA-A*02:12': ['A2', 'A*02:12', 'MHCI'],
            'HLA-A*02:13': ['A2', 'A*02:13', 'MHCI'],
            'HLA-A*02:14': ['A2', 'A*02:14', 'MHCI'],
            'HLA-A*02:15': ['A2', 'A*02:15', 'MHCI'],
            'HLA-A*02:16': ['A2', 'A*02:16', 'MHCI'],
            'HLA-A*02:17': ['A2', 'A*02:17', 'MHCI'],
            'HLA-DR1': ['DR1', np.nan, 'MHCII'],
            'HLA-DR11': ['DR11', np.nan, 'MHCII'],
            'HLA-DR15': ['DR15', np.nan, 'MHCII'],
            'HLA-DR5': [None, np.nan, 'MHCII'],  # not clear: expert -> DR11, DR8 DR12 !?
            'HLA-DR4': ['DR4', np.nan, 'MHCII'],
            'HLA-DQ2.5': [None, np.nan, 'MHCII'],  # DQ2?? not clear
            'HLA-C*07:02': ['Cw7', 'C*07:02', 'MHCI'],
            'HLA-DRB1': ['DR1', np.nan, 'MHCII'],
            'DRB1*04-01': ['DR4', 'DRB1*04:01', 'MHCII'],
            '\xa0HLA-DQB1*06:02': ['DQ6', 'DQB1*06:02', 'MHCII'],  # what is this
            'HLA-A2:01': ['A2', 'A*02:01', 'MHCI'],
            'HLA-B*08': ['B8', np.nan, 'MHCI'],
            'HLA-A*03': ['A3', np.nan, 'MHCI'],
            'HLA-B*27': ['B27', np.nan, 'MHCI'],
            'HLA-B*57': ['B57', np.nan, 'MHCI'],
            'HLA-B*44:03:08': ['B44', 'B*44:03', 'MHCI'],
            'HLA-DRA*01': [None, np.nan, 'MHCII'],  # not in database
            'HLA-A*03:01': ['A3', 'A*03:01', 'MHCI'],
            'HLA-A*11:01': ['A11', 'A*11:01', 'MHCI'],
            'HLA-A*24:02': ['A24', 'A*24:02', 'MHCI'],
            'HLA-B*08:01': ['B8', 'B*08:01', 'MHCI'],
            'HLA-B*35:01': ['B35', 'B*35:01', 'MHCI'],
            'HLA-B*58': ['B58', np.nan, 'MHCI'],
            'HLA-B*42': ['B42', np.nan, 'MHCI'],
            'HLA-B*53': ['B53', np.nan, 'MHCI'],
            'HLA-B*27:05': ['B27', 'B*27:05', 'MHCI'],
            'HLA-B*35:08': ['B35', 'B*35:08', 'MHCI'],
            'HLA-DRA*01:01': [None, None, 'MHCII'],  # not in database
            'HLA-B*81:01': ['B81', 'B*81:01', 'MHCI'],
            'HLA-B*42:01': ['B42', 'B*42:01', 'MHCI'],
            'HLA-B*15': [None, np.nan, 'MHCI'],  # B62, 75, 72, 62, 70,71, 75, 77, 63 etc
            'HLA-B*51:01': ['B51', 'B*51:01', 'MHCI'],
            'HLA-B*07': ['B7', np.nan, 'MHCI'],
            'HLA-DQA1*01:02': [None, None, 'MHCII'],  # not in database
            'HLA-DPA1*02:01': [None, None, 'MHCII'],  # not in database
            'HLA-A*02:01:59': ['A2', 'A*02:01', 'MHCI'],
            'HLA-A*02:01:48': ['A2', 'A*02:01', 'MHCI'],
            'HLA-DRA*01:02:03': [None, None, 'MHCII'],  # not in database
            'HLA-B*08:01:29': ['B8', 'B*08:01', 'MHCI'],
            'HLA-A*02:256': ['A2', 'A*02:256', 'MHCI'],
            'HLA-B*35:08:01': ['B35', 'B*35:08', 'MHCI'],
            'HLA-E*01:01:01:03': [None, None, 'MHCI'],  # not in database
            'HLA-DRA*01:01:02': [None, None, 'MHCII'],  # not in database
            'HLA-B*35:42:01': ['B35', 'B*35:42', 'MHCI'],
            'HLA-B*57:06': ['B57', 'B*57:06', 'MHCI'],
            'HLA-B*44:05:01': ['B44', 'B*44:05:01', 'MHCI'],
            'HLA-A*02:01:98': ['A2', 'A*02:01', 'MHCI'],
            'HLA-A*24:02:84': ['A24', 'A*24:02', 'MHCI'],
            'HLA-B*27:05:31': ['B27', 'B*27:05', 'MHCI'],
            'HLA-DQA1*03:01:01': [None, None, 'MHCII'],  # not in database
            'HLA-B*51:193': ['B51', 'B*51:193', 'MHCI'],
            'HLA-DQA1*05:01:01:02': [None, None, 'MHCII'],  # not in database
            'HLA-A*02:01:110': ['A2', 'A*02:01', 'MHCI'],
            'HLA-B*35': ['B35', np.nan, 'MHCI'],
            'HLA-A*11': ['A11', np.nan, 'MHCI'],
            'HLA-B*12': [None, np.nan, 'MHCI'],  # not found. Could it be B44, B45 or B47?
            'HLA-A*68:01': ['A68', 'A*68:01', 'MHCI'],
            'HLA-B*18': ['B18', np.nan, 'MHCI'],
            'HLA class I': [np.nan, np.nan, 'MHCI'],
            'HLA-B*37:01': ['B37', 'B*37:01', 'MHCI'],
            'HLA-A*80:01': ['A80', 'A*80:01', 'MHCI'],
            'HLA-B27': ['B27', np.nan, 'MHCI'],
            'HLA-A11': ['A11', np.nan, 'MHCI'],
            'HLA-B57': ['B57', np.nan, 'MHCI'],
            'HLA-B*44:03': ['B44', 'B*44:03', 'MHCI'],
            'HLA-B8': ['B8', np.nan, 'MHCI'],
            'HLA-B*35:01, HLA-B*35:08': ['B35', np.nan, 'MHCI'],
            'HLA-A*30:02': ['A30', 'A*30:02', 'MHCI'],
            'HLA-B*38:01': ['B38', 'B*38:01', 'MHCI'],
            'HLA-B*44:02': ['B44', 'B*44:02', 'MHCI'],
            'HLA-B35': ['B35', np.nan, 'MHCI'],
            'HLA-B18': ['B18', np.nan, 'MHCI'],
            'HLA-B*15:01': ['B15', np.nan, 'MHCI'],
            'HLA-B*57:01, HLA-B*57:03': ['B57', np.nan, 'MHCI'],
            'HLA-C*14:02': [np.nan, np.nan, 'MHCI'],  # What is C14?
            'HLA-A*29:02': ['A29', 'A*29:02', 'MHCI'],
            'HLA-Cw3': [np.nan, np.nan, 'MHCI'],  # cannot find in the database, but does exist on wikipedia
            'HLA-A*02:01 K66A, E63Q mutant': ['A2', 'A*02:01', 'MHCI'],
            'HLA-B*18:01': ['B18', 'B18*:01', 'MHCI'],
            'HLA-B*40:01': ['B40', np.nan, 'MHCI'],
            'HLA-B*44:02, HLA-B*44:05': ['B44', np.nan, 'MHCI'],
            'HLA-B*44:02, HLA-B*44:03, HLA-B*44:05': ['B44', np.nan, 'MHCI'],
            'HLA-A*02:01, HLA-A*02:01 A150P mutant': ['A2', 'A*02:01', 'MHCI'],
            'HLA-A*30:01': ['A30', 'A*30:01', 'MHCI'],
            'HLA-A*02:01, HLA-A*02:01 K66A, E63Q mutant': ['A2', 'A*02:01', 'MHCI'],
            'HLA-B*35:08 Q65A, T69A, Q155A mutant': ['B35', 'B*35:08', 'MHCI'],
            'HLA-A*02:01 K66A mutant': ['A2', 'A*02:01', 'MHCI'],
            'HLA-C*06:02': [np.nan, np.nan, 'MHCI'],  # have to look into what / how to annotate the HLAC
            'HLA-C*07:01': [np.nan, np.nan, 'MHCI'],
            'HLA-C*05:01': [np.nan, np.nan, 'MHCI'],
            'HLA-B*35:08, HLA-B*35:08 Q65A, T69A, Q155A mutant': ['B35', 'B*35:08', 'MHCI'],
            'HLA-B*35:08 Q155A mutant': ['B35', 'B*35:08', 'MHCI'],
            'HLA-A*02:01, HLA-A*02:01 A69G mutant, HLA-A*02:01 E166A mutant, HLA-A*02:01 K66A mutant, HLA-A*02:01 R65A mutant, HLA-A*02:01 T163A mutant, HLA-A*02:01 W167A mutant': [
                'A2', 'A*02:01', 'MHCI'],
            'HLA-B*44:03, HLA-B*44:08': ['B44', 'B*44:08', 'MHCI'],
            'HLA-C*08:02': [np.nan, np.nan, 'MHCI']}

VDJDB_METHODS_DICT = {'tetramer-sort': 'multimer',
                      'tetramer-sort,antigen-loaded-targed': 'multimer',
                      'dextramer-sort': 'multimer',
                      'antigen-loaded-targets,cultured-T-cells': 'functional',
                      'dextramer-sort,cultured-T-cells': 'multimer',
                      'pentamer-sort': 'multimer',
                      'cultured-T-cells,beads,tetramer-sort': 'multimer',
                      'pelimer-sort': 'multimer',
                      'tetramer-magnetic-selection': 'multimer',
                      'multimer-sort': 'multimer',
                      'antigen-loaded-targets,dextramer-sort': 'multimer',
                      'cultured-T-cells,tetramer-sort': 'multimer',
                      'dextramer-sort,beads': 'multimer',
                      'antigen-expressing-targets': 'functional',
                      'MHC-peptide-beads,cultivated-t-cells': 'multimer',
                      'streptamer-sort': 'multimer',
                      'tetramer-sort,cultured-T-cells': 'multimer',
                      'antigen-loaded-targets': 'functional',
                      'cla': 'functional',
                      'peptide-restimulation,tetramer-sort': 'multimer',
                      'tetramer-sort,beads': 'multimer',
                      'tetramer-sort,limiting-dilution-cloning': 'multimer',
                      'limiting-dilution-cloning,tetramer-sort': 'multimer',
                      'antigen-loaded-targets, pentamer-sort': 'functional',
                      'cd8null-tetramer-sort': 'multimer',
                      'cultured-T-cells': 'functional',
                      'antigen-loaded-targets, IFNg capture assay': 'functional',
                      'antigen loaded targets, IFNg capture assay': 'functional',
                      'antigen-loaded-targets, tetramer-sort': 'multimer',
                      'nan': np.nan,
                      'cultured-T-cells,antigen-expressing-targets': 'functional',
                      'limiting-diffusion-cloning,antigen-expressing-targets': 'functional',
                      'antigen-loaded-target, CTL culture': 'functional',
                      'Tetramer-sort': 'multimer',
                      'limiting-dilution-cloning': np.nan,
                      'tetramer-sort,antigen-expressing cells,cultured-T-cells': 'multimer',
                      'antigen-loaded-target,limiting-dilution-cloning': 'functional',
                      'beads,tetramer-sort': 'multimer',
                      'limited-dilution-cloning,antigen-expressing-targets': 'functional',
                      'antigen-loaded-targets,pentamer-sort': 'multimer',
                      'CTL clone': np.nan,
                      'antigen-loaded-targets,cloning-by-limiting-dilution': 'functional'
                      }

VDJDB_SINGLE_CELL_DICT = {'nan': np.nan,
                          'yes': 'Yes',
                          'tetramer-sort': 'No'}

MCPAS_ASSAY_DICT = {'1.0': 'multimer',
                    '2.0': 'functional',
                    '2.1': 'functional',
                    '2.2': 'functional',
                    '2.3': 'functional',
                    '2.4': 'functional',
                    '2.5': 'functional',
                    '3.0': 'remove',
                    '4.0': 'ngs',
                    'nan': np.nan}

EPITOPE_MHC_DICT_10x = {'nan': np.nan,
                        'VTEHDTLLY': 'HLA-A*01:01',
                        "CYTWNQMNL": 'HLA-A*24:02',
                        "AYAQKIFKI": "HLA-A*24:02",
                        "QYDPVAALF": "HLA-A*24:02",
                        "KLGGALQAK": "HLA-A*03:01",
                        "RLRAEAQVK": "HLA-A*03:01",
                        "RIAAWMATY": "HLA-A*03:01",
                        "IVTDFSVIK": "HLA-A*11:01",
                        "AVFDRKSDAK": "HLA-A*11:01",
                        "IPSINVHHY": "HLA-B*35:01",
                        "QPRAPIRPI": "HLA-B*07:02",
                        "TPRVTGGGAM": "HLA-B*07:02",
                        "RPPIFIRRL": "HLA-B*07:02",
                        "RPHERNGFTVL": "HLA-B*07:02",
                        "RAKFKQLL": "HLA-B*08:01",
                        "ELRRKMMYM": "HLA-B*08:01",
                        "FLRGRAYGL": "HLA-B*08:01"}

# The TRAJ allele can be imputed from the CDR3A region in 2 cases
TRAJ_IMPUTATION_DICT = {'TRAJ37': {'GNTGKLIF': ['TRAJ37*01', 'GSGNTGKLIFGQGTTLQVKP'],  # GSGNTGKLI is allele *01
                                   'SNTGKLIF': ['TRAJ37*02', 'GSSNTGKLIFGQGTTLQVKP'],  # XXSXXXXXX is allele *02
                                   },
                        'TRAJ24': {'FEF': ['TRAJ24*01', 'TTDSWGKFEFGAGTQVVVTP'],  # TTDSWGKFE is allele *01
                                   'LQF': ['TRAJ24*02', 'TTDSWGKLQFGAGTQVVVTP'],  # XXXXXXXLQ is allele *02
                                   'FQF': ['TRAJ24*03', 'TTDSWGKFQFGAGTQVVVTP'],  # XXXXXXXFQ is allele *03
                                   }}

# these are most options of VDJ allele notations that are found in the database created by parser.py
# TODO: MAKE IT ABLE TO ADD NEW ONES: using self.TRAV etc
TRAV = ['nan', 'TRAV38-2DV8', 'TRAV29DV5', 'TRAV13-1', 'TRAV1-1', 'TRAV5', 'TRAV14DV4', 'TRAV1-2', 'TRAV8-2', 'TRAV24',
        'TRAV41', 'TRAV19', 'TRAV8-1', 'TRAV27', 'TRAV17', 'TRAV35', 'TRAV22', 'TRAV12-3', 'TRAV13-2', 'TRAV8-3',
        'TRAV10', 'TRAV39', 'TRAV9-2', 'TRAV2', 'TRAV21', 'TRAV3', 'TRAV12-2', 'TRAV16', 'TRAV26-1', 'TRAV38-1',
        'TRAV6', 'TRAV12-1', 'TRAV4', 'TRAV25', 'TRAV36DV7', 'TRAV8-6', 'TRAV8-4', 'TRAV26-2', 'TRAV30', 'TRAV34',
        'TRAV23DV6', 'TRAV40', 'TRAV20', 'TRAV18', 'TRAV9-1', 'TRAV9', 'TRAV14-1', 'TRAV1', 'TRAV8', 'TRAV1-4',
        'TRAV28', 'TRAV2-1', 'TRAV26', 'TRAV15', 'TRAV2-3', 'TRAV4-1', 'TRAV15-1', 'TRAV23', 'TRAV2-2', 'TRAV12',
        'TRAV14', 'TRAV12-4', 'TRAV38-2', 'TRAV29', 'TRAV31', 'TRAV11', 'TRAV1-01', 'TRAV251', 'TRAV21-1', 'TRAV21-2',
        'TRAV3-1', 'TRAV35-1', 'TRAV29/DV5', 'TRAV14/DV4', 'TRAV23/DV6', 'TRAV12-2, TRAV21', 'TRAV2-01', 'TRAV35:01',
        'TRAV12-3:01', 'TRAV25:01', 'TRAV3-01', 'TRAV10-01', 'TRAV1-2:01', 'TRAV14/DV4:01', 'TRAV29/DV5:01',
        'TRAV1-1:01', 'TRAV6-01', 'TRAV38-2/DV8:01', 'TRAV13-2:01', 'TRAV5-01', 'TRAV23/DV6:01', 'TRAV9-2:01',
        'TRAV41:01', 'TRAV38-1:01', 'TRAV14/DV4:02', 'TRAV8-1:01', 'TRAV12-2:01', 'TRAV22:01', 'TRAV8-6:02',
        'TRAV8-6:01', 'TRAV8-3:01', 'TRAV13-1:01', 'TRAV41-01', 'TRAV8-2:01', 'TRAV21-01', 'TRAV26-2:01', 'TRAV8-4:01',
        'TRAV4:01', 'TRAV16:01', 'TRAV36/DV7:01', 'TRAV12-1:01', 'TRAV26-1:01', 'TRAV16-01', 'TRAV38-2/DV8',
        'TRAV24:01', 'TRAV5:01', 'TRAV3:01', 'TRAV20:02', 'TRAV6:02', 'TRAV19:01', 'TRDV1', 'TRAV8-2/8-4',
        'TRAV36/DV7', 'TRAV7', 'TRAV17:01', 'TRAV29/DV5:02', 'TRAV29/DV5:05', 'TRAV6:03', 'TRAV27:01', 'TRAV26-1:03',
        'TRAV2:02', 'TRAV8-4:03', 'TRAV78', 'TRAV36/DV7:02', 'TRAV24:01/TRAV24:02', 'TRAV24:01,', 'TRAV24:02',
        'TRAV3 F', 'TRAV1-2\xa0', 'TRAV22 F', 'TRAV09-2', 'TRAV29/DV5 F', 'TRAV8-1 F', 'TRAV25 F', 'TRAV38-2/DV8 F',
        'TRAV13-2/13-2:02', 'TRAV41 F', 'TRAV121', 'TRAV-2', 'TRAV08-3', 'TRAV08-2', 'TRAV08-4', 'TRAV01-2',
        'TRAV08-1', 'TRAV38-2/DV8*01', 'TRAV26-1*01', 'TRAV26-2*01', 'TRAV24*01', 'TRAV39*01', 'TRAV4*01',
        'TRAV12-2*01', 'TRAV17*01', 'TRAV29/DV5*01', 'TRAV38-1*01', 'TRAV3*01 F', 'TRAV41*01', 'TRAV8-6*01',
        'TRAV9-2*01', 'TRAV1-2*01', 'TRAV14/DV4*01', 'TRAV13-2*01', 'TRAV16*01', 'TRAV53', 'TRAV56', 'TRAV33',
        'TRAV13', 'TRAV32', 'TRAV36', 'TRAV37', 'TRAV38', 'TRAV42', 'TRAV43', 'TRAV44', 'TRAV45', 'TRAV46', 'TRAV47',
        'TRAV48', 'TRAV49', 'TRAV50', 'TRAV51', 'TRAV52', 'TRAV54', 'TRAV55', 'TRAV14DV4.', 'TRAV29DV5.', 'TRAV23DV6.',
        'TRA1-2', 'TRAV38-2/DV8-1', 'TRAV22-1', '36/DV7*01', 'TRAV8-3*01', 'TRAV20-1', 'TRAV24-1', 'TRAV14/DV4*03',
        'TRAV10-2', 'TRAV6-1', 'mTRDV2-2', 'mTRAV14D-1', 'TRAV29-1', 'TRAV12-02', 'TRAV25-01*01', 'TRAV19-01*01',
        'TRAV39-01*01', 'TRAV20-01', 'TRAV13-01', 'TRAV14-01', 'TRAV29-01', 'TRAV12-01', 'TRAV27-01*01', 'TRDV02-01',
        'TRAV24-01*01', 'TRAV26-01', 'TRAV20*01', 'TRAV13-1*01', 'TRAV36/DV7*01', 'TRAV19*01', 'TRAV8-2*01',
        'TRAV27*01', 'TRAV35*01', 'TRAV23/DV6*01', 'TRAV18*01', 'TRAV22*01', 'TRAV21*01', 'TRAV8-4*01', 'TRAV5*01',
        'TRAV12-1*01', 'TRAV12-3*01', 'TRAV3*01', 'TRAV34*01', 'TRAV8-1*01', 'TRAV1-1*01', 'TRAV2*01', 'TRAV30*01',
        'TRAV6*01', 'TRAV10*01', 'TRAV25*01', 'TRAV40*01', 'TRAV9-1*01', 'TRAV7*01', 'TRAV14/DV4*02', 'TRAV8-6*02']
TRAJ = ['nan', 'TRAJ34', 'TRAJ5', 'TRAJ48', 'TRAJ53', 'TRAJ6', 'TRAJ37', 'TRAJ49', 'TRAJ9', 'TRAJ23', 'TRAJ31',
        'TRAJ36', 'TRAJ54', 'TRAJ40', 'TRAJ42', 'TRAJ26', 'TRAJ28', 'TRAJ8', 'TRAJ32', 'TRAJ13', 'TRAJ11', 'TRAJ58',
        'TRAJ10', 'TRAJ4', 'TRAJ33', 'TRAJ17', 'TRAJ21', 'TRAJ15', 'TRAJ52', 'TRAJ24', 'TRAJ22', 'TRAJ43', 'TRAJ35',
        'TRAJ18', 'TRAJ12', 'TRAJ30', 'TRAJ20', 'TRAJ7', 'TRAJ50', 'TRAJ39', 'TRAJ56', 'TRAJ3', 'TRAJ57', 'TRAJ47',
        'TRAJ44', 'TRAJ41', 'TRAJ16', 'TRAJ45', 'TRAJ27', 'TRAJ29', 'TRAJ46', 'TRAJ38', 'TRAJ14', 'TRAJ25', 'TRAJ1',
        'TRAJ5-1', 'TRAJ2', 'TRAJ1-3', 'TRAJ9-1', 'TRAJ16-5', 'TRAJ1-8', 'TRAJ10-1', 'TRAJ3-2', 'TRAJ9-14', 'TRAJ14-1',
        'TRAJ16-1', 'TRAJ37-2', 'TRAJ4-01', 'TRAJ5-01', 'TRAJ3-01', 'TRAJ1-01', 'TRAJ24:02', 'TRAJ2-01', 'TRAJ9-01',
        'TRAJ6-01', 'TRAJ8-01', 'TRAJ13:02', 'TRAJ49:01', 'TRAJ20:01', 'TRAJ42:01', 'TRAJ24:01', 'TRAJ21:01',
        'TRAJ7:01', 'TRAJ39:01', 'TRAJ36:01', 'TRAJ54:01', 'TRAJ16:01', 'TRAJ58:01', 'TRAJ19', 'TRAJ53:01',
        'TRAJ17:01', 'TRAJ', 'TRAJ53:02', 'TRAJ53:05', 'TRAJ43:01', 'TRAJ2-1', 'TRAJ41:01', 'TRAJ3-1', 'TRAJ34:01',
        'TRAJ31:01', 'TRAJ44:01', 'TRAJ38:01', 'TRAJ37:01', 'TRAJ48:01', 'TRAJ9:01', 'TRAJ18:01', 'TRAJ26:01',
        'TRAJ52:01', 'TRAJ29:01', 'TRAJ56:01', 'TRAJ45:01', 'TRAJ3:01', 'TRAJ11:01', 'TRAJ27:01', 'TRAJ33:01',
        'TRAJ37:02', 'TRAJ23:01', 'TRAJF', 'TRAJ36/DV7:02', 'TRAJ57:01', 'TRAJ6:01', 'TRAJ32:02', 'TRAJ22:01',
        'TRAJ10:01', 'TRAJ32:01', 'TRAJ5:01', 'TRAJ2-2', 'TRAJ3:58 0', 'TRAJ6: 56 0', 'TRAJ52 67 0', 'TRAJ26:51',
        'TRAJ9:53', 'TRAJ4 ', 'TRAJ49 ', 'TRAJ36 ', 'TRAJ20 ', 'TRAJ5 58', 'TRAJ5 58 ', 'TRAJ28 ', 'TRAJ13 ',
        'TRAJ21 ', 'TRAJ47-1', 'TRAJ20-1', 'TRAJ24-1', 'TRAJ32-1', 'TRAJ44-1', 'TRAJ4-1', 'TRAJ37-1', 'TRAJ53-1',
        'TRAJ34-1', 'TRAJ21-1', 'TRAJ33-1', 'TRAJ45-1', 'TRAJ54-1', 'TRAJ17-1', 'TRAJ42-1', 'TRAJ28-1', 'TRAJ40-1',
        'TRAJ22-5', 'TRAJ57/58', 'TRAJ21/', 'TRAJ65', 'TRAJ48/', 'TRAJ34/58', 'TRAJ1-16', 'TRAJ2-5 ', 'TRAJ38/57 ',
        'TRAJ3 49 ', 'TRAJ61', 'TRAJ34/43', 'TRAJ42/55', 'TRAJ22-6', 'TRAJ42/', 'TRAJ33/51', 'TRAJ52/64', 'TRAJ49/54',
        'TRAJ31-1', 'TRAJ45 61 0', 'TRAJ5 60 0', 'TRAJ53 63 0', 'TRAJ40 59 0', 'TRAJ23 60 0', 'TRAJ49 48 0',
        'TRAJ53 58 0', 'TRAJ49 54 0', 'TRAJ24 53 0', 'TRAJ56 59 0', 'TRAJ45 66 0', 'TRAJ10 59 0', 'TRAJ44 58 0',
        'TRAJ54 56 0', 'TRAJ4 62 0', 'TRAJ29 60 0', 'TRAJ54 54 0', 'TRAJ5 56 0', 'TRAJ15 60 0', 'TRAJ49 50 0',
        'TRAJ42 56 0', 'TRAJ32 57 0', 'TRAJ39 57 0', 'TRAJ58 60 0', 'TRAJ37 59 0', 'TRAJ37 60 0', 'TRAJ31 52 0',
        'TRAJ39 50 0', 'TRAJ49 53 0', 'TRAJ12 50 0', 'TRAJ10 60 0', 'TRAJ52 66 0', 'TRAJ10 56 0', 'TRAJ39 58 0',
        'TRAJ27 56 0', 'TRAJ50 56 0', 'TRAJ27 55 0', 'TRAJ16 55 1', 'TRAJ49 55 0', 'TRAJ49 52 0', 'TRAJ16 57 1',
        'TRAJ22 52 0', 'TRAJ12 55 0', 'TRAJ4 58 0', 'TRAJ6 59 0', 'TRAJ47 57 0', 'TRAJ20 54 0', 'TRAJ9 56 0',
        'TRAJ17 62 0', 'TRAJ29 56 0', 'TRAJ6 58 0', 'TRAJ28 60 0', 'TRAJ32 59 0', 'TRAJ11 60 0', 'TRAJ54 53 0',
        'TRAJ54 60 0', 'TRAJ20 53 0', 'TRAJ28 56 0', 'TRAJ20 39 0', 'TRAJ4 57 0', 'TRAJ5 59 1', 'TRAJ43 46 0',
        'TRAJ8 51 0', 'TRAJ30 57 0', 'TRAJ9 57 0', 'TRAJ23 45 0', 'TRAJ49 56 0', 'TRAJ17 60 0', 'TRAJ33 57 0',
        'TRAJ15 55 0', 'TRAJ8 56 0', 'TRAJ6 57 0', 'TRAJ18 62 0', 'TRAJ20 50 0', 'TRAJ42 65 0', 'TRAJ53 60 0',
        'TRAJ52 69 0', 'TRAJ22 51 0', 'TRAJ11 56 0', 'TRAJ7 56 0', 'TRAJ17 55 0', 'TRAJ45 56 0', 'TRAJ8 54 0',
        'TRAJ15 57 0', 'TRAJ13 58 0', 'TRAJ41 55 0', 'TRAJ13 59 0', 'TRAJ7 52 0', 'TRAJ27 58 0', 'TRAJ39 56 0',
        'TRAJ21 54 0', 'TRAJ7 55 0', 'TRAJ18 61 0', 'TRAJ8 53 0', 'TRAJ40 54 0', 'TRAJ52 65 0', 'TRAJ58 59 0',
        'TRAJ40 48 0', 'TRAJ9 55 0', 'TRAJ34 44 0', 'TRAJ43 54 1', 'TRAJ31 46 0', 'TRAJ22 58 0', 'TRAJ29 57 0',
        'TRAJ22 59 0', 'TRAJ15 58 0', 'TRAJ53 61 0', 'TRAJ22 61 0', 'TRAJ37 58 0', 'TRAJ23 58 0', 'TRAJ48 56 0',
        'TRAJ33 52 0', 'TRAJ26 56 0', 'TRAJ34 50 0', 'TRAJ54 57 0', 'TRAJ21 55 0', 'TRAJ57 58 0', 'TRAJ45 58 0',
        'TRAJ3 56 0', 'TRAJ57 57 0', 'TRAJ35 56 0', 'TRAJ23 59 0', 'TRAJ38 61 0', 'TRAJ57 61 0', 'TRAJ28 63 0',
        'TRAJ27 57 0', 'TRAJ48 46 0', 'TRAJ50 60 0', 'TRAJ54 55 0', 'TRAJ40 57 0', 'TRAJ39 52 0', 'TRAJ58 58 0',
        'TRAJ17 53 0', 'TRAJ8 60 0', 'TRAJ5 58 1', 'TRAJ42 58 0', 'TRAJ11 53 0', 'TRAJ29 59 0', 'TRAJ27 59 0',
        'TRAJ37 61 0', 'TRAJ38 62 0', 'TRAJ43 47 0', 'TRAJ56 58 0', 'TRAJ56 56 0', 'TRAJ11 59 0', 'TRAJ29 54 0',
        'TRAJ20 52 0', 'TRAJ44 58 1', 'TRAJ39 59 0', 'CATSESSGQTYEQYF', 'TRAJ43-01*01', 'TRAJ17-01*01', 'TRAJ27-01*01',
        'TRAJ42-01*01', 'TRAJ39-01*01', 'TRAJ53-01*01', 'TRAJ52-01*01', 'TRAJ45-01*01', 'TRAJ50-01*01', 'TRAJ26-01*01',
        'TRAJ20-01*01', 'TRDJ01-01*01', 'TRAJ47-01*01', 'TRAJ49-01*01', 'TRAJ40-01*01', 'TRAJ22-01*01', 'TRAJ41-01*01',
        'TRAJ37-01', 'TRAJ18-01*01', 'TRAJ30-01*01', 'TRAJ29-01*01', 'TRAJ54-01*01', 'TRAJ06-01*01', 'TRAJ34-01*01',
        'TRDJ1', 'TRAJ39*01', 'TRAJ10*01', 'TRAJ42*01', 'TRAJ32*01', 'TRAJ13*01', 'TRAJ8*01', 'TRAJ47*01', 'TRAJ29*01',
        'TRAJ17*01', 'TRAJ38*01', 'TRAJ44*01', 'TRAJ20*01', 'TRAJ49*01', 'TRAJ26*01', 'TRAJ34*01', 'TRAJ40*01',
        'TRAJ33*01', 'TRAJ23*01', 'TRAJ57*01', 'TRAJ43*01', 'TRAJ6*01', 'TRAJ30*01', 'TRAJ24*01', 'TRAJ48*01',
        'TRAJ31*01', 'TRAJ45*01', 'TRAJ28*01', 'TRAJ21*01', 'TRAJ11*01', 'TRAJ15*01', 'TRAJ3*01', 'TRAJ53*01',
        'TRAJ7*01', 'TRAJ9*01', 'TRAJ22*01', 'TRAJ37*01', 'TRAJ54*01', 'TRAJ52*01', 'TRAJ50*01', 'TRAJ12*01',
        'TRAJ27*01', 'TRAJ41*01', 'TRAJ4*01', 'TRAJ18*01', 'TRAJ36*01', 'TRAJ56*01', 'TRAJ5*01', 'TRAJ46*01',
        'TRAJ13*02', 'TRAJ24*02', 'TRAJ14*01']
TRBV = ['nan', 'TRBV2', 'TRBV5-1', 'TRBV4-3', 'TRBV15', 'TRBV3-1', 'TRBV7-9', 'TRBV12-4', 'TRBV20-1', 'TRBV7-7',
        'TRBV13', 'TRBV18', 'TRBV7-2', 'TRBV19', 'TRBV11-2', 'TRBV30', 'TRBV4-2', 'TRBV6-3', 'TRBV9', 'TRBV10-2',
        'TRBV6-5', 'TRBV5-6', 'TRBV6-4', 'TRBV27', 'TRBV12-3', 'TRBV28', 'TRBV11-3', 'TRBV6-6', 'TRBV10-3', 'TRBV14',
        'TRBV4-1', 'TRBV24-1', 'TRBV6-1', 'TRBV5-4', 'TRBV25-1', 'TRBV7-3', 'TRBV29-1', 'TRBV7-8', 'TRBV5-5',
        'TRBV7-6', 'TRBV12-5', 'TRBV16', 'TRBV11-1', 'TRBV5-8', 'TRBV7-4', 'TRBV21-1', 'TRBV6-7', 'TRBV10-1',
        'TRBV5-7', 'TRBV23-1', 'TRBV6-8', 'TRBV6-9', 'TRBV5-3', 'TRBV22', 'TRBV6', 'TRBV7', 'TRBV13-2', 'TRBV2-01',
        'TRBV4', 'TRBV12', 'TRBV11-01', 'TRBV7-08', 'TRBV25-01', 'TRBV20-01', 'TRBV21', 'TRBV6-05', 'TRBV5',
        'TRBV6-03', 'TRBV8', 'TRBV3', 'TRBV23', 'TRBV17', 'TRBV29-01', 'TRBV12-01', 'TRBV3-01', 'TRBV6-07', 'TRBV5-01',
        'TRBV9-01', 'TRBV8-2', 'TRBV4-01', 'TRBV6-01', 'TRBV12-03', 'TRBV21-03', 'TRBV6-04', 'TRBV18-01', 'TRBV7-03',
        'TRBV16-01', 'TRBV1-01', 'TRBV11', 'TRBV1', 'TRBV13-3', 'TRBV30-05', 'TRBV12-04', 'TRBV28-01', 'TRBV9-03',
        'TRBV7-09', 'TRBV11-02', 'TRBV4-02', 'TRBV11-03', 'TRBV5-05', 'TRBV5-04', 'TRBV5-08', 'TRBV4-03', 'TRBV16-03',
        'TRBV7-02', 'TRBV24-01', 'TRBV20', 'TRBV6-02', 'TRBV10-02', 'TRBV7-06', 'TRBV6-02/3', 'TRBV5-06', 'TRBV6-06',
        'TRBV24', 'TRBV10-03', 'TRBV25', 'TRBV10', 'TRBV29', 'TRBV5-07', 'TRBV10-01', 'TRBV3-02', 'TRBV1-02',
        'TRBV13-06', 'TRBV6-09', 'TRBV12-05', 'TRBV21-01', 'TRBV7-04', 'TRBV7-07', 'TRBV1-05', 'TRBV6-02/TRBV6-03',
        'TRBV12-02', 'TRBV6-02/6-03', 'TRBV23-01', 'TRBV19-01', 'TRBV2-03', 'TRBV12-03/4 ', 'TRBV12-03/4', 'TRBV14-01',
        'TRBV2-05', 'TRBV27-01', 'TRBV13-01', 'TRBV-012', 'TRBV15-01', 'TRBV30-01', 'TRBV15-02', 'TRBV7-08/TRBV7-02',
        'TRBV12-03/TRBV12-04', 'TRBV12-03,TRBV12-04', 'TRBV7-08, TRBV7-02', 'TRBV6-01,TRBV6-06',
        'TRBV6-01,TRBV6-05,TRBV6-06', 'TRBV6-02,TRBV6-03', 'TRBV12-04 ', 'TRBV6-03 ', ' TRBV6-02', ' TRBV6-01',
        ' TRBV6-06', 'TRBV19:02', 'TRBV7-01', 'TRBV9-02', 'TRBV6-4:02', 'TRBV6-2', 'TRBV27:01', 'TRBV20-1:01',
        'TRBV25-1:01', 'TRBV7-3:01', 'TRBV12-4:01', 'TRBV2-1', 'TRBV7-1', 'TRBV27-1', 'TRBV13-1', 'TRBV28-1',
        'TRBV5-1:01', 'TRBV18:01', 'TRBV20-1:03-07', 'TRBV6-2/3', 'TRBV201-1', 'TRBV12-3/4', 'TRBV12-3/12-4',
        'TRBV3-2', 'TRBV6-2/6-3', 'TRBV3-1:01', 'TRBV1-1', 'TRBV24-1:01', 'TRBV6-2:01', 'TRBV19:01', 'TRBV28:01',
        'TRBV6-5:01', 'TRBV29/DV5:01', 'TRBV25:01', 'TRBV5-8:01', 'TRBV6-1:01', 'TRBV4-2:01', 'TRBV7-8:01',
        'TRBV12-3:01', 'TRBV26-1:03', 'TRBV29-1:01', 'TRBV5-6:01', 'TRBV6-6:01', 'TRBV4:01', 'TRBV26-1:01',
        'TRBV17:01', 'TRBV8-4:05', 'TRBV1-5:01', 'TRBV2-2:01', 'TRBV2-1:01', 'TRBV1-2:01', 'TRBV1-1:01', 'TRBV1-4:01',
        'TRBV2-6:01', 'TRBV3:01', 'TRBV4-1:01', 'TRBV5-5:01', 'TRBV8-1:01', 'TRBV8-2:01', 'TRBV7-9:01', 'TRBV2-7:01',
        'TRBV14:01', 'TRBV', 'TRBV4-3:01', 'TRBV9:01', 'TRBV5-4:01', 'TRBV19:03', 'TRBV19:04', 'TRBV19:06',
        'TRBV22:01', 'TRBV26-2:01', 'TRBV23/DV6:01', 'TRBV7-8:03', 'TRBV19:05', 'TRBV7-2:01', 'TRBV2-1/TRBV2-3',
        'TRBV2-1,', 'TRBV2-1/TRBV2-2/TRBV2-3', 'TRBV2-2', 'TRBV2-1/TRBV2-2', 'TRBV6-2,6-3', 'TRBV12-3,12-4',
        'TRBV7-9*01', 'TRBV3-1/3-2', 'TRBV6-1/6-5/6-6', 'TRBV12-2', 'TRBV-2', 'TRBV-3', 'TRBV-6', 'TRBV-5', 'TRBV-4',
        'TRBV-7', 'TRBV-9', 'TRBV8-3', 'TRBV18-1', 'TRBV9-1', 'TRBV15-1*01', 'TRBV22-1*01', 'TRBV11-1*01',
        'TRBV7-7*01', 'TRBV7-1*01', 'TRBV28-1*01', 'TRBV6-9*01', 'TRBV27-1*01', 'TRBV21-1*01', 'TRBV26-1*01',
        'TRBV13-1*01', 'TRBV5-6*01', 'TRBV5-1*01', 'TRBV5-5*01', 'TRBV5-2*01', 'TRBV18-1*01', 'TRBV1-1*01',
        'TRBV5-4*01', 'TRBV12-1*01', 'TRBV23-1*01', 'TRBV6-1*01', 'TRBV12-5*01', 'TRBV10-3*01', 'TRBV11-2*02',
        'TRBV7-2*01', 'TRBV8-2*01', 'TRBV10-2*01', 'TRBV6-8*01', 'TRBV7-3*01', 'TRBV17-1*01', 'TRBV6-7*01',
        'TRBV14-1*01', 'TRBV7-5*01', 'TRBV25-1*01', 'TRBV4-2*01', 'TRBV5-7*01', 'TRBV2-1*01', 'TRBV16-1', 'TRBV19-1',
        'TRBV4-3*01', 'TRBV9-02*01', 'TRBV11-3*01', 'TRBV12-2*01', 'TRBV7-8*01', 'TRBV7-6*01', 'TRBV30-1*01',
        'TRBV7-4*01', 'TRBV4-1*01', 'TRBV29-1*01', 'TRBV7-2*02', 'TRBV20-1*01', 'TRBV1-2', 'TRBV1-3', 'TRBV7-08*01',
        'TRBJ2-7', 'TRBV05', 'TRBV05-06*01', 'TRBV04-03*01', 'TRBV15-01*01', 'TRBV11-02*02', 'TRBV07-09', 'TRBV06',
        'TRBV05-04*01', 'TRBV07-02*01', 'TRBV04-01*01', 'TRBV10-03*01', 'TRBV03', 'TRBV28-01*01', 'TRBV05-01*01',
        'TRBV30-01*01', 'TRBV09-01', 'TRBV09', 'TRBV07', 'TRBV02', 'TRBV04', 'TRBV19*01', 'TRBV7-2*03', 'TRBV28*01',
        'TRBV27*01', 'TRBV6-3*01', 'TRBV6-5*01', 'TRBV6-6*01', 'TRBV6-2*01', 'TRBV30*01', 'TRBV9*01', 'TRBV11-2*01',
        'TRBV12-3*01', 'TRBV2*01', 'TRBV24-1*01', 'TRBV18*01', 'TRBV15*01', 'TRBV3-1*01', 'TRBV12-4*01', 'TRBV13*01',
        'TRBV14*01', 'TRBV16*01', 'TRBV10-1*01', 'TRBV6-4*01', 'TRBV5-8*01', 'TRBV13S1', 'TRBV20S1', 'TRBV6-6*02',
        'TRBV19*02', 'TRBV9*02', 'TRBV6-4*02', 'TRBV7-9*03']
TRBD = ['nan', 'None', 'TRBD2', 'TRBD1', 'TRBD2-1', 'TRBD1-1', 'TRBD2-2', 'TRBD1-2', 'TRBD1, TRBD2', 'TRBD1:01',
        'TRBD2:01', 'TRBDTRBJ2-6', 'TRBD2:02', 'TRBD', 'TRBD2-2:01', 'TRBD2-1*02', 'TRBD1-1*01', 'TRBD2-1*01',
        'TRBD1*01', 'TRBD2*02', 'TRBD2*01', 'TRBD01-01*01', 'TRBD02-01*02', 'unknown', 'TRBD02-01', 'TRBD02-01*01',
        'na', 'TRBD1,TRBD2']
TRBJ = ['nan', 'TRBJ2-1', 'TRBJ2-3', 'TRBJ2-7', 'TRBJ1-3', 'TRBJ1-2', 'TRBJ2-2', 'TRBJ1-1', 'TRBJ2-6', 'TRBJ2-5',
        'TRBJ1-6', 'TRBJ1-5', 'TRBJ1-4', 'TRBJ2-4', 'TRBJ3-1', 'TRBJ2', 'TRBJ1-7', 'TRBJ5-6', 'TRBJ2-7:01',
        'TRBJ1-2:01', 'TRBJ2-1:01', 'TRBJ2-2:01', 'TRBJ2-5:01', 'TRBJ2-3:01', 'TRBJ1-5:01', 'TRBJ1-1:01', 'TRBJ2-4:01',
        'TRBJ1-3:01', 'TRBJ1-4:01', 'TRBJ2-6:01', 'TRAJ30', 'TRAJ26', 'TRBJ1-6:02', 'TRBJ53:01', 'TRBJTRBJ2-7',
        'TRBJ54:01', 'TRBJ3-1:01', 'TRBJ42:01', 'TRBJ39:01', 'TRBJ45:01', 'TRBJ2:01', 'TRBJ1:01', 'TRBJ31:01',
        'TRBJ43:01', 'TRBJ33:01', 'TRBJ', 'TRBJ27:01', 'TRBJ37:02', 'TRBJ8:01', 'TRBJ52:01', 'TRBJ9:01', 'TRBJ36:01',
        'TRBJ2-7:02', 'TRBJ1-246', 'TRBJ1-203', 'TRBJ1-59', 'TRBJ1-258', 'TRBJ1-164', 'TRBJ1-304', 'TRBJ1-89',
        'TRBJ1-229', 'TRBJ1-223', 'TRBJ1-305', 'TRBJ1-90', 'TRBJ1-91', 'TRBJ1-92', 'TRBJ1-93', 'TRBJ1-94', 'TRBJ1-95',
        'TRBJ1-96', 'TRBJ1-97', 'TRBJ1-98', 'TRBJ1-99', 'TRBJ1-100', 'TRBJ1-101', 'TRBJ1-102', 'TRBJ1-103',
        'TRBJ1-104', 'TRBJ1-105', 'TRBJ1-106', 'TRBJ1-107', 'TRBJ1-108', 'TRBJ1-109', 'TRBJ1-110', 'TRBJ1-111',
        'TRBJ1-112', 'TRBJ1-113', 'TRBJ1-114', 'TRBJ1-115', 'TRBJ1-130', 'TRBJ1-314', 'TRBJ1-205', 'TRBJ1-234',
        'TRBJ1-292', 'TRBJ1-116', 'TRBJ1-117', 'TRBJ1-49', 'TRBJ1-8', 'TRBJ1-9', 'TRBJ1-118', 'TRBJ1-220', 'TRBJ1-238',
        'TRBJ1-241', 'TRBJ1-274', 'TRBJ1-306', 'TRBJ1-119', 'TRBJ1-44', 'TRBJ1-216', 'TRBJ1-302', 'TRBJ1-163',
        'TRBJ1-12', 'TRBJ1-26', 'TRBJ1-31', 'TRBJ1-40', 'TRBJ1-41', 'TRBJ1-63', 'TRBJ1-120', 'TRBJ1-121', 'TRBJ1-219',
        'TRBJ1-259', 'TRBJ1-122', 'TRBJ1-123', 'TRBJ1-235', 'TRBJ1-32', 'TRBJ1-124', 'TRBJ1-14', 'TRBJ1-33',
        'TRBJ1-125', 'TRBJ1-126', 'TRBJ1-195', 'TRBJ1-196', 'TRBJ1-197', 'TRBJ1-198', 'TRBJ1-232', 'TRBJ1-233',
        'TRBJ1-287', 'TRBJ1-268', 'TRBJ1-87', 'TRBJ1-127', 'TRBJ1-208', 'TRBJ1-129', 'TRBJ1-168', 'TRBJ1-209',
        'TRBJ1-278', 'TRBJ1-279', 'TRBJ1-131', 'TRBJ1-132', 'TRBJ1-133', 'TRBJ1-20', 'TRBJ1-21', 'TRBJ1-22',
        'TRBJ1-23', 'TRBJ1-134', 'TRBJ1-135', 'TRBJ2-6*01', 'TRBJ2-3*01', 'TRBJ2-1*01', 'TRBJ1-5*01', 'TRBJ2-2*01',
        'TRBJ1-6*01', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-1*01', 'TRBJ2-7*01', 'TRBJ2-4*01', 'TRBJ1-4*01', 'TRBJ2-5*01',
        'TRBJ1-6*02', 'Donor 13', 'Negative', 'TRBJ01-02*01', 'TRBJ01-04*01', 'TRBJ02-01*01', 'TRBJ02-05*01',
        'TRBJ02-07*01', 'TRBJ01-03*01', 'TRBJ01-01*01', 'TRBJ01-06*01', 'TRBJ02-02*01', 'TRBJ01-05*01', 'TRBJ02-03*01',
        'TRBJ02-06*01', 'TRBJ1S2']
