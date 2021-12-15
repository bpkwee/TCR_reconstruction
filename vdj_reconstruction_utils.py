"""
This script is imputes the amino acid (aa) and nucleotide (nt) sequences from IMGT.
The constants TRAV, TRAJ, TRBV, TRBV and TRBJ are all types of notations of VDJ segments found in the database.
These entries are added to a translation dict and checked for direct matches in the IMGT fasta files.
Hereafter a series of mutations to the entries are done to see if the database entries will match the IMGT standard format.
The resulting translation dictionaries translate database V, D and J entries to aa and nt sequence.
"""

from collections import defaultdict
import os
from pysam import FastaFile
import re
from datetime import datetime
import json
import numpy as np
import logging
from collections import Counter

from TCR_reconstruction.constants import START_CODON, TRANSLATION_TABLE, AMBIGUOUS_AA_CODE, VDJ_CDR3_ALL_END_MOTIFS, \
    TRAJ_IMPUTATION_DICT, VDJ_NAME, TRAV, TRAJ, TRBV, TRBD, TRBJ

from TCR_reconstruction.logger.logger import init_logger

logger = logging.getLogger(__name__)


def reconstruct_full_tcr(v_column_nt, v_column_aa, j_column_nt, j_column_aa, cdr3_column_aa, include_leader=False,
                         exlude_C_F_W=False):
    # TODO: C and F/W checks do currently not allow to reconstruct the full seq
    # TODO: add allele inference of J (nt level) segment from CDR3 region (?)
    """
    Reconstructs the full TCR sequence from the V, J and CDR3 columns.

    Parameters
    ----------
    v_column_nt : list
        List of the V gene nucleotide sequences.
    v_column_aa : list
        List of the V gene amino acid sequences.
    j_column_nt : list
        List of the J gene nucleotide sequences.
    j_column_aa : list
        List of the J gene amino acid sequences.
    cdr3_column_aa : list
        List of the CDR3 amino acid sequences.
    include_leader : bool, optional
        If True, the leader sequence will be included in the full sequence.
        The default is True.

    Returns
    -------
    full_sequence_column : list
        List of the full TCR sequences.

    """
    full_sequence_column = []
    missing_c_fw_counter = Counter()
    # only if all columns contain a value
    for i, cdr3 in enumerate(cdr3_column_aa):
        try:
            if v_column_nt[i] is not None and j_column_nt[i] is not None and cdr3_column_aa[i] is not None and \
                    j_column_aa[i] is not None:
                if include_leader:
                    # find C104 and translate the leader+v
                    trxv_until_c_104 = find_c104(v_column_nt[i].upper(), nt_or_aa='nt')
                    trxv_aa = translate_nt_to_aa(trxv_until_c_104)
                else:
                    # the amino acid sequence is not including the leader
                    trxv_aa = find_c104(v_column_aa[i].upper(), nt_or_aa='aa')

                try:
                    assert cdr3[0] == 'C'
                except AssertionError:
                    # print('BE AWARE: the cdr3: {0} at index {1} did not start with a \'C\''.format(cdr3, i))
                    logger.debug(
                        'BE AWARE: the cdr3: {0} at index {1} did not start with a \'C\'. The full sequence will not be reconstructed!'.format(
                            cdr3, i))
                    missing_c_fw_counter.update(['no_C_at_start'])
                    if exlude_C_F_W:
                        full_sequence_column.append(np.nan)
                        continue

                try:
                    assert cdr3[-1] == 'F' or cdr3[-1] == 'W'
                except AssertionError:
                    # print('BE AWARE: the cdr3: {0} at index {1} did not end with a \'F\' or \'W\''.format(cdr3, i))
                    logger.debug(
                        'BE AWARE: the cdr3: {0} at index {1} did not end with a \'F\' or \'W\'. this is a warning, but the sequence will be reconstructed.'.format(
                            cdr3, i))
                    missing_c_fw_counter.update(['did_not_end_with_F_or_W'])
                    if exlude_C_F_W:
                        full_sequence_column.append(np.nan)
                        continue

                # add trbv and cdr3 together
                trxv_plus_cdr3 = trxv_aa + cdr3

                # Cannot use the script to match CDR3 end motifs in J segments, because the
                # "TRBJ2-4*01" and "TRBJ1-3*01" have multiple matches in multiple reading frames, both coding for
                # valid amino acid sequences. therefore the aa seq will be used to find the pattern
                trxj_aa = j_column_aa[i]
                index_j = find_cdr3_end_motif(trxj_aa, VDJ_CDR3_ALL_END_MOTIFS)
                full_seq = trxv_plus_cdr3 + trxj_aa[index_j + 1:]

                # full_seq = trxv_plus_cdr3
                full_sequence_column.append(full_seq)
            else:
                full_sequence_column.append(np.nan)


        except KeyError as key_error:
            logger.debug(
                'KEYERROR: {0}, at index: {1}. adding np.nan instead of full seq for cdr3 {2}'.format(key_error, i,
                                                                                                      cdr3))
            full_sequence_column.append(np.nan)
        except AttributeError as att_error:
            logger.debug('ATTRIBUTE_ERROR: {0}, at index: {1}. adding np.nan instead of full seq'.format(att_error, i))
            full_sequence_column.append(np.nan)
        except TypeError as type_error:
            logger.debug('TYPE_ERROR: {0}, at index: {1}. adding np.nan instead of full seq'.format(type_error, i))
            full_sequence_column.append(np.nan)
    logger.warning('did not reconstruct (reason:count) : {} '.format(missing_c_fw_counter))

    return full_sequence_column


def translate_nt_to_aa(dna):
    """
    Translates a DNA sequence into an amino acid sequence.

    Parameters
    ----------
    dna : str
        The DNA sequence to translate.

    Returns
    -------
    str
        The translated amino acid sequence.
    """
    protein = []
    end = len(dna) - (len(dna) % 3) - 1
    dna = dna.upper()

    for i in range(0, end, 3):
        codon = dna[i:i + 3]

        aminoacid = TRANSLATION_TABLE.get(codon, AMBIGUOUS_AA_CODE)
        protein.append(aminoacid)

    return "".join(protein)


def find_c104(v_seq, nt_or_aa):
    # TODO: This C is important for the cloning. Discuss what C to take & where in the script it is done. (nt and aa level)
    """
    Finds the first C in the V gene sequence.

    Parameters
    ----------
    v_seq : str
        The V gene sequence.
    nt_or_aa : str
        The sequence is in nucleotides or amino acids ('aa' or 'nt').

    Returns
    -------
    str
        The sequence until the first C in the V gene sequence.
    """

    if nt_or_aa == 'nt':
        # assert if the sequence starts with ATG,
        assert v_seq[0:3].upper() == START_CODON

        # the reading frame is determined to remove the last nts
        residual_nts = len(v_seq) % 3
        reading_frame_v = v_seq[:len(v_seq) - residual_nts]

        # find the first C going from the last codon to the left
        for codon_end in range(len(reading_frame_v), 0, -3):
            if TRANSLATION_TABLE[v_seq[codon_end - 3: codon_end]] == 'C':
                # return the sequence without the C, as this is the starting AA from the CDR3
                return v_seq[:codon_end - 3]

    elif nt_or_aa == 'aa':
        # find the first C going from the last codon to the left
        for index_c104 in range(len(v_seq) - 1, -1, -1):
            if v_seq[index_c104] == 'C':
                # return the sequence without the C, as this is the starting AA from the CDR3
                return v_seq[:index_c104]


def find_cdr3_end_motif(j_amino_acids, allowed_end_motifs):
    """
    Search the CDR3 end motif in the amino acid sequence of J

    Parameters
    ----------
    j_amino_acids : str
        The amino acids of the J gene.
    allowed_end_motifs : list
        The list of allowed end motifs.

    Returns
    -------
    int
        The index of the end of the CDR3 region.
    """

    for motif in allowed_end_motifs:
        assert motif in VDJ_CDR3_ALL_END_MOTIFS
        for idx in range(len(j_amino_acids) - len(motif) + 1):
            valid_end = True
            for i, aa in enumerate(motif):
                valid_end = valid_end and (j_amino_acids[idx + i] == aa or aa == 'X')
            if valid_end:
                if idx != 0:  # Fix to prevent the first hit to mess up the allignment  because of double end motifs in TRAJ35 and TRAJ16
                    return idx
                else:
                    continue
    return None


def find_start_codons(sequence_nt):
    """

    Finds the indices of all start codons in a given sequence.

    Parameters
    ----------
    sequence_nt : str
        The nucleotide sequence to be searched.

    Returns
    -------
    start_codon_list : list
        A list of the indices of all start codons in the sequence.

    Examples
    --------
    >>> find_start_codons('ATGTTATGAGGGTCAATGCTGCC')
    [0, 5, 15]
    """
    start_codon_list = []
    for idx in range(len(sequence_nt) - len(START_CODON) + 1):
        if sequence_nt[idx:idx + 3] == START_CODON:
            start_codon_list.append(idx)

    return start_codon_list


def reconstruct_vdj(vdj_column, vdj_type, translation_dict_path_nt, translation_dict_path_aa):
    """
    This function is used in the refinement of the database in the DataParser class (refine_dataset()).
    It matches the raw vdj allele names to the IMGT database using a translation
    dictionary (constructed in vdj_construction_utilt.py). This dictionary should be kept up to date if new data is added.
    :param vdj_column: column containing the v, d and j segments and v and j segments
    :param vdj_type: TRAV, TRAJ, TRBV, TRBD or TRBJ as a string
    :return: IMGT standardized name for aa and nt and
    """

    # RECONSTRUCTION FULL SEQ test script path
    nt_path = translation_dict_path_nt
    aa_path = translation_dict_path_aa

    with open(aa_path) as aa_dict:
        translation_dict_aa = json.load(aa_dict)
    with open(nt_path) as nt_dict:
        translation_dict_nt = json.load(nt_dict)

    parsed_vdj_column_aa = []
    parsed_vdj_column_nt = []
    full_seq_vdj_column_aa = []
    full_seq_vdj_column_nt = []

    for vdj_value in vdj_column[vdj_type]:
        try:
            parsed_vdj_column_aa.append(translation_dict_aa[vdj_type][vdj_value][0])
            full_seq_vdj_column_aa.append(translation_dict_aa[vdj_type][vdj_value][1])
            # logger.debug('MATCHED {} in the amino acid dictionary'.format(vdj_value))
        except (TypeError, KeyError):
            # logger.debug('could not match {} in the amino acid dictionary'.format(vdj_value))
            parsed_vdj_column_aa.append(np.nan)
            full_seq_vdj_column_aa.append(np.nan)

        try:
            parsed_vdj_column_nt.append(translation_dict_nt[vdj_type][vdj_value][0])
            full_seq_vdj_column_nt.append(translation_dict_nt[vdj_type][vdj_value][1])
            # logger.debug('MATCHED {} in the nucleotide dictionary'.format(vdj_value))
        except (TypeError, KeyError):
            # logger.debug('could not match {} in the nucleotide dictionary'.format(vdj_value))
            parsed_vdj_column_nt.append(np.nan)
            full_seq_vdj_column_nt.append(np.nan)

    return parsed_vdj_column_aa, full_seq_vdj_column_aa, parsed_vdj_column_nt, full_seq_vdj_column_nt


def IMGT_to_vdj_dict(path):
    """
    Initializes the dictionary used for translating to IMGT standard notation (uses downloaded fasta files).
    Returns a dictionary of IMGT-formatted FASTA sequences, where the key is the VDJ allele name (e.g. TRAV1-1*01)
    and the value is a nested dictionary of sequences, where the key is the sequence name (e.g. TRAV1-1*01|TRAJ9*01)
    and the value is the sequence.

    :param path: Path to the directory containing the IMGT-formatted FASTA files.
    :type path: str
    :return: A dictionary of IMGT-formatted FASTA sequences.
    :rtype: dict
    """
    vdj_seqs = defaultdict(dict)

    # loop over all files and fetch sequences for each vdj allele (aa or nt)
    for filename in os.listdir(path):
        if filename.endswith('.fasta'):
            vdj_name_1 = filename[0:4] + '_' + filename[-8:-6]  # 0:4 is the TRAV, -8:-6 is nt or aa
            vdj_seqs[vdj_name_1] = {}

            # add to nested dict (ie: dict[TRAJnt][TRAJ9*01]: sequence)
            with FastaFile(filename) as fasta_file:
                for ref in fasta_file.references:
                    vdj_seqs[vdj_name_1][ref.split('|')[1]] = fasta_file.fetch(
                        reference=ref)  # middle of ref is TRAJ9*01 etc

    # revert to normal dict, to avoid accidentally adding empty entries that are not from IMGT during imputing
    return dict(vdj_seqs)


def make_vdj_segment_translation_dict():
    # TODO: make able to input self. TRAJ etc fromt the parser.
    """
    This function creates a dictionary of dictionaries for the V, D, and J gene segment annotations found in the database.
    :return: A dictionary of dictionaries for the V, D, and J gene segments.
    """
    VDJ_DICT = defaultdict(dict)
    for i, VDJ in enumerate([TRAV, TRAJ, TRBV, TRBD, TRBJ]):
        VDJ_DICT[VDJ_NAME[i]] = defaultdict(list)
        for entry in VDJ:
            entry = entry.strip()  # remove the spaces around the entries
            VDJ_DICT[VDJ_NAME[i]][entry] = None

    return VDJ_DICT


def impute_missing(vdj_dict, aa_nt, imgt):
    # TODO: discuss with wouter and marius if this is okay (assume that you can add this): conclusion: compare the effect on the final result (accuracy)
    #  if entry+*01 exists but not entry+*02 --> add the *01 to the entry and match to IMGT
    # TODO: does it work without continue?
    """
    Impute the missing annotations for the V D and J.
    :param vdj_dict: annotations as keys and sequence as values
    :param aa_nt: amino acid or nucleotide dictionary
    :param imgt: imgt gold standard dict with annotations as keys and sequence as values
    :return: vdj translation dict with imputed annotations and corresponding sequences
    """
    # for TRAV, TRAJ, TRBV, TRBD, TRBJ:
    for vdj_ab in VDJ_NAME:
        logger.info(
            'Imputing missing keys for: {0}{1} #########################################\n\n'.format(vdj_ab, aa_nt))
        for entry in vdj_dict[vdj_ab]:

            # if not none: the entry is already present in the correct IMGT format
            if vdj_dict[vdj_ab][entry] is not None:
                logger.debug(' entry exists for: {}'.format(entry))
                continue

            else:
                if not re.match(vdj_ab, entry):
                    logger.debug(
                        '\tentry: {0} did not start with: {1}, adding None to entry in dict.'.format(entry, vdj_ab))
                    vdj_dict[vdj_ab][
                        entry] = None

                else:
                    suffix = entry[4:]

                    # if first 2 numbers match to a single key in IMGT, match them:
                    matches = [[key, val] for key, val in imgt[vdj_ab + aa_nt].items() if
                               key.startswith(vdj_ab + suffix[0:2])]

                    if len(matches) == 1:
                        # un-nest NESTED dict because of list comprehension
                        vdj_dict[vdj_ab][entry] = [matches[0][0], matches[0][1]]
                        logger.debug('\t{0} first 2 digits exactly 1 key: {1}'.format(entry, matches[0][0]))
                        continue

                    else:
                        # match entry + *01 and if match return key[:-3], which gives the key without *01 info.
                        with_01 = [key[:-3] for key, val in imgt[vdj_ab + aa_nt].items() if key == entry + '*01']
                        # with_02 = [key[:-3] for key, val in imgt[vdj_ab + aa_nt].items() if key == entry + '*02']
                        # hits_with_01_not_02 = [entry_with_01 for entry_with_01 in with_01 if entry_with_01 not in with_02]
                        if len(with_01) > 0:
                            fix = entry + '*01'
                            vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                            logger.debug('\t added *01 to {0} and found an exact match for: {1}'.format(entry, fix))
                            continue

                        # some entries have -0 when this can never be the case (ie should be TRAV1-01 istead of TRAV1-1)
                        if re.search('-0', entry):
                            fix = entry.replace('-0', '-')
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug('\t replaced -0 for - for entry: {0}, to {1}'.format(entry, fix))
                                continue
                            except KeyError:
                                try:
                                    fix = entry.replace('-', '*')
                                    vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                    logger.debug('\t replaced - for * for entry: {0}, to {1}'.format(entry, fix))
                                    continue
                                except KeyError:
                                    try:
                                        fix = entry.replace('0', '', 1) + '*01'
                                        vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                        logger.debug(
                                            '\t replaced -0 for - , and added *01 for entry: {0}, to {1}'.format(entry,
                                                                                                                 fix))
                                        continue
                                    except KeyError:
                                        try:
                                            fix = entry.replace('0', '', 2)
                                            vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                            logger.debug(
                                                'removed first and 2nd 0 for entry: {0}, fix: {1}'.format(entry, fix))
                                            continue
                                        except KeyError:
                                            logger.debug('Could not fix -0 for {}'.format(entry))
                                            continue

                        # an entry cannot begin with TRAV01, as it should be TRAV1.
                        if re.search('[TRAV]0\d', entry):
                            fix = entry.replace('0', '', 1)
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug('removed 0 for {0}, fix: {1}'.format(entry, fix))
                                continue
                            except KeyError:
                                try:
                                    fix = entry.replace('0', '', 1) + '*01'
                                    vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                    logger.debug('removed 0 and added *01 for {0}, fix: {1}'.format(entry, fix))
                                    continue
                                except KeyError:
                                    try:
                                        fix = entry.replace('0', '', 2)
                                        vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                        logger.debug('removed first and 2nd 0 {0}, fix: {1}'.format(entry, fix))
                                        continue
                                    except KeyError:
                                        try:
                                            fix = entry.replace('0', '', 2) + '*01'
                                            vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                            logger.debug(
                                                'removed first and 2nd 0 and added *01 {0}, fix: {1}'.format(entry,
                                                                                                             fix))
                                            continue
                                        except KeyError:
                                            logger.debug('adding zeros did not work for {0}'.format(entry))

                        if re.search(r":", suffix):  # entries with a ':' are changed to a *
                            suffix = suffix.replace(':', '*')
                            suffix = suffix.strip(',')
                            fix = vdj_ab + suffix
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug('\treplacing \':\' for \'*\' success for: {0} -> {1}'.format(entry, fix))
                                continue
                            except KeyError:
                                logger.debug('\tmissed exception replacing : for * {0} -> {1}'.format(entry, fix))

                        # find the /DV exceptions
                        if re.search('/DV', entry):
                            if re.search(' F', entry):
                                try:
                                    fix = entry.strip(' F') + '*01'
                                    vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                    logger.debug('\tstripping " F" successfull for: {0}, {1}'.format(entry, fix))
                                    continue
                                except KeyError as key_error:
                                    logger.debug('\tfailed stripping " F"  for {0}, {1}'.format(entry, key_error))

                        if re.search('DV', entry):
                            if re.search('\.', entry):
                                fixing = entry.strip('.')
                                logger.debug('\tthis entry: {} ,has a dot, stripping it now'.format(entry))
                            else:
                                fixing = entry

                            try:  # without (also adding) *01 had no successes
                                fix = fixing.replace('DV', '/DV') + '*01'
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug('\t replacing DV for /DV succsesfully for : {0}, {1}'.format(entry, fix))
                                continue
                            except KeyError as key_error:
                                logger.debug('failed replacing DV for /DV for: {0}, {1}'.format(entry, key_error))

                        # replace - for *
                        if bool(re.search(r"-", suffix)):
                            suffix = suffix.replace('-', '*')
                            fix = vdj_ab + suffix
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug('\treplacing - for * success for: {0} -> {1}'.format(entry, fix))
                                continue
                            except KeyError:

                                logger.debug('\treplacing - for * FAILED for: {0} -> {1}'.format(entry, fix))

                    if len(matches) > 1:
                        logger.debug('\t{} first 2 digits matches more than one key'.format(entry))
                        # continue
                    if len(matches) == 0:
                        logger.debug('\t{} first 2 digits matches no keys'.format(entry))
                        # continue

    return vdj_dict


def match_IMGT_to_translation_dict(IMGT, translation_dict, aa_nt):
    """
    Match all the annotations in the database to IMGT gold standard anntotations and extract corresponding sequences.
    This dictionary will be expanded on in the impute_missing() function)

    :param IMGT: IMGT dictionary with VDJ gene annotations as keys and sequence as value
    :param translation_dict: annotations from database as keys and empty values
    :param aa_nt: amino acid or nucleotide
    :return: translation dict, where exact matches to IMGT contain the corresponding sequence
    """
    logger.info(' Making {} translation dict: ##############################################\n\n'.format(aa_nt))
    for vdj in VDJ_NAME:

        logger.info('\t making translation dict for: {}'.format(vdj))
        logger.debug(
            '\tlength translation dict {0} before adding new keys: {1}'.format(vdj, len(translation_dict[vdj])))
        for entry in IMGT[vdj + aa_nt]:
            if entry != 'nan':
                try:
                    translation_dict[vdj][entry] = [entry, IMGT[vdj + aa_nt][entry]]

                except KeyError as key_error:  # should not be needed, because it uses a default dict!
                    logger.debug('KEYERROR: {} adding as a new entry to the translation_dict'.format(key_error))
                    translation_dict[vdj][entry] = None
        logger.debug(
            '\tlength translation dict {0} after adding new keys: {1}\n'.format(vdj, len(translation_dict[vdj])))

        # return to a normal dict for the next step
        translation_dict[vdj] = dict(translation_dict[vdj])

    return dict(translation_dict)


def impute_traj_column(cdr3a_column, traj_column):
    # TODO: This function is not used anymore, as the aa level is not informative for imputing J columns
    #   it might be possible for nt level, but was not yet tried.
    imputed_traj_list = []
    imputed_traj_aa_list = []
    for i, cdr3a in enumerate(cdr3a_column):
        try:
            traj_variant = traj_column[i][0:6]
            if traj_variant == 'TRAJ24':
                motif = cdr3a[-3:]
                if motif in TRAJ_IMPUTATION_DICT['TRAJ24'].keys():
                    imputed_traj, imputed_traj_aa = TRAJ_IMPUTATION_DICT['TRAJ24'][motif]  # return the correct
                    logger.debug('IMPUTED TRAJ24 : {0} for index {1}'.format(imputed_traj, i))
                    # print('IMPUTED {}'.format(imputed_traj))
                else:
                    imputed_traj, imputed_traj_aa = np.nan, np.nan
            elif traj_variant == 'TRAJ37':
                motif = cdr3a[-8:]
                if motif in TRAJ_IMPUTATION_DICT['TRAJ37'].keys():
                    imputed_traj, imputed_traj_aa = TRAJ_IMPUTATION_DICT['TRAJ37'][motif]  # return the correct
                    logger.debug('IMPUTED TRAJ37 : {0} for index {1}'.format(imputed_traj, i))
                    # print('IMPUTED {}'.format(imputed_traj))
                else:
                    imputed_traj, imputed_traj_aa = np.nan, np.nan
            else:
                imputed_traj, imputed_traj_aa = np.nan, np.nan

        # occurs when an empty traj is encountered or the annotation is not 6 characters long.
        except (KeyError, TypeError) as error_msg:
            logger.info(
                'ERROR: {0}, {1} could not impute the traj {2}, for index {3}'.format(error_msg, cdr3a, traj_column[i],
                                                                                      i))
            imputed_traj, imputed_traj_aa = np.nan, np.nan

        imputed_traj_list.append(imputed_traj)
        imputed_traj_aa_list.append(imputed_traj_aa)
    return imputed_traj_list, imputed_traj_aa_list


if __name__ == "__main__":
    """
    When running this script (instead of importing) it will create the translation dictionaries and impute missing annotations.
    """

    current_datetime = str(datetime.now())[0:13].replace(' ', '_') + 'h_'  # only take date and hour
    current_datetime = current_datetime.replace(':', '-')

    # # save imgt from fasta to dict in a json file

    # quick hack: need to be in the correct directory for the FastaFile function
    os.chdir('VDJ_gene_sequences/IMGT_versions/after_benchmark/functional_with_L-PART1+V-EXON_after_benchmark')
    # imgt = IMGT_to_vdj_dict(path='')
    # with open('imgt_to_py_dict_' + current_datetime + '.json', 'w', encoding='utf-8') as file_imgt:
    #     json.dump(imgt, file_imgt, ensure_ascii=False, indent=4)
    with open('imgt_to_py_dict_after_benchmark.json') as imgt_file:
        imgt = json.load(imgt_file)

    # make the aa translation dictionary
    logger = init_logger('vdj_parsing_aa.log')
    vdj_translation_dict_before_aa = make_vdj_segment_translation_dict()
    vdj_translation_dict_aa = match_IMGT_to_translation_dict(imgt, vdj_translation_dict_before_aa, '_aa')
    vdj_translation_dict_aa = impute_missing(vdj_translation_dict_aa, '_aa', imgt=imgt)
    with open('vdj_translation_dict_aa_' + current_datetime + '.json', 'w', encoding='utf-8') as vdj_transl:
        json.dump(vdj_translation_dict_aa, vdj_transl, ensure_ascii=False, indent=4)

    # make the nt translation dictionary
    logger = init_logger('vdj_parsing_nt.log')
    vdj_translation_dict_before_nt = make_vdj_segment_translation_dict()
    vdj_translation_dict_nt = match_IMGT_to_translation_dict(imgt, vdj_translation_dict_before_nt, '_nt')
    vdj_translation_dict_nt = impute_missing(vdj_translation_dict_nt, '_nt', imgt=imgt)
    with open('vdj_translation_dict_nt_' + current_datetime + '.json', 'w',
              encoding='utf-8') as vdj_transl:
        json.dump(vdj_translation_dict_nt, vdj_transl, ensure_ascii=False, indent=4)
