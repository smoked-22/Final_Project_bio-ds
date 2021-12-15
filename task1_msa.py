# Copyright (c) 2021 Hyejin Lee (KAIST)

import sys
import datetime
import argparse

sys.setrecursionlimit(5000)


class ScoreParam():
    def __init__(self, template_seq, match=7, mismatch=-2, indel=-10):
        ############################################
        # implement function
        # Input: a list of strings (aligned sequences)
        ############################################
        self.match = match
        self.mismatch = mismatch
        self.indel = indel
        self.temp_seq = template_seq[1]

    def score(self, pos, nucleotide):
        ############################################
        # implement this function
        # Input: self, pos (int), nucleotide (str)
        # Output : a score for a specific nucleotide in a specific position
        ############################################
        if self.temp_seq[pos] == nucleotide:
            return self.match
        else:
            return self.mismatch


def greedy_align(template_seq, seq, score_dict):
    """
    aligns a sequence (str) to template_seq (list of str), which is MSA.
    same as two-way alignment, but has a difference in scoring function.
    """

    score = ScoreParam(template_seq, **score_dict)
    len_s1 = len(template_seq[1])
    len_s2 = len(seq)
    print("Len1 : {} Len2 : {}".format(len_s1, len_s2))
    scores = [[float("-inf") for _ in range(len_s2 + 1)] for _ in range(len_s1 + 1)]
    backtrack = [[None for _ in range(len_s2 + 1)] for _ in range(len_s1 + 1)]
    scores[0][0] = 0
    # total_len = (len_s1 + 1) * (len_s2 + 1)
    for i in range(len_s1 + 1):
        # large_loop = i * (len_s2 + 1)
        for j in range(len_s2 + 1):
            # print("Loop at = {}/{}\r".format(large_loop+j, total_len), end='')
            if i > 0 and j > 0:
                n_score = score.score(i - 1, seq[j - 1])
                if scores[i - 1][j - 1] + n_score > scores[i][j]:
                    scores[i][j] = scores[i - 1][j - 1] + n_score
                    backtrack[i][j] = "m"
            if i > 0 and scores[i - 1][j] + score.indel > scores[i][j]:
                scores[i][j] = scores[i - 1][j] + score.indel
                backtrack[i][j] = "d"
            if j > 0 and scores[i][j - 1] + score.indel > scores[i][j]:
                scores[i][j] = scores[i][j - 1] + score.indel
                backtrack[i][j] = "r"
    print("")
    return helper_backtrack(backtrack, template_seq, seq, len_s1, len_s2, [], "")


def helper_backtrack(backtrack_arr, template_seq, seq, i, j, template_alignment, seq2_alignment):
    """
    backtracking using the score matrix
    """
    if i == 0 and j == 0: return template_alignment + [seq2_alignment]
    if backtrack_arr[i][j] == "d":
        return helper_backtrack(backtrack_arr, template_seq, seq, i - 1, j,
                                add_aa(template_seq, template_alignment, i - 1), "-" + seq2_alignment)
    elif backtrack_arr[i][j] == "r":
        return helper_backtrack(backtrack_arr, template_seq, seq, i, j - 1,
                                add_aa(template_seq, template_alignment, "-"), seq[j - 1] + seq2_alignment)
    elif backtrack_arr[i][j] == "m":
        return helper_backtrack(backtrack_arr, template_seq, seq, i - 1, j - 1,
                                add_aa(template_seq, template_alignment, i - 1), seq[j - 1] + seq2_alignment)


def add_aa(template_seq, seq_list, pos):
    """
    helper function to modify a list of str (MSA; template_seq)
    """
    if len(seq_list) == 0:
        if pos == "-":
            return ["-" for seq in template_seq]
        return [seq[pos] for seq in template_seq]
    else:
        if pos == "-":
            return ["-" + seq for i, seq in enumerate(seq_list)]
        return [template_seq[i][pos] + seq for i, seq in enumerate(seq_list)]


def make_msa_from_file(filename):
    """
    input : fasta file
    output : list of sequences
    """
    return_lst = []
    msa_template = open(filename, "r")
    for line in msa_template:
        if '>' not in line:
            return_lst.append(line.strip())
    return return_lst


def save_list_as_file(lst, filename):
    f = open(filename, "w")
    for seq in lst:
        f.write(seq + '\n')
    f.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("match")
    parser.add_argument("mismatch")
    parser.add_argument("indel")
    args = parser.parse_args()
    score_dict = {"match": args.match, "mismatch": args.mismatch,
                  "indel": args.indel}
    print(">> Opening fasta files\n")
    print(datetime.datetime.now())

    template_msa = make_msa_from_file("spike_msa_template.fasta")
    files = ['HCOV19-ENGLAND-2021-04-19.fasta']\
        # ,\
        #      'HCOV19-ENGLAND-2021-05-03.fasta',\
    #          'HCOV19-ENGLAND-2021-05-17.fasta',\
    #          'HCOV19-ENGLAND-2021-05-31.fasta',\
    #          'HCOV19-ENGLAND-2021-06-14.fasta',\
    #          'HCOV19-ENGLAND-2021-06-28.fasta']
    lengths = [len(template_msa)]

    print(">> Running MSA algorithm\n")
    start_idx = lengths[0]
    for filenum, filename in enumerate(files):
        # Iteratively conduct sequence alignment to template MSA
        print(">> making MSA of file ", filename)
        print(datetime.datetime.now())

        seqlist = make_msa_from_file(filename)
        lengths.append(len(seqlist))

        for i, seq in enumerate(seqlist):
            print(">> sequence #{}".format(i))
            template_msa = greedy_align(template_msa, seq, score_dict)

        # save each MSA after alignment of each fasta files
        save_list_as_file(template_msa[start_idx: (start_idx + lengths[filenum + 1])], 'A-' + filename)
        start_idx += lengths[filenum]

    # re-save the final MSA to each file 
    # there is a possibility that the final MSA is different from the previous ones
    print(">> Saving to files\n")
    save_list_as_file(template_msa, "full_MSA.temp")
    start_idx = lengths[0]
    for filenum, filename in enumerate(files):
        save_list_as_file(template_msa[start_idx: (start_idx + lengths[filenum + 1])], 'A-' + filename)
        start_idx += lengths[filenum]


if __name__ == "__main__":
    main()
