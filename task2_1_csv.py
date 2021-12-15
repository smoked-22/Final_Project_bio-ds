# Copyright (c) 2021 Hyejin Lee (KAIST)

def main():
    template_file = open("gene.fna", "r")
    template_file.readline()
    template_gene = ""
    for line in template_file:
        template_gene += line.strip()
    template_gene = "----" + template_gene + "ACGA"

    nt_list = ["A", "T", "G", "C", "-"]

    files = ['A-HCOV19-ENGLAND-2021-04-19.fasta']
    # ,\
    #         'A-HCOV19-ENGLAND-2021-05-03.fasta',\
    #          'A-HCOV19-ENGLAND-2021-05-17.fasta',\
    #          'A-HCOV19-ENGLAND-2021-05-31.fasta',\
    #          'A-HCOV19-ENGLAND-2021-06-14.fasta',\
    #          'A-HCOV19-ENGLAND-2021-06-28.fasta']

    mutation_counter = [{}]*5
    lengths = []
    for filename in files:
        date = filename[17:27]
        lst = make_list_from_file(filename)
        lengths.append(len(lst))
        pos_dict = make_pos_dict(lst)
        
        for pos in range(5,len(template_gene)-5):
            mutation_counter.append({})
            template_nt = template_gene[pos]
            for nt in nt_list:
                if nt != template_nt:
                    if nt not in mutation_counter[pos]:
                        mutation_counter[pos][nt] = []
                    mutation_counter[pos][nt].append(pos_dict[pos].count(nt))

    full_list = []
    for pos in range(5, len(template_gene)-6):
        counter = mutation_counter[pos]
        for nt in counter :
            if counter[nt] != [0,0,0,0,0,0]:
                full_list.append((pos, template_gene[pos], nt, counter[nt]))
    
    f = open("mutation_list.csv", "w")
    f.write("mutation,04-19,05-03,05-17,05-31,06-14,06-28\n")
    for pos, original, changed, counter in full_list:
        if changed == '-' :
            mutation = str(pos)+ "del"+ original
        else:
            mutation = str(pos) + original + ">" + changed
        f.write(f'{mutation},{counter[0]/lengths[0]*100},{counter[1]/lengths[1]*100},\
                {counter[2]/lengths[2]*100},{counter[3]/lengths[3]*100},{counter[4]/lengths[4]*100},{counter[5]/lengths[5]*100}\n')
    f.close()

    f = open("mutation_cumulative_list.csv", "w")
    f.write("mutation,04-19,05-03,05-17,05-31,06-14,06-28\n")
    ln = 0
    for pos, original, changed, counter in full_list:
        if changed == '-' :
            mutation = str(pos)+ "del"+ original
        else:
            mutation = str(pos) + original + ">" + changed
        f.write(f'{mutation},{counter[0]/lengths[0]*100},{(counter[0]+counter[1])/(lengths[0]+lengths[1])*100},\
                {(counter[0]+counter[1]+counter[2])/(lengths[0]+lengths[1]+lengths[2])*100},\
                {(counter[0]+counter[1]+counter[2]+counter[3])/(lengths[0]+lengths[1]+lengths[2]+lengths[3])*100},\
                {(counter[0]+counter[1]+counter[2]+counter[3]+counter[4])/(lengths[0]+lengths[1]+lengths[2]+lengths[3]+lengths[4])*100},\
                {(counter[0]+counter[1]+counter[2]+counter[3]+counter[4]+counter[5])/(lengths[0]+lengths[1]+lengths[2]+lengths[3]+lengths[4]+lengths[5])*100}\n')
    f.close()
            
        
def make_list_from_file(filename):
    f = open(filename, "r")
    lst = []
    for line in f:
        lst.append(line.strip())
    f.close()
    return lst

def make_pos_dict(seq_list):
    pos_dict = []
    pos_dict = [[seq[pos] for seq in seq_list] for pos in range(len(seq_list[0]))]
    return pos_dict

if __name__ == "__main__":
    main()
    
