
if __name__ == '__main__':

    # file="/home/ubuntu/trinotate/Trinotate_0301.xls"

    file = "/home/ubuntu/trinotate_E/Trinotate_0312.xls"

    dict = {}
    gene_list=["id\tGO\tfunction"]
    isoform_list=["id\tGO\tfunction"]
    with open(file, 'r') as fin:
        next(fin)
        visited = set()
        for lines in fin.readlines():
            gene_id = lines.strip().split('\t')[0]
            iso_id = lines.strip().split('\t')[1]
            go = lines.strip().split("\t")[14]
            go_function = go.split("`")[0]
            go_only = go_function.split("^")[0]
            key = gene_id
            if key not in dict:
                dict[key] = 0
            if key not in visited:
                dict[key] = go
                # make list
                gene_list.append("{}\t{}\t{}".format(gene_id, go_only, go_function.replace("^", ":")))
                isoform_list.append("{}\t{}\t{}".format(iso_id, go_only, go_function.replace("^", ":")))
                visited.add(key)


    output_gene = "/home/ubuntu/onlyGO_0312_gene.txt"
    output_isoform = "/home/ubuntu/onlyGO_0312_isoform.txt"

    with open(output_gene, 'w') as fout:
        new_go_list = '\n'.join(gene_list)
        fout.write(new_go_list)

    with open(output_isoform, 'w') as fout:
        new_go_list = '\n'.join(isoform_list)
        fout.write(new_go_list)










