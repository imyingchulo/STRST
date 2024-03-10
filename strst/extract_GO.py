
if __name__ == '__main__':

    file="/home/ubuntu/trinotate/Trinotate_0301.xls"

    dict = {}
    go_list=["id\tGO\tfunction"]
    with open(file, 'r') as fin:
        next(fin)
        for lines in fin.readlines():
            id = lines.strip().split('\t')[0]
            go = lines.strip().split("\t")[14]
            go_function = go.split("`")[0]
            go_only = go_function.split("^")[0]
            key = id
            if key not in dict:
                dict[key] = 0
            dict[key] = go

            # make list
            go_list.append("{}\t{}\t{}".format(id ,go_only, go_function.replace("^", ":")))

    output = "/home/ubuntu/onlyGO_0303.txt"

    with open(output, 'w') as fout:
        new_go_list = '\n'.join(go_list)
        fout.write(new_go_list)










