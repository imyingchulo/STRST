


if __name__ == '__main__':

    number_dict={}
    data_dict={}
    sample = "/home/ubuntu/sample_input.txt"
    final_marker_list = []
    with open(sample, 'r') as fin:
        line = fin.readlines()
        visited =set()
        for i in line:
            raw_marker = i.strip().split("\t")[0]
            final_marker = i.strip().split("_")[1]

            if final_marker not in visited:
                final_marker_list.append(final_marker)
                visited.add(final_marker)

            key = final_marker
            # Check if key exists in number_dict, if not, initialize it with 0
            if key not in number_dict:
                number_dict[key] = 0

            number_dict[key] += 1

            if key not in data_dict:
                data_dict[key] = []

            data_dict[key].append("{}\t{}_rep{}\t{}".format(i.strip().split("\t")[0], i.strip().split("\t")[1], number_dict[key], i.strip().split("\t")[2]))

    raw =[]
    out = "/home/ubuntu/new_sample_list.txt"
    for i in final_marker_list:
        if i in data_dict:
            raw.append('\n'.join(data_dict[i]))

    clean_raw = '\n'.join(raw)
    with open(out, 'w') as fout:
        # print(clean_raw)
        fout.write(clean_raw)



