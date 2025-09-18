
def read_fasta(path, name):
    data = []
    f_loc = path + '/' + name
    try:
        with open(f_loc, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if line[0] != '>':
                    data.append(line)
        # print(name, len(data), len(data[0]))
        return data
    except Exception as e:
        print(f"Error reading {f_loc}: {e}")
        return []



def read_data(path):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1 + short2 + long1