


def parse_slurm_tasks_per_node(string):
    per_node = string.split(',')
    result = []
    for node in per_node:
        parts = node.split('(')
        count = parts[0]
        if len(parts) == 2:
            nodes = parts[1]
        else:
            nodes = None
            
        try:
            count = int(count)
        except:
            count = 0 # unparsable number
        if nodes:
            nodes = nodes[1:-1] # skip the 'x' character and the closing ')' character
            try:
                nodes = int(nodes)
            except:
                nodes = 1 # unparsable number, assume 1
            for _ in range(nodes):
                result.append(count)
        else:
            result.append(count)
    return result
    
def parse_slurm_nodelist(string):
    result = []
    
    name_characters = []
    position = 0
    while position < len(string):
        char = string[position]
        if char == '[':
            name = ''.join(name_characters)
            ids, position = parse_ids(string, position)
            for x in ids:
                result.append(name + x)
            name_characters = []
        elif char == ',':
            name = ''.join(name_characters)
            result.append(name)
            name_characters = []
            position += 1
        else:
            name_characters.append(char)
            position += 1
    if len(name_characters) > 0:
        name = ''.join(name_characters)
        result.append(name)
        name_characters = []
    return result
    
def parse_ids(string, position):
    result = []
    end = string.index(']',position)
    count_ranges = string[position+1:end]
    for count_range in count_ranges.split(','):
        if '-' in count_range:
            from_id, to_id = count_range.split('-')
            for number in range(int(from_id), int(to_id) + 1):
                result.append(str(number))
        else:
            result.append(count_range)
    return result, end+1
        
if __name__ == "__main__":
    print(parse_slurm_tasks_per_node("10(x4),3"))
    print(parse_slurm_nodelist("tcn[595,597-598,600-606],tcn100"))
