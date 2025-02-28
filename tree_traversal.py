from ete3 import PhyloTree
import math

def removeDuplicates(myDict):
    result = {}
    temp = list(myDict.keys())
    for i in range(len(temp) - 1): 
        if myDict[temp[i]] != 0: 
            result[temp[i]] = myDict[temp[i]]
    result[temp[-1]] = myDict[temp[-1]]
    return result

def flip_dict(dic_to_flip):
    currKeys = list(dic_to_flip.keys())
    currValues = list(dic_to_flip.values()) #list of lists
    temp = []
    output = {}
    for sublst in currValues:
        for val in sublst:
            temp.append(val)
    keyMax = max(temp) # Get key values for dictionary
    for i in range(keyMax + 1):
        store = []
        for entry in currKeys: 
            if i in dic_to_flip[entry]:
                store.append(entry)
        output[i] = store
    return output

def parseNewick(inTree):
    distances = {}
    for node in inTree.traverse("postorder"): 
        if node.is_leaf(): 
            distances[node.name] = node.dist
        else: 
            curr_lst = node.children
            running_str = "("
            for child in curr_lst:
                if child.name:  # Check if the child has a name
                    running_str += child.name + ","
            running_str = running_str.rstrip(",") + ")"  # Remove trailing comma
            node.name = running_str  
            distances[running_str] = node.dist
    return removeDuplicates(distances)

# Takes in dictionary with node/leaf and distances 
def break_interval(dist_dict): 
    interval = 0.2
    output = {}
    key_values = list(dist_dict.keys())
    for i in range(len(key_values)):
        if ("(" in key_values[i]): # This is an internal node or root
            previous_stop = output[key_values[i - 1]][-1] # access previous leaf/node from output 
            new_stop = previous_stop + (math.ceil(dist_dict[key_values[i]] / interval)) # Get new stop 
            output[key_values[i]] = [i for i in range(previous_stop, new_stop + 1)]
        else: # This is a leaf 
            output[key_values[i]] = [i for i in range(math.ceil(dist_dict[key_values[i]] / interval))]
    return output

# For Wednesday - merge dictionaries 
final = {}

# Doesn't return anything, just updates final
def mergeDict(input_dict, outputDict): 
    global final
    if len(outputDict) == 0: 
        final = input_dict
    else:
        for element in list(input_dict.keys()):
            if element in list(outputDict.keys()):
                mergedlist = outputDict[element] + input_dict[element]
                final.update({element:mergedlist})
            else:
                final.update({element:input_dict[element]})

for i in range(65000):
    print(i) 
    curr_str = "rep_" + str(i) + ".tre"
    try:
        curr_file = open(curr_str, 'r')
        this_tree = PhyloTree(curr_file.readline())
        mergeDict(flip_dict(break_interval(parseNewick(this_tree))), final)
    except: 
        continue

def count_repeated(in_list):
    return_dict = {}
    store_vals = []
    for topology in in_list:
        if topology not in store_vals: 
            store_vals.append(topology)
            curr_sum = in_list.count(topology)
            return_dict[topology] = curr_sum
    return return_dict

def countTopologies():
    global final 
    outdict = {}
    for key, value in final.items():
        outdict[key] = count_repeated(value)
    return outdict

def findProportion(topology_count_dict): # Takes in a dictionary 
    temp_dict = {}
    totalSum = sum(topology_count_dict.values())
    for key, value in topology_count_dict.items():
        temp_dict[key] = value
    return temp_dict

hold = countTopologies()
final_output = {}

for key, value in hold.items(): 
    final_output[key] = findProportion(value)

open("myDictionary.txt", 'w').write(str(final_output))
