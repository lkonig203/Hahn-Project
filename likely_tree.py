from collections import Counter
from itertools import permutations

# Brings trees to simplest form 
def hide_integers(inlst):
    newlst = inlst[0].strip()
    output = ""
    i = 0 
    while i < len(newlst): 
        if (newlst[i] != '_') and not newlst[i].isdigit() and (newlst[i] != '.') and (newlst[i] != ':'):
            output = output + newlst[i]
        i = i + 1
    return output 

def assignValues(inputString):
    temp = {}
    output = ""
    for char in inputString: 
        if char.isalpha() and char not in temp.keys(): 
            temp[char] = 1
    for char in inputString: 
        if not char.isalpha(): # if it's not a letter
            output = output + char 
        else: # if it is a letter
            output = output + char # add letter 
            output = output + "_" + str(temp[char])
            temp[char] = temp[char] + 1
    return output

def grouped_permutations(lst):
    groups = {}
    for item in lst:
        letter = item.split('_')[0]
        groups.setdefault(letter, []).append(item)

    group_permutations = {key: list(permutations(value)) for key, value in groups.items()}

    group_keys = list(groups.keys())
    group_orders = permutations(group_keys)

    result = []
    for group_order in group_orders:
        for perm_combination in permutations([group_permutations[key] for key in group_order]):
            for combination in zip(*perm_combination):
                result.append([item for sublist in combination for item in sublist])
    
    return result

def replacenodes(inTree, newordering): 
    currindex = 0 
    i = 0
    letters = [item.split('_')[0] for item in newordering]
    output = ""
    while i < len(inTree): 
        if inTree[i].isalpha(): # if current index is a letter
            if inTree[i] in letters: 
                output = output + newordering[currindex]
                currindex = currindex + 1 # go to next item in newordering 
                i = i + 3
            else: 
                output = output + inTree[i]
                i = i + 1
        else: 
            output = output + inTree[i]
            i = i + 1
    return output 

def generate_alt(sampleTree): #output = list containing alternative trees
    print(sampleTree)
    letter_list = [char for char in sampleTree if char.isalpha()]
    letter_counts = Counter(letter_list)
    duplicatedletters = []
    result = []
    if len(letter_list) <= 2 or max(list(letter_counts.values())) < 2: 
        result.append(sampleTree)
        return result 
    for index, char in enumerate(sampleTree): 
        if char.isalpha(): 
            if letter_counts[char] > 1: #Letter appears more than once
                duplicatedletters.append(char + sampleTree[index+1] + sampleTree[index+2])
    possible_ordering = grouped_permutations(duplicatedletters)
    for ordering in possible_ordering: 
        result.append(replacenodes(sampleTree, ordering))
    return result

for i in range(1000):
    curr = "rep_" + str(i) + ".tre" 
    try:
        file = open(curr, 'r')
        cleaned_tree = generate_alt(assignValues(hide_integers(file.readlines())))
        print(cleaned_tree)
    except Exception as e:
        print(e)



