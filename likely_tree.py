from collections import Counter
from itertools import permutations
from ete3 import Tree
import math 
import ast
import pandas

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

def removeDuplicateTrees(treelst): 
    if len(treelst) == 1: 
        return treelst
    output = []
    visited = []
    for tree in treelst: 
        if tree not in visited: 
            sublst = [tree]
            tree1 = Tree(tree)
            visited.append(tree)  # Don't want to visit this tree anymore
            for subtree in treelst: 
                tree2 = Tree(subtree)
                rf, _, _, _, _, _, _ = tree1.robinson_foulds(tree2)
                if rf == 0 and subtree not in visited:  # Ensure we're not adding already visited trees
                    sublst.append(subtree)
                    visited.append(subtree)
            output.append(sublst)
    return [elem[0] for elem in output]

def generate_alt(sampleTree): #output = list containing alternative trees
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
    possible_ordering = permutations(duplicatedletters)
    for ordering in possible_ordering: 
        result.append(replacenodes(sampleTree, ordering))
    return result

def find_valid(alt_tree_lst, fileTree): 
    outlst = []
    for tree in alt_tree_lst: 
        new_tree = hide_integers([tree])
        original_tree = hide_integers(fileTree)
        if original_tree == new_tree: 
            outlst.append(tree)
    return outlst 

def isequal(tree1, tree2):
    if tree1[-1] != ";":
        tree1 = tree1 + ";"
    if tree2[-1] != ";": 
        tree2 = tree2 + ";"
    firstTree = Tree(tree1)
    secondTree = Tree(tree2)
    if len(tree1) == 2 or len(tree2) == 2: 
        if tree1 == tree2: 
            return True
        else:
            return False 
    else: 
        rf, _, _, _, _, _, _ = firstTree.robinson_foulds(secondTree)
        if rf == 0: 
            return True 
        else: 
            return False 

savedValues = {}

for i in range(5):
    try:
        curr = "rep_" + str(i) + ".tre" 
        file = open(curr, 'r').readlines()
        currentTree = Tree(file[0])
        treeHeight = currentTree.get_farthest_node()[1]
        interval_quantity = math.ceil(treeHeight/0.001)
        file2 = ast.literal_eval(open('myDictionary.txt', 'r').read())
        cleaned_tree = removeDuplicateTrees(find_valid(generate_alt(assignValues(hide_integers(file))), file))
        resultDictionary = {}
        for xtree in cleaned_tree:
            for ytree in file2[interval_quantity]:
                if isequal(xtree, ytree): 
                    resultDictionary[xtree] = file2[interval_quantity][ytree]
        print(file2[interval_quantity])
        print("True Tree: ", file[0]) 
        print("Predicted Tree: ",max(resultDictionary, key = resultDictionary.get))
        print("------------------------------------------------------")
        savedValues[file[0]] = max(resultDictionary, key = resultDictionary.get)
    except Exception as e: 
        print(e)

print(savedValues)