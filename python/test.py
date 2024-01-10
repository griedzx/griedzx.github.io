import cmath

z = 1 + 2j
sin_z = cmath.sin(z)
print(f"The sin of {z} is {sin_z}")

#极坐标转换表示
r,theta = cmath.polar(z)
print (f"The polar coordinates of {z} are (r={r}, theta={theta})")
print (f"The polar coordinates of {z} are ({r:.2f},{theta:.2f})")

from math import factorial
factorial(70)

import math
math.factorial(70)

a_list = [1,2,3,4,5]
alist = list(range(1,6))
blst = [1, 'a', 3.6, 2+5j]

alist + blst
blst * 3

a = {}
a = dict()
b = {'x':3, 'y':4}
c = dict(uid=105, login='Lumberjack', name='Michael Palin')

cset = set('anc')

# read(size)
# readline(size)
# readlines()

# for eachline in file:
#     print (eachline)

import sys

DNAstr = 'ATGCTGATCG'

A_num = C_num = G_num = T_num = 0

for i in range(len(DNAstr)):
    if DNAstr[i]=='A':
        A_num = A_num + 1
    elif DNAstr[i]=='C':
        C_num = C_num + 1
    elif DNAstr[i]=='G':
        G_num += 1 #another way
    elif DNAstr[i]=='T':
        T_num += 1 #another way
        
print("A_num =", A_num)
print("C_num =", C_num)
print("G_num =", G_num)
print("T_num =", T_num)

total_num = len(DNAstr)
print("A_frq =", A_num/total_num)
print("C_frq =", C_num/total_num)
print("G_frq =", G_num/total_num)
print("T_frq =", T_num/total_num)

dna_dic = {'A': 0, 'C': 0, 'G': 0, 'T': 0} #the dict for DNA numbers
for n in DNAstr:
    dna_dic[n] += 1
print(dna_dic) #print out the DNA numbers in a dict format

dna_frq_dic = {}
DNAs = 'ACGT' #the 4 kinds of nucleiotides
for d in DNAs: #calculate the frequencies
    dna_frq_dic[d] = dna_dic[d]/total_num
print(dna_frq_dic) #print out the DNA frequencies in a dict format

output_file_name = 'D:/frq.txt'
output_file = open(output_file_name, 'wt') #open file in a writing and text model
output_string = '\t'.join(DNAs)
output_file.write(output_string + '\n')
output_string = ''
for d in DNAs:
    output_string += str(dna_dic[d]) + '\t'
output_string = output_string.strip()
output_string += '\n'
output_file.write(output_string)
olst = [] #create a list for DNA frequencies
for d in DNAs:
    olst.append(str(dna_frq_dic[d]))
output_string = '\t'.join(olst) + '\n'
output_file.write(output_string)
output_file.close() #please check the file D:\frq.txt to see what you got

class Dog:
    def __init__(self):
        self.mouth = "big"
    def bark(self):
        print("Woof!")
        
class Dog:
    jaw = ["sharp", 32]
    paw = ["cute", 4]
    def __init__(self):
        self.mouth = "big"
        
print(Dog.jaw)
xiaobai = Dog()
print(xiaobai.paw)