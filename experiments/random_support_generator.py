import random
import os, sys
from math import *

#########################################
# functions

def inList(list, k):
	found = 0
	for l in list:
		if l == k:
			found = 1
			break
	return found

def generate_dense(d,mi,degree):
	supports = []
	random.seed()
	for i in range(d+1):
		support = []
		while len(support) < mi[i]:
			point=[]
                        p_sum=0
                        for di in range(d):
                                coord = random.randint(0, degree)
                                point.append(coord)
			        p_sum += coord
                        if p_sum <= degree and inList(support,point)==0:
			        support.append(point)
		supports.extend(support)
	return supports

def generate_sparse(d,mi,degree):
	supports = []
        degree = int(floor(degree / 2))
	#print degree
        random.seed()
	for i in range(d+1):
		support = []
		while len(support) < mi[i]:
			point=[]
                        p_sum=0
                        for di in range(d):
                                coord = random.randint(0, degree)
                                point.append(coord)
			        p_sum += coord
                        if inList(support,point)==0:
			        support.append(point)
		supports.extend(support)
	return supports

def generate_u_dense(d,mi,degree):
	supports = []
	#print degree
	random.seed()
	for i in range(d):
		support = []
		while len(support) < mi[i]:
			point=[]
			p_sum=0
			for di in range(d):
				coord = random.randint(0, degree)
				point.append(coord)
				p_sum += coord
			if p_sum <= degree and inList(support,point)==0:
				support.append(point)
		supports.extend(support)
	# construct u polynomial
	u_support=[]
	u=[0]*d
	u_support.append(u)
	for i in range(d):
		u=[0]*d
		u[i]=1
		u_support.append(u)
	supports.extend(u_support)
	return supports

def generate_u_sparse(d,mi,degree):
	supports = []
	degree = int(floor(degree / 2))
	#print degree
	random.seed()
	for i in range(d):
		support = []
		while len(support) < mi[i]:
			point=[]
			p_sum=0
			for di in range(d):
				coord = random.randint(0, degree)
				point.append(coord)
				p_sum += coord
			if inList(support,point)==0:
				support.append(point)
		supports.extend(support)
	# construct u polynomial
	u_support=[]
	u=[0]*d
	u_support.append(u)
	for i in range(d):
		u=[0]*d
		u[i]=1
		u_support.append(u)
	supports.extend(u_support)
	return supports

def write2file(d,mi,supports,count):
	sum_m=0
	for m in mi:
		sum_m+=m
	fres = open("inputs/input"+str(sum_m)+"_"+str(count), 'w')
	fres.write(str(d)+"\n")
	for m in mi:
		fres.write(str(m)+" ")
	fres.write("\n")
	fres.write(str(supports))
	fres.write("\n")

# ATTENTION 
# dimension = number of variables 
# !!!!!!!!!!!
d=6
# !!!!!!!!!!!

############################################
# 

supports=[]
# generate
for count in range(10):
        for n in range(12,30):
                mi=[n/(d+1)]*(d+1)
                for i in range(n%(d+1)):
                        mi[i]+=1
                #degree = int(floor(pow(10,log(n)/d)))
                degree = n
                print mi
                supports = generate_sparse(d,mi,degree)
                print supports
                write2file(d,mi,supports,count)


############################################
# u-resultants
'''
supports=[]
# generate
for count in range(10):
        for n in range(10,50):
                mi=[n/(d)]*(d)
                for i in range(n%(d)):
                        mi[i]+=1
                mi.append(d+1)
                degree = int(floor(pow(10,log(n)/d)))
                #degree = n
                #print degree
                supports = generate_u_sparse(d,mi,degree)
                #print supports
                write2file(d,mi,supports,count)
'''
