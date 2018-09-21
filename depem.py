#!/usr/bin/python

from sympy.matrices import *
from sympy import Symbol, nsimplify

import numpy as np

import sys

if ( len(sys.argv) != 2 ):
    print 'usage:', sys.argv[0], '<script name>'
    exit()

with open(sys.argv[1], 'r') as f:
    data = [line.split() for line in f if line[0] != '%']

mbrsa=data[0]
mbrsb=data[1]

var=[]
for i in range(len(data[2])):
    var.append(Symbol(data[2][i]))

cvars=Matrix([data[3]])
c=Matrix(data[4:4+len(mbrsa)])
cp=Matrix(data[4+len(mbrsa):5+2*len(mbrsa)])


vGa=[]
vGb=[]
vpa=[]
Wa=[]
Wb=[]
for i in range(len(mbrsa)):
    vGa.append('G(' + mbrsa[i] + ')')
    vpa.append('p(' + mbrsa[i] + ')')
    vGb.append('G(' + mbrsb[i] + ')')
    Wtmpa=[]
    Wtmpb=[]
    for j in range(len(mbrsa)):
        if (i<j):
            Wtmpa.append('W(' + mbrsa[i] + ',' + mbrsa[j] + ')')
            Wtmpb.append('W(' + mbrsb[i] + ',' + mbrsb[j] + ')')
        else:
            Wtmpa.append(0)
            Wtmpb.append(0)
    Wa.append(Wtmpa)
    Wb.append(Wtmpb)

Wa=Matrix(Wa)

vGa=Matrix([vGa]).T
vGb=Matrix([vGb]).T
vpa=Matrix([vpa]).T

'''
print 'MINERAL PROPORTIONS (SET A)'
cva=c.col_join(-cvars).T

ns=cva.nullspace(simplified=True)
ns=np.transpose(ns)

for i in range(len(mbrsa)):
    print 'p(' + mbrsa[i] + ') =', ns[i][0]
print ''

print 'MINERAL PROPORTIONS (SET B)'
cvb=cp.col_join(-cvars).T
ns=cvb.nullspace(simplified=True)
ns=transpose(ns)
for i in range(len(mbrsa)):
    print 'p(' + mbrsa[i] + ') =', ns[i][0]
print ''
'''


def simplify_matrix(arr):
    def f(i,j):
        return nsimplify(arr[i][j], tolerance=1.e-9)
    return Matrix( len(arr), len(arr[0]), f )


A=simplify_matrix(np.linalg.lstsq(c.T, cp.T, rcond=None)[0].round(10))
pp=simplify_matrix(np.linalg.lstsq(cp.T, c.T, rcond=None)[0].round(10))*vpa

print 'A -> B'
for i in range(len(mbrsb)):
    print 'p(' + mbrsb[i] + ') =', pp[i]
print ''

Q=A.T*Wa*A

U=np.copy(Q)
L=np.copy(Q)
for i in range(len(mbrsb)):
    for j in range(len(mbrsb)):
        if (i>=j):
            U[i,j]=0
            L[j,i]=0

D=Matrix(np.diag(Q)).T

Qp3=Matrix([[0]*len(U)]*len(U))
for i in range(len(U)):
    if D[i]!=0:
        for j in range(len(U)):    
            if (j<i):
                Qp3[j,i]=Qp3[j,i]+D[i]
            else:
                if (j>i):
                    Qp3[i,j]=Qp3[i,j]+D[i]

Wbans=U + L.T - Qp3
print 'W\''
for i in range(len(U)):
    for j in range(len(U)):
        if (Wb[i][j] != 0):
            print Wb[i][j],  '=' , Wbans[i,j]
print ''

print 'G\''
vGbans=A.T*vGa +  D.T
for i in range(len(vGb)):
    print vGb[i],  '=' , vGbans[i]
