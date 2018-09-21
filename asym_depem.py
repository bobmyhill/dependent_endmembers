#!/usr/bin/python

from sympy.matrices import Matrix
from sympy import Symbol, nsimplify
import numpy as np

import sys

if ( len(sys.argv) != 2 ):
    print 'usage:', sys.argv[0], '<script name>'
    exit()

with open(sys.argv[1], 'r') as f:
    data = [line.split() for line in f if line[0] != '%']

mbrsa=data[0]
alpha=data[1]
mbrsb=data[2]


dalpha=Matrix([[0]*len(alpha)]*len(alpha))
valpha=Matrix([0]*len(alpha))
for i in range(len(alpha)):
    dalpha[i,i]=nsimplify(float(alpha[i]))
    valpha[i]=nsimplify(float(alpha[i]))

var=[]
for i in range(len(data[3])):
    var.append(Symbol(data[4][i]))

cvars=Matrix([data[4]])
c=Matrix(data[5:5+len(mbrsa)])
cp=Matrix(data[5+len(mbrsa):6+2*len(mbrsa)])

##############################################

nmbrs = len(mbrsa)
vGa=Matrix([['G({0})'.format(mbrsa[i]) for i in range(nmbrs)]]).T
vGb=Matrix([['p({0})'.format(mbrsa[i]) for i in range(nmbrs)]]).T
vpa=Matrix([['G({0})'.format(mbrsb[i]) for i in range(nmbrs)]]).T
vpb=Matrix([['p({0})'.format(mbrsb[i]) for i in range(nmbrs)]]).T
Wa=Matrix([['W({0},{1})'.format(mbrsa[i],mbrsa[j]) if i<j else 0 for j in range(nmbrs)]
           for i in range(nmbrs)])
Wb=Matrix([['W({0},{1})'.format(mbrsb[i],mbrsb[j]) if i<j else 0 for j in range(nmbrs)]
           for i in range(nmbrs)])


for i in range(len(mbrsa)):
    for j in range(len(mbrsa)):
        if (i<j):
            Wa[i,j]=nsimplify(2./(dalpha[i,i]+dalpha[j,j]))*Wa[i,j]
            
'''
print 'MINERAL PROPORTIONS (SET A)'
cva=c.col_join(-cvars).T

ns=cva.nullspace(simplified=True)
ns=transpose(ns)

for i in range(len(mbrsa)):
    print 'p(' + mbrsa[i] + ') =', ns[i][0]
print ''

print 'MINERAL PROPORTIONS (SET B)'
cvb=cp.col_join(-cvars).T
ns=cvb.nullspace(simplified=True)
ns=transpose(ns)
for i in range(len(mbrsa)):
    print 'p(' + mbrsa[i] + ') =', ns[i][0]
print 
'''
def symplify(f):
    return nsimplify(f.round(10), tolerance=1.e-9)

def simplify_matrix(arr):
    def f(i,j):
        return symplify(arr[i][j])
    return Matrix( len(arr), len(arr[0]), f )


A=simplify_matrix(np.linalg.lstsq(c.T, cp.T, rcond=None)[0])

alphap=Matrix([nsimplify(a) for a in A.T*valpha])

dalphap=Matrix([[0]*len(alpha)]*len(alpha))
dinvalphap=Matrix([[0]*len(alpha)]*len(alpha))
for i in range(len(alphap)):
    dalphap[i,i]=symplify(alphap[i])
    dinvalphap[i,i]=symplify(1./alphap[i])

pp=simplify_matrix(np.linalg.lstsq(cp.T, c.T, rcond=None)[0])*vpa

print('A -> B')
for i in range(len(mbrsb)):
    print('p({0}) = {1}'.format(mbrsb[i], pp[i]))
print('')

print('New alphas')
for i in range(len(mbrsb)):
    print('a({0}) = {1}'.format(mbrsb[i], alphap[i]))
print('')

B=dalpha*A*dinvalphap
Q=B.T*Wa*B

U=np.copy(Q)
L=np.copy(Q)
for i in range(len(mbrsb)):
    for j in range(len(mbrsb)):
        if (i>=j):
            U[i,j]=0
            L[j,i]=0

D=Matrix(np.diag(Q)).T

Q3=Matrix([[0]*len(D)]*len(D))
for i in range(len(D)):
    Q3[i,i]=D[i]

Qp3=Matrix([[0]*len(U)]*len(U)) # N.B Qp3=-Qf
for i in range(len(U)):
    if D[i]!=0:
        for j in range(len(U)):    
            if (j<i):
                Qp3[j,i]+=Qp3[j,i]+D[i]
            else:
                if (j>i):
                    Qp3[i,j]+=Qp3[i,j]+D[i]

Wbans=U + L.T - Qp3

print('W\'')
for i in range(len(U)):
    for j in range(len(U)):
        if (Wb[i,j] != 0):
            print('{0} = {1}'.format(Wb[i,j], nsimplify((alphap[i] + alphap[j])/2.)*Wbans[i,j]))
            
print('')

print('G\'')
vGbans=A.T*vGa + Q3*alphap
for i in range(len(vGb)):
    print('{0} = {1}'.format(vGb[i], vGbans[i]))
