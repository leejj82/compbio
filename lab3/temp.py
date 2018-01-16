#!/usr/bin/env python
import sys
from datetime import datetime, date, time

#A=0,C=1,G=2,T=3
for i in range(0,10):
    for j in range(11,20):
        if j==15:
            print j
            break
        
student_tuples = [('john', 'A', 15), ('jane', 'B', 12),('dave', 'B', 10),]
print sorted(student_tuples, key=lambda student: student[2]) 