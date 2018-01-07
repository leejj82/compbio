#!/usr/bin/env python
import sys
from datetime import datetime, date, time

def func1():
    return

def f1(arr):
    for i in range(0,10):
        f2(i)
        print i
                

def f2(j): 
    if j==3:
        return
    else:
        f1(j)
    
f1(1)
                           
    