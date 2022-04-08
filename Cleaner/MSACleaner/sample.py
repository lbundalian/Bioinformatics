# -*- coding: utf-8 -*-
"""

Java OOP vs Python OOP


"""


import sys
import getopt
import os


x = sys.argv[1]
y = sys.argv[2]


class Animal():
    
    age = 0
    name = None
    breed = None
    color = None

    # 5 properties    
    
    def __init__(self,_name,_age, _breed, _color):
        
        self.name = _name
        self.age = _age
        self.breed = _breed
        self.color = _color
        self.introduce()
        
    # 3 methods
    
    def introduce(self):
        
        
        
        print("My name is {0}, I am a {1}. I am {2} and my color is {3}".format(self.name,self.breed,self.age,self.color))
        

    def bark(self,times : int):
        
        for i in range(0, times):
            print("Aww")



if __name__ == '__main__':

    dog = Animal("sam",1,"beagle","brown")
    dog.bark(10)