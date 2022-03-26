# For Loop Sum
# For Loop sum
# FOR loop to 20 where every second index is added up to a double integer value. 
# Whats the result?
n = 20
sum_of_values = 0
for i in range(0,n,2):
    sum_of_values = sum_of_values+i
    
print(sum_of_values)        


#Average
#Average value of 10 array values (2,5,8,11,14,17,20,23,26,29)
val = [2,5,8,11,14,17,20,23,26,29]
print(sum(val)/len(val))


#Traffic Light
#Implement a state machine for a classic traffic light:

#5 seconds red
#2 seconds red and orange
#5 seconds green
#3 seconds green blinking
#2 seconds orange
# Start again
from enum import Enum
import time

class state(Enum):
    RED = 1
    GREEN = 2
    ORANGE = 3
    
traffic_light = state.RED

print(time.time())    
    
