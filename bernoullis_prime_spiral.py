# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:15:20 2023

@author: Indranil

Plotting (r, theta) as (integer, integer) : polar co-ordinates
(Needs optimization in terms of checking for a prime)
"""

import matplotlib.pyplot as plt

# Un-comment for Bernoullis Spiral

x_lim = 1000

theta = [i for i in range(x_lim)]
r = [i for i in range(x_lim)]

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, r, 'k.')

ax.grid(True)
plt.show()


# # Un-comment for primes
# x_lim = 10000
# primes = []

# for i in range(2,x_lim):
#     prime_check = True
        
#     for j in range(2,i-1):
#         if i/j - i//j == 0:
#             prime_check = False
#             break
    
#     if prime_check:
#         primes.append(i)        

# theta = [primes[i] for i in range(len(primes))]
# r = [primes[i] for i in range(len(primes))]

# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# ax.plot(theta, r, 'k.')

# ax.grid(True)
# plt.show()