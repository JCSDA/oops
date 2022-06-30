#!/usr/bin/env python3
"""
@author: Benjamin Menetrier
@description: call generic cost function plotting
"""

def func(args):
    # Load module and run function
    module = __import__("cost")
    subfunc = getattr(module, "subfunc")
    subfunc(args)
