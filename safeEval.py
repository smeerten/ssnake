import math
import re

#saf(ish) parsing of math statements
def safeEval(inp):
    env = vars(math).copy()
    env["locals"]   = None
    env["globals"]  = None
    env["__name__"] = None
    env["__file__"] = None
    env["__builtins__"] = None
    inp =  re.sub('([0-9]+)[k,K]','\g<1>*1024',inp) #WF: allow 'K' input
    try:
        return eval(inp,env)
    except:
        print("Could not interpret input")
        return 0
