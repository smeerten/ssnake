import numpy as np

def euro(val, num):
    firstDigit = '%.0e' % val
    firstDigit = int(firstDigit[0])
    order = int(np.log10(val))
    if num < 1:
        return
    numStep = int(num)//3 + 1
    if firstDigit == 1:
        subset = [1,2,5]
    elif firstDigit == 2:
        subset = [2,5,10]
    elif firstDigit==5:
        subset = [5,10,20]
    else:
        return
    returnVal = np.tile(subset,numStep)
    orderArray = np.repeat(range(numStep),3)+order
    returnVal = returnVal*10.0**orderArray
    return returnVal[:num]
