from numpy import *
from scipy.sparse import csr_matrix
import time
def elapsed_time(start_time):
        """Returns elapsed time given a start time"""
        now = time.clock()
        elapsed = now - start_time
        secs = float((int(10000*elapsed)))/10000
        mins = int(elapsed / 60)
        hours = int(mins / 60)
        mins = mins - hours * 60
        secs = secs - mins * 60 - hours * 3600
        if mins == 0:
            print('That took ' + str(secs) + ' seconds.\n')
        elif hours == 0:
            print('That took ' + str(mins) + ' minutes and ' +
                              str(secs) + ' seconds.\n')
        else:
            print('That took ' + str(hours) + ' hours ' +
                              str(mins) + ' minutes and ' + str(secs) +
                              ' seconds.\n')
        return now

test=100* random.rand(10000,1000)
print test
test[0:500]=nan
st=time.clock()
test2=csr_matrix(test)
st=elapsed_time(st)

# print test
# test [0]=nan
# test2=test-1
# print test2
duh2 = test2-test2
st=elapsed_time(st)
duh=test-test
st=elapsed_time(st)


st=time.clock()
duh = where(test2>0,test2,1)
st=elapsed_time(st)

st=time.clock()
duh = where(test>0,test,1)
st=elapsed_time(st)

st=time.clock()
duh2 = where(test2>0,test2,1)
st=elapsed_time(st)
