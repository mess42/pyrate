import pprint, pickle, pickletools
import pylab as plt
import jsonpickle

#with open('optical_sys.pkl', 'wb') as output:
#    pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)

with open('optical_sys.pkl', 'rb') as input:
    s = pickle.load(input)

s2 = jsonpickle.decode()

pprint.pprint(s)



fig = plt.figure(1)
ax = fig.add_subplot(111)
s2.draw2d(ax, color="red")

plt.show()
