import pprint, pickle
import pylab as plt

#with open('optical_sys.pkl', 'wb') as output:
#    pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)

with open('optical_sys.pkl', 'rb') as input:
    s = pickle.load(input)

pprint.pprint(s)

fig = plt.figure(1)
ax = fig.add_subplot(111)
s.draw2d(ax, color="red")

plt.show()
