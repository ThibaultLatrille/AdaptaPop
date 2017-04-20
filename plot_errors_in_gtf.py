import pickle as pickle
import matplotlib.pyplot as plt
# data_path = "/pandata/tlatrill/AdaptaPop/data"
data_path = "/mnt/sda1/AdaptaPop/data"

f = open("{0}/gtf_errors.p".format(data_path), "rb")
data_frame = pickle.load(f)
f.close()

mk_file = "79_mk_test.out"

RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"

my_dpi = 96
plt.figure(figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
print([annot for name, annot, seq, tot in data_frame])
annots = [(i, annot) for i, (name, annot, seq, tot) in enumerate(data_frame) if annot > 0]
plt.plot([i for i, annot in annots], [annot for i, annot in annots],
            color=YELLOW, lw=1)
plt.scatter([i for i, annot in annots], [annot for i, annot in annots],
            label="CSD with unknown start-end", c=YELLOW, s=70)
plt.plot(range(len(data_frame)), [seq for name, annot, seq, tot in data_frame],
            color=RED, lw=1)
plt.scatter(range(len(data_frame)), [seq for name, annot, seq, tot in data_frame],
            label="CDS with sequence lenght not multiple of 3", c=RED, s=70)
plt.plot(range(len(data_frame)), [tot for name, annot, seq, tot in data_frame],
            color=BLUE, lw=1)
plt.scatter(range(len(data_frame)), [tot for name, annot, seq, tot in data_frame],
            label="Total number of CDS", c=BLUE, s=70)
plt.legend()
plt.xlim(0, len(data_frame)-1)
plt.ylim((1, max([tot for name, annot, seq, tot in data_frame])))
plt.yscale('log')
plt.legend(loc=3, fontsize=18)
plt.savefig("errors_in_cds.svg", format="svg")
plt.show()


print('Test completed')

