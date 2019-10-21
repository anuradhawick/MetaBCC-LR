import h5py
from multiprocessing import Pool
import sys
import os


dsk_output = sys.argv[1]
output = sys.argv[2]
threads = int(sys.argv[3])

input_file = h5py.File(dsk_output, "r")
keys = list(map(str, input_file["dsk/solid"].keys()))
input_file.close()

def write(key):
    try:
        cmd = """h5dump -d /dsk/solid/"""+ str(key) +""" -y --ddl """+ dsk_output +""" | awk '{if($1=="{" || $1 == "}," || $1=="}") { printf "" } else if (match($1, /,$/) || $1 == "") { printf $1}  else print $1}' > """ + output + "/" + str(key) +""".chunk"""
        os.system(cmd)
    except Exception as e:
        print("error", e)
    print("Finished chunk - " + str(key))


p = Pool(threads)

p.map(write, keys)
p.close()
p.join()