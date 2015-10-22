import os
import os.path
import subprocess
import sys

if __name__ == "__main__":
    dirname = sys.argv[1]
    for x in os.listdir(dirname):
        if x.endswith('.crt'):
            try:
                filename = os.path.join(dirname, x)
                filehash = subprocess.check_output(['openssl', 'x509', '-noout', '-hash', '-in', filename]).strip()
                filehash += '.0'
                hash_filename = os.path.join(dirname, filehash)
                if os.path.exists(hash_filename):
                    print x, filehash
                    os.remove(hash_filename)
                os.symlink(x, hash_filename)
            except:
                print "error in handling file:", filename
