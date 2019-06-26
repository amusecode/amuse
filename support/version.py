major = 12
minor = 0
micro = 0
pre = ""
post = ""
version = "%s.%s.%s" % (major, minor, micro)
if pre != "":
    version += ".%s" % pre
if post != "":
    version += ".%s" % post

def main():
    print("%s" % version)

if __name__ == "__main__":
    main()
