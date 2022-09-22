from setuptools_scm import get_version
version = get_version()


def main():
    print(("%s" % version))


if __name__ == "__main__":
    main()
