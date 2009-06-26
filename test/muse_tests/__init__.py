import sys,os

name_of_the_musedir_variable = 'MUSEDIR'

if not name_of_the_musedir_variable in os.environ:
	message = """
	No MUSEDIR environment variable found
	Please set MUSEDIR to the root directory of the muse distribution
	"""
	print message
	raise Exception(message)

sys.path.insert(0,os.environ['MUSEDIR'])

