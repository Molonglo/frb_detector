from distutils.core import setup
from distutils.extension import Extension
#from Cython.Build import cythonize
from Cython.Distutils import build_ext

"""
ext_modules = [
    Extension("name", ["name.pyx"]),
	]
	"""
ext_modules = [
		Extension("nsmd", ["nsmd.pyx"],
			extra_objects=["/home/dada/linux_64/lib/libmopsr.a"],
			include_dirs=["/home/dada/opt/sofa/include/",
				"/home/dada/linux_64/include/",
				"/home/dada/psrdada/src/"],
			libraries = ["mopsr","psrdada",
				"cpgplot","pgplot","gfortran",
				"X11","png","sofa_c"],
			library_dirs = ["/home/dada/psrdada/src/",
				"/usr/lib64",
				"/home/dada/linux_64/pgplot/",
				"/home/dada/opt/sofa/lib/"]
				)
		]
#/home/dada/linux_64/lib/
#extra_objects=["/home/dada/linux_64/lib/libmopsr.a"]

setup(
		name ="nsmd",
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
		)



#ext_modules = cythonize(ext_modules),



#$python setup.py build_ext --inplace
#--inplace is for the
