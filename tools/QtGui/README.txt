For building : 

* Install LAPACK and BLAS development packages.
* Install Qt 4.x SDK (from website or with you packages manager). 
* Depending of your OS or compiler, Run :
	'qmake-qt4 -spec linux-g++-64'
	'qmake-qt4 -spec linux-g++-32'
	'qmake-qt4 -spec linux-icc-64'
	'qmake-qt4 -spec win32-g++'
	... a lot of options are available
* Run 'make'.

The executable is './crystal_to_box_gui'

