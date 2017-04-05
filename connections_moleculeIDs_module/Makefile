
all: c

clean:
	rm -rf build nm_conn_molID*.so nm_conn_molID.cpp

c:
	rm -rf build nm_conn_molID*.so nm_conn_molID.cpp build
	python3 setup.py build_ext --inplace
	rm -rf build
