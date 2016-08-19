all:
	cd ./test/unit; make all_unit_tests

time_test:
	cd ./src; make full_time_test
