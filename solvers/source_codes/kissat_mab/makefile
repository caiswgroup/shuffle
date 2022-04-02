all:
	$(MAKE) -C "/pub/data/zxd/shuffle/solvers/kissat_mab/build"
kissat:
	$(MAKE) -C "/pub/data/zxd/shuffle/solvers/kissat_mab/build" kissat
tissat:
	$(MAKE) -C "/pub/data/zxd/shuffle/solvers/kissat_mab/build" tissat
clean:
	rm -f "/pub/data/zxd/shuffle/solvers/kissat_mab"/makefile
	-$(MAKE) -C "/pub/data/zxd/shuffle/solvers/kissat_mab/build" clean
	rm -rf "/pub/data/zxd/shuffle/solvers/kissat_mab/build"
coverage:
	$(MAKE) -C "/pub/data/zxd/shuffle/solvers/kissat_mab/build" coverage
indent:
	$(MAKE) -C "/pub/data/zxd/shuffle/solvers/kissat_mab/build" indent
test:
	$(MAKE) -C "/pub/data/zxd/shuffle/solvers/kissat_mab/build" test
.PHONY: all clean coverage indent kissat test tissat
