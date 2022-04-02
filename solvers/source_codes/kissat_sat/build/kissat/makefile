all:
	$(MAKE) -C "/home/zhangxd/solvers/rdpaper/kissat_sat/build/kissat/build"
clean:
	rm -f "/home/zhangxd/solvers/rdpaper/kissat_sat/build/kissat"/makefile
	-$(MAKE) -C "/home/zhangxd/solvers/rdpaper/kissat_sat/build/kissat/build" clean
	rm -rf "/home/zhangxd/solvers/rdpaper/kissat_sat/build/kissat/build"
coverage:
	$(MAKE) -C "/home/zhangxd/solvers/rdpaper/kissat_sat/build/kissat/build" coverage
indent:
	$(MAKE) -C "/home/zhangxd/solvers/rdpaper/kissat_sat/build/kissat/build" indent
test:
	$(MAKE) -C "/home/zhangxd/solvers/rdpaper/kissat_sat/build/kissat/build" test
.PHONY: all clean coverage indent test
