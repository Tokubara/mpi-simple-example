main: main.c
	mpicc $^ -O0 -g -o $@ -D M=${M}
run: main
	mpirun -n 4 --oversubscribe ./main
clean:
	rm -rf ./main
debug:
	objdump -d -S -l main > objdump_log.txt
