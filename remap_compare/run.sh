icpc binary_cmp.cxx -o binary_cmp
/bin/rm -f binary_cmp.o
./binary_cmp weights.scrip weights.out 1 > weights.log
./binary_cmp initial_weights.scrip initial_weights.out 1 > initial_weights.log
./binary_cmp final_weights.scrip final_weights.out 1 > final_weights.log
./binary_cmp grid_info.scrip grid_info.out 2 > grid_info.log
diff -q grid_info.scrip grid_info.out
diff -q initial_weights.scrip initial_weights.out
diff -q final_weights.scrip final_weights.out
diff -q weights.scrip weights.out
