#
# Output the component name followed by the error and convergence rate (if
# possible) for each resolution
#
# The input is assumed to be one line containing the component name followed
# by the error at each resolution (which all differ by a factor of 2).
#

# There should only be one line
NR == 1 {
  printf("%-16s  ",$1);
  for (i = 2; i <= NF; i++) {
    printf("%7.1e ",$i);
    if (i == 2) {
      printf("(n/a)");
    } else {
      if ($i == 0) {
        printf("(----)");
      } else {
        printf("(%4.2f)",log($(i-1)/$i)/log(2.0));
      }
    }
    if (i < NF) {
      printf("  ");
    }
  }
  printf("\n");
}
