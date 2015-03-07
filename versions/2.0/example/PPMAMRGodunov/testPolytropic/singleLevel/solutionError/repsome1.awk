#
# Print out the column headers given one line with the minimum and maximum
# resolution
#

# Look for line the one line
NR == 1 && NF == 2 {
  # Get the minimum nad maximum resolutions
  minres = $1;
  maxres = $2;

  # Skip the first column
  printf("%16s","");
  printf("  ");

  for (res = minres; 2*res <= maxres; res *= 2) {
    header=sprintf("%d (%d)",res,2*res);

    if (res == minres) {
      spaceLeft = 13 - length(header);
    } else {
      spaceLeft = 14 - length(header);
    }

    for (i = 0; i < spaceLeft / 2; i++) {
      printf(" ");
    }

    printf("%s",header);

    for ( ; i < spaceLeft; i++) {
      printf(" ");
    }

    if (2*res < maxres) {
      printf("  ");
    }
  }

  printf("\n");
}
